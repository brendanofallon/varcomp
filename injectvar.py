__author__ = 'bofallon'

import pysam
import os
import subprocess
import json
import ConfigParser as cp
import argparse
import read_simulator as rs
import callers
import time
import comparators

def gen_alt_genome(variant, orig_genome_path, dest_filename, overwrite=False, window_size=2000):
    """
    Generate a new genome fasta that contains the given variant
    :param variant: Tuple of (chr, pos, ref, alt)
    :param orig_genome_path:
    :param dest_filename:
    :param overwrite: OK to overwrite existing dest_filename
    :param window_size: Number of ref genome bases to include
    :return: Length of genome created
    """
    if os.path.exists(dest_filename) and not overwrite:
        raise ValueError("Destination " + dest_filename + " exists and overwrite is set to False")

    mod_var_start = variant.start + 1
    ref = pysam.FastaFile(orig_genome_path)
    seq = ref.fetch(variant.chrom, mod_var_start-window_size/2, mod_var_start+window_size/2)
    alt = seq[0:window_size/2-1] + variant.alts[0] + seq[window_size/2+len(variant.ref)-1:]
    #print "Ref coords of injection window: " + str( variant[1]-window_size/2) + "-" + str( variant[1]+window_size/2)
    #print "Replacing " + seq[window_size/2-5:window_size/2+5] + " with " + alt[len(alt)/2-5:len(alt)/2+5]
    dest = open(dest_filename, "w")
    dest.write(">" + variant.chrom + "\n")
    dest.write(alt)
    dest.close()

    #Generate fai
    dest_index = open(dest_filename + ".fai", "w")
    dest_index.write(str(variant.chrom) + "\t" + str(len(alt)) + "\t" + str(len(variant.chrom)+1) + "\t" + str(len(alt)+1) + "\t" + str(len(alt)+1) + "\n")
    dest_index.close()
    return len(alt)

def generate_reads(alt_genome_path, chr, pos, mean_coverage=250, prefix="test-reads"):
    """
    Generate reads in fastq format from the altered genome, return paths to the files generated
    :param alt_genome_path:
    :param mean_coverage: Approx (very) mean coverage
    :param prefix: filename prefix for output files
    :return:
    """
    generator = rs.ReadSimulator(alt_genome_path, chr, pos)
    r1_filename = prefix + "_R1.fastq"
    r2_filename = prefix + "_R2.fastq"
    a_output = open(r1_filename, "w")
    b_output = open(r2_filename, "w")
    for x in range(mean_coverage):
        (a, b) = generator.gen_read_pair()
        a_output.write(a + "\n")
        b_output.write(b + "\n")
    a_output.close()
    b_output.close()
    return (r1_filename, r2_filename)


def create_bam(ref_genome, reads1, reads2, bwapath, samtoolspath):
    dest = reads1.replace("_1.fq", "") + ".bam"
    cmd = bwapath + " mem " + " -R \"" + "\\t".join(['@RG', 'ID:test', 'SM:sample', 'PL:Illumina']) + "\" " + ref_genome + " " + reads1 + " " + reads2 + " 2> /dev/null | " + samtoolspath + " sort -T sorttmp -O bam - > " + dest + "\n" + samtoolspath + " index " + dest + "\n"
    script_path = "align.sh"
    with open(script_path, "w") as script_fh:
        script_fh.write(cmd)
    os.chmod(script_path, 0755)
    subprocess.check_call(["bash", script_path])
    return dest


def write_vcf(variant, filename, conf, gt="1/1"):
    fh = open(filename, "w")
    fh.write("##fileformat=VCFv4.1\n")
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n')
    fh.write(variant.chrom + "\t")
    fh.write(str(variant.start+1) + "\t") #Remember - internally 0-based coords, but in vcf 1-based
    fh.write("." + "\t")
    fh.write(variant.ref + "\t")
    fh.write(variant.alts[0] + "\t")
    fh.write("100" + "\t")
    fh.write("PASS" + "\t")
    fh.write("." + "\t")
    fh.write("GT" + "\t")
    fh.write(gt + "\n")
    fh.close()
    return callers.compress_vcf(filename, conf)


def process_variant(variant, results, conf):
    """
    Process the given variant, update results dict
    :param variant:
    :param results_by_caller:
    :param results_by_method:
    :param conf:
    :return:
    """

    tmpdir = "tmp-working" + str(time.time())[-7:].replace(".", "")
    try:
        os.mkdir(tmpdir)
    except:
        pass
    os.chdir(tmpdir)

    if len(variant.alts) != 1:
        raise ValueError("Only one alt per variant now")
    orig_vcf = write_vcf(variant, "test_input.vcf", conf)
    ref_path = conf.get('main', 'ref_genome')

    alt_genome_path = 'alt_genome.fa'
    alt_genome_size = gen_alt_genome(variant, ref_path, alt_genome_path, overwrite=True)


    (reads1, reads2) = generate_reads(alt_genome_path, variant.chrom, alt_genome_size/2, prefix='inputvar')
    bam = create_bam(ref_path, reads1, reads2, conf.get('main', 'bwa_path'), conf.get('main', 'samtools_path'))

    variant_callers = callers.get_callers()
    variants = {}
    for caller in variant_callers:
        vars = variant_callers[caller](bam, ref_path, variant.chrom, variant.start-250, variant.start+250, conf)
        variants[caller] = vars

    for method_name, comp in comparators.get_comparators().iteritems():
        for caller, vars in variants.iteritems():
            result = comp(orig_vcf, vars, conf)
            print "Result for " + " ".join( str(variant).split()[0:5]) + ": " + caller + " / " + method_name + ": " + result
            results[caller][method_name][result] += 1

    os.chdir("..")
    os.system("rm -rf " + tmpdir)

def process_vcf(input_vcf, conf):
    """
    Iterate over entire vcf file, processing each variant individually and collecting results
    :param input_vcf:
    :param conf:
    :return:
    """
    input_vcf = pysam.VariantFile(input_vcf)

    #Initialize results structure
    all_results = {}
    for caller_name in callers.get_callers():
        all_results[caller_name] = {}
        for method_name in comparators.get_comparators():
            all_results[caller_name][method_name] = {
                comparators.NO_VARS_FOUND_RESULT: 0,
                comparators.NO_MATCH_RESULT: 0,
                comparators.MATCH_RESULT: 0,
                comparators.PARTIAL_MATCH: 0
            }


    for input_var in input_vcf:
        try:
            process_variant(input_var, all_results, conf)
        except Exception as ex:
            print "Error processing variant " + str(input_var) + " : " + str(ex)

    for caller in all_results:
        print "Caller: " + caller
        for method in all_results[caller]:
            print "\t" + method
            for result, count in all_results[caller][method].iteritems():
                print "\t\t" + result + "\t:\t" + str(count)


if __name__=="__main__":
    parser = argparse.ArgumentParser("Inject, simulate, call, compare")
    parser.add_argument("-c", "--conf", help="Path to configuration file", default="./comp.conf")
    parser.add_argument("-v", "--vcf", help="Input vcf file")
    args = parser.parse_args()

    conf = cp.SafeConfigParser()
    conf.read(args.conf)

    process_vcf(args.vcf, conf)
