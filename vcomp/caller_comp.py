import ConfigParser as cp
import argparse
import os
import random
import string

import pysam

import injectvar, bam_simulation, callers, comparators


def process_variant(variant_batch, results, conf, batchnum, homs):
    """
    Process the given variant, update results dict
    :param variant:
    :param results_by_caller:
    :param results_by_method:
    :param conf:
    :return:
    """

    tmpdir = "tmp-varcomp-" + str(batchnum) + "-" + "".join([random.choice(string.ascii_uppercase + string.ascii_lowercase) for _ in range(8)])
    try:
        os.mkdir(tmpdir)
    except:
        pass
    os.chdir(tmpdir)

    ref_path = conf.get('main', 'ref_genome')

    bed = callers.vars_to_bed(variant_batch)
    bam = bam_simulation.gen_alt_bam(ref_path, variant_batch, conf, homs)

    variant_callers = callers.get_callers()
    variants = {}
    for caller in variant_callers:
        vars = variant_callers[caller](bam, ref_path, bed, conf)
        variants[caller] = vars

    for variant in variant_batch:
        orig_vcf = bam_simulation.write_vcf(variant, "test_input.vcf", conf)
        bed = callers.vars_to_bed([variant])
        vgraph_comp = comparators.get_comparators()['vgraph']
        for caller, vars in variants.iteritems():

            result = comparators.compare_genotype(orig_vcf, vars, bed)
            if result is None:
                result = vgraph_comp(orig_vcf, vars, conf, bed)
            print "Result for " + " ".join( str(variant).split()[0:5]) + ": " + caller + ": " + result
            results[caller][result] += 1

    os.chdir("..")
    os.system("rm -rf " + tmpdir)

def var_sort(a, b):
    if a.chrom == b.chrom:
        return a.start - b.start
    if a.chrom > b.chrom:
        return -1
    else:
        return 1

def process_vcf(input_vcf, homs, conf):
    """
    Iterate over entire vcf file, processing each variant individually and collecting results
    :param input_vcf:
    :param conf:
    :return:
    """


    batches = injectvar.batch_variants(pysam.VariantFile(input_vcf))
    assert sum([len(b) for b in batches]) == len(list(pysam.VariantFile(input_vcf)))
    #for batch in batches:
    #    print "\n Batch of size " + str(len(batch))
    #    for v in sorted(batch, cmp=var_sort):
    #        print "  " + v.chrom + "\t" + str(v.start) + "\t" + str(v.ref) + "\t" + str(v.alts[0])
    for num, batch in enumerate(batches):
        injectvar.process_batch(sorted(batch, cmp=var_sort),  str(num), conf, homs, disable_flagging=True, read_depth=40)




if __name__=="__main__":
    parser = argparse.ArgumentParser("Test variant callers")
    parser.add_argument("-c", "--conf", help="Path to configuration file", default="./comp.conf")
    parser.add_argument("-v", "--vcf", help="Input vcf file")
    parser.add_argument("--het", help="Run all variants as hets (default false, run everything as homs)", action='store_true')
    args = parser.parse_args()

    conf = cp.SafeConfigParser()
    conf.read(args.conf)

    process_vcf(args.vcf, not args.het, conf)
