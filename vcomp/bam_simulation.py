import os
import random
import string
import subprocess

import pysam

import callers
import read_simulator as rs


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

def generate_reads(alt_genome_path, chr, pos, read_count=250, prefix="test-reads", read1_fh=None, read2_fh=None):
    """
    Generate reads in fastq format from the altered genome, return paths to the files generated
    :param alt_genome_path:
    :param read_count: Total number of read pairs to generate
    :param prefix: filename prefix for output files
    :return:
    """
    generator = rs.ReadSimulator(alt_genome_path, chr, pos)
    r1_filename = prefix + "_R1.fastq"
    r2_filename = prefix + "_R2.fastq"
    close = False
    if read1_fh is None:
        read1_fh = open(r1_filename, "w")
        close = True
    if read2_fh is None:
        read2_fh = open(r2_filename, "w")
        close = True

    for x in range(read_count):
        (a, b) = generator.gen_read_pair()
        read1_fh.write(a + "\n")
        read2_fh.write(b + "\n")
    if close:
        read1_fh.close()
        read2_fh.close()
    return (r1_filename, r2_filename)


def create_bam(ref_genome, reads1, reads2, bwapath, samtoolspath):
    dest = reads1.replace("_1.fq", "") + ".bam"
    cmd = bwapath + " mem " + " -R \'" + "\t".join(['@RG', 'ID:test', 'SM:sample', 'PL:Illumina']) + "\' " + ref_genome + " " + reads1 + " " + reads2 + " | " + samtoolspath + " sort -T sorttmp -O bam - > " + dest + "\n" + samtoolspath + " index " + dest + "\n"
    script_path = "./align.sh"
    with open(script_path, "w") as script_fh:
        script_fh.write(cmd)
    os.chmod(script_path, 0755)
    subprocess.check_call(script_path, shell=True)
    return dest

def gen_alt_fq(ref_path, variants, homs=True, depth=250):
    reads1 = "input_r1.fq"
    reads2 = "input_r2.fq"
    read1_fh = open(reads1, "w")
    read2_fh = open(reads2, "w")
    for variant in variants:
        alt_genome_path = 'alt_genome' + "".join([random.choice(string.ascii_lowercase + string.ascii_uppercase) for _ in range(10)]) + '.fa'
        alt_genome_size = gen_alt_genome(variant, ref_path, alt_genome_path, overwrite=True)

        if homs:
            generate_reads(alt_genome_path, variant.chrom, alt_genome_size/2,read_count=depth, read1_fh=read1_fh, read2_fh=read2_fh)
        else:
            generate_reads(ref_path, variant.chrom, variant.start, read_count=depth/2, read1_fh=read1_fh, read2_fh=read2_fh)
            generate_reads(alt_genome_path, variant.chrom, alt_genome_size/2, read_count=depth/2, read1_fh=read1_fh, read2_fh=read2_fh)

    read1_fh.close()
    read2_fh.close()
    return (reads1, reads2)

def gen_alt_bam(ref_path, conf, reads):
    reads1, reads2 = reads
    bam = create_bam(ref_path, reads1, reads2, conf.get('main', 'bwa_path'), conf.get('main', 'samtools_path'))
    verify_reads(reads1, reads2, bam, conf)
    return bam

def verify_reads(fq1, fq2, bam, conf):
    """
    Verify that all reads in the input fastq file are present in the bam file
    """
    r1 = len(list([line for line in open(fq1, "r") if line.strip()=='+']))
    r2 = len(list([line for line in open(fq2, "r") if line.strip()=='+']))
    cmd = conf.get('main', 'samtools_path') + " flagstat " + bam
    info = subprocess.check_output(cmd, shell=True, executable="/bin/bash")
    tot_line = info.split('\n')[0]
    bc = int(tot_line.split(' ')[0])
    if (r1 + r2) > bc:
        raise ValueError("BAM does not have same number of reads as input fastqs")

