import os
import subprocess
from collections import defaultdict

import pysam

import read_simulator as rs
import vcomp.util as util

ALL_HETS="all hets"
CIS = "cis"
TRANS = "trans"
ALL_HOMS="all homs"

def gen_alt_genome(chrom, pvars, orig_genome_path, dest_filename, overwrite=False, window_size=2000):
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

    pvars = sorted(pvars, key=lambda x: x[0], reverse=True)
    ref_genome = pysam.FastaFile(orig_genome_path)
    mod_var_start = pvars[0][0] + 1 #Start position of first variant
    window_start = mod_var_start-window_size/2
    window_end = mod_var_start+window_size/2
    seq = ref_genome.fetch(chrom, window_start, window_end)
    newseq = seq

    for start, ref, alt in pvars:
        dstart = start - window_start
        newseq = newseq[0:dstart] + alt + newseq[dstart+len(ref):]

    dest = open(dest_filename, "w")
    dest.write(">" + chrom + "\n")
    dest.write(newseq)
    dest.close()

    #Generate fai
    dest_index = open(dest_filename + ".fai", "w")
    dest_index.write(str(chrom) + "\t" + str(len(newseq)) + "\t" + str(len(chrom)+1) + "\t" + str(len(newseq)+1) + "\t" + str(len(newseq)+1) + "\n")
    dest_index.close()
    return len(newseq)

def generate_reads(alt_genome_path, chr, pos, read_count=250, prefix="test-reads", read1_fh=None, read2_fh=None):
    """
    Generate reads in fastq format from the altered genome, return paths to the files generated
    :param alt_genome_path:
    :param read_count: Total number of read pairs to generate
    :param prefix: filename prefix for output files
    :return: Paths to two fastq files containing reads
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
    cmd = bwapath + " mem " + " -I 250.0,50,500 -R \'" + "\t".join(['@RG', 'ID:test', 'SM:sample', 'PL:Illumina']) + "\' " + ref_genome + " " + reads1 + " " + reads2 + " | " + samtoolspath + " sort -T sorttmp -O bam - > " + dest + "\n" + samtoolspath + " index " + dest + "\n"
    script_path = "./align.sh"
    with open(script_path, "w") as script_fh:
        script_fh.write(cmd)
    os.chmod(script_path, 0755)
    subprocess.check_call(script_path, shell=True)
    return dest

def gen_bam_stats(bamfile, region=None):
    """
    Calculate various metrics, like total reads, reads with MQ > X, perhaps some measure of reference bias, etc.
    :param bamfile:
    :return:
    """
    stats = defaultdict(int)
    af = pysam.AlignmentFile(bamfile)
    for read in af.fetch(region[0], region[1], region[2]):
        stats['total_reads'] += 1
        if read.is_proper_pair:
            stats['properpair'] += 1
        if read.mapping_quality > 20:
            stats['mq20'] += 1
        if read.mapping_quality > 40:
            stats['mq40'] += 1
        softclipped = 0
        try:
            for item in read.cigartuples:
                if item[0] == 4:
                    softclipped += item[1]
        except:
            pass
        stats['softclipped_bases'] += softclipped
        if softclipped > 0:
            stats['softclipped_reads'] += 1
    return stats

def collect_alts(vset):
    """
    Given a list of variants and a 'policy' describing how to arrange them, construct two haplotypes representing
    phased versions of the variants. Each 'haplotype' is a list of a start position, reference allele, and a
    single alt allele (these lists are then used to construct a fasta genome from which we can simulate reads).
    Note that CIS and TRANS are not very well defined if there are multi-alt variants in the vars
    :param vset:
    :return: Tuple of two 'haplotypes'
    """
    hap1 = []
    hap2 = []
    policy = vset['policy']
    for i,var in enumerate(vset['vars']):
        if policy == ALL_HOMS:
            if len(var.alts)>1:
                raise ValueError("Can't use ALL_HOMS policy for variants with multiple alts")
            hap1.append( (var.start, var.ref, var.alts[0]) )
            hap2.append( (var.start, var.ref, var.alts[0]) )
        elif policy == TRANS:
            first = i%2==0
            for alt in var.alts:
                if first:
                    hap1.append( (var.start, var.ref, alt) )
                    hap2.append( (var.start, var.ref, var.ref) ) #Required to add non-alt allele so we can still generate an alt genome
                else:
                    hap2.append( (var.start, var.ref, alt) )
                    hap1.append( (var.start, var.ref, var.ref) )

                first = not first
        elif policy == CIS:
            first = True
            for alt in var.alts:
                if first:
                    hap1.append( (var.start, var.ref, alt) )
                    hap2.append( (var.start, var.ref, var.ref) )
                else:
                    hap2.append( (var.start, var.ref, alt) )
                    hap1.append( (var.start, var.ref, var.ref) )
                first = not first
        elif policy == ALL_HETS:
            if len(var.alts)==1:
                hap1.append( (var.start, var.ref, var.alts[0]) )
                hap2.append( (var.start, var.ref, var.ref) )
            elif len(var.alts)==2:
                hap1.append( (var.start, var.ref, var.alts[0]) )
                hap2.append( (var.start, var.ref, var.alts[1]) )
            else:
                raise ValueError('Cant handle triploid / polyploid variants')
        else:
            raise ValueError('Unrecognized phasing policy: ' + policy)

    return hap1, hap2

def gen_alt_fq(ref_path, variant_sets, read_count, dest_prefix="input"):
    """
    Generate a batch of simulated reads independently for the variants in each variant_set
    Each set contains a list of variants and a policy describing cis / trans configuration
    :param ref_path: Path to reference fasta
    :param variant_sets: List of sets of variants to inject into reference
    :param read_count:
    :return:
    """
    reads1 = dest_prefix + "_R1.fq"
    reads2 = dest_prefix + "_R2.fq"
    read1_fh = open(reads1, "w")
    read2_fh = open(reads2, "w")
    read_count = int(read_count * 1.25)  # Make sure depth in IGV is about whatever read_count is
    for vset in variant_sets:
        chrom = vset['vars'][0].chrom
        hap1, hap2 = collect_alts(vset)

        alt_genome_path = 'alt_genome' + util.randstr() + '.fa'
        alt_genome_size = gen_alt_genome(chrom, hap1, ref_path, alt_genome_path, overwrite=True)
        generate_reads(alt_genome_path, chrom, alt_genome_size / 2, read_count=read_count / 2, read1_fh=read1_fh, read2_fh=read2_fh)
        os.remove(alt_genome_path)
        os.remove(alt_genome_path + ".fai")

        alt_genome_path = 'alt_genome' + util.randstr() + '.fa'
        alt_genome_size = gen_alt_genome(chrom, hap2, ref_path, alt_genome_path, overwrite=True)
        generate_reads(alt_genome_path, chrom, alt_genome_size / 2, read_count=read_count / 2, read1_fh=read1_fh, read2_fh=read2_fh)
        os.remove(alt_genome_path)
        os.remove(alt_genome_path + ".fai")

    read1_fh.close()
    read2_fh.close()
    return (reads1, reads2)

def gen_alt_bam(ref_path, conf, reads):
    """
    Align reads to reference, sort them, and generate an indexed .bam file. This assumes
    BWA, but we should allow this to be defined in a configuration.
    :param ref_path: Path to reference genome
    :param conf: Configuration containing paths to BWA, samtools, etc
    :param reads: Paths to reads to align (assumes paired-end)
    :return: Path to bam file
    """
    #TODO: Allow different alignment tools
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

