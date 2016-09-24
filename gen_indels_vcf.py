import sys
import csv
import argparse
import random

from collections import defaultdict

import pysam

from kevinlib.structs import interval_tree

bases = ['A', 'C', 'T', 'G']
region_margin = 20


def read_regions(filename):
    for line in open(filename):
        toks  = line.split()
        chrom = toks[0]
        start = int(toks[1])+region_margin
        stop  = int(toks[2])-region_margin
        if stop > start:
            yield chrom, start, stop


def gen_insertion(ref, chrom, location, size):
    insertion = "".join(random.choice(bases) for _ in range(size))
    ref_base = ref.fetch(chrom, location, location+1)
    return chrom, str(location+1), ".", ref_base, ref_base + insertion


def gen_deletion(ref, chrom, location, size):
    ref_bases = ref.fetch(chrom, location, location+size+1)
    return chrom, str(location+1), ".", ref_bases, ref_bases[0]


def gen_duplication(ref, chrom, location, size):
    ref_bases = ref.fetch(chrom, location, location+size)
    return chrom, str(location+1), ".", ref_bases, ref_bases + ref_bases


def gen_inverse_dup(ref, chrom, location, size):
    # FIXME: The interval arithmetic looks wrong.
    #        fetch should take half open intervals, so the +1 looks spurious.
    ref_bases = ref.fetch(chrom, location, location+size+1)
    return chrom, str(location+1), ".", ref_bases, ref_bases + revcomp(ref_bases)


def gen_snp(ref, chrom, location, size):
    ref_bases = ref.fetch(chrom, location, location+1)
    alt = ref_bases
    while alt == ref_bases:
        alt = random.choice(bases)
    return chrom, str(location+1), ".", ref_bases, alt


def gen_blocksub(ref, chrom, location, size):
    """
    Length-preserving variant
    :param ref:
    :param chrom:
    :param location:
    :param size:
    :return:
    """
    # FIXME: Each alt base may not match the corresponding ref base.  This
    # makes easier variants than if inserting truly random sequence, since
    # the alleles are position-wise distinct.  Perhaps make sure that only
    # the first and last base are distinct from ref?
    ref_bases = ref.fetch(chrom, location, location+size)
    alt = ""
    for b in ref_bases:
        a = b
        while a==b:
            a = random.choice(bases)
        alt += a
    return chrom, str(location+1), ".", ref_bases, alt


def gen_ins_mnp(ref, chrom, location, size):
    insertion = "".join(random.choice(bases) for _ in range(size))
    ref_base = ref.fetch(chrom, location, location+1)
    alt = ref_base
    while alt == ref_base:
        alt = random.choice(bases)
    return chrom, str(location+1), ".", ref_base, ref_base + insertion + "," + ref_base + alt


def gen_del_snp_mnp(ref, chrom, location, size):
    ref_bases = ref.fetch(chrom, location, location+size+1)
    del_alt = ref_bases[0]
    okbases = ['A', 'C', 'T', 'G']
    okbases.remove(ref_bases[1])
    alt = ref_bases[0] + random.choice(okbases) + ref_bases[2:]
    return chrom, str(location+1), ".", ref_bases, alt + "," + del_alt


def pick_location(regions, blacklist, min_safe_dist=2000):
    chrom, start, stop = random.choice(regions)
    loc = random.randint(start, stop)

    while (loc-min_safe_dist, loc+min_safe_dist) in blacklist[chrom]:
        chrom, start, stop = random.choice(regions)
        loc = random.randint(start, stop)
    blacklist[chrom].insert(loc, loc+1)

    return chrom, loc


def generate(ref, regions, output):
    reps_per_size = 10
    blacklist = defaultdict(interval_tree.IntervalTree)
    output.write("##fileformat=VCFv4.1\n")
    output.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    output.flush()
    out = csv.writer(output, delimiter='\t')
    out.writerow(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample'])

    gen_vars = [gen_deletion, gen_insertion, gen_duplication, gen_blocksub]
    # gen_del_snp_mnp, gen_snp, gen_inverse_dup

    for rep in range(0, reps_per_size):
        for size in [1,2,3,4,5] + range(10, 150, 10):
            for gen_var in gen_vars:
                loc = pick_location(regions, blacklist)
                var = gen_var(ref, loc[0], loc[1], size)
                out.writerow(var + ('.', '.', '.'))


def revcomp(bases):
    """
    Return reverse-complemented bases, uses the revcomp_lookup global lookup dictionary
    :return: Reverse-complemented bases as a string
    """
    revcomp_lookup={'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G' }
    return ''.join(revcomp_lookup[b] for b in bases[::-1])


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Create a VCF file with many simulated indels")
    parser.add_argument("-f", "--fasta", help="Reference fasta", required=True)
    parser.add_argument("-r", "--regions", help="Input BED file with target regions", required=True)
    args = parser.parse_args()

    regions = list(read_regions(args.regions))
    ref = pysam.FastaFile(args.fasta)

    generate(ref, regions, sys.stdout)
