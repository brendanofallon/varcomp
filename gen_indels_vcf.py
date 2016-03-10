
import argparse
import intervaltree
import pysam
import random
import sys
from collections import defaultdict

bases = ['A', 'C', 'T', 'G']
region_margin = 20

def to_region(region_str):
    toks = region_str.split()
    start = int(toks[1])+region_margin
    end = int(toks[2])-region_margin
    if (end > start):
        return (toks[0], start, end)
    else:
        return None

def read_regions(path):
    regions = [to_region(r) for r in open(path, "r") if r[0] != '#']
    return [reg for reg in regions if reg is not None]

def gen_insertion(ref, chr, location, size):
    insertion = "".join([random.choice(bases) for _ in range(size)])
    ref_base = ref.fetch(chr, location, location+1)
    return "\t".join([chr, str(location+1), ".", ref_base, ref_base + insertion])

def gen_deletion(ref, chr, location, size):
    ref_bases = ref.fetch(chr, location, location+size+1)
    alt = ref_bases[0]
    return "\t".join([chr, str(location+1), ".", ref_bases, alt])

def gen_duplication(ref, chr, location, size):
    ref_bases = ref.fetch(chr, location, location+size+1)
    alt = ref_bases + ref_bases
    return "\t".join([chr, str(location+1), ".", ref_bases, alt])

def gen_inverse_dup(ref, chr, location, size):
    ref_bases = ref.fetch(chr, location, location+size+1)
    alt = ref_bases + revcomp(ref_bases)
    return "\t".join([chr, str(location+1), ".", ref_bases, alt])

def gen_snp(ref, chr, location, size):
    ref_bases = ref.fetch(chr, location, location+1)
    alt = random.choice(bases)
    while alt == ref_bases:
        alt = random.choice(bases)
    return "\t".join([chr, str(location+1), ".", ref_bases, alt])


def gen_blocksub(ref, chr, location, size):
    """
    Length-preserving variant
    :param ref:
    :param chr:
    :param location:
    :param size:
    :return:
    """
    ref_bases = ref.fetch(chr, location, location+size+1)
    alt = ""
    for b in ref_bases:
        a = random.choice(bases)
        while a==b:
            a = random.choice(bases)
        alt = alt + a
    return "\t".join([chr, str(location+1), ".", ref_bases, alt])

def gen_ins_mnp(ref, chr, location, size):
    insertion = "".join([random.choice(bases) for _ in range(size)])
    ref_base = ref.fetch(chr, location, location+1)
    alt = random.choice(bases)
    while alt == ref_base:
        alt = random.choice(bases)
    return "\t".join([chr, str(location+1), ".", ref_base, ref_base + insertion + "," + ref_base + alt])

def gen_del_snp_mnp(ref, chr, location, size):
    ref_bases = ref.fetch(chr, location, location+size+1)
    del_alt = ref_bases[0]
    okbases = ['A', 'C', 'T', 'G']
    okbases.remove(ref_bases[1])
    alt = ref_bases[0] + random.choice(okbases) + ref_bases[2:]
    return "\t".join([chr, str(location+1), ".", ref_bases, alt + "," + del_alt ])

def pick_location(regions, blacklist=None, min_safe_dist = 1000):
    region = random.choice(regions)
    loc = random.randint(region[1], region[2])
    if blacklist is not None:
        hits = blacklist[region[0]][loc-min_safe_dist:loc+min_safe_dist]
        while len(hits)>0:
            region = random.choice(regions)
            loc = random.randint(region[1], region[2])
            hits = blacklist[region[0]][loc-min_safe_dist:loc+min_safe_dist]
        blacklist[region[0]].addi(loc, loc+1)


    return (region[0], loc)



def generate_all(ref, regions, output):
    reps_per_size = 20
    repeats = 1
    blacklist = defaultdict(intervaltree.IntervalTree)
    output.write("##fileformat=VCFv4.1\n")
    output.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    output.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n')
    for rep in range(0, reps_per_size):
        for size in range(1, 121, 10):
            loc = pick_location(regions, blacklist)
            #var = gen_deletion(ref, loc[0], loc[1], size)
            #var = gen_del_snp_mnp(ref, loc[0], loc[1], size)
            #var = gen_insertion(ref, loc[0], loc[1], size)
            # var = gen_snp(ref, loc[0], loc[1], size)
            #var = gen_duplication(ref, loc[0], loc[1], size)
            var = gen_blocksub(ref, loc[0], loc[1], size)
            #var = gen_inverse_dup(ref, loc[0], loc[1], size)
            for _ in range(repeats):
                output.write(var + "\t" + "\t".join(['.', '.', '.']) + "\n")


def revcomp(bases):
    """
    Return reverse-complemented bases, uses the revcomp_lookup global lookup dictionary
    :return: Reverse-complemented bases as a string
    """
    result = []
    revcomp_lookup={'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G' }
    for b in bases[::-1]:
        result.append(revcomp_lookup[b])
    return "".join(result)

if __name__=="__main__":
    parser = argparse.ArgumentParser("Create a VCF file with many simulated indels")
    parser.add_argument("-f", "--fasta", help="Reference fasta")
    parser.add_argument("-r", "--regions", help="Input BED file with target regions")
    args = parser.parse_args()

    regions = read_regions(args.regions)
    ref = pysam.FastaFile(args.fasta)

    generate_all(ref, regions, sys.stdout)


