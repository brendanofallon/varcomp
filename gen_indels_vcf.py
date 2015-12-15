
import argparse
import pysam
import random
import sys

bases = ['A', 'C', 'T', 'G']
region_margin = 20

def to_region(region_str):
    toks = region_str.split()
    return (toks[0], int(toks[1])+region_margin, int(toks[2])-region_margin)

def read_regions(path):
    return [to_region(r) for r in open(path, "r") if r[0] != '#']

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

def pick_location(regions):
    region = random.choice(regions)
    loc = random.randint(region[1], region[2])
    return (region[0], loc)

def generate_all(ref, regions, output):
    reps_per_size = 100
    output.write("##fileformat=VCFv4.1\n")
    output.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    output.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n')
    for rep in range(0, reps_per_size):
        loc = pick_location(regions)
        for size in range(1, 150, 10):
            #var = gen_deletion(ref, loc[0], loc[1], size)
            var = gen_insertion(ref, loc[0], loc[1], size)
            #var = gen_duplication(ref, loc[0], loc[1], size)
            #var = gen_inverse_dup(ref, loc[0], loc[1], size)
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


