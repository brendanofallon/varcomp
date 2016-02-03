
import subprocess
import gzip
import random
import string
from collections import namedtuple

HOM_REF_GT = "Hom ref."
HET_GT = "Het"
HOM_ALT_GT = "Hom alt."
HET_NONREF = "Het, non ref"
HET_WITHREF = "Triploid / complex"

ALL_HOMREF_GTS = ["0/0", "0|0"]
ALL_HET_GTS = ["0/1", "1/0", "1|0", "0|1"]
ALL_HOMALT_GTS = ["1/1", "1|1"]

DEFAULT_CONTIG_ORDER=['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6', '7', '8','9', 'MT', 'X','Y']

def var_comp(v1, v2):
    """
    Comparator for two variants - chromosome first (according to DEFAULT_CONTIG_ORDER), then position
    :param v1:
    :param v2:
    :return:
    """
    v1c = DEFAULT_CONTIG_ORDER.index(v1[0])
    v2c = DEFAULT_CONTIG_ORDER.index(v2[0])
    if v1c == v2c:
        return int(v1[1]) - int(v2[1])
    else:
        return v1c - v2c

def sort_vcf(vcf, conf):

    tmpfile = vcf.replace(".vcf", ".sort" + randstr() + ".vcf").replace(".gz", "")
    vars = []
    ofh = open(tmpfile, "w")
    if vcf.endswith(".gz"):
        fh = gzip.open(vcf)
    else:
        fh = open(vcf)

    for line in fh.readlines():
         if len(line)>0 and line[0]=='#':
             ofh.write(line)
         else:
             vars.append(line.split('\t'))
    for var in sorted(vars, cmp=var_comp):
        ofh.write('\t'.join(var))

    ofh.close()
    fh.close()
    return bgz_tabix(tmpfile, conf)

def randstr(length=8):
    return "".join([random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(length)])

def bgz_tabix(path, conf):
    """
    If the path does not end in .gz bgzip the file, then index with tabix and return the potentially modified filename
    :return: Filename of compressed file
    """

    if not path.endswith(".gz"):
        cmd = conf.get('main', 'bgzip_path') + " " + path
        subprocess.check_call(cmd.split())
        path = path + ".gz"
    cmd = conf.get('main', 'tabix_path') + " -f " + path
    subprocess.check_call(cmd.split())

    return path



def set_genotypes(orig_vcf, newGT, region, conf):
    """
    Create a new VCF file that is identical to the given VCF, except that all GT info fields are set to 'newGT'
    """
    fh = None
    if orig_vcf.endswith(".gz"):
        fh = gzip.open(orig_vcf, "r")
    else:
        fh = open(orig_vcf, "r")

    newvcf = orig_vcf.replace(".vcf", ".gtmod" + randstr() + ".vcf").replace(".gz", "")
    ofh =open(newvcf, "w")
    for line in fh.readlines():
        if len(line)==0 or line[0]=='#':
            ofh.write(line)
        else:
            toks = line.split('\t')
            chr = toks[0]
            start = int(toks[1])

            if region is not None and (chr != region.chr or start<region.start or start>=region.end):
                ofh.write(line)
                continue

            if len(toks)<10:
                ofh.write(line)
            else:
                if "," in toks[4]:
                    raise GTModException('Cant set GT for multi-alt variants.')
                infoitems = [newGT]
                if ':' in toks[9]:
                    infoitems.extend(toks[9].strip().split(':')[1:] )
                newinfo = ":".join(infoitems)
                ofh.write('\t'.join(toks[0:9] + [newinfo]) + "\n")
    fh.close()
    ofh.close()
    bgz_vcf = bgz_tabix(newvcf, conf)
    return bgz_vcf

class GTModException(Exception):

    def __init__(self, msg=None):
        self.msg = msg

    def __str__(self):
        return "GT modification exception: " + msg

def region_to_bedfile(region):
    """
    Write the given region to its own one-line bed file, return the filename
    :param region:
    :return:
    """
    filename = "tmpbed-" + randstr() + ".bed"
    with open(filename, "w") as fh:
        fh.write("\t".join([region.chr, str(region.start), str(region.end)]) + "\n")
    return filename

def vars_to_bed(variants, window=500):
    """
    Generate a bed file containing regions that span each variant, each region is centered
    on the variant start position and extends 'window' bp in each direction
    The resulting file is NOT sorted by
    :param variants: List of variants to create a bed file for
    :param window:Number of bases to extend (from variant start position) in each direction
    :return: Name of bed file created
    """
    bedfilename = "var_regions" + randstr() + ".bed"
    with open(bedfilename, "w") as bfh:
        for var in variants:
            bfh.write("\t".join([var.chrom, str(var.start-window), str(var.start+window)]) + "\n")

    return bedfilename

def read_regions(bedfile):
    """
    Generator for iterating over
    :param bedfile:
    :return:
    """
    Region = namedtuple('Region', ['chr', 'start', 'end'])
    for line in open(bedfile).readlines():
        if len(line)==0 or line[0]=='#':
            continue
        toks = line.split('\t')
        yield Region(toks[0], int(toks[1]), int(toks[2]))

def find_matching_var(vars, region):
    """
    Collect a list of the variants in the vars list whose start position is contained in the given region
    :param vars:
    :param region:
    :param gt:
    :param conf:
    :return:
    """
    matches = [var for var in vars if var.chrom==region.chr and var.start >= region.start and var.start <= region.end]
    return matches


def write_vcf(variants, filename, conf, gt="1/1"):
    """
    Write the variants in the list to a vcf file. Doesn't do any sorting. By default genotype fields are hom alt (1/1)
    :param variants: List of input variants to write
    :param filename: Destination filename of vcf to write
    :param conf: Configuration (needed for paths to tabix, bgzip)
    :param gt: Genotype field
    :return:
    """
    fh = open(filename, "w")
    fh.write("##fileformat=VCFv4.1\n")
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n')
    for variant in variants:
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
    return compress_vcf(filename, conf)


def compress_vcf(input_vcf, conf):
    """
    If the input vcf's filename does not end with .gz, compress and index it with bgzip / tabix
    :param input_vcf:
    :param conf:
    :return: Name of compressed vcf file, typically input_vcf + '.gz'
    """
    if not input_vcf.endswith(".gz"):
        cmd = conf.get('main', 'bgzip_path') + " -f " + input_vcf
        subprocess.check_call(cmd, shell=True)
        cmd = conf.get('main', 'tabix_path') + " -f " + input_vcf + ".gz"
        subprocess.check_call(cmd, shell=True)
        input_vcf = input_vcf + '.gz'
    return input_vcf


def get_first_gt(var):
    """
    Returns string version of GT field.. Hack until we can get pysam to work..
    :param var:
    """
    sample = var.samples[0]
    if 'GT' not in sample:
        return None
    gts = sample['GT']
    refcount = len([g for g in gts if g==0])
    altcount = len(gts)-refcount
    alt_types = set([g for g in gts if g!=0])

    if len(alt_types)==0: #GT is like 0
        return HOM_REF_GT

    if len(alt_types)==1: #GT is like 0/1 or 1/1 or 0/2
        if altcount==0:
            return HOM_REF_GT
        if refcount==1 and altcount==1:
            return HET_GT
        return HOM_ALT_GT

    if len(alt_types)>1: #GT is like 1/2 or 0/1/2...
        if refcount==0:
            return HET_NONREF
        else:
            return HET_WITHREF


def canadd(var, batch, max_batch_size, min_safe_dist=2000):
    """
    Helper for variant batching function
    :param var: Single variant
    :param batch: A batch to consider adding the variant to
    :param max_batch_size: Maximum size of a batch
    :param min_safe_dist: Minimum distance required between variants in the batch
    :return:
    """
    if len(batch)>=max_batch_size:
        return False
    for b in batch:
        if var.chrom == b.chrom and abs(b.start - var.start)<min_safe_dist:
            return False
    return True

def batch_variants(vars, max_batch_size=50, min_safe_dist=2000):
    """
    Given a list of variants, group them into batches such that no batch contains two variants
    whose start positions are within min_safe_dist bases of each other
    :param vars: List of variants
    :param max_batch_size: Maximum number of variants per batch
    :param min_safe_dist: Min permissible distance between two variants in batch
    :return: List of batches containing variants
    """
    batches = []
    vars = list(vars)
    while len(vars)>0:
        var = vars.pop(0)
        unfilled_batches = [b for b in batches if len(b)<max_batch_size]
        found = False
        for b in unfilled_batches:
            if canadd(var, b, max_batch_size, min_safe_dist=min_safe_dist):
                b.append(var)
                found = True
                break

        if not found:
            #Need to make a new batch for this variant, it doesn't fit anywhere
            batch = []
            batch.append(var)
            batches.append(batch)

    return batches