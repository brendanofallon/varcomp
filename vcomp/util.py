import os
import random
import string
import subprocess
import traceback as tb

from collections import namedtuple

import pysam

from kevinlib.core.smartfile import smartfile

ALL_HOMREF_GTS = ["0/0", "0|0"]
ALL_HET_GTS = ["0/1", "1/0", "1|0", "0|1", "1/2", "2/1", "1|2", "2|1"]
ALL_HOMALT_GTS = ["1/1", "1|1", "2/2", "2|2"]

VCF_MISSING = "."

DEFAULT_CONTIG_ORDER = ['chr1',  'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
                        'chr17', 'chr18', 'chr19', 'chr2',  'chr20', 'chr21', 'chr22', 'chr3',
                        'chr4',  'chr5',  'chr6',  'chr7',  'chr8',  'chr9',  'chrX',  'chrY', 'chrM']

Region = namedtuple('Region', ['chrom', 'start', 'stop'])
Variant = namedtuple('Variant', ['chrom', 'start', 'ref', 'alts', 'gt'])
ErrorVariant = namedtuple('ErrorVariant', ['chrom', 'start', 'msg'])


def strip_extensions(filename, exts):
    '''
    Strip listed file extensions from filename from the left.  This process
    Stops at the first unrecognized extension.

    >>> strip_extensions('/usr/local/foo.bar', 'bar')
    '/usr/local/foo'
    >>> strip_extensions('foo.g.vcf.gz', ['g','vcf','gz'])
    'foo'
    >>> strip_extensions('/a/b/c/foo.g.vcf.gz', ['g','vcf','gz'])
    '/a/b/c/foo'
    >>> strip_extensions('/a/b/c/foo.Z.g.vcf.gz', ['Z'])
    '/a/b/c/foo.Z.g.vcf.gz'
    '''
    if isinstance(exts, basestring):
        exts = [exts]
    prefix,suffix = os.path.split(filename)
    parts = suffix.split('.')
    while parts and parts[-1] in exts:
        parts.pop()
    return os.path.join(prefix, '.'.join(parts))


def var_comp(v1, v2, contig_list=DEFAULT_CONTIG_ORDER):
    """
    Comparator for two tokenized VCF lines - chromosome first (according to DEFAULT_CONTIG_ORDER), then position
    :param v1: First variant
    :param v2: Second variant
    :return: Sorting key
    """
    v1c = contig_list.index(v1[0])
    v2c = contig_list.index(v2[0])
    if v1c == v2c:
        return int(v1[1]) - int(v2[1])
    else:
        return v1c - v2c


def variant_comp(v1, v2, contig_list=DEFAULT_CONTIG_ORDER):
    """
    Adapter to make var_comp sorter compatible with actual VCF variants
    :param v1: First variant
    :param v2: Second variant
    :return: Sorting key
    """
    return var_comp((v1.chrom, v1.start), (v2.chrom, v2.start), contig_list=contig_list)


def sort_vcf(vcf, conf):
    tmpfile = vcf.replace(".vcf", ".sort" + randstr() + ".vcf").replace(".gz", "")
    vars = []
    ofh = open(tmpfile, "w")

    for line in smartfile(vcf):
        if line.startswith('#'):
            ofh.write(line)
        else:
            vars.append(line.split('\t'))

    for var in sorted(vars, cmp=var_comp):
        ofh.write('\t'.join(var))

    ofh.close()

    return bgz_tabix(tmpfile, conf)


def randstr(length=8):
    return "".join(
        [random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(length)])


def bgz_tabix(filename, conf):
    """
    If the filename does not end in .gz bgzip the file, then index with tabix and return the potentially modified filename
    :return: Filename of compressed file
    """
    if not filename.endswith(".gz"):
        cmd = conf.get('main', 'bgzip') + " -f " + filename
        subprocess.check_call(cmd, shell=True)
        filename += ".gz"
    cmd = conf.get('main', 'tabix') + " -f " + filename
    subprocess.check_call(cmd, shell=True)
    return filename


def pysamVar_to_Variant(pvar, default_gt):
    if 'GT' in pvar.format:
        gt = '/'.join(pvar.samples[0]['GT'].allele_indices)
    else:
        if default_gt is None:
            raise ValueError('No default GT specified, and variant does not contain GT information: ' + str(pvar))

        # Special case, if default_gt is "0/1" or "0|1", and there are multiple alts, assume we want
        # het-alt case
        # FIXME: Why is this special case needed?
        if default_gt in ALL_HET_GTS and len(pvar.alts) == 2:
            gt = "1|2"
        else:
            gt = default_gt

    return Variant(pvar.chrom, pvar.start, pvar.ref, pvar.alts, gt)


def set_genotypes(orig_vcf, newGT, region, conf, compress_result=True):
    """
    Create a new VCF file that is identical to the given VCF, except that all GT info fields are set to 'newGT'
    """
    # FIXME: Update to use pysam.VariantFile
    newvcf = '{}.gtmod.{}.vcf'.format(strip_extensions(orig_vcf, ['gz','vcf']), randstr())
    ofh = open(newvcf, 'w')
    for line in smartfile(orig_vcf):
        if line.startswith('#'):
            ofh.write(line)
        else:
            toks  = line.split('\t')
            chrom = toks[0]
            start = int(toks[1])

            # FIXME: Interval intersection logic off?  stop=start+rlen is not taken into consideration.
            if region is not None and (chrom != region.chrom or start < region.start or start >= region.stop):
                ofh.write(line)
                continue

            if len(toks) < 10:
                ofh.write(line)
            else:
                if ',' in toks[4]:
                    raise GTModException('Cant set GT for multi-alt variants.')
                if toks[4] == VCF_MISSING:
                    infoitems = ['0']
                else:
                    infoitems = [newGT]
                if ':' in toks[9]:
                    infoitems.extend(toks[9].strip().split(':')[1:])
                newinfo = ':'.join(infoitems)
                ofh.write('\t'.join(toks[0:9] + [newinfo]) + '\n')

    ofh.close()
    if compress_result:
        newvcf = bgz_tabix(newvcf, conf)

    return newvcf


class GTModException(Exception):
    def __init__(self, msg=None):
        self.msg = msg

    def __str__(self):
        return "GT modification exception: " + self.msg


def region_to_bedfile(region):
    """
    Write the given region to its own one-line bed file, return the filename
    :param region:
    :return:
    """
    filename = "tmpbed-" + randstr() + ".bed"
    with open(filename, "w") as fh:
        fh.write("\t".join([region.chrom, str(region.start), str(region.stop)]) + "\n")
    return filename


def vars_to_bed(variants, window=500):
    """
    Generate a bed file containing regions that span each variant, each region is centered
    on the variant start position and extends 'window' bp in each direction
    Nearby variants are merged into the same region.
    :param variants: List of variants to create a bed file for
    :param window:Number of bases to extend (from variant start position) in each direction
    :return: Name of bed file created
    """
    bedfilename = "var_regions" + randstr() + ".bed"
    cur_chrom = None
    interval_start = -1
    interval_stop = -1
    with open(bedfilename, "w") as bfh:
        for item in variants:
            if type(item) == dict:
                var = item['vars'][0]
            else:
                var = item
            if var.chrom == cur_chrom and var.start >= interval_start and (var.start + len(var.ref)) < interval_stop:
                interval_stop = var.start + window
            else:
                if cur_chrom is not None:
                    bfh.write("\t".join([cur_chrom, str(interval_start), str(interval_stop)]) + "\n")
                cur_chrom = var.chrom
                interval_start = var.start - window
                interval_stop = var.start + window

        if cur_chrom is not None:
            bfh.write("\t".join([cur_chrom, str(interval_start), str(interval_stop)]) + "\n")

    return bedfilename


def remove_halfcalls(vcf):
    """
    Create a new VCF file by removing all half-called variants from the original vcf file
    """
    newvcf = os.path.basename(vcf).replace(".vcf", ".nohalfcalls.vcf").replace(".gz", "")
    cmd = "zgrep '^#' " + vcf + " > " + newvcf
    subprocess.check_call(cmd, shell=True)
    with open(newvcf, "aw") as ofh:
        for var in pysam.VariantFile(vcf):
            if any(None in s.allele_indices for s in var.samples):
                continue
            ofh.write(str(var))
    return newvcf


def read_regions(bedfile):
    """
    Generator for creating / iterating over Regions in a BED file
    :param bedfile: Path to .BED formatted file to read
    :return: generator yielding Regions for each line in input file
    """
    for line in open(bedfile):
        if line.startswith('#'):
            continue
        toks = line.split('\t')
        yield Region(toks[0], int(toks[1]), int(toks[2]))


def is_empty(vcf):
    """
    Return True if the vcf file contains at least one variant, this doesn't use pysam
    :param vcf: Path to possibly-gzipped vcf file
    """
    for line in smartfile(vcf):
        if line and not line.startswith('#'):
            return False
    return True


def find_matching_var(vars, region):
    """
    Collect a list of the variants in the vars list whose start position is contained in the given region
    :param vars:
    :param region:
    :param conf:
    :return: List of variants in region
    """
    if isinstance(vars, str):
        with pysam.VariantFile(vars) as vfile:
            matches = [var for var in vfile.fetch(region.chrom, region.start, region.stop)
                       if var.start >= region.start and var.start <= region.stop]
    else:
        matches = [var for var in vars
                   if var.chrom == region.chrom and var.start >= region.start and var.start <= region.stop]
    return matches


def gen_snp(chrom, pos, gt, ref_genome):
    """
    Generate a randomly chosen SNP variant at the given position
    """
    currentbase = ref_genome.fetch(chrom, pos, pos + 1)
    bases = ['A', 'C', 'G', 'T']
    bases.remove(currentbase)
    newbase = random.choice(bases)
    newvar = Variant(chrom, pos, currentbase, (newbase,), gt)
    return newvar


def write_vcf(variants, filename, conf):
    """
    Write the variants in the list to a vcf file. Doesn't do any sorting. By default genotype fields are hom alt (1/1)
    :param variants: List of input variants to write
    :param filename: Destination filename of vcf to write
    :param conf: Configuration (needed for paths to tabix, bgzip)
    :param gt: Genotype field
    """
    fh = open(filename, "w")
    fh.write("##fileformat=VCFv4.1\n")
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n')
    for variant in sorted(variants, cmp=variant_comp):
        fh.write(variant.chrom + "\t")
        fh.write(str(variant.start + 1) + "\t")  # Remember - internally 0-based coords, but in vcf 1-based
        fh.write("." + "\t")
        fh.write(variant.ref + "\t")
        fh.write(",".join(variant.alts) + "\t")
        fh.write("100" + "\t")
        fh.write("PASS" + "\t")
        fh.write("." + "\t")
        fh.write("GT" + "\t")
        fh.write(variant.gt + "\n")
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
        cmd = conf.get('main', 'bgzip') + " -f " + input_vcf
        subprocess.check_call(cmd, shell=True)
        cmd = conf.get('main', 'tabix') + " -f " + input_vcf + ".gz"
        subprocess.check_call(cmd, shell=True)
        input_vcf = input_vcf + '.gz'
    return input_vcf


def canadd(var, batch, max_batch_size, min_safe_dist=2000):
    """
    Helper for variant batching function
    :param var: Single variant
    :param batch: A batch to consider adding the variant to
    :param max_batch_size: Maximum size of a batch
    :param min_safe_dist: Minimum distance required between variants in the batch
    :return:
    """
    if len(batch) >= max_batch_size:
        return False
    for b in batch:
        if var.chrom == b.chrom and abs(b.start - var.start) < min_safe_dist:
            return False
    return True


def batch_variants(vcf, max_batch_size=1000, min_safe_dist=2000):
    """
    Given a list of variants, group them into batches such that no batch contains two variants
    whose start positions are within min_safe_dist bases of each other
    :param vcf: Filename of VCF file containing variants to batch
    :param max_batch_size: Maximum number of variants per batch
    :param min_safe_dist: Min permissible distance between two variants in batch
    :return: List of VCF files containing subsets of variants
    """

    batches = []
    header = []

    for x in smartfile(vcf):
        if x.startswith('#'):
            header.append(x)
        else:
            break

    name = strip_extensions(os.path.basename(vcf), ['gz','vcf'])

    vars = list(pysam.VariantFile(vcf))

    while vars:
        var = vars.pop(0)
        unfilled_batches = [b for b in batches if len(b) < max_batch_size]
        found = False
        for b in unfilled_batches:
            if canadd(var, b, max_batch_size, min_safe_dist=min_safe_dist):
                b.append(var)
                found = True
                break

        if not found:
            # Need to make a new batch for this variant, it doesn't fit anywhere
            batches.append([var])

    # return batches
    files = []
    for i, batch in enumerate(batches, 1):
        batchname = '{}.batch{:03d}.'.format(name, i) + randstr() + ".vcf"
        with open(batchname, 'w') as out:
            for x in header:
                out.write(x)
            for x in batch:
                out.write(str(x))
        files.append(batchname)
    return files


class TempDir(object):
    """
    On entry, create a new temporary directory and chdir into it. Move back to the original dir on exit.
    If any exceptions are encountered, try to write them to a file in the temporary directory.
    """

    ALWAYS_DELETE = "always delete"
    DELETE_NO_EXCEPTION = "delete no exception"
    NEVER_DELETE = "never delete"

    DELETION_POLICIES = [ALWAYS_DELETE, DELETE_NO_EXCEPTION, NEVER_DELETE]

    def __init__(self, dirname=None, del_policy=DELETE_NO_EXCEPTION):
        if dirname is None:
            self.dirname = "tmp-working-" + randstr()
        else:
            self.dirname = dirname
        self.origdir = None
        self.deletion_policy = del_policy
        if not del_policy in TempDir.DELETION_POLICIES:
            raise ValueError('Unrecognized deletion policy: ' + str(del_policy))

    def __enter__(self):
        self.origdir = os.getcwd()
        try:
            os.mkdir(self.dirname)
        except:
            pass
        os.chdir(self.dirname)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            try:
                with open("exception.info.txt", "a") as fh:
                    tb.print_tb(exc_tb, file=fh)
                    fh.write(str(exc_type) + " : " + str(exc_val) + "\n")
            except:
                pass
        os.chdir(self.origdir)
        if self.deletion_policy == TempDir.ALWAYS_DELETE \
                or (self.deletion_policy == TempDir.DELETE_NO_EXCEPTION and exc_val is None):
            os.system("rm -rf " + self.dirname)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
