
import pysam
import subprocess
import sys
import os
import gzip
import time
import ConfigParser as cp


ZYGOSITY_MATCH="Zygosity match"
ZYGOSITY_EXTRA_ALLELE="Extra allele"
ZYGOSITY_MISSING_ALLELE="Missing allele"
ZYGOSITY_MISSING_TWO_ALLELES="Missing two alleles!"

ALLELE_MATCH="Alleles matched"
ALLELE_MISMATCH="Alleles did not match"
ALLELE_EXTRA="Additional variants identified"

NO_VARS_FOUND_RESULT="No variants identified"
MATCH_RESULT="Variants matched"
INCORRECT_GENOTYPE_RESULT="Genotype mismatch"
NO_MATCH_RESULT="Variants did not match"
PARTIAL_MATCH="Partial variant match"

all_result_types = (NO_VARS_FOUND_RESULT, MATCH_RESULT, INCORRECT_GENOTYPE_RESULT, NO_MATCH_RESULT)
genotype_result_types = (NO_VARS_FOUND_RESULT, INCORRECT_GENOTYPE_RESULT)
var_result_types = set(all_result_types)-set(genotype_result_types)


HOM_REF_GT = "Hom ref."
HET_GT = "Het"
HOM_ALT_GT = "Hom alt."

hom_ref_gts = ["0/0", "0|0"]
het_gts = ["0/1", "1/0", "1|0", "0|1"]
hom_alt_gts = ["1/1", "1|1"]



class VariantCompResult(object):

    def __init__(self, orig_vars, called_vars):
        self.orig_vars = orig_vars
        self.called_vars = called_vars
        #Try to find matches among the orig_vars and the called vars, and keep track of which vars
        #haven't been matched
        (orig_unmatched, matches, called_unmatched) = self._match_vars()
        self.orig_unmatched = orig_unmatched
        self.matches = matches
        self.called_unmatched = called_unmatched

    @staticmethod
    def from_raw(orig_unmatched, matched, called_unmatched):
        """
        Generate a new VariantCompResult object from variants that have already been matched - this is used
        by vgraph and vcfeval callers which handle their own variant matching
        :param orig_unmatched:
        :param matched:
        :param called_unmatched:
        :return:
        """
        res = VariantCompResult([], [])
        res.orig_unmatched = orig_unmatched
        res.matches = matched
        res.called_unmatched = called_unmatched
        return res

    def _match_vars(self):
        matches = []
        unmatched_orig = []
        unmatched_caller = [v for v in self.called_vars]
        for ovar in self.orig_vars:
            match_found = False
            for cvar in self.called_vars:
                match_found = test_var_equiv(ovar, cvar)
                if match_found:
                    break
            if match_found:
                matches.append( (ovar, cvar) )
                unmatched_caller.remove(cvar)
            else:
                unmatched_orig.append(ovar)

        return (unmatched_orig, matches, unmatched_caller)

    def match_count(self):
        return len(self.matches)

    def zygosity_calc(self):
        return [comp_zyg(a,b) for a,b in self.matches]

    def has_unmatched_orig(self):
        return len(self.orig_unmatched)>0

    def has_unmatched_caller(self):
        return len(self.called_unmatched)>0

    def __str__(self):
        return "[" + str(len(self.orig_unmatched)) + ", " + str(len(self.matches)) + ", " + str(len(self.called_unmatched)) + "; " +  ",".join(self.zygosity_calc())  + "]"

    def short(self):
        if len(self.matches)>0 and len(self.orig_unmatched)==0 and len(self.called_unmatched)==0:
            return MATCH_RESULT
        if len(self.matches)==0 and len(self.orig_unmatched)>0:
            return NO_MATCH_RESULT
        if len(self.matches)==0 and len(self.orig_unmatched)==0 and len(self.called_unmatched)==0:
            return NO_VARS_FOUND_RESULT
        return PARTIAL_MATCH

def comp_zyg(tvar, qvar):
    """
    Compare the two variants for genotype (GT format field) and return a ZYGOSITY result of some sort
    :param tvar: True variant
    :param qvar: Query variant
    :return:
    """
    t = get_first_gt(tvar)
    q = get_first_gt(qvar)
    if t == HET_GT:
        if q == HET_GT:
            return ZYGOSITY_MATCH
        if q == HOM_ALT_GT:
            return ZYGOSITY_EXTRA_ALLELE
        if q == HOM_REF_GT:
            return ZYGOSITY_MISSING_ALLELE
    if t == HOM_ALT_GT:
        if q == HET_GT:
            return ZYGOSITY_MISSING_ALLELE
        if q == HOM_ALT_GT:
            return ZYGOSITY_MATCH
        if q == HOM_REF_GT:
            return ZYGOSITY_MISSING_TWO_ALLELES
    raise ValueError('Should never get here - was true var Hom Ref?')


def read_all_vars(vcf, bed=None):
    """
    Try to read all the variants from th egiven vcf file into a list. If there's an error
    reading the vcf, return an empty list
    :param vcf:VCF to read variants from
    :param bed:If not none, only read variants in these regions
    :return:
    """
    vars = []
    try:
        if bed is None:
            vars = list(pysam.VariantFile(vcf))
        else:
            vfh = pysam.VariantFile(vcf)
            for line in open(bed, "r"):
                line = line.strip()
                if len(line)==0 or line[0]=='#':
                    continue
                toks = line.split()
                chr = toks[0]
                reg_start = int(toks[1])
                reg_end = int(toks[2])
                vars.extend(list(vfh.fetch(chr, reg_start, reg_end)))
    except Exception as ex:
        #VCF files with zero variants cause the pysam.VariantFile to throw an exception. We ignore
        #these exceptions and return an empty list
        #print "Wha? " + str(ex)
        pass
    return vars

def test_var_equiv(var1, var2):
    alts1 = str([str(a) + "," for a in var1.alts])
    alts2 = str([str(a) +","  for a in var2.alts])
    gt1 = get_first_gt(var1)
    gt2 = get_first_gt(var2)
    return var1.chrom == var2.chrom and var1.start == var2.start and var1.ref == var2.ref and  alts1==alts2 and gt1 == gt2

def get_first_gt(var):
    """
    Hack until we can get pysam to work..
    :param var:
    :return:
    """
    toks = str(var).split()
    if len(toks) <= 9:
        return None

    gt = toks[9].split(":")[0]
    if gt in hom_ref_gts:
        return HOM_REF_GT
    if gt in hom_alt_gts:
        return HOM_ALT_GT
    if gt in het_gts:
        return HET_GT
    return gt

def compare_nonorm(orig_vcf, caller_vcf, bed, conf):
    """
    Compare variants with absolutely no normalization
    :return:
    """
    orig_vars = read_all_vars(orig_vcf, bed)
    caller_vars = read_all_vars(caller_vcf, bed)
    return VariantCompResult(orig_vars, caller_vars)


def compare_vt(orig_vcf, caller_vcf, bed, conf):
    """
    Compare variants after using vt to normalize
    :param conf: configuration object with path to reference genome and vt binary
    :return: String describing variant matching result
    """
    err = open("/dev/null")
    norm_orig_vcf = orig_vcf.replace(".vcf", ".norm.vt.vcf")
    norm_orig_cmd = conf.get('main', 'vt_path') + " normalize " + " -r " + conf.get('main', 'ref_genome') + " " + orig_vcf + " -o " + norm_orig_vcf
    subprocess.check_output(norm_orig_cmd.split(), stderr=err)
    norm_orig_vcf = bgz_tabix(norm_orig_vcf, conf)
    normed_orig_vars = read_all_vars(norm_orig_vcf, bed)

    norm_caller_vcf = caller_vcf.replace(".vcf", ".norm.vt.vcf")
    norm_caller_cmd = conf.get('main', 'vt_path') + " normalize " + " -r " + conf.get('main', 'ref_genome') + " " + caller_vcf + " -o " + norm_caller_vcf
    subprocess.check_output(norm_caller_cmd.split(), stderr=err)
    norm_caller_vcf = bgz_tabix(norm_caller_vcf, conf)
    normed_caller_vars = read_all_vars(norm_caller_vcf, bed)
    err.close()
    return VariantCompResult(normed_orig_vars, normed_caller_vars)

def compare_bcftools(orig_vcf, caller_vcf, bed, conf):
    norm_orig_vcf = orig_vcf.replace(".vcf.gz", ".norm.bcftools" + str(time.time())[-6:].replace(".", "") + ".vcf")
    norm_orig_cmd = conf.get('main', 'bcftools_path') + " norm " + " -c w -f " + conf.get('main', 'ref_genome') + " " + orig_vcf + " -o " + norm_orig_vcf
    subprocess.check_call(norm_orig_cmd.split())
    norm_orig_vcf = bgz_tabix(norm_orig_vcf, conf)
    normed_orig_vars = read_all_vars(norm_orig_vcf, bed)

    norm_caller_vcf = caller_vcf.replace(".vcf.gz", ".norm.bcftools.vcf")
    norm_caller_cmd = conf.get('main', 'bcftools_path') + " norm " + " -c w -f " + conf.get('main', 'ref_genome') + " " + caller_vcf + " -o " + norm_caller_vcf
    subprocess.check_output(norm_caller_cmd.split(), stderr=None)
    norm_caller_vcf = bgz_tabix(norm_caller_vcf, conf)
    normed_caller_vars = read_all_vars(norm_caller_vcf, bed)
    return VariantCompResult(normed_orig_vars, normed_caller_vars)


def compare_vgraph(orig_vcf, caller_vcf, bed, conf):
    #Slightly tricky - must make sure executed python script is executed in a workable virtualenv. For now we
    #execute it in whatever virtual env this script is being processed in.

    #First, test to see if there are any vars present...
    orig_vars = read_all_vars(orig_vcf, bed)
    orig_var_count = len(orig_vars)
    caller_vars = read_all_vars(caller_vcf, bed)
    caller_var_count = len(caller_vars)
    if orig_var_count > 0 and caller_var_count == 0:
        return NO_VARS_FOUND_RESULT

    bedcmd = ""
    if bed is not None:
        bedcmd = " -i " + bed
    vg_cmd = conf.get('main', 'vgraph_path') + " --reference " + conf.get('main', 'ref_genome') + bedcmd + " " + orig_vcf + " " + caller_vcf
    cmd = [sys.executable]
    cmd.extend(vg_cmd.split())
    result = subprocess.check_output(cmd, env=os.environ.copy())
    if "MISMATCH!" in result:
        return VariantCompResult.from_raw(orig_vars, [], caller_vars)
    elif "MATCH!" in result:
        return VariantCompResult.from_raw([], zip(orig_vars, orig_vars), [])
    return VariantCompResult.from_raw(orig_vars, [], caller_vars)

def bgz_tabix(path, conf):
    """
    If the path does not end in .gz bgzip the file, then index with tabix and return the potentially modified filename
    :return: Filename of compressed file
    """
    try:
        if not path.endswith(".gz"):
            cmd = conf.get('main', 'bgzip_path') + " " + path
            subprocess.check_call(cmd.split())
            path = path + ".gz"

        cmd = conf.get('main', 'tabix_path') + " -f " + path
        subprocess.check_call(cmd.split())
    except Exception as ex:
        raise ex


    return path

def compare_vcfeval(orig_vcf, caller_vcf, bed, conf):
    output_dir = "vcfeval-output" + str(time.time())[-6:].replace(".", "")
    cmd = "java -Xmx2g -jar " + conf.get('main', 'rtg_jar') + " vcfeval -t " + conf.get('main', 'rtg_ref_sdf') + " -o " + output_dir + " -b " + orig_vcf + " -c " + caller_vcf
    if bed is not None:
        cmd = cmd + " --bed-regions " + bed
    subprocess.check_output(cmd, shell=True, executable="/bin/bash")
    orig_vars = read_all_vars(orig_vcf, bed)
    tp_vars = read_all_vars(output_dir + "/tp.vcf.gz")
    fp_vars = read_all_vars(output_dir + "/fp.vcf.gz")
    fn_vars = read_all_vars(output_dir + "/fn.vcf.gz")
    return VariantCompResult.from_raw(fn_vars, zip(orig_vars, tp_vars), fp_vars)



def set_genotypes(orig_vcf, newGT, conf):
    """
    Modify all GT info fields to be whatever 'newGT' is
    """
    fh = None
    if orig_vcf.endswith(".gz"):
        fh = gzip.open(orig_vcf, "r")
    else:
        fh = open(orig_vcf, "r")

    newvcf = orig_vcf.replace(".vcf", ".gtmod" + str(time.time())[-6:].replace(".", "") + ".vcf").replace(".gz", "")
    ofh =open(newvcf, "w")
    for line in fh.readlines():
        if len(line)==0 or line[0]=='#':
            ofh.write(line)
        else:
            toks = line.split('\t')
            if len(toks)<10:
                ofh.write(line)
            else:
                infoitems = [newGT]
                if ':' in toks[9]:
                    infoitems.extend(toks[9].strip().split(':')[1:] )
                newinfo = ":".join(infoitems)
                ofh.write('\t'.join(toks[0:9] + [newinfo]) + "\n")
    fh.close()
    ofh.close()
    bgz_vcf = bgz_tabix(newvcf, conf)
    return bgz_vcf

def get_comparators():
    return {
        "nothing": compare_nonorm,
        "vt": compare_vt,
        "bcftools": compare_bcftools,
        "vgraph": compare_vgraph,
        "vcfeval": compare_vcfeval
    }


if __name__=="__main__":
    conf = cp.SafeConfigParser()
    conf.read("comp.conf")

    vcf = "test_gt.vcf.gz"
    newvcf = set_genotypes(vcf, "1/1", conf)
    print "New vcf file is: " + newvcf

