
import pysam
import subprocess
import sys
import os
import ConfigParser as cp

NO_VARS_FOUND_RESULT="No variants identified"
MATCH_RESULT="Variants matched"
INCORRECT_GENOTYPE_RESULT="Genotype mismatch"
NO_MATCH_RESULT="Variants did not match"
PARTIAL_MATCH="Partial variant match"

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
    except:
        pass
    return vars

def test_var_equiv(var1, var2):
    alts1 = str([str(a) + "," for a in var1.alts])
    alts2 = str([str(a) +","  for a in var2.alts])
    return var1.chrom == var2.chrom and var1.start == var2.start and var1.ref == var2.ref and  alts1==alts2

def compare_raw(orig_vars, caller_vars):
    """
    Compare all original variant to all caller vars to see if there are any exact matches

    :param orig_vars:
    :param caller_vars:
    :return:
    """
    if len(orig_vars)>0 and len(caller_vars) == 0:
        return NO_VARS_FOUND_RESULT

    if len(caller_vars)>0 and len(orig_vars)>0:
        caller_gt = get_first_gt(caller_vars[0])
        orig_gt = get_first_gt(orig_vars[0])
        if caller_gt != orig_gt:
            return INCORRECT_GENOTYPE_RESULT

    #There are multiple input and/or output vars
    #For every input var, is there at least one match?
    matches = []
    unmatched = []
    for ovar in orig_vars:
        match_found = False
        for cvar in caller_vars:
            match_found = test_var_equiv(ovar, cvar)
            if match_found:
                break
        if match_found:
            matches.append( (ovar, cvar) )
        else:
            unmatched.append(ovar)

    if len(matches)>0 and len(unmatched)==0:
        return MATCH_RESULT
    if len(matches)==0 and len(unmatched)>0:
        return NO_MATCH_RESULT
    if len(matches)>0 and len(unmatched)>0:
        return PARTIAL_MATCH

def get_first_gt(var):
    """
    Hack until we can get pysam to work..
    :param var:
    :return:
    """
    toks = str(var).split()
    if len(toks) <= 9:
        return None

    return toks[9].split(":")[0]

def compare_nonorm(orig_vcf, caller_vcf, conf=None):
    """
    Compare variants with absolutely no normalization
    :return:
    """
    orig_vars = read_all_vars(orig_vcf)
    caller_vars = read_all_vars(caller_vcf)
    return compare_raw(orig_vars, caller_vars)


def compare_vt(orig_vcf, caller_vcf, conf):
    """
    Compare variants after using vt to normalize
    :param conf: configuration object with path to reference genome and vt binary
    :return: String describing variant matching result
    """
    err = open("/dev/null")
    norm_orig_vcf = orig_vcf.replace(".vcf", ".norm.vt.vcf")
    norm_orig_cmd = conf.get('main', 'vt_path') + " normalize " + " -r " + conf.get('main', 'ref_genome') + " " + orig_vcf + " -o " + norm_orig_vcf
    subprocess.check_output(norm_orig_cmd.split(), stderr=err)
    normed_orig_vars = read_all_vars(norm_orig_vcf)

    norm_caller_vcf = caller_vcf.replace(".vcf", ".norm.vt.vcf")
    norm_caller_cmd = conf.get('main', 'vt_path') + " normalize " + " -r " + conf.get('main', 'ref_genome') + " " + caller_vcf + " -o " + norm_caller_vcf
    subprocess.check_output(norm_caller_cmd.split(), stderr=err)
    normed_caller_vars = read_all_vars(norm_caller_vcf)
    err.close()
    return compare_raw(normed_orig_vars, normed_caller_vars)

def compare_bcftools(orig_vcf, caller_vcf, conf):
    norm_orig_vcf = orig_vcf.replace(".vcf", ".norm.bcftools.vcf")
    norm_orig_cmd = conf.get('main', 'bcftools_path') + " norm " + " -c w -f " + conf.get('main', 'ref_genome') + " " + orig_vcf + " -o " + norm_orig_vcf
    subprocess.check_call(norm_orig_cmd.split())
    normed_orig_vars = read_all_vars(norm_orig_vcf)

    norm_caller_vcf = caller_vcf.replace(".vcf", ".norm.bcftools.vcf")
    norm_caller_cmd = conf.get('main', 'bcftools_path') + " norm " + " -c w -f " + conf.get('main', 'ref_genome') + " " + caller_vcf + " -o " + norm_caller_vcf
    subprocess.check_output(norm_caller_cmd.split(), stderr=None)
    normed_caller_vars = read_all_vars(norm_caller_vcf)
    return compare_raw(normed_orig_vars, normed_caller_vars)


def compare_vgraph(orig_vcf, caller_vcf, conf, bed=None):
    #Slightly tricky - must make sure executed python script is executed in a workable virtualenv. For now we
    #execute it in whatever virtual env this script is being processed in.

    #First, test to see if there are any vars present...
    orig_var_count = len(read_all_vars(orig_vcf, bed))
    caller_vars = read_all_vars(caller_vcf, bed)
    caller_var_count = len(caller_vars)
    if orig_var_count > 0 and caller_var_count == 0:
        return NO_VARS_FOUND_RESULT

    for var in caller_vars:
        gt = get_first_gt(var)
        if gt != "1/1":
            return INCORRECT_GENOTYPE_RESULT

    bedcmd = ""
    if bed is not None:
        bedcmd = " -i " + bed
    vg_cmd = conf.get('main', 'vgraph_path') + " --reference " + conf.get('main', 'ref_genome') + bedcmd + " " + orig_vcf + " " + caller_vcf
    cmd = [sys.executable]
    cmd.extend(vg_cmd.split())
    result = subprocess.check_output(cmd, env=os.environ.copy())
    if "MISMATCH!" in result:
        return NO_MATCH_RESULT
    elif "MATCH!" in result:
        return MATCH_RESULT
    return NO_MATCH_RESULT




def get_comparators():
    return {
        "nothing": compare_nonorm,
        "vt": compare_vt,
        "bcftools": compare_bcftools,
        "vgraph": compare_vgraph
    }


if __name__=="__main__":
    conf = cp.SafeConfigParser()
    conf.read("comp.conf")
    compare_vgraph("tmp-workingRPVTKsVf/test2.vcf.gz", "tmp-workingRPVTKsVf/output-fb.vcf.gz", conf, "tmp-workingRPVTKsVf/var_regionsFoTUZUhjvJ.bed")