
import pysam
import subprocess
import sys
import os
import gzip

NO_VARS_FOUND_RESULT="No variants identified"
MATCH_RESULT="Variants matched"
NO_MATCH_RESULT="Variants did not match"
PARTIAL_MATCH="Partial variant match"

def read_all_vars(vcf):
    """
    Try to read all the variants from th egiven vcf file into a list. If there's an error
    reading the vcf, return an empty list
    :param vcf:
    :return:
    """
    vars = []
    try:
        vars = list(pysam.VariantFile(vcf))
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


def compare_vgraph(orig_vcf, caller_vcf, conf):
    #Slightly tricky - must make sure executed python script is executed in a workable virtualenv. For now we
    #execute it in whatever virtual env this script is being processed in.

    #First, test to see if there are any vars present...
    orig_var_count = len(read_all_vars(orig_vcf))
    caller_var_count = len(read_all_vars(caller_vcf))
    if orig_var_count > 0 and caller_var_count == 0:
        return NO_VARS_FOUND_RESULT

    vg_cmd = conf.get('main', 'vgraph_path') + " --reference " + conf.get('main', 'ref_genome') + " " + orig_vcf + " " + caller_vcf
    cmd = [sys.executable]
    cmd.extend(vg_cmd.split())
    result = subprocess.check_output(cmd, env=os.environ.copy())
    if "MATCH!" in result:
        return MATCH_RESULT
    else:
        return NO_MATCH_RESULT



def get_comparators():
    return {
        "nothing": compare_nonorm,
        "vt": compare_vt,
        "bcftools": compare_bcftools,
        "vgraph": compare_vgraph
    }