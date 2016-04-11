
import pysam
import subprocess
import os
import vcomp.util


ALLELE_MATCH="Alleles matched"
ALLELE_MISMATCH="Alleles did not match"
ALLELE_EXTRA="Additional variants identified"


def get_comparators():
    return {
        #"raw": compare_raw,
        "vgraph": compare_vgraph,
        #"vcfeval": compare_vcfeval,
        #"happy": compare_happy
    }

def compare_raw(orig_vcf, caller_vcf, bed, conf):
    """
    Compare variants with absolutely no normalization. This still requires a little work since there
    are no guarantees about variants being in the right order or the positions of false pos / false neg vars
    :return:
    """
    orig_vars = read_all_vars(orig_vcf, bed)
    caller_vars = read_all_vars(caller_vcf, bed)

    matches = []
    unmatched_orig = []
    unmatched_caller = [v for v in caller_vars]
    for ovar in orig_vars:
        match_found = False
        for cvar in caller_vars:
            match_found = test_var_equiv(ovar, cvar)
            if match_found:
                break
        if match_found:
            matches.append( (ovar, cvar) )
            unmatched_caller.remove(cvar)
        else:
            unmatched_orig.append(ovar)

    return (unmatched_orig, matches, unmatched_caller)



def compare_vgraph(orig_vcf, caller_vcf, bed, conf):
    #Slightly tricky - must make sure executed python script is executed in a workable virtualenv. For now we
    #execute it in whatever virtual env this script is being processed in.

    #First, test to see if there are any vars present...
    orig_vars = read_all_vars(orig_vcf, bed)
    caller_vars = read_all_vars(caller_vcf, bed)
    if len(caller_vars)==0:
        return (read_all_vars(orig_vcf, bed), [], caller_vars)

    orig_out = "vgout-orig." + vcomp.util.randstr() + ".vcf"
    caller_out = "vgout-caller." + vcomp.util.randstr() + ".vcf"
    bedcmd = ""
    if bed is not None:
        bedcmd = " --include-regions " + bed
    vg_cmd = conf.get('main', 'vgraph_path') + " --out1 " + orig_out + " --out2 " + caller_out + " --reference " + conf.get('main', 'ref_genome') + bedcmd + " " + orig_vcf + " " + caller_vcf
    ignored = subprocess.check_output(vg_cmd, env=os.environ.copy(), shell=True)

    unmatched_orig = []
    matches = []
    unmatched_caller = []

    for ovar in pysam.VariantFile(orig_out):
        bd = ovar.samples[0]['BD']
        bk = ovar.samples[0]['BK']

        if bd == 'X':
            unmatched_orig.append(ovar)
        if bd == 'N':
            unmatched_orig.append(
                vcomp.util.ErrorVariant(chrom=ovar.chrom, start=ovar.start, msg="vgraph error code N"))

    for cvar in pysam.VariantFile(caller_out):
        bd = cvar.samples[0]['BD']
        if bd == '=':
            matches.append((cvar, cvar))
        if bd == 'X':
            unmatched_caller.append(cvar)
        if bd == 'N':
            unmatched_caller.append(vcomp.util.ErrorVariant(chrom=cvar.chrom, start=cvar.start, msg="vgraph error code N"))

    return (unmatched_orig, matches, unmatched_caller)



def compare_vcfeval(orig_vcf, caller_vcf, bed, conf):

    #vcfeval will throw an exception if caller vcf is empty, so catch that case here
    caller_vars = read_all_vars(caller_vcf, bed)
    if len(caller_vars)==0:
        return (read_all_vars(orig_vcf, bed), [], caller_vars)

    output_dir = "vcfeval-output" + vcomp.util.randstr()
    cmd = "java -Djava.io.tmpdir=. -Xmx4g -jar " + conf.get('main', 'rtg_jar') + " vcfeval -t " + conf.get('main', 'rtg_ref_sdf') + " --all-records -o " + output_dir + " -b " + orig_vcf + " -c " + caller_vcf
    if bed is not None:
        cmd = cmd + " --bed-regions " + bed
    subprocess.check_output(cmd, shell=True, executable="/bin/bash")
    # orig_vars = read_all_vars(orig_vcf, bed)
    tp_vars = read_all_vars(output_dir + "/tp.vcf.gz")
    fp_vars = read_all_vars(output_dir + "/fp.vcf.gz")
    fn_vars = read_all_vars(output_dir + "/fn.vcf.gz")
    return (fn_vars, zip(tp_vars, tp_vars), fp_vars)


def compare_happy(orig_vcf, caller_vcf, bed, conf):
    orig_vars = read_all_vars(orig_vcf, bed)
    caller_vars = read_all_vars(caller_vcf, bed)
    if len(caller_vars)==0:
        return (read_all_vars(orig_vcf, bed), [], caller_vars)
    output_prefix = "happyoutput-" + vcomp.util.randstr(6)

    all_chrs = set([v.chrom for v in orig_vars])
    all_chrs.update([v.chrom for v in caller_vars])
    all_chrs = list(all_chrs)

    bedarg = ""
    if bed is not None:
        bedarg = " -T " + bed
    cmd = conf.get('main', 'happy_path') + " " + orig_vcf + " " + caller_vcf + " " + bedarg + " -o " + output_prefix + " --scratch-prefix=. --include-nonpass -r " + conf.get('main', 'ref_genome') + " -l " + ",".join(all_chrs) + " --no-fixchr-truth --no-fixchr-query -V"
    ignored = subprocess.check_output(cmd, shell=True)

    orig_unmatched = []
    matches = []
    caller_unmatched = []

    for var in pysam.VariantFile(output_prefix + ".vcf.gz"):
        vstr = str(var)
        if "type=FP" in vstr:
            caller_unmatched.append(var)
            if "kind=gtmismatch" in vstr or "kind=alpartial" in vstr:
                orig_unmatched.append(var)
        elif "type=TP" in vstr:
            matches.append( (var, var) )
        elif "type=FN" in vstr:
            orig_unmatched.append(var)
        else:
            raise ValueError('Unable to parse hap.py variant type: ' + vstr)
    return (orig_unmatched, matches, caller_unmatched)


def read_all_vars(vcf, bed=None):
    """
    Try to read all the variants from the given vcf file into a list. If there's an error
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
    """
    Compare two vars for equality of chrom, pos, ref, and alt and return True if everything is equal.
    :param var1:
    :param var2:
    :return: True if both variant records are identical (contain same alts with same GT fields)
    """
    v1alts = var1.alts if var1.alts is not None else [] #I guess these can sometimes be None?
    v2alts = var2.alts if var2.alts is not None else []

    alts1 = str([str(a) + "," for a in v1alts])
    alts2 = str([str(a) +","  for a in v2alts])
    gt1 = vcomp.util.get_first_gt(var1)
    gt2 = vcomp.util.get_first_gt(var2)
    return var1.chrom == var2.chrom and var1.start == var2.start and var1.ref == var2.ref and  alts1==alts2 and gt1 == gt2





