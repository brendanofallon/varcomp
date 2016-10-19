import os
import time
import subprocess
import gzip
import pysam

import vcomp.util


ALLELE_MATCH="Alleles matched"
ALLELE_MISMATCH="Alleles did not match"
ALLELE_EXTRA="Additional variants identified"


def get_comparators():
    return {
        #"raw": compare_raw,
        #"vgraph": compare_vgraph,
        "vcfeval": compare_vcfeval,
        #"happy": compare_happy
    }

def compare_raw(orig_vcf, caller_vcf, bed, conf):
    """
    Compare variants with absolutely no normalization.  This still requires
    a little work since there are no guarantees about variants being in the
    right order or the positions of false pos / false neg vars
    :return:

    FIXME: Should update to use vgraph's parallel iterator and only apply
           the trivial variant matcher.  This will be much more elegant.  A
           quicker/dirtier fix would be to hash the called vars, rather than
           relying on a quadratic search.
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
    # Slightly tricky - must make sure executed python script is executed in
    # a workable virtualenv.  For now we execute it in whatever virtualenv
    # this script is being processed in.

    # KBJ: sub-shells should automatically inherit the current processes's
    # environment.  Not sure why the virtualenv is tricky unless vgraph is
    # submitted via SGE/OGE qsub.

    #First, test to see if there are any vars present...
    #orig_vars = read_all_vars(orig_vcf, bed)
    caller_vars = read_all_vars(caller_vcf, bed)
    if not caller_vars:
        return (read_all_vars(orig_vcf, bed), [], caller_vars)

    orig_out = "vgout-orig." + vcomp.util.randstr() + ".vcf"
    caller_out = "vgout-caller." + vcomp.util.randstr() + ".vcf"
    bedcmd = ("--include-regions " + bed) if bed else ''

    cmd = '{vgraph} --out1 {orig_out} --out2 {caller_out} --reference {genome} {bedcmd} {orig_vcf} {caller_vcf} > /dev/null'
    cmd = cmd.format(vgraph=conf.get('main', 'vgraph'), orig_out=orig_out, caller_out=caller_out,
                     genome=conf.get('main', 'ref_genome'), bedcmd=bedcmd,
                     orig_vcf=orig_vcf, caller_vcf=caller_vcf)

    time.sleep(1)
    subprocess.check_call(cmd, shell=True)

    unmatched_orig = []
    matches = []
    unmatched_caller = []

    for ovar in pysam.VariantFile(orig_out):
        bd = ovar.samples[0]['BD']
        #bk = ovar.samples[0]['BK']

        if bd == 'X':
            unmatched_orig.append(ovar)
        if bd == 'N':
            unmatched_orig.append(vcomp.util.ErrorVariant(chrom=ovar.chrom, start=ovar.start, msg="vgraph error code N"))

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
    if not caller_vars:
        return (read_all_vars(orig_vcf, bed), [], caller_vars)

    output_dir = "vcfeval-output" + vcomp.util.randstr()
    cmd = "java -Djava.io.tmpdir=. -Xmx4g -jar " + conf.get('main', 'rtg_jar') + " vcfeval -t " + conf.get('main', 'rtg_ref_sdf') + " --all-records -o " + output_dir + " -b " + orig_vcf + " -c " + caller_vcf
    if bed is not None:
        cmd = cmd + " --bed-regions " + bed
    subprocess.check_output(cmd, shell=True)
    # orig_vars = read_all_vars(orig_vcf, bed)
    tp_vars = read_all_vars(output_dir + "/tp.vcf.gz")
    fp_vars = read_all_vars(output_dir + "/fp.vcf.gz")
    fn_vars = read_all_vars(output_dir + "/fn.vcf.gz")
    return (fn_vars, zip(tp_vars, tp_vars), fp_vars)


def compare_happy(orig_vcf, caller_vcf, bed, conf):
    orig_vars = read_all_vars(orig_vcf, bed)
    caller_vars = read_all_vars(caller_vcf, bed)
    if not caller_vars:
        return (read_all_vars(orig_vcf, bed), [], caller_vars)
    output_prefix = "happyoutput-" + vcomp.util.randstr(6)

    all_chrs = set(v.chrom for v in orig_vars)
    all_chrs.update(v.chrom for v in caller_vars)
    all_chrs = list(all_chrs)

    bedarg = ""
    if bed is not None:
        bedarg = " -T " + bed
    cmd = conf.get('main', 'happy') + " " + orig_vcf + " " + caller_vcf + " " + bedarg + " -o " + output_prefix + " --scratch-prefix=. --include-nonpass -r " + conf.get('main', 'ref_genome') + " -l " + ",".join(all_chrs) + " --no-fixchr-truth --no-fixchr-query -V"
    subprocess.check_call(cmd, shell=True, stdout=open('/dev/null'))

    orig_unmatched = []
    matches = []
    caller_unmatched = []

    for var in pysam.VariantFile(output_prefix + ".vcf.gz"):
        # FIXME: use Variant API
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


def vcf_has_no_vars(vcf):
    """ Returns true if the vcf does not contain any variants. This does not use pysam. """
    if vcf.endswith(".gz"):
        fh = gzip.open(vcf)
    else:
        fh = open(vcf)

    for line in fh:
        if len(line)>1 and line[0] != '#':
            fh.close()
            return True

    fh.close()
    return False

def read_all_vars(vcf, bed=None):
    """
    Try to read all the variants from the given vcf file into a list. If there's an error
    reading the vcf, return an empty list
    :param vcf:VCF to read variants from
    :param bed:If not none, only read variants in these regions
    :return:
    """
    if vcf_has_no_vars(vcf):
        return []
    if bed is None:
        vars = list(pysam.VariantFile(vcf))
    else:
        # FIXME: Use bed reader
        vfh = pysam.VariantFile(vcf)
        vars = []
        for line in open(bed, "r"):
            line = line.strip()
            if not line or line[0]=='#':
                continue
            toks = line.split()
            chrom = toks[0]
            start = int(toks[1])
            stop = int(toks[2])
            vars.extend(list(vfh.fetch(chrom, start, stop)))
    return vars


def test_var_equiv(var1, var2):
    """
    Compare two VariantRecords for genotype equality
    :param var1:
    :param var2:
    :return: True if both variant records have identical genotype calls
             (treated as unphased).  No calls and partial calls always return False.
    """
    a1 = sorted(var1.samples[0].alleles) if 'GT' in var1 else [None]
    a2 = sorted(var2.samples[0].alleles) if 'GT' in var2 else [None]
    if None in (a1+a2):
        return False
    v1 = [var1.chrom, var1.start, var1.stop] + a1
    v2 = [var2.chrom, var2.start, var2.stop] + a2
    return v1 == v2
