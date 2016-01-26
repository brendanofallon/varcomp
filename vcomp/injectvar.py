__author__ = 'bofallon'

import ConfigParser as cp
import argparse
import os
import sys
import traceback as tb
import util
import pysam
import time
from collections import namedtuple
import bam_simulation, callers, comparators, normalizers


NO_VARS_FOUND_RESULT="No variants identified"
MATCH_RESULT="Variants matched"
NO_MATCH_RESULT="Variants did not match"
MATCH_WITH_EXTRA_RESULT= "Additional false variants present"

ZYGOSITY_MATCH="Zygosity match"
ZYGOSITY_EXTRA_ALLELE="Extra allele"
ZYGOSITY_MISSING_ALLELE="Missing allele"
ZYGOSITY_MISSING_TWO_ALLELES="Missing two alleles!"

all_result_types = (MATCH_RESULT, NO_MATCH_RESULT, NO_VARS_FOUND_RESULT, MATCH_WITH_EXTRA_RESULT, ZYGOSITY_MISSING_ALLELE, ZYGOSITY_EXTRA_ALLELE)

def result_from_tuple(tup):
    """
    Comparators return a tuple consisting of unmatched_orig, matches, and unmatched_caller variants
    This examines the tuple and returns a short descriptive string of the result
    :param tup:
    :return:
    """
    unmatched_orig = tup[0]
    matches = tup[1]
    unmatched_caller = tup[2]

    #ONLY return a match if everything matches perfectly
    if len(unmatched_orig)==0 and len(unmatched_caller)==0 and len(matches)>0:
        return MATCH_RESULT

    if len(unmatched_orig)==0 and len(matches)>0 and len(unmatched_caller)>0:
        return MATCH_WITH_EXTRA_RESULT

    if len(unmatched_orig)>0 and len(matches)==0:
        if len(unmatched_caller)>0:
            return NO_MATCH_RESULT
        else:
            return NO_VARS_FOUND_RESULT

    return NO_MATCH_RESULT


def should_keep_dir(var_res):
    """
    Decide whether or not to flag the analysis dir for this variant for non-deletion (usually, we delete all tmp dirs)
    :param var_res: 3 layer dict of the form [caller][normalizer][comparator] containing results strings
    :param var: variant
    :return: Tuple of (boolean, suffix, comment), boolean indicates keep or not, suffix is applied to the tmpdir, comment is written to a file in the dir
    """
    #Want to flag following situations:
    #vgraph / vcfeval / happy disagree on anything
    #variant match with nonorm, but mismatch with vapleft or vt
    #variant mismatch with vapleft / raw but correctly matched with vgraph / vcfeval / etc

    keep = False
    suffix = None
    comment = None

    for caller in var_res:
        for norm in var_res[caller]:

            #vgraph_result = var_res[caller][norm]["vgraph"]
            vcfeval_result = var_res[caller][norm]["vcfeval"]
            happy_result = var_res[caller][norm]["happy"]

            if vcfeval_result != happy_result:
            # if vgraph_result != vcfeval_result or vcfeval_result != happy_result:
                keep = True
                suffix = "comp-mismatch"
                comment = "\n".join(["caller: " + caller, "norm:" + norm, "vgraph:" + vgraph_result, "vcfeval:" + vcfeval_result, "happy:"+ happy_result])

        nonorm_vcfeval_result = var_res[caller]["nonorm"]["vcfeval"]
        vap_vcfeval_result = var_res[caller]["vapleft"]["vcfeval"]
        vap_raw_result = var_res[caller]["vapleft"]["raw"]

        if nonorm_vcfeval_result == MATCH_RESULT and vap_vcfeval_result != MATCH_RESULT:
            keep = True
            suffix = "vapleft-broken"
            comment = "\n".join(["caller: " + caller, "nonorm / vcfeval:" + nonorm_vcfeval_result, "vapleft / vcfeval:" + vap_vcfeval_result])

        if vap_raw_result != MATCH_RESULT and nonorm_vcfeval_result == MATCH_RESULT:
            keep = True
            suffix = "vapleft-failed"
            comment = "\n".join(["caller: " + caller, "vapleft / raw:" + vap_raw_result, "nonorm / vcfeval:" + nonorm_vcfeval_result])

    return (keep, suffix, comment)

def compare_single_var(result, bedregion, orig_vars, caller_vars, comparator, inputgt, conf):
    """
    Determine a result string for the given result tuple. Not trivial since for NO_MATCH_RESULTS we need to
     determine if a simple genotype change will produce a match
    :param result: Result tuple produced by a comparator
    :param bedregion: Genomic region containing result
    :param orig_vars: 'original' (truth) variant set
    :param caller_vars: Variants produced by caller
    :param comparator: Comparator function
    :param inputgt: True variant genotype
    :param conf: Configuration
    :return:
    """
    result_str = result_from_tuple(result)
    remove_tmpdir = True
    tmpdir_suffix = ""
    if result_str == NO_MATCH_RESULT:
        gt_mod_vars = util.set_genotypes(caller_vars, inputgt, bedregion, conf)
        bedfile = util.region_to_bedfile(bedregion)
        gt_mod_result = comparator(orig_vars, gt_mod_vars, bedfile, conf)
        if result_from_tuple(gt_mod_result) == MATCH_RESULT:
            if inputgt in "0/1":
                result_str = ZYGOSITY_EXTRA_ALLELE
                remove_tmpdir = False
                tmpdir_suffix = "-extra-allele"
            else:
                result_str = ZYGOSITY_MISSING_ALLELE
                remove_tmpdir = False
                tmpdir_suffix = "-missing-allele"

    return result_str, remove_tmpdir, tmpdir_suffix

def split_results(allresults, bed):
    """
    all_results is the result of a call to a comparator, so it's a
    tuple of (unmatched_orig (FN), matches, unmatched_caller (FP)). This function
    breaks the variants all_results into the same structure, but with separate
    entries for each region in the bed file.
    :param allresults: Tuple containing results from a single comparator call
    :param bed: BED file to split regions by
    :return: List of tuples containing same data as allresults, but organized by bed region
    """
    reg_results = []
    for region in util.read_regions(bed):
        fns = [v for v in allresults[0] if v.chrom==region.chr and v.start >= region.start and v.start < region.end]
        matches = [v for v in allresults[1] if v[0].chrom==region.chr and v[0].start >= region.start and v[0].start < region.end]
        fps = [v for v in allresults[2] if v.chrom==region.chr and v.start >= region.start and v.start < region.end]
        reg_results.append( (fns, matches, fps) )
        #Sanity check...
        if len(fns)+len(matches)==0:
            raise ValueError('Uh oh, did not find any matching original vars for region!')

    return reg_results

def process_batch(variant_batch, conf, homs, keep_tmpdir=False):
    """
    Process the given variant by creating a fake 'genome' with the variant, simulating reads from it,
     aligning the reads to make a bam file, then using different callers, variant normalizers, and variant
     comparison methods to generate results. The results are just written to a big text file, which needs to
     be parsed by a separate utility to generate anything readable.
    :param variant_batch: pysam.Variant object to simulate
    :param conf: Configuration containing paths to all required binaries / executables / genomes, etc.
    :param homs: Boolean indicating whether variants should be simulated as hets or homs
    :return:
    """

    tmpdir = "tmp-working-" + util.randstr()
    try:
        os.mkdir(tmpdir)
    except:
        pass
    os.chdir(tmpdir)

    #Allow single variants to be passed in as well for easier debugging
    if type(variant_batch)!=list:
        variant_batch = [variant_batch]


    #The GT field to use in the true input VCF
    true_gt = "0/1"
    if homs:
        true_gt = "1/1"

    ttime = time.time()
    # print "Begin : " + str(time.time() - ttime)
    try:
        orig_vcf = util.write_vcf(variant_batch, "test_input.vcf", conf, true_gt)
        ref_path = conf.get('main', 'ref_genome')

        bed = util.vars_to_bed(variant_batch)
        reads = bam_simulation.gen_alt_fq(ref_path, variant_batch, homs, depth=250)
        bam = bam_simulation.gen_alt_bam(ref_path, conf, reads)

        # print "Done with bam gen : " + str(time.time() - ttime)

        var_results = {}
        variant_callers = callers.get_callers()
        variants = {}

        for caller in variant_callers:
            vars = variant_callers[caller](bam, ref_path, bed, conf)
            variants[caller] = vars


        # print "Done with calling : " + str(time.time() - ttime)



        remove_tmpdir = not keep_tmpdir
        tmpdir_suffix = ""
        for normalizer_name, normalizer in normalizers.get_normalizers().iteritems():
            normed_orig_vcf = normalizer(orig_vcf, conf)

            # print "Normalizing " + normalizer_name + " : " + str(time.time() - ttime)

            for caller in variants:
                normed_caller_vcf = normalizer(variants[caller], conf)

                for comparator_name, comparator in comparators.get_comparators().iteritems():
                    all_results = comparator(normed_orig_vcf, normed_caller_vcf, None, conf)

                    single_results = split_results(all_results, bed)
                    for region, result in zip(util.read_regions(bed), single_results):
                        result, remove_tmpdir, tmpdir_suffix = compare_single_var(result, region, normed_orig_vcf, normed_caller_vcf, comparator, true_gt, conf)
                        match_vars = util.find_matching_var( pysam.VariantFile(orig_vcf), region)
                        if len(match_vars)!=1:
                            raise ValueError('Unable to find original variant from region!')
                        match_var = " ".join(str(match_vars[0]).split()[0:5])
                        if match_var not in var_results:
                            var_results[match_var] = {}
                        if caller not in var_results[match_var]:
                            var_results[match_var][caller] = {}
                        if normalizer_name not in var_results[match_var][caller]:
                            var_results[match_var][caller][normalizer_name] = {}

                        var_results[match_var][caller][normalizer_name][comparator_name] = result
                        print "Result for " + " ".join(str(match_var).split()[0:5]) + ": " + caller + " / " + normalizer_name + " / " + comparator_name + ": " + result
                    # print "All comps for " + caller + " / " + normalizer_name + " : " + str(time.time() - ttime)

        for origvar in var_results.keys():
            keep, suffix, comment = should_keep_dir(var_results[origvar])
            if keep:
                remove_tmpdir = False
                tmpdir_suffix = suffix
                with open("flag.info.txt", "a") as fh:
                    fh.write(str(comment) + "\n")
                with open("var.info.txt", "w") as fh:
                    fh.write("Causative variant: " + str(origvar) + "\n")

    except Exception as ex:
        print "Error processing variant " + str(variant_batch) + " : " + str(ex)
        tb.print_exc(file=sys.stdout)
        remove_tmpdir = False
        tmpdir_suffix = "error"
        try:
            with open("exception.info.txt", "a") as fh:
                fh.write(str(ex) + "\n")
        except:
            #we tried...
            pass

    os.chdir("..")
    if remove_tmpdir:
        os.system("rm -rf " + tmpdir)
    else:
        dirname = "tmpdir" + tmpdir_suffix
        os.system("mv " + tmpdir + " " + dirname)

def process_vcf(input_vcf, homs, conf, keep_tmpdir=False):
    """
    Iterate over entire vcf file, processing each variant individually and collecting results
    :param input_vcf:
    :param conf:
    :return:
    """
    input_vcf = pysam.VariantFile(input_vcf)
    process_batch(list(input_vcf), conf, homs, keep_tmpdir)
    # for input_var in input_vcf:
    #     process_variant(input_var, conf, homs, keep_tmpdir)


if __name__=="__main__":
    parser = argparse.ArgumentParser("Inject, simulate, call, compare")
    parser.add_argument("-c", "--conf", help="Path to configuration file", default="./comp.conf")
    parser.add_argument("-v", "--vcf", help="Input vcf file")
    parser.add_argument("-k", "--keep", help="Dont delete temporary directories", action='store_true')
    parser.add_argument("--het", help="Run all variants as hets (default false, run everything as homs)", action='store_true')
    args = parser.parse_args()

    conf = cp.SafeConfigParser()
    conf.read(args.conf)
    process_vcf(args.vcf, not args.het, conf, args.keep)
