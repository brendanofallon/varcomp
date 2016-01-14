__author__ = 'bofallon'

import ConfigParser as cp
import argparse
import os
import sys
import time
import traceback as tb
import util
import pysam

import bam_simulation, callers, comparators, normalizers


NO_VARS_FOUND_RESULT="No variants identified"
MATCH_RESULT="Variants matched"
NO_MATCH_RESULT="Variants did not match"

ZYGOSITY_MATCH="Zygosity match"
ZYGOSITY_EXTRA_ALLELE="Extra allele"
ZYGOSITY_MISSING_ALLELE="Missing allele"
ZYGOSITY_MISSING_TWO_ALLELES="Missing two alleles!"

all_result_types = (MATCH_RESULT, NO_MATCH_RESULT, NO_VARS_FOUND_RESULT, ZYGOSITY_MISSING_ALLELE, ZYGOSITY_EXTRA_ALLELE)

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
    if len(unmatched_orig)>0 and len(matches)==0:
        if len(unmatched_caller)>0:
            return NO_MATCH_RESULT
        else:
            return NO_VARS_FOUND_RESULT
    if len(unmatched_orig)==0 and len(matches)>0:
        return MATCH_RESULT
    return NO_MATCH_RESULT



def process_variant(variant, conf, homs, keep_tmpdir=False):
    """
    Process the given variant by creating a fake 'genome' with the variant, simulating reads from it,
     aligning the reads to make a bam file, then using different callers, variant normalizers, and variant
     comparison methods to generate results. The results are just written to a big text file, which needs to
     be parsed by a separate utility to generate anything readable.
    :param variant: pysam.Variant object to simulate
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

    vcf_gt = "0/1"
    if homs:
        vcf_gt = "1/1"
    try:
        orig_vcf = bam_simulation.write_vcf(variant, "test_input.vcf", conf, vcf_gt)
        ref_path = conf.get('main', 'ref_genome')

        bed = callers.vars_to_bed([variant])
        bam = bam_simulation.gen_alt_bam(ref_path, [variant], conf, homs)

        variant_callers = callers.get_callers()
        variants = {}
        for caller in variant_callers:
            vars = variant_callers[caller](bam, ref_path, bed, conf)
            variants[caller] = vars


        remove_tmpdir = not keep_tmpdir
        tmpdir_suffix = ""
        for normalizer_name, normalizer in normalizers.get_normalizers().iteritems():
            normed_orig_vars = normalizer(orig_vcf, conf)

            for caller in variants:
                normed_caller_vars = normalizer(variants[caller], conf)

                vcfeval_result = None #Special case code here for flagging results where vgraph and vcfeval differ
                vgraph_result = None

                for comparator_name, comparator in comparators.get_comparators().iteritems():
                    result = comparator(normed_orig_vars, normed_caller_vars, bed, conf)
                    result_str = result_from_tuple(result)
                    if result_str == NO_MATCH_RESULT:
                        gt_mod_vars = util.set_genotypes(normed_caller_vars, vcf_gt, conf)
                        gt_mod_result = comparator(normed_orig_vars, gt_mod_vars, bed, conf)
                        if result_from_tuple(gt_mod_result) == MATCH_RESULT:
                            if vcf_gt in "0/1":
                                result_str = ZYGOSITY_EXTRA_ALLELE
                                remove_tmpdir = False
                                tmpdir_suffix = caller + "-" + normalizer_name + "-extra-allele"
                            else:
                                result_str = ZYGOSITY_MISSING_ALLELE
                                remove_tmpdir = False
                                tmpdir_suffix = caller + "-" + normalizer_name + "-missing-allele"
                    if comparator_name == "vgraph":
                        vgraph_result = result_str
                    if comparator_name == "vcfeval":
                        vcfeval_result = result_str
                    if vcfeval_result is not None and vgraph_result is not None and vcfeval_result != vgraph_result:
                        remove_tmpdir = False
                        tmpdir_suffix = caller + "-" + normalizer_name + "-comp-conflict"
                        with open("conflict.info.txt", "a") as fh:
                            fh.write(caller + "\t" + normalizer_name + "\tvgraph: " + vgraph_result + "\tvcfeval: " + vcfeval_result + "\n")

                    print "Result for " + " ".join( str(variant).split()[0:5]) + ": " + caller + " / " + normalizer_name + " / " + comparator_name + ": " + result_str
    except Exception as ex:
        print "Error processing variant " + str(variant) + " : " + str(ex)
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
        dirname = "tmpdir-chr" + variant.chrom + "-" + str(variant.start) + "-" + tmpdir_suffix
        os.system("mv " + tmpdir + " " + dirname)

def process_vcf(input_vcf, homs, conf, keep_tmpdir=False):
    """
    Iterate over entire vcf file, processing each variant individually and collecting results
    :param input_vcf:
    :param conf:
    :return:
    """
    input_vcf = pysam.VariantFile(input_vcf)
    for input_var in input_vcf:
        process_variant(input_var, conf, homs, keep_tmpdir)


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
