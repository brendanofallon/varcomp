__author__ = 'bofallon'

import pysam
import os
import sys
import ConfigParser as cp
import argparse
import callers
import time
import bam_simulation
import comparators
import traceback as tb


def process_variant(variant, results, conf, homs):
    """
    Process the given variant, update results dict
    :param variant:
    :param results_by_caller:
    :param results_by_method:
    :param conf:
    :return:
    """

    tmpdir = "tmp-working" + str(time.time())[-7:].replace(".", "")
    try:
        os.mkdir(tmpdir)
    except:
        pass
    os.chdir(tmpdir)

    vcf_gt = "0/1"
    if homs:
        vcf_gt = "1/1"
    orig_vcf = bam_simulation.write_vcf(variant, "test_input.vcf", conf, vcf_gt)
    ref_path = conf.get('main', 'ref_genome')

    bed = callers.vars_to_bed([variant])
    bam = bam_simulation.gen_alt_bam(ref_path, [variant], conf, homs)

    variant_callers = callers.get_callers()
    variants = {}
    for caller in variant_callers:
        vars = variant_callers[caller](bam, ref_path, bed, conf)
        variants[caller] = vars

    for method_name, comp in comparators.get_comparators().iteritems():
        for caller, vars in variants.iteritems():
            result = comp(orig_vcf, vars, bed, conf)
            result_str = result.short()
            if result_str == comparators.NO_MATCH_RESULT or result_str == comparators.PARTIAL_MATCH:
                gt_mod_vars = comparators.set_genotypes(vars, vcf_gt, conf)
                gt_mod_result = comp(orig_vcf, gt_mod_vars, bed, conf)
                if gt_mod_result.short() == comparators.MATCH_RESULT:
                    if vcf_gt in "0/1":
                        result_str = comparators.ZYGOSITY_EXTRA_ALLELE
                    else:
                        result_str = comparators.ZYGOSITY_MISSING_ALLELE
            print "Result for " + " ".join( str(variant).split()[0:5]) + ": " + caller + " / " + method_name + ": " + result_str
            results[caller][method_name][result_str] += 1

    os.chdir("..")
    os.system("rm -rf " + tmpdir)

def process_vcf(input_vcf, homs, conf):
    """
    Iterate over entire vcf file, processing each variant individually and collecting results
    :param input_vcf:
    :param conf:
    :return:
    """
    input_vcf = pysam.VariantFile(input_vcf)

    #Initialize results structure
    all_results = {}
    for caller_name in callers.get_callers():
        all_results[caller_name] = {}
        for method_name in comparators.get_comparators():
            all_results[caller_name][method_name] = {
                comparators.NO_VARS_FOUND_RESULT: 0,
                comparators.NO_MATCH_RESULT: 0,
                comparators.MATCH_RESULT: 0,
                comparators.ZYGOSITY_MISSING_ALLELE: 0,
                comparators.ZYGOSITY_MISSING_TWO_ALLELES: 0,
                comparators.ZYGOSITY_EXTRA_ALLELE: 0,
            }


    for input_var in input_vcf:
        try:
            process_variant(input_var, all_results, conf, homs)
        except Exception as ex:
            print "Error processing variant " + str(input_var) + " : " + str(ex)
            tb.print_exc(file=sys.stdout)

    for caller in all_results:
        print "Caller: " + caller
        for method in all_results[caller]:
            print "\t" + method
            for result, count in all_results[caller][method].iteritems():
                print "\t\t" + result + "\t:\t" + str(count)


if __name__=="__main__":
    parser = argparse.ArgumentParser("Inject, simulate, call, compare")
    parser.add_argument("-c", "--conf", help="Path to configuration file", default="./comp.conf")
    parser.add_argument("-v", "--vcf", help="Input vcf file")
    parser.add_argument("--het", help="Run all variants as hets (default false, run everything as homs)", action='store_true')
    args = parser.parse_args()

    conf = cp.SafeConfigParser()
    conf.read(args.conf)
    process_vcf(args.vcf, not args.het, conf)
