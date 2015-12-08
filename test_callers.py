
import argparse
import ConfigParser as cp
import callers
import pysam
import os
import bam_simulation
import time
import comparators
import traceback as tb

def process_variant(variant, results, conf):
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

    if len(variant.alts) != 1:
        raise ValueError("Only one alt per variant now")
    orig_vcf = bam_simulation.write_vcf(variant, "test_input.vcf", conf)
    ref_path = conf.get('main', 'ref_genome')

    bam = bam_simulation.gen_alt_bam(ref_path, variant, conf)

    variant_callers = callers.get_callers()
    variants = {}
    for caller in variant_callers:
        vars = variant_callers[caller](bam, ref_path, variant.chrom, variant.start-250, variant.start+250, conf)
        variants[caller] = vars

    vgraph_comp = comparators.get_comparators()['vgraph']
    for caller, vars in variants.iteritems():
        result = vgraph_comp(orig_vcf, vars, conf)
        print "Result for " + " ".join( str(variant).split()[0:5]) + ": " + caller + ": " + result
        results[caller][result] += 1

    os.chdir("..")
    os.system("rm -rf " + tmpdir)


def process_vcf(input_vcf, conf):
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
        all_results[caller_name] = {
                comparators.NO_VARS_FOUND_RESULT: 0,
                comparators.NO_MATCH_RESULT: 0,
                comparators.MATCH_RESULT: 0,
                comparators.PARTIAL_MATCH: 0
            }


    for input_var in input_vcf:
        try:
            process_variant(input_var, all_results, conf)
        except Exception as ex:
            print "Error processing variant " + str(input_var) + " : " + str(ex)
            tb.print_exc()

    for caller in all_results:
        print "Caller: " + caller
        for result, count in all_results[caller].iteritems():
            print "\t\t" + result + "\t:\t" + str(count)





if __name__=="__main__":
    parser = argparse.ArgumentParser("Test variant callers")
    parser.add_argument("-c", "--conf", help="Path to configuration file", default="./comp.conf")
    parser.add_argument("-v", "--vcf", help="Input vcf file")
    args = parser.parse_args()

    conf = cp.SafeConfigParser()
    conf.read(args.conf)

    process_vcf(args.vcf, conf)