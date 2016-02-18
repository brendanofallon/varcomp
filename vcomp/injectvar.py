__author__ = 'bofallon'

import ConfigParser as cp
import argparse
import os
import sys
import traceback as tb
import util
import pysam
import logging
import random
import json
import bam_simulation, callers, comparators, normalizers
from collections import namedtuple, defaultdict

NO_VARS_FOUND_RESULT="No variants identified"
MATCH_RESULT="Variants matched"
NO_MATCH_RESULT="Variants did not match"
MATCH_WITH_EXTRA_RESULT= "Additional false variants present"
ERROR_RESULT="Error"

ZYGOSITY_MATCH="Zygosity match"
ZYGOSITY_EXTRA_ALLELE="Extra allele"
ZYGOSITY_MISSING_ALLELE="Missing allele"
ZYGOSITY_MISSING_TWO_ALLELES="Missing two alleles!"

all_result_types = (MATCH_RESULT, NO_MATCH_RESULT, NO_VARS_FOUND_RESULT, MATCH_WITH_EXTRA_RESULT, ZYGOSITY_MISSING_ALLELE, ZYGOSITY_EXTRA_ALLELE, ERROR_RESULT)

ExSNPInfo = namedtuple('ExSNPInfo', ['policy', 'dist'])

def result_from_tuple(tup):
    """
    Comparators return a tuple consisting of unmatched_orig, matches, and unmatched_caller variants
    This examines the tuple and returns a short descriptive string of the result
    :param tup:
    :return: Short, human readable string describing result
    """
    unmatched_orig = tup[0]
    matches = tup[1]
    unmatched_caller = tup[2]

    for v in unmatched_orig:
        if type(v) is util.ErrorVariant:
            return ERROR_RESULT
    for v in unmatched_caller:
        if type(v) is util.ErrorVariant:
            return ERROR_RESULT

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
    if result_str == NO_MATCH_RESULT  and (inputgt in util.ALL_HET_GTS or inputgt in util.ALL_HOMALT_GTS):
        try:
            gt_mod_vars = util.set_genotypes(caller_vars, inputgt, bedregion, conf)
            bedfile = util.region_to_bedfile(bedregion)
            gt_mod_result = comparator(orig_vars, gt_mod_vars, bedfile, conf)
            if result_from_tuple(gt_mod_result) == MATCH_RESULT:
                if inputgt in util.ALL_HET_GTS:
                    result_str = ZYGOSITY_EXTRA_ALLELE
                else:
                    result_str = ZYGOSITY_MISSING_ALLELE
        except util.GTModException as ex:
            logging.warning('Exception while attempting to modify GTs: ' + str(ex))

    return result_str

def split_results(allresults, bed):
    """
    allresults is the result of a call to a comparator, so it's a
    tuple of (unmatched_orig (FN), matches, unmatched_caller (FP)). This function
    breaks the allresults into a list with separate
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

    return reg_results

def create_variant_sets(vars, ex_snp_info, default_policy, ref_genome):
    """
    Create a list of variant 'sets', where each set is a list of possibly-phased variants to add to two
    alternate reference genomes. If ex_snp_info isn't None, extra SNPs are added to each input variant, otherwise,
    every variant ends up in its own unique set.
    Default-pol
    :param vars:
    :param ex_snp_info:
    :param default_policy: Genotype policy for original variants - either None (read GT from sample field), ALL_HETS, or ALL_HOMS
    :param ref_genome:
    :return:
    """
    sets = []
    if ex_snp_info is not None and ex_snp_info.dist >=0:
        raise ValueError('Positive snp distances are currently not supported')

    default_gt = None
    if default_policy ==  bam_simulation.ALL_HOMS:
        default_gt = "1|1"
    if default_policy ==  bam_simulation.ALL_HETS:
        default_gt = "0|1"

    for var in vars:
        vset = {}
        if ex_snp_info is None:
            vset['policy'] = default_policy
            vset['vars'] = [ util.pysamVar_to_Variant(var, default_gt) ]
        else:
            vset['policy'] = ex_snp_info.policy

            if ex_snp_info.policy == bam_simulation.ALL_HOMS:
                snp_gt = "1/1"
                if default_gt is None:
                    var_gt = "1/1"
                else:
                    var_gt = default_gt
            elif ex_snp_info.policy == bam_simulation.CIS:
                snp_gt = "0|1"
                if default_gt is None:
                    var_gt = "0|1"
                else:
                    var_gt = default_gt
            elif ex_snp_info.policy == bam_simulation.TRANS:
                snp_gt = "0|1"
                if default_gt is None or default_policy == bam_simulation.ALL_HETS:
                    var_gt = "1|0"
                else:
                    var_gt = default_gt
            newsnp = util.gen_snp(var.chrom, var.start + ex_snp_info.dist, snp_gt, ref_genome)
            vset['vars'] = [newsnp]
            vset['vars'].append(util.pysamVar_to_Variant(var, var_gt))
        sets.append(vset)

    return sets

def emit_batch_output(results, bamstats, output):
    """
    Write output for a batch of input variants (with individual entries for each caller/normalizer/comparator
      combination) to the given output handle.
    :param results: Four-level deep dict containing results for [input variant string][caller][normalizer][comparator]
    :param bamstats: Dictionary containing statistics for bam file
    :param output: File-like object to which formatted output will be written
    """
    for var, vresults in results.iteritems():
        json.dump({
            "variant":var,
            "bamstats": bamstats[var],
            "results": vresults
        }, output)
        output.write("\n")


def process_batch(variant_batch, batchname, conf, gt_policy, ex_snp=None, output=sys.stdout, keep_tmpdir=False, disable_flagging=False, read_depth=250):
    """
    Process the given batch of variants by creating a fake 'genome' with the variants, simulating reads from it,
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

    ref_path = conf.get('main', 'ref_genome')
    bam_stats = defaultdict(dict)
    remove_tmpdir = not keep_tmpdir
    try:

        variant_sets = create_variant_sets(variant_batch, ex_snp, gt_policy, pysam.FastaFile(ref_path))
        allvars = []
        for vset in variant_sets:
            allvars.extend( vset['vars'] )
        variant_batch = sorted(allvars, cmp=util.variant_comp)
        orig_vcf = util.write_vcf(variant_batch, "test_input.vcf", conf)

        bed = util.vars_to_bed(variant_sets)
        reads = bam_simulation.gen_alt_fq(ref_path, variant_sets, read_depth)
        bam = bam_simulation.gen_alt_bam(ref_path, conf, reads)

        var_results = defaultdict(dict)
        variant_callers = callers.get_callers()
        variants = {}

        #Call all variants using each variant caller, store in variants dict
        for caller in variant_callers:
            logging.info("Running variant caller " + caller)
            vars = variant_callers[caller](bam, ref_path, bed, conf)
            variants[caller] = vars

        #Compute bam statistics separately for each region, and store them in a dictionary indexed
        #by the same key used to store individual varian results
        for region in util.read_regions(bed):
            match_vars = util.find_matching_var( pysam.VariantFile(orig_vcf), region)
            match_var = "/".join([" ".join(str(mvar).split()[0:5]) for mvar in match_vars])
            bam_stats[match_var] = bam_simulation.gen_bam_stats(bam, region)


        for normalizer_name, normalizer in normalizers.get_normalizers().iteritems():
            logging.info("Running normalizer " + normalizer_name)
            normed_orig_vcf = normalizer(orig_vcf, conf)

            for caller in variants:
                normed_caller_vcf = normalizer(variants[caller], conf)

                for comparator_name, comparator in comparators.get_comparators().iteritems():
                    all_results = comparator(normed_orig_vcf, normed_caller_vcf, None, conf)
                    single_results = split_results(all_results, bed)
                    logging.info("Running comparator " + comparator_name)
                    for region, result in zip(util.read_regions(bed), single_results):
                        match_vars = util.find_matching_var( pysam.VariantFile(orig_vcf), region)
                        if len(match_vars)==0:
                            raise ValueError('Unable to find original variant from region!')

                        result = compare_single_var(result, region, normed_orig_vcf, normed_caller_vcf, comparator, "/".join([str(i) for i in match_vars[0].samples[0]['GT']]), conf)

                        match_var = "/".join([" ".join(str(mvar).split()[0:5]) for mvar in match_vars])

                        if caller not in var_results[match_var]:
                            var_results[match_var][caller] = defaultdict(dict)

                        var_results[match_var][caller][normalizer_name][comparator_name] = result



        #Iterate over all results and write to standard output. We do this here instead of within the loops above
        #because it keeps results organized by variant, which makes them easier to look at
        emit_batch_output(var_results, bam_stats, output)

    except Exception as ex:
        logging.error("Error processing variant batch " + batchname + " : " + str(ex))
        tb.print_exc(file=sys.stderr)
        remove_tmpdir = False
        try:
            with open("exception.info.txt", "a") as fh:
                fh.write(str(ex) + "\n")
        except:
            #we tried...
            pass

    os.chdir("..")
    if remove_tmpdir:
        os.system("rm -rf " + tmpdir)


def process_vcf(vcf, gt_default, conf, output, snp_info=None, single_batch=False, keep_tmpdir=False, disable_flagging=False, read_depth=250):
    """
    Perform analyses for each variant in the VCF file.
    :param input_vcf:
    :param single_batch: Assume all variants in VCF are part of one batch and process them all simultaneously
    :param keep_tmpdir: Preserve tmpdirs created (otherwise delete them, unless they are flagged)
    :param conf:
    """

    input_vars = pysam.VariantFile(vcf)
    logging.info("Processing variants in file " + vcf)
    if single_batch:
        logging.info("Processing all variants as one batch")
        process_batch(list(input_vars), vcf.replace(".vcf", "-tmpfiles"), conf, gt_default, ex_snp=snp_info, output=output, keep_tmpdir=keep_tmpdir, read_depth=read_depth, disable_flagging=disable_flagging)
    else:
        for batchnum, batch in enumerate(util.batch_variants(input_vars, max_batch_size=1000, min_safe_dist=2000)):
            logging.info("Processing batch #" + str(batchnum) + " containing " + str(len(batch)) + " variants")
            process_batch(batch, vcf.replace(".vcf", "-tmpfiles"), conf, gt_default, ex_snp=snp_info, output=output, keep_tmpdir=keep_tmpdir, read_depth=read_depth, disable_flagging=disable_flagging)



if __name__=="__main__":
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    parser = argparse.ArgumentParser("Inject, simulate, call, compare")
    parser.add_argument("-c", "--conf", help="Path to configuration file", default="./comp.conf")
    parser.add_argument("-v", "--vcf", help="Input vcf file(s)", nargs="+")
    parser.add_argument("-k", "--keep", help="Dont delete temporary directories", action='store_true')
    parser.add_argument("-b", "--batch", help="Treat each input VCF file as a single batch (default False)", action='store_true')
    parser.add_argument("-s", "--seed", help="Random seed", default=None)
    parser.add_argument("-o", "--output", help="Output destination", default=sys.stdout)
    parser.add_argument("-r", "--readdepth", help="Number of reads to generate per variant", default=250, type=int)
    parser.add_argument("-d", "--dontflag", help="Disable temp dir flagging", action='store_true')
    parser.add_argument("-a", "--addsnp", help="Add a SNP upstream of each variant", action='store_true')
    parser.add_argument("-t", "--trans", help="If SNP is added, add it in trans (default cis)", action='store_true')
    parser.add_argument("--snphom", help="Added SNPs are homozygous (default het)", action='store_true')
    parser.add_argument("--het", help="Force all simulated variants to be hets", action='store_true')
    parser.add_argument("--hom", help="Force all simulated variants to be homozygotes", action='store_true')
    args = parser.parse_args()

    conf = cp.SafeConfigParser()
    conf.read(args.conf)

    if type(args.output) is str:
        args.output = open(args.output, "w")

    if args.seed is not None:
        random.seed(args.seed)

    if args.het and args.hom:
        raise ValueError('Specify just one of --het or --hom')

    gt_default = None
    if args.het:
        gt_default = bam_simulation.ALL_HETS
    if args.hom:
        gt_default = bam_simulation.ALL_HOMS

    snp_inf = None
    if args.addsnp:
        if args.trans:
            snp_policy = bam_simulation.TRANS
        else:
            snp_policy = bam_simulation.CIS
        if args.snphom:
            snp_policy = bam_simulation.ALL_HOMS
        snp_inf = ExSNPInfo(policy=snp_policy, dist=-4)


    for vcf in args.vcf:
        process_vcf(vcf, gt_default, conf, args.output, snp_info=snp_inf, single_batch=args.batch, keep_tmpdir=args.keep, read_depth=args.readdepth, disable_flagging=args.dontflag)

    try:
        args.output.close()
    except:
        pass
