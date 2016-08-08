import os
import pysam
import util
import logging
from collections import defaultdict
import traceback as tb
import sys
from sim import bam_simulation

NO_VARS_FOUND_RESULT="No variants identified"
MATCH_RESULT="Variants matched"
NO_MATCH_RESULT="Variants did not match"
MATCH_WITH_EXTRA_RESULT= "Additional false variants present"
ERROR_RESULT="Error"

MISSING_QUAL = -1

ZYGOSITY_MATCH="Zygosity match"
ZYGOSITY_EXTRA_ALLELE="Extra allele"
ZYGOSITY_MISSING_ALLELE="Missing allele"
ZYGOSITY_MISSING_TWO_ALLELES="Missing two alleles!"

class VariantProcessor(object):

    def __init__(self, variant_callers, normalizers, comparators, output_reporter, conf):
        self.callers = variant_callers
        self.normalizers = normalizers
        self.comparators = comparators
        self.reporter = output_reporter
        self.conf = conf

    def create_input_vcf(self, raw_vars, ex_snp, gt_policy):
        ref_path = self.conf.get('main', 'ref_genome')
        variant_sets = create_variant_sets(raw_vars, ex_snp, gt_policy, pysam.FastaFile(ref_path))
        allvars = []
        for vset in variant_sets:
            allvars.extend(vset['vars'])
        variant_batch = sorted(allvars, cmp=util.variant_comp)
        orig_vcf = util.write_vcf(variant_batch, "test_input.vcf", self.conf)
        return orig_vcf, variant_sets

    def call_variants(self, bam, bed):
        """
        Run each of the variant callers on the given bam file, restricted to the given input regions
        """
        ref_path = self.conf.get('main', 'ref_genome')
        called_vars = {}
        for caller in self.callers:
            logging.info("Running variant caller " + caller)
            vars = self.callers[caller](bam, ref_path, bed, self.conf)
            called_vars[caller] = vars
        return called_vars



    def collect_bam_stats(self, bam, bed, orig_vcf):
        """
        For each bed region compute some bam-level stats and return them in a dict
        The key of the dict is the var_key of the matching original (input) variant
        """
        bam_stats = defaultdict(dict)
        for region in util.read_regions(bed):
            key = var_key(util.find_matching_var(orig_vcf, region))
            bam_stats[key] = bam_simulation.gen_bam_stats(bam, region)
        return bam_stats

    def collect_var_quals(self, caller_vars, bed, orig_vcf):
        """
        For each call and input variant, find the quality of the matching called variant, if there is one
        Return a dict[variant key][caller] for each input variant
        Missing variants (ref calls) are assigned a quality of MISSING_QUAL
        """
        var_quals = defaultdict(dict)
        for region in util.read_regions(bed):
            key = var_key(util.find_matching_var(orig_vcf, region))
            for caller in caller_vars:
                # Avoid pysam bug with empty vcfs
                if util.is_empty(caller_vars[caller]):
                    var_quals[key][caller] = MISSING_QUAL
                else:
                    with pysam.VariantFile(caller_vars[caller]) as cvars:
                        cvar = util.find_matching_var(cvars, region)
                        var_quals[key][caller] = find_qual(cvar)
        return var_quals

    def process_batch(self, vcf, batchname, gt_policy, ex_snp=None, keep_tmpdir=False, read_depth=250, reads=None):
        """
        Process the given batch of variants by creating a fake 'genome' with the variants, simulating reads from it,
         aligning the reads to make a bam file, then using different callers, variant normalizers, and variant
         comparison methods to generate results. The results are just written to a big text file, which needs to
         be parsed by a separate utility to generate anything readable.
        :param vcf: .vcf file containing variants to simulate
        :param conf: Configuration containing paths to all required binaries / executables / genomes, etc.
        :param homs: Boolean indicating whether variants should be simulated as hets or homs
        :return:
        """
        raw_vars = list(pysam.VariantFile(vcf))

        tmpdir_del_policy = util.TempDir.DELETE_NO_EXCEPTION
        if keep_tmpdir:
            tmpdir_del_policy = util.TempDir.NEVER_DELETE

        tmp_dirname = batchname + "-" + util.randstr()
        with util.TempDir(dirname=tmp_dirname, del_policy=tmpdir_del_policy):
            ref_path = self.conf.get('main', 'ref_genome')
            var_results = defaultdict(dict)

            orig_vcf, variant_sets = self.create_input_vcf(raw_vars, ex_snp, gt_policy)
            bed = util.vars_to_bed(variant_sets)
            if reads is None:
                reads = bam_simulation.gen_alt_fq(ref_path, variant_sets, read_depth)
            bam = bam_simulation.gen_alt_bam(ref_path, self.conf, reads)

            caller_variants = self.call_variants(bam, bed)
            bam_stats = self.collect_bam_stats(bam, bed, orig_vcf)
            var_quals = self.collect_var_quals(caller_variants, bed, orig_vcf)
            for normalizer_name, normalizer in self.normalizers.iteritems():
                logging.info("Running normalizer " + normalizer_name)
                normed_orig_vcf = normalizer(orig_vcf, self.conf)

                for caller in caller_variants:
                    normed_caller_vcf = normalizer(caller_variants[caller], self.conf)

                    for comparator_name, comparator in self.comparators.iteritems():
                        logging.info("Running comparator " + comparator_name + " (normalizer " + normalizer_name + ")")
                        all_results = comparator(normed_orig_vcf, normed_caller_vcf, None, self.conf)
                        single_results = split_results(all_results, bed)
                        for region, result in zip(util.read_regions(bed), single_results):
                            match_vars = util.find_matching_var(orig_vcf, region)
                            if len(match_vars) == 0:
                                raise ValueError('Unable to find original variant from region ' + str(region))
                            result = compare_single_var(result,
                                                        region,
                                                        normed_orig_vcf,
                                                        normed_caller_vcf,
                                                        comparator,
                                                        "/".join([str(i) for i in match_vars[0].samples[0]['GT']]),
                                                        self.conf)
                            key = var_key(match_vars)
                            if caller not in var_results[key]:
                                var_results[key][caller] = defaultdict(dict)
                            var_results[key][caller][normalizer_name][comparator_name] = result
            #Iterate over all results and write to standard output. We do this here instead of within the loops above
            #because it keeps results organized by variant, which makes them easier to look at
            self.reporter.write_output(var_results, var_quals, bam_stats)

    def compare_test_vcf(self, raw_orig_vcf, raw_test_vcf):

        raw_orig_vcf = os.path.abspath(raw_orig_vcf)
        raw_test_vcf = os.path.abspath(raw_test_vcf)
        orig_vars = list(pysam.VariantFile(raw_orig_vcf))
        tmp_dirname = raw_test_vcf.replace(".gz", "").replace(".vcf", "") + "-vcomp-" + util.randstr()
        with util.TempDir(dirname=tmp_dirname):
            orig_vcf = util.bgz_tabix(raw_orig_vcf, self.conf)

            test_vcf = util.remove_halfcalls(raw_test_vcf)
            test_vcf = util.bgz_tabix(test_vcf, self.conf)
            caller_name = test_vcf.replace(".gz", "").replace(".vcf", "")
            bed = util.vars_to_bed(orig_vars)
            var_results = defaultdict(dict)
            var_quals = self.collect_var_quals({caller_name: test_vcf}, bed, orig_vcf)
            bamstats = defaultdict(dict)


            for normalizer_name, normalizer in self.normalizers.iteritems():
                logging.info("Running normalizer " + normalizer_name)
                normed_orig_vcf = normalizer(orig_vcf, self.conf)
                normed_caller_vcf = normalizer(test_vcf, self.conf)

                for comparator_name, comparator in self.comparators.iteritems():
                    logging.info("Running comparator " + comparator_name + " (normalizer " + normalizer_name + ")")
                    all_results = comparator(normed_orig_vcf, normed_caller_vcf, None, self.conf)
                    single_results = split_results(all_results, bed)
                    for region, result in zip(util.read_regions(bed), single_results):
                        match_vars = util.find_matching_var(orig_vcf, region)
                        if len(match_vars) == 0:
                            raise ValueError('Unable to find original variant from region ' + str(region))
                        result = compare_single_var(result,
                                                    region,
                                                    normed_orig_vcf,
                                                    normed_caller_vcf,
                                                    comparator,
                                                    "/".join([str(i) for i in match_vars[0].samples[0]['GT']]),
                                                    self.conf)
                        key = var_key(match_vars)
                        if caller_name not in var_results[key]:
                            var_results[key][caller_name] = defaultdict(dict)
                        var_results[key][caller_name][normalizer_name][comparator_name] = result
                        bamstats[key] = {}
        # Iterate over all results and write to standard output. We do this here instead of within the loops above
        # because it keeps results organized by variant, which makes them easier to look at
        self.reporter.write_output(var_results, var_quals, bamstats)


def find_qual(vars):
    qual = MISSING_QUAL
    if vars is None or len(vars)==0:
        return qual

    var = vars[0]
    if var.qual is None:
        if 'GQ' in var.samples[0]:
            qual = var.samples[0]['GQ']
    else:
        qual = var.qual

    return qual

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
    if result_str == NO_MATCH_RESULT and (inputgt in util.ALL_HET_GTS or inputgt in util.ALL_HOMALT_GTS):
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
    :param vars: List of variants (not a VCF)
    :param ex_snp_info: Information describing additional SNPs to add to input variants
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


def var_key(vars):
    """
    Unique str for one or more variants
    """
    if type(vars) is not list:
        vars = [vars]
    return "/".join([" ".join(str(var).split()[0:5]) for var in vars])
