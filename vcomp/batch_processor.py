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

    def __init__(self, variant_callers, normalizers, comparators, output_reporter):
        self.callers = variant_callers
        self.normalizers = normalizers
        self.comparators = comparators
        # self.read_simulator = read_simulator
        self.reporter = output_reporter


    def process_batch(self, vcf, batchname, conf, gt_policy, ex_snp=None, keep_tmpdir=False, read_depth=250, reads=None):
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
        variant_batch = list(pysam.VariantFile(vcf))
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
                allvars.extend(vset['vars'])
            variant_batch = sorted(allvars, cmp=util.variant_comp)
            orig_vcf = util.write_vcf(variant_batch, "test_input.vcf", conf)

            bed = util.vars_to_bed(variant_sets)
            if reads is None:
                reads = bam_simulation.gen_alt_fq(ref_path, variant_sets, read_depth)

            bam = bam_simulation.gen_alt_bam(ref_path, conf, reads)

            var_results = defaultdict(dict)
            variants = {}
            var_quals = defaultdict(dict)

            #Call all variants using each variant caller, store in variants dict
            for caller in self.callers:
                logging.info("Running variant caller " + caller)
                vars = self.callers[caller](bam, ref_path, bed, conf)
                variants[caller] = vars

            #Compute bam statistics separately for each region, and store them in a dictionary indexed
            #by the same key used to store individual varian results

            # origvars = pysam.VariantFile(orig_vcf, 'rb')
            for region in util.read_regions(bed):
                match_vars = util.find_matching_var( orig_vcf, region)
                match_var = "/".join([" ".join(str(mvar).split()[0:5]) for mvar in match_vars])
                bam_stats[match_var] = bam_simulation.gen_bam_stats(bam, region)
                for caller in self.callers:
                    #Avoid pysam bug with empty vcfs
                    if util.is_empty(variants[caller]):
                        var_quals[match_var][caller] = MISSING_QUAL
                    else:
                        with pysam.VariantFile(variants[caller]) as cvars:
                            cvar = util.find_matching_var(cvars, region)
                            var_quals[match_var][caller] = find_qual(cvar)


            for normalizer_name, normalizer in self.normalizers.iteritems():
                logging.info("Running normalizer " + normalizer_name)
                normed_orig_vcf = normalizer(orig_vcf, conf)

                for caller in variants:
                    normed_caller_vcf = normalizer(variants[caller], conf)

                    for comparator_name, comparator in self.comparators.iteritems():
                        all_results = comparator(normed_orig_vcf, normed_caller_vcf, None, conf)
                        single_results = split_results(all_results, bed)
                        logging.info("Running comparator " + comparator_name + " (normalizer " + normalizer_name + ")")
                        for region, result in zip(util.read_regions(bed), single_results):
                            match_vars = util.find_matching_var(orig_vcf, region)
                            if len(match_vars) == 0:
                                raise ValueError('Unable to find original variant from region ' + str(region))

                            result = compare_single_var(result, region, normed_orig_vcf, normed_caller_vcf, comparator, "/".join([str(i) for i in match_vars[0].samples[0]['GT']]), conf)

                            match_var = "/".join([" ".join(str(mvar).split()[0:5]) for mvar in match_vars])

                            if caller not in var_results[match_var]:
                                var_results[match_var][caller] = defaultdict(dict)

                            var_results[match_var][caller][normalizer_name][comparator_name] = result

            #Iterate over all results and write to standard output. We do this here instead of within the loops above
            #because it keeps results organized by variant, which makes them easier to look at
            self.reporter.write_output(var_results, var_quals, bam_stats)

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