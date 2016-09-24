import os
import sys
import imp
import json
import logging
import random
import argparse

import ConfigParser as cp
from collections import namedtuple

import pysam


import vcomp.util as util
import vcomp.batch_processor as bp
from   vcomp.sim import bam_simulation

from vcomp.plugins import core_callers, normalizers as core_norms, comparators as core_comps


all_result_types = (bp.MATCH_RESULT, bp.NO_MATCH_RESULT, bp.NO_VARS_FOUND_RESULT, bp.MATCH_WITH_EXTRA_RESULT,
                    bp.ZYGOSITY_MISSING_ALLELE, bp.ZYGOSITY_EXTRA_ALLELE, bp.ERROR_RESULT)

ExSNPInfo = namedtuple('ExSNPInfo', ['policy', 'dist'])


class JsonReporter(object):
    """
    An output reporter that writes a little json object for each input variant
    """

    def __init__(self, outputfile=sys.stdout):
        self.output = outputfile

    def write_output(self, results, quals, bamstats):
        """
        Write output for a batch of input variants (with individual entries for each caller/normalizer/comparator
          combination) to the given output handle.
        :param results: Four-level deep dict containing [input variant string][caller][normalizer][comparator]
        :param bamstats: Dictionary containing statistics for bam file
        :param output: File-like object to which formatted output will be written
        """
        for var, vresults in results.iteritems():
            json.dump({
                "variant":var,
                "caller_quals": quals[var],
                "bamstats": bamstats[var],
                "results": vresults
            }, self.output)
            self.output.write("\n")


def gen_reads(vcf, dest_vcf, dest_fq_prefix, ex_snp, gt_policy, read_depth, conf):
    """
    Generate fastqs for the given set of input variants. This code is fired when the user supplies the --generate-fqs
    arg, and closely mimics the fastq generation code in VariantProcessor
    :param vars: List of variants
    :param dest_vcf: Destination filename for final VCF (may be gzipped)
    :param dest_fq_prefix: Destination prefix for fastq files
    :param ex_snp: Info for extra SNP addition
    :param gt_policy: Policy describing genotype (hets, homs, from file, etc.)
    :param read_depth:
    :param conf:
    """

    #First, make sure there aren't variants that are too close to process independently...
    batches = util.batch_variants(vcf, max_batch_size=1e9)
    if len(list(batches))>1:
        raise ValueError('The VCF file ' + vcf + ' contains variants that are too close to include in a single set of fastqs, please ensure no two variants are within 2kb of each other')
    vars = list(pysam.VariantFile(vcf))
    variant_sets = bp.create_variant_sets(vars, ex_snp, gt_policy, pysam.FastaFile( conf.get('main', 'ref_genome')))
    allvars = []
    for vset in variant_sets:
        allvars.extend(vset['vars'])
    variant_batch = sorted(allvars, cmp=util.variant_comp)
    final_vcf = util.write_vcf(variant_batch, dest_vcf, conf)
    logging.info("Writing full VCF to " + final_vcf)
    reads = bam_simulation.gen_alt_fq(conf.get('main', 'ref_genome'), variant_sets, read_depth, dest_prefix=dest_fq_prefix)
    logging.info("Writing fastqs to " + reads[0] + ", " + reads[1])



def load_components(conf, section, callable_name):
    """
    Create a dict of string -> callables by examining the configuration object, loading any modules
    defined in the given section, and then calling the 'callable_name' function in the module. That function should
    return a dict, which we add to the dict to return.
    :param conf:  Configuration object
    :param section:  Section to examine in configuration (e.g. 'callers')
    :param callable_name: Name of callable to execute to get new components (e.g. 'get_callers')
    :return: Dict containing loaded components, mapping string -> callable
    """
    components = {}

    #Allow missing section
    try:
        items = conf.items(section)
    except:
        return components

    for item in items:
        try:
            if not os.path.isabs(item[1]):
                pdir = os.path.split(__file__)[0]
                mod_path = os.path.split( pdir )[0] + "/" + item[1]
            else:
                mod_path = item[1]
            logging.info("Loading plugins from module " + mod_path)
            mod = imp.load_source(item[0], mod_path)
            if not callable_name in dir(mod):
                raise ImportError('Module ' + item[1] + ' does not define a function called ' + callable_name)
            result = mod.__dict__[callable_name]()
            components.update(result)

        except ImportError as e:
            logging.fatal("Could not import caller module " + item[1] + ": " + str(e))
            raise e

    return components


def process_vcf(vcf, gt_default, conf, output, callers, fqs=None, snp_info=None, single_batch=False, keep_tmpdir=False, read_depth=250):
    """
    Perform analyses for each variant in the VCF file.
    :param input_vcf: Path to vcf file containing variants to process
    :param single_batch: Assume all variants in VCF are part of one batch and process them all simultaneously
    :param keep_tmpdir: Preserve tmpdirs created (otherwise delete them, unless they are flagged)
    :param conf: Configuration object
    """

    variant_callers = core_callers.get_callers()
    variant_callers.update(load_components(conf, 'callers', 'get_callers'))

    normalizers = core_norms.get_normalizers()
    normalizers.update(load_components(conf, 'normalizers', 'get_normalizers'))

    comparators = core_comps.get_comparators()
    comparators.update(load_components(conf, 'comparators', 'get_comparators'))

    if callers is not None and len(callers)>0:
        callers_to_use = {}
        for caller in callers:
            if caller not in variant_callers:
                raise KeyError('No variant caller ' + caller + ' found in callers')
            callers_to_use[caller] = variant_callers[caller]
        variant_callers = callers_to_use

    if fqs is not None:
        fqs = [os.path.abspath(fq) for fq in fqs]

    processor = bp.VariantProcessor(variant_callers, normalizers, comparators, JsonReporter(output), conf)
    logging.info("Processing variants in file " + vcf)
    if single_batch:
        logging.info("Processing all variants as one batch")
        tmp_dir = '{}-tmp'.format(util.strip_extensions(os.path.basename(vcf), ['vcf','gz']))
        processor.process_batch(vcf, tmp_dir, gt_default, ex_snp=snp_info, keep_tmpdir=keep_tmpdir, read_depth=read_depth, reads=fqs)
    else:
        batches = util.batch_variants(vcf, max_batch_size=1000, min_safe_dist=2000)
        for batchnum, batch_vcf in enumerate(batches, 1):
            logging.info('Processing batch #{} of {}'.format(batchnum, len(batches)))
            tmp_dir = '{}-batch{:03d}-tmp'.format(util.strip_extensions(os.path.basename(vcf), ['vcf','gz']), batchnum)
            processor.process_batch(batch_vcf, tmp_dir, gt_default, ex_snp=snp_info, keep_tmpdir=keep_tmpdir, read_depth=read_depth, reads=fqs)
            os.remove(batch_vcf)

def process_test_vcf(orig_vcf, test_vcf, output, conf):
    normalizers = core_norms.get_normalizers()
    normalizers.update(load_components(conf, 'normalizers', 'get_normalizers'))

    comparators = core_comps.get_comparators()
    comparators.update(load_components(conf, 'comparators', 'get_comparators'))
    processor = bp.VariantProcessor({}, normalizers, comparators, JsonReporter(output), conf)
    processor.compare_test_vcf(orig_vcf, test_vcf)

def main(args):
    """
    Respond to command line args, check for basic config errors, and perform analyses
    :param args:
    """
    conf = cp.SafeConfigParser()
    conf.read(args.conf)

    output = sys.stdout
    if isinstance(args.output, str):
        output = open(args.output, "w")

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

    if args.generate_fqs:
        if len(args.vcf)>1:
            raise ValueError('Only one VCF supported for now')
        vcf = args.vcf[0]
        fastq_prefix = util.strip_extensions(vcf, ['vcf','gz'])
        logging.info("Generating reads for vcf file " + vcf)
        suffix = ".d" + str(args.readdepth)
        if gt_default == bam_simulation.ALL_HETS:
            suffix = suffix + ".het"
        elif gt_default == bam_simulation.ALL_HOMS:
            suffix = suffix + ".hom"
        fastq_prefix += suffix
        gen_reads(vcf, fastq_prefix + ".truth.vcf", fastq_prefix, snp_inf, gt_default, args.readdepth, conf)
        return

    if args.test_vcf:
        if len(args.vcf) > 1:
            raise ValueError('Only one VCF supported for now')
        process_test_vcf(args.vcf[0], args.test_vcf, output, conf)
        return

    for vcf in args.vcf:
        logging.info("Processing vcf file " + vcf)
        process_vcf(vcf, gt_default, conf, output, args.callers, fqs=args.fqs, snp_info=snp_inf,
                         single_batch=args.batch, keep_tmpdir=args.keep, read_depth=args.readdepth)


if __name__=="__main__":
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    parser = argparse.ArgumentParser("Inject, simulate, call, compare")
    parser.add_argument("-c", "--conf", help="Path to configuration file", default="./comp.conf")
    parser.add_argument("-v", "--vcf", help="Input vcf file(s)", nargs="+", required=True)
    parser.add_argument("-k", "--keep", help="Dont delete temporary directories", action='store_true')
    parser.add_argument("-b", "--batch", help="Treat each input VCF file as a single batch (default False)", action='store_true')
    parser.add_argument("-s", "--seed", help="Random seed", default=None)
    parser.add_argument("-o", "--output", help="Output destination", default=sys.stdout)
    parser.add_argument("-r", "--readdepth", help="Number of reads to generate per variant", default=200, type=int)
    parser.add_argument("-a", "--addsnp", help="Add a SNP upstream of each variant", action='store_true')
    parser.add_argument("-t", "--trans", help="If SNP is added, add it in trans (default cis)", action='store_true')
    parser.add_argument("--snphom", help="Added SNPs are homozygous (default het)", action='store_true')
    parser.add_argument("--het", help="Force all simulated variants to be hets", action='store_true')
    parser.add_argument("--hom", help="Force all simulated variants to be homozygotes", action='store_true')
    parser.add_argument("--callers", help="Comma separated list of variant callers to use (default: use all)", action='append')
    parser.add_argument("--fqs", help="Dont generate fastqs, use these instead (two entries expected)", action='append')
    parser.add_argument("--generate-fqs", help="Generate fastqs only, do not perform any variant calling or comparison", action='store_true')
    parser.add_argument("--test-vcf", help="Compare an already-called VCF file the input VCF")
    args = parser.parse_args()

    main(args)
