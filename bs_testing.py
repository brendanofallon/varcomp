
from vcomp import bam_simulation as bs
from vcomp import injectvar as ij
import ConfigParser as cp
from collections import namedtuple
import pysam

Var = namedtuple('Var', ['chrom', 'start', 'ref', 'alts'])
VarSet = namedtuple('VarSet', ['policy', 'vars'])
ExSNPInfo = namedtuple('ExSNPInfo', ['policy', 'dist'])
conf = cp.SafeConfigParser()
conf.read("comp.conf")

vars = pysam.VariantFile('test_single.vcf')
# v = [Var('9', 127818330, 'GCTG', ('G',)), Var('9', 127818335, 'C', ('A',))]
# vars = [VarSet(bs.TRANS, v)]

snp_inf = ExSNPInfo(policy=bs.TRANS, dist=-1)
sets = ij.create_variant_sets(vars, snp_inf, bs.ALL_HOMS, pysam.FastaFile(conf.get('main', 'ref_genome')))
reads = bs.gen_alt_fq(conf.get('main', 'ref_genome'), sets, 200)
bam = bs.gen_alt_bam(conf.get('main', 'ref_genome'), conf, reads)