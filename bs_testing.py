
from vcomp import bam_simulation as bs
import ConfigParser as cp
import pysam

conf = cp.SafeConfigParser()
conf.read("comp.conf")

vars = pysam.VariantFile('test_mnp.vcf')
reads = bs.gen_alt_fq(conf.get('main', 'ref_genome'), vars, 200, policy=bs.EACH_ALT)
bam = bs.gen_alt_bam(conf.get('main', 'ref_genome'), conf, reads)