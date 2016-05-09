
import unittest
import os
from collections import namedtuple
from vcomp import util

MockVariant = namedtuple('MockVariant', ['chrom', 'start', 'ref', 'alleles', 'samples'])

class TestUtils(unittest.TestCase):

    DATA_DIR = "test_data"
    DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), DATA_DIR)
    HALFCALL_VCF = "test_halfcalls.vcf"
    TEST_BED = "test.bed"
    EMPTY_VCF = "empty.vcf"
    NORMAL_VCF = "normal.vcf"
    COMPLEX_VCF = "morecomplex.vcf"

    def test_rm_halfcalls(self):

        input_vcf = os.path.join(TestUtils.DATA_PATH, TestUtils.HALFCALL_VCF)
        try:
            no_halfcalls = util.remove_halfcalls(input_vcf)
            calls = {}
            with open(no_halfcalls) as fh:
                for line in fh:
                    if line[0] == '#':
                        continue
                    toks = line.split()
                    calls[toks[1]] = toks[9]

                self.assertTrue(len(calls) == 3)
                for pos, gt in calls.iteritems():
                    self.assertTrue('.' not in gt.split(':')[0])

        except Exception as ex:
            print str(ex)

        finally:
            os.remove(no_halfcalls)

    def test_read_bedfile(self):
        bed_path = os.path.join(TestUtils.DATA_PATH, TestUtils.TEST_BED)
        regions = [r for r in util.read_regions(bed_path)]
        self.assertTrue(len(regions) == 3)
        self.assertTrue(regions[0].chr == '1')
        self.assertTrue(regions[0].start == 10)
        self.assertTrue(regions[0].end == 20)

        self.assertTrue(regions[1].chr == '1')
        self.assertTrue(regions[1].start == 100)
        self.assertTrue(regions[1].end == 200)

        self.assertTrue(regions[2].chr == '10')
        self.assertTrue(regions[2].start == 55)
        self.assertTrue(regions[2].end == 77)

    def test_vars_to_bed(self):
        vars = [
                MockVariant(chrom='1', start=1000, ref='A', alleles=('G', ), samples=None),
                MockVariant(chrom='1', start=5000,  ref='A', alleles=('G', ), samples=None),
                MockVariant(chrom='1', start=5050,  ref='T', alleles=('G', ), samples=None),
                MockVariant(chrom='2', start=2000, ref='T', alleles=('G',), samples=None),
                MockVariant(chrom='2', start=2010, ref='T', alleles=('G',), samples=None),
                MockVariant(chrom='2', start=3000, ref='T', alleles=('G',), samples=None),
            ]

        bed = util.vars_to_bed(vars, window=500)
        regions = [r for r in util.read_regions(bed)]
        self.assertTrue(len(regions)==4)
        self.assertTrue(regions[0].chr == '1')
        self.assertTrue(regions[0].start == 500)
        self.assertTrue(regions[0].end == 1500)

        self.assertTrue(regions[1].chr == '1')
        self.assertTrue(regions[1].start == 4500)
        self.assertTrue(regions[1].end == 5550)

        self.assertTrue(regions[2].chr == '2')
        self.assertTrue(regions[2].start == 1500)
        self.assertTrue(regions[2].end == 2510)

        self.assertTrue(regions[3].chr == '2')
        self.assertTrue(regions[3].start == 2500)
        self.assertTrue(regions[3].end == 3500)

        os.remove(bed)


    def test_empty_vcf(self):
        empty_vcf = os.path.join(TestUtils.DATA_PATH, TestUtils.EMPTY_VCF)
        non_empty_vcf = os.path.join(TestUtils.DATA_PATH, TestUtils.HALFCALL_VCF)

        self.assertTrue(util.is_empty(empty_vcf))
        self.assertFalse(util.is_empty(non_empty_vcf))

    # def test_find_regions(self):
    #
    #     vars = [
    #         MockVariant(chrom='1', start=10, end=20, ref='A', alt='G'),
    #         MockVariant(chrom='1', start=10, end=20, ref='A', alt='G'),
    #         MockVariant(chrom='1', start=10, end=20, ref='A', alt='G'),
    #     ]

    def test_gtmod(self):
        orig_vcf = os.path.join(TestUtils.DATA_PATH, TestUtils.COMPLEX_VCF)

        orig_vars = ["-".join(line.split('\t')[0:5])
                     for line in open(orig_vcf)
                     if line[0] != '#']

        mod_vcf = util.set_genotypes(orig_vcf, "0/1", None, None, compress_result=False)
        new_vars = ["-".join(line.split('\t')[0:5])
                     for line in open(mod_vcf)
                     if line[0] != '#']

        self.assertListEqual(orig_vars, new_vars)

        for line in open(mod_vcf):
            if line[0]=='#':
                continue
            else:
                toks = line.split('\t')
                alt = toks[4]
                gt = toks[9].split(":")[0]
                if alt != util.VCF_MISSING:
                    self.assertTrue(gt == "0/1")
                else:
                    self.assertTrue(gt == "0")

        os.remove(mod_vcf)


if __name__=="__main__":
    unittest.main()