
import unittest
import os
from collections import namedtuple
from vcomp import util

MockVariant = namedtuple('MockVariant', ['chrom', 'start', 'end', 'ref', 'alleles', 'samples'])

class TestUtils(unittest.TestCase):

    DATA_DIR = "test_data"
    DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), DATA_DIR)
    HALFCALL_VCF = "test_halfcalls.vcf"
    TEST_BED = "test.bed"
    EMPTY_VCF = "empty.vcf"
    NORMAL_VCF = "normal.vcf"

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


    def test_empty_vcf(self):
        empty_vcf = os.path.join(TestUtils.DATA_PATH, TestUtils.EMPTY_VCF)
        non_empty_vcf = os.path.join(TestUtils.DATA_PATH, TestUtils.HALFCALL_VCF)

        self.assertTrue(util.is_empty(empty_vcf))
        self.assertFalse(util.is_empty(non_empty_vcf))

    def test_find_regions(self):

        vars = [
            MockVariant(chrom='1', start=10, end=20, ref='A', alt='G'),
            MockVariant(chrom='1', start=10, end=20, ref='A', alt='G'),
            MockVariant(chrom='1', start=10, end=20, ref='A', alt='G'),
        ]

if __name__=="__main__":
    unittest.main()