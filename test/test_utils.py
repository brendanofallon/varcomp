
import unittest
import os
from vcomp import util

class TestUtils(unittest.TestCase):

    DATA_DIR = "test_data"
    DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), DATA_DIR)
    HALFCALL_VCF = "test_halfcalls.vcf"

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


if __name__=="__main__":
    unittest.main()