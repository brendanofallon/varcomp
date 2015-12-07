"""
Very simple tool for simulating paired end ngs reads. No base calling or pcr errors supported.
"""

import random
import pysam

revcomp_lookup={
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}
class ReadSimulator(object):

    def __init__(self, ref_genome, target_chr, target_pos, mean_template_size=250, stdev_template_size=50, read_len=100):
        self.target_chr = target_chr
        self.target_pos = target_pos
        self.target_pos_stdev = 50
        self.mean_template_size = mean_template_size
        self.stdev_template_size = stdev_template_size
        self.flanking_bases=2000
        self.read_len=read_len
        ref = pysam.FastaFile(ref_genome)
        self.seq = ref.fetch(target_chr, max(0, target_pos-self.flanking_bases), target_pos+self.flanking_bases)
        self.counter = 0
        self.quals = "".join(['~' for x in range(self.read_len)])

    def gen_read_pair(self):
        """
        Generate a read pair and return a tuple of fastq-formatted strings suitable for writing to output files
        """

        template_pos = int(random.gauss(self.target_pos, self.target_pos_stdev ))
        #template_pos = int(random.gauss(len(self.seq)/2, self.target_pos_stdev ))
        template_size = int(random.gauss(self.mean_template_size, self.stdev_template_size))
        templ_seq = self.seq[template_pos-template_size/2:template_pos+template_size/2]
        first_read = templ_seq[0:self.read_len]
        second_read = revcomp(templ_seq[-self.read_len:])
        ref_templ_mid = self.target_pos-self.flanking_bases+template_pos
        ref_read_start = ref_templ_mid-template_size/2
        # print "Mid: " + str(template_pos) + " length: " + str(template_size)
        # print "Ref templ midpoint: " + str(ref_templ_mid) + " Read 1 start:" + str(ref_read_start)
        # print "Seq: " + templ_seq
        # print "1st: " + first_read
        # print "2nd: " + second_read
        first_read_name = '@' + str(self.counter) + ":" + self.target_chr + ":" + str(ref_read_start) + "/1"
        second_read_name = '@' + str(self.counter) + ":" + self.target_chr + ":" + str(ref_read_start) + "/2"
        first_fq = first_read_name + "\n" + first_read + "\n+\n" + self.quals
        second_fq = second_read_name + "\n" + second_read + "\n+\n" + self.quals
        self.counter += 1
        return (first_fq, second_fq)


def revcomp(bases):
    """
    Return reverse-complemented bases, uses the revcomp_lookup global lookup dictionary
    :return: Reverse-complemented bases as a string
    """
    result = []
    for b in bases[::-1]:
        result.append(revcomp_lookup[b])
    return "".join(result)

if __name__=="__main__":
    gen = ReadSimulator("/Users/bofallon/resources/GRCh38.p2/GRCh38.p2.fa", "9", 127830920)

    prefix = "test"
    a_output = open(prefix + "_R1.fastq", "w")
    b_output = open(prefix + "_R2.fastq", "w")
    for x in range(100):
        (a, b) = gen.gen_read_pair()
        a_output.write(a + "\n")
        b_output.write(b + "\n")
    a_output.close()
    b_output.close()



