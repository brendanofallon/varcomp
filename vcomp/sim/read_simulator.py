"""
Very simple tool for simulating paired end ngs reads. No base calling or pcr errors supported.
"""

import random
import pysam
import string

revcomp_lookup={
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

class ReadSimulator(object):
    """
    ReadSimulators generate a pairs of fastq reads from a fasta file (typically, a modified reference genome)
    This version generates 'perfect' (base-calling error free) reads, also assumes that template midpoints are
    Normally (Gaussian-ly) distributed. Template sizes also vary in a Gaussian fashion, with mean and stdev given as
     keywords args to the constructor. This creates coverage histograms that look like hybrid capture data.
    """

    def __init__(self, ref_genome, target_chr, target_mid, mean_template_size=250, stdev_template_size=50, read_len=100):
        """
        Create a new ReadSimulator
        :param ref_genome: Path to fasta file containing 'reference' to simulate from
        :param target_chr: Contig / chr name in fasta file
        :param target_mid: Base offset in reference contig describing the
        :param mean_template_size: Mean size in bp of templates
        :param stdev_template_size: Stdev of template size
        :param read_len: Read length
        """
        self.target_chr = target_chr
        self.target_pos = target_mid #Mean midpoint of templates, typically a simulated variant is close to here
        self.target_pos_stdev = 50
        self.mean_template_size = mean_template_size
        self.stdev_template_size = stdev_template_size
        self.flanking_bases=2000
        self.read_len=read_len
        ref = pysam.FastaFile(ref_genome)
        self.seq_start = max(0, target_mid - self.flanking_bases)
        self.seq_end = target_mid + self.flanking_bases
        self.seq = ref.fetch(target_chr, self.seq_start, self.seq_end)
        self.counter = 0
        self.quals = "".join(['Z' for x in range(self.read_len)])

    def gen_read_pair(self):
        """
        Generate a single read pair as a tuple of fastq-formatted strings suitable for writing to output files
        """
        template_pos = int(random.gauss(len(self.seq)/2, self.target_pos_stdev ))
        template_size = int(random.gauss(self.mean_template_size, self.stdev_template_size))
        templ_seq = self.seq[template_pos-template_size/2:template_pos+template_size/2]
        first_read = templ_seq[0:self.read_len]
        second_read = revcomp(templ_seq[-self.read_len:])
        ref_templ_mid = self.target_pos-self.flanking_bases+template_pos
        ref_read_start = ref_templ_mid-template_size/2
        rnd = "".join([ random.choice(string.ascii_lowercase + string.ascii_uppercase) for _ in range(8)])
        first_read_name = '@' + str(self.counter) + ":" + self.target_chr + ":" + rnd + ":" + str(ref_read_start)
        second_read_name = '@' + str(self.counter) + ":" + self.target_chr + ":" + rnd + ":" + str(ref_read_start)
        first_fq = first_read_name + "\n" + first_read + "\n+\n" + self.quals[0:len(first_read)]
        second_fq = second_read_name + "\n" + second_read + "\n+\n" + self.quals[0:len(second_read)]
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




