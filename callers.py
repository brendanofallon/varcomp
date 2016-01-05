
import subprocess
import time
import string
import random

def compress_vcf(input_vcf, conf):
    """
    If the input vcf's filename does not end with .gz, compress and index it with bgzip / tabix
    :param input_vcf:
    :param conf:
    :return: Name of compressed vcf file, typically input_vcf + '.gz'
    """
    if not input_vcf.endswith(".gz"):
        cmd = conf.get('main', 'bgzip_path') + " -f " + input_vcf
        subprocess.check_call(cmd, shell=True)
        cmd = conf.get('main', 'tabix_path') + " -f " + input_vcf + ".gz"
        subprocess.check_call(cmd, shell=True)
        input_vcf = input_vcf + '.gz'
    return input_vcf

def vars_to_bed(variants, window=250):
    """
    Generate a bed file containing regions that span each variant, each region is centered on the variant start position and extends
    'window' bp in each direction
    :param variants:
    :param window:
    :return:
    """
    bedfilename = "var_regions" + "".join([random.choice(string.ascii_lowercase + string.ascii_uppercase) for _ in range(10)]) + ".bed"
    bfh = open(bedfilename, "w")
    for var in variants:
        bfh.write("\t".join([var.chrom, str(var.start-window), str(var.start+window)]) + "\n")
    bfh.close()
    return bedfilename

def call_variant_fb(bam, orig_genome_path, bed, conf=None):
    vcfoutput = "output-fb.vcf"
    cmd=[conf.get('main', 'freebayes_path'), "-f", orig_genome_path, "-t", bed, "-b", bam, "-v", vcfoutput]
    #print "Executing " + " ".join(cmd)
    subprocess.check_output(cmd)
    return compress_vcf(vcfoutput, conf)

def call_variant_platypus(bam, orig_genome_path, bed, conf=None):
    vcfoutput = "output-platypus.vcf"
    #err = open("/dev/null")
    #cmd=["python", conf.get('main', 'platypus_path'), "callVariants", "--refFile", orig_genome_path, "--bamFiles", bam, "--regions", contig+":" + str(start) + "-" + str(end), "-o", vcfoutput]
    cmd= "python " + conf.get('main', 'platypus_path') + " callVariants --refFile " + orig_genome_path + " --bamFiles " + bam + " --regions " + bed + " -o " + vcfoutput
    subprocess.check_call(cmd, shell=True)
    #err.close()
    return compress_vcf(vcfoutput, conf)

def call_wecall(bam, orig_genome_path, bed, conf=None):
    vcfoutput = "output-wc.vcf"
    cmd=conf.get('main', 'wecall_path') + " --refFile " + orig_genome_path + " --inputs " + bam + " --regions " + bed + " --output " + vcfoutput
    subprocess.check_call(cmd, shell=True)
    return compress_vcf(vcfoutput, conf)

def call_variant_gatk_hc(bam, orig_genome_path, bed, conf=None):
    vcfoutput = "output-hc.vcf"
    err = open("/dev/null")
    cmd="java -Xmx1g -jar " + conf.get('main', 'gatk_path') + " -T HaplotypeCaller -R " + orig_genome_path +" -I " + bam + " -L " + bed + " -o " + vcfoutput
    #print "Executing " + cmd
    subprocess.check_output(cmd, shell=True, stderr=err)
    err.close()
    return compress_vcf(vcfoutput, conf)



def call_variant_rtg(bam, orig_genome_path, bed, conf=None):
    output_dir = "rtg-output-" + str(time.time()).replace(".", "")[-7:]
    vcfoutput = output_dir + "/snps.vcf.gz"
    cmd=["java", "-jar", conf.get('main', 'rtg_jar'), "snp", "-t", conf.get('main', 'rtg_ref_sdf'), "--bed-regions", bed, "-o", output_dir, bam]
    subprocess.check_output(cmd)
    #vars = pysam.VariantFile(vcfoutput, "rb")
    return vcfoutput


def get_callers():
    return {
        "freebayes": call_variant_fb,
        "platypus": call_variant_platypus,
        "rtg": call_variant_rtg,
        "gatk-hc": call_variant_gatk_hc
        #"wecall": call_wecall
    }
