
import subprocess
import pysam
import time

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

def call_variant_fb(bam, orig_genome_path, contig, start, end, conf=None):
    vcfoutput = "output-fb.vcf"
    cmd=[conf.get('main', 'freebayes_path'), "-f", orig_genome_path, "-b", bam, "-v", vcfoutput]
    subprocess.check_output(cmd)
    return compress_vcf(vcfoutput, conf)

def call_variant_platypus(bam, orig_genome_path, contig, start, end, conf=None):
    vcfoutput = "output-platypus.vcf"
    #err = open("/dev/null")
    #cmd=["python", conf.get('main', 'platypus_path'), "callVariants", "--refFile", orig_genome_path, "--bamFiles", bam, "--regions", contig+":" + str(start) + "-" + str(end), "-o", vcfoutput]
    cmd= "python " + conf.get('main', 'platypus_path') + " callVariants --refFile " + orig_genome_path + " --bamFiles " + bam + " --regions " + contig + ":" + str(start) + "-" + str(end) + " -o " + vcfoutput
    subprocess.check_call(cmd, shell=True)
    #err.close()
    return compress_vcf(vcfoutput, conf)


def call_variant_gatk_hc(bam, orig_genome_path, contig, start, end, conf=None):
    vcfoutput = "output-hc.vcf"
    err = open("/dev/null")
    cmd="java -Xmx1g -jar " + conf.get('main', 'gatk_path') + " -T HaplotypeCaller -R " + orig_genome_path +" -I " + bam + " -L " + contig + ":" + str(start) + "-" + str(end) + " -o " + vcfoutput
    #print "Executing " + cmd
    subprocess.check_output(cmd, shell=True)
    err.close()
    return compress_vcf(vcfoutput, conf)



def call_variant_rtg(bam, orig_genome_path, contig, start, end, conf=None):
    output_dir = "rtg-output-" + str(time.time()).replace(".", "")[-7:]
    vcfoutput = output_dir + "/snps.vcf.gz"
    cmd=["java", "-jar", conf.get('main', 'rtg_jar'), "snp", "-t", conf.get('main', 'rtg_ref_sdf'), "-o", output_dir, bam]
    subprocess.check_output(cmd)
    #vars = pysam.VariantFile(vcfoutput, "rb")
    return vcfoutput


def get_callers():
    return {
        "freebayes": call_variant_fb,
        "platypus": call_variant_platypus,
        "rtg": call_variant_rtg,
        "gatk-hc": call_variant_gatk_hc
    }