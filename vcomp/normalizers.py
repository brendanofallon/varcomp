
import subprocess
import random
import string
import util

def normalize_nothing(orig_vcf, conf):
    """
    Just copy the original vcf to a new, identical vcf file.
    """
    newvcf_name = orig_vcf.replace(".vcf", ".nonorm.vcf")
    cmd = "cp " + orig_vcf + " " + newvcf_name
    subprocess.check_call(cmd, shell=True)

    newvcf_name = util.bgz_tabix(newvcf_name, conf)
    return newvcf_name

def normalize_vap_leftalign(orig_vcf, conf):
    err = open("/dev/null")

    orig_vcf = util.sort_vcf(orig_vcf, conf)
    tmp_vcf = orig_vcf.replace(".vcf", ".vap.tmp.vcf").replace(".gz", "")
    final_vcf = orig_vcf.replace(".vcf", ".vap.leftaligned.vcf")
    norm_orig_cmd = conf.get('main', 'vcfallelicprimitives_path') + " " + orig_vcf
    tmp_output=subprocess.check_output(norm_orig_cmd, shell=True)
    with open(tmp_vcf, "w") as fh:
        fh.write(tmp_output)


    cmd = "java -Djava.io.tmpdir=. -Xmx1g -jar " + conf.get('main', 'gatk_path') + " -T LeftAlignAndTrimVariants -R " + conf.get('main', 'ref_genome') + " -V " + tmp_vcf + " -o " + final_vcf
    subprocess.check_output(cmd, shell=True)
    err.close()

    return util.bgz_tabix(final_vcf, conf)

def normalize_vt(orig_vcf, conf):
    """
    Use vt to normalize
    :param conf: configuration object with path to reference genome and vt binary
    :return: String describing variant matching result
    """
    err = open("/dev/null")
    norm_orig_vcf = orig_vcf.replace(".vcf", ".norm.vt.vcf")
    norm_orig_cmd = conf.get('main', 'vt_path') + " normalize " + " -r " + conf.get('main', 'ref_genome') + " " + orig_vcf + " -o " + norm_orig_vcf
    subprocess.check_output(norm_orig_cmd.split(), stderr=err)
    norm_orig_vcf = util.bgz_tabix(norm_orig_vcf, conf)
    err.close()
    return norm_orig_vcf

def normalize_bcftools(orig_vcf, conf):
    """
    Use bcftools to normalize.
    :param orig_vcf:
    :param conf:
    :return:
    """
    norm_orig_vcf = orig_vcf.replace(".vcf.gz", ".norm.bcftools" + rndstr(6) + ".vcf")
    norm_orig_cmd = conf.get('main', 'bcftools_path') + " norm " + " -c w -f " + conf.get('main', 'ref_genome') + " " + orig_vcf + " -o " + norm_orig_vcf
    subprocess.check_call(norm_orig_cmd.split())
    norm_orig_vcf = util.bgz_tabix(norm_orig_vcf, conf)
    return norm_orig_vcf


def get_normalizers():
    return {
        'vapleft': normalize_vap_leftalign,
        'nonorm': normalize_nothing,
        'vt': normalize_vt,
        # 'bcftools': normalize_bcftools
    }


def rndstr(length):
    return "".join([random.choice(string.ascii_lowercase+string.ascii_uppercase) for _ in range(length)])