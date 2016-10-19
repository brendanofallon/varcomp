import subprocess

from vcomp import util


def get_normalizers():
    return {
        # 'vapleft': normalize_vap_leftalign,
        'nonorm': normalize_nothing,
        # 'vt': normalize_vt,
        # 'bcftools': normalize_bcftools
    }


def normalize_nothing(orig_vcf, conf):
    """
    Just copy the original vcf to a new, identical vcf file.
    """
    new_vcf = orig_vcf.replace(".vcf", ".nonorm.vcf")
    cmd = 'cp {} {}'.format(orig_vcf, new_vcf)
    subprocess.check_call(cmd, shell=True)
    return util.bgz_tabix(new_vcf, conf)


def normalize_vap_leftalign(orig_vcf, conf):
    orig_vcf = util.sort_vcf(orig_vcf, conf)
    tmp_vcf = orig_vcf.replace(".vcf", ".vap.tmp.vcf").replace(".gz", "")
    final_vcf = orig_vcf.replace(".vcf", ".vap.leftaligned.vcf")
    norm_orig_cmd = conf.get('main', 'vcfallelicprimitives') + " " + orig_vcf
    subprocess.check_call(norm_orig_cmd, shell=True, stdout=file(tmp_vcf, 'w'))

    no_et = ""
    try:
        no_et = " -et NO_ET -K " + conf.get('main', 'gatk_no_et')
    except:
        pass

    cmd = "java -Djava.io.tmpdir=. -Xmx1g -jar " + conf.get('main', 'gatk') + " -T LeftAlignAndTrimVariants " + no_et + " -R " + conf.get('main', 'ref_genome') + " -U ALLOW_SEQ_DICT_INCOMPATIBILITY -V " + tmp_vcf + " -o " + final_vcf
    subprocess.check_call(cmd, shell=True)

    return util.bgz_tabix(final_vcf, conf)


def normalize_vt(orig_vcf, conf):
    """
    Use vt to normalize
    :param conf: configuration object with path to reference genome and vt binary
    :return: String describing variant matching result
    """
    norm_orig_vcf = orig_vcf.replace(".vcf", ".norm.vt.vcf")
    norm_orig_cmd = conf.get('main', 'vt') + " normalize " + " -r " + conf.get('main', 'ref_genome') + " " + orig_vcf + " -o " + norm_orig_vcf
    subprocess.check_call(norm_orig_cmd.split(), stderr=open("/dev/null"))
    return util.bgz_tabix(norm_orig_vcf, conf)


def normalize_bcftools(orig_vcf, conf):
    """
    Use bcftools to normalize.
    :param orig_vcf:
    :param conf:
    :return:
    """
    norm_orig_vcf = orig_vcf.replace(".vcf.gz", ".norm.bcftools" + util.randstr() + ".vcf")
    norm_orig_cmd = conf.get('main', 'bcftools') + " norm " + " -c w -f " + conf.get('main', 'ref_genome') + " " + orig_vcf + " -o " + norm_orig_vcf
    subprocess.check_call(norm_orig_cmd.split())
    return util.bgz_tabix(norm_orig_vcf, conf)
