import time
import subprocess

from vcomp import util


def get_callers():
    return {
        #"freebayes": call_variant_fb,
        #"samtools": call_variant_mp_bcf,
        #"platypus": call_variant_platypus,
        #"varscan": call_variant_varscan,
        #platypus-asm": call_variant_platypus_asm,
        #"rtg": call_variant_rtg,
        "gatk-hc": call_variant_gatk_hc,
        "gatk-hc-somatic": call_variant_gatk_hc_somatic,
        #"sentieon-hc": call_variant_sentieon_hc,
        #"somaticppl": call_variant_somaticppl,
		#"wecall": call_wecall,
        #"gatk-ug": call_variant_gatk_ug
        #"wecall": call_wecall
        #"freebayes-mre": call_variant_fb_minrepeatentropy,

        #"gatk-hc": call_variant_gatk_hc_emit_all,
        #"gatk-ug": call_variant_gatk_ug_emit_all,
        #"varscan": call_variant_varscan_emit_all,
    }


def call_variant_platypus_asm(bam, genome, bed, conf=None):
    vcfoutput = "output-platypus.vcf"
    cmd= "python " + conf.get('main', 'platypus') + " callVariants --assemble=1 --assembleBadReads=1 --refFile " + genome + " --bamFiles " + bam + " --regions " + bed + " -o " + vcfoutput
    subprocess.check_call(cmd)
    return util.compress_vcf(vcfoutput, conf)


def call_variant_fb(bam, genome, bed, conf=None):
    vcfoutput = "output-fb.vcf"
    cmd=[conf.get('main', 'freebayes'), "-f", genome, "-t", bed, "-b", bam, "-v", vcfoutput]
    subprocess.check_call(cmd, stdout=open('/dev/null'))
    sorted_vcf = util.sort_vcf(vcfoutput, conf)
    return util.compress_vcf(sorted_vcf, conf)


def call_variant_fb_minrepeatentropy(bam, genome, bed, conf=None):
    vcfoutput = "output-fb.vcf"
    cmd=[conf.get('main', 'freebayes'), "-f", genome, "--no-partial-observations", "--min-repeat-entropy", "1", "-t", bed, "-b", bam, "-v", vcfoutput]
    subprocess.check_call(cmd, stdout=open('/dev/null'))
    sorted_vcf = util.sort_vcf(vcfoutput, conf)
    return util.compress_vcf(sorted_vcf, conf)


def call_variant_platypus(bam, genome, bed, conf=None):
    vcfoutput = "output-platypus.vcf"
    cmd= "python " + conf.get('main', 'platypus') + " callVariants --refFile " + genome + " --bamFiles " + bam + " --regions " + bed + " -o " + vcfoutput
    subprocess.check_call(cmd, stdout=open('/dev/null'))
    return util.compress_vcf(vcfoutput, conf)


def call_wecall(bam, genome, bed, conf=None):
    vcfoutput = "output-wc.vcf"
    cmd=conf.get('main', 'wecall') + " --refFile " + genome + " --inputs " + bam + " --regions " + bed + " --output " + vcfoutput
    subprocess.check_call(cmd, stdout=open('/dev/null'))
    return util.compress_vcf(vcfoutput, conf)


def call_variant_gatk_hc(bam, genome, bed, conf=None):
    vcfoutput = "output-hc.vcf"

    cmd='''java -Xmx1g -Djava.io.tmpdir=. -jar {gatk} -T HaplotypeCaller -R {genome} -I {bam} -U ALLOW_SEQ_DICT_INCOMPATIBILITY \
           -L {bed} -o {vcfoutput}'''
    cmd = cmd.format(gatk=conf.get('main', 'gatk'), genome=genome, bam=bam, bed=bed, vcfoutput=vcfoutput)

    subprocess.check_call(cmd, stdout=open('/dev/null'), stderr=subprocess.STDOUT, shell=True)
    return util.compress_vcf(vcfoutput, conf)

def call_variant_gatk_hc_somatic(bam, genome, bed, conf=None):
    vcfoutput = "output-hc-somatic.vcf"

    cmd='''java -Xmx1g -Djava.io.tmpdir=. -jar {gatk} -T HaplotypeCaller -R {genome} -I {bam} -U ALLOW_SEQ_DICT_INCOMPATIBILITY \
           -L {bed} -o {vcfoutput} -ploidy 20'''
    cmd = cmd.format(gatk=conf.get('main', 'gatk'), genome=genome, bam=bam, bed=bed, vcfoutput=vcfoutput)

    subprocess.check_call(cmd, stdout=open('/dev/null'), stderr=subprocess.STDOUT, shell=True)

    return util.set_genotypes(vcfoutput, "0/1", None, conf)


def call_variant_sentieon_hc(bam, genome, bed, conf=None):
    vcfoutput = "output-sentieonhc.vcf"
    cmd='{sentieon} driver -r {genome} -i {bam} --algo Haplotyper {vcfoutput}'
    cmd=cmd.format(sentieon=conf.get('main', 'sentieon'), genome=genome, bam=bam, vcfoutput=vcfoutput)
    subprocess.check_call(cmd, stdout=open('/dev/null'), stderr=subprocess.STDOUT, shell=True)
    return util.compress_vcf(vcfoutput, conf)

def call_variant_somaticppl(bam, genome, bed, conf=None):
	vcfoutput = "output-somaticppl.vcf"
	cmd = '''{somaticcaller} {genome} {bed} {bam} {output}'''.format(somaticcaller=conf.get('main', 'somaticppl_path'), bam=bam, genome=genome, bed=bed, output=vcfoutput)
	subprocess.check_call(cmd, shell=True)
	return util.compress_vcf(vcfoutput, conf)

def call_variant_gatk_ug(bam, genome, bed, conf=None):
    vcfoutput = "output-ug.vcf"
    no_et = ""
    try:
        no_et = " -et NO_ET -K " + conf.get('main', 'gatk_no_et')
    except:
        pass
    cmd="java -Xmx1g -Djava.io.tmpdir=. -jar " + conf.get('main', 'gatk') + " -T UnifiedGenotyper -glm BOTH " + no_et + " -R " + genome +" -I " + bam + " -L " + bed + " -o " + vcfoutput
    subprocess.check_call(cmd, stdout=open('/dev/null'), stderr=subprocess.STDOUT)
    return util.compress_vcf(vcfoutput, conf)



def call_variant_gatk_hc_emit_all(bam, genome, bed, conf=None):
    vcfoutput = "output-hc.vcf"
    no_et = ""
    try:
        no_et = " -et NO_ET -K " + conf.get('main', 'gatk_no_et')
    except:
        pass
    cmd="java -Xmx1g -Djava.io.tmpdir=. -jar " + conf.get('main', 'gatk') + " -T HaplotypeCaller -stand_emit_conf 1.0 " + no_et + " -R " + genome +" -I " + bam + " -U ALLOW_SEQ_DICT_INCOMPATIBILITY -L " + bed + " -o " + vcfoutput
    subprocess.check_call(cmd, stdout=open('/dev/null'), stderr=subprocess.STDOUT)
    return util.compress_vcf(vcfoutput, conf)


def call_variant_gatk_ug_emit_all(bam, genome, bed, conf=None):
    vcfoutput = "output-ug.vcf"
    no_et = ""
    try:
        no_et = " -et NO_ET -K " + conf.get('main', 'gatk_no_et')
    except:
        pass
    cmd="java -Xmx1g -Djava.io.tmpdir=. -jar " + conf.get('main', 'gatk') + " -T UnifiedGenotyper -glm BOTH -stand_emit_conf 1.0 " + no_et + " -R " + genome +" -I " + bam + " -L " + bed + " -o " + vcfoutput
    subprocess.check_call(cmd, stdout=open('/dev/null'), stderr=subprocess.STDOUT)
    return util.compress_vcf(vcfoutput, conf)


def call_variant_rtg(bam, genome, bed, conf):
    output_dir = "rtg-output-" + util.randstr()
    vcfoutput = output_dir + "/snps.vcf.gz"
    cmd=["java", "-Djava.io.tmpdir=.", "-jar", conf.get('main', 'rtg_jar'), "snp", "-t", conf.get('main', 'rtg_ref_sdf'), "--bed-regions", bed, "-o", output_dir, bam]
    subprocess.check_call(cmd, stdout=open('/dev/null'), stderr=subprocess.STDOUT)
    return vcfoutput


def call_variant_varscan(bam, genome, bed, conf):
    pre_output = "varscan." + util.randstr() + ".mpileup"
    vcfoutput = "output-vs." + util.randstr() + ".vcf"
    bedarg = ""
    if bed is not None:
        bedarg = " -l " + bed
    cmd = conf.get('main','samtools') + ' mpileup ' + ' -f ' + genome + " -o " + pre_output + " " + bedarg + " " + bam
    subprocess.check_call(cmd)
    cmd2 = "java -Xmx2g -jar " + conf.get('main', 'varscan') + ' mpileup2cns ' + pre_output + ' --variants --output-vcf 1 --output-file ' + vcfoutput
    subprocess.check_call(cmd2, stdout=open(vcfoutput, 'w'))
    return util.bgz_tabix(vcfoutput, conf)


def call_variant_varscan_emit_all(bam, genome, bed, conf):
    pre_output = "varscan." + util.randstr() + ".mpileup"
    vcfoutput = "output-vs." + util.randstr() + ".vcf"
    bedarg = ""
    if bed is not None:
        bedarg = " -l " + bed
    cmd = conf.get('main','samtools') + ' mpileup ' + ' -f ' + genome + " -o " + pre_output + " " + bedarg + " " + bam
    subprocess.check_call(cmd)
    cmd2 = "java -Xmx2g -jar " + conf.get('main', 'varscan') + ' mpileup2cns ' + pre_output + '  --p-value 0.5 --variants --output-vcf 1 --output-file ' + vcfoutput
    subprocess.check_call(cmd2, stdout=open(vcfoutput, 'w'))
    return util.bgz_tabix(vcfoutput, conf)


def call_variant_mp_bcf(bam, genome, bed, conf):
    pre_output = "mpileup." + util.randstr() + ".vcf"
    vcfoutput = "output-mp." + util.randstr() + ".vcf"
    bedarg = ""
    if bed is not None:
        bedarg = " -l " + bed
    cmd = conf.get('main','samtools') + ' mpileup ' + ' -f ' + genome + " -uv " + " -o " + pre_output + " " + bedarg + " " + bam
    subprocess.check_call(cmd)
    cmd2 = conf.get('main', 'bcftools') + ' call ' + ' -mv ' + ' -o ' + vcfoutput + " " + pre_output
    subprocess.check_call(cmd2)
    return util.bgz_tabix(vcfoutput, conf)
