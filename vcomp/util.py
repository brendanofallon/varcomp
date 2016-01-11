
import subprocess
import gzip
import time

def bgz_tabix(path, conf):
    """
    If the path does not end in .gz bgzip the file, then index with tabix and return the potentially modified filename
    :return: Filename of compressed file
    """
    try:
        if not path.endswith(".gz"):
            cmd = conf.get('main', 'bgzip_path') + " " + path
            subprocess.check_call(cmd.split())
            path = path + ".gz"

        cmd = conf.get('main', 'tabix_path') + " -f " + path
        subprocess.check_call(cmd.split())
    except Exception as ex:
        raise ex


    return path



def set_genotypes(orig_vcf, newGT, conf):
    """
    Create a new VCF file that is identical to the given VCF, except that all GT info fields are set to 'newGT'
    """
    fh = None
    if orig_vcf.endswith(".gz"):
        fh = gzip.open(orig_vcf, "r")
    else:
        fh = open(orig_vcf, "r")

    newvcf = orig_vcf.replace(".vcf", ".gtmod" + str(time.time())[-6:].replace(".", "") + ".vcf").replace(".gz", "")
    ofh =open(newvcf, "w")
    for line in fh.readlines():
        if len(line)==0 or line[0]=='#':
            ofh.write(line)
        else:
            toks = line.split('\t')
            if len(toks)<10:
                ofh.write(line)
            else:
                infoitems = [newGT]
                if ':' in toks[9]:
                    infoitems.extend(toks[9].strip().split(':')[1:] )
                newinfo = ":".join(infoitems)
                ofh.write('\t'.join(toks[0:9] + [newinfo]) + "\n")
    fh.close()
    ofh.close()
    bgz_vcf = bgz_tabix(newvcf, conf)
    return bgz_vcf