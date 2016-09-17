# varcomp
Tools for calling and comparing variants from simulated read sets. The primary function is to take a list of variants (in VCF format), simulate fastq reads from them, align those reads to a reference, and then run multiple variant callers, normalization tools, and comparison tools to see which callers produced results that matched the input. 

##Read simulation

Varcomp has the ability to take an input VCF file, insert the variants into a modified reference genome, and simulate fastq-formatted reads from the new genome. This is useful for generating static read sets that can  be shared, run through multiple different pipelines, etc. Basic usage looks like:

    python vcomp/injectvar.py -v my_variants.vcf --generate-fqs --het 

This creates several output files:

    my_variants_final.vcf.gz
    my_variants_R1.fastq
    my_variants_R2.fastq
    
The first file contains the set of 'final' variants, which may not be identical to the input variants (for instance, if a zygosity argument like --het or --hom was supplied, or additional snps added via the --addsnp command). The other two files contain the first and second paired-end reads.

##Variant caller benchmarking

Varcomp can also take an input variant list, simulate reads (or use existing simulated fastqs), and run several variant callers, and compare the results to the input set of variants to see which caller had the higher accuracy. Basic usage looks like:

    python vcomp/injectvar.py -v my_variants.vcf > my_output.txt
    
If you've already created an input set of fastqs, you can use them by: 

    python vcomp/injectvar.py -v my_variants.vcf --fqs my_variants_R1.fastq --fqs my_variants_R2.fastq > my_output.txt
    
If the input VCF doesn't not contain any sample or genotype information, then you must specify one of --het or --hom to tell varcomp which genotype you'd like the variants to be. For instance:

    python vcomp/injectvar.py -v my_variants.vcf --hom > my_hom_output.txt

    
###Modifying input variants

Given an input VCF, varcomp can perform several modifications of the variants prior to simulating reads. In particular, zygosities can be forced to het or hom, and extra SNPs can be added upstream of each input variant, either in cis or trans. 

    --het        Force all input variants to be heterozygous
    --hom        Force all input variants to be homozygous
    --addsnp     Add upstream a single SNP upstream of each input variant (default, variants added in cis)
    --addsnp --trans Add upstream SNP in trans
   
One of --het or --hom is *required* if the input VCF does not contain a GT format entry for each variant.


##Configuration
 
 varcomp makes extensive use of multiple external applications (callers, normalizers, aligners, comparison tools, etc). The paths to these applications must be defined in a standard python configuration file. By default, `injectvar.py` looks for a file called `comp.conf` in the active directory, but you can specify a path to it by using the `-c /path/to/configuration/file` argument. The format of the file is standard python configuration file, essentially just a list of key=value pairs where the values are the paths to various executables, with a few section headers (see example below). 
 
 A minimal configuration file must include a path to the reference genome, the path to samtools, tabix, bgzip, and bwa executables. You'll likely want at least one caller and one comparator tool to be able to do anything. Here's a quick example: 
 
    [main]
    ref_genome=/Users/example/resources/GRCh37.p5/GRCh37.p5.fa
    samtools_path=/Users/bexample/tools/samtools
    bgzip_path=/Users/example/tools/bgzip
    tabix_path=/Users/example/tools/tabix
    bwa_path=/Users/example/tools/bwa
 
    [callers]
    example_caller=/path/to/caller_module.py

    [comparators]
    example_comparator=/path/to/comparator_module.py
   
###Built-in callers and comparators

Many commonly used variant callers are comparators are 'built in' and do not require any custom code. However, they still require an entry in the configuration file to specify the path to the executable. As of 9/16/2016, varcomp contains built in support for the following callers:

  - GATK UnifiedGenotyper v3.6 (requires config entry `gatk_path=/path/to/GenomeAnalysisTK.jar`)
  - GATK HaplotypeCaller  v3.6 (requires config entry `gatk_path=/path/to/GenomeAnalysisTK.jar`)
  - Freebayes v1.0.1 (requires config entry `freebayes_path=/path/to/freebayes`)
  - bcftools 1.3.1 (requires config entry `bcftools_path=/path/to/bcftools`)
  - VarScan 2.4.1 (requires config entry `varscan_path=/path/to/VarScan.jar`)
  - RTG 3.6.1 (requires config entry `rtg_jar=/path/to/RTG.jar` AND `rtg_ref_sdf=/path/to/ref_genome.sdf` for SDF-formatted ref genome)
  - Platypus 0.8. (requires config entry `platypus_path=/path/to/Platypus.py`)
   
Similarly, the following variant comparison tools are all supported

  - [vcfeval](http://realtimegenomics.com/products/rtg-tools/) (requires config entry `rtg_jar=/path/to/RTG.jar` AND `rtg_ref_sdf=/path/to/ref_genome.sdf` for SDF-formatted ref genome)
  - [hap.py](https://github.com/Illumina/hap.py) (requires config entry `happy_path=/slc-ngs/opt/hap.py-0.2.7/bin/hap.py`)
  - [vgraph](https://github.com/bioinformed/vgraph) (requires config entry `vgraph_path=vgraph`)

###Adding new callers, normalizers, etc

Adding new callers, normalization tools, or comparators is easy and doesn't require a rebuild. To add a new caller, simply create a new python file containing a function that that executes the caller on a given BAM file and returns the VCF. That file also must define a function called `get_callers` that returns a dict mapping a unique string to the caller's function. Here's a quick example:

    import vcomp.util 
    
    def call_variant_freebayes(bam, ref_fasta, bed, conf):
        vcfoutput = "output-fb.vcf"
        cmd=["/path/to/freebayes", "-f", ref_fasta, "-t", bed, "-b", bam, "-v", vcfoutput]
        subprocess.check_output(cmd)
        return util.sort_vcf(vcfoutput, conf) # Automatically sort & bgzip & tabix index 
    
    def get_callers():
        return { "freebayes": call_variant_freebayes }

The new function (`call_variant_freebayes` in the example above) should have the same signature as in the example. the `bam`, `ref_fasta` and `bed` args are paths to those objects on the filesystem. The `conf` argument is an instance of a standard python [ConfigParser](https://docs.python.org/2/library/configparser.html) that contains information about the current configuration file. 
Lastly, add a line to the main configuration file under the 'callers' section stating where to find your Python file, like so:

    [callers]
    fb_caller=/path/to/freebayes_caller.py

 That should be it. By default a typical `varcomp` run will find and execute your newly added caller.  
 
 ##Controlling callers
 
 By default varcomp will execute every caller it finds. If you'd like to run some subset of the callers, use the --callers argument, like so:
 
     python vcomp/injectvar.py -v my_vars.vcf --callers freebayes
     


## Docker based setup

There are three Docker files present in this repository:
1. Dockerfile-vcomp-base which sets up all of the third-party tools used in the benchmakrs by default
2. Dockerfile-vcomp which adds the varcomp code onto the image produced Dockerfile-vcomp-base
3. Dockerfile which builds the entire image on its own (i.e. Dockerfile-vcomp-base + Dockerfile-vcomp)

First and second Dockerfiles are intended for development usage, as building the image with the third-party tools takes a bit of time.
Thus, to build the docker image start by cloning this repository:

        git clone https://github.com/goranrakocevic/varcomp.git
        cd varcomp

Building with the default set of options requires access to GATK JAR file, and the Platypus source folder (by default, they will be sought in the varcomp folder, so you should either copy/link them there, or modify the appropriate dockerfile to point to the right path)

You can then build the image either by calling:

        docker build -t vcomp-base:v1 -f Dockerfile-vcomp-base .
        docker build -t vcomp:v1 -f Dockerfile-vcomp .

Or:
    
        docker build -t vcomp:v1 -f Dockerfile .

File comp-docker.conf holds the configurations for all of the tools, while you may need to adjust the input files.
You should be ready  start varcomp:

        docker run -v /path/to/data:/data   \
                   -v /path/to/out/dir:/out \
                   -w /out                  \
                   vcomp:v1                 \
                   python /opt/varcomp/vcomp/injectvar.py -v /data/test.single.vcf --het -c /data/comp-docker.conf


