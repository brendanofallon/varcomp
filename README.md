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

Due to a current limitation with the simulation procedure, the minimum distance between any two variants in the input VCF cannot be less than 2kb. Otherwise, nearby variants might interfere in an unpredictable manner. 

##Variant caller benchmarking

Varcomp can also take an input variant list, simulate reads (or use existing simulated fastqs), and run several variant callers, and compare the results to the input set of variants to see which caller had the higher accuracy. Basic usage looks like:

    python vcomp/injectvar.py -v my_variants.vcf > my_output.txt
    
If you've already created an input set of fastqs, you can use them by: 

    python vcomp/injectvar.py -v my_variants.vcf --fqs my_variants_R1.fastq --fqs my_variants_R2.fastq > my_output.txt
    
###Modifying input variants

Given an input VCF, varcomp can perform several modifications of the variants prior to simulating reads. In particular, zygosities can be forced to het or hom, and extra SNPs can be added upstream of each input variant, either in cis or trans. 

    --het        Force all input variants to be heterozygous
    --hom        Force all input variants to be homozygous
    --addsnp     Add upstream a single SNP upstream of each input variant (default, variants added in cis)
    --addsnp --trans Add upstream SNP in trans
   
One of --het or --hom is *required* if the input VCF does not contain a GT format entry for each variant.


##Configuration
 
 varcomp makes extensive use of multiple external applications (callers, normalizers, aligners, comparison tools, etc). The paths to these applications must be defined in a standard python configuration file. By default, `injectvar.py` looks for a file called `comp.conf` in the active directory, but you can specify a path to it by using the `-c /path/to/configuration/file` argument. The format of the file is pretty straightforward, just a list of key=value pairs where the values are the paths to various executables. 
 
 TODO: Exactly what is required  in a minimal configuration file, samtools, tabix, and bgzip? 
 
##Adding new callers, normalizers, etc

Adding new callers, normalization tools, or comparators is easy and doesn't require a rebuild. To add a new caller, simply create a new python file containing a function that that executes the caller on a given BAM file and returns the VCF. That file also must define a function called `get_callers` that returns a dict mapping a unique string to the caller's function. Here's a quick example:

    import vcomp.util 
    
    def call_variant_freebayes(bam, ref_fasta, bed, conf):
        vcfoutput = "output-fb.vcf"
        cmd=["/path/to/freebayes", "-f", ref_fasta, "-t", bed, "-b", bam, "-v", vcfoutput]
        subprocess.check_output(cmd)
        return util.sort_vcf(vcfoutput, conf)
    
    def get_callers():
        return { "freebayes": call_variant_freebayes }

The new function (`call_variant_freebayes` in the example above) should have the same signature as in the example. the `bam`, `ref_fasta` and `bed` args are paths to those objects on the filesystem. The `conf` argument is an instance of a standard python ConfigParser[https://docs.python.org/2/library/configparser.html] that contains information about the current configuration file. 
Lastly, add a line to the main configuration file under the 'callers' section stating where to find your Python file, like so:

    [callers]
    fb_caller=/path/to/freebayes_caller.py

 That should be it. By default a typical `varcomp` run will find and execute your newly added caller.  



