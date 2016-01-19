# varcomp
Tools for calling and comparing variants from simulated read sets

This tool consists of two main modules. The first, `vcomp/injectvar.py` processes an input VCF of variants and writes results to standard output (you'll probably want to redirect it to a file). The second, `vcomp/parse_injectvar_output.py` parses the results file and emits aggregated results. Typical usage looks something like:

    python vcomp/injectvar.py -v my_variants.vcf > my_output.txt
  
Each line in `my_output.txt` is the result for a single caller / normalizer / comparator combo. The format should be pretty self-explanatory, although the list of all possible result types is steadily changing. 

 To process the output to produce some aggregated results, do something like:
 
     python vcomp/parse_injectvar_output.py my_output.txt
     
 Various aggregate statistics will be written to standard output. Exactly what is computed is in flux and is likely to remain so for some time. 
 
 varcomp makes extensive use of multiple external applications (callers, normalizers, aligners, comparison tools, etc). The paths to these applications must be defined in a configuration file called `comp.conf`. By default, `injectvar.py` looks for a file called `comp.conf` in the active directory, but you can specify a path to it by using the `-c /path/to/configuration/file` argument. The format of the file is pretty straightforward, just a list of key=value pairs where the values are the paths to various executables. 
 
# Under the hood
 
 While the module `injectvar` contains basic processing structure, the callers, normalizers, and comparators are defined in the `callers.py`, `normalizers.py` and `comparators.py` modules. Each caller / normalizer / comparator is simply a python callable, and these are stored in a dict that associates each callable with a name. In each module, there's a function called something like `get_callers()` that returns this dictionary. There's no command line argument control over callers or anything - so adjustments to the callers / normalizers / etc. must be done by commenting out lines within the `get_callers()` etc functions. 
 
# Temporary and flagged dirs

 For each variant, `injectvar` creates a temporary working directory and keeps all of the VCFs, BAMs, etc in there. These are named something like `tmp-working129873` (the numbers may change). These are usually deleted when the variant is done being processed, but in some cases it's nice to leave them there for closer investigation. For instance, if vgraph / vcfeval / hap.py disagree, we want to know exactly what happened and be able to examine each VCF. To facilitate this, certain types of results can trigger a temporary directory to be 'flagged' and preserved. In this case, the name of the temporary directory is changed to something like `tmpdir-chr10-128376-comp-conflict` and a text file called `flag.info.txt` is created in the directory that summarizes what happened. (This also happens for errors, but in this case the file is called `exception.info.txt`.)
