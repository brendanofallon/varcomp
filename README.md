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
 
  
