import sys

from vcomp import comparators

results_by_method = {}
results = {}
vgraph_vcfeval_mismatches = []
all_methods = []

def formatted(results, tot, result):
    if result in results:
        return "{:.5}".format(100.0 * results[result] / tot) + "\t"
    else:
        return "0.0\t"

def find_vgraph_vcfeval_mismatches(var_results, var):
    mismatches = []
    for caller in var_results:
        if "vgraph:" in var_results[caller] and "vcfeval:" in var_results[caller]:
            if var_results[caller]["vgraph:"] != var_results[caller]["vcfeval:"]:
                mismatches.append( (caller, var, var_results[caller]["vgraph:"], var_results[caller]["vcfeval:"]) )
    return mismatches

#Build big dict with all results, organized by caller and then method and then result type

var_results = {} #Stores all results for a single variant, helps comparing results for a single var
prev_var = "-1"
for line in open(sys.argv[1], "r"):
    if line.startswith("Result for"):

        toks=line.strip().split(' ', 10)
        if len(toks)<10:
            continue
        ref = toks[5]
        alt = toks[6].replace(":", "")
        caller = toks[7]
        method = toks[9]
        result = toks[10]
        thisvar = toks[2] + ":" + toks[3] + ":" + ref + ":" + alt

        if thisvar != prev_var:
            vgraph_vcfeval_mismatches.extend( find_vgraph_vcfeval_mismatches(var_results, prev_var))
            var_results = {}
            prev_var = thisvar

        if not caller in var_results:
            var_results[caller] = {}

        var_results[caller][method] = result

        if not caller in results:
            results[caller] = {}

        if not method in results[caller]:
            results[caller][method] = {}
        resultmap = results[caller][method]
        if not method in all_methods:
            all_methods.append(method)

        if result in resultmap:
            resultmap[result] += 1
        else:
            resultmap[result] = 1


for caller in results:
    print "\n\t" + caller
    for result in sorted(comparators.all_result_types):
        print "\t" + result,
    print ""
    for method in sorted(all_methods):
        print method + "\t",
        tot = 0.0
        for result in sorted(comparators.all_result_types):
            if result in results[caller][method]:
                tot += results[caller][method][result]
        for result in sorted(comparators.all_result_types):
            print formatted(results[caller][method], tot,  result),
        print "\n",


print "\nCaller comparison: (=matches / (matches + not matched + no vars found + incorrect genotypes))"

method = "vgraph:"
print "\t" + method
for caller in results:
    print caller + "\t",
    tot = 0.0
    method = "vgraph:"
    for result in sorted(comparators.all_result_types):
        if result in results[caller][method]:
            tot += results[caller][method][result]
    print formatted(results[caller][method], tot, comparators.MATCH_RESULT),
    print "\n",

print "\nCaller err breakdown: mismatch / incorrect genotype / no variant found"

method = "vgraph:"
print "\t" + method
for caller in results:
    print caller + "\t",
    tot = 0.0
    method = "vgraph:"
    for result in sorted(comparators.all_result_types):
        if result in results[caller][method]:
            tot += results[caller][method][result]
    print formatted(results[caller][method], tot, comparators.NO_MATCH_RESULT),
    print formatted(results[caller][method], tot, comparators.INCORRECT_GENOTYPE_RESULT),
    print formatted(results[caller][method], tot, comparators.NO_VARS_FOUND_RESULT),
    print "\n",

print "\n\nMatcher comparison: (=matches / (matches + not matched))"
for method in sorted(all_methods):
    print "\t" + method,
print "\n",

for caller in results:
    print caller + "\t",
    for method in sorted(all_methods):
        tot = 0.0
        for result in sorted(comparators.var_result_types):
            if result in results[caller][method]:
                tot += results[caller][method][result]
        print formatted(results[caller][method], tot, comparators.MATCH_RESULT),
    print "\n",


print "\n\n vgraph / vcfeval mismatches:"
for mismatch in vgraph_vcfeval_mismatches:
    print mismatch[1] + "\t" + mismatch[0] + "\t" + " vgraph: " + mismatch[2] + "  vcfeval: " + mismatch[3]