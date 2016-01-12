import sys

import injectvar

results_by_method = {}
results = {}
vgraph_vcfeval_mismatches = []
all_methods = []
all_comp_methods = []

def formatted(results, tot, result):
    if tot==0.0:
        return "NA\t"
    if result in results:
        if tot==1:
            return str(results[result]) + "\t"
        else:
            return "{:.5}".format(100.0 * results[result] / tot) + "\t"
    else:
        return "0\t"

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

        toks=line.strip().split(' ', 12)
        if len(toks)<10:
            continue
        ref = toks[5]
        alt = toks[6].replace(":", "")
        caller = toks[7]
        norm_method = toks[9]
        comp_method = toks[11]
        result = toks[12]
        thisvar = toks[2] + ":" + toks[3] + ":" + ref + ":" + alt

        if thisvar != prev_var:
            vgraph_vcfeval_mismatches.extend( find_vgraph_vcfeval_mismatches(var_results, prev_var))
            var_results = {}
            prev_var = thisvar

        if not norm_method in all_methods:
            all_methods.append(norm_method)

        if not comp_method in all_comp_methods:
            all_comp_methods.append(comp_method)

        if not caller in var_results:
            var_results[caller] = {}

        if not norm_method in var_results[caller]:
            var_results[caller][norm_method] = {}

        var_results[caller][norm_method][comp_method] = result

        if not caller in results:
            results[caller] = {}

        if not norm_method in results[caller]:
            results[caller][norm_method] = {}

        if not comp_method in results[caller][norm_method]:
            results[caller][norm_method][comp_method] = {}

        resultmap = results[caller][norm_method][comp_method]

        if result in resultmap:
            resultmap[result] += 1
        else:
            resultmap[result] = 1


for caller in results:
    print "\n" + caller
    for result in sorted(injectvar.all_result_types):
        print "\t" + result,
    print ""
    for norm_method in sorted(all_methods):
        for comp_method in sorted(all_comp_methods):
            print norm_method + "/" + comp_method + "\t",
            # tot = 0.0
            # for result in sorted(injectvar.all_result_types):
            #     if result in results[caller][norm_method][comp_method]:
            #         tot += results[caller][norm_method][comp_method][result]
            for result in sorted(injectvar.all_result_types):
                print formatted(results[caller][norm_method][comp_method], 1, result),
            print "\n",


print "\nCaller comparison: (=matches / (matches + not matched + no vars found + incorrect genotypes))"

comp_methods = ("vgraph:", "vcfeval:")
norm_method = "nonorm"
for comp_method in comp_methods:
    print "\t" + comp_method,
print ""

for caller in results:
    print caller + "\t",
    for comp_method in comp_methods:
        # tot = 0.0
        # for result in sorted(injectvar.all_result_types):
        #     if result in results[caller][norm_method][comp_method]:
        #         tot += results[caller][norm_method][comp_method][result]
        print formatted(results[caller][norm_method][comp_method], 1, injectvar.MATCH_RESULT),
    print "\n",

print "\nCaller err breakdown: mismatch / extra allele / missing allele / no variant found"

norm_method = "nonorm"
comp_method = "vgraph:"
print "\t" + norm_method
for caller in results:
    print caller + "\t",
    # tot = 0.0
    # for result in sorted(injectvar.all_result_types):
    #     if result in results[caller][norm_method][comp_method]:
    #         tot += results[caller][norm_method][comp_method][result]
    print formatted(results[caller][norm_method][comp_method], 1, injectvar.NO_MATCH_RESULT),
    print formatted(results[caller][norm_method][comp_method], 1, injectvar.ZYGOSITY_EXTRA_ALLELE),
    print formatted(results[caller][norm_method][comp_method], 1, injectvar.ZYGOSITY_MISSING_ALLELE),
    print formatted(results[caller][norm_method][comp_method], 1, injectvar.NO_VARS_FOUND_RESULT),
    print "\n",

# print "\n\nMatcher comparison: (=matches / (matches + not matched))"
# for norm_method in sorted(all_methods):
#     print "\t" + norm_method,
# print "\n",
#
# for caller in results:
#     print caller + "\t",
#     for norm_method in sorted(all_methods):
#         tot = 0.0
#         for result in sorted(injectvar.var_result_types):
#             if result in results[caller][norm_method]:
#                 tot += results[caller][norm_method][result]
#         print formatted(results[caller][norm_method], 1, injectvar.MATCH_RESULT),
#     print "\n",


print "\n\n vgraph / vcfeval mismatches:"
for mismatch in vgraph_vcfeval_mismatches:
    print mismatch[1] + "\t" + mismatch[0] + "\t" + " vgraph: " + mismatch[2] + "  vcfeval: " + mismatch[3]