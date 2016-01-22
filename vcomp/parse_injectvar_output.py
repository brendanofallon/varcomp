import sys

import injectvar

results_by_method = {}
results = {}
vgraph_vcfeval_mismatches = []
normalizer_issues = []
vap_fails = []
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
        for norm_method in var_results[caller]:
            if "vgraph:" in var_results[caller][norm_method] and "vcfeval:" in var_results[caller][norm_method] and var_results[caller][norm_method]["vcfeval:"] != injectvar.NO_VARS_FOUND_RESULT:
                if var_results[caller][norm_method]["vgraph:"] != var_results[caller][norm_method]["vcfeval:"]:
                    mismatches.append( (caller, norm_method, var, var_results[caller][norm_method]["vgraph:"], var_results[caller][norm_method]["vcfeval:"]) )
    return mismatches


def find_normalizer_breaks(var_results, var):
    breaks = []
    comp_method = "vgraph:"
    for caller in var_results:
        try:
            res = None
            res_method = None
            for norm_method in var_results[caller]:
                if res is None:
                    res = var_results[caller][norm_method][comp_method]
                    res_method = norm_method
                if var_results[caller][norm_method][comp_method] != res:
                    breaks.append( (caller, norm_method, var, var_results[caller][norm_method][comp_method], res_method + "=" +res) )
        except:
            pass
    return breaks


def find_vapfails_vgraph_hits(var_results, var):
    fixes = []
    for caller in var_results:
        try:
            vap_result = var_results[caller]["vapleft"]["raw:"]
            vcfeval_result = var_results[caller]["nonorm"]["vcfeval:"]
            if vap_result != vcfeval_result:
                fixes.append( (caller, var, "vap/leftalign:" + vap_result, "vcfeval:" + vcfeval_result))
        except KeyError:
            pass
    return fixes




#Build big dict with all results, organized by caller and then method and then result type

var_results = {} #Stores all results for a single variant, helps comparing results for a single var
prev_var = "-1"
for line in open(sys.argv[1], "r"):
    if line.startswith("Result for"):

        toks=line.strip().split(' ', 12)
        if len(toks)<12:
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
            normalizer_issues.extend( find_normalizer_breaks(var_results, prev_var))
            vap_fails.extend( find_vapfails_vgraph_hits(var_results, prev_var))
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
                try:
                    print formatted(results[caller][norm_method][comp_method], 1, result),
                except:
                    pass
            print "\n",


print "\nCaller comparison: (=matches / (matches + not matched + no vars found + incorrect genotypes))"

comp_methods = ("vcfeval:",)
norm_method = "nonorm"
for comp_method in comp_methods:
    print "\t" + comp_method,
print ""

for caller in results:
    print caller + "\t",
    for comp_method in comp_methods:
        tot = 0.0
        for result in sorted(injectvar.all_result_types):
            if result in results[caller][norm_method][comp_method]:
                tot += results[caller][norm_method][comp_method][result]
        print formatted(results[caller][norm_method][comp_method], tot, injectvar.MATCH_RESULT),
    print "\n",

print "\nCaller err breakdown: mismatch / extra allele / missing allele / no variant found"

norm_method = "nonorm"
comp_method = "vcfeval:"
print "\t" + norm_method
for caller in results:
    print caller + "\t",
    tot = 0.0
    for result in sorted(injectvar.all_result_types):
        if result in results[caller][norm_method][comp_method]:
            tot += results[caller][norm_method][comp_method][result]
    print formatted(results[caller][norm_method][comp_method], tot, injectvar.NO_MATCH_RESULT),
    print formatted(results[caller][norm_method][comp_method], tot, injectvar.ZYGOSITY_EXTRA_ALLELE),
    print formatted(results[caller][norm_method][comp_method], tot, injectvar.ZYGOSITY_MISSING_ALLELE),
    print formatted(results[caller][norm_method][comp_method], tot, injectvar.NO_VARS_FOUND_RESULT),
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
    print mismatch[1] + "\t" + mismatch[0] + "\t" + mismatch[2] + " vgraph: " + mismatch[3] + "  vcfeval: " + mismatch[4]

print "\n\n Normalizer changing vgraph results:"
for case in normalizer_issues:
    print "\t".join(list(case))

print "\n\n VAP/leftalign fails, fixed by vcfeval:"
for case in vap_fails:
    print "\t".join(list(case))
