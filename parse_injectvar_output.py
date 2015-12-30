
import sys
import comparators

results_by_method = {}
results = {}
all_methods = []

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

def formatted(results, tot, result):
    if result in results:
        return "{:.5}".format(100.0 * results[result] / tot) + "\t"
    else:
        return "0.0\t"

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