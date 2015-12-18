
import sys

results_by_method = {}
results = {}
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

        if result in resultmap:
            resultmap[result] += 1
        else:
            resultmap[result] = 1


for caller in results:
    print "Caller: " + caller
    for method in results[caller]:
        print "\t" + method
        for result, count in results[caller][method].iteritems():
            print "\t\t" + result + "\t:\t" + str(count)