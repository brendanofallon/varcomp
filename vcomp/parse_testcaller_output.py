
import sys

results_by_size = {}
results = {}
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
        if not caller in results:
            results[caller] = [0, 0]
        counts = results[caller]


        size = len(alt) - len(ref)
        if not size in results_by_size:
            results_by_size[size] = {}
        if not caller in results_by_size[size]:
            results_by_size[size][caller] = [0, 0]
        size_counts = results_by_size[size][caller]

        if 'Variants matched' in result:
            counts[0] += 1
            size_counts[0] += 1
        else:
            counts[1] += 1
            size_counts[1] += 1

#print "Totals by caller: "
for caller in results:
    counts = results[caller]
    percent = float(counts[0]) / (float(counts[0] + counts[1]))*100
 #   print caller + "\t:\t" + str(counts[0]) + " / " + str(counts[0] + counts[1]) + "\t" + str(percent)[0:6] + "%"

#print "Totals by size / caller"

print ",".join(["size", "caller", "hits", "attempts"])
for size in sorted(results_by_size.keys()):
    for caller in results_by_size[size]:
        counts = results_by_size[size][caller]
        percent = float(counts[0]) / (float(counts[0] + counts[1]))*100
        #print "" + caller + "".join([' ' for _ in range(max(0, 12-len(caller)))]) +  "\t:\t" + str(counts[0]) + " / " + str(counts[0] + counts[1]) + "\t" + str(percent)[0:6] + "%"
        print ",".join([str(size), caller, str(counts[0]), str(counts[0] + counts[1])])

