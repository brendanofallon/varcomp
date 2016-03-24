
import json
import injectvar, batch_processor
import sys
from collections import defaultdict
# import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import itertools

VARIANT = "variant"
QUALS="caller_quals"
BAMSTATS = "bamstats"
RESULTS = "results"

GATK_HC="gatk-hc"
RTG="rtg"
FREEBAYES="freebayes"
VGRAPH = "vgraph"
HAPPY = "happy"
VAP_LEFT = "vapleft"
VT = "vt"
VCFEVAL = "vcfeval"
NO_NORM = "nonorm"
RAW_COMP = "raw"

LABEL_SUBS = {
    batch_processor.NO_VARS_FOUND_RESULT: "Ref. call",
    batch_processor.ZYGOSITY_MISSING_ALLELE: "Zygosity error",
    batch_processor.NO_MATCH_RESULT: "Other mismatch",
    batch_processor.MATCH_RESULT: "Match",
}

class Tabelize(object):
    """
    Writes tab-delimited results in table form to sys.stdout
    """
    def __init__(self):
        self.first = True

    def perform_op(self, results):
        if self.first:
            print "\t".join(["variant", "quality", "totalreads", "properpairs", "mq20", "mq40", "softclipped_reads", "softclipped_bases", "caller", "normalizer", "comparator", "result"])
            self.first = False
        var = results[VARIANT]
        for caller, caller_results in results[RESULTS].iteritems():
            qual = results[QUALS][caller]
            for normalizer, norm_results in caller_results.iteritems():
                for comparator, comp_results in norm_results.iteritems():
                    totreads = results[BAMSTATS]['total_reads']
                    pp = results[BAMSTATS]['properpair']
                    if 'mq20' in results[BAMSTATS]:
                        mq20 = results[BAMSTATS]['mq20']
                    else:
                        mq20 = 0
                    if 'mq40' in results[BAMSTATS]:
                        mq40 = results[BAMSTATS]['mq40']
                    else:
                        mq40 = 0
                    if 'softclipped_reads' in results[BAMSTATS]:
                        screads = results[BAMSTATS]['softclipped_reads']
                    else:
                        screads = 0

                    scbases = results[BAMSTATS]['softclipped_bases']

                    print "\t".join([var, str(qual), str(totreads), str(pp), str(mq20), str(mq40), str(screads), str(scbases), caller, normalizer, comparator, comp_results])

    def finalize(self):
        pass


class NormBreakFinder(object):
    """
    Finds results in which use of VT or V.A.P / leftalign caused vgraph to mismatch, but 'nonorm' was match
    """

    def __init__(self):
        self.breaks = {}
        self.comp_method = VGRAPH
        self.norm_method1 = VAP_LEFT
        self.nonorm_method = NO_NORM

    def perform_op(self, results):
        for caller in results[RESULTS]:
            try:
                 if results[RESULTS][caller][self.nonorm_method][self.comp_method] == batch_processor.MATCH_RESULT and results[RESULTS][caller][self.norm_method1][self.omp_method] != batch_processor.MATCH_RESULT:
                     self.breaks[results[VARIANT] + "-" + caller] = self.norm_method1 + ": " + results[RESULTS][caller][self.norm_method1][self.comp_method] + "  " + self.nonorm_method + ": " + results[RESULTS][caller][self.nonorm_method][self.comp_method]
            except:
                pass

    def finalize(self):
        print "Normalization (" + self.norm_method1 + ") causing comparator (" + self.comp_method + ") mismatches:"
        if len(self.breaks)==0:
            print "\tNo normalizer breaking issues detected"
        else:
             for var, result in self.breaks:
                 print var + "\t" + result


class VAPFailsVgraphHits(object):

    def __init__(self):
        self.vap_hits = defaultdict(dict)
        self.vt_hits = defaultdict(dict)
        self.naive_hits = defaultdict(dict)
        self.vgraph = VGRAPH
        self.vapleft = VAP_LEFT
        self.vt = VT
        self.nonorm = NO_NORM
        self.rawcomp = RAW_COMP
        self.tot = defaultdict(int)
        self.callers = [GATK_HC]

    def perform_op(self, results):
        var_results = results[RESULTS]
        vartype = get_vartype(results[VARIANT])

        caller = self.callers[0]
        try:
            vt_result = var_results[caller][self.vt][self.rawcomp]
            vap_result = var_results[caller][self.vapleft][self.rawcomp]
            graph_result = var_results[caller][self.nonorm][self.vgraph]
            naive_result = var_results[caller][self.nonorm][self.rawcomp]

            # if graph_result != batch_processor.MATCH_RESULT:
            #     return
            # if vap_result != batch_processor.MATCH_RESULT:
            #     return

            if vt_result != graph_result:
                self.vt_hits[vartype][results[VARIANT]] = caller + " " + self.vt + ": " + vt_result + "  " + self.vgraph + ": " + graph_result
            if vap_result != graph_result:
                self.vap_hits[vartype][results[VARIANT]] = caller + " " + self.vapleft + ": " + vap_result + "  " + self.vgraph + ": " + graph_result
            if naive_result != graph_result:
                self.naive_hits[vartype][results[VARIANT]] = caller + " " + self.nonorm + ": " + naive_result + "  " + self.vgraph + ": " + graph_result
            self.tot[vartype] += 1
        except KeyError:
            pass

    def finalize(self):
        print "VAP / raw match fails, matched by " + self.vgraph

        tothits = 0.0
        tottot = 0.0
        for vartype in sorted(self.tot):
            tothits += len(self.vap_hits[vartype])
            tottot += self.tot[vartype]
            print vartype + ": ",
            for var, result in self.vap_hits[vartype].iteritems():
                print var + "\t" + result
            print vartype + "\t" + str(len(self.vap_hits[vartype])/float(self.tot[vartype])) #+ "\t" + str(len(self.vap_hits[vartype])) + "\t" + str(self.tot[vartype])
        print "VAP / LeftAlign overall:\t" + str(tothits / tottot) + "\t" + str(tottot)

        tothits = 0.0
        tottot = 0.0
        print "\n\nVT / raw match fails, matched by " + self.vgraph
        for vartype in sorted(self.tot):
            tothits += len(self.vt_hits[vartype])
            tottot += self.tot[vartype]
            print vartype + ": ",
            for var, result in self.vt_hits[vartype].iteritems():
               print var + "\t" + result
            print vartype + "\t" + str(len(self.vt_hits[vartype])/float(self.tot[vartype])) #+ "\t" + str(len(self.vt_hits[vartype])) + "\t" + str(self.tot[vartype])
        print "VT:\t" + str(tothits / tottot) + "\t" + str(tottot)

        tothits = 0.0
        tottot = 0.0
        print "\n\nNaive match fails, matched by " + self.vgraph
        for vartype in sorted(self.tot):
            tothits += len(self.naive_hits[vartype])
            tottot += self.tot[vartype]
            #print vartype + ": ",
            # for var, result in self.naive_hits[vartype].iteritems():
            #     print var + "\t" + result
            print vartype + "\t" + str(len(self.naive_hits[vartype])/float(self.tot[vartype])) #+ "\t" + str(len(self.naive_hits[vartype])) + "\t" + str(self.tot[vartype])
        print "Naive:\t" + str(tothits / tottot)+  "\t" + str(tottot)

        print "\n\nTotal variants by type:"
        for vartype, tot in self.tot.iteritems():
            print vartype + " : " + str(tot)

class CallerSummary(object):

    def __init__(self):
        self.summary = {}
        self.comparator = VGRAPH
        self.normalizer = NO_NORM


    def perform_op(self, results):
        var_results = results[RESULTS]
        for caller in var_results:
            if caller not in self.summary:
                self.summary[caller] = defaultdict(int)
            cresult = var_results[caller][self.normalizer][self.comparator]
            self.summary[caller][cresult] += 1

    def finalize(self):
        print "Caller summary:"
        print "caller\t" + "\t".join(injectvar.all_result_types)
        for caller in self.summary:
            tot = 0
            for result in self.summary[caller]:
                tot += self.summary[caller][result]
            print caller,
            for res in injectvar.all_result_types:
                print "\t{:.5}".format(100.0*float(self.summary[caller][res])/float(tot)),
            print ""

        # plot_callers( ([self.summary], ), ["Overall"])
        plot_results(self.summary)

class CallerSummaryBySize(object):

    def __init__(self, breaks=None):
        self.comparator = VGRAPH
        self.normalizer = NO_NORM
        if breaks is None:
            self.breaks = [10, 25, 50, 1000]
        else:
            self.breaks = breaks
        self.ins_summary = [{} for _ in range(len(self.breaks))]
        self.del_summary = [{} for _ in range(len(self.breaks))]

    def _index(self, size):
        for i, el in enumerate(self.breaks):
            if size < el:
                return i
        return None

    def perform_op(self, results):
        var_results = results[RESULTS]

        var = results[VARIANT].split("/")[-1]
        toks = var.split()
        ref = toks[3]
        alt = toks[4]
        if len(ref) == len(alt):
            idx = self._index(len(ref))
        else:
            idx = self._index( abs( len(ref)-len(alt)))
        if idx is None:
            return

        if len(alt) > len(ref):
            summary = self.ins_summary[idx]
        else:
            summary = self.del_summary[idx]

        for caller in var_results:
            if caller not in summary:
                summary[caller] = defaultdict(int)
            cresult = var_results[caller][self.normalizer][self.comparator]
            summary[caller][cresult] += 1

    def _emit_summary(self, summary):
        print "caller\t" + "\t".join(injectvar.all_result_types)
        for caller in summary:
            tot = 0
            for result in summary[caller]:
                tot += summary[caller][result]
            print caller + "(" + str(tot) + ")",
            for res in injectvar.all_result_types:
                print "\t{:.5}".format(100.0*float(summary[caller][res])/float(tot)),
            print ""

    def finalize(self):
        print "Caller summary:"
        bstrs = [str(i) + "-" + str(j) for i,j in zip([0] + self.breaks[0:-1], self.breaks)]
        for i, bstr in enumerate(bstrs):
            print "\nInsertions, Size range: " + bstr
            self._emit_summary(self.ins_summary[i])

        for i, bstr in enumerate(bstrs):
            print "\nDeletions, Size range: " + bstr
            self._emit_summary(self.del_summary[i])

        #plot_callers( (self.ins_summary, self.del_summary), bstrs)
        # if len(self.ins_summary)>0:
        #      plot_callers_line(self.ins_summary, [str(x) for x in self.breaks + [str(self.breaks[-1]) + "+" ]])
        # else:

        plot_callers_line(self.ins_summary, [str(x) for x in self.breaks + [str(self.breaks[-1]) + "+" ]])


class AccuracyBySoftclip(object):
    """
    Compute match percentage as a function of number of softclipped reads
    """

    def __init__(self):
        self.stat = "softclipped_reads"
        self.bins = [x/20.0 for x in range(10)]
        self.matches = {}
        self.totals = {}
        self.norm = NO_NORM
        self.comp = VGRAPH

    def perform_op(self, results):
        var_results = results[RESULTS]
        vartype = get_vartype(results[VARIANT])
        if vartype != "Insertion (10-20)":
            return

        try:
            stat = float(results[BAMSTATS][self.stat]) / float(results[BAMSTATS]['total_reads'])
        except KeyError:
            stat = 0
        bin = self._index(stat)
        if bin is None:
            return
        for caller in var_results:
            if caller not in self.totals:
                self.totals[caller] = [0 for _ in range(len(self.bins))]
                self.matches[caller] = [0 for _ in range(len(self.bins))]

            self.totals[caller][bin] += 1

            res1 = var_results[caller][self.norm][self.comp]
            if res1 == batch_processor.MATCH_RESULT:
                self.matches[caller][bin] += 1


    def finalize(self):
        for caller in self.totals:
            print caller + "\t",
        print ""
        for i,l,u in zip(range(len(self.bins)), [0] + self.bins[0:-1], self.bins):
            print str(u),

            for caller in self.totals:
                mval = float(self.matches[caller][i])
                totval = float(self.totals[caller][i])
                if mval==0.0:
                    print "\t0.0",
                else:
                    print "\t" + "{:3g}".format(mval/totval),
            print ""


    def _index(self, val):
        for j, bin in enumerate(self.bins):
            if val < bin:
                return j
        return None




class GraphCompMismatches(object):

    def __init__(self):
        self.mismatches = {}
        self.comp1 = VGRAPH
        self.comp2 = HAPPY
        self.comp3 = VCFEVAL

    def perform_op(self, results):
        var_results = results[RESULTS]
        for caller in var_results:
            for norm_method in var_results[caller]:
                res1 = var_results[caller][norm_method][self.comp1]
                res2 = var_results[caller][norm_method][self.comp2]
                res3 = var_results[caller][norm_method][self.comp3]

                if res1 != res2 or res2 != res3:
                    self.mismatches[results[VARIANT]] = caller + "/ " + norm_method + ": " + self.comp1 + ":" +res1 + "\t" + self.comp2 + ": " +res2 + "\t" + self.comp3 + ": " + res3

    def finalize(self):
        print "Graph comparator mismatches:"
        if len(self.mismatches)==0:
            print "\tNone detected"
        else:
            for k,v in self.mismatches.iteritems():
                print k + "\t" + v



def plot_results(data):
    fig = plt.figure()
    fig.patch.set_facecolor('white')
    # sns.set(style="white")
    xmod = 0.5
    xlocs = []
    labels = []
    ax = plt.subplot(211)
    ax.yaxis.grid(True)
    vals = {}

    results = [batch_processor.NO_MATCH_RESULT, batch_processor.ZYGOSITY_MISSING_ALLELE, batch_processor.NO_VARS_FOUND_RESULT]

    colors = ('green', 'orange')
    ymax = 0

    tots = defaultdict(int)
    for caller in data:
        for res in [batch_processor.MATCH_RESULT, batch_processor.NO_VARS_FOUND_RESULT, batch_processor.NO_MATCH_RESULT, batch_processor.ZYGOSITY_MISSING_ALLELE]:
            tots[caller] += data[caller][res]

    #First, plot the matches on the first plot
    v = []
    labels = []
    for caller in data:
        v.append(data[caller][batch_processor.MATCH_RESULT] / float(tots[caller]))
        labels.append(caller)
    vals[caller] = v
    prev = v
    bars = ax.bar(range(len(v)), v, color=colors[0], label="Correct call")
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width()/3.5, 0.925*bar.get_height(), "{:.3n}".format(bar.get_height()), size='small', color='black')
    v = []
    labels = []
    for caller in data:
        v.append(sum([data[caller][batch_processor.NO_MATCH_RESULT], data[caller][batch_processor.ZYGOSITY_MISSING_ALLELE], data[caller][batch_processor.NO_VARS_FOUND_RESULT]]))
        labels.append(caller)
    vals[caller] = v

    ax.bar(range(len(v)), v, bottom=prev, color=colors[1], label="Incorrect call")


    ax.set_ylim([0,1.0])
    ax.set_ylabel("Fraction of correct calls")
    ax.set_xticklabels(labels, rotation=45, horizontalalignment='left')
     # ax.set_xticks([x+0.3 for x in xlocs])
    ax.legend(loc='lower right', fontsize='small')


    colors = ('red', 'orange', 'yellow')
    prev = None
    ax = plt.subplot(212)
    ax.yaxis.grid(True)
    tots = defaultdict(int)
    for caller in data:
        for res in [batch_processor.NO_VARS_FOUND_RESULT, batch_processor.NO_MATCH_RESULT, batch_processor.ZYGOSITY_MISSING_ALLELE]:
            tots[caller] += data[caller][res]

    for j, res in enumerate([batch_processor.NO_VARS_FOUND_RESULT, batch_processor.NO_MATCH_RESULT, batch_processor.ZYGOSITY_MISSING_ALLELE]):
        v = []
        labels = []
        for i, caller in enumerate(data):
            v.append(data[caller][res] / float(tots[caller]))
            labels.append(caller)
        vals[caller] = v

        label = res
        if label in LABEL_SUBS:
            label = LABEL_SUBS[res]

        bars = ax.bar(range(len(v)), v, bottom=prev, color=colors[j], label=label)
        for bar in bars:
            if bar.get_height() > 0:
                height_offset = bar.get_height()-0.068 if bar.get_height() > 0.06 else bar.get_height()-0.05
                if prev is not None:
                    height_offset += prev[bars.index(bar)]
                    if prev[bars.index(bar)] > 0.95:
                        height_offset = 1.02
                    # if height_offset-prev[bars.index(bar)] < 0.1:
                    #     height_offset = prev[bars.index(bar)]+0.02
                # height_offset = min(0.925, height_offset)
                ax.text(bar.get_x() + bar.get_width()/3.5, height_offset, "{:.2g}".format(bar.get_height()), size='small', color='black')
        if prev is None:
            prev = v
        else:
            prev = [a+b for a,b in zip(prev, v)]

    ax.set_ylim([0,1.0])
    ax.set_ylabel("Fraction of erroneous calls by type")
    ax.set_xticklabels(labels, rotation=45, horizontalalignment='left')
     # ax.set_xticks([x+0.3 for x in xlocs])
    ax.legend(loc='right', fontsize='small')

    plt.tight_layout()
    plt.savefig("calls.png", dpi=200)
    plt.show()

def plot_callers_line(data, titles):
    fig = plt.figure()
    fig.patch.set_facecolor('white')
    cmap = cm.get_cmap('CMRmap')
    markers = ['.', 'o', 'v', '+', 's', '*', '>', 'x', 'D', '<', '^']
    target_result = batch_processor.MATCH_RESULT


    #Need to accumulate a single list of accuracies for a given caller across sizes
    accuracies = defaultdict(list)
    for i, res_set in enumerate(data):
        for caller in res_set:
            tot = sum(res_set[caller].values())
            accuracies[caller].append(float(res_set[caller][target_result]) / float(tot))

    ax = plt.subplot(111)
    c = 0.0
    xlocs = range(len(titles))
    for caller, vals in accuracies.iteritems():
        plt.plot(xlocs[0:len(vals)], vals, linewidth=2.0, color=cmap(c/float(len(accuracies))), marker=markers[int(c%len(markers))], label=caller)
        c += 1


    ax.legend(loc='lower right')
    ax.set_xticks(range(len(titles)))
    ax.set_xticklabels(titles)
    ax.set_ylabel("Accuracy (matches / total assessments)")
    ax.set_xlabel("Size")
    plt.show()



def plot_callers(data, titles):
    fig = plt.figure()
    fig.patch.set_facecolor('white')
    rows = len(data)
    cols = len(data[0])
    xmod = 0.5
    ymax = 0
    results = [batch_processor.MATCH_RESULT, batch_processor.NO_MATCH_RESULT,batch_processor.ZYGOSITY_MISSING_ALLELE, batch_processor.NO_VARS_FOUND_RESULT]
    for c, ins in enumerate( itertools.chain(*data)):
        labels = []
        ax = plt.subplot(rows, cols, c+1)
        prev = None
        colors = ('green', 'yellow', 'orange', 'red',)
        for j, res in enumerate(results):
            v = []
            labels = []
            for i, caller in enumerate(ins):
                v.append(ins[caller][res])
                labels.append(caller)

            label = res
            if label in LABEL_SUBS:
                label = LABEL_SUBS[res]
            vmax = max(v) if len(v)>0 else 0
            if vmax > ymax:
                ymax = vmax
            ax.bar(range(len(v)), v, width=0.5, bottom=prev, color=colors[j], label=label)
            if prev is None:
                prev = v
            else:
                prev = [a+b for a,b in zip(prev, v)]

        ax.set_ylim([0,ymax])
        ax.yaxis.grid(True)
        # bars = ax.bar(xlocs, vals, 0.8)
        ax.set_xticklabels(labels, size='small', rotation=45, horizontalalignment='left')
        # ax.set_xticks([x+0.3 for x in xlocs])
        # for bar in bars:
        #     ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(), str(int(bar.get_height())), ha='center', va='bottom', size='small')

        ax.set_title(titles[c % cols])

    plt.subplots_adjust(hspace=0.3)
    plt.show()



def parseline(line):
    """
    Convert the line of input into a results dict
    :param line:
    :return:
    """
    return json.loads(line)


def sizebin(size):
    bins = [-30, -20, -10, 0, 5, 10, 20, 50]
    if size < bins[0]:
        return "<" + str(bins[0])
    for l, u in zip(bins[:-1], bins[1:]):
        if size>l and size<=u:
            return str(l) + "-" + str(u)
    return ">" + str(bins[-1])

def get_vartype(varstr):
    vars = varstr.split("/")
    toks = vars[-1].split()
    ref = toks[3]
    alt = toks[4]
    if len(ref)==1 and len(alt)==1:
        return "SNP"
    if len(ref)<2 and len(alt)>1:
        return "Insertion (" + sizebin(len(alt)-len(ref)) + ")"
    if len(alt)<2 and len(ref)>1:
        return "Deletion (" + sizebin(len(ref)-len(alt)) + ")"
    return "MNP " + sizebin( len(ref))

def main(path, operations=[]):
    line_num = 0
    with open(path) as fh:
        for line in fh.readlines():
            line_num += 1
            if len(line)==0 or line[0] == '#':
                continue
            try:
                results = parseline(line)
            except Exception as ex:
                sys.stderr.write("Error parsing line #" + str(line_num) + ": " + str(ex))
                continue

            for op in operations:
                op.perform_op(results)


    for op in operations:
        print "\n"
        op.finalize()


if __name__=="__main__":
    if len(sys.argv)<2:
        print "Please enter the name of the results file to parse"
        exit(1)

    ops = [

        Tabelize(),
        #NormBreakFinder(),
        #VAPFailsVgraphHits(),
        #GraphCompMismatches(),
        #CallerSummary(),
        #CallerSummaryBySize(range(1, 151, 10)),
        #AccuracyBySoftclip()
    ]
    main(sys.argv[1], ops)