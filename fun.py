
import sys
import pandas as pd
from vcomp import batch_processor as bp
import numpy as np
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
from collections import defaultdict
import itertools as itools

ALL_CALLERS = ['samtools', 'varscan', 'platypus', 'rtg', 'gatk-hc', 'gatk-ug', 'freebayes', 'freebayes-mre']
SOME_CALLERS = ['gatk-hc', 'clcbio']
CLCBIO_CALLER = ['clcbio']
ALL_RESULTS = (bp.MATCH_RESULT, bp.NO_MATCH_RESULT, bp.NO_VARS_FOUND_RESULT, bp.MATCH_WITH_EXTRA_RESULT, bp.ZYGOSITY_MISSING_ALLELE, bp.ZYGOSITY_EXTRA_ALLELE, bp.ERROR_RESULT)

CALLER_SUBS = {
    'freebayes': 'Freebayes 1.0.1',
    'gatk-ug': 'GATK UG 3.5.0',
    'gatk-hc': 'GATK HC 3.5.0',
    'varscan': 'VarScan 2.4.1',
    'rtg': 'RTG 3.6',
    'freebayes-mre': 'Freebayes 1.0.2',
    'samtools': 'BCFTools 1.3.1',
    'platypus': 'Platypus 0.8',
    'clcbio': "CLCBio"
}

# A simple color palette
colors = [
    'dodgerblue',
    'orangered',
    'limegreen',
    'gold',
    'darkmagenta',
    'darkslateblue',
    'mediumorchid',
    'silver',
    'darkviolet',
    'cyan',
    'mediumspringgreen',
]

CLCBIO="CLCBio"
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

def freq_matches(d, callers=ALL_CALLERS):
    if callers is None:
        matches =  (d['result'] == bp.MATCH_RESULT).sum()
        total = float(len(d))
        return matches / total
    else:
        return [((d['caller']==c) & (d['result']==bp.MATCH_RESULT)).sum() / float((d['caller']==c).sum()) for c in callers]


def find_caller(name):
    if name in CALLER_SUBS:
        return CALLER_SUBS[name]
    else:
        if "clc" in name:
            return CLCBIO

def filename2label(fname):
    zyg = "?"
    depth = "?"
    if "het" in fname:
        zyg = "Heterozygous"
    else:
        zyg = "Homozygous"

    if "highdepth" in fname or "d200" in fname:
        depth = "depth 200"
    if "d30" in fname:
        depth = "depth 30"

    if "cis" in fname:
        return "SNP in cis"
    if "trans" in fname:
        return "SNP in trans"

    vtype = ""
    if "dels" in fname:
        vtype = "Deletions,"
    if "mnp" in fname:
        vtype = "MNVs,"
    if "ins" in fname:
        vtype = "Insertions,"
    if "dup" in fname:
        vtype = "Tandem Duplications,"

    return vtype + " " + zyg + ", " + depth

def filename2vartype(fname):
    vtype = ""
    if "dels" in fname:
        vtype = "Deletions"
    if "mnp" in fname:
        vtype = "MNVs"
    if "ins" in fname:
        vtype = "Insertions"
    if "dup" in fname:
        vtype = "Tandem Duplications"
    return vtype

def find_max_F1(data, caller):
    normalizer = NO_NORM
    comparator = VGRAPH
    qs = (data['caller'] == caller) & (data['normalizer'] == normalizer) & (data['comparator'] == comparator)
    matches = qs & (data['result'] == bp.MATCH_RESULT)
    mismatches = qs & (data['result'] == bp.NO_MATCH_RESULT)
    missed = qs & (data['result'] == bp.NO_VARS_FOUND_RESULT)

    max_F1 = 0.0
    max_qual = 0.0
    cutoffs = sorted( data['quality'][qs].unique() )

    cutoffs = cutoffs[1:25]
    max_prec = 0.0
    max_recall = 0.0
    for qual_cutoff in cutoffs:
        gt_qual = data['quality'] > qual_cutoff
        tps = (matches & gt_qual).sum()
        fps = (mismatches & gt_qual).sum()
        fns = (missed | (matches & (~gt_qual))).sum()
        if tps + fps == 0:
            continue

        if tps + fns == 0:
            continue

        precision = float(tps) / float(tps + fps)
        recall = float(tps) / float(tps + fns)
        F1 = 2.0 * precision * recall / (precision + recall)
        # print "Cutoff: {} tp: {} fp: {} fn: {} F1: {}".format(qual_cutoff, tps, fps, fns, F1)
        if F1 >= max_F1:
            max_F1 = F1
            max_qual = qual_cutoff
            max_prec = precision
            max_recall = recall
    return max_F1, max_qual, max_prec, max_recall

def plot_quals(data, callers=ALL_CALLERS):
    normalizer = NO_NORM
    comparator = VGRAPH
    if len(callers) == 1:
        cols = 1
        rows = 1
    else:
        cols = 2
        rows = max(1, len(callers) / cols + 1)

    max_var_size = 85
    sizes = data['variant'].map(var_size)
    for idx in range(rows * cols):
        if idx >= len(callers):
            break
        caller = callers[idx]
        ax = plt.subplot(rows, cols, idx + 1)
        qs = (data['caller']==caller) & (data['normalizer']==normalizer) & (data['comparator']==comparator) & (sizes < max_var_size)
        tot = qs.sum()
        matches = qs & (data['result']==bp.MATCH_RESULT)
        mismatches = qs & (data['result']==bp.NO_MATCH_RESULT)
        match_quals = data['quality'][matches]
        mismatch_quals = data['quality'][mismatches]
        # results[caller] = (match_quals, mismatch_quals)
        if len(match_quals)>5:
            sns.distplot(match_quals, ax=ax, label='Correct', kde_kws={'bw': 0.25})
        if len(mismatch_quals) > 5:
            sns.distplot(mismatch_quals, ax=ax, label='Incorrect', color='red', kde_kws={'bw': 0.25})
        ax.set_title(CALLER_SUBS[caller])
        maxqual = max(match_quals)
        if maxqual>1000:
            ax.set_xlim([0, 1000])
        else:
            ax.set_xlim([0, int( (maxqual*110.0)/100.0)])
        plt.locator_params(axis='y', nbins=4)
        if idx%2 == 0:
            ax.set_ylabel("Frequency")
        if idx > 5:
            ax.set_xlabel("Quality")
        else:
            ax.set_xlabel("")
        if idx==1:
            plt.legend(loc='upper right', fontsize='medium')

        if len(mismatch_quals)>10:
            uval, pval = stats.mannwhitneyu(match_quals, mismatch_quals)
        else:
            pval = 0.0
        max_F1, qual, prec, recall = find_max_F1(data, caller)
        print "\t".join([CALLER_SUBS[caller], str(tot), str(len(match_quals)), str(len(mismatch_quals)), "{:.3}".format(qual), "{:.3}".format(prec), "{:.3}".format(recall), "{:.4}".format(max_F1), "{:.5}".format(pval)])
        # plt.axvline(opt_q, linewidth=2.0, color='green')
    plt.tight_layout(h_pad=0.25, w_pad=1.0)
    plt.show()


def result_freq(data, callers=ALL_CALLERS):
    """
    Return proportion of
    :param data: DataFrame containing tabelized raw data
    :return:
    """
    normalizer = NO_NORM
    comparator = VGRAPH
    result = defaultdict(dict)
    for caller in callers:
        qs = (data['caller']==caller) & (data['normalizer']==normalizer) & (data['comparator']==comparator)
        qsum = float(qs.sum())
        for res_val in ALL_RESULTS:
            rsum = (qs & (data['result'] == res_val)).sum()
            result[caller][res_val] = rsum / qsum
    return result


def size_bin_str(idx, bins):
    if idx+1 >= len(bins):
        raise ValueError('nope')
    return "{}".format( (bins[idx]+bins[idx+1])/2 )



def size_bin_index(size, bins):
    if size < bins[0]:
        return None
    for l, u in zip(bins[:-1], bins[1:]):
        if size>l and size<=u:
            return bins.index(l)
    return None

def make_bin_indexer(bins):
    def index(b):
        return size_bin_index(b, bins)
    return index

def var_size(varstr):
    vars = varstr.split("/")
    toks = vars[-1].split()
    if len(toks)<5:
        raise ValueError("Could not parse variant tokens : " + varstr)

    ref = toks[3]
    alt = toks[4]
    if len(ref)==len(alt):
        return len(ref)
    else:
        return abs(len(alt) - len(ref))

def accuracy_by_size(data, callers=ALL_CALLERS, bins=range(1,150,20), result_type=bp.MATCH_RESULT):
    normalizer = NO_NORM
    comparator = VGRAPH
    result = defaultdict(list)
    sizes = data["variant"].map(var_size)
    size_bins = sizes.map(make_bin_indexer(bins))
    bin_strs = [size_bin_str(i, bins) for i in range(len(bins)-1)]
    for caller in callers:
        qs = (data['caller'] == caller) & (data['normalizer'] == normalizer) & (data['comparator'] == comparator)
        matches = qs & (data['result'] == result_type)
        cresults = []
        for bin in range(len(bins)):
            bq = size_bins == bin
            hits = (bq & matches).sum()
            all = (bq & qs).sum()
            if all == 0.0:
                cresults.append(np.nan)
            else:
                cresults.append(float(hits) / float(all))
        result[caller] = (bin_strs, cresults)
    return result


def plot_errors(ax, data, callers=ALL_CALLERS, colors=('red', 'orange', 'yellow',), result_types=(bp.NO_VARS_FOUND_RESULT, bp.NO_MATCH_RESULT, bp.ZYGOSITY_MISSING_ALLELE,), title=None, show_legend=True, ymax=None, show_ylab=True):
    totwidth = 0.8
    width = totwidth/1.0
    offset = -0.05
    prev = None
    results = result_freq(data, callers=callers)
    for j, res in enumerate(result_types):
        vals = [results[caller][res] for caller in callers]
        bars =plt.bar(np.arange(0-offset, len(callers)-offset), vals, bottom=prev, width=width, color=colors[j], label=res)
        if prev is None:
            prev = vals
        else:
            prev = [a+b for a,b in zip(prev, vals)]

    for bar in bars:
        ax.text(bar.get_x() + 0.1, prev[bars.index(bar)]+0.005, "{:1.2f}".format(prev[bars.index(bar)]), size='small', color='black')
    offset += totwidth/len(data)

    caller_labels = [CALLER_SUBS[caller]
                     if caller in CALLER_SUBS
                     else caller
                     for caller in callers]
    ax.set_xticklabels(caller_labels, rotation=75, horizontalalignment='left', size='small')
    if show_ylab:
        ax.set_ylabel("Fraction of all calls")
    ax.set_title(title)
    if ymax is not None:
        ax.set_ylim([0,ymax])
    if show_legend:
        ax.legend(loc=(1.0, 0.7), fontsize='medium')

def plot_vartypes(data, callers=ALL_CALLERS):
    fig = plt.figure()
    if len(data) == 1:
        cols = 1
        rows = 1
    else:
        cols = 2
        rows = len(data)/2
    callers = list(callers)
    result_types = (bp.MATCH_RESULT, bp.ZYGOSITY_MISSING_ALLELE, bp.NO_MATCH_RESULT, bp.NO_VARS_FOUND_RESULT)
    # result_types = (bp.MATCH_RESULT, )
    linestyles = ['-', '--', '-.', ':']
    for idx in range(rows*cols):
        name = data.keys()[idx]
        bins = range(0, 120, 10)
        ax = plt.subplot(rows, cols, idx + 1)

        for rtype in result_types:
            res = accuracy_by_size(data[data.keys()[idx]], bins=bins, callers=callers, result_type=rtype)

            for caller, acc in res.iteritems():
                ax.plot(acc[1],
                        label=find_caller(caller) + " - " + rtype,
                        linestyle=linestyles[result_types.index(rtype)],
                        color=colors[callers.index(caller)])

            print "Plotting caller {} -> {}".format(caller, find_caller(caller))
            ax.set_xticklabels(acc[0], fontsize='small')
            ax.set_title(filename2label(data.keys()[idx]))

        if idx in [0, 2]:
            ax.set_ylabel("Fraction of correct calls")

        ax.set_xlabel("Variant size")
        plt.locator_params(axis='x', nbins=len(acc[0]))


    legend = ax.legend(loc=(0.6, 0.65), fontsize='medium', frameon=True, framealpha=1.0)
    legend.get_frame().set_facecolor('white')

    plt.tight_layout()
    plt.show()


def plot_overall(names, data, callers=ALL_CALLERS):
    fig = plt.figure()
    cols = min(len(names), 2)
    rows = max(1, len(names)/2)
    for idx, fname in enumerate(names):
        ax = plt.subplot(rows, cols, idx+1)
        plot_errors(ax, data[fname], callers=callers, title=filename2label(fname), show_legend=(idx==-1), ymax=0.25, show_ylab=(idx in [0,2]))

    plt.tight_layout()
    plt.show()

def collect_callers(data):
    return data['caller'].unique()

def compute_ms90(data, callers=ALL_CALLERS):
    fig = plt.figure()
    fig.set_size_inches(6, 5)
    sns.set(font_scale=1.1)
    if len(data) == 1:
        cols = 1
        rows = 1
    else:
        cols = 2
        rows = len(data) / 2

    inputfiles = list(sorted(data.keys(), key=lambda x: filename2vartype(x)))
    callers = list(callers)

    bins = range(0, 150, 10)
    accuracy = 0.90

    bar_group_width = 0.8
    bar_width = bar_group_width / len(inputfiles)

    rtype = bp.MATCH_RESULT

    vtype_sums = defaultdict(int)
    for idx, inputfile in enumerate(inputfiles):
        result_data = []
        ax = plt.subplot("111")
        label = filename2vartype(inputfile)
        offset = bar_group_width*(float(idx)/len(inputfiles)-0.5)
        print "\n" + filename2label(inputfile)
        res = accuracy_by_size(data[inputfile], bins=bins, callers=callers, result_type=rtype)

        for caller in res:
            bin_strs, cresults = res[caller]
            vals = list(itools.takewhile(lambda x: x[1]>accuracy, enumerate(cresults)))
            ms90 = int(bin_strs[len(vals)])
            result_data.append(ms90)
            vtype_sums[label] += ms90
            print "{},{}".format(CALLER_SUBS[caller], bin_strs[len(vals)])

        xleft = [1.0 + x+offset for x in range(len(callers))]

        ax.bar(xleft,
               result_data,
               bar_width,
               color=colors[inputfiles.index(inputfile)],
               label=label)
        ax.set_xticklabels([""] + map(lambda i: CALLER_SUBS[i] if i in CALLER_SUBS else i, callers), rotation=40, horizontalalignment='right', va='top')
        ax.set_ylim([0, 140])
        if idx%2==0:
            ax.set_ylabel("MS90 (bp)")
        leg = ax.legend(loc='upper left', fontsize='medium', frameon=True, framealpha=1.0)
        leg.get_frame().set_edgecolor('b')
        leg.get_frame().set_facecolor('w')

    for key, val in vtype_sums.iteritems():
        print "{},{}".format(key, float(val)/len(callers))

    plt.tight_layout()
    plt.show()




def compute_venn(data, caller, normalizer_a, normalizer_b, comparator_a, comparator_b):
    all_a = (data['caller'] == caller) & (data['normalizer'] == normalizer_a) & (data['comparator'] == comparator_a)
    match_a = (all_a & (data['result'] == bp.MATCH_RESULT))
    tot_a_vars = set(data['variant'][all_a])
    match_a_vars = data['variant'][match_a]
    avars = set(match_a_vars)

    # for v in match_a_vars:
    #     x = match_a_vars == v
    #     if x.sum()>1:
    #         print "whoa, {} hits for {}".format(x.sum(), v)


    all_b = (data['caller'] == caller) & (data['normalizer'] == normalizer_b) & (data['comparator'] == comparator_b)
    tot_b_vars = set(data['variant'][all_b])
    match_b = (all_b & (data['result'] == bp.MATCH_RESULT))
    match_b_vars = data['variant'][match_b]
    bvars = set(match_b_vars)

    print "A: Found {} total and {} matches".format(len(tot_a_vars), len(avars))
    print "B: Found {} total and {} matches".format(len(tot_b_vars), len(bvars))
    isect = avars.intersection(bvars)
    uniq_a = avars - isect
    uniq_b = bvars - isect
    int_size = len(isect)
    print "Match intersection: {}".format(int_size)
    print "Unique A: {}".format(len(uniq_a))
    print "Unique B: {}".format(len(uniq_b))


def main(args):
    data = {}
    for arg in args:
        data[arg] = pd.read_csv(arg, sep='\t')
    # compute_venn(data[args[0]], 'gatk-ug', NO_NORM, NO_NORM, VCFEVAL, VGRAPH)
    # compute_ms90(data, callers=collect_callers(data[args[0]]))
    plot_vartypes(data, callers=collect_callers(data[args[0]]))
    # plot_overall(args, data, callers=collect_callers(data[args[0]])) # SOME_CALLERS)
    # plot_quals(data[args[0]])

if __name__=="__main__":
    main(sys.argv[1:])
