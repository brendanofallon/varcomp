import sys

from   collections import deque
from   itertools   import groupby
from   operator    import attrgetter
from   heapq       import heappop, heappush

import pysam

DEFAULT_CONTIG_ORDER=['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20',
                      '21', '22', '3', '4', '5', '6', '7', '8','9', 'MT', 'X','Y']

def var_comp(v1, v2):
    v1c = DEFAULT_CONTIG_ORDER.index(v1[0])
    v2c = DEFAULT_CONTIG_ORDER.index(v2[0])
    if v1c == v2c:
        return v1[1] - v2[1]
    else:
        return v1c - v2c


def window(seq, n=2):
    '''windowed iterator

    >>> map(tuple, window('ABCD'))
    [('A', 'B'), ('B', 'C'), ('C', 'D')]
    '''
    it = iter(seq)
    win = deque((next(it, None) for _ in xrange(n)), maxlen=n)
    yield win
    append = win.append
    for e in it:
        append(e)
        yield win


def decorate_vars(vars, d=1000):
    '''Generator to decorate VCF records vars, yielding:
       (contig, max(0, start-d), stop+d, var)
    '''
    for var in vars:
        yield var.contig, max(0, var.start - d), var.stop + 1000, var



def collect_solitary(items):
    '''
    Collect solitary (non-overlapping) decorated variants.  These are
    extracted before partitioning in order to later stuff into the smallest
    partitions to balance sizes
    '''
    solitary, popular = [], []

    nxt = None
    last_stop = -1


    for here, nxt in window(items):
        if nxt is None or (here[1] >= last_stop and here[2] <= nxt[1]):
            solitary.append(here)
        else:
            popular.append(here)
        last_stop = here[2]

    # Handle last item
    if nxt is not None:
        if nxt[1] >= last_stop:
            solitary.append(nxt)
        else:
            popular.append(nxt)

    return solitary, popular


def partition_intervals(items):
    '''
    Partition decorated variants into a minimal number of non-overlapping sets
    '''
    parts = [ (-1, []) ]

    for item in items:
        contig, start, stop, var = item

        if parts[0][0] <= start:
            _, p = heappop(parts)
            p.append(item)
        else:
            p = [item]

        heappush(parts, (stop, p))

    return [ p for (_, p) in parts ]


def merge_partitions(parts, solitary):
    '''
    Merge partitioned and solitary variants.  Partitioned variants are
    merged using a min-size greedy algorithm.  Solitary variants are added
    sequentially to the minimum sized partitions.
    '''

    # Start with an empty set of partitions.  For each contig, add each new
    # partition to the current set of partitions such that len(combined) =
    # max( len(old), len(new) ) in descending size order of new partitions
    # and selecting the smallest old partition.  These two greedy choices
    # are not optimal, since partitions are added contig by contig and not
    # globally.

    heap = []
    for part in parts:
        next_heap = []

        part.sort(key=len, reverse=True)

        for vars in part:
            if heap:
                _, p = heappop(heap)
                p.extend(vars)
            else:
                p = vars
            heappush(next_heap, (len(p), p) )

        for item in heap:
            heappush(next_heap, item)

        heap = next_heap

    # Stuff in solitary elements to the smallest partitions even out sizes
    for s in solitary:
        _, p = heappop(heap)
        p.append(s)
        heappush(heap, (len(p), p) )

    # Undecorate heap to obtain final partitions
    return [p for _,p in heap]


def main(args):
    filename = args[1]
    vars = pysam.VariantFile(filename)

    solitary = []
    parts = []

    for contig, contig_vars in groupby(vars, attrgetter('contig')):
        contig_vars = decorate_vars(contig_vars)
        contig_solitary, contig_popular = collect_solitary(contig_vars)
        contig_parts = partition_intervals(contig_popular)

        solitary += contig_solitary
        parts.append(contig_parts)

    parts = merge_partitions(parts, solitary)

    if filename == '-':
        basename = 'vars_part{:03d}.vcf'
    elif filename.endswith('.vcf') or filename.endswith('.bcf'):
        basename = filename[:-4] + '_part{:03d}' + filename[-4:]
    elif filename.endswith('.vcf.gz'):
        basename = filename[:-7] + '_part{:03d}' + filename[-7:]

    for i,part in enumerate(parts, 1):
        out_vars = pysam.VariantFile(basename.format(i), 'w', header=vars.header)
        for var in sorted(list(part), cmp=var_comp):
            out_vars.write(var[3])
        out_vars.close()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main(sys.argv)
