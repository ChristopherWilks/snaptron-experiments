#!/bin/env python2.7
#calculate basic statistics from base raw coverage results from the Snaptron client
import sys
import argparse
from collections import Counter


def mean(vals):
    s = sum(vals)
    count = len(vals)
    return s/float(count)

def median(vals):
    count = len(vals)
    if count == 1:
        return float(vals[0])
    vals = sorted(vals)
    idx1 = int(count/2)
    idx2 = idx1 + 1
    #are we odd?
    if float(idx1) != count/2.0:
        idx1 = idx2
    return (vals[idx1]+vals[idx2])/2.0


ops={'sum':sum, 'mean':mean, 'median':median,'none':lambda x: x}
reverse_ops={sum:'sum', mean:'mean', median:'medians', ops['none']:'none'}
def main(args):
    #which column do the base coverge values start in?
    #if doing a sample subset, there's only one initial non-value column
    starting_col = args.base_start_col
    sample_ops = [ops[op] for op in args.sample_stat.split(',')]
    if args.subset:
        starting_col = 1
    header = not args.noheader
    base_counts = []
    no_sample_stats = args.sample_stat == 'none'
    for line in sys.stdin:
        #skip header
        if header:
            header = False
            continue
        fields = line.rstrip().split('\t')
        if no_sample_stats:
            base_counts.append([float(x) for x in fields[starting_col:]])
            continue
        stats = [sample_op([float(x) for x in fields[starting_col:]]) for sample_op in sample_ops]
        if not args.suppress_rows:
            if not args.row_labels:
                sys.stdout.write("sample %s\t%s\n" % (args.sample_stat, '\t'.join([str(stat) for stat in stats])))
            else:
                sys.stdout.write("%s\t%s\n" % ('\t'.join(fields[0:starting_col]), '\t'.join([str(stat) for stat in stats])))
        base_counts.append(stats)
    if not args.suppress_final and len(base_counts) > 0:
        if no_sample_stats:
            sys.stdout.write("%s\t" % (args.label))
            sums = [0]*len(base_counts[0])
            counts = Counter()
            for row in base_counts:
                counts.update({i:int(y) for (i,y) in enumerate(row)})
            sys.stdout.write("\t".join(str(counts[y]) for y in sorted(counts.keys()))+"\n")
        else:
            base_ops = [ops[op] for op in args.base_stat.split(',')]
            for op in base_ops:
                for (i,_) in enumerate(sample_ops):
                    stat = op([b[i] for b in base_counts])
                    sys.stdout.write("%s of %ss\t%s\n" % (str(reverse_ops[op]), str(reverse_ops[sample_ops[i]]), str(stat)))

def create_parser(disable_header=False):
    parser = argparse.ArgumentParser(description='Base coverage statistics')
    parser.add_argument('--subset', action='store_const', const=True, default=False, help='use if Snaptron client was called with --samples ...')
    parser.add_argument('--sample-stat', metavar='cross-sample summary(s) operation to perform', type=str, default='mean', help='how to summarize the data *per row*, "mean", "median", "sum", or 2 or more of them (e.g. median,mean,sum,none)?')
    parser.add_argument('--base-stat', metavar='cross-base summary(s) operation to perform', type=str, default='mean', help='how to summarize across all rows')
    parser.add_argument('--suppress-rows', action='store_const', const=True, default=False, help='should we suppress output of individual row numbers?')
    parser.add_argument('--suppress-final', action='store_const', const=True, default=False, help='should we suppress the final output of aggregate stats?')
    parser.add_argument('--row-labels', action='store_const', const=True, default=False, help='should we include the labels at the start of each row?')
    parser.add_argument('--noheader', action='store_const', const=True, default=False, help='if there is no header in the in put data')
    parser.add_argument('--base-start-col', metavar='column index starting from 0', type=int, default=4, help='column index where base values start')
    parser.add_argument('--label', metavar='label to preface the output stats', type=str, default="", help='first field to be output to describe the range in input')
    return parser

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    main(args)
    
