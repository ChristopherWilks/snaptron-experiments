#!/bin/env python2.7
#calculate basic statistics from base raw coverage results from the Snaptron client
import sys
import argparse


def mean(vals):
    s = sum([float(x) for x in vals])
    count = len(vals)
    return s/float(count)

def median(vals):
    count = len(vals)
    if count == 1:
        return float(vals[0])
    vals = sorted([float(x) for x in vals])
    idx1 = int(count/2)
    idx2 = idx1 + 1
    #are we odd?
    if float(idx1) != count/2.0:
        idx1 = idx2
    return (vals[idx1]+vals[idx2])/2.0


ops={'mean':mean,'median':median}
def main(args):
    #which column do the base coverge values start in?
    #if doing a sample subset, there's only one initial non-value column
    starting_col = 4
    sample_op = ops[args.sample_stat]
    if args.subset:
        starting_col = 1
    header = True
    base_counts = []
    for line in sys.stdin:
        #skip header
        if header:
            header = False
            continue
        fields = line.rstrip().split('\t')
        stat = sample_op(fields[starting_col:])
        if not args.suppress_rows:
            sys.stdout.write("sample %s\t%s\n" % (args.sample_stat, str(stat)))
        base_counts.append(stat)
    if len(base_counts) > 0:
        stat = ops[args.base_stat](base_counts)
        sys.stdout.write("aggregate coordinate range %s of cross-sample per-base coverage %s\t%s\n" % (args.base_stat, args.sample_stat, str(stat)))


def create_parser(disable_header=False):
    parser = argparse.ArgumentParser(description='Base coverage statistics')
    parser.add_argument('--subset', action='store_const', const=True, default=False, help='use if Snaptron client was called with --samples ...')
    parser.add_argument('--sample-stat', metavar='cross-sample summary operation to perform', type=str, default='median', help='how to summarize the data *per row*, "mean" or "median"?')
    parser.add_argument('--base-stat', metavar='cross-base summary operation to perform', type=str, default='median', help='how to summarize across all rows')
    parser.add_argument('--suppress-rows', action='store_const', const=True, default=False, help='should we suppress output of individual row numbers?')
    return parser

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    main(args)
    
