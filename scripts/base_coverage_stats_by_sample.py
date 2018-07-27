#!/bin/env python2.7
#calculate basic statistics from base raw coverage results from the Snaptron client
import sys
import argparse


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

ops={'sum':sum, 'mean':mean, 'median':median}
reverse_ops={sum:'sum', mean:'mean', median:'median'}
def main(args):
    #which column do the base coverge values start in?
    #if doing a sample subset, there's only one initial non-value column
    starting_col = args.base_start_col
    if args.subset:
        starting_col = 1
    header = not args.noheader
    final_per_sample = []
    num_rows=0.0
    for line in sys.stdin:
        fields = line.rstrip().split('\t')
        #skip header
        if header:
            header = False
            sys.stdout.write('operation\t'+'\t'.join(fields[args.base_start_col:])+'\n')
            continue
        num_rows += 1
        #init 
        if len(final_per_sample) == 0:
            final_per_sample = [0]*len(fields)
        final_per_sample = [sum([float(x),final_per_sample[j]]) for (j,x) in enumerate(fields[starting_col:])]
   
    
    if 'sum' in args.final_stat:
        if args.label is not None:
            sys.stdout.write("%s\t" % args.label)
        sys.stdout.write('sum\t'+'\t'.join([str(x) for x in final_per_sample])+'\n')
    if 'mean' in args.final_stat:
        if args.label is not None:
            sys.stdout.write("%s\t" % args.label)
        sys.stdout.write('mean\t'+'\t'.join([str(x/num_rows) for x in final_per_sample])+'\n')

def create_parser(disable_header=False):
    parser = argparse.ArgumentParser(description='Base coverage statistics')
    parser.add_argument('--subset', action='store_const', const=True, default=False, help='use if Snaptron client was called with --samples ...')
    parser.add_argument('--final-stat', metavar='cross-sample summary(s) operation to perform', type=str, default='sum', help='how to summarize the data *per row*, "mean", "median", "sum", or 2 or more of them (e.g. mean,sum, or both)?')
    parser.add_argument('--noheader', action='store_const', const=True, default=False, help='if there is no header in the in put data')
    parser.add_argument('--base-start-col', metavar='column index starting from 0', type=int, default=4, help='column index where base values start')
    parser.add_argument('--label', metavar='label to preface the output stats', type=str, default=None, help='first field to be output to describe the range in input')
    return parser

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    main(args)
    
