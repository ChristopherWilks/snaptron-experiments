#!/bin/env python3
#calculate basic statistics from base raw coverage results from the Snaptron client non-streaming using NP and related
import sys
import argparse
import numpy as np
from scipy import stats

def slope(indices, rows):
    lm = stats.linregress(indices, rows)
    return lm.slope

ops={'sum':np.sum, 'mean':np.mean, 'median':np.median, 'slope':slope}
def main(args):
    #which column do the base coverge values start in?
    #if doing a sample subset, there's only one initial non-value column
    starting_col = args.base_start_col
    if args.subset:
        starting_col = 1
    header = not args.noheader
    rows = []
    for line in sys.stdin:
        #skip header
        if header:
            header = False
            continue
        fields = line.rstrip().split('\t')
        row = np.array([int(n) for n in fields[starting_col:]])
        rows.append(row)
    if len(rows) > 0:
        b = np.array(rows)
        indices = range(0,len(b))
        row_length = len(b[0])
        stats = []
        if args.base_stat == 'slope':
            stats = [slope(indices, b[:,i]) for i in range(0,row_length)]
            max_val = np.max(b, axis=0)
            if args.normalize:
                stats = np.divide(np.array(stats), max_val)
            if args.include_max:
                sys.stderr.write("%s\t%s\n" % (args.label,'\t'.join([str(x) for x in max_val]))) 
        else:
            stats = ops[args.base_stat](b, axis=0)
        sys.stdout.write("%s\t%s\n" % (args.label,'\t'.join([str(x) for x in stats]))) 


def create_parser(disable_header=False):
    parser = argparse.ArgumentParser(description='Base coverage statistics')
    parser.add_argument('--subset', action='store_const', const=True, default=False, help='use if Snaptron client was called with --samples ...')
    parser.add_argument('--base-stat', metavar='cross-base summary(s) operation to perform', type=str, default='mean', help='how to summarize across all rows')
    parser.add_argument('--base-start-col', metavar='column index starting from 0', type=int, default=4, help='column index where base values start')
    parser.add_argument('--label', metavar='label to preface the output stats', type=str, default="", help='first field to be output to describe the range in input')
    parser.add_argument('--noheader', action='store_const', const=True, default=False, help='if there is no header in the in put data')
    parser.add_argument('--normalize', action='store_const', const=True, default=False, help='additionally compute factor(s) to make cross-sample comparisons')
    parser.add_argument('--include-max', action='store_const', const=True, default=False, help='additionally write maximum values per sample across rows to STDERR (useful for post-facto normalizion)')
    return parser

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    main(args)
    
