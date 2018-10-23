#!/bin/env python2.7
import sys
import re
import argparse

START_OF_COUNTS_COL=5
MIN_EXON_SZ=1000

def load_annotations(args):
    annots = set()
    with open(args.annotation_file,"rb") as fin:
        for line in fin:
            (chrm,start,end,gene) = line.rstrip().split('\t')[:4]
            annots.add("%s:%s:%s:%s" % (gene,chrm,start,end))
    return annots


def process_counts(args, annots):
    with open(args.counts_file,"rb") as fin:
        #fairly big file, do one line at a time
        (pg1,pbp,pchrm,pstart,pend) = (None,0,None,0,0)
        counts = []
        for line in fin:
            line = line.rstrip()
            #remove decimals which are useless for counts
            if args.remove_decimals:
                line = line.replace('.0\t','\t')[:-2]
            fields = line.split('\t')
            (g,bp,chrm,start,end) = fields[:START_OF_COUNTS_COL]
            g1=g.split(args.delimiter)[0]
            if args.join_coords_for_label:
                g = ':'.join([g,chrm,start,end])
            if start == 'start':
                sys.stdout.write(line+"\n")
                continue
            #if this is not a mini-exon, continue as normal
            if g in annots:
                if pg1 is not None:
                    sys.stdout.write("%s\t%d\t%s\t%d\t%d\t" % (pg1,pbp,pchrm,pstart,pend))
                    sys.stdout.write("\t".join([str(x) for x in counts])+"\n")
                sys.stdout.write(line+"\n")
                (pg1,pbp,pchrm,pstart,pend) = (None,0,None,0,0)
                counts = []
                continue
            #mini-exon but we have a previous split gene to finish
            elif pstart == 0 or pg1 != g1:
                if pg1 is not None:
                    sys.stdout.write("%s\t%d\t%s\t%d\t%d\t" % (pg1,pbp,pchrm,pstart,pend))
                    sys.stdout.write("\t".join([str(x) for x in counts])+"\n")
                pstart = int(start)
                pend = 0
                pbp = 0
                counts = [0]*len(fields[START_OF_COUNTS_COL:])
            pg1=g1
            pbp+=int(bp)
            pchrm=chrm
            #in case the coordinates are not in order for these splits
            #this can happen since the sort is non-numeric
            start = int(start)
            if start < pstart:
                pstart=start
            end = int(end)
            if end > pend:
                pend=end
            counts = [counts[i]+int(x) for (i,x) in enumerate(fields[START_OF_COUNTS_COL:])] 
                
        if pg1 is not None:
            sys.stdout.write("%s\t%d\t%s\t%d\t%d\t" % (pg1,pbp,pchrm,pstart,pend))
            sys.stdout.write("\t".join([str(x) for x in counts])+"\n")

def create_parser(disable_header=False):
    parser = argparse.ArgumentParser(description='Fix up artificially split intervals for large coverage jobs')
    parser.add_argument('--annotation-file', metavar='/path/to/annotation_file', type=str, default=None, help='file contains list of chrom,start,end,labels for intervals')
    parser.add_argument('--counts-file', metavar='/path/to/counts_file', type=str, default=None, help='file contains group,length,chrom,start,end,count0,count1,...countN')
    parser.add_argument('--join-coords-for-label', action='store_const', const=True, default=False, help='join the label,chrom,start,end to match annotation label')
    parser.add_argument('--remove-decimals', action='store_const', const=True, default=False, help='remove decimal and right hand 0s')
    parser.add_argument('--delimiter', metavar=':', type=str, default=':', help='delimiter of label + split index (default ":")')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    annots = load_annotations(args)

    process_counts(args, annots)

if __name__ == '__main__':
    main()
