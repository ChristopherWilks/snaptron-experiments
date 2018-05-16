#!/bin/env python2.7
import sys

START_OF_COUNTS_COL=5
MIN_EXON_SZ=1000

def load_annotations(f):
    annots = set()
    with open(f,"rb") as fin:
        for line in fin:
            (chrm,start,end,gene) = line.rstrip().split('\t')[:4]
            #annots["%s:%s:%d-%d" % (g,chrm,start,end)] = 1
            annots.add("%s:%s:%s:%s" % (gene,chrm,start,end))
    return annots


def process_counts(annots, f):
    with open(f,"rb") as fin:
        #fairly big file, do one line at a time
        (pg,pg1,pbp,pchrm,pstart,pend) = (None,None,0,None,0,0)
        counts = []
        for line in fin:
            line = line.rstrip()
            fields = line.split('\t')
            (g,bp,chrm,start,end) = fields[:START_OF_COUNTS_COL]
            g1=g.split(':')[0]
            if start == 'start':
                sys.stdout.write(line+"\n")
                continue
            if g in annots:
                if pg is not None:
                    sys.stdout.write("%s:%s:%d:%d\t%d\t%s\t%d\t%d\t" % (pg1,pchrm,pstart,pend,pbp,pchrm,pstart,pend))
                    sys.stdout.write("\t".join([str(x) for x in counts])+"\n")
                sys.stdout.write(line+"\n")
                (pg,pg1,pbp,pchrm,pstart,pend) = (None,None,0,None,0,0)
                counts = []
                continue
            elif pstart == 0 or pg1 != g1:
                if pg is not None:
                    sys.stdout.write("%s:%s:%d:%d\t%d\t%s\t%d\t%d\t" % (pg1,pchrm,pstart,pend,pbp,pchrm,pstart,pend))
                    sys.stdout.write("\t".join([str(x) for x in counts])+"\n")
                pstart = int(start)
                pend = 0
                pbp = 0
                counts = [0.0]*len(fields[START_OF_COUNTS_COL:])
            pg=g
            pg1=g1
            pbp+=int(bp)
            pchrm=chrm
            pend=int(end)
            counts = [counts[i]+float(x) for (i,x) in enumerate(fields[START_OF_COUNTS_COL:])] 
            #reached final mini-exon in this real-exon group of lines
            if int(bp) != MIN_EXON_SZ:
                sys.stdout.write("%s:%s:%d:%d\t%d\t%s\t%d\t%d\t" % (pg1,pchrm,pstart,pend,pbp,pchrm,pstart,pend))
                sys.stdout.write("\t".join([str(x) for x in counts])+"\n")
                (pg,pg1,pbp,pchrm,pstart,pend) = (None,None,0,None,0,0)
                counts = []
                
        if pg is not None:
            sys.stdout.write("%s:%s:%d:%d\t%d\t%s\t%d\t%d\t" % (pg1,pchrm,pstart,pend,pbp,pchrm,pstart,pend))
            sys.stdout.write("\t".join([str(x) for x in counts])+"\n")

def main():
    gene_annotationF = sys.argv[1]
    gene_exon_countsF = sys.argv[2]

    annots = load_annotations(gene_annotationF)

    process_counts(annots, gene_exon_countsF)

if __name__ == '__main__':
    main()
