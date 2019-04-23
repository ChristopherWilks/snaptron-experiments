#!/usr/bin/env python2
import sys

GROUP_FIELD = 0
GENE_FIELD = 12
GENE_NAME_SUBFIELD = 1
CHRM_FIELD = 3
END_FIELD = 5
STRAND_FIELD = 7
ANNOTATION_FIELD = 8
#coverage_sum, per sample group per exon/jx
TOTAL_COUNTS_FIELD = 15 

groups_ = sys.argv[1]
groups = groups_.split(',')
group_set = set(groups)
groups_fh = {}
for (i,g) in enumerate(groups):
    groups_fh[g] = [open(g+".counts.tsv","w"),i] 

gene = ""
idx = 0
gchrm = ""
#gstart = sys.maxint
#gend = -1
gstart = "-1"
gend = "-1"
gstrand = "?"

map2group = {}
#track type status of a mapped element (exon/annotation jx/novel jx)
map2type = {}

for line in sys.stdin:
    line = line.rstrip()
    fields = line.split('\t')
    group = fields[GROUP_FIELD]
    (chrm, start, end) = fields[CHRM_FIELD:END_FIELD+1]
    strand = fields[STRAND_FIELD]
    key = "\t".join([chrm,start,end,strand])
    skip_gff = False
    if key in map2group:
        idx_str = map2group[key][1]
    #output the coordinates to the GFF
    else:
        idx+=1
        idx_str = str(idx)
        if idx < 10:
            idx_str = "0"+idx_str
        map2group[key] = [set(),idx_str]
        #exon line
        ftype = "E"
        if fields[1][-1] == 'E':
            if len(gene) == 0:
                gene = fields[GENE_FIELD].split(':')[GENE_NAME_SUBFIELD]
                gchrm = chrm
                gstart = start
                gstrand = strand
            sys.stdout.write("\t".join([chrm,"ScalaUtils","exonic_part",start,end,".",strand,".","gene_id "+gene+"; tx_set tx1; num 0"+idx_str+"; gene_set "+gene])+"\n")
            gend = end
        else:
            #junction line (novel or annotated)
            ftype = "J"
            jx_type = "splice_site"
            tx_set = "tx1"
            gene_set = gene
            if fields[ANNOTATION_FIELD] == "0":
                ftype = "N"
                jx_type = "novel_" + jx_type 
                tx_set = "UNKNOWN_TX"
                gene_set = "UNKNOWN_GENE_SET"
            sys.stdout.write("\t".join([chrm,"ScalaUtils",jx_type,start,end,".",strand,".","gene_id "+gene+"; tx_set "+tx_set+"; num 0"+idx_str+"; gene_set "+gene_set])+"\n")
    #now output counts
    (fh,gidx) = groups_fh[group]
    fh.write(gene+":"+ftype+"0"+idx_str+"\t"+fields[TOTAL_COUNTS_FIELD]+"\n")
    #track which counts groups we're writing to for this jx/exon
    map2group[key][0].add(group)
    map2type[key] = ftype
#finall write out the whole gene line
sys.stdout.write("\t".join([gchrm,"ScalaUtils","aggregate_gene",gstart,gend,".",gstrand,".","gene_id "+gene+"; tx_set tx1; num 000; gene_set "+gene])+"\n")

#write 0's for the exons/jxs that didn't have a count for a particular sample group file
for key in map2group.keys():
    (group_set_, idx_str) = map2group[key]
    for g in (group_set - group_set_):
        (fh,i) = groups_fh[g]
        fh.write(gene+":"+ftype+"0"+idx_str+"\t0\n")

for (fh,i) in groups_fh.values():
    #TODO: write gene counts as 0 for now (change this to actually query the gene from snaptron)
    fh.write(gene+":A000\t0\n")
    fh.close()
