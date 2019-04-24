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
SAMPLES_FIELD = 13

#hugo name of gene to match with exon & gene query results
#e.g. CD99
gene_to_find = sys.argv[1]
decoderF = open(gene_to_find+".decoder.txt","w")
#this is the pseudo # of samples
#represents the semi-random partitioning 
#of the true samples into bins which will
#be given to JS as the "samples"
num_bins = int(sys.argv[2])
#get groups if passed in
groups_ = sys.argv[3]
groups_fh = {}
group_set = set()
groups = groups_.split(',')
#map groups to bin idx
group2bin = {}
for (i,g) in enumerate(groups):
    group2bin[g] = 0
    for j in range(0,num_bins):
        kg = g+str(j)
        groups_fh[kg] = [open(kg+".counts.tsv","w"),i] 
        decoderF.write(kg+"\t"+g+"\n")

gene = ""
idx = 0

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
    ftype = "A"
    if key in map2group:
        idx_str = map2group[key][1]
        ftype = map2type[key]
    #output the coordinates to the GFF
    else:
        idx+=1
        idx_str = str(idx)
        if idx < 10:
            idx_str = "0"+idx_str
        map2group[key] = [set(),idx_str]
        #exon line
        ftype = fields[1][-1];
        #gene or exon
        if ftype == 'G' or ftype == 'E':
            gene = fields[GENE_FIELD].split(':')[GENE_NAME_SUBFIELD]
            if gene != gene_to_find:
                del map2group[key]
                continue
            str_type = "exonic_part"
            if ftype == 'G':
                ftype = 'A'
                str_type = "aggregate_gene"
                #reset idx
                idx -= 1
                idx_str = "00"
                map2group[key][1] = idx_str
            sys.stdout.write("\t".join([chrm,"ScalaUtils",str_type,start,end,".",strand,".","gene_id "+gene+"; tx_set tx1; num 0"+idx_str+"; gene_set "+gene])+"\n")
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
    bin_idx = 0
    bin_sums = [0]*num_bins
    for (sid,cov) in [f.split(':') for f in fields[SAMPLES_FIELD].split(',')[1:]]:
        bin_sums[bin_idx] += int(cov)
        #track which counts groups we're writing to for this jx/exon
        bin_idx += 1
        if bin_idx >= num_bins:
            bin_idx = 0
    for (b,s) in enumerate(bin_sums):
        gk = group+str(b)
        fh = groups_fh[gk][0]
        fh.write(gene+":"+ftype+"0"+idx_str+"\t"+str(s)+"\n")
        map2group[key][0].add(gk)
        group_set.add(gk)
    map2type[key] = ftype

#write 0's for the exons/jxs that didn't have a count for a particular sample group file
for key in map2group.keys():
    (group_set_, idx_str) = map2group[key]
    ftype = map2type[key]
    for g in (group_set - group_set_):
        (fh,i) = groups_fh[g]
        fh.write(gene+":"+ftype+"0"+idx_str+"\t0\n")

for (fh,i) in groups_fh.values():
    fh.close()
decoderF.close()
