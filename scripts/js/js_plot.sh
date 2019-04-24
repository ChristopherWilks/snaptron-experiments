#!/usr/bin/env bash
#set -eox

SNAPTRON='http://snaptron.cs.jhu.edu'

#compilation
#e.g. "gtex"
comp=$1
#normal WSI snaptron query syntax:
#novel exon from ASCOT for gene IMPDH1 (exon boundary)
#e.g. chr7:128393029-128394277
region=$2
#e.g. "Pituitary,Brain,Spleen"
#must be at least 2 groups
groups=$3
#gene gname to match exon/gene queries
gname=$4
#e.g. rfilter=coverage_sum>:1000&contains=1
filters=$5

if [[ -n $filters ]]; then
    filters="&${filters}"
fi

echo -n "" > ${gname}.jxs
echo -n "" > ${gname}.exons
echo -n "" > ${gname}.genes
for g in `cat <(echo $groups) | tr , " "`; do
    curl "$SNAPTRON/${comp}/snaptron?regions=${region}${filters}&sfilter=SMTS:${g}&header=0" | perl -ne 'print "'$g'\t$_";' >> ${gname}.jxs
    curl "$SNAPTRON/${comp}/exons?regions=${region}&sfilter=SMTS:${g}&header=0" | perl -ne 'print "'$g'\t$_";' >> ${gname}.exons
    curl "$SNAPTRON/${comp}/genes?regions=${region}&sfilter=SMTS:${g}&header=0" | perl -ne 'print "'$g'\t$_";' >> ${gname}.genes
done

cat ${gname}.genes ${gname}.exons ${gname}.jxs | python2 format_for_js.py $gname 2 $groups > ${gname}.gff
cat <(echo 'sample.ID	group.ID') <(sort IMPDH1.decoder.txt) > IMPDH1.decoder.txt.sorted
