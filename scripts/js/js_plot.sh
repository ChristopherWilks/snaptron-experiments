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
#e.g. rfilter=coverage_sum>:1000&contains=1
filters=$4
#to use when labeling files
name=$5

if [[ -z $name ]]; then
    name=$comp
fi

if [[ -n $filters ]]; then
    filters="&${filters}"
fi

echo -n "" > ${name}.jxs
echo -n "" > ${name}.exons
for g in `cat <(echo $groups) | tr , " "`; do
    curl "$SNAPTRON/${comp}/snaptron?regions=${region}${filters}&sfilter=SMTS:${g}&header=0" | perl -ne 'print "'$g'\t$_";' >> ${name}.jxs
    curl "$SNAPTRON/${comp}/exons?regions=${region}&header=0" | perl -ne 'print "'$g'\t$_";' >> ${name}.exons
done

cat ${name}.exons ${name}.jxs | python2 format_for_js.py $groups > ${name}.gff
