#!/usr/bin/env bash

#e.g. srav2
compilation=$1
#e.g. coordinates of cassette exon: chr8:117146819-117146848
exon=$2
#e.g. coordinates of junction spanning/excluding cassette exon: chr8:117135399-117146953
exclusion=$3
#need strand for splicing: e.g. "+"
strand=$4
min_count=$5

#first determine 2 inclusion junctions and write them out to a snaptron query file
outf=`perl -e '$o="'$strand'"; $comp="'$compilation'"; $inc="'$exon'"; $ex="'$exclusion'"; ($c,$s,$e)=split(/[:-]/,$inc); $s--; $e++; ($c2,$s2,$e2)=split(/[:-]/,$ex); open(OUT,">$comp.$c.$s.$e.psi.snap.tsv"); print OUT "region\texact\tfilters\tgroup\n$c:$s2-$s\t1\tstrand=$o\tINCLUSION\n$c:$e-$e2\t1\tstrand=$o\tINCLUSION\n$ex\t1\tstrand=$o\tEXCLUSION\n"; close(OUT); print "$comp.$c.$s.$e.psi.snap.tsv";'`
#extract the 2 inclusion junction coordinates
inc1=`grep INCLUSION $outf | cut -f 1 | head -1`
inc2=`grep INCLUSION $outf | cut -f 1 | tail -n1`

#get raw stats from Snaptron per junction
cat <(echo -n "INCLUSION_1\t") <(curl "http://snaptron.cs.jhu.edu/${compilation}/snaptron?regions=${inc1}&exact=1&rfilter=strand:+" | cut -f 1-12,14-) > data/diabetes.${compilation}.tsv
cat <(echo -n "INCLUSION_2\t") <(curl "http://snaptron.cs.jhu.edu/${compilation}/snaptron?regions=${inc2}&exact=1&rfilter=strand:+&header=0" | cut -f 1-12,14-) >> data/diabetes.${compilation}.tsv
cat <(echo -n "EXCLUSION\t") <(curl "http://snaptron.cs.jhu.edu/${compilation}/snaptron?regions=${exclusion}&exact=1&rfilter=strand:+&header=0" | cut -f 1-12,14-) >> data/diabetes.${compilation}.tsv

#now get PSI for all samples and summarize across all samples
python2 client/query_snaptron.py --query-file $outf --function psi --datasrc ${compilation} --min-count $min_count --summarize > ${outf}.psi_out.tsv 2> ${outf}.psi_summary.tsv

head -1 $outf > ${outf}.jir
#only want one inclusion junction for the JIR
#avoids double counting
tail -n2 $outf >> ${outf}.jir

#now get JIR for all 
python2 client/query_snaptron.py --query-file ${outf}.jir --function jir --datasrc ${compilation} --min-count-jir $min_count --summarize > ${outf}.jir_out.tsv 2> ${outf}.jir_summary.tsv
