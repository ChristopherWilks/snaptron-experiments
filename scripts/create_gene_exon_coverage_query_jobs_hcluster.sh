#!/bin/bash
#creates list of parallel jobs to construct gene/exon coverage sums per all samples
#in a Snaptron compilation using a specific annotation

#1) parameter is the TSV of disjoint exon coordinates, format(tab-delimted):
#chrm,start,end,gene_name
#2) parameter is the Snaptron compilation name string (e.g. encode1159)
#3) parameter is the number of splits you want

#get this script's path
p=`perl -e '@f=split("/","'$0'"); pop(@f); print "".join("/",@f)."\n";'`

#echo $p

#with auto determined ~balanced splits keeping genes together w/ individual headers on separate filesystems, assumes 1-base input; used for both human/mouse
cat <(echo "$1") $1 | perl -ne 'chomp; $s=$_; if($s=~/\.tsv$/) { print "$s\n"; next; } @f=split(/\t/,$s); ($c,$s,$e,$g)=split(/\t/,$s); print "$c:$s-$e\t1\t$g\n";' | perl -ne 'BEGIN { $g="NA"; open(OUT,">dummyfile"); $i=0; $i2=-1; $A=66; } chomp; $ls=$_; if($ls=~/\.tsv$/) { $fn=$ls; $wc=`wc -l $fn`; @fn=split(/\//,$fn); $fn=pop(@fn); chomp($wc); $nl=int($wc/'${3}'); $mnl='${3}'*$nl; next; } ($c,$m,$g)=split(/\t/,$ls); $SPLIT=1 if($i++ % $nl == 0 && $i < $mnl); if($SPLIT && $g ne $pg) { $SPLIT=undef; $i2++; close(OUT); $fs="/data".chr($A++); $nfn="$fs/'${2}'.$i2"; print "python client/bulk_base_intervals.py --bulk-query-file $nfn --endpoint bases --datasrc '${2}' > $nfn.run 2>&1\n"; open(OUT,">$nfn"); print OUT "region\tcontains\tgroup\n";} $pg=$g; if($c !~ /region/) {$c1=$c; $c1=~s/-/:/; $g.=":".$c1;} print OUT "$c\t$m\t$g\n"; END { `rm dummyfile`;}'
