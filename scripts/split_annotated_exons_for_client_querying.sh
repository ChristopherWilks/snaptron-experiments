#!/bin/bash
#takes a tab-delimited file of annotated exons of format:
#chr	start	end	gene_id	exon_length	strand
#and splits them into no more than 1000 bp chunks
#this enables them to be queries for base coverage
#using the Snaptron client w/o disconnects
#usage: cat recount_hg38_gencodev25_disjoint_exons.tsv | split_annotated_exons_for_client_querying.sh > recount_hg38_gencodev25_disjoint_exons.1k_splits.tsv

cat /dev/stdin | perl -ne 'BEGIN { $MIN_SZ=1000; $N=1000; } chomp; $s=$_; @f=split(/\t/,$s); ($chrm,$start,$end,$g,$sz,$strand)=split(/\t/,$s); if($sz <= $MIN_SZ) { print "$s\n"; next; } $n=$sz/$N; $n1=int($n); $n2=$sz-($N*$n1); $i=$start; $idx=0; while(($i+$N)-1 <= $end) { $idx_=$idx; $idx_="0".$idx if($idx<10); $nend = $i+$N-1; $nsz = ($nend - $i) + 1; print "$chrm\t$i\t$nend\t$g\t$nsz\t$strand\n"; $i=$i+$N; $idx++; } $idx_= $idx; $idx_="0".$idx if($idx<10); $nend = $i+$n2-1; $nsz = ($nend - $i) + 1; print "$chrm\t$i\t$nend\t$g\t$nsz\t$strand\n" if($nsz > 0);' | sort -t'	' -k4,4 -k1,1 -k2,2n -k3,3n
