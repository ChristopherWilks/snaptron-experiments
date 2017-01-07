#!/bin/sh

cd client/
#SRAv1 (gives results consistent with paper)
#Nellore et al. further filtered their results by requiring total coverage 
#be >= 50.  This is done here via Perl
python query_snaptron.py --query-file ../data/alk_alt_tss.hg19.snap.tsv --function jir --datasrc srav1 | perl -ne 'chomp; $s=$_; if($s=~/jir_score/) { print "$s\n"; next}; @f=split(/\t/,$_); next if($f[1]+$f[2] < 50); print "".join("\t",@f)."\n";' > ../alk_alt_tss.hg19.srav1.jir_results.tsv

#GTEx
python query_snaptron.py --query-file ../data/alk_alt_tss.hg38.snap.tsv --function jir --datasrc gtex | perl -ne 'chomp; $s=$_; if($s=~/jir_score/) { print "$s\n"; next}; @f=split(/\t/,$_); next if($f[1]+$f[2] < 50); print "".join("\t",@f)."\n";' > ../alk_alt_tss.hg38.gtex.jir_results.tsv

#TCGA
python query_snaptron.py --query-file ../data/alk_alt_tss.hg38.tcga.snap.tsv --function jir --datasrc tcga | perl -ne 'chomp; $s=$_; if($s=~/jir_score/) { print "$s\n"; next}; @f=split(/\t/,$_); next if($f[1]+$f[2] < 50); print "".join("\t",@f)."\n";' > ../alk_alt_tss.hg38.tcga.jir_results.tsv
