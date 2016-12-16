#!/bin/sh

cd client/
#SRAv1 (gives results consist with paper)
python query_snaptron.py --query-file ../data/alk_alt_tss.hg19.snap.tsv --function jir --datasrc srav1 > ../data/alk_alt_tss.hg19.srav1.jir_results.tsv

#GTEx
python query_snaptron.py --query-file ../data/alk_alt_tss.hg38.snap.tsv --function jir --datasrc gtex > ../data/alk_alt_tss.hg38.gtex.jir_results.tsv

#TCGA
python query_snaptron.py --query-file ../data/alk_alt_tss.hg38.tcga.snap.tsv --function jir --datasrc tcga > ../data/alk_alt_tss.hg38.tcga.jir_results.tsv
