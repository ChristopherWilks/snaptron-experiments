#!/bin/sh

cd client/
python query_snaptron.py --query-file ../data/alk_alt_tss.hg19.snap.tsv --function jir --datasrc srav1 > ../data/alk_alt_tss.hg19.jir_results.tsv

python query_snaptron.py --query-file ../data/alk_alt_tss.hg38.snap.tsv --function jir --datasrc tcga > ../data/alk_alt_tss.hg38.jir_results.tsv

