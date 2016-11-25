#!/bin/sh

cd client/
python query_snaptron.py --query-file ../data/alk_alt_tss.hg19.tsv --function jir --datasrc srav1 > ../data/alk_alt_tss.hg19.jir_results.tsv
