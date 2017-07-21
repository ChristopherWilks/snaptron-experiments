#!/bin/bash
python client/query_snaptron.py --query-file data/alk_alt_tss.hg19.snap.tsv --datasrc srav1 --function jir | wc -l
python client/query_snaptron.py --region "chr2:29446395-30142858" --contains 1 --filters "samples_count>=100&annotated=1" | wc -l
#python client/query_snaptron.py --region "chr11:82970135-82997450" --contains 1 --filters "samples_count>=100&coverage_sum>=1000" --metadata "description:cortex" | wc -l
python client/query_snaptron.py --region "chr11:82970135-82997450" --contains 1 --filters "samples_count>=7&coverage_sum>=10" --metadata "description:cortex" | wc -l
python client/query_snaptron.py --region "chr2:29446395-30142858" --contains 1 --function exon --exon-length 57-57 | wc -l
python client/query_snaptron.py --query-file data/small_ssc_test.snap.tsv --function shared --datasrc gtex
python client/query_snaptron.py --metadata "description:cortex" >cortex_samples.tsv
python client/query_snaptron.py --query-file data/sox9_retained.hg38.tcga.snap.tsv --datasrc tcga --function jir > sox9_jir_samples.tsv
python client/query_snaptron.py --query-file data/test_psi_abcd3.snap.tsv --function psi

python client/query_snaptron.py --bulk-query-file data/rel_splices.hg38.snap.tsv

python client/query_snaptron.py --func intersection --query-file data/intersect_test.snap.tsv
