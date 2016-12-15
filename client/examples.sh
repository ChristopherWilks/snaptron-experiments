#!/bin/bash -x
python query_snaptron.py --query-file ../data/alk_alt_tss.hg19.snap.tsv --datasrc srav1 --function jir | wc -l
python query_snaptron.py --region "chr2:29446395-30142858" --contains 1 --thresholds "samples_count>=100&annotated=1" | wc -l
python query_snaptron.py --region "chr11:82970135-82997450" --contains 1 --thresholds "samples_count>=100&coverage_sum>=1000" --filters "description:cortex" | wc -l
python query_snaptron.py --region "chr2:29446395-30142858" --contains 1 --function exon --exon-length 57-57 | wc -l
python query_snaptron.py --query-file ../data/small_ssc_test.snap.tsv --function shared --datasrc gtex
