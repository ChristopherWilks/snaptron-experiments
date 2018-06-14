#!/bin/bash

cd client/
python query_snaptron.py --query-file ../data/rel_splices.hg38.snap.tsv --function ts --datasrc gtex | tee ../rel_ts_list.tsv | ../scripts/adjust_kruskal_wallis_input.sh | tee ../rel_ts_list.tsv.adjusted | Rscript ../scripts/tissue_specificty_testing.R
