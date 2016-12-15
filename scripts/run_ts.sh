#!/bin/bash

cd client/
python query_snaptron.py --query-file ../data/rel_splices.hg38.snap.tsv --function ts --datasrc gtex | tee ../rel_ts_list.tsv | Rscript ../scripts/tissue_specificty_testing.R
