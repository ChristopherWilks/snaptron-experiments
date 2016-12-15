#!/bin/bash

cd client/
python query_snaptron.py --query-file ../data/rel_splices.hg38.tsv --function ts --datasrc gtex | tee ../ts_list.tsv | Rscript ../scripts/tissue_specificty_testing.R
