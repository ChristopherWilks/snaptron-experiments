#!/bin/bash

cd client/
python query_snaptron.py --query-file ../data/rel_splices.hg38.tsv --function ts --datasrc gtex | Rscript ../scripts/tissue_specificty_testing.R
