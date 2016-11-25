#!/bin/bash

cat data/novel_exons.raw.hg38.bed | python scripts/bed2snaptron.py > ./data/novel_exons.hg38.snap.tsv

cd client/
python query_snaptron.py --query-file ../data/novel_exons.hg38.snap.tsv --function shared --datasrc srav2 > ../data/novel_exons.hg38.ssc_results.tsv
