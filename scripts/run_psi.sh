#!/bin/sh

cd client/

python query_snaptron.py --query-file ../data/test_psi_abcd3.snap.tsv --function psi > ../test_psi_abcd3.snap.samples.tsv
