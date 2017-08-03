#!/bin/bash

cd client/
python query_snaptron.py --func intersection --query-file ../data/test_intersection2.snap.tsv > ../test_intersection2.snap.junctions.tsv
