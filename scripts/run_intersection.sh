#!/bin/bash

cd client/
python query_snaptron.py --func intersection --query-file ../data/test_intersection.snap.tsv > ../test_intersection.snap.junctions.tsv
