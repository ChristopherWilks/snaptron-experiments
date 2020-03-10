#similar to example.sh, but specifically geared to be exactly like the Snapcount vignette for smoke testing

echo "## Basic Queries"
python client/query_snaptron.py --region "CD99" --datasrc gtex --endpoint genes | tail -n+2 | wc -l
python client/query_snaptron.py --region "CD99" --datasrc gtex | tail -n+2 > snapcount_jxs1
wc -l snapcount_jxs1
cut -f 13 snapcount_jxs1 | cut -d',' -f 2- | tr "," \\n | cut -d':' -f 1 | sort -u | wc -l
python client/query_snaptron.py --region "CD99" --metadata "SMTS:Brain" --datasrc gtex | tail -n+2 | wc -l
python client/query_snaptron.py --region "CD99" --metadata "SMTS:Brain" --datasrc gtex --endpoint exons | tail -n+2 | wc -l
python client/query_snaptron.py --region "CD99" --metadata "SMTS:Brain" --datasrc gtex --filters "annotated=1" | tail -n+2 > snapcount_jxs2
wc -l snapcount_jxs2
cut -f 13 snapcount_jxs2 | cut -d',' -f 2- | tr "," \\n | cut -d':' -f 1 | sort -u | wc -l

## PSI
echo "## PSI"
python client/query_snaptron.py --query-file data/test_psi_snapcount.snap.tsv --function psi --datasrc gtex | tail -n+2 > snapcount_psi.test
cut -f 1-6 snapcount_psi.test | sort -k1,1nr -k6,6 | head -6

## JIR
echo "## JIR"
python client/query_snaptron.py --query-file data/test_jir_snapcount.snap.tsv --function jir --datasrc srav2 | tail -n+2 > snapcount_jir.test
cut -f 1-6 snapcount_jir.test | sort -k1,1nr -k6,6 | head -6

## SSC
echo "## SSC"
python client/query_snaptron.py --query-file data/small_ssc_test.snap.tsv --function shared --datasrc gtex

## TS (single)
echo "## TS (single)"
python client/query_snaptron.py --query-file data/test_ts1_snapcount.snap.tsv --function ts --datasrc gtex | sort -k2,2 | head -6

## TS (double)
echo "## TS (double)"
python client/query_snaptron.py --query-file data/test_ts2_snapcount.snap.tsv --function ts --datasrc gtex | sort -k2,2 | head -6

## Union (Merge)
echo "## Union (merge)"
python client/query_snaptron.py --datasrc gtex,tcga --bulk-query-file data/test_union_snapcount.snap.tsv --bulk-query-stdout | tail -n+2 > snapcount_union.test
wc -l snapcount_union.test
cut -f 13 snapcount_union.test | cut -d',' -f 2- | tr "," \\n | cut -d':' -f 1 | sort -u | wc -l

## Intersection  (check)
# can't directly do this via a query, but we can find it easily from the results of the Union/merge query above
echo "## Intersection (check)"
cut -f 18 snapcount_union.test | fgrep "," | wc -l
