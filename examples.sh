#!/bin/bash
python client/query_snaptron.py --query-file data/alk_alt_tss.hg19.snap.tsv --datasrc srav1 --function jir | wc -l
python client/query_snaptron.py --region "chr2:29446395-30142858" --contains 1 --filters "samples_count>=100&annotated=1" | wc -l
#python client/query_snaptron.py --region "chr11:82970135-82997450" --contains 1 --filters "samples_count>=100&coverage_sum>=1000" --metadata "description:cortex" | wc -l
python client/query_snaptron.py --region "chr11:82970135-82997450" --contains 1 --filters "samples_count>=7&coverage_sum>=10" --metadata "description:cortex" | wc -l
python client/query_snaptron.py --region "chr2:29446395-30142858" --contains 1 --function exon --exon-length 57-57 | wc -l
python client/query_snaptron.py --query-file data/small_ssc_test.snap.tsv --function shared --datasrc gtex
python client/query_snaptron.py --metadata "description:cortex" | wc -l
python client/query_snaptron.py --query-file data/sox9_retained.hg38.tcga.snap.tsv --datasrc tcga --function jir | wc -l
python client/query_snaptron.py --query-file data/test_psi_abcd3.snap.tsv --function psi | wc -l

python client/query_snaptron.py --bulk-query-file data/rel_splices.hg38.snap.tsv --bulk-query-stdout | wc -l

python client/query_snaptron.py --func intersection --query-file data/test_intersection.snap.tsv | wc -l
python client/query_snaptron.py --func intersection --query-file data/test_intersection2.snap.tsv | wc -l

python client/query_snaptron.py --region "chr2:29446395-30142858" --contains 1 --filters "samples_count>=100&annotated=1" --normalize jxcov | wc -l

python client/query_snaptron.py --region "chr2:29446395-30142858" --contains 1 --filters "samples_count>=100&annotated=1" --endpoint genes --normalize recount --datasrc supermouse | grep "ENSMUSG00000039844.19"| cut -f 13 | perl -ne 'chomp; @f=split(/,/,$_); shift(@f); print "".join("\t",(map { ($a,$b)=split(/:/); $b;} @f))."\n";' > g2a
diff tests/test_gene_coverages.tsv g2a

#donor/acceptor demonstration for doing splices mates score (mouse gene AGRN)
python client/query_snaptron.py --datasrc supermouse --region "chr4:156177442-156177442" --acceptor - --min-count 1 | cut -f 1,2,3,4,5-12,232
python client/query_snaptron.py --datasrc supermouse --bulk-query-file ./tests/agrn_acceptor.tsv --min-count 1 --bulk-iterator |  cut -f 1,2,3,4,5-12,232
