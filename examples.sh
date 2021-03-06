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

tests/gtex_recount_normalization_check.91.sh

#test bulk exon/gene coverage derivation script
python client/bulk_base_intervals.py --bulk-query-file tests/gtex.test.gene.exon.coverage.tsv --endpoint bases --datasrc gtex --summary gene_exon
diff tests/gtex.test.genes.coverage.expected tests/gtex.test.gene.exon.coverage.tsv.snapout.genes.tsv
diff tests/gtex.test.exons.coverage.expected tests/gtex.test.gene.exon.coverage.tsv.snapout.exons.tsv

#donor/acceptor demonstration for doing splices mates score (mouse gene AGRN)
python client/query_snaptron.py --datasrc supermouse --region "chr4:156177442-156177442" --acceptor - --min-count 1 | cut -f 1,2,3,4,5-12,190
python client/query_snaptron.py --datasrc gtex,tcga --bulk-query-file data/alk_alt_tss.hg38.snap.tsv --bulk-query-stdout | cut -f1-10 | grep ":" | cut -f 1,3-5,7 | sort -u > both.all.shared.coords.test

#APSI GTEx query
python client/query_snaptron.py --query-file data/gtex.apsi.8.snapin.tsv --function apsi --datasrc gtex --sample-group-file data/gtex.samples.groups.tsv > gtex.apsi.8.snapin.tsv.output
diff tests/test.gtexv1.apsi.8.tsv gtex.apsi.8.snapin.tsv.output
