#!/bin/bash

#included for this analysis: an example of converrting BED formatted *exons* to Snaptron formatted junction finding input file
cat data/novel_exons.raw.hg38.bed | python scripts/bed2snaptron.py > ./data/novel_exons.hg38.snap.tsv

#run the actual analysis using the SRAv2 data source
cd client/
python query_snaptron.py --query-file ../data/novel_exons.hg38.snap.tsv --function shared --datasrc srav2 | perl -ne '$s=$_; if(!$p) { print "group\tvalidation\tshared_sample_counts\n"; $p=1; next} $s=~s/ /\t/; print "$s";'  > ../data/novel_exons.hg38.ssc_results.tsv

cd ../
cat ./data/novel_exons.hg38.ssc_results.tsv | Rscript ./scripts/validation_boxplot.R
