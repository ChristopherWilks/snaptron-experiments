#!/bin/bash
#assumes Rscript is in the path

#included for this analysis: an example of converting BED formatted *exons* to Snaptron formatted, junction-finding input file
cat data/novel_exons.raw.hg38.bed | python scripts/bed2snaptron.py > ./data/novel_exons.hg38.snap.tsv

cd client/
#run the actual analysis using the SRAv2 data source
#a slight re-formatting is needed to rename the header to be specific for the validation analysis
python query_snaptron.py --query-file ../data/novel_exons.hg38.snap.tsv --function shared --datasrc srav2 | perl -ne '$s=$_; if(!$p) { print "group\tvalidation\tshared_sample_counts\n"; $p=1; next} $s=~s/ /\t/; print "$s";'  > ../data/novel_exons.hg38.ssc_results_srav2.tsv

#run the actual analysis using the GTEx data source
python query_snaptron.py --query-file ../data/novel_exons.hg38.snap.tsv --function shared --datasrc gtex | perl -ne '$s=$_; if(!$p) { print "group\tvalidation\tshared_sample_counts\n"; $p=1; next} $s=~s/ /\t/; print "$s";'  > ../data/novel_exons.hg38.ssc_results_gtex.tsv

cd ../
cat ./data/novel_exons.hg38.ssc_results_srav2.tsv | Rscript ./scripts/validation_boxplot.R
mv shared_sample_counts.pdf shared_sample_counts_srav2.pdf

cat ./data/novel_exons.hg38.ssc_results_gtex.tsv | Rscript ./scripts/validation_boxplot.R
mv shared_sample_counts.pdf shared_sample_counts_gtex.pdf
