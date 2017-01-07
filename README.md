# snaptron-experiments
contains code and scripts to re-create analysis function experiments from the Snaptron paper

Requirements:

Python 2.7

R with ggplot2 


Intermediate results will still be produced even if Rscript
is not found in the path.

## Analyses

1. Shared sample count (SSC) 
[Script](scripts/run_ssc.sh)

   [HG38 Input file](data/novel_exons.raw.hg38.bed)

2. Tissue specificity (TS)
[Script](scripts/run_ts.sh)

   [HG38 Input file](data/rel_splices.hg38.snap.tsv)

3. Junction Inclusion Ratio (JIR)
[Script](scripts/run_jir.sh)

   [HG19 Input file](data/alk_alt_tss.hg19.snap.tsv)

   [GTEx HG38 Input file](data/alk_alt_tss.hg38.snap.tsv)

   [TCGA HG38 Input file](data/alk_alt_tss.hg38.tcga.snap.tsv)


## General Snaptron Client

   [query_snaptron.py](client/query_snaptron.py)

## Examples

Run the [examples script](examples.sh) to see some of the Snaptron client's various options.
