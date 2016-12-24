# snaptron-experiments
contains code and scripts to re-create analysis function experiments from the Snaptron paper


## Analyses

1. Shared sample count (SSC)

[Script](scripts/run_ssc.sh)

[Input file](data/novel_exons.raw.hg38.bed)

2. Tissue specificity (TS)

[Script](scripts/run_ts.sh)

[Input file](data/rel_splices.hg38.snap.tsv)

3. Junction Inclusion Ratio (JIR)

[Script](scripts/run_jir.sh)

[Input file](data/alk_alt_tss.hg19.snap.tsv)

[Input file](data/alk_alt_tss.hg38.snap.tsv)

[Input file](data/alk_alt_tss.hg38.tcga.snap.tsv)


## General Snpatron Client

[query_snaptron.py](client/query_snaptron.py)
