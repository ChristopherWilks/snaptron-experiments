# snaptron-experiments
Contains code and scripts to re-create analyses from the Snaptron paper.
Also contains a general purpose client for querying the Snaptron web services.

Requirements:

* Python 2.7

* R with ggplot2 


Intermediate results will still be produced even if Rscript
is not found in the path.

**Ask questions in the project's**

[![Join the chat at https://gitter.im/snaptron/Lobby](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/snaptron/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## Analyses

You can run all three analyses from the paper via this [script](./run_all.sh):
	
	./run_all.sh > run_all.out 2>&1

There will be some delay (typically a few minutes)
when a Snaptron compilation (e.g. GTEx)
is accessed for the first time as the whole of the
sample metadata will be downloaded and cached locally.

Output from the scripts is dumped in the working directory.

Intermediate data downloaded from the Snaptron web services is stored in the `snaptron_tmp` directory.


The SSC analysis takes several minutes to complete (~15) as it 
has to make more than 200 queries to both the GTEx and SRAv2
compilations.

The other two should complete within a minute.


1. Shared sample count (SSC) 
[Script](scripts/run_ssc.sh)

  Input:

	[HG38 Input file](data/novel_exons.raw.hg38.bed)

  Output:

	* novel_exons.hg38.ssc_results_srav2.tsv

	* novel_exons.hg38.ssc_results_gtex.tsv

	* shared_sample_counts.pdf

	* The following is written to standard error at the end of each compilation's run:

		* Exons with 0 SSC

		* \# of exons with > 0 SSC

		* \# of exons with > 0 SSC which are fully annotated


2. Tissue specificity (TS)
[Script](scripts/run_ts.sh)

  Input:

	[HG38 Input file](data/rel_splices.hg38.snap.tsv)

  Output:
	* P-values writen to standard out

	* rel_ts_list.tsv


3. Junction Inclusion Ratio (JIR)
[Script](scripts/run_jir.sh)

  Input:

	[HG19 Input file](data/alk_alt_tss.hg19.snap.tsv)

	[GTEx HG38 Input file](data/alk_alt_tss.hg38.snap.tsv)

	[TCGA HG38 Input file](data/alk_alt_tss.hg38.tcga.snap.tsv)

  Output:
	* alk_alt_tss.hg19.srav1.jir_results.tsv

	* alk_alt_tss.hg38.gtex.jir_results.tsv

	* alk_alt_tss.hg38.tcga.jir_results.tsv

4. Percent Spliced In (PSI)
[Script](scripts/run_psi.sh)

  Input: 
  	[HG38 Example ABCD3 Cassette Exon](data/test_psi_abcd3.snap.tsv)
	
  Output:
  	* test_psi_abcd3.snap.samples.tsv


## General Snaptron Client

   [query_snaptron.py](client/query_snaptron.py)

## Examples

Run the [examples script](examples.sh) to see some of the Snaptron client's various options.

