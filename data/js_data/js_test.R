#!/usr/bin/env Rscript
library(JunctionSeq)
countData<-lapply(list.files(pattern="*.counts.tsv"), read.delim,header=FALSE)
gff.file<-"IMPDH1.gff"
design <- data.frame(condition = factor(c("Pituitary","Brain","Spleen")))
jscs = readJunctionSeqCounts(countdata = countData,samplenames = c("SAMP1","SAMP2","SAMP3"),design = design,flat.gff.file = gff.file)
plotJunctionSeqResultsForGene(geneID="IMPDH1",jscs=jscs,plot.type="rawCounts",displayTranscripts = TRUE, plot.exon.results = TRUE, plot.junction.results = TRUE, plot.novel.junction.results = TRUE, plot.untestable.results = TRUE)



