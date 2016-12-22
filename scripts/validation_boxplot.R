#!/usr/bin/env Rscript

library('ggplot2')

# Read tables
args = commandArgs(trailingOnly=TRUE)
m_s <- read.table(args[1], sep='\t', quote='', comment.char='', header=T)
m_g <- read.table(args[2], sep='\t', quote='', comment.char='', header=T)
m_s$dataset <- 'SRAv2'
m_g$dataset <- 'GTEx'
m <- rbind(m_s, m_g)

# Plot
pdf(file='shared_sample_counts.pdf', width=4, height=4)
ggplot(m, aes(factor(dataset), shared_sample_counts)) +
	  geom_boxplot(aes(fill=factor(validation, labels=c('Failed', 'Passed')))) +
	  theme_bw() +  # less gray
	  labs(x='\nData compilation', y='Shared sample count (SSC)\n') +  # change axis labels
	  scale_fill_discrete(name="Validation") +  # change legend label
	  theme(legend.justification=c(0,1), legend.position=c(0,1))  # put legend in upper left
dev.off()
