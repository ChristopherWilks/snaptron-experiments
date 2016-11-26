#!/usr/bin/env Rscript
infile <- file("stdin")
m<-read.delim(infile)
mx<-max(subset(m,validation == 'validated', select=c(shared_sample_counts)))
m$group <- factor(m$validation,labels=c('Did not validate', 'Validated'))

pdf("./shared_sample_counts.pdf")
boxplot(m$shared_sample_counts ~ m$group, ylab="Shared-sample count")
tg <- wilcox.test(m$shared_sample_counts ~ m$group, alternative="less")
tgp <- paste(c("Wilcox pval",round(tg$p.value,digits=6)), collapse=" ")
text(1.45, mx+10, tgp)
dev.off()
