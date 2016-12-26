#!/usr/bin/env Rscript

library(dplyr)

m <- read.table(file("stdin"), sep="\t", header=T)
groups <- levels(unique(m$group))
for(grp in groups) {
	  print(grp)
	  print(m %>% filter(group == grp) %>%
                      group_by(tissue) %>%
                      summarize(present_frac=mean(shared), present=sum(shared == 1), absent=sum(shared == 0), n=n()) %>%
                      arrange(desc(present_frac)))
}
