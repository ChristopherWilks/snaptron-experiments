




#SRAv2's sample ID range is: 0-50097
#some will always have 0 coverage

library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)

#this is lifted from Leo Collado's snaptron_query(...) method in the recount2 package:
#https://github.com/leekgroup/recount/blob/master/R/snaptron_query.R
## Helper function for some special variables
to_chr_list <- function(x) {
  r <- strsplit(x, ',')
  i <- which(sapply(r, function(y) {
    y[[1]] == '0'
  }))
  if (length(i) > 0)
    r[i] <- NA
  return(CharacterList(r))
}

getSnaptronDataFrame <- function(url, verbose = FALSE)
{
  s1 <- Sys.time()
  query_data <- read_tsv(url)
  e1 <- Sys.time()
  tt1 <- e1 - s1
  
  if (verbose)
    message(paste(
      tt1,
      "seconds to query Snaptron and download",
      dim(query_data)[1],
      "junction rows"
    ))
  
  samples <- str_split(query_data$samples, ",")
  t <-
    tibble(
      snaptron_id = rep(query_data$snaptron_id, sapply(samples, length)),
      sample = unlist(samples)
    )
  query_data %>% select(-samples) %>% left_join(t, by = "snaptron_id") %>%
    filter(sample != "") %>% separate(sample, into = c("sample", "coverage"))
  query_data
}

url <-
  'http://snaptron.cs.jhu.edu/srav2/snaptron?regions=chr6:1-514015&rfilter=samples_count>:100&header=1'
sresults <- getSnaptronDataFrame(url, verbose = TRUE)
# str(sresults['ranges'])
# str(sresults['sample_counts'])
