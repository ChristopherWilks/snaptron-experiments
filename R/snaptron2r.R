
#SRAv2's sample ID range is: 0-50097
#some will always have 0 coverage

library(RCurl)
library(GenomicRanges)

#this is lifted from Leo Collado's snaptron_query(...) method in the recount2 package:
#https://github.com/leekgroup/recount/blob/master/R/snaptron_query.R
## Helper function for some special variables
to_chr_list <- function(x) {
  r <- strsplit(x, ',')
  i <- which(sapply(r, function(y) { y[[1]] == '0' }))
  if(length(i) > 0) r[i] <- NA
  return(CharacterList(r))
}

getSnaptronDataFrame<-function(url, verbose=FALSE)
{
  s1<-Sys.time()
  res<-read.csv(textConnection(getURI(url)), sep='\t')
  e1<-Sys.time()
  tt1<-e1-s1
  if(verbose) message(paste(tt1,"seconds to query Snaptron and download",dim(res)[1],"junction rows"))

  #this block is from Leo as well
  colnames(res) <- c('type', 'snaptron_id', 'chromosome', 'start', 'end',
                     'length', 'strand', 'annotated', 'left_motif', 'right_motif',
                     'left_annotated', 'right_annotated', 'samples',
                     'samples_count', 'coverage_sum', 'coverage_avg', 'coverage_median',
                     'source_dataset_id')
  result <- GRanges(seqnames = res[, 'chromosome'], 
                    IRanges(as.numeric(res[, 'start']), as.numeric(res[, 'end'])),
                    strand = res[, 'strand'])
  result$type <- as.factor(res[, 'type'])
  result$snaptron_id <- as.integer(res[, 'snaptron_id'])
  result$annotated <- to_chr_list(as.character(res[, 'annotated']))
  result$left_motif <- res[, 'left_motif']
  result$right_motif <- res[, 'right_motif']
  result$left_annotated <- to_chr_list(as.character(res[, 'left_annotated']))
  result$right_annotated <- to_chr_list(as.character(res[, 'right_annotated']))
  result$samples_count <- as.integer(res[, 'samples_count'])
  result$coverage_sum <- as.integer(res[, 'coverage_sum'])
  result$coverage_avg <- as.numeric(res[, 'coverage_avg'])
  result$coverage_median <- as.numeric(res[, 'coverage_median'])
  result$source_dataset_id <- as.integer(res[, 'source_dataset_id'])

  #start of cwilks' code; Leo might not want to own this part :)
  nr<-nrow(res)
  nc<-ncol(res)
  rows<-NULL
  cols<-NULL
  vals<-NULL
  
  s1<-Sys.time()
  i<-1
  while(i<=nr)
  {
    c2<-strsplit(as.character(res$samples[i]),",")
    aa<-unlist(lapply(c2[[1]],function(x) { strsplit(x,':')[[1]][1] }))
    ab<-unlist(lapply(c2[[1]],function(x) { strsplit(x,':')[[1]][2] }))
    
    rows<-c(rows,rep(as.integer(res[i, 'snaptron_id']),length(aa[2:length(aa)])))
    cols<-c(cols,aa[2:length(aa)])
    vals<-c(vals,ab[2:length(ab)])
    i<-i+1
  }
  e1<-Sys.time()
  tt3<-e1-s1
  if(verbose) message(paste(tt3,"seconds to setup sample count matrix"))
  #start at 0 since the sample ID's do
  #convert to a sparse matrix representation in a data frame
  sample_counts<-data.frame(rows,as.numeric(cols),as.numeric(vals),stringsAsFactors=FALSE)
  colnames(sample_counts)<-c("junction","sample","coverage")
  return(list("ranges"=result,"sample_counts"=sample_counts))
}

url<-'http://snaptron.cs.jhu.edu/srav2/snaptron?regions=chr6:1-514015&rfilter=samples_count>:100&header=0'
sresults<-getSnaptronDataFrame(url,verbose=TRUE)
str(sresults['ranges'])
str(sresults['sample_counts'])
