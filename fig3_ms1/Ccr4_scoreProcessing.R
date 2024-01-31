library(stats4)
library(stringr)

#format counts
format_counts = function(files){
  table_ls <- lapply(files, function(x){
    f = read.table(x, header=T,sep = '')
    colnames(f) = c(x,'Name')
    #delete the unassign read (row 1)
    f = f[2:length(f[,1]),]
    return(f)
  })
  counts = Reduce(function(...) merge(..., all = TRUE, by = "Name"), table_ls)
  counts[is.na(counts)] <- 0
  return(counts)
}

#load rep1 and rep2 YFP files and iRFP files
setwd("~/post-transcriptional-idrs/raw_counts/Ccr4_sortSeq/Ccr4_Rep1/")
ccr4_rep1 = format_counts(list.files(pattern ='.txt'))
colnames(ccr4_rep1)[2:ncol(ccr4_rep1)] = sapply(strsplit(colnames(ccr4_rep1)[2:ncol(ccr4_rep1)],'_'), 
                                                function(x){x[5]})
setwd("~/post-transcriptional-idrs/raw_counts/Ccr4_sortSeq/Ccr4_Rep2/")
ccr4_rep2 = format_counts(list.files(pattern ='.txt'))
colnames(ccr4_rep2)[2:ncol(ccr4_rep2)] = sapply(strsplit(colnames(ccr4_rep2)[2:ncol(ccr4_rep2)],'_'), 
                                                function(x){x[5]})

#cutoff peps that have less than 25 reads in unsorted, or less than 5 reads in sorted
ccr4_rep1 = ccr4_rep1[ccr4_rep1$`US-1` >= 25 & ccr4_rep1$`US-2` >= 25 & rowSums(ccr4_rep1[,2:5]) >=5,]
ccr4_rep2 = ccr4_rep2[ccr4_rep2$`US-2` >= 25 & ccr4_rep2$`US-2` >= 25 & rowSums(ccr4_rep2[,2:5]) >=5,] 

#prepare pseudocounts for each replicate
pseudo_ccr4_r1 <- data.frame(FL = ccr4_rep1$FL / sum(ccr4_rep1$FL),
                             L = ccr4_rep1$L / sum(ccr4_rep1$L),
                             R = ccr4_rep1$R / sum(ccr4_rep1$R),
                             FR = ccr4_rep1$FR/ sum(ccr4_rep1$FR))

pseudo_ccr4_r2 <- data.frame(FL = ccr4_rep2$FL / sum(ccr4_rep2$FL),
                             L = ccr4_rep2$L / sum(ccr4_rep2$L),
                             R = ccr4_rep2$R / sum(ccr4_rep2$R),
                             FR = ccr4_rep2$FR/ sum(ccr4_rep2$FR))

# bounds are c(LOWER, UPPER). #bins correspond to flow gates (Left to Right)
loglikeBin <- function(bounds, peak) {
  log(pnorm(bounds[[2]], mean=peak) - pnorm(bounds[[1]], mean=peak))
}

binBounds = list(FL=qnorm(c(0.001, 0.24999)),
                 L=qnorm(c(0.25,0.49999)),
                 R=qnorm(c(0.50,0.749999)),
                 FR=qnorm(c(0.75,0.99999)))

mlePeak <- function(pseudos, bounds) {
  peak0 <- 0
  
  nll <- function(peak) {
    ll <- c(sum(pseudos$FL) * loglikeBin(bounds$FL, peak),
            sum(pseudos$L)  * loglikeBin(bounds$L,  peak),
            sum(pseudos$R)  * loglikeBin(bounds$R,  peak),
            sum(pseudos$FR) * loglikeBin(bounds$FR, peak))
    ll[!is.finite(ll)] <- -36
    res <- -sum(1e8*ll, na.rm=TRUE)
    res
  }
  
  mle(minuslogl = nll,
      start=list(peak=peak0),
      method="L-BFGS-B",
      lower=-3, upper=3,
      control=(REPORT=1))
}

#run yfp peak for REp1 adn rep2
activity_score_r1 <- lapply(seq(1,nrow(pseudo_ccr4_r1)), function(i) { mlePeak(pseudo_ccr4_r1[i,], binBounds)})
ccr4_rep1$activity_score_r1 <- sapply(activity_score_r1,coef) 
activity_score_r2 <- lapply(seq(1,nrow(pseudo_ccr4_r2)), function(i) { mlePeak(pseudo_ccr4_r2[i,], binBounds)})
ccr4_rep2$activity_score_r2 <- sapply(activity_score_r2,coef)

all_counts = merge(ccr4_rep1, ccr4_rep2, by='Name') #.x = rep 1, .y = rep2
scores = all_counts[,c('Name','activity_score_r1','activity_score_r2')]
scores$avg_activity = rowMeans(scores[c('activity_score_r1','activity_score_r2')])
scores$activity_sd = apply(scores[c('activity_score_r1','activity_score_r2')],1,sd)

setwd('~/post-transcriptional-idrs/processed_scores/')
write.csv(all_counts, "ccr4_allCounts.csv",row.names = F)
write.csv(scores, 'ccr4_scores.csv',row.names = F)
