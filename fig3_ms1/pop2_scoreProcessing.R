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
setwd("~/yeast-idr-analysis-main/raw_counts/Pop2_sortSeq/Pop2_sortSeq_Rep1/")
pop2_rep1 = format_counts(list.files(pattern ='.txt'))
colnames(pop2_rep1)[2:ncol(pop2_rep1)] = sapply(strsplit(colnames(pop2_rep1)[2:ncol(pop2_rep1)],'_'), 
                                                function(x){x[4]})
setwd("~/yeast-idr-analysis-main/raw_counts/Pop2_sortSeq/Pop2_sortSeq_Rep2/")
pop2_rep2 = list.files(pattern ='.txt')
pop2_rep2 = format_counts(list.files(pattern ='.txt'))
colnames(pop2_rep2)[2:ncol(pop2_rep2)] = sapply(strsplit(colnames(pop2_rep2)[2:ncol(pop2_rep2)],'_'), 
                                                function(x){x[4]})

#cutoff peps that have less than 25 reads in unsorted, or less than 5 reads in sorted
pop2_rep1 = pop2_rep1[pop2_rep1$`US-1` >= 25 & pop2_rep1$`US-2` >= 25 & rowSums(pop2_rep1[,2:5]) >=5,]
pop2_rep2 = pop2_rep2[pop2_rep2$`US-2` >= 25 & rowSums(pop2_rep2[,2:5]) >=5,] #US-1 had low counts, only filtered on US-2

#prepare pseudocounts for each replicate
pseudo_pop2_r1 <- data.frame(FL = pop2_rep1$FL / sum(pop2_rep1$FL),
                            L = pop2_rep1$L / sum(pop2_rep1$L),
                            R = pop2_rep1$R / sum(pop2_rep1$R),
                            FR = pop2_rep1$FR/ sum(pop2_rep1$FR))

pseudo_pop2_r2 <- data.frame(FL = pop2_rep2$FL / sum(pop2_rep2$FL),
                            L = pop2_rep2$L / sum(pop2_rep2$L),
                            R = pop2_rep2$R / sum(pop2_rep2$R),
                            FR = pop2_rep2$FR/ sum(pop2_rep2$FR))

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
activity_score_r1 <- lapply(seq(1,nrow(pseudo_pop2_r1)), function(i) { mlePeak(pseudo_pop2_r1[i,], binBounds)})
pop2_rep1$activity_score_r1 <- sapply(activity_score_r1,coef) 
activity_score_r2 <- lapply(seq(1,nrow(pseudo_pop2_r2)), function(i) { mlePeak(pseudo_pop2_r2[i,], binBounds)})
pop2_rep2$activity_score_r2 <- sapply(activity_score_r2,coef)

all_counts = merge(pop2_rep1, pop2_rep2, by='Name') #.x = rep 1, .y = rep2
scores = all_counts[,c('Name','activity_score_r1','activity_score_r2')]
scores$avg_activity = rowMeans(scores[c('activity_score_r1','activity_score_r2')])
scores$activity_sd = apply(scores[c('activity_score_r1','activity_score_r2')],1,sd)

setwd('~/yeast-idr-analysis-main/processed_scores/')
write.csv(all_counts, "pop2_allCounts.csv",row.names = F)
write.csv(scores, 'pop2_scores.csv',row.names = F)
