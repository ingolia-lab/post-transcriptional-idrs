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
setwd("~/post-transcriptional-idrs/raw_counts/dms_scan_sortSeq/mutscan_rep1/")
dms_yfp_R1 = format_counts(list.files(pattern ='.txt')[!grepl('iRFP',list.files(pattern='.txt'))])
colnames(dms_yfp_R1) = c('Name','US-1','US-2','FL','FR','L','R')
dms_irfp_R1 = format_counts(list.files(pattern ='.txt')[!grepl('YFP',list.files(pattern='.txt'))])
colnames(dms_irfp_R1) = c('Name','FL','FR','L','R','US-1','US-2')

setwd("~/post-transcriptional-idrs/raw_counts/dms_scan_sortSeq/mutscan_rep2/")
dms_yfp_R2 = format_counts(list.files(pattern ='.txt')[!grepl('iRFP',list.files(pattern='.txt'))])
colnames(dms_yfp_R2) = c('Name','US','FL','FR','L','R')
dms_irfp_R2 = format_counts(list.files(pattern ='.txt')[!grepl('YFP',list.files(pattern='.txt'))])
colnames(dms_irfp_R2) = c('Name','FL','FR','L','R','US-1','US-2')

#cutoff peps that have less than 25 reads in unsorted, or less than 5 reads in sorted
dms_yfp_R1 = dms_yfp_R1[dms_yfp_R1$`US-1` >= 25 & dms_yfp_R1$`US-2` >= 25 & rowSums(dms_yfp_R1[,4:7]) >=5,]
dms_yfp_R2 = dms_yfp_R2[dms_yfp_R2$US >= 25 & rowSums(dms_yfp_R2[,3:6]) >=5,] 
dms_irfp_R1 = dms_irfp_R1[dms_irfp_R1$`US-1` >= 25 & dms_irfp_R1$`US-2` >= 25 & rowSums(dms_irfp_R1[,2:5]) >=5,]
dms_irfp_R2 = dms_irfp_R2[dms_irfp_R2$`US-1` >= 25 & dms_irfp_R2$`US-2` >= 25 & rowSums(dms_irfp_R2[,2:5]) >=5,]

#prepare pseudocounts for each replicate
pseudo_yfp_r1 <- data.frame(FL = dms_yfp_R1$FL / sum(dms_yfp_R1$FL),
                             L = dms_yfp_R1$L / sum(dms_yfp_R1$L),
                             R = dms_yfp_R1$R / sum(dms_yfp_R1$R),
                             FR = dms_yfp_R1$FR/ sum(dms_yfp_R1$FR))

pseudo_yfp_r2 <- data.frame(FL = dms_yfp_R2$FL / sum(dms_yfp_R2$FL),
                             L = dms_yfp_R2$L / sum(dms_yfp_R2$L),
                             R = dms_yfp_R2$R / sum(dms_yfp_R2$R),
                             FR = dms_yfp_R2$FR/ sum(dms_yfp_R2$FR))

pseudo_irfp_r1 <- data.frame(FL = dms_irfp_R1$FL / sum(dms_irfp_R1$FL),
                            L = dms_irfp_R1$L / sum(dms_irfp_R1$L),
                            R = dms_irfp_R1$R / sum(dms_irfp_R1$R),
                            FR = dms_irfp_R1$FR/ sum(dms_irfp_R1$FR))

pseudo_irfp_r2 <- data.frame(FL = dms_irfp_R2$FL / sum(dms_irfp_R2$FL),
                            L = dms_irfp_R2$L / sum(dms_irfp_R2$L),
                            R = dms_irfp_R2$R / sum(dms_irfp_R2$R),
                            FR = dms_irfp_R2$FR/ sum(dms_irfp_R2$FR))

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

#run yfp and iRFP for both rep 1 and rep 2
activity_score_r1 <- lapply(seq(1,nrow(dms_yfp_R1)), function(i) { mlePeak(dms_yfp_R1[i,], binBounds)})
dms_yfp_R1$activity_score_r1 <- sapply(activity_score_r1,coef) 

activity_score_r2 <- lapply(seq(1,nrow(dms_yfp_R2)), function(i) { mlePeak(dms_yfp_R2[i,], binBounds)})
dms_yfp_R2$activity_score_r2 <- sapply(activity_score_r2,coef)

stability_score_r1 <- lapply(seq(1,nrow(dms_irfp_R1)), function(i) { mlePeak(dms_irfp_R1[i,], binBounds)})
dms_irfp_R1$stability_score_r1 <- sapply(stability_score_r1,coef) 

stability_score_r2 <- lapply(seq(1,nrow(dms_irfp_R2)), function(i) { mlePeak(dms_irfp_R2[i,], binBounds)})
dms_irfp_R2$stability_score_r2 <- sapply(stability_score_r2,coef) 

#make counts and score data and save in processed reads
colnames(dms_yfp_R1)[2:7] = paste0(colnames(dms_yfp_R1)[2:7],'_YFP_','R1')
colnames(dms_yfp_R2)[2:6] = paste0(colnames(dms_yfp_R2)[2:6],'_YFP_','R2')
colnames(dms_irfp_R1)[2:7] = paste0(colnames(dms_irfp_R1)[2:7],'_irfp_','R1')
colnames(dms_irfp_R2)[2:7] = paste0(colnames(dms_irfp_R2)[2:7],'_irfp_','R2')
all_counts = merge(merge(merge(dms_yfp_R1, dms_yfp_R2, by='Name'), dms_irfp_R1, by='Name'), dms_irfp_R2, by='Name')
scores = all_counts[,c('Name','activity_score_r1','activity_score_r2','stability_score_r1','stability_score_r2')]
scores$avg_activity = rowMeans(scores[,c("activity_score_r1",'activity_score_r2')])
scores$activity_sd = apply(scores[,c("activity_score_r1",'activity_score_r2')],1,sd)
scores$avg_stability = rowMeans(scores[,c("stability_score_r1",'stability_score_r2')])
scores$stability_sd = apply(scores[,c("activity_score_r1",'activity_score_r2')],1,sd)

setwd("~/post-transcriptional-idrs/processed_scores/")
write.csv(all_counts,"dms_counts.csv")
write.csv(scores, 'dms_scores.csv')




