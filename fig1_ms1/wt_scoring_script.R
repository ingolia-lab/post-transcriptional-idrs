library(stats4)


#load rep1 and rep2 YFP files and iRFP files
setwd("~/post-transcriptional-idrs/raw_counts/WT_sortSeq/WT_sortSeq_Rep1/")
wt_yfp_R1 = list.files(pattern ='.txt')[!grepl('iRFP',list.files(pattern='.txt'))]
wt_irfp_R1 = list.files(pattern ='.txt')[!grepl('YFP',list.files(pattern='.txt'))]

setwd("~/post-transcriptional-idrs/raw_counts/WT_sortSeq/WT_sortSeq_Rep2/")
wt_yfp_R2 = list.files(pattern ='.txt')[!grepl('iRFP',list.files(pattern='.txt'))]
wt_irfp_R2 = list.files(pattern ='.txt')[!grepl('YFP',list.files(pattern='.txt'))]


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

#format column names
setwd("~/post-transcriptional-idrs/raw_counts/WT_sortSeq/WT_sortSeq_Rep1/")
yfp_counts_r1 = format_counts(wt_yfp_R1)
irfp_counts_r1 = format_counts(wt_irfp_R1)
colnames(yfp_counts_r1) = c("Name",'US_1','US_2','FL','FR','L','R')
colnames(irfp_counts_r1) = c('Name','FL','FR','L','R','US_1','US_2')
setwd("~/post-transcriptional-idrs/raw_counts/WT_sortSeq/WT_sortSeq_Rep2/")
yfp_counts_r2 = format_counts(wt_yfp_R2)
irfp_counts_r2 = format_counts(wt_irfp_R2)
colnames(yfp_counts_r2) = c("Name",'US_1','US_2','FL','FR','L','R')
colnames(irfp_counts_r2) = c('Name','FL','FR','L','R','US_1','US_2')


#cutoff peps that have less than 25 reads in unsorted, or less than 5 reads in sorted
yfp_counts_r1 = yfp_counts_r1[yfp_counts_r1$US_1 >= 25 & rowSums(yfp_counts_r1[,4:7]) >=5,] #low reads on US_2, so only filtered on US1
yfp_counts_r2 = yfp_counts_r2[yfp_counts_r2$US_1 >= 25  & yfp_counts_r2$US_2 >= 25 & rowSums(yfp_counts_r2[,4:7]) >=5,]
irfp_counts_r1 = irfp_counts_r1[irfp_counts_r1$US_1 >= 25 & rowSums(irfp_counts_r1[,2:5]) >= 5, ]
irfp_counts_r2 = irfp_counts_r2[irfp_counts_r2$US_1 >= 25 & irfp_counts_r2$US_2 >= 25 & rowSums(irfp_counts_r2[,2:5]) >= 5,]

#prepare pseudocounts for yfp and irfp
pseudo_yfp_r1 <- data.frame(FL = yfp_counts_r1$FL / sum(yfp_counts_r1$FL),
                            L = yfp_counts_r1$L / sum(yfp_counts_r1$L),
                            R = yfp_counts_r1$R / sum(yfp_counts_r1$R),
                            FR = yfp_counts_r1$FR/ sum(yfp_counts_r1$FR))

pseudo_yfp_r2 <- data.frame(FL = yfp_counts_r2$FL / sum(yfp_counts_r2$FL),
                          L = yfp_counts_r2$L / sum(yfp_counts_r2$L),
                          R = yfp_counts_r2$R / sum(yfp_counts_r2$R),
                          FR = yfp_counts_r2$FR/ sum(yfp_counts_r2$FR))

pseudo_irfp_r1 <- data.frame(FL = irfp_counts_r1$FL / sum(irfp_counts_r1$FL),
                            L = irfp_counts_r1$L / sum(irfp_counts_r1$L),
                            R = irfp_counts_r1$R / sum(irfp_counts_r1$R),
                            FR = irfp_counts_r1$FR/ sum(irfp_counts_r1$FR))

pseudo_irfp_r2 <- data.frame(FL = irfp_counts_r2$FL / sum(irfp_counts_r2$FL),
                            L = irfp_counts_r2$L / sum(irfp_counts_r2$L),
                            R = irfp_counts_r2$R / sum(irfp_counts_r2$R),
                            FR = irfp_counts_r2$FR/ sum(irfp_counts_r2$FR))


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
activity_score_r1 <- lapply(seq(1,nrow(pseudo_yfp_r1)), function(i) { mlePeak(pseudo_yfp_r1[i,], binBounds)})
yfp_counts_r1$activity_score_r1 <- sapply(activity_score_r1,coef) 
activity_score_r2 <- lapply(seq(1,nrow(pseudo_yfp_r2)), function(i) { mlePeak(pseudo_yfp_r2[i,], binBounds)})
yfp_counts_r2$activity_score_r2 <- sapply(activity_score_r2,coef)

#do stability scores
stability_score_r1 = lapply(seq(1,nrow(pseudo_irfp_r1)), function(i) { mlePeak(pseudo_irfp_r1[i,], binBounds)})
irfp_counts_r1$stability_score_r1 = sapply(stability_score_r1,coef) 
stability_score_r2 = lapply(seq(1,nrow(pseudo_irfp_r2)), function(i) { mlePeak(pseudo_irfp_r2[i,], binBounds)})
irfp_counts_r2$stability_score_r2 = sapply(stability_score_r2, coef)

#merge files into master dataset that includes raw counts. This is a raw file; cleanup up one below
#change col_names to get reps, remove the unsortered from the iRFP, since they're duplicated in the yfp/irfp datasets
colnames(yfp_counts_r1) = c('Name','US_1_1','US_1_2','yfp_FL_1','yfp_FR_1','yfp_L_1','yfp_R_1','activity_score_r1')
colnames(yfp_counts_r2) = c('Name','US_2_1','US_2_2','yfp_FL_2','yfp_FR_2','yfp_L_2','yfp_R_2','activity_score_r2')
irfp_counts_r1 = irfp_counts_r1[,c(1:5,8)]
colnames(irfp_counts_r1) = c('Name','irfp_FL_1','irfp_FR_1','irfp_L_1','irfp_R_1', 'stability_score_r1')
irfp_counts_r2 = irfp_counts_r2[,c(1:5,8)]
colnames(irfp_counts_r2) = c('Name','irfp_FL_2','irfp_FR_2','irfp_L_2','irfp_R_2', 'stability_score_r2')

#make datafile that has the counts and only the raw scores (easiest to just use raw scores going forward)
all_counts = merge(merge(merge(yfp_counts_r1, yfp_counts_r2, by='Name'), irfp_counts_r1, by='Name'), irfp_counts_r2, by='Name')
scores = all_counts[,c("Name",'activity_score_r1','activity_score_r2','stability_score_r1','stability_score_r2')]
scores$avg_activity = rowMeans(scores[,c("activity_score_r1",'activity_score_r2')])
scores$activity_sd = apply(scores[,c("activity_score_r1",'activity_score_r2')],1,sd)
scores$avg_stability = rowMeans(scores[,c("stability_score_r1",'stability_score_r2')])
scores$stability_sd = apply(scores[,c("activity_score_r1",'activity_score_r2')],1,sd)
setwd("~/post-transcriptional-idrs/processed_scores/")

write.csv(all_counts,"wt_yfp_iRFP_counts.csv")
write.csv(scores, 'wt_yfp_irfp_scores.csv')




