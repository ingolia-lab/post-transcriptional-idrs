library(ggplot2)
setwd("~/yeast-idr-analysis-main/processed_scores/")

sgd = read.delim("SGD_features.tab", header=FALSE, quote="",
                 col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                             "parent", "sgdid2", "chrom", "start", "end",
                             "strand", "genpos", "cver", "sver", "desc"))
scores = read.csv("wt_yfp_irfp_scores.csv")
scores$yorf = sapply(scores$Name, function(x){sub("^(.*?)_.*", "\\1", x)})
scores$yorf = toupper(substr(scores$yorf, 1, nchar(scores$yorf)-1))
#get systematic names
scores$sysName = ifelse(scores$yorf %in% sgd$gene,
                        sgd[match(scores$yorf,sgd$gene),'name'],scores$yorf)

setwd("~/yeast-idr-analysis-main/suppFigs/suppFig_2/rbp_agreement/")
rbp = read.csv("yeast_RBP_dataset.csv")
#sum up datasets that ID the protein as an RBP
rbp$occurrences = sapply(seq(1,nrow(rbp)), function(x){sum(grepl('YES',rbp[x,5:12]))})

#get the strongest repressor (lowest score) per orf
get_max_score <-function(df){
  genes = unique(df$yorf)
  master = data.frame(matrix(nrow=length(genes), ncol=ncol(df)))
  for(i in 1:length(genes)){
    temp = df[df$yorf == genes[i],]
    master[i,] = temp[temp$avg_activity == min(temp$avg_activity),]
  }
  colnames(master) = colnames(df)
  return(master)
}

most_repress = get_max_score(scores)
most_repress$occurrences = ifelse(most_repress$sysName %in% rbp$ID,
                                 rbp[match(most_repress$sysName,rbp$ID),'occurrences'],0)

blues = c('#f1eef6','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#034e7b')

rbp_cdf = ggplot(most_repress, aes(x = avg_activity, color = factor(occurrences))) +
  stat_ecdf(geom = "step", size=1.5) +
  labs(x = "Activity Scores", y = "Cumulative Probability") +
  scale_color_manual(values = blues, name = "Occurrences") +
  scale_y_continuous(expand = c(0,0.008))  + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        panel.border = element_blank())

if(!file.exists('rbp_cdf.pdf')){
  ggsave('rbp_cdf.pdf',rbp_cdf,height = 8, width = 14, units = 'cm')
}


