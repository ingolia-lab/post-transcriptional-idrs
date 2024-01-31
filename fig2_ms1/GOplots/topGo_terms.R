library(ggplot2)
library(viridis)
setwd("~/post-transcriptional-idrs/processed_scores/")

scores = read.csv("wt_yfp_irfp_scores.csv")
scores$yorf = sapply(scores$Name, function(x){sub("^(.*?)_.*", "\\1", x)})
scores$yorf = substr(scores$yorf, 1, nchar(scores$yorf)-1)

repressors = scores[scores$avg_activity <= -1,]
repressors = unique(repressors$yorf)


#for GO
setwd("~/post-transcriptional-idrs/fig2_ms1/GOplots/")
if(!(file.exists('repressive_orfs.txt')) & !(file.exists('all_orfs.txt'))){
  write.table(repressors, file = 'repressive_orfs.txt',sep='\t', row.names = F, quote = F)
  write.table(unique(scores$yorf), file='all_orfs.txt', sep='\t', row.names=F, quote= F)
}


if(!file.exists('pantherDB_GOterms.csv')){
  print('Run pantherDB to get GO terms')
}else{
  go_terms = read.csv("pantherDB_GOterms.csv")
}

colnames(go_terms) = c('GO_process', 'all_orfs.txt - REFLIST','repressive_orfs_found',
                       'repressive_orfs_expected','repressive_orfs.txt_(over/under)','fold_enrichment','rawPval','FDR')
go_terms = go_terms[7:nrow(go_terms),]



#plot select ones with GO terms > 3. 
top_terms = go_terms[go_terms$fold_enrichment > 3, ]
top_terms$GO_process = gsub(',','\n',top_terms$GO_process)
top_terms$GO_process = gsub("\\(","\n(",top_terms$GO_process)
toShow = top_terms[c(1,3,4,5,11),]
toShow$FDR = as.numeric(toShow$FDR)
toShow$repressive_orfs_found = as.numeric(toShow$repressive_orfs_found)
toShow$GO_process <- factor(toShow$GO_process, 
                           levels = toShow$GO_process[order(toShow$fold_enrichment)])

topTerms = ggplot(toShow, aes(x = as.numeric(fold_enrichment), y =GO_process, color=FDR,size=repressive_orfs_found)) + 
  geom_point() + 
  scale_size_continuous(range = c(13/27*7, 27/27*7)) + 
  labs(x = "Fold Enrichment", y = "GO Term") + 
  xlim(c(3.5,4.5)) + 
  theme_minimal()  + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        aspect.ratio = 1.5,
        axis.ticks.length = unit(0.1,'cm'),
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))


  

if(!file.exists('topGO_terms.pdf')){
  ggsave('topGo_terms.pdf',topTerms,height = 3.5, width = 10, units='in')
}



###generate full table for supplement
go_terms$FDR = as.numeric(go_terms$FDR)
go_terms$repressive_orfs_found = as.numeric(go_terms$repressive_orfs_found)
go_terms$GO_process <- factor(go_terms$GO_process, 
                            levels = go_terms$GO_process[order(go_terms$fold_enrichment)])
go_terms$fold_enrichment = as.numeric(go_terms$fold_enrichment)



full_terms = ggplot(go_terms, aes(x = as.numeric(fold_enrichment), y =GO_process, color=FDR,size=repressive_orfs_found)) + 
  geom_point() + 
  scale_size_continuous(range = c(7/251*89, 251/251*9)) + 
  labs(x = "Fold Enrichment", y = "GO Term") + guides(size = guide_legend(title = "Repress Found")) + 
  theme_minimal()  +  theme(axis.text.y = element_text(size = 12))

setwd("~/post-transcriptional-idrs/suppFigs/suppFig_2/")
if(!file.exists('full_GOtermsV2.pdf')){
  ggsave('full_GOtermsV2.pdf',full_terms,height = 22, width = 13, units='in')
}



