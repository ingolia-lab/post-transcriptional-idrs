#cascade plots for supp
library(ggplot2)
library(ggrastr)

setwd("~/yeast-idr-analysis-main/processed_scores/")
scores = read.csv('wt_yfp_irfp_scores.csv')
protLen = read.csv("protein_len.csv")

#format Scores data table (add yorf, start, stop)
scores$yorf = sapply(scores$Name, function(x){sub("^(.*?)_.*$", "\\1", x)})
scores$start = sapply(scores$Name, function(x){sub("^.*?_(.*?)_.*$", "\\1", x)})
scores$stop = sapply(scores$Name, function(x){sub("^.*_(.*)$", "\\1", x)})
scores$len = sapply(scores$yorf, function(x){protLen[match(x, protLen$Standard.Name),'Length']})

repressive_act <- aggregate(avg_activity ~ yorf, scores,function(x){min(x)})
repressive_act$rank = rank(repressive_act$avg_activity)
pal = colorRampPalette(c('#b2182b','#ef8a62','#fddbc7','#d1e5f0','#67a9cf','#2166ac'))
repressive_act$col = pal(nrow(repressive_act))
ordered_scores <- merge(scores, repressive_act, by = "yorf")


activities_plot = ggplot(ordered_scores, aes(x=rank, y=avg_activity.x,color=rank)) + 
  geom_point(size=0.2) + scale_color_gradientn(colors = pal(nrow(ordered_scores))) + 
  scale_x_continuous(expand = c(0,0), limits = c(-25,length(unique(scores$yorf))+25))+ 
  theme_classic() + ylab('Activity Score') + xlab('Rank') + geom_hline(yintercept = -1, linetype='dashed') + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())


#now do the same for stability scores
unstable <- aggregate(avg_stability ~ yorf, scores,function(x){min(x)})
unstable$rank = rank(unstable$avg_stability)
pal = colorRampPalette(c('#b2182b','#ef8a62','#fddbc7','#d1e5f0','#67a9cf','#2166ac'))
unstable$col = pal(nrow(unstable))
ordered_stab <- merge(scores, unstable, by = "yorf")

stabilities_plot = ggplot(ordered_stab, aes(x=rank, y=avg_stability.x,color=rank)) + 
  geom_point(size=0.2) + scale_color_gradientn(colors = pal(nrow(ordered_stab))) + 
  scale_x_continuous(expand = c(0,0), limits = c(-25,length(unique(scores$yorf))+25))+ 
  theme_classic() + ylab('Stability Score') + xlab('Rank') + geom_hline(yintercept = -1, linetype='dashed') + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/yeast-idr-analysis-main/suppFigs/suppFig_2/cascade_plots/")
if(!file.exists('stability_cascade.pdf') & !file.exists('activity_cascade.pdf')){
  ggsave('activity_cascade.pdf', rasterize(activities_plot, layers = 'Point',dpi=300),width=12, height=8, units='cm')
  ggsave('stability_cascade.pdf', rasterize(stabilities_plot, layers = 'Point',dpi=300),width=12, height=8, units='cm')
}






setwd("~/yeast-idr-analysis-main/suppFigs/suppFig_2/")
ggsave('orfs_peps_rank.pdf',rasterize(cascade_plot, layers='Point',dpi=300), width = 12, height=8, units='cm')

#######??????
#perform GO enrichment on stable peptides
stable_peps = scores[scores$avg_stability >= -1, ]
#get unique orfs from repressors
repressive_orfs = unique(stable_peps[stable_peps$avg_activity <= -1, 'yorf'])
repressive_orfs = substr(repressive_orfs,1,nchar(repressive_orfs)-1)
inactive_orfs = unique(stable_peps[stable_peps$avg_activity > -1, 'yorf'])
inactive_orfs = substr(inactive_orfs,1,nchar(inactive_orfs)-1)
write.table(repressive_orfs, file = 'repressive_orfs.txt',sep='\t', row.names = F, quote=F,col.names = F)
write.table(inactive_orfs, file = 'inactive_orfs.txt',sep='\t', row.names = F,quote = F,col.names = F)
#input into PanterDB, use statistical over-represntation and fisher's Exact Test and bonferroni hypothesis correction
if(file.exists('pantherGO_terms.txt')){
  print('using pantherGO_terms.txt')
  go_terms = read.table('pantherGO_terms.txt',sep='\t',col.names = paste0("V",seq_len(7)), fill=T)
}else{
  print('no GOterm file found')
}