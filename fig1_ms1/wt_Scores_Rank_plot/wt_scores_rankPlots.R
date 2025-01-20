library(ggplot2)
library(ggrastr)

setwd("~/post-transcriptional-idrs/processed_scores/")
scores = read.csv('wt_yfp_irfp_scores.csv')
protLen = read.csv("protein_len.csv")

#format Scores data table (add yorf, start, stop)
scores$yorf = sapply(scores$Name, function(x){sub("^(.*?)_.*$", "\\1", x)})
scores$start = sapply(scores$Name, function(x){sub("^.*?_(.*?)_.*$", "\\1", x)})
scores$stop = sapply(scores$Name, function(x){sub("^.*_(.*)$", "\\1", x)})
scores$len = sapply(scores$yorf, function(x){protLen[match(x, protLen$Standard.Name),'Length']})

#rank order of scores that are stable
df = data.frame(ranked = rank(scores$avg_activity), 
                avg_act=scores$avg_activity,
                id = ifelse(scores$avg_activity <= -1, 'repressor','inact'),
                col= ifelse(scores$avg_activity <= -1, '#f03b20','#000000')) 
  

rank_order_plot = ggplot(df, aes(x=ranked,y=avg_act)) +
  geom_line(aes(group=id,color=col)) + scale_color_manual(values=c('#000000','#f03b20')) +
  ylab('Activity Score') + xlab('Rank') + geom_hline(yintercept = -1) + 
  scale_y_continuous(limits=c(-1.9,1.97),expand=c(0,0)) + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())



  
setwd("~/post-transcriptional-idrs/fig1_ms1/wt_Scores_Rank_plot/")

if(!file.exists('rank_plot.pdf')){
  ggsave('rank_plot.pdf',rank_order_plot,height = 5, width = 10, units='cm')
}












