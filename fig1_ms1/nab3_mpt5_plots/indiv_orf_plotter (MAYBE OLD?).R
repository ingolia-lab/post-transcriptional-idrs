library(ggplot2)
library(gridExtra)

setwd("~/post-transcriptional-idrs/processed_scores/")
scores = read.csv('wt_yfp_irfp_scores.csv')
protLen = read.csv("protein_len.csv")

#format Scores data table (add yorf, start, stop)
scores$yorf = sapply(scores$Name, function(x){sub("^(.*?)_.*$", "\\1", x)})
scores$start = as.numeric(sapply(scores$Name, function(x){sub("^.*?_(.*?)_.*$", "\\1", x)}))
scores$stop = as.numeric(sapply(scores$Name, function(x){sub("^.*_(.*)$", "\\1", x)}))
scores$len = as.numeric(sapply(scores$yorf, function(x){protLen[match(x, protLen$Standard.Name),'Length']}))
scores$cols <- ifelse(scores$avg_stability >= -1, "#67a9cf", "#ef8a62")


nab3 = scores[scores$yorf == 'Nab3p',]
nab3 = nab3[order(nab3$start),]
nab3$group = cumsum(c(nab3$start[1], diff(nab3$start) > 50))


nab3_plot <- ggplot(nab3, aes(x = start, y = avg_activity)) +
  geom_segment(aes(x=start,xend = stop, y= avg_activity, yend = avg_activity,color=cols), size=1) +
  geom_line(aes(x = start+24, y = avg_activity, group=group),color='#bdbdbd')+ 
  geom_errorbar(aes(x = start+24, ymin = avg_activity - activity_sd, ymax = avg_activity + activity_sd), width = 10, color='black',size=0.2) +
  theme_classic() + scale_color_manual(values = c("#67a9cf", "#ef8a62"))+
  xlim(0,nab3$len[1])+ ylim(c(-1.7,1.2)) +  ylab('Activity Score') + xlab('Residue') + 
  geom_vline(xintercept = c(323,405,445,557)) + #domain boundries from RRM/ABD
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          panel.border = element_blank())


Mpt5 = scores[scores$yorf == 'Mpt5p',]
Mpt5 = Mpt5[order(Mpt5$start),]
Mpt5$group = cumsum(c(Mpt5$start[1], diff(Mpt5$start) > 50))



Mpt5_plot <- ggplot(Mpt5, aes(x = start, y = avg_activity)) +
  geom_segment(aes(x=start,xend = stop, y= avg_activity, yend = avg_activity,color=cols), size=1) +
  geom_line(aes(x = start+24, y = avg_activity, group=group),color='#bdbdbd')+ 
  geom_errorbar(aes(x = start+24, ymin = avg_activity - activity_sd, ymax = avg_activity + activity_sd), width = 10, color='black',size=0.2) +
  theme_classic() + scale_color_manual(values = c("#67a9cf", "#ef8a62"))+ ylab('Activity Score') + xlab('Residue') + 
  xlim(0,Mpt5$len[1])+ ylim(c(-1.7,1)) + 
  geom_vline(xintercept = c(141,540)) + #boundries for armadillo like helical
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())


setwd("~/post-transcriptional-idrs/fig1_ms1/nab3_mpt5_plots/")
if(!file.exists('nab3_plot.pdf') & !file.exists('mpt5_plots.pdf')){
  ggsave('nab3_plot.pdf', plot=nab3_plot, height = 5, width = 10, units='cm')
  ggsave('mpt5_plot.pdf', plot=Mpt5_plot, height = 5, width = 10, units='cm')
}



###for supplemental 5, Eap1 


eap1 = scores[scores$yorf == 'Eap1p',]
eap1 = eap1[order(eap1$start),]
eap1$group = cumsum(c(eap1$start[1], diff(eap1$start) > 50))



eap1_plot <- ggplot(eap1, aes(x = start, y = avg_activity)) +
  geom_segment(aes(x=start,xend = stop, y= avg_activity, yend = avg_activity,color=cols), size=1) +
  geom_line(aes(x = start+24, y = avg_activity, group=group),color='#bdbdbd')+ 
  geom_errorbar(aes(x = start+24, ymin = avg_activity - activity_sd, ymax = avg_activity + activity_sd), width = 10, color='black',size=0.2) +
  theme_classic() + scale_color_manual(values = c("#67a9cf", "#ef8a62"))+ ylab('Activity Score') + xlab('Residue') + 
  xlim(0,eap1$len[1])+ ylim(c(-1.7,1)) + 
  geom_vline(xintercept = c(109,115,299,418)) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

setwd('~/post-transcriptional-idrs/suppFigs/suppFig_5/')
if(!file.exists('eap1_fl.pdf')){
  ggsave('eap1_fl.pdf',eap1_plot, height = 5, width = 10, units='cm')
}