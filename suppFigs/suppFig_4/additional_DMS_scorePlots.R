library(readxl)
library(ggplot2) 
library(viridis)
library(ggpointdensity)
library(ggrastr)
library(flowCore)
setwd("~/post-transcriptional-idrs/")
load(file = 'gateFlow.R')
setwd('~/post-transcriptional-idrs/suppFigs/suppFig_4/')
dms_flow = format_flow_data(data.frame(exprs(read.FCS('dms_yfpRFP_sort.fcs',transformation = FALSE))),
                      fsca_low = 75000,fsca_high = 245000, rat = 1.4)
dms_flow = dms_flow[dms_flow$good_ah == T & dms_flow$good_ssc==T, ]
div = quantile(dms_flow$Ratio..FITC.A.PE.Texas.Red.A, probs = c(0.25, 0.5, 0.75)) #get quartiles

#generate flow density
dms_density = ggplot(data = dms_flow, aes(x=Ratio..FITC.A.PE.Texas.Red.A)) +
  geom_density() +
  geom_vline(xintercept =  div, linetype='dashed') + 
  scale_y_continuous(expand = c(0,0))  + 
  xlab('YFP/RFP Ratio') + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.border = element_blank())

#save and set color gradient in illustrator
ggsave('dms_yfpRFP.pdf', dms_density, width = 4.8, height = 4.8,units='cm')

setwd("~/post-transcriptional-idrs/processed_scores/")
dms_scores = read.csv('dms_scores.csv')
wt_scores = read.csv("wt_scores_sortSeq.csv")
#get repressor and inactive_control frag list
setwd("~/post-transcriptional-idrs/screening_lib_seqs/")
active_frags = read_xlsx('dms_library_seqs.xlsx', sheet = 'active_DMS')
inact_frags = read_xlsx('dms_library_seqs.xlsx', sheet = 'inactive_DMS')

active_frags = unique(sub("_[^_]+_[^_]+$", "", active_frags$Peptide))
inact_frags = unique(sub("_[^_]+_[^_]+$", "", inact_frags$Peptide))

setwd("~/post-transcriptional-idrs/suppFigs/suppFig_4/")

#comparison between replicates for activity and stability scores
rep1_v_rep2_act = ggplot(data = dms_scores, mapping = aes(x=activity_score_r1, y = activity_score_r2)) +
  geom_pointdensity(show.legend = FALSE, size=0.5) + scale_color_viridis() + 
  xlab('Activity Score, Rep-1') + ylab('Activity Score, Rep-2') + 
  coord_fixed(ratio = 1) + 
  theme_classic() + theme(panel.grid = element_blank(),
                            axis.line = element_line(color = "black"),
                            axis.ticks = element_line(color = "black",),
                            axis.ticks.length = unit(0.1,'cm'),
                            legend.position = "none",
                            panel.border = element_blank())
rep1_v_rep2_stab = ggplot(data = dms_scores, mapping = aes(x=stability_score_r1, y = stability_score_r2)) +
  geom_pointdensity(show.legend = FALSE, size=0.5) + scale_color_viridis() + 
  xlab('Stability Score, Rep-1') + ylab('Stability Score, Rep-2') + 
  coord_fixed(ratio = 1) + 
  theme_classic() + theme(panel.grid = element_blank(),
                          axis.line = element_line(color = "black"),
                          axis.ticks = element_line(color = "black",),
                          axis.ticks.length = unit(0.1,'cm'),
                          legend.position = "none",
                          panel.border = element_blank())

#print corrs (entered in illustrator manually), rasterize points. 
print(cor(dms_scores$activity_score_r1, dms_scores$activity_score_r2))
print(cor(dms_scores$stability_score_r1, dms_scores$stability_score_r2))


ggsave('dms_activity_cor.pdf',rasterize(rep1_v_rep2_act, layers='Point',dpi=300), height = 6, width = 6, units='cm')
ggsave('dms_stability_cor.pdf',rasterize(rep1_v_rep2_stab, layers='Point',dpi=300), height = 6, width = 6, units='cm')


#plot agreement between activity scores in original WT vs the WT frags activity in the DMS screen
dms_wt_peps = dms_scores[grepl('WT', dms_scores$Name),]
dms_wt_peps$base_frag = sub("_[^_]+_[^_]+$", "", dms_wt_peps$Name)
dms_wt_peps$wt_score = wt_scores[match(dms_wt_peps$base_frag, wt_scores$Peptide),'YFP_mean']
dms_wt_peps$inact = ifelse(dms_wt_peps$base_frag %in% inact_frags, T, F)

dms_v_wt_score = ggplot(dms_wt_peps, aes(x=wt_score,y=avg_activity,color=inact)) +
  geom_point(size=2,alpha=0.6,stroke=NA) + 
  scale_color_manual(values = c("FALSE" = 'black', "TRUE" = "#f03b20")) + 
  theme_classic() +   coord_fixed(ratio = 1) + 
  xlab('Activity Score in WT') + ylab('Activity Score in DMS') + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/post-transcriptional-idrs/suppFigs/suppFig_4/")
ggsave('dms_vs_wt_actScore.pdf',dms_v_wt_score, height = 6, width = 6, units = 'cm')



