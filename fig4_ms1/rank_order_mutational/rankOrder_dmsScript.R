library(ggplot2)
library(readxl)

setwd("~/yeast-idr-analysis-main/processed_scores/")
dms_scores = read.csv('dms_scores.csv')
wt_scores = read.csv("wt_scores_sortSeq.csv")
#get repressor and inactive_control frag list
setwd("~/yeast-idr-analysis-main/screening_lib_seqs/")
active_frags = read_xlsx('dms_library_seqs.xlsx', sheet = 'active_DMS')
inact_frags = read_xlsx('dms_library_seqs.xlsx', sheet = 'inactive_DMS')

active_frags = unique(sub("_[^_]+_[^_]+$", "", active_frags$Peptide))
inact_frags = unique(sub("_[^_]+_[^_]+$", "", inact_frags$Peptide))

#rank order plot
dms_wt_peps = dms_scores[grepl('WT', dms_scores$Name),]
dms_wt_peps$base_frag = sub("_[^_]+_[^_]+$", "", dms_wt_peps$Name)
dms_wt_peps$inact = ifelse(dms_wt_peps$base_frag %in% inact_frags, T, F)
dms_wt_peps$rank = rank(dms_wt_peps$avg_activity)
dms_act = dms_wt_peps[dms_wt_peps$inact == F, ]
dms_inact = dms_wt_peps[dms_wt_peps$inact == T, ]

rank_order_plot = ggplot(dms_act, aes(x=rank,y=avg_activity))+
  geom_point(size=2, color='#f03b20') + 
  geom_point(data = dms_inact, aes(x=rank, y=avg_activity), size=2, color='black') + 
  theme_classic() + xlab('Rank') + ylab('Average Activity') + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

setwd('~/yeast-idr-analysis-main/fig4_ms1/rank_order_mutational/')
if(!file.exists('rankOrder_dms.pdf')){
  ggsave(filename = 'rankOrder_dms.pdf',plot = rank_order_plot, height=4.5, width = 6,units='cm')
}

