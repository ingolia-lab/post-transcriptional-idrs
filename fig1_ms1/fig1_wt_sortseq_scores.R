#plotting activity scores from the WT_sort_Seq for paper. Stability and activty scores

setwd("~/post-transcriptional-idrs/processed_scores/")
wt_scores = read.csv("wt_yfp_irfp_scores.csv")

library(ggplot2)
library(viridis)
library(ggpointdensity)
library(ggrastr)

#bplot correlation between activity scores
p = ggplot(data = wt_scores, mapping = aes(x=activity_score_r1, y = activity_score_r2)) +
geom_pointdensity(show.legend = FALSE) + scale_color_viridis() + 
  theme_classic() + theme(axis.ticks.length=unit(.25, "cm")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks.length=unit(0.1086, "cm")) + 
        xlim(-2,2) + ylim(-2,2) + coord_fixed()


rasterize(p, layers='Point', dpi=300)
print(cor(wt_scores$activity_score_r1,wt_scores$activity_score_r2))
setwd("~/post-transcriptional-idrs/fig1_ms1/")
ggsave('wt_sortSeq_rep1_rep2.pdf',height = 8, width = 8)

#plot correlation between stability scores
stability = ggplot(data = wt_scores, mapping = aes(x=stability_score_r1, y = stability_score_r2)) +
  geom_pointdensity(show.legend = FALSE) + scale_color_viridis() + 
  theme_classic() + theme(axis.ticks.length=unit(.25, "cm")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks.length=unit(0.1086, "cm")) + 
  xlim(-2,2) + ylim(-2,2) + coord_fixed()

rasterize(stability, layers='Point', dpi=300)
print(cor(wt_scores$stability_score_r1,wt_scores$stability_score_r2))
setwd("~/post-transcriptional-idrs/suppFigs/suppFig_1/")
ggsave('wt_stability_rep1rep2.pdf',height = 8, width = 8)


#plot activity vs stability
actStab = ggplot(data = wt_scores, mapping = aes(x=avg_activity, y = avg_stability)) +
  geom_pointdensity(show.legend = FALSE) + scale_color_viridis() + 
  theme_classic() + theme(axis.ticks.length=unit(.25, "cm")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks.length=unit(0.1086, "cm")) + 
  xlim(-2,2) + ylim(-2,2) + coord_fixed()

rasterize(actStab, layers='Point', dpi=300)
print(cor(wt_scores$avg_activity,wt_scores$avg_stability))
ggsave('act_stability_comp.pdf',height = 8, width = 8)








