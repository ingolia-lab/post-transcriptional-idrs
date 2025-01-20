library(ggplot2)
library(reshape2)
#load precision recall values
compo_s_pr = read.csv("~/post-transcriptional-idrs/fig5_ms1/idr_logit_models/composition_model_logit_singles/pr_compo_single.csv")
compo_dipep_pr = read.csv('~/post-transcriptional-idrs/fig5_ms1/idr_logit_models/composition_model_logit_withDipeps/pr_compo_dipep.csv')
motif_pr = read.csv('~/post-transcriptional-idrs/fig5_ms1/idr_logit_models/motif_model_logit/pr_motif.csv')
la_pr = read.csv("~/post-transcriptional-idrs/fig5_ms1/idr_logit_models/LightAttention_model/LightAttention_PRC.csv")

compo_s_pr$group = paste0('compo_single ', round(compo_s_pr$AUPRC[1],4))
compo_dipep_pr$group = paste0('compo_dipep ', round(compo_dipep_pr$AUPRC[1],4))
motif_pr$group =paste0('motif ', round(motif_pr$AUPRC[1],4))
la_pr$group = paste0('la ', round(la_pr$AUPRC[1],4))

colors = c('#e41a1c','#377eb8','#4daf4a','#984ea3')

all_pr = rbind(compo_s_pr,compo_dipep_pr,motif_pr,la_pr)


auc_plot = ggplot(all_pr, aes(x = recall, y = precision, group = group, color = group)) + 
  geom_path(size=1.5) + theme_classic() + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(expand = c(0,0.01),limits = c(0,1)) + 
  scale_x_continuous(expand = c(0,0.01),limits = c(0,1)) + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        panel.border = element_blank())




setwd("~/post-transcriptional-idrs/fig5_ms1/auprc_curves/")
if(!file.exists('auprc_curves.pdf')){
  ggsave('auprc_curves.pdf', auc_plot,width = 18,height =18, units='cm' )
}
  











