library(ggplot2)

setwd('~/post-transcriptional-idrs/fig5_ms1/idr_logit_models/composition_model_logit_singles/')
compo_s_coef = read.csv("compo_single_coef.csv")

setwd('~/post-transcriptional-idrs/fig5_ms1/idr_logit_models/motif_model_logit/')
motif_coef = read.csv('motif_coef.csv')

coef_df = data.frame('aas' = motif_coef$aas, 'motif' = motif_coef$coef,'compo' = compo_s_coef$coef,
                     colors =  c(rep('#ca0020',5),rep('#f4a582',2),rep('#f7f7f7',4),rep('#0571b0',7), rep('#92c5de',2)))

coef_plot = ggplot(coef_df, aes(x=motif,y=compo,label=aas,color=colors)) + 
  geom_text(size=4,fontface='bold') + 
  scale_color_identity() +
  theme_minimal() +
  ylab('Composition Coefficients') + xlab('Motif Coefficients') + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())



