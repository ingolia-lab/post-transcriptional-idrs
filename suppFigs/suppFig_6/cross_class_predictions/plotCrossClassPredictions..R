library(ggplot2)
setwd("~/post-transcriptional-idrs/suppFigs/suppFig_6/cross_class_predictions/")

compo = read.csv("compo_pred_motif_model.csv")
motif = read.csv("motif_pred_compo_model.csv")

df = data.frame('pred'=c(compo$motif_prob, motif$compo_model), repressor=c(rep('compo',nrow(compo)),
                                                                         rep('motif',nrow(motif))))

df$col = ifelse(df$repressor == 'compo', '#d7191c','#2c7bb6')

cc_plot = ggplot(df, aes(x=pred, fill= repressor, color=repressor)) + 
  geom_histogram(position='identity', alpha=0.5) + 
  scale_color_manual(values=c("#d7191c", "#2c7bb6")) +
  scale_fill_manual(values=c("#d7191c", "#2c7bb6")) +
  scale_y_continuous(expand = c(0,3)) + 
  scale_x_continuous(expand = c(0.01,0.01)) + 
  ylab('Counts') + xlab('Probability') + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

if(!file.exists("cross_class_plot.pdf")){
  ggsave("cross_class_plot.pdf",cc_plot, height = 7, width = 7, units='cm')
}
