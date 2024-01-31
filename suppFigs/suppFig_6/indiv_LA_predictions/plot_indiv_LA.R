library(ggplot2)
setwd("~/yeast-idr-analysis-main/suppFigs/suppFig_6/indiv_LA_predictions/")

sgn = read.csv("sgn1_151_200_pred.csv")
mrn = read.csv("mrn1_61_110.csv")
sgn_scram = read.csv("sgn_scram5_pred.csv")


sgn_plot = ggplot(sgn, aes(x=X, y=pred))+
  geom_point(color = '#377eb8') + 
  scale_x_continuous(expand=c(0.01,0.01), breaks = sgn$X, labels=sgn$aa) + 
  geom_line(color = '#377eb8') + 
  labs(x = "Amino Acid Sequence", y = "Prediction") + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())


scram_plot = ggplot(sgn_scram, aes(x=X, y=pred))+
  geom_point(color = '#377eb8') + 
  scale_x_continuous(expand=c(0.01,0.01), breaks = sgn_scram$X, labels=sgn_scram$aa) + 
  geom_line(color = '#377eb8') + 
  labs(x = "Amino Acid Sequence", y = "Prediction") + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank()) 

mrn_plot = ggplot(mrn, aes(x=X, y=pred))+
  geom_point(color = '#377eb8') + 
  scale_x_continuous(expand=c(0.01,0.01), breaks = mrn$X, labels=mrn$aa) + 
  geom_line(color = '#377eb8') + 
  labs(x = "Amino Acid Sequence", y = "Prediction") + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

if(!file.exists('sgn_indiv.pdf') & !(file.exists('mrn_indiv.pdf'))){
  ggsave('sgn_indiv.pdf', sgn_plot, height = 6, width =12, units='cm')
  ggsave('mrn_indiv.pdf', mrn_plot, height = 6, width =12, units='cm')
}

if(!file.exists("sgn_scram_indiv.pdf")){
  ggsave('sgn_scram_indiv.pdf', scram_plot, height = 6, width =12, units='cm')
}
  

 
