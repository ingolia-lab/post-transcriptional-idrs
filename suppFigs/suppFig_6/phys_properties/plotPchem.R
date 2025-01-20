library(ggplot2)
library(ggpubr)
library(ggpointdensity)
library(viridis)
library(stringr)

setwd("~/post-transcriptional-idrs/suppFigs/suppFig_6/phys_properties/")

omega = read.csv("omega_vals.csv")
omega_plot = ggplot(omega, aes(x=activity,y=omega_aro)) + 
  geom_pointdensity(size=1) + scale_color_viridis() + 
  xlab('Activity Score') + ylab('Omega, aromatic') +  
  stat_cor(method="pearson") + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        aspect.ratio = 1,
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

if(!file.exists("omegas.pdf")){
  ggsave('omegas.pdf', omega_plot, height = 7, width = 7, units = 'cm')
}

#plot correlation with Y/F/W
hydro_df = read.csv("kd_hydro.csv")
hydro_df$aros = rowSums(hydro_df[,c('Y','F','W')])
hydro_df$num_neg = sapply(hydro_df$seq, function(x){str_count(x, 'D') + str_count(x,'E')})


#plot correlations with amino acid types and hydroscore

kd_hydro_plot = ggplot(hydro_df, aes(x=KD_hydro,y=activity_score)) + 
  geom_pointdensity() + scale_color_viridis() + 
  stat_cor(method="pearson") + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        aspect.ratio = 1,
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))


aro_plot = ggplot(hydro_df, aes(x=aros,y=activity_score))+
  geom_pointdensity() + scale_color_viridis() + 
  stat_cor(method="pearson") + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        aspect.ratio = 1,
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

neg_plot = ggplot(hydro_df, aes(x=num_neg,y=activity_score))+
  geom_pointdensity() + scale_color_viridis() + 
  stat_cor(method="pearson") + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        aspect.ratio = 1,
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

  

kd_hydro_plot = rasterize(kd_hydro_plot, layers='Point',dpi=300)
aro_plot = rasterize(aro_plot, layers='Point',dpi=300)
neg_plot = rasterize(neg_plot, layers='Point',dpi=300)

if(!file.exists('kd_hyd_plot.pdf') &
   !file.exists('numAros.pdf') & 
   !file.exists('numNegs.pdf')){
  ggsave('kd_hyd_plot.pdf', kd_hydro_plot, height = 10,width = 10, units='cm')
  ggsave('numAros.pdf', aro_plot, height = 10,width = 10, units='cm')
  ggsave('numNegs.pdf', neg_plot, height = 10,width = 10, units='cm')
}










