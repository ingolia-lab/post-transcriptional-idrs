library(ggplot2)
library(viridis)
library(ggrastr)
library(ggpubr)


setwd("~/yeast-idr-analysis-main/processed_scores/")

ccr4 = read.csv("ccr4_scores.csv")
pop2 = read.csv('pop2_scores.csv')
dcp2 = read.csv('dcp2_scores.csv')

ccr4_cor = ggplot(ccr4, aes(x=activity_score_r1,y=activity_score_r2))+
  geom_pointdensity(size = 0.5, show.legend = FALSE) + scale_color_viridis() + 
  stat_cor(method="pearson") + 
  scale_y_continuous( expand=c(0,0), limits = c(-2,2)) + 
  scale_x_continuous(expand=c(0,0), limits = c(-2,2)) + coord_fixed() + 
  theme_classic() + ylab('-Ccr4 Activity Score, Rep2') + xlab('-Ccr4 Activity Score, Rep1') +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        aspect.ratio = 1,
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

pop2_cor = ggplot(pop2, aes(x=activity_score_r1,y=activity_score_r2))+
  geom_pointdensity(size = 0.5, show.legend = FALSE) + scale_color_viridis() + 
  stat_cor(method="pearson") + 
  scale_y_continuous( expand=c(0,0), limits = c(-2,2)) + 
  scale_x_continuous(expand=c(0,0), limits = c(-2,2)) + coord_fixed() + 
  theme_classic() + ylab('-Pop2 Activity Score, Rep2') + xlab('-Pop2 Activity Score, Rep1') +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        aspect.ratio = 1,
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

dcp2_cor = ggplot(dcp2, aes(x=activity_score_r1,y=activity_score_r2))+
  geom_pointdensity(size = 0.5, show.legend = FALSE) + scale_color_viridis() + 
  stat_cor(method="pearson") + 
  scale_y_continuous( expand=c(0,0), limits = c(-2,2)) + 
  scale_x_continuous(expand=c(0,0), limits = c(-2,2)) + coord_fixed() + 
  theme_classic() + ylab('-Dcp2 Activity Score, Rep2') + xlab('-Dcp2 Activity Score, Rep1') +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        aspect.ratio = 1,
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

ccr4_cor = rasterize(ccr4_cor, layers='Point',dpi=300)
pop2_cor = rasterize(pop2_cor, layers='Point',dpi=300)
dcp2_cor = rasterize(dcp2_cor, layers='Point',dpi=300)

composite_plot <- ccr4_cor + pop2_cor + dcp2_cor
composite_plot <- composite_plot + plot_layout(ncol = 3)

setwd("~/yeast-idr-analysis-main/suppFigs/suppFig_3/")
if(!file.exists('degron_score_cors.pdf')){
  ggsave('degron_score_cors.pdf',composite_plot, height = 8, width = 32, units = 'cm')
}
