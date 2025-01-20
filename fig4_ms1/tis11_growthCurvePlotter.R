library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(minpack.lm)

setwd("~/post-transcriptional-idrs/fig5_ms1/")


gc = read.csv("tis11_ironGrowth.csv")
gc$time..s. = gc$time..s./3600 #convert time to hours
#set minimum point to 0
mins = apply(gc[2:ncol(gc)],2,min)
gc[,2:ncol(gc)] = sweep(gc[,2:ncol(gc)],MARGIN = 2,STATS = mins,FUN = "-")

sc_df = gc[,c(T,grepl('sc',colnames(gc[2:ncol(gc)])))]
ferro_df = gc[,c(T,grepl('ferrozine',colnames(gc[2:ncol(gc)])))]
iron_add_df = gc[,c(T,grepl('iron',colnames(gc[2:ncol(gc)])))]

sc_df = melt(sc_df, id.vars = 'time..s.', variable.name = 'strain', value.name = "value")
ferro_df = melt(ferro_df, id.vars = 'time..s.', variable.name = 'strain', value.name = "value")
iron_add_df = melt(iron_add_df, id.vars = 'time..s.', variable.name = 'strain', value.name = "value")

purp = brewer.pal(3,'Purples')
org = brewer.pal(3,'Oranges')
grns = brewer.pal(3,'Greens')

ferro_plot <- ggplot(ferro_df, aes(x = time..s., y = value, color=strain)) +
  geom_line(size=1.25) + scale_color_manual(values = c(purp,org,grns)) + 
  labs(x = "Time (h)", y = "OD600") + theme_classic() + 
  scale_y_continuous(expand = c(0,0.005), limits = c(0,0.23)) + 
  scale_x_continuous(expand= c(0,0.005),limits=c(0,24)) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

sc_plot <- ggplot(sc_df, aes(x = time..s., y = value, color=strain)) +
  geom_line(size=1.25) + scale_color_manual(values = c(purp,org,grns)) + 
  labs(x = "Time (h)", y = "OD600") + theme_classic() + 
  scale_y_continuous(expand = c(0,0.02), limits = c(0,1.25)) + 
  scale_x_continuous(expand= c(0,0.005),limits=c(0,24)) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())


iron_add <- ggplot(iron_add_df, aes(x = time..s., y = value, color=strain)) +
  geom_line(size=1.25) + scale_color_manual(values = c(purp,org,grns)) + 
  labs(x = "Time (h)", y = "OD600") + theme_classic() + 
  scale_y_continuous(expand = c(0,0.02), limits = c(0,1.25)) + 
  scale_x_continuous(expand= c(0,0.005),limits=c(0,24)) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

#write files to appropriate locations
if(!file.exists('ferrozine_gc.pdf')){
  ggsave('ferrozine_gc.pdf',ferro_plot, height = 4, width=6.5, units='cm')
  
}


if(!file.exists('~/post-transcriptional-idrs/suppFigs/suppFig_5/tis11_gc_sc.pdf')){
  setwd("~/post-transcriptional-idrs/suppFigs/suppFig_5/")
  ggsave('tis11_gc_sc.pdf',sc_plot, height = 4, width=6.5, units='cm')
}

if(!file.exists('~/post-transcriptional-idrs/suppFigs/suppFig_5/tis11_gc_addBack.pdf')){
  setwd("~/post-transcriptional-idrs/suppFigs/suppFig_5/")
  ggsave('tis11_gc_addBack.pdf',iron_add, height = 4, width=6.5, units='cm')
}

