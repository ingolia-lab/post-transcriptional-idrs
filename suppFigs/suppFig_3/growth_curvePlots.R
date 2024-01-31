setwd("~/post-transcriptional-idrs/suppFigs/suppFig_3/")
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(minpack.lm)


growth = read.csv("degronStrain_growthCurves.csv")
growth$time..s. = growth$time..s./3600 #convert to hours
tir = growth[, c(T,grepl("tir", names(growth[2:ncol(growth)])))]
ccr4 = growth[, c(T,grepl("ccr4", names(growth[2:ncol(growth)])))]
pop2 = growth[, c(T,grepl("pop2", names(growth[2:ncol(growth)])))]
dcp2 = growth[, c(T,grepl("dcp2", names(growth[2:ncol(growth)])))]

#separate datasets
tir <- melt(tir, id.vars = "time..s.", variable.name = "condition", value.name = "value")
ccr4 <- melt(ccr4, id.vars = "time..s.", variable.name = "condition", value.name = "value")
pop2 <- melt(pop2, id.vars = "time..s.", variable.name = "condition", value.name = "value")
dcp2 <- melt(dcp2, id.vars = "time..s.", variable.name = "condition", value.name = "value")

blues = brewer.pal(3,'Blues')
reds = brewer.pal(3,'Reds')

#plot representative trace for each one

tir_plot <- ggplot(tir, aes(x = time..s., y = value, color=condition)) +
  geom_path() + scale_color_manual(values = c(blues,reds)) + 
  labs(x = "Time (h)", y = "OD600") + theme_classic() + 
  scale_x_continuous(expand=c(0,0)) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

ccr4_plot <- ggplot(ccr4, aes(x = time..s., y = value, color=condition)) +
  geom_path() + scale_color_manual(values = c(blues,reds)) + 
  labs(x = "Time (h)", y = "OD600") + theme_classic() + 
  scale_x_continuous(expand=c(0,0)) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

pop2_plot <- ggplot(pop2, aes(x = time..s., y = value, color=condition)) +
  geom_path() + scale_color_manual(values = c(blues,reds)) + 
  labs(x = "Time (h)", y = "OD600") + theme_classic() + 
  scale_x_continuous(expand=c(0,0)) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

dcp2_plot <- ggplot(dcp2, aes(x = time..s., y = value, color=condition)) +
  geom_path() + scale_color_manual(values = c(blues,reds)) + 
  labs(x = "Time (h)", y = "OD600") + theme_classic() + 
  scale_x_continuous(expand=c(0,0)) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())


combined_plot <- tir_plot + ccr4_plot + pop2_plot + dcp2_plot +
  plot_layout(ncol = 4, widths = c(1, 1, 1, 1))

if(!file.exists('raw_growth.pdf')){
  ggsave('raw_growth.pdf', combined_plot, height = 4, width = 16, units='cm')
}

#fit growth curves, k = max, m = min, r = growth rate
fit_growth <- function(y, x) {
  lin = nlsLM(y ~ k / (1 + ((k - m) / m) * exp(-r * x)),
        start = list(k = 1, r = 0.2, m = 0.1))
  return(coef(lin))
}
fits <- sapply(seq(2, ncol(growth)), function(x) {
  fit_growth(growth[, x], growth[, 1])
})
colnames(fits) = colnames(growth)[2:ncol(growth)]
df_fits = data.frame(rate = fits['r',], strain = sub("_.*", "", colnames(fits)),
                     chem = sub(".*_", "", colnames(fits)))

df_fits$strain <- factor(df_fits$strain, levels = c('tir1','ccr4','pop2','dcp2'))


rates = ggplot(df_fits, aes(x=strain, y=rate, fill=chem)) + 
  geom_bar(position = "dodge", stat = "summary",fun = "mean", color="black", alpha=0.7)+
  geom_jitter(aes(x = strain), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1),size=2)+
  geom_errorbar(stat = 'summary', position =position_dodge(width = 0.9), width = 0.2,color='#525252') +
  theme_classic()+  scale_y_continuous(expand = c(0,0), limits = c(0,0.5))  + 
  theme(text = element_text(size=20), legend.position = "none") + 
  scale_fill_manual(values = rep(c('#67a9cf','#ef8a62'),4)) + 
  ylab('Rate (1/h') + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        aspect.ratio = 1,
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

if(!file.exists('growthRates.pdf')){
  ggsave(filename = 'growthRates.pdf', rates, height = 16, width = 16, units='cm')
}






