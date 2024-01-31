setwd("~/post-transcriptional-idrs/")
load(file = 'gateFlow.R')
setwd("~/post-transcriptional-idrs/fig1_ms1/ebs1_ded1_data/")

library(flowCore)
library(ggpubr)
library(ggplot2)

file_names = list.files(pattern ='fcs')
list_tables <- lapply(file_names, 
                      function(x){
                        #applies gating scheme
                        format_flow_data(data.frame(exprs(read.FCS(x,transformation = FALSE))),
                                         fsca_low = 75000,fsca_high = 245000, rat = 1.4)
                      })

#QC to look at gating
#lapply(list_tables, function(x){plot(x$FSC.A,x$FSC.H, col=ifelse(x$good_ah,'#8856a7','#bdbdbd'))})
#lapply(list_tables, function(x){plot(x$FSC.A,x$SSC.A, col=ifelse(x$good_ssc,'#8856a7','#bdbdbd'))})

#keep only gated values and FITC & RFP values above 0
gated_tables = lapply(list_tables, function(x){x[x$good_ah == TRUE & x$good_ssc == TRUE &
                                                   x$FITC.A >= 0 & x$PE.Texas.Red.A >= 0, ]})
gated_tables =  lapply(gated_tables, function(x){x$YFPratio = x$FITC.A/x$PE.Texas.Red.A;return(x)})
gated_tables =  lapply(gated_tables, function(x){x$IRFPratio = x$APC.Cy7.A/x$PE.Texas.Red.A;return(x)})


yfp_ratios = sapply(gated_tables, function(x){mean(x[,'YFPratio'])})

#yJL004 is no tether, yJL021 is N-iRFP+GGS
df = data.frame(val = yfp_ratios, name = c(rep('Ded1',3),rep('Ebs1',3),rep('Empty',3),rep('No Tether',3)))

frags = ggplot(df, aes(x=name, y=val, fill='black')) + 
  geom_bar(position='dodge', stat='summary',fun.y='mean', fill='#2171b5') +
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.2,color='#525252') +
  geom_point(aes(x = name), shape = 16, size = 0.5, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,6))  + 
  theme_classic() + ylab('YFP/RFP Ratio') + xlab('') +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/post-transcriptional-idrs/fig1_ms1/")
ggsave('ebs_ded_barplot.pdf',plot=frags, width = 4.95, height = 3.35, units = 'cm')

###below all goes to supp
#plot iRFP/RFP ratio to supp
#yJL004 is no tether, yJL021 is N-iRFP+GGS
irfp_ratios = sapply(gated_tables, function(x){mean(x[,'IRFPratio'])})
df_irfp = data.frame(val = irfp_ratios, name = c(rep('Ded1',3),rep('Ebs1',3),rep('Empty',3),rep('No Tether',3)))

irfp_frags = ggplot(df_irfp, aes(x=name, y=val, fill='black')) + 
  geom_bar(position='dodge', stat='summary',fun.y='mean', fill='#2171b5') +
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.2,color='#525252') +
  geom_point(aes(x = name), shape = 16, size = 0.5, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.045))  + 
  theme_classic() + ylab('iRFP/RFP Ratio') + xlab('') +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/post-transcriptional-idrs/suppFigs/suppFig_1/")
ggsave('iRFP_barplot.pdf',plot=irfp_frags, width = 4.5, height = 4, units = 'cm')


combo = rbind(transform(gated_tables[[1]], group='ded1'), 
              transform(gated_tables[[4]], group='ebs1'))

#individual flow density for supp
indiv_plots = ggplot(combo, aes(x=PE.Texas.Red.A)) + 
  geom_density(color ='#bdbdbd', fill ='#bdbdbd', alpha=0.5) + facet_grid(group ~.) +
  geom_density(aes(x=FITC.A), color='#8856a7', fill ='#8856a7', alpha=0.5) + scale_x_log10() + 
  scale_y_continuous(expand = c(0,0), limits = c(0,2))  + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        strip.text.y = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())


ggsave('indiv_plots.pdf',plot=indiv_plots, width = 4.95, height = 7, units = 'cm')


combo_irfp = rbind(transform(gated_tables[[4]], group='ded1'), 
              transform(gated_tables[[7]], group='None'))

#remove negative iRFP values
combo_irfp = combo_irfp[combo_irfp$APC.Cy7.A > 0 , ]



indiv_irfp = ggplot(combo_irfp, aes(x=PE.Texas.Red.A)) + 
  geom_density(color ='#bdbdbd', fill ='#bdbdbd', alpha=0.5) + facet_grid(group ~.) +
  geom_density(aes(x=APC.Cy7.A), color='#8856a7', fill ='#8856a7', alpha=0.5) + scale_x_log10() + 
  scale_y_continuous(expand = c(0,0), limits = c(0,2.2))  + 
  theme_classic() + 
  labs(x='Flourescence (A.U)') + 
  theme(panel.grid = element_blank(),
        strip.text.y = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

ggsave('irfp_density_plot.pdf',plot=indiv_irfp, width = 4.95, height = 7, units = 'cm')









