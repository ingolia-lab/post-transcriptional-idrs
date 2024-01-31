setwd("~/post-transcriptional-idrs/")
load(file = 'gateFlow.R')
setwd('~/post-transcriptional-idrs/fig1_ms1/iJL039_sorts/')

library(flowCore)
library(ggplot2)
library(viridis)
library(ggpointdensity)
library(ggrastr)

#for YFP plotting
yfp_sort = data.frame(exprs(read.FCS('20200718_Yeast GFP RFP Ratio SortK562 d1 AD 8_H2_lib_d6.fcs', transformation = F)))
yfp_sort = format_flow_data(df=yfp_sort,fsca_low = 75000,fsca_high = 245000)

#plot YFP/RFP ratio of sorted data
yfp_sort = yfp_sort[yfp_sort$good_ah == T & yfp_sort$good_ssc == T, ]
div = quantile(yfp_sort$Ratio..FITC.A.PE.Texas.Red.A, probs = c(0.25, 0.5, 0.75)) #get quartiles
         
rat_density = ggplot(data = yfp_sort, aes(x=Ratio..FITC.A.PE.Texas.Red.A)) +
  geom_density() +
  geom_vline(xintercept =  div, linetype='dashed') + 
  scale_y_continuous(expand = c(0,0))  + 
  xlab('YFP/RFP Ratio') + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

#color gradient in illustartor
setwd("~/post-transcriptional-idrs/fig1_ms1/")
ggsave('gated_yfpRFP.pdf', rat_density, width = 4.8, height = 4.1,units='cm')

#iRFP plots, need to use ggplot
setwd('~/post-transcriptional-idrs/fig1_ms1/iJL039_sorts/')
iRFP_sort = data.frame(exprs(read.FCS('20200718_Yeast APCCy7,2f, RFP Ratio Sort Yeast.fcs', transformation = F)))

iRFP_sort = format_flow_data(df=iRFP_sort,fsca_low = 75000,fsca_high = 245000)
iRFP_sort = iRFP_sort[iRFP_sort$good_ah == T & iRFP_sort$good_ssc == T, ]


irfp_plot = ggplot(data = iRFP_sort, mapping = aes(x=PE.Texas.Red.A , y = APC.Cy7.A)) +
  geom_pointdensity(show.legend = FALSE) + scale_color_viridis() + 
  theme_classic() + theme(axis.ticks.length=unit(.25, "cm")) +
  scale_x_log10(limits=c(1000,200000)) + scale_y_log10() + xlab("RFP") + ylab('iRFP') +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks.length=unit(0.1086, "cm")) 
    

    
setwd("~/post-transcriptional-idrs/fig1_ms1/")
rasterize(irfp_plot, layers='Point', dpi=300)
ggsave('irfp_RFP_gating.pdf',height = 9, width = 9, units = 'cm')








