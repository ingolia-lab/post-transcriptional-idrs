setwd('~/yeast-idr-analysis-main/fig1_ms1/iJL039_sorts/')

library(flowCore)
library(ggplot2)
library(viridis)
library(ggpointdensity)
library(ggrastr)

#for YFP plotting
yfp_sort = data.frame(exprs(read.FCS('20200718_Yeast GFP RFP Ratio SortK562 d1 AD 8_H2_lib_d6.fcs', transformation = F)))

format_flow_data <-function(df, fsca_low=75000, fsca_high=245000,rat = 1.4){
  #remove saturated values
  df = df[df[,'FSC.A'] < max(df[,'FSC.A']) & df[,'SSC.A'] < max(df[,'SSC.A']),]
  
  #for fsca/h, want to get everything above a certain ratio and trim the extreme high/lows
  df$good_ah = ifelse(df[,'FSC.A']/df[,'FSC.H'] < rat & df[,'FSC.A']>fsca_low & df[,'FSC.A']<fsca_high,TRUE, FALSE) 
  
  #ssc/fsca gated within 2 sigma of a linear fit and extreme high/lows trimmer
  ssc_lm = lm(df$SSC.A~df$FSC.A)
  df$good_ssc = ifelse(df$SSC.A < ssc_lm$coefficients[2]* df$FSC.A + ssc_lm$coefficients[1]+(2*sigma(ssc_lm)) &
                         df$SSC.A > ssc_lm$coefficients[2]* df$FSC.A + ssc_lm$coefficients[1]-(2*sigma(ssc_lm)) &
                         df$FSC.A > fsca_low & df$FSC.A < fsca_high, TRUE, FALSE)
  return(df)
}

yfp_sort = format_flow_data(df=yfp_sort,fsca_low = 75000,fsca_high = 245000)

fsc_plot = ggplot(yfp_sort, aes(FSC.A,FSC.H, color = good_ah)) +
  rasterize(geom_point(),dpi=300) + scale_color_manual(values = c('#bdbdbd','#8856a7')) + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

ssc_plot = ggplot(yfp_sort, aes(FSC.A,SSC.A, color = good_ssc)) +
  rasterize(geom_point(), dpi=300) + scale_color_manual(values = c('#bdbdbd','#8856a7')) + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

setwd('~/yeast-idr-analysis-main/suppFigs/suppFig_1/')
ggsave(filename = 'fsc_plot.pdf',plot = fsc_plot,height = 8, width = 8)
ggsave(filename = 'ssc_plot.pdf',plot = ssc_plot,height = 8, width = 8)



