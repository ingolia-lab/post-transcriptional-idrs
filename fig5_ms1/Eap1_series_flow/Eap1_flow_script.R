setwd("~/post-transcriptional-idrs/")
load(file = 'gateFlow.R')
library(flowCore)
library(ggpubr)
library(ggplot2)

setwd("~/post-transcriptional-idrs/fig5_ms1/Eap1_series_flow/")


file_names = list.files(pattern ='fcs')
list_tables <- lapply(file_names, 
                      function(x){
                        #applies gating scheme
                        format_flow_data(data.frame(exprs(read.FCS(x,transformation = FALSE))),
                                         fsca_low = 75000,fsca_high = 245000, rat = 1.4)
                      })

names = substr(file_names,1,nchar(file_names)-4)
gated_tables = lapply(list_tables, function(x){x[x$good_ah == TRUE & x$good_ssc == TRUE &
                                                   x$FITC.A >= 0 & x$PE.Texas.Red.A >= 0, ]})
gated_tables =  lapply(gated_tables, function(x){x$YFPratio = x$FITC.A/x$PE.Texas.Red.A;return(x)})
yfp_ratios = sapply(gated_tables, function(x){mean(x[,'YFPratio'])})

df = data.frame(val = yfp_ratios, name = names, variant = sapply(strsplit(names,'_'),'[[',1))

df$val = df$val/max(df$val)

df$variant = factor(df$variant, levels=c('wt','yxl','m1','m2','300s'))


eap1_flow = ggplot(df, aes(x=variant, y=val, fill='black')) + 
  geom_bar(position='dodge', stat='summary',fun.y='mean', fill='#67a9cf') +
  geom_jitter(aes(x = variant), shape = 16, size = 2, 
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1)) + 
  geom_errorbar(stat = 'summary', position = 'dodge', width = 0.2,size=1,color='#525252') +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.05))  + 
  theme_classic() + ylab('YFP/RFP Ratio') + xlab('') +
  theme(panel.grid = element_blank(),
        axis.text=element_text(size=12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/post-transcriptional-idrs/fig5_ms1/Eap1_series_flow/")
if(!file.exists('eap1_flow_series.pdf')){
  ggsave('eap1_flow_series.pdf', eap1_flow,height = 7, width=10, units='cm')
}
