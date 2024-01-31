setwd("~/post-transcriptional-idrs/")
load(file = 'gateFlow.R')
library(flowCore)
library(ggpubr)
library(ggplot2)

setwd("~/post-transcriptional-idrs/fig4_ms1/select_mut_flow/") 

file_names = list.files(pattern ='fcs')
list_tables <- lapply(file_names, 
                      function(x){
                        #applies gating scheme
                        format_flow_data(data.frame(exprs(read.FCS(x,transformation = FALSE))),
                                         fsca_low = 75000,fsca_high = 245000, rat = 1.4)
                      })

tether = sub(".*?_([^_]+_[^_]+)_.*", "\\1", file_names)

gated_tables = lapply(list_tables, function(x){x[x$good_ah == TRUE & x$good_ssc == TRUE &
                                                   x$FITC.A >= 0 & x$PE.Texas.Red.A >= 0, ]})
gated_tables =  lapply(gated_tables, function(x){x$YFPratio = x$FITC.A/x$PE.Texas.Red.A;return(x)})
yfp_ratios = sapply(gated_tables, function(x){mean(x[,'YFPratio'])})

df = data.frame(val = yfp_ratios, name = tether)
df$gene = sapply( strsplit(df$name,'_'),'[[',1)
df$cond = sapply( strsplit(df$name,'_'),'[[',2)
#normalize yfp ratio to 1
df$val = df$val/max(df$val)

df$gene <- factor(df$gene, levels = c('Pho92','Caf20', 'Cth1'))
df$cond <- factor(df$cond, levels = c('WT','Mut'))



ratios = ggplot(df, aes(x=gene, y=val, fill=cond)) + 
  geom_bar(position = "dodge", stat = "summary",fun = "mean", color="black", alpha=0.7) +
  geom_jitter(aes(x = gene), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1),size=2) +
  geom_errorbar(stat = 'summary', position =position_dodge(width = 0.9), width = 0.2,color='#525252') +
  theme_classic()+  scale_y_continuous(expand = c(0,0), limits = c(0,1.1)) + 
  theme(text = element_text(size=20), legend.position = "none") + 
  scale_fill_manual(values = rep(c('#67a9cf','#ef8a62'),4)) + 
  ylab('YFP/RFP Ratio') + xlab('') + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/post-transcriptional-idrs/fig4_ms1/")
ggsave("mut_ratios.pdf",ratios, height = 7, width=7, units='cm')