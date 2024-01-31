setwd("~/post-transcriptional-idrs/")
load(file = 'gateFlow.R')
setwd('~/post-transcriptional-idrs/fig1_ms1/frag_vs_sortseq_analysis/Joseph T/')
library(flowCore)
library(minpack.lm)
library(ggplot2)

#note, indiv genes in same order as flow scripts
indiv_genes = read.csv("iJL051_IndividualGenesV2.csv")

file_names = list.files(pattern ='fcs')
list_tables <- lapply(file_names, 
                      function(x){
                        #applies gating scheme
                        format_flow_data(data.frame(exprs(read.FCS(x,transformation = FALSE))),
                                         fsca_low = 75000,fsca_high = 245000, rat = 1.4)
                      })

#keep only gated values and FITC & RFP values above 0
gated_tables = lapply(list_tables, function(x){x[x$good_ah == TRUE & x$good_ssc == TRUE &
                                                   x$FITC.A >= 0 & x$PE.Texas.Red.A >= 0, ]})
#format dataframe with yfp/rfp ratios, sds, sort_Seq scores and sortseq sd
gated_tables =  lapply(gated_tables, function(x){x$YFPratio = x$FITC.A/x$PE.Texas.Red.A;return(x)})
yfp_ratios = sapply(gated_tables, function(x){mean(x[,'YFPratio'])})
avg_ratios = sapply(seq(1,length(yfp_ratios)-1, by=2), function(x){mean(yfp_ratios[x:(x+1)])})
sds = sapply(seq(1,length(yfp_ratios)-1, by=2), function(x){sd(yfp_ratios[x:(x+1)])})
df = data.frame(names = unique(sub("^[^_]*_([^_]+)_.*$", "\\1", file_names)), ratio = avg_ratios, sd = sds)
df$sort_score = indiv_genes[match(df$names, indiv_genes$yJLstrain),'YFP_mean']
df$sort_sd = indiv_genes[match(df$names, indiv_genes$yJLstrain),'YFP_sd']
df = df[order(df$sort_score),]
colors = colorRampPalette(c('#f03b20','#ffeda0'))(16)
df$colors = colors


flow_vs_seq = ggplot(df, aes(x = sort_score, y = ratio, color=colors)) +
  geom_point(size=2) +
  scale_color_manual(values = df$colors) + 
  geom_errorbar(aes(ymin = ratio - sd, ymax = ratio + sd), width = 0.00, color='black') +
  geom_errorbarh(aes(xmin = sort_score - sort_sd, xmax = sort_score + sort_sd), height = 0.00, color='black') +
  theme_classic() + ylab('YFP/RFP Ratio') + xlab('Activity Score') +
  geom_smooth(method = "nls", formula = y ~ A * exp(-k * x) + B, method.args = list(start = list(A=1,B=1,k=1)), 
              se = FALSE, color = "black", fullrange=T) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/post-transcriptional-idrs/fig1_ms1/")
ggsave('frag_vs_sortSeq.pdf',plot=flow_vs_seq, width = 4.95, height = 3.35, units = 'cm')


