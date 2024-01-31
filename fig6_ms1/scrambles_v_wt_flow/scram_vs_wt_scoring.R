library(ggplot2)
library(flowCore)

setwd("~/post-transcriptional-idrs/")
load(file = 'gateFlow.R')
setwd("~/post-transcriptional-idrs/fig6_ms1/scrambles_v_wt_flow/")

file_names = list.files(pattern ='fcs')
list_tables <- lapply(file_names, 
                      function(x){
                        #applies gating scheme
                        format_flow_data(data.frame(exprs(read.FCS(x,transformation = FALSE))),
                                         fsca_low = 75000,fsca_high = 245000, rat = 1.4)
                      })


#keep only gated values and FITC & RFP values above 0
gated_tables = lapply(list_tables, function(x){x[x$good_ah == TRUE & x$good_ssc == TRUE &
                                                   x$FITC.A > 0 & x$PE.Texas.Red.A > 0, ]})
gated_tables =  lapply(gated_tables, function(x){x$YFPratio = x$FITC.A/x$PE.Texas.Red.A;return(x)})


yfp_ratios = sapply(gated_tables, function(x){mean(x[,'YFPratio'])})
yfp_means = sapply(seq(1,length(yfp_ratios)-2,by=3), function(x){mean(yfp_ratios[x:(x+2)])})
yfp_sd = sapply(seq(1,length(yfp_ratios)-2,by=3), function(x){sd(yfp_ratios[x:(x+2)])})





#load activity scores from dms data
setwd("~/post-transcriptional-idrs/processed_scores/")
dms_scores = read.csv("dms_scores.csv")

#scrambles for the mutants
names = c('Mrn1p_61_110_scramble_1',
          'Mrn1p_61_110_WT_0',
          'Ngr1p_510_559_scramble_3',
          'Ngr1p_510_559_WT_0',
          'Rim4p_602_651_scramble_2',
          'Rim4p_602_651_WT_0',
          'Sgn1p_151_200_scramble_5',
          'Sgn1p_151_200_WT_0',
          'Sup35p_41_90_scramble_4',
          'Sup35p_41_90_WT_0')

dms_sub = dms_scores[dms_scores$Name %in% names,]


df = data.frame('pep' = names, 
                'yfp_mean' = yfp_means, 
                'yfp_sd' = yfp_sd, 
                'dms_act' = dms_sub$avg_activity, 
                'dms_sd' = dms_sub$activity_sd, 
                'type' = rep(c('scram','wt'),5),
                'parent' = c(rep('mrn1',2), rep('ngr1',2), rep('rim4',2), rep('sgn1',2), rep('sup35',2)))

#use black or colors for error bars? 
indiv_scram = ggplot(df, aes(x=yfp_means, y=dms_act, color=parent, shape=type)) + 
  geom_point(size=5) + 
  scale_shape_manual(values=c('wt' = 19, 'scram' = 17)) + 
  scale_color_manual(values = c('mrn1' = '#e41a1c','ngr1'='#377eb8','rim4' = '#4daf4a','sgn1' = '#984ea3','sup35' = '#ff7f00')) + 
  scale_fill_manual(values = c('mrn1' = '#e41a1c','ngr1'='#377eb8','rim4' = '#4daf4a','sgn1' = '#984ea3','sup35' = '#ff7f00')) + 
  geom_errorbar(aes(x = yfp_mean, ymin = dms_act - dms_sd, ymax = dms_act + dms_sd), width = 0.05, color='black',size=0.25) +
  geom_errorbarh(aes(y = dms_act, xmin = yfp_mean - yfp_sd, xmax = yfp_mean + yfp_sd), height = 0.05, color='black', size=0.25) +
  geom_smooth(method='lm', formula = y~x, se=FALSE, aes(group=1), color='black', linetype='dashed') + 
  xlab("YFP/RFP ratio") + ylab("Activity Score") + theme_classic() + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/post-transcriptional-idrs/fig6_ms1/")
if(!file.exists('scram_vs_wt.pdf')){
  ggsave('scram_vs_wt.pdf',plot = indiv_scram, width = 8, height = 8, units='cm')
}

