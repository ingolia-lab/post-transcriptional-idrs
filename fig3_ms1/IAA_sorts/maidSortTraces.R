setwd("~/post-transcriptional-idrs/")
load(file = 'gateFlow.R')
library(ggplot2)
library(gridExtra)

setwd("~/post-transcriptional-idrs/fig1_ms1/iJL039_sorts/")
wt = format_flow_data(data.frame(exprs(read.FCS('20200718_Yeast GFP RFP Ratio SortK562 d1 AD 8_H2_lib_d6.fcs',transformation = FALSE))),
                      fsca_low = 75000,fsca_high = 245000, rat = 1.4)
wt = wt[wt$good_ah == T & wt$good_ssc==T, ]
setwd("~/post-transcriptional-idrs/fig3_ms1/IAA_sorts/")

ccr4 = format_flow_data(data.frame(exprs(read.FCS('ccr4maid_flow.fcs',transformation = FALSE))),
                        fsca_low = 75000,fsca_high = 245000, rat = 1.4)
ccr4 = ccr4[ccr4$good_ah ==T & ccr4$good_ssc ==T, ]
pop2 = format_flow_data(data.frame(exprs(read.FCS('pop2maid.fcs',transformation = FALSE))),
                       fsca_low = 75000,fsca_high = 245000, rat = 1.4)
pop2 = pop2[pop2$good_ah ==T & pop2$good_ssc ==T, ]
dcp2 = format_flow_data(data.frame(exprs(read.FCS('dcp2maid_flow.fcs',transformation = FALSE))),
                        fsca_low = 75000,fsca_high = 245000, rat = 1.4)
dcp2 = dcp2[dcp2$good_ah==T & dcp2$good_ssc==T, ]

pal = c('#66c2a5','#fc8d62','#8da0cb')

flowTraces_raw = ggplot() +
  geom_density(data = wt, aes(x=Ratio..FITC.A.PE.Texas.Red.A),color='#000000', size=2) +
  geom_density(data = ccr4, aes(x=Ratio..FITC.A.PE.Texas.Red.A),color=pal[1], size=2) + 
  geom_density(data = pop2, aes(x=Ratio..FITC.A.PE.Texas.Red.A), color=pal[2],size=2) + 
  geom_density(data = dcp2, aes(x=Ratio..FITC.A.PE.Texas.Red.A), color=pal[3],size=2) + 
  scale_y_continuous(expand = c(0,0),limits=c(0,1.6e-05))  + 
  xlab('YFP/RFP Ratio') + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1)) + 
  labs(color = "Depletion")

if(!file.exists('maid_flowDensity.pdf')){
  ggsave('maid_flowDensity.pdf',flowTraces_raw, height=9,width=9,units='cm')
}


#normalize density to peak at 1
norm_density = function(flow){
  temp = density(flow$Ratio..FITC.A.PE.Texas.Red.A)
  temp$y = temp$y/max(temp$y)
  temp = data.frame(x = temp$x, y = temp$y)
  return(temp)
}

wt = norm_density(wt)
ccr4 = norm_density(ccr4)
pop2 = norm_density(pop2)
dcp2 = norm_density(dcp2)

flowTraces_norm = ggplot(data = wt, aes(x=x,y=y)) +
  geom_line(data = wt, aes(x=x,y=y,color='None'), size=2) +
  geom_line(data = ccr4, aes(x=x,y=y,color='Ccr4'), size=2) + 
  geom_line(data = pop2, aes(x=x,y=y,color='Pop2'),size=2) + 
  geom_line(data = dcp2, aes(x=x,y=y,color='Dcp2'),size=2) + 
  scale_y_continuous(limits = c(0,1.1), expand = c(0,0))  + 
  xlab('YFP/RFP Ratio') + 
  theme_classic() + 
  scale_color_manual(values = c('#000000', '#377eb8', '#ff7f00', '#4daf4a'))+ 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1)) + 
  labs(color = "Depletion")



