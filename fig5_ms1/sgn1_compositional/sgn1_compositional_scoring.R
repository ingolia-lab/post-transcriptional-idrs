library(ggplot2)
library(flowCore)
library(tidyr)

setwd("~/post-transcriptional-idrs/")
load(file = 'gateFlow.R')
setwd("~/post-transcriptional-idrs/fig5_ms1/sgn1_compositional/")

file_names = list.files(pattern ='fcs')
list_tables <- lapply(file_names, 
                      function(x){
                        #applies gating scheme
                        format_flow_data(data.frame(exprs(read.FCS(x,transformation = FALSE))),
                                         fsca_low = 30000,fsca_high = 80000, rat = 1.4)
                      })

#keep only gated values and FITC & RFP values above 0. PE.A on this cytometer
gated_tables = lapply(list_tables, function(x){x[x$good_ah == TRUE & x$good_ssc == TRUE &
                                                   x$FITC.A > 0 & x$PE.A > 0, ]})
gated_tables =  lapply(gated_tables, function(x){x$YFPratio = x$FITC.A/x$PE.A;return(x)})

#add names
ids = sapply(strsplit(file_names,'\\.'),'[[',1)

for(i in 1:length(gated_tables)){
  gated_tables[[i]]$name = ids[i]
}


#format and get mean/sd
master = do.call(rbind, gated_tables)
master_long = gather(master, key='channel', value='value','FITC.A','PE.A','APC.A','YFPratio')
master_long = master_long[,c('name','channel','value')]

#this if for individualpoints. get mean values of distribution, then average teh replicates
df = master_long %>% group_by(name, channel) %>% summarise(mean_val = mean(value,na.rm=TRUE))
df$tether = paste0(sapply(strsplit(df$name,'_'),'[[',1),'_',sapply(strsplit(df$name,'_'),'[[',2))
df$replicate = sapply(strsplit(df$name,'_'),'[[',3)
df = subset(df,channel=='YFPratio')

#normalize to highest empty and take avg/sd
df$mean_val = df$mean_val/(max(df[df$tether=='empty_empty','mean_val']))
df = df %>% group_by(tether, channel) %>% summarize(avg = mean(mean_val,na.rm=TRUE), sd = sd(mean_val, na.rm=TRUE))

#remove empty for comparing to predict score, and assemble df
df = subset(df, tether != 'empty_empty')
preds = read.csv("sgn1_predictiveScores.csv")
df$pred = preds$score[c(2,11,3:10,1,12)]

sgn1_cor = ggplot(df,aes(x=pred,y=avg)) + 
  geom_point(size=5,color='#2171b5') + 
  geom_errorbar(aes(x = pred, ymin = avg - sd, ymax = avg + sd), width = 0.3, color='black',size=0.25) +
  ylab("YFP/RFP ratio") + xlab("Predicted Score") + theme_classic() + 
  geom_smooth(method='lm', formula = y~x, se=FALSE, aes(group=1), color='black', linetype='dashed') + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

if(!file.exists('sgn1_compositional.pdf')){
  ggsave('sgn1_compositional.pdf',plot = sgn1_cor, width = 8, height = 8, units='cm')
}
  






