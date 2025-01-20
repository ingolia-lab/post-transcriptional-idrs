library(ggplot2)
library(dplyr)

setwd("~/post-transcriptional-idrs/")
load('gateFlow.R')
setwd('~/post-transcriptional-idrs/suppFigs/suppFig_6/mrn1_compo_flow/')



file_names = list.files(pattern ='fcs')
list_tables <- lapply(file_names, 
                      function(x){
                        #applies gating scheme
                        format_flow_data(data.frame(exprs(read.FCS(x,transformation = FALSE))),
                                         fsca_low = 30000,fsca_high = 80000, rat = 1.4)
                      })



gated_tables = lapply(list_tables, function(x){x[x$good_ah == TRUE & x$good_ssc == TRUE &
                                                   x$FITC.A > 0 & x$APC.A > 0 & x$PE.A >0,]})

ids = paste0(sapply(strsplit(file_names,'_'),'[[',3),'_',sapply(strsplit(file_names,'_'),'[[',4))

for(i in 1:length(gated_tables)){
  gated_tables[[i]]$name = ids[i]
}

#make YFP and iRFP ratio
gated_tables =  lapply(gated_tables, function(x){x$YFPratio = x$FITC.A/x$PE.A;return(x)})
gated_tables =  lapply(gated_tables, function(x){x$iRFPratio = x$APC.A/x$PE.A;return(x)})


#organize as longform and summarize table
master = do.call(rbind, gated_tables)
master_long = gather(master, key='channel', value='value','FITC.A','PE.A','APC.A','YFPratio','iRFPratio')
master_long = master_long[,c('name','channel','value')]

df = master_long %>% group_by(name, channel) %>% summarise(mean_val = mean(value,na.rm=TRUE))
df$tether = sapply(strsplit(df$name,'_'),'[[',1)
df$replicate = sapply(strsplit(df$name,'_'),'[[',2)


#now plot measure vs prediction
preds = read.csv('~/post-transcriptional-idrs/suppFigs/suppFig_6/mrn1_compo_flow/mrn1_variants_NtoC.csv')
yfp_rat = subset(df,channel=='YFPratio')
#normalize to the highest empty
yfp_rat$mean_val = yfp_rat$mean_val/(max(yfp_rat[yfp_rat$tether=='empty','mean_val']))


output = yfp_rat %>% group_by(tether)  %>% summarize(avg = mean(mean_val,na.rm=TRUE), sd = sd(mean_val,na.rm=TRUE))

#add predicted score and remove the empty for plotting
output$pred_score = preds[match(output$tether,preds$id),'pred_score']
output = output[!is.na(output$pred_score),]


mrn1_cor = ggplot(output, aes(pred_score,avg)) + 
  geom_point(size=5,color='#67a9cf') + 
  geom_errorbar(aes(x = pred_score, ymin = avg - sd, ymax = avg + sd), width = 0.3, color='black',size=0.25) +
  ylab("YFP/RFP ratio") + xlab("Predicted Score") + theme_classic() + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank())

if(!file.exists('mrn1_compositional.pdf')){
  ggsave('mrn1_compositional.pdf',plot = mrn1_cor, width = 8, height = 8, units='cm')
}


