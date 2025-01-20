setwd("~/post-transcriptional-idrs/")
load(file = 'gateFlow.R')
setwd('~/post-transcriptional-idrs/fig1_ms1/frag_vs_sortseq_analysis/Joseph T/')
library(flowCore)
library(minpack.lm)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

#note, indiv genes in same order as flow scripts
indiv_genes = read.csv("iJL051_IndividualGenesV2.csv")
indiv_genes = indiv_genes %>% arrange(YFP_mean)
#format names
indiv_genes$Peptide = paste0(sapply(strsplit(indiv_genes$Peptide,'_'),'[[',1),'\n',
                             '(',sapply(strsplit(indiv_genes$Peptide,'_'),'[[',2),'-',
                             sapply(strsplit(indiv_genes$Peptide,'_'),'[[',3),')')

indiv_genes$color = pal = colorRampPalette(c('#f03b20','#ffeda0'))(16)


file_names = list.files(pattern ='fcs')
list_tables <- lapply(file_names, 
                      function(x){
                        #applies gating scheme
                        format_flow_data(data.frame(exprs(read.FCS(x,transformation = FALSE))),
                                         fsca_low = 75000,fsca_high = 245000, rat = 1.4)
                      })

#format dataframe with yfp/rfp ratios, sds, sort_Seq scores and sortseq sd
gated_tables =  lapply(gated_tables, function(x){x$YFPratio = x$FITC.A/x$PE.Texas.Red.A;return(x)})

gated_tables = lapply(list_tables, function(x){x[x$good_ah == TRUE & x$good_ssc == TRUE &
                                                   x$FITC.A >= 0 & x$PE.Texas.Red.A >= 0, ]})

ids = paste0(sapply(strsplit(file_names,'_'),'[[',2),rep('_',length(file_names)), sapply(strsplit(file_names,'_'),'[[',3))
ids = substr(ids,1,(nchar(ids)-4))

for(i in 1:length(gated_tables)){
  gated_tables[[i]]$name = ids[i]
  gated_tables[[i]]$strain = strsplit(ids[i],'_')[[1]][1]
  gated_tables[[i]]$replicate = strsplit(ids[i],'_')[[1]][2]
}

master = do.call(rbind, gated_tables)
master_long = gather(master, key='name', value='value','PE.Texas.Red.A')
master_long$id = paste0(master_long$strain,'_',master_long$replicate)
master_long$pep = indiv_genes[match(master_long$strain, indiv_genes$yJLstrain),'Peptide']


master_long$pep = factor(master_long$pep,levels=unique(indiv_genes$Peptide))


p = ggplot(master_long,aes(x=pep,y=value, group=id,fill=strain)) + 
  geom_violin(position="dodge", alpha=0.5) + scale_fill_manual(values = rep('#bdbdbd',16)) +
  scale_y_log10() + theme_classic() + ylab('mScarlet fluorescence') + 
  xlab('tethering construct') + 
  theme(legend.position = 'none',
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))


setwd("~/post-transcriptional-idrs/suppFigs/suppFig_1/")
if(!file.exists('rfp_distributions.pdf')){
  ggsave('rfp_distributions.pdf',p,height = 5.94,width =12.84,units = 'cm')
}
