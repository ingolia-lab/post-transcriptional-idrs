library(ggplot2)

setwd("~/post-transcriptional-idrs/fig6_ms1/idr_logit_models/composition_model_logit_withDipeps/")
dipep_coef = read.csv('compo_dipep_coef.csv', na.strings = "")


#remove singles
dipep_coef = dipep_coef[nchar(dipep_coef$dipep) ==2 ,]

df = data.frame('first_aa' = substr(dipep_coef$dipep,1,1), 
                'second_aa' = substr(dipep_coef$dipep,2,2),
                'coef' = dipep_coef$coef)

df$first_aa = factor(df$first_aa, levels=unique(df$first_aa))
df$second_aa = factor(df$second_aa, levels=unique(df$second_aa))


dipep_coef = ggplot(df, aes(first_aa,second_aa, fill=coef)) + geom_tile()  + 
  scale_fill_gradient2(low = "#d7191c", mid = "white", high = "#2c7bb6", midpoint = 0) +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        panel.border = element_blank())

setwd("~/post-transcriptional-idrs/suppFigs/suppFig_6/dipep_coefs_heatmap/")
if(!file.exists('dipep_coefs.pdf')){
  ggsave('dipep_coefs.pdf', plot = dipep_coef, height = 7, width = 7, units = 'cm')
}



