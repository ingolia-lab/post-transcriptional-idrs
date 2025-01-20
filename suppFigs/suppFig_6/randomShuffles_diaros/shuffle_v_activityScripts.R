library(ggplot2)
library(ggpubr)
library(viridis)
library(ggpointdensity)
library(ggrastr)
library(stringr)
library(cowplot)
#shuffled seqs are Sgn1_151_200, Mrn1p_61_110, Ngr1p_510_559, Rim4p_602_651,Sup35p_41_90

setwd("~/post-transcriptional-idrs/suppFigs/suppFig_6/randomShuffles_diaros/")
df = read.csv('scoredShuffles.csv')

#count occurrences of diaromatics . str_count will return 1 for FWY, where you'd want 2, 
#these functions deal with correctly counting dipeps for this purpose (eg. FW + WY)
target_dipeps = function(resis = c('W','F','Y')){
  dipeps = expand.grid(resis,resis)
  dipeps = paste0(dipeps$Var1, dipeps$Var2)
  return(dipeps) #pattern matching later
}

count_dipeps = function(pep, diAros = target_dipeps()){
  pattern_len = nchar(diAros[1]) #get pattern length (2 for dipep), all grids are same
  dipeps = c()
  for(i in 1:(nchar(pep)-pattern_len+1)){
    dipeps[i] = substr(pep,i,(i+pattern_len-1))
  }
  counts = sum(str_count(dipeps,paste0(diAros,collapse='|'))) 
  return(counts)
}

df$diaro_count = sapply(df$seq,count_dipeps)

diAro_simulated_scatter = ggplot(df, aes(x=score,y=diaro_count)) + geom_pointdensity(size=1) + 
  scale_color_viridis() + 
  facet_wrap(~id,ncol = 5) + 
  xlab('Predicted Score (linear)') + ylab('Number di-aromatics') +  
  stat_cor(method="pearson") + 
  scale_y_continuous(breaks = round(seq(min(df$diaro_count), max(df$diaro_count), by = 3),1)) + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        aspect.ratio = 1,
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

#add histogram above to better show point density (all still predicted repressive if using logit)
#adjust in illustrator to get final plot. 
diAro_simulated_hist = ggplot(df, aes(x=score)) + 
  geom_density() +
  facet_wrap(~id,ncol = 5) + 
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

final_plot <- plot_grid(diAro_simulated_hist, diAro_simulated_scatter,
                        ncol = 1, align = "v", rel_heights = c(0.3, 1))


if(!file.exists('shuffle_v_dipep.pdf')){
  ggsave('shuffle_v_dipep.pdf',final_plot, height = 8, width = 8,units='in')
}


#compare dipep of compostional repressors and their scrambles
dms = read.csv('~/post-transcriptional-idrs/processed_scores/dms_scores.csv')
compo = read.csv('~/post-transcriptional-idrs/HMM_analysis/composition.csv')
dms$parent = paste0(sapply(strsplit(dms$Name,'_'),'[[',1),'_',
                    sapply(strsplit(dms$Name,'_'),'[[',2),'_',
                    sapply(strsplit(dms$Name,'_'),'[[',3))

dms_seqs = read.csv('~/post-transcriptional-idrs/sequence_tables_dna/Table_S4_Mut.csv')

dms = dms[dms$parent %in% compo$pep & grepl('WT|scramble',dms$Name),]
dms = dms[dms$avg_stability >= -1,]
dms$pepSeq = dms_seqs$pepSeq[match(dms$Name, dms_seqs$name)]

dms$dipeps = sapply(dms$pepSeq,count_dipeps)

diAro_measured = ggplot(dms, aes(x=avg_activity,y=dipeps)) + geom_pointdensity(size=2) + 
  scale_color_viridis() + 
  xlab('Activity Score') + ylab('Number di-aromatics') +  
  stat_cor(method="pearson") + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        aspect.ratio = 1,
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

diAro_measured = rasterize(diAro_measured, layers='Point',dpi=600)

if(!file.exists('measured_vs_diaro.pdf')){
  ggsave('measured_vs_diaro.pdf',diAro_measured, height = 2.5,width = 2.5, units = 'in')
}


