library(seqinr)
library(ggplot2)

setwd("~/post-transcriptional-idrs/fig4_ms1/")

#must run after doing HMM analysis

submotifs = read.csv("subMotifs_hmm.csv")

#remove submotifs less than 3 AA (455/471 pass)
submotifs = submotifs[nchar(submotifs$seq) > 3, ]

actives = submotifs[submotifs$state == 'active',]
inacts = submotifs[submotifs$state == 'inactive',]

setwd("~/post-transcriptional-idrs/fig4_ms1/motif_enrichment/")
if(!file.exists('active_subregions.fa') & !file.exists('inact_subregions.fa')){
  write.fasta(sequences = as.list(actives$seq), names = as.list(actives$pep),file.out = 'active_subregions.fa')
  write.fasta(sequences = as.list(inacts$seq), names = as.list(inacts$pep),file.out = 'inact_subregions.fa')
}

#########################################################
#use composition profiler to get values/significance####
#########################################################



#use composition profiler to get values/significance, save in folder
if(!file.exists('motif_enrichment_analysis.csv')){
  print('use composition profiler to get motif enrichments')
}else{
  enrich = read.csv('motif_enrichment_analysis.csv')
  aa_order = c('D','E','K','R','H','N','Q','S','T','C','A','V','M','I','L','F','Y','W','P','G')
  enrich$AA = factor(enrich$AA, levels=aa_order)
  enrich$significant = ifelse(enrich$pval <= 0.05,'sig','insig')
  
  
  aa_plot = ggplot(enrich, aes(x = AA, y = mean, fill = significant)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev), width = 0.2,size=0.5,  position = position_dodge(0.9)) +
    labs(x = "Amino Acid", y = "Mean Value") +
    scale_fill_manual(values = c("sig" = "#2171b5", "insig" = "#bdbdbd")) +
    geom_hline(yintercept = 0, col = "black") + 
    theme_classic() + ylab('Enrichment (motifs/inactive regions)') + xlab('') +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black",),
          axis.ticks.length = unit(0.1,'cm'),
          legend.position = "none",
          panel.border = element_blank())
  
  ggsave('motif_inact_enrichment.pdf',plot=aa_plot, height = 8, width = 13.5, units = 'cm')
}
