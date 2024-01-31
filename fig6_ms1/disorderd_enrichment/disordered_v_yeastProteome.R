#enrichment analysis for motif+compo vs yeast proteome 
library(stringr)
library(seqinr)

setwd("~/post-transcriptional-idrs/screening_lib_seqs/")
idr_lib = read.csv('wt_library_translation.csv')
idr_lib$yorf = sapply(idr_lib$Name, function(x){strsplit(x,'_')[[1]][1]})
idr_lib$start = sapply(idr_lib$Name, function(x){as.numeric(strsplit(x,'_')[[1]][2])})
idr_lib$stop = sapply(idr_lib$Name, function(x){as.numeric(strsplit(x,'_')[[1]][3])})

#load yeast proteome
setwd("~/post-transcriptional-idrs/disordered_library_generation/")
yeast_proteome = read.csv("sgd_nonDubious_Proteins.csv")
#remove stopcodons
yeast_proteome$Sequence = substr(yeast_proteome$Sequence,1,nchar(yeast_proteome$Sequence)-1)

#get groups of overlapping peptides
get_sections = function(df = idr_lib){
  master = data.frame('yorf' = NA, 'sect_start' = NA, 'sect_stop' = NA)
  for(orf in unique(df$yorf)){
    temp = df[df$yorf == orf,]
    temp = temp[order(temp$start),]
    temp$group = cumsum(c(temp$start[1], diff(temp$start) > 50))
    for(group_num in unique(temp$group)){
      master= rbind(master, c(orf,
                              min(temp[temp$group == group_num, 'start']),
                              max(temp[temp$group == group_num, 'stop'])
      ))
    }
  }
  return(master[complete.cases(master),])
}

#take groups and yeast, proteome_extract disorderd regions
get_disorderd_regions = function(sections, yp = yeast_proteome){
  aas = c()
  for(i in 1:nrow(sections)){
    orf_idx = match(sections$yorf[i], yp$Standard.Name)
    prot_seq = yp[orf_idx,'Sequence']
    aas = c(aas, substr(prot_seq,sections$sect_start[i],sections$sect_stop[i]))
  }
  sections$disordered_seq = aas
  return(sections)
}

#run both functions to get disorderd regions, wrapped in one
main = function(df = idr_lib, ref = yeast_proteome){
  sect = get_sections()
  disordered = get_disorderd_regions(sections = sect,yp = ref)
  disordered$yorf = paste0(disordered$yorf,'_disorder') #add to make names distinguishable
  return(disordered)
}


disorder = main()

#write out fastas
setwd("~/post-transcriptional-idrs/fig6_ms1/disorderd_enrichment/")
if(!file.exists('disordered_aas.fa') & !file.exists('yorfs.fa')){
  write.fasta(sequences = as.list(disorder$disordered_seq),names = as.list(disorder$yorf), file.out = 'disordered_aas.fa')
  write.fasta(sequences = as.list(yeast_proteome$Sequence),names = as.list(yeast_proteome$Standard.Name), file.out = 'yorfs.fa')
}


#####use composition profiler to get enrichment, download results as text and copy to csv#######

#use composition profiler to get values/significance, save in folder
if(!file.exists('disordered_composition_enrichment.csv')){
  print('use composition profiler to get enrichment of disorderd AAs')
}else{
  composition = read.csv('disordered_composition_enrichment.csv')
  aa_order = c('D','E','K','R','H','N','Q','S','T','C','A','V','M','I','L','F','Y','W','P','G')
  composition$AA = factor(composition$AA, levels=aa_order)
  composition$significant = ifelse(composition$pval <= 0.05,'sig','insig')
  
  
  aa_plot = ggplot(composition, aes(x = AA, y = mean, fill = significant)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev), width = 0.2,size=0.5,  position = position_dodge(0.9)) +
    labs(x = "Amino Acid", y = "Mean Value") +
    scale_fill_manual(values = c("sig" = "#2171b5", "insig" = "#bdbdbd")) +
    geom_hline(yintercept = 0, col = "black") + 
    theme_classic() + ylab('Enrichment (Disordered Regions/Yeast Proteome)') + xlab('') +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black",),
          axis.ticks.length = unit(0.1,'cm'),
          legend.position = "none",
          panel.border = element_blank())
  
  if(!file.exists('disordered_enrichment.pdf')){
    ggsave('disordered_enrichment.pdf',plot=aa_plot, height = 8, width = 13.5, units = 'cm')
  }
}

