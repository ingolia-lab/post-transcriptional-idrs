library(readxl)
library(depmixS4)
library(ggplot2)
library(qvalue)

setwd("~/post-transcriptional-idrs/processed_scores/")
dms_scores = read.csv("dms_scores.csv")
setwd("~/post-transcriptional-idrs/screening_lib_seqs/")
pep_seqs = read.csv("mutational_library_translation.csv")
repressors = read_xlsx('dms_library_seqs.xlsx', sheet = 'active_DMS')
repressors = unique(sub("_[^_]+_[^_]+$", "", repressors$Peptide))
setwd("~/post-transcriptional-idrs/processed_scores/")
dms_scores = read.csv("dms_scores.csv")

if(!(file.exists('~/post-transcriptional-idrs/processed_scores/hmm_scores.csv'))){
  #format hmmdata
  format_data <-function(df=dms_scores,pep, stability_cutoff = -1){
    foo = df[grepl(paste0(pep,'_Mut'), df$Name),]
    foo$start = sapply(seq(1,nrow(foo)),function(x){as.numeric(strsplit(foo$Name[x],'_')[[1]][5])}) 
    foo$cutoff_activity = ifelse(foo$avg_stability <= stability_cutoff, NA,foo$avg_activity) #replace activity scores of unstable with NA
    #if missing peptides, fill in with Na
    missing = which(!(seq(1,46) %in% foo$start)) 
    empty = data.frame(matrix(ncol=ncol(foo),nrow=length(missing)))
    colnames(empty) = colnames(foo)
    empty$start = missing
    foo = rbind(foo,empty)
    return(foo[order(foo$start),])
  } 
  
  #two peps didn't fit well to 2state model; ignoring
  problems = c("Ded1p_500_549","Nba1p_391_440",'Tbs1p_892_941')
  reps_2state = repressors[!(repressors %in% problems)]
  
  #fit each pep to HMM and extrac tthe lower state as a proxy for activity score; include stdev of both states
  get_hmm_params <- function(df=dms_scores, peps=reps_2state, stability_cutoff = -1,states=2){
    activities = data.frame(name = peps,
                            act_state=rep(0,length(peps)),
                            act_sd=rep(0,length(peps)),
                            inact_state=rep(0,length(peps)),
                            inact_sd=rep(0,length(peps)))
    for(pep in peps){
      print(pep)
      temp = format_data(df=df, pep=pep, stability_cutoff = stability_cutoff)
      hmm = fit(depmix(temp$avg_activity~1, data=temp, nstates = states))
      if(summary(hmm)[1] > summary(hmm)[2]){
        activities[activities$name == pep,] = c(pep, summary(hmm)[2],summary(hmm)[4],summary(hmm)[1],summary(hmm)[3])
      }else{
        activities[activities$name == pep,] = c(pep, summary(hmm)[1],summary(hmm)[3],summary(hmm)[2],summary(hmm)[4])
      }
    }
    return(activities)
  }
  
  #write file, this takes a bit
  hmm_activities = get_hmm_params()
  write.csv(hmm_activities, file='~/post-transcriptional-idrs/processed_scores/hmm_scores.csv',row.names = F)
}else{
  hmm_scores = read.csv("~/post-transcriptional-idrs/processed_scores/hmm_scores.csv", header=T)
}

#use the HMM fit as the mean activity (more accurate). compare to average of scrambles
concensus_scram = function(raw_scores = dms_scores, hmm = hmm_scores, stability_cutoff=-1){
  df = data.frame(name = hmm$name,
                  num_scrams = rep(0,nrow(hmm)),
                  mean_scram_score = rep(0,nrow(hmm)),
                  scram_sd = rep(0,nrow(hmm)))
  for(pep in hmm$name){
    temp = raw_scores[grepl(paste0(pep,'_scramble'), raw_scores$Name) & raw_scores$avg_stability >= stability_cutoff,]
    df[df$name == pep,'num_scrams'] = nrow(temp) 
    df[df$name == pep,'mean_scram_score'] = mean(temp$avg_activity)
    df[df$name == pep,'scram_sd'] = sd(temp$avg_activity)
  }
  return(df)
}

scram_df = concensus_scram()
#require at least 2 scrams
scram_df = scram_df[scram_df$num_scrams >= 2, ]

master = merge(hmm_scores, scram_df, by='name')
master$z_raw = sapply(seq(1,nrow(master)), function(x){
  2*pnorm(-abs(master$mean_scram_score[x] - master$act_state[x])/mean(dms_scores$activity_sd)) #includes sd for unstable peps, 
})
master$z_adj = p.adjust(master$z_raw, 'BH')

#get differences between scramble vs WT score and HMM inact vs HMM act
master$scram_diff = master$mean_scram_score - master$act_state
master$hmm_diff = master$inact_state-master$act_state

num_sig = nrow(master[master$z_adj <= 0.05 & master$scram_diff > 0,])
annot = paste0('Padj < 0.05: ', num_sig,'/',nrow(master))

scram_plot = ggplot(master, aes(x=act_state,y=mean_scram_score))+
  geom_point(aes(color = z_adj < 0.05), size=1) + theme_classic() + 
  xlab('HMM Activity Score') + ylab('Mean Scramble Score') +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#f03b20")) + 
  geom_abline(slope=1, intercept = 0, size=1.5) + 
  annotate(geom = 'text', label = annot, x = -Inf, y = Inf, hjust = 0, vjust = 1) +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))


hmm_v_scram = ggplot(master, aes(x=scram_diff,y=hmm_diff))+
  geom_point(aes(color = z_adj < 0.05), size=1) + theme_classic() + 
  xlab('Scram-WT score') + ylab('inact-act score') +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#f03b20")) + 
  coord_fixed() + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 1, 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

#write figures
setwd("~/post-transcriptional-idrs/fig6_ms1/")
if(!(file.exists('~/post-transcriptional-idrs/fig6_ms1/scram_diff.pdf'))){
  ggsave('scram_diff.pdf', scram_plot,height = 8, width = 8, units  = 'cm')
}
if(!(file.exists('~/post-transcriptional-idrs/fig6_ms1/hmm_diff_plot.pdf'))){
  ggsave('hmm_diff_plot.pdf', hmm_v_scram, height = 8, width = 8, units = 'cm')
}



#used from some of the computational modelling
if(!(file.exists('~/post-transcriptional-idrs/processed_scores/scram_vs_hmm.csv'))){
  write.csv(master, file = 'scram_vs_hmm.csv',row.names = F)
}
  

