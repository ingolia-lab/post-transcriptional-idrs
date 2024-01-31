library(depmixS4)

setwd("~/post-transcriptional-idrs/processed_scores/")
dms_scores = read.csv("dms_scores.csv")
setwd("~/post-transcriptional-idrs/screening_lib_seqs/")
pep_seqs = read.csv("mutational_library_translation.csv")

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

format_hmm <- function(df=dms_scores, pep, stability_cutoff = -1,states=2){
  print(pep)
  temp = format_data(df=df, pep=pep, stability_cutoff = stability_cutoff) #format data
  hmm = fit(depmix(temp$avg_activity~1, data=temp, nstates = states)) 
  post = posterior(hmm, type='viterbi')[1]
  print(post)
  temp$state = post$state
  s1 = summary(hmm, which='all')[1] 
  s2 = summary(hmm, which='all')[2]
  temp$hmm_act = sapply(temp$state, function(x){ifelse(x == 1, s1,s2)}) #match activity to state
  temp$end = temp$start + 4 #for downstream plotting
  return(temp[complete.cases(temp),]) 
}
##############################################################
###collection of individual HMMs used in different figures#####
##############################################################

#hmm plots used in figures
pho92_df = format_hmm(pep='Pho92p_17_66')

pho92_plot <- ggplot(pho92_df, aes(x = start, y = cutoff_activity)) +
  geom_segment(aes(x=start,xend = end, y= cutoff_activity, yend = cutoff_activity), color = "#ef8a62", size=2) +
  geom_line(aes(x=start+2,y=cutoff_activity), size=1,color='#b2182b') + 
  geom_line(aes(x=start+2,y=hmm_act), size=1.5,color='#67a9cf') + 
  geom_errorbar(aes(x = start+2, ymin = cutoff_activity - activity_sd, ymax = cutoff_activity + activity_sd), size=0.25) + 
  xlab("Position") +ylab("Activity Score") + ylim(-2.1,1.3) + theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

caf20_df = format_hmm(pep='Caf20p_27_76')

caf20_plot <- ggplot(caf20_df, aes(x = start, y = cutoff_activity)) +
  geom_segment(aes(x=start,xend = end, y= cutoff_activity, yend = cutoff_activity), color = "#ef8a62", size=2) +
  geom_line(aes(x=start+2,y=cutoff_activity), size=1,color='#b2182b') + 
  geom_line(aes(x=start+2,y=hmm_act), size=1.5,color='#67a9cf') + 
  geom_errorbar(aes(x = start+2, ymin = cutoff_activity - activity_sd, ymax = cutoff_activity + activity_sd), size=0.25) + 
  xlab("Position") +ylab("Activity Score") + ylim(-2.1,1.3) + theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/post-transcriptional-idrs/fig4_ms1/")
if(!file.exists('pho92_hmm_plot.pdf') & !file.exists('caf20_hmm_plot.pdf')){
  ggsave('pho92_hmm_plot.pdf', pho92_plot, height = 4, width=6.5, units='cm')
  ggsave('caf20_hmm_plot.pdf', caf20_plot, height = 4, width=6.5, units='cm')
}


###for tis11 and Eap1, fig 5

#plot tis11 pep
tis11_df = format_hmm(pep='Tis11p_28_77')

tis11_plot <- ggplot(tis11_df, aes(x = start, y = cutoff_activity)) +
  geom_segment(aes(x=start,xend = end, y= cutoff_activity, yend = cutoff_activity), color = "#ef8a62", size=2) +
  geom_line(aes(x=start+2,y=cutoff_activity), size=1,color='#b2182b') + 
  geom_line(aes(x=start+2,y=hmm_act), size=1.5,color='#67a9cf') + 
  geom_errorbar(aes(x = start+2, ymin = cutoff_activity - activity_sd, ymax = cutoff_activity + activity_sd), size=0.25) + 
  xlab("Position") +ylab("Activity") + ylim(-2.1,1.1) + theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/post-transcriptional-idrs/fig4_ms1/")
if(!file.exists('tis11_hmm_plot.pdf')){
  ggsave('tis11_hmm_plot.pdf', tis11_plot, height = 4, width=6.5, units='cm')
}

eap1_m1_df = format_hmm(pep='Eap1p_299_348')

eap1_m1_plot <- ggplot(eap1_m1_df, aes(x = start, y = cutoff_activity)) +
  geom_segment(aes(x=start,xend = end, y= cutoff_activity, yend = cutoff_activity), color = "#ef8a62", size=2) +
  geom_line(aes(x=start+2,y=cutoff_activity), size=1,color='#b2182b') + 
  geom_line(aes(x=start+2,y=hmm_act), size=1.5,color='#67a9cf') + 
  geom_errorbar(aes(x = start+2, ymin = cutoff_activity - activity_sd, ymax = cutoff_activity + activity_sd), size=0.25) + 
  xlab("Position") +ylab("Activity Score") + ylim(-1,1.5) + theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

eap1_m2_df = format_hmm(pep='Eap1p_369_418')

eap1_m2_plot <- ggplot(eap1_m2_df, aes(x = start, y = cutoff_activity)) +
  geom_segment(aes(x=start,xend = end, y= cutoff_activity, yend = cutoff_activity), color = "#ef8a62", size=2) +
  geom_line(aes(x=start+2,y=cutoff_activity), size=1,color='#b2182b') + 
  geom_line(aes(x=start+2,y=hmm_act), size=1.5,color='#67a9cf') + 
  geom_errorbar(aes(x = start+2, ymin = cutoff_activity - activity_sd, ymax = cutoff_activity + activity_sd), size=0.25) + 
  xlab("Position") +ylab("Activity Score") + ylim(-1,1.5) + theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())

setwd("~/post-transcriptional-idrs/fig5_ms1/eap1_motif_hmms/")
if(!file.exists('eap1_motif_1.pdf') & !file.exists('eap1_motif_2.pdf')){
  ggsave('eap1_motif_1.pdf', eap1_m1_plot, height = 4, width=6.5, units='cm')
  ggsave('eap1_motif_2.pdf', eap1_m2_plot, height = 4, width=6.5, units='cm')
}


#plot sgn1_151_200
sgn1_df = format_hmm(pep='Sgn1p_151_200')


sgn1_plot <- ggplot(sgn1_df, aes(x = start, y = cutoff_activity)) +
  geom_segment(aes(x=start,xend = end, y= cutoff_activity, yend = cutoff_activity), color = "#ef8a62", size=2) +
  geom_line(aes(x=start+2,y=cutoff_activity), size=1,color='#b2182b') + 
  geom_line(aes(x=start+2,y=hmm_act), size=1.5,color='#67a9cf') + 
  geom_errorbar(aes(x = start+2, ymin = cutoff_activity - activity_sd, ymax = cutoff_activity + activity_sd), size=0.25) + 
  xlab("Position") +ylab("Activity Score") + ylim(-2.1,1.5) + theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.border = element_blank())


setwd("~/post-transcriptional-idrs/fig6_ms1/")
ggsave('sgn1_hmm_plot.pdf',sgn1_plot, height = 6.5, width=6.5, units='cm') #save as square



