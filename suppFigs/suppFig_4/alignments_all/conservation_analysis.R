library(msa)
library(ggplot2)

#load the blosum_matrix and sgd features
data(BLOSUM62)
sgd <- read.delim("~/post-transcriptional-idrs/processed_scores/SGD_features.tab", header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))



#load score, format yorf and match sgid to orf
setwd("~/yeast-idr-analysis-main/processed_scores/")
scores = read.csv("wt_scores_sortSeq.csv")

#create column that has sgdid;; a bit annoying and slow
scores$yorf = toupper(substr(scores$yorf,1,nchar(scores$yorf)-1)) #remove the 'p' at the end
scores$sgid = rep(0,nrow(scores))
for(i in 1:nrow(scores)){
  if(scores$yorf[i] %in% sgd$gene){
    scores$sgid[i] = sgd[sgd$gene == scores$yorf[i],'name']
  }else{
    scores$sgid[i] = scores$yorf[i] 
  }
}

setwd("~/post-transcriptional-idrs/suppFigs/suppFig_4/alignments_all/ungapped_alignments/")

#load alignment, make MSA, get conservation, map to Scer, remove gaps, return scores for Scer
make_cons_dfs = function(f_name, species='Scer'){
  seqs = readAAStringSet(f_name)
  align = msa(seqs)
  cons_score = msaConservationScore(align, BLOSUM62)
  #makes df of seq, conservation score for the Scer 
  temp = data.frame('seq' = unlist(strsplit(as.character(align@unmasked[species]),'')),
                    'score' = cons_score)
  #remove gaps
  out = temp[temp$seq != '-',]
  out$idx = seq(1,nrow(out))
  return(out)
}

#generate all the conservation scores for what we have, this takes quite a bit of time
if(!(dir.exists('conservation_scores'))){
  dir.create('conservation_scores')
  f_names = list.files(pattern='*aa.fasta')
  for(f in f_names){
    x = make_cons_dfs(f)
    orf = sapply(strsplit(f,'\\.'),'[[',2)
    write.csv(x,paste0('conservation_scores/',orf,'_conScore.csv'))
  }
}


#get conservation score for each peptide, using mean
get_cons_scores = function(df=scores){
  out = df
  out$cons_score = rep(0,nrow(out))
  counter = 1
  for(sgid in unique(df$sgid)){
    print(sgid)
    #subset the data for sgid and load conservation scores, avoids slow loading/subsetting
    temp = df[df$sgid == sgid,]
    con_score = read.csv(paste0("conservation_scores/",sgid,'_conScore.csv'))

    for(i in 1:nrow(temp)){
      #average conservation score for a peptide, add to df using idx
      seq_score = mean(con_score[temp[i,]$start:temp[i,]$stop,'score'])
      out$cons_score[counter] = seq_score
      counter = counter + 1
    }
  }
  return(out)
}




#only keep scores that have an alignment
yorfs_present = sapply(strsplit(list.files(path='conservation_scores/',pattern='.csv'),'_'),'[[',1)
scores = scores[scores$sgid %in% yorfs_present,]

#get mean score of fragment and then assign as repressor
scores_cons = get_cons_scores(scores)
scores_cons$repressor = ifelse(scores_cons$YFP_mean <= -1, 'repressive_frag', 'inactive_frag')


#analyze residues within and outside motifs, check if within motif is more conserved
motifs = read.csv("~/post-transcriptional-idrs//HMM_analysis/subMotifs_hmm.csv")
#replace pep name with start stop. the last number is the 'start' position within the pep, and the length can be used
motifs$pep_start = as.numeric(sapply(strsplit(motifs$pep,'_'),'[[',2))
motifs$motif_start = motifs$pep_start + as.numeric(sapply(strsplit(motifs$pep,'_'),'[[',4)) -1
motifs$motif_stop = motifs$motif_start + motifs$len - 1
motifs$yorf = sapply(strsplit(motifs$pep,'_'),'[[',1)
motifs$yorf = toupper(substr(motifs$yorf,1,nchar(motifs$yorf)-1))
#exclude motifs < 3 aa, and then format the name. adds a sgid name for matching. remove NA (40 inact/act subfrags)
motifs = subset(motifs, len >= 3)
motifs$id = paste0(motifs$yorf,'_',motifs$motif_start,'_',motifs$motif_stop)

#remove motifs that don't have an sgid or an alignment
motifs$sgid = rep(0,nrow(motifs))
for(i in 1:nrow(motifs)){
  if(motifs$yorf[i] %in% sgd$gene){
    motifs$sgid[i] = sgd[sgd$gene == motifs$yorf[i],'name']
  }else{
    motifs$sgid[i] = NA
  }
}

#keeps 350/473
alignment_yorfs = sapply(strsplit(list.files('conservation_scores/',pattern = '.csv'),'_'),'[[',1)
motifs = motifs[complete.cases(motifs$sgid) & motifs$sgid %in% alignment_yorfs,]


#annoyingly rewritten to handle the motif dataframes, but works similar to get_cons_scores above.
#cleaner too, since not operating on all orfs at once. meant to be called in sapply
get_frag_cons = function(frag_name, df){
  temp = df[df$pep == frag_name,]
  conservations = read.csv(paste0("conservation_scores/",temp$sgid[1],'_conScore.csv'))
  score = mean(conservations[temp$motif_start[1]:temp$motif_stop[1],'score'])
  return(score)
}

#assign mean conservation score to each subfrag
motifs$cons_score = sapply(motifs$pep,function(x){get_frag_cons(x,motifs)})
motifs$state = ifelse(motifs$state=='inactive','outside_motif','within_motif')


#for plotting: make dataframe of full repressors and motifs/inact regions, plot all on same
full_df = data.frame(name=as.character(),cons_score=as.numeric(),id=as.character())
full_df = rbind(full_df,scores_cons[,c('Peptide','cons_score','repressor')])
temp = motifs[,c('pep','cons_score','state')]
colnames(temp) = c('Peptide','cons_score','repressor')
full_df = rbind(full_df,temp)

cons_plot = ggplot(full_df,aes(x=cons_score,color=repressor)) + 
  stat_ecdf(size=1) + scale_color_manual(values =c("#e41a1c",'#377eb8','#4daf4a','#984ea3')) + 
  scale_y_continuous(expand=c(0,0),limits=c(-0.01,1.01)) + theme_classic() + 
  xlab('Conservation Score') + ylab('Cumulative Probability') + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black",),
        axis.ticks.length = unit(0.1,'cm'),
        aspect.ratio = 1.0,
        panel.border = element_blank(), 
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
  )

if(!file.exists('../cons_ecdf.pdf')){
  ggsave('../cons_ecdf.pdf', cons_plot, height = 6,width = 6, units='in')
}


#stats....only within_motif different from others
combos = combn(unique(full_df$repressor),2)

ks_df <- data.frame(comparison=character(6), D=numeric(6), p=numeric(6), stringsAsFactors=F)

for(i in 1:6){  
  k <- ks.test(full_df[full_df$repressor==combos[1,i],'cons_score'],full_df[full_df$repressor==combos[2,i],'cons_score'])
  ks_df$comparison[i] = paste0(combos[1,i],'_',combos[2,i])
  ks_df$D[i] = k$statistic
  ks_df$p[i] = k$p.value
}















