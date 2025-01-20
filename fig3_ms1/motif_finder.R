library(depmixS4)


setwd("~/post-transcriptional-idrs/processed_scores/")
scores = read.csv("dms_scores.csv")
hmm_acts = read.csv("scram_vs_hmm.csv") #if this doesn't exist, run script in fig6 (scramble_assignment.R)

setwd("~/post-transcriptional-idrs/screening_lib_seqs/")
aa_seq = read.csv("mutational_library_translation.csv")

#get motifs as ones that are sig diff in scramble and have negative scram vs hmm scores
motif = hmm_acts[hmm_acts$scram_diff >= 0 & hmm_acts$z_adj <= 0.05,] #excludes peps with more acitivty when scrambled, small num
composition = hmm_acts[hmm_acts$z_adj > 0.05, ]

format_data <-function(df=scores,pep, stability_cutoff = -1){
  foo = df[grepl(paste0(pep,'_Mut'), df$Name),]
  foo$start = sapply(seq(1,nrow(foo)),function(x){as.numeric(strsplit(foo$Name[x],'_')[[1]][5])}) 
  foo$cutoff_activity = ifelse(foo$avg_stability <= stability_cutoff, NA,foo$avg_activity) #replace activity scores of unstable with NA
  #if   missing peptides, fill in with Na
  missing = which(!(seq(1,46) %in% foo$start)) 
  empty = data.frame(matrix(ncol=ncol(foo),nrow=length(missing)))
  colnames(empty) = colnames(foo)
  empty$start = missing
  foo = rbind(foo,empty)
  return(foo[order(foo$start),])
} 

fit_hmm_with_aa <- function(df=scores, pep, aas=aa_seq, stability_cutoff = -1,states=2){
  print(pep)
  temp = format_data(df=df, pep=pep, stability_cutoff = stability_cutoff) #format data
  hmm = fit(depmix(temp$avg_activity~1, data=temp, nstates = states)) #fit to hmm, gets state values
  post = posterior(hmm, type='viterbi')[1] #get which position is in each state
  temp$state = post$state
  s1 = summary(hmm, which='all')[1] 
  s2 = summary(hmm, which='all')[2]
  temp$hmm_act = sapply(temp$state, function(x){ifelse(x == 1, s1,s2)}) #match activity to state
  temp$end = temp$start + 4 #for downstream plotting
  temp$aa = strsplit(aas[aas$name == paste0(pep,'_WT_0'),'seq'],'')[[1]][1:nrow(temp)]
  return(temp[complete.cases(temp),]) 
}

#build df with transition numbers and state len
get_HMM_params <- function(df=scores, pep_ls, aas = aa_seq, stability_cutoff = -1, states=2){
  params = data.frame(pep = rep(0,length(pep_ls)),
                      num_transitions = rep(0,length(pep_ls)),
                      min_s1_len = rep(0,length(pep_ls)),
                      min_s2_len = rep(0,length(pep_ls)))
  
  for(i in 1:length(pep_ls)){
    hmm = fit_hmm_with_aa(df=df, pep=pep_ls[i], aas=aa_seq, stability_cutoff = stability_cutoff,states = states)
    state_lens = with(rle(hmm$state),data.frame(values,lengths))
    params[i,] = list((pep_ls[i]), 
                sum(diff(hmm$state) != 0), #get number of transitions
                min(state_lens$lengths[state_lens$values == 1]), #get shortest state 1 
                min(state_lens$lengths[state_lens$values == 2])) #get shortest state 2
    
  }
  return(params)
}

motif_params = get_HMM_params(pep_ls=motif$name)
motif_params$size_crit = ifelse(motif_params$min_s1_len >= 2 & motif_params$min_s2_len >= 2, TRUE,FALSE )
motif_pass = motif_params[motif_params$num_transitions <= 5 & motif_params$size_crit,]


#this  divides the HMM fits into active/inactive regions, aware that the first residue of a functional motif
#starts at the last residue of the first inactivating mutation. 
find_motif_V4 <- function(df=scores, pep_ls,aas = aa_seq, stability_cutoff = -1,states=2,minlen=4){
  master = data.frame('pep' = NA, 'activity' = NA, 'seq' = NA, 'state'=NA)
  for(pep in pep_ls){
    #get the aa_sequence and initialize start at 1 
    pepseq = strsplit(aas[grepl(paste0(pep,'_WT_0'),aas$name),'seq'],'')[[1]]
    pos = 1

    #fit hmm, divide into regions within same state
    hmm = fit_hmm_with_aa(df=df, pep=pep, aas=aa_seq, stability_cutoff = stability_cutoff,states = states)
    #divide the hmm into list of separate regions of continuous activity. 
    hmm$region = c(0, cumsum(abs(diff(hmm$state)) > 0))
    regions = split(hmm, hmm$region)
    #here, multiple if loops evaluate different scenarios. 
    for(i in 1:length(regions)){
      #1: only evaluate if regions is longer that min specified
      if(nrow(regions[[i]]) > minlen){
        #for inert mutation (i.e is in the more repressive activity), collect all from start to stop +4. If its a deleterious, collect from start +4
        if(regions[[i]]$hmm_act[1] == min(hmm$hmm_act) & i != length(regions)){
          subseq = paste0(pepseq[pos:(regions[[i]]$start[nrow(regions[[i]])]+4)],collapse = '')
          master[nrow(master)+1,] = list(paste0(pep,'_',pos),
                                                min(hmm$hmm_act),
                                                subseq,'inactive')
          #reset the new starting position for the next motif, which must be a deleterious mutation
          pos = (regions[[i]]$start[nrow(regions[[i]])] + 5)
        }      
        #check if region consists of disruptive mutations. For deleterious mutations, begin counting from the 5th position
        #the last residue mutated of the first window is responsible for activity. 
        if(regions[[i]]$hmm_act[1] == max(hmm$hmm_act) & i != length(regions)){
          subseq = paste0(pepseq[pos:regions[[i]]$start[nrow(regions[[i]])]], collapse = '')
          master[nrow(master)+1,] = list(paste0(pep,'_',pos),
                                      max(hmm$hmm_act),
                                      subseq, 'active')
          pos = (regions[[i]]$start[nrow(regions[[i]])]+1)

        }
        #Assumes all residues in the last region, from x to 50, are in the same state
        if(i == length(regions)){
          subseq = paste0(pepseq[pos:50],collapse = '') #hardcode the 50, as all pepseqs are 50aa
          master[nrow(master)+1,] = list(paste0(pep,'_',pos),
                                      regions[[i]]$hmm_act[1],
                                      subseq,
                                      ifelse(regions[[i]]$hmm_act[1] == max(hmm$hmm_act),'active','inactive'))
        }
      }
    }
  }
  return(master[2:nrow(master),])
}
  
#active state means it has higher activity score (in this case, activity when mutated).Converse for inactive
submotifs = find_motif_V4(pep_ls = motif_pass$pep)
submotifs$len = nchar(submotifs$seq)


#qc = data.frame(pep = unique(sub("^(([^_]*_){2}[^_]*)_.*$", "\\1",submotifs$pep)), 
 #               assembled = rep(0,length(unique(sub("^(([^_]*_){2}[^_]*)_.*$", "\\1",submotifs$pep)))))
#qc$assembled = sapply(qc$pep,function(x){paste0(submotifs[grepl(x,submotifs$pep),'seq'],collapse = '')})
#qc$len = nchar(qc$assembled) #check that there aren't any super long submotifs, looks good




setwd("~/post-transcriptional-idrs/fig3_ms1/")

write.csv(submotifs,'subMotifs_hmm.csv',row.names = F)
write.csv(motif, 'motifs_full.csv',row.names = F) #also write the full motif data

#get the composition seq

compo_output = data.frame(pep=composition$name,
                          activity = composition$act_state,
                          seq = aa_seq[aa_seq$name %in% paste0(composition$name,'_WT_0'),'seq'],
                          state = rep('active',nrow(composition)),
                          len = rep(50,nrow(composition)))

if(!(file.exists('~/post-transcriptional-idrs/fig5_ms1/'))){
  write.csv(compo_output,'~/post-transcriptional-idrs/composition.csv',row.names = F)
}


