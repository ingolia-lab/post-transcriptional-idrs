#load composition peps
import pandas as pd
import numpy as np 
import pickle
import os
import warnings
warnings.filterwarnings(action='ignore')


#load composition peps and models
home_dir = os.path.expanduser('~')
path = home_dir + '/post-transcriptional-idrs/fig6_ms1/idr_logit_models/composition_model_logit_singles/compo_single_model_trained.pkl'
compo_model = pickle.load(open(path,'rb'))
compo_peps = pd.read_csv("~/post-transcriptional-idrs/HMM_analysis/composition.csv")


#load motif peps and models
path = home_dir + '/post-transcriptional-idrs/fig6_ms1/idr_logit_models/motif_model_logit/motif_model_trained.pkl'
motif_model = pickle.load(open(path, 'rb'))
submotifs = pd.read_csv('~/post-transcriptional-idrs/HMM_analysis/subMotifs_hmm.csv')     
active_motifs = submotifs[(submotifs['state'] == 'active') & (submotifs['seq'].str.len() > 2)]

#add AAs
aas = ['D','E','K','R','H','Q','N','C','S','T','A','V','M','L','I','F','Y','W','G','P']

#count aas, convert into percentages for motifs
for aa in aas:
    compo_peps[aa] = compo_peps['seq'].str.count(f'(?=({aa}))')
    active_motifs[aa] = active_motifs['seq'].str.count(f'(?=({aa}))')


#run predictions using opposite model
compo_probs_motif_model = motif_model.predict_proba(compo_peps[aas])
motif_probs_compo_model = compo_model.predict_proba(active_motifs[aas])

compo_peps['motif_prob'] = compo_probs_motif_model[:,1]
active_motifs['compo_model'] = motif_probs_compo_model[:,1]


#write predictions out
path = home_dir + '/post-transcriptional-idrs/suppFigs/suppFig_6/cross_class_predictions/'
compo_path = os.path.join(path,'compo_pred_motif_model.csv')
motif_path = os.path.join(path,'motif_pred_compo_model.csv')

compo_peps.to_csv(compo_path, index=False)
active_motifs.to_csv(motif_path, index=False)
