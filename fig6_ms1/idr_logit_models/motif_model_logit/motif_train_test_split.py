import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import pickle
import warnings
warnings.filterwarnings(action='ignore')


'''
This script takes the motif-dependent peptides and splits them into train/test. Same concepts as the composition
Does some pruning to remove inactives that are 'borderline'
1. removes YFP_mean activity scores between -1 and -0.75
2. removes inacts that are unstable (iRPF_score = -1)
3. removes inacts that overlap with repressors
4. inacts are highly overlapping; to balance dataset more effectively, remove inacts that overlap by > 30 with previous inact in an yorf
5. for motif dependent peptides, use fraction of AA rather than absolute number; helps normalize for length. In the composition comparison, everything is 50 AA
'''

#load WT scores and motif_dependent peptides
wt_scores = pd.read_csv("~/post-transcriptional-idrs/HMM_analysis/wt_scores_seqs.csv")
#load the 'motifs_full' to help eliminate fragments that overlap with the motif-dependent peps
motifs_full = pd.read_csv("~/post-transcriptional-idrs/HMM_analysis/motifs_full.csv")
wt_scores['repressor'] = [1 if score <= -1 else 0 for score in wt_scores['YFP_mean']]


#remove edge cases (-0.999 to -0.75), select repressors. removes the composition_dependnet peps 
wt_scores = wt_scores[(wt_scores['YFP_mean'] > -0.75) | (wt_scores['Peptide'].isin(motifs_full['name']))]
#remove unstable ones if they're inactive
wt_scores = wt_scores[(wt_scores['iRFP_mean'] >= -1) | wt_scores['Peptide'].isin(motifs_full['name'])]
repressors = wt_scores[wt_scores['repressor'] == 1]
inacts = wt_scores[wt_scores['repressor'] == 0]

#make sure inacts is sorted alphabetically by peptide
inacts = inacts.sort_values(by='Peptide')
inacts = inacts.reset_index(drop=True)

inacts_overlaps = []

#eliminate peptides that overlap with previous pep by > 30 .
for idx, row in inacts.iloc[:-1].iterrows():
    current_start = row['start']
    current_yorf = row['yorf']
    if idx < len(inacts)-2:
        next_row = inacts.loc[idx+1]
        two_row = inacts.loc[idx+2]

    if row['Peptide'] in inacts_overlaps:
        continue

    if row['yorf'] == next_row['yorf'] and next_row['start'] < (row['start'] + 19) and two_row['start'] < (next_row['start'] + 19):
        inacts_overlaps.append(next_row['Peptide'])


#remove peptides that overlap with both motif and composition dependent petides
overlap_with_repressor = []
for pep in tqdm(repressors['Peptide']):
    rep = pep.split('_')
    for inact in inacts['Peptide']:
        nonfunc = inact.split('_')
        if (rep[0] == nonfunc[0]): #check if same orf
            #check if the inact's start/stop are within the boundries of the repressors
            if (int(rep[2]) > int(nonfunc[1]) > int(rep[1])) or (int(rep[1]) < int(nonfunc[2]) < int(rep[2])):
                overlap_with_repressor.append(nonfunc)


#there will be duplicate here; remove
overlap_with_repressor = list(set(['_'.join(over) for over in overlap_with_repressor]))

#there are peptides that are both in inactivte overlaps and overlap with repressor. remove these
print(f'{len(overlap_with_repressor)} unique overlap with repressor, {len(inacts_overlaps)} unique inact overlaps')
inacts_to_remove = set(overlap_with_repressor + inacts_overlaps)
print(f'{len(inacts_to_remove)} total inacts to remove')

#remove overlaps and peps taht overlap with the motif-dependnet peps. load the submotifs and extract only the mutationally-sensitive regions to use as the repressors 
inacts_filtered = inacts.loc[~(inacts['Peptide'].isin(inacts_to_remove))] #filter inacts
submotifs = pd.read_csv('~/post-transcriptional-idrs/HMM_analysis/subMotifs_hmm.csv')	#load submotifs, select only the mutationally sensitive regions 

active_motifs = submotifs[(submotifs['state'] == 'active') & (submotifs['seq'].str.len() > 2)] #get active motifs and ensure > 2 residues
active_motifs['repressor'] = np.ones(len(active_motifs))

inact_df = pd.DataFrame({'pep':inacts_filtered['Peptide'], 'repressor':inacts_filtered['repressor'], 'seq':inacts_filtered['seq']})
active_df = pd.DataFrame({'pep':active_motifs['pep'], 'repressor':active_motifs['repressor'], 'seq': active_motifs['seq']})


filtered_df = pd.concat([inact_df,active_df], axis=0)

aas = ['D','E','K','R','H','Q','N','C','S','T','A','V','M','L','I','F','Y','W','G','P']

#count aas, convert into percentages for motifs
for aa in tqdm(aas):
    filtered_df[aa] = filtered_df['seq'].str.count(f'(?=({aa}))')

for aa in tqdm(aas):
    filtered_df[aa] /= filtered_df['seq'].apply(len)  

print(f'do all sum to 1 after norming?: {all(filtered_df[aas].sum(axis=1)-1 < 10e-10)}')

x_train, x_test, y_train, y_test = train_test_split(filtered_df[aas], filtered_df['repressor'],
                                                    test_size=0.20, random_state=18)

with open('x_train_motif.pkl', 'wb') as file:
    pickle.dump(x_train, file)

with open('x_test_motif.pkl', 'wb') as file:
    pickle.dump(x_test, file)

with open('y_train_motif.pkl', 'wb') as file:
    pickle.dump(y_train, file)

with open('y_test_motif.pkl', 'wb') as file:
    pickle.dump(y_test, file)
