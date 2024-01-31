import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
import itertools
from tqdm import tqdm
import pickle
from collections import Counter
import os
import warnings
warnings.filterwarnings(action='ignore')

'''
This script takes the composition-dependent peptides and splits them into train/test
Does some pruning to remove inactives that are 'borderline' 
1. removes YFP_mean activity scores between -1 and -0.75
2. removes inacts that are unstable (iRPF_score = -1)
3. removes inacts that overlap with repressors
4. inacts are highly overlapping; to balance dataset more effectively, remove inacts that overlap by > 30 with previous inact in an yorf
'''

#some formatting, remove ones that overalp with repressors
wt_scores = pd.read_csv("~/post-transcriptional-idrs/HMM_analysis/wt_scores_seqs.csv")
compo= pd.read_csv("~/post-transcriptional-idrs/HMM_analysis/composition.csv")
wt_scores['repressor'] = [1 if score <= -1 else 0 for score in wt_scores['YFP_mean']] 

#remove edge cases (-0.999 to -0.75), select repressors
wt_scores = wt_scores[(wt_scores['YFP_mean'] > -0.75) | (wt_scores['Peptide'].isin(compo['pep']))]
#remove unstable ones if they're inactive
wt_scores = wt_scores[(wt_scores['iRFP_mean'] >= -1) | wt_scores['Peptide'].isin(compo['pep'])]
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

#there will be duplicate values here if compositions share. Remove duplicates
overlap_with_repressor = list(set(['_'.join(over) for over in overlap_with_repressor]))

#there are peptides that are both in inactivte overlaps and overlap with repressor. remove these
print(f'{len(overlap_with_repressor)} unique overlap with repressor, {len(inacts_overlaps)} unique inact overlaps')
inacts_to_remove = set(overlap_with_repressor + inacts_overlaps)
print(f'{len(inacts_to_remove)} total inacts to remove')


#remove overlaps. wt_scores already filtered for composition peps
filtered = wt_scores.loc[~(wt_scores['Peptide'].isin(inacts_to_remove))] 
print(f'{len(filtered)} peps to use, repressors: {sum(filtered["repressor"])}, inacts: {sum(filtered["repressor"] == 0)}')



#enumerate dipeptides and single AAs
aas = ['D','E','K','R','H','Q','N','C','S','T','A','V','M','L','I','F','Y','W','G','P']
dipeps = [''.join(comb) for comb in itertools.product(aas, repeat=2)]
labels = dipeps + aas



#count occurrences of each single AA and dipep
for lab in tqdm(labels):
    filtered[lab] = filtered['seq'].str.count(f'(?=({lab}))') 

print(f'99 dipeps+singles for all seqs?: {(filtered[labels].sum(axis=1) == 99).all()}')


x_train, x_test, y_train, y_test = train_test_split(filtered[labels], filtered['repressor'], 
                                                    test_size=0.20, random_state=18)

with open('x_train_compo_dipep.pkl', 'wb') as file:
    pickle.dump(x_train, file)

with open('x_test_compo_dipep.pkl', 'wb') as file:
    pickle.dump(x_test, file)

with open('y_train_compo_dipep.pkl', 'wb') as file:
    pickle.dump(y_train, file)

with open('y_test_compo_dipep.pkl', 'wb') as file:
    pickle.dump(y_test, file)
