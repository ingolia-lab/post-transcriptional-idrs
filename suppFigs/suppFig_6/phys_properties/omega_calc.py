import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings(action='ignore')




dms_scores = pd.read_csv("~/yeast-idr-analysis-main/processed_scores/dms_scores.csv")
compo= pd.read_csv("~/yeast-idr-analysis-main/HMM_analysis/composition.csv")
dms_trans = pd.read_csv('~/yeast-idr-analysis-main/screening_lib_seqs/mutational_library_translation.csv')

#merge dataset, select peptides that are stable, and only include WT and scrambles of hte composition peps
dms_scores = pd.merge(dms_scores, dms_trans, left_on='Name', right_on='name', how='inner')
dms_scores = dms_scores.drop(['Unnamed: 0','name'], axis=1)



dms_scores['parent_pep'] = dms_scores['Name'].str.split('_').str[:3].str.join('_') #get parent name
compos_only = dms_scores[dms_scores['parent_pep'].isin(compo['pep'])]
compos_only = compos_only[compos_only['avg_stability'] >= -1] #filter for stability
compos_only = compos_only[compos_only['Name'].str.contains('scramble') | 
                          compos_only['Name'].str.contains('WT')]

compos_only['mutation_type'] = ['WT' if 'WT' in name else 'scramble' for name in compos_only['Name']]

#first, calculate the number of aromatics
compos_only['num_aromatics'] = [sum(seq.count(aro) for aro in ('F','Y','W')) for seq in compos_only['seq']]
#eliminate ones with fewer than 5 aromatics (i.e remove with less than 10%)
filtered = compos_only[compos_only['num_aromatics'] >=5]

def aro_deviation(df, blob_size = 5, get_max=False):
    devs_all = []
    for idx, row in df.iterrows():
        sigma_local = []
        #if get_max, organizes seq so that all the aromatics are at the end. otherwise, use the seq
        if get_max:
            seq = ''.join(sorted(row['seq'], key=lambda x: (x in 'YFW', x)))
        else:
            seq = row['seq']

        #calculate the global sigma of a sequence. constant irrespective of aromatic location
        sigma_total = ((row['num_aromatics']/len(seq))- (1-(row['num_aromatics']/len(seq))))**2
        
        #calc average local sigma (F_aro in blob - F_other in blob)^2
        for j in range(0,(len(seq)-blob_size+1)):
            #calc number of aros in local blob
            aro_local = sum(seq[j:(j+blob_size)].count(aro) for aro in ('F','Y','W'))
            #calc local sigma for blob (frac aro_local - frac_other_local)^2. append to list
            sigma_local.append(((aro_local/blob_size)-(1-(aro_local/blob_size)))**2)
            
        #get deviation. sum of (sig_local-sig_total)^2 / N_blobs, for each sig_local 
        devs_all.append(sum([(local - sigma_total)**2 for local in sigma_local])/(len(seq)-blob_size+1))
        sigma_local=[]
    return(devs_all)
        
        
        
filtered['sigma_seq'] = aro_deviation(filtered)  
filtered['sigma_max'] = aro_deviation(filtered, get_max=True)
filtered['omega_aro'] = filtered['sigma_seq']/filtered['sigma_max']

output = pd.DataFrame({'Name': filtered['Name'],'activity':filtered['avg_activity'],'omega_aro':filtered['omega_aro']})

#output df
path = os.path.join(os.getcwd(),'omega_vals.csv')
output.to_csv(path,index=False)
