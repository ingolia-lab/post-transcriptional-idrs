import os
import pandas as pd
import numpy as np
import localcider
from localcider.sequenceParameters import SequenceParameters

peps = pd.read_csv("~/post-transcriptional-idrs/processed_scores/wt_yfp_irfp_scores.csv")
#remove inactives that are unstable 
peps =  peps[(peps['avg_stability'] > -1) | (peps['avg_activity'] <= -1)]


seqs = pd.read_csv("~/post-transcriptional-idrs/screening_lib_seqs/wt_library_translation.csv")

peps =pd.merge(peps, seqs, on='Name', how='inner')
peps['param_seq'] = [SequenceParameters(x) for x in list(peps['seq'])]
peps['KD_hydropathy'] = [x.get_mean_hydropathy() for x in list(peps['param_seq'])]
peps['Y'] = [seq.count('Y') for seq in peps['seq']]
peps['F'] = [seq.count('F') for seq in peps['seq']]
peps['W'] = [seq.count('W') for seq in peps['seq']]

df = pd.DataFrame({'peps':peps['Name'], 'seq':peps['seq'],'activity_score' : peps['avg_activity'] ,'KD_hydro': peps['KD_hydropathy'], 'Y': peps['Y'], 'F': peps['F'], 'W': peps['W']})

path = os.path.join(os.getcwd(),'kd_hydro.csv')

df.to_csv(path, index=False)


