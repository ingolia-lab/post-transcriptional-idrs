import pickle
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
import itertools
from sklearn.metrics import precision_recall_curve,auc
import os

#load the gridsearchvals and test sets. 
#fit model, export the PRscores as a file to then use in R for plotting
with open('composition_dipep_gridSearch.pkl', 'rb') as file:
    grid_compo_dipep = pickle.load(file)

with open('x_train_compo_dipep.pkl', 'rb') as file:
    x_train = pickle.load(file)

with open('y_train_compo_dipep.pkl', 'rb') as file:
    y_train = pickle.load(file)

with open('x_test_compo_dipep.pkl', 'rb') as file:
    x_test = pickle.load(file)

with open('y_test_compo_dipep.pkl', 'rb') as file:
    y_test = pickle.load(file)




#get aas
aas = ['D','E','K','R','H','Q','N','C','S','T','A','V','M','L','I','F','Y','W','G','P']
dipeps = [''.join(comb) for comb in itertools.product(aas, repeat=2)]
labels = dipeps + aas

best_vals = grid_compo_dipep.best_params_

#get nbest params for model, fit to training, predict values and write to file for plotting in R
logreg_compo_dipep =  LogisticRegression(penalty = best_vals['penalty'], solver=best_vals['solver'], 
                                 class_weight= best_vals['class_weight'], 
                                 C = best_vals['C'],
                                 fit_intercept=best_vals['fit_intercept'], 
                                 max_iter=best_vals['max_iter'])

logreg_compo_dipep.fit(x_train[labels], y_train)
y_prob = logreg_compo_dipep.predict_proba(x_test[labels])

precision, recall, _ = precision_recall_curve(y_test, y_prob[:,1])

#save to dataframe for plotting in R, write model out too. 
df = pd.DataFrame({'precision': precision, 'recall':recall, 'AUPRC' : auc(recall, precision)})

f_path = os.path.join(os.getcwd(),'pr_compo_dipep.csv')

df.to_csv(f_path, index=False)

with open('compo_dipep_model_trained.pkl','wb') as file:
	pickle.dump(logreg_compo_dipep, file)

coef_double = pd.DataFrame({'dipep' : labels, 'coef' : logreg_compo_dipep.coef_[0]})

f_path = os.path.join(os.getcwd(),'compo_dipep_coef.csv')
coef_double.to_csv(f_path, index=False)


