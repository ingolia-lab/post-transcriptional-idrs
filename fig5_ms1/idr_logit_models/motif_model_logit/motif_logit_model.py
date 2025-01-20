import pickle
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve,auc
import os

#load the gridsearchvals and test sets.
#fit model, export the PRscores as a file to then use in R for plotting
with open('motif_gridSearch.pkl', 'rb') as file:
    grid_motif = pickle.load(file)

with open('x_train_motif.pkl', 'rb') as file:
    x_train = pickle.load(file)

with open('y_train_motif.pkl', 'rb') as file:
    y_train = pickle.load(file)

with open('x_test_motif.pkl', 'rb') as file:
    x_test = pickle.load(file)

with open('y_test_motif.pkl', 'rb') as file:
    y_test = pickle.load(file)

aas = ['D','E','K','R','H','Q','N','C','S','T','A','V','M','L','I','F','Y','W','G','P']

best_vals = grid_motif.best_params_

#get best params for model, fit to training, predict values and write to file for plotting in R
logreg_motif =  LogisticRegression(penalty = best_vals['penalty'], solver=best_vals['solver'],
                                 class_weight= best_vals['class_weight'],
                                 C = best_vals['C'],
                                 fit_intercept=best_vals['fit_intercept'],
                                 max_iter=best_vals['max_iter'])


logreg_motif.fit(x_train[aas], y_train)
y_prob = logreg_motif.predict_proba(x_test[aas])

precision, recall, _ = precision_recall_curve(y_test, y_prob[:,1])
#write these out to dataframe for plotting in R. Then save the model to pickle
df = pd.DataFrame({'precision': precision, 'recall':recall, 'AUPRC' : auc(recall, precision)})

f_path = os.path.join(os.getcwd(),'pr_motif.csv')

df.to_csv(f_path, index=False)

#save model
with open('motif_model_trained.pkl','wb') as file:
	pickle.dump(logreg_motif, file)

coef_motif = pd.DataFrame({'aas' : aas, 'coef' : logreg_motif.coef_[0]})

f_path = os.path.join(os.getcwd(),'motif_coef.csv')
coef_motif.to_csv(f_path, index=False)
