import numpy as np
import pandas as pd
import sys 
import os 
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import average_precision_score, make_scorer
import pickle
import itertools

with open('x_train_compo_dipep.pkl', 'rb') as file:
	x_train = pickle.load(file)

with open('y_train_compo_dipep.pkl', 'rb') as file:
	y_train = pickle.load(file)


logistic_model = LogisticRegression()
scorer = make_scorer(average_precision_score)

#enumerate single and dipep labels
aas = ['D','E','K','R','H','Q','N','C','S','T','A','V','M','L','I','F','Y','W','G','P']
dipeps = [''.join(comb) for comb in itertools.product(aas, repeat=2)]
labels = dipeps + aas


logistic_model = LogisticRegression()
scorer = make_scorer(average_precision_score)

param_grid = {
    'C': np.logspace(-5,-1,20),
    'class_weight': [{0: 1, 1: weight} for weight in np.logspace(-2,2,20)],
    'penalty': ['l2'],
    'solver': ['liblinear'],
    'fit_intercept' : [True],
    'max_iter' : [500]
}

# create gridsearch model
grid_search = GridSearchCV(logistic_model, param_grid, cv=5, scoring=scorer, verbose=2, n_jobs=-1)

# fit gridsearch
grid_search.fit(x_train[labels], y_train)

# Print the best parameters and the corresponding AUC-PR score
print("Best parameters found:", grid_search.best_params_)
print("Best AUC-PR score:", grid_search.best_score_)

with open('composition_dipep_gridSearch.pkl', 'wb') as file:
    pickle.dump(grid_search, file)










