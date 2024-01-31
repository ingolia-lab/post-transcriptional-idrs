import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import os
from sklearn.metrics import auc, precision_recall_curve
import pickle

#define the models
class logit_attentionV3(torch.nn.Module):
    # build the constructor
    def __init__(self, n_inputs, n_outputs):
        super(logit_attentionV3, self).__init__()
        self.weights = nn.Parameter(torch.rand(1,n_inputs)) #initialize random tensor, helps model converge faster 
        self.linear = torch.nn.Linear(n_inputs, n_outputs, bias=True)

    # make predictions
    #here, applying the weighting vector to the L2 normalized matrix
    def forward(self, x):
        weighting_vec = self.weights*x #multiply each vector by the trainable weight, 
        weighting_vec_sum = torch.sum(weighting_vec, dim=-1) #gets dot prod
        norm_weights = torch.nn.functional.softmax(weighting_vec_sum, dim=1)  # Normalize weights, use correct dimensions
        norm_weights_unsque = norm_weights.unsqueeze(-1) #format 
        weighted_x =  x*norm_weights_unsque #multiply matrix by weights
        weighted_avg = weighted_x.sum(dim=1) #get mean across column
        y_pred = torch.sigmoid(self.linear(weighted_avg))
        return y_pred

#user can load model and parameters and run on own device. 
