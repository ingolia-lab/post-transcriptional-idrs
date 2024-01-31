import torch
import torch
import torch.nn as nn

#make model for predicting each AA, 
#uses same parameters, but need to split the operations to make a per-residue prediction, applying logit to each row of weidghted matrix

#note to user: load weights from mainFig_6 for LA, then run predictions. For figure with individual plots, were run on remote server with different dependencies. 

class logit_ORF(torch.nn.Module):
    # build the constructor
    def __init__(self, n_inputs, n_outputs):
        super(logit_ORF, self).__init__()
        self.weights = nn.Parameter(torch.rand(1,n_inputs)) #initialize random tensor, helps model converge faster 
        self.linear = torch.nn.Linear(n_inputs, n_outputs, bias=True)
        
    def weighted_embedding(self, x):
        weighting_vec = self.weights*x #multiply each vector by the trainable weight, 
        weighting_vec_sum = torch.sum(weighting_vec, dim=-1) #gets dot prod
        norm_weights = torch.nn.functional.softmax(weighting_vec_sum, dim=1)  # Normalize weights, use correct dimensions
        norm_weights_unsque = norm_weights.unsqueeze(-1) #format 
        weighted_x =  x*norm_weights_unsque #multiply matrix by weights
        return weighted_x         #no need to average here, since you want the Nx1280 represetionat
    
    def score_weighted(self, weighted_x):
        y_pred = torch.sigmoid(self.linear(weighted_x))
        return y_pred

    def forward(self, x):
        weighted_x = self.weighted_embedding(x)
        return self.score_weighted(weighted_x)
