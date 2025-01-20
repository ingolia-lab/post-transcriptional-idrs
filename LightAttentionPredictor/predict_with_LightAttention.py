import torch
import torch.nn as nn
import torch.optim as optim
import esm 

'''
script with example to run prediction on your_seq

all models use same learned weights but with slightly different outputs based on user's end goal, applying a sigmoid or linear activation
function to an ESM representation of a given protein sequence

logit_attentionV3 returns logit prediction of a specific sequence being a repressor (closer to 1  = repressor, closer to 0 = inactivate)
linear_attentionV3 returns a linear prediction of the repressive sequence. Helpful for generating mutational series
logit_perResi returns logit prediction for each residue to a given sequence. Helpful to examine which residues contribute to repressor function 

user can define which model to use as predictor in predict_effect
outputs (my_seq, tensor_of_prediction). tensor can be parsed later for downstream analysis. 
'''


#using Sgn1_151_200 as example, replace with your sequence of interest
my_seq = 'NSQNQYFQQWQWNYPLMAYPNPDTFPYYPPYPPNQSPNQNFGYNKNNYYR'

#check if cuda is support and use cuda
print('IS CUDA supported?:' , torch.cuda.is_available())
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f'using {device} ...')

# Load esm model and predictor model. Model trained with esm1b_t33_650M_UR50S 
esm_model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
esm_model.to(device) #transfer model to GPU

batch_converter = alphabet.get_batch_converter()
esm_model.eval()  # disables dropout for deterministic results

#returns single value predicting if a sequence is a repressor or not
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


#returns single value predicting if sequence is a repressor, using linear activation function. More helpful for finding
#mutations that reduce, but not eliminate, activity.
class linear_attentionV3(torch.nn.Module):
    # build the constructor
    def __init__(self, n_inputs, n_outputs):
        super(linear_attentionV3, self).__init__()
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
        #use linear to understand variants that reduce function but stil be active enought to  otherwise be called repressor
        y_pred = self.linear(weighted_avg) #replace sigmoid with linear acitvation function.
        return y_pred

#returns a n x 1 array where each n corresponds to a logit prediction of each residue's probability of being a repressive
class logit_perResi(torch.nn.Module):
    # build the constructor
    def __init__(self, n_inputs, n_outputs):
        super(logit_perResi, self).__init__()
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


###############
#specify which model you want to use for the predictor. Setting to logit_attentionV3
predictor = linear_attentionV3(n_inputs=1280, n_outputs=1)
#########################


#load weights, assumes weights in same folder but can change paths
predictor.load_state_dict(torch.load('LA_model_wScrams.pth',map_location=torch.device(device)))
if not next(predictor.parameters()).is_cuda:
    predictor.to(device)
    print(f'predictor transferred to {device}')


#to predict score: 
def predict_effect(my_seq,predictor=predictor, esm_model=esm_model):
    #format to make prediction
    batch_labels, batch_strs, batch_tokens = batch_converter([('',my_seq)])
    #make representation
    with torch.no_grad():
        results = esm_model(batch_tokens.to(device), repr_layers=[33]) #transfer to device
    token_representations = results["representations"][33]
    #note, representations have an extra column at start/end for padding, remove these. 
    token_representations = token_representations[:,1:-1,:]
    #score representation
    y = predictor(token_representations)
    return(my_seq,y) #note, tensor object. can parse later


#EXAMPLE:
my_prediction = predict_effect(my_seq)
print(my_prediction)
