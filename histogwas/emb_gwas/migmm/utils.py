import torch
import torch.nn as nn
import torch.nn.functional as F


class WPlatSoftmax(nn.Module):
    def __init__(self, dim=1):
        super(WPlatSoftmax, self).__init__()
        self.dim = dim
        self.shift = 0. #nn.Parameter(torch.zeros(1))

    def forward(self, x):
        sigmoid_x = self.raw_weights(x)
        sigmoid_sum = torch.sum(sigmoid_x, dim=self.dim, keepdim=True)
        output = sigmoid_x / sigmoid_sum
        return output

    def raw_weights(self, x):
        shifted_x = x + self.shift
        #sigmoid_x = 0.5 * (torch.tanh(shifted_x) + 1)
        #sigmoid_x = torch.sigmoid(shifted_x)
        #sigmoid_x = 0.5 * (torch.tanh(torch.exp(shifted_x)) + 1)
        sigmoid_x = F.softplus(shifted_x)
        return sigmoid_x


class WSoftmax(nn.Module):
    def __init__(self, dim=1):
        super(WSoftmax, self).__init__()
        self.softmax = nn.Softmax(dim=dim)
        self.shift = 0 #nn.Parameter(torch.zeros(1))

    def forward(self, x):
        shifted_x = x + self.shift
        return self.softmax(shifted_x)

    def raw_weights(self, x):
        shifted_x = x + self.shift
        return torch.exp(shifted_x)
        

def get_weight_act(activation, dim=1):
    if activation == 'softmax':
        return WSoftmax(dim=dim)
    elif activation == 'platsoftmax':
        return WPlatSoftmax(dim=dim)
    else:
        raise ValueError(f"Unsupported activation: {activation}")

