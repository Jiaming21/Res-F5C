#!/usr/bin/env python
# coding: utf-8

import torch
import torch.nn as nn
import pandas as pd
import numpy as np

target_dir2 = '/home/jiaming/webserver/webserver_home' 
target_dir1 = pd.read_table(target_dir2+'/target_dir1.txt')

class ConvNet(torch.nn.Module):

    def __init__(self, num_classes):
        super(ConvNet, self).__init__()
        
        #########################
        ### 1st residual block
        #########################
        
        self.block_1 = torch.nn.Sequential(
                torch.nn.Conv1d(in_channels=1,
                                out_channels=4,
                                kernel_size=1,
                                stride=1,
                                padding=0),
                torch.nn.BatchNorm1d(4),
                torch.nn.ReLU(inplace=True),
                torch.nn.Conv1d(in_channels=4,
                                out_channels=1,
                                kernel_size=3,
                                stride=1,
                                padding=1),
                torch.nn.BatchNorm1d(1)
        )
        
        self.block_2 = torch.nn.Sequential(
                torch.nn.Conv1d(in_channels=1,
                                out_channels=4,
                                kernel_size=1,
                                stride=1,
                                padding=0),
                torch.nn.BatchNorm1d(4),
                torch.nn.ReLU(inplace=True),
                torch.nn.Conv1d(in_channels=4,
                                out_channels=1,
                                kernel_size=3,
                                stride=1,
                                padding=1),
                torch.nn.BatchNorm1d(1)
        )

        #########################
        ### Fully connected
        #########################        
        self.linear_1 = torch.nn.Linear(245, num_classes) 

        
    def forward(self, x):
        
        #########################
        ### 1st residual block
        #########################
        shortcut = x
        x = self.block_1(x)
        x = torch.nn.functional.relu(x + shortcut)
        
        #########################
        ### 2nd residual block
        #########################
        shortcut = x
        x = self.block_2(x)
        x = torch.nn.functional.relu(x + shortcut)
        
        #########################
        ### Fully connected
        #########################
        logits = self.linear_1(x.view(-1, 245))
        return logits

device = "cpu"

domain_feature = pd.read_csv('/home/jiaming/webserver/webserver_home/webserver_domain_encoding.csv').fillna(0)

sequence_feature = pd.read_csv('/home/jiaming/webserver/webserver_home/webserver_sequence_encoding_OH_ND.csv').fillna(0)

features = pd.merge(domain_feature, sequence_feature)
features = features.iloc[:,1:]
features = np.expand_dims(features, axis = -2)
features = torch.from_numpy(features.astype(np.float32)).to(device).type(torch.float)
print(features)

net = ConvNet(num_classes=2).to(device)
net.load_state_dict(torch.load('/home/jiaming/webserver/webserver_home/ResNet_weights.pth', map_location=lambda storage, loc: storage))

net.eval() 
with torch.no_grad():
    preds = net(features).detach().cpu().numpy()
    
print(preds)

preds = preds[:,1]
preds = pd.DataFrame(preds)

preds.to_csv('/home/jiaming/webserver/webserver_home/preds.csv', index=False)
print("preds has been saved to '/home/jiaming/webserver/webserver_home/preds.csv'")
