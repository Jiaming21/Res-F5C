#!/usr/bin/env python
# coding: utf-8

# In[18]:


#读取模型参数
import torch
import torch.nn as nn
import pandas as pd
import numpy as np


# In[19]:


class Net_conv(nn.Module):
    def __init__(self):
        super(Net_conv, self).__init__()
        self.conv1 = nn.Conv1d(1, 16, kernel_size=3, stride=1, padding=1)
        self.dropout1 = nn.Dropout(0.2)  # 添加第一个Dropout层
        self.conv2 = nn.Conv1d(16, 32, kernel_size=3, stride=1, padding=1)
        self.dropout2 = nn.Dropout(0.5)  # 添加第二个Dropout层
        self.fc1 = nn.Linear(7840, 64)
        self.dropout3 = nn.Dropout(0.5)  # 添加第三个Dropout层
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, 1)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.relu(self.conv1(x))
        x = self.dropout1(x)  # 在第一个卷积层后应用第一个Dropout层
        x = self.relu(self.conv2(x))
        x = self.dropout2(x)  # 在第二个卷积层后应用第二个Dropout层
        x = x.view(x.size(0), -1)
        x = self.relu(self.fc1(x)) 
        x = self.dropout3(x)  # 在第一个全连接层后应用第三个Dropout层
        x = self.relu(self.fc2(x))
        x = self.sigmoid(self.fc3(x))
        return x


# In[20]:


device = "cpu"


# In[51]:


domain_feature = pd.read_csv('/Users/jiaming/Desktop/f5c/webserver/webserver_domain_encoding.csv').fillna(0)

sequence_feature = pd.read_csv('/Users/jiaming/Desktop/f5c/webserver/webserver_sequence_encoding_OH_ND.csv') 

features = pd.merge(domain_feature, sequence_feature)
features = features.iloc[:,1:]
features = np.expand_dims(features, axis = -2)
features = torch.from_numpy(features.astype(np.float32)).to(device).type(torch.float)


# In[52]:


#加载模型参数
net = Net_conv().to(device)
net.load_state_dict(torch.load('/Users/jiaming/Desktop/f5c/webserver/model_parameters.pth'))

with torch.no_grad():
    net.eval()           
    preds = net(features)
    preds = preds.detach().cpu().numpy() 


# In[53]:


preds = pd.DataFrame(preds, columns=['preds'])
preds


# In[54]:


preds.to_csv('/Users/jiaming/Desktop/f5c/webserver/preds.csv', index=False)
print("preds has been saved to '/Users/jiaming/Desktop/f5c/webserver/preds.csv'")


# In[ ]:





# In[ ]:





# In[ ]:




