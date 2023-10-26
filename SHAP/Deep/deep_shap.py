#!/usr/bin/env python
# coding: utf-8

# In[5]:


#import torch
#import torch.nn as nn
#import torch.optim as optim
#import torchvision
#import torchvision.transforms as transforms
#import shap
#
#torch.manual_seed(0)
#
#class SimpleNet(nn.Module):
#    def __init__(self):
#        super(SimpleNet, self).__init__()
#        self.fc1 = nn.Linear(28*28, 500)
#        self.fc2 = nn.Linear(500, 10)
#
#    def forward(self, x):
#        x = x.view(-1, 28*28)
#        x = torch.relu(self.fc1(x))
#        x = self.fc2(x)
#        return x
#
#transform = transforms.Compose([transforms.ToTensor(), transforms.Normalize((0.5,), (0.5,))])
#
#train_dataset = torchvision.datasets.FashionMNIST(root='./data', train=True, download=True, transform=transform)
#train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=32, shuffle=True)
#
#model = SimpleNet()
#criterion = nn.CrossEntropyLoss()
#optimizer = optim.Adam(model.parameters(), lr=0.001)
#
#for epoch in range(5):  # loop over the dataset multiple times
#    for inputs, labels in train_loader:
#        optimizer.zero_grad()
#        outputs = model(inputs)
#        loss = criterion(outputs, labels)
#        loss.backward()
#        optimizer.step()
#
#torch.save(model.state_dict(), "model.pth")
#
#loaded_model = SimpleNet()
#loaded_model.load_state_dict(torch.load("model.pth"))
#loaded_model.eval()
#
#background = train_dataset.data[:100].float().view(-1, 28*28) / 255.0  
#explainer = shap.DeepExplainer(loaded_model, background)
#test_data = train_dataset.data[100:110].float().view(-1, 28*28) / 255.0  
#shap_values = explainer.shap_values(test_data)
#
## 5. 使用shap.summary_plot
#shap.summary_plot(shap_values, test_data.detach().numpy())
#


# In[ ]:


import torch.nn as nn
import torch.optim as optim
import torch
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

class Net_conv(nn.Module):
    def __init__(self):
        super(Net_conv, self).__init__()
        self.conv1 = nn.Conv1d(1, 16, kernel_size=3, stride=1, padding=1)
        self.dropout1 = nn.Dropout(0.2) 
        self.conv2 = nn.Conv1d(16, 32, kernel_size=3, stride=1, padding=1)
        self.dropout2 = nn.Dropout(0.2) 
        self.fc1 = nn.Linear(1280, 64)
        self.dropout3 = nn.Dropout(0.5) 
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, 1)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.relu(self.conv1(x))
        x = self.dropout1(x) 
        x = self.relu(self.conv2(x))
        x = self.dropout2(x) 
        x = x.view(x.size(0), -1)
        x = self.relu(self.fc1(x)) 
        x = self.dropout3(x) 
        x = self.relu(self.fc2(x))
        x = self.sigmoid(self.fc3(x))
        return x

# Load data
pos = pd.read_csv('/Users/jiaming/Desktop/f5c/pos_domain_encoding.csv')
pos = pos.iloc[:, 1:]
neg = pd.read_csv('/Users/jiaming/Desktop/f5c/neg_domain_encoding.csv')
neg = neg.iloc[:, 1:]

raw_datas = np.concatenate((pos, neg), axis=0)
raw_labels = np.concatenate(([1] * pos.shape[0], [0] * neg.shape[0]), axis=0)

np.random.seed(1)
indices = np.random.permutation(raw_labels.shape[0])

X = raw_datas[indices, :]
y = raw_labels[indices]

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

X_train = torch.from_numpy(np.expand_dims(X_train.astype(np.float32), axis = -2)).to("mps")
X_test = torch.from_numpy(np.expand_dims(X_test.astype(np.float32), axis = -2)).to("mps")
y_train = torch.from_numpy(np.expand_dims(y_train.astype(np.float32), axis = -1)).to("mps")
y_test = torch.from_numpy(np.expand_dims(y_test.astype(np.float32), axis = -1)).to("mps")

model = Net_conv().to("mps")

# Hyperparameters
num_epochs = 10

criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

for epoch in range(num_epochs):
    total_loss = 0

    outputs = model(X_train)
    loss = criterion(outputs, y_train)

    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    
    total_loss += loss.item()
    
    print(f"Epoch [{epoch + 1}/{num_epochs}], Loss: {total_loss:.4f}")

import shap

# Create a SHAP explainer object
explainer = shap.GradientExplainer(model, X_train[:30]) # using a subset for efficiency

# Compute SHAP values for a sample of the test set
shap_values = explainer.shap_values(X_train[:30])
shap_values = np.squeeze(shap_values, axis=1)

print(shap_values.shape)
print(X_train[:30].squeeze(1).shape)

# Plot the SHAP values
fig1 = shap.summary_plot(shap_values, X_train[:30].squeeze(1), feature_names=pos.columns, show=False)
plt.savefig('fig1.pdf', bbox_inches='tight', format='pdf')
plt.close(fig1) 

fig2 = shap.summary_plot(shap_values, X_train[:30].squeeze(1), feature_names=pos.columns, plot_type="bar", show=False)
plt.savefig('fig2.pdf', bbox_inches='tight', format='pdf')
plt.close(fig2) 


# In[ ]:




