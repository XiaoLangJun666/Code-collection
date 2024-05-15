import pandas as pd
from scipy.stats import pearsonr
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import LabelEncoder
import torch.optim as optim
from sklearn.ensemble import AdaBoostClassifier
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn import linear_model
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import accuracy_score
import math
from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict
from sklearn.metrics import accuracy_score, confusion_matrix
import pickle
import shap
from sklearn.preprocessing import OneHotEncoder
import  seaborn as sns
import matplotlib.pyplot as plt
import catboost as cb
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score

np.set_printoptions(threshold=np.inf)

df=pd.read_csv('H:/FYP/New_Data/Chayi_gene/Gene_matrix_40', index_col=[0], sep='\t')
df=df.T

x=df.drop(["Label"],axis=1)
df1=np.zeros((x.shape[0],x.shape[1]),dtype=np.float32)
for a in range(0,x.shape[0]-1):
    for b in range(0,x.shape[1]):
        if float(x.iloc[a,b])<=1:
            x.iloc[a,b]=1
        df1[a,b]=math.log2(float(x.iloc[a,b]))

y=df['Label']
# df2=pd.DataFrame(df1,columns=x.columns,index=x.index)

label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y)
print(y)
print(y_encoded)

df2 = y_encoded.astype(np.int64)

X_train,X_test,y_train,y_test=train_test_split(df1,df2,test_size=0.2,random_state=42)





#MLP


train_X=torch.tensor(X_train)
test_X=torch.tensor(X_test)


Y_train = torch.tensor(y_train.squeeze())
Y_test = torch.tensor(y_test.squeeze())

dataset_train=TensorDataset(train_X,Y_train)
dataset_test=TensorDataset(test_X,Y_test)

train_data=DataLoader(dataset_train,batch_size=64)
test_data=DataLoader(dataset_test,batch_size=64)




# Define MLP model
class MLP(nn.Module):
    def __init__(self):
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(807, 128)
        self.fc2 = nn.Linear(128, 33)

    def forward(self, x):
        x=self.fc1(x)
        x = nn.functional.sigmoid(x)
        x = self.fc2(x)
        return x

# Set device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Define hyperparameters
learning_rate = 0.001
weight_decay = 0.001
num_epochs = 50

# Define model, loss function, and optimizer
model = MLP().to(device)
criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate, weight_decay=weight_decay)



train_acc=[]
train_loss=[]
test_Acc=[]
test_Loss=[]
# Train the model
for epoch in range(num_epochs):
    running_loss = 0.0
    total_correct = 0
    total_samples = 0
    for i, (inputs, labels) in enumerate(train_data):
        inputs = inputs.to(device)
        labels = labels.to(device)

        # Forward pass
        outputs = model(inputs.float())
        loss = criterion(outputs, labels)

        # Backward pass and optimize
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        # Compute running loss and accuracy
        running_loss += loss.item() * inputs.size(0)
        _, predicted = torch.max(outputs.data, 1)
        total_correct += (predicted == labels).sum().item()
        total_samples += labels.size(0)

    # Compute and print epoch statistics
    epoch_loss = running_loss / len(train_data.dataset)
    epoch_acc = 100.0 * total_correct / total_samples
    train_acc.append(epoch_acc)
    train_loss.append(epoch_loss)
    print(f"Epoch [{epoch+1}/{num_epochs}], Train Loss: {epoch_loss:.4f}, Train Acc: {epoch_acc:.2f}%")

# Test the model
model.eval()
with torch.no_grad():
    running_loss = 0.0
    total_correct = 0
    total_samples = 0
    for inputs, labels in test_data:
        inputs = inputs.to(device)
        labels = labels.to(device)

        # Forward pass
        outputs = model(inputs.float())
        loss = criterion(outputs, labels)

        # Compute running loss and accuracy
        running_loss += loss.item() * inputs.size(0)
        _, predicted = torch.max(outputs.data, 1)
        total_correct += (predicted == labels).sum().item()
        total_samples += labels.size(0)

    # Compute and print test set statistics
    test_loss = running_loss / len(test_data.dataset)
    test_acc = 100.0 * total_correct / total_samples

    test_Acc.append(test_acc)
    test_Loss.append(test_loss)
    print(f"Test Loss: {test_loss:.4f}, Test Acc: {test_acc:.2f}%")

model.eval()

# Make predictions on the test set
y_pred = []
y_true = []
with torch.no_grad():
    for x, y in test_data:
        x = x.to(device)
        y = y.to(device)

        output = model(x)
        pred = output.argmax(dim=1, keepdim=True)
        y_pred.extend(pred.squeeze().tolist())
        y_true.extend(y.squeeze().tolist())

# Compute the confusion matrix
cm = confusion_matrix(y_true, y_pred)
print(cm)
sensitivity_scores = recall_score(y_true, y_pred, average=None)
print(sensitivity_scores)
precision = precision_score(y_true, y_pred, average=None)
print(precision)

