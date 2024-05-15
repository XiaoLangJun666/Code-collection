import pandas as pd
import numpy as np
import json
from collections import OrderedDict
# a=pd.read_csv('finalBlosum62/' + "O00311" + '.csv',sep=',')
#
# a=a.values
# b=a[:,1:]
#
# print(b)
# print(b.shape)





#
#     proteins = json.load(open('data/kiba/' + "proteins.txt"), object_pairs_hook=OrderedDict)
#     for t in proteins.keys():
#
#             if prot == proteins[t]:
#                 data = pd.read_csv('finalBlosum62/' + t + '.csv',sep=',')
#                 data = data.values
#                 matrix=data[:20,1:]

# matrix = np.zeros((20, 20),dtype="float64")
# proteins = json.load(open('data/kiba/' + "proteins.txt"), object_pairs_hook=OrderedDict)
# for t in proteins.keys():
#
#             if prot == proteins[t]:
#                 data = pd.read_csv('finalBlosum62/' + t + '.csv',sep=',')
#                 data = data.values
#                 #matrix = np.zeros((5000, 20))
#                 index = 0
#                 protein_len = []
#                 for str in data:
#                     if index<=19:
#                         a = str[0].split(' ')
#                         a = a[0]
#                         array = [float(x) for x in str[1:] if x != '']
#                         array = np.array(array[:20])
#                         matrix[index, :] = array
#                         index = index + 1

# data=pd.read_csv('finalBlosum62/'+"O00141"+".csv",sep=',')
# data=data.values
# index=0
# while index<=19:
#         matrix[index,:]=data[index,1:]
#         index=index+1




# proteins = json.load(open('data/kiba/' + "proteins.txt"), object_pairs_hook=OrderedDict)
# pt={}
# dt={}
# for t in proteins.keys():
#     pt[proteins[t]]=t
#     dt[t]=np.zeros([20,20])
#
#
# def seq_cat(prot):
#
#     na=pt[prot]
#     matrix=dt[na]
#     matrix=matrix[0:20,:]
#
#     return matrix
#
# df = pd.read_csv('data/' + 'kiba' + '_train.csv')
# test_prots=list(df['target_sequence'])
# XT = [seq_cat(t) for t in test_prots]

import numpy as np
import pandas as pd
import sys, os
from random import shuffle
import torch
import torch.nn as nn
from models.gat import GATNet
from models.gat_gcn import GAT_GCN
from models.gcn import GCNNet
from models.ginconv import GINConvNet
from utils import *
from torch_geometric.data import DataLoader
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import roc_curve, confusion_matrix, precision_score, recall_score, auc
import time

da = pd.read_csv("data/kinase/kinase_grayscale_cold_protein_balance_train.csv", sep=',')
da = da.values
dt = {}
b = 0
for i in da[:, 6]:
    if i not in dt.keys():
        dt[i] = da[:, 2][b]
    b = b + 1
df = pd.read_csv('data/' + "kinase" + '/kinase_grayscale_cold_protein_balance_test.csv')
df = df.values
c = 0
for k in df[:, 2]:
     if k not in dt.keys():
         dt[k] = df[:, 3][c]
     c = c + 1

dh = pd.read_csv('data/' + "kinase" + '/kinase_grayscale_cold_protein_balance_valid.csv')
dh = dh.values
h = 0
for j in dh[:, 6]:
     if j not in dt.keys():
         dt[j] = dh[:, 2][h]
     h = h + 1



seq_voc = "ABCDEFGHIKLMNOPQRSTUVWXYZ"
seq_dict = {v:(i+1) for i,v in enumerate(seq_voc)}
seq_dict_len = len(seq_dict)



dc={}

for t in dt.keys():



          data1 = pd.read_csv('FinalBLOSUM/' + t + ".csv", sep=',')
          data1 = data1.values
          matrix1=data1[:len(dt[t]),1:]



          data2 = pd.read_csv('pssm/' + t + '.fasta.pssm', sep='\t', header=None)
          data2 = data2.values
          matrix2 = np.zeros((data2.shape[0], 20), dtype="float64")
          index = 0
          for ds in data2:
                array = [float(x) for x in ds[1:] if x != '']
                array = np.array(array)
                matrix2[index, :] = array
                index = index + 1



          data3 = pd.read_csv("psi-pred/" + t + ".ss2")
          data3 = data3.values
          matrix3 = np.zeros([data3.shape[0], 3])
          x = 0
          y = 0
          for j in data3:
            for i in range(0, len(str(j))):
                if str(j)[i] == ".":
                    s = str(j)[i - 1:i + 4]
                    matrix3[x, y] = s
                    y = y + 1
            y = 0
            x = x + 1



          data4 = pd.read_csv('ASAquick/asaq.' + t + ".fasta/asaq.pred", sep=',', header=None)
          data4 = data4.values
          matrix4 = np.zeros(([data4.shape[0] - 1, 1]))
          for i in range(data4.shape[0] - 1):
              matrix4[i] = float(str(data4[i]).split(' ')[2])

          matrix6 = np.zeros((len(dt[t]),1),dtype='float64')
          for i in range (0,len(dt[t])):
              matrix6[i]=seq_dict[dt[t][i]]


          dc[t] = np.concatenate((matrix2,matrix6), 1)
          print(dc[t].shape)
#
#
def seq_cat(prot):
    matrix = []
    for t in dt.keys():
        if prot == dt[t]:
            if dc[t].shape[0] > 20:
                matrix = dc[t][:20, :]
            else:
                matrix = np.zeros((500, dc[t].shape[1]))
                matrix[:dc[t].shape[0],:]=dc[t]

    return matrix




df = pd.read_csv('data/' + "kinase" + '/kinase_grayscale_cold_protein_balance_train.csv')
train_prots= list(df['sequence.y'])
XT = [seq_cat(t) for t in train_prots]

print(XT)
# XT = []
# for t in train_prots:
#     XT0 = []
#     XT0.append(seq_cat1(t))
#     XT0.append(seq_cat2(t))
#     XT.append(XT0)
train_prots= np.asarray(XT)



train_prots=np.asarray(XT)
#
#
#
#
#
#
#

# data4=pd.read_csv('ASAquick/asaq.'+"O00141"+".fasta/asaq.pred",sep=',',header=None)
# data4=data4.values
# matrix4=np.zeros(([data4.shape[0]-1, 1]))
# for i in range(data4.shape[0]-1):
#     matrix4[i]=float(str(data4[i]).split(' ')[2])
# print(matrix4)
# print(data4)
#
# print(str(data4[2]).split(' ')[2])

# da=pd.read_csv("data/kinase/kinase_grayscale_cold_protein_balance_train.csv",sep=',')
# da=da.values
# dt={}
# b=0
# for i in da[:,6]:
#      if i not in dt.keys():
#          dt[i]=da[:,2][b]
#      b=b+1
# dc = pd.read_csv('data/' + "kinase" + '/kinase_grayscale_cold_protein_balance_test.csv',sep=',')
# dc = dc.values
# m = 0
# for e in dc[:, 2]:
#  if e not in dt.keys():
#      dt[e] = dc[:, 3][m]
#  m = m + 1
#
#
# df = pd.read_csv('data/' + "kinase" + '/kinase_grayscale_cold_protein_balance_test.csv',sep=',')
# df=df.values
# dn={}
# c=0
# for k in df[:,2]:
#      if k not in dn.keys():
#          dn[k]=df[:,3][c]
#      c=c+1
# a=True
# for t in dn.keys():
#     if dt[t]==dn[t]:
#         print(True)
#     else:
#         a=False
#         print(t)


# a=np.zeros(20)
# b=np.ones((20,5))
# print(b[:,:-1])
# print(b[:,-1])
# for k in range(0,b.shape[0]):
#     a[k]=b[k]
# print(a.shape,b.shape)
# print(a)
# print(a.shape)