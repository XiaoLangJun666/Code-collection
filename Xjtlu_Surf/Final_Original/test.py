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

# print(TestbedDataset(root='data/processed', dataset="kinase" + '_grayscale_cold_protein_balance_train'))



da=pd.read_csv("data/kinase/kinase_grayscale_cold_protein_balance_train.csv",sep=',')
da=da.values
dt={}
b=0
for i in da[:,6]:
     if i not in dt.keys():
         dt[i]=da[:,2][b]
     b=b+1

# for t in dt.keys():
#     df = pd.read_csv("psi-pred/" + t + ".ss2")
#     d = df.values
#     n = np.zeros([d.shape[0], 3])
#     x = 0
#     y = 0
#     for j in d:
#         for i in range(0, len(str(j))):
#             if str(j)[i] == ".":
#                 s = str(j)[i - 1:i + 4]
#                 n[x, y] = s
#                 y = y + 1
#         y = 0
#         x = x + 1
dc={}

for t in dt.keys():
        data = pd.read_csv('pssm/' + t + '.fasta.pssm', sep='\t', header=None)
        data = data.values
        matrix = np.zeros((data.shape[0], 20), dtype="float64")
        index = 0
        for ds in data:
            array = [float(x) for x in ds[1:] if x != '']
            array = np.array(array)
            matrix[index, :] = array
            index = index + 1

        df = pd.read_csv("psi-pred/" + t + ".ss2")
        d = df.values
        n = np.zeros([d.shape[0], 3])
        x = 0
        y = 0
        for j in d:
            for i in range(0, len(str(j))):
                if str(j)[i] == ".":
                    s = str(j)[i - 1:i + 4]
                    n[x, y] = s
                    y = y + 1
            y = 0
            x = x + 1

        dc[t] = np.concatenate((matrix, n), 1)
print(dc)
# def seq_cat(prot):
#     matrix = []
#     for t in dt.keys():
#         print(dt[t])
#          # if prot == dt[t]:
#          #     print(1)
#     #         matrix = dc[t][0:20, :]
#     #         print(1)
#     return matrix
# df = pd.read_csv('data/' + "kinase" + '/kinase_grayscale_cold_protein_balance_train.csv')
# train_prots=list(df['sequence.y'])
# XT = [seq_cat(t) for t in train_prots]
# train_prots=np.asarray(XT)





