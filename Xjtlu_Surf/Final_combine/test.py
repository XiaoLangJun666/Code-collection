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

# da=pd.read_csv("data/kinase/kinase_grayscale_cold_protein_balance_train.csv",sep=',')
# da=da.values
# dt={}
# b=0
# for i in da[:,6]:
#      if i not in dt.keys():
#          dt[i]=da[:,2][b]
#      b=b+1
# df = pd.read_csv('data/' + "kinase" + '/kinase_grayscale_cold_protein_balance_test.csv')
# df = df.values
# c = 0
# for k in df[:, 2]:
#      if k not in dt.keys():
#          dt[k] = df[:, 3][c]
#      c = c + 1
#
#
# dc={}
#
# for t in dt.keys():
#
#
#
#           data1 = pd.read_csv('FinalBLOSUM/' + t + ".csv", sep=',')
#           data1 = data1.values
#           matrix1=data1[:len(dt[t]),1:]
#
#
#
#           data2 = pd.read_csv('pssm/' + t + '.fasta.pssm', sep='\t', header=None)
#           data2 = data2.values
#           matrix2 = np.zeros((data2.shape[0], 20), dtype="float64")
#           index = 0
#           for ds in data2:
#                 array = [float(x) for x in ds[1:] if x != '']
#                 array = np.array(array)
#                 matrix2[index, :] = array
#                 index = index + 1
#
#
#
#           data3 = pd.read_csv("psi-pred/" + t + ".ss2")
#           data3 = data3.values
#           matrix3 = np.zeros([data3.shape[0], 3])
#           x = 0
#           y = 0
#           for j in data3:
#             for i in range(0, len(str(j))):
#                 if str(j)[i] == ".":
#                     s = str(j)[i - 1:i + 4]
#                     matrix3[x, y] = s
#                     y = y + 1
#             y = 0
#             x = x + 1
#
#
#
#           data4 = pd.read_csv('ASAquick/asaq.' + t + ".fasta/asaq.pred", sep=',', header=None)
#           data4 = data4.values
#           matrix4 = np.zeros(([data4.shape[0] - 1, 1]))
#           for i in range(data4.shape[0] - 1):
#               matrix4[i] = float(str(data4[i]).split(' ')[2])
#
#           dc[t] = np.concatenate((matrix1, matrix2, matrix3,matrix4), 1)

#
#
# def seq_cat(prot):
#     matrix = []
#     for t in dt.keys():
#          if prot == dt[t] :
#              matrix = dc[t][0:20, :]
#     if matrix==[]:
#             print(prot)
#
#
#     return matrix
#
# #
# df = pd.read_csv('data/' + "kinase" + '/kinase_grayscale_cold_protein_balance_test.csv')
# train_prots=list(df['sequence.y'])
# XT = [seq_cat(t) for t in train_prots]
#
# print(dc["O00141"])

# train_prots=np.asarray(XT)
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
#设计一个包含所有蛋白质名字和序列的字典
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

#建立一个表格A代表domain名字，B代表Accession，C代表蛋白质编号，D和E分别是域的from和to
dr=pd.read_csv('final-Pfam/result.dom',sep=',')
dr=dr.values
dp=pd.DataFrame(np.arange(1770).reshape(354,5),columns=list('ABCDE'))

for i in range(2,356):
    a=dr[i]
    b = 0
    for x in str(a).split(' '):
        if x!='':
            if b==0:
               x=x[2:len(x)]
               dp["A"][i-2]=x
            elif b==1:
                dp["B"][i-2]=x
            elif b==3:
                dp["C"][i-2]=x
            elif b==19:
                dp["D"][i-2]=x
            elif b==20:
                dp["E"][i-2]=x
            b=b+1


#字典dc代表每一个蛋白质编码所包含的accession
dc={}
for t in dt.keys():
    for i in range(0,354):
        if dp["C"][i]==t and t not in dc.keys():
             dc[t]=dp["B"][i]
        elif dp["C"][i]==t:
             dc[t] = dc[t]+","+dp["B"][i]

#建立字典cc查看每一个accession出现的频率
cc={}
for a in dp["B"]:
    if a not in cc.keys():
        cc[a]=1
    else:
        cc[a]=cc[a]+1

dc = {}

for t in dt.keys():

    data1 = pd.read_csv('FinalBLOSUM/' + t + ".csv", sep=',')
    data1 = data1.values
    matrix1 = data1[:len(dt[t]), 1:]

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

    # 建立一个dataframe，列名为所有出现的accession，共19种，行数和蛋白质序列长度保持一致
    data5 = pd.DataFrame(np.zeros(len(dt[t]) * 19).reshape(len(dt[t]), 19), columns=list(cc.keys()))

    for i in range(0, 354):

        if dp["C"][i] == t:
            a = dp["B"][i]
            b = int(dp["D"][i])
            c = int(dp["E"][i])
            for j in range(b, c + 1):
                data5[a][j] = 1
    matrix5 = data5.to_numpy()

    dc[t] = np.concatenate((matrix2, matrix3, matrix4,matrix5), 1)


#建立一个dataframe，列名为所有出现的accession，共19种，行数和蛋白质序列长度保持一致
# for u in dt.keys():
#     if u=="O00311":
#       nm=pd.DataFrame(np.zeros(len(dt[u])*19).reshape(len(dt[u]),19),columns=list(cc.keys()))
#
#       for i in range(0,354):
#
#           if dp["C"][i]==u:
#               a=dp["B"][i]
#               b=int(dp["D"][i])
#               c=int(dp["E"][i])
#               for j in range(b,c+1):
#                    nm[a][j]=1
#       matrix4=nm.to_numpy()
#       for i in range(0,matrix4.shape[0]):
#           print(matrix4[i])










# dp=pd.DataFrame(np.arange(1770).reshape(354,5),columns=list('ABSCE'))
# print(dp["A"][0])
# print(dp["A"][1])
# print(dp["A"][2])
# a="['PK_Tyr_Ser-Thr"
# print(a[2:len(a)])

