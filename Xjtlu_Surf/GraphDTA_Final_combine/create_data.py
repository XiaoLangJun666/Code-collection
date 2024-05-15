import pandas as pd
import numpy as np
import os
import json,pickle
from collections import OrderedDict
from rdkit import Chem
from rdkit.Chem import MolFromSmiles
import networkx as nx
from utils import *

def atom_features(atom):
    return np.array(one_of_k_encoding_unk(atom.GetSymbol(),['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na','Ca', 'Fe', 'As', 'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb','Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti', 'Zn', 'H','Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr','Cr', 'Pt', 'Hg', 'Pb', 'Unknown']) +
                    one_of_k_encoding(atom.GetDegree(), [0, 1, 2, 3, 4, 5, 6,7,8,9,10]) +
                    one_of_k_encoding_unk(atom.GetTotalNumHs(), [0, 1, 2, 3, 4, 5, 6,7,8,9,10]) +
                    one_of_k_encoding_unk(atom.GetImplicitValence(), [0, 1, 2, 3, 4, 5, 6,7,8,9,10]) +
                    [atom.GetIsAromatic()])

def one_of_k_encoding(x, allowable_set):
    if x not in allowable_set:
        raise Exception("input {0} not in allowable set{1}:".format(x, allowable_set))
    return list(map(lambda s: x == s, allowable_set))

def one_of_k_encoding_unk(x, allowable_set):
    """Maps inputs not in the allowable set to the last element."""
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))

def smile_to_graph(smile):
    mol = Chem.MolFromSmiles(smile)
    
    c_size = mol.GetNumAtoms()
    
    features = []
    for atom in mol.GetAtoms():
        feature = atom_features(atom)
        features.append( feature / sum(feature) )

    edges = []
    for bond in mol.GetBonds():
        edges.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
    g = nx.Graph(edges).to_directed()
    edge_index = []
    for e1, e2 in g.edges:
        edge_index.append([e1, e2])
        
    return c_size, features, edge_index


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
"""
以上这部分建立一个字典，字典的key是蛋白质标签，内容是其氨基酸序列
"""
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
"""
以上这部分是pfam的基础读取的代码
"""


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



    dc[t] = np.concatenate((matrix1,matrix2,matrix3,matrix5), 1)

"""
建立一个新的字典，字典的key是蛋白质的标签，以上的matrix1，2，3，4，5分别代表 blo62，pssm，psi-pred，asaquick和pfam.可以按照需要在在最后一句代码中更改所需要的特征组合
"""


def seq_cat(prot):
    matrix = []
    for t in dt.keys():
        if prot == dt[t]:
            if dc[t].shape[0] >750:
                matrix = dc[t][:750, :]
            else:
                matrix = np.zeros((750, dc[t].shape[1]))
                matrix[:dc[t].shape[0], :] = dc[t]

    return matrix




compound_iso_smiles = []
for dt_name in ['kinase']:
    opts = ['train', 'valid', 'test']
    for opt in opts:
        df = pd.read_csv('data/' + dt_name + '/kinase_grayscale_cold_protein_balance_' + opt + '.csv')
        compound_iso_smiles += list( df['canonical_smiles.y'] )
compound_iso_smiles = set(compound_iso_smiles)
smile_graph = {}
for smile in compound_iso_smiles:
    g = smile_to_graph(smile)
    smile_graph[smile] = g

datasets = ['kinase']
# convert to PyTorch data format
for dataset in datasets:
    processed_data_file_train = 'data/processed/' + dataset + '_grayscale_cold_protein_balance_train.pt'
    processed_data_file_valid = 'data/processed/' + dataset + '_grayscale_cold_protein_balance_valid.pt'
    processed_data_file_test = 'data/processed/' + dataset + '_grayscale_cold_protein_balance_test.pt'
    if ((not os.path.isfile(processed_data_file_train)) or (not os.path.isfile(processed_data_file_valid)) or (not os.path.isfile(processed_data_file_test))):
    #if (not os.path.isfile(processed_data_file_test)):   
        df = pd.read_csv('data/' + dataset + '/kinase_grayscale_cold_protein_balance_train.csv')
        train_drugs, train_prots,  train_Y = list(df['canonical_smiles.y']),list(df['sequence.y']),list(df['label.y'])
        XT = [seq_cat(t) for t in train_prots]
        train_drugs, train_prots, train_Y = np.asarray(train_drugs), np.asarray(XT), np.asarray(train_Y)
       	df = pd.read_csv('data/' + dataset + '/kinase_grayscale_cold_protein_balance_valid.csv')
        valid_drugs, valid_prots, valid_Y = list(df['canonical_smiles.y']),list(df['sequence.y']),list(df['label.y'])
        XT = [seq_cat(t) for t in valid_prots]
        valid_drugs, valid_prots,  valid_Y = np.asarray(valid_drugs), np.asarray(XT), np.asarray(valid_Y)
        df = pd.read_csv('data/' + dataset + '/kinase_grayscale_cold_protein_balance_test.csv')
        test_drugs, test_prots,  test_Y = list(df['canonical_smiles.y']),list(df['sequence.y']),list(df['label.y'])
        XT = [seq_cat(t) for t in test_prots]
        test_drugs, test_prots,  test_Y = np.asarray(test_drugs), np.asarray(XT), np.asarray(test_Y)

        # make data PyTorch Geometric ready
        print('preparing ', dataset + '_grayscale_cold_protein_balance_train.pt in pytorch format!')
        train_data = TestbedDataset(root='data/processed', dataset=dataset+'_grayscale_cold_protein_balance_train', xd=train_drugs, xt=train_prots, y=train_Y,smile_graph=smile_graph)
        print('preparing ', dataset + '_grayscale_cold_protein_balance_valid.pt in pytorch format!')
        valid_data = TestbedDataset(root='data/processed', dataset=dataset+'_grayscale_cold_protein_balance_valid', xd=valid_drugs, xt=valid_prots, y=valid_Y,smile_graph=smile_graph)
        print('preparing ', dataset + '_grayscale_cold_protein_balance_test.pt in pytorch format!')
        test_data = TestbedDataset(root='data/processed', dataset=dataset+'_grayscale_cold_protein_balance_test', xd=test_drugs, xt=test_prots, y=test_Y,smile_graph=smile_graph)
        print(processed_data_file_train, ' and ', processed_data_file_valid, 'and', processed_data_file_test, ' have been created')
        #print(processed_data_file_test, ' have been created')        
    else:
        print(processed_data_file_train, ' and ', processed_data_file_valid, 'and', processed_data_file_test, ' are already created')
        print(processed_data_file_test, ' have been created')  
