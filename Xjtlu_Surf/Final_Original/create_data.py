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

def seq_cat(prot):
    x = np.zeros(max_seq_len)
    for i, ch in enumerate(prot[:max_seq_len]): 
        x[i] = seq_dict[ch]
    return x 

seq_voc = "ABCDEFGHIKLMNOPQRSTUVWXYZ"
seq_dict = {v:(i+1) for i,v in enumerate(seq_voc)}
seq_dict_len = len(seq_dict)
max_seq_len = 1000
'''
此处的max_seq_len因为我电脑性能的问题改成了100，可以根据需要改为1000，后续网络模型相应的更改会标注在GCN模型上作为代表，其余网络模型只需要按照GCN一次更改即可
'''

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
