import pandas as pd
import numpy as np
import math
from collections import Counter
nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
    "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
    "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]

# nm=['ACC']
# for t in nm:
#     # df=pd.read_csv("Data_Pre_Primary/" + t + ".RSEM_genes_normalized__data_Level_3", index_col=[0], sep='\t')
#     df=pd.read_csv("/data/haochun/fyp/Data_Pre_Primary/" + t + ".RSEM_genes_normalized__data_Level_3", index_col=[0], sep='\t')
#     df.drop(labels="Label", axis=0, inplace=True)
#     df = df.astype(float)
#     round=True
#     for s in df.columns:
#         for j in df.index:
#             if df.loc[j,s] != 0 and round==True :
#                 mini = df.loc[j, s]
#                 round=False
#             elif df.loc[j,s] !=0 and round ==False and df.loc[j,s]< mini:
#                 mini=df.loc[j,s]
#     df.replace(to_replace=0, value=mini, inplace=True)
#     # for s in df.columns:
#     #     for j in df.index:
#     #         if df.loc[j,s] == 0:
#     #             mini=df.loc[j,s]
#
#
#     # df.loc["Label"] = t
#     # outputpath = 'test'+t+"test"
#     # df.to_csv(outputpath, sep='\t', index=True, header=True)
#
#     df.loc["Label"] = t
#     outputpath = '/data/haochun/fyp/Data_None_Zero/'+t+'.RSEM_genes_normalized__data_Level_3'
#     df.to_csv(outputpath, sep='\t', index=True, header=True)

# round=True
# for t in nm:
#     df=pd.read_csv("H:/FYP/Data_Pre_Primary/" + t + ".RSEM_genes_normalized__data_Level_3", index_col=[0], sep='\t')
#     if round==True:
#         sz=df.columns
#         round=False
#     elif round==False:
#         sz=sz.append(df.columns)
#
# dict1=dict(Counter(sz))
#
# dict2={key:value for key, value in dict1.items() if value>1}
#
# print(dict2)

