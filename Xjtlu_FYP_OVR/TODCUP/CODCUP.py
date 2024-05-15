import pandas as pd
import numpy as np
import math


#TOD_CUP code ,TCGA_matrix是已经合并的所有TCGA数据20532个gene。9158个样本
# df1=pd.read_csv("H:/FYP/summer/Data/wts/TCGA_matrix",sep='\t',index_col=[0])
#
# dic1={}
# R=True
# for i in df1.index:
#     if i == "ZNRF1"  and R ==True:
#         R=False
#         for j in df1.columns:
#             df1.loc[i, j] = float(df1.loc[i, j])
#     elif R ==False and i !="Label":
#         for j in df1.columns:
#             df1.loc[i, j] = float(df1.loc[i, j])
#
# for k in df1.index:
#     if k != "Label" :
#         print(k)
#         point = np.median(df1.loc[k, :], axis=0)
#         lst=[]
#         for s in df1.columns:
#             lst.append(abs(df1.loc[k,s]-point))
#
#     dic1[k]=np.median(lst)
#
# dic2=sorted(dic1.items(),key=lambda item: item[1],reverse=True)
#
# lst2=[]
# for key in range(0,5000):
#     lst2.append(dic2[key][0])
#
# df3=pd.DataFrame(df1,index=lst2)
#
# lst2.append('Label')
#
# df2=pd.DataFrame(df1,index=lst2)
#
# outputpath='H:/FYP/summer/Data/TOD_CUP/MAD_5000_Label'
# df2.to_csv(outputpath,sep='\t',index=True,header=True)
#
# outputpath1='H:/FYP/summer/Data/TOD_CUP/MAD_5000'
# df3.to_csv(outputpath1,sep='\t',index=True,header=True)



'''
wts2样本筛选5000个 TOD_CUP的基因
'''

# df1=pd.read_csv('H:/FYP/summer/Data/wts2/Sample_matrix_fill_log_max',sep='\t',index_col=[0])
#
# df2=pd.read_csv('H:/FYP/summer/Data/TOD_CUP/MAD_5000',sep='\t',index_col=[0])
#
#
# sz1=df2.index
#
# PD1=pd.DataFrame(df1,index=sz1)

# PD2=PD1
#
# for k in PD2.columns:
#     PD2.loc['Label',k]="BRCA"
#
#
#
# outputpath='H:/FYP/summer/Data/TOD_CUP/MAD_5000_wts2_Label'
# PD2.to_csv(outputpath,sep='\t',index=True,header=True)

# outputpath1='H:/FYP/summer/Data/TOD_CUP/MAD_5000_wts2'
# PD1.to_csv(outputpath1,sep='\t',index=True,header=True)


"""
TCGA Metastatic Data
"""


# df=pd.read_csv('H:/FYP/summer/Data/GEO/THCA/GSE202413_RNA_FPKM_Clean.txt',index_col=[0],sep='\t')
# for s in df.index:
#     if df[df.index==s].shape[0] !=1:
#         df4=df[df.index==s]
#         df.drop(index=s,inplace=True)
#         for j in df.columns:
#             df.loc[s,j]=df4[j].mean()
# df2=pd.read_csv("H:/FYP/summer/Data/TOD_CUP/MAD_5000",index_col=[0],sep='\t')
# df3=pd.read_csv("H:/FYP/summer/Data/wts/TCGA_matrix",index_col=[0],sep='\t')
# # df2.drop(index='Label',inplace=True)
# k1=df.index
# sz=[]
# for k in df2.index:
#     name=k.split('|')
#     sz.append(name[0])
# sz1=[]
# sz2=[]
# for t in sz:
#     if t in k1:
#         sz1.append(t)
#     else:
#         sz2.append(t)
#
# bg1=pd.DataFrame(df,index=sz1)
#
# for s in sz2:
#     sum=0
#     for t in df3.loc[s,:]:
#         sum=sum+float(t)
#         if sum>= 9158:
#             bg1.loc[s] = math.log2(sum / 9158)
#         else:
#             bg1.loc[s] = 0
#
#
#
# bg2=pd.DataFrame(bg1,index=sz)
# outputpath="H:/FYP/summer/Data/TOD_CUP/GEO_THCA"
# bg2.to_csv(outputpath, sep='\t', index=True, header=True)








