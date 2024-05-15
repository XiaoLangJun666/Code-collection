import pandas as pd
import pandas as pd
import pickle
import numpy as np
import math
import  seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict
from sklearn.metrics import accuracy_score, confusion_matrix
#将样本合并，生成矩阵
# nm=['2200136','2200454','2200640','2200682','2200777','2200778','2200874','2201042','2300110','2300675']
# round=True
# for n in nm:
#     df=pd.read_csv('H:/FYP/summer/Data/wts2/20230801-SZ'+n+'.genes.results', index_col=[0], sep='\t')
#     df=df[['FPKM']]
#     df=df.iloc[:-4,:]
#     name=n+"_FPKM"
#     df.rename(columns={"FPKM":name}, inplace=True)
#     if round==True:
#         DF=df
#         round=False
#     else:
#         DF=pd.concat([DF,df],axis=1)
# print(DF.index)
# df2=pd.read_csv('H:/FYP/summer/Data/wts2/Gene.table',header=None,index_col=[0],sep='\t')
# # df2=df2.iloc[:-4,:]
# # df2.rename(columns={1:"gene_name"},inplace=True)
#
# for k in DF.index:
#     re=df2.loc[k,1]
#     DF.rename(index={k:re},inplace=True)
#
#
#
# outputpath="H:/FYP/summer/Data/wts2/Sample_matrix"
# DF.to_csv(outputpath, sep='\t', index=True, header=True)




#去重复（选取表达量最大的）
# df1=pd.read_csv('H:/FYP/summer/Data/wts2/Sample_matrix',sep='\t',index_col=[0])
# print(df1.shape)
#
# for k in df1.index:
#     if df1[df1.index==k].shape !=(1,10):
#         print(k)
#         df2=df1[df1.index==k]
#         dic1={}
#         for s in range(0,len(df2.index)):
#             dic1[s]=sum(df2.iloc[s,:])
#         Largest= max(dic1, key=lambda i: dic1[i])
#         df1.drop(index=k,inplace=True)
#         df1.loc[k,:]=df2.iloc[Largest,:]
# outputpath="H:/FYP/summer/Data/wts2/Sample_matrix_max"
# df1.to_csv(outputpath, sep='\t', index=True, header=True)

#去重复（选取平均值）
# df1=pd.read_csv('H:/FYP/summer/Data/wts2/Sample_matrix',sep='\t',index_col=[0])
# print(df1.shape)
#
# for k in df1.index:
#     if df1[df1.index==k].shape !=(1,10):
#         df2=df1[df1.index==k]
#         df1.drop(index=k,inplace=True)
#         for s in df1.columns:
#             df1.loc[k,s]=df2[s].mean()
#
#
# outputpath="H:/FYP/summer/Data/wts2/Sample_matrix_average"
# df1.to_csv(outputpath, sep='\t', index=True, header=True)
#去重复（和）
# df1=pd.read_csv('H:/FYP/summer/Data/wts2/Sample_matrix',sep='\t',index_col=[0])
# print(df1.shape)
#
# for k in df1.index:
#     if df1[df1.index==k].shape !=(1,10):
#         df2=df1[df1.index==k]
#         df1.drop(index=k,inplace=True)
#         for s in df1.columns:
#             df1.loc[k,s]=df2[s].sum()
#
# print(df1.shape)
# outputpath="H:/FYP/summer/Data/wts2/Sample_matrix_Sum"
# df1.to_csv(outputpath, sep='\t', index=True, header=True)




#去重复（第一个）
# df1=pd.read_csv('H:/FYP/summer/Data/wts/Sampe_matrix',sep='\t',index_col=[0])
# print(df1.shape)
# sz=[]
# for k in df1.index:
#     if k not in sz:
#         df1.loc[k,"Bo"]=1
#         sz.append(k)
#     else:
#         df1.loc[k,"Bo"]=0
# outputpath="H:/FYP/summer/Data/wts/Sample_matrix_boolean"
# df1.to_csv(outputpath, sep='\t', index=True, header=True)




#合并TCGA的数据
import pandas as pd
import numpy as np
import math
# nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
#     "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
#     "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]
#
# round=True
# for s in nm:
#     df1=pd.read_csv('H:/FYP/Data_Pre/'+s+'.RSEM_genes_normalized__data_Level_3',sep='\t',index_col=[0])
#     if round==True:
#         DF = df1
#         round=False
#     else:
#         DF=pd.concat([DF,df1],axis=1)
#
# for k in DF.index:
#     name=k.split('|')[0]
#     DF.rename(index={k:name},inplace=True)
#
# outputpath="H:/FYP/summer/Data/wts/TCGA_matrix"
# DF.to_csv(outputpath, sep='\t', index=True, header=True)



#将样本的基因与TCGA的基因overlap，并填补空缺的基因（使用TCGA中所有样本的均值）

# df1=pd.read_csv('H:/FYP/summer/Data/wts2/Sample_matrix_average',sep='\t',index_col=[0])
# df2=pd.read_csv('H:/FYP/summer/Data/wts/TCGA_matrix',sep='\t',index_col=[0])
# df2.drop(index='Label',inplace=True)
#
# #统计个数
# # sz1=df1.index
# # sz2=df2.index
# #
# # sum=0
# # for k in sz2:
# #     if k in sz1:
# #         sum+=1
# # print(sum,sum/len(sz1),sum/len(sz2))
# #
# DF=pd.DataFrame(columns=['2200136_FPKM','2200454_FPKM','2200640_FPKM','2200682_FPKM','2200777_FPKM','2200778_FPKM','2200874_FPKM','2201042_FPKM','2300110_FPKM','2300675_FPKM'])
# sz1=df1.index
# sz2=df2.index
# for k in sz2:
#     if k in sz1:
#         DF.loc[k,:]=df1.loc[k,:]
#     else:
#         sum=0
#         for t in df2.loc[k,:]:
#             sum=sum+float(t)
#         if sum >=9158:
#             DF.loc[k]=math.log2(sum/9158)
#         else:
#             DF.loc[k]=0
#
# DF=pd.DataFrame(DF,index=sz2)
#
# outputpath="H:/FYP/summer/Data/wts2/Sample_matrix_fill_log_average"
# DF.to_csv(outputpath, sep='\t', index=True, header=True)










#生成对应的测试样本矩阵
# nm=['50','100','150','200','250','300']
# for s in nm:
#     df1=pd.read_csv('H:/FYP/summer/Data/wts2/Sample_matrix_fill_log_sum',sep='\t',index_col=[0])
#     df2=pd.read_csv("H:/FYP/summer/Data/Chayi/Gene_matrix_"+s,sep='\t',index_col=[0])
#     df2.drop(index="Label",inplace=True)
#     sz=[]
#     for k in df2.index:
#         name=k.split('|')
#         sz.append(name[0])
#     df3=pd.DataFrame(df1,index=sz)
#
#     outputpath='H:/FYP/summer/Data/wts2/log_clinical/sum/Chayi_'+s
#     df3.to_csv(outputpath,sep='\t',index=True,header=True)



#样本结果预测
# time=['50','100','150','200','250','300']
# for p in time:
#     df1=pd.read_csv('H:/FYP/summer/Data/wts2/log_clinical/sum/Chayi_'+p,sep='\t',index_col=[0])
#     for i in range(0,df1.shape[0]):
#         for j in range(0,df1.shape[1]):
#             if float(df1.iloc[i,j])==0:
#                 df1.iloc[i,j]=1
#             df1.iloc[i,j]=math.log2(abs(float(df1.iloc[i,j])))
#     df1=df1.T
#     loaded_model=pickle.load(open("H:/FYP/summer/Model/Chayi_2_asag/Gene_Model_"+p+'.dat',"rb"))
#     pre_y=loaded_model.predict(df1)
#     print(pre_y)


#寻找频次出现最多的方法
df1=pd.read_csv('H:/FYP/summer/Data/wts2/Result/clinical_result_2_logfill.txt',sep='\t',index_col=[0])
#统计每个样本出现的癌症预测的类别的分布
yb={}
for k in range(5,15):
    single={}
    for s in range(len(df1.index)):
        pre=df1.iloc[s,k]
        if pre not in single.keys():
            single[pre]=1
        else:
            single[pre]=single[pre]+1

    name=df1.columns[k]

    single=sorted(single.items(),key=lambda item:item[1],reverse=True)
    yb[name]=single
print(yb)
#统计表现最好的方法
# std=['STAD','SARC','OV','STAD','KIRC','KIRC','STAD','BLCA','MESO','KIRC']
#
# for s in range(len(df1.index)):
#     sum = 0
#     for k in range(4, 14):
#         pre=df1.iloc[s,k]
#         if k==4:
#             if pre=="STAD" or pre=="KIRC":
#                 sum+=1
#         elif k==6:
#             if pre=="OV" or pre =="BRCA":
#                 sum+=1
#         else:
#             if pre==std[k-4]:
#                 sum+=1
#     index=df1.index[s]
#     df1.loc[index,'accuracy']=sum/10
#
#
# outputpath='H:/FYP/summer/clinical_result_Fill_accuracy'
# df1.to_csv(outputpath,sep='\t',index=True,header=True)


