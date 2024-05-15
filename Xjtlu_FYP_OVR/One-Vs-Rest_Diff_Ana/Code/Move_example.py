import pandas as pd
import pickle
import numpy as np
import math
import  seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict
from sklearn.metrics import accuracy_score, confusion_matrix
nm=["BRCA","CESC","ESCA","HNSC","PAAD","PCPG","PRAD","SARC","SKCM","THCA"]
#对转移癌选择我们模型中使用的特征基因
df1=pd.read_csv('H:/FYP/summer/Data/TOD_CUP/MAD_5000',sep='\t',index_col=[0])
sz=df1.index
round=True
for t in nm:
    df2=pd.read_csv('H:/FYP/Data_Pre_move/' + t + ".RSEM_genes_normalized__data_Level_3",sep='\t',index_col=[0])
    for k in df2.index:
        if k !='Label':
            name=k.split("|")[0]
            df2.rename({k:name},inplace=True)
    if round==True:
        pd1=pd.DataFrame(df2,index=sz)
        round=False
    elif round==False:
        pd2=pd.DataFrame(df2,index=sz)
        pd1=pd.concat([pd1,pd2],axis=1)

# for n in pd1.index:
#     re=n.replace('-','.')
#     pd1.rename(index={n:re},inplace=True)

outputpath1='H:/FYP/summer/Data/TOD_CUP/MAD_5000_TSP_Metastatic'
pd1.to_csv(outputpath1,sep='\t',index=True,header=True)

#预测模型在转移癌上的预测能力
# number=['50','100','150','200','250','300']
#
# tst=pd.read_csv("H:/FYP/summer/Data/Meta_TSP/Gene_matrix_Join",index_col=[0],sep='\t')
# tst=tst.T
# tst_x=tst.drop(["Label"],axis=1)
# df2=np.zeros((tst_x.shape[0],tst_x.shape[1]),dtype=np.float32)
# for a in range(0,tst_x.shape[0]):
#     for b in range(0,tst_x.shape[1]):
#         if float(tst_x.iloc[a,b])<=1:
#             tst_x.iloc[a,b]=1
#         df2[a,b]=math.log2(float(tst_x.iloc[a,b]))
# tst_y=tst['Label']
#
# loaded_model=pickle.load(open("H:/FYP/summer/Model/TSP/Join_Saga.dat","rb"))
#
#
# # loaded_model.ptrdict(tst_x)
#
# tst_accuracy=accuracy_score(tst_y,loaded_model.predict(df2))
# print("Training accuracy: ", tst_accuracy)
# tst_train_pred = cross_val_predict(loaded_model, df2, tst_y, cv=10)
# tst_conf_matrix = confusion_matrix(tst_y, tst_train_pred)
# print("Training confusion matrix:\n", tst_conf_matrix)
#
# y_pred=loaded_model.predict(df2)
#
#
#
# test_conf_matrix=confusion_matrix(tst_y,y_pred,labels=nm)
# fn=[]
# for i in range(len(nm)):
#     fn.append(sum(test_conf_matrix[i,:]) - test_conf_matrix[i,i])
#
# fnr = []
# for i in range(len(nm)):
#     fnr.append(fn[i] / sum(test_conf_matrix[i,:]))
# print(fnr)
# overall_fnr=sum(fn)/sum(sum(test_conf_matrix))
#
# print(overall_fnr)
# sns.heatmap(test_conf_matrix, annot=False, cmap='RdYlBu')
# plt.xlabel('Predicted Labels')
# plt.ylabel('True Labels')
# plt.title('Confusion Matrix')
# plt.show()