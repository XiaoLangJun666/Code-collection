#Lung Cancer SRP058626
import pandas as pd
import math
import pickle
import sklearn
print(sklearn.__version__)


# df=pd.read_csv("H:/FYP/New_Data/Out_Data/Lung_FPKM.txt",sep='\t',index_col=[0])
# for k in df.columns:
#     a=k.split('_')
#     if a[1] !='HCC' :
#         df.drop(labels=k,axis=1,inplace=True)
# df1=pd.read_csv("H:/FYP/New_Data/Chayi_gene/Gene_matrix_150",sep='\t',index_col=[0])
# df1.drop(labels='Label',axis=0,inplace=True)
# for s in df1.index:
#     b=s.split('|')
#     df1.rename(index={s: b[0]}, inplace=True)
#
# c=df1.index
# df2=pd.DataFrame(df,index=c)
# df2.fillna(value=1)
# for i in range(0,df2.shape[0]):
#     for j in range(0,df2.shape[1]):
#         if float(df2.iloc[i,j])==0:
#             df2.iloc[i,j]=1
#         df2.iloc[i,j]=math.log2(float(df2.iloc[i,j]))
#
# outputpath1='H:/FYP/New_Data/Out_Data/Lung'
# df2.to_csv(outputpath1,sep='\t',index=True,header=True)
# df2=pd.read_csv("H:/FYP/New_Data/Out_Data/Lung",sep='\t',index_col=[0])
# df2=df2.T
# loaded_model=pickle.load(open("H:/FYP/New_Data/Chayi_gene/Gene_model_150.dat","rb"))
# pre_y=loaded_model.predict(df2)
# print(pre_y)


#Breast  SRP043470
# df=pd.read_csv("H:/FYP/New_Data/Out_Data/Breast/Breast_Counts.txt",sep='\t',index_col=[0])
#
# df1=pd.read_csv("H:/FYP/New_Data/Chayi_gene/Gene_matrix_150",sep='\t',index_col=[0])
# df1.drop(labels='Label',axis=0,inplace=True)
# for s in df1.index:
#     b=s.split('|')
#     df1.rename(index={s: b[0]}, inplace=True)
# c=df1.index
# df=pd.DataFrame(df,index=c)
#
# for i in range(0,df.shape[0]):
#     for j in range(0,df.shape[1]):
#         if float(df.iloc[i,j])==0:
#             df.iloc[i,j]=1
#         df.iloc[i,j]=math.log2(float(df.iloc[i,j]))
#
#
#
# outputpath='H:/FYP/New_Data/Out_Data/Breast/Breast_matrix'
# df.to_csv(outputpath,sep='\t',index=True,header=True)
# df=pd.read_csv("H:/FYP/New_Data/Out_Data/Breast/Breast_matrix",sep='\t',index_col=[0])
# df=df.T
# loaded_model=pickle.load(open("H:/FYP/New_Data/Chayi_gene/Gene_model_150.dat","rb"))
# pre_y=loaded_model.predict(df)
# print(pre_y)


#Kid
# df=pd.read_csv("H:/FYP/New_Data/Out_Data/KID/Kid_matrix.txt",sep='\t',index_col=[0])
#
# df1=pd.read_csv("H:/FYP/New_Data/Chayi_gene/Gene_matrix_150",sep='\t',index_col=[0])
# df1.drop(labels='Label',axis=0,inplace=True)
# for s in df1.index:
#     b=s.split('|')
#     df1.rename(index={s: b[0]}, inplace=True)
# c=df1.index
#
# df=df[~df.index.duplicated()]
# df2=pd.DataFrame(df,index=c)
# for i in range(0,df2.shape[0]):
#     for j in range(0,df2.shape[1]):
#         if float(df2.iloc[i,j])==0:
#             df2.iloc[i,j]=1
#         df2.iloc[i,j]=math.log2(float(df2.iloc[i,j]))
# outputpath='H:/FYP/New_Data/Out_Data/KID/KID_matrix'
# df2.to_csv(outputpath,sep='\t',index=True,header=True)

# df=pd.read_csv("H:/FYP/New_Data/Out_Data/KID/KID_matrix",sep='\t',index_col=[0])
# df=df.T
# loaded_model=pickle.load(open("H:/FYP/New_Data/Chayi_gene/Gene_model_150.dat","rb"))
# pre_y=loaded_model.predict(df)
# print(pre_y)


#LUNG_Brest  GSE215384
# df=pd.read_csv("H:/FYP/New_Data/Out_Data/Breast_Liver/Liver_Breast.txt",sep='\t',index_col=[0])
# df=df[~df.index.duplicated()]
# df1=pd.read_csv("H:/FYP/New_Data/model/TSP_UNION_JOIN/Union_Join_TSP_JOIN_matrix",sep='\t',index_col=[0])
# df1.drop(labels='Label',axis=0,inplace=True)
# for s in df1.index:
#     b=s.split('|')
#     df1.rename(index={s: b[0]}, inplace=True)
# c=df1.index
# df2=pd.DataFrame(df,index=c)
# for i in range(0,df2.shape[0]):
#     for j in range(0,df2.shape[1]):
#         if float(df2.iloc[i,j])==0:
#             df2.iloc[i,j]=1
#         df2.iloc[i,j]=math.log2(float(df2.iloc[i,j]))
# outputpath='H:/FYP/New_Data/Out_Data/Breast_Liver/matrix_TSP'
# df2.to_csv(outputpath,sep='\t',index=True,header=True)
# df=pd.read_csv("H:/FYP/New_Data/Out_Data/Breast_Liver/matrix_TSP",sep='\t',index_col=[0])
# df=df.T
# loaded_model=pickle.load(open("H:/FYP/New_Data/model/TSP_UNION_JOIN/UNION_JOIN_TSP_JOIN_model.dat","rb"))
# pre_y=loaded_model.predict(df)
# print(pre_y)



#BLCA GSE215947
# df=pd.read_csv("H:/FYP/New_Data/Out_Data/BLCA/GSE215947_FPKMs_allsamples.txt",sep='\t',index_col=[0])
# df=df[~df.index.duplicated()]
# df1=pd.read_csv("H:/FYP/New_Data/model/TSP_UNION_JOIN/Union_Join_TSP_JOIN_matrix",sep='\t',index_col=[0])
# df1.drop(labels='Label',axis=0,inplace=True)
# for s in df1.index:
#     b=s.split('|')
#     df1.rename(index={s: b[0]}, inplace=True)
# c=df1.index
# df2=pd.DataFrame(df,index=c)
# for i in range(0,df2.shape[0]):
#     for j in range(0,df2.shape[1]):
#         if float(df2.iloc[i,j])==0:
#             df2.iloc[i,j]=1
#         df2.iloc[i,j]=math.log2(float(df2.iloc[i,j]))
# outputpath='H:/FYP/New_Data/Out_Data/BLCA/matrix_TSP'
# df2.to_csv(outputpath,sep='\t',index=True,header=True)

# df=pd.read_csv("H:/FYP/New_Data/Out_Data/BLCA/matrix_TSP",sep='\t',index_col=[0])
# df=df.T
# loaded_model=pickle.load(open("H:/FYP/New_Data/model/TSP_UNION_JOIN/UNION_JOIN_TSP_JOIN_model.dat","rb"))
# pre_y=loaded_model.predict(df)
# print(pre_y)

###GSE164862
# df=pd.read_csv("H:/FYP/New_Data/Out_Data/BLCA/GSE164862_Experiment_processed_data.txt",sep='\t',index_col=[0])
# df=df[~df.index.duplicated()]
# df1=pd.read_csv("H:/FYP/New_Data/model/TSP_UNION_JOIN/Union_Join_TSP_JOIN_matrix",sep='\t',index_col=[0])
# df1.drop(labels='Label',axis=0,inplace=True)
# for s in df1.index:
#     b=s.split('|')
#     df1.rename(index={s: b[0]}, inplace=True)
# c=df1.index
# df2=pd.DataFrame(df,index=c)
# for i in range(0,df2.shape[0]):
#     for j in range(0,df2.shape[1]):
#         if float(df2.iloc[i,j])==0:
#             df2.iloc[i,j]=1
#         df2.iloc[i,j]=math.log2(float(df2.iloc[i,j]))
# outputpath='H:/FYP/New_Data/Out_Data/BLCA/2_matrix_TSP'
# df2.to_csv(outputpath,sep='\t',index=True,header=True)

# df=pd.read_csv("H:/FYP/New_Data/Out_Data/BLCA/2_matrix_TSP",sep='\t',index_col=[0])
# df=df.T
# loaded_model=pickle.load(open("H:/FYP/New_Data/model/TSP_UNION_JOIN/UNION_JOIN_TSP_JOIN_model.dat","rb"))
# pre_y=loaded_model.predict(df)
# print(pre_y)


#COAD

# df=pd.read_csv("H:/FYP/New_Data/Out_Data/COAD/COAD.txt",sep='\t',index_col=[0])
# df=df[~df.index.duplicated()]
# df1=pd.read_csv("H:/FYP/New_Data/model/TSP_UNION_JOIN/Union_Join_TSP_JOIN_matrix",sep='\t',index_col=[0])
# df1.drop(labels='Label',axis=0,inplace=True)
# for s in df1.index:
#     b=s.split('|')
#     df1.rename(index={s: b[0]}, inplace=True)
# c=df1.index
# df2=pd.DataFrame(df,index=c)
# for i in range(0,df2.shape[0]):
#     for j in range(0,df2.shape[1]):
#         if float(df2.iloc[i,j])==0:
#             df2.iloc[i,j]=1
#         df2.iloc[i,j]=math.log2(float(df2.iloc[i,j]))
# outputpath='H:/FYP/New_Data/Out_Data/COAD/COAD_matrix_TSP'
# df2.to_csv(outputpath,sep='\t',index=True,header=True)

# df=pd.read_csv("H:/FYP/New_Data/Out_Data/COAD/COAD_matrix_TSP",sep='\t',index_col=[0])
# df=df.T
# loaded_model=pickle.load(open("H:/FYP/New_Data/model/TSP_UNION_JOIN/UNION_JOIN_TSP_JOIN_model.dat","rb"))
# pre_y=loaded_model.predict(df)
# print(pre_y)


#HNSC
# df=pd.read_csv("H:/FYP/New_Data/Out_Data/HNSC/GSE207182_lcpm_rnaSeq_pdxHNtumor.txt",sep='\t',index_col=[0])
# df=df[~df.index.duplicated()]
# df1=pd.read_csv("H:/FYP/New_Data/Chayi_gene/Gene_matrix_200",sep='\t',index_col=[0])
# df1.drop(labels='Label',axis=0,inplace=True)
# for s in df1.index:
#     b=s.split('|')
#     df1.rename(index={s: b[0]}, inplace=True)
# c=df1.index
# df2=pd.DataFrame(df,index=c)
# for i in range(0,df2.shape[0]):
#     for j in range(0,df2.shape[1]):
#         if float(df2.iloc[i,j])==0:
#             df2.iloc[i,j]=1
#         df2.iloc[i,j]=math.log2(float(df2.iloc[i,j]))
# outputpath='H:/FYP/New_Data/Out_Data/HNSC/HNSC_matrix'
# df2.to_csv(outputpath,sep='\t',index=True,header=True)

# df=pd.read_csv("H:/FYP/New_Data/Out_Data/HNSC/HNSC_matrix",sep='\t',index_col=[0])
# for i in range(0,df.shape[0]):
#     for j in range(0,df.shape[1]):
#         if float(df.iloc[i,j])==0:
#             df.iloc[i,j]=1
#         df.iloc[i,j]=math.log2(abs(float(df.iloc[i,j])))
# df=df.T
# loaded_model=pickle.load(open("H:/FYP/New_Data/Chayi_gene/Gene_model_200.dat","rb"))
# pre_y=loaded_model.predict(df)
# print(pre_y)

#Breast
df=pd.read_csv("H:/FYP/New_Data/Out_Data/HNSC/GSE207182_lcpm_rnaSeq_pdxHNtumor.txt",sep='\t',index_col=[0])
df=df[~df.index.duplicated()]
df1=pd.read_csv("H:/FYP/New_Data/Chayi_gene/Gene_matrix_100",sep='\t',index_col=[0])
df1.drop(labels='Label',axis=0,inplace=True)
# for k in df.index:
#     a=k.split('|')
#     df.rename(index={k: a[1]}, inplace=True)
# df=df[~df.index.duplicated()]
for s in df1.index:
    b=s.split('|')
    df1.rename(index={s: b[0]}, inplace=True)
c=df1.index
df2=pd.DataFrame(df,index=c)
# for i in range(0,df2.shape[0]):
#     for j in range(0,df2.shape[1]):
#         if float(df2.iloc[i,j])==0:
#             df2.iloc[i,j]=1
#         df2.iloc[i,j]=math.log2(float(df2.iloc[i,j]))
outputpath='H:/FYP/New_Data/Out_Data/HNSC/matrix_100'
df2.to_csv(outputpath,sep='\t',index=True,header=True)
# df1.drop(labels='Label',axis=0,inplace=True)
# for k in df.index:
#     a = k.split('|')
#     for s in df1.index:
#         b = s.split('|')
#         if a[1]==b[0]:
#             df.rename(index={k:s},inplace=True)
#
# df=df[~df.index.duplicated()]
#
# c=df1.index
# df2=pd.DataFrame(df,index=c)
# for i in range(0,df2.shape[0]):
#     for j in range(0,df2.shape[1]):
#         if float(df2.iloc[i,j])==0:
#             df2.iloc[i,j]=1
#         df2.iloc[i,j]=math.log2(float(df2.iloc[i,j]))
# outputpath='H:/FYP/New_Data/Out_Data/KID/KID_Matrix_1'
# df2.to_csv(outputpath,sep='\t',index=True,header=True)

df=pd.read_csv("H:/FYP/New_Data/Out_Data/HNSC/matrix_100_1.txt",sep='\t',index_col=[0])

df=df.T
loaded_model=pickle.load(open("H:/FYP/New_Data/Chayi_gene/Gene_model_100.dat","rb"))
pre_y=loaded_model.predict(df)
print(pre_y)