import pandas as pd
from scipy.stats import pearsonr
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


#生成test文件
# df=pd.read_csv('FYP_R/TCGA_Test/Data_Test_Label', index_col=[0], sep='\t')
# DF1=df.iloc[0:20,0:5]
# DF2=df.iloc[0:20,10:15]
# DF = pd.concat([DF1, DF2], axis=1)
# outputpath='H:/FYP/PER_test'
# DF.to_csv(outputpath,sep='\t',index=True,header=True)


#
nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
    "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
    "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]

# nm=["ACC"]


#pearsonr 生成基因对的得分>0.9
# for t in nm:
#     # df=pd.read_csv("/data/haochun/fyp/Data_None_Zero/" + t + ".RSEM_genes_normalized__data_Level_3", index_col=[0], sep='\t')
#     df=pd.read_csv("/data/haochun/fyp/Data_Pre_Primary/" + t + ".RSEM_genes_normalized__data_Level_3", index_col=[0], sep='\t')
#     df.drop(labels="Label", axis=0, inplace=True)
#     df=df.T
#     df=df.astype(float)
#     a=df.corr(method="pearson")
#     dic={}
#     for n in a.columns:
#         for m in a.index:
#             if n!=m:
#                 val=a.loc[m,n]
#                 if abs(val)>=0.9:
#                     # if n not in ge:
#                     #     ge.append(str(n))
#                     # if m not in ge:
#                     #     ge.append(str(m))
#                     na=str(n)+","+str(m)
#                     dic[na]=abs(val)
#
#     dic=sorted(dic.items(),key=lambda item:item[1],reverse=True)
#     gene_df=pd.DataFrame(dic)
#     rownames=gene_df.iloc[:,0]
#     gene_df=gene_df.iloc[:,1]
#     gene_df.index=rownames
#     outputpath = '/data/haochun/fyp/New/Pearson/Pearson_'+t+'_gene'
#     # outputpath = '/data/haochun/fyp/Data_Zero_Pearson/Pearson_'+t+'_gene'
#     gene_df.to_csv(outputpath, sep='\t', index=True, header=True)

    # print(gene_df)
    # dic=dic[0:200]
    # print(dic[0:200])
    # print(len(dic[0:200]))
    # ge=[]
    # for k in range(0,len(dic)):
    #     name=dic[k][0].split(',')
    #     if name[0] not in ge:
    #         ge.append(name[0])
    #     if name[1] not in ge:
    #         ge.append(name[1])
    # print(ge)
    # print(len(ge))



#找pearson系数下的特意基因，第一种筛选基因的方法（选择一个癌种的基因与剩下所有的合并的基因的差集）

# new_data=pd.DataFrame(columns=["gene","length"])
# for t in nm:
#     print(t)
#     diff_gene=[]
#     nm2 = ["ACC","BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH",
#           "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG",
#           "PRAD", "READ", "SARC", "SKCM", "STAD",  "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
#     total_gene = []
#     for s in nm2:
#         if s!=t:
#             df1=pd.read_csv('H:/FYP/New_Data/Pearson/Pearson_'+s+"_gene", index_col=[0], sep='\t')
#             for ge in df1.index:
#                 name=ge.split(',')
#                 if name[0] not in total_gene:
#                     total_gene.append(name[0])
#                 if name[1] not in total_gene:
#                     total_gene.append(name[1])
#
#     df2=pd.read_csv('H:/FYP/New_Data/Pearson/Pearson_'+t+"_gene", index_col=[0], sep='\t')
#     for wei in df2.index:
#         ni=wei.split(',')
#         if ni[0] not in total_gene and ni[0] not in diff_gene:
#             diff_gene.append(ni[0])
#         if ni[1] not in total_gene and ni[1] not in diff_gene:
#             diff_gene.append(ni[1])
#
#     new_data.loc[t,"gene"]=str(diff_gene)
#     new_data.loc[t,"Length"]=len(diff_gene)
#
# outputpath='H:/FYP/New_Data/diff_gene_Pearson'
# new_data.to_csv(outputpath,sep='\t',index=True,header=True)


#第二种筛选基因的方法，根据person系数排序，每个基因都选最大的10个
# total_gene=[]
# new_data=pd.DataFrame(columns=["gene","length"])
# for t in nm:
#     print(t)
#     df = pd.read_csv('H:/FYP/New_Data/Pearson/Pearson_' + t + "_gene", index_col=[0], sep='\t')
#     dic={}
#     for n in df.index:
#         dic[n]=df.loc[n,"1"]
#     dic=sorted(dic.items(),key=lambda item:item[1],reverse=True)
#     se_ge=[]
#     for k in range(0,len(dic)):
#         if len(se_ge) < 10:
#             name=dic[k][0].split(',')
#             if name[0] not in se_ge:
#                 se_ge.append(name[0])
#             if name[1] not in se_ge:
#                 se_ge.append(name[1])
#     new_data.loc[t, "gene"] = str(se_ge)
#     new_data.loc[t,"Length"]=len(se_ge)
#     for s in se_ge:
#         if s not in total_gene:
#             total_gene.append(s)
#
# new_data.loc["total", "gene"] = str(total_gene)
# new_data.loc["total", "Length"] = len(total_gene)
# outputpath='H:/FYP/New_Data/top_gene_Pearson'
# new_data.to_csv(outputpath,sep='\t',index=True,header=True)




# 生成矩阵top

# df1=pd.read_csv('H:/FYP/New_Data/top_gene_Pearson', index_col=[0], sep='\t')
# gene=df1.loc["total",'gene']
# gene=gene.split('\', \'')
# gene[0]=str(gene[0])[2:len(gene[0])]
# gene[301]=str(gene[301])[0:len(gene[301])-2]
# gene.append("Label")
# nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
#     "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
#     "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]
#
# round=True
# for t in nm:
#     df2=pd.read_csv("H:/FYP/Data_Pre/" + t + ".RSEM_genes_normalized__data_Level_3", index_col=[0], sep='\t')
#     if round==True:
#         DF=pd.DataFrame(df2,index=gene)
#         round=False
#     else:
#         DF1=pd.DataFrame(df2,index=gene)
#         DF=pd.concat([DF,DF1],axis=1)
#
#
#
#
# outputpath='H:/FYP/New_Data/Pearson_top_matrix'
# DF.to_csv(outputpath,sep='\t',index=True,header=True)



#生成矩阵，diff
nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
    "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
    "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]

df1=pd.read_csv('H:/FYP/New_Data/diff_gene_Pearson', index_col=[0], sep='\t')
total_gene=[]
for k in nm:
    gene=df1.loc[k,'gene']
    gene=gene.split('\', \'')
    gene[0]=str(gene[0])[2:len(gene[0])]
    s=len(gene)-1
    gene[s]=str(gene[s])[0:len(gene[s])-2]
    for p in gene:
        if p not in total_gene and p != "":
            total_gene.append(p)
total_gene.append("Label")

round=True
for t in nm:
    df2=pd.read_csv("H:/FYP/Data_Pre/" + t + ".RSEM_genes_normalized__data_Level_3", index_col=[0], sep='\t')
    if round==True:
        DF=pd.DataFrame(df2,index=total_gene)
        round=False
    else:
        DF1=pd.DataFrame(df2,index=total_gene)
        DF=pd.concat([DF,DF1],axis=1)

outputpath='H:/FYP/New_Data/Pearson_diff_matrix'
DF.to_csv(outputpath,sep='\t',index=True,header=True)











