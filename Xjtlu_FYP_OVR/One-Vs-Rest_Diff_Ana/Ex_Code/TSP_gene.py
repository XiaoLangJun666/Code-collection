import pandas as pd

df=pd.read_csv('H:/FYP/New_Data/Union/Rule_Gene.csv', sep='\t')

gene=[]
for k in range(0,df.shape[0]):
    if k % 2== 0:
        raft=df.iloc[k,0]
        raft=raft.split('\"')
        raft=raft[1]
        if raft not in gene:
            gene.append(raft)

df1=pd.read_csv('H:/FYP/New_Data/Join/Rule_Gene.csv', sep='\t')

gene1=[]
for k in range(0,df1.shape[0]):
    if k % 2== 0:
        raft=df1.iloc[k,0]
        raft=raft.split('\"')
        raft=raft[1]
        if raft not in gene1:
            gene1.append(raft)

gene2=[val for val in gene1 if val in gene]


print(len(gene),len(gene1),len(gene2))
# nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
#     "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
#     "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]
# df1=pd.read_csv('H:/FYP/New_Data/diff_gene_Pearson', index_col=[0], sep='\t')
# total_gene=[]
# for k in nm:
#     gene1=df1.loc[k,'gene']
#     gene1=gene1.split('\', \'')
#     gene1[0]=str(gene1[0])[2:len(gene1[0])]
#     s=len(gene1)-1
#     gene1[s]=str(gene1[s])[0:len(gene1[s])-2]
#     for p in gene1:
#         if p not in total_gene and p != "":
#             total_gene.append(p)
# print(len(total_gene))
# print(len(gene))
# real_gene=[]
# for j in gene:
#     if j in total_gene:
#         real_gene.append(j)
# print(len(real_gene))



# group=[("COAD","READ"),("LUSC","LUAD"),("ESCA","STAD"),("CHOL","LIHC")]
# #group=[("ESCA","STAD")]
# sel=[]
# for s in group:
#     name1=s[0]
#     name2=s[1]
#     df2=pd.read_csv('H:/FYP/New_Data/Select_Diff/'+name1+"_"+name2+'_result.csv',sep=',',index_col=[0])
#     ix=df2.index
#     print(len(ix))
#     for k in ix:
#         if k not in sel:
#             sel.append(k)

# for t in sel:
#     if t not in gene2:
#         gene2.append(t)

# gene2.append("Label")
#
# nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
#     "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
#     "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]
#
# round=True
# for j in nm:
#     df2=pd.read_csv("H:/FYP/Data_Pre/" + j + ".RSEM_genes_normalized__data_Level_3", index_col=[0], sep='\t')
#     if round==True:
#         DF=pd.DataFrame(df2,index=gene2)
#         round=False
#     else:
#         DF1=pd.DataFrame(df2,index=gene2)
#         DF=pd.concat([DF,DF1],axis=1)
#
# for j in DF.index:
#         c=0
#         if j!="Label":
#             for k in DF.columns:
#                 if float(DF.loc[j][k])==0:
#                     c = c + 1
#             if c > int(DF.shape[1]) * 0.25:
#                 DF.drop(labels=j, axis=0, inplace=True)
#
# outputpath='H:/FYP/New_Data/TSP_UNION_JOIN_JOIN_matrix'
# DF.to_csv(outputpath,sep='\t',index=True,header=True)



