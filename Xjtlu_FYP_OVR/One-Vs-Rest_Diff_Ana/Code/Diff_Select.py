import pandas as pd

#生成基因表格
nm = ["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
    "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
    "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]


dic={}
new_data=pd.DataFrame(columns=["gene","length"])

for t in nm:
    df=pd.read_csv('H:/FYP/summer/Data/Chayi_reduced/Reduced_result/Result_'+t+"_Others",sep=',',index_col=[0])
    dic1={}
    lst=[]
    for k in df.index:
        dic1[k]=abs(df.loc[k,"log2FoldChange"])
    dic1=sorted(dic1.items(),key=lambda item:item[1],reverse=True)
    dic_final=dic1[0:500]  #这一步可以选择我们的H的值，最大的前500个这里是
    for s in dic_final:
        lst.append(s[0])
    dic[t]=lst


for n in nm:
    gene1=dic[n]
    gene_rest=[]
    gene_select=[]
    nm2 = ["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
    "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
    "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]
    for m in nm2:
        if n!=m:
            gene2=dic[m]
            for j in gene2:
                if j not in gene_rest:
                    gene_rest.append(j)

    for i in gene1:
        if i not in gene_rest:
            gene_select.append(i)

    new_data.loc[n,"gene"]=str(gene_select)
    new_data.loc[n,'length']=len(gene_select)

outputpath="H:/FYP/summer/Data/Chayi_reduced/Reduced_Gene_list/Gene_list_500"
new_data.to_csv(outputpath, sep='\t', index=True, header=True)


#生成训练矩阵
nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH",
    "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
    "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"]

number=['450','500']

for q in number:
    print(q)
    df1=pd.read_csv('H:/FYP/summer/Data/Chayi_reduced/Reduced_Gene_list/Gene_list_'+q, index_col=[0], sep='\t')
    total_gene=[]
    for k in nm:
        gene1=df1.loc[k,'gene']
        gene1=gene1.split('\', \'')
        gene1[0]=str(gene1[0])[2:len(gene1[0])]
        s=len(gene1)-1
        gene1[s]=str(gene1[s])[0:len(gene1[s])-2]
        for p in gene1:
            if p not in total_gene and p != "":
                total_gene.append(p)


    total_gene.append("Label")
    round=True

    for t in nm:
        mt=pd.read_csv("H:/FYP/summer/Data/Chayi_reduced/Data_Pre/" + t + ".RSEM_genes_normalized__data_Level_3", index_col=[0], sep='\t')
        if round ==True:
            pd1=pd.DataFrame(mt,index=total_gene)
            round=False
        else:
            pd2=pd.DataFrame(mt,index=total_gene)
            pd1=pd.concat([pd1,pd2],axis=1)

    for k in pd1:
        name=k.split('|')[0]
        pd1.rename({k:name},inplace=True)

    outputpath="H:/FYP/summer/Data/Chayi_reduced/Reduced_Gene_list/Gene_matrix_"+q
    pd1.to_csv(outputpath, sep='\t', index=True, header=True)









