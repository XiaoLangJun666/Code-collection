import pandas as pd

nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","GBMLGG","HNSC","KICH",
    "KIPAN","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
    "PRAD","READ","SARC","SKCM","STAD","STES","TGCT","THCA","THYM","UCEC","UCS","UVM"]
for t in nm:
    ### 读取文件，因为这里有两种文件名所以用了try
    try:
        df=pd.read_csv("H:/FYP/Data/"+t+".RSEM_genes_normalized__data.Level_3/"+t+".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",index_col=[0],sep='\t')
    except:
        df=pd.read_csv("H:/FYP/Data/"+t+".RSEM_genes_normalized__data.Level_3/"+t+".normalized__data.txt",index_col=[0],sep='\t')

    ### 删掉gene_id因为我觉得没有用
    df.drop(labels='gene_id',axis=0,inplace=True)

    ### 修改gene的名字，因为一些缺失
    gid=["LOC100130426|100130426","UBE2Q2P3|100133144","UBE2Q2P2|100134869",
          "HMGB1P1|10357",    "TIMM23|10431",    "MOXD2|136542",
          "LOC155060|155060",   "RNU12-2P|26823",    "SSX9|280660",
          "LOC317712|317712",   "CXorf67|340602",   "EFCAB8|388795",
          "SRP14P1|390284",   "LOC391343|391343",   "TRIM75P|391714",
          "SPATA31B1P|404770",   "REXO1L6P|441362",   "SDR16C6P|442388",
          "LOC553137|553137",   "KIAA1618|57714",    "LOC645851|645851",
          "RGPD7|652919",   "HSPB1P1|653553",   "PPBPP1|728045",
          "FRMPD2P2|728603",   "ANKRD20A20P|728788",   "TMPRSS11E2|729884",
          "GTPBP6|8225",     "EFCAB12|90288"]
    for i in range (0,29):
        df.rename(index={df.index[i]:gid[i]},inplace=True)
    df.rename(index={"SLC35E2|728661":"SLC35E2B|728661"},inplace=True)

    ###通过TCGA的名字给sample打标签，因为要做肿瘤的分型，所以删除掉所有normal,新建Label行作为标签，
    for n in df.columns:
        a = n.split('-')
        if int(a[3][0:2]) != 1 and int(a[3][0:2]) != 3: #选择tcga barcode中为1和3的作为原位癌的数据
            df.drop(labels=n, axis=1, inplace=True)

    # 删除超过25%表达值为0的基因
    for j in df.index:
        c=0
        for k in df.columns:
            if float(df.loc[j][k])==0:
                c = c + 1
        if c > int(df.shape[1]) * 0.25:
            df.drop(labels=j, axis=0, inplace=True)



    df.loc["Label"] =t
    '''
    Data_Pre 输出的是没有对基因数量经过clean的数据（删除25%以上为0的基因）。
    Data_Pre_Primary是经过处理clean的数据（删掉25%以上为0的基因）。
    '''
    # outputpath='H:/FYP/Data_Pre_Total/'+t+".RSEM_genes_normalized__data_Level_3"
    outputpath='H:/FYP/Data_Pre/'+t+".RSEM_genes_normalized__data_Level_3"
    df.to_csv(outputpath,sep='\t',index=True,header=True)


