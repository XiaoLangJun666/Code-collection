import pandas as pd
from statsmodels import robust
import numpy as np

# t='BLCA'
# try:
#     df = pd.read_csv(
#         "Data/" + t + ".RSEM_genes_normalized__data.Level_3/" + t + ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
#         index_col=[0], sep='\t')
# except:
#     df = pd.read_csv("Data/" + t + ".RSEM_genes_normalized__data.Level_3/" + t + ".normalized__data.txt", index_col=[0],
#                      sep='\t')
# df.drop(labels='gene_id', axis=0, inplace=True)
# gid=["LOC100130426|100130426","UBE2Q2P3|100133144","UBE2Q2P2|100134869",
#           "HMGB1P1|10357",    "TIMM23|10431",    "MOXD2|136542",
#           "LOC155060|155060",   "RNU12-2P|26823",    "SSX9|280660",
#           "LOC317712|317712",   "CXorf67|340602",   "EFCAB8|388795",
#           "SRP14P1|390284",   "LOC391343|391343",   "TRIM75P|391714",
#           "SPATA31B1P|404770",   "REXO1L6P|441362",   "SDR16C6P|442388",
#           "LOC553137|553137",   "KIAA1618|57714",    "LOC645851|645851",
#           "RGPD7|652919",   "HSPB1P1|653553",   "PPBPP1|728045",
#           "FRMPD2P2|728603",   "ANKRD20A20P|728788",   "TMPRSS11E2|729884",
#           "GTPBP6|8225",     "EFCAB12|90288"]
# for i in range(0, 29):
#     df.rename(index={df.index[i]: gid[i]}, inplace=True)
# df.rename(index={"SLC35E2|728661": "SLC35E2B|728661"}, inplace=True)
# for t in df.columns:
#     if int(t[13:15]) >= 10 and int(t[13:15]) < 19:
#         df.drop(labels=t, axis=1, inplace=True)
# df.loc["Label"] = ("ACC")
# print(df.loc['Label'])


# df=pd.read_csv("Data_Pre/ACC.RSEM_genes_normalized__data_Level_3",index_col=[0],sep='\t')
# df.drop(labels='Label', axis=0, inplace=True)
# def get_median(data):
#     data.sort()
#     half = len(data) // 2
#     return (data[half] + data[~half]) / 2
# for i in df.columns:
#     median_sum=[]
#     median_1=df[i].median()
#     for j in range(0,len(df.index)):
#         median_sum.append(abs(float(df.iloc[j][i])-median_1))
#     median_2=get_median(median_sum)
#     print(median_2)

#看看是不是都是肿瘤
nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","GBMLGG","HNSC","KICH",
    "KIPAN","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
    "PRAD","READ","SARC","SKCM","STAD","STES","TGCT","THCA","THYM","UCEC","UCS","UVM"]
# nm=["SKCM"]

#找转移癌症
tj={}
for k in nm:
    tj[k]=0
for t in nm:
        try:
            df = pd.read_csv("Data/" + t + ".RSEM_genes_normalized__data.Level_3/" + t + ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",index_col=[0], sep='\t')
        except:
            df = pd.read_csv("Data/" + t + ".RSEM_genes_normalized__data.Level_3/" + t + ".normalized__data.txt",index_col=[0], sep='\t')
        for s in df.columns:
            a=s.split('-')
            if int(a[3][0:2])==6 or int(a[3][0:2])==7:
                tj[t]=tj[t]+1


print(tj)


# tj={}
# bh=range(1,10)
# for k in bh:
#     tj[k]=0
# c=0
# for t in nm:
#         df=pd.read_csv("Data_Pre/"+t+".RSEM_genes_normalized__data_Level_3",index_col=[0],sep='\t')
#         for k in df.columns:
#             a=k.split('-')
#             if int(a[3][0:2])==6:
#                 c=c+1
#
#         for k in bh:
#             if int(a[3][0:2])==k:
#                 tj[k]=tj[k]+1
# print(tj)

###读取临床数据
# for t in nm:
#     df2=pd.read_csv("tcga_clinical_level4/gdac.broadinstitute.org_"+t+".Clinical_Pick_Tier1.Level_4.2016012800.0.0/"+t+".clin.merged.picked.txt",index_col=[0],sep='\t')
#     for  k in df2.columns:
#         a=k.split('-')
#         if a[2][0:2]
#     print(t)
# for t in nm:
#         df=pd.read_csv("Data_Pre/"+t+".RSEM_genes_normalized__data_Level_3",index_col=[0],sep='\t')
#         for t in df.columns:
#             new=t.replace('-','_')
#             df.rename(columns={t:new},inplace=True)
#         print(df.columns)



