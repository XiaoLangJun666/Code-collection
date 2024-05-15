import pandas as pd
from statsmodels import robust
import numpy as np


nm=["ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","HNSC","KICH",
    "KIPAN","KIRC","KIRP","LIHC","LUAD","LUSC","MESO","OV","PAAD",
    "PRAD","READ","SKCM","STAD","STES","TGCT","THCA","THYM","UCEC","UCS","UVM"]
# nm=["BLCA"]   #GBM，GBML，LAML，LGG,PCPG，SARC 这六个没有临床数据


# ##读取临床数据中可判断的原位癌的数量
# dc={}
# for k in nm:
#     dc[k]=0
# total=0
# sum=0

# for t in nm:
#     cal=0
#     #df1 = pd.read_csv("Data_Pre/" + t + ".RSEM_genes_normalized__data_Level_3", index_col=[0], sep='\t')
#     df2=pd.read_csv("tcga_clinical_level4/gdac.broadinstitute.org_"+t+".Clinical_Pick_Tier1.Level_4.2016012800.0.0/"+t+".clin.merged.picked.txt",index_col=[0],sep='\t')
#     for s in df2.columns:
#         c = df2.loc['pathology_M_stage',s]
#         b = df2.loc['pathology_T_stage',s]
#         if c is not np.nan and b is not np.nan:
#             k=str(c).lower()
#             j=str(b).lower()
#             if k=="m0" and j is not "t0" and j is not 'tx':
#                 cal+=1
#     total=total+len(df2.columns)
#     sum=sum+cal
#     dc[t]=(cal,len(df2.columns))
#
# print(dc)
# print(total)
# print(sum)
# print(sum/total)



dc={}
for k in nm:
    dc[k]=0

dui_total=0
fu_total=0
for t in nm:

        dui = 0
        fu = 0
        df1 = pd.read_csv("Data_Pre/" + t + ".RSEM_genes_nor                           malized__data_Level_3", index_col=[0], sep='\t')
        df2=pd.read_csv("tcga_clinical_level4/gdac.broadinstitute.org_"+t+".Clinical_Pick_Tier1.Level_4.2016012800.0.0/"+t+".clin.merged.picked.txt",index_col=[0],sep='\t')
        for  k in df2.columns:
            a=k.split('-')
            a=a[2]  #临床数据的tcga barcode中的序号
            a=a.upper()
            round=False

            for d in df1.columns:
                if round==False:
                    b=d.split('-')
                    c=b[3][0:2]
                    b=b[2]
                    if a==b:
                        round=True
                        dui=dui+1
                        if int(c)==1 or int(c)==3:
                            pm=df2.loc["pathology_M_stage",k]
                            pt=df2.loc["pathology_T_stage",k]
                            if pm is not np.nan and pt is not np.nan:
                                if pm == "m0" and pt is not "t0" and pt is not 'tx':
                                    fu=fu+1
        dc[t]=(fu,dui,fu/dui)
        dui_total=dui_total+dui
        fu_total=fu_total+fu

print(dc)
print(fu_total,dui_total,fu_total/dui_total)




