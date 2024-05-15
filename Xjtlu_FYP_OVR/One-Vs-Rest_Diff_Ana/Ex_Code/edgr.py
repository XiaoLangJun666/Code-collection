import pandas as pd
import numpy as np
import math

df = pd.read_csv("/data/haochun/fyp/TCGA_union/Data_TCGA_union_Label", index_col=[0], sep='\t')


for c in range(0,len(df.columns)):
    for s in range(0,len(df.index)-1) :
        df.iloc[s,c]=int(float(df.iloc[s,c]))
        s=s+1
    c=c+1
for k in df.columns:
    k1=df.loc["Label",k]
    k2=k+"-"+k1
    df.rename(columns={k:k2},inplace=True)#如果先重命名，那在下面就找不到列名

outputpath='/data/haochun/fyp/TCGA_union/diff'
df.to_csv(outputpath,sep='\t',index=True,header=True)


