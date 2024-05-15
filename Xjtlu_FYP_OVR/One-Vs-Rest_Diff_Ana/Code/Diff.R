library(DESeq2)
library(dplyr)
# 胶质瘤组：GBM，LGG，GBMLGG。食管胃癌组：STAD，STES，ESCA
#读取表达
df<-read.csv('/data/haochun/fyp/TCGA_diff_int/ESCA_STAD_gene',sep='\t',header=TRUE,row.names=1)
#读取标签
lab<-read.csv('/data/haochun/fyp/TCGA_diff_int/ESCA_STAD',sep='\t',header=TRUE,row.names=1)
#数据过滤
df=df[rowMeans(df)>1,]
all(rownames(lab)%in%colnames(df))
all(rownames(lab)==colnames(df))
#制作差异矩阵
dds <-  DESeqDataSetFromMatrix(countData = df,colData = lab,design = ~ condition) 

dim(dds)
#过滤
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)
#差异比较
dep <- DESeq(dds)
res <- results(dep)
diff = res
diff <- na.omit(diff)  ## 去除NA
dim(diff)
foldChange = 2
padj = 0.05
#
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)
rownames(diffsig)
write.csv(diffsig, "/data/haochun/fyp/TCGA_diff_int/GBM_LGG_diffsig.csv")