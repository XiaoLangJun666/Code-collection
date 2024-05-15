if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
library('GEOquery')
gse <- getGEO(filename="GSE64810_series_matrix.txt.gz",GSEMatrix=TRUE)
gse$submission_date
gse
show(pData(phenoData(gse[[1]]))[1:5,c(1,6,8)])
cd<-pData(phenoData(gse))
class(cd)
round(sd(as.numeric(gse$`age of onset:ch1`),na.rm=TRUE),2)
na.pass()


gse$`diagnosis:ch1`
gse$`age of death:ch1`

plot_data <- data.frame(
  diagnosis = as.factor(gse$`diagnosis:ch1`),
  y_value = as.numeric(as.character(gse$`age of death:ch1`))
)
plot_data <- na.omit(plot_data)
ggplot(plot_data, aes(x = diagnosis, y = y_value,fill=diagnosis)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Age of Death by Diagnosis", 
       x = "Diagnosis", 
       y = "Age of Death")

ggplot(plot_data, aes(x = y_value, fill = diagnosis)) +
  geom_histogram(bins = 10, color = "black") + # 设置binwidth，您可以根据需要调整
  labs(title = "Histogram of H-V Cortical Scores by Diagnosis",
       x = "H-V Cortical Score",
       y = "Frequency") +
  theme_minimal()
ggplot(plot_data, aes(x = y_value, fill = diagnosis)) +
  geom_density(alpha = 0.5) + # 设置透明度，以便看到重叠部分
  scale_fill_brewer(palette = "Set1") + # 使用预定义的颜色板
  labs(title = "Density Plot of Y Values by Diagnosis",
       x = "Y Value",
       y = "Density") +
  theme_minimal()


#counts
cm<-read.csv('GSE64810_mlhd_DESeq2_norm_counts_adjust.txt',sep='\t')
var<-apply(cm[-1],1,var)
var_th<-quantile(var,probs=0.95)
var_th
cm1<-cm[apply(cm[-1],1,var)>=var_th,]

cm1<-cm1[apply(cm1[-1],1,function(x) sum(x!=0))>60,]

nrow(cm1)

scatter_data <- data.frame(gene_name = cm[1])
scatter_data$median_count<-apply(cm[-1],1,median)
scatter_data$filtered <- cm[[1]] %in% cm1[[1]]
scatter_data$variance <- apply(cm[-1], 1, var)
scatter_data$num_zeros <- apply(cm[-1], 1, function(x) sum(x == 0))

sum(scatter_data$filtered==TRUE)


p1 <- ggplot(scatter_data, aes(x = median_count, y = variance, color = filtered)) +
  geom_point() +
  scale_color_manual(values = c("darkblue","lightgray")) +
  scale_y_log10() +
  scale_x_log10()+
  ggtitle("Median Count vs Variance")


rownames(cm2)<-cm2[[1]]

cm2<-cm2[-1]
pca_results <- prcomp(scale(t(cm2)), center=FALSE, scale=FALSE)
summary(pca_results)

PCA_tibble<-as_tibble(pca_results$x)%>%
  mutate(Sample=ifelse(substr(rownames(pca_results$x),1,1)=='H','Huntington',ifelse(substr(rownames(pca_results$x),1,1)=='C','Control',NA)))%>%
  select(Sample,PC1,PC2)
  
ggplot(PCA_tibble,aes(x=PC1,y=PC2,color=Sample))+
  labs(color = "SixSubtypesClassification")+
  geom_point()+
  theme_classic()

cm1[-1]<-lapply(names(cm1[-1]),function(x) cm1[-1][[x]]*10^6/vapply(cm1[-1],sum,1)[[x]])



count_mat<-as.matrix(cm1[-1])
rownames(count_mat)<-cm1[[1]]
dds<-DESeqDataSetFromMatrix(
  countData=round(count_mat),
  colData=tibble(sample_name=colnames(cm1[-1])),
  design=~1
)
dds<-estimateSizeFactors(dds)
deseq_norm_counts<-as_tibble(counts(dds,normalized=TRUE))

cm2<-deseq_norm_counts


log_transformed_counts <- log1p(cm1[-1])
pheatmap(log_transformed_counts, 
         color = colorRampPalette(c("blue", "white", "red"))(50), # 颜色渐变
         scale = "row", # 按行标准化
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         legend = TRUE)


pca_results <- prcomp(t(cm1[-1]),center=FALSE, scale=TRUE)
PCA_tibble<-as_tibble(pca_results$x)%>%
  mutate(Sample=ifelse(substr(rownames(pca_results$x),1,1)=='H','Huntington',ifelse(substr(rownames(pca_results$x),1,1)=='C','Control',NA)))
prop_var <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100

ggplot(PCA_tibble, aes(x = PC1, y = PC2, color = Sample)) +
  geom_point() +
  labs(x = paste0("PC1: ", round(prop_var[1]), "% variance"),
       y = paste0("PC2: ", round(prop_var[2]), "% variance"))+
  theme_classic()

ggplot(PCA_tibble, aes_string(x = "PC1", y = paste0("PC", 1))) +
  geom_beeswarm() +
  theme_minimal() +
  labs(title = paste("Beeswarm Plot of Top", num_pcs, "Principal Components"))


PCA_tibble%>%
  pivot_longer(PC1:PC20,names_to="PC",values_to="projection")%>%
  mutate(PC=fct_relevel(PC,str_c("PC",1:20))) %>%
  ggplot(aes(x=PC,y=projection,color=Sample)) +
  geom_beeswarm() + labs(title="PCA Projection Plot") +
  theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))
ncol(PCA_tibble)


pc1=1
pc2=2


pc_1<-paste0('PC',1)
pc_2<-paste0('PC',2)

ggplot(PCA_tibble, aes(x = .data[[pc_1]], y = .data[[pc_2]], color = Sample)) +
  geom_point() +
  labs(x = paste0(pc_1,': ', round(prop_var[pc1]), "% variance"),
       y = paste0(pc_2,': ', round(prop_var[pc2]), "% variance"))+
  theme_classic()


DE<-read.csv('GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz',sep='\t')
DE<-as_tibble(DE)%>%
  arrange(desc(log2FoldChange))%>%
  select(symbol,log2FoldChange)
library('fgsea')
rnk_list<-setNames(DE$log2FoldChange,DE$symbol)
head(rnk_list)
tail(rnk_list,10)
hallmark_pathways_fgsea <- fgsea::gmtPathways('msigdb_v2023.2.Hs_GMTs/c5.go.bp.v2023.2.Hs.symbols.gmt')

fgsea_results <- fgsea(hallmark_pathways_fgsea, rnk_list, minSize =1, maxSize=1500)

fgsea_results <- fgsea_results %>% as_tibble()


fgsea_results %>% arrange(padj)



neg_nes<-fgsea_results[fgsea_results$NES<0,]
neg_nes_sorted<-neg_nes[order(neg_nes$NES,decreasing=FALSE),]
pos_nes<-fgsea_results[fgsea_results$NES>0,]
pos_nes_sorted<-pos_nes[order(pos_nes$NES,decreasing=TRUE),]
top_NES<-rbind(head(pos_nes_sorted,10),
               head(neg_nes_sorted,10)
)
top_NES$color<-ifelse(top_NES$NES>0,"Positive","Negative")
top_NES$pathway <- str_sub(top_NES$pathway, 1, 40)

top_NES$short_pathway <- ifelse(nchar(top_NES$pathway) > 40, 
                             substr(top_NES$pathway, 1, 40), 
                             top_NES$pathway)



# top_NES<-top_NES%>%
#   arrange(desc(NES))
p<-ggplot(top_NES,aes(x=NES,y = reorder(short_pathway, NES),fill=color))+
  geom_bar(stat='identity')+
  scale_fill_manual(values = c("Positive" = "red", "Negative" = "blue"))+
  theme_minimal() +
  labs(title='fgsea results for Hallmark MSlgDB gene se',x='Normalized Enrichment Score (NES)')

ggplotly(p, tooltip = "y") %>%
  layout(hoverlabel = list(bgcolor = "white"),
         yaxis = list(title = "Pathway")) %>%
  style(hoverinfo = "text") %>%
  add_trace(text = ~top_NES$pathway)




library(plotly)
top_positive_nes <- fgsea_results %>%
  filter(padj < .25 & abs(NES) > 0.5) %>%
  slice_max(NES, n=5)


top_positive_nes %>%
  mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot() +
  geom_bar(aes(x=pathway, y=NES, fill = padj < .25), stat='identity') +
  scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
  theme_minimal() +
  ggtitle('fgsea results for Hallmark MSigDB gene sets') +
  ylab('Normalized Enrichment Score (NES)') +
  xlab('') +
  coord_flip()
