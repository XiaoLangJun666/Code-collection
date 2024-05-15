#Imports
library(tidyverse)
library(DESeq2)

#' Load a tsv located at specific location `filename` into a tibble
#'
#'
#' @param filename (str): the path to a specific file (ie 'file/path/to/file.tsv')
#'
#' @return tibble: a (g x 1+m) tibble with a 'gene' column followed by
#' sample names as column names.
#'
#' @note Column 'gene' should be first and the only column to contain strings.
#' Data in sample_name columns CANNOT be strings
#'
#' @example `verse_counts <- read_data('verse_counts.tsv')`

read_data <- function(filename){
  mt<-read_tsv(filename)  
  return(mt)
}


#' Filter out genes with zero variance
#'
#'
#' @param verse_counts tibble: a (g x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns with sample names as column names.
#'
#' @return tibble: a (n x 1+m) tibble with a 'gene' column followed by m columns
#' of raw counts with genes that have zero variance across samples removed
#'
#' @note (g >= n)
#'
#' @example `filtered_counts <- filter_zero_var_genes(verse_counts)`

filter_zero_var_genes <- function(verse_counts) {
    return(verse_counts[apply(verse_counts[-1],1,var)!=0,])
}


#' Extract time point information from sample name
#'
#'
#' @param str string: sample name from count data.
#'
#' @return string: string character representing sample time point
#'
#' @example `timepoint_from_sample("vAd_1")`
#' output:`"Ad"`

timepoint_from_sample <- function(x) {
    former<-strsplit(x,'_')[[1]][1]
    timepoint<-strsplit(former,'v')[[1]][2]
    return(timepoint)
}


#' Grab sample replicate number from sample name
#'
#'
#' @param str  string: sample name from count data.
#'
#' @return string: string character represent sample replicate number
#'
#' @example `sample_replicate("vAd_1")`
#' output: `"1"`

sample_replicate <- function(x) {
    repli<-strsplit(x,'_')[[1]][2]
    return(repli)
}


#' Generate sample-level metadata from sample names.
#'
#' Will include columns named "sample", "timepoint", and "replicate" that store
#' sample names, sample time points, and sample replicate, respectively.
#'
#'
#' @param sample_names vector: character vector of length (_S_) consisting of sample
#' names from count data.
#'
#' @return tibble: a (_S_ x 3) tibble with column names "sample",
#' "timepoint", and "replicate". "sample"holds sample_names; "timepoint"
#' stores sample time points; and "replicate" stores sample replicate
#'
#' @note _S_ < m
#'
#' @example `meta <- meta_info_from_labels(colnames(count_data)[colnames(count_data)!='gene'])`

meta_info_from_labels <- function(sample_names) {
  timepoint<-vapply(sample_names,timepoint_from_sample,'a')
  repli<-vapply(sample_names,sample_replicate,'a')
  tb<-tibble('sample'=sample_names,
             'timepoint'=timepoint,
             'replicate'=repli)
  return(tb)
}


#' Calculate total read counts for each sample in a count data.
#'
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @return tibble or named vector of read totals from each sample. Vectors must
#' be length `_S_ `, a tibble can be `(1 x _S_)` with sample names as columns
#' names OR `(_S_ x 2)` with columns ("sample", "value")
#'
#' @examples `get_library_size(count_data)`

get_library_size <- function(count_data) {
    return(vapply(count_data[-1],sum,1))
}


#' Normalize raw count data to counts per million WITH pseudocounts using the
#' following formula:
#'     count / (sample_library_size/10^6)
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m columns of cpm normalized read counts
#'
#' @examples
#' `normalize_by_cpm(count_data)`

normalize_by_cpm <- function(count_data) {
    count_data[-1]<-lapply(names(count_data[-1]),function(x) count_data[-1][[x]]*10^6/get_library_size(count_data)[[x]])

    return(count_data)
}

#' Normalize raw count matrix using DESeq2
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts

#' @param meta_data tibble: sample-level information tibble corresponding to the
#' count matrix columns
#'
#' @return tibble: DESeq2 normalized count matrix
#' @export
#'
#' @examples
#' `deseq_normalize(count_data, meta_data)`
deseq_normalize <- function(count_data, meta_data) {
    count_mat<-as.matrix(count_data[-1])
    row.names(count_mat)<-count_data$gene
    dds<-DESeqDataSetFromMatrix(
      countData=count_mat,
      colData=tibble(sample_name=colnames(count_data[-1])),
      design=~1
    )
    dds<-estimateSizeFactors(dds)
    deseq_norm_counts<-as_tibble(counts(dds,normalized=TRUE))%>%
      mutate(gene=count_data$gene)%>%
      relocate(gene)
    
    return(deseq_norm_counts)
}


#' Perform and plot PCA using processed data.
#'
#' PCA is performed over genes, and samples should be colored by time point.
#' Both `y` and `x` axis should have percent of explained variance included.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param meta tibble: sample-level meta information (_S_ x 3)
#' @param title string: title for plot
#'
#' @return ggplot: scatter plot showing each sample in the first two PCs.
#'
#' @examples
#' `plot_pca(data, meta, "Raw Count PCA")`

plot_pca <- function(data, meta, title="") {
    pca_result <- prcomp(t(data))
    pca_data <- as.data.frame(pca_result$x[, 1:2])
    colnames(pca_data) <- c("PC1", "PC2")
    pca_data$sample=rownames(pca_data)
    pca_data<-left_join(pca_data,meta,by='sample')
    prop_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
    p<-ggplot(pca_data, aes(x = PC1, y = PC2, color = timepoint)) +
        geom_point() +
        labs(x = paste0("PC1: ", round(prop_var[1]), "% variance"),
             y = paste0("PC2: ", round(prop_var[2]), "% variance"),title=title)
    return(p)
}


#' Plot gene count distributions for each sample using boxplots.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale the `y` axis to log10 values.
#' Default is FALSE, and y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: boxplot show gene count distributions for each sample
#'
#' @example `plot_sample_distributions(data, scale_y_axis=TRUE, title='Raw Count Distributions')`

plot_sample_distributions <- function(data, scale_y_axis=FALSE, title="") {
    x<-pivot_longer(
        data,
        c(colnames(data)),
        names_to='sample',
        values_to='counts'
      )%>%
        mutate(sample=factor(sample,levels=colnames(data)))
    
    p<-ggplot(x,aes(x=sample,y=counts,fill=sample,color=sample))
    
    if (scale_y_axis==TRUE){
      return(p+geom_boxplot(fill='white')+labs(title = title)+scale_y_log10())
    }
    else{
      return(p+geom_boxplot(fill='white')+labs(title = title))
    }
    
}


#' Plot relationship between mean read counts and variability over all genes.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale to y-axis to log10 values. Default
#' is false, and the y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: A scatter plot where the x-axis is the rank of gene ordered by mean
#' count over all samples, and the y-axis is the observed variance of the
#' given gene. Each dot should have their transparency increased. The scatter
#' plot should also be accompanied by a line representing the average mean and
#' variance values.
#'
#' @example `plot_variance_vs_mean(data, scale_y_axis=TRUE, title='variance vs mean (raw counts)')`

plot_variance_vs_mean <- function(data, scale_y_axis=FALSE, title="") {
    gene_stats<-tibble(
      mean_count=rowMeans(data),
      variance=apply(data,1,var)
    )
    gene_stats <- gene_stats %>%
      arrange(mean_count) %>%
      mutate(rank = rank(mean_count,ties.method = 'min'))
    p<-ggplot(gene_stats, aes(x = rank, y = variance)) +
      geom_point(alpha=0.5,stat='identity') +
      geom_smooth(method = "auto", se = FALSE, color = "blue") +  
      labs(x = "Rank(Mean)", y = "Variance",title=title)
    if (scale_y_axis==TRUE){
      p<-p+scale_y_continuous(trans = "log10")
    }
    
    return(p)
}

