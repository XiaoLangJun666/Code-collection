library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
#BiocManager::install("fgsea")
library('fgsea')


#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
    counts<-read_tsv(counts_csv)
    metadata<-read.csv(metafile_csv)
    counts<-counts%>%
      dplyr::select('gene',starts_with(selected_times))
    counts_matrix<-as.matrix(counts[-1])
    row.names(counts_matrix)<-counts$gene
    metadata<-metadata%>%
      filter(timepoint%in%selected_times)%>%
      dplyr::select(samplename,timepoint)%>%
      mutate(timepoint=factor(timepoint,levels=selected_times))
    
    dds<-SummarizedExperiment(
      assays=list(counts=counts_matrix),
      colData=metadata,
      rowData=rownames(counts_matrix)
    )
    return(dds)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
    design_formula<-formula(design)
    dds<-DESeqDataSet(se,design=design_formula)
    dds<-DESeq(dds)
    res<-as_tibble(results(dds))%>%
      mutate(genes=rownames(results(dds)))%>%
      relocate(genes)%>%
      arrange(pvalue)
    
    results<-list(dds,res)
    
    return(results)
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
    labeled_results<-deseq2_res%>%
      mutate(volc_plot_status=ifelse(padj<padj_threshold & log2FoldChange>0,"UP",ifelse(padj<padj_threshold &log2FoldChange<0,'DOWN',"NS")))
    return( as_tibble(labeled_results))
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
    pvalue_plot<-ggplot(labeled_results,aes(x=pvalue))+
      geom_histogram(bins=50,fill='skyblue',color='black')+
      labs(title='Histogram of raw pvalues obtained from DE analysis (vP0 vs. vAd)',
           x='pvalue',
           y='count')+
      theme(plot.title = element_text(hjust = 0.5))+
      theme_minimal() 
    return(pvalue_plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
    labeled_results_filtered<-labeled_results%>%
      filter(padj<padj_threshold)
    log2FoldChange_plot<-ggplot(labeled_results_filtered,aes(x=log2FoldChange))+
      geom_histogram(bins=100,fill='skyblue',color='black')+
      labs(title='Histogram of Log2FoldChanges obtained from DE analysis (vP0 vs. vAd)',
           x='log2FoldChange',
           y='count')+
      theme(plot.title = element_text(hjust = 0.5))+
      theme_minimal() 
    return(log2FoldChange_plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
    top_genes<-labeled_results[order(labeled_results$padj),]
    top_genes<-top_genes$genes[1:num_genes]
    normalized_counts_selected<-as_tibble(counts(dds_obj,normalized=TRUE)[top_genes,])
    normalized_counts_selected<-log10(normalized_counts_selected)
    normalized_counts_selected$genes<-top_genes
    flap_counts<-pivot_longer(normalized_counts_selected,
                              c(colnames(normalized_counts_selected[,1:4])),
                              names_to='samplenames',
                              values_to='counts')
    p<-ggplot(flap_counts,aes(x=genes,y=counts,color=samplenames))+
      geom_jitter()+
      theme_minimal()+
      labs(title='Plot of Log10(normalized counts) for top ten DE genes',
           x=NULL,
           y='log10(norm_counts)')+
      theme(axis.text.x=element_text(angle=90))
    return(p)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
    labeled_results_add<-mutate(labeled_results,
                                '-log10(padj)'=-log10(padj),
    )%>%
      filter(!is.na(volc_plot_status))
    p<-ggplot(labeled_results_add,aes(x=log2FoldChange,y=-log10(padj),color=volc_plot_status))+
      geom_point()+
      theme_minimal()+
      labs(title='Volcano plot of DESeq2 differential expression results (vP0 vs. vAd)')
    return(p)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
    id2gene<-read.csv(id2gene_path,sep='\t',header=FALSE)
    labeled_results<-labeled_results
    results_symbols<-left_join(labeled_results,id2gene,by=c(genes='V1'))
    result_tibble<-results_symbols%>%
      arrange(desc(log2FoldChange))%>%
      rename('V2'='symbols')%>%
      dplyr::select(symbols,log2FoldChange)%>%
      filter(!is.na(log2FoldChange))
    
    rnk_list<-setNames(result_tibble$log2FoldChange,result_tibble$symbols)
    return(rnk_list)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
    gene_set_collection<-gmtPathways(gmt_file_path)
    fgsea_result<-fgsea(pathways=gene_set_collection,
                        stats=rnk_list,
                        minSize=min_size,
                        maxSize=max_size)
    return(as_tibble(fgsea_result))
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){

    return(NULL)
}

