source('main.R')
library('testthat')


test_names <- c('ENSMUSG00000102693.2', 
                'ENSMUSG00000064842.3', 
                'ENSMUSG00000051951.6', 
                'ENSMUSG00000102851.2',
                'ENSMUSG00000103377.2',
                'ENSMUSG00000104017.2')

test_data <- data.frame(
  log2FoldChange = c(8, -8, 3, -4, 5, -6),
  pval = rep(c(.001, .01, .1), 2),
  padj = rep(c(.01, .10, .15), 2))

row.names(test_data) <- test_names

test_tib <- test_data %>% 
  as_tibble(rownames='genes') %>%
  mutate(volc_plot_status = c('UP', rep('NS', 2), 'DOWN', rep('NS', 2)))

describe('make_se()', {
  se <- make_se('data/verse_counts.tsv', 'data/sample_metadata.csv', c('vP0', 'vAd'))
  
  it("returns a SummarizedExperiment object", {
    expect_true(class(se) == "SummarizedExperiment")
  })
  it("returns the right matrix dimensions, 55416 rows and 4 columns", {
    expect_equal(dim(assays(se)$counts), c(55416, 4))
  })
  it("returns the colnames in coldata in the specified, correct order", {
    expect_true(all.equal(colnames(se), c('vP0_1', 'vP0_2', 'vAd_1', 'vAd_2')))
  })
  it("stores the counts as a matrix for DESeq2 in the SE object", {
    expect_true(is.matrix(assays(se)$counts))
  })
})

describe("return_deseq_res()", {
  se <- make_se('data/verse_counts.tsv', 'data/sample_metadata.csv', c('vP0', 'vAd'))
  deseq2_results <- return_deseq_res(se, ~timepoint)

  it("returns a list of the DESeq2 results and the dds object", {
    expect_equal(class(deseq2_results), 'list')
  })
  it("returns both the results and the dds object", {
    expect_equal(length(deseq2_results), 2)
  })
})

describe("label_res()", {
  test_subset <- test_tib %>% dplyr::select(log2FoldChange, padj, volc_plot_status)

  func_labels <- label_res(test_data, .10)
  func_cols <- func_labels %>% dplyr::select(log2FoldChange, padj, volc_plot_status)

  it("correctly labels each row by correct padj threshold and log2fc sign", {
    expect_equal(func_cols, test_subset)
  })
})

describe("DESeq2 results return the same directionality for most significant genes", {
  se <- make_se('data/verse_counts.tsv', 'data/sample_metadata.csv', c('vP0', 'vAd'))
  dds <- DESeqDataSet(se, design = ~timepoint)
  dds <- DESeq(dds)
  res <- results(dds) %>% as_tibble(rownames='genes')

  fc_sign_neg <- res %>%
    filter(genes == 'ENSMUSG00000026418.17') %>%
    dplyr::select(log2FoldChange) %>%
    pull()

  fc_sign_pos <- res %>%
    filter(genes == 'ENSMUSG00000002500.16') %>%
    dplyr::select(log2FoldChange) %>%
    pull()

  it("returns ENSMUSG00000026418.17 with a negative log2FC, if this fails,
     please check your factor reference levels", {
    expect_true(fc_sign_neg < 0)
  })
  it("returns ENSMUSG00000002500.16 with a positive log2FC, if this fails,
     please check your factor reference levels", {
    expect_true(fc_sign_pos > 0)
  })
})

describe("make_ranked_log2fc()", {
  results <- make_ranked_log2fc(test_tib, 'data/id2gene.txt')

  symbols <- c('4933401J01Rik', 'Gm37180', 'Xkr4', 'Gm18956', 'Gm37363', 'Gm26206')
  log2fc <- c(8,5,3,-4,-6,-8)
  test_vec <- setNames(log2fc, symbols)

  it("should return a named vector matching the known results", {
    expect_mapequal(results, test_vec)
  })
})

describe("run_fgsea()", {
  se <- make_se('data/verse_counts.tsv', 'data/sample_metadata.csv', c('vP0', 'vAd'))
  dds <- DESeqDataSet(se, design = ~timepoint)
  dds <- DESeq(dds)
  res <- results(dds) %>% as_tibble(rownames='genes')
  rnk_list <- make_ranked_log2fc(res, 'data/id2gene.txt')

  fgsea_res <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)

  it("returns a tibble of the fgsea output", {
    expect_true(is_tibble(fgsea_res))
  })
  it("has the right columns", {
    expect_true(all(c('NES', 'ES', 'pval', 'padj', 'log2err', 'size', 'leadingEdge') %in% names(fgsea_res)))
  })
  it("was run with the correct minSize", {
    expect_true(min(fgsea_res$size) >= 15)
  })
  it("was run with the correct maxSize", {
    expect_true(max(fgsea_res$size) <= 500)
  })
})

