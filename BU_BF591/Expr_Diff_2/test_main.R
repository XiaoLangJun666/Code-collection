#!/usr/bin/Rscript
source("main.R")
library(testthat)

csv <- paste0("data/verse_counts.tsv")

describe("Data Loading and Reduction", {
  
  test_df <- load_n_trim(csv)
  
  it("should have the correct column names", {
    expect_equal(names(test_df), c("vP0_1", "vP0_2", "vAd_1", "vAd_2"))
  })
  
  it("should have the correct number of rows and columns", {
    expect_equal(dim(test_df), c(55416, 4))
  })
  
  it("should be a data.frame", {
    expect_equal(class(test_df), "data.frame")
  })
})

describe("DESeq2 Functionality", {
  
  load("data/mock_counts_df.RData") # loads the counts_df object into env
  coldata <- data.frame(condition = rep(c("day4", "day7"), each=2),
                        type="paired-end")
  row.names(coldata) <- c("vP4_1", "vP4_2", "vP7_1", "vP7_2")
  deseq <- run_deseq(counts_df, coldata, 10, "condition_day7_vs_day4")
  
  it("should return a DESeqResults object", {
    expect_equal(class(deseq)[1], "DESeqResults")
  })
  
  it("should have the correct dimensions", {
    expect_equal(dim(deseq), c(19127, 6))
  })
  
  it("should have pvalue and padj columns", {
    expect_equal(c("pvalue", "padj") %in% names(deseq), c(TRUE, TRUE))
  })
})

describe("edgeR Functionality", {
  
  load("data/mock_counts_df.RData")
  group <- factor(rep(c(1,2), each=2))
  edger_res <- run_edger(counts_df, group)
  
  it("should have the correct column names", {
    expect_equal(names(edger_res), c("logFC", "logCPM", "PValue"))
  })
  
  it("should have the correct dimensions", {
    expect_equal(dim(edger_res), c(15026, 3))
  })
})

describe("Limma + Voom Integration", {
  
  load("data/mock_counts_df.RData")
  group <- factor(rep(c(1,2), each=2))
  design <- data.frame(day4=1, day4vsday7=c(0, 0, 1, 1))
  row.names(design) <- c("vP4_1", "vP4_2", "vP7_1", "vP7_2")
  voom_res <- run_limma(counts_df, design, group)
  
  it("should have the correct dimensions", {
    expect_equal(dim(voom_res), c(15026, 6))
  })
   
  it("should have the correct column names", {
    expect_equal(names(voom_res), c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))
  })
})
