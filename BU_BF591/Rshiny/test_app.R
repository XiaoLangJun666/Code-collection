library(shiny)
library(testthat)

source("app.R")

testServer(app = server, {
  # load_data()
  session$setInputs(file = list(datapath = "data/deseq_res.csv"))
  # expect_equal(dim(load_data()), c(19407, 7))

  # test volcano_plot() function
  read.csv("data/deseq_res.csv")
    .[order(.$padj), ] %>%
    head() %>%
    volcano_plot(., "log2FoldChange", "padj",-200, "blue", "green") -> plotted
  expect_equal(class(plotted$layers[[1]]$geom)[1], "GeomPoint")

  # # test volcano_plot() output
  # suppressWarnings(
  #   session$setInputs(
  #     file = list(datapath = "data/deseq_res.csv"),
  #     x_axis = "log2FoldChange",
  #     y_axis = "padj",
  #     slider = -200,
  #     base = "blue",
  #     highlight = "green"
  #   )
  # )
  # expect_true(grepl("data:image", output$volcano[["src"]]))
  # 
  # # test table
  # table_obj <- draw_table(read.csv("data/deseq_res.csv"), -210)
  # expect_equal(dim(table_obj), c(2, 7))
  # expect_true("ENSMUSG00000026418.17" %in% table_obj[,1])
  # 
  # # test table output
  # suppressWarnings(session$setInputs(
  #   file = list(datapath = "data/deseq_res.csv"),
  #   slider = -210,
  # ))
  # expect_true(grepl("ENSMUSG00000026418.17", output$table))
  # expect_true(grepl("table shiny-table", output$table))
})
