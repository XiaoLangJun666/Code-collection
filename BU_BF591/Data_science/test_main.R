#!/usr/bin/Rscript
source("main.R")
library(testthat)

#generate consistent sample data for all tests using it
set.seed(42)
fx_intensity_df <- data.frame(GSM971958 = runif(1000, 0, 15),
                              GSM971959 = runif(1000, 0, 15),
                              GSM971967 = runif(1000, 0, 15),
                              GSM971968 = runif(1000, 0, 15))

row.names(fx_intensity_df) <- paste0(1:1000, "_s_at")
fx_intensity_mat <- as.matrix(fx_intensity_df)


test_that("read_data loads correctly", {
  #checks dimensions and is.data.frame()
  
  res <- read_data("data/example_intensity_data.csv", " ")
  
  expect_equal(dim(res), c(54675, 35))
  expect_true(is.data.frame(res))
  
})

test_that("calculate_variance squares each element in a vector and divides by sum", {
  #mimic the object by using a named list so input arg is the same

  test_obj <- list(sdev = c(8, 4, 2))
  test_res <- c(0.76190476, 0.19047619, 0.04761905)

  expect_true(all(dplyr::near(calculate_variance_explained(test_obj), test_res, tol=.1)))
})

test_that("make_variance_tibble", {
  # Checks if length of tibble is the same as number of PCs, verifies it's a tibble,
  # and ensures correct calculation of cumulative sum for a specific column.

  test_pca_res <- prcomp(scale(t(fx_intensity_mat)), scale=FALSE, center=FALSE)
  test_ve <- c(1, 1, 1, 1)

  function_tib <- make_variance_tibble(test_ve, test_pca_res)

  expect_equal(dim(function_tib), c(length(test_ve), 3))
  expect_true(is_tibble(function_tib))

  # New check for correct calculation of cumulative sum on specified column (e.g., "cumulative_variance")
  expected_cumsum <- cumsum(test_ve)
  expect_equal(function_tib$cumulative, expected_cumsum)
})


test_that("list_significant_probes returns the probeids with padj < threshold", {
  # Defining the test_tibble within the test block
  test_tibble <- tibble(
    probeid = c("1_s_at", "2_s_at", "3_s_at", "4_s_at"),
    padj = c(0.005, 0.008, 0.01, 0.015)
  )

  expect_equal(list_significant_probes(test_tibble, .01), c('1_s_at', '2_s_at'))
})

test_that("return_de_intensity returns a matrix and the selected probes", {
  #tests is.matrix() and since function subsets the matrix should be identical
  #no need to worry about rounding, dimensions, etc. ?

  test_mat <- fx_intensity_mat[c('1_s_at', '2_s_at'),]
  function_mat <- return_de_intensity(fx_intensity_df, c('1_s_at', '2_s_at'))

  expect_true(is.matrix(function_mat))
  expect_true(identical(function_mat, test_mat))
})

