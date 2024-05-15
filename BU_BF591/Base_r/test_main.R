#!/usr/bin/Rscript

library(testthat)

# if you change the name of your script, this line must be changed as well
source("main.R")


test_that("less_than_zero() tests", {
  expect_equal(less_than_zero(-1), TRUE)
  expect_equal(less_than_zero(0), FALSE)
  expect_equal(less_than_zero(1), FALSE)
  expect_equal(
    less_than_zero(c(-1,0,1,10)),
    c(TRUE,FALSE,FALSE,FALSE),
    desc="less_than_zero() should be able to accept a numeric vector and return a logical vector"
  )
  expect_equal(
    less_than_zero(matrix(-2:1, nrow=2, byrow=T)),
    matrix(c(TRUE, TRUE, FALSE, FALSE), nrow=2, byrow=T),
    desc="less_than_zero() should be able to accept a numeric matrix and return a logical matrix"
  )
})

test_that("is_between() tests", {
  expect_equal(is_between(1,0,5), TRUE)
  expect_equal(is_between(-1,0,5), FALSE)
  expect_equal(
    is_between(c(-1,1,2),0,5),
    c(FALSE, TRUE, TRUE),
    desc="is_between() should be able to accept a numeric vector and return a logical vector"
  )
  expect_equal(
    is_between(matrix(1:4, nrow=2, byrow=T),0,3),
    matrix(c(TRUE, TRUE, FALSE, FALSE), nrow=2, byrow=T),
    desc="is_between() should be able to accept a numeric matrix and return a logical matrix"
  )
})

test_that("rm_na() tests", {
  expect_equal(rm_na(3), 3)
  expect_equal(rm_na(c(1,2,3)), c(1,2,3))
  expect_equal(rm_na(c(1,2,NA)), c(1,2))
})

test_that("row_medians() tests", {
  m <- matrix(1:9, nrow=3, byrow=T)
  expect_equal(row_medians(m), c(2,5,8))
})

test_that("summarize_rows() produces correct result when computing min()", {
  m <- matrix(1:9, nrow=3, byrow=T)
  expect_equal(!!summarize_rows(m, min), c(1,4,7))
})
test_that("summarize_rows() produces correct result when computing mean()", {
  m <- matrix(1:9, nrow=3, byrow=T)
  expect_equal(summarize_rows(m, mean), c(2,5,8))
})

test_that("summarize_matrix() tests", {
  m <- matrix(1:9, nrow=3, byrow=T)
  m_summary <- summarize_matrix(m, na.rm=FALSE)
  expect_true(is.data.frame(m_summary))
  expect_equal(
    colnames(m_summary),
    c("mean", "stdev", "median", "min", "max", "num_lt_0", "num_btw_1_and_5", "num_na"),
  )
  expect_equal(m_summary$mean, c(2,5,8))
  expect_equal(m_summary$stdev, c(1, 1, 1))
  expect_equal(m_summary$median, c(2, 5, 8))
  expect_equal(m_summary$min, c(1, 4, 7))
  expect_equal(m_summary$max, c(3, 6, 9))
  expect_equal(m_summary$num_lt_0, c(0, 0, 0))
  expect_equal(m_summary$num_btw_1_and_5, c(2, 1, 0))
  expect_equal(m_summary$num_na, c(0, 0, 0))
})

# these tests are bonus - if you want to try the challenge, uncomment!
test_that("Bonus tests! Make your code work when there are NAs!", {
 m <- matrix(1:9, nrow=3, byrow=T)
 m[2,1] <- NA
 expect_equal(
   summarize_rows(m, mean, na.rm=FALSE),
   c(2,NA,8),
 )
 expect_equal(
   summarize_rows(m, mean, na.rm=TRUE),
   c(2,5.5,8),
 )
 expect_equal(
   summarize_rows(m, min, na.rm=FALSE),
   c(1,NA,7),
 )
 expect_equal(
   summarize_rows(m, min, na.rm=TRUE),
   c(1,5,7),
 )

 m_summary <- summarize_matrix(m, na.rm=FALSE)
 expect_equal(m_summary$mean, c(2,NA,8))
 expect_equal(m_summary$stdev, c(1, NA, 1))
 expect_equal(m_summary$median, c(2, NA, 8))
 expect_equal(m_summary$min, c(1, NA, 7))
 expect_equal(m_summary$max, c(3, NA, 9))
 expect_equal(m_summary$num_lt_0, c(0, NA, 0))
 expect_equal(m_summary$num_btw_1_and_5, c(2, NA, 0))
 expect_equal(m_summary$num_na, c(0, 1, 0))

 m_summary <- summarize_matrix(m, na.rm=TRUE)
 expect_equal(m_summary$mean, c(2,5.5,8))
 expect_equal(m_summary$stdev, c(1, sd(c(5,6)), 1))
 expect_equal(m_summary$median, c(2, 5.5, 8))
 expect_equal(m_summary$min, c(1, 5, 7))
 expect_equal(m_summary$max, c(3, 6, 9))
 expect_equal(m_summary$num_lt_0, c(0, 0, 0))
 expect_equal(m_summary$num_btw_1_and_5, c(2, 0, 0))
 expect_equal(m_summary$num_na, c(0, 0, 0))
})