#!/usr/bin/Rscript
source("main.R")
library(testthat)
library(tidyverse)


# Test for read_expression_table function
test_that("read_expression_table returns tibble with subject_id column", {
  result <- read_expression_table('data/example_intensity_data_subset.csv')
  
  # Check for subject_id in column names
  expect_true("subject_id" %in% names(result), 
              info = "subject_id column is missing from the returned tibble. Make sure to include this column in your function output.")
  
  # Check if result is a tibble
  expect_true(is_tibble(result), 
              info = "The returned result is not a tibble. Ensure your function returns a tibble.")
})



test_that("load_metadata function loads metadata correctly", {
  # The expected output using readr's read_csv function
  expected_metadata <- readr::read_csv("data/proj_metadata.csv")

  # The output from the developer's function
  calculated_metadata <- load_metadata("data/proj_metadata.csv")

  # Test: Compare the content of the expected and calculated metadata
  expect_equal(expected_metadata, calculated_metadata,
               info = "The loaded metadata does not match the expected content.")
})




# Test for period_to_underscore function
test_that("period_to_underscore replaces period with underscore", {
  result <- period_to_underscore("foo.bar.baz")

  # Check string manipulation
  expect_equal(result, "foo_bar_baz",
               info = "The function did not replace periods with underscores correctly. Check the string manipulation in your function.")
})


# Test for rename_and_select function
test_that("rename_and_select renames and selects the correct columns", {

  fake_tibble <- tibble(Age_at_diagnosis = rep(10, 10), SixSubtypesClassification = rep('C4', 10), normalizationcombatbatch=rep('batch1', 10), Sex = rep(10, 10), TNM_Stage = rep(10, 10), Tumor_Location = rep(10, 10), geo_accession = rep(10, 10), KRAS_Mutation = rep(10, 10), extra = rep(10, 10), extra2 = rep(10, 10))

  result <- rename_and_select(fake_tibble)

  expected_names <- c("Sex", "Age", "TNM_Stage", "Tumor_Location", "geo_accession", "KRAS_Mutation", "Subtype", "Batch")
  expect_identical(names(result), expected_names, info = "The returned tibble does not have the expected column names. Ensure your function renames and selects columns correctly.")
})




test_that("stage_as_factor adds a Stage column with the correct format", {

  # Define the fake_tibble inside this test function
  fake_tibble <- tibble::tibble(
    Age_at_diagnosis = rep(10, 10),
    SixSubtypesClassification = rep('C4', 10),
    normalizationcombatbatch = rep('batch1', 10),
    Sex = rep(10, 10),
    TNM_Stage = rep(10, 10), # We're using the same value 10 for this example. Ideally, this would be different stages.
    Tumor_Location = rep(10, 10),
    geo_accession = rep(10, 10),
    KRAS_Mutation = rep(10, 10),
    extra = rep(10, 10),
    extra2 = rep(10, 10)
  )

  result <- stage_as_factor(fake_tibble)

  expect_true("Stage" %in% names(result), info = "Stage column is missing from the returned tibble. Ensure your function adds this column.")

  expect_true(is.factor(result$Stage), info = "Stage column should be of type factor. Make sure you've set it as such.")

  expect_true(all(grepl("^stage ", result$Stage)), info = "Entries in the Stage column do not follow the 'stage x' format. Double-check the string manipulation in your function.")

  # Since our fake tibble has 10 as the TNM_Stage for all rows, let's also check if the Stage is correctly formatted as "stage 10"
  expect_true(all(result$Stage == "stage 10"), info = "Stage column does not have the expected values. Check the transformation in your function.")
})

test_that("mean_age_by_sex calculates correct mean age for given sex", {
  # Create dummy data
  set.seed(123)  # Setting seed for reproducibility
  data <- tibble(
    Sex = sample(c("M", "F"), 100, replace = TRUE),
    Age = sample(15:90, 100, replace = TRUE)
  )

  for (sex in c("M", "F")) {
    mean_age <- data %>%
      filter(Sex == sex) %>%
      summarize(mean_age = mean(Age, na.rm = TRUE)) %>%
      pull(mean_age)

    result <- mean_age_by_sex(data, sex)

    expect_equal(as.double(result), as.double(mean_age),
                     info = paste("The calculated mean age for", sex, "does not match the expected value. Check the aggregation logic in your function."))
  }
})



test_that("age_by_stage calculates average age per stage correctly", {
  # Create mock test data
  test_data <- tibble(
    Stage = c(rep('stage 1', 4), rep('stage 2', 3), rep('stage 3', 2), rep('stage 4', 1)),
    Age = c(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)
  )

  # Calculate expected result
  expected_result <- test_data %>%
    group_by(Stage) %>%
    summarize(mean_avg = mean(Age)) %>%
    ungroup() %>%
    arrange(Stage)

  # Apply function on mock test data
  calculated_result <- age_by_stage(test_data)

  # Ensure calculated_result has the required columns
  expect_true(all(c("Stage", "mean_avg") %in% colnames(calculated_result)),
              info = "The result should have 'Stage' and 'mean_avg' columns.")

  # Compare the calculated result to the expected result
  expect_equivalent(
    calculated_result,
    expected_result,
    info = "The calculated result does not match the expected result. Please check the function logic."
  )
})



# Test for subtype_stage_cross_tab function
test_that("subtype_stage_cross_tab returns correct cross-tabulated table", {
  # Create test data
  test_stage <- tibble(Stage = c(rep('stage 4', 2), rep('stage 3', 1), rep('stage 1', 4), rep('stage 2', 3)),
                       Subtype = c(rep('C4', 3), rep('C3', 7)))

  # Create expected cross-tabulated result
  test_crosstab <- tribble(
    ~Stage, ~C3, ~C4,
    "stage 1", 4, 0,
    "stage 2", 3, 0,
    "stage 3", 0, 1,
    "stage 4", 0, 2
  )

  # Calculate cross-tabulated result using the function
  calculated_crosstab <- subtype_stage_cross_tab(test_stage)

  # Ensure that there are no NA values in the result
  expect_false(any(is.na(calculated_crosstab)), info = "Your table should not have any NA values. Fill missing pairs with zeros.")

  # Ensure 'Stage' is present as a column
  expect_true("Stage" %in% colnames(calculated_crosstab), info = "'Stage' should be present as a column in the calculated table.")

  # Check if the calculated result matches the expected result
  expect_equal(test_crosstab, calculated_crosstab, info = "The calculated cross-tabulated table does not match the expected result.")
})


# Test for summarize_expression function
test_that("summarize_expression behaves as expected", {

  # Create a synthetic tibble
  exprs <- tibble(
    subject_id = c("A", "B"),
    probe1 = c(1, 2),
    probe2 = c(2, 7),
    probe3 = c(3, 15)
  )


  result <- summarize_expression(exprs)

  # Check that result is a tibble
  expect_true(is_tibble(result), info = "The returned result is not a tibble. Ensure your function returns a tibble.")

  # Check that result has the expected columns
  expect_identical(names(result), c("mean_exp", "variance", "probe"), info = "The returned tibble does not have the expected column names.")

  # Check number of rows
  expect_equal(nrow(result), 3, info = "Mismatch in the number of rows of the result compared to the number of columns in the expression matrix.")

  # Check values in mean_exp column
  values_mean <- c(1.5, 4.5, 9)
  names_mean <- c('probe1', 'probe2', 'probe3')
  test_means <- setNames(values_mean, names_mean)
  expect_identical(result %>% pull(mean_exp), test_means, info = "Mismatch between the expected means and the means in the mean_exp column.")

  # Check values in variance column
  values_var <- c(0.5, 12.5, 72)
  names_var <- c('probe1', 'probe2', 'probe3')
  test_var <- setNames(values_var, names_var)
  expect_identical(result %>% pull(variance), test_var, info = "Mismatch between the expected variances and the variances in the variance column.")

})




