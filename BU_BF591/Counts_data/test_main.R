#!/usr/bin/Rscript
source("main.R")
library(testthat)

test_data <- as_tibble(read.csv("data/test.tsv", sep=""))
td2 <- test_data[1:1000,]

####Testing read_data()####
test_that("Function read_data() loads properly", {
  data <- read_data("data/verse_counts.tsv")
  type <- unlist(data %>% summarise_all(class))
  
  #tests dimensions
  expect_equal(dim(data), c(55416, 9)) 
  
  #tests to ensure column names
  expect_equal(colnames(data), c("gene",  "vP0_1", "vP0_2", "vP4_1", "vP4_2", "vP7_1", "vP7_2", "vAd_1", "vAd_2")) 
  
  #tests to ensure that gene is a character type
  expect_equal(type[['gene']], "character")
  
  #tests to ensure that the other columns are all numeric
  for(name in names(type)[names(type)!='gene']){
    expect_setequal(type[[name]], "numeric")
  }
  
  expect_false(has_rownames(data), "Your read_data() function returns something that has rownames. You should not be returning anything with rownames as tibbles should not have them")
  
  expect_true(is_tibble(data), "Your read_data() function does not return a tibble. Please make sure that it does")
})

####Testing filter_zero_var_genes()####
test_that("Function filter_zero_var_genes() filters properly", {

  #expected rows and columns
  col <- length(colnames(test_data))
  row <- 999

  #run function
  zero_var_filtered <- filter_zero_var_genes(test_data)
  #check for tibble
  expect(is_tibble(zero_var_filtered), "filter_zero_var_genes() does not output a tibble")

  if(length(zero_var_filtered)==row){
    #Notify students that they may have transposed
    fail("Check to make sure you didn't transpose where you shouldn't've; you have way too many columns")
  }else if(isFALSE(is.list(zero_var_filtered))){
    fail("Data is in the wrong format or NULL. Check your function's output")
  }else{
  #checks for dimensions
  dims <- dim(zero_var_filtered)
  expect(dims[1]==row, paste0("Number of rows is incorrect nrow(user_function) - nrow(expected) == ", dims[1]-row))
  expect(dims[2]==col, paste0("Number of columns is incorrect ncol(user_function) - ncol(expected) == ", dims[2]-col))

  #checks for 'gene' column as column 1
  expect(colnames(zero_var_filtered)[1] == "gene", paste0("filter_zero_var_genes()'s output should have a first column named 'gene.' Returned first column named: ", colnames(zero_var_filtered)[1]))

  #checks for colname order
  expect_named(zero_var_filtered,  names(test_data))

  #check if expected genes were removed
  user_removed <- setdiff(test_data$gene, zero_var_filtered$gene)
  expected_removed <- c("CJNN53716", "WHGN65413")
  expect(all(user_removed %in% expected_removed), paste0(length(setdiff(user_removed, expected_removed)), " extra gene(s) were filtered out"))
  expect(all(expected_removed %in% user_removed), paste0(length(setdiff(expected_removed, user_removed)), " gene(s) should have been filtered out but weren't"))


  #Makes sure output has no rownames
  expect(isFALSE(has_rownames(zero_var_filtered)), "Your output should not have row names but it does")
  }

  #Make sure all data in a row is the same
  expect(isTRUE(all_equal(test_data[!(test_data$gene %in% user_removed),], zero_var_filtered)),
         "If you see this and your filter_zero_var_genes() returns all of its columns, then your gene data is somehow being shuffled;

Gene names in your output are not correlating with their true counts")


})

####Testing timepoint_from_sample()####
test_that("timepoint_from_sample() returns correctly",{

  #test for base case
  expect_equal(timepoint_from_sample("vAd_1"), "Ad")

  #test for letter/number variance
  expect_equal(timepoint_from_sample("vP7_2"), "P7")

  #Cases outside of the sample data but fit within the requirements as defined in my outline
  expect_equal(timepoint_from_sample("vP7_8"), "P7") # last number not a 1 or 2
  expect_equal(timepoint_from_sample("vBd_1"), "Bd") # Capital letter not A or P
  expect_equal(timepoint_from_sample("vBd_7"), "Bd") # Both the above cases
  expect_equal(timepoint_from_sample("vJ1_9"), "J1") # 1 as 3rd character
  expect_equal(timepoint_from_sample("vJ1_1"), "J1") # 1 as 3rd char w/ a second 1 as replicate number


})

####Testing sample_replicate()####
test_that("sample_replicate() output is correct", {

  #basic instance, checks for a string
  expect_equal(sample_replicate("vAd_1"), "1")

  #Ensures that the correct number get returned (ignores the "4")
  expect_equal(sample_replicate("vP4_1"), "1")

  #Number variations
  expect_equal(sample_replicate("vP3_2"), "2")
  expect_equal(sample_replicate("vG8_0"), "0")
  expect_equal(sample_replicate("vPk_8"), "8")

  # #Cases that don't have to work, correct functions may return differently depending on method
  # sample_replicate("vPk_82") # ("82", "8", "2", NULL, ...) # multiple numbers after the `_`
  # sample_replicate("v8k_8g") # ("8g", "8", "g", NA, ...) # It's okay if there are letters after `_`
  # sample_replicate("vLj_1_19") # ("1", "19", "9", ... ) # Returning "1_19" instead could be valid, just use `expect_equal(sample_replicate("vLj_1_19"), "1_19")`
  # sample_replicate("v1_An") # ("An", "A", NULL, ...) # doesn't have to be numbers after the `_` at all
  # sample_replicate("v_<3") #("<3", "3", NULL, ...) #non alphanumerical for kicks

})

####Testing meta_info_from_labels()####
test_that("meta_info_from_labels() is in the correct format", {
  #t_data_gene <- c("vK3_1", "vMr_4", "vAd_2", "vLj_9")

  desired_output <- tibble(sample=c("vK3_2", "vMr_1","vK3_1", "vMr_2"), timepoint=c("K3", "Mr", "K3", "Mr"), replicate=c("2","1","1","2"))

  t_tib <- meta_info_from_labels(colnames(test_data)[-1])%>%
    mutate(across(everything(), as.character)) #Did not have a test case that required this but going to do it anyways

  #If someone added named vectors to the tibble, then they would fail some of these tests. Auto unnaming columns for them
  for(i in 1:ncol(t_tib)){
    t_tib[[i]] <- unname(t_tib[[i]])
  }#don't need to do this for rows

  #check if output is a tibble
  expect(is_tibble(t_tib), "meta_info_from_labels() does not output a tibble")

  #check tibble dimensions
  row <- nrow(desired_output)
  col <- ncol(desired_output)
  dims <- dim(t_tib)
  expect(dims[1]==row, paste0("Number of rows is incorrect nrow(user_function) - nrow(expected) == ", dims[1]-row))
  expect(dims[2]==col, paste0("Number of columns is incorrect ncol(user_function) - ncol(expected) == ", dims[2]-col))


  #check column names
  expect_named(t_tib, c('sample', 'timepoint', 'replicate'))

  #check row order
  expect_equal(list(t_tib$sample), list(desired_output$sample))

  #check if timepoint and replicate points are aligned with their proper samples
  expect(all(t(t_tib) %in% t(desired_output)), "Data is not consistent by rows. If you have the correct columns, then there might be some data shuffling going on in one or more of your functions")

})

####Testing get_library_size()####
test_that("get_library_size() sums up sample reads correctly", {

  expected <- c("vK3_1" = 254130, "vMr_2" = 284628,"vK3_2" = 267519, "vMr_1" = 283515 )

   lib_size <- get_library_size(test_data)

  if(is_tibble(lib_size)){
    if(ncol(lib_size)<=2){ #long tibble
    #Test for exact
    lib_size <- deframe(lib_size)

    }else if(dim(lib_size)[1]==1){ #wide tibble
      #Test for exact
      lib_size <- unlist(lib_size)
    }else {
      fail("Tibble output has too many rows and/or columns")
    }
  }

  #check for length
  expect(length(lib_size)==length(expected), paste0("Unexpected number of samples returned. Expected 4, got", length(lib_size)))
  #check for names
  expect_named(lib_size, names(expected), ignore.order = TRUE)

  #check for name values
  for(name in names(expected)){
    expect(lib_size[[name]]==expected[[name]], paste0("Library sum value for ", name, " is incorrect"))
  }

})

####Testing normalize_by_cpm####
test_that("normalize_by_cpm() has the correct output",{

  cpm <- normalize_by_cpm(td2)%>%
    mutate(across(where(is.numeric), round, 3))

  sample_answers <-tibble(gene = c("SZUG00090", "JIZU74666", "AWDL96666", "HDLU58875"),
                          vK3_2 = c(44.8566270059323, 59.8088360079097, 59.8088360079097, 112.141567514831),
                          vMr_1 = c(59.9615540623953, 56.4344038234309, 42.3258028675731, 88.1787559741107),
                          vK3_1 = c(55.0899146106324, 59.0249085113918, 39.3499390075945, 78.6998780151891),
                          vMr_2 = c(49.1870090082494, 59.7270823671599, 31.6202200767317, 91.3473024438917))%>%
    mutate(across(where(is.numeric), round, 3))
  #Testing for tibble
  expect(is_tibble(cpm), "normalize_by_cpm() does not output a tibble")

  #Testing dimensions
  dims <- dim(cpm)
  row <- nrow(td2)
  col <- ncol(td2)
  expect(dims[1]==row, paste0("Number of rows returned by `normalize_by_cpm()` is incorrect; nrow(user_function) - nrow(expected) == ", dims[1]-row))
  expect(dims[2]==col, paste0("Number of columns returned by `normalize_by_cpm()` is incorrect; ncol(user_function) - ncol(expected) == ", dims[2]-col))

  #Testing colnames, must be in the same order
  expect_named(cpm, colnames(td2))

  #testing for rownames
  expect(isFALSE(has_rownames(cpm)), "Your normalize_by_cpm() output has rownames, do not put rownames on your tibbles")

  #testing for genes
  expect(all(td2$gene %in% cpm$gene), "Somehow your normalize_by_cpm() output is not returning gene names properly or is filtering out genes")

  #Testing the whole thing, checking for equation execution and row consistency
  for(g in sample_answers$gene){
    expect_equal(cpm[cpm$gene==g,], sample_answers[sample_answers$gene==g,], .03)
  }
})

####Testing deseq_normalize####
test_that("deseq_normalize() has correct output",{
  meta_data <- tibble(sample=c("vK3_2", "vMr_1","vK3_1", "vMr_2"), timepoint=c("K3", "Mr", "K3", "Mr"), replicate=c("2","1","1","2"))

  normalized <- suppressWarnings(deseq_normalize(td2, meta_data) %>%
    mutate(across(where(is.numeric), round, 3)), class='warning' )

  sample_answers <-tibble(gene = c("DYOG91349", "MCAB52159", "QUIS14366", "DRWM24408"),
                          vK3_2 =  c(23.1379797208499, 13.0779885378717, 25.1499779574455, 456.723599707211),
                          vMr_1 = c(16.5231110249736, 9.71947707351388, 30.130378927893, 489.861644505099),
                          vK3_1 = c(17.4371217965184, 9.23141742168622, 25.6428261713506, 413.362357882172),
                          vMr_2 =  c(15.6092727506343, 15.6092727506343, 40.9743409704149, 458.522387049882)) %>%
    mutate(across(where(is.numeric), round, 3))

  #Testing for tibble
  expect(is_tibble(normalized), "deseq_normalized() does not output a tibble")

  #Testing dimensions
  dims <- dim(normalized)
  row <- nrow(td2)
  col <- ncol(td2)
  expect(dims[1]==row, paste0("Number of rows returned by `deseq_normalized()` is incorrect; nrow(user_function) - nrow(expected) == ", dims[1]-row))
  expect(dims[2]==col, paste0("Number of columns returned by `deseq_normalized()` is incorrect; ncol(user_function) - ncol(expected) == ", dims[2]-col))

  #Testing colnames, must be in the same order
  expect_named(normalized, colnames(td2))

  #testing for rownames
  expect(isFALSE(has_rownames(normalized)), "Your deseq_normalized() output has rownames, do not put rownames on your tibbles")

  #testing for genes
  expect(all(td2$gene %in% normalized$gene), "Somehow your deseq_normalized() output is not returning gene names properly or is filtering out genes")

  #Testing the whole thing, checking for equation execution and row consistency
  for(g in sample_answers$gene){
     expect_equal(normalized[normalized$gene==g,], sample_answers[sample_answers$gene==g,], .03)
  }

})

####Testing plot_sample_distributions() ####
test_that("plot_sample_distributions is performing correctly()",{
  test_title <- "A Sample Title"
  row <- (ncol(td2)-1) * nrow(td2)
  col <- 2
  ggplot <- plot_sample_distributions(td2[-1], FALSE,test_title)
  dims <- dim(ggplot$data)
  #testing that the data got used

  expect(dims[1]==row, paste0("Number of rows by plot_sample_distributions() is incorrect; nrow(user_function) - nrow(expected) == ", dims[1]-row))
  expect(dims[2]==col, paste0("Number of columns by plot_sample_distributions() is incorrect; ncol(user_function) - ncol(expected) == ", dims[2]-col))


  #testing the data itself
  data <- ggplot$plot_env$data
  for(n in colnames(td2)[-1]){
    expect(all(td2[[n]][order(td2[[n]])]==data[[n]][order(data[[n]])]), paste0("Column ", n, "had incorrect data transferred"))
  }


  #Testing data log10 transformation
  ggplot_log <- plot_sample_distributions(td2[-1], TRUE, "Title Here")
  if(!ggplot_log$plot_env$scale_y_axis){
    log_test <- td2[-1] %>%
    mutate(across(where(is.numeric), log10))
  }else{
    log_test <-td2[-1]
  }

  log_data <- ggplot_log$plot_env$data

  for(n in colnames(td2)[-1]){
    expect(all(log_test[[n]][order(log_test[[n]])]==log_data[[n]][order(log_data[[n]]) ] ), paste0("Column ", n, " has incorrect data processed for log10() transformation"))
  }
})

####Testing plot_variance_vs_mean() ####
test_that("that plot_variance_vs_mean is working properly",{
  test_title <- "A Sample Title"

  row <-  nrow(td2)
  col <- 3
  ggplot <- plot_variance_vs_mean(td2[-1], FALSE,test_title)
  dims <- dim(ggplot$data)

  #testing that the data got used
  expect(dims[1]==row, paste0("Number of rows by plot_variance_vs_mean() is incorrect; nrow(user_function) - nrow(expected) == ", dims[1]-row))
  expect(dims[2]==col, paste0("Number of columns by plot_variance_vs_mean() is incorrect; ncol(user_function) - ncol(expected) == ", dims[2]-col))

  #testing the data itself
  data <- ggplot$plot_env$data
  for(n in colnames(td2)[-1]){
    expect(all(td2[[n]][order(td2[[n]])]==data[[n]][order(data[[n]])]), paste0("Column ", n, "had incorrect data transferred"))
  }

  #Testing data log10 transformation
  ggplot_log <- plot_variance_vs_mean(td2[-1], TRUE, "Title Here")
  if(!ggplot_log$plot_env$scale_y_axis){
    log_test <- td2[-1] %>%
      mutate(across(where(is.numeric), log10))
  }else{
    log_test <-td2[-1]
  }

  log_data <- ggplot_log$plot_env$data

  for(n in colnames(td2)[-1]){
    expect(all(log_test[[n]][order(log_test[[n]])]==log_data[[n]][order(log_data[[n]]) ] ), paste0("Column ", n, " has incorrect data processed for log10() transformation"))
  }
})


####Testing plot_pca() ####
test_that("plot_pca",{
  test_title <- "A Sample Title"

  meta_data <- tibble::tibble(sample=c("vK3_2", "vMr_1","vK3_1", "vMr_2"), timepoint=c("K3", "Mr", "K3", "Mr"), replicate=c("2","1","1","2"))
  row <- 4
  col <- 5
  ggplot <- plot_pca(td2[-1], meta_data ,test_title)
  dims <- dim(ggplot$data)

  #testing that the data got used
  expect(dims[1]==row, paste0("Number of rows by plot_pca() is incorrect; nrow(user_function) - nrow(expected) == ", dims[1]-row))
  expect(dims[2]==col, paste0("Number of columns by plot_pca() is incorrect; ncol(user_function) - ncol(expected) == ", dims[2]-col))

  #Testing the PCA components
  pca <- prcomp(t(td2[-1]))
  expect(all(ggplot$data$PC1==pca$x[,1]), "Expected PC1 data not found. Did you use the correct PC?")
  expect(all(ggplot$data$PC2==pca$x[,2]), "Expected PC2 data not found. Did you use the correct PC?")

  #testing the data itself
  data <- ggplot$plot_env$data
  for(n in colnames(td2)[-1]){
    expect(all(td2[[n]][order(td2[[n]])]==data[[n]][order(data[[n]])]), paste0("Column ", n, "had incorrect data transferred"))
  }
})

