library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
# library(purrr)


# ----------------------- Helper Functions to Implement ------------------------

#' Read the expression data "csv" file.
#'
#' Function to read microarray expression data stored in a csv file. The
#' function should return a sample x gene tibble, with an extra column named
#' "subject_id" that contains the geo accession ids for each subject.
#'
#' @param filename (str): the file to read.
#'
#' @return
#' @export
#'
#' @examples expr_mat <- read_expression_table('example_intensity_data_subset.csv')
read_expression_table <- function(filename) {
  expre<-readr::read_csv(filename,col_names=FALSE)
  expre_1<-t(expre)
  expre_2<-as_tibble(expre_1)
  colnames(expre_2)<-expre_2[1,]
  colnames(expre_2)[1]<-'subject_id'
  expre_2<-expre_2[2:nrow(expre_2),]
  return (expre_2)
  
  
  
  # expr_mat<-read.table(
  #   filename,
  #   header=TRUE,
  #   sep=' '
  # )
  #   
  # expre_mat<-tibble::as_tibble(
  #   t(expr_mat),
  #   rownames='subject_id'
  # )

  
  #use pivot_longer and pivot_wider
}



#' Load Metadata from Specified CSV File
#'
#' This function reads the provided CSV file into a dataframe.
#'
#' @param filepath (character) The path to the CSV file.(data/proj_metadata.csv)
#'
#' @return A dataframe containing the loaded metadata.
#'

load_metadata <- function(filepath) {
  
  # TODO: Use appropriate function to read in the CSV file
  #metadata <-   
    
    # Return the loaded metadata
    #return(metadata)
  meta<-readr::read_csv(filepath)
  return(meta)
}






#' Replaces all '.' in a string with '_'
#'
#' @param str String to operate upon.
#'
#' @return reformatted string.
#' @export
#'
#' @examples
#' period_to_underscore("foo.bar")
#' "foo_bar"
period_to_underscore <- function(str) {
  new_string<-str_replace_all(str,'\\.','_')
  #[.] also can be used to replace \\.
  return (new_string)
}


# rename variables:
# Age_at_diagnosis to Age
# SixSubtypesClassification to Subtype
# normalizationcombatbatch to Batch

#' Rename and select specified columns.
#'
#' Function to rename Age_at_diagnosis, SixSubtypesClassification, and
#' normalizationcombatbatch columns to Age, Subtype, and Batch, respectively. A
#' subset of the data should be returned, containing only the Sex, Age, TNM_Stage,
#' Tumor_Location, geo_accession, KRAS_Mutation, Subtype, and Batch columns.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) renamed and subsetted tibble
#' @export
#'
#' @examples rename_and_select(metadata)
#' 
#' 
rename_and_select <- function(data) {
  renamed_tibble<-rename(data,
            Age=Age_at_diagnosis  ,
           Subtype=SixSubtypesClassification  ,
            Batch=normalizationcombatbatch )%>%
  dplyr::select(Sex, Age, TNM_Stage, Tumor_Location,
                                 geo_accession, KRAS_Mutation, Subtype, Batch)%>%
    
  return (selected_tibble)
}


#' Create new "Stage" column containing "stage " prefix.
#'
#' Creates a new column "Stage" with elements following a "stage x" format, where
#' `x` is the cancer stage data held in the existing TNM_Stage column. Stage
#' should have a factor data type.
#'
#' @param data  (tibble) metadata information for each sample
#'
#' @return (tibble) updated metadata with "Stage" column
#' @export
#'
#' @examples metadata <- stage_as_factor(metadata)
stage_as_factor <- function(data) {
  dt<-mutate(data,Stage=factor(paste0('stage ',TNM_Stage)))
  return (dt)
}


#' Calculate age of samples from a specified sex.
#'
#' @param data (tibble) metadata information for each sample
#' @param sex (str) which sex to calculate mean age. Possible values are "M"
#' and "F"
#'
#' @return (float) mean age of specified samples
#' @export
#'
#' @examples mean_age_by_sex(metadata, "F")
mean_age_by_sex <- function(data, sex) {
  rst<-mean(data$Age[data$Sex==sex])
  return (rst)
}


#' Calculate average age of samples within each cancer stage. Stages should be
#' from the newly created "Stage" column.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) summarized tibble containing average age for all samples from
#' each stage. Name the newly created column containing the average, 'mean_avg'
#' @export
#'
#' @examples age_by_stage(data)
age_by_stage <- function(data) {
  summary_tb<-dplyr::group_by(data,
                              Stage)%>%
    # dplyr::summarize('\'(mean(Age))\''=round(mean(Age),1))
    dplyr::summarize(mean_avg=round(mean(Age),1))
  return (summary_tb)
}

#' Create a cross tabulated table for Subtype and Stage using dplyr methods.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) table where rows are the cancer stage of each sample, and the
#' columns are each cancer subtype. Elements represent the number of samples from
#' the corresponding stage and subtype. If no instances of a specific pair are
#' observed, a zero entry is expected.
#' @export
#'
#' @examples cross_tab <- dplyr_cross_tab(metadata)
subtype_stage_cross_tab <- function(data) {
  summary_tb<-dplyr::group_by(data,
                              Stage)%>%
    dplyr::summarize('C3'=sum(Subtype=="C3"),"C4"=sum(Subtype=="C4"))
  return (summary_tb)
}

#' Summarize average expression and probe variability over expression matrix.
#'
#' @param exprs An (n x p) expression matrix, where n is the number of samples,
#' and p is the number of probes.
#'
#' @return A summarized tibble containing `main_exp`, `variance`, and `probe`
#' columns documenting average expression, probe variability, and probe ids,
#' respectively.

summarize_expression <- function(exprs) {
  
  exprs<-exprs[,2:ncol(exprs)]%>%
    mutate_all(as.numeric)
  
  new_tb<-tibble()
  
  for (colnum in 1:ncol(exprs)){
    new_tb[colnum,'mean_exp']=colMeans(exprs[,colnum])
    new_tb[colnum,'variance']=var(exprs[,colnum])
    new_tb[colnum,'probe']=colnames(exprs)[colnum]
  }
  
  return (new_tb)
}
