

rm(list=ls())
# Install the released version from CRAN using
if (!requireNamespace("multiclassPairs", quietly = TRUE)) {
  install.packages("multiclassPairs")
}

# Or install the dev version from GitHub using
# if (!requireNamespace("multiclassPairs", quietly = TRUE)) {
#  if (!requireNamespace("devtools", quietly = TRUE)) {
#    install.packages("devtools")
#  }
#  library(devtools) # this package is needed to install from GitHub
#  install_github("NourMarzouka/multiclassPairs", build_vignettes = TRUE)
#}

# Install the dependencies from Bioconductor
# BiocManager, Biobase, and switchBox packages from Bioconductor are needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("Biobase", quietly = TRUE)) {
  BiocManager::install("Biobase")
}
if (!requireNamespace("switchBox", quietly = TRUE)) {
  BiocManager::install("switchBox")
}

# load multiclassPairs library
library(multiclassPairs)




df2<-read.csv('H:/FYP/summer/Data/TOD_CUP/MAD_5000_Label',sep='\t',header=TRUE,row.names=1)#label
df1<-read.csv("H:/FYP/summer/Data/TOD_CUP/MAD_5000",sep='\t',header=TRUE,row.names=1)
n<-ncol(df1)
set.seed(1234)
training_samples<-sample(1:n,size=n*0.7)
train<-df1[,training_samples]
test<-df1[,-training_samples]
L1<-df2[nrow(df2),training_samples]
L1<-as.character(L1)
L2<-df2[nrow(df2),-training_samples]
L2<-as.character(L2)
object<-ReadData(Data=train,
                 Labels=L1,
                 Platform=NULL,
                 verbose=FALSE)

object
filtered_genes <- filter_genes_TSP(data_object = object,
                                   filter = "one_vs_one",
                                   platform_wise = FALSE,
                                   featureNo = 5000,
                                   UpDown = TRUE,
                                   verbose = TRUE)
filtered_genes