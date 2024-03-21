# initialize.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# The idea of this script is to initialize folder structure to store artefacts made
# in this project and to load necessary R packages.



# Script structure
# >>>>>>>>>>>>>>>>
# 1. Directory settings: Directory structure for result storage.
# 2. Library: Loading of depended R packages and local functions
# 3. Summary of the purpose of each scripts and their ordinal positions




# 1. Directory settings: Directory structure for result storage.
# ==============================================================================

# Output directory
if(!dir.exists("results/tables")){
  dir.create(path = "results/tables", recursive = TRUE, mode = "0700")
}
if(!dir.exists("results/figures")){
  dir.create(path = "results/figures", recursive = TRUE, mode = "0700")
}
if(!dir.exists("results/data")){
  dir.create(path = "results/data", recursive = TRUE, mode = "0700")
}
if(!dir.exists("results/documents")){
  dir.create(path = "results/documents", recursive = TRUE, mode = "0700")
}


out_tables <- "results/tables/"
out_figures <- "results/figures/"
out_data <- "results/data/"
out_documents <- "results/documents/"


#
# ==============================================================================




# 2. Library: Loading of depended R packages and local functions
# ==============================================================================

# Public

library(tidyverse)
library(sva) # ComBat
library(pROC) # roc, auc
library(irr) # kappa2: Cohen's kappa
library(hablar) # retype: Transforms all elements into simple classes
library(DescTools) # RobScale: Robust Scaling With Median and Mad
library(rmarkdown) # render: To render the Rmarkdown file into specific output format
library(curatedBreastData) # data("clinicalData"): to comapre with present study
library(xlsx) # write.xlsx wrirte xlsx files in multiple sheets

# library(openxlsx) #
# library(writexl) # write_xlsx: wrirte xlsx files in multiple sheets
# library(readxl) # read_excel: Read xls and xlsx files
# library(altmeta) # metahet: Meta-Analysis Heterogeneity Measures


# Private

source("code/r/functions.R")


#
# ==============================================================================




# 3. Summary of the purpose of each scripts and their ordinal positions
# ==============================================================================

# Script structure:
# 1. initialize.R
#       a. Initialize folder structure to store artefacts fo the analysis
#       b. Load necessary R packages
# 2. functions.R
#       a. Supporting functions for other scripts.
# 3. s1_geo_search_result_curation.R
#       a. Compile the GEO search results.
# 4. s2_geo_dataset_extraction.R
#       a. Extract the sample characteristics and expression matrix from the series matrix.
#       b. Format sample characteristics both computationally and manually.
# 5. s3_sample_characteristic_curation.R
#       a. Manual curation and annotation of sample characteristics names, and re-coding of characteristics values for consistency.
# 6. s4_expression_data_integration.R
#       a. NA filtering, converting to log2 space, maximum probeset collapse, and data integration.
# 7. s5_clinical_data_integration.R
#       a. Selection of relevant clinical-pathological variables, their integration, and the curation of treatment regimens.
# 8. s6_process_and_validate_integrated_expression_matrix.R
#       a. Comparison of methods to correct dataset effect from the pooled gene expression data
#           i. Non-corrected (log2)
#           ii. Combat (log2 + combat)
#           iii. Quantile (log2 + per dataset quantile gene-rescaling)
#           iv. Madquantile (log2 + MAD sample-rescaling + per dataset quantile gene-rescaling)
#       b. The comparison is made by comparing
#           i. RLE plots
#           ii. IHC vs. PAM50 molecular subtype agreement (Cohen’s Kappa)
#           iii. The ability of  the ER/HER2 gene to predict the respective IHC status (AUC).
# 9. s7_figures_tables_data.R
#       a. Generate journal-specific figures and tables.


#
# ==============================================================================

