# main.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# The idea of this main script is to manage all other scrips from here.
# However, in this data pooling project no script can run without manual
# intervention. Scripts which needs manual interventions are prefixed with
# keyword "s*_manual". Hence in effect this main script functions as a place
# holder where in the details of other sripts and how they are connected are
# documented.


# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Directory settings: Directory structure for result storage.
# 2. Library: Loading of depended R packages and local functions
# 3. Summary of the pupose of each scripts and their ordinal positions
# 4. Details of each scripts and the main Robjects/results created within




# 1. Directory settings: Directory structure for result storage.
# ==============================================================================

# Output directory
if(!dir.exists("results/main")){
  dir.create(path = "results/main", recursive = TRUE, mode = "0700")
}
outdir <- "results/main/"

# # Tmp directory
# if(!dir.exists("results/tmp")){
#   dir.create(path = "results/tmp", recursive = TRUE, mode = "0700")
# }
# tmpdir <- "results/tmp/"

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
# library(readxl) # read_excel: Read xls and xlsx files
# library(altmeta) # metahet: Meta-Analysis Heterogeneity Measures


# Private
source("R/functions.R")

#
# ==============================================================================




# 3. Summary of the pupose of each scripts and their ordinal positions
# ==============================================================================

  # Note:
  # The keyword "manual" in script filename implies that some part of the script
  # requires manual cleaning/formatting/curation etc.

  # main.R:
  #       This script !!!!!!!
  #       Main script from which other scripts are managed.
  # functions.R:
  #       Supporting functions for main and other scripts.
  #       Sourced under "Library" section.
  # s1_manual_geo_search_result_curation.R:
  #       Compile geo search results.
  # s2_manual_geo_dataset_extraction.R:
  #       Extract series, sample and expression data from each series matrix.
  #       Format sample characteristics both computationally and manually.
  # s3_manual_sample_characterisitic_curation.R:
  #       Manual curation and annotation of sample characteristics names and
  #         recodeing of characterisitcs value for consistancy.
  # s4_manual_expression_data_integration.R:
  #       NA filtering, convert to log2 space, max-var collapsing, data integration.
  # s5_manual_clinical_data_integration.R:
  #       Selecting relevant clincal-pahological variables, integration,
  #         treatment regimen curation.
  # s6_process_and_validate_integrated_expression_matrix.R:
  #       Dataset effect corretion method comparison with respect to
  #       1) Relative Log Expression (RLE) plots
  #       2) IHC vs PAM50 molecular subtype agreement using Cohen's Kappa coeffcient
  #       3) Ability of ER/HER2 gene in predicting respective IHC status
  #         measured using AUC (Area under the curve of ROC analysis).
  #       The dataset effect correction methods compared were:
  #       1) Non-corrected (original log2)
  #       2) Combat (log2 + combat)
  #       3) Quantile (log2 + per dataset quantile gene scaling)
  #       4) Madquantile (log2 + MAD sample scaling + per dataset quantile gene scaling)
  # s7_figures_tables_data.R:
  #       Generate journal specific figures, tables and dataset.

#
# ==============================================================================




# 4. Details of each scripts and the main Robjects/results created within
# ==============================================================================


# s1_manual_geo_search_result_curation.R
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Do not run !!!
# To run, go to each script and run code chunk by chunk !!!
# source("R/s1_manual_geo_search_result_curation.R")

# Robject created: 1) geo_series_summary
# Relevant Robjects: 1) geo_series_summary




# s2_manual_geo_dataset_extraction.R
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Do not run !!!
# To run, go to each script and run code chunk by chunk !!!
# source("R/s2_manual_geo_dataset_extraction.R")

# Robject created: 1) geo (44 series matrices, 40 datasets)
# Relevant Robjects: 1) geo_series_summary, 2) geo





# s3_manual_sample_characterisitic_curation.R
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Do not run !!!
# To run, go to each script and run code chunk by chunk !!!
# source("R/s3_manual_sample_characterisitic_curation.R")

# Robject updated: geo ( 1) 5 series matrices filtered out, which gives final
#                           39 series matrices, 36 datasets.
#                        2) Updated with curated and annotated clinical data.)
# Relevant Robjects: 1) geo_series_summary, 2) geo




# s4_manual_expression_data_integration.R
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Do not run !!!
# To run, go to each script and run code chunk by chunk !!!
# source("R/s4_manual_expression_data_integration.R")

# Robject updated: geo ( 1) Removed GSE4056 due to large no.of NA expression,
#                           which gives final 38 series matrices, 35 datasets.
#                        2) Removed 12 samples from GSE69031.
#                        3) Removed NA genes from 10 datasets. )
# Robject created: gpl (Cleaned gpl annotation list.)
#                  geo_maxvar_annot (Max-var annotation list per series matrix in geo.)
#                  geo_tidy (Simialr to geo, but hold only cleaned max-var
#                             collapsed expression and cleaned clinical data.
#                             Cleaned clinical data will be interated in the following
#                             scirpt. Note that geo_tidy also contains
#                             38 series matrices, 35 datasets)
#                  geo_expr (Integrated gene X sample expression matrix, non-corrected.
#                           Note that 9 additional seris matrices were filtered out
#                           from geo_tidy due to poor genome representation of arrays
#                           or expression missingness. Hence, geo_expr contains only
#                           29 sereis matrices, 28 datasets.)
#                  geo_expr_meta (Series matrix -  sample geo accession mapping)
# Relevant Robjects:  1) geo_series_summary, 2) geo, 3) gpl, 4) geo_maxvar_annot
#                     5) geo_tidy, 6) geo_expr, 7) geo_expr_meta,




# s5_manual_clinical_data_integration.R
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Do not run !!!
# To run, go to each script and run code chunk by chunk !!!
# source("R/s5_manual_clinical_data_integration.R")

# Robject updated: geo_tidy (Updated with cleaned clinical data.)
# Robject created: geo_series_matrix_summary (Integrated series and series_matrix summary.)
#                  geo_clin (Integrated clinical data with regimen classification.
#                            Note that geo_clin is generated from geo.RData and contains
#                            38 series matrices, 35 datasets.)
# Relevant Robjects:  1) geo_series_summary, 2) geo, 3) gpl, 4) geo_maxvar_annot
#                     5) geo_tidy, 6) geo_expr, 7) geo_expr_meta,
#                     8) geo_series_matrix_summary, 9) geo_clin




# s6_process_and_validate_integrated_expression_matrix.R
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Do not run !!!
# To run, go to each script and run code chunk by chunk !!!
# source("R/s6_process_and_validate_integrated_expression_matrix.R")

# Robject created: geo_expr_list (List of dataset effect correted pooled data.
#                                 Contains 1) Non-corrected, 2) Combat corrected,
#                                 3) Quantile corrected, and 4) Madquantile corrected
#                                 pooled datasets. All veriosn of dataets
#                                 contains 3736 samples, except Madquntile corrected
#                                 dataset which contins 3734 samples. The two samples
#                                 were filtered out due to zero MAD.
#                  geo_clin_list (Simialr strucutre to geo_expr_list.
#                                 But contains clinical data with respective
#                                 method's qc metrics (RLE))
# Relevant Robjects:  1) geo_series_summary, 2) geo, 3) gpl, 4) geo_maxvar_annot
#                     5) geo_tidy, 6) geo_expr, 7) geo_expr_meta,
#                     8) geo_series_matrix_summary, 9) geo_clin,
#                     10) geo_expr_list, 11) geo_clin_list




# s7_figures_tables_data.R
# >>>>>>>>>>>>>>>>>>>>>>>>

# Do not run !!!
# To run, go to each script and run code chunk by chunk !!!
# source("R/s7_figures_tables_data.R")

# Output: Figures, tables, text files.
# Relevant Robjects:  1) geo_series_summary, 2) geo, 3) gpl, 4) geo_maxvar_annot
#                     5) geo_tidy, 6) geo_expr, 7) geo_expr_meta,
#                     8) geo_series_matrix_summary, 9) geo_clin,
#                     10) geo_expr_list, 11) geo_clin_list


#
# ==============================================================================


