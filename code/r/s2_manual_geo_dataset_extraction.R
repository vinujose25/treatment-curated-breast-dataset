# s2_manual_geo_dataset_extraction.R

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# 1. Compuational extraction of series, sample and expression data from each
# geo series matrices.
# 2. For sample characterisitcs which are failed to format(curate) compuateionally,
# do a manual curation.
# 3. Manual formatting performed:
# 3.1. Changing charcteristic names (column names)
# 3.2. Correcting mis-alinged samples_characterisitcs spreadsheet-cells due to
#     missingness.
# 3.3. Removing the concatenated characterisitcs name to each value.


# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Compuational extraction of series, sample and expression data from each
# geo series matrices.
# 2. Manual formatting of var_qc_failed_sample_characteristics
# 3. Append the formatted(curated) sample characterisitics to geo object
# 4. Backup: Unexpected behaviours or  bugs in format_geo() which are resolved.




# 1. Compuational extraction of series, sample and expression data from each
# geo series matrices.
# =============================================================================


files = c(
  list.files("data/gse", full.names = TRUE)
)
geo <- format_geo(files = files)

names(geo) <-  str_split_fixed(names(geo),"/", 3)[, 3] %>%
  str_replace("_series_matrix.txt.gz", "") %>%
  str_replace("-", "_")

#
# ==============================================================================



# 2. Manual formatting of var_qc_failed_sample_characteristics
# ==============================================================================


# Identify series matrices with issues with sample characterisitics
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

table(purrr::map_chr(geo, ~(.x$dataset_issues)))
# No Yes
# 24  20
nme <- names(geo)[purrr::map_chr(geo, ~(.x$dataset_issues)) == "Yes"]
# [1] "GSE114403"        "GSE16446"         "GSE20194"
# [4] "GSE20271"         "GSE20685"         "GSE21997_GPL5325"
# [7] "GSE22093"         "GSE22226_GPL1708" "GSE22226_GPL4133"
# [10] "GSE22358"         "GSE25066"         "GSE28844"
# [13] "GSE31863"         "GSE41998"         "GSE50948"
# [16] "GSE66999"         "GSE6861"          "GSE75678"
# [19] "GSE8465_GPL1390"  "GSE8465_GPL887"



# Delete comment !!!!!!!!!!!!!!!!!!!!!!!!!!
# # For the below series, already processed data present in
# # "results/geo/var_qc_failed_sample_characteristics" is used
#
# # "GSE20194" "GSE20271"           "GSE41998"           "GSE50948"
# # "GSE6861"  "GSE25066"
# End delete comment !!!!!!!!!!!!!!!!!!!!!!



# Write out sample characterisitcs along with complete sample data
# for manuall curation.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


for(inme in nme){
  # Full sample data is written out for reference to identify problem in
  # sample characterisitcs if necessary.
  write_tsv(
    geo[[inme]]$sample,
    path = str_c(out_data, inme, "_sample.tsv")
  )

  # Sample characteristics data of series matrices for which atleast one
  # sample characterisitc has doubious value.
  write_tsv(
    geo[[inme]]$sample_characteristics,
    path = str_c(out_data, inme, "_sample_characteristics.tsv")
  )
}
# The above files are kept in results/data/var_qc_failed_sample_characteristics/
# The corrected files have the suffix "_cleaned"

# Manual formatting performed:
# 1. Changing charcteristic names (column names)
# 2. Correcting mis-alinged samples_characterisitcs spreadsheet-cells due to
#     missingness.
# 3. Removing the concatenated characterisitcs name to each value.

#
# ==============================================================================



# 3.Append the formatted(curated) sample characterisitics to geo object
# ==============================================================================

files <- list.files("results/data/var_qc_failed_sample_characteristics", full.names = TRUE)
files <- files[str_detect(files, "_cleaned")]
id <- str_split_fixed(files,"/",4)[,4]
id <- str_split_fixed(id, "_sample_characteristics_cleaned.tsv", 2)[, 1]
names(files) = id

# files = files[c("GSE8465_GPL1390","GSE8465_GPL887")]


for(nme in names(files)){

  print(nme)
  ifile = files[nme]

  geo[[nme]]$sample_characteristics <-
    read_tsv(file = ifile, col_types = cols())
}
length(geo) # 44

# Checking appended data
# geo$GSE8465_GPL1390$sample_characteristics
# geo$GSE8465_GPL887$sample_characteristics
# geo$GSE143846$sample_characteristics


# Saving
# save(geo, file = str_c(out_data, "geo.RData")) # updated geo

#
# ==============================================================================



# 4. Backup: Unexpected behaviours or  bugs in format_geo() which are resolved.
# Kept it as a back up
# ==============================================================================

# GSE6861
# >>>>>>>>>>
# Issue: Sample charactristic column doesn't follow the usual format of
# "characteristic name : characteristic value". Hence the column name in
# GSE6861$sample_characteristics is the first value of the
# characterisitic value vector".
# !!!!! This scenario is a bug, fix it, if no ":", generate custom column names, v1,v2 etc

# GSE20271
# >>>>>>>>>>>
# Issue1: The sample characteristic colum with values "dlda30 pred (1=pcr, 0=rd): 1"
# did not parse correctly. The values are empty in the parsed column.
# This is due to the fact that the parenthesis in the column name are special characters
# The split pattern uses column name to split. This is clearly a bug, and can be fixed by
# use only the ":" as pattern for split.

# Issue2: The sample characteristic colum with values "post chemo +ln/total: 0/12"
# although would parse correctly, the forward slash (or a  dash) would be interpreted as
#  date by any spreadsheet like sotware.
#  Its better to not parse it and genrate a var_qc varning instead.

# Note !!!!!!!!!
# The column "post chemo +ln/total" requires computational cleaning
# !!!!!!!!!!!!!!


# GSE25066
# >>>>>>>>>
# Comment !!!
# Introduced a new column "type_taxane". Which is hide under "drfs_even_time_years"
# column

# GSE41998
# >>>>>>>>
# Comment !!!!
# If NAs present in sample characteristics column, var_qc will fail.
# NA values should present as "charactristic name : NA", to pass var_qc check.

# GSE75685
# >>>>>>>>>>>
# The coulmns "rna 260/280 (nanodrop)"	and "rna 260/230 (nanodrop)" hase the issue
# same as issue1 of GSE20271.
#
# End of backup notes

#
# ==============================================================================















