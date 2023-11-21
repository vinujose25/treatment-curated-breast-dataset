 # s4_manual_expression_data_integration.R

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# 1. Cleaning: removing genes with at leans one NA, log space converion
# 2. Max-var collapsing: Multiple probes representing same gene is collpased
# by using the maximum varient probe#
# 3. Integrate multiple dataset using common genes


# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Prepare data

# !!! Part 1 (cleaning) !!!!!!!!!!!!(Part 1 filtered one dataset out)!!!!!!!!!!!
# 2. Remove NA expression values
# 3. Convert expression data to log space

# !!! Part 2 (maxvar collapsing + integration) !!!(Part 2 filtered 7 datsets out)
# 4. Do max-var collapsing
# 5. Genrate geo_tidy, a clean version of geo with only expression and clinical
# slots to hold max-var collaped expression data and cleaned clinical data.
# 6. Expression data integration



# 1. Prepare data
# ==============================================================================

# load(str_c(out_data, "geo.RData"))
length(geo) # 39 series matrices

#
# ==============================================================================




# Part 1 (cleaning) !!!!!!!!!!!!(Part 1 filtered one dataset out)!!!!!!!!!!!!!!!
# 2) Remove NA expression values
# 3) Convert expression data to log space
# ==============================================================================


# Identify class of all expression profiles
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

for(inme in names(geo)){
  xx <- geo[[inme]]$expression[, -1]
  cls <- purrr::map_chr(xx,~(class(.x))) %>% unique()
  print(str_c(inme, cls, sep = " "))
}

# All numeric except GSE31863 and GSE69031
# !!!!!!!!!!!!!!!!!!!!
# "GSE31863 character" !!!! contains "null" for missing values
# "GSE69031 numeric" "GSE69031 logical" !!!! Few samples got missing expression for all probesets


# GSE31863: Char to Numeric
geo$GSE31863$expression[,-1] <- purrr::map(
  geo$GSE31863$expression[,-1],
  ~(as.numeric(.x))) %>%
  bind_cols()

# GSE69031: identify probsets with NAs for all samples
idx <- purrr::map_lgl(
  geo$GSE69031$expression,
  ~(all(is.na(.x))))
table(idx)
# FALSE  TRUE
# 119    12
# This will be removed later in the script,
#  while exploring NAs in other datasets !!!!!!!!!!!!!!!!!!!!!!!!!!




# Explore scale of expression data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

for(inme in names(geo)){

  xx <- geo[[inme]]$expression[, -1]
  if (nrow(xx) > 0) {
    print(paste(inme, ", ", str_c(dim(xx), collapse = " x "), ", ",
                str_c(range(xx), collapse = " - ")))
  } else {
    print(paste(inme, ", ", str_c(dim(xx), collapse = " x ")))
  }
}

# [1] "GSE109710 ,  815 x 173 ,  -13.2930747796743 - 10.1827742985521"
# [1] "GSE114403 ,  784 x 100 ,  -8.549153547 - 47081"
# [1] "GSE130786 ,  41000 x 110 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE143222 ,  579 x 55 ,  0 - 78830.53087"
# [1] "GSE143846 ,  42545 x 44 ,  -8.78502 - 12.247937"
# [1] "GSE16391 ,  54675 x 55 ,  2.224229557 - 14.93376538"
# [1] "GSE16446 ,  54675 x 120 ,  1.979789284 - 14.59179297"
# [1] "GSE18728 ,  54675 x 61 ,  5.643856049 - 15.07926178"
# [1] "GSE19615 ,  54675 x 115 ,  0.402 - 14222.174"
# [1] "GSE20194 ,  22283 x 278 ,  -3.1553 - 20.0716"
# [1] "GSE20271 ,  22283 x 178 ,  0 - 31168.7"
# [1] "GSE20685 ,  54627 x 327 ,  0 - 15.881485921769"
# [1] "GSE21974 ,  41000 x 57 ,  0.6029321 - 18.005314"
# [1] "GSE21997_GPL1390 ,  22575 x 35 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE21997_GPL5325 ,  44290 x 32 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE21997_GPL7504 ,  45220 x 31 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE22093 ,  22283 x 103 ,  -2.994992002 - 18.61908686"
# [1] "GSE22219 ,  24332 x 216 ,  5.970680292 - 15.39221539"
# [1] "GSE22226_GPL1708 ,  44290 x 130 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE22226_GPL4133 ,  45220 x 20 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE22358 ,  44290 x 158 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE23988 ,  22283 x 61 ,  -2.434076364 - 18.70437604"
# [1] "GSE25066 ,  22283 x 508 ,  -5.026351893 - 21.0321122"
# [1] "GSE28844 ,  54675 x 61 ,  2.1109 - 14.6917"
# [1] "GSE31863 ,  26766 x 143 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE32603 ,  35069 x 248 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE32646 ,  54613 x 115 ,  -5.254537741 - 16.53461248"
# [1] "GSE34138 ,  27506 x 178 ,  -1.887129192 - 16.25492742"
# [1] "GSE4056 ,  21888 x 190 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE41998 ,  22277 x 279 ,  2.498 - 14.773"
# [1] "GSE42822 ,  22283 x 91 ,  2.807297 - 14.3512"
# [1] "GSE45255 ,  22268 x 139 ,  -2.7055 - 17.3699"
# [1] "GSE4779 ,  61297 x 102 ,  2.6593512253927 - 14.8591598999343"
# [1] "GSE50948 ,  54675 x 156 ,  0.514264471 - 13.80251035"
# [1] "GSE66305 ,  54675 x 88 ,  1.73888 - 14.8914"
# [1] "GSE66999 ,  41028 x 76 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE6861 ,  61359 x 161 ,  3.26114481513878 - 14.8676734288273"
# [1] "GSE69031 ,  22215 x 130 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE76360 ,  44810 x 100 ,  5.830292236 - 14.15377231"





# Explore expression datasets with NA expression values and filter NAs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# [1] "GSE130786 ,  41000 x 110 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE21997_GPL1390 ,  22575 x 35 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE21997_GPL5325 ,  44290 x 32 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE21997_GPL7504 ,  45220 x 31 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE22226_GPL1708 ,  44290 x 130 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE22226_GPL4133 ,  45220 x 20 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE22358 ,  44290 x 158 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE31863 ,  26766 x 143 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE32603 ,  35069 x 248 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE4056 ,  21888 x 190 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE66999 ,  41028 x 76 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE69031 ,  22215 x 130 ,  NA" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



explore_na <- function(xx){
  #xx expresison matrix tibble; first column ID_REF contains probeset ids
  print("sample")
  idx_sample <- purrr::map_lgl(xx, ~(any(is.na(.x))))
  print("gene")
  idx_gene <- purrr::map_lgl(t(xx) %>% as_tibble(), ~(any(is.na(.x))))

  return(list(na_sample = idx_sample, na_gene = idx_gene))
}


# Names of expression data with NAs
nme <- c("GSE130786", "GSE21997_GPL1390", "GSE21997_GPL5325", "GSE21997_GPL7504",
         "GSE22226_GPL1708", "GSE22226_GPL4133", "GSE22358", "GSE31863",
         "GSE32603", "GSE4056", "GSE66999", "GSE69031")
na_index <- purrr::map(geo[nme],~(explore_na(.x$expression)))

purrr::map_lgl(na_index,~(any(.x$na_sample))) # all dataset has NAs


for(inme in nme){
  ixx = na_index[[inme]]
  print(inme)
  print("sample")
  print(table(ixx$na_sample))
  print("gene")
  print(table(ixx$na_gene))
  print("--------------------")
}

# [1] "GSE130786" #### remove genes
# [1] "sample"
# FALSE  TRUE
# 1   110
# [1] "gene"
# FALSE  TRUE
# 36286  4714

# [1] "GSE21997_GPL1390" ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 1    35
# [1] "gene"
# FALSE  TRUE
# 21631   944

# [1] "GSE21997_GPL5325" ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 2    31
# [1] "gene"
# FALSE  TRUE
# 43898   392


# [1] "GSE21997_GPL7504"  ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 5    27
# [1] "gene"
# FALSE  TRUE
# 44940   280

# [1] "GSE22226_GPL1708"  ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 7   124
# [1] "gene"
# FALSE  TRUE
# 28349 15941

# [1] "GSE22226_GPL4133"  ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 2    19
# [1] "gene"
# FALSE  TRUE
# 44795   425

# [1] "GSE22358"  ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 3   156
# [1] "gene"
# FALSE  TRUE
# 41810  2480

# [1] "GSE31863"  ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 1   143
# [1] "gene"
# FALSE  TRUE
# 11851 14915

# [1] "GSE32603"  ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 1   248
# [1] "gene"
# FALSE  TRUE
# 12901 22168

# [1] "GSE4056"  ##### remove dataset !!!!! (more than 90% genes have atleast one NA value)
# [1] "sample"
# FALSE  TRUE
# 1   190
# [1] "gene"
# FALSE  TRUE
# 149 21739

# [1] "GSE66999" ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 2    75
# [1] "gene"
# FALSE  TRUE
# 38706  2322

# [1] "GSE69031" ##### remove samples (12 samples has no expression)
# [1] "sample"
# FALSE  TRUE
# 119    12
# [1] "gene"
# TRUE
# 22215 n



# !!!!!!!!! Removing GSE4056 !!!!!!!!!!!!!!!!!!!!
#
geo <- geo[-which(names(geo) == "GSE4056")]




# !!!!!!!!! Removing 12 samples of GSE69031 with no expression !!!!!!!!!!!!!!!!
#
names(geo$GSE69031)
# [1] "series"                 "sample"                 "expression"
# [4] "sample_characteristics" "dataset_issues"         "clinical"

# Remove 12 samples from  sample, expression, sample_characterisitcs and clinical
nme <- names(na_index$GSE69031$na_sample)[na_index$GSE69031$na_sample]

# clinical, sample,sample_characteristics are of identical row order
identical(geo$GSE69031$sample$Sample_geo_accession, geo$GSE69031$clinical$Sample_geo_accession)
# TRUE
idx = which(geo$GSE69031$sample$Sample_geo_accession %in% nme)
identical(nme, geo$GSE69031$sample$Sample_geo_accession[idx]) # TRUE

geo$GSE69031$sample <- geo$GSE69031$sample[-idx, ]
geo$GSE69031$sample_characteristics <- geo$GSE69031$sample_characteristics[-idx, ]
geo$GSE69031$clinical <- geo$GSE69031$clinical[-idx, ]

identical(names(geo$GSE69031$expression), names(na_index$GSE69031$na_sample)) # TRUE
geo$GSE69031$expression <- geo$GSE69031$expression[, -(which(na_index$GSE69031$na_sample))]

head2(geo$GSE69031$sample) # 118  56
head2(geo$GSE69031$sample_characteristics) # 118  25
head2(geo$GSE69031$clinical) # 118  50
head2(geo$GSE69031$expression) # 22215   119

identical(geo$GSE69031$sample$Sample_geo_accession, geo$GSE69031$clinical$Sample_geo_accession)
identical(geo$GSE69031$sample$Sample_geo_accession, names(geo$GSE69031$expression)[-1])
# Both TRUE



# !!!!!!!!! Removing NA genes from the rest !!!!!!!!!!!!!!!!
na_index <- na_index[!(names(na_index) %in% c("GSE4056", "GSE69031"))]

for(inme in names(na_index)){
  ixx = na_index[[inme]]
  geo[[inme]]$expression <- geo[[inme]]$expression[!ixx$na_gene, ]
}



# Checking if NAs still exist in geo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Names of expression data with NAs
nme <- c("GSE130786", "GSE21997_GPL1390", "GSE21997_GPL5325", "GSE21997_GPL7504",
         "GSE22226_GPL1708", "GSE22226_GPL4133", "GSE22358", "GSE31863",
         "GSE32603", "GSE66999", "GSE69031") # "GSE4056" is removed from geo
na_index <- purrr::map(geo[nme],~(explore_na(.x$expression)))

purrr::map_lgl(na_index,~(any(.x$na_sample))) # all dataset has no more NAs


identical(geo$GSE69031$sample$Sample_geo_accession,
          geo$GSE69031$clinical$Sample_geo_accession) # TRUE





# Explore scale of expression data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

for(inme in names(geo)){

  xx <- geo[[inme]]$expression[, -1]
  if (nrow(xx) > 0) {
    print(paste(inme, ", ", str_c(dim(xx), collapse = " x "), ", ",
                str_c(range(xx), collapse = " - ")))
  } else {
    print(paste(inme, ", ", str_c(dim(xx), collapse = " x ")))
  }
}

# [1] "GSE109710 ,  815 x 173 ,  -13.2930747796743 - 10.1827742985521"
# [1] "GSE114403 ,  784 x 100 ,  -8.549153547 - 47081" !!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE130786 ,  36286 x 110 ,  -2.85191 - 2.40524"
# [1] "GSE143222 ,  579 x 55 ,  0 - 78830.53087" !!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE143846 ,  42545 x 44 ,  -8.78502 - 12.247937"
# [1] "GSE16391 ,  54675 x 55 ,  2.224229557 - 14.93376538"
# [1] "GSE16446 ,  54675 x 120 ,  1.979789284 - 14.59179297"
# [1] "GSE18728 ,  54675 x 61 ,  5.643856049 - 15.07926178"
# [1] "GSE19615 ,  54675 x 115 ,  0.402 - 14222.174" !!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE20194 ,  22283 x 278 ,  -3.1553 - 20.0716"
# [1] "GSE20271 ,  22283 x 178 ,  0 - 31168.7" !!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE20685 ,  54627 x 327 ,  0 - 15.881485921769"
# [1] "GSE21974 ,  41000 x 57 ,  0.6029321 - 18.005314"
# [1] "GSE21997_GPL1390 ,  21631 x 35 ,  -10.643 - 10.834"
# [1] "GSE21997_GPL5325 ,  43898 x 32 ,  -8.33 - 10.436"
# [1] "GSE21997_GPL7504 ,  44940 x 31 ,  -11.729 - 9.953"
# [1] "GSE22093 ,  22283 x 103 ,  -2.994992002 - 18.61908686"
# [1] "GSE22219 ,  24332 x 216 ,  5.970680292 - 15.39221539"
# [1] "GSE22226_GPL1708 ,  28349 x 130 ,  -10.505 - 11.122"
# [1] "GSE22226_GPL4133 ,  44795 x 20 ,  -10.104 - 10.005"
# [1] "GSE22358 ,  41810 x 158 ,  -9.931 - 10.699"
# [1] "GSE23988 ,  22283 x 61 ,  -2.434076364 - 18.70437604"
# [1] "GSE25066 ,  22283 x 508 ,  -5.026351893 - 21.0321122"
# [1] "GSE28844 ,  54675 x 61 ,  2.1109 - 14.6917"
# [1] "GSE31863 ,  11851 x 143 ,  -12.07781529 - 10.14477606"
# [1] "GSE32603 ,  12901 x 248 ,  -10.20233089 - 10.34187329"
# [1] "GSE32646 ,  54613 x 115 ,  -5.254537741 - 16.53461248"
# [1] "GSE34138 ,  27506 x 178 ,  -1.887129192 - 16.25492742"
# [1] "GSE41998 ,  22277 x 279 ,  2.498 - 14.773"
# [1] "GSE42822 ,  22283 x 91 ,  2.807297 - 14.3512"
# [1] "GSE45255 ,  22268 x 139 ,  -2.7055 - 17.3699"
# [1] "GSE4779 ,  61297 x 102 ,  2.6593512253927 - 14.8591598999343"
# [1] "GSE50948 ,  54675 x 156 ,  0.514264471 - 13.80251035"
# [1] "GSE66305 ,  54675 x 88 ,  1.73888 - 14.8914"
# [1] "GSE66999 ,  38706 x 76 ,  3.321928024 - 17.30794716"
# [1] "GSE6861 ,  61359 x 161 ,  3.26114481513878 - 14.8676734288273"
# [1] "GSE69031 ,  22215 x 118 ,  3.758497 - 15.08126"
# [1] "GSE76360 ,  44810 x 100 ,  5.830292236 - 14.15377231"


# Non-log2 datasets
# [1] "GSE114403 ,  784 x 100 ,  -8.549153547 - 47081" !!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE143222 ,  579 x 55 ,  0 - 78830.53087" !!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE19615 ,  54675 x 115 ,  0.402 - 14222.174" !!!!!!!!!!!!!!!!!!!!!!!!!
# [1] "GSE20271 ,  22283 x 178 ,  0 - 31168.7" !!!!!!!!!!!!!!!!!!!!!!!!!


summary(geo$GSE114403$expression[, -1] %>% unlist())
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -8.55    -2.49    -0.84    56.18     0.95 47081.00
# Nanostring data; filter NEG and POS controls, everythig else is in log2 scale

summary(geo$GSE143222$expression[, -1] %>% unlist())
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00     1.63     5.60   110.07    23.04 78830.53
# Nanostring data; No control present, convert to log2 scale

summary(geo$GSE19615$expression[, -1] %>% unlist())
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.402    24.704    53.993   163.672   135.482 14222.174
# u133plus2 convert to log2 scale

summary(geo$GSE20271$expression[, -1] %>% unlist())
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0   112.1   194.6   399.9   342.5 31168.7
# u133a convert to log2 scale



# !!!!!!!!!!!! GSE114403 filter neg and pos control probes !!!!!!!!!!!!!!!!!!!!!
xx <- geo$GSE114403$expression %>%
  dplyr::filter(!str_detect(ID_REF,"POS_")) %>%  # 6
  dplyr::filter(!str_detect(ID_REF,"NEG_")) # 8
range(xx[,-1])
# [1] -8.549154 13.329376 # log2 scale
dim(xx) # 770 101
dim(geo$GSE114403$expression) # 784 101

geo$GSE114403$expression <- xx


# !!!!!!!!!!!! GSE143222 convert to log2 scale !!!!!!!!!!!!!!!!!!!!!
xx <- geo$GSE143222$expression
range(xx[,-1])
# 0.00 78830.53
xx[,-1] <- log2(xx[,-1] + (min(xx[,-1])+0.1))
range(xx[,-1])
# -3.321928 16.266469

geo$GSE143222$expression <- xx


# !!!!!!!!!!!! GSE19615 convert to log2 scale !!!!!!!!!!!!!!!!!!!!!
xx <- geo$GSE19615$expression
range(xx[,-1])
# 0.402 14222.174
xx[,-1] <- log2(xx[,-1] + (min(xx[,-1])+0.1))
range(xx[,-1])
# -0.1456053 13.7959053

geo$GSE19615$expression <- xx



# !!!!!!!!!!!! GSE20271 convert to log2 scale !!!!!!!!!!!!!!!!!!!!!
xx <- geo$GSE20271$expression
range(xx[,-1])
# 0.0 31168.7
xx[,-1] <- log2(xx[,-1] + (min(xx[,-1])+0.1))
range(xx[,-1])
# -3.321928 14.927815

geo$GSE20271$expression <- xx



# Double checking scale of expression matrices

# Scale
for(inme in names(geo)){

  xx <- geo[[inme]]$expression[, -1]
  if (nrow(xx) > 0) {
    print(paste(inme, ", ", str_c(dim(xx), collapse = " x "), ", ",
                str_c(range(xx), collapse = " - ")))
  } else {
    print(paste(inme, ", ", str_c(dim(xx), collapse = " x ")))
  }
}
# [1] "GSE109710 ,  815 x 173 ,  -13.2930747796743 - 10.1827742985521"
# [1] "GSE114403 ,  770 x 100 ,  -8.549153547 - 13.32937591"
# [1] "GSE130786 ,  36286 x 110 ,  -2.85191 - 2.40524"
# [1] "GSE143222 ,  579 x 55 ,  -3.32192809488736 - 16.2664686998096"
# [1] "GSE143846 ,  42545 x 44 ,  -8.78502 - 12.247937"
# [1] "GSE16391 ,  54675 x 55 ,  2.224229557 - 14.93376538"
# [1] "GSE16446 ,  54675 x 120 ,  1.979789284 - 14.59179297"
# [1] "GSE18728 ,  54675 x 61 ,  5.643856049 - 15.07926178"
# [1] "GSE19615 ,  54675 x 115 ,  -0.145605322246899 - 13.7959053134741"
# [1] "GSE20194 ,  22283 x 278 ,  -3.1553 - 20.0716"
# [1] "GSE20271 ,  22283 x 178 ,  -3.32192809488736 - 14.9278149917673"
# [1] "GSE20685 ,  54627 x 327 ,  0 - 15.881485921769"
# [1] "GSE21974 ,  41000 x 57 ,  0.6029321 - 18.005314"
# [1] "GSE21997_GPL1390 ,  21631 x 35 ,  -10.643 - 10.834"
# [1] "GSE21997_GPL5325 ,  43898 x 32 ,  -8.33 - 10.436"
# [1] "GSE21997_GPL7504 ,  44940 x 31 ,  -11.729 - 9.953"
# [1] "GSE22093 ,  22283 x 103 ,  -2.994992002 - 18.61908686"
# [1] "GSE22219 ,  24332 x 216 ,  5.970680292 - 15.39221539"
# [1] "GSE22226_GPL1708 ,  28349 x 130 ,  -10.505 - 11.122"
# [1] "GSE22226_GPL4133 ,  44795 x 20 ,  -10.104 - 10.005"
# [1] "GSE22358 ,  41810 x 158 ,  -9.931 - 10.699"
# [1] "GSE23988 ,  22283 x 61 ,  -2.434076364 - 18.70437604"
# [1] "GSE25066 ,  22283 x 508 ,  -5.026351893 - 21.0321122"
# [1] "GSE28844 ,  54675 x 61 ,  2.1109 - 14.6917"
# [1] "GSE31863 ,  11851 x 143 ,  -12.07781529 - 10.14477606"
# [1] "GSE32603 ,  12901 x 248 ,  -10.20233089 - 10.34187329"
# [1] "GSE32646 ,  54613 x 115 ,  -5.254537741 - 16.53461248"
# [1] "GSE34138 ,  27506 x 178 ,  -1.887129192 - 16.25492742"
# [1] "GSE41998 ,  22277 x 279 ,  2.498 - 14.773"
# [1] "GSE42822 ,  22283 x 91 ,  2.807297 - 14.3512"
# [1] "GSE45255 ,  22268 x 139 ,  -2.7055 - 17.3699"
# [1] "GSE4779 ,  61297 x 102 ,  2.6593512253927 - 14.8591598999343"
# [1] "GSE50948 ,  54675 x 156 ,  0.514264471 - 13.80251035"
# [1] "GSE66305 ,  54675 x 88 ,  1.73888 - 14.8914"
# [1] "GSE66999 ,  38706 x 76 ,  3.321928024 - 17.30794716"
# [1] "GSE6861 ,  61359 x 161 ,  3.26114481513878 - 14.8676734288273"
# [1] "GSE69031 ,  22215 x 118 ,  3.758497 - 15.08126"
# [1] "GSE76360 ,  44810 x 100 ,  5.830292236 - 14.15377231"

# All in log scale !!!!!!!!!!!!!!!!!!!!!


# Convert ID_REF to character (Some expression$ID_REF are numerics)
geo <- purrr::map(
  geo,
  function(x){
    print(class(x$expression$ID_REF))
    x$expression$ID_REF <- as.character(x$expression$ID_REF)
    x
  }
)

# Saving
# save(geo, file = str_c(out_data,"geo.RData"))

#
# ==============================================================================




# Part 2 (maxvar collapsing + integration) !!!(Part 2 filtered 7 datsets out)!!!
# 4. Do max-var collapsing
# ==============================================================================


# load(str_c(out_data,"geo.RData"))


table(purrr::map_chr(geo,~class(.x$expression$ID_REF)))
# All character
purrr::map(geo,~(.x$expression[1:3,1:5]))
# ID_REF char, all others double


# Preparing gse - gpl(pltform) mapping
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

gse_gpl_map <- tibble(
  Series_accession = str_split_fixed(names(geo), "_", 2)[, 1],
  Series_matrix_accession = names(geo),
  Series_matrix_platform_id = purrr::map_chr(geo, ~(unique(.x$clinical$Sample_platform_id)))
)

table(gse_gpl_map$Series_matrix_platform_id)
# GPL1352  GPL1390 GPL14374 GPL14550 GPL14668  GPL1708  GPL18649 GPL24546 GPL24993
# 2        1        1        1        1        1        1        1        1
# GPL4133  GPL5325   GPL570   GPL571  GPL6098  GPL6480  GPL6884  GPL6947  GPL7504
# 1        2        9        2        1        3        1        1        1
# GPL96
# 7


# Manually download and clean annotations
# Must present column ID_REF, Ncbi_gene_id, Hugo_gene_symbol
# Downloaded gpl files are in the folder data/gpl/


# Accession number formats
# https://www.ncbi.nlm.nih.gov/Sequin/acc.html
# http://www.nslc.wustl.edu/elgin/genomics/bio4342/1archives/2006/AccReference.pdf


# Load annotations
# >>>>>>>>>>>>>>>>

files = list.files("data/gpl", full.names = TRUE)
files = files[str_detect(files, "_cleaned")]
gpl <- purrr::map(
  files,
  function(x){
    xx <- read_tsv(file = x,
                   col_types = cols(.default = col_character()),
                   comment = "#",
                   guess_max = 5000)
    names(xx) <- names(xx) %>%
      str_replace_all(" ","_") %>%
      str_to_sentence()
    xx
  }
)
names(gpl) <- files

names(gpl) <-  str_split_fixed(names(gpl),"/", 3)[, 3]
names(gpl) <-  str_split_fixed(names(gpl),"-", 2)[, 1]


# columns in gpl
purrr::map(gpl, names) %>% unlist() %>% table()
# Ensembl_gene_id Ensembl_transcript_id        Entrez_gene_id
# 2                     5                    13
# Genbank_accession      Hugo_gene_symbol             Id_ref
# 18                    17                    19
# Illumina_gene   Illumina_transcript            Probe_name
# 2                     2                     2
# Refseq_accession
# 12


# Checking congruence between geo expression matrix and gpl platform
gse_gpl_cong <- purrr::map2_dfr(
  gse_gpl_map$Series_matrix_accession,
  gse_gpl_map$Series_matrix_platform_id,
  function(nme_gse, nme_gpl, geo, gpl){

    xgse <- geo[[nme_gse]]$expression
    xgpl <- gpl[[nme_gpl]] %>% dplyr::rename(ID_REF = "Id_ref")

    list(
      # eq:equal, ex:expression, an:annotation
      n_eq = nrow(xgse) == nrow(xgpl),
      n_ex =nrow(xgse),
      n_ex2an = sum(xgse$ID_REF %in% xgpl$ID_REF),
      # proportion of expressed probes(non-NA) with annotations
      prop_ex2n = sum(xgse$ID_REF %in% xgpl$ID_REF)/nrow(xgse) %>% round(digits = 2),
      n_an = nrow(xgpl),
      n_an2ex = sum(xgpl$ID_REF %in% xgse$ID_REF),
      # proportion of annotated genes with (non-NA) expression value
      prop_an2ex = sum(xgpl$ID_REF %in% xgse$ID_REF)/nrow(xgpl) %>% round(digits = 2),
      gse = nme_gse,
      gpl = nme_gpl
    )
  },
  geo,
  gpl
)
# Due to removal of genes which has at least one NA in the dataset causes
# non-congruense between expression and annotation.
# 1) prop_ex2n < 1 means annotation issue, which is hard to correct.
# 2) prop_an2ex < 1 means missignness in expression values for annotated genes, mainly due to
# removal of probes with atleast one NA expression values per dataset.
# Large missingness in expression values is problematic.

gse_gpl_cong %>% as.data.frame()
#     n_eq  n_ex n_ex2an prop_ex2n  n_an n_an2ex prop_an2ex              gse      gpl
# 1  FALSE   815     815         1   819     815  0.9951160        GSE109710 GPL24546
# 2  FALSE   770     770         1   788     770  0.9771574        GSE114403 GPL24993
# 3  FALSE 36286   36286         1 41108   36286  0.8826992        GSE130786  GPL6480
# 4  FALSE   579     579         1   598     579  0.9682274        GSE143222 GPL18649
# 5   TRUE 42545   42545         1 42545   42545  1.0000000        GSE143846 GPL14550
# 6   TRUE 54675   54675         1 54675   54675  1.0000000         GSE16391   GPL570
# 7   TRUE 54675   54675         1 54675   54675  1.0000000         GSE16446   GPL570
# 8   TRUE 54675   54675         1 54675   54675  1.0000000         GSE18728   GPL570
# 9   TRUE 54675   54675         1 54675   54675  1.0000000         GSE19615   GPL570
# 10  TRUE 22283   22283         1 22283   22283  1.0000000         GSE20194    GPL96
# 11  TRUE 22283   22283         1 22283   22283  1.0000000         GSE20271    GPL96
# 12 FALSE 54627   54627         1 54675   54627  0.9991221         GSE20685   GPL570
# 13 FALSE 41000   41000         1 41108   41000  0.9973728         GSE21974  GPL6480
# 14 FALSE 21631   21631         1 22579   21631  0.9580141 GSE21997_GPL1390  GPL1390
# 15 FALSE 43898   43898         1 44294   43898  0.9910597 GSE21997_GPL5325  GPL5325
# 16 FALSE 44940   44940         1 45220   44940  0.9938080 GSE21997_GPL7504  GPL7504
# 17  TRUE 22283   22283         1 22283   22283  1.0000000         GSE22093    GPL96
# 18 FALSE 24332   24332         1 24385   24332  0.9978265         GSE22219  GPL6098
# 19 FALSE 28349   28349         1 44290   28349  0.6400768 GSE22226_GPL1708  GPL1708 ! probe missing
# 20 FALSE 44795   44795         1 45220   44795  0.9906015 GSE22226_GPL4133  GPL4133
# 21 FALSE 41810   41810         1 44294   41810  0.9439202         GSE22358  GPL5325
# 22  TRUE 22283   22283         1 22283   22283  1.0000000         GSE23988    GPL96
# 23  TRUE 22283   22283         1 22283   22283  1.0000000         GSE25066    GPL96
# 24  TRUE 54675   54675         1 54675   54675  1.0000000         GSE28844   GPL570
# 25 FALSE 11851   11851         1 26823   11851  0.4418223         GSE31863 GPL14374 ! probe missing
# 26 FALSE 12901   12901         1 35073   12901  0.3678328         GSE32603 GPL14668 ! probe missing
# 27 FALSE 54613   54613         1 54675   54613  0.9988660         GSE32646   GPL570
# ?????????????????????
# 28 FALSE 27506   27506         1 48803   27506  0.5636129         GSE34138  GPL6884 ! probe missing
# ?????????????? Probes are missing in the original series matrix (not due to NAs)
# 29  TRUE 22277   22277         1 22277   22277  1.0000000         GSE41998   GPL571
# 30  TRUE 22283   22283         1 22283   22283  1.0000000         GSE42822    GPL96
# 31 FALSE 22268   22268         1 22283   22268  0.9993268         GSE45255    GPL96
# 32 FALSE 61297   61297         1 61359   61297  0.9989896          GSE4779  GPL1352
# 33  TRUE 54675   54675         1 54675   54675  1.0000000         GSE50948   GPL570
# 34  TRUE 54675   54675         1 54675   54675  1.0000000         GSE66305   GPL570
# 35 FALSE 38706   38706         1 41108   38706  0.9415686         GSE66999  GPL6480
# 36  TRUE 61359   61359         1 61359   61359  1.0000000          GSE6861  GPL1352
# 37 FALSE 22215   22215         1 22277   22215  0.9972169         GSE69031   GPL571
# 38 FALSE 44810   44810         1 49576   44810  0.9038648         GSE76360  GPL6947


# id mapping
load("data/ncbi_gene_annot/Homo_sapiens.gene_info.RData")
annot <- Homo_sapiens.gene_info %>%
  dplyr::rename(Ncbi_gene_id = "GeneID",
                Hugo_gene_symbol = "Symbol",
                Ensembl_gene_id = "EnsemblGeneID") %>%
  dplyr::mutate(Ncbi_gene_id = as.character(Ncbi_gene_id))



# Checking if atleast Ncbi_gene_id, Hugo_gene_symbol, or Ensembl_gene_id
# is present in every annotation
table(purrr::map_lgl(
  gpl,
  ~(any(c("Entrez_gene_id", "Hugo_gene_symbol", "Ensembl_gene_id") %in% names(.x))))
)
# TRUE
# 19
# i.e all annotation contains atleast one of the tested columns.



# Clean annotations gpl
# >>>>>>>>>>>>>>>>>>>>>

gpl <- purrr::map(
  gpl,
  function(x, annot){

    if ("Entrez_gene_id" %in% names(x)) {

      x <- x %>%
        dplyr::filter( !(is.na(Entrez_gene_id) | str_detect(Entrez_gene_id, "///") | str_detect(Entrez_gene_id, "\\|") | str_detect(Entrez_gene_id, ";") | str_detect(Entrez_gene_id, ",")) ) %>%
        dplyr::rename(Ncbi_gene_id = "Entrez_gene_id",
                      ID_REF = "Id_ref") %>%
        dplyr::select(ID_REF, Ncbi_gene_id) %>%
        left_join(annot, by = "Ncbi_gene_id")

    } else if ("Ensembl_gene_id" %in% names(x)) {

      x <- x %>%
        dplyr::filter( !(is.na(Ensembl_gene_id) | str_detect(Ensembl_gene_id, "///") | str_detect(Ensembl_gene_id, "\\|") | str_detect(Ensembl_gene_id, ";") | str_detect(Ensembl_gene_id, ",")) ) %>%
        dplyr::rename(ID_REF = "Id_ref") %>%
        dplyr::select(ID_REF, Ensembl_gene_id) %>%
        left_join(annot, by = "Ensembl_gene_id")

    } else if ("Hugo_gene_symbol" %in% names(x)) {


      x <- x %>%
        dplyr::filter( !(is.na(Hugo_gene_symbol) | str_detect(Hugo_gene_symbol, "///") | str_detect(Hugo_gene_symbol, "\\|") | str_detect(Hugo_gene_symbol, ";") | str_detect(Hugo_gene_symbol, ",")) ) %>%
        dplyr::rename(ID_REF = "Id_ref") %>%
        dplyr::select(ID_REF, Hugo_gene_symbol) %>%
        left_join(annot, by = "Hugo_gene_symbol")

    }


    x %>%
      dplyr::filter( !(is.na(Ncbi_gene_id) | str_detect(Ncbi_gene_id, "///") | str_detect(Ncbi_gene_id, "\\|") | str_detect(Ncbi_gene_id, ";") | str_detect(Ncbi_gene_id, ",")) ) %>%
      dplyr::filter( !(is.na(Ensembl_gene_id) | str_detect(Ensembl_gene_id, "///") | str_detect(Ensembl_gene_id, "\\|") | str_detect(Ensembl_gene_id, ";") | str_detect(Ensembl_gene_id, ",")) ) %>%
      dplyr::filter( !(is.na(Hugo_gene_symbol) | str_detect(Hugo_gene_symbol, "///") | str_detect(Hugo_gene_symbol, "\\|") | str_detect(Hugo_gene_symbol, ";") | str_detect(Hugo_gene_symbol, ",")) ) %>%
      dplyr::select(ID_REF,Ncbi_gene_id, Hugo_gene_symbol, Ensembl_gene_id)


  },
  annot
)


# columns in gpl
purrr::map(gpl, names) %>% unlist() %>% table()
# Ensembl_gene_id Hugo_gene_symbol           ID_REF     Ncbi_gene_id
#              19               19               19               19


# Exploring platform size(No.of probes/genes)
xx <- data.frame(
  gpl = names(gpl),
  Size = purrr::map_int(gpl, nrow),
  Unique_probes = purrr::map_int(gpl, ~(length(unique(.x$ID_REF)))),
  Unique_genes = purrr::map_int(gpl, ~(length(unique(.x$Ncbi_gene_id))))
)
xx
#   gpl       Size Unique_probes Unique_genes
# 1 GPL1352  40817         40817        18076
# 2 GPL1390  10234         10234         9924
# 3 GPL14374  8737          8737         8623
# 4 GPL14550 27630         27630        19917
# 5 GPL14668 17912         17912        11483
# 6 GPL1708  27531         27531        17650
# 7 GPL18649   559           559          559
# 8 GPL24546   791           791          786
# 9 GPL24993   733           733          733
# 10 GPL4133  30844         30844        18272
# 11 GPL5325  13537         13537         9992
# 12 GPL570   39766         39766        18943
# 13 GPL571   19334         19334        12147
# 14 GPL6098  11777         11777         9601
# 15 GPL6480  29504         29504        18447
# 16 GPL6884  29098         29098        19179
# 17 GPL6947  29072         29072        18925
# 18 GPL7504  34633         34633        17056
# 19 GPL96    19334         19334        12147

gse_gpl_cong %>%
  left_join(xx, by = "gpl") %>%
  as.data.frame()
#     n_eq  n_ex n_ex2an prop_ex2n  n_an n_an2ex prop_an2ex              gse      gpl  Size Unique_probes Unique_genes
# 1  FALSE   815     815         1   819     815  0.9951160        GSE109710 GPL24546   791           791          786 !!!
# 2  FALSE   770     770         1   788     770  0.9771574        GSE114403 GPL24993   733           733          733 !!!
# 3  FALSE 36286   36286         1 41108   36286  0.8826992        GSE130786  GPL6480 29504         29504        18447
# 4  FALSE   579     579         1   598     579  0.9682274        GSE143222 GPL18649   559           559          559 !!!
# 5   TRUE 42545   42545         1 42545   42545  1.0000000        GSE143846 GPL14550 27630         27630        19917
# 6   TRUE 54675   54675         1 54675   54675  1.0000000         GSE16391   GPL570 39766         39766        18943
# 7   TRUE 54675   54675         1 54675   54675  1.0000000         GSE16446   GPL570 39766         39766        18943
# 8   TRUE 54675   54675         1 54675   54675  1.0000000         GSE18728   GPL570 39766         39766        18943
# 9   TRUE 54675   54675         1 54675   54675  1.0000000         GSE19615   GPL570 39766         39766        18943
# 10  TRUE 22283   22283         1 22283   22283  1.0000000         GSE20194    GPL96 19334         19334        12147
# 11  TRUE 22283   22283         1 22283   22283  1.0000000         GSE20271    GPL96 19334         19334        12147
# 12 FALSE 54627   54627         1 54675   54627  0.9991221         GSE20685   GPL570 39766         39766        18943
# 13 FALSE 41000   41000         1 41108   41000  0.9973728         GSE21974  GPL6480 29504         29504        18447
# 14 FALSE 21631   21631         1 22579   21631  0.9580141 GSE21997_GPL1390  GPL1390 10234         10234         9924 !!!
# 15 FALSE 43898   43898         1 44294   43898  0.9910597 GSE21997_GPL5325  GPL5325 32443         32443        16614
# 16 FALSE 44940   44940         1 45220   44940  0.9938080 GSE21997_GPL7504  GPL7504 34633         34633        17056
# 17  TRUE 22283   22283         1 22283   22283  1.0000000         GSE22093    GPL96 19334         19334        12147
# 18 FALSE 24332   24332         1 24385   24332  0.9978265         GSE22219  GPL6098 11777         11777         9601 !!!
# 19 FALSE 28349   28349         1 44290   28349  0.6400768 GSE22226_GPL1708  GPL1708 27531         27531        17650
# 20 FALSE 44795   44795         1 45220   44795  0.9906015 GSE22226_GPL4133  GPL4133 30844         30844        18272
# 21 FALSE 41810   41810         1 44294   41810  0.9439202         GSE22358  GPL5325 32443         32443        16614
# 22  TRUE 22283   22283         1 22283   22283  1.0000000         GSE23988    GPL96 19334         19334        12147
# 23  TRUE 22283   22283         1 22283   22283  1.0000000         GSE25066    GPL96 19334         19334        12147
# 24  TRUE 54675   54675         1 54675   54675  1.0000000         GSE28844   GPL570 39766         39766        18943
# 25 FALSE 11851   11851         1 26823   11851  0.4418223         GSE31863 GPL14374 14841         14841        13788
# 26 FALSE 12901   12901         1 35073   12901  0.3678328         GSE32603 GPL14668 17912         17912        11483
# 27 FALSE 54613   54613         1 54675   54613  0.9988660         GSE32646   GPL570 39766         39766        18943
# 28 FALSE 27506   27506         1 48803   27506  0.5636129         GSE34138  GPL6884 29098         29098        19179
# 29  TRUE 22277   22277         1 22277   22277  1.0000000         GSE41998   GPL571 19334         19334        12147
# 30  TRUE 22283   22283         1 22283   22283  1.0000000         GSE42822    GPL96 19334         19334        12147
# 31 FALSE 22268   22268         1 22283   22268  0.9993268         GSE45255    GPL96 19334         19334        12147
# 32 FALSE 61297   61297         1 61359   61297  0.9989896          GSE4779  GPL1352 40817         40817        18076
# 33  TRUE 54675   54675         1 54675   54675  1.0000000         GSE50948   GPL570 39766         39766        18943
# 34  TRUE 54675   54675         1 54675   54675  1.0000000         GSE66305   GPL570 39766         39766        18943
# 35 FALSE 38706   38706         1 41108   38706  0.9415686         GSE66999  GPL6480 29504         29504        18447
# 36  TRUE 61359   61359         1 61359   61359  1.0000000          GSE6861  GPL1352 40817         40817        18076
# 37 FALSE 22215   22215         1 22277   22215  0.9972169         GSE69031   GPL571 19334         19334        12147
# 38 FALSE 44810   44810         1 49576   44810  0.9038648         GSE76360  GPL6947 29072         29072        18925




# Save gpl !!!!!!!!
#
# save(gpl, file = str_c(out_data, "gpl.RData"))
# load(str_c(out_data, "gpl.RData"))




# Maxvar collapseing
# >>>>>>>>>>>>>>>>>>>

# geo: selected datasets
# gse_gpl_map: series matrix to annotation dtaframe mapping
# gpl: cleaned gpl

geo_maxvar_annot <- purrr::map(
  names(geo),
  function(xnme_gse, xgse_list, xgse_gpl_map, xgpl_list){

    xnme_gpl <- xgse_gpl_map %>%
      dplyr::filter(Series_matrix_accession == xnme_gse) %>%
      dplyr::select(Series_matrix_platform_id) %>%
      tibble::deframe()

    xgse <- xgse_list[[xnme_gse]]$expression
    xgpl <- xgpl_list[[xnme_gpl]]

    print(str_c(xnme_gse, " - ", xnme_gpl,"; ", class(xgse$ID_REF), class(xgpl$ID_REF), class(xgpl$Ncbi_gene_id)))

    get_max_var_annot(ge = xgse, xannot = xgpl,
                      pset_colname = "ID_REF", gene_colname = "Ncbi_gene_id")

  },
  xgse_list = geo,
  xgse_gpl_map = gse_gpl_map,
  xgpl_list = gpl
)
names(geo_maxvar_annot) <- names(geo)


# Save RData !!!!!!!
#
# save(geo_maxvar_annot, file = str_c(out_data, "geo_maxvar_annot.RData"))
# load(str_c(out_data, "geo_maxvar_annot.RData"))



# Checking congruence between geo maxvar-expression matrix and gpl platform
gse_gpl_cong2 <- purrr::map2_dfr(
  gse_gpl_map$Series_matrix_accession,
  gse_gpl_map$Series_matrix_platform_id,
  function(nme_gse, nme_gpl, geo, gpl, geo_maxvar_annot){

    xgse <- geo[[nme_gse]]$expression
    xgpl <- gpl[[nme_gpl]]
    xannot <-geo_maxvar_annot[[nme_gse]]

    list(
      #ex:expression, an:annotation
      n_ex =nrow(xgse),
      n_exmv = nrow(xannot),
      prop_exmv2ex = sum(xannot$ID_REF %in% xgse$ID_REF)/nrow(xannot) %>% round(digits = 2),
      n_an = nrow(xgpl),
      n_anuq = xgpl$Ncbi_gene_id %>% unique() %>% length(),
      n_anup2exmv = sum((xgpl$Ncbi_gene_id %>% unique()) %in% xannot$Ncbi_gene_id),
      prop_anup2exmv = sum((xgpl$Ncbi_gene_id %>% unique()) %in% xannot$Ncbi_gene_id) / (xgpl$Ncbi_gene_id %>% unique() %>% length()) %>% round(digits = 2),
      gse = nme_gse,
      gpl = nme_gpl
    )
  },
  geo,
  gpl,
  geo_maxvar_annot
)

gse_gpl_cong2 %>% as.data.frame()
# Expectation:
# prop_exmv2ex muat be 1, ie. all maxvar genes must present in full expr matrix
# prop_anup2exmv close to 1, ie. proportion of unique genes in gpl present in gse.
#  Note1: genes are removed if any sample has got NA expression,
#  hence prop_anup2exmv can be less han 1. More missingness could be a problem.
#  Note2: gse34138(GPL6884) probes are missing in the original series matrix.

#     n_ex n_exmv prop_exmv2ex  n_an n_anuq n_anup2exmv prop_anup2exmv              gse      gpl
# 1    815    786            1   791    786         786      1.0000000        GSE109710 GPL24546 ! nanostring ! discarded, see below
# 2    770    733            1   733    733         733      1.0000000        GSE114403 GPL24993 ! nanostring ! discarded, see below
# 3  36286  17146            1 29504  18447       17146      0.9294736        GSE130786  GPL6480
# 4    579    544            1   559    559         544      0.9731664        GSE143222 GPL18649 ! nanostring ! discarded, see below
# 5  42545  19917            1 27630  19917       19917      1.0000000        GSE143846 GPL14550
# 6  54675  18943            1 39766  18943       18943      1.0000000         GSE16391   GPL570
# 7  54675  18943            1 39766  18943       18943      1.0000000         GSE16446   GPL570
# 8  54675  18943            1 39766  18943       18943      1.0000000         GSE18728   GPL570
# 9  54675  18943            1 39766  18943       18943      1.0000000         GSE19615   GPL570
# 10 22283  12147            1 19334  12147       12147      1.0000000         GSE20194    GPL96
# 11 22283  12147            1 19334  12147       12147      1.0000000         GSE20271    GPL96
# 12 54627  18943            1 39766  18943       18943      1.0000000         GSE20685   GPL570
# 13 41000  18447            1 29504  18447       18447      1.0000000         GSE21974  GPL6480
# 14 21631   9560            1 10234   9924        9560      0.9633212 GSE21997_GPL1390  GPL1390 ! discarded, see below
# 15 43898  16562            1 32443  16614       16562      0.9968701 GSE21997_GPL5325  GPL5325
# 16 44940  17016            1 34633  17056       17016      0.9976548 GSE21997_GPL7504  GPL7504
# 17 22283  12147            1 19334  12147       12147      1.0000000         GSE22093    GPL96
# 18 24332   9593            1 11777   9601        9593      0.9991668         GSE22219  GPL6098 ! discarded, see below
# 19 28349  13042            1 27531  17650       13042      0.7389235 GSE22226_GPL1708  GPL1708 ! moderate missingness ! discarded, see below
# 20 44795  18182            1 30844  18272       18182      0.9950744 GSE22226_GPL4133  GPL4133
# 21 41810  16193            1 32443  16614       16193      0.9746599         GSE22358  GPL5325
# 22 22283  12147            1 19334  12147       12147      1.0000000         GSE23988    GPL96
# 23 22283  12147            1 19334  12147       12147      1.0000000         GSE25066    GPL96
# 24 54675  18943            1 39766  18943       18943      1.0000000         GSE28844   GPL570
# 25 11851   7263            1 14841  13788        7263      0.5267624         GSE31863 GPL14374 ! lagre missingness ! discarded, see below
# 26 12901   5934            1 17912  11483        5934      0.5167639         GSE32603 GPL14668 ! lagre missingness ! discarded, see below
# 27 54613  18943            1 39766  18943       18943      1.0000000         GSE32646   GPL570
# 28 27506  16071            1 29098  19179       16071      0.8379478         GSE34138  GPL6884 ! moderate missingness ! discarded, see below
# 29 22277  12147            1 19334  12147       12147      1.0000000         GSE41998   GPL571
# 30 22283  12147            1 19334  12147       12147      1.0000000         GSE42822    GPL96
# 31 22268  12147            1 19334  12147       12147      1.0000000         GSE45255    GPL96
# 32 61297  18076            1 40817  18076       18076      1.0000000          GSE4779  GPL1352
# 33 54675  18943            1 39766  18943       18943      1.0000000         GSE50948   GPL570
# 34 54675  18943            1 39766  18943       18943      1.0000000         GSE66305   GPL570
# 35 38706  17829            1 29504  18447       17829      0.9664986         GSE66999  GPL6480
# 36 61359  18076            1 40817  18076       18076      1.0000000          GSE6861  GPL1352
# 37 22215  12147            1 19334  12147       12147      1.0000000         GSE69031   GPL571
# 38 44810  17885            1 29072  18925       17885      0.9450462         GSE76360  GPL6947

identical(names(geo), names(geo_maxvar_annot)) # TRUE

#
# ==============================================================================




# 5. Genrate geo_tidy, a clean version of geo with only expression and clinical
# slots to hold max-var collaped expression data and cleaned clinical data.
# ==============================================================================

# maxvar expression matrix list
geo_tidy <- purrr::map2(
  geo,
  geo_maxvar_annot,
  function(x, y){

    y1 <- y %>%
      dplyr::mutate_all(~(as.character(.x))) %>%
      dplyr::left_join(
        x$expression %>% dplyr::mutate_at("ID_REF",~(as.character(.x))),
        by = "ID_REF") %>%
      dplyr::select(-c("ID_REF","Hugo_gene_symbol","Ensembl_gene_id")) %>%
      dplyr::mutate(
        Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id)
      )
    print(paste(nrow(x$expression), nrow(y), identical(y1$Ncbi_gene_id %>% str_replace("ncbi_", ""),  y$Ncbi_gene_id)))

    list(expression = y1, clinical = NA)
  }
)


# Save RData !!!!!!!!!!
#
# save(geo_tidy, file = str_c(out_data, "geo_tidy.RData"))
# load(str_c(out_data, "geo_tidy.RData"))

#
# ==============================================================================




# 6. Expression data integration
# ==============================================================================

# Nanostring series
#     n_ex n_exmv prop_exmv2ex  n_an n_anuq n_anup2exmv prop_anup2exmv              gse      gpl
# 1    815    786            1   791    786         786      1.0000000        GSE109710 GPL24546
# 2    770    733            1   733    733         733      1.0000000        GSE114403 GPL24993
# 4    579    544            1   559    559         544      0.9731664        GSE143222 GPL18649


# Series with large missingness
# 25 11851   7263            1 14841  13788        7263      0.5267624         GSE31863 GPL14374 ! lagre missingness
# 26 12901   5934            1 17912  11483        5934      0.5167639         GSE32603 GPL14668 ! lagre missingness
# "GPL14374" ~ "Swegene///SWEGENE Human_v2.1.1 55K_condensed///26819"
# "GPL14668" ~ "UCSF_DeRisi_Lab///UCSF/HAQQLAB_Human_40986_ISPY///35069"



# Ideintifying common genes
# >>>>>>>>>>>>>>>>>>>>>>>>>>

# Common genes in all, except nanostring series and low quality datasets
# Low quality is defined as poor representation of transcriptome
# due to gene missingness.
uid = geo_maxvar_annot$GSE130786$Ncbi_gene_id # non-naostring series
for(inme in names(geo_maxvar_annot)){
  if(!(inme %in% c("GSE109710", "GSE114403", "GSE143222", # Nanostring
                   "GSE21997_GPL1390", # nearly halved unique genes; Agilent (2005 array, 22575 probes)
                   "GSE22219", # nearly halved unique genes; GPL6098 Illumina (2007 array, 24385 probes)
                   "GSE22226_GPL1708", # moderate missingness; Agilent (2004 array, 44290 probes)
                   "GSE31863", # large missingness; GPL14374 Swegene (2011 array, 26819 probes)
                   "GSE32603", # large missingness; GPL14668 UCSF_DeRisi_Lab (2012 array, 35069 probes)
                   "GSE34138" # small missingness; GPL6884 Illumina (2008 array, 48803 probes)
  )))
    uid <- intersect(uid, geo_maxvar_annot[[inme]]$Ncbi_gene_id)
  print(paste(inme, length(uid)))
}
uid_geo = uid #  9184

# [1] "GSE109710 17146"
# [1] "GSE114403 17146"
# [1] "GSE130786 17146"
# [1] "GSE143222 17146"
# [1] "GSE143846 16897"
# [1] "GSE16391 15925"
# [1] "GSE16446 15925"
# [1] "GSE18728 15925"
# [1] "GSE19615 15925"
# [1] "GSE20194 11097"
# [1] "GSE20271 11097"
# [1] "GSE20685 11097"
# [1] "GSE21974 11097"
# [1] "GSE21997_GPL1390 11097"
# [1] "GSE21997_GPL5325 10586"
# [1] "GSE21997_GPL7504 10402"
# [1] "GSE22093 10402"
# [1] "GSE22219 10402"
# [1] "GSE22226_GPL1708 10402"
# [1] "GSE22226_GPL4133 10351"
# [1] "GSE22358 10150"
# [1] "GSE23988 10150"
# [1] "GSE25066 10150"
# [1] "GSE28844 10150"
# [1] "GSE31863 10150"
# [1] "GSE32603 10150"
# [1] "GSE32646 10150"
# [1] "GSE34138 10150"
# [1] "GSE41998 10150"
# [1] "GSE42822 10150"
# [1] "GSE45255 10150"
# [1] "GSE4779 10049"
# [1] "GSE50948 10049"
# [1] "GSE66305 10049"
# [1] "GSE66999 9798"
# [1] "GSE6861 9798"
# [1] "GSE69031 9798"
# [1] "GSE76360 9184"



# Does er/pr/her2/pam50 genes present in uid_geo ?
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# ESR1 "2099" %in% uid_geo # TRUE
# ERBB2 "2064" %in% uid_geo # TRUE
# PGR "5241" %in% uid_geo # TRUE

# PAM5050 genes
pam50 <- read_tsv("data/pam50/PAM50/bioclassifier_R/pam50_annotation.txt")
table(pam50$EntrezGene %>% as.character() %in% uid_geo)
# FALSE  TRUE
# 14    36




# Consolidate geo expression
# >>>>>>>>>>>>>>>>>>>>>>>>>>

# Filter sub standard datasets >>>>>>>>>
nme <- c("GSE109710", "GSE114403", "GSE143222", # Nanostring
         "GSE21997_GPL1390", # nearly halved unique genes; Agilent
         "GSE22219", # nearly halved unique genes; GPL6098 Illumina
         "GSE22226_GPL1708", # moderate missingness; Agilent
         "GSE31863", # large missingness; GPL14374 Swegene
         "GSE32603", # large missingness; GPL14668 UCSF_DeRisi_Lab
         "GSE34138" # small missingness; GPL6884 Illumina
)

# "GPL24546" ~ "Nanostring///NanoString nCounter human Gene Expression Custom CodeSet///815",
# "GPL24993" ~ "Nanostring///Nanostring Immune Oncology 360 platform///784",
# "GPL18649" ~ "Nanostring///NanoString nCounter GX Human Immunology v2///594",
#
# "GPL1390" ~ "Agilent///Agilent Human 1A Oligo UNC custom Microarrays///22575",
# "GPL6098" ~ "Illumina///Illumina humanRef-8 v1.0 expression beadchip///24385"
# "GPL1708" ~ "Agilent///Agilent-012391 Whole Human Genome Oligo Microarray G4112A (Feature Number version)///44290"
# "GPL14374" ~ "Swegene///SWEGENE Human_v2.1.1 55K_condensed///26819",
# "GPL14668" ~ "UCSF_DeRisi_Lab///UCSF/HAQQLAB_Human_40986_ISPY///35069",
# "GPL6884" ~ "Illumina///Illumina HumanWG-6 v3.0 expression beadchip///48803"



geo_expr <- purrr::map(
  names(geo_tidy)[!(names(geo_tidy) %in% nme)],
  function(nme_gse, geo_tidy, genes){
    print(nme_gse)
    geo_tidy[[nme_gse]]$expression %>%
      dplyr::filter(Ncbi_gene_id %in% all_of(genes)) %>%
      t_tibble(names_x_desc = "Sample_geo_accession") %>% #to facilitate merging
      dplyr::mutate(Series_matrix_accession = nme_gse) %>%
      dplyr::select("Series_matrix_accession", "Sample_geo_accession", all_of(genes))

  },
  geo_tidy,
  genes = str_c("ncbi_", uid_geo)
)
# list

# Checking whther samples per dataset are unique
table(purrr::map_lgl(
  geo_expr,~(nrow(.x) == (.x$Sample_geo_accession %>% unique() %>% length()))
)
)
# TRUE
# 29
# Per dataset samples are unique


# integration
geo_expr <- bind_rows(geo_expr) # merging dataframes
head2(geo_expr)
# $head
#   Series_matrix_accession Sample_geo_accession ncbi_400451 ncbi_57099 ncbi_23339
# 1 GSE130786               GSM3753582                0.340      0.238     -0.102
# 2 GSE130786               GSM3753583               -0.262      0.0634    -0.204
# 3 GSE130786               GSM3753584                0.0394     0.135      0.160
# 4 GSE130786               GSM3753585               -0.571      0.109     -0.261
# 5 GSE130786               GSM3753586               -0.179      0.135     -0.0628
# $dim
# [1] 3744 9186; (samples x genes)




# Checking whther samples in pooled dataset are unique
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

nrow(geo_expr) == (geo_expr$Sample_geo_accession %>% unique() %>% length()) # FALSE
nrow(geo_expr) # 3744
geo_expr$Sample_geo_accession %>% unique() %>% length() # 3740

# Sample redundancy @ Sample_geo_accession !!!!!!!!!!!!!!!!!!
# Discard all samples from geo_expr for which redundancy exists !!!!!!!!


idx <- purrr::map(
  geo_expr$Sample_geo_accession %>% unique(),
  function(x, id){
    which(x == id)
  },
  id = geo_expr$Sample_geo_accession
)
idx[purrr::map_int(idx,length) >1]
# [[1]]
# [1] 1374 1686
# [[2]]
# [1] 1375 1687
# [[3]]
# [1] 1376 1688
# [[4]]
# [1] 1377 1689


# Exploring expression values of redundant samples !!!!!!!
geo_expr[idx[purrr::map_int(idx,length) >1] %>% unlist(),]
#   Series_matrix_a… Sample_geo_acce… ncbi_400451 ncbi_57099 ncbi_23339 ncbi_11261 ncbi_57082
# 1 GSE21997_GPL5325 GSM556167             -0.663      0.357      0.008      0.588      -4.63
# 2 GSE22358         GSM556167             -0.663      0.357      0.008      0.588      -4.63
# 3 GSE21997_GPL5325 GSM556168              1.24      -0.199     -0.228      0.188      -3.54
# 4 GSE22358         GSM556168              1.24      -0.199     -0.228      0.188      -3.54
# 5 GSE21997_GPL5325 GSM556169              0.131      0.072      0.103      0.269      -3.63
# 6 GSE22358         GSM556169              0.131      0.072      0.103      0.269      -3.63
# 7 GSE21997_GPL5325 GSM556170              0.328      0.336     -0.068      0.523      -2.97
# 8 GSE22358         GSM556170              0.328      0.336     -0.068      0.523      -2.97


# Exploring Arms of redundat samples !!!!!!!
geo$GSE21997_GPL5325$clinical$Arm_neoadj %>% unique()
# [1] "Docetaxel"   "Doxorubicin" NA
geo$GSE22358$clinical$Arm_neoadj %>% unique()
# [1] "Docetaxel+Capecitabine"             "Docetaxel+Capecitabine+Trastuzumab"
# [3] NA
# No Doxorubicine info in GSE22358 !!!!!
# Incinsistant arms !!!!!!!!!!!, discard all redundat samples


# Redundant sample filtering
# Note: all redundant samples are filtered as it is unsure which arms are correct
dim(geo_expr) # [1] 3744 9186
geo_expr <- geo_expr[-(idx[purrr::map_int(idx,length) >1] %>% unlist()),]
dim(geo_expr) # [1] 3736 9186


# Note:
# Some studies may have used same samples and may exist in GEO with different
# GSM sample accessions. These samples could be identified using co-correlation
# matrix, as redundant samples might have highly correlated.
# But we are ignoring this step, as we don't expect it may change results.
# However the above samples are removed because of identical Sample_geo_accession.



# Wide to long format ((samples x genes) to (genes x samples))
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Mapping b/w series matrix and sample
geo_expr_meta <- geo_expr %>%
  dplyr::select(Series_matrix_accession, Sample_geo_accession)

# Gene X Sample integrated expression matrix
geo_expr <- t_tibble(geo_expr %>% dplyr::select(-Series_matrix_accession),
                     names_x_desc = "Ncbi_gene_id")
head2(geo_expr)
# $head
#   Ncbi_gene_id GSM3753582 GSM3753583 GSM3753584 GSM3753585
# 1 ncbi_400451      0.340     -0.262      0.0394     -0.571
# 2 ncbi_57099       0.238      0.0634     0.135       0.109
# 3 ncbi_23339      -0.102     -0.204      0.160      -0.261
# 4 ncbi_11261      -0.0819    -0.0967     0.156      -0.163
# 5 ncbi_57082      -0.259     -0.0952    -0.145      -0.257
# $dim
# [1] 9184 3737; genes x samples


# Save RData !!!!!
# save(geo_expr, file = str_c(out_data, "geo_expr.RData"))
# save(geo_expr_meta, file = str_c(out_data, "geo_expr_meta.RData"))


# Write_out annotation data
# >>>>>>>>>>>>>>>>>>>>>>>>>

x <- geo_expr %>%
  dplyr::select(Ncbi_gene_id) %>%
  dplyr::mutate(Ncbi_gene_id = str_replace(Ncbi_gene_id, "ncbi_", ""))
x <- x %>%
  dplyr::left_join(annot, by = "Ncbi_gene_id") %>%
  dplyr::mutate(ID_REF = str_c("ncbi_", Ncbi_gene_id)) %>%
  dplyr::select(ID_REF, Ncbi_gene_id, Hugo_gene_symbol, Ensembl_gene_id)

identical(x$ID_REF, geo_expr$Ncbi_gene_id) # TRUE

write_tsv(
  x = x,
  path = str_c(out_data, "geo_expr_annot.tsv")
)

#
# ==============================================================================

