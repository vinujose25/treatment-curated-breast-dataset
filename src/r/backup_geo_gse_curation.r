# geo.R


# Directory settings
# ==============================================================================
#
# Working directory
# setwd("/home/josev/proj/on/TIL.IMM.FIB/")

# Output directory
if(!dir.exists("results/geo_gse_curation")){
  dir.create(path = "results/geo_gse_curation", recursive = TRUE, mode = "0700")
}
outdir <- "results/geo_gse_curation/"
#
# ==============================================================================




# Public library
# ==============================================================================
#
library(tidyverse)
library(GEOquery)
library(sva)
library(readxl)
#
# ==============================================================================




# Local library
# ==============================================================================
#
source("R/geo_gse_curation_fun.R")
#
# ==============================================================================





# Load and format geo data
# =============================================================================
#

files = c(
  list.files("data/geo_data3/gse", full.names = TRUE)
  )
geo <- format_geo(files = files)

names(geo) <-  str_split_fixed(names(geo),"/", 4)[, 4] %>%
  str_replace("_series_matrix.txt.gz", "") %>%
  str_replace("-", "_")


# # # >>>>>>>>>>>>>>>>>>>>>>>
# # GSE8465_GPL1390
# # GSE8465_GPL887
# # GSE143846
# xfiles = c(
#   list.files("data/geo_data3/gse", full.names = TRUE)
# )
# xfiles <- xfiles[c(str_which(xfiles, "GSE8465"),str_which(xfiles, "GSE143846"))]
# xgeo <- format_geo(files = xfiles)
#
# names(xgeo) <-  str_split_fixed(names(xgeo),"/", 4)[, 4] %>%
#   str_replace("_series_matrix.txt.gz", "") %>%
#   str_replace("-", "_")
# xgeo$GSE8465_GPL1390$dataset_issues # Yes
# xgeo$GSE8465_GPL887$dataset_issues # Yes
# xgeo$GSE143846$dataset_issues # No
# glimpse(xgeo$GSE8465_GPL1390$sample_characteristics)
# glimpse(xgeo$GSE8465_GPL887$sample_characteristics)
# glimpse(xgeo$GSE143846$sample_characteristics)
# # # >>>>>>>>>>>>>>>>>>>>>>>>>>

# Save RData !!!!!!!!
# save(geo, file = str_c(outdir,"geo.RData"))
# load(str_c(outdir,"geo.RData"))



# Backup notes: !!!!!!!!!!!!!!!!!!!!!!!!!

# Unexpected behaviours /  bugs in format_geo()
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
# End of backup notes !!!!!!!!!!!!!!!!!!!!!!!!!!
#
# ==============================================================================





# Manual correction of var_qc_failed_sample_characteristics#
# ==============================================================================

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

# For the below series, already processed data present in
# "results/geo/var_qc_failed_sample_characteristics" is used

# "GSE20194" "GSE20271"           "GSE41998"           "GSE50948"
# "GSE6861"  "GSE25066"


for(inme in nme){
  write_tsv(
    geo[[inme]]$sample,
    path = str_c(outdir, inme, "_sample.tsv")
  )

  write_tsv(
    geo[[inme]]$sample_characteristics,
    path = str_c(outdir, inme, "_sample_characteristics.tsv")
  )
}
# The above files are kept in results/geo_gse_curation/var_qc_failed_sample_characteristics/
# The corrected files have the suffix "_cleaned"

# Manual cleaning performed:
# 1. Changing charcteristic names (column names)
# 2. Correcting mis-alinged samples_characterisitcs spreadsheet-cells due to
# missingness.
# 3. Removing the concatenated characterisitcs name to each value.



# Append cleaned sample characteristics to geo
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
files <- list.files("results/geo_gse_curation/var_qc_failed_sample_characteristics", full.names = TRUE)
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

# geo$GSE8465_GPL1390$sample_characteristics
# geo$GSE8465_GPL887$sample_characteristics
# geo$GSE143846$sample_characteristics

# Save RData !!!!!!!!
# save(geo, file = str_c(outdir,"geo.RData"))
# load(str_c(outdir,"geo.RData"))
#
# ==============================================================================





# Manually rename and format columns of sample_characteristics#
# ==============================================================================


# all sample characteristic names
xx <- purrr::map(geo, function(x){names(x$sample_characteristics)}) %>%
  unlist()
names(xx) <- NULL
xx <- sort(xx) %>% tibble::enframe(name = NULL)
write_tsv(xx, str_c(outdir, "all_column_names.txt"))


# purrr::map(geo,~((.x$sample %>% names())))
c("Sample_title", "Sample_geo_accession",
  "Sample_type", "Sample_channel_count",
  "Sample_source_name_ch1", "Sample_organism_ch1",
  "Sample_platform_id")
# inme = "GSE32603"
# inme = "GSE31863"
# inme = "GSE8465_GPL1390"
# inme = "GSE8465_GPL887"
# inme = "GSE143846"


for(inme in names(geo)){

  igeo <- geo[[inme]]

  if (unique(igeo$sample$Sample_channel_count) == 1) {
    colnme <- c("Sample_title", "Sample_geo_accession",
                "Sample_type", "Sample_channel_count",
                "Sample_source_name_ch1", "Sample_organism_ch1",
                "Sample_platform_id")

  } else {
    colnme <- c("Sample_title", "Sample_geo_accession",
                "Sample_type", "Sample_channel_count",
                "Sample_source_name_ch1", "Sample_organism_ch1",
                "Sample_source_name_ch2", "Sample_organism_ch2",
                "Sample_platform_id")

  }

  xx <- bind_cols(
    igeo$sample[, colnme],
    # Append the above columns from sapme to sample _characteristics
    # This will make complete clinical data which cna be mapped to expression data
    igeo$sample_characteristics
  )

  write_tsv(x = xx, path = str_c(outdir,inme,".tsv"))
}
# The above files are stored in results/geo/reformatted_clinical_data/
# The reformatted files have the name suffix, "_cleaned.tsv"



# Reformatting performed:
# 1. Keep all original sample_characterisitcs
# 2. Extract hidden sample characterisitcs from text and append columns.
# 3. Whereever possible stick to the following format
#    Er. Pr, Her2, Node: pos vs neg
#    Age: years
#    Grade: 1-3
#    Size: cm
#    Histology: as is
#    Ethinicity: as is
#    Time_os/rfs/ddfs/organ_specific_metastasis: years
#    Event_os/rfs/ddfs/organ_specific_metastasis event: 1/0 event/noevent
#
# Note !!!!!!!!!!!!!!!
# Take dfs: disease free survival as the merged endpoint
# Within each merged dataset map the most releveant endpoint to dfs
# Eg. map dfs, dmfs, mfs, rfs to dfs
#
#    Response: pCR/nopCR (neoadj endpointr)
#    Arm: radio_chemo_hormone_her2therapy (Treatment arm)
#    Arm_description: therapy descriptoin
#    Sampe_id1/2/3: if additional sample ids are present other than Sample_title
#
#    Additional potential columns:
#    Node_collected: no.of nodes ivestigated
#    Node_pos: no.of nodes positive
#    Ki67: 0-100
#    Archive_method: ffpe, freezing, rna_later
#    Til: percentage (average multiple counts)
#    Leukocyte_infiltration: percentage
#    Regimen: neoadj/adj
#    Therapy_radio: yes, no
#    Therapy_hormone: yes, no
#    Therapy_her2: yes, no
#    Therapy_chemo: yes, no
#    Therapy_hormone_class: therapy given seperated by ///
#    Therapy_her2_class: therapy given seperated by ///
#    Therapy_chemo_class: therapy given seperated by ///
#    Timepoint: for paired data




# List of geo datasets which have processed previously
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nme_old <- list.files("results/geo/reformatted_clinical_data/",
                      full.names = FALSE)
nme_new <- list.files("results/geo_gse_curation/reformatted_clinical_data/",
                      full.names = FALSE)

nme_old <- nme_old[str_detect(nme_old, ".tsv")]
nme_old <- nme_old[!str_detect(nme_old, "cleaned")]

nme_new <- nme_new[str_detect(nme_new, ".tsv")]
nme_new <- nme_new[!str_detect(nme_new, "cleaned")]

length(nme_old) # 48
length(nme_new) # 44

length(intersect(nme_old, nme_new)) # 13

intersect(nme_old, nme_new)
# [1] "GSE114403.tsv" "GSE19615.tsv"  "GSE20194.tsv"  "GSE20271.tsv"  "GSE25066.tsv"
# [6] "GSE32646.tsv"  "GSE4056.tsv"   "GSE41998.tsv"  "GSE45255.tsv"  "GSE4779.tsv"
# [11] "GSE50948.tsv"  "GSE6861.tsv"   "GSE76360.tsv"
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>





# Computational formatting of sample characterisitics of the 40 datasets selected
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# T stage
# TX Primary tumor cannot be assessed
# T0 No evidence of primary tumor
# T1 size < 2cm
# T2 2cm < size < 5cm
# T3 size > 5cm
# T4 Tumor of any size with direct extensionto the chest wall and/or to the skin
# (ulceration or skin nodules).


clean_names <- function(x){
  x %>%
    str_trim() %>%
    str_replace_all(" ", "_") %>%
    str_to_sentence()
}




# GSE25066
# >>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE25066.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx)

table(xx$Grade)
# 1               2               3   4=Indeterminate
# 32             180             259              15
table(xx$Clinical_nodal_status)
# N0  N1  N2  N3
# 157 244  66  41
table(xx$Clinical_t_stage)
# T0  T1  T2  T3  T4
# 3  30 255 145  75
table(xx$Er_status_ihc)
# I   N   P
# 4 205 297
table(xx$Pr_status_ihc)
# I   N   P
# 5 258 243
table(xx$Her2_status)
# I   N   P
# 4 485   6
table(xx$Drfs_1_event_0_censored)
# 0   1
# 397 111
table(xx$Pathologic_response_pcr_rd)
# pCR  RD
# 99 389
table(xx$Tissue)
# breast cancer tumor
# 284
table(xx$Type_taxane)
# Taxol Taxotere
# 106       92

xx %>%
  dplyr::group_by(Type_taxane, Pathologic_response_pcr_rd, Drfs_1_event_0_censored) %>%
  dplyr::summarise(N = n())
#   Type_taxane Pathologic_response_pcr_rd Drfs_1_event_0_censored     N
# 1 Taxol       pCR                                              0    19
# 2 Taxol       RD                                               0    51
# 3 Taxol       RD                                               1    22
# 4 Taxol       NA                                               0    13 # adj
# 5 Taxol       NA                                               1     1 # adj
# 6 Taxotere    pCR                                              0    20 # +FEC
# 7 Taxotere    pCR                                              1     3 # +FEC
# 8 Taxotere    RD                                               0    50 # +FEC
# 9 Taxotere    RD                                               1    17 # +FEC
# 10 Taxotere    NA                                               1     2 # adj
# 11 NA          pCR                                              0    53
# 12 NA          pCR                                              1     4
# 13 NA          RD                                               0   189
# 14 NA          RD                                               1    60
# 15 NA          NA                                               0     2
# 16 NA          NA                                               1     2


xx %>%
  dplyr::group_by(Type_taxane, Source) %>%
  dplyr::summarise(N = n())
#   Type_taxane Source         N
# 1 Taxol       LBJ/IN/GEI    56 ([Paclitaxel+FAC] | [FAC|FEC]) + [Hormone] ???
# 2 Taxol       MDACC         50 ([Paclitaxel+FAC] | [FAC|FEC]) + [Hormone] ???
# 3 Taxotere    LBJ/IN/GEI     2 Docetaxel+Capecitabin+FEC + [Hormone]
# 4 Taxotere    MDACC         36 Docetaxel+Capecitabin+FEC + [Hormone]
# 5 Taxotere    USO           54 Docetaxel+Capecitabin+FEC + [Hormone]
# 6 NA          ISPY          83 Taxane+AC + [Hormone]
# 7 NA          MDACC        227 T+FAC + [Hormone]



xx1 <- xx %>%
  dplyr::mutate(
    Age = Age_years,
    Grade = if_else(Grade == "4=Indeterminate", NA_character_, Grade) %>%
      purrr::map_chr(~(str_c("G", .x))),
    Node_cat = Clinical_nodal_status,
    Node_bin = if_else(Clinical_nodal_status == "N0", "neg", "pos"),
    Size_cat = Clinical_t_stage,
    Size_bin = if_else( str_replace(Clinical_t_stage, "T", "") %>% as.integer() <= 1,
                        "samall", "large"), # small = T0/T1 large = >T1
    Er = if_else(Er_status_ihc == "I", NA_character_, Er_status_ihc) %>%
      purrr::map_chr(~(if_else(.x == "P", "pos", "neg"))),
    Pr = if_else(Pr_status_ihc == "I", NA_character_, Pr_status_ihc) %>%
      purrr::map_chr(~(if_else(.x == "P", "pos", "neg"))),
    Her2 = if_else(Her2_status == "I", NA_character_, Her2_status) %>%
      purrr::map_chr(~(if_else(.x == "P", "pos", "neg"))),
    Event_dfs = Drfs_1_event_0_censored %>% as.integer(),
    Time_dfs = Drfs_even_time_years,
    Response_pathological = if_else(Pathologic_response_pcr_rd == "RD",
                                    "npCR", "pCR"),

    # Arm = purrr::map_chr(
    # Type_taxane,
    # ~(case_when(.x == "Taxol" ~ "Paclitaxel+FAC|FEC",
    #             .x == "Taxotere" ~ "Docetaxel+FEC",
    #             TRUE ~ "Taxane+AC|FAC"))

    # From text:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5638042/
    # "In the discovery cohort, there were 227 FNAs (MDACC) and 83 CBXs (I-SPY), and all chemotherapy was administered as neoadjuvant treatment. In the validation cohort, there were 157 FNAs (MDACC, Peru, USO) and 41 CBXs (MDACC, LBJ, Spain), and 165 of 198 patients received all chemotherapy as neoadjuvant treatment."

    # Adjuvat endocrine!!!!!!!!!
    # "Patients with any nuclear immunostaining for ER in the tumor cells were considered as eligible for adjuvant endocrine therapy."

    # Chemotherapy Regimen		    Discovery Validation
    # Entirely Neoadjuvant        (310)     (198)
    # T × 12 → FAC × 4 → Sx	      227	      73  Taxol
    # AC × 4 → T/Tx × 4 → Sx	    83*	      –
    # TxX × 4 → FEC × 4 → Sx	    –	        92
    #
    # Partial Neoadjuvant
    # FAC/FEC × 6 → Sx → T × 12	  –	        18  Taxol
    #
    # Entirely Adjuvant
    # Sx → T × 12 → FAC/FEC × 4	  –	        12  Taxol
    # Sx → TxX × 4 → FEC × 4	    –	        2
    # Sx → Tx × 4 → FEC × 4	      – 	      1

    # In the validation cohort
    # 1) the 227 samples with T+FAC comes from MDACC and
    # 2) the 83 samples with Taxane+AC comes from ISPY
    #
    # In the validation cohort
    # 1) only the docetaxel arm is selected (comes mainly from MDACC ans USO).
    # Docetaxel arm is purely neoadjuvant(!!! ambiguity in 3 cases,
    # where it can be potentialy purely adjuvant. howevre 3 ambigous samples will
    # not affect results, and can be considered for analysis)
    # 2) !!!! discarded the paclitaxel arms, as it potentially contain regimens spanning
    # neoadjuvat/adjuvant/both regimens.

    # Note!!! Hormone therapy is given only in adjuvant regimen

    Regimen = "neoadj_adj",
    Arm_neoadj = purrr::map2_chr(
      Type_taxane, Source,
      ~(case_when(
        (is.na(.x) & .y == "MDACC") ~ "Paclitaxel+FAC",
        (is.na(.x) & .y == "ISPY") ~ "Taxane+AC",
        (.x == "Taxotere") ~ "Docetaxel+Capecitabin+FEC"
      ))
    ),

    Arm_adj = if_else(Er == "pos", "Hormone", "No_hormone"),

    Arm_detailed = "In the neoadj regimen \"Docetaxel+Capecitabin+FEC\", three samples may come from entirely adjuvant setting."
  )


table(xx1$Arm_neoadj)
# Docetaxel+Capecitabin+FEC            Paclitaxel+FAC                 Taxane+AC
# 92                       227                        83
table(xx1$Arm_adj)
# Hormone No_hormone
# 297        205

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE25066_cleaned.tsv"
)




# GSE20685
# >>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE20685.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx)

table(xx$Event_death)
# 0   1
# 244  83
table(xx$Event_metastasis)
# 0   1
# 244  83
table(xx$Regional_relapse)
# 0   1
# 282  25

xx %>%
  dplyr::group_by(Event_metastasis, Regional_relapse) %>%
  dplyr::summarise(N = n())
#   Event_metastasis Regional_relapse     N
# 1                0                0   224
# 2                0                1     8
# 3                0               NA    12
# 4                1                0    58
# 5                1                1    17
# 6                1               NA     8

table(is.na(xx$`Follow_up_duration_(years)`))
# FALSE
# 327
table(is.na(xx$`Time_to_metastasis_(years)`))
# FALSE  TRUE
# 83   244
table(is.na(xx$`Time_to_relapse_(years)`))
# FALSE  TRUE
# 25   302

table(xx$T_stage)
# 1  1c   2   3   4
# 99   2 188  26  12
table(xx$N_stage)
# 0   1   2   3
# 137  87  63  40

table(xx$Adjuvant_chemotherapy)
# no unknown     yes
# 54       5     268
table(xx$`Regimen_(caf_vs_cmf)`)
# CAF     CMF      no  others unknown
# 88      61      54     119       5
table(xx$Neoadjuvant_chemotherapy)
# yes
# 30


xx %>%
  dplyr::group_by(Adjuvant_chemotherapy, Neoadjuvant_chemotherapy, `Regimen_(caf_vs_cmf)`) %>%
  dplyr::summarise(N = n())
#   Adjuvant_chemotherapy Neoadjuvant_chemotherapy `Regimen_(caf_vs_cmf)`     N
# 1 no                    NA                       no                        54 !+- Hormone
# 2 unknown               yes                      unknown                    1
# 3 unknown               NA                       unknown                    4
# 4 yes                   yes                      CAF                        1 !+- Hormone
# 5 yes                   yes                      others                    28 !Taxane +- Hormone
# 6 yes                   NA                       CAF                       87 !+- Hormone
# 7 yes                   NA                       CMF                       61 !+- Hormone
# 8 yes                   NA                       others                    91 !Taxane +- Hormone

# The regimen others can be mapped to taxane containing regiement based on the
# following text.
# "The selected tissue samples spanned the major transition periods of adjuvant chemotherapy from CMF (cyclophosphamide, methotrexate and fluorouracil) to CAF (cyclophosphamide, doxorubicin, fluorouracil) and to taxane-based regimens."
# "All patients were treated by a multidisciplinary team according to the guidelines consistent with the National Comprehensive Cancer Network [18]. Following modified radical mastectomy or breast-conserving surgery plus dissection of axillary nodes, patients received radiotherapy, adjuvant chemotherapy, and/or hormonal therapy, if indicated. Neoadjuvant chemotherapy was administered to patients with locally advanced disease."

# Filter samples treated with neoadj chemotherapy


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age_at_diagnosis,
    Node_cat = str_c("N",N_stage),
    Node_bin = if_else(N_stage < 1, "neg", "pos"),
    Size_cat = str_c("T",T_stage),
    Size_bin = if_else( str_replace(T_stage, "c", "") %>% as.integer() <= 1,
                        "samall", "large"), # small = T1/T1c large = >T1
    Event_dfs = Event_metastasis %>% as.integer(),
    Time_dfs = purrr::map2_dbl(
      `Time_to_metastasis_(years)`,`Follow_up_duration_(years)`,
      ~(if_else(is.na(.x), .y, .x))
    ),

    Regimen = "adj",
    Arm_adj = purrr::map_chr(
      `Regimen_(caf_vs_cmf)`,
      ~(case_when(.x == "CAF" ~ "FAC[+Hormone]",
                  .x == "CMF" ~ "CMF[+Hormone]",
                  # .x == "others" ~ "Taxane_based_regimen[+Hormone]",
                  # Others amy include hormone-only therapy too, hence set as NA
                  .x == "no" ~ "No_systemic_therapy[+Hormone]"))
    ),
    # Setting Arm_adj = NA for neoadjuvant treated samples, to exclude them from analysis
    Arm_adj = if_else(is.na(Neoadjuvant_chemotherapy), Arm_adj, NA_character_)
  )

table(xx1$Arm_adj)
# CMF[+Hormone]                 FAC[+Hormone] No_systemic_therapy[+Hormone]
# 61                            87                            54

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE20685_cleaned.tsv"
)




# GSE41998
# >>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE41998.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx)

table(xx$Gender)
# female
# 279
table(xx$Treatment_arm)
# Ixabepilone        none  Paclitaxel
# 138          14         127
# "Patients received sequential neoadjuvant therapy starting with 4 cycles of AC (doxorubicin 60 mg/m2 intravenously and cyclophosphamide 600 mg/m2 intravenously) given every 3 weeks, followed by 1:1 randomization to either ixabepilone (40 mg/m2 3-hour infusion) every 3 weeks for 4 cycles, or paclitaxel (80 mg/m2 1-hour infusion) weekly for 12 weeks."
table(xx$Ac_response)
# 0   complete response    partial response progressive disease      stable disease
# 3                  40                 161                   5                  64
# unable to determine
# 6
table(xx$Basetms)
# < 2 cm       > 5 cm     2 - 5 cm not reported
# 3           99          174            3
table(xx$Er)
# negative positive
# 171      108
table(xx$Pr)
# negative positive  UNKNOWN
# 179       99        1
table(xx$Her2stat)
# negative positive
# 251       28
table(xx$Her2)
# other positive
# 251       28
table(xx$Mnssl)
# not reported peri-menopausal post-menopausal  pre-menopausal
# 6              13             124             136
table(xx$Pcr)
# 0  No Yes
# 20 184  69
# "Surgical specimens were evaluated by a staff pathologist at each study site; no central pathology review was conducted. The pCR rate was evaluated as the primary endpoint, with pCR defined by no histologic evidence of residual invasive adenocarcinoma in the breast and axillary lymph nodes, with or without the presence of ductal carcinoma in situ."
table(xx$Pcrrcb1)
# 0  No Yes
# 20 167  86


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Size_cat = purrr::map_chr(
      Basetms,
      ~(case_when(.x == "< 2 cm" ~ "T1",
                  .x == "2 - 5 cm" ~ "T2",
                  .x == "> 5 cm" ~ "T3"))
    ),
    Size_bin = purrr::map_chr(
      Basetms,
      ~(case_when(.x == "< 2 cm" ~ "small",
                  .x == "2 - 5 cm" ~ "large",
                  .x == "> 5 cm" ~ "large"))
    ),
    Er = if_else(Er == "negative", "neg", "pos"),
    Pr = if_else(Pr == "UNKNOWN", NA_character_, Pr) %>%
      purrr::map_chr(~(if_else(.x == "positive", "pos", "neg"))),
    Her2 = if_else(Her2stat == "negative", "neg", "pos"),
    Menopause =  purrr::map_chr(
      Mnssl,
      ~(case_when(.x == "peri-menopausal" ~ "peri",
                  .x == "post-menopausal" ~ "post",
                  .x == "per-menopausal" ~ "pre"))
    ),
    Response_pathological = purrr::map_chr(
      Pcr,
      ~(case_when(.x == "No" ~ "npCR",
                  .x == "Yes" ~ "pCR"))
    ),
    Response =  Response_pathological,
    Regimen = "neoadj",
    Arm_neoadj = purrr::map_chr(
      Treatment_arm,
      ~(case_when(.x == "Ixabepilone" ~ "Ixabepilone+AC",
                  .x == "Paclitaxel" ~ "Paclitaxel+AC",
                  .x == "none" ~ "Untreated"))
    )
  )

table(xx1$Arm_neoadj)
# Ixabepilone+AC  Paclitaxel+AC      Untreated
# 138            127             14

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE41998_cleaned.tsv"
)




# GSE20194
# >>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE20194.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx)

table(xx$Tbefore)
# 0   1   2   3   4
# 3  23 147  50  53
table(xx$Nbefore)
# 0   1   2   3
# 79 125  31  42
table(xx$Bmngrd)
# 1   2   3
# 13 104 150
table(xx$Er_status)
# N   P
# 114 164
table(xx$Pr_status)
# N   P
# 157 121
table(xx$Her2_status)
# N   P
# 219  59
table(xx$Histology)
# IC/DCIS                     IDC              IDC + DCIS IDC, medullary features
# 1                     216                       4                       1
# IDC/anaplastic                IDC/DCIS                 IDC/IBC                 IDC/ILC
# 1                      20                       2                       7
# IDC/ILC/DCIS            IDC/ILC/LCIS                 IDC/LDC                IDC+DCIS
# 1                       1                       1                       3
# IDC+ILC                     ILC                 ILC/IDC                     IMC
# 3                      10                       4                       1
# IMC/IDC
# 1
table(xx$Race)
# asian    black hispanic    mixed    white
# 18       29       42        3      176
table(xx$Pcr_vs_rd)
# pCR  RD
# 56 222
table(xx$Treatment_code)
# FAC       FACT FACT+XRT/X        FEC       FECT       TFAC    TFAC/HT       TFEC     TH/FAC
# 3          1          1          1          2        210          1         33          6
# TH/FEC      Tonly      TXFAC
# 2          1          9


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = str_c("G", Bmngrd),
    Node_cat = str_c("N", Nbefore),
    Node_bin = if_else(Nbefore == 0, "neg", "pos"),
    Size_cat = str_c("T", Tbefore),
    Size_bin = if_else(Tbefore <= 1, "small", "large"),
    Er = if_else(Er_status == "N", "neg", "pos"),
    Pr = if_else(Pr_status == "N", "neg", "pos"),
    Her2 = if_else(Her2_status == "N", "neg", "pos"),
    Ethnicity = Race,
    Histology = Histology,
    Response_pathological = if_else(Pcr_vs_rd == "RD", "npCR", "pCR"),
    Response =  Response_pathological,

    Regimen = "neoadj",
    Arm_neoadj = Treatment_code
    # Arm:
    # There exists missmatch between treatemnt_code and tratment_comments
    # Do manual formatting
  )

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2880423/
# Patients received 6 months of preoperative (neoadjuvant) chemotherapy including paclitaxel, 5-fluorouracil, cyclophosphamide, and doxorubicin, followed by surgical resection of the cancer. Response to preoperative chemotherapy was categorized as a pathologic complete response (pCR = no residual invasive cancer in the breast or lymph nodes) or residual invasive cancer (RD)


# old version of formatted
xx2 <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE20194_cleaned_old.tsv")
identical(xx1$Sample_geo_accession, xx2$Sample_geo_accession) # TRUE

table(xx1$Arm_neoadj)
table(xx2$Arm)


xx2$Arm <- xx2$Arm %>%
  str_replace("anastrozole",  "Anastrozole") %>%
  str_replace("Xelonda",  "Capecitabine") %>%
  str_replace("Gleevec",  "Imatinib") %>%
  str_replace("_H",  "_Trastuzumab") %>%
  str_replace("_Tam",  "_Tamoxifen") %>%
  str_replace("TF",  "Paclitaxel_F") %>%
  purrr::map_chr(~(if_else(.x == "T", "Paclitaxel", .x)))


identical(xx1$Sample_geo_accession, xx2$Sample_geo_accession) # TRUE
xx1$Arm_neoadj <- xx2$Arm

table(xx1$Arm_neoadj)
# FAC                                            FEC
# 3                                              1
# Paclitaxel                                 Paclitaxel_FAC
# 5                                            181
# Paclitaxel_FAC_Anastrozole                    Paclitaxel_FAC_Capecitabine
# 1                                              1
# Paclitaxel_FAC_Capecitabine_Letrozole_Imatinib                       Paclitaxel_FAC_Tamoxifen
# 1                                              1
# Paclitaxel_FAC_Trastuzumab                                 Paclitaxel_FAX
# 4                                              2
# Paclitaxel_FEC                     Paclitaxel_FEC_Trastuzumab
# 73                                              4


write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE20194_cleaned.tsv"
)




# GSE32603
# >>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE32603.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx)

table(xx$Tissue)
# Breast
# 248
table(xx$Treatment)
# F1  F2  F3  F4  FS
# 130  11  40   5  62
table(xx$Biopsy_time)
# T1  T2  TS
# 141  45  62

xx %>%
  dplyr::group_by(Biopsy_time, Treatment) %>%
  dplyr::summarise(N = n())
#   Biopsy_time Treatment     N
# 1 T1          F1          130 !!! AC+T, imputed from Table 1 of PMID:22198468
# 2 T1          F2           11 !!! discard, ambigious arm
# 3 T2          F3           40 !!! discard, only posttreatment samples
# 4 T2          F4            5 !!! discard, only postttreatment sample
# 5 TS          FS           62 !!! discard, only posttreatment sample

table(xx$Recurrence_yes_is_1)
# 0   1
# 169  79
table(xx$`Hr-positive_yes_is_1`)
# 0   1
# 91 154
table(xx$`Her2-positive_yes_is_1`)
# 0   1
# 186  55
table(xx$Pcr_yes_is_1)
# 0   1
# 193  51
table(xx$Rcb_class)
# 0  1  2  3
# 50 14 99 65


xx1 <- xx %>%
  dplyr::mutate(
    Size_bin = rep("large", nrow(xx)), # minimum tumor size >= 3cm
    Hr = if_else(`Hr-positive_yes_is_1` == 0, "neg", "pos"),
    Her2 = if_else(`Her2-positive_yes_is_1` == 0, "neg", "pos"),
    Response_pathological = if_else(Pcr_yes_is_1 == 0, "npCR", "pCR"),
    Response =  Response_pathological,
    Event_dfs = Recurrence_yes_is_1 %>% as.integer(),
    Time_dfs =  Rfs_days/365, # years

    Regimen = "neoadj", # adj treatments are ambigous hense discarded

    Arm_neoadj = if_else(Treatment == "F1", "Taxane+AC", NA_character_),
    # Largest arm in the trial is AC + T, all other arms are discarded
    Timepoint = purrr::map_chr(
      Biopsy_time,
      ~(case_when(.x == "T1" ~ "Pre_treatment",
                  .x == "T2" ~ "On_treatment",
                  .x == "TS" ~ "Post_treatment"))
    )
  )


xx1 %>%
  dplyr::group_by(Timepoint, Arm_neoadj) %>%
  dplyr::summarise(N = n())
#   Timepoint      Arm           N
# 1 On_treatment   NA           45
# 2 Post_treatment NA           62
# 3 Pre_treatment  Taxane+AC   130
# 4 Pre_treatment  NA           11


table(xx1$Arm_neoadj)
# Taxane+AC
# 130

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE32603_cleaned.tsv"
)




# GSE22219
# >>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE22219.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx)

table(xx$Tissue)
# Breast Cancer
# 216
table(xx$Er_status)
# 0   1
# 82 134
table(xx$Tumour_grade)
# .  1  2  3
# 25 41 87 63
table(xx$`Distant-relapse_event`)
# 0   1
# 134  82



xx1 <- xx %>%
  dplyr::mutate(
    Age = Patient_age,

    # T0 No evidence of primary tumor
    # T1 size < 2cm
    # T2 2cm < size < 5cm
    # T3 size > 5cm
    Size = Tumour_size,
    Size_bin = if_else(Tumour_size < 2, "small", "large"),
    Size_cat = purrr::map_chr(
      Tumour_size,
      ~(case_when(.x == 0 ~ "T0",
                  .x <= 2 ~ "T1",
                  .x <= 5 ~ "T2",
                  .x > 5 ~ "T3"))
    ),

    # pNX: Regional lymph nodes cannot be assessed (e.g., previously removed, or not removed for pathologic study)
    # pN0: No regional lymph node metastasis identified histologically
    # pN1: Micrometastases; or metastases in 1–3 axillary lymph nodes; and/or in internal mam mary nodes with metastases detected by sentinel lymph node biopsy but not clinically detected
    # pN2: Metastases in 4–9 axillary lymph nodes; or in clinically detected. internal mammary lymph nodes in the absence of axillary lymph node metastases
    # pN3: Metastases in ten or more axillary lymphnodes; or in infraclavicular (level III axillary) lymph nodes; or in clinically detected.
    Node = Nodes_involved,
    Node_bin = if_else( Nodes_involved == 0, "neg", "pos"),
    Node_cat = purrr::map_chr(
      Nodes_involved,
      ~(case_when(.x == 0 ~ "N0",
                  .x <= 3 ~ "N1",
                  .x <= 9 ~ "N2",
                  .x > 9 ~ "N3"))
    ),
    Grade = if_else(Tumour_grade == ".", NA_character_, str_c("G", Tumour_grade)),
    Er = if_else(Er_status == 0, "neg", "pos"),
    Event_dfs = `Distant-relapse_event` %>% as.integer(),
    Time_dfs = `Distant-relapse_free_survival`, # years
    # Arm imputed from supplementary methods of PMID: 21737487
    # "Patients received surgery followed by adjuvant chemotherapy, adjuvant hormone therapy, both of these therapies, or no adjuvant treatment. Tamoxifen was used as endocrine therapy for 5 years in all ER+ BCs. Patients who were <50 years of age, with lymph node positive tumors, or ER– and/or >3 cm in diameter, received adjuvant cyclophosphamide, methotrexate, and 5-fluorouracil (CMF). Patients >50 years of age with ER–, lymph node–positive tumors also received CMF for six cycles, at a thrice weekly intravenous regimen"

    Regimen = "adj",

    Arm_adj = if_else(
      Age < 50,
      if_else(Node > 0 | (Er == "neg" & Size > 3), "CMF", "noCMF"),
      if_else(Node > 0 & Er == "neg", "CMF", "noCMF")
    ),
    Arm_adj = if_else(Er == "pos", str_c(Arm_adj,"_Tamoxifen"), Arm_adj),
    Arm_adj = if_else(Arm_adj == "noCMF", "No_systemic_therapy", Arm_adj),
    Arm_adj = if_else(Arm_adj == "noCMF_Tamoxifen", "Tamoxifen", Arm_adj)
  )


table(xx1$Arm_adj)
# CMF       CMF_Tamoxifen No_systemic_therapy           Tamoxifen
# 38                  14                  44                 120
#
# CMF = 38+14 = 52 (There is mismatch in Er pos in text and geo-data, further 3 samples are missing )
# Tamoxifen = 14+120 = 134 (There is mismatch in Er pos in text and geo-data)




write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE22219_cleaned.tsv"
)




# GSE4056
# >>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE4056.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx)

xx[,"Grading"] <- str_c("Grading: ", xx[,"Grading"] %>% tibble::deframe())

yy <- str_split_fixed(xx[,"Grading"] %>% deframe(), ",", 8) %>% as_tibble()
names(yy) <- get_var_name(yy) %>% str_trim()
xx <- bind_cols(xx[, 1:9], yy)
names(xx) <- clean_names(names(xx))
glimpse(xx) # change Ypn to ypN # post therapy node status

names(xx) <- str_replace(names(xx), "Ypn", "ypN")



glimpse(xx %>% dplyr::filter(str_detect(Sample_title, "Cy3")))
# All characteristics values are NAs
glimpse(xx %>% dplyr::filter(str_detect(Sample_title, "Cy5")))
# All characteristics values are present

# Discard Cy3 samples, as they are dye swaping experiment of Cy5
# Take only Cy5

xx1 <- xx %>% dplyr::filter(str_detect(Sample_title, "Cy5"))
xx2 <- xx %>% dplyr::filter(str_detect(Sample_title, "Cy3"))

glimpse(xx1)
glimpse(xx2)

table(xx1$Grading)
# Grading: G1      Grading: G2      Grading: G3 Grading: no data
# 4               40               46                4
table(xx1$Her2)
# HER2: 0  HER2: 1  HER2: 2  HER2: 3
# 68        4        2       20
table(xx1$Er)
# ER: 0        ER: 1       ER: 12        ER: 2        ER: 3        ER: 4        ER: 5        ER: 6
# 29            1           20            5            1            7            1           14
# ER: 8        ER: 9  ER: no data
# 10            5            1
table(xx1$Pr)
# PR: 0   PR: 1  PR: 12   PR: 2   PR: 4   PR: 6   PR: 8   PR: 9
# 34       1      23       3       9      11       7       6
table(xx1$Therapy_status_pathological)
# Therapy_Status_Pathological: pCR  Therapy_Status_Pathological: pNC  Therapy_Status_Pathological: pPR
# 24                                 9                                61
table(xx1$ypN)
# ypN: N+  ypN: N0  ypN: ND
# 28       65        1


xx3 <- xx1 %>%
  dplyr::rename(
    Grade = "Grading",
    Her2_score = "Her2",
    Er_score = "Er",
    Pr_score = "Pr",
    Response_pathological = "Therapy_status_pathological"
  ) %>%
  dplyr::mutate(
    Grade = str_replace(Grade, "Grading: ", ""),
    Grade = str_replace(Grade, "no data", NA_character_),
    Her2_score = str_replace(Her2_score, " HER2: ", "") %>% as.numeric(),
    Her2 = if_else(Her2_score > 0, "pos", "neg"),
    Er_score = str_replace(Er_score, " ER: ", "") %>% as.numeric(),
    Er = if_else(Er_score > 0, "pos", "neg"),
    Pr_score = str_replace(Pr_score, " PR: ", "") %>% as.numeric(),
    Pr = if_else(Pr_score > 0, "pos", "neg"),
    Response_pathological = str_replace(Response_pathological %>% str_trim(),
                                        "Therapy_Status_Pathological: ", ""),
    Response = if_else(Response_pathological == "pCR", "pCR", "npCR"), # integrated
    ypN = str_trim(ypN) %>% str_replace("ypN: ", "") %>% str_replace("ND", NA_character_),
    Ct_max = str_replace(Ct_max %>% str_trim(), "cT_max: ", "") %>% as.numeric(),
    Size = Ct_max, # in cm
    Size_bin = if_else(Size < 2, "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when(.x == 0 ~ "T0",
                  .x <= 2 ~ "T1",
                  .x <= 5 ~ "T2",
                  .x > 5 ~ "T3"))
    ),
    Bcl2_score = str_replace(Bcl2_score %>% str_trim(), "bcl2_Score: ", "") %>%
      as.numeric(),

    Regimen = "neoadj",
    Arm_neoadj = "Docetaxel+Gemcitabine+Epirubicin"
    # "Within two sequential phase I and II studies, patients received one of the following two regimens: (1) five cycles of E (90 to 100 mg/m2) in combination with G (1,250 mg/m2) repeated every 2 weeks, sequentially followed by four cycles Doc (80 to 100 mg/m2) every 2 weeks (GEsDoc) with prophylactic pegfilgrastim; or (2) six cycles of G (800 mg/m2), E (60 to 90 mg/m2), and Doc (60 to 75 mg/m2) every 3 weeks (GEDoc) with prophylactic filgrastim. G (800 mg/m2) was repeated on day 8 of each cycle."
  )


table(xx3$Arm_neoadj)
# Docetaxel+Gemcitabine+Epirubicin
# 94

write_tsv(
  x = xx3,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE4056_cleaned.tsv"
)




# GSE34138
# >>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE34138.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx)


table(xx$Response_breast)
# pCR+npCR mamma          PR+NR
# 58            119
table(xx$Response)
# CR CR+kl    NE  NOCR    PD
# 28    10     1   137     2
table(xx$Treatment)
# AC
# 178
table(xx$Clinical_subtype)
# LUM  TN
# 120  58
table(xx$Tissue)
# breast cancer tumor
# 178


xx1 <- xx %>%
  dplyr::mutate(
    Er_score = Ihc_er_percentage,
    Er = if_else(Er_score <= 10, "neg", "pos"),
    Pr_score = Ihc_pr_percentage,
    Pr = if_else(Pr_score <= 10, "neg", "pos"),
    Her2 = "neg", # imputed from publication and trial design
    Response_clinical = Response,
    Response = if_else(Response_clinical == "CR", "pCR", "npCR"),

    Regimen = "neoadj",
    Arm_neoadj = Treatment
  )


xx1 %>%
  dplyr::group_by(str_c(Er,Her2),Clinical_subtype) %>%
  dplyr::summarise(N = n())
#   `str_c(Er, Her2)` Clinical_subtype     N
# 1 negneg            LUM                  1
# 2 negneg            TN                  57
# 3 posneg            LUM                118
# 4 NA                LUM                  1
# 5 NA                TN                   1

table(xx1$Arm_neoadj)
# AC
# 178


write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE34138_cleaned.tsv"
)



# GSE20271
# >>>>>>>>>>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE20271.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx)


table(xx$Pcr_or_rd)
# pCR  RD
# 26 152
table(xx$Race)
# A  B  H  W
# 1 13 83 81
table(xx$Histology)
# IDC                             IDC/IBC
# 127                                   2
# IDC/ILC                                 ILC
# 3                                   4
# ILC/IDC       infiltrating Ductal Carcinoma
# 1                                   1
# Infiltrating Ductal Carcinoma      Infiltrating Lobular Carcinoma
# 33                                   1
# infiltrating Mucosecretor Carcinoma      Inflitrating Lobular Carcinoma
# 1                                   2
# Mixed; Ductal and Lobular                   Only Infiltrating
# 2                                   1
table(xx$Prechemo_t)
# 0   1   2   3   4 N/A
# 2  11  76  37  51   1
table(xx$Prechemo_n)
# 0   1   2   3 N/A
# 59  71  38   9   1
table(xx$Bmn_grade)
# 1   2   3 N/A
# 15  61  72  30
table(xx$Er_status)
# N  P
# 80 98
table(xx$Pr_status)
# N  P
# 95 83
table(xx$Her_2_status)
# N   P
# 152  26
length(table(xx$Preoperative_treatment))
# 62 arms; needs manual cleaning !!!!!!!!!!!!!!! use previous version !!!!
table(xx$`Treatment_received_(1=fac,_2=t/fac)`)
# 1  2
# 87 91
table(xx$`Randomized_(1=fac,_2=t/fac)`)
# 1  2
# 94 84


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Ethnicity = Race,
    Size_cat = str_c("T",as.integer(Prechemo_t)),
    Size_bin = if_else(as.integer(Prechemo_t) <= 1, "small", "large"),
    Node_cat = str_c("N",as.integer(Prechemo_n)),
    Node_bin = if_else(as.integer(Prechemo_n) == 0, "neg", "pos"),
    Grade =  str_c("G",as.integer(Bmn_grade)),
    Er_score = `Er%_positive`,
    Er = if_else(Er_status == "N", "neg", "pos"),
    Pr_score = `Pr%_positive`,
    Pr = if_else(Pr_status == "N", "neg", "pos"),
    Her2_score.ihc = Her_2_ihc,
    Her2_score.fish = Her_2_fish,
    Her2 = if_else(Her_2_status == "N", "neg", "pos"),
    Response_pathological = Pcr_or_rd,
    Response = if_else(Pcr_or_rd == "pCR", "pCR", "npCR"),

    Regimen = "neoadj",
    Arm_neoadj = rep(NA, nrow(xx))
)

xx2 <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE20271_cleaned_old.tsv")
names(xx2) <- clean_names(names(xx2))
glimpse(xx2) # non standard character
glimpse(xx2 %>% dplyr::select(-c("Sample_procurement_details", "Archive_details")))
table(xx2$Arm) # n = 178
# FAC  FEC TFAC TFEC
# 49   38   46   45

identical(xx2$Sample_geo_accession, xx1$Sample_geo_accession) # TRUE
xx1$Arm_neoadj <- xx2$Arm

table(xx1$Arm_neoadj)
# FAC  FEC TFAC TFEC
# 49   38   46   45

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE20271_cleaned.tsv"
)



# GSE109710
# >>>>>>>>>>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE109710.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx)


table(xx$Gender)
# Female
# 173
table(xx$Tissue)
# breast cancer tumor
# 173
table(xx$Tumor_type)
# HER2 Breast Cancer
# 173
table(xx$Estrogen_receptor_status)
# NO YES
# 94  79
table(xx$Grade)
# 0  1  2  3
# 33  6 68 66
table(xx$Progesterone_receptor_status)
# NO YES
# 110  63
table(xx$Pcr)
# 0   1
# 70 103
table(xx$Dfs.event)
# .   0   1
# 11 136  26

xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade =  str_c("G",as.integer(Grade)),
    Er = if_else(Estrogen_receptor_status == "NO", "neg", "pos"),
    Pr = if_else(Progesterone_receptor_status == "NO", "neg", "pos"),
    Her2 = "pos",
    Response_pathological = Pcr,
    Response = if_else(Pcr == 1, "pCR", "npCR"),
    Event_dfs = Dfs.event %>% as.integer(),
    Time_dfs = xx$Dfs.time %>% as.numeric() / 365,

    Regimen = "neoadj_adj",
    Arm_neoadj = "Docetaxel+Trastuzumab+Pertuzumab+[FEC|Carboplatin]",
    Arm_adj = "[Trastuzumab+Chemo+Hormone+Radio]"
  )

table(xx1$Arm_neoadj)
# Docetaxel+Trastuzumab+Pertuzumab+[FEC|Carboplatin]
# 173
table(xx1$Arm_adj)
# [Trastuzumab+Chemo+Hormone+Radio]
# 173

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE109710_cleaned.tsv"
)




# GSE6861
# >>>>>>>>>>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE6861.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 161

yy <- read_csv(file = "data/geo_data3/gse_supplementary/GSE6861_demo_from_paper_mmc5.csv")
names(yy) <- clean_names(names(yy))
glimpse(yy) # 125

xx <- xx %>%
  dplyr::mutate(
    Patient_id = str_replace(Sample_title, "HB", "") %>%
      str_replace("bis", "") %>% as.numeric()
  ) %>%
  left_join(yy , by = "Patient_id")
glimpse(xx) # 161

table(xx$Treatment[str_detect(xx$Sample_title, "bis")])
# FEC
# 66
# Sample_title with pattern "bis" is FEC arm !!!!!!!!!!!!

table(xx$Pcr)
# npCR  pCR
# 95   66
table(xx$Response)
# No pCR    pCR
# 70     55
table(xx$Treatment)
# FEC TET
# 66  59
table(xx$Size)
# T1 T2 T3
# 3 76 46
table(xx$Node)
# N0 N1 N2
# 49 66 10
table(xx$Pr)
# neg pos
# 122   3
table(xx$Grade)
# 1  2  3
# 2 37 70
table(xx$Histology)
# Ductal lobular Lobular   Other
# 113       2       4       6


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade =  str_c("G",as.integer(Grade)),
    Size_cat = Size,
    Size_bin = if_else(Size == "T1", "small", "large"),
    Node_cat = Node,
    Node_bin = if_else(Node == "N0", "neg", "pos"),
    Er = "neg",
    Pr = Pr,
    Response_pathological = Pcr,
    Response = Pcr,

    Regimen = "neoadj",
    Arm_neoadj = if_else(str_detect(Sample_title, "bis"), "FEC", "Docetaxel+Epirubicin")
  )

table(xx1$Arm_neoadj)
# Docetaxel+Epirubicin                  FEC
#                   59                  102
write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE6861_cleaned.tsv"
)



# GSE22358
# >>>>>>>>>>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE22358.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 161

table(xx$Tissue)
# breast cancer normal breast !!!!!!!!!!!! discard normal breast(set Arm = NA)
# 154             4
table(xx$Study)
# Xena
# 154
table(xx$P53_status_by_amplichip)
# ---   MUTANT WILDTYPE
# 8       73       73
table(xx$P53_measured_by_ihc)
# ---     NEGATIVE     POSITIVE UNDETERMINED
# 1           53           79           21
table(xx$Grade)
# ---   1   2   3
# 15  19  50  70
table(xx$`Er_(0_=_negative;_1_=_positive)`)
# ---   0   1
# 1  71  82
table(xx$`Pgr_(0_=_negative;_1_=_positive)`)
# ---   0   1
# 1  98  55
table(xx$`Her2_(0_=_negative;_1_=_positive)`)
# 0   1
# 120  34
table(xx$Response)
# ---      COMPLETE RESPONSE NEAR COMPLETE RESPONSE            NO RESPONSE
# 32                     20                     10                      6
# PARTIAL RESPONSE
# 86
table(xx$Neoadjuvant_chemotherapy)
# (docetaxel-Capecitabine)+ Trastuzumab                docetaxel-Capecitabine
# 34                                   120

xx %>%
  dplyr::group_by(`Her2_(0_=_negative;_1_=_positive)`,Neoadjuvant_chemotherapy, Tissue) %>%
  dplyr::summarise(N = n())
#   `Her2_(0_=_negative;_1_=_positive)` Neoadjuvant_chemotherapy              Tissue            N
# 1                                   0 docetaxel-Capecitabine                breast cancer   120
# 2                                   1 (docetaxel-Capecitabine)+ Trastuzumab breast cancer    34
# 3                                  NA NA                                    normal breast     4

xx1 <- xx %>%
  dplyr::mutate(
    Mut_p53 = purrr::map_chr(
      P53_status_by_amplichip,
      ~(case_when(.x == "MUTANT" ~ "Mut",
                  .x == "WILDTYPE" ~ "WT"))
    ),
    Ihc_p53 = purrr::map_chr(
      P53_measured_by_ihc,
      ~(case_when(.x == "NEGATIVE" ~ "neg",
                  .x == "POSITIVE" ~ "pos"))
    ),
    Size = `Tumor_size_(cm)`,
    Size_bin = if_else(Size <=2, "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when(.x == 0 ~ "T0",
                  .x <=2 ~ "T1",
                  .x <=5 ~ "T2",
                  .x >5 ~ "T3"))
    ),
    Grade =  str_c("G",as.integer(Grade)),
    Er = if_else( `Er_(0_=_negative;_1_=_positive)` %>% as.integer() == 0, "neg", "pos"),
    Pr = if_else(`Pgr_(0_=_negative;_1_=_positive)` %>% as.integer() == 0, "neg", "pos"),
    Her2 = if_else(`Her2_(0_=_negative;_1_=_positive)` %>% as.integer() == 0, "neg", "pos"),
    Response_pathological = if_else(Response == "---", NA_character_, Response),
    Response = if_else(
      (Response_pathological == "COMPLETE RESPONSE" |
        Response_pathological == "NEAR COMPLETE RESPONSE"),
      "pCR", "npCR"),

    Regimen = "neoadj",
    Arm_neoadj = if_else(Her2 == "neg",
                  "Docetaxel+Capecitabine", "Docetaxel+Capecitabine+Trastuzumab")
  )

table(xx1$Arm_neoadj)
# Docetaxel+Capecitabine Docetaxel+Capecitabine+Trastuzumab
# 120                                 34
write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE22358_cleaned.tsv"
)


# GSE50948
# >>>>>>>>>>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE50948.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 156


table(xx$Invasive_tumor_grade)
# 2  3
# 68 86
table(xx$Race)
# Black Caucasian  Oriental
# 1       154         1
table(xx$Er)
# ER- ER+
#   104  52
table(xx$Pr)
# PR- PR+
#   121  35
table(xx$Her2)
# HER2- HER2+
#   42   114
table(xx$Arm)
# HER2- CT   HER2+ CT HER2+ CT H
# 42         51         63
table(xx$Treatment)
# neoadjuvant doxorubicin/paclitaxel (AT) followed by cyclophosphamide/methotrexate/fluorouracil (CMF)
# 93
# neoadjuvant doxorubicin/paclitaxel (AT) followed by cyclophosphamide/methotrexate/fluorouracil (CMF) + Trastuzumab for 1 year
# 63
table(xx$Menopausal.status)
# post  pre
# 84   72
table(xx$Pcr)
# pCR  RD
# 53 103


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Til = `Inflammatory_cells_[%]`,
    Grade = str_c("G",Invasive_tumor_grade),
    Ethnicity = Race,
    Er = if_else(Er == "ER-", "neg", "pos"),
    Pr = if_else(Pr == "PR-", "neg", "pos"),
    Her2 = if_else(Her2 == "HER2-", "neg", "pos"),
    Menopause = Menopausal.status,
    Response_pathological = Pcr,
    Response = if_else(Pcr == "RD", "npCR", "pCR"),

    Regimen = "neoadj",
    Arm_neoadj = if_else(
      Treatment == "neoadjuvant doxorubicin/paclitaxel (AT) followed by cyclophosphamide/methotrexate/fluorouracil (CMF)",
      "Paclitaxel+Doxorubicin+CMF",
      "Paclitaxel+Doxorubicin+Trastuzumab+CMF"
    )
  )

table(xx1$Arm_neoadj)
# Paclitaxel+Doxorubicin+CMF Paclitaxel+Doxorubicin+Trastuzumab+CMF
# 93                                     63

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE50948_cleaned.tsv"
)



# GSE22226_GPL4133
# >>>>>>>>>>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE22226_GPL4133.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 20

table(xx$`Histologic_grade_(1=grade_i_(low);__2=_grade_ii_(intermediate);_3=_grade_iii_(high);_4=indeterminate)`)
# 1  2  3
# 3 11  6
table(xx$`Histology_(1=necrosis;_2=ductal_carcinoma;_3=lobular;_4=mixed_ductal/lobular_carcinoma;_5=other_6=no_invasive_tumor_present;)`)
# 2  3  4
# 16  3  1
table(xx$`Clinical_t_stage_(pre-chemotherapy)_(1_=_<=_2cm;_2_=_>2-5_cm;_3_=_>_5cm)`)
# 1  2  3
# 1 12  7
table(xx$`Er_(0=negative;_1=positive;)`)
# 0  1
# 5 14
table(xx$`Pgr__(0=negative;_1=positive;)`)
# 0  1
# 6 13
table(xx$`Her2_(0=negative;_1=positive;)`)
# 0  1
# 15  4
table(xx$ `Neoadjuvant_chemotherapy*_(1_=_ac_only;_2_=_ac_+_t;_3_=_ac_+_t_+_herceptin;_4_=_ac_+_t_+_other;_[a_=_doxorubicin,_c_=_cyclophosphamide,_t_=_taxane])` )
# 2  3
# 17  3
table(xx$`Pathological_complete_response_(pcr)` )
# No Yes
# 16   4
table(xx$Rcb_class)
# 0=RCB index 0          I=RCB index less than or equal to 1.36
# 4                                               3
# II=RCB index greater than 1.36 or equal to 3.28                 III=RCB index greater than 3.28
# 9                                               3
# Unavailable
# 1
table(xx$`Relapse-free_survival_indicator_(1=event;_local_or_distant_progression_or_death,_0=censor_at_last_follow-up)`)
# 0  1
# 16  4
table(xx$`Survival_status_(7=alive;_8=dead,_9=lost)`)
# 7  8
# 19  1



xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = str_c("G", `Histologic_grade_(1=grade_i_(low);__2=_grade_ii_(intermediate);_3=_grade_iii_(high);_4=indeterminate)`),
    Grade = if_else(Grade == "G4", NA_character_, Grade),
    Histology = purrr::map_chr(
      `Histology_(1=necrosis;_2=ductal_carcinoma;_3=lobular;_4=mixed_ductal/lobular_carcinoma;_5=other_6=no_invasive_tumor_present;)`,
      ~(case_when( .x == 1 ~ "necrosis",
                   .x == 2 ~ "ductal_carcinoma",
                   .x == 3 ~ "lobular",
                   .x == 4 ~ "mixed_ductal/lobular_carcinoma",
                   .x == 5 ~ "other",
                   .x == 6 ~ "no_invasive_tumor_present"))
    ),
    Size = `Clinical_tumor_size_(pre-chemotherapy),_cm_(-1_=_n/a)`, # cm
    Size = if_else(Size == -1, NA_real_, Size),
    Size_bin = if_else(Size <=2 , "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when( .x == 0 ~ "T0",
                   .x <= 2 ~ "T1",
                   .x <= 5 ~ "T2",
                   .x > 4 ~ "T3"))
    ),
    Er = if_else(`Er_(0=negative;_1=positive;)` == 0, "neg", "pos"),
    Pr = if_else(`Pgr__(0=negative;_1=positive;)` == 0, "neg", "pos"),
    Her2 = if_else(`Her2_(0=negative;_1=positive;)` == 0, "neg", "pos"),
    Response_pathological = purrr::map_chr(
      `Pathological_complete_response_(pcr)`,
      ~(case_when(.x == "No" ~ "npCR",
                  .x == "Yes" ~ "pCR"))
    ),
    Response = Response_pathological,
    Event_dfs = `Relapse-free_survival_indicator_(1=event;_local_or_distant_progression_or_death,_0=censor_at_last_follow-up)`,
    Time_dfs = `Relapse-free_survival_time_–_time_from_chemo_start_date_until_earliest_[local_or_distant_progression_or_death_(time_unit_is_days)]` / 365,

    Regimen = "neoadj_adj",
    Arm_neoadj = purrr::map_chr(
      `Neoadjuvant_chemotherapy*_(1_=_ac_only;_2_=_ac_+_t;_3_=_ac_+_t_+_herceptin;_4_=_ac_+_t_+_other;_[a_=_doxorubicin,_c_=_cyclophosphamide,_t_=_taxane])`,
      ~(case_when(.x == 1 ~ "AC",
                  .x == 2 ~ "Taxane+AC",
                  .x == 3 ~ "Taxane+AC+Trastuzumab",
                  .x == 4 ~ "Taxane+AC+Other"))
    ),
    Arm_adj = "[±Chemo±Hormone±Trastuzumab]"
    # Chemo in Arm_adj is imputed from related publication, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3434983/, under "study treatment and procedures subheading"

  )

table(xx1$Arm_neoadj)
# axane+AC Taxane+AC+Trastuzumab
# 17                     3
table(xx1$Arm_adj)
# [±Chemo±Hormone±Trastuzumab]
# 20

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE22226_GPL4133_cleaned.tsv"
)


# GSE22226_GPL1708
# >>>>>>>>>>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE22226_GPL1708.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 130


table(xx$`Histologic_grade_(1=grade_i_(low);__2=_grade_ii_(intermediate);_3=_grade_iii_(high);_4=indeterminate)`)
# 1  2  3  4
# 7 52 69  1
table(xx$`Histology_(1=necrosis;_2=ductal_carcinoma;_3=lobular;_4=mixed_ductal/lobular_carcinoma;_5=other_6=no_invasive_tumor_present;)`)
# 1   2   3   4   5   6
# 1 105   9   2   9   1
table(xx$`Clinical_t_stage_(pre-chemotherapy)_(1_=_<=_2cm;_2_=_>2-5_cm;_3_=_>_5cm)`)
# 1  2  3
# 3 51 71
table(xx$`Er_(0=negative;_1=positive;)`)
# 0  1
# 58 65
table(xx$`Pgr__(0=negative;_1=positive;)`)
# 0  1
# 70 51
table(xx$`Her2_(0=negative;_1=positive;)`)
# 0  1
# 73 39
table(xx$ `Neoadjuvant_chemotherapy*_(1_=_ac_only;_2_=_ac_+_t;_3_=_ac_+_t_+_herceptin;_4_=_ac_+_t_+_other;_[a_=_doxorubicin,_c_=_cyclophosphamide,_t_=_taxane])` )
# 1   2   3   4
# 4 112  11   2
table(xx$`Pathological_complete_response_(pcr)` )
# N/A  No Yes
# 5  92  32
table(xx$Rcb_class)
# 0=RCB index 0          I=RCB index less than or equal to 1.36
# 32                                               9
# II=RCB index greater than 1.36 or equal to 3.28                 III=RCB index greater than 3.28
# 50                                              23
# Unavailable
# 13
table(xx$`Relapse-free_survival_indicator_(1=event;_local_or_distant_progression_or_death,_0=censor_at_last_follow-up)`)
# 0  1
# 89 40
table(xx$`Survival_status_(7=alive;_8=dead,_9=lost)`)
# 7  8  9
# 98 27  4



xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = str_c("G", `Histologic_grade_(1=grade_i_(low);__2=_grade_ii_(intermediate);_3=_grade_iii_(high);_4=indeterminate)`),
    Grade = if_else(Grade == "G4", NA_character_, Grade),
    Histology = purrr::map_chr(
      `Histology_(1=necrosis;_2=ductal_carcinoma;_3=lobular;_4=mixed_ductal/lobular_carcinoma;_5=other_6=no_invasive_tumor_present;)`,
      ~(case_when( .x == 1 ~ "necrosis",
                   .x == 2 ~ "ductal_carcinoma",
                   .x == 3 ~ "lobular",
                   .x == 4 ~ "mixed_ductal/lobular_carcinoma",
                   .x == 5 ~ "other",
                   .x == 6 ~ "no_invasive_tumor_present"))
    ),
    Size = `Clinical_tumor_size_(pre-chemotherapy),_cm_(-1_=_n/a)`, # cm
    Size = if_else(Size == -1, NA_real_, Size),
    Size_bin = if_else(Size <=2 , "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when( .x == 0 ~ "T0",
                   .x <= 2 ~ "T1",
                   .x <= 5 ~ "T2",
                   .x > 4 ~ "T3"))
    ),
    Er = if_else(`Er_(0=negative;_1=positive;)` == 0, "neg", "pos"),
    Pr = if_else(`Pgr__(0=negative;_1=positive;)` == 0, "neg", "pos"),
    Her2 = if_else(`Her2_(0=negative;_1=positive;)` == 0, "neg", "pos"),
    Response_pathological = purrr::map_chr(
      `Pathological_complete_response_(pcr)`,
      ~(case_when(.x == "No" ~ "npCR",
                  .x == "Yes" ~ "pCR"))
      ),
    Response = Response_pathological,
    Event_dfs = `Relapse-free_survival_indicator_(1=event;_local_or_distant_progression_or_death,_0=censor_at_last_follow-up)`,
    Time_dfs = `Relapse-free_survival_time_–_time_from_chemo_start_date_until_earliest_[local_or_distant_progression_or_death_(time_unit_is_days)]` / 365,

    Regimen  = "neoadj_adj",
    Arm_neoadj = purrr::map_chr(
      `Neoadjuvant_chemotherapy*_(1_=_ac_only;_2_=_ac_+_t;_3_=_ac_+_t_+_herceptin;_4_=_ac_+_t_+_other;_[a_=_doxorubicin,_c_=_cyclophosphamide,_t_=_taxane])`,
      ~(case_when(.x == 1 ~ "AC",
                  .x == 2 ~ "Taxane+AC",
                  .x == 3 ~ "Taxane+AC+Trastuzumab",
                  .x == 4 ~ "Taxane+AC+Other"))
    ),
    Arm_adj = "[±Chemo±Hormone±Trastuzumab]"
    # Chemo in Arm_adj is imputed from related publication, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3434983/, under "study treatment and procedures subheading"

  )

table(xx1$Arm_neoadj)
# AC             Taxane+AC       Taxane+AC+Other Taxane+AC+Trastuzumab
# 4                   112                     2                    11
table(xx1$Arm_adj)
# [±Chemo±Hormone±Trastuzumab]
# 130


write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE22226_GPL1708_cleaned.tsv"
)



# GSE31863
# >>>>>>>>>>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE31863.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 143

# DNA index and S-phase fraction (SPF)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1971595/
# S_dnaind1/2/3/4
# S_spf_nbr_pop/S_spf_percent/S_spf_status_0low_1high


table(xx$`Series,_1_=_lr+/rt+,_2_=_lr-/rt+,_3_=_lr+/rt-,_4_=_lr-/rt-`)
# 1  2  3  4
# 30 47 22 44
table(xx$P_treat_chemo_0no_1yes)
# 0   1
# 140   3
table(xx$P_treat_horm_0no_1yes)
# 0   1
# 129  14
table(xx$P_treat_radio_0no_1yes)
# 0  1
# 66 77
table(xx$`Characteristics:_s_found_by_screening_0no_1yes:_ch1`)
# 0  1
# 64 57
table(xx$S_histology_type)
# Ductal_medullary    Ductal_mucous   Ductal_tubular       Ductal_UNS          Lobular
# 5                2                2               68                6
# n/a          Tubular     Tubuloductal    Tubulolobular
# 33               18                8                1
table(xx$S_meno)
# post Post  pre  Pre
# 2   84    9   37
table(xx$S_surgery_desc)
# 1   3
# 119   2
table(xx$S_surgery_extended_0no_1yes)
# 0   1
# 112   9
table(xx$S_er_pos_neg)
# er_neg er_pos
# 29    114
table(xx$S_pgr_pos_neg)
# pgr_neg pgr_pos
# 47      88
table(xx$P_cause_of_death)
# alive          brca         other other_disease
# 101            13             1             3
table(xx$P_dead_followup_0no_1yes)
# 0   1
# 111  14
table(xx$P_local_recurrence_0no_1yes)
# 0  1
# 91 52
table(xx$S_recurrence_0no_1yes)
# 0   1
# 140   3

table(xx$P_treat_radio_0no_1yes, xx$P_local_recurrence_0no_1yes)
#   0  1
# 0 44 22
# 1 47 30

xx %>%
  dplyr::group_by(P_treat_radio_0no_1yes, P_local_recurrence_0no_1yes,
                  P_treat_horm_0no_1yes, P_treat_chemo_0no_1yes) %>%
  dplyr::summarise(N = n())
#   P_treat_radio_0no_… P_local_recurrence_0… P_treat_horm_0no_… P_treat_chemo_0no…     N
# 1                   0                     0                  0                  0    39
# 2                   0                     0                  0                  1     1
# 3                   0                     0                  1                  0     4
# 4                   0                     1                  0                  0    18
# 5                   0                     1                  0                  1     1
# 6                   0                     1                  1                  0     2
# 7                   0                     1                  1                  1     1
# 8                   1                     0                  0                  0    46
# 9                   1                     0                  1                  0     1
# 10                   1                     1                  0                  0    24
# 11                   1                     1                  1                  0     6

xx1 <- xx %>%
  dplyr::mutate(
    Age = S_age_year,
    Size = S_tumor_size_mm %>% as.numeric()  / 10, # cm
    Size = if_else(Size == -1, NA_real_, Size),
    Size_bin = if_else(Size <=2 , "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when( .x == 0 ~ "T0",
                   .x <= 2 ~ "T1",
                   .x <= 5 ~ "T2",
                   .x > 4 ~ "T3"))
    ),
    Er = S_er_pos_neg %>% str_replace("er_", ""),
    Pr = S_pgr_pos_neg %>% str_replace("pgr_", ""),
    Histology = xx$S_histology_type %>% str_replace("n/a", NA_character_),
    Menopause = S_meno,
    Therapy_chemo = if_else(P_treat_chemo_0no_1yes == 0, "no", "yes"),
    Therapy_hormone = if_else(P_treat_horm_0no_1yes == 0, "no", "yes"),
    Therapy_radio = if_else(P_treat_radio_0no_1yes == 0, "no", "yes"),
    Event_dfs = P_local_recurrence_0no_1yes,
    Time_dfs = if_else(P_days_to_recurrence == -1,  P_days_followup, P_days_to_recurrence),
    Time_dfs = Time_dfs / 365, # years

    Regimen = "adj",
    Arm_adj = if_else(Therapy_radio == "no", "Untreated", "Radiotherapy"),
    Arm_adj = if_else(Therapy_hormone == "no", Arm_adj, NA_character_), # filtering hormone treated
    Arm_adj = if_else(Therapy_chemo == "no", Arm_adj, NA_character_) # filtering chemo treated
  )

table(xx1$Arm_adj)

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE31863_cleaned.tsv"
)



# GSE45255
# >>>>>>>>>>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE45255.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 139

table(xx$Ln_status)
# LN- LN+
#   94  45
table(xx$Er_status)
# ER- ER+
#   48  89
table(xx$Pgr_status)
# PgR- PgR+
#   60   72
table(xx$Her2_status)
# He- He+
#   69  62
table(xx$Histological_grade)
# G1 G2 G3
# 17 52 67
table(xx$Histology)
# invasive ductal and lobular carcinoma                         invasive ductal carcinoma
# 1                                                85
# invasive ductal carcinoma and papillary carcinoma                        invasive lobular carcinoma
# 1                                                 4
# medullary carcinoma                                mucinous carcinoma
# 2                                                 1
# papillary carcinoma
# 1
table(xx$`Adjuvant_treated?_(0=no,_1=yes)`)
# 0   1
# 10 121
table(xx$Characteristics)
# endocrine? (0=no, 1=yes): 0  endocrine? (0=no, 1=yes): 1 endocrine? (0=no, 1=yes): NA
# 54                           77                            8
table(xx$`Chemo?_(0=no,_1=yes)`)
# 0  1
# 55 76
table(xx$Treatment_type)
# AC                ACx3/CMFx1                      ACx4               ACx4 cycles
# 2                         1                         2                        18
# ACx4; Arimidex/Tamoxifen           ACx4; Tamoxifen           ACx6; Tamoxifen        ACx8 cycles; taxol
# 1                         1                         1                         1
# anastrozole;taxol        Arimidex/Tamoxifen     Arimidex+CMFx2 cycles                       CAF
# 1                         1                         1                         9
# CAFx6 cycles                       CMF          CMFx2; Tamoxifen                     CMFx6
# 2                         4                         1                         4
# CMFx6; Arimidex/Tamoxifen          CMFx6; Tamoxifen           ECx4; Tamoxifen                      none
# 1                         3                         1                        10
# Tam+AC           Tam+ACx4 cycles            Tam+ACx8;taxol          Tam+CAFx2 cycles
# 1                         8                         1                         1
# Tam+CAFx6 cycles          Tam+CMFx6 cycles                 tamoxifen                 Tamoxifen
# 1                         7                        23                        20
# tamoxifen;goserelin Tamoxifen+Zoledronic acid
# 1                         3
table(xx$`Dfs_event_(defined_as_any_type_of_recurrence_or_death_from_breast_cancer)`)
# 0  1
# 73 22
table(xx$`Dmfs_event_(defined_as_distant_metastasis_or_death_from_breast_cancer)`)
# 0   1
# 104  32
table(xx$`Dss_event_(defined_as_death_from_breast_cancer)`)
# 0   1
# 116  18



xx1 <- xx %>%
  dplyr::mutate(
    Node_bin = if_else(Ln_status == "LN-", "neg", "pos"),
    Er =  if_else(Er_status == "ER-", "neg", "pos"),
    Pr =  if_else(Pgr_status == "PgR-", "neg", "pos"),
    Her2 =  if_else(Her2_status == "He-", "neg", "pos"),
    Grade = Histological_grade,
    Size =`Size_(mm)`  / 10, # cm
    Size_bin = if_else(Size <=2 , "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when( .x == 0 ~ "T0",
                   .x <= 2 ~ "T1",
                   .x <= 5 ~ "T2",
                   .x > 4 ~ "T3"))
    ),
    Age = Patient_age,
    Histology = Histology,
    Therapy_adjuvant = if_else(`Adjuvant_treated?_(0=no,_1=yes)` == 0, "no", "yes"),
    Therapy_hormone = purrr::map_chr(
      Characteristics,
      ~(case_when(.x == "endocrine? (0=no, 1=yes): 0" ~ "no",
                  .x == "endocrine? (0=no, 1=yes): 1" ~ "yes"))
    ),
    Therapy_chemo = if_else(`Chemo?_(0=no,_1=yes)` == 0, "no", "yes"),
    Event_dfs = `Dmfs_event_(defined_as_distant_metastasis_or_death_from_breast_cancer)`,
    Time_dfs = Dmfs_time, # years

    Regimen = "adj",
    Arm_adj = Treatment_type
  )


yy <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE45255_cleaned_old.tsv")
names(yy) <- clean_names(names(yy))
glimpse(yy) # 139

identical(xx1$Sample_geo_accession, yy$Sample_geo_accession) # TRUE

length(table(xx1$Arm_adj)) # 30
length(table(yy$Arm)) # 18

# resetting arm
xx1$Arm_adj <- yy$Arm %>% str_replace("adjuvant_","")

table(xx1$Arm_adj)
# AC                 AC_CMF          AC_Paclitaxel        Anastrozole_CMF
# 22                      1                      1                      1
# Anastrozole_Paclitaxel                    CAF                    CMF    No_systemic_therapy
# 1                     11                      8                     10
# Tam                 Tam_AC        Tam_Anastrozole     Tam_Anastrozole_AC
# 43                     12                      1                      1
# Tam_Anastrozole_CMF                Tam_CAF                Tam_CMF                 Tam_EC
# 1                      2                     11                      1
# Tam_Goserelin     Tam_ZoledronicAcid
# 1                      3

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE45255_cleaned.tsv"
)


# GSE69031
# >>>>>>>>>>>>>>>>

xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE69031.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 139


table(xx$Sex)
# female
# 130
table(xx$Organism_part)
# mammary gland
# 130
table(xx$Ethnicity )
# African American        Caucasian            Other
# 21               94               11
table(xx$Tumor_grading)
# 1  2  3
# 14 46 65
table(xx$Node_positive)
# no yes
# 59  71

table(xx$Estrogen_receptor_status)
# negative positive
# 46       84
table(xx$Progesterone_receptor_status)
# n/a negative positive
# 2       54       74
table(xx$`Erbb2_positive_(ihc)`)
# n/a  no yes
# 41  78  11
table(xx$`P53_positive_(ihc)`)
# n/a  no yes
# 39  67  24

table(xx$Recurrence)
# n/a  no yes
# 1  90  39
table(xx$Distal_recurrence)
# n/a  no yes
# 1 102  27
table(xx$Dead_of_disease)
# n/a  no yes
# 1 100  29
table(xx$Alive_at_endpoint)
# n/a  no yes
# 1  84  45

table(xx$Disease_state)
# breast cancer
# 130
table(xx$Hormonal_therapy)
# n/a  no yes
# 3  53  74
table(xx$Radiation_treatment)
# n/a  no yes
# 2  61  67
table(xx$Chemotherapy_treatment)
# n/a  no yes
# 2  60  68



xx1 <- xx %>%
  dplyr::mutate(
    Age = Age_at_diagnosis,
    Grade = str_c("G", Tumor_grading),
    Node_bin = purrr::map_chr(
      Node_positive,
      ~(case_when(.x == "no" ~ "neg",
                  .x == "yes" ~ "pos"))
    ),
    Size = `Tumor_size_(mm)` %>% as.numeric() / 10,
    Size_bin = if_else(Size <=2 , "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when( .x == 0 ~ "T0",
                   .x <= 2 ~ "T1",
                   .x <= 5 ~ "T2",
                   .x > 4 ~ "T3"))
    ),


    Er = purrr::map_chr(
      Estrogen_receptor_status,
      ~(case_when(.x == "negative" ~ "neg",
                  .x == "positive" ~ "pos"))
    ),
    Pr = purrr::map_chr(
      Progesterone_receptor_status,
      ~(case_when(.x == "negative" ~ "neg",
                  .x == "positive" ~ "pos"))
    ),
    Her2 = purrr::map_chr(
      `Erbb2_positive_(ihc)`,
      ~(case_when(.x == "no" ~ "neg",
                  .x == "yes" ~ "pos"))
    ),
    Ihc_p53 = purrr::map_chr(
      `P53_positive_(ihc)`,
      ~(case_when(.x == "no" ~ "neg",
                  .x == "yes" ~ "pos"))
    ),

    Gender = Sex,
    Ethnicity = Ethnicity,

    Event_dfs = purrr::map_dbl(
      Recurrence,
      ~(case_when(.x == "no" ~ 0,
                  .x == "yes" ~ 1))
    ),
    Time_dfs = Recurrence_time %>% as.numeric(),

    Therapy_radio = if_else(Radiation_treatment == "n/a",
                            NA_character_, Radiation_treatment),
    Therapy_hormone = if_else(Hormonal_therapy == "n/a",
                            NA_character_, Hormonal_therapy),
    Therapy_chemo = if_else(Chemotherapy_treatment == "n/a",
                            NA_character_, Chemotherapy_treatment),

    # Methods
    # "About half of the tumors were node positive, 67% were estrogen receptor positive, 60% received tamoxifen, and half received adjuvant chemotherapy (typically adriamycin and cytoxan). "

    Regimen = "adj",
    Arm_adj = str_c(Therapy_radio, Therapy_hormone, Therapy_chemo) %>%
      purrr::map_chr(
        ~(case_when(
          .x == "nonono" ~ "Untreated",
          .x == "nonoyes" ~ "AC",
          .x == "noyesno" ~ "Tamoxifen",
          .x == "noyesyes" ~ "Tamoxifen+AC",
          .x == "yesnono" ~ "Radio",
          .x == "yesnoyes" ~ "Radio+AC",
          .x == "yesyesno" ~ "Radio+Tamoxifen",
          .x == "yesyesyes" ~ "Radio+Tamoxifen+AC"
        ))
      )
  )

table(xx1$Arm_adj)
# AC              Radio           Radio+AC    Radio+Tamoxifen Radio+Tamoxifen+AC
# 16                  8                 22                 21                 15
# Tamoxifen       Tamoxifen+AC          Untreated
# 24                 14                  7
write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE69031_cleaned.tsv"
)



# GSE16446 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE16446.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 139

table(xx$Agebin)
# 0  1
# 70 50
table(xx$`T`)
# T1 T2 T3 T4
# 17 83  5 15
table(xx$N)
# N0 N1 N2 N3
# 55 60  3  2
table(xx$Grade)
# 1  2  3
# 2 20 92
table(xx$Her2fishbin)
# 0  1
# 62 31
table(xx$Top2atri)
# -1  0  1
# 11 69 12
table(xx$Topoihc)
# <10    0    1   10   15 17.5    2   20   25   30   40    5   50   60   70   80   90
# 2    3    1   18    8    1    5   12    6    7    4    9    4    1    2    3    2
table(xx$Pcr)
# 0  1
# 98 16
table(xx$Dmfs_event)
# 0  1
# 89 25
table(xx$Os_event)
# 0  1
# 98 16


xx1 <- xx %>%
  dplyr::mutate(
    Age_bin = if_else(Agebin == 0, "young", "old"),
    Age_detailed = "young:<=50,old>50",
    Grade = str_c("G", Grade),
    Node_cat = N,
    Node_bin = if_else(Node_cat == "N0", "neg", "pos"),
    Size_cat = `T`,
    Size_bin = if_else(Size_cat == "T1" , "small", "large"),

    Er = "neg", # from lterature
    Her2 = if_else(Her2fishbin == 0, "neg", "pos"),
    Her2_score.fish = Her2fish,
    Top2a = purrr::map_chr(
      Top2atri,
      ~(case_when(.x == -1 ~ "Deletion",
                  .x ==  0 ~ "Normal",
                  .x ==  1 ~ "Amplification"))
    ),
    Topo_score = Topoihc,

    Event_dfs = Dmfs_event,
    Time_dfs = Dmfs_time / 365, # years

    Response_pathological = if_else(Pcr == 0, "npCR", "pCR"),
    Response = Response_pathological,

    # "The neoadjuvant Trial of Principle (TOP) study, in which patients with estrogen receptor (ER) –negative tumors were treated with anthracycline (epirubicin) monotherapy, was specifically designed to evaluate the predictive value of topoisomerase II-α (TOP2A) and develop a gene expression signature to identify those patients who do not benefit from anthracyclines."

    Regimen = "neoadj",
    Arm_neoadj = "Epirubicin"
  )

table(xx1$Arm_neoadj)
# Epirubicin
# 120

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE16446_cleaned.tsv"
)



# GSE32646 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE32646.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 115

table(xx$Clinical_t_stage)
# 1  2  3  4
# 5 87 18  5
table(xx$Lymph_node_status)
# negative positive
# 32       83
table(xx$Histological_grade)
# 1  2  3
# 16 78 21
table(xx$Er_status_ihc)
# negative positive
# 44       71
table(xx$Pr_status_ihc)
# negative positive
# 70       45
table(xx$Her2_status_fish)
# negative positive
# 81       34
table(xx$Pathologic_response_pcr_ncr)
# nCR pCR
# 88  27


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = str_c("G", Histological_grade),
    Node_bin = if_else(Lymph_node_status == "negative", "neg", "pos"),
    Size_cat = str_c("T", Clinical_t_stage),
    Size_bin = if_else(Size_cat == "T1" , "small", "large"),

    Er = if_else(Er_status_ihc == "negative", "neg", "pos"),
    Pr = if_else(Pr_status_ihc == "negative", "neg", "pos"),
    Her2 = if_else(Her2_status_fish == "negative", "neg", "pos"),

    Response_pathological = if_else(Pathologic_response_pcr_ncr == "nCR", "npCR", "pCR"),
    Response = Response_pathological,


    # "Primary breast cancer patients (n  =  123, T1‐4b N0‐1 M0) who were consecutively recruited for the present study had been treated with NAC consisting of paclitaxel (80 mg/m2) weekly for 12 cycles followed by 5‐FU (500 mg/m2), epirubicin (75 mg/m2) and cyclophosphamide (500 mg/m2) every 3 weeks for four cycles (paclitaxel followed by 5‐fluorouracil/epirubicin/cyclophosphamide [P‐FEC]) at Osaka University Hospital between 2004 and 2010. "

    Regimen = "neoadj",
    Arm_neoadj = "Paclitaxel+FEC"
  )

table(xx1$Arm_neoadj)
# Paclitaxel+FEC
# 115

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE32646_cleaned.tsv"
)



# GSE19615 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE19615.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 115

table(xx$`Tumor_recurrence_(36mo)`)
# no yes
# 100  15
table(xx$Histology_type)
# Ductal Lobular   Mixed
# 89      14      12
table(xx$`Grade_(modified,_bloom,_richardson)`)
# I  II III
# 23  28  64
table(xx$Er)
# neg     pos pos-low
# 45      66       4
table(xx$Pr)
# neg     pos pos-low
# 51      52      12
table(xx$Her.2)
# low pos (2+)          neg     pos (3+)
# 6           79           30
table(xx$Lymph_nodes)
# negative pos.micromet     positive
# 62            2           51
table(xx$Adjuvant_chemotherapy)
# AC x 2, Taxol                      AC x 4                    AC-taxol
# 1                          39                          25
# AC-taxol.taxotere                AC-taxol(x1)          AC-taxol/Herceptin
# 1                           1                           3
# AC-taxotere                         CAF                         CMF
# 3                           3                           4
# Epirubicin, Cytoxan, Xeloda                        none                     Unknown
# 1                          28                           6
table(xx$Chemo_class)
# Anthracycline-based                none               Other         Trastuzumab             Unknown
# 74                  28                   4                   3                   6
table(xx$Hormonal_rx)
# Arimidex     none      Tam  unknown
# 2       47       62        4
table(xx$`Distant_recur_(yn)`)
# N   Y
# 101  14


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = purrr::map_chr(
      `Grade_(modified,_bloom,_richardson)`,
      ~(case_when(.x == "I" ~ "G1",
                  .x == "II" ~ "G2",
                  .x == "III" ~ "G3"))
    ),
    Node_bin = if_else(Lymph_nodes == "negative", "neg", "pos"),
    Size = `Tumor_size_(cm)`,
    Size_bin = if_else(Size <=2 , "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when(.x == 0 ~ "T0",
                  .x <= 2 ~ "T1",
                  .x <= 5 ~ "T2",
                  .x > 5 ~ "T3"))
    ),
    Er = if_else(Er == "neg", "neg", "pos"),
    Pr = if_else(Pr == "neg", "neg", "pos"),
    Her2 = if_else(Her.2 == "neg", "neg", "pos"),

    Event_dfs = if_else(`Distant_recur_(yn)` == "N", 0 ,1),
    Time_dfs =`Distant_recurrence_free_survival_(mo)` / 12, # years

    Histology = Histology_type,


    Therapy_chemo_name = str_c(
      if_else(str_detect(Adjuvant_chemotherapy, "Taxol"), "Paclitaxel", ""),
      if_else(str_detect(Adjuvant_chemotherapy, "taxol"), "Paclitaxel", ""),
      if_else(str_detect(Adjuvant_chemotherapy, "taxotere"), "Docetaxel", ""),
      if_else(str_detect(Adjuvant_chemotherapy, "Herceptin"), "Trastuzumab", ""),
      if_else(str_detect(Adjuvant_chemotherapy, "Epirubicin, Cytoxan, Xeloda"), "Capecitabine+EC", ""),
      if_else(str_detect(Adjuvant_chemotherapy, "AC"), "AC", ""),
      if_else(str_detect(Adjuvant_chemotherapy, "CAF"), "CAF", ""),
      if_else(str_detect(Adjuvant_chemotherapy, "CMF"), "CMF", ""),
      if_else(Adjuvant_chemotherapy == "none", "no_therapy", ""),
      if_else(str_detect(Adjuvant_chemotherapy, "Unknown"), NA_character_, "")
    ) %>%
      str_replace("Paclitaxel","Paclitaxel+") %>%
      str_replace("Docetaxel","Docetaxel+") %>%
      str_replace("Trastuzumab","Trastuzumab+"),
    Therapy_hormone_name = purrr::map_chr(
      Hormonal_rx,
      ~(case_when(.x == "Arimidex" ~ "Anastrozole",
                  .x == "Tam" ~ "Tamoxifen",
                  .x == "unknown" ~ NA_character_,
                  .x == "none" ~ "no_therapy"))
    ),

    Regimen = "adj",
    Arm_adj = str_c(Therapy_chemo_name, Therapy_hormone_name, sep=" ") %>%
      str_replace(" no_therapy","") %>%
      str_replace("no_therapy ","") %>%
      str_trim() %>%
      str_replace(" ","+"),
    Arm_adj = if_else(Arm_adj == "no_therapy", "No_systemic_therapy", Arm_adj)

  )


table(xx1$Arm_adj)
# AC                      AC+Anastrozole
# 22                                   1
# AC+Tamoxifen                         Anastrozole
# 16                                   1
# CAF                       CAF+Tamoxifen
# 1                                   2
# Capecitabine+EC                       CMF+Tamoxifen
# 1                                   4
# Docetaxel+AC              Docetaxel+AC+Tamoxifen
# 1                                   1
# No_systemic_therapy                       Paclitaxel+AC
# 7                                  12
# Paclitaxel+AC+Tamoxifen   Paclitaxel+Docetaxel+AC+Tamoxifen
# 15                                   1
# Paclitaxel+Trastuzumab+AC Paclitaxel+Trastuzumab+AC+Tamoxifen
# 1                                   2
# Tamoxifen
# 20

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE19615_cleaned.tsv"
)




# GSE130786 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE130786.tsv")
names(xx) <- clean_names(names(xx))
names(xx)[10] <- "Sample_type.1" # redundant coullmn name
glimpse(xx) # 110

# Ref: https://cancerres.aacrjournals.org/content/73/24_Supplement/S1-02

table(xx$`Er_status_(ihc_staining_results)`)
# Neg Pos
# 48  62
table(xx$`Pr_status_(ihc_staining_results)`)
# Neg Pos
# 64  46
table(xx$Treatment_group)
# TCH TCHTy  TCTy
# 31    50    29
table(xx$Drug_response)
# PCR  RD
# 47  63



xx1 <- xx %>%
  dplyr::mutate(
    Er = if_else(`Er_status_(ihc_staining_results)` == "Neg", "neg", "pos"),
    Pr = if_else(`Pr_status_(ihc_staining_results)` == "Neg", "neg", "pos"),
    Her2 = "pos",

    Response_pathological = if_else(Drug_response == "RD", "npCR", "pCR"),
    Response = Response_pathological,

    Regimen = "neoadj",
    Arm_neoadj = purrr::map_chr(
      Treatment_group,
      ~(case_when(.x == "TCH" ~ "Docetaxel+Carboplatin+Trastuzumab",
                  .x == "TCHTy" ~ "Docetaxel+Carboplatin+Trastuzumab+Lapatinib",
                  .x == "TCTy" ~ "Docetaxel+Carboplatin+Lapatinib"))
    )
  )

table(xx1$Arm_neoadj)
# Docetaxel+Carboplatin+Lapatinib           Docetaxel+Carboplatin+Trastuzumab
# 29                                          31
# Docetaxel+Carboplatin+Trastuzumab+Lapatinib
# 50


write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE130786_cleaned.tsv"
)



# GSE22093 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE22093.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 110

table(xx$Tissue)
# breast cancer FNA biopsy
# 103
table(xx$Pcr.v.rd)
# pCR  RD
# 28  69
table(xx$P53_status)
# MUT n/a  WT
# 58   3  42
table(xx$Er_positive_vs_negative_by_immunohistochemistry)
# ERneg ERpos
# 56    42
table(xx$Prechemo_t)
# 0  1  2  3  4
# 1  2 51 26 18
table(xx$Prechemo_n)
# 0  1  2  3
# 21 16 10  4
table(xx$Bmn.grade)
# 1  2  3
# 3 29 47


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = str_c("G", Bmn.grade),
    Node_bin = if_else(Prechemo_n == 0, "neg", "pos"),
    Size_bin = if_else(Prechemo_t <= 1, "small", "large"),
    Size_cat = str_c("T", Prechemo_t),

    Er = if_else(Er_positive_vs_negative_by_immunohistochemistry == "ERneg", "neg", "pos"),
    Her2 = "neg",
    Mut_p53 = purrr::map_chr(
      P53_status,
      ~(case_when(.x == "MUT" ~ "Mut",
                  .x == "WT" ~ "Wt"))
    ),
    Response_pathological = if_else(Pcr.v.rd == "RD", "npCR", "pCR"),
    Response = Response_pathological,

    Regimen = "neoadj",
    Arm_neoadj = "FAC|FEC"
  )

table(xx1$Arm_neoadj)

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE22093_cleaned.tsv"
)



# GSE4779 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE4779.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 102

table(xx$Label)
# npCR  pCR
# 63   39
table(xx$Estrogen_receptor)
# En Ep
# 65 37
table(xx$Progesterone_receptor)
# P- Pn Pp
# 2 75 25
table(xx$Size)
# T- T1 T2 T3
# 3  2 63 34
table(xx$Node)
# N- N0 N1 N2
# 3 37 55  7


xx1 <- xx %>%
  dplyr::mutate(
    Node_cat = if_else(Node == "N-", NA_character_, Node),
    Node_bin = if_else(Node_cat == "N1", "neg", "pos"),
    Size_cat = if_else(Size == "T-", NA_character_, Size),
    Size_bin = if_else(Size_cat == "T1", "small", "large"),

    Er = if_else(Estrogen_receptor == "En", "neg", "pos"),
    Pr = if_else(Progesterone_receptor == "P-", NA_character_,
                 if_else(Progesterone_receptor == "Pn", "neg", "pos")),

    Response_pathological = Label,
    Response = Response_pathological,
# "The samples were taken from the FEC arm (5-fluorouracil, epirubicin, cyclophosphamide) of the EORTC 10994 trial. EORTC 10994 is a phase III clinical trial comparing FEC with ET (epirubicin, docetaxel) in patients with large operable, locally advanced or inflammatory breast cancer (www.eortc.be)."
    Regimen = "neoadj",
    Arm_neoadj = "FEC"
  )

table(xx1$Arm_neoadj)


write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE4779_cleaned.tsv"
)





# GSE114403 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE114403.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 100

table(xx$`Er/pr_status_(1_if_er+_or_pr+;_0_if_er-_and_pr-)`)
# 0  1
# 26 74
table(xx$`Pcr_status_(1_if_pcr;_0_if_residual_disease_or_progressed_before_surgery)`)
# 0  1
# 74 26
table(xx$`Treatment_arm_(1_if_assigned_to_bevacizumab_(arm_a);_0_if_assigned_to_control_(arm_b_+_c))`)
# 0  1
# 55 45

table(xx$Sample_source_name_ch1)
# post-treatment breast cancer  pre-treatment breast cancer
# 50                           50

summary(xx$Til.count.reader.1)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#    0.00    6.25   15.00   22.35   30.00   90.00       2
summary(xx$Til.count.reader.2)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#    0.00    0.00    5.00   14.67   20.00   90.00       2
summary(xx$Pd.l1.ihc.tumor....r.1)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#   0.000   0.000   0.000   1.614   0.000  20.000      17
summary(xx$Pd.l1.ihc.tumor....r.2)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#   0.000   0.000   0.000   1.506   0.000  30.000      17
summary(xx$Pd.l1.ihc.ic....r.1)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#  0.000   0.000   0.000   4.639   5.000  40.000      17
summary(xx$Pd.l1.ihc.ic....r.1.1)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#   0.000   0.000   0.000   2.482   1.000  40.000      17

xx1 <- xx %>%
  dplyr::mutate(

    Hr = if_else(`Er/pr_status_(1_if_er+_or_pr+;_0_if_er-_and_pr-)` == 0, "neg", "pos"),
    Her2 = "neg", # from text

    Response_pathological = if_else(`Pcr_status_(1_if_pcr;_0_if_residual_disease_or_progressed_before_surgery)` == 0, "npCR", "pCR"),
    Response = Response_pathological,

    Timepoint = if_else(Sample_source_name_ch1 == "pre-treatment breast cancer",
                        "Pre_treatment", "Post_treatment"),

    Til_stromal = (xx$Til.count.reader.1 + xx$Til.count.reader.2) / 2,
    Pdl1_stromal = (Pd.l1.ihc.tumor....r.1 + Pd.l1.ihc.tumor....r.2) /2,
    Pdl1_tumor = (Pd.l1.ihc.ic....r.1 + Pd.l1.ihc.ic....r.1.1) / 2,

    # Introduction
    # "S0800 (NCT00856492) was a randomized, 3-arm, Phase II trial that tested if inclusion of bevacizumab with neoadjuvant chemotherapy could improve pCR rates in HER2-negative, locally advanced, or inflammatory breast cancer (IBC). The three arms of the trial were weekly nab-paclitaxel and bevacizumab followed by dose-dense doxorubicin/cyclophosphamide (ddAC) (Arm A), nab-paclitaxel followed by ddAC, (Arm B), and ddAC followed by nab-paclitaxel (Arm C). Patients were randomly allocated in 2:1:1 ratio to arms A, B and C, respectively. For the primary efficacy analysis, the arms B and C were combined."

    Regimen = "neoadj",
    Arm_neoadj = if_else(`Treatment_arm_(1_if_assigned_to_bevacizumab_(arm_a);_0_if_assigned_to_control_(arm_b_+_c))` == 0,
                         "Nab_paclitaxel+AC", "Nab_paclitaxel+Bevacizumab+AC"),
    Arm_neoadj = if_else(Timepoint == "Post_treatment", NA_character_, Arm_neoadj), # filtering post treatment profiles

  )

table(xx1$Arm_neoadj)
# Nab_paclitaxel+AC Nab_paclitaxel+Bevacizumab+AC
# 27                            23

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE114403_cleaned.tsv"
)


# GSE76360 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE76360.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 100

table(xx$Patient_status)
# HER2+ Breast Cancer
# 100
table(xx$Timepoint)
# baseline     post
# 50       50
table(xx$Response_at_surgery)
# NOR OBJR  pCR
# 12   66   18
table(xx$Er_status)
# Neg Pos
# 56  44
table(xx$Pr_status)
# Neg Pos
# 72  28


xx1 <- xx %>%
  dplyr::mutate(

    Er = if_else(Er_status == "Neg", "neg", "pos"),
    Pr = if_else(Pr_status == "Neg", "neg", "pos"),
    Her2 = "pos",

    Response_pathological = if_else(Response_at_surgery == "pCR", "pCR", "npCR"),
    Response = Response_pathological,

    Timepoint = if_else(Timepoint == "baseline",
                        "Pre_treatment", "Post_treatment"),

    # Methods
    # "Tissue samples were collected from patients enrolled on two independent multicenter phase II neoadjuvant trials (03-311, NCT00148668; and BrUOG 211B, NCT00617942) for clinical stage II–III HER2+ breast cancer, with HER2 positivity defined as either overexpression by immunohistochemical (IHC) stain (3+) or a FISH ratio for HER/CEP17 of >2.0."
    # 03-311:
    # "The 03-311 trial was conducted at the Dana-Farber Cancer Institute and the Yale Comprehensive Cancer Center (YCCC). Once a patient signed informed consent, pretreatment tumor biopsies were obtained, and then she received a loading dose of trastuzumab (8 mg/kg). Repeat tumor biopsies were obtained approximately 14 days later, after which the patient started treatment with either vinorelbine 25 mg/m2 weekly × 12 weeks (NH) or docetaxel 75 mg/m2 and carboplatin AUC 6 every 3 weeks (TCH) with trastuzumab 2 mg/kg weekly for 12 weeks (TCH); the choice of chemotherapy regimen and subsequent surgical management of the breast and axilla were at the discretion of the treating physicians.
    # BrUOG 211B:
    # "BrUOG 211B (211B) was conducted by the Brown University Oncology Group (BrUOG) at its participating hospitals and at the YCCC and the City of Hope Comprehensive Cancer Center (COHCCC). Enrolled patients underwent baseline tumor biopsies, then received “run-in” treatment determined by their institution, with patients enrolled at a BrUOG institution receiving a loading dose (6 mg/kg) of trastuzumab while patients enrolled at YCCC or COHCCC received two weekly doses of nab-paclitaxel 100 mg/m2. Repeat tumor biopsies were obtained between 10 and 14 days after trastuzumab or the first dose of nab-paclitaxel. All patients were then treated with carboplatin AUC 6 every 3 weeks and nab-paclitaxel 100 mg/m2 and trastuzumab 2 mg/kg weekly for 18 weeks, followed by surgical management of the breast and axilla at the discretion of their treating physicians. "

    Regimen = "neoadj",
    Arm_neoadj = "Trastuzumab",
    Arm_neoadj = if_else(Timepoint == "Post_treatment", NA_character_, Arm_neoadj) # filtering
  )

table(xx1$Arm_neoadj)
# Trastuzumab
# 50

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE76360_cleaned.tsv"
)





# GSE21997_GPL1390 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE21997_GPL1390.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 35


table(xx$Study)
# Madrid
# 35
table(xx$`Histologic_grade_(1=grade_i_(low);__2=_grade_ii_(intermediate);_3=_grade_iii_(high);_9_=_n/a)`)
# 1  2  3
# 2 24  9
table(xx$`Histology_(1=necrosis;_2=ductal_carcinoma;_3=lobular;_4=mixed_ductal/lobular_carcinoma;_5=other__6=no_invasive_tumor_present;_9_=_n/a_)`)
# 2  3  5
# 30  3  2
table(xx$`Tumor_size_(-1_=_n/a;__1_=_<=_2cm;_2_=_2-5_cm;_3_=_>_5cm)`)
# 1  2  3
# 1 10 24
table(xx$`Er_(0=negative;_1=positive;_9_=_n/a)`)
# 0  1
# 15 20
table(xx$`Pgr__(0=negative;_1=positive;_9_=_n/a)`)
# 0  1
# 21 14
table(xx$`Her2_(0=negative;_1=positive;_9_=_n/a)`)
# 0  1
# 26  9
table(xx$Neoadjuvant_chemotherapy)
# Doxorubicin      Taxane
# 23          12
table(xx$Pcr)
# No Yes
# 30   5
table(xx$Residual_cancer_burden_index_class)
# 0   I  II III
# 5   3  13  14



xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = str_c("G", `Histologic_grade_(1=grade_i_(low);__2=_grade_ii_(intermediate);_3=_grade_iii_(high);_9_=_n/a)`),
    Size = if_else(`Tumor_size_(pre-chemotherapy),_cm_(-1_=_n/a)` == -1,
                   NA_real_, `Tumor_size_(pre-chemotherapy),_cm_(-1_=_n/a)`),
    Size_bin = if_else(Size <=2, "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when(.x == 0 ~ "T0",
                  .x <= 2 ~ "T1",
                  .x <= 5 ~ "T2",
                  .x > 5 ~ "T3"))
    ),

    Er = if_else(`Er_(0=negative;_1=positive;_9_=_n/a)` == 0, "neg", "pos"),
    Pr = if_else(`Pgr__(0=negative;_1=positive;_9_=_n/a)` == 0, "neg", "pos"),
    Her2 = if_else(`Her2_(0=negative;_1=positive;_9_=_n/a)` == 0, "neg", "pos"),

    Response_pathological = if_else(Pcr == "No", "npCR", "pCR"),
    Response = Response_pathological,



    Histology = purrr::map_chr(
      `Histology_(1=necrosis;_2=ductal_carcinoma;_3=lobular;_4=mixed_ductal/lobular_carcinoma;_5=other__6=no_invasive_tumor_present;_9_=_n/a_)`,
      ~(case_when(.x == 1 ~ "necrosis",
                  .x == 2 ~ "ductal_carcinoma",
                  .x == 3 ~ "lobular",
                  .x == 4 ~ "mixed_ductal/lobular_carcinoma",
                  .x == 5 ~ "other",
                  .x == 6 ~ "no_invasive_tumor_present"))
    ),

    # Abstract
    # "Four cycles of doxorubicin (75 mg/m2) or docetaxel (100 mg/m2) were compared as presurgical chemotherapy for breast cancer"

    Regimen = "neoadj",
    Arm_neoadj = if_else(Neoadjuvant_chemotherapy == "Taxane", "Docetaxel", "Doxorubicin")

  )


table(xx1$Arm_neoadj)
# Docetaxel Doxorubicin
# 12          23

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE21997_GPL1390_cleaned.tsv"
)





# GSE21997_GPL5325 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE21997_GPL5325.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 32


table(xx$Study)
# Madrid
# 28
table(xx$`Histologic_grade_(1=grade_i_(low);__2=_grade_ii_(intermediate);_3=_grade_iii_(high);_9_=_n/a)`)
# 1  2  3
# 1 12 15
table(xx$`Histology_(1=necrosis;_2=ductal_carcinoma;_3=lobular;_4=mixed_ductal/lobular_carcinoma;_5=other__6=no_invasive_tumor_present;_9_=_n/a_)`)
# 2  3  5
# 21  5  2
table(xx$`Tumor_size_(-1_=_n/a;__1_=_<=_2cm;_2_=_2-5_cm;_3_=_>_5cm)`)
# 2  3
# 16 12
table(xx$`Er_(0=negative;_1=positive;_9_=_n/a)`)
# 0  1
# 14 14
table(xx$`Pgr__(0=negative;_1=positive;_9_=_n/a)`)
# 0  1
# 12 16
table(xx$`Her2_(0=negative;_1=positive;_9_=_n/a)`)
# 0  1
# 21  7
table(xx$Neoadjuvant_chemotherapy)
# Doxorubicin      Taxane
# 14          14
table(xx$Pcr)
# No Yes
# 27   1
table(xx$Residual_cancer_burden_index_class)
# 0   I  II III
# 1   4   8  15


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = str_c("G", `Histologic_grade_(1=grade_i_(low);__2=_grade_ii_(intermediate);_3=_grade_iii_(high);_9_=_n/a)`),
    Size = if_else(`Tumor_size_(pre-chemotherapy),_cm_(-1_=_n/a)` == -1,
                   NA_real_, `Tumor_size_(pre-chemotherapy),_cm_(-1_=_n/a)`),
    Size_bin = if_else(Size <=2, "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when(.x == 0 ~ "T0",
                  .x <= 2 ~ "T1",
                  .x <= 5 ~ "T2",
                  .x > 5 ~ "T3"))
    ),

    Er = if_else(`Er_(0=negative;_1=positive;_9_=_n/a)` == 0, "neg", "pos"),
    Pr = if_else(`Pgr__(0=negative;_1=positive;_9_=_n/a)` == 0, "neg", "pos"),
    Her2 = if_else(`Her2_(0=negative;_1=positive;_9_=_n/a)` == 0, "neg", "pos"),

    Response_pathological = if_else(Pcr == "No", "npCR", "pCR"),
    Response = Response_pathological,

    Histology = purrr::map_chr(
      `Histology_(1=necrosis;_2=ductal_carcinoma;_3=lobular;_4=mixed_ductal/lobular_carcinoma;_5=other__6=no_invasive_tumor_present;_9_=_n/a_)`,
      ~(case_when(.x == 1 ~ "necrosis",
                  .x == 2 ~ "ductal_carcinoma",
                  .x == 3 ~ "lobular",
                  .x == 4 ~ "mixed_ductal/lobular_carcinoma",
                  .x == 5 ~ "other",
                  .x == 6 ~ "no_invasive_tumor_present"))
    ),

    # Abstract
    # "Four cycles of doxorubicin (75 mg/m2) or docetaxel (100 mg/m2) were compared as presurgical chemotherapy for breast cancer"

    Regimen = "neoadj",
    Arm_neoadj = if_else(Neoadjuvant_chemotherapy == "Taxane", "Docetaxel", "Doxorubicin")
  )

table(xx1$Arm_neoadj)
# Docetaxel Doxorubicin
# 14          14


write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE21997_GPL5325_cleaned.tsv"
)



# GSE21997_GPL7504 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE21997_GPL7504.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 31


table(xx$Study)
# Madrid
# 31
table(xx$`Histologic_grade_(1=grade_i_(low);__2=_grade_ii_(intermediate);_3=_grade_iii_(high);_9_=_n/a)`)
# 2  3
# 16 15
table(xx$`Histology_(1=necrosis;_2=ductal_carcinoma;_3=lobular;_4=mixed_ductal/lobular_carcinoma;_5=other__6=no_invasive_tumor_present;_9_=_n/a_)`)
# 2  3  4
# 27  3  1
table(xx$`Tumor_size_(-1_=_n/a;__1_=_<=_2cm;_2_=_2-5_cm;_3_=_>_5cm)`)
# 2  3
# 14 17
table(xx$`Er_(0=negative;_1=positive;_9_=_n/a)`)
# 0  1
# 10 21
table(xx$`Pgr__(0=negative;_1=positive;_9_=_n/a)`)
# 0  1
# 9 22
table(xx$`Her2_(0=negative;_1=positive;_9_=_n/a)`)
# 0  1
# 21 10
table(xx$Neoadjuvant_chemotherapy)
# Doxorubicin      Taxane
# 17          14
table(xx$Pcr)
# No Yes
# 26   5
table(xx$Residual_cancer_burden_index_class)
# 0  II III
# 5  15  11

xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = str_c("G", `Histologic_grade_(1=grade_i_(low);__2=_grade_ii_(intermediate);_3=_grade_iii_(high);_9_=_n/a)`),
    Size = if_else(`Tumor_size_(pre-chemotherapy),_cm_(-1_=_n/a)` == -1,
                   NA_real_, `Tumor_size_(pre-chemotherapy),_cm_(-1_=_n/a)`),
    Size_bin = if_else(Size <=2, "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when(.x == 0 ~ "T0",
                  .x <= 2 ~ "T1",
                  .x <= 5 ~ "T2",
                  .x > 5 ~ "T3"))
    ),

    Er = if_else(`Er_(0=negative;_1=positive;_9_=_n/a)` == 0, "neg", "pos"),
    Pr = if_else(`Pgr__(0=negative;_1=positive;_9_=_n/a)` == 0, "neg", "pos"),
    Her2 = if_else(`Her2_(0=negative;_1=positive;_9_=_n/a)` == 0, "neg", "pos"),

    Response_pathological = if_else(Pcr == "No", "npCR", "pCR"),
    Response = Response_pathological,

    Histology = purrr::map_chr(
      `Histology_(1=necrosis;_2=ductal_carcinoma;_3=lobular;_4=mixed_ductal/lobular_carcinoma;_5=other__6=no_invasive_tumor_present;_9_=_n/a_)`,
      ~(case_when(.x == 1 ~ "necrosis",
                  .x == 2 ~ "ductal_carcinoma",
                  .x == 3 ~ "lobular",
                  .x == 4 ~ "mixed_ductal/lobular_carcinoma",
                  .x == 5 ~ "other",
                  .x == 6 ~ "no_invasive_tumor_present"))
    ),

    # Abstract
    # "Four cycles of doxorubicin (75 mg/m2) or docetaxel (100 mg/m2) were compared as presurgical chemotherapy for breast cancer"

    Regimen = "neoadj",
    Arm_neoadj = if_else(Neoadjuvant_chemotherapy == "Taxane", "Docetaxel", "Doxorubicin")
  )

table(xx1$Arm_neoadj)
# Docetaxel Doxorubicin
# 14          17

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE21997_GPL7504_cleaned.tsv"
)




# GSE42822 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE42822.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 91

table(xx$`Pcr_(1)_vs_rd_(0)`)
# 0  1
# 54 37
table(xx$Prechemo_t)
# 1  2  3  4
# 1 34 47  9
table(xx$Prechemo_n)
# N0 N1 N2 N3
# 29 45  9  5
table(xx$`Prechemo_n.status_(1`)
# positive, 0: negative): 0  positive, 0: negative): 1 positive, 0: negative): NA
# 29                         59                          3
table(xx$Bmn.grade)
# G1/2   G3
# 23   53
table(xx$`Er.status_(1`)
# positive, 0: negative): 0  positive, 0: negative): 1 positive, 0: negative): NA
# 52                         38                          1
table(xx$`Pr.status_(1`)
# positive, 0: negative): 0  positive, 0: negative): 1 positive, 0: negative): NA
# 50                         40                          1
table(xx$`Her2.status_(1`)
# positive, 0: negative): 0  positive, 0: negative): 1 positive, 0: negative): NA
# 54                         34                          3
table(xx$Treatment.type)
# FEC/TX FEC/TX+H
# 66       25



xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = if_else(Bmn.grade == "G1/2","G1_2", "G3"),
    Node_cat = Prechemo_n,
    Node_bin = if_else(Node_cat == "N0", "neg", "pos"),
    Size_cat = str_c("T", Prechemo_t),
    Size_bin = if_else(Size_cat == "T1", "small", "large"),


    Er = `Er.status_(1` %>% str_replace("positive, 0: negative\\): ", "") %>% as.integer(),
    Pr = `Pr.status_(1` %>% str_replace("positive, 0: negative\\): ", "") %>% as.integer(),
    Her2 = `Her2.status_(1` %>% str_replace("positive, 0: negative\\): ", "") %>% as.integer(),
    Er = if_else(Er == 0, "neg", "pos"),
    Pr = if_else(Pr == 0, "neg", "pos"),
    Her2 = if_else(Her2 == 0, "neg", "pos"),

    Response_pathological = if_else(`Pcr_(1)_vs_rd_(0)` == 0, "npCR", "pCR"),
    Response = Response_pathological,

    # Abstract
    # "Tumor samples were obtained from patients with stage II-III breast cancer before starting neoadjuvant chemotherapy with four cycles of 5-fluorouracil/epirubicin/cyclophosphamide (FEC) followed by four cycles of docetaxel/capecitabine (TX) on US Oncology clinical trial 02-103. Most patients with HER-2-positive cancer also received trastuzumab (H). "

    Regimen = "neoadj",
    Arm_neoadj = if_else(Treatment.type == "FEC/TX",
                  "Docetaxel+Capecitabine+FEC",
                  "Docetaxel+Capecitabine+Trastuzumab+FEC")
  )

table(xx1$Arm_neoadj)
# Docetaxel+Capecitabine+FEC Docetaxel+Capecitabine+Trastuzumab+FEC
# 66                                     25

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE42822_cleaned.tsv"
)




# GSE66305 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE66305.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 88


table(xx$Disease_state)
# HER2-positive breast cancer
# 88
table(xx$Tissue)
# breast cancer
# 88
table(xx$`Cher-lob_arm`)
# A  B  C
# 23 31 34
table(xx$Arm_description)
# chemotherapy+lapatinib           chemotherapy+trastuzumab chemotherapy+trastuzumab+lapatinib
# 31                                 23                                 34
table(xx$`Pcr_(1=yes)`)
# 0  1
# 61 27


xx1 <- xx %>%
  dplyr::mutate(
    Her2 = "pos",

    Response_pathological = if_else(`Pcr_(1=yes)` == 0, "npCR", "pCR"),
    Response = Response_pathological,

    # "CHER-LOB is a phase II randomized multicenter trial in which 121 patients with primary HER2-positive breast cancer were randomized to receive preoperative chemotherapy with weekly paclitaxel for 12 weeks followed by 4 weekly courses over 3 weeks of the FEC regimen (fluorouracil, epirubicin, and cyclophosphamide) plus either trastuzumab (arm A), lapatinib (arm B), or the combination of trastuzumab and lapatinib (arm C)."
    Regimen = "neoadj",
    Arm_neoadj = purrr::map_chr(
      Arm_description,
      ~(case_when(.x == "chemotherapy+lapatinib" ~ "Paclitaxel+Lapatinib+FEC",
                  .x == "chemotherapy+trastuzumab" ~ "Paclitaxel+Trastuzumab+FEC",
                  .x == "chemotherapy+trastuzumab+lapatinib" ~ "Paclitaxel+Lapatinib+Trastuzumab+FEC"))
    )
  )


table(xx1$Arm_neoadj)
# Paclitaxel+Lapatinib+FEC Paclitaxel+Lapatinib+Trastuzumab+FEC
# 31                                   34
# Paclitaxel+Trastuzumab+FEC
# 23


write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE66305_cleaned.tsv"
)





# GSE66999 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE66999.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 76


table(xx$Tissue)
# human breast tumor core biopsy
# 76
table(xx$Treatment)
# pretreatment
# 76
table(xx$Schedule)
# Schedule A Schedule B
# 41         35
table(xx$Clinical_response)
# Complete Response (CR)    Partial Response (PR) Progressive Disease (PD)
# 13                       57                        1
# Stable Disease (SD)
# 5
table(xx$`Estrogen-receptor`)
# negative positive
# 35       38
table(xx$Progesterone_receptor)
# negative positive
# 46       27
table(xx$Her2)
# negative positive
# 55       18
table(xx$Topo2)
# negative positive
# 14       58
table(xx$Experimental_batches)
# A  B  C  D  E  F  G
# 19  7  4 18 11 16  1


xx1 <- xx %>%
  dplyr::mutate(

    Er = if_else(`Estrogen-receptor` == "negative", "neg", "pos"),
    Pr = if_else(Progesterone_receptor == "negative", "neg", "pos"),
    Her2 = if_else(Her2 == "negative", "neg", "pos"),

    Response_clinical = purrr::map_chr(
      Clinical_response,
      ~(case_when(.x == "Complete Response (CR)" ~ "CR",
                  .x == "Partial Response (PR)" ~ "PR",
                  .x == "Stable Disease (SD)" ~ "SD",
                  .x == "Progressive Disease (PD)" ~ "PD"))
    ),
    Response = if_else(Response_clinical == "CR", "pCR", "npCR"),
    Topo2 = if_else(Topo2 == "negative", "neg", "pos"),

    # "Various doses of epirubicin and docetaxel were administered to patients in either a standard q3 weekly (Schedule A) or dose dense q2 weekly (Schedule B) regimen. Doses for Schedule A were 75 mg/m2 IV of docetaxel and 75, 90, 105, or 120 mg/m2 IV of epirubicin (with 6 mg of pegfilgrastim per cycle on day 2 to prevent neutropenia). Doses for both docetaxel and epirubicin in Schedule B were 50, 60, and 70 mg/m2 IV (with identical pegfilgrastim support). For each schedule, phase I was dose finding for phase II. Patients were allocated to the various phases and dosing regimens of the trial without randomization. None of the patients received trastuzumab in the initial years of the study and HER2+ patients were not enrolled on study, once trastuzumab funding became available. "

    Regimen = "neoadj",
    Arm_neoadj = "Docetaxel+Epirubicin"
  )


table(xx1$Arm_neoadj)
# Docetaxel+Epirubicin
# 76

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE66999_cleaned.tsv"
)






# GSE28844 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE28844.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 61

table(xx$Gender)
# female
# 61
table(xx$Tissue)
# breast tumor
# 61
table(xx$Trucut_biopsy)
# 3.5 cm x 1.75 mm
# 32
table(xx$Tumor_biopsy_acquisition)
# Such biopsy was taken at diagnosis
# 32
# Such biopsy was taken from surgery specimen and selected by the pathologist in charge within 30 min. after tumor removal
# 29
table(xx$Treatment)
# A                      B                      C none (Prechemotherapy)
# 17                      5                      7                     18
# none(Prechemotherapy)
# 14
table(xx$`Pathologic_response_to_chemotherapy_(miller_&_payne_grade)`)
# 1  2  3  4  5
# 4 24 20  7  6

# Arm
# Treatment A: Epirubicin 90 mg/m2-Cyclophosphamide 600 mg/m2, 3 cycles bi-weekly and Taxol 150 mg/m2-Gemcitabine 2500 mg/m2, 6 cycles bi-weekly ± weekly Herceptin 4 mg/Kg during the first week, 2 mg/Kg for the remaining 11 cycles; Treatment B: Doxorubicin 60 mg/m2-Pemetrexed 500 mg/m2, 4 cycles tri-weekly and Taxotere 100 mg/m2, 4 cycles tri-weekly; Treatment C: Doxorubicin 60 mg/m2-Cyclophosphamide 600 mg/m2, 4 cycles tri-weekly and Taxotere 100 mg/m2, 4 cycles tri-weekly

# Treatment A: Epirubicin+Cyclophosphamide+Taxol+Gemcitabine±Herceptin
# Treatment B: Doxorubicin+Pemetrexed+Taxotere
# Treatment C: Doxorubicin+Cyclophosphamide+Taxotere
# Anthracyclin+Taxane+-Antimetabolite+-Alkalyting_agent+-Trastuzumab


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,

    Response_clinical = str_c("MPG", # Miller & Payne Grade
      `Pathologic_response_to_chemotherapy_(miller_&_payne_grade)`),
    Response = if_else(Response_clinical == "MPG4" | Response_clinical == "MPG5",
                       "pCR", "npCR"),
    Timepoint = if_else(Tumor_biopsy_acquisition == "Such biopsy was taken at diagnosis",
                        "Pre_treatment", "Post_treatment"),


    # GEO description:
    # "After informed consent, patients with a histologically confirmed diagnosis of breast cancer and scheduled chemotherapy treatment based on Anthracyclines and Taxanes (Treatment A: Epirubicin 90 mg/m2-Cyclophosphamide 600 mg/m2, 3 cycles bi-weekly and Paclitaxel 150 mg/m2-Gemcitabine 2500 mg/m2, 6 cycles bi-weekly ± weekly Herceptin 4 mg/Kg during the first week, 2 mg/Kg for the remaining 11 cycles; Treatment B: Doxorubicin 60 mg/m2-Pemetrexed 500 mg/m2, 4 cycles tri-weekly and Docetaxel 100 mg/m2, 4 cycles tri-weekly; Treatment C: Doxorubicin 60 mg/m2-Cyclophosphamide 600 mg/m2, 4 cycles tri-weekly and Docetaxel 100 mg/m2, 4 cycles tri-weekly ) were recruited for this study. Pre-chemotherapy and post-chemotherapy biopsies were examined by a pathologist who determined the Miller & Payne grade for each patient. Matching pairs of pre-chemotherapy and post-chemotherpy samples were divided into 3 groups according to Miller & Payne grade: group of bad response (Miller & Payne grades 1 and 2), group of mid response (Miller & Payne grade 3) and group of good response (Miller & Payne grades 4 and 5)."

    # Taxanes: Docetaxel|Paclitaxel
    # Anthracyclines: Doxorubicin|Epirubicin
    # Anti-metabolite: Gemcitabine|Pemetrexed
    # Alkalytin_agent: Cyclophosphamide
    # Heceptin
    Regimen = "neoadj",
    Arm_neoadj = "Anthracyclin+Taxane[±Antimetabolite±Alkalyting_agent±Trastuzumab]",
    Arm_neoadj = if_else(Timepoint == "Post_treatment", NA_character_, Arm_neoadj) # filtering
  )


table(xx1$Arm_neoadj)
# Anthracyclin+Taxane[±Antimetabolite±Alkalyting_agent±Trastuzumab]
# 32

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE28844_cleaned.tsv"
)




# GSE23988 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE23988.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 61

table(xx$Pcr.v.rd)
# pCR  RD
# 20  41
table(xx$Er_positive_vs_negative)
# ERneg ERpos
# 29    32
table(xx$Grade)
# 1  2  3
# 1 19 37
table(xx$Prechemo_t_stage)
# 1  2  3
# 1 20 40
table(xx$Prechemo_nodal_status)
# 0  1  2  3
# 21 32  5  3


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = str_c("G", Grade),
    Node_bin = if_else(Prechemo_nodal_status == 0, "neg", "pos"),
    Node_cat = str_c("N", Prechemo_nodal_status),
    Size = Prechemo_tumor_size,
    Size_bin = if_else(Prechemo_t_stage <= 1, "small", "large"),
    Size_cat = str_c("T", Prechemo_t_stage),

    Er = str_replace(Er_positive_vs_negative, "ER", ""),
    Her2 = "neg", # from text

    Response_pathological = if_else(Pcr.v.rd == "RD", "npCR", "pCR"),
    Response = Response_pathological,

    Regimen = "neoadj",
    Arm_neoadj = "Docetaxel+Capecitabine+FEC"

  )


table(xx1$Arm_neoadj)

# Docetaxel+Capecitabine+FEC
# 61


write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE23988_cleaned.tsv"
)





# GSE18728 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE18728.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 61

table(xx$Tissue_type)
# breast tumor biopsies
# 61
table(xx$Time_point)
# bl c2 C2 or
# 21 17  1 22
table(xx$Response_category)
# NR  R
# 38 23
table(xx$Er_original)
# focally pos         neg         pos
# 2          29          30
table(xx$Pr)
# focally pos         neg         pos
# 3          36          22
table(xx$Her2_summary)
# neg pos
# 44  17
table(xx$Upin)
# 1  3  4  5  6  7  8 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
# 1  2  1  1  3  2  3  1  1  3  3  2  3  2  3  3  3  3  1  2  2  2  2  3  2  2  3  2



xx1 <- xx %>%
  dplyr::mutate(

    Er = if_else(Er_original == "neg", "neg", "pos"),
    Pr = if_else(Pr == "neg", "neg", "pos"),
    Her2 = Her2_summary,

    Response_clinical = if_else(Response_category == "R", "Response", "No_response"),
    Response = if_else(Response_clinical == "No_response", "npCR", "pCR"),


    Timepoint = purrr::map_chr(
      Time_point,
      ~(case_when(.x == "bl" ~ "Pre_treatment",
                  .x == "c2" ~ "Mid_treatment",
                  .x == "C2" ~ "Mid_treatment",
                  .x == "or" ~ "Post_treatment"))
    ),


    # "All patients received four cycles of docetaxel and capecitabine administered every 21 days. Patients were initially treated with docetaxel (75 mg/m2 i.v.) on day 1 and capecitabine (1,000 mg/m2 p.o.) twice daily on days 2–15 every 21 days for four cycles. Due to excessive toxicity in the first 10 patients treated on protocol, both agents were dose reduced (docetaxel to 60 mg/m2 and capecitabine to 937.5 mg/m2 twice daily). After completing four cycles of docetaxel/capecitabine (TX), all patients received four cycles of adriamycin (60 mg/m2) and cyclophosphamide (600 mg/m2) on day 1 and every 21 days; six patients with a poor response to the initial chemotherapy received adriamycin and cyclophosphamide (AC) prior to surgery, the remainder had this additional therapy post-operatively. "
    #
    # "We classified patients as either responders or non-responders based on change in tumor size by clinical exam and pathologic response (see Table 1). For patients treated pre-operatively with AC, response was assessed at the completion of TX."

    Regimen = "neoadj",
    Arm_neoadj = "Docetaxel+Capecitabine",
    Arm_neoadj = if_else(Timepoint == "Pre_treatment", Arm_neoadj, NA_character_)
  )


table(xx1$Arm_neoadj)
# Docetaxel+Capecitabine
# 21

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE18728_cleaned.tsv"
)


# GSE21974 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE21974.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 61

table(xx$`Before/after_chemotherapy`)
# after before
# 25     32
table(xx$Response)
# NR  R
# 39 18
table(xx$Sinn)
# 1  2  3  4
# 25 14 10  8
table(xx$Ypt)
# 0    1   1a   1b   1c 1mic    2   is
# 12    2    9    9   13    2    6    4
table(xx$Ypn)
# 0  1 1a  2 2a  3  x
# 41  1  7  2  3  2  1
table(xx$Grade)
# 1  2  3
# 2 27 28
table(xx$`Er(%)`)
# 0 10 15 40 50 60 70 75 80 90
# 24  2  2  1  8  5  2  2  5  6
table(xx$`Pr(%)`)
# 0  3 10 20 25 30 40 50 60 80 90 95
# 24  2  5  1  2  2  2  4  4  4  4  3
table(xx$Her2neu_score)
# 0  1  2  3
# 18 24  2 13


xx1 <- xx %>%
  dplyr::mutate(
    Grade = str_c("G", Grade),
    Size = Initial_cm,
    Size_bin = if_else(Size <=2, "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when(.x == 0 ~ "T0",
                  .x <= 2 ~ "T1",
                  .x <= 5 ~ "T2",
                  .x > 5 ~ "T3"))
    ),

    Er_score = `Er(%)`,
    Pr_score = `Pr(%)`,
    Her2_score = Her2neu_score,

    Er = if_else(Er_score <=10, "neg", "pos"),
    Pr = if_else(Pr_score <=10, "neg", "pos"),
    Her2 = if_else(Her2_score <= 1, "neg", "pos"),
    #Her2 cutoff:  https://www.breastcancer.org/symptoms/diagnosis/her2

    Response_pathological = if_else(Response == "R", "pCR", "npCR"),
    # imputed from Table1; Response and pCR numbers were matching
    Response = Response_pathological,

    Timepoint = if_else(`Before/after_chemotherapy` == "before",
                        "Pre_treatment", "Post_treatment"),


    # After establishing the diagnosis of invasive breast cancer, all women underwent sequential NAC with 4 cycles of epirubicine 90 mg/m2 and cyclophosphamide 600 mg/m2 every 3  weeks, followed by 4 cycles of docetaxel 100 mg/m2 every 3  weeks. At the time of diagnosis and at 4 time points during NAC and after completion of therapy, lesions were ultrasono-graphically measured to assess response. After completion of 4  cycles of epirubicine/cyclophosphamide, a second two-pass high speed core biopsy of the primary tumor was performed and tissue was stored for analysis. All patients underwent surgery after completion of chemotherapy.

    Regimen = "neoadj",
    Arm_neoadj = "Docetaxel+EC",
    Arm_neoadj = if_else(Timepoint == "Pre_treatment", Arm_neoadj, NA_character_)
  )


table(xx1$Arm_neoadj)
# Docetaxel+EC
# 32

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE21974_cleaned.tsv"
)


# GSE143222 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE143222.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 61

table(xx$Tissue)
# Triple Negative Breast cancer
# 55
table(xx$Status_after_5_years)
# Alive  Dead
# 34    21
table(xx$Relapse)
# non reccurence     reccurence
# 35             20
table(xx$Pcr)
# complete remission non complete remission
# 18                     37
table(xx$`Pre/post`)
# Post  Pre
# 14   41


xx1 <- xx %>%
  dplyr::mutate(
    Er = "neg",
    Pr = "neg",
    Her2 = "neg",

    Event_os = if_else(Status_after_5_years == "Dead", 1, 0),
    Event_dfs = if_else(Relapse == "reccurence", 1, 0),
    Response_pathological = if_else(Pcr == "complete remission", "pCR", "npCR"),
    Response = Response_pathological,


    Timepoint = if_else(`Pre/post` == "Pre",
                        "Pre_treatment", "Post_treatment"),

    # All patients received anthracy-cline and taxane-based regimens that included four cyclesof  60 mg/m2adriamycin  and  600 mg/m2cyclophos-phamide followed by four cycles of 75 mg/m2docetaxel.
    # Ref: PMID: 26006068
    # "Prognostic and predictive value of NanoString-based immune-related gene signatures in a neoadjuvant setting of triple-negativebreast cancer: relationship to tumor-infiltrating lymphocytes"
    Regimen = "neoadj",
    Arm_neoadj = "Docetaxel+AC",
    Arm_neoadj = if_else(Timepoint == "Pre_treatment", Arm_neoadj, NA_character_)
  )


table(xx1$Arm_neoadj)
# Docetaxel+AC
# 41

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE143222_cleaned.tsv"
)




# GSE16391 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE16391.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 55

table(xx$`Case_control_0-_case,_1_-_control`)
# 0  1
# 10 38
table(xx$Post_menopausal_status)
# 1  2
# 48  7
table(xx$Er_pgr)
# 1  2
# 41 14
table(xx$Adj_neoadj_chemotherapy_received)
# 0  1
# 35 20
table(xx$Grade)
# 1  2  3
# 2 35 18
table(xx$Local_therapy)
# 1  2  3  4
# 34  1 14  6
table(xx$Node)
# 0  1
# 22 33
table(xx$Treatment)
# 0  1
# 32 23
table(xx$Her2_status)
# 0  1
# 42  3
table(xx$Size)
# 1  2
# 20 35
table(xx$Tissue)
# primary breast tumor
# 55



xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = str_c("G", Grade),
    Node_bin = if_else(Node == 0, "neg", "pos"),
    Size_bin = if_else(Size == 1, "small", "large"),


    Er = "pos",
    Pr = if_else(Er_pgr == 1, "pos", "neg"),
    Her2 = if_else(Her2_status == 1, "pos", "neg"),

    Event_dfs = if_else(`Case_control_0-_case,_1_-_control` == 0, 1, 0), # case =0
    Time_dfs = T_rfs_months / 12,

    Therapy_chemo = if_else(Adj_neoadj_chemotherapy_received == 1, "yes", "no"),
    Therapy_radio = purrr::map_chr(
      Local_therapy,
      ~(case_when(.x == 1 ~ "yes",
                  .x == 2 ~ "no",
                  .x == 3 ~ "yes",
                  .x == 4 ~ "no"))
    ),
    Therapy_surgery = "yes",
    Therapy_surgery_type = purrr::map_chr(
      Local_therapy,
      ~(case_when(.x == 1 ~ "Breast_conserving_surgery",
                  .x == 2 ~ "Breast_conserving_surgery",
                  .x == 3 ~ "Mastectomy",
                  .x == 4 ~ "Mastectomy"))
    ),
    Therapy_hormone = "yes",
    Therapy_hormone_name = if_else(Treatment == 0, "Letrozole", "Tamoxifen"),


    # The following variables were imputed from Table 2 of associated publication.
    # Size_bin,Er,Pr,Event_dfs,Therapy_radio,Therapy_surgery,Therapy_surgery_type,
    # Therapy_hormone, Therapy_hormone_name, Menopause

    Menopause = if_else(Post_menopausal_status == 2, "pre", "post"),
    # Menopause achived after therapy is set to "pre"
    # Consider only pre-treatment scenario.

    Regimen = "adj",
    Arm_adj = str_c(
      if_else(Therapy_radio == "yes", "Radio+",""),
      if_else(Therapy_chemo == "yes", "Chemo+",""),
      Therapy_hormone_name
    )

  )


table(xx1$Arm_adj)
# Chemo+Tamoxifen             Letrozole Radio+Chemo+Letrozole Radio+Chemo+Tamoxifen
# 1                     3                    10                     9
# Radio+Letrozole       Radio+Tamoxifen             Tamoxifen
# 19                    10                     3

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE16391_cleaned.tsv"
)




# GSE75678 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE75678.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 54


# Note !!!!!!!!!!!!!!!!!!!!!!
# This dataset is discarded as regimen is not explisitly mentioned in GEO descriptions or in published text. Published. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5469719/
# !!!!!!!!!!!!!!!!!!!!!!!!!!!

table(xx$Histological_grade)
# 1   2   3 N/A
# 2  23  27   2
table(xx$Estrogen_receptor_immunohistochemestry)
# Neg Pos
# 25  29
table(xx$Progesteron_receptor_immunohistochemestry)
# Neg Pos POs
# 28  25   1
table(xx$Her2_immunohistochemestry)
# N/A Neg Pos
# 3  37  14
table(xx$Her2_fish)
# N/A Neg Pos
# 51   2   1
table(xx$P63_ihc)
# Absent    N/A    Neg    Pos
# 19     17     17      1
table(xx$`T`)
# 1  2  3  4
# 5 16 22 11
table(xx$N)
# 0  1  2  3
# 8 26 19  1
table(xx$Response_to_treatment)
# Complete    Partial       Poor Progresive
# 6         15         22         11
table(xx$Chemotherapy1)
# Adriamicin Doxorubicin  Epirubicin         N/A  Paclitaxel   Tamoxifen       Taxol Trastuzumab
# 14           1           2           2          29           1           2           3
table(xx$Chemotherapy2)
# Adriamicin Cyclophosphamide      Doxorubicin      Epirrubicin              N/A            Taxol
# 2               21               11                1               15                2
# Trastazumab
# 2
table(xx$Chemotherapy3)
# Adramicin       Adriamicin Cyclophosphamide              N/A       Paclitaxel            Taxol
# 1                5               13               26                4                2
# Trastuzumab
# 3
table(xx$Chemotherapy4)
# N/A       Taxol Trastuzumab
# 50           1           3
table(xx$Radiotherapy)
# N/A  No Yes
# 3  24  27
table(xx$Place_to_metastasis)
# Liver  Lung   N/A
# 1     3    50

# Study published Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5469719/
# Title: "A New Gene Expression Signature for Triple-Negative Breast Cancer Using Frozen Fresh Tissue before Neoadjuvant Chemotherapy"



xx1 <- xx %>%
  dplyr::mutate(
    Age = Age_at_diagnosis,
    Grade = str_c("G", Histological_grade %>% as.integer()),
    Node_cat = str_c("N", N %>% as.integer()),
    Node_bin = if_else(Node_cat == "N0", "neg", "pos"),
    Size = Initial_highest_diameter_in_mm_by_image / 10, # cm
    Size_bin = if_else(Size <=2, "small", "large"),
    Size_cat = purrr::map_chr(
      Size,
      ~(case_when(.x == 0 ~ "T0",
                  .x <= 2 ~ "T1",
                  .x <= 5 ~ "T2",
                  .x > 5 ~ "T3"))
    ),

    Er = if_else(Estrogen_receptor_immunohistochemestry == "Neg", "neg", "pos"),
    Pr = if_else(Progesteron_receptor_immunohistochemestry == "Neg", "neg", "pos"),
    Her2 = if_else(Her2_immunohistochemestry == "N/A", # if NA the HER2.fish
                   Her2_fish, Her2_immunohistochemestry),
    Her2 = str_replace(Her2, "N/A", NA_character_) %>% str_to_lower(),


    Response_clinical = purrr::map_chr(
      Response_to_treatment ,
      ~(case_when(.x == "Complete" ~ "CR",
                  .x == "Partial" ~ "PR",
                  .x == "Poor" ~ "SD",
                  .x == "Progresive" ~ "PR",))
    ),
    Response = if_else(Response_clinical == "CR" | Response_clinical == "PR",
                       "pCR", "npCR"),

    Therapy_radio = if_else(Radiotherapy == "N/A", NA_character_, Radiotherapy) %>%
      str_to_lower(),

    Arm = str_c(Chemotherapy1, Chemotherapy2, Chemotherapy3, Chemotherapy4, sep = " ") %>%
      str_replace_all("N/A", "") %>%
      str_replace_all("Adramicin", "Doxorubicin") %>%
      str_replace_all("Adriamicin", "Doxorubicin") %>%
      str_replace_all("Epirrubicin", "Epirubicin") %>%
      str_replace_all("Taxol", "Paclitaxel") %>%
      str_replace_all("Trastazumab", "Trastuzumab") %>%
      str_trim(),
    Arm = str_c(
      if_else(str_detect(Arm, "Paclitaxel"), "Paclitaxel ", ""),
      if_else(str_detect(Arm, "Doxorubicin"), "Doxorubicin ", ""),
      if_else(str_detect(Arm, "Epirubicin"), "Epirubicin ", ""),
      if_else(str_detect(Arm, "Cyclophosphamide"), "Cyclophosphamide ", ""),
      if_else(str_detect(Arm, "Trastuzumab"), "Trastuzumab ", ""),
      if_else(str_detect(Arm, "Tamoxifen"), "Tamoxifen ", "")
    ) %>%
      str_trim() %>%
      str_replace_all(" ", "+"),
    Arm = if_else(Arm == "", NA_character_, Arm)
  )

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE75678_cleaned.tsv"
)


# To aid in Arm formatting !!!!!!!!!!
yy <- names(table(xx1$Arm %>% str_replace_all("N/A", "") %>% str_trim()))
yy <- purrr::map(yy,~(str_split(string = .x, pattern = " ")))
names(table(unlist(yy)))
# [1] ""                 "Adramicin"        "Adriamicin"       "Cyclophosphamide"
# [5] "Doxorubicin"      "Epirrubicin"      "Epirubicin"       "Paclitaxel"
# [9] "Tamoxifen"        "Taxol"            "Trastazumab"      "Trastuzumab"
# Cleaned
# [1] ""                 "Cyclophosphamide" "Doxorubicin"      "Epirubicin"
# [5] "Paclitaxel"       "Tamoxifen"        "Trastuzumab"



# GSE55348 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE55348.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 53

# Note !!!!!!!!!!!!!!!!!!!!!!
# This dataset is discarded as the samples set includes a mixture of neoadj/adj treatment (GHEA cohort), and the regimen is not explicitly mentioned.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!

table(xx$`T`)
# T1 T2 T3 T4
# 21 19  5  8
table(xx$N)
# neg pos
# 7  46
table(xx$Grade)
# II III
# 14  39
table(xx$Er)
# neg pos
# 21  32
table(xx$Pgr)
# neg pos
# 29  24
table(xx$Her2)
# neg pos
# 2  51
table(xx$Event)
# 0  1
# 30 23

summary(xx$Time_to_relapse/325)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.06769 1.53846 2.03077 2.38752 3.33539 5.01539


xx1 <- xx %>%
  dplyr::mutate(
    Age = Age,
    Grade = if_else(Grade == "II", "G2", "G3"),
    Node_bin = N,
    Size_cat = `T`,
    Size_bin = if_else(Size_cat == "T1", "small", "large"),

    Er = Er,
    Pr = Pgr,
    Her2 = Her2,

    Event_dfs = Event,
    Time_dfs = Time_to_relapse / 325, # years

    Arm = "Paclitaxel+Doxorubicin+Trastuzumab+CMF+-Hormone"
  )

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE55348_cleaned.tsv"
)





# GSE8465_GPL1390 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE8465_GPL1390.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 2


# Note !!!!!!!!!!!!!!!!!!!!!!
# This dataset is discarded as only post_treatment profiles are present (imputed from text).
# !!!!!!!!!!!!!!!!!!!!!!!!!!!



table(xx$Sample_characteristics_ch1)
# Stratagene Human Universal Reference that contained 1/10 added MCF7 and ME16C RNAs
# 2


xx <- xx %>%
  dplyr::mutate(`Pretreatment_breast_tumor_from_neoadjuvant_trial;treatment_cohort` = str_c("treatment_cohort:",`Pretreatment_breast_tumor_from_neoadjuvant_trial;treatment_cohort`)) %>%
  dplyr::rename(Sample_characteristics_ch2 = "Pretreatment_breast_tumor_from_neoadjuvant_trial;treatment_cohort")



yy <- str_split_fixed(xx$Sample_characteristics_ch2, ";", 11) %>% as_tibble()

var_name <- get_var_name(yy)
var_qc <- get_var_qc(yy)

names(yy) <- var_name
yy <- purrr::pmap_dfr(

  list(x = yy,
       y = var_name,
       z = var_qc),

  function(x, y, z){

    if (z) {
      # str_split_fixed(x, str_c(y,": "), 2)[, 2]
      # If the column name "y" contains special character, eg "(",
      # the above code will not parse correctly

      str_split_fixed(x, ":", 2)[, 2] %>% str_trim()
    } else {
      x
    }

  }

)
glimpse(yy)

xx <- xx %>%
  dplyr::select(-Sample_characteristics_ch2) %>%
  bind_cols(yy)
names(xx) <- clean_names(names(xx))
glimpse(xx)

xx1 <- xx %>%
  dplyr::mutate(
    Grade = purrr::map_chr(
      Tumor_grade,
      ~(case_when(.x == "moderately differentiated" ~ "G2"))
    ),
    # Ref: https://www.mdanderson.org/patients-family/diagnosis-treatment/a-new-diagnosis/cancer-grade-vs--cancer-stage.html
    Node_cat = str_split_fixed(Staging_tnm, "N", 2)[, 2],
    Node_cat = str_c("N", Node_cat),
    Node_bin = if_else(Node_cat == "N0", "neg", "pos"),
    Size_cat = str_split_fixed(Staging_tnm, "N", 2)[, 1],
    Size_bin = if_else(Size_cat == "T0" | Size_cat == "T1", "small", "large"),

    Er = if_else(Er_status == "negative", "neg", "pos"),
    Pr = if_else(Pr_status == "negative", "neg", "pos"),
    Her2 = purrr::map_chr(
      Her2_status,
      ~(case_when(.x == "1+" ~ "neg"))
    ),

    Response_pathological =  purrr::map_chr(
      Pathological_complete_response,
      ~(case_when(.x == "no" ~ "npCR",
                  .x == "yes" ~ "pCR"))
    ),
    Response_clinical = if_else(Overall_best_response == "Not Assessed",
                                NA_character_, Overall_best_response) %>%
      str_replace("\\(As assessed by SWOG criteria, Green et al., Invest New Drugs 10:239-253 \\(1992\\)\\)",
                  "") %>%
      str_trim(),
    Response = Response_pathological,

    Arm = Treatment_cohort, # post therapy profiles
    Histology = Pathological_diagnosis
  )

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE8465_GPL1390_cleaned.tsv"
)





# GSE8465_GPL887 >>>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE8465_GPL887.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 44


# Note !!!!!!!!!!!!!!!!!!!!!!
# This dataset is discarded as only post_treatment profiles are present (imputed from text).
# !!!!!!!!!!!!!!!!!!!!!!!!!!!


table(xx$Sample_characteristics_ch1)
# Stratagene Human Universal Reference that contained 1/10 added MCF7 and ME16C RNAs
# 44


xx <- xx %>%
  dplyr::mutate(`Pretreatment_breast_tumor_from_neoadjuvant_trial;treatment_cohort` = str_c("treatment_cohort:",`Pretreatment_breast_tumor_from_neoadjuvant_trial;treatment_cohort`)) %>%
  dplyr::rename(Sample_characteristics_ch2 = "Pretreatment_breast_tumor_from_neoadjuvant_trial;treatment_cohort")



yy <- str_split_fixed(xx$Sample_characteristics_ch2, ";", 11) %>% as_tibble()

var_name <- get_var_name(yy)
var_qc <- get_var_qc(yy)

names(yy) <- var_name
yy <- purrr::pmap_dfr(

  list(x = yy,
       y = var_name,
       z = var_qc),

  function(x, y, z){

    if (z) {
      # str_split_fixed(x, str_c(y,": "), 2)[, 2]
      # If the column name "y" contains special character, eg "(",
      # the above code will not parse correctly

      str_split_fixed(x, ":", 2)[, 2] %>% str_trim()
    } else {
      x
    }

  }

)
glimpse(yy)

xx <- xx %>%
  dplyr::select(-Sample_characteristics_ch2) %>%
  bind_cols(yy)
names(xx) <- clean_names(names(xx))
glimpse(xx)

table(xx$Treatment_cohort)
# Doxorubicin                  Gemcitabine Gemcitabine plus Doxorubicin
# 20                           12                           12
table(xx$Overall_best_response)
# CR (As assessed by SWOG criteria, Green et al., Invest New Drugs 10:239-253 (1992))
# 15
# Not Assessed
# 3
# PD (As assessed by SWOG criteria, Green et al., Invest New Drugs 10:239-253 (1992))
# 5
# PR (As assessed by SWOG criteria, Green et al., Invest New Drugs 10:239-253 (1992))
# 21
table(xx$Completed_therapy)
# CR before last cycle                   no         Not Assessed                  yes
# 6                    9                    2                   27
table(xx$Pathological_complete_response)
# no Not Assessed          yes
# 23            8           13
table(xx$Staging_tnm)
# T2N0 T2N1 T3N0 T3N1 T4N0 T4N1 T4N2
# 2    3   16   11    6    4    2
table(xx$Pathological_diagnosis)
# Ductal Mucinous breast carcinoma                     Other
# 40                         1                         3
table(xx$Tumor_grade)
# moderately differentiated     poorly differentiated                   unknown       well differentiated
# 22                         4                        17                         1
table(xx$Er_status)
# negative Not Assessed     positive
# 21            1           22
table(xx$Pr_status)
# negative Not Assessed     positive
# 23            1           20
table(xx$Her2_status)
# 0           1+           2+           3+ Not Assessed Not detected
# 19            5            2           12            1            5



xx1 <- xx %>%
  dplyr::mutate(
    Grade = purrr::map_chr(
      Tumor_grade,
      ~(case_when(.x == "well differentiated" ~ "G1",
                  .x == "moderately differentiated" ~ "G2",
                  .x == "poorly differentiated" ~ "G3"))
    ),
    # Ref: https://www.mdanderson.org/patients-family/diagnosis-treatment/a-new-diagnosis/cancer-grade-vs--cancer-stage.html
    Node_cat = str_split_fixed(Staging_tnm, "N", 2)[, 2],
    Node_cat = str_c("N", Node_cat),
    Node_bin = if_else(Node_cat == "N0", "neg", "pos"),
    Size_cat = str_split_fixed(Staging_tnm, "N", 2)[, 1],
    Size_bin = if_else(Size_cat == "T0" | Size_cat == "T1", "small", "large"),

    Er = if_else(Er_status == "Not Assessed",
                 NA_character_,
                 if_else(Er_status == "negative", "neg", "pos")),
    Pr = if_else(Pr_status == "Not Assessed",
                 NA_character_,
                 if_else(Pr_status == "negative", "neg", "pos")),
    Her2 = purrr::map_chr(
      Her2_status,
      ~(case_when(.x == "0" ~ "neg",
                  .x == "1+" ~ "neg",
                  .x == "2+" ~ "pos",
                  .x == "3+" ~ "pos"))
    ),

    Response_pathological =  purrr::map_chr(
      Pathological_complete_response,
      ~(case_when(.x == "no" ~ "npCR",
                  .x == "yes" ~ "pCR"))
    ),
    Response_clinical = if_else(Overall_best_response == "Not Assessed",
                                NA_character_, Overall_best_response) %>%
      str_replace("\\(As assessed by SWOG criteria, Green et al., Invest New Drugs 10:239-253 \\(1992\\)\\)",
                  "") %>%
      str_trim(),
    Response = Response_pathological,

    Arm = Treatment_cohort %>% str_replace(" plus ", "+"), # post therapy profiles
    Histology = Pathological_diagnosis
  )



write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE8465_GPL887_cleaned.tsv"
)



# GSE143846 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE143846.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 44


table(str_detect(xx$Sample_source_name_ch1,"Non-Downstage"))
# FALSE  TRUE
# 9    35
table(xx$Age)
# >40 ≤40
# 27  17
table(xx$Histologic_type)
# D/IS     D/L   D/oth  Ductal Lobular
# 1       1       1      38       3
table(xx$Tumor_subtype )
# HER2         Luminal       Luminal A       Luminal B triple negative Triple negative
# 4               1               7              23               1               8
table(xx$Reference)
# Strategene Universal Human Reference RNA
# 44


xx1 <- xx %>%
  dplyr::mutate(
    Age_bin = if_else(Age == "≤40", "young", "old"),
    Age_detailed = "Young:age<=40;Old:age>40",

    Response_clinical = if_else(str_detect(Sample_source_name_ch1,"Non-Downstage"), "non_downstage", "downstage"),
    Response = if_else(Response_clinical == "downstage", "pCR", "npCR"),

    Patient_id = str_split_fixed(Sample_source_name_ch1, "Downstage_", 2)[, 2],
    Histology = Histologic_type,


    # Stromal compartment expression profiles !!!!!!!!!!!!!!!!!!!!!!!!!1
    # "Neoadjuvant chemotherapy followed the hospital treatment protocol, consisting of 4 cycles of doxorubicin 60 mg/m2 and cyclophosphamide 600 mg/m2, every 21 days, followed by 4 cycles of paclitaxel 174 mg/m2 every 21 days (or 80 mg/m2 weekly for 12 weeks). Response was defined as pathological complete response (PCR) or downstaging to maximum ypT1a-b/ypN0, after chemotherapy."

    Regimen = "neoadj",
    Arm_neoadj = "Paclitaxel+AC"
  )

table(xx1$Arm_neoadj)
# Paclitaxel+AC
# 44

write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE143846_cleaned.tsv"
)




# GSE55374 >>>>>>>>>>>>>>>>
xx <- read_tsv(file = "results/geo_gse_curation/reformatted_clinical_data/GSE55374.tsv")
names(xx) <- clean_names(names(xx))
glimpse(xx) # 36


# Note !!!!!!!!!!!!!!!!!!!!!!
# This dataset is discarded as all patient responded to treatment(aired ER +ve).
# !!!!!!!!!!!!!!!!!!!!!!!!!!!

table(xx$Tissue)
# breat tumor IDC breat tumor ILC
# 18              18
table(xx$Treatment)
# 2 wk endocrine therapy 3 mo endocrine therapy          pre-treatment
# 10                     13                     13
table(xx$Clinical_response)
# Responder
# 36




xx1 <- xx %>%
  dplyr::mutate(
    Er = "pos", # from text

    Response_clinical = "response",
    Response = "pCR",

    Arm = "Letrozol",

    Histology = str_replace(Tissue, "breat tumor ", ""),
    Timepoint = purrr::map_chr(
      Treatment,
      ~(case_when(.x == "pre-treatment" ~ "pre_treatment",
                  .x == "2 wk endocrine therapy" ~ "mid_treatment_2weeks",
                  .x == "3 mo endocrine therapy" ~ "mid_treatment_3months",))
    )
  )



write_tsv(
  x = xx1,
  path = "results/geo_gse_curation/reformatted_clinical_data/GSE55374_cleaned.tsv"
)

#
# ==============================================================================




# geo_clean: append reformatted clinical and expression(log scale) data to geo and subset it
# ==============================================================================


# Clinical
# >>>>>>>>>>

files <- list.files("results/geo_gse_curation/reformatted_clinical_data", full.names = TRUE)
files <- files[str_detect(files, "_cleaned\\.tsv")]
id <- str_split_fixed(files,"/",4)[,4]
id <- str_split_fixed(id, "_cleaned.tsv", 2)[, 1]
names(files) = id


# Filter formatted clinical data of sereis matrix  discrded during reformatting step
nme = c("GSE75678","GSE55348","GSE8465_GPL1390","GSE8465_GPL887","GSE55374")
table(names(files) %in% nme)
table(names(geo) %in% nme)
# FALSE  TRUE
# 39     5

files <- files[!(names(files) %in% nme)]
geo_clean <- geo[!(names(geo) %in% nme)]

# rm(geo)
# gc()

for(nme in names(files)){

  print(nme)
  ifile = files[nme]

  geo_clean[[nme]]$clinical <-
    read_tsv(file = ifile, col_types = cols(.default = col_character()))
}
length(geo_clean) # 39


# Regimen column must present in all clinical dataframes
table(purrr::map_lgl(geo_clean,~("Regimen" %in% names(.x$clinical))))
# TRUE
# 39

# Arm column should not present in all clinical dataframes
table(purrr::map_lgl(geo_clean,~("Arm" %in% names(.x$clinical))))
# FALSE  TRUE
# 38     1
purrr::map_lgl(geo_clean,~("Arm" %in% names(.x$clinical)))
# GSE50948 = TRUE, This "Arm" column come from original sample_characteristics



# Expression
# >>>>>>>>>>


# identify class of all expression profiles
for(inme in names(geo_clean)){
  xx <- geo_clean[[inme]]$expression[, -1]
  cls <- purrr::map_chr(xx,~(class(.x))) %>% unique()
  print(str_c(inme, cls, sep = " "))
}

# All numeric except GSE31863 and GSE69031
# !!!!!!!!!!!!!!!!!!!!
# "GSE31863 character" !!!! contains "null" for missing values
# "GSE69031 numeric" "GSE69031 logical" !!!! Few samples got missing expression for all probesets


geo_clean$GSE31863$expression[,-1] <- purrr::map(
  geo_clean$GSE31863$expression[,-1],
  ~(as.numeric(.x))) %>%
  bind_cols()


idx <- purrr::map_lgl(
  geo_clean$GSE69031$expression,
  ~(all(is.na(.x))))
table(idx)
# FALSE  TRUE
# 119    12
# This will be removed later in the script, while exploring NAs in other datasets


# scale
for(inme in names(geo_clean)){

  xx <- geo_clean[[inme]]$expression[, -1]
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



# Explore expression datasets with NA expression values
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



# Filter NAs  from geo_clean
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
na_index <- purrr::map(geo_clean[nme],~(explore_na(.x$expression)))

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
# TRUE
# 110
# [1] "gene"
# FALSE  TRUE
# 36286  4714

# [1] "GSE21997_GPL1390" ##### remove genes
# [1] "sample"
# TRUE
# 35
# [1] "gene"
# FALSE  TRUE
# 21631   944

# [1] "GSE21997_GPL5325" ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 1    31
# [1] "gene"
# FALSE  TRUE
# 43898   392


# [1] "GSE21997_GPL7504"  ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 4    27
# [1] "gene"
# FALSE  TRUE
# 44940   280

# [1] "GSE22226_GPL1708"  ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 6   124
# [1] "gene"
# FALSE  TRUE
# 28349 15941

# [1] "GSE22226_GPL4133"  ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 1    19
# [1] "gene"
# FALSE  TRUE
# 44795   425

# [1] "GSE22358"  ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 2   156
# [1] "gene"
# FALSE  TRUE
# 41810  2480

# [1] "GSE31863"  ##### remove genes
# [1] "sample"
# TRUE
# 143
# [1] "gene"
# FALSE  TRUE
# 11851 14915

# [1] "GSE32603"  ##### remove genes
# [1] "sample"
# TRUE
# 248
# [1] "gene"
# FALSE  TRUE
# 12901 22168

# [1] "GSE4056"  ##### remove dataset!!!! more than 90% genes have atleast one NA value
# [1] "sample"
# TRUE
# 190
# [1] "gene"
# FALSE  TRUE
# 149 21739

# [1] "GSE66999" ##### remove genes
# [1] "sample"
# FALSE  TRUE
# 1    75
# [1] "gene"
# FALSE  TRUE
# 38706  2322

# [1] "GSE69031" ##### remove samples(12 samples has no expression)
# [1] "sample"
# FALSE  TRUE
# 118    12
# [1] "gene"
# TRUE
# 22215



# !!!!!!!!! Removing GSE4056 !!!!!!!!!!!!!!!!!!!!
#
geo_clean <- geo_clean[-which(names(geo_clean) == "GSE4056")]




# !!!!!!!!! Removing 12 samples of GSE69031 with no expression !!!!!!!!!!!!!!!!
#
names(geo_clean$GSE69031)
# [1] "series"                 "sample"                 "expression"
# [4] "sample_characteristics" "dataset_issues"         "clinical"

# Remove 12 samples from  sample, expression, sample_characterisitcs and clinical
nme <- names(xx$GSE69031$na_sample)[xx$GSE69031$na_sample]

# clinical, sample,sample_characteristics are of identical row order
identical(geo_clean$GSE69031$sample$Sample_geo_accession, geo_clean$GSE69031$clinical$Sample_geo_accession)
# TRUE
idx = which(geo_clean$GSE69031$sample$Sample_geo_accession %in% nme)
identical(nme, geo_clean$GSE69031$sample$Sample_geo_accession[idx]) # TRUE

geo_clean$GSE69031$sample <- geo_clean$GSE69031$sample[-idx, ]
geo_clean$GSE69031$sample_characteristics <- geo_clean$GSE69031$sample_characteristics[-idx, ]
geo_clean$GSE69031$clinical <- geo_clean$GSE69031$clinical[-idx, ]

identical(names(geo_clean$GSE69031$expression), names(xx$GSE69031$na_sample)) # TRUE
geo_clean$GSE69031$expression <- geo_clean$GSE69031$expression[, -(which(xx$GSE69031$na_sample))]



# !!!!!!!!! Removing NA genes from the rest !!!!!!!!!!!!!!!!
na_index <- na_index[!(names(na_index) %in% c("GSE4056", "GSE69031"))]

for(inme in names(na_index)){
  ixx = na_index[[inme]]
  geo_clean[[inme]]$expression <- geo_clean[[inme]]$expression[!ixx$na_gene, ]
}



# Checking if NAs still exist in geo_clean !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Names of expression data with NAs
nme <- c("GSE130786", "GSE21997_GPL1390", "GSE21997_GPL5325", "GSE21997_GPL7504",
         "GSE22226_GPL1708", "GSE22226_GPL4133", "GSE22358", "GSE31863",
         "GSE32603", "GSE66999", "GSE69031") # "GSE4056" is removed from geo_clean
na_index <- purrr::map(geo_clean[nme],~(explore_na(.x$expression)))

purrr::map_lgl(na_index,~(any(.x$na_sample))) # all dataset has no more NAs


identical(geo_clean$GSE69031$sample$Sample_geo_accession, geo_clean$GSE69031$clinical$Sample_geo_accession) # TRUE

# End of exploring/Filtering NAs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# scale
for(inme in names(geo_clean)){

  xx <- geo_clean[[inme]]$expression[, -1]
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


summary(geo_clean$GSE114403$expression[, -1] %>% unlist())
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -8.55    -2.49    -0.84    56.18     0.95 47081.00
# Nanostring data; filter NEG and POS controls, everythig else is in log2 scale

summary(geo_clean$GSE143222$expression[, -1] %>% unlist())
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00     1.63     5.60   110.07    23.04 78830.53
# Nanostring data; No control present, convert to log2 scale

summary(geo_clean$GSE19615$expression[, -1] %>% unlist())
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.402    24.704    53.993   163.672   135.482 14222.174
# u133plus2 convert to log2 scale

summary(geo_clean$GSE20271$expression[, -1] %>% unlist())
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0   112.1   194.6   399.9   342.5 31168.7
# u133a convert to log2 scale



# !!!!!!!!!!!! GSE114403 filter neg and pos control probes !!!!!!!!!!!!!!!!!!!!!
xx <- geo_clean$GSE114403$expression %>%
  dplyr::filter(!str_detect(ID_REF,"POS_")) %>%  # 6
  dplyr::filter(!str_detect(ID_REF,"NEG_")) # 8
range(xx[,-1])
# [1] -8.549154 13.329376 # log2 scale
dim(xx) # 770 101
dim(geo_clean$GSE114403$expression) # 784 101

geo_clean$GSE114403$expression <- xx


# !!!!!!!!!!!! GSE143222 convert to log2 scale !!!!!!!!!!!!!!!!!!!!!
xx <- geo_clean$GSE143222$expression
range(xx[,-1])
# 0.00 78830.53
xx[,-1] <- log2(xx[,-1] + (min(xx[,-1])+0.1))
range(xx[,-1])
# -3.321928 16.266469

geo_clean$GSE143222$expression <- xx


# !!!!!!!!!!!! GSE19615 convert to log2 scale !!!!!!!!!!!!!!!!!!!!!
xx <- geo_clean$GSE19615$expression
range(xx[,-1])
# 0.402 14222.174
xx[,-1] <- log2(xx[,-1] + (min(xx[,-1])+0.1))
range(xx[,-1])
# -0.1456053 13.7959053

geo_clean$GSE19615$expression <- xx



# !!!!!!!!!!!! GSE20271 convert to log2 scale !!!!!!!!!!!!!!!!!!!!!
xx <- geo_clean$GSE20271$expression
range(xx[,-1])
# 0.0 31168.7
xx[,-1] <- log2(xx[,-1] + (min(xx[,-1])+0.1))
range(xx[,-1])
# -3.321928 14.927815

geo_clean$GSE20271$expression <- xx



# Double checking scale of xpression matrices

# scale
for(inme in names(geo_clean)){

  xx <- geo_clean[[inme]]$expression[, -1]
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
geo_clean <- purrr::map(
  geo_clean,
  function(x){
    print(class(x$expression$ID_REF))
    x$expression$ID_REF <- as.character(x$expression$ID_REF)
    x
  }
)



# # Save RData !!!!!!!!
# save(geo_clean, file = str_c(outdir,"geo_clean.RData"))
# load(str_c(outdir,"geo_clean.RData"))
#
# ==============================================================================




# Integrate geo expression data
# ==============================================================================

load(str_c(outdir,"geo_clean.RData"))


table(purrr::map_chr(geo_clean,~class(.x$expression$ID_REF)))
# All character
purrr::map(geo_clean,~(.x$expression[1:3,1:5]))
# ID_REF char, all others double



# Maxvar collapsing of expresion values
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

gse_gpl_map <- tibble(
  Series_accession = str_split_fixed(names(geo_clean), "_", 2)[, 1],
  Series_matrix_accession = names(geo_clean),
  Series_matrix_platform_id = purrr::map_chr(geo_clean, ~(unique(.x$clinical$Sample_platform_id)))
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


# Accession number formats
# https://www.ncbi.nlm.nih.gov/Sequin/acc.html
# http://www.nslc.wustl.edu/elgin/genomics/bio4342/1archives/2006/AccReference.pdf


# Load annotations
files = list.files("data/geo_data3/gpl", full.names = TRUE)
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

names(gpl) <-  str_split_fixed(names(gpl),"/", 4)[, 4]
names(gpl) <-  str_split_fixed(names(gpl),"-", 2)[, 1]


# columns in gpl
purrr::map(gpl, names) %>% unlist() %>% table()
# Ensembl_gene_id Ensembl_transcript_id        Entrez_gene_id
# 2                     5                    13
# Genbank_accession      Hugo_gene_symbol                Id_ref
# 18                    17                    19
# Illumina_gene   Illumina_transcript            Probe_name
# 2                     2                     2
# Refseq_accession
# 12


# Checking congruence between geo expression matrix and gpl platform
gse_gpl_cong <- purrr::map2_dfr(
  gse_gpl_map$Series_matrix_accession,
  gse_gpl_map$Series_matrix_platform_id,
  function(nme_gse, nme_gpl, geo_clean, gpl){

    xgse <- geo_clean[[nme_gse]]$expression
    xgpl <- gpl[[nme_gpl]] %>% dplyr::rename(ID_REF = "Id_ref")

    list(
      # eq:equal, ex:expression, an:annotation
      n_eq = nrow(xgse) == nrow(xgpl),
      n_ex =nrow(xgse),
      n_ex2an = sum(xgse$ID_REF %in% xgpl$ID_REF),
      prop_ex2n = sum(xgse$ID_REF %in% xgpl$ID_REF)/nrow(xgse) %>% round(digits = 2),
      n_an = nrow(xgpl),
      n_an2ex = sum(xgpl$ID_REF %in% xgse$ID_REF),
      prop_an2ex = sum(xgpl$ID_REF %in% xgse$ID_REF)/nrow(xgpl) %>% round(digits = 2),
      gse = nme_gse,
      gpl = nme_gpl
    )
  },
  geo_clean,
  gpl
)


# Due to removal of genes which has at least one NA in the dataset causes
# non-congruense between expression and annotation.
# Focus on platforms/series with missing annotations(ie. prop_ex2n < 1; annotation issue).
# Ignore missing probes in expression data (ie. prop_an2ex < 1; NA expression issue), as it
# is due to removal of probes with NA expression values.

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




gpl_clean <- purrr::map(
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


# columns in gpl_clean
purrr::map(gpl_clean, names) %>% unlist() %>% table()
# Ensembl_gene_id Hugo_gene_symbol           ID_REF     Ncbi_gene_id
#              19               19               19               19



xx <- tibble(
  gpl = names(gpl_clean),
  Size = purrr::map_int(gpl_clean, nrow),
  Unique_probes = purrr::map_int(gpl_clean, ~(length(unique(.x$ID_REF)))),
  Unique_genes = purrr::map_int(gpl_clean, ~(length(unique(.x$Ncbi_gene_id)))),
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
# 1  FALSE   815     815         1   819     815  0.9951160        GSE109710 GPL24546   791           791          786
# 2  FALSE   770     770         1   788     770  0.9771574        GSE114403 GPL24993   733           733          733
# 3  FALSE 36286   36286         1 41108   36286  0.8826992        GSE130786  GPL6480 29504         29504        18447
# 4  FALSE   579     579         1   598     579  0.9682274        GSE143222 GPL18649   559           559          559
# 5   TRUE 42545   42545         1 42545   42545  1.0000000        GSE143846 GPL14550 27630         27630        19917
# 6   TRUE 54675   54675         1 54675   54675  1.0000000         GSE16391   GPL570 39766         39766        18943
# 7   TRUE 54675   54675         1 54675   54675  1.0000000         GSE16446   GPL570 39766         39766        18943
# 8   TRUE 54675   54675         1 54675   54675  1.0000000         GSE18728   GPL570 39766         39766        18943
# 9   TRUE 54675   54675         1 54675   54675  1.0000000         GSE19615   GPL570 39766         39766        18943
# 10  TRUE 22283   22283         1 22283   22283  1.0000000         GSE20194    GPL96 19334         19334        12147
# 11  TRUE 22283   22283         1 22283   22283  1.0000000         GSE20271    GPL96 19334         19334        12147
# 12 FALSE 54627   54627         1 54675   54627  0.9991221         GSE20685   GPL570 39766         39766        18943
# 13 FALSE 41000   41000         1 41108   41000  0.9973728         GSE21974  GPL6480 29504         29504        18447
# 14 FALSE 21631   21631         1 22579   21631  0.9580141 GSE21997_GPL1390  GPL1390 10234         10234         9924
# 15 FALSE 43898   43898         1 44294   43898  0.9910597 GSE21997_GPL5325  GPL5325 32443         32443        16614
# 16 FALSE 44940   44940         1 45220   44940  0.9938080 GSE21997_GPL7504  GPL7504 34633         34633        17056
# 17  TRUE 22283   22283         1 22283   22283  1.0000000         GSE22093    GPL96 19334         19334        12147
# 18 FALSE 24332   24332         1 24385   24332  0.9978265         GSE22219  GPL6098 11777         11777         9601
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




# Save gpl_clean !!!!!!!!
#
# save(gpl_clean, file = str_c(outdir, "gpl_clean.RData"))
# load(str_c(outdir, "gpl_clean.RData"))




# Maxvar collapseing >>>>>>>>>>>>>>>>>>>>>>>

# geo_clean: selected datasets
# gse_gpl_map: series matrix to annotation dtaframe mapping
# gpl_clean: cleaned gpl

geo_clean_maxvar_annot <- purrr::map(
  names(geo_clean),
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
  xgse_list = geo_clean,
  xgse_gpl_map = gse_gpl_map,
  xgpl_list = gpl_clean
)
names(geo_clean_maxvar_annot) <- names(geo_clean)


# Save RData !!!!!!!
#
# save(geo_clean_maxvar_annot, file = str_c(outdir, "geo_clean_maxvar_annot.RData"))
# load(str_c(outdir, "geo_clean_maxvar_annot.RData"))




# Checking congruence between geo maxvar-expression matrix and gpl-clean platform
 gse_gpl_cong2 <- purrr::map2_dfr(
  gse_gpl_map$Series_matrix_accession,
  gse_gpl_map$Series_matrix_platform_id,
  function(nme_gse, nme_gpl, geo_clean, gpl_clean, geo_clean_maxvar_annot){

    xgse <- geo_clean[[nme_gse]]$expression
    xgpl <- gpl_clean[[nme_gpl]]
    xannot <-geo_clean_maxvar_annot[[nme_gse]]

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
  geo_clean,
  gpl_clean,
  geo_clean_maxvar_annot
)

 gse_gpl_cong2 %>% as.data.frame()
 # Expectation:
 # prop_exmv2ex muat be 1, ie. all maxvar genes must present in full expr matrix
 # prop_anup2exmv close to 1, ie. proportion of unique genes in gpl present in gse.
 #  Note1: genes are removed if any sample has got NA expression,
 #  hence prop_anup2exmv can be less han 1. More missingness could be a problem.
 #  Note2: gse34138(GPL6884) probes are missing in the original series matrix.

 #     n_ex n_exmv prop_exmv2ex  n_an n_anuq n_anup2exmv prop_anup2exmv              gse      gpl
 # 1    815    786            1   791    786         786      1.0000000        GSE109710 GPL24546
 # 2    770    733            1   733    733         733      1.0000000        GSE114403 GPL24993
 # 3  36286  17146            1 29504  18447       17146      0.9294736        GSE130786  GPL6480
 # 4    579    544            1   559    559         544      0.9731664        GSE143222 GPL18649
 # 5  42545  19917            1 27630  19917       19917      1.0000000        GSE143846 GPL14550
 # 6  54675  18943            1 39766  18943       18943      1.0000000         GSE16391   GPL570
 # 7  54675  18943            1 39766  18943       18943      1.0000000         GSE16446   GPL570
 # 8  54675  18943            1 39766  18943       18943      1.0000000         GSE18728   GPL570
 # 9  54675  18943            1 39766  18943       18943      1.0000000         GSE19615   GPL570
 # 10 22283  12147            1 19334  12147       12147      1.0000000         GSE20194    GPL96
 # 11 22283  12147            1 19334  12147       12147      1.0000000         GSE20271    GPL96
 # 12 54627  18943            1 39766  18943       18943      1.0000000         GSE20685   GPL570
 # 13 41000  18447            1 29504  18447       18447      1.0000000         GSE21974  GPL6480
 # 14 21631   9560            1 10234   9924        9560      0.9633212 GSE21997_GPL1390  GPL1390
 # 15 43898  16562            1 32443  16614       16562      0.9968701 GSE21997_GPL5325  GPL5325
 # 16 44940  17016            1 34633  17056       17016      0.9976548 GSE21997_GPL7504  GPL7504
 # 17 22283  12147            1 19334  12147       12147      1.0000000         GSE22093    GPL96
 # 18 24332   9593            1 11777   9601        9593      0.9991668         GSE22219  GPL6098
 # 19 28349  13042            1 27531  17650       13042      0.7389235 GSE22226_GPL1708  GPL1708 ! moderate missingness
 # 20 44795  18182            1 30844  18272       18182      0.9950744 GSE22226_GPL4133  GPL4133
 # 21 41810  16193            1 32443  16614       16193      0.9746599         GSE22358  GPL5325
 # 22 22283  12147            1 19334  12147       12147      1.0000000         GSE23988    GPL96
 # 23 22283  12147            1 19334  12147       12147      1.0000000         GSE25066    GPL96
 # 24 54675  18943            1 39766  18943       18943      1.0000000         GSE28844   GPL570
 # 25 11851   7263            1 14841  13788        7263      0.5267624         GSE31863 GPL14374 ! lagre missingness
 # 26 12901   5934            1 17912  11483        5934      0.5167639         GSE32603 GPL14668 ! lagre missingness
 # 27 54613  18943            1 39766  18943       18943      1.0000000         GSE32646   GPL570
 # 28 27506  16071            1 29098  19179       16071      0.8379478         GSE34138  GPL6884
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



identical(names(geo_clean), names(geo_clean_maxvar_annot)) # TRUE

# maxvar expression matrix list
geo_clean_expr <- purrr::map2(
  geo_clean,
  geo_clean_maxvar_annot,
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
    y1
  }
)


# Save RData !!!!!!!!!!
#
# save(geo_clean_expr, file = str_c(outdir, "geo_clean_expr.RData"))
# load(str_c(outdir, "geo_clean_expr.RData"))




# Select dataset and merge them
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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



# Common genes in all, except nanostring series and low quality datasets
# Low quality is defined as poor representation of transcriptome
# due to gene missingness.
uid = geo_clean_maxvar_annot$GSE130786$Ncbi_gene_id # non-naostring series
for(inme in names(geo_clean_maxvar_annot)){
  if(!(inme %in% c("GSE109710", "GSE114403", "GSE143222", # Nanostring
                   "GSE21997_GPL1390", # nearly halved unique genes; Agilent
                   "GSE22219", # nearly halved unique genes; GPL6098 Illumina
                   "GSE22226_GPL1708", # moderate missingness; Agilent
                   "GSE31863", # large missingness; GPL14374 Swegene
                   "GSE32603", # large missingness; GPL14668 UCSF_DeRisi_Lab
                   "GSE34138" # small missingness; GPL6884 Illumina
                   )))
    uid <- intersect(uid, geo_clean_maxvar_annot[[inme]]$Ncbi_gene_id)
  print(paste(inme, length(uid)))
}
uid_geo = uid #  9184


# # gpl manufacturer map >>>>>>>>>>>>>>
#
# platform <- tibble(
#   Series_matrix_platform_id = gse_gpl_map$Series_matrix_platform_id %>% unique(),
#   Platform_manufacturer = NA,
#   Platform_probe_number = NA,
#   Platform_title = NA
# ) %>%
#   dplyr::mutate(
#     Platform_detailed = purrr::map_chr(
#       Series_matrix_platform_id,
#       ~(case_when(
#         .x == "GPL24546" ~ "Nanostring///NanoString nCounter human Gene Expression Custom CodeSet///815",
#         .x == "GPL24993" ~ "Nanostring///Nanostring Immune Oncology 360 platform///784",
#         .x == "GPL6480" ~ "Agilent///Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)///41108",
#         .x == "GPL18649" ~ "Nanostring///NanoString nCounter GX Human Immunology v2///594",
#         .x == "GPL14550" ~ "Agilent///Agilent-028004 SurePrint G3 Human GE 8x60K Microarray (Probe Name Version)///42545",
#         .x == "GPL570" ~ "Affymetrix///[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array///54675",
#         .x == "GPL96" ~ "Affymetrix///[HG-U133A] Affymetrix Human Genome U133A Array///22283",
#         .x == "GPL1390" ~ "Agilent///Agilent Human 1A Oligo UNC custom Microarrays///22575",
#         .x == "GPL5325" ~ "Agilent///Agilent UNC Perou Lab Homo sapiens 1X44K Custom Array///44290",
#         .x == "GPL7504" ~ "Agilent///Agilent Axon scanner UNC custom 4X44K without Virus///45220",
#         .x == "GPL6098" ~ "Illumina///Illumina humanRef-8 v1.0 expression beadchip///24385",
#         .x == "GPL1708" ~ "Agilent///Agilent-012391 Whole Human Genome Oligo Microarray G4112A (Feature Number version)///44290", # !!!!!!!!!
#         .x == "GPL4133" ~ "Agilent///Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Feature Number version)///45220",
#         .x == "GPL14374" ~ "Swegene///SWEGENE Human_v2.1.1 55K_condensed///26819", #!!!!!!!!
#         .x == "GPL14668" ~ "UCSF_DeRisi_Lab///UCSF/HAQQLAB_Human_40986_ISPY///35069", #!!!!!!!!!!!
#         .x == "GPL6884" ~ "Illumina///Illumina HumanWG-6 v3.0 expression beadchip///48803", # !!!!!!!!!!
#         .x == "GPL571" ~ "Affymetrix///[HG-U133A_2] Affymetrix Human Genome U133A 2.0 Array///22277",
#         .x == "GPL1352" ~ "Affymetrix///[U133_X3P] Affymetrix Human X3P Array///61359",
#         .x == "GPL6947" ~ "Illumina///Illumina HumanHT-12 V3.0 expression beadchip///49576"
#       ))
#     ),
#     Platform_manufacturer = str_split_fixed(Platform_detailed, "///", 3)[, 1],
#     Platform_probe_number = str_split_fixed(Platform_detailed, "///", 3)[, 3],
#     Platform_title = str_split_fixed(Platform_detailed, "///", 3)[, 2]
#   )
#
# table(platform$Platform_manufacturer)
# # Affymetrix         Agilent        Illumina      Nanostring         Swegene UCSF_DeRisi_Lab
# # 4               7               3               3               1               1
#
#
#
# # Manufacturerer based datset grouping >>>>>>>>>>>>>>>>
#
# # Affymetrix
# uid = gpl_clean$GPL570$Ncbi_gene_id
# for(inme in platform$Series_matrix_platform_id[platform$Platform_manufacturer == "Affymetrix"]){
#   uid <- intersect(uid, gpl_clean[[inme]]$Ncbi_gene_id)
# }
# length(uid) # 11898
# uid_affy <- uid
#
#
# # Agilent
# uid = gpl_clean$GPL6480$Ncbi_gene_id
# for(inme in platform$Series_matrix_platform_id[platform$Platform_manufacturer == "Agilent"]){
#   uid <- intersect(uid, gpl_clean[[inme]]$Ncbi_gene_id)
# }
# length(uid) # 9342
# uid_agi <- uid
#
#
# # Illumina
# uid = gpl_clean$GPL6098$Ncbi_gene_id
# for(inme in platform$Series_matrix_platform_id[platform$Platform_manufacturer == "Illumina"]){
#   uid <- intersect(uid, gpl_clean[[inme]]$Ncbi_gene_id)
# }
# length(uid) # 9487
# uid_illu <- uid
#
#
# # Nanostring
# uid = gpl_clean$GPL24546$Ncbi_gene_id
# for(inme in platform$Series_matrix_platform_id[platform$Platform_manufacturer == "Nanostring"]){
#   uid <- intersect(uid, gpl_clean[[inme]]$Ncbi_gene_id)
# }
# length(uid) # 96
# uid_nano <- uid
#
#
# # Others
# uid_other <- intersect(gpl_clean$GPL14374$Ncbi_gene_id,
#                        gpl_clean$GPL14668$Ncbi_gene_id)
# length(uid_other) # 9129
#
#
# # affy-agi-illu platforms
# length(intersect(uid_affy,uid_illu)) # 7563
# length(intersect(uid_affy,uid_agi)) # 4700
# length(intersect(uid_agi,uid_illu)) # 3645
# uid_aai = intersect(uid_agi, intersect(uid_affy, uid_illu))
# length(uid_aai) # 3161
#
# uid_aai_v1 = uid_aai
#
#
#
# # all platforms
# length(intersect(gpl_clean$GPL14374$Ncbi_gene_id, uid_aai)) # 2114
# length(intersect(gpl_clean$GPL14668$Ncbi_gene_id,uid_aai)) # 2544
# uid_all = intersect(uid_other, uid_aai)
# length(uid_all) # 1768
#
# uid_all_v1 = uid_all
#
#
# length(intersect(uid_nano, uid_all)) # 15 :)
#
#
#
# # Proportion of module genes present in uid_all_v1/uid_aai_v1
# load("results/scratch/module_list.RData")
#
# uid_module <- purrr::map(module_list,~(.x$Ncbi_gene_id)) %>% unlist() %>% unique() %>%
#   str_replace("ncbi_", "")
# length(uid_module) # 1425
#
# length(intersect(uid_module, uid_all))/length(uid_module) # 0.18 !!!!!!!!!
# length(intersect(uid_module, uid_aai))/length(uid_module) # 0.30 !!!!!!!!!
#
# length(intersect(uid_module, uid_affy))/length(uid_module) # 0.85
# length(intersect(uid_module, uid_agi))/length(uid_module) # 0.39 !!!!!!!!!
# length(intersect(uid_module, uid_illu))/length(uid_module) # 0.0.77
#
# # Note:
# # Some agilent arrays are non-standard and missing many module genes !!!!!!!!!
#
#
#
# # Identify and remove non-standard agilent array >>>>>>>>>>>>>>>>>>>
#
# # Agilent
# uid = gpl_clean$GPL6480$Ncbi_gene_id
# for(inme in platform$Series_matrix_platform_id[platform$Platform_manufacturer == "Agilent"]){
#   if(!(inme %in% c("GPL1390","GPL5325")))
#     uid <- intersect(uid, gpl_clean[[inme]]$Ncbi_gene_id)
#   print(paste(inme, length(uid)))
# }
# uid_agi = uid
# # ignoring GPL1390 and GPL5325 !!!!!!!!
# # [1] "GPL6480 18447"
# # [1] "GPL14550 18174"
# # [1] "GPL1390 18174"
# # [1] "GPL5325 18174"
# # [1] "GPL7504 16625"
# # [1] "GPL1708 16299"
# # [1] "GPL4133 16267"
#
# # ignoring GPL1390 !!!!!!!!
# # [1] "GPL6480 18447"
# # [1] "GPL14550 18174"
# # [1] "GPL1390 18174"
# # [1] "GPL5325 9820"
# # [1] "GPL7504 9646"
# # [1] "GPL1708 9568"
# # [1] "GPL4133 9556"
#
# # all
# # [1] "GPL6480 18447"
# # [1] "GPL14550 18174"
# # [1] "GPL1390 9778"
# # [1] "GPL5325 5915"
# # [1] "GPL7504 5840"
# # [1] "GPL1708 5824"
# # [1] "GPL4133 5823"
#
#
# # Note:
# # Discarding the agilent platforms GPL1390  and GPL5325 increases the number
# # of common genes in the remaining agilent platforms.
#
# # Proportion of module genes present in updated uid_agi
# length(intersect(uid_module, uid_agi))/length(uid_module) # 0.93 (old:0.39)
# # !!!!!!!! improved module gene presence in aglent arrays
#
#
#
#
# # Updated within and across platform genes >>>>>>>>>>>>>>>>>>>>
#
# length(intersect(uid_affy,uid_illu)) # 7563
# length(intersect(uid_affy,uid_agi)) # 11320 !!! updated
# length(intersect(uid_agi,uid_illu)) # 9014 !!!updated
# uid_aai = intersect(uid_agi, intersect(uid_affy, uid_illu)) #!!! updated
# length(uid_aai) # 7350
#
# uid_aai_v2 = uid_aai
#
#
# length(intersect(gpl_clean$GPL14374$Ncbi_gene_id, uid_aai)) # 4487
# length(intersect(gpl_clean$GPL14668$Ncbi_gene_id,uid_aai)) # 5698
# uid_all = intersect(uid_other, uid_aai)
# length(uid_all) # 3670
#
# uid_all_v2 = uid_all
#
# length(intersect(uid_nano, uid_all)) # 39 :)
#
#
#
# # Proportion of module genes present in updated common genesets
#
# load("results/scratch/module_list.RData")
#
# uid_module <- purrr::map(module_list,~(.x$Ncbi_gene_id)) %>% unlist() %>% unique() %>%
#   str_replace("ncbi_", "")
# length(uid_module) # 1425
#
# length(intersect(uid_module, uid_all))/length(uid_module) # 0.37 !!! updated
# length(intersect(uid_module, uid_aai))/length(uid_module) # 0.68 !!! updated
#
# length(intersect(uid_module, uid_affy))/length(uid_module) # 0.85
# length(intersect(uid_module, uid_agi))/length(uid_module) # 0.93 !!! updated
# length(intersect(uid_module, uid_illu))/length(uid_module) # 0.0.77
#
#
# # Potential dataset combinations
# length(uid_all_v1) # 1768 # all non-nanostring dataset
# length(uid_all_v2) # 3670 # all non-nanostring dataset except GPL1390, GPL5325 agi platforms
# length(uid_aai_v1) # 3161 # affy-agi-illu datasets
# length(uid_aai_v2) # 7350 # affy-agi-illu datasets except GPL1390, GPL5325 agi platforms




# Is er/pr/her2/pam50 genes present in uid_geo#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# ESR1 "2099" %in% uid_geo # TRUE
# ERBB2 "2064" %in% uid_geo # TRUE
# PGR "5241" %in% uid_geo # TRUE


# pam50 genes in
pam50 <- read_tsv("data/pam50/PAM50/bioclassifier_R/pam50_annotation.txt")
table(pam50$EntrezGene %>% as.character() %in% uid_geo)
# FALSE  TRUE
# 14    36



# Validate with external data whether signture subsets derived from
# uid_geo are congruent with the full signatures.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# Module list
load("/home/josev/proj/on/TIL_DTX_NVB/data/sig/module_list.RData")
module_list <- purrr::map(
  module_list,
  function(x){
    x %>%
      dplyr::mutate_at("Ncbi_gene_id", ~str_c("ncbi_", .x))
  }
)



names(module_list)
# [1] "Desmedt08_immune_response"     "Teschendorff07_immune2_erneg"
# [3] "Hamy16_immunity"               "Farmer09_interferon_mx1"
# [5] "Hamy16_interferon"             "Bergamaschi08_ecm1"
# [7] "Hamy16_matrix"                 "Desmedt08_proliferation"
# [9] "Farmer09_proliferation_tpx2"   "Dai2005_prolif"
# [11] "Lundberg2017_prolif"           "Ehmsen2019_choles"
# [13] "Kimbung2016_choles"            "Kuzhu2016_choles"
# [15] "Sorrentino2014_choles"         "Simigdala2016_choles"
# [17] "Triulzi2013_ecm"               "Nirmal2018_b_cells"
# [19] "Nirmal2018_t_cells"            "Nirmal2018_macrophages"
# [21] "Nirmal2018_monocytes"          "Nirmal2018_neutrophils"
# [23] "Nirmal2018_nk_cells"           "Nirmal2018_plasma_cells"
# [25] "Nirmal2018_interferon"         "Nirmal2018_proliferation"
# [27] "Nirmal2018_translation"        "Tamborero2018_b_cells"
# [29] "Tamborero2018_eosinophils"     "Tamborero2018_macrophages"
# [31] "Tamborero2018_mast_cells"      "Tamborero2018_nkbrigh"
# [33] "Tamborero2018_nkdim"           "Tamborero2018_neutrophils"
# [35] "Tamborero2018_th"              "Tamborero2018_tcm"
# [37] "Tamborero2018_tem"             "Tamborero2018_tfh"
# [39] "Tamborero2018_adc"             "Tamborero2018_idc"
# [41] "Tamborero2018_cd8_t"           "Tamborero2018_tgd"
# [43] "Tamborero2018_treg"            "Tamborero2018_cytotoxic_cells"
# [45] "Truntzer2018_lymphoids"        "Truntzer2018_myeloids"
# [47] "Truntzer2018_cancer"           "Truntzer2018_stroma"
# [49] "Naba2014_ecm_core"             "Naba2014_ecm_low_metastasis"
# [51] "Naba2014_ecm_high_metastasis"  "Yang2018_immune"
# [53] "Yang2018_prolif"               "Gruosso2019_Immune"
# [55] "Gruosso2019_Fibrosis"          "Gruosso2019_Interferon"
# [57] "Gruosso2019_Cholesterol"


# Extracted from "TIL_DTX_NVB/R/til_dtx_nvb_figures.R"
#
# dplyr::rename(
#   TIL_HE = "StrLy_Mean",
#
#   TILsig = "tilsig",
#   Immune_TILsig = "tilsig_up",
#   Fibrosis_TILsig = "tilsig_down",
#
#   Immune1 = "Hamy16_immunity",
#   Immune2 = "Gruosso2019_Immune",
#   Immune3 = "Desmedt08_immune_response",
#   Immune4 = "Yang2018_immune",
#
#   Interferon1 = "Hamy16_interferon",
#   Interferon2 = "Gruosso2019_Interferon",
#   Interferon3 = "Nirmal2018_interferon",
#   Interferon4 = "Farmer09_interferon_mx1",
#
#   Cholesterol1 = "Kimbung2016_choles",
#   Cholesterol2 = "Sorrentino2014_choles",
#
#   Fibrosis1 = "Hamy16_matrix",
#   Fibrosis2 = "Gruosso2019_Fibrosis",
#   Fibrosis3 = "Triulzi2013_ecm",
#   Fibrosis4 = "Naba2014_ecm_core",
#   Fibrosis5 = "Truntzer2018_stroma",
#
#
#   Proliferation1 = "Desmedt08_proliferation",
#   Proliferation2 = "Yang2018_prolif",
#   Proliferation3 = "Nirmal2018_proliferation",
#   Proliferation4 = "Farmer09_proliferation_tpx2",
#   Proliferation5 = "Dai2005_prolif",
#
#
#   T_cells = "Nirmal2018_t_cells",
#   Macrophages = "Nirmal2018_macrophages",
#   Plasma_cells = "Nirmal2018_plasma_cells"
# )


# Tcga (Discovery)
# >>>>>>>>>>>>>
load("results/clean_data_2/tcga_voom.RData")

tcga_clin <- as_tibble(tcga_voom$targets) %>%
  dplyr::filter(!is.na(Til_percentage))
# Note xtil is now tcga_clin !!!!!!!!!!

# Ref: https://doi.org/10.1016/j.cell.2018.02.052
# An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics
tcga_cdr <- read_excel(
  path = "data/liu2018_tcga_cdr_NIHMS978596-supplement-1.xlsx",
  sheet =  "TCGA-CDR",
  na = c("=#N/A",""),
  guess_max = 5000
)[, -1] %>% # removing row number
  dplyr::rename_all(~(str_c("Cdr_",.x) %>% str_to_title())) %>%
  dplyr::rename(Case_id = "Cdr_bcr_patient_barcode")

tcga_clin <- tcga_clin %>%
  left_join(tcga_cdr, by = "Case_id")

tcga_expr <- as_tibble(tcga_voom$E, rownames = "Ncbi_gene_id") %>%
  t_tibble(names_x_desc = "Case_id") %>%
  dplyr::filter(Case_id %in% tcga_clin$Case_id)
# Note xdata is now tcga_expr !!!!!!!!!!


identical(tcga_clin$Case_id, tcga_expr$Case_id) # TRUE




# Rody TNBC dataset (Validation)
# >>>>>>>>>>>>>
rody_clin <- read_tsv(
  file = "../TIL_DTX_NVB/data/rody/GSE31519_TNBC_SampleInfo_BCR.txt.gz",
  col_names = FALSE) %>%
  t() %>%
  as_tibble(.name_repair = "minimal")

names(rody_clin) <- rody_clin[1,] %>% as.character()
rody_clin <- rody_clin[-1,] %>%
  dplyr::mutate(
    fu_120 = as.numeric(fu_120),
    ev120 = as.integer(ev120),
    age = as.numeric(age)
  )
glimpse(rody_clin)


load("../TIL_DTX_NVB/data/rody/rody_annot_max_var.RData")
load("../TIL_DTX_NVB/data/rody/rody.RData")
rody_expr <- rody_annot_max_var %>%
  inner_join(rody, by = "Probeset_Id") %>%
  dplyr::select(-c("Probeset_Id", "Gene_Symbol")) %>%
  dplyr::rename(Ncbi_gene_id = "Entrez_Gene") %>%
  dplyr::mutate(Ncbi_gene_id = str_c("ncbi_",Ncbi_gene_id)) %>%
  t_tibble(names_x_desc = "GSM_accession") %>%
  dplyr::filter(GSM_accession %in% rody_clin$GSM_accession)

rody_clin <- rody_clin %>%
  dplyr::filter(GSM_accession %in% rody_expr$GSM_accession) %>%
  dplyr::arrange(GSM_accession)

rody_expr <- rody_expr %>%
  dplyr::arrange(GSM_accession)

identical(rody_clin$GSM_accession, rody_expr$GSM_accession) # TRUE




# Metabric (Validation)
# >>>>>>>>>>>>>
load("/home/josev/proj/off/supp/datasets/adjuvant/METABRIC/EXP/METABRIC.RData")
mbric_clin <- as_tibble(demo) %>% dplyr::mutate(samplename = as.character(samplename))
mbric_expr <- as_tibble(data, rownames = "samplename")

identical(as.character(annot$EntrezGene.ID),
          str_replace(as.character(annot$probe),"geneid.","")) # TRUE
rm(annot, data, demo)

names(mbric_expr) <- str_replace(names(mbric_expr),"geneid.","ncbi_")

# removing NAs
index <- purrr::map_lgl(mbric_expr,~(is.na(.x) %>% any()))
table(index)
# FALSE  TRUE
# 24914    11

mbric_expr <- mbric_expr[, !index]

# index <- purrr::map_lgl(t(mbric_expr) %>% as_tibble(), ~(is.na(.x) %>% any()))
# table(index) # all false

dim(mbric_clin)
# [1] 1992 53
dim(mbric_expr)
# [1]   1992 24914

identical(mbric_clin$samplename, mbric_expr$samplename) # TRUE





# Module list ful and subset versions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>


module_list2 <- list(
  full = module_list,
  geo = purrr::map(module_list,
                      function(sig,id){sig %>% dplyr::filter(Ncbi_gene_id %in% id)},
                      id = str_c("ncbi_", uid_geo))
)

size_module <- tibble(
  all = purrr::map_int(module_list2$full, nrow),
  geo = purrr::map_int(module_list2$geo, nrow),
  prop_in_geo = geo/all
  )
summary(size_module$prop_in_geo)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1538  0.6129  0.7083  0.6757  0.7647  0.8889
names(module_list)[size_module$prop_in_geo<.5]
# [1] "Nirmal2018_nk_cells"     "Nirmal2018_plasma_cells" "Nirmal2018_translation"
# [4] "Tamborero2018_nkdim"





# Tcga module score correlations
# >>>>>>>>>>>>>
scorelist <- purrr::map(
  module_list2,
  function(xsiglist, xdata){
    get_module_score(x = xdata, module_list = xsiglist, by = "Ncbi_gene_id")
  },
  xdata = tcga_expr
)

xx <- tibble(
  signame = names(module_list2$full),
  pearson = diag(cor(scorelist$full[,-1],scorelist$geo[,-1], method = "pearson")),
  spearman = diag(cor(scorelist$full[,-1],scorelist$geo[,-1], method = "spearman"))
)
summary(xx$pearson)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.5933  0.9700  0.9892  0.9677  0.9967  1.0000
summary(xx$spearman)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.5801  0.9667  0.9883  0.9651  0.9963  1.0000

xx %>%
  dplyr::filter(pearson < 0.95 | spearman < 0.9)
#   signame                     pearson spearman
# 1 Kimbung2016_choles            0.851    0.824
# 2 Nirmal2018_neutrophils        0.936    0.929
# 3 Tamborero2018_eosinophils     0.869    0.844
# 4 Tamborero2018_nkbrigh         0.593    0.580
# 5 Tamborero2018_nkdim           0.897    0.899
# 6 Tamborero2018_th              0.937    0.934
# 7 Tamborero2018_tem             0.851    0.836
# 8 Tamborero2018_tgd             0.926    0.924
# 9 Naba2014_ecm_low_metastasis   0.888    0.867




# Metabric module score correlations
# >>>>>>>>>>>>.
scorelist <- purrr::map(
  module_list2,
  function(xsiglist, xdata){
    get_module_score(x = xdata, module_list = xsiglist, by = "Ncbi_gene_id")
  },
  xdata = mbric_expr
)

xx <- tibble(
  signame = names(module_list2$full),
  pearson = diag(cor(scorelist$full[,-1],scorelist$geo[,-1], method = "pearson")),
  spearman = diag(cor(scorelist$full[,-1],scorelist$geo[,-1], method = "spearman"))
)
summary(xx$pearson)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7884  0.9501  0.9869  0.9660  0.9958  0.9987
summary(xx$spearman)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7632  0.9371  0.9842  0.9596  0.9952  0.9987

xx %>%
  dplyr::filter(pearson < 0.9 | spearman < 0.9)
#   signame                     pearson spearman
# 1 Kimbung2016_choles            0.889    0.864
# 2 Nirmal2018_neutrophils        0.888    0.879
# 3 Tamborero2018_nkbrigh         0.788    0.763
# 4 Tamborero2018_nkdim           0.906    0.900
# 5 Tamborero2018_tem             0.891    0.876
# 6 Tamborero2018_tgd             0.912    0.893
# 7 Naba2014_ecm_low_metastasis   0.830    0.814




# Rody module score correlations
# >>>>>>>>>>>>>
scorelist <- purrr::map(
  module_list2,
  function(xsiglist, xdata){
    get_module_score(x = xdata, module_list = xsiglist, by = "Ncbi_gene_id")
  },
  xdata = rody_expr
)


xx <- tibble(
  signame = names(module_list2$full),
  pearson = diag(cor(scorelist$full[,-1],scorelist$geo[,-1], method = "pearson")),
  spearman = diag(cor(scorelist$full[,-1],scorelist$geo[,-1], method = "spearman"))
)
summary(xx$pearson)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.6480  0.9532  0.9828  0.9632  0.9944  0.9987
summary(xx$spearman)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.6717  0.9407  0.9783  0.9559  0.9929  0.9982

xx %>%
  dplyr::filter(pearson < 0.9 | spearman < 0.9)
#   signame                     pearson spearman
# 1 Nirmal2018_nk_cells           0.884    0.868
# 2 Tamborero2018_eosinophils     0.928    0.889
# 3 Tamborero2018_nkbrigh         0.797    0.773
# 4 Tamborero2018_nkdim           0.648    0.672
# 5 Tamborero2018_tem             0.889    0.868
# 6 Naba2014_ecm_low_metastasis   0.885    0.858


# dplyr::rename(
#   TIL_HE = "StrLy_Mean",
#
#   TILsig = "tilsig",
#   Immune_TILsig = "tilsig_up",
#   Fibrosis_TILsig = "tilsig_down",
#
#   Immune1 = "Hamy16_immunity",
#   Immune2 = "Gruosso2019_Immune",
#   Immune3 = "Desmedt08_immune_response",
#   Immune4 = "Yang2018_immune",
#
#   Interferon1 = "Hamy16_interferon",
#   Interferon2 = "Gruosso2019_Interferon",
#   Interferon3 = "Nirmal2018_interferon",
#   Interferon4 = "Farmer09_interferon_mx1",
#
#   Cholesterol1 = "Kimbung2016_choles",
#   Cholesterol2 = "Sorrentino2014_choles",
#
#   Fibrosis1 = "Hamy16_matrix",
#   Fibrosis2 = "Gruosso2019_Fibrosis",
#   Fibrosis3 = "Triulzi2013_ecm",
#   Fibrosis4 = "Naba2014_ecm_core",
#   Fibrosis5 = "Truntzer2018_stroma",
#
#
#   Proliferation1 = "Desmedt08_proliferation",
#   Proliferation2 = "Yang2018_prolif",
#   Proliferation3 = "Nirmal2018_proliferation",
#   Proliferation4 = "Farmer09_proliferation_tpx2",
#   Proliferation5 = "Dai2005_prolif",
#
#
#   T_cells = "Nirmal2018_t_cells",
#   Macrophages = "Nirmal2018_macrophages",
#   Plasma_cells = "Nirmal2018_plasma_cells"
# )


# Summary !!!!!!!!!!!!!!!!!!!
# The TIME signatures are significantly correlated to between full and
# reduced versions.
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




# Consolidate geo_clean expression
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Do two versions of analysis
# 1) pooled using uid_geo and by discarding sub-standard datasets, including nanostring
# 2) per dataset based

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
  names(geo_clean_expr)[!(names(geo_clean_expr) %in% nme)],
  function(nme_gse, geo_clean_expr, genes){
    print(nme_gse)
    geo_clean_expr[[nme_gse]] %>%
      dplyr::filter(Ncbi_gene_id %in% all_of(genes)) %>%
      t_tibble(names_x_desc = "Sample_geo_accession") %>%
      dplyr::mutate(Series_matrix_accession = nme_gse) %>%
      dplyr::select("Series_matrix_accession", "Sample_geo_accession", genes)

  },
  geo_clean_expr,
  genes = str_c("ncbi_", uid_geo)
)


# Checking whther samples per dataset are unique
table(purrr::map_lgl(
  geo_expr,~(nrow(.x) == (.x$Sample_geo_accession %>% unique() %>% length()))
)
)
# TRUE
# 29
# Per dataset samples are unique


# samples x genes
geo_expr <- bind_rows(geo_expr)
head2(geo_expr)
# $head
#   Series_matrix_accession Sample_geo_accession ncbi_400451 ncbi_57099 ncbi_23339
# 1 GSE130786               GSM3753582                0.340      0.238     -0.102
# 2 GSE130786               GSM3753583               -0.262      0.0634    -0.204
# 3 GSE130786               GSM3753584                0.0394     0.135      0.160
# 4 GSE130786               GSM3753585               -0.571      0.109     -0.261
# 5 GSE130786               GSM3753586               -0.179      0.135     -0.0628
# $dim
# [1] 3744 9186; samples x genes


# Checking whther samples in pooled dataset are unique
nrow(geo_expr) == (geo_expr$Sample_geo_accession %>% unique() %>% length()) # FALSE
nrow(geo_expr) # 3744
geo_expr$Sample_geo_accession %>% unique() %>% length() # 3740

# Sample redundancy @ Sample_geo_accession!!!!!!!!!!!!!!!!!!
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


# # Exploring Arms of redundat samples !!!!!!!
# geo_clin %>%
#   dplyr::filter(Series_matrix_accession == "GSE21997_GPL5325") %>%
#   dplyr::distinct(Arm)
# #   Arm
# # 1 Docetaxel
# # 2 Doxorubicin
# # 3 NA
# geo_clin %>%
#  dplyr::filter(Series_matrix_accession == "GSE22358") %>%
#  dplyr::distinct(Arm)
# #   Arm
# # 1 Docetaxel+Capecitabine
# # 2 Docetaxel+Capecitabine+Trastuzumab
# # 3 NA



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



# genes x samples
geo_expr_meta <- geo_expr %>%
  dplyr::select(Series_matrix_accession, Sample_geo_accession)
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



# Data transformation comparison >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# original
purrr::map_dbl(geo_expr[, -1],min) %>% summary()
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -11.7290  -0.6769   0.4419   0.4081   2.7142   5.9298
purrr::map_dbl(geo_expr[, -1],max) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.8559 13.7909 14.3345 13.9696 15.6581 20.9422
purrr::map_dbl(geo_expr[, -1],IQR) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1329  2.5338  2.9128  2.7778  3.3356  5.3239


# Rescaled
xx <- geo_expr[, -1] %>% # -1 for Ncbi_gene_id
  purrr::map(~(genefu::rescale(x = .x, q=0.05)) ) %>%
  bind_cols()
head2(xx)

purrr::map_dbl(as_tibble(xx),min) %>% summary()
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -2.78398 -0.42307 -0.20386 -0.35183 -0.09191  0.00000
purrr::map_dbl(as_tibble(xx),max) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.198   1.366   1.481   1.542   1.670   3.377
purrr::map_dbl(as_tibble(xx),IQR) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1852  0.3408  0.3634  0.3617  0.3872  0.5987


# comabt (original)
xx1 <- sva::ComBat(as.matrix(geo_expr[,-1]), #geo_expr[,-1]),
                  batch = geo_expr_meta$Series_matrix_accession)

head2(xx1)
purrr::map_dbl(as_tibble(xx1),min) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -6.3175  0.1916  0.9357  0.8026  1.6216  2.6951
purrr::map_dbl(as_tibble(xx1),max) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 10.91   12.60   12.96   13.28   13.57   20.30
purrr::map_dbl(as_tibble(xx1),IQR) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.498   2.155   2.260   2.275   2.383   4.026


# rescaled + combat
xx2 <- sva::ComBat(as.matrix(xx), #geo_expr[,-1]),
                   batch = geo_expr_meta$Series_matrix_accession)

head2(xx2)
purrr::map_dbl(as_tibble(xx2),min) %>% summary()
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -1.000092 -0.303376 -0.213211 -0.226384 -0.122613  0.001606
purrr::map_dbl(as_tibble(xx2),max) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.146   1.301   1.359   1.399   1.454   2.323
purrr::map_dbl(as_tibble(xx2),IQR) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2129  0.2760  0.2858  0.2883  0.2979  0.5476


# rescaled + combat(meanonly)
xx3 <- sva::ComBat(as.matrix(xx), #geo_expr[,-1]),
                   mean.only = TRUE,
                   batch = geo_expr_meta$Series_matrix_accession)

head2(xx3)
purrr::map_dbl(as_tibble(xx3),min) %>% summary()
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -1.97651 -0.38806 -0.26940 -0.26111 -0.05613  0.07242
purrr::map_dbl(as_tibble(xx3),max) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.136   1.267   1.345   1.395   1.460   3.093
purrr::map_dbl(as_tibble(xx3),IQR) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2128  0.2775  0.2999  0.2982  0.3167  0.4377


# correlation
cor(geo_expr$GSM3753582, xx[, "GSM3753582"]) # rescaled
# [1] 1
cor(geo_expr$GSM3753582, xx1[, "GSM3753582"]) # combat(centering+scaling)
# [1] 0.2633442
cor(geo_expr$GSM3753582, xx2[, "GSM3753582"]) # rescaled + combat(centering+scaling)
# [1] 0.3323717
cor(geo_expr$GSM3753582, xx3[, "GSM3753582"]) # rescaled + combat(centering)
# [1] 0.6831184


# IQR
purrr::map_dbl(geo_expr[, -1],IQR) %>% summary() # original
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1329  2.5338  2.9128  2.7778  3.3356  5.3239
purrr::map_dbl(as_tibble(xx),IQR) %>% summary() # rescaled
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1852  0.3408  0.3634  0.3617  0.3872  0.5987
purrr::map_dbl(as_tibble(xx1),IQR) %>% summary() # combat(centering+scaling)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.498   2.155   2.260   2.275   2.383   4.026
purrr::map_dbl(as_tibble(xx2),IQR) %>% summary() # rescaled + combat(centering+scaling)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2129  0.2760  0.2858  0.2883  0.2979  0.5476
purrr::map_dbl(as_tibble(xx3),IQR) %>% summary() # rescaled + combat(centering)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2128  0.2775  0.2999  0.2982  0.3167  0.4377



# Scaling/Batch correction summary >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# A simple rescaling keep the original structure of the data and
# will give comparable IQR values, simialr to that obtained from combat !!!!

# Rescaling samples
geo_expr[, -1] <- geo_expr[, -1] %>% # -1 for Ncbi_gene_id
  purrr::map(~(genefu::rescale(x = .x, q=0.05)) ) %>%
  bind_cols()
head2(geo_expr)

purrr::map_dbl(geo_expr[, -1],IQR) %>% summary() # original
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1852  0.3409  0.3636  0.3619  0.3872  0.5987



# Save RData !!!!!!!!!!!!
#
# save(geo_expr, file = str_c(outdir, "geo_expr.RData"))
# save(geo_expr_meta, file = str_c(outdir, "geo_expr_meta.RData")) # batch(series information)
# load(str_c(outdir, "geo_expr.RData"))
# load(str_c(outdir, "geo_expr_meta.RData"))
# ==============================================================================




# Integrate geo clinical data
# ==============================================================================


# Geo series and series-matrix summary
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# Sereis summary >>>>>>>>>>>>
# "geo_search_taxane.antracyclin.breast_31july2020_SearchResult_annotated.tsv"
# is the annotated version of GEO search results.
# Annotation include, regimen, sample_procrument, archive_method, end point etc
geo_series_summary <- read_tsv(file = "results/geo_search_result_curation/geo_search_taxane.antracyclin.breast_31july2020_SearchResult_annotated.tsv") %>%
  # contains invalid multibyte strings hence removed
  dplyr::select(-c("Sample_procurement_details", "Archive_details")) %>%
  dplyr::rename_all(~(str_c("Series_", .x) %>%
                        str_to_sentence() %>%
                        str_replace("Series_series","Series")))


glimpse(geo_series_summary)
# Note:
# names(geo_series_summary)[10:ncol(geo_series_summary)] are manually curated to
# get an overview of the dataset.
# Most accurate Arm and Regimen information can be found in the clinical data
# associated with the  R object "geo_clean".
# However "Series_sample_procurement" and "Series_archive_method" can be used.
#
# [1] "Series_arm"                "Series_arm_description"    "Series_survival"
# [4] "Series_regimen"            "Series_sample_procurement" "Series_archive_method"
# [7] "Series_selected"           "Series_comment"


# Series-matrix summary (extracted from geo list; contains only) >>>>>>>
geo_series_matrix_summary <- tibble(
  Series_accession = str_split_fixed(names(geo_clean), "_", 2)[, 1],
  Series_matrix_accession = names(geo_clean),
  Series_matrix_probe_size = purrr::map_int(geo_clean,~nrow(.x$expression)),
  Series_matrix_gene_size = purrr::map_int(geo_clean_maxvar_annot,~nrow(.x)),
  Series_matrix_platform = purrr::map_chr(geo_clean,~(unique(.x$clinical$Sample_platform_id))),
  Series_matrix_regimen =  purrr::map_chr(geo_clean,~(unique(.x$clinical$Regimen))),
  Series_matrix_sample_size = purrr::map_int(geo_clean,~nrow(.x$clinical)),
  # Samples with treatment info
  Series_matrix_regimen_size = purrr::map2_int(
    Series_matrix_regimen, geo_clean,
    ~(if_else(str_detect(.x, "neoadj"),
              sum(!is.na(.y$clinical$Arm_neoadj)),
              sum(!is.na(.y$clinical$Arm_adj))))),
  Series_matrix_has_neoadj = purrr::map_lgl(
    geo_clean, ~(any(str_detect(names(.x$clinical), "_neoadj")))),
  Series_matrix_has_adj = purrr::map_lgl(
    geo_clean, ~(any(str_detect(names(.x$clinical), "_adj")))),
  Series_matrix_selected = NA,
  Series_matrix_comment = NA
)

glimpse(geo_series_matrix_summary)


# Checking congruence between Series_matrix_regimen, Series_matrix_has_neoadj,
#  and Series_matrix_has_adj
geo_series_matrix_summary %>% select(Series_matrix_accession,Series_matrix_regimen,Series_matrix_has_neoadj,Series_matrix_has_adj) %>% as.data.frame()
# expectation: adj regimen has no neoadj arm and vice-versa
# neoadj_adj has both neoadj and adj arm.
# Note that GSE16391 adj dataset seems failed expectations, but this is due to
# an original column wth pattern "_neoadj"


geo_summary <- full_join(
  geo_series_summary, geo_series_matrix_summary,
  by = "Series_accession") %>%
  dplyr::mutate(
    Series_matrix_selected = if_else(is.na(Series_matrix_accession), "no", "yes")
  )
glimpse(geo_summary)

table(geo_summary$Series_matrix_selected)
# no yes
# 137  38
table(geo_summary$Series_regimen) # irrelevant
# adj  cell line     neoadj neoadj/adj
# 15          7         42          2
table(geo_summary$Series_matrix_regimen)
# adj     neoadj neoadj_adj
# 7         27          4

geo_summary %>%
  dplyr::group_by(Series_matrix_selected, Series_matrix_regimen) %>%
  dplyr::summarise(Sum = sum(Series_matrix_sample_size),
                   Sum_regimen = sum(Series_matrix_regimen_size))
#   Series_matrix_selected Series_matrix_regimen   Sum Sum_regimen
# 1 no                     NA                       NA          NA
# 2 yes                    adj                    1113         954
# 3 yes                    neoadj                 3078        2743
# 4 yes                    neoadj_adj              831         724


table(geo_summary$Series_sample_procurement)
# biopsy
# 3
# Biopsy
# 1
# core biopsy
# 1
# core needle biopsy
# 14
# Core needle biopsy
# 2
# core needle biopsy/surgery
# 1
# fine needle biopsy
# 1
# Fine needle biopsy
# 1
# FNA
# 3
# FNA/CBX
# 3
# incisional biopsy
# 1
# radical mastectomy or breast-conserving surgery
# 1
# surgery
# 7
table(geo_summary$Series_archive_method)
# ffpe          frozen        RNAlater RNAlater/Frozen
# 5              28               6               3


geo_summary <- geo_summary %>%
  dplyr::mutate(

    Series_sample_procurement = purrr::map_chr(
      Series_sample_procurement %>% str_to_lower(),
      ~(case_when(
        .x == "biopsy" ~ "needle biopsy",
        .x == "fna/cbx" ~ "needle biopsy",
        .x == "core biopsy" ~ "needle biopsy",
        .x == "incisional biopsy" ~ "needle biopsy",
        .x == "radical mastectomy or breast-conserving surgery" ~ "surgery",
        TRUE ~ .x
      )
      )
    ),

    Series_sample_procurement = Series_sample_procurement %>%
      str_to_lower() %>%
      str_replace("\\/", "+") %>%
      str_replace("fna", "fine needle biopsy") %>%
      str_replace_all(" ", "_") %>%
      str_to_sentence() %>%
      str_replace_all("surgery", "Surgery") %>%
      str_replace_all("core", "Core"),

    Series_archive_method = Series_archive_method %>%
      str_replace("\\/", "|") %>%
      str_replace("ffpe", "FFPE") %>%
      str_replace("frozen", "Frozen")

  )


geo_summary %>%
  dplyr::group_by(Series_matrix_selected,Series_sample_procurement) %>%
  dplyr::summarise(N = n())
#   Series_matrix_selected Series_sample_procurement      N
# 1 no                     Core_needle_biopsy             2
# 2 no                     Needle_biopsy                  4
# 3 no                     Surgery                        5
# 4 no                     NA                           126
# 5 yes                    Core_needle_biopsy            14
# 6 yes                    Core_needle_biopsy+Surgery     1
# 7 yes                    Fine_needle_biopsy             5
# 8 yes                    Needle_biopsy                  5
# 9 yes                    Surgery                        3
# 10 yes                    NA                            10

geo_summary %>%
  dplyr::group_by(Series_matrix_selected,Series_archive_method) %>%
  dplyr::summarise(N = n())
#   Series_matrix_selected Series_archive_method     N
# 1 no                     FFPE                      2
# 2 no                     Frozen                    7
# 3 no                     RNAlater                  2
# 4 no                     RNAlater|Frozen           2
# 5 no                     NA                      124
# 6 yes                    FFPE                      3
# 7 yes                    Frozen                   21
# 8 yes                    RNAlater                  4
# 9 yes                    RNAlater|Frozen           1
# 10 yes                    NA                        9


#
# # Save geo_summary !!!!!!!!!!!!!!!!!!!
# #
# # save(geo_summary, file = str_c(outdir,"geo_summary.RData"))
# # load(str_c(outdir,"geo_summary.RData"))



# Consolidate clinical data of selected sereis matrices
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# adj platform
table(geo_summary %>%
        dplyr::filter(Series_matrix_selected == "yes" & Series_matrix_regimen == "adj") %>%
        dplyr::select(Series_matrix_platform) %>%
        tibble::deframe())
# GPL14374   GPL570   GPL571  GPL6098    GPL96
# SWEGENE    u133p2   u133a2  illumina   u133a
#        1        3        1        1        1
# affy: 5
# illu: 1
# custom: 1


# neoadj platform
table(geo_summary %>%
        dplyr::filter(Series_matrix_selected == "yes" & str_detect(Series_regimen,"neoadj")) %>%
        dplyr::select(Series_matrix_platform) %>%
        tibble::deframe())
# GPL1352  GPL1390 GPL14550 GPL14668  GPL1708 GPL18649 GPL24546 GPL24993
# x3p      Agi     Agi      UCSF      Agi     Nstring  Nstring  Nstring
# 2        1        1        1        1        1        1        1
# GPL4133  GPL5325   GPL570   GPL571  GPL6480  GPL6884  GPL6947  GPL7504
# Agi      Agi       u133p2   u133a2  Agi      illu     illu     Agi
# 1        2        6        1        3        1        1        1
# GPL96
# u133a
# 6

# affy: 13, affyexon:2
# nano: 3,
# agil: 10
# illu: 2
# custom: 1


# For neoadjuvant endpoint is Response
# For adjuvant endpoint id Time_dfs, Event_dfs
nme <- geo_summary %>%
  dplyr::filter(Series_matrix_selected == "yes") %>%
  dplyr::select(Series_matrix_accession) %>%
  tibble::deframe()
# n = 38


geo_clin <- NULL
for(inme in nme){
  print(inme)
  xx <- geo_clean[[inme]]$clinical
  xx$Series_matrix_accession <- inme
  names(xx) <- str_to_sentence(names(xx))
  geo_clin <- bind_rows(geo_clin, xx)
}

dim(geo_clin) # 5022 429
write_tsv(x = tibble(Column_name = sort(names(geo_clin))),
          path = str_c(outdir, "column_name_geo_clin.tsv"))

nme <- c(
  "Series_matrix_accession",
  "Sample_title", "Sample_geo_accession",
  "Sample_channel_count", "Sample_platform_id",
  "Sample_organism_ch1", "Sample_organism_ch2",
  "Sample_source_name_ch1", "Sample_source_name_ch2",
  "Sample_type","Sample_type_ch2",

  "Age", "Age_bin", "Age_detailed",
  "Grade",
  "Node", "Node_bin", "Node_cat",
  "Size", "Size_bin", "Size_cat",

  "Hr", "Er", "Pr", "Her2",

  "Regimen", "Arm_neoadj", "Arm_adj", "Arm_detailed",
  "Response", "Response_clinical", "Response_pathological",
  "Event_dfs", "Time_dfs",

  "Therapy_adjuvant",
  "Therapy_chemo", "Therapy_chemo_name",
  "Therapy_hormone", "Therapy_hormone_name",  #!!!!!!!!!!!!
  "Therapy_radio", "Therapy_surgery", "Therapy_surgery_type",

  "Timepoint",
  "Histology", "Gender", "Ethnicity", "Menopause",
  "Til", "Til_stromal", "Pdl1_stromal", "Pdl1_tumor"
  # Er_score
  # Pr_score
  # Her2_score
  # Her2_score.ihc
  # Her2_score.fish
  # Bcl2_score
  #
  # Mut_p53
  # Ihc_p53
  #
  # Top2a
  # Topo2
  # Topo_score
)


geo_clin <- geo_clin %>% dplyr::select(all_of(nme))
dim(geo_clin) # 5128 51
glimpse(geo_clin)

geo_clin %>%
  dplyr::group_by(Arm_neoadj = !is.na(geo_clin$Arm_neoadj),
                  Arm_adj = !is.na(geo_clin$Arm_adj),
                  dfs = !is.na(geo_clin$Time_dfs),
                  Pcr = !is.na(geo_clin$Response)) %>%
  dplyr::summarise(N=n())
#   Arm_neoadj Arm_adj dfs   Pcr       N
# 1 FALSE      FALSE   FALSE FALSE    19
# 2 FALSE      FALSE   FALSE TRUE    207
# 3 FALSE      FALSE   TRUE  FALSE   151
# 4 FALSE      FALSE   TRUE  TRUE    117
# 5 FALSE      TRUE    FALSE FALSE   108
# 6 FALSE      TRUE    TRUE  FALSE   953 ! adj only
# 7 TRUE       FALSE   FALSE FALSE    72
# 8 TRUE       FALSE   FALSE TRUE   2434 ! neoadj only
# 9 TRUE       FALSE   TRUE  FALSE     9
# 10 TRUE       FALSE   TRUE  TRUE    234
# 11 TRUE       TRUE    FALSE TRUE     11 ! neoadj only
# 12 TRUE       TRUE    TRUE  FALSE   401 ! adj only
# 13 TRUE       TRUE    TRUE  TRUE    306 ! adj + neoadj





# Clean Arm_neoadj and Arm_adj
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>

# 1) In the combined neoadj_adj regimens consider only neoadj regimen.
# 2) Merge Arm_neoadj and Arm_adj to Arm
# 3) Ignore local therapies; surgery and radiation
# 4) Clean ambigous arms
# 5) Clean relevant arms

geo_clin <- geo_clin %>%
  dplyr::mutate(
    # Clean regimen
    Regimen_updated = if_else( Regimen == "neoadj_adj", "neoadj", Regimen),
    # Merge neoadj and adj arms
    Arm = if_else( Regimen_updated == "neoadj", Arm_neoadj, Arm_adj)
    )
glimpse(geo_clin)



table(geo_clin$Regimen)
# adj     neoadj neoadj_adj
# 1113       3078        831
table(geo_clin$Regimen_updated)
# adj neoadj
# 1113   3909

geo_clin$Arm_neoadj %>% unique() %>% length() # 49
geo_clin$Arm_adj %>% unique() %>% length() # 53
geo_clin$Arm %>% unique() %>% length() # 93



# Explore arms >>>>>>>>>>>>
geo_clin$Arm %>% unique()


# Ambigous arm definitions:
#
# [+Hormone] > means hormone therapy if er/pr positive; if er/pr missing set as +Other
# [FEC|Carboplatin] > +Alkalyting_agent(Carboplatin/Cyclophosphamide)+Other
# Other
# FAC|FEC
# Anthracyclin+Taxane[±Antimetabolite±Alkalyting_agent±Trastuzumab] > Anthracyclin+Taxane+other
# Chemo > keep it intact

# [+Hormone] exploring
#
# Note: Hormone therapy assignment based on pam50 subtype may be misleading.!!!!
# Note: Pam50 and er/pr/her2 subtypes are not 100% congruent, hence may lead to wrong assignment.!!!!
# table(geo_clin$Er[str_detect(geo_clin$Arm, "\\[\\+Hormone\\]")]) # all NA
# table(geo_clin$Pr[str_detect(geo_clin$Arm, "\\[\\+Hormone\\]")]) # all NA
# table(geo_clin$Her2[str_detect(geo_clin$Arm, "\\[\\+Hormone\\]")]) # all NA
#
# table(geo_clin$Series_matrix_accession[str_detect(geo_clin$Arm, "\\[\\+Hormone\\]")])
# # GSE20685
# # 202

# Radio therapy exploring
#
# sum(geo_clin$Arm == "Radiotherapy", na.rm = T) # 70
# sum(geo_clin$Arm == "Radio", na.rm = T) # 7
# sum(str_detect(geo_clin$Arm, "Radiotherapy"), na.rm = T) # 70
# sum(str_detect(geo_clin$Arm, "Radio\\+"), na.rm = T) # 100
# sum(str_detect(geo_clin$Arm, "Radio"), na.rm = T) # 177 (7 Radio, 100 Radio+, and 70 Radiotherapy)


geo_clin <- geo_clin %>%
  dplyr::mutate(

    Arm_updated = purrr::map_chr(
      Arm %>% str_trim(),
      ~(case_when(
        # removing local therapy
        .x == "Radiotherapy" ~ "No_systemic_therapy",
        .x == "Radio" ~ "No_systemic_therapy",
        # cleaning ambigous arms
        .x == "No_systemic_therapy[+Hormone]" ~ NA_character_, #"Other",
        .x == "Anthracyclin+Taxane[±Antimetabolite±Alkalyting_agent±Trastuzumab]" ~ "Anthracyclin+Taxane+Other",
        # cleaning relevant arms
        .x == "TFEC" ~ "Paclitaxel_FEC",
        .x == "TFAC" ~ "Paclitaxel_FAC",
        .x == "Tam" ~ "Tamoxifen",
        TRUE ~ .x)
      )
    )  %>%
      # removing local therapy
      str_replace("Radio\\+", "") %>%

      # cleaning ambigous arms
      str_replace("\\[\\+Hormone\\]","+Other") %>%
      str_replace("\\[FEC\\|Carboplatin\\]","Other") %>%
      str_replace("Chemo","Other") %>%

      # cleaning relevant arms (ignore supportive therapy)
      str_replace("CAF","FAC") %>%
      str_replace("Tam_", "Tamoxifen_") %>%
      str_replace("_ZoledronicAcid",""),



    Arm_expanded = purrr::map_chr(
      Arm_updated,
      ~(.x %>%
          str_replace("No_systemic_therapy", "NoSystemicTherapy") %>%
          str_replace("Nab_paclitaxel","Paclitaxel") %>%
          str_replace_all("_","+") %>%
          str_replace("FAC","Fluorouracil+Doxorubicin+Cyclophosphamide") %>%
          str_replace("FEC","Fluorouracil+Epirubicin+Cyclophosphamide") %>%
          str_replace("CMF","Cyclophosphamide+Methotrexate+Fluorouracil") %>%
          str_replace("FAX","Fluorouracil+Doxorubicin+Other") %>%
          str_replace("AC","Doxorubicin+Cyclophosphamide") %>% # if do it first unexpected behaviour happens
          str_replace("EC","Epirubicin+Cyclophosphamide")

      )
    )
  )
geo_clin$Arm_updated %>% unique() %>% sort() # 79
geo_clin$Arm_expanded %>% unique() %>% sort() # 74

drugs <- geo_clin$Arm_expanded %>%
  unique() %>%
  str_split("\\+") %>%
  unlist() %>%
  str_split("\\|") %>%
  unlist() %>%
  unique()
# [1] "Taxane"            "Doxorubicin"       "Cyclophosphamide"
# [4] "Paclitaxel"        "Fluorouracil"      "Docetaxel"
# [7] "Capecitabin"       "Epirubicin"        NA
# [10] "Other"             "Methotrexate"      "Untreated"
# [13] "Ixabepilone"       "Trastuzumab"       "Anastrozole"
# [16] "Capecitabine"      "Tamoxifen"         "Letrozole"
# [19] "Imatinib"          "NoSystemicTherapy" "Pertuzumab"
# [22] "Goserelin"         "Carboplatin"       "Lapatinib"
# [25] "Bevacizumab"       "Anthracyclin"

# write_tsv(x=tibble(Drug_name = drugs), path = str_c(outdir, "drugs.tsv"))
drugs <- read_tsv(file = str_c(outdir, "drugs_cleaned.tsv"))
drugs <- drugs %>%
  dplyr::mutate(
    Drug_class = str_replace_all(Drug_class, " ", "_") %>%
      str_to_lower() %>%
      str_replace("epothilone_b","taxane"), # due to function similar to taxane
    Drug_class_active = if_else(Drug_class == "targetted", Drug_target, Drug_class) %>%
      str_replace_all(" ", "_") %>%
      str_to_lower()
  )

table(drugs$Drug_class)
# alkylating_agent   anthracycline   antimetabolite        targetted
#                2               3                4                9
# taxane
#      4
table(drugs$Drug_class_active)
# alkylating_agent   anthracycline   antimetabolite               er
#                2               3                4                4
# her2           taxane  tyrosine_kinase             vegf
#    3                4                1                1
drugs %>% dplyr::select(Drug_name, Drug_class_active) %>% as.data.frame()
#            Drug_name Drug_class_active
# 1             Taxane            taxane
# 2        Doxorubicin     anthracycline
# 3   Cyclophosphamide  alkylating_agent
# 4         Paclitaxel            taxane
# 5       Fluorouracil    antimetabolite
# 6          Docetaxel            taxane
# 7        Capecitabin    antimetabolite
# 8         Epirubicin     anthracycline
# 9               <NA>              <NA>
#   10             Other              <NA>
#   11      Methotrexate    antimetabolite
# 12         Untreated              <NA>
#   13       Ixabepilone            taxane
# 14       Trastuzumab              her2
# 15       Anastrozole                er
# 16      Capecitabine    antimetabolite
# 17         Tamoxifen                er
# 18         Letrozole                er
# 19          Imatinib   tyrosine_kinase
# 20 NoSystemicTherapy              <NA>
#   21        Pertuzumab              her2
# 22         Goserelin                er
# 23       Carboplatin  alkylating_agent
# 24         Lapatinib              her2
# 25       Bevacizumab              vegf
# 26      Anthracyclin     anthracycline
glimpse(geo_clin)


geo_clin <- geo_clin %>%
  dplyr::mutate(

    Has_anthracycline = if_else(
      str_detect(Arm_expanded, "Doxorubicin") |
        str_detect(Arm_expanded, "Epirubicin") |
        str_detect(Arm_expanded, "Anthracyclin"),
      "yes",
      "no"),
    Has_antimetabolite = if_else(
      str_detect(Arm_expanded, "Fluorouracil") |
        str_detect(Arm_expanded, "Capecitabin") | # also detects Capecitabine
        str_detect(Arm_expanded, "Methotrexate"),
      "yes",
      "no"),
    Has_alkylating_agent = if_else(
      str_detect(Arm_expanded, "Cyclophosphamide") |
        str_detect(Arm_expanded, "Carboplatin"),
      "yes",
      "no"),
    Has_aaa = if_else(str_c(Has_anthracycline, Has_antimetabolite,Has_alkylating_agent) == "yesyesyes",
                      "yes",
                      "no"),
    Has_fac = if_else(
      str_detect(Arm_updated, "FAC"),
      "yes",
      "no"),
    Has_fec = if_else(
      str_detect(Arm_updated, "FEC"),
      "yes",
      "no"),
    Has_cmf = if_else(
      str_detect(Arm_updated, "CMF"),
      "yes",
      "no"),
    Has_taxane = if_else(
      str_detect(Arm_expanded, "Docetaxel") |
        str_detect(Arm_expanded, "Paclitaxel") | # includes Nab-paclitaxel
        str_detect(Arm_expanded, "Ixabepilone") | # epithelone_B, considered as taxane due to mtubule stabilization
        str_detect(Arm_expanded, "Taxane"),
      "yes",
      "no"),
    Has_hormone = if_else(
      str_detect(Arm_expanded, "Anastrozole") |
        str_detect(Arm_expanded, "Tamoxifen") |
        str_detect(Arm_expanded, "Letrozole") |
        str_detect(Arm_expanded, "Goserelin"),
      "yes",
      "no"),
    Has_her2_agent = if_else(
      str_detect(Arm_expanded, "Trastuzumab") |
        str_detect(Arm_expanded, "Pertuzumab") |
        str_detect(Arm_expanded, "Lapatinib"),
      "yes",
      "no"),
    Has_other = if_else(
      str_detect(Arm_expanded, "Imatinib") |
        str_detect(Arm_expanded, "Bevacizumab") |
        str_detect(Arm_expanded, "Other"),
      "yes",
      "no"),
    Has_no_treatment = if_else(
      str_detect(Arm_expanded, "Untreated") |
        str_detect(Arm_expanded, "NoSystemicTherapy"),
      "yes",
      "no"),


    Name_anthracycline = paste(
      if_else(str_detect(Arm_expanded, "Doxorubicin"), "Doxorubicin", ""),
      if_else(str_detect(Arm_expanded, "Epirubicin"), "Epirubicin", ""),
      if_else(str_detect(Arm_expanded, "Anthracyclin"), "", "")
    ),
    Name_antimetabolite = paste(
      if_else(str_detect(Arm_expanded, "Fluorouracil"), "Fluorouracil",""),
      if_else(str_detect(Arm_expanded, "Capecitabin"), "Capecitabine", ""),
      if_else(str_detect(Arm_expanded, "Methotrexate"), "Methotrexate", "")
    ),
    Name_alkylating_agent = paste(
      if_else(str_detect(Arm_expanded, "Cyclophosphamide"), "Cyclophosphamide",""),
      if_else(str_detect(Arm_expanded, "Carboplatin"), "Carboplatin", "")
    ),
    Name_taxane = paste(
      if_else(str_detect(Arm_expanded, "Docetaxel"), "Docetaxel", ""),
      if_else(str_detect(Arm_expanded, "Paclitaxel"), "Paclitaxel", ""),
      if_else(str_detect(Arm_expanded, "Ixabepilone"), "Ixabepilone", ""),
      if_else(str_detect(Arm_expanded, "Taxane"), "", "")
    ),
    Name_hormone = paste(
      if_else(str_detect(Arm_expanded, "Anastrozole"), "Anastrozole", ""),
      if_else(str_detect(Arm_expanded, "Tamoxifen"), "Tamoxifen", ""),
      if_else(str_detect(Arm_expanded, "Letrozole"), "Letrozole", ""),
      if_else(str_detect(Arm_expanded, "Goserelin"), "Goserelin", "")
    ),
    Name_her2_agent = paste(
      if_else(str_detect(Arm_expanded, "Trastuzumab"), "Trastuzumab", ""),
      if_else(str_detect(Arm_expanded, "Pertuzumab"), "Pertuzumab", ""),
      if_else(str_detect(Arm_expanded, "Lapatinib"), "Lapatinib", "")
    ),
    Name_other = paste(
      if_else(str_detect(Arm_expanded, "Imatinib"), "Imatinib", ""),
      if_else(str_detect(Arm_expanded, "Bevacizumab"), "Bevacizumab", ""),
      if_else(str_detect(Arm_expanded, "Other"), "", "")
    )
  ) %>%
  dplyr::mutate_at(
    vars(starts_with("Name_")),
    ~(.x %>%  str_replace_all("NA","") %>%
        str_trim() %>% str_replace_all("  "," ") %>% str_replace_all(" ","+"))
    # "The .vars argument lets you specify columns in the same way that you would
    # specify columns in select(), provided you put that specification inside the function vars()."
    # https://stackoverflow.com/questions/45634566/using-mutate-at-from-dplyr
  )
glimpse(geo_clin)

table(geo_clin$Has_aaa)
table(geo_clin$Has_fac)
table(geo_clin$Has_fec)
table(geo_clin$Has_cmf)


geo_clin <- geo_clin %>%
  dplyr::mutate(
    Arm_anthracycline = paste(Has_anthracycline, Has_taxane, Has_no_treatment) %>%
      purrr::map_chr(~(case_when(.x == "no no yes" ~ "No_systemic_treatment",
                                 .x == "no yes no" ~ "Taxane",
                                 .x == "yes no no" ~ "Anthracycline",
                                 .x == "yes yes no" ~ "Anthracycline+Taxane"))),
    Arm_aaa = paste(Has_aaa, Has_taxane, Has_no_treatment) %>%
      purrr::map_chr(~(case_when(.x == "no no yes" ~ "No_systemic_treatment",
                                 .x == "no yes no" ~ "Taxane",
                                 .x == "yes no no" ~ "AAA",
                                 .x == "yes yes no" ~ "AAA+Taxane"))),
    Arm_fac = paste(Has_fac, Has_taxane, Has_no_treatment) %>%
      purrr::map_chr(~(case_when(.x == "no no yes" ~ "No_systemic_treatment",
                                 .x == "no yes no" ~ "Taxane",
                                 .x == "yes no no" ~ "FAC",
                                 .x == "yes yes no" ~ "FAC+Taxane"))),
    Arm_fec = paste(Has_fec, Has_taxane, Has_no_treatment) %>%
      purrr::map_chr(~(case_when(.x == "no no yes" ~ "No_systemic_treatment",
                                 .x == "no yes no" ~ "Taxane",
                                 .x == "yes no no" ~ "FEC",
                                 .x == "yes yes no" ~ "FEC+Taxane"))),
    Arm_cmf = paste(Has_cmf, Has_taxane, Has_no_treatment) %>%
      purrr::map_chr(~(case_when(.x == "no no yes" ~ "No_systemic_treatment",
                                 .x == "no yes no" ~ "Taxane",
                                 .x == "yes no no" ~ "CMF",
                                 .x == "yes yes no" ~ "CMF+Taxane")))
  )


geo_clin <- geo_clin %>%
  dplyr::mutate(
    Arm_fac = if_else(Arm_updated == "FAC|FEC", NA_character_, Arm_fac),
    Arm_fec = if_else(Arm_updated == "FAC|FEC", NA_character_, Arm_fec),
  )


table(geo_clin$Arm_anthracycline)
# Anthracycline  Anthracycline+Taxane No_systemic_treatment      Taxane
# 995                  2180                   216                   504
table(geo_clin$Arm_aaa)
# AAA            AAA+Taxane No_systemic_treatment                Taxane
# 503                  1187                   216                  1497
table(geo_clin$Arm_fac)
# FAC            FAC+Taxane No_systemic_treatment                Taxane
# 155                   462                   216                  2222
table(geo_clin$Arm_fec)
# FEC            FEC+Taxane No_systemic_treatment                Taxane
# 243                   569                   216                  2115
table(geo_clin$Arm_cmf)
# CMF            CMF+Taxane No_systemic_treatment                Taxane
# 139                   156                   216                  2528

glimpse(geo_clin)
dim(geo_clin) # 5022 79

table(is.na(geo_clin$Arm_updated))
# FALSE  TRUE
# 4367   655

geo_clin %>%
  dplyr::group_by(NA_arm = is.na(Arm_updated), Arm_anthracycline) %>%
  dplyr::summarise(N = n())
#   NA_arm Arm_anthracycline         N
# 1 FALSE  Anthracycline           995
# 2 FALSE  Anthracycline+Taxane   2180
# 3 FALSE  No_systemic_treatment   216
# 4 FALSE  Taxane                  504
# 5 FALSE  NA                      472
# 6 TRUE   NA                      655
# 995+2180+216+504 = 3895



# Update geo_clin with subtyping and cellluar estimates
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


identical(geo_clin$Sample_geo_accession, names(geo_expr)[-1]) # FALSE
dim(geo_clin) #[1] 5022   79
dim(geo_expr) #[1] 9184 3737 # not all dataset is integrated


geo_clin <- geo_clin %>%
  dplyr::mutate(
    Has_expr_pooled = if_else(Sample_geo_accession %in% names(geo_expr)[-1],
                              "yes", "no")
  )



# Do ihc subtyping for all samples ingeo_clin which has er(hr) and her2
# Do molecular subtyping for samples in geo_expr


# Subtyping IHC based >>>>>>>>>>>>>>>>>>>>>

table(geo_clin$Hr)
# neg pos
# 117 228
table(is.na(geo_clin$Hr))
# FALSE  TRUE
# 345  4677

geo_clin <- geo_clin %>%
  dplyr::mutate(
    ErPr = paste(Er, Pr), # temp column
    Hr = if_else(
      is.na(Hr),
       case_when(ErPr == "neg NA" ~ "neg",
                 ErPr == "neg neg" ~ "neg",
                 ErPr == "neg pos" ~ "pos",
                 ErPr == "pos NA" ~ "pos",
                 ErPr == "pos neg" ~ "pos",
                 ErPr == "pos pos" ~ "pos"),
      Hr),
    Subtype_ihc = purrr::map_chr(
      str_c(Hr, Her2),
      ~case_when(.x == "negneg" ~ "TN",
                 .x == "negpos" ~ "HER2",
                 .x == "posneg" ~ "HR",
                 .x == "pospos" ~ "HER2")
    )
  ) %>%
  dplyr::select(-ErPr)

table(geo_clin$Hr)
# neg  pos
# 2017 2446
table(is.na(geo_clin$Hr))
# FALSE  TRUE
# 4463   559

table(geo_clin$Subtype_ihc)
# HER2   HR   TN
# 1031 1603 1094
table(is.na(geo_clin$Subtype_ihc))
# FALSE  TRUE
# 3728  1294



# Subtyping based on expression bimodality >>>>>>>>>>>>>>>>>>>>>>

identical(geo_clin$Sample_geo_accession, names(geo_expr)[-1]) # FALSE
dim(geo_clin) #[1] 5022   79
dim(geo_expr) #[1] 9184 3737 # not all dataset is integrated

table(names(geo_expr)[-1] %in% geo_clin$Sample_geo_accession)
# TRUE
# 3736

# Note:
# Generate Subtype_expr from geo_expr and merge to geo_clin.

"ncbi_2099" %in% geo_expr$Ncbi_gene_id # ESR1 present
"ncbi_5241" %in% geo_expr$Ncbi_gene_id # PGR present
"ncbi_2064" %in% geo_expr$Ncbi_gene_id # ERBB2 present


dim(geo_expr) # 9184 3737

xdata <- geo_expr[, -1] %>% t()
colnames(xdata) <- geo_expr$Ncbi_gene_id
rownames(xdata) <- names(geo_expr)[-1]
dim(xdata) # 3736 9184


esr1 <- genefu::bimod(
  x = tibble(probe ="ncbi_2099", EntrezGene.ID = 2099,coefficient = 1),
  data = xdata,
  annot = tibble(EntrezGene.ID = str_replace(rownames(xdata), "ncbi_","")),
  do.mapping = FALSE, do.scale = TRUE, verbose = TRUE
)
# Warning message:
#   In genefu::bimod(x = tibble(probe = "ncbi_2099", EntrezGene.ID = 2099,  :
# 7 is the most likely number of Gaussians!

pgr <- genefu::bimod(
  x = tibble(probe = "ncbi_5241", EntrezGene.ID = 5241, coefficient = 1),
  data = xdata,
  annot = tibble(EntrezGene.ID = str_replace(rownames(xdata), "ncbi_","")),
  do.mapping = FALSE, do.scale = TRUE, verbose = TRUE
)
# Warning message:
# In genefu::bimod(x = tibble(probe = "ncbi_5241", EntrezGene.ID = 5241,  :
#                               7 is the most likely number of Gaussians!

erbb2 <- genefu::bimod(
  x = tibble(probe = "ncbi_2064", EntrezGene.ID = 2064, coefficient = 1),
  data = xdata,
  annot = tibble(EntrezGene.ID = str_replace(rownames(xdata), "ncbi_","")),
  do.mapping = FALSE, do.scale = TRUE, verbose = TRUE
)
# Warning message:
# In genefu::bimod(x = tibble(probe = "ncbi_2064", EntrezGene.ID = 2064,  :
#                               6 is the most likely number of Gaussians!


table(esr1$status) # very wrong !!!!!!!!!
# 0    1
# 78  3658
table(pgr$status) # very likely !!!!!!!!!
# 0    1
# 2916 820
table(erbb2$status) # very likely !!!!!!!!!!
#  0    1
# 3198  538


# biomodality checking
x <- geo_expr %>% filter(Ncbi_gene_id == "ncbi_2099") %>% as.numeric()
plot(density(x[-1]), main = "ESR1") # bimod

x <- geo_expr %>% filter(Ncbi_gene_id == "ncbi_5241") %>% as.numeric()
plot(density(x[-1]), main = "PGR") # trimod ?

x <- geo_expr %>% filter(Ncbi_gene_id == "ncbi_2064") %>% as.numeric()
plot(density(x[-1]), main = "ERBB2") # bimod


xx <- tibble(
  Sample_geo_accession = rownames(xdata),
  Esr1_expr = if_else(esr1$status == 1, "pos", "neg"),
  Pgr_expr = if_else(pgr$status == 1, "pos", "neg"),
  Her2_expr = if_else(erbb2$status == 1, "pos", "neg"),

  Hr_expr = purrr::map_chr(
    str_c(Esr1_expr, Pgr_expr, sep = "_"),
    function(x){
      case_when(
        x == "neg_neg" ~ "neg",
        x == "neg_pos" ~ "pos",
        x == "pos_neg" ~ "pos",
        x == "pos_pos" ~ "pos",
        TRUE ~ NA_character_
      )
    }
  ),

  Subtype_expr = purrr::map_chr(
    str_c(Hr_expr,Her2_expr, sep = "_"),
    function(x){
      case_when(
        x == "neg_neg" ~ "TN",
        x == "neg_pos" ~ "HER2",
        x == "pos_neg" ~ "HR",
        x == "pos_pos" ~ "HER2",
        TRUE ~ NA_character_
      )
    }
  )
)

dim(xx) # 3736    6
table(xx$Sample_geo_accession %in% geo_clin$Sample_geo_accession)
# TRUE
# 3736

geo_clin <- geo_clin %>%
  left_join(xx, by = "Sample_geo_accession")




# Subtyping pam50 based >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pam50 <- list()
pam50[["centroid"]] <- read_tsv(file = "data/pam50/PAM50/bioclassifier_R/pam50_centroids.txt")
pam50[["annot"]] <- read_tsv(file = "data/pam50/PAM50/bioclassifier_R/pam50_annotation.txt")

pam50 <- pam50$centroid %>%
  dplyr::rename(pcrID = "X1") %>%
  dplyr::left_join(pam50$annot, by = "pcrID") %>%
  dplyr::mutate(EntrezGene = str_c("ncbi_", EntrezGene)) %>%
  dplyr::filter(EntrezGene %in% geo_expr$Ncbi_gene_id) %>%
  dplyr::rename(Ncbi_gene_id = "EntrezGene")
nrow(pam50) # 36


# expression
xx <-  pam50 %>%
  dplyr::select(Ncbi_gene_id) %>%
  dplyr::left_join(geo_expr, by = "Ncbi_gene_id")

xx[, -1] <- xx[, -1] %>%
  apply(1,function(x){x-median(x)}) %>% # gene centering, the resut is transposed
  t() %>% as.data.frame()


# pam50 subtyping
xx_corr <- cor(pam50[, c("Basal", "Her2", "LumA", "LumB", "Normal")],
               xx[,-1],
               method = "pearson")
# Median centering is valid only with pearson corelation.
# Foe spearman correlation, median centering has no effect.
head(xx_corr[,1:10])

pam50_subtype <- purrr::map_dfr(
  colnames(xx_corr),
  function(sample, xx_corr, subtype){

    corr <- xx_corr[, sample]

    list(
      Sample_geo_accession = sample,
      Subtype_pam50 = subtype[which(corr == max(corr))],
      Pam50_centroid_corr = corr[which(corr == max(corr))]
    )

  },
  xx_corr,
  subtype = rownames(xx_corr)
)


pam50_subtype <- pam50_subtype %>%
  dplyr::mutate(
    Subtype_pam50_2 = purrr::map_chr(
      Subtype_pam50,
      function(x){
        case_when(
          x == "Basal" ~ "TN",
          x == "Her2" ~ "HER2",
          x == "LumA" ~ "HR",
          x == "LumB" ~ "HR",
          TRUE ~ NA_character_
        )
      }
    )

  )

dim(geo_clin) # 5022   85
dim(pam50_subtype) # 3736 4
table(pam50_subtype$Sample_geo_accession %in% geo_clin$Sample_geo_accession)
# TRUE
# 3736

geo_clin <- geo_clin %>%
  left_join(pam50_subtype, by = "Sample_geo_accession")

dim(geo_clin) # 5022   88

table(geo_clin$Subtype_pam50)
# Basal   Her2   LumA   LumB Normal
# 1079    625   1085    572    375
table(geo_clin$Subtype_pam50_2)
# HER2   HR   TN
# 625 1657 1079


# Subtyping Kappa agreement >>>>>>>>>>>>>>>>>>>>>>>>>>>

# Kappa,κ ranges from 0 to 1, with 0 indicating no relation and 1 indicating
# a perfect concordance. Typically qualitative descriptions are associated
# with intervals [κ ≤ 0.20, slight agreement; 0.20 < κ ≤ 0.40, fair agreement;
# 0.40 < κ ≤ 0.60, moderate agreement, 0.60 < κ ≤ 0.80, substantial agreement;
# and 0.80 < κ ≤ 0.99, almost perfect agreement
# See Benjamins Three gene model paper


# IHC vs Expression(Pam50 and bimodal) >>>>>>>>>>>>>>>>>>>>

rating <- geo_clin %>%
  dplyr::select(Er, Pr, Hr, Her2, Esr1_expr, Pgr_expr, Her2_expr,
                Subtype_expr, Subtype_ihc, Subtype_pam50, Subtype_pam50_2)
dim(rating) # 5022 11
rating <- na.omit(rating)
dim(rating) # 2325 11

table(rating$Subtype_ihc)
# HER2   HR   TN
# 653 1061  611
table(rating$Subtype_expr)
# HER2   HR   TN
# 408 1856   61
table(rating$Subtype_pam50_2)
# HER2   HR   TN
# 446 1176  703

# Hormone receptor binary status agreement
xx <- irr::kappa2(ratings = rating %>% na.omit() %>% dplyr::select(Er, Esr1_expr))
xx$value # 0.04214222, poor
xx <- irr::kappa2(ratings = rating %>% na.omit() %>% dplyr::select(Pr, Pgr_expr))
xx$value # 0.296443, fair
xx <- irr::kappa2(ratings = rating %>% na.omit() %>% dplyr::select(Her2, Her2_expr))
xx$value # 0.5227322, moderate


# Subtype  agreement
xx <- irr::kappa2(ratings = rating %>% na.omit() %>% dplyr::select(Subtype_ihc, Subtype_expr))
xx$value #  0.2645115, fair
xx <- irr::kappa2(ratings = rating %>% na.omit() %>% dplyr::select(Subtype_ihc, Subtype_pam50_2))
xx$value #  0.6543398, substantial
xx <- irr::kappa2(ratings = rating %>% na.omit() %>% dplyr::select(Subtype_expr, Subtype_pam50_2))
xx$value # 0.2764711, fair


# Note:
# For the downstream analysis, use IHC and pam50 subtypes





# Integrate sample procrument and preservation methods to geo_clin >>>>>>>>>>>>>

glimpse(geo_summary)

table(geo_summary$Series_matrix_accession %in% geo_clin$Series_matrix_accession)
# FALSE  TRUE
# 137    38


geo_clin <- geo_clin %>%
  left_join(geo_summary %>%
              dplyr::select(Series_matrix_accession,
                            Series_sample_procurement,
                            Series_archive_method),
            by = "Series_matrix_accession")

geo_clin %>%
  dplyr::group_by(Has_expr_pooled,Series_archive_method) %>%
  dplyr::summarise(N = n())
#   Has_expr_pooled Series_archive_method     N
# 1 no              FFPE                    155
# 2 no              Frozen                 1088
# 3 no              NA                       43
# 4 yes             FFPE                    156
# 5 yes             Frozen                 1524
# 6 yes             RNAlater                826
# 7 yes             RNAlater|Frozen         508
# 8 yes             NA                      722



# # Save RData !!!!!!!!!!!!!!!!!
# #
# save(geo_clin, file = str_c(outdir, "geo_clin.RData"))
# load(str_c(outdir, "geo_clin.RData"))


# geo_clean_clin (similar to geo_clean_expr list)
# Useful for perdataset analysis(meta)
#
geo_clean_clin <- list()
for(inme in names(geo_clean_expr)){

  geo_clean_clin[[inme]] <- geo_clin %>%
    dplyr::filter(Series_matrix_accession == inme)

  print(
    paste(
      inme,
      nrow(geo_clean_clin[[inme]]) == (ncol(geo_clean_expr[[inme]])-1),
      identical(geo_clean_clin[[inme]]$Sample_geo_accession, names(geo_clean_expr[[inme]])[-1])
    )
  )
  # all TRUE; geo_clean_clin and geo_clean_expr are congruent
}


purrr::map2_lgl(
  geo_clean_clin,
  geo_clean_expr,
  ~(identical(.x$Sample_geo_accession, names(.y)[-1]))
) %>% table()
# TRUE all identical
# 38


# # Save RData !!!!!!!!!!!!!!!!!
# #
# save(geo_clean_clin, file = str_c(outdir, "geo_clean_clin.RData"))
# load(str_c(outdir, "geo_clean_clin.RData"))

#
# ==============================================================================


