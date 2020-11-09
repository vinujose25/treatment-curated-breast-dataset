# s5_manual_clinical_data_integration.R

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Integrate clinical data avaialable from geo.RData
# Treatment regimen clssification
# Then split out the integrated clinical data and
# populate the clinical slot of geo_tidy object.


# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Prepare data
# 2. Generate series matrix summary and integrate series summary into it
# (This is to update integrated clinical data with sample procrument and
# preservation info later)
# 3. Integrate clinical data with relevant characteristics
# 4. Dissect treatment regiment and reclassify regimen in a consistant way
# 5. Update integrated clinical data with sample preocrument and preservation details
# 6. Split out integrated clinical data and populate clinical slot of geo_tidy



# 1. Prepare data
# ==============================================================================

load("results/main/geo.RData")
length(geo) #38 series matrices

load("results/main/geo_series_summary.RData")
glimpse(geo_series_summary)

#
# ==============================================================================




# 2. Generate series matrix summary and integrate series summary into it
# (This is to update integrated clinical data with sample procrument and
# preservation info later)
# ==============================================================================


# Genarate GEO series matrix summary
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Series-matrix summary (extracted from geo list; contains only) >>>>>>>
geo_series_matrix_summary <- tibble(
  Series_accession = str_split_fixed(names(geo), "_", 2)[, 1],
  Series_matrix_accession = names(geo),
  Series_matrix_probe_size = purrr::map_int(geo,~nrow(.x$expression)),
  Series_matrix_gene_size = purrr::map_int(geo_maxvar_annot,~nrow(.x)),
  Series_matrix_platform = purrr::map_chr(geo,~(unique(.x$clinical$Sample_platform_id))),
  Series_matrix_regimen =  purrr::map_chr(geo,~(unique(.x$clinical$Regimen))),
  Series_matrix_sample_size = purrr::map_int(geo,~nrow(.x$clinical)),
  # Samples with treatment info
  Series_matrix_regimen_size = purrr::map2_int(
    Series_matrix_regimen, geo,
    ~(if_else(str_detect(.x, "neoadj"),
              sum(!is.na(.y$clinical$Arm_neoadj)),
              sum(!is.na(.y$clinical$Arm_adj))))),
  Series_matrix_has_neoadj = purrr::map_lgl(
    geo, ~(any(str_detect(names(.x$clinical), "_neoadj")))),
  Series_matrix_has_adj = purrr::map_lgl(
    geo, ~(any(str_detect(names(.x$clinical), "_adj")))),
  Series_matrix_selected = NA,
  Series_matrix_comment = NA
)

glimpse(geo_series_matrix_summary)


# Checking congruence between Series_matrix_regimen, Series_matrix_has_neoadj,
#  and Series_matrix_has_adj
geo_series_matrix_summary %>%
  select(Series_matrix_accession, Series_matrix_regimen,
         Series_matrix_has_neoadj, Series_matrix_has_adj) %>%
  as.data.frame()
# expectation: adj regimen has no neoadj arm and vice-versa
# neoadj_adj has both neoadj and adj arm.
# Note that GSE16391 adj dataset seems failed expectations, but this is due to
# an original column wth pattern "_neoadj"



# Update geo_series_matrix_summary with geo_series_summary
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

geo_series_matrix_summary <- full_join(
  geo_series_summary, geo_series_matrix_summary,
  by = "Series_accession") %>%
  dplyr::mutate(
    Series_matrix_selected = if_else(is.na(Series_matrix_accession), "no", "yes")
  )
glimpse(geo_series_matrix_summary)

table(geo_series_matrix_summary$Series_matrix_selected)
# no yes
# 137  38
table(geo_series_matrix_summary$Series_regimen) # irrelevant
# adj  cell line     neoadj neoadj/adj
# 15          7         42          2
table(geo_series_matrix_summary$Series_matrix_regimen)
# adj     neoadj neoadj_adj
# 7         27          4

geo_series_matrix_summary %>%
  dplyr::group_by(Series_matrix_selected, Series_matrix_regimen) %>%
  dplyr::summarise(Sum = sum(Series_matrix_sample_size),
                   Sum_regimen = sum(Series_matrix_regimen_size))
#   Series_matrix_selected Series_matrix_regimen   Sum Sum_regimen
# 1 no                     NA                       NA          NA
# 2 yes                    adj                    1113         954
# 3 yes                    neoadj                 3078        2743
# 4 yes                    neoadj_adj              831         724


table(geo_series_matrix_summary$Series_sample_procurement)
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
table(geo_series_matrix_summary$Series_archive_method)
# ffpe          frozen        RNAlater RNAlater/Frozen
# 5              28               6               3


geo_series_matrix_summary <- geo_series_matrix_summary %>%
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


geo_series_matrix_summary %>%
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

geo_series_matrix_summary %>%
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


# Save geo_series_matrix_summary !!!!!!!!!!!!!!!!!!!
#
# save(geo_series_matrix_summary, file = str_c(outdir,"geo_series_matrix_summary.RData"))
# load(str_c(outdir,"geo_series_matrix_summary.RData"))

#
# ==============================================================================




# 3. Integrate clinical data with relevant characteristics
# ==============================================================================


# Consolidate clinical data of selected sereis matrices
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# adj platform
table(geo_series_matrix_summary %>%
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
table(geo_series_matrix_summary %>%
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


# For neoadjuvant, endpoint is Response
# For adjuvant, endpoint id Time_dfs, Event_dfs
nme <- geo_series_matrix_summary %>%
  dplyr::filter(Series_matrix_selected == "yes") %>%
  dplyr::select(Series_matrix_accession) %>%
  tibble::deframe()
# n = 38


geo_clin <- NULL
for(inme in nme){
  print(inme)
  xx <- geo[[inme]]$clinical
  xx$Series_matrix_accession <- inme
  names(xx) <- str_to_sentence(names(xx))
  geo_clin <- bind_rows(geo_clin, xx)
}
dim(geo_clin) # 5022 429

write_tsv(x = tibble(Column_name = sort(names(geo_clin))),
          path = str_c(outdir, "geo_clin_column_name.tsv"))


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


  # "Therapy_adjuvant",     # Redundant info, already encoded in Arm* !!!!!!!!
  # "Therapy_chemo", "Therapy_chemo_name",
  # "Therapy_hormone", "Therapy_hormone_name",  #!!!!!!!!!!!!
  # "Therapy_radio", "Therapy_surgery", "Therapy_surgery_type",

  "Timepoint",
  "Histology", "Gender", "Ethnicity", "Menopause",
  "Til", "Til_stromal", "Pdl1_stromal", "Pdl1_tumor"

  # Irrelevant info !!!!!!!!!!!!!!!!

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
dim(geo_clin) # 5022 43
glimpse(geo_clin)

#
# ==============================================================================




# 4. Dissect treatment regiment and reclassify regimen in a consistant way
# ==============================================================================

# Clean Arm_neoadj and Arm_adj
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>

# 1) In the combined neoadj_adj regimens consider only neoadj regimen.
# 2) Merge Arm_neoadj and Arm_adj to Arm
# 3) Ignore local therapies; surgery and radiation
# 4) Clean ambigous arms
# 5) Clean relevant arms

geo_clin <- geo_clin %>%
  dplyr::mutate(
    # 1) Clean regimen
    Regimen_updated = if_else( Regimen == "neoadj_adj", "neoadj", Regimen),
    # 2) Merge neoadj and adj arms
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



# Explore arms
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



# Ignore local therapies(surgery and radiation) and clean ambigous arms
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


# Write out drugs.tsv for manual annoattion !!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# write_tsv(x=tibble(Drug_name = drugs), path = str_c(outdir, "drugs.tsv"))


# Loading annoatated and cleaned drugs.tsv !!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
drugs <- read_tsv(file = str_c(outdir, "drugs_cleaned.tsv"))

drugs <- drugs %>%
  dplyr::mutate(
    Drug_class = str_replace_all(Drug_class, " ", "_") %>%
      str_to_lower() %>%
      str_replace("epothilone_b","taxane"), # due to function similar to taxane
    # !!!!!!!!!! epothilone_b  as taxane

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



# Dissecting the original treatment regimen
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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

    Has_taxane = if_else(
      str_detect(Arm_expanded, "Docetaxel") |
        str_detect(Arm_expanded, "Paclitaxel") | # includes Nab-paclitaxel
        str_detect(Arm_expanded, "Ixabepilone") | # epithelone_B, considered as taxane due to mtubule stabilization
        str_detect(Arm_expanded, "Taxane"),
      "yes",
      "no"),

    Has_chemo = paste(
      Has_anthracycline, Has_antimetabolite,
      Has_alkylating_agent,
      # Has_aaa, Has_fac, Has_fec, Has_cmf, # redundant
      Has_taxane
    ) %>%
      purrr::map_chr(
        ~(
          ifelse(str_detect(.x, "yes"),
                 "yes", # any yes is 100% chemo
                 ifelse(str_detect(.x, "NA"),
                        NA_character_, # no & NA is not 100% no, hence NA
                        "no")
          )
        )
      ),

    # new
    Class_chemo = paste(Has_anthracycline, Has_antimetabolite, Has_alkylating_agent) %>%
      purrr::map_chr(~(case_when(.x == "no no no" ~ "000",
                                 .x == "yes no no" ~ "A00",
                                 .x == "no yes no" ~ "0A0",
                                 .x == "no no yes" ~ "00A",
                                 .x == "yes yes no" ~ "AA0",
                                 .x == "yes no yes" ~ "A0A",
                                 .x == "no yes yes" ~ "0AA",
                                 .x == "yes yes yes" ~ "AAA"))),
    Class_chemo = str_c(Class_chemo, Has_taxane) %>%
      str_replace("yes", "+Taxane") %>%
      str_replace("no", "+noTaxane"),


    Has_aaa = if_else(str_c(Has_anthracycline, Has_antimetabolite, Has_alkylating_agent) == "yesyesyes",
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

    # Note !!!!!!!!!!!!
    # Without is.na(Has_*) check, NAs won't propogate in the following code.

    Name_anthracycline = if_else(
      is.na(Has_anthracycline),
      NA_character_,
      paste(
        if_else(str_detect(Arm_expanded, "Doxorubicin"), "Doxorubicin", ""),
        if_else(str_detect(Arm_expanded, "Epirubicin"), "Epirubicin", ""),
        if_else(str_detect(Arm_expanded, "Anthracyclin"), "UnknownAnthracyclin", "")
      )
    ),
    Name_antimetabolite = if_else(
      is.na(Has_antimetabolite),
      NA_character_,
      paste(
        if_else(str_detect(Arm_expanded, "Fluorouracil"), "Fluorouracil",""),
        if_else(str_detect(Arm_expanded, "Capecitabin"), "Capecitabine", ""),
        if_else(str_detect(Arm_expanded, "Methotrexate"), "Methotrexate", "")
      )
    ),
    Name_alkylating_agent = if_else(
      is.na(Has_alkylating_agent),
      NA_character_,
      paste(
        if_else(str_detect(Arm_expanded, "Cyclophosphamide"), "Cyclophosphamide",""),
        if_else(str_detect(Arm_expanded, "Carboplatin"), "Carboplatin", "")
      )
    ),
    Name_taxane = if_else(
      is.na(Has_taxane),
      NA_character_,
      paste(
        if_else(str_detect(Arm_expanded, "Docetaxel"), "Docetaxel", ""),
        if_else(str_detect(Arm_expanded, "Paclitaxel"), "Paclitaxel", ""),
        if_else(str_detect(Arm_expanded, "Ixabepilone"), "Ixabepilone", ""),
        if_else(str_detect(Arm_expanded, "Taxane"), "UnknownTaxane", "")
      )
    ),
    Name_hormone = if_else(
      is.na(Has_hormone),
      NA_character_,
      paste(
        if_else(str_detect(Arm_expanded, "Anastrozole"), "Anastrozole", ""),
        if_else(str_detect(Arm_expanded, "Tamoxifen"), "Tamoxifen", ""),
        if_else(str_detect(Arm_expanded, "Letrozole"), "Letrozole", ""),
        if_else(str_detect(Arm_expanded, "Goserelin"), "Goserelin", "")
      )
    ),
    Name_her2_agent =  if_else(
      is.na(Has_her2_agent),
      NA_character_,
      paste(
        if_else(str_detect(Arm_expanded, "Trastuzumab"), "Trastuzumab", ""),
        if_else(str_detect(Arm_expanded, "Pertuzumab"), "Pertuzumab", ""),
        if_else(str_detect(Arm_expanded, "Lapatinib"), "Lapatinib", "")
      )
    ),
    Name_other = if_else(
      is.na(Has_other),
      NA_character_,
      paste(
        if_else(str_detect(Arm_expanded, "Imatinib"), "Imatinib", ""),
        if_else(str_detect(Arm_expanded, "Bevacizumab"), "Bevacizumab", ""),
        if_else(str_detect(Arm_expanded, "Other"), "UnknownTherapy", "")
      )
    )



  ) %>% # end of mutate
  dplyr::mutate_at(
    vars(starts_with("Name_")),
    ~(.x %>%
        str_replace_all("NA","") %>%
        str_trim() %>%
        str_replace_all("  "," ") %>%
        str_replace_all(" ","+"))
    # "The .vars argument lets you specify columns in the same way that you would
    # specify columns in select(), provided you put that specification inside the function vars()."
    # https://stackoverflow.com/questions/45634566/using-mutate-at-from-dplyr
  )
glimpse(geo_clin)


geo_clin$Has_anthracycline %>% is.na() %>% table()
# FALSE  TRUE
# 4367   655
geo_clin$Name_anthracycline %>% is.na() %>% table()
# FALSE  TRUE
# 4367   655


table(geo_clin$Has_chemo)
# no  yes
# 550 3817
table(geo_clin$Class_chemo)
# 000+noTaxane   000+Taxane   00A+Taxane   0A0+Taxane 0AA+noTaxane A00+noTaxane
# 550          219          110          175          138          174
# A00+Taxane A0A+noTaxane   A0A+Taxane   AA0+Taxane AAA+noTaxane   AAA+Taxane
# 167          318          824            2          503         1187
table(geo_clin$Has_hormone)
# no  yes
# 3969  398
table(geo_clin$Has_her2_agent)
# no  yes
# 3799  568
table(geo_clin$Has_other)
# no  yes
# 3966  401
table(geo_clin$Has_no_treatment)
# no  yes
# 4151  216


table(geo_clin$Has_aaa)
# no  yes
# 2677 1690
table(geo_clin$Has_fac)
# no  yes
# 3647  720
table(geo_clin$Has_fec)
# no  yes
# 3452  915
table(geo_clin$Has_cmf)
# no  yes
# 4072  295



# Reclassifying treatment regimen
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

geo_clin <- geo_clin %>%
  dplyr::mutate(

    # In the literature focus is on Anthracycline +/- taxane, hence Taxane
    #   is included in every arm defition.
    # Also the following arm classification ignores treatments given outside
    # its definition(except taxanes); such as hormone, her2, other.

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
                                 .x == "yes yes no" ~ "CMF+Taxane"))),

    Arm_chemo = paste(Class_chemo, Has_no_treatment) %>%
      str_replace(" no", "") %>%
      purrr::map_chr(~(case_when(.x == "000+noTaxane yes" ~ "No_systemic_treatment",
                                 .x == "NA NA" ~ NA_character_,
                                 TRUE ~ .x))),
    Arm_her2 = if_else(Has_her2_agent == "no",
                       "No_her2_agent",
                       Name_her2_agent),
    Arm_hormone = if_else(Has_hormone == "no",
                          "No_hormone_therapy",
                          Name_hormone),
    Arm_other = if_else(Has_other == "no",
                        "No_other_therapy",
                        Name_other),

    # integrating has_no_treatment
    Arm_her2 = if_else(Has_no_treatment == "yes",
                       "No_systemic_treatment",
                       Arm_her2),
    Arm_hormone = if_else(Has_no_treatment == "yes",
                          "No_systemic_treatment",
                          Arm_hormone),
    Arm_other = if_else(Has_no_treatment == "yes",
                        "No_systemic_treatment",
                        Arm_other)

  )

geo_clin <- geo_clin %>%
  dplyr::mutate(
    # Discarding ambigous arm
    Arm_fac = if_else(Arm_updated == "FAC|FEC", NA_character_, Arm_fac),
    Arm_fec = if_else(Arm_updated == "FAC|FEC", NA_character_, Arm_fec),
  )


table(is.na(geo_clin$Arm_updated))
# FALSE  TRUE
# 4367   655

table(geo_clin$Arm_chemo %>% is.na())
# FALSE  TRUE
# 4367   655
table(geo_clin$Arm_her2 %>% is.na())
# FALSE  TRUE
# 4367   655
table(geo_clin$Arm_hormone %>% is.na())
# FALSE  TRUE
# 4367   655
table(geo_clin$Arm_other %>% is.na())
# FALSE  TRUE
# 4367   655


table(geo_clin$Arm_her2)
# Lapatinib          No_her2_agent  No_systemic_treatment
# 60                   3583                    216
# Trastuzumab  Trastuzumab+Lapatinib Trastuzumab+Pertuzumab
# 251                     84                    173
table(geo_clin$Arm_hormone)
# Anastrozole Anastrozole+Tamoxifen             Letrozole
# 5                     3                    33
# No_hormone_therapy No_systemic_treatment             Tamoxifen
# 3753                   216                   356
# Tamoxifen+Goserelin
# 1
table(geo_clin$Arm_other)
# Bevacizumab              Imatinib      No_other_therapy
# 23                     1                  3750
# No_systemic_treatment        UnknownTherapy
# 216                   377
table(geo_clin$Arm_chemo)
# 000+noTaxane            000+Taxane            00A+Taxane
# 334                   219                   110
# 0A0+Taxane          0AA+noTaxane          A00+noTaxane
# 175                   138                   174
# A00+Taxane          A0A+noTaxane            A0A+Taxane
# 167                   318                   824
# AA0+Taxane          AAA+noTaxane            AAA+Taxane
# 2                   503                  1187
# No_systemic_treatment
# 216
table(geo_clin$Class_chemo)
# 000+noTaxane   000+Taxane   00A+Taxane   0A0+Taxane 0AA+noTaxane A00+noTaxane
# 550          219          110          175          138          174
# A00+Taxane A0A+noTaxane   A0A+Taxane   AA0+Taxane AAA+noTaxane   AAA+Taxane
# 167          318          824            2          503         1187


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
dim(geo_clin) # 5022 77


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

#
# ==============================================================================




# 5. Update integrated clinical data with sample preocrument and preservation details
# ==============================================================================

# Integrate sample procrument and preservation methods to geo_clin
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

glimpse(geo_series_matrix_summary)

table(geo_series_matrix_summary$Series_matrix_accession %in% geo_clin$Series_matrix_accession)
# FALSE  TRUE
# 137    38


geo_clin <- geo_clin %>%
  left_join(geo_series_matrix_summary %>%
              dplyr::select(Series_matrix_accession,
                            Series_sample_procurement,
                            Series_archive_method),
            by = "Series_matrix_accession") %>%
  dplyr::mutate(
    Has_expr_pooled = if_else(Sample_geo_accession %in% geo_expr_meta$Sample_geo_accession,
                              "yes", "no")
  )

length(geo_clin$Series_matrix_accession %>% unique()) # 38
length(geo_expr_meta$Series_matrix_accession %>% unique()) # 29


geo_clin %>%
  dplyr::group_by(Has_expr_pooled, Series_archive_method) %>%
  dplyr::summarise(N = n())
#   Has_expr_pooled Series_archive_method     N
# 1 FALSE           FFPE                    155
# 2 FALSE           Frozen                 1088
# 3 FALSE           NA                       43
# 4 TRUE            FFPE                    156
# 5 TRUE            Frozen                 1524
# 6 TRUE            RNAlater                826
# 7 TRUE            RNAlater|Frozen         508
# 8 TRUE            NA                      722



head2(geo_clin)
#   Series_matrix_accessi… Sample_title Sample_geo_accession Sample_channel_cou… Sample_platform_…
# 1 GSE25066               1002         GSM615096            1                   GPL96
# 2 GSE25066               1005         GSM615097            1                   GPL96
# 3 GSE25066               1009         GSM615098            1                   GPL96
# 4 GSE25066               1011         GSM615099            1                   GPL96
# 5 GSE25066               1016         GSM615100            1                   GPL96
# $dim
# [1] 5022   80


# # Save geo_clin.RData !!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# #
# save(geo_clin, file = str_c(outdir, "geo_clin.RData"))
# load(str_c(outdir, "geo_clin.RData"))

#
# ==============================================================================




# 6. Split out integrated clinical data and populate clinical slot of geo_tidy
# ==============================================================================

# Integrate geo_clin to geo_tidy
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# geo_tidy (similar to geo), but contained cleaned max-var collapsed expression
#   and curated/annotated clinical data.
# Useful for per dataset analysis(meta)


for(inme in names(geo_tidy)){

  geo_tidy[[inme]]$clinical <- geo_clin %>%
    dplyr::filter(Series_matrix_accession == inme)

  print(
    paste(
      inme,
      nrow(geo_tidy[[inme]]$clinical) == (ncol(geo_tidy[[inme]]$expression)-1),
      identical(geo_tidy[[inme]]$clinical$Sample_geo_accession,
                names(geo_tidy[[inme]]$expression)[-1])
    )
  )
  # all TRUE; geo_tidy[[]]$clinical and geo_tidy[[]]$expression are congruent
}


purrr::map_lgl(
  geo_tidy,
  ~(identical(.x$clinical$Sample_geo_accession, names(.x$expression)[-1]))
) %>% table()
# TRUE all identical
# 38


# # Save RData !!!!!!!!!!!!!!!!!
# #
# save(geo_tidy, file = str_c(outdir, "geo_tidy.RData"))
# load(str_c(outdir, "geo_tidy.RData"))

#
# ==============================================================================



