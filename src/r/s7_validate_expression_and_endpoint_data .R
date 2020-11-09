# 7_validate_expression_and_endpoint_data.R

# 1) log+scaling+centering of pooled data
# 2) RLE plots
# 3) ER/HER2 IHC status prediction ability (AUC)
# 4) PAM50 vs Subtype-IHC Kappa agreement
# 5) ER - immune response association (save plots)
# 6) Heterogenity test per arm

# load("results/main/geo_expr.RData")
# load("results/main/geo_expr_meta.RData")
# load("results/main/geo_clin.RData")


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

table(geo_clin$Subtype_ihc)
# HER2   HR   TN
# 1031 1603 1094
table(is.na(geo_clin$Subtype_ihc))
# FALSE  TRUE
# 3728  1294



geo_clin %>%
  dplyr::filter(Regimen_updated == "neoadj") %>%
  dplyr::group_by(Response, Response_clinical, Response_pathological, Event_dfs) %>%
  dplyr::summarise(N = n()) %>%
  as.data.frame()
#    Response Response_clinical  Response_pathological Event_dfs    N
# 1      npCR             CR+kl                   <NA>        NA   10
# 2      npCR              MPG1                   <NA>        NA    4
# 3      npCR              MPG2                   <NA>        NA   24
# 4      npCR              MPG3                   <NA>        NA   20
# 5      npCR                NE                   <NA>        NA    1
# 6      npCR       No_response                   <NA>        NA   38
# 7      npCR              NOCR                   <NA>        NA  137
# 8      npCR     non_downstage                   <NA>        NA   35
# 9      npCR                PD                   <NA>        NA    3
# 10     npCR                PR                   <NA>        NA   57
# 11     npCR                SD                   <NA>        NA    5
# 12     npCR              <NA>                      0         0   46
# 13     npCR              <NA>                      0         1   15
# 14     npCR              <NA>                      0        NA    9
# 15     npCR              <NA>            NO RESPONSE        NA    6
# 16     npCR              <NA>                   npCR         0  288
# 17     npCR              <NA>                   npCR         1  148
# 18     npCR              <NA>                   npCR        NA 1214
# 19     npCR              <NA>       PARTIAL RESPONSE        NA   86
# 20     npCR              <NA>                     RD        NA  255
# 21      pCR                CR                   <NA>        NA   41
# 22      pCR         downstage                   <NA>        NA    9
# 23      pCR              MPG4                   <NA>        NA    7
# 24      pCR              MPG5                   <NA>        NA    6
# 25      pCR          Response                   <NA>        NA   23
# 26      pCR              <NA>                      1         0   90
# 27      pCR              <NA>                      1         1   11
# 28      pCR              <NA>                      1        NA    2
# 29      pCR              <NA>      COMPLETE RESPONSE        NA   20
# 30      pCR              <NA> NEAR COMPLETE RESPONSE        NA   10
# 31      pCR              <NA>                    pCR         0  107
# 32      pCR              <NA>                    pCR         1   14
# 33      pCR              <NA>                    pCR        NA  568
# 34     <NA>              <NA>                   npCR         0  290
# 35     <NA>              <NA>                   npCR         1   99
# 36     <NA>              <NA>                    pCR         0   92
# 37     <NA>              <NA>                    pCR         1    7
# 38     <NA>              <NA>                   <NA>         0   18
# 39     <NA>              <NA>                   <NA>         1   11
# 40     <NA>              <NA>                   <NA>        NA   83

# Set if Respose as Response_pathological if Response == NA

geo_clin <- geo_clin %>%
  dplyr::mutate(Response = if_else(is.na(Response),
                                    Response_pathological,
                                    Response))

geo_clin %>%
  dplyr::filter(Regimen_updated == "adj") %>%
  # dplyr::filter(Selected_for_analysis & Regimen_updated == "adj") %>%
  dplyr::group_by(Response, Response_clinical, Response_pathological, Event_dfs) %>%
  dplyr::summarise(N = n())
#   Response Response_clinical Response_pathological Event_dfs     N
# 1 NA       NA                NA                            0   793
# 2 NA       NA                NA                            1   309
# 3 NA       NA                NA                           NA    11








# 6) Heterogenity test per arm
# ==============================================================================

glimpse(geo_clin %>% dplyr::select(Arm_chemo, Arm_her2, Arm_hormone, Arm_other))

xx <- geo_clin %>%
  dplyr::group_by(Arm_chemo, Arm_her2, Arm_hormone, Arm_other, Regimen_updated, Subtype_ihc, Series_matrix_accession) %>%
  # , Series_matrix_accession
  dplyr::summarise(N = n(), N_dataset = unique(Series_matrix_accession) %>% length()) %>%
  as.data.frame()

# Selected xx arms with minimum 2 datasets to test for heterogeneity
# Arm and dataset selection summary
# Arm selection criteria:
# 1) No.of datasets per arm per subtype >= 2,
# 2) No.of samples per arm per subtype >= 50, and
# 3) No.of samples per dataset per arm per subtype >= 10



# Arm only summary
xx <- geo_clin %>%
  dplyr::group_by(Arm_chemo, Arm_her2, Arm_hormone, Arm_other) %>%
  dplyr::summarise(N_dataset = n_distinct(Series_matrix_accession), N_sample = n()) %>%
  as.data.frame()
write_tsv(x =xx, path = str_c(outdir,"heterogeneity_arm_selection_level1.tsv"))

# Arm + Subtype summary
xx <- geo_clin %>%
  dplyr::group_by(Arm_chemo, Arm_her2, Arm_hormone, Arm_other, Subtype_ihc) %>%
  dplyr::summarise(N_dataset = n_distinct(Series_matrix_accession), N_sample = n()) %>%
  as.data.frame()
write_tsv(x =xx, path = str_c(outdir,"heterogeneity_arm_selection_level2.tsv"))

# Arm + Subtype + Regimen summary
xx <- geo_clin %>%
  dplyr::group_by(Arm_chemo, Arm_her2, Arm_hormone, Arm_other, Subtype_ihc, Regimen_updated) %>%
  dplyr::summarise(N_dataset = n_distinct(Series_matrix_accession), N_sample = n()) %>%
  as.data.frame()
write_tsv(x =xx, path = str_c(outdir,"heterogeneity_arm_selection_level3.tsv"))


# Arm + Subtype + Regimen + Dataset summary
xx <- geo_clin %>%
  dplyr::group_by(Arm_chemo, Arm_her2, Arm_hormone, Arm_other, Subtype_ihc, Regimen_updated, Series_matrix_accession) %>%
  dplyr::summarise(N_dataset = n_distinct(Series_matrix_accession), N_sample = n()) %>%
  as.data.frame()
write_tsv(x =xx, path = str_c(outdir,"heterogeneity_arm_selection_level4.tsv"))



# Select samples for analysis
geo_clin <- geo_clin %>%
  dplyr::mutate(
    Selected_for_analysis =

      # Arm: Tamoxifen only
      (Arm_chemo == "000+noTaxane" & Arm_her2 == "No_her2_agent" &
         Arm_hormone == "Tamoxifen" & Arm_other == "No_other_therapy" &
         Subtype_ihc == "HR" & Series_matrix_accession %in% c("GSE16391",
                                                              "GSE19615",
                                                              "GSE45255",
                                                              "GSE69031")) |

      # Arm: A0A+Taxane
      (Arm_chemo == "A0A+Taxane" & Arm_her2 == "No_her2_agent" &
         Arm_hormone == "No_hormone_therapy" & Arm_other == "No_other_therapy" &
         Subtype_ihc == "HR" & Series_matrix_accession %in% c("GSE114403",
                                                              "GSE21974",
                                                              "GSE22226_GPL1708",
                                                              "GSE22226_GPL4133",
                                                              "GSE25066",
                                                              "GSE32603",
                                                              "GSE41998")) |
      (Arm_chemo == "A0A+Taxane" & Arm_her2 == "No_her2_agent" &
         Arm_hormone == "No_hormone_therapy" & Arm_other == "No_other_therapy" &
         Subtype_ihc == "HER2" & Series_matrix_accession %in% c("GSE22226_GPL1708",
                                                                "GSE32603",
                                                                "GSE41998"
         )) |
      (Arm_chemo == "A0A+Taxane" & Arm_her2 == "No_her2_agent" &
         Arm_hormone == "No_hormone_therapy" & Arm_other == "No_other_therapy" &
         Subtype_ihc == "TN" & Series_matrix_accession %in% c("GSE143222",
                                                              "GSE22226_GPL1708",
                                                              "GSE25066",
                                                              "GSE32603",
                                                              "GSE41998")) |

    # Arm: AAA+noTaxane
      (Arm_chemo == "AAA+noTaxane" & Arm_her2 == "No_her2_agent" &
         Arm_hormone == "No_hormone_therapy" & Arm_other == "No_other_therapy" &
         Subtype_ihc == "HR" & Series_matrix_accession %in% c("GSE20271",
                                                              "GSE22093")) |
      (Arm_chemo == "AAA+noTaxane" & Arm_her2 == "No_her2_agent" &
         Arm_hormone == "No_hormone_therapy" & Arm_other == "No_other_therapy" &
         Subtype_ihc == "TN" & Series_matrix_accession %in% c("GSE20271",
                                                              "GSE22093")) |

      # Arm: AAA+Taxane
      (Arm_chemo == "AAA+Taxane" & Arm_her2 == "No_her2_agent" &
         Arm_hormone == "No_hormone_therapy" & Arm_other == "No_other_therapy" &
         Subtype_ihc == "HR" & Series_matrix_accession %in% c("GSE20194",
                                                              "GSE20271",
                                                              "GSE23988",
                                                              "GSE25066",
                                                              "GSE32646",
                                                              "GSE42822",
                                                              "GSE50948")) |
      (Arm_chemo == "AAA+Taxane" & Arm_her2 == "No_her2_agent" &
         Arm_hormone == "No_hormone_therapy" & Arm_other == "No_other_therapy" &
         Subtype_ihc == "HER2" & Series_matrix_accession %in% c("GSE20194",
                                                                "GSE20271",
                                                                "GSE32646",
                                                                "GSE42822",
                                                                "GSE50948")) |
      (Arm_chemo == "AAA+Taxane" & Arm_her2 == "No_her2_agent" &
         Arm_hormone == "No_hormone_therapy" & Arm_other == "No_other_therapy" &
         Subtype_ihc == "TN" & Series_matrix_accession %in% c("GSE20194",
                                                              "GSE20271",
                                                              "GSE23988",
                                                              "GSE25066",
                                                              "GSE32646",
                                                              "GSE42822",
                                                              "GSE50948")) |

      # Arm: AAA+Taxane+Trastuzumab
      (Arm_chemo == "AAA+Taxane" & Arm_her2 == "Trastuzumab" &
         Arm_hormone == "No_hormone_therapy" & Arm_other == "No_other_therapy" &
         Subtype_ihc == "HER2" & Series_matrix_accession %in% c("GSE42822",
                                                                "GSE50948"))
  )

table(geo_clin$Selected_for_analysis)
# FALSE  TRUE
# 2695  1967

xx <- geo_clin %>%
  dplyr::filter(Selected_for_analysis) %>%
  dplyr::group_by(Arm_chemo, Arm_her2, Arm_hormone, Arm_other, Subtype_ihc,
                  Regimen_updated, Series_matrix_accession, is.na(Response)) %>%
  dplyr::summarise(N_dataset = n_distinct(Series_matrix_accession), N_sample = n()) %>%
  as.data.frame()
sum(xx$N_sample) # 1967
# Identical to heterogeneity_arm_selection_level1-2-3-4.ods(sheet: selected_arms)
# NA's doesn't breake the minimum sample criteria of datasets(Minimum 10 samples),
# Hence the no.of arms and dataset dataset remains intact, but some sample will
#   be discarded due to "NA" response data.

head2(geo_clin) # 5022   82


# >>>
# Clinical retyping
geo_clin <- geo_clin %>%
  dplyr::mutate_all(~(hablar::retype(.x)))
# >>>


geo_clin %>%
  dplyr::filter(Selected_for_analysis & Regimen_updated == "neoadj") %>%
  dplyr::group_by(Response, Response_clinical, Response_pathological) %>%
  dplyr::summarise(N = n())
#   Response Response_clinical Response_pathological     N
# 1 npCR     NA                npCR                   1167
# 2 npCR     NA                RD                      245
# 3 pCR      NA                pCR                     455
# 4 NA       NA                NA                       26

geo_clin %>%
  dplyr::filter(Selected_for_analysis & Regimen_updated == "adj") %>%
  dplyr::group_by(Response, Response_clinical, Response_pathological, Event_dfs) %>%
  dplyr::summarise(N = n())
# Response Response_clinical Response_pathological Event_dfs     N
# 1 NA       NA                NA                            0    62
# 2 NA       NA                NA                            1    10
# 3 NA       NA                NA                           NA     2





xx <- geo_clin %>%
  dplyr::filter(Selected_for_analysis & Regimen_updated == "neoadj") %>%
  dplyr::group_by(Arm_chemo, Arm_her2, Arm_hormone, Arm_other, Subtype_ihc, Regimen_updated, Series_matrix_accession) %>%
  dplyr::summarise(N_dataset = n_distinct(Series_matrix_accession),
                   N_sample = n(),
                   OR = ,
                   SE_OR) %>%
  as.data.frame()
sum(xx$N_sample) # 1967

#   Arm_chemo  Arm_her2  Arm_hormone Arm_other Subtype_ihc Regimen_updated Series_matrix_a… N_dataset N_sample
# 1 A0A+Taxane No_her2_… No_hormone… No_other… HER2        neoadj          GSE22226_GPL1708         1       27
# 2 A0A+Taxane No_her2_… No_hormone… No_other… HER2        neoadj          GSE32603                 1       28
# 3 A0A+Taxane No_her2_… No_hormone… No_other… HER2        neoadj          GSE41998                 1       27
# 4 A0A+Taxane No_her2_… No_hormone… No_other… HR          neoadj          GSE114403                1       20
# 5 A0A+Taxane No_her2_… No_hormone… No_other… HR          neoadj          GSE21974                 1       15
# 6 A0A+Taxane No_her2_… No_hormone… No_other… HR          neoadj          GSE22226_GPL1708         1       40
# 7 A0A+Taxane No_her2_… No_hormone… No_other… HR          neoadj          GSE22226_GPL4133         1       13
# 8 A0A+Taxane No_her2_… No_hormone… No_other… HR          neoadj          GSE25066                 1       36
# 9 A0A+Taxane No_her2_… No_hormone… No_other… HR          neoadj          GSE32603                 1       62
# 10 A0A+Taxane No_her2_… No_hormone… No_other… HR          neoadj          GSE41998                 1      106

xx <-  geo_clin %>%
  dplyr::filter(Selected_for_analysis & Regimen_updated == "neoadj" &
                  Arm_chemo == "A0A+Taxane" & Subtype_ihc == "HER2" &
                  Series_matrix_accession == "GSE22226_GPL1708") %>%
  dplyr::mutate(Response = if_else(Response == "npCR", 0, 1))

# 1264
fit1 <- glm(
  formula = Response ~ 1,
  data = xx,
  family = "binomial"
)
fit1$coefficients
# (Intercept)
# -0.3101549
fit_sum <- summary(fit1)


xx <-  geo_clin %>%
  dplyr::filter(Selected_for_analysis & Regimen_updated == "neoadj" &
                  Arm_chemo == "A0A+Taxane") %>%
  dplyr::mutate(Response = if_else(Response == "npCR", 0, 1))

# 1264
fit1 <- glm(
  formula = Response ~ Subtype_ihc,
  data = xx,
  family = "binomial"
)
fit_sum <- summary(fit1)
fit1$coefficients
# (Intercept) Subtype_ihcHR Subtype_ihcTN
# -0.55961579   -1.62216445   -0.01092907
