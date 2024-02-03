# s7_figures_tables_data.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# 1. Write out dataset effect corrected data in tsv format.
# 2. Write out matching clinical data with qc metrics and clinico-patho variables.
# 3. Generate RLE plots, AUC table, Kappa table related to dataset effect comparison
# 4. Generate patient characteristics and drug-class summary.
# 5. Note that the study flow chart and the integrated script structure + study flow chart
#    figures are created in "LibreOffice Draw"


# Script structure
# >>>>>>>>>>>>>>>>
#
# 1. Prepare data
# 2. Write out data as tsv files
# 3. Generate RLE plots
# 4. Compute Cohen's Kappa and generate table plot
# 5. Compute ROC AUC and generate table plot
# 6. Summarize patient characteristics and generate table plot
# 7. Summarize drug classes and generate table plot
# 8. Summarize non relevant patient characteristics
# 9. Generate a table containing GEO search with series details,
#    - selected/rejected flag, reason for rejection.
# 10. Compare curatedBreastData package with present study



# 1. Prepare data
# ==============================================================================

# To generate qc related to the reliability of expression data matrix
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

load(str_c(out_data,"geo_expr_list.RData"))
load(str_c(out_data,"geo_clin_list.RData"))

names(geo_expr_list)[1] # Original
names(geo_expr_list)[1] = "Non-corrected"

names(geo_clin_list)[1] # Original
names(geo_clin_list)[1] = "Non-corrected"


names(geo_expr_list)[4] # Madquantile
names(geo_expr_list)[4] = "MAD-quantile"

names(geo_clin_list)[4] # Madquantile
names(geo_clin_list)[4] = "MAD-quantile"



# For comaprison with curatedBreastData R package.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

data("clinicalData")
planey <- clinicalData$clinicalTable %>% as_tibble()
planey_annot <- clinicalData$clinicalVarDef %>% as_tibble()
rm(clinicalData)



# To generate a table with original serach result and result of dataset selection
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

load(str_c(out_data, "geo_series_summary.RData"))
load(str_c(out_data, "geo_series_matrix_summary.RData"))


#
# ==============================================================================



# 2. Write out data as tsv files
# ==============================================================================

purrr::map(
  names(geo_expr_list),
  function(nme, lst){
    x <- lst[[nme]]
    x[,-1] <- round(x[,-1], digits = 3)
    nme <- str_replace(nme,"-","") # Non-corrected to Noncorrected
    write_delim(x = x, path = str_c(out_data,"geo_expr_", nme, ".tsv"), delim = "\t")
    TRUE
  },
  geo_expr_list
)


purrr::map(
  names(geo_clin_list),
  function(nme, lst){
    x <- lst[[nme]]
    nme <- str_replace(nme,"-","") # Non-corrected to Noncorrected
    write_delim(x = x, path = str_c(out_data,"geo_clin_", nme, ".tsv"), delim = "\t")
    TRUE
  },
  lst = geo_clin_list
)


rm(geo_expr_list)
gc()

#
# ==============================================================================



# 3. Genrate RLE plots
# ==============================================================================

xx <- purrr::map(
  names(geo_clin_list),
  function(nme, geo_clin_list){
    bind_cols(
      tibble(Batch_correction_methode = nme),
      geo_clin_list[[nme]] %>%
        dplyr::select(Series_matrix_accession, Sample_geo_accession,
                      "Rle_min", "Rle_first_qu", "Rle_median", "Rle_mean",
                      "Rle_third_qu", "Rle_max")
    )
  },
  geo_clin_list
) %>%
  bind_rows()

p <- xx %>%
  tidyr::gather("Rle_class", "Rle_value", "Rle_first_qu","Rle_median", "Rle_third_qu") %>%
  dplyr::mutate(Series_matrix_accession = factor(Series_matrix_accession,
                                                 levels = unique(Series_matrix_accession)),
                Sample_geo_accession = factor(Sample_geo_accession,
                                              levels = unique(Sample_geo_accession)),
                Batch_correction_methode = Batch_correction_methode %>%
                  str_replace("MAD-quantile", "MAD-quantile*"), # to highlight filtering of 2 samples
                Batch_correction_methode = factor(Batch_correction_methode, levels = unique(Batch_correction_methode))) %>%
  ggplot(aes(x = Sample_geo_accession, y = Rle_value, group = Rle_class)) +
  geom_line(aes(color = Series_matrix_accession)) +
  guides(color = "none") +
  labs(x = "Samples", y = "RLE") +
  facet_wrap(facets = ~ Batch_correction_methode, scales = "free_y", ncol = 1) +
  theme(legend.position = "bottom",
        axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


pdf2(file = str_c(out_figures, "figure_rle.pdf"), width = 7.5, height = 4.5)
print(p)
dev.off()

#
# ==============================================================================



# 4. Compute Cohen's Kappa and generate table plot
# ==============================================================================


# Cohen's kappa
# >>>>>>>>>>>>>

# Kappa,κ ranges from 0 to 1, with 0 indicating no relation and 1 indicating
# a perfect concordance. Typically qualitative descriptions are associated
# with intervals [κ ≤ 0.20, slight agreement; 0.20 < κ ≤ 0.40, fair agreement;
# 0.40 < κ ≤ 0.60, moderate agreement, 0.60 < κ ≤ 0.80, substantial agreement;
# and 0.80 < κ ≤ 0.99, almost perfect agreement
# See Benjamins Three gene model paper

table(geo_clin_list$`Non-corrected`$Subtype_ihc) # Reference
# HER2   HR   TN
# 756  1232  850
table(geo_clin_list$`Non-corrected`$Subtype_pam50)
# Basal         Her2         LumA         LumB       Normal Unclassified
# 1089          639         1140          492          341           35
table(geo_clin_list$`Non-corrected`$Subtype_pam50_2) # Test subtype
# HER2   HR   TN
# 639 1632 1089


xx <- purrr::map_dfr(
  names(geo_clin_list),
  function(nme, xx){
    x = xx[[nme]]
    list(method = nme,
         kappa = irr::kappa2(ratings = x %>%
                               dplyr::select(Subtype_ihc, Subtype_pam50_2))$value
    )
  },
  geo_clin_list
)
#   method        kappa
# 1 Non-corrected 0.621
# 2 Combat        0.618
# 3 Quantile      0.634
# 4 MAD-quantile  0.628  !!! filtered 2 samples with MAD = 0


write_delim(x = xx, path = str_c(out_tables,"table_kappa.tsv"), delim = "\t")


p <- xx %>%
  dplyr::mutate(
    y = method,
    kappa = round(kappa, digits = 3)
  ) %>%
  tidyr::gather("key", "value", -y) %>%
  dplyr::mutate(
    key = key %>%
      str_to_title() %>%
      factor(levels = c("Method","Kappa")),
    y = factor(y, levels = y %>% unique() %>% rev()),

    ) %>%
  ggplot(aes(x = key, y = y, label = value)) +
  geom_text(size = 3.5, hjust = 0.5) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black"),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        panel.grid = element_blank()
  )


pdf2(file = str_c(out_tables, "table_kappa.pdf"), height = 1.5, width = 2)
print(p)
dev.off()

#
# ==============================================================================



# 5. Compute ROC AUC and generate table plot
# ==============================================================================

# AUC
# >>>>

xx <- purrr::map(
  geo_clin_list,
  function(x){
    list(er_gene = pROC::roc(formula = Er  ~ Er_gene, data = x),
         her2_gene = pROC::roc(formula = Her2 ~ Her2_gene, data = x))
  }
)



xx <- purrr::map_dfr(
  names(xx),
  function(nme, xx){
    x = xx[[nme]]
    list(method = nme,
         # er_sig = auc(x$er_sig) %>% as.numeric(),
         # her2_sig = auc(x$her2_sig) %>% as.numeric(),
         er_gene = auc(x$er_gene) %>% as.numeric(),
         her2_gene = auc(x$her2_gene) %>% as.numeric())
  }, xx)
#   method          er_gene her2_gene
# 1 Non-corrected   0.791     0.606
# 2 Combat          0.880     0.729
# 3 Quantile        0.895     0.776
# 4 Madquantile     0.892     0.770 !!! filtered 2 samples with MAD = 0


write_delim(x = xx, path = str_c(out_tables,"table_auc.tsv"), delim = "\t")



p <- xx %>%
  dplyr::mutate(
    y = method,
    er_gene = round(er_gene, digits = 3),
    her2_gene = round(her2_gene, digits = 3)
  ) %>%
  tidyr::gather("key", "value", -y) %>%
  dplyr::mutate(
    key = key %>%
      str_replace("method", "Method") %>%
      str_replace("er_gene", "ER-gene\n(AUC)") %>%
      str_replace("her2_gene", "HER2-gene\n(AUC)") %>%
      factor(levels = c("Method", "ER-gene\n(AUC)", "HER2-gene\n(AUC)")),
    y = factor(y, levels = y %>% unique() %>% rev()),

  ) %>%
  ggplot(aes(x = key, y = y, label = value)) +
  geom_text(size = 3.5, hjust = 0.5) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black"),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        panel.grid = element_blank()
  )


pdf2(file = str_c(out_tables, "table_auc.pdf"), height = 1.5, width = 3)
print(p)
dev.off()

#
# ==============================================================================



# 6. Summarize patient characteristics and generate table plot
# ==============================================================================

clin <- geo_clin_list$`Non-corrected`
dim(clin) # 3736 92


xx <- hablar::retype(clin)
dim(xx) # 3736 92


# Exploring relevant characterisitics
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# NA percentage
xx$Age %>% is.na() %>% sum()/3736
# [1] 0.2438437
xx$Age_bin %>% is.na() %>%  sum()/3736
# [1] 0.9561028
xx$Grade %>% is.na() %>% sum()/3736
# [1] 0.388651
xx$Node %>% is.na() %>% sum()/3736
# [1] 0.9245182
xx$Node_bin %>% is.na() %>% sum()/3736
# [1] 0.3642934
xx$Node_cat %>% is.na() %>% sum()/3736
# [1] 0.5230193
xx$Size %>% is.na() %>% sum()/3736
# [1] 0.731531
xx$Size_bin %>% is.na() %>% sum()/3736
# [1] 0.2002141
xx$Size_cat %>% is.na() %>% sum()/3736
# [1] 0.2149358

xx$Er %>% is.na() %>% sum()/3736
# [1] 0.1442719
xx$Pr %>% is.na() %>% sum()/3736
# [1] 0.2312634
xx$Hr %>% is.na() %>% sum()/3736
# [1] 0.1442719
xx$Her2 %>% is.na() %>% sum()/3736
# [1] 0.2149358


xx$Subtype_ihc %>% is.na() %>% sum()/3736
# 0.240364
xx$Subtype_pam50 %>% is.na() %>% sum()/3736
# 0

xx$Response_pathological %>% is.na() %>% sum()/3736
# [1] 0.2917559
xx$Response_clinical %>% is.na() %>% sum()/3736
# [1] 0.9352248
xx$Response %>% is.na() %>% sum()/3736
# [1] 0.2269807

xx$Event_dfs %>% is.na() %>% sum()/3736
# [1] 0.6292827
xx$Time_dfs %>% is.na() %>% sum()/3736
# [1] 0.6563169

xx$Arm_chemo %>% is.na() %>% sum()/3736
# [1] 0.1199143
xx$Arm_her2 %>% is.na() %>% sum()/3736
# [1] 0.1199143
xx$Arm_hormone %>% is.na() %>% sum()/3736
# [1] 0.1199143
xx$Arm_other %>% is.na() %>% sum()/3736
# [1] 0.1199143


xx$Arm_chemo %>% table() %>% as_tibble()
# 1 000+noTaxane            214
# 2 000+Taxane               34
# 3 00A+Taxane              110
# 4 0A0+Taxane              175
# 5 0AA+noTaxane             86
# 6 A00+noTaxane            151
# 7 A00+Taxane              167
# 8 A0A+noTaxane            136
# 9 A0A+Taxane              478
# 10 AA0+Taxane                2
# 11 AAA+noTaxane            503
# 12 AAA+Taxane             1187
# 13 No_systemic_treatment    45
xx$Arm_her2 %>% table() %>% as_tibble()
# 1 Lapatinib                60
# 2 No_her2_agent          2859
# 3 No_systemic_treatment    45
# 4 Trastuzumab             240
# 5 Trastuzumab+Lapatinib    84
xx$Arm_hormone %>% table() %>% as_tibble()
# 1 Anastrozole               5
# 2 Anastrozole+Tamoxifen     3
# 3 Letrozole                33
# 4 No_hormone_therapy     2979
# 5 No_systemic_treatment    45
# 6 Tamoxifen               222
# 7 Tamoxifen+Goserelin       1
xx$Arm_other %>% table() %>% as_tibble()
# 1 Imatinib                  1
# 2 No_other_therapy       3040
# 3 No_systemic_treatment    45
# 4 UnknownTherapy          202


xx %>%
  group_by(Regimen_updated) %>%
  summarise(N_dataset = unique(Series_matrix_accession) %>% length(),
            N_samples = n())
# Regimen_updated N_dataset N_samples
# adj                     5       754
# neoadj                 24      2982



# Get patient summary for neoadj and adj seperatly
# and display them side by side
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pt_sum_adj <- get_patient_summary(
  xx1 = clin %>% filter(Regimen_updated == "adj")
)

pt_sum_neoadj <- get_patient_summary(
  xx1 = clin %>% filter(Regimen_updated == "neoadj")
)

# Introducing "key" variable to position characterisitcs correctly in the plot
pt_sum_adj <- pt_sum_adj %>% dplyr::mutate(key = str_c(variable, cat))
pt_sum_neoadj <- pt_sum_neoadj %>% dplyr::mutate(key = str_c(variable, cat))


# Not all key (e.g. Arm classes) are present in neoadj/adj data seperatley.
# Hence first create a list of comprehensive list of keys, variables and value.
# Then merge individual neoadj/adj summary to it.
k <- full_join(
  pt_sum_neoadj %>% select(key, variable, cat),
  pt_sum_adj %>% select(key, variable, cat),
  by = "key"
) %>%
  dplyr::mutate(
    key = factor(
      key,
      levels = c(
        "Age<=50FALSE",                     "Age<=50TRUE",
        "Age<=50NA",                        "GradeG1",
        "GradeG1/2",                        "GradeG2",
        "GradeG3",                          "GradeNA",
        "Nodeneg",                          "Nodepos",
        "NodeNA",                           "SizeT0",
        "SizeT1",                           "SizeT2",
        "SizeT3",                           "SizeT4",
        "SizeNA",                           "ERneg",
        "ERpos",                            "ERNA",
        "PRneg",                            "PRpos",
        "PRNA",                             "HRneg",
        "HRpos",                            "HRNA",
        "HER2neg",                          "HER2pos",
        "HER2NA",                           "Subtype IHCTN",
        "Subtype IHCHER2",                  "Subtype IHCHR",
        "Subtype IHCNA",                    "Subtype PAM50Basal",
        "Subtype PAM50Her2",                "Subtype PAM50LumB",
        "Subtype PAM50LumA",                "Subtype PAM50Normal",
        "Subtype PAM50Unclassified",        "Subtype PAM50NA",
        "pCRnpCR",                          "pCRpCR",
        "pCRNA",                            "DFSno-event",
        "DFSevent",                         "DFSNA",

        "Arm_chemo000+noTaxane",            "Arm_chemo000+Taxane",
        "Arm_chemo00A+Taxane",              "Arm_chemo0A0+Taxane",
        "Arm_chemoA00+noTaxane",            "Arm_chemoA00+Taxane",
        "Arm_chemoA0A+noTaxane",            "Arm_chemoA0A+Taxane",
        "Arm_chemoAA0+Taxane",              "Arm_chemo0AA+noTaxane",
        "Arm_chemoAAA+noTaxane",            "Arm_chemoAAA+Taxane",
        "Arm_chemoNo_systemic_treatment",   "Arm_chemoNA",

        "Arm_her2Lapatinib",                "Arm_her2Trastuzumab",
        "Arm_her2Trastuzumab+Lapatinib",    "Arm_her2No_her2_agent",
        "Arm_her2No_systemic_treatment",    "Arm_her2NA",

        "Arm_hormoneAnastrozole",           "Arm_hormoneLetrozole",
        "Arm_hormoneTamoxifen",             "Arm_hormoneAnastrozole+Tamoxifen",
        "Arm_hormoneTamoxifen+Goserelin",   "Arm_hormoneNo_hormone_therapy",
        "Arm_hormoneNo_systemic_treatment", "Arm_hormoneNA",

        "Arm_otherImatinib",                "Arm_otherUnknownTherapy",
        "Arm_otherNo_other_therapy",        "Arm_otherNo_systemic_treatment",
        "Arm_otherNA"
      )
    )
  ) %>%
  dplyr::arrange(key)


# Cleaning k (comprehensive list of keys, variables and values)
k <- k %>%
  dplyr::mutate(
    variable.x = ifelse(is.na(variable.x), variable.y, variable.x),
    cat.x = ifelse(is.na(cat.x), cat.y, cat.x)
  ) %>%
  dplyr::rename(
    variable = "variable.x",
    cat = "cat.x"
  ) %>%
  dplyr::select(key, variable, cat)



# Merging neoadj and adj summary to k
pt_sum <- full_join(k,
                    pt_sum_neoadj %>%
                      dplyr::select(key, n) %>%
                      dplyr::rename(neoadj = "n"),
                    by = "key")
pt_sum <- full_join(pt_sum,
                    pt_sum_adj %>%
                      dplyr::select(key, n) %>%
                      dplyr::rename(adj = "n"),
                    by = "key")

pt_sum <- pt_sum %>%
  dplyr::mutate(
    # setting NA counts/frequencies as zero
    neoadj = ifelse(is.na(neoadj), 0, neoadj),
    adj = ifelse(is.na(adj), 0, adj)
  )

glimpse(pt_sum)


# Writing out patient summary !!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write_csv(x = pt_sum, path = str_c(out_tables, "table_patient_summary.csv"))


# Remove repition of variable names in the table for reomving clutter
idx <- purrr::map(
  unique(pt_sum$variable),
  function(x,id){
    which(x == id)[-1] # to keep the first occurance of variable
  },
  id=pt_sum$variable
) %>% unlist()

pt_sum$variable[idx] = ""
pt_sum %>% as.data.frame()


# Table plot
p <- pt_sum %>%
  dplyr::rename("value" = "cat",
                "neoadj-count" = "neoadj",
                "adj-count" = "adj") %>%
  dplyr::rename_all(~str_to_title(.x)) %>%
  dplyr::mutate(
    Value = Value %>%
      str_replace("treatment", "treat.") %>%
      str_replace("therapy", "treat.") %>%
      str_replace("UnknownTherapy", "Unknown_treat.") %>%
      str_replace("Trastuzumab", "Tra.") %>%
      str_replace("Tamoxifen", "Tam."),
    Variable = Variable %>%
      str_replace("Arm_hormone", "Arm_hr") %>%
      str_replace("pCR", "Response")
  ) %>%
  tidyr::gather(key = "X", value = "Label", -1) %>%
  dplyr::mutate(
    X = factor(
      X,
      levels = c("Variable", "Value", "Neoadj-Count", "Adj-Count")
    ),
    Key = factor(
      Key, levels = Key %>% unique() %>% rev()
    ),
    Key_class = ifelse(str_detect(Key, "Arm"), "Treatment", "Clinical-Pathological")
  ) %>%
  ggplot(aes(x = X, y = Key, label = Label)) +
  geom_text(size = 3.5, hjust = -Inf) +
  facet_wrap(facets = ~Key_class, nrow = 1, scales = "free") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black"),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        panel.grid = element_blank()
  )


pdf2(file = str_c(out_tables, "table_patient_summary.pdf"), height = 8, width = 7.5)
print(p)
dev.off()

#
# ==============================================================================



# 7. Summarize drug classes and generate table plot
# ==============================================================================

drug_sum <- read_tsv(file = str_c(out_data, "drugs_cleaned.tsv")) %>%
  na.omit()

drug_sum <- drug_sum %>%
  dplyr::arrange(Drug_class, Drug_target, Drug_name) %>%
  dplyr::select(Drug_class, Drug_target, Drug_name, Drug_function) %>%
  dplyr::mutate(Key = str_c(Drug_class, ".", Drug_name),
                Key = factor(Key, levels = Key),
                Drug_function = ifelse(Key == "Targetted.Imatinib",
                                       "Bcr-Abl pathway",
                                       Drug_function),
                Drug_target = Drug_target %>%
                  str_replace("Er", "ER signalling") %>%
                  str_replace("Her2", "HER2 signalling") %>%
                  str_replace("VEGF", "VEGF signalling"),
                Drug_target = ifelse(
                  Drug_name == "Lapatinib",
                  "HER2/EGFR signalling",
                  Drug_target
                )
  )

write_delim(x = drug_sum, path = str_c(out_tables,"table_drug_class_summary.tsv"), delim = "\t")


# TO removing replication of drug class in the table plot (removes clutter)
idx <- purrr::map(
  unique(drug_sum$Drug_class),
  function(x,id){
    which(x == id)[-1]
  },
  id=drug_sum$Drug_class
) %>% unlist()

drug_sum$Drug_class[idx] = ""
drug_sum%>% as.data.frame()


# Table plot
p <- drug_sum %>%
  tidyr::gather(key = "X", value = "Label", -Key, -Drug_function) %>%
  dplyr::mutate(
    X = factor(
      X %>% str_replace("_"," ") %>% str_to_title(),
      levels = c("Drug Class", "Drug Name", "Drug Target", "Drug Function")
    ),
    Key = factor(
      Key, levels = Key %>% unique() %>% rev()
    )
  ) %>%
  ggplot(aes(x = X, y = Key, label = Label)) +
  geom_text(size = 3, hjust = 0.5) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black"),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        panel.grid = element_blank()
  )


pdf2(file = str_c(out_tables, "table_drug_class_summary.pdf"), height = 5, width = 3.75)
print(p)
dev.off()

#
# ==============================================================================



# 8. Summarize non relevant patient characterisitics
# ==============================================================================

clin <- geo_clin_list$`Non-corrected`
dim(clin) # 3736 92

clin %>%
  group_by(Regimen_updated) %>%
  summarise(N_dataset = Series_matrix_accession %>% unique() %>% length(),
            N_sample = n())
#   Regimen_updated N_dataset N_sample
# 1 adj                     5      754
# 2 neoadj                 24     2982

clin %>%
  filter(Regimen_updated == "neoadj") %>%
  select(Series_matrix_accession) %>%
  deframe() %>%
  unique()
# [1] "GSE25066"         "GSE41998"         "GSE20194"
# [4] "GSE20271"         "GSE6861"          "GSE22358"
# [7] "GSE50948"         "GSE22226_GPL4133" "GSE16446"
# [10] "GSE32646"         "GSE130786"        "GSE22093"
# [13] "GSE4779"          "GSE76360"         "GSE21997_GPL5325"
# [16] "GSE21997_GPL7504" "GSE42822"         "GSE66305"
# [19] "GSE66999"         "GSE28844"         "GSE23988"
# [22] "GSE18728"         "GSE21974"         "GSE143846"
# 24 series matrices from 23 datasets

clin %>%
  filter(Regimen_updated == "adj") %>%
  select(Series_matrix_accession) %>%
  deframe() %>%
  unique()
# [1] "GSE20685" "GSE45255" "GSE69031" "GSE19615" "GSE16391"
# 5 series matrices from 5 datasets


table(clin$Ethnicity) %>% sum()
# [1] 716
table(clin$Histology) %>% sum()
# [1] 913
table(clin$Menopause) %>% sum()
# [1] 348
table(clin$Til) %>% sum()
# [1] 154
table(clin$Til_stromal) %>% sum()
# [1] 0
table(clin$Pdl1_stromal) %>% sum()
# [1] 0
table(clin$Pdl1_tumor) %>% sum()
# [1] 0
table(clin$Timepoint) %>% sum()
# [1] 279

table(clin$Timepoint)
# Mid_treatment Post_treatment  Pre_treatment
# 18            126            135
table(clin$Timepoint %>% is.na())
# FALSE  TRUE
# 279  3457



clin %>%
  filter(!is.na(Timepoint)) %>%
  select(Series_matrix_accession) %>% deframe() %>% unique()
# [1] "GSE76360" "GSE28844" "GSE18728" "GSE21974"

clin %>%
  filter(Timepoint=="Mid_treatment") %>%
  select(Series_matrix_accession) %>% deframe() %>% unique()
# [1] "GSE18728"
clin %>%
  filter(Timepoint=="Post_treatment") %>%
  select(Series_matrix_accession) %>% deframe() %>% unique()
# [1] "GSE76360" "GSE28844" "GSE18728" "GSE21974"
clin %>%
  filter(Timepoint=="Pre_treatment") %>%
  select(Series_matrix_accession) %>% deframe() %>% unique()
# [1] "GSE76360" "GSE28844" "GSE18728" "GSE21974"

clin %>%
  group_by(Timepoint, Regimen) %>%
  summarise(N = n())
#   Timepoint      Regimen        N
# 1 Mid_treatment  neoadj        18
# 2 Post_treatment neoadj       126
# 3 Pre_treatment  neoadj       135
# 4 NA             adj          754
# 5 NA             neoadj      2175
# 6 NA             neoadj_adj   528

# Timepoint present only in neoadjuvant datasets
# NOTE: Remove mid/post treatment biopsies before any analysis !!!!!!!!!


clin %>%
  filter(Regimen == "neoadj") %>%
  group_by(Response, Response_clinical) %>%
  summarise(N = n())
#   Response Response_clinical     N
# 1 npCR     MPG1                  4
# 2 npCR     MPG2                 24
# 3 npCR     MPG3                 20
# 4 npCR     No_response          38
# 5 npCR     non_downstage        35
# 6 npCR     PD                    1
# 7 npCR     PR                   57
# 8 npCR     SD                    5
# 9 npCR     NA                 1555
# 10 pCR      CR                   13
# 11 pCR      downstage             9
# 12 pCR      MPG4                  7
# 13 pCR      MPG5                  6
# 14 pCR      Response             23
# 15 pCR      NA                  583
# 16 NA       NA                   74

# Clinical responses were mapped to pathological response


# Identify other relevant characteristics to summarise
names(clin)
# [1] "Series_matrix_accession"   "Sample_title"
# [3] "Sample_geo_accession"      "Sample_channel_count"
# [5] "Sample_platform_id"        "Sample_organism_ch1"
# [7] "Sample_organism_ch2"       "Sample_source_name_ch1"
# [9] "Sample_source_name_ch2"    "Sample_type"
# [11] "Sample_type_ch2"           "Age"
# [13] "Age_bin"                   "Age_detailed"
# [15] "Grade"                     "Node"
# [17] "Node_bin"                  "Node_cat"
# [19] "Size"                      "Size_bin"
# [21] "Size_cat"                  "Hr"
# [23] "Er"                        "Pr"
# [25] "Her2"                      "Regimen"
# [27] "Arm_neoadj"                "Arm_adj"
# [29] "Arm_detailed"              "Response"
# [31] "Response_clinical"         "Response_pathological"
# [33] "Event_dfs"                 "Time_dfs"
# [35] "Timepoint"                 "Histology"
# [37] "Gender"                    "Ethnicity"
# [39] "Menopause"                 "Til"
# [41] "Til_stromal"               "Pdl1_stromal"
# [43] "Pdl1_tumor"                "Regimen_updated"
# [45] "Arm"                       "Arm_updated"
# [47] "Arm_expanded"              "Has_anthracycline"
# [49] "Has_antimetabolite"        "Has_alkylating_agent"
# [51] "Has_taxane"                "Has_chemo"
# [53] "Class_chemo"               "Has_aaa"
# [55] "Has_fac"                   "Has_fec"
# [57] "Has_cmf"                   "Has_hormone"
# [59] "Has_her2_agent"            "Has_other"
# [61] "Has_no_treatment"          "Name_anthracycline"
# [63] "Name_antimetabolite"       "Name_alkylating_agent"
# [65] "Name_taxane"               "Name_hormone"
# [67] "Name_her2_agent"           "Name_other"
# [69] "Arm_anthracycline"         "Arm_aaa"
# [71] "Arm_fac"                   "Arm_fec"
# [73] "Arm_cmf"                   "Arm_chemo"
# [75] "Arm_her2"                  "Arm_hormone"
# [77] "Arm_other"                 "Series_sample_procurement"
# [79] "Series_archive_method"     "Has_expr_pooled"
# [81] "Subtype_ihc"               "Rle_min"
# [83] "Rle_first_qu"              "Rle_median"
# [85] "Rle_mean"                  "Rle_third_qu"
# [87] "Rle_max"                   "Er_gene"
# [89] "Her2_gene"                 "Subtype_pam50"
# [91] "Pam50_centroid_corr"       "Subtype_pam50_2"

#
# ==============================================================================



# 9. Genarate a table containing GEO search with series details,
# selected/rejected flag, reason for rejection.
# ==============================================================================

dim(geo_series_summary) # 172 19
dim(geo_series_matrix_summary) # 175 30

glimpse(geo_series_summary)
glimpse(geo_series_matrix_summary)

geo_series_matrix_summary$Series_matrix_selected %>% table()
# no yes
# 137  38
geo_series_matrix_summary$Series_selected %>% table()
# no yes
# 132  43

geo_series_summary$Series_selected %>% table()
# no  yes
# 132  40


# Genrate series summary from series_matrix_summary and then collapse selected
# multi-series-matrix per geo series if present.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

dataset_summary <- geo_series_matrix_summary %>%
  left_join(geo_clin_list$`Non-corrected` %>%
              dplyr::select(Series_matrix_accession) %>%
              dplyr::mutate(Final_filtering = "yes") %>%
              dplyr::distinct(Series_matrix_accession, .keep_all =TRUE),
            by = "Series_matrix_accession") %>%
  dplyr::select(
    "Series_id",
    "Series_accession",
    "Series_matrix_accession",
    "Series_selected",
    "Series_matrix_selected",
    "Final_filtering",
    "Series_comment",
    "Series_matrix_comment",
    "Series_title",
    "Series_description",
    "Series_organism",
    "Series_type",
    "Series_platform",
    "Series_sample_size",
    "Series_ftp",
    "Series_arm",
    "Series_arm_description",
    "Series_survival",
    "Series_regimen",
    "Series_sample_procurement",
    "Series_sample_procurement_details",
    "Series_archive_method",
    "Series_archive_details"
  ) %>%
  dplyr::mutate(
    Series_id = as.character(Series_id),
    Series_status = paste(
      #!! str_c() wiil propogate NAs
      Series_selected, Series_matrix_selected, Final_filtering)
  ) %>%
  dplyr::select(
    -c("Series_selected",
       "Series_matrix_selected",
       "Final_filtering"
    )) %>%
  dplyr::mutate(

    Series_status = purrr::map_chr(
      Series_status,
      ~switch (
        .x,
        "no no NA" = "non_compactable",
        "yes no NA" = "non_compactable",
        "yes yes NA" = "poor_genome_representation",
        "yes yes yes" = "selected",
        "non_compactable"
      )
      ),

    Series_status = ifelse(
      Series_accession == "GSE4056",
      "missingness",
      Series_status),

    Series_selected = ifelse(
      Series_status == "selected",
      "yes", "no")
  )



dataset_summary$Series_status %>% table()
# missingness            non_compactable poor_genome_representation
# 1                        136                          9
# selected
# 29
dataset_summary$Series_selected %>% table()
# no yes
# 146  29




dim(dataset_summary) # 175 22 !!!
dataset_summary$Series_accession %>% unique() %>% length() # 172 !!!

# The difference is due to multiple series matrices per geo series.
# Collpase the series matrices manually, to make the summary sereis centric !!!!!
# Only collapse the selected series matrices, to avoid ambigiuty.

idx <- purrr::map(dataset_summary$Series_accession %>% unique(),
                  function(x, id){which(x==id)},
                  dataset_summary$Series_accession)
names(idx) = dataset_summary$Series_accession %>% unique()
idx <- idx[purrr::map_int(idx, length) >1]
# $GSE22226
# [1] 59 60
# $GSE21997
# [1] 94 95 96

dataset_summary[unlist(idx), ] %>%
  dplyr::select(Series_accession,Series_matrix_accession, Series_status)
#   Series_accession Series_matrix_accession Series_status
# 1 GSE22226         GSE22226_GPL1708        poor_genome_representation !!! discard
# 2 GSE22226         GSE22226_GPL4133        selected

# 3 GSE21997         GSE21997_GPL1390        poor_genome_representation !!! discard
# 4 GSE21997         GSE21997_GPL5325        selected
# 5 GSE21997         GSE21997_GPL7504        selected

# Removing non-selected seires matrices from serieses containing multiple series matrices
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
purrr::map_int(idx, ~(.x[1]))
# GSE22226 GSE21997
# 59       94
dataset_summary <- dataset_summary[-purrr::map_int(idx, ~(.x[1])),]


# Repeat the collapsing process (as indices changed due to the above filtering)
# !!!!!!!!!!!!!!!!!!!!!!!1

idx <- purrr::map(dataset_summary$Series_accession %>% unique(),
                  function(x, id){which(x==id)},
                  dataset_summary$Series_accession)
names(idx) = dataset_summary$Series_accession %>% unique()
idx <- idx[purrr::map_int(idx, length) >1]
# $GSE21997
# [1] 93 94

dataset_summary[unlist(idx), ] %>%
  dplyr::select(Series_accession,Series_matrix_accession, Series_status)
#   Series_accession Series_matrix_accession Series_status
# 1 GSE21997         GSE21997_GPL5325        selected
# 2 GSE21997         GSE21997_GPL7504        selected


dataset_summary$Series_matrix_accession[idx$GSE21997[1]] =
  dataset_summary$Series_matrix_accession[idx$GSE21997] %>% str_c(collapse = "///")

idx <- purrr::map(idx,~(.x <- .x[-1])) %>% unlist() # keeping only the updated row
dataset_summary <- dataset_summary[-idx,]

dim(dataset_summary) # 172 22


dataset_summary$Series_status %>% table()
# missingness            non_compactable poor_genome_representation
# 1                        136                          7
# selected
# 28
dataset_summary$Series_selected %>% table()
# no yes
# 144  28
# In agreement with supplementary figure 1: Integrated script structure and
# study flow chart.


write_tsv(x = dataset_summary,
          path = str_c(out_tables, "table_geo_series_summary.tsv"))

#
# ==============================================================================



# 10. Comapre curatedBreastData package with present study
# ==============================================================================


# a. Dataset's presence or absence in planey and present dataset
# b. Detected_and_integrated_in_both
# c. Detected_and_integrated_only_in_planey
# d. Detected_and_integrated_only_in_present
# e. Detected_only_in_present
# f. Detected_in_both_integrated_only_in_planey



# Prepare dataset summaries
# >>>>>>>>>>>>>>>>>>>>>>>>>


# Planey et.al R package "curatedBreastData"

glimpse(planey)
# relevant variables
# $ study_ID                                       <int> 32646, 32646, 32646, 32646, 32…
# $ patient_ID                                     <int> 809184, 809185, 809186, 809187…
# $ GEO_GSMID                                      <int> 809184, 809185, 809186, 809187…

planey_annot[1:5,]
# 1 dbUniquePatient… Unique patient id created for this database.
# 2 study_ID         GEO GSE study ID.
# 3 patient_ID       patient ID used in phenoData, expression matrix.
# 4 GEO_GSMID        GEO GSM patient (sample) ID. Often the same as patient ID, most samp…
# 5 platform_ID      general platform ID.

identical(planey$patient_ID, planey$GEO_GSMID) # TRUE

dim(planey) # 2719 139

planey <- planey %>%
  dplyr::mutate(Series_accession = str_c("GSE", study_ID), # used for merging
                Series_accession_curatedBreastData = Series_accession,
                Series_selected_curatedBreastData = "yes") %>%
  dplyr::select(Series_accession, Series_accession_curatedBreastData,
                Series_selected_curatedBreastData) %>%
  dplyr::distinct(Series_accession, .keep_all = TRUE)

dim(planey) # 24 3



# Present study

dim(dataset_summary) # 172 22


# Merged Planey and present study

merged_dataset_summary <- planey %>%
  dplyr::full_join(dataset_summary, by = "Series_accession") %>%
  dplyr::rename_all(~str_replace_all(.x,"curatedBreastData","planey_dataset")) %>%
  dplyr::rename("Series_selected_present_dataset" = "Series_selected")



dataset_comparison <- list()

# a. Datasets in planey and present dataset
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dataset_comparison[["Planey_and_present_study"]] <- merged_dataset_summary %>%
  dplyr::group_by(Series_selected_present_dataset,
                  Series_selected_planey_dataset) %>%
  dplyr::summarise(Dataset_count = n()) %>%
  dplyr::rename_all(~str_replace_all(.x,"Series_selected_","")) %>%
  dplyr::rename_all(~str_to_sentence(.x))
# NA: means not detected according to dataset search criteria.
# Note that for Planey_dataset, the list of dataset detected according to
# search criteria was not available.
# Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3814460/pdf/amia_tbi_2013_138.pdf
# no: means the dataset detected by search criteria was not selected for
# further processing according to dataset filtering criteria.
# yes: means the dataset detected by search criteria was selected for
# further processing.

# "Detect_integrate_both"
# "Detect_integrate_planey"
# "Detect_integrate_present"
# "Detect_only_in_present"
# "Detect_both_integrae_planey"


# b. Detected_and_integrated_in_both
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dataset_comparison[["Both_study_detect_integrate"]] <- merged_dataset_summary %>%
  dplyr::filter(Series_selected_present_dataset == "yes" &
                  Series_selected_planey_dataset == "yes")



# c. Detected_and_integrated_only_in_planey
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dataset_comparison[["Planey_only_detect_integrate"]] <- merged_dataset_summary %>%
  dplyr::filter(Series_selected_present_dataset %>% is.na()  &
                  Series_selected_planey_dataset == "yes")



# d. Detected_and_integrated_only_in_present
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dataset_comparison[["Present_only_detect_integrate"]] <- merged_dataset_summary %>%
  dplyr::filter(Series_selected_present_dataset == "yes" &
                  Series_selected_planey_dataset %>% is.na())




# e. Detected_only_in_present
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dataset_comparison[["Present_only_detect"]] <- merged_dataset_summary %>%
  dplyr::filter(Series_selected_present_dataset == "no" &
                  Series_selected_planey_dataset %>% is.na())



# f. Detected_in_both_integrated_only_in_planey
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dataset_comparison[["Both_detect_planey_only_integrate"]] <- merged_dataset_summary %>%
  dplyr::filter(Series_selected_present_dataset == "no" &
                  Series_selected_planey_dataset == "yes")



write.xlsx(x = dataset_comparison$Planey_and_present_study %>% as.data.frame(),
           file = str_c(out_tables, "table_dataset_comparison_Planey_vs_Prenset.xlsx"),
           sheetName="Planey_and_present_study",
           append = FALSE)
write.xlsx(x = dataset_comparison$Both_study_detect_integrate %>% as.data.frame(),
           file = str_c(out_tables, "table_dataset_comparison_Planey_vs_Prenset.xlsx"),
           sheetName="Both_study_detect_integrate",
           append = TRUE)
write.xlsx(x = dataset_comparison$Planey_only_detect_integrate %>% as.data.frame(),
           file = str_c(out_tables, "table_dataset_comparison_Planey_vs_Prenset.xlsx"),
           sheetName="Planey_only_detect_integrate",
           append = TRUE)
write.xlsx(x = dataset_comparison$Present_only_detect_integrate %>% as.data.frame(),
           file = str_c(out_tables, "table_dataset_comparison_Planey_vs_Prenset.xlsx"),
           sheetName="Present_only_detect_integrate",
           append = TRUE)
write.xlsx(x = dataset_comparison$Present_only_detect %>% as.data.frame(),
           file = str_c(out_tables, "table_dataset_comparison_Planey_vs_Prenset.xlsx"),
           sheetName="Present_only_detect",
           append = TRUE)
write.xlsx(x = dataset_comparison$Both_detect_planey_only_integrate %>% as.data.frame(),
           file = str_c(out_tables, "table_dataset_comparison_Planey_vs_Prenset.xlsx"),
           sheetName="Both_detect_planey_only_integrate",
           append = TRUE)


#
# ==============================================================================

