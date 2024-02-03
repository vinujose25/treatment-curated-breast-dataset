# s6_process_and_validate_integrated_expression_matrix


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# 1. Dataset effect correction method comparison.
#
#    Compare the following dataset effect correction methods
#    a. original (log2) vs
#    b. combat(log2 + combat) vs
#    c. quantile (log2 + per dataset quantile gene scaling) vs
#    d. mad_quantile (log2 + mad sample scaling + per dataset quantile gene scaling)
#
#    with respect to the  following metrics from each version of corrected dataset
#
#    a. RLE plots (save plots)
#    b. PAM50 vs Subtype-IHC Kappa agreement
#    c. ER/HER2 IHC status prediction ability (AUC) (save plots)


# Script structure
# >>>>>>>>>>>>>>>>
#
# 1. Prepare data
# 2. Compute list of dataset effect corrected data
# 3. Compute list of clinical data with respective qc metrics from corrected data
# 4. generate and display rle plots a  nd other tables



# 1. Prepare data
# ==============================================================================

load("results/data/geo_expr.RData")
load("results/data/geo_clin.RData")


# Data-metadata congruence
# >>>>>>>>>>>>>>>>>>>>>>>

clin <- geo_clin %>% dplyr::filter(Has_expr_pooled == "yes")
expr <- geo_expr %>% dplyr::select(1, all_of(clin$Sample_geo_accession))

identical(names(expr)[-1], clin$Sample_geo_accession) #TRUE

dim(clin) # 3736 80
dim(expr) # 9184 3737


clin %>%
  group_by(Regimen_updated) %>%
  summarise(N_dataset = Series_matrix_accession %>% unique() %>% length(),
            N_sample = n())
#   Regimen_updated N_dataset N_sample
# 1 adj                     5      754
# 2 neoadj                 24     2982

clin %>%
  dplyr::filter(Regimen_updated == "neoadj") %>%
  dplyr::select(Series_matrix_accession) %>%
  deframe() %>% unique()
# [1] "GSE25066"         "GSE41998"         "GSE20194"
# [4] "GSE20271"         "GSE6861"          "GSE22358"
# [7] "GSE50948"         "GSE22226_GPL4133" "GSE16446"
# [10] "GSE32646"         "GSE130786"        "GSE22093"
# [13] "GSE4779"          "GSE76360"         "GSE21997_GPL5325"
# [16] "GSE21997_GPL7504" "GSE42822"         "GSE66305"
# [19] "GSE66999"         "GSE28844"         "GSE23988"
# [22] "GSE18728"         "GSE21974"         "GSE143846"
# 24 series matrices representing 23 neoadj datasets !!!!!!!!!!

clin %>%
  dplyr::filter(Regimen_updated == "adj") %>%
  dplyr::select(Series_matrix_accession) %>%
  deframe() %>% unique()
# [1] "GSE20685" "GSE45255" "GSE69031" "GSE19615" "GSE16391"
# 5 series matrices represnting 5 adj datasets !!!!!!!!!!


# pryr::object_size(expr[,-1]) # 275 MB


# Create a list of expr and clin data to hold processed expr and expr qc metrics

geo_expr_list <- list(
  # original (log2)
  Original = NA,
  # combat (log2 + combat)
  Combat = NA,
  # per dataset quantile gene scaling (log2 + quantile)
  Quantile = NA,
  # mad sample scaling + per dataset quantile gene scaling (log2 + mad_quantile)
  # filtered 2 samples due to zero MAD value !!!!!!!!!!
  Madquantile = NA
)


geo_clin_list <- list(
  # original (log2)
  Original = NA,
  # combat (log2 + combat)
  Combat = NA,
  # per dataset quantile gene scaling (log2 + quantile)
  Quantile = NA,
  # mad sample scaling + per dataset quantile gene scaling (log2 + mad_quantile)
  # filtered 2 samples due to zero MAD value !!!!!!!!!!
  Madquantile = NA
)

#
# ==============================================================================




# 2. Compute list of dataset effect corrected data
# ==============================================================================

# Originl
geo_expr_list$Original <- expr

# Combat
geo_expr_list$Combat <- bind_cols(
  expr %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = FALSE # default is FALSE
  ) %>%
    as_tibble()
)
# Found29batches
# Adjusting for0covariate(s) or covariate level(s)
# Standardizing Data across genes
# Fitting L/S model and finding priors
# Finding parametric adjustments
# Adjusting the Data


# Quantile: Per dataset gene scaling
geo_expr_list$Quantile <- expr %>%
  dataset_gene_scaling(dataset_map = clin)


# Mad_quantile: MAD sample scaling + per dataset gene scaling
geo_expr_list$Madquantile <- bind_cols(
  expr %>% dplyr::select(1),
  expr[, -1] %>%
    purrr::map(~(DescTools::RobScale(x = .x))) %>%
    bind_cols()
) %>%
  dataset_gene_scaling(dataset_map = clin)
# Error message !!!
# Error in quantile.default(x, probs = 1 - (q/2), na.rm = na.rm) :
#  missing values and NaN's not allowed if 'na.rm' is FALSE

# Reason for Mad_quantile error !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# MAD sample scaling made two expression profiles invalid due to unknown systematic issue,
# probabily a sample processing issue in the original study.
xx <- expr[, -1] %>%
  purrr::map(~(DescTools::RobScale(x = .x))) %>%
  bind_cols()
purrr::map(xx, median) %>% unlist() %>% is.na() %>% table()
# FALSE  TRUE
# 3734     2
purrr::map(xx, median) %>% unlist() %>% is.na() %>% which()
# GSM1635819 GSM1635820
# 3377       3378

expr$GSM1635819 %>% range() # 3.818736 15.139882
expr$GSM1635819 %>% median() # 3.818736 !!! Median is the minimum value
expr$GSM1635819 %>% mad() # 0
(expr$GSM1635819 == expr$GSM1635819 %>% median()) %>% table()
(expr$GSM1635819 == expr$GSM1635819 %>% min()) %>% table()
# FALSE  TRUE
# 4508  4676
# More than half of the genes have minimum expression value and are identical.
# This implies that MAD(Median absolute deviation) will be zero for this sample.
# Hence MAD rescaled values will be NAs; (x-median)/MAD = NA/Inf
# Roughly half the expression got lowest expression value
# Seems like a systaematic issue, most probably an array processing issue.
# Discard !!!!!!!!!!!!!


expr$GSM1635820 %>% range() # 4.101677 15.445718
expr$GSM1635820 %>% median() # 4.101677 !!! Median is the minimum value
expr$GSM1635820 %>% mad() # 0
(expr$GSM1635820 == expr$GSM1635820 %>% median()) %>% table()
(expr$GSM1635820 == expr$GSM1635820 %>% min()) %>% table()
# FALSE  TRUE
# 3526  5658
# More than half of the genes have minimum expression value and are identical.
# This implies that MAD(Median absolute deviation) will be zero for this sample.
# Hence MAD rescaled values will be NAs; (x-median)/MAD = NA/Inf
# Roughly half the expression got lowest expression value
# Seems like a systaematic issue, most probably an array processing issue.
# Discard !!!!!!!!!!!!!


# MAD for all samples
xx <- expr[, -1] %>%
purrr::map(~(mad(x = .x))) %>%
  unlist()
table(xx == 0) # no.of samples with MAD == 0
# FALSE  TRUE
# 3734     2
summary(xx)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.000   1.867   2.117   2.033   2.440   4.383


clin %>%
  dplyr::filter(Sample_geo_accession %in% all_of(c("GSM1635819", "GSM1635820"))) %>%
  glimpse()
# These two samples are from a single study(Series), GSE66999.
# Both are ER and PR positive and HER2 negative.
# Neoadjuvant regimen.
# Treatment: Docetaxel + Epirubicine
# Core needle biopsy, preserved by snap freezing
# Profiled by Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (year 2008)
#   having 41108 probes


# Remove "GSM1635819" and  "GSM1635820", and redo Mad_quantile correction !!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
head2(clin) # 3736 80
head2(expr) # 9184 3737
clin2 <- clin %>%
  dplyr::filter(!(Sample_geo_accession %in% all_of(c("GSM1635819", "GSM1635820"))))
expr2 <- expr %>% dplyr::select(1, all_of(clin2$Sample_geo_accession))
head2(clin2) # 3734 81
head2(expr2) # 9184 3735

geo_expr_list$Madquantile <- bind_cols(
  expr2 %>% dplyr::select(1),
  expr2[, -1] %>%
    purrr::map(~(DescTools::RobScale(x = .x))) %>%
    bind_cols()
) %>%
  dataset_gene_scaling(dataset_map = clin2)


geo_expr_list %>% purrr::map(dim)
# $Original
# [1] 9184 3737
# $Combat
# [1] 9184 3737
# $Quantile
# [1] 9184 3737
# $Madquantile
# [1] 9184 3735 # 2 samples were filtered out due to zero MAD


# Save dataset effect corrected data list !!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# save(geo_expr_list, file = str_c(out_data, "geo_expr_list.RData"))

#
# ==============================================================================




# 3. Compute list of clinical data with respective qc metrics from corrected data
# ==============================================================================


# IHC subtyping, not depend on molecular data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

table(is.na(clin$Hr)) # all NAs
# TRUE
# 3736

clin <- clin %>%
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

table(clin$Subtype_ihc)
# HER2   HR   TN
# 756 1232  850
table(is.na(clin$Subtype_ihc))
# FALSE  TRUE
# 2838   898




# Module list
# >>>>>>>>>>>>

# ER/HER2/Immune gene
module_list <- list(
                   Er_gene = tibble(
                     Ncbi_gene_id = "ncbi_2099",
                     Hugo_gene_symbol = "ESR1",
                     Direction = 1
                   ),
                   Her2_gene = tibble(
                     Ncbi_gene_id = "ncbi_2064",
                     Hugo_gene_symbol = "ERBB2",
                     Direction = 1
                   )
                 )



# PAM50 centroids
pam50 <- list()
pam50[["centroid"]] <- read_tsv(file = "data/pam50/PAM50/bioclassifier_R/pam50_centroids.txt")
pam50[["annot"]] <- read_tsv(file = "data/pam50/PAM50/bioclassifier_R/pam50_annotation.txt")

pam50 <- pam50$centroid %>%
  dplyr::rename(pcrID = "...1") %>%
  dplyr::left_join(pam50$annot, by = "pcrID") %>%
  dplyr::mutate(EntrezGene = str_c("ncbi_", EntrezGene)) %>%
  dplyr::filter(EntrezGene %in% expr$Ncbi_gene_id) %>%
  dplyr::rename(Ncbi_gene_id = "EntrezGene")
nrow(pam50) # 36




# Function to compute qc metrics, RLE and others
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

update_clin <- function(clin, expr, module_list, pam50){


  # RLE
  # >>>

  cat("RLE...")

  x <- expr

  x[, -1] <- x[, -1] %>%
    apply(1,function(xx){xx-median(xx)}) %>% # gene centering, the resut is transposed
    t() %>% as_tibble()

  rle <- bind_cols(
    tibble(Sample_geo_accession = names(x)[-1]),
    x[, -1] %>% purrr::map_dfr(summary) %>% purrr::map_dfr(~(as.numeric(.x)))
  ) %>%
    dplyr::rename(First_Qu = "1st Qu.",
                  Third_Qu = "3rd Qu.") %>%
    dplyr::rename_all(~(str_c("Rle_", .x) %>% str_remove_all("\\."))) %>%
    dplyr::rename(Sample_geo_accession = "Rle_Sample_geo_accession") %>%
    dplyr::rename_all(~(str_to_title(.x)))

  cat(identical(rle$Sample_geo_accession, clin$Sample_geo_accession))
  cat("\n")

  clin <- clin %>% left_join(rle, by = "Sample_geo_accession")



  # Relevant genes subsetting
  # >>>>>>>>>>>>>>>>>>>>>>>>>

  g <- purrr::map(module_list,~(.x$Ncbi_gene_id)) %>%
    unlist() %>%
    c(pam50$Ncbi_gene_id) %>%
    unique()

  expr <- expr %>% dplyr::filter(Ncbi_gene_id %in% g)



  # Module score
  # >>>>>>>>>>>>

  cat("Score...")

  # converting x to genes on column and samples on rows
  x = t_tibble(expr , names_x_desc = "Sample_geo_accession")


  # Module score as weighted average
  score <- get_module_score(x = x, module_list = module_list, by = "Ncbi_gene_id")

  cat(identical(score$Sample_geo_accession, clin$Sample_geo_accession))
  cat("\n")

  clin <- clin %>% left_join(score, by = "Sample_geo_accession")



  # PAM50 subtyping
  # >>>>>>>>>>>>>>>

  cat("PAM50...")

  x <- pam50 %>%
    dplyr::select(Ncbi_gene_id) %>%
    dplyr::left_join(expr, by = "Ncbi_gene_id")


  # Median gene centering recommanded
  x[, -1] <- x[, -1] %>%
    apply(1,function(xx){xx-median(xx)}) %>% # gene centering, the resut is transposed
    t() %>% as.data.frame()


  # Pearson correlation
  # (Median centering is valid only with pearson corelation.
  # Foe spearman correlation, median centering has no effect).
  x_corr <- cor(pam50[, c("Basal", "Her2", "LumA", "LumB", "Normal")],
                x[,-1],
                method = "pearson")


  # Max correlated centroid as subtype, but if max corr <.1 the tumor is unclassified
  subtype <- purrr::map_dfr(
    colnames(x_corr),
    function(sample, x_corr, subtype){

      corr <- x_corr[, sample]

      list(
        Sample_geo_accession = sample,
        Subtype_pam50 = if_else(max(corr) < .1,
                                "Unclassified",
                                subtype[which(corr == max(corr))]),
        Pam50_centroid_corr = corr[which(corr == max(corr))]
      )

    },
    x_corr,
    subtype = rownames(x_corr)
  ) %>%
    dplyr::mutate(
      Subtype_pam50_2 = purrr::map_chr(
        Subtype_pam50,
        function(x){
          case_when(
            x == "Basal" ~ "TN",
            x == "Her2" ~ "HER2",
            x == "LumA" ~ "HR",
            x == "LumB" ~ "HR",
            TRUE ~ NA_character_ # Normals will be set as NA
          )
        }
      )

    )

  cat(identical(subtype$Sample_geo_accession, clin$Sample_geo_accession))
  cat("\n")

  clin <- clin %>% left_join(subtype, by = "Sample_geo_accession")


  clin
}





# Populate geo_clin_list with respective qc metrics
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# geo_clin_list <- list(
#   # original (log2)
#   Original = NA,
#   # combat (log2 + combat)
#   Combat = NA,
#   # per dataset quantile gene scaling (log2 + quantile)
#   Quantile = NA,
#   # mad sample scaling + per dataset quantile gene scaling (log2 + mad_quantile)
#   # filtered 2 samples due to zero MAD value !!!!!!!!!!
#   Madquantile = NA
# )



# original (log2)
# x <- expr
geo_clin_list$Original <- update_clin(
  clin, geo_expr_list$Original, module_list, pam50
)

# combat (log2 + combat)
# x <- expr_combat
# geo_clin_list$Combat <- update_clin(clin, x, module_list, pam50)
geo_clin_list$Combat <- update_clin(
  clin, geo_expr_list$Combat, module_list, pam50
)


# per dataset quantile gene scaling (log2 + quantile)
# x <- expr_quantile
# geo_clin_list$Quantile <- update_clin(clin, x, module_list, pam50)
geo_clin_list$Quantile <- update_clin(
  clin, geo_expr_list$Quantile, module_list, pam50
)


# mad sample scaling + per dataset quantile gene scaling (log2 + mad_quantile)
identical(
  clin %>%
    dplyr::filter(Sample_geo_accession %in% names(geo_expr_list$Madquantile)) %>%
    dplyr::select(Sample_geo_accession) %>%
    deframe(),
  names(geo_expr_list$Madquantile)[-1]
)
identical(
  clin2$Sample_geo_accession,
  names(geo_expr_list$Madquantile)[-1]
)
# TRUE
# x <- expr_madquantile
geo_clin_list$Madquantile <- update_clin(
  clin = clin %>%
    dplyr::filter(Sample_geo_accession %in% names(geo_expr_list$Madquantile)),
  geo_expr_list$Madquantile, module_list, pam50
)

geo_clin_list %>% purrr::map(dim)
# $Original
# [1] 3736   92
# $Combat
# [1] 3736   92
# $Quantile
# [1] 3736   93 # additional column for NAs in RLE ("Rle_na's")
# $Madquantile
# [1] 3734   92


# Save qc metrics associated with dataset effect correction methods !!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# save(geo_clin_list, file = str_c(out_data, "geo_clin_list.RData"))

#
# ==============================================================================




# 4. generate and disply rle plots and other tables
# ==============================================================================


# RLE plots
# >>>>>>>>>

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
                  str_replace("Mad_quantile", "Mad_quantile*"), # to highlight filtering of 2 samples
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

print(p)




# PAM50 vs Subtype-IHC
# >>>>>>>>>>>>>>>>>>>>

# Kappa,κ ranges from 0 to 1, with 0 indicating no relation and 1 indicating
# a perfect concordance. Typically qualitative descriptions are associated
# with intervals [κ ≤ 0.20, slight agreement; 0.20 < κ ≤ 0.40, fair agreement;
# 0.40 < κ ≤ 0.60, moderate agreement, 0.60 < κ ≤ 0.80, substantial agreement;
# and 0.80 < κ ≤ 0.99, almost perfect agreement
# See Benjamins Three gene model paper

table(geo_clin_list$Original$Subtype_ihc) # Reference
# HER2   HR   TN
# 756  1232  850
table(geo_clin_list$Original$Subtype_pam50)
# Basal         Her2         LumA         LumB       Normal Unclassified
# 1089          639         1140          492          341           35
table(geo_clin_list$Original$Subtype_pam50_2) # Test subtype
# HER2   HR   TN
# 639 1632 1089


xx <- geo_clin_list
purrr::map_dfr(
  names(geo_clin_list),
  function(nme, xx){
    x = xx[[nme]]
    list(type = nme,
         kappa = irr::kappa2(ratings = x %>%
                               dplyr::select(Subtype_ihc, Subtype_pam50_2))$value
    )
  }, xx = geo_clin_list) %>%
  as.data.frame()
#           type     kappa
# 1     Original 0.6210252
# 2       Combat 0.6180238
# 3     Quantile 0.6342430
# 4 Mad_quantile 0.6276156 !!! filtered 2 samples with MAD = 0




# ER/HER2 IHC status prdiction
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>

xx <- purrr::map(
  geo_clin_list,
  function(x){
    list(er_gene = pROC::roc(formula = Er  ~ Er_gene, data = x),
         her2_gene = pROC::roc(formula = Her2 ~ Her2_gene, data = x))
  }
)



purrr::map_dfr(
  names(xx),
  function(nme, xx){
    x = xx[[nme]]
    list(type = nme,
         # er_sig = auc(x$er_sig) %>% as.numeric(),
         # her2_sig = auc(x$her2_sig) %>% as.numeric(),
         er_gene = auc(x$er_gene) %>% as.numeric(),
         her2_gene = auc(x$her2_gene) %>% as.numeric())
  }, xx) %>%
  as.data.frame()
#           type   er_gene her2_gene
# 1     Original 0.7912363 0.6055798
# 2       Combat 0.8803596 0.7290365
# 3     Quantile 0.8947382 0.7762824
# 4 Mad_quantile 0.8921448 0.7704629 !!! filtered 2 samples with MAD = 0


#
# ==============================================================================








