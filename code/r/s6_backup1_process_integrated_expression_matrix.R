# s6_backup1_process_integrated_expression_matrix

# A) log2 vs log2+scaling vs log2+combat correction of expression data w.r.t
#   1) RLE plots (save plots)
#   2) ER/HER2 IHC status prediction ability (AUC) (save plots)
#   3) PAM50 vs Subtype-IHC Kappa agreement
#   4) ER - immune response association (save plots)
# B) Gruosso gene-module computation (Computed in A, while computing module scores)


load("results/main/geo_expr.RData")
load("results/main/geo_clin.RData")


# A) Log2 vs Log2+Scale vs Log2+ConBAT
# ==============================================================================


# data/metadata list
# >>>>>>>>>>>>>>>>>>

# 1) congruent data-metadata
# 2) IHC subtyping, not depend on molecular data
# 3) data list
# 4) metadata list
# 5) Module list
# 6) Module score
# 7) Pam50



# 1) congruent data-metadata

clin <- geo_clin %>% dplyr::filter(Has_expr_pooled == "yes")
expr <- geo_expr %>% dplyr::select(1, all_of(clin$Sample_geo_accession))

identical(names(expr)[-1], clin$Sample_geo_accession) #TRUE

pryr::object_size(expr[,-1]) # 275 MB




# Internal testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# RLE: log2 vs log2+combat_m vs log2 combat

# Rationale of this testing: Previous analysis showed log2+combat option is not
# suffice to correct batch effect, but in the below analysis combat successfully
# eliminated batch effect. Hence it is presumed that in the previous analysis,
# combat's mean_only option is set true while correcting batch effect.


dlist <- list( # template !!!!!!!!!!

  # original
  log2 = expr,

  # combat sample scaling and mean centering
  log2_combat =  bind_cols(
    expr %>% dplyr::select(1),
    sva::ComBat(
      dat = as.matrix(expr[, -1]),
      batch = clin$Series_matrix_accession
    ) %>%
      as_tibble()
  ),

  # combat mean centering
  log2_combat_m =  bind_cols(
    expr %>% dplyr::select(1),
    sva::ComBat(
      dat = as.matrix(expr[, -1]),
      batch = clin$Series_matrix_accession,
      mean.only = TRUE # default is FALSE
    ) %>%
      as_tibble()
  )
)


rle <- function(x){

x[, -1] <- x[, -1] %>%
  apply(1,function(xx){xx-median(xx)}) %>% # gene centering, the resut is transposed
  t() %>% as_tibble()

rle <- bind_cols(
  tibble(Sample_geo_accession = names(x)[-1]),
  x[, -1] %>% purrr::map_dfr(summary)
) %>%
  dplyr::rename(First_Qu = "1st Qu.",
                Third_Qu = "3rd Qu.") %>%
  dplyr::rename_all(~(str_c("Rle_", .x) %>% str_remove_all("\\."))) %>%
  dplyr::rename(Sample_geo_accession = "Rle_Sample_geo_accession") %>%
  dplyr::rename_all(~(str_to_title(.x)))

rle
}


mdlist <- purrr::map(
  dlist,
  function(x, clin){
    xrle <- rle(x)
    print(identical(xrle$Sample_geo_accession, clin$Sample_geo_accession))

    clin <- clin %>%
      left_join(xrle, by = "Sample_geo_accession")
  },
  clin
)

# RLE plots: ALL
xx <- purrr::map(
  names(mdlist),
  function(nme, mdlist){
    bind_cols(
      tibble(Batch_correction_methode = nme),
      mdlist[[nme]] %>%
        dplyr::select(Series_matrix_accession, Sample_geo_accession,
                      "Rle_min", "Rle_first_qu", "Rle_median", "Rle_mean",
                      "Rle_third_qu", "Rle_max")
    )
  },
  mdlist
) %>%
  bind_rows()

p <- xx %>%
  tidyr::gather("Rle_class", "Rle_value", "Rle_first_qu","Rle_median", "Rle_third_qu") %>%
  dplyr::mutate(Series_matrix_accession = factor(Series_matrix_accession,
                                                 levels = unique(Series_matrix_accession)),
                Sample_geo_accession = factor(Sample_geo_accession,
                                              levels = unique(Sample_geo_accession)),
                Batch_correction_methode = factor(Batch_correction_methode, levels = unique(Batch_correction_methode))) %>%
  ggplot(aes(x = Sample_geo_accession, y = Rle_value, group = Rle_class)) +
  geom_line(aes(color = Series_matrix_accession)) +
  guides(color = guide_legend(title = "Batch (Series_matrix_accession)",
                              title.position = "top")) +
  labs(x = "Samples", y = "RLE") +
  facet_wrap(facets = ~ Batch_correction_methode, scales = "free_y", ncol = 1) +
  theme(legend.position = "bottom",
        axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

  pdf2(file = str_c(outdir, "rle_internal_test.pdf"))
  print(p)
  dev.off()

rm(dlist, mdlist, p)
gc()

#
# End of Internal testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






# Sample filtering for reliable comaprison !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# MAD sample scaling made two expression profiles invalid due to unknown systematic issue,
# probabily a sample processing issue in the original study.
# This was revealed during processing, hence a reanlysis is performed after
# discarding the two invalid samples.

xx <- expr[, -1] %>%
  purrr::map(~(DescTools::RobScale(x = .x))) %>%
  bind_cols()
purrr::map(xx, median) %>% unlist() %>% is.na() %>% table()
# FALSE  TRUE
# 3734     2
purrr::map(xx, median) %>% unlist() %>% is.na() %>% which()
# GSM1635819 GSM1635820
# 3377       3378
clin %>%
  dplyr::filter(Sample_geo_accession %in% all_of(c("GSM1635819", "GSM1635820"))) %>%
  glimpse()
# These two samples are froma single study(Series), GSE66999.
# Both are ER and PR positive and HER2 negative.
# Neoadjuvant regimen.
# Treatment: Docetaxel + Epirubicine
# Core needle biopsy, preserved by snap freezing
# Profiled by Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (year 2008)
#   having 41108 probes

# Remove these samples from expr and clin !!!!!!!!!!!!!!!!!!

head2(clin) # 3736 80
head2(expr) # 9184 3737

clin <- clin %>%
  dplyr::filter(!(Sample_geo_accession %in% all_of(c("GSM1635819", "GSM1635820"))))
expr <- expr %>% dplyr::select(1, all_of(clin$Sample_geo_accession))

head2(clin) # 3734 81
head2(expr) # 9184 3735

#
# End of sample filtering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





# # mad sample scaling
# dlist$log2_mscale = bind_cols(
#   dlist$log2 %>% dplyr::select(1),
#   dlist$log2[, -1] %>%
#     purrr::map(~(DescTools::RobScale(x = .x))) %>%
#     bind_cols()
# )
# # mad sample scaling + global sample median centering
# purrr::map(dlist$log2_mscale[, -1], median) %>% unlist() %>% is.na() %>% table()
# # FALSE  TRUE
# # 3734     2
# purrr::map(dlist$log2_mscale[, -1], median) %>% unlist() %>% is.na() %>% which()
# # GSM1635819 GSM1635820
# # 3377       3378
# dlist$log2_mscale$GSM1635819 %>% is.na() %>% table()
# dlist$log2_mscale$GSM1635819 %>% is.infinite() %>% table()
# dlist$log2_mscale$GSM1635820 %>% is.na() %>% table()
# dlist$log2_mscale$GSM1635820 %>% is.infinite() %>% table()
# # both samples got all expression as either Inf or NaN
# mad(dlist$log2$GSM1635819) # 0
# mad(dlist$log2$GSM1635820) # 0
# range(dlist$log2$GSM1635819) # 3.818736 15.139882
# range(dlist$log2$GSM1635820) # 4.101677 15.445718
# median(dlist$log2$GSM1635819) # 3.818736 # Median is the lowest value !!!
# median(dlist$log2$GSM1635820) # 4.101677 # Median is the lowest value !!!
#
# (dlist$log2$GSM1635819 == median(dlist$log2$GSM1635819)) %>% table()
# # FALSE  TRUE
# # 4508  4676
# (dlist$log2$GSM1635820 == median(dlist$log2$GSM1635820)) %>% table()
# # FALSE  TRUE
# # 3526  5658
# # Roughly half the expression got lowest expression value
# # Seems like a systaematic issue, most probably an array processing issue.
# # Discard !!!!!!!!!!!!!




# 2) IHC subtyping, not depend on molecular data

table(is.na(clin$Hr)) # all NAs
# TRUE
# 3734

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
# 756 1230  850
table(is.na(clin$Subtype_ihc))
# FALSE  TRUE
# 2836   898





# 3) Module list

# ER/HER2/Immune signature
xx1 <- read_csv(file = "data/gene_module/Desmedt2008_Supplementary_Table_S1.csv")
module_list <- list()
for(i in  unique(xx1$Biological_process)){
  module_list[[str_c("desmedt2008_",i)]] <- xx1 %>%
    dplyr::filter(Biological_process == i) %>%
    dplyr::rename(Ncbi_gene_id = "EntrezGene.ID",
                  Hugo_gene_symbol = "HUGO.gene.symbol") %>%
    dplyr::mutate(Direction = if_else(coefficient < 0, -1, 1),
                  Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id),
                  Hugo_gene_symbol = as.character(Hugo_gene_symbol)) %>%
    dplyr::select(Ncbi_gene_id, Hugo_gene_symbol, Direction)
}
names(module_list) <- str_to_title(names(module_list))

# ER/HER2/Immune gene
module_list <- c(module_list,
                 list(
                   Er_gene = module_list$Desmedt2008_er[1, ],
                   Her2_gene = module_list$Desmedt2008_her2[1, ],
                   Immune_gene = module_list$Desmedt2008_immune_response[1, ]
                 ))



# PAM50 centroids
pam50 <- list()
pam50[["centroid"]] <- read_tsv(file = "data/pam50/PAM50/bioclassifier_R/pam50_centroids.txt")
pam50[["annot"]] <- read_tsv(file = "data/pam50/PAM50/bioclassifier_R/pam50_annotation.txt")

pam50 <- pam50$centroid %>%
  dplyr::rename(pcrID = "X1") %>%
  dplyr::left_join(pam50$annot, by = "pcrID") %>%
  dplyr::mutate(EntrezGene = str_c("ncbi_", EntrezGene)) %>%
  dplyr::filter(EntrezGene %in% expr$Ncbi_gene_id) %>%
  dplyr::rename(Ncbi_gene_id = "EntrezGene")
nrow(pam50) # 36


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
    x[, -1] %>% purrr::map_dfr(summary)
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







# 3) mdata list

mdlist <- list( # template !!!!!!!!!!

  # original
  log2 = NA,
  # combat sample scaling and mean centering
  log2_combat = NA,
  # benjamin: per dataset gene centering
  log2_center = NA,


  # genefu::rescale (quantile sample scaling)
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # quantile sample scaling
  log2_qscale = NA,
  # quantile sample scaling + global sample median centering
  log2_qscale_center = NA,
  # quantile sample scaling + combat sample centering
  log2_qscale_combat_m = NA,
  # quantile sample scaling + combat sample scaling and centering
  log2_qscale_combat = NA,


  # + per dataset gene rescaling using genefu (dgc)
  # quantile sample scaling
  log2_qscale_dgc = NA,
  # quantile sample scaling + global sample median centering
  log2_qscale_center_dgc = NA,
  # quantile sample scaling + combat sample centering
  log2_qscale_combat_m_dgc = NA,
  # quantile sample scaling + combat sample scaling and centering
  log2_qscale_combat_dgc = NA,


  # + pooled dataset gene rescaling using genefu (ggc)
  # quantile sample scaling
  log2_qscale_ggc = NA,
  # quantile sample scaling + global sample median centering
  log2_qscale_center_ggc = NA,
  # quantile sample scaling + combat sample centering
  log2_qscale_combat_m_ggc = NA,
  # quantile sample scaling + combat sample scaling and centering
  log2_qscale_combat_ggc = NA,


  # MAD sample rescaling
  # >>>>>>>>>>>>>>>>>>>>

  # NOTE: mad sample scaling + global sample median centering is irrelavant as
  # mad sample scaling will make all sample medians as zero.

  # mad sample scaling
  log2_mscale = NA,
  # mad sample scaling + combat sample centering
  log2_mscale_combat_m = NA,
  # mad sample scaling + combat sample scaling and centering
  log2_mscale_combat = NA,


  # + per dataset gene rescaling using genefu (dgc)
  # mad sample scaling
  log2_mscale_dgc = NA,
  # mad sample scaling + combat sample centering
  log2_mscale_combat_m_dgc = NA,
  # mad sample scaling + combat sample scaling and centering
  log2_mscale_combat_dgc = NA,


  # + pooled dataset gene rescaling using genefu (ggc)
  # mad sample scaling
  log2_mscale_ggc = NA,
  # mad sample scaling + combat sample centering
  log2_mscale_combat_m_ggc = NA,
  # mad sample scaling + combat sample scaling and centering
  log2_mscale_combat_ggc = NA
)


# log2 (original)
# expr

# qscale (genefu quantile sample scaling)
expr_qscale <- bind_cols(
  expr %>% dplyr::select(1),
  expr[, -1] %>%
    purrr::map(~(genefu::rescale(x = .x, q = 0.05))) %>%
    bind_cols()
)

# mscale (mad sample scaling)
expr_mscale <- bind_cols(
  expr %>% dplyr::select(1),
  expr[, -1] %>%
    purrr::map(~(DescTools::RobScale(x = .x))) %>%
    bind_cols()
)


# head2(expr)
# head2(expr_qscale)
# head2(expr_mscale)


# log2
x <- expr
mdlist$log2 <- update_clin(clin, x, module_list, pam50)

# log2 + combat sample scaling and mean centering
x <- bind_cols(
  expr %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = FALSE # default is FALSE
  ) %>%
    as_tibble()
)
mdlist$log2_combat <- update_clin(clin, x, module_list, pam50)

# benjamin: per dataset gene centering
x <- expr %>%
  dataset_gene_scaling(dataset_map = clin)
mdlist$log2_center <- update_clin(clin, x, module_list, pam50)




# genefu::rescale (quantile sample scaling)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# quantile sample scaling
x <- expr_qscale
mdlist$log2_qscale <- update_clin(clin, x, module_list, pam50)

# quantile sample scaling + global sample median centering
gm <- purrr::map(expr_qscale[, -1], median) %>% unlist() %>% median() # gm = 0.5269001
x <- bind_cols(
  expr_qscale %>% dplyr::select(1),
  expr_qscale[, -1] %>%
    purrr::map(function(x, gm){ x + (gm - median(x)) }, gm) %>%
    bind_cols()
)
mdlist$log2_qscale_center <- update_clin(clin, x, module_list, pam50)

# quantile sample scaling + combat sample centering
x <- bind_cols(
  expr_qscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_qscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = TRUE # default is FALSE
  ) %>%
    as_tibble()
)
mdlist$log2_qscale_combat_m <- update_clin(clin, x, module_list, pam50)

# quantile sample scaling + combat sample scaling and centering
x <- bind_cols(
  expr_qscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_qscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = FALSE # default is FALSE
  ) %>%
    as_tibble()
)
mdlist$log2_qscale_combat <- update_clin(clin, x, module_list, pam50)


# + per dataset gene rescaling using genefu (dgc)
# quantile sample scaling
x <- expr_qscale %>%
  dataset_gene_scaling(dataset_map = clin)
mdlist$log2_qscale_dgc <- update_clin(clin, x, module_list, pam50)

# quantile sample scaling + global sample median centering
gm <- purrr::map(expr_qscale[, -1], median) %>% unlist() %>% median() # gm = 0.5269001
x <- bind_cols(
  expr_qscale %>% dplyr::select(1),
  expr_qscale[, -1] %>%
    purrr::map(function(x, gm){ x + (gm - median(x)) }, gm) %>%
    bind_cols()
) %>%
  dataset_gene_scaling(dataset_map = clin)
mdlist$log2_qscale_center_dgc <- update_clin(clin, x, module_list, pam50)

# quantile sample scaling + combat sample centering
x <- bind_cols(
  expr_qscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_qscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = TRUE # default is FALSE
  ) %>%
    as_tibble()
) %>%
  dataset_gene_scaling(dataset_map = clin)
mdlist$log2_qscale_combat_m_dgc <- update_clin(clin, x, module_list, pam50)

# quantile sample scaling + combat sample scaling and centering
x <- bind_cols(
  expr_qscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_qscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = FALSE # default is FALSE
  ) %>%
    as_tibble()
) %>%
  dataset_gene_scaling(dataset_map = clin)
mdlist$log2_qscale_combat_dgc <- update_clin(clin, x, module_list, pam50)



# + pooled dataset gene rescaling using genefu (ggc)
# # quantile sample scaling
x <- expr_qscale %>%
  global_gene_scaling()
mdlist$log2_qscale_ggc <- update_clin(clin, x, module_list, pam50)

# quantile sample scaling + global sample median centering
gm <- purrr::map(expr_qscale[, -1], median) %>% unlist() %>% median() # gm = 0.5107782
x <- bind_cols(
  expr_qscale %>% dplyr::select(1),
  expr_qscale[, -1] %>%
    purrr::map(function(x, gm){ x + (gm - median(x)) }, gm) %>%
    bind_cols()
) %>%
  global_gene_scaling()
mdlist$log2_qscale_center_ggc <- update_clin(clin, x, module_list, pam50)

# quantile sample scaling + combat sample centering
x <- bind_cols(
  expr_qscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_qscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = TRUE # default is FALSE
  ) %>%
    as_tibble()
) %>%
  global_gene_scaling()
mdlist$log2_qscale_combat_m_ggc <- update_clin(clin, x, module_list, pam50)

# quantile sample scaling + combat sample scaling and centering
x <- bind_cols(
  expr_qscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_qscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = FALSE # default is FALSE
  ) %>%
    as_tibble()
) %>%
  global_gene_scaling()
mdlist$log2_qscale_combat_ggc <- update_clin(clin, x, module_list, pam50)


# MAD sample rescaling
# >>>>>>>>>>>>>>>>>>>>

# mad sample scaling
x <- expr_mscale
mdlist$log2_mscale = update_clin(clin, x, module_list, pam50)

# mad sample scaling + combat sample centering
x = bind_cols(
  expr_mscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_mscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = TRUE # default is FALSE
  ) %>%
    as_tibble()
)
mdlist$log2_mscale_combat_m = update_clin(clin, x, module_list, pam50)

# mad sample scaling + combat sample scaling and centering
x = bind_cols(
  expr_mscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_mscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = FALSE # default is FALSE
  ) %>%
    as_tibble()
)
mdlist$log2_mscale_combat = update_clin(clin, x, module_list, pam50)



# + per dataset gene rescaling using genefu (dgc)
# mad sample scaling
x <- expr_mscale %>%
  dataset_gene_scaling(dataset_map = clin)
mdlist$log2_mscale_dgc <- update_clin(clin, x, module_list, pam50)

# mad sample scaling + combat sample centering
x <- bind_cols(
  expr_mscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_mscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = TRUE # default is FALSE
  ) %>%
    as_tibble()
) %>%
  dataset_gene_scaling(dataset_map = clin)
mdlist$log2_mscale_combat_m_dgc <- update_clin(clin, x, module_list, pam50)

# mad sample scaling + combat sample scaling and centering
x <- bind_cols(
  expr_mscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_mscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = FALSE # default is FALSE
  ) %>%
    as_tibble()
) %>%
  dataset_gene_scaling(dataset_map = clin)
mdlist$log2_mscale_combat_dgc <- update_clin(clin, x, module_list, pam50)


# + pooled dataset gene rescaling using genefu (ggc)
# mad sample scaling
x <- expr_mscale %>%
  global_gene_scaling()
mdlist$log2_mscale_ggc <- update_clin(clin, x, module_list, pam50)

# mad sample scaling + combat sample centering
x <- bind_cols(
  expr_mscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_mscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = TRUE # default is FALSE
  ) %>%
    as_tibble()
) %>%
  global_gene_scaling()
mdlist$log2_mscale_combat_m_ggc <- update_clin(clin, x, module_list, pam50)

# mad sample scaling + combat sample scaling and centering
x <- bind_cols(
  expr_mscale %>% dplyr::select(1),
  sva::ComBat(
    dat = as.matrix(expr_mscale[, -1]),
    batch = clin$Series_matrix_accession,
    mean.only = FALSE # default is FALSE
  ) %>%
    as_tibble()
) %>%
  global_gene_scaling()
mdlist$log2_mscale_combat_ggc <- update_clin(clin, x, module_list, pam50)


# Selected preprocessings !!!!!!!!!!
# names(mdlist)
relevant_comparison = c(
  "log2", "log2_combat",  "log2_center", "log2_mscale_dgc"
)



# RLE plots: ALL
xx <- purrr::map(
  names(mdlist),
  function(nme, mdlist){
    bind_cols(
      tibble(Batch_correction_methode = nme),
      mdlist[[nme]] %>%
        dplyr::select(Series_matrix_accession, Sample_geo_accession,
                      "Rle_min", "Rle_first_qu", "Rle_median", "Rle_mean",
                      "Rle_third_qu", "Rle_max")
    )
  },
  mdlist
) %>%
  bind_rows()

p <- xx %>%
  tidyr::gather("Rle_class", "Rle_value", "Rle_first_qu","Rle_median", "Rle_third_qu") %>%
  dplyr::mutate(Series_matrix_accession = factor(Series_matrix_accession,
                                                 levels = unique(Series_matrix_accession)),
                Sample_geo_accession = factor(Sample_geo_accession,
                                              levels = unique(Sample_geo_accession)),
                Batch_correction_methode = factor(Batch_correction_methode, levels = unique(Batch_correction_methode))) %>%
  ggplot(aes(x = Sample_geo_accession, y = Rle_value, group = Rle_class)) +
  geom_line(aes(color = Series_matrix_accession)) +
  guides(color = "none") +
  labs(x = "Samples", y = "RLE") +
  facet_grid(facets = Batch_correction_methode ~ 1, scales = "free_y") +
  theme(legend.position = "bottom",
        axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, #size = 5,
                                    margin = margin(0,.75,0,.75, "cm")))

pdf2(file = str_c(outdir, "rle_all.pdf"))
print(p)
dev.off()


# RLE plots: Relevant
xx <- purrr::map(
  relevant_comparison,
  function(nme, mdlist){
    bind_cols(
      tibble(Batch_correction_methode = nme),
      mdlist[[nme]] %>%
        dplyr::select(Series_matrix_accession, Sample_geo_accession,
                      "Rle_min", "Rle_first_qu", "Rle_median", "Rle_mean",
                      "Rle_third_qu", "Rle_max")
    )
  },
  mdlist
) %>%
  bind_rows()

p <- xx %>%
  tidyr::gather("Rle_class", "Rle_value", "Rle_first_qu","Rle_median", "Rle_third_qu") %>%
  dplyr::mutate(Series_matrix_accession = factor(Series_matrix_accession,
                                                 levels = unique(Series_matrix_accession)),
                Sample_geo_accession = factor(Sample_geo_accession,
                                              levels = unique(Sample_geo_accession)),
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

pdf2(file = str_c(outdir, "rle_relevant.pdf"))
print(p)
dev.off()




# ER/HER2 IHC status prdiction
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>
xx <- purrr::map(
  mdlist,
  function(x){
    list(er_sig = pROC::roc(formula = Er ~ Desmedt2008_er, data = x),
         her2_sig = pROC::roc(formula = Her2 ~ Desmedt2008_her2, data = x),
         er_gene = pROC::roc(formula = Er  ~ Er_gene, data = x),
         her2_gene = pROC::roc(formula = Her2 ~ Her2_gene, data = x))
  }
)



# All comparison
purrr::map_dfr(
  names(xx),
  function(nme, xx){
    x = xx[[nme]]
    list(type = nme,
         er_sig = auc(x$er_sig) %>% as.numeric(),
         her2_sig = auc(x$her2_sig) %>% as.numeric(),
         er_gene = auc(x$er_gene) %>% as.numeric(),
         her2_gene = auc(x$her2_gene) %>% as.numeric())
  }, xx) %>%
  as.data.frame()
#                        type    er_sig  her2_sig   er_gene her2_gene
# 1                      log2 0.8266707 0.5158213 0.7912057 0.6059456
# 2               log2_combat 0.8659527 0.6787362 0.8803074 0.7290693
# 3               log2_center 0.8713231 0.6335382 0.8946938 0.7762707
# 4               log2_qscale 0.8897351 0.7185182 0.8704406 0.7826218
# 5        log2_qscale_center 0.8883630 0.8015831 0.8785903 0.7857920
# 6      log2_qscale_combat_m 0.8852806 0.7362623 0.8931958 0.7714326
# 7        log2_qscale_combat 0.8670841 0.6764160 0.8801898 0.7298033
# 8           log2_qscale_dgc 0.8699530 0.6338772 0.8941928 0.7858419
# 9    log2_qscale_center_dgc 0.8735654 0.6347407 0.8927380 0.7755157
# 10 log2_qscale_combat_m_dgc 0.8699530 0.6338772 0.8941928 0.7858419
# 11   log2_qscale_combat_dgc 0.8699530 0.6338772 0.8941928 0.7858419
# 12          log2_qscale_ggc 0.8864809 0.7078944 0.8704406 0.7826218
# 13   log2_qscale_center_ggc 0.8862289 0.7862115 0.8785903 0.7857920
# 14 log2_qscale_combat_m_ggc 0.8826360 0.7189139 0.8931958 0.7714326
# 15   log2_qscale_combat_ggc 0.8642486 0.6661578 0.8801898 0.7298033
# 16              log2_mscale 0.8852414 0.7986424 0.8789322 0.7860389
# 17     log2_mscale_combat_m 0.8854597 0.7162002 0.8947043 0.7628577
# 18       log2_mscale_combat 0.8663730 0.6792137 0.8795716 0.7284482
# 19          log2_mscale_dgc 0.8730009 0.6366874 0.8921448 0.7704629
# 20 log2_mscale_combat_m_dgc 0.8730009 0.6366874 0.8921448 0.7704629
# 21   log2_mscale_combat_dgc 0.8730009 0.6366874 0.8921448 0.7704629
# 22          log2_mscale_ggc 0.8825407 0.7825287 0.8789322 0.7860389
# 23 log2_mscale_combat_m_ggc 0.8821656 0.7033049 0.8947043 0.7628577
# 24   log2_mscale_combat_ggc 0.8626688 0.6693808 0.8795716 0.7284482



# Relevant comaprison
purrr::map_dfr(
  relevant_comparison,
  function(nme, xx){
    x = xx[[nme]]
    list(type = nme,
         er_sig = auc(x$er_sig) %>% as.numeric(),
         her2_sig = auc(x$her2_sig) %>% as.numeric(),
         er_gene = auc(x$er_gene) %>% as.numeric(),
         her2_gene = auc(x$her2_gene) %>% as.numeric())
  }, xx) %>%
  as.data.frame()
#              type    er_sig  her2_sig   er_gene her2_gene
# 1            log2 0.8266707 0.5158213 0.7912057 0.6059456
# 2     log2_combat 0.8659527 0.6787362 0.8803074 0.7290693
# 3     log2_center 0.8713231 0.6335382 0.8946938 0.7762707
# 4 log2_mscale_dgc 0.8730009 0.6366874 0.8921448 0.7704629





# PAM50 vs Subtype-IHC
# >>>>>>>>>>>>>>>>>>>>

# Kappa,κ ranges from 0 to 1, with 0 indicating no relation and 1 indicating
# a perfect concordance. Typically qualitative descriptions are associated
# with intervals [κ ≤ 0.20, slight agreement; 0.20 < κ ≤ 0.40, fair agreement;
# 0.40 < κ ≤ 0.60, moderate agreement, 0.60 < κ ≤ 0.80, substantial agreement;
# and 0.80 < κ ≤ 0.99, almost perfect agreement
# See Benjamins Three gene model paper

table(mdlist$log2$Subtype_ihc) # Reference
# HER2   HR   TN
# 756  1230  850
table(mdlist$log2$Subtype_pam50)
# Basal         Her2         LumA         LumB       Normal Unclassified
# 1089          639         1139          491          341           35
table(mdlist$log2$Subtype_pam50_2) # Test subtype
# HER2   HR   TN
# 639 1630 1089



# All comparison
xx <- mdlist
purrr::map_dfr(
  names(xx),
  function(nme, xx){
    x = xx[[nme]]
    list(type = nme,
         kappa = irr::kappa2(ratings = x %>%
                               dplyr::select(Subtype_ihc, Subtype_pam50_2))$value
    )
  }, xx = mdlist) %>%
  as.data.frame()
#                        type     kappa
# 1                      log2 0.6208237
# 2               log2_combat 0.6173302
# 3               log2_center 0.6340461
# 4               log2_qscale 0.6470775
# 5        log2_qscale_center 0.6428843
# 6      log2_qscale_combat_m 0.6255837
# 7        log2_qscale_combat 0.6184341
# 8           log2_qscale_dgc 0.6396011
# 9    log2_qscale_center_dgc 0.6306087
# 10 log2_qscale_combat_m_dgc 0.6396011
# 11   log2_qscale_combat_dgc 0.6396011
# 12          log2_qscale_ggc 0.6372117
# 13   log2_qscale_center_ggc 0.6254652
# 14 log2_qscale_combat_m_ggc 0.6335823
# 15   log2_qscale_combat_ggc 0.6123502
# 16              log2_mscale 0.6355030
# 17     log2_mscale_combat_m 0.6302831
# 18       log2_mscale_combat 0.6137476
# 19          log2_mscale_dgc 0.6276156
# 20 log2_mscale_combat_m_dgc 0.6276156
# 21   log2_mscale_combat_dgc 0.6276156
# 22          log2_mscale_ggc 0.6232505
# 23 log2_mscale_combat_m_ggc 0.6366000
# 24   log2_mscale_combat_ggc 0.6135366



# Relevant comparison
xx <- mdlist
purrr::map_dfr(
  relevant_comparison,
  function(nme, xx){
    x = xx[[nme]]
    list(type = nme,
         kappa = irr::kappa2(ratings = x %>%
                               dplyr::select(Subtype_ihc, Subtype_pam50_2))$value
    )
  }, xx = mdlist) %>%
  as.data.frame()
#              type     kappa
# 1            log2 0.6208237
# 2     log2_combat 0.6173302
# 3     log2_center 0.6340461
# 4 log2_mscale_dgc 0.6276156



# ER - immune response association boxplot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# All comparison
xx <- mdlist
xx_stat <- purrr::map_dfr(
  names(xx),
  function(nme, xx){

    x = xx[[nme]]
    wilcox_sig <- wilcox.test(x$Desmedt2008_immune_response ~ factor(x$Er))
    wilcox_gene <- wilcox.test(x$Immune_gene ~ factor(x$Er))

    list(type = nme,
         p_sig = wilcox_sig$p.value,
         p_gene = wilcox_gene$p.value,
         label_sig = str_c("p ", format(wilcox_sig$p.value, digits = 2, scientific = T)),
         label_gene = str_c("p ", format(wilcox_gene$p.value, digits = 2, scientific = T))
    )
  }, xx = mdlist)

xx_dat <- purrr::map(
  names(xx),
  function(nme, xx){

    x = xx[[nme]]

    bind_cols(
      tibble(type = nme),
      x %>% dplyr::select(Desmedt2008_immune_response, Immune_gene, Er)
    )
  }, xx = mdlist) %>%
  bind_rows()

p <- xx_dat %>%
  dplyr::mutate(type = factor(type, levels = unique(type))) %>%
  ggplot(aes(x = Er, y = Desmedt2008_immune_response)) +
  geom_boxplot() +
  geom_text(data = xx_stat, aes(x = 2, y = +Inf, label = label_sig), vjust = 1) +
  labs(x = "ER status (IHC)", y = "Immune module score") +
  facet_wrap(facets = ~ type, scales = "free_y", ncol = 4) +
  theme(
    strip.text = element_text(size = 8)
  )

pdf2(file = str_c(outdir, "er_immsig_boxplot_all.pdf"))
print(p)
dev.off()


p <- xx_dat %>%
  dplyr::mutate(type = factor(type, levels = unique(type))) %>%
  ggplot(aes(x = Er, y = Immune_gene)) +
  geom_boxplot() +
  geom_text(data = xx_stat, aes(x = 2, y = +Inf, label = label_gene), vjust = 1) +
  labs(x = "ER status (IHC)", y = "Immune gene (STAT1)") +
  facet_wrap(facets = ~ type, scales = "free_y", ncol = 4) +
  theme(
    strip.text = element_text(size = 8)
  )

pdf2(file = str_c(outdir, "er_immgene_boxplot_all.pdf"))
print(p)
dev.off()



# Relevant comparison
xx <- mdlist
xx_stat <- purrr::map_dfr(
  relevant_comparison,
  function(nme, xx){

    x = xx[[nme]]
    wilcox_sig <- wilcox.test(x$Desmedt2008_immune_response ~ factor(x$Er))
    wilcox_gene <- wilcox.test(x$Immune_gene ~ factor(x$Er))

    list(type = nme,
         p_sig = wilcox_sig$p.value,
         p_gene = wilcox_gene$p.value,
         label_sig = str_c("p ", format(wilcox_sig$p.value, digits = 2, scientific = T)),
         label_gene = str_c("p ", format(wilcox_gene$p.value, digits = 2, scientific = T))
    )
  }, xx = mdlist)

xx_dat <- purrr::map(
  relevant_comparison,
  function(nme, xx){

    x = xx[[nme]]

    bind_cols(
      tibble(type = nme),
      x %>% dplyr::select(Desmedt2008_immune_response, Immune_gene, Er)
    )
  }, xx = mdlist) %>%
  bind_rows()

p <- xx_dat %>%
  dplyr::mutate(type = factor(type, levels = unique(type))) %>%
  ggplot(aes(x = Er, y = Desmedt2008_immune_response)) +
  geom_boxplot() +
  geom_text(data = xx_stat, aes(x = 2, y = +Inf, label = label_sig), vjust = 1) +
  labs(x = "ER status (IHC)", y = "Immune module score") +
  facet_wrap(facets = ~ type, scales = "free_y", ncol = 4) +
  theme(
    strip.text = element_text(size = 8)
  )

pdf2(file = str_c(outdir, "er_immsig_boxplot_relevant.pdf"))
print(p)
dev.off()


p <- xx_dat %>%
  dplyr::mutate(type = factor(type, levels = unique(type))) %>%
  ggplot(aes(x = Er, y = Immune_gene)) +
  geom_boxplot() +
  geom_text(data = xx_stat, aes(x = 2, y = +Inf, label = label_gene), vjust = 1) +
  labs(x = "ER status (IHC)", y = "Immune gene (STAT1)") +
  facet_wrap(facets = ~ type, scales = "free_y", ncol = 4) +
  theme(
    strip.text = element_text(size = 8)
  )

pdf2(file = str_c(outdir, "er_immgene_boxplot_relevant.pdf"))
print(p)
dev.off()


#
# ==============================================================================




# Backup code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 5) Module list

# # ER/HER2(Desmedt), Immune/Ifn(Gruosso)
# xx1 <- read_csv(file = "data/gene_module/Desmedt2008_Supplementary_Table_S1.csv")
# xx2 <- read_csv(file = "data/gene_module/Gruosso2019_TIME.csv")
#
# module_list <- list()
# for(i in  unique(xx1$Biological_process)){
#   module_list[[str_c("desmedt2008_",i)]] <- xx1 %>%
#     dplyr::filter(Biological_process == i) %>%
#     dplyr::rename(Ncbi_gene_id = "EntrezGene.ID",
#                   Hugo_gene_symbol = "HUGO.gene.symbol") %>%
#     dplyr::mutate(Direction = if_else(coefficient < 0, -1, 1),
#                   Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id),
#                   Hugo_gene_symbol = as.character(Hugo_gene_symbol)) %>%
#     dplyr::select(Ncbi_gene_id, Hugo_gene_symbol, Direction)
# }
# for(i in  unique(xx2$Module_name)){
#   module_list[[str_c("gruosso2019_", i)]] <- xx2 %>%
#     dplyr::filter(Module_name == i) %>%
#     dplyr::mutate(Direction = if_else(Direction < 0, -1, 1),
#                   Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id),
#                   Hugo_gene_symbol = as.character(Hugo_gene_symbol)) %>%
#     dplyr::select(Ncbi_gene_id, Hugo_gene_symbol, Direction)
# }
# names(module_list) <- str_to_title(names(module_list))
#





