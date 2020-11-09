# s1_manual_geo_search_result_curation.R

# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
# Curate the search result text file from geo computationally and
# convert it into a table strucutre(geo_series_summary) for easy manual annotation.
# Manually annotate "geo_series_summary" with information related to
# regimen/arm, sample procrument method, sample archival method,
# whether selected for further processing etc. for each series matrix.


# Script strucutre
# >>>>>>>>>>>>>>>>
# 1. Curate computationally, the geo search results exported from geo as a text file.
# 2. Write out the curated geo_search_result for the manual annotation with regimen,
# sample procrument method, archival method etc.
# 3. Load and save the curated and annotated series matrix summary.




# 1. Curate computationally, the geo search results exported from geo as a text file.
# ==============================================================================

xx <- read_delim(
  file = "data/geo_search_taxane.antracyclin.breast_31july2020_SearchResult.txt",
  delim = "^",
  col_names = FALSE)

geo_series_summary <- tibble(
  Title = NA,
  Description = xx$X1[str_detect(xx$X1, "(Submitter supplied)")],
  Organism = xx$X1[str_detect(xx$X1, "Organism:\t")],
  Type = xx$X1[str_detect(xx$X1, "Type:")],
  Platform = xx$X1[str_detect(xx$X1, "Platform") & str_detect(xx$X1, "Samples")],
  Sample_size = NA,
  Ftp = xx$X1[str_detect(xx$X1, "FTP download:")],
  Series_accession = xx$X1[str_detect(xx$X1, "Series\t\tAccession:")],
  Id = NA
)


# xx$X1 = "1. Drug-induced change in gene expression across NCI-60 cell line"
# The following code will identify the no.of serieses identified by geo search
yy <- str_split_fixed(xx$X1, pattern = "\\. ", n = 2)[,1]
yy <- as.integer(yy)
table(is.na(yy))
# FALSE  TRUE
# 172  1032
geo_series_summary$Title <- xx$X1[!is.na(yy)]


geo_series_summary <- geo_series_summary %>%
  dplyr::mutate(
    Title = str_split_fixed(Title, "\\. ", n = 2)[,2] %>% str_trim(),
    Description = str_replace(Description, "\\(Submitter supplied\\)", "") %>% str_trim(),
    Organism = str_replace(Organism, "Organism:", "") %>% str_trim(),
    Type = str_replace(Type, "Type:", "") %>% str_trim(),
    Platform = str_replace(Platform, "Platform:", "") %>% str_trim(),
    Platform = str_replace(Platform, "Platforms:", "") %>% str_trim(),
    Sample_size = str_split(Platform, " ") %>%
      purrr::map_int(~(.x[length(.x)-1] %>% as.integer())),
    Platform = purrr::map2_chr(
      Platform,
      Sample_size,
      ~(str_split_fixed(.x, str_c(.y, " S"), 2)[, 1])) %>% str_trim(),
    Ftp = str_split_fixed(Ftp, "ftp:", 2)[,2] %>% str_trim(),
    Ftp = str_c("ftp:", Ftp),
    Id = str_split_fixed(Series_accession, "ID:", 2)[,2] %>% str_trim(),
    Series_accession = str_split_fixed(Series_accession, "ID:", 2)[,1] %>% str_trim(),
    Series_accession = str_split_fixed(Series_accession, "Accession:", 2)[,2] %>% str_trim(),
  )

glimpse(geo_series_summary)
sum(geo_series_summary$Sample_size) # 37497

#
# ==============================================================================




# 2. Write out the curated geo_search_result for the manual annotation with regimen,
# sample procrument method, archival method etc.
# ==============================================================================

# Write out "geo_series_summary" for manual annotation !!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write_tsv(
  x = geo_series_summary,
  path = str_c(
    outdir,
    "geo_search_taxane.antracyclin.breast_31july2020_SearchResult_curated.tsv"
  )
)
# !!!!!! Manual annotation !!!!!!!!!!!!!!
# Annotate "geo_series_summary" manually with information related to
# regimen/arm, sample procrument method, sample archival method,
# whether selected for further processing etc. Save the contents as
# geo_series_summary_taxane.antracyclin.breast_31july2020_SearchResult .tsv(.ods)

#
# ==============================================================================




# 3. Load and save the curated and annotated series matrix summary.
# ==============================================================================

# Load annotated geo_series_summary !!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
geo_series_summary <- read_tsv(
  file = "results/main/geo_search_taxane.antracyclin.breast_31july2020_SearchResult_annotated.tsv"
) %>%
    dplyr::rename_all(~(str_c("Series_", .x) %>%
                          str_to_sentence() %>%
                          str_replace("Series_series","Series")))


glimpse(geo_series_summary)

# Saving
# save(geo_series_summary , file = str_c(outdir, "geo_series_summary.RData"))

#
# ==============================================================================







