# functions.R

unique_names <- function(x){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # The functions ensures that the extracted sample characteristics names
  # are unique. Redundancy is corrected by adding 1,2,3 etc as suffix to
  # redundant characteristics names.

  # input
  # >>>>>
  #
  # x: Character vector of column names (can contain redundant column name)

  # output
  # >>>>>>
  #
  # A character vector of unique column name

  x_count <- table(x)

  for(ix in names(x_count)){

    if (x_count[ix] == 1) {
      next
    } else {
      x[which(x == ix)] <- str_c(ix, seq(1,x_count[ix]),sep = ".")
    }

  }

  return(x)
}



get_row_range <- function(file, n_max = 250){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # The series matrix consolidates multilevel data into a single file.
  # The multiple levels were series, sample and expression data.
  # The expression data is readily distinguishable by data read functions
  # (such as "read_csv") when the comment option is used, as all other metadata
  # (series and sample) in the series matrix start with the character "!".
  # Hence by setting comment == "!", the expression data can be read automatically.
  # However there is no way to distinguish series and sample data, except
  # a blank line that separates them. Hence to automate reading in series and
  # sample data, the skip and n_max options of data read functions (eg. read_csv)
  # must be utilized. This script identifies the correct row ranges (i.e. no.of
  # rows to skip and maximum no of rows to read) of series and
  # sample data, by exploiting the blank line that separates them.

  # input
  # >>>>>
  #
  # file: Character path to compressed series matrix.
  #       To feed into read_tsv() function.
  # n_max: An increment to be feed into read_tsv() function to read the
  #       compressed series matrix. The idea of n_max is to ensure that
  #       minimal data is read into R (and no expression data) to  extract
  #       row indices of series and sample data. Reading of entire series matrix
  #       costs performance issue and is not necessary to read in all expression
  #       data to extract row ranges of series or sample data.
  #       If the row ranges cannot be found in first n_max lines
  #       of data, then the algorithm repeats by doubling n_max. This will go on
  #       until n_max reaches 2000 lines, a hard cutoff to prevent indefinite
  #       looping. It is highly unlikely that a series matrix may contain >1000
  #       sample characteristics (+series data).

  # output
  # >>>>>>
  #
  # List of no.of lines to skip ("skip" option) and maximum line to read
  # ("n_max" option) to extract series and sample data.




  # While loop is to accommodate the unexpected scenario of
  # -more than 250 lines of series and sample characteristics.

  while (n_max < 2000) {

    suppressWarnings(
      xx <- read_tsv(file = file,
                     col_names = FALSE,
                     col_types = "c",
                     comment = "",
                     skip = 0,
                     n_max = n_max,
                     skip_empty_rows = FALSE)
    )


    # To get row indices to distinguish series and sample data

    index <- c(
      which(is.na(xx$X1)), # first empty raw is read as NA
      which(xx$X1 == "!series_matrix_table_begin")
    )

    if (length(index) == 2) {
      break
    } else {
      n_max <- n_max + 250
    }

  }


  if (length(index) == 2) {

    skip <- c(
      "series" = 0, # series data starts from row 1
      "sample" = index[1] # sample data starts (index[1] + 1)
    )

    n_max <- c(
      "series" = index[1] - 1,
      "sample" = (index[2] - 1) - index[1]
    )


  } else {

    # if indices couldn't find after increasing n_max

    warning("\nRow indices to distinguish series and sample data could not find. Try reviewing the respective series matrix manually. Not reading the series and sample data.")

    skip <- c(
      "series" = 0,
      "sample" = 0
    )

    n_max <- c(
      "series" = 0,
      "sample" = 0
    )

  }

  return(list(skip = skip, n_max = n_max))

}


get_var_name <- function(x){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # The sample characteristics extracted from series matrix usually have
  # the column name "Sample_characteristics_ch1" (or"_ch2" if dual channel).
  # The actual characteristic name is encoded in the value as
  # "characteristic-name:characteristics-value". This function extracts
  # the characteristics name. The function exploits the format of colon
  # separated characteristics name and value pair for extracting characteristics name.



  # input
  # >>>>>
  #
  # x: A tibble of sample characteristics to extract the characteristics names

  # output
  # >>>>>>
  #
  # A character vector of extracted names


  nme <- x %>%
    purrr::map_chr(
      function(xx){

        # colon separate characteristic name and value
        no_colon <- xx %>% str_detect(":", negate = TRUE) %>% all()
        if(is.na(no_colon)){ # empty cell value generate NA
          no_colon = TRUE
        }

        if (no_colon) {

          # unknown column name
          # characteristic-value seperatoer is missing
          return("X")

        } else {

          xx <- str_split_fixed(xx, ":", 2)[, 1]
          # above code return "" for missing/NA cell values
          xx <- xx[!xx == ""]

          # The first occuring characteristic name as column name
          return(unique(xx)[1])

          # Error log
          # >>>>>>>>>

          # # The logic of picking the max occurred characteristic name as column name,
          # # seems complicated; hence removed.
          # xx <- table(xx)
          # # expecting the max occurrence as column name
          # # This is to accommodate inconsistency arise due to missingness in some columns,
          # # as this will cause the column values to shift left and ending up in entirely
          # #  different column.
          # return(names(xx)[xx == max(xx)])

          # End of error log
          # >>>>>>>>>>>>>>>>

        }
      })

  # # making unknown column names unique
  # idx <- which(nme == "X")
  # if(length(idx) > 1){
  #   nme[idx] <- str_c("X",1:length(idx))
  # }

  return(unique_names(nme)) # unique names append .123 for redundant names
}


get_var_qc <- function(x){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # The sample characteristics extracted from series matrix usually have
  # the column name "Sample_characteristics_ch1" (or"_ch2" if dual channel).
  # The actual characteristic name is encoded in the value as
  # "characteristic-name:characteristics-value". This function checks
  # the expected behaviors and if it fails raises a flag for that
  # characteristics column. The checks performed were:
  # 1) The presence of ":" that separates characteristics name and value.
  # 2) The presence of unique characteristics name.
  # The reason for the second check is that if any of the value for a
  # sample characteristic is missing and is not encoded as NAs, the entire
  # row (samples on rows and characteristic on columns) towards the right of the
  # missing value shifts to the left n no.of times, where n is the no.of missing
  #  values that is not encoded as NAs. This problem creates curation issues and
  #  consumes lot of time if do it manually for all series matrices.
  #  With the inclusion of variable qc checks, the curator can focus only on
  #  the flagged datasets/columns.


  # input
  # >>>>>
  #
  # x: A tibble of sample characteristics to extract the characteristics names

  # output
  # >>>>>>
  #
  # A logical vector of characteristics QC values.



  qc <- x %>%
    purrr::map_lgl(
      function(xx){

        # colon seperate characteristic name and value
        no_colon <- xx %>% str_detect(":", negate = TRUE) %>% all()
        if(is.na(no_colon)){ # empty cell value generate NA
          no_colon = TRUE
        }

        if (no_colon) {

          # unknown column name
          # characteristic-value seperatoer is missing
          # qc failed
          return(FALSE)

        } else {

          return(
            str_split_fixed(xx, ":", 2)[, 1] %>%
              unique() %>%
              length() == 1
          )
        }

      })

  return(qc)
}


format_geo <- function(files){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # This script extracts series, sample and expression data from series matrix file.
  # Curate the sample characteristics computationally and unexpected behaviors
  # in sample characteristics will be flagged for manual intervention.


  # input
  # >>>>>
  #
  # file: Character path to compressed series matrix. To feed into read_tsv() function.


  # output
  # >>>>>>
  #
  # Output is a list containing series, sample, expression,
  # curated sample characteristics, and a flag denoting whether there are
  # issues in sample characteristics (Yes means there are issues).



  data_list <- purrr::map(

    files,

    function(file){

      print(file)

      # Get the no.of lines to skip and read to extract series and sample data
      row_range <- get_row_range(file = file)

      data <- list()


      cat("Getting series data.\n")

      data[["series"]] <- read_tsv(

        file = file,
        col_names = FALSE,
        col_types = cols(),
        skip = row_range$skip["series"],
        n_max = row_range$n_max["series"]

      ) %>%
        # Series and sample data are in wide format
        t() %>%
        as.data.frame() %>%
        as_tibble()

      names(data$series) <- data$series[1, ] %>%
        t() %>%
        str_replace(., "!", "")

      data$series <- as_tibble(data$series[-1, ],
                               .name_repair = "unique") %>%
        dplyr::mutate_if(is.factor, ~as.character(.x))


      cat("Getting sample data.\n")

      data[["sample"]] <- read_tsv(

        file = file,
        col_names = FALSE,
        col_types = cols(),
        skip = row_range$skip["sample"],
        n_max = row_range$n_max["sample"]

      ) %>%
        # Series and sample data are in wide format
        t() %>%
        as.data.frame() %>%
        as_tibble()

      names(data$sample) <- data$sample[1, ] %>%
        t() %>%
        str_replace(., "!", "")

      data$sample <- as_tibble(data$sample[-1, ],
                               .name_repair = "unique") %>%
        dplyr::mutate_if(is.factor, ~as.character(.x))



      cat("Getting expression data.\n")

      data[["expression"]] <- read_tsv(

        file = file,
        col_names = TRUE,
        col_types = cols(),
        comment = "!"

      )



      cat("Formating sample characteristics.\n")

      # Error log !!!!!!!!!!!
      # >>>>>>>>>>>>>>>>>>>>>

      # Old code (version:0)
      # Error:
      # Doesn't accomodate dual channel array.
      #
      # data[["sample_characteristics"]] <- data$sample  %>%
      #   dplyr::select(str_which(names(data$sample),
      #                           "Sample_characteristics_ch1"))

      # Old code (version:1)
      # Error:
      # For dual channel array, it is not sure that which channel contains test
      # samples.
      # i.e ch2 may contains reference sample and ch1 test sample, and vice-versa
      #
      # # New code accomodating dual channel arrays
      # if(any(str_detect(names(data$sample), "Sample_characteristics_ch2"))){
      #   # accomodating dual channel
      #   data[["sample_characteristics"]] <- data$sample  %>%
      #     dplyr::select(str_which(names(data$sample),
      #                             "Sample_characteristics_ch2"))
      # } else {
      #   data[["sample_characteristics"]] <- data$sample  %>%
      #     dplyr::select(str_which(names(data$sample),
      #                             "Sample_characteristics_ch1"))
      # }

      # End of error log !!!!!!!!!!!
      # >>>>>>>>>>>>>>>>>>>>>>>>>>>>


      # New code extracteing sample characteristics from all channels in
      # dual channel arrays
      data[["sample_characteristics"]] <- data$sample  %>%
        dplyr::select(str_which(names(data$sample),
                                "Sample_characteristics"))


      # The sample characterisitcs extracted from series matrix usually have
      # the column name "Sample_characterisitcs_ch1" (or"_ch2" if dual channel).
      # The actual charcaterisitc name is encoded in the value as
      # "characterisitc-name:characterisicts-value".
      # The get_var_name() function extracts the characterisitics name.

      var_name <- get_var_name(data$sample_characteristics)

      # To check whether manual intervention is needed for sample characterisitcs
      # curation.
      var_qc <- get_var_qc(data$sample_characteristics)

      names(data$sample_characteristics) <- var_name


      data[["sample_characteristics"]] <- purrr::pmap_dfr(

        list(x = data[["sample_characteristics"]],
             y = var_name,
             z = var_qc),

        function(x, y, z){

          if (z) {

            # Log of previously fixed error !!!
            # str_split_fixed(x, str_c(y,": "), 2)[, 2]
            # If the column name "y" contains special character, eg "(",
            # the above code will not parse correctly

            # Extracting sample cgaracterisitc value
            # Each column in x is exprected to be of the strucutre
            # "chracteristic-name:characterisitc-value"

            str_split_fixed(x, ":", 2)[, 2] %>% str_trim()

          } else {
            # If variable qc is failed, no curation is performed.
            x
          }

        }

      )



      data[["dataset_issues"]] <- if_else(any(!var_qc), "Yes", "No")
      # Yes - there are issues

      data

    } # end map function

  ) # end map

  names(data_list) <- files
  data_list

} # end main function


clean_names <- function(x){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # Cleans the input character vector, by trimming and removing space.

  # input
  # >>>>>
  #
  # x: Character vector of column names or others

  # output
  # >>>>>>
  #
  # A character vector of cleaned names

  x %>%
    str_trim() %>%
    str_replace_all(" ", "_") %>%
    str_to_sentence()
}


# updated t_tibble()
t_tibble <- function(x, names_x_desc = NULL){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # Tranposes a gene expression tibble.
  # By definition, row names are not recommended for tibble, instead any
  # rowname must be included as a column in the tibble. For instance the
  # as_tibble() function has "rownames" option to give names for existing rownames
  # so as to include the rownames as the first column of a tibble. However,
  # rownames are conventional and usually included in a tibble as the first column.
  # In this case, the default transpose function may not work properly as
  # transposing a tibble may cause data type conflicts and may convert the
  # entire column of a transposed tibble as the with the type of first column of
  # non-transposed tibble. This will be case with gene expression tibbles with
  # first column as gene symbols, and column names are sample names.
  # Note that this function is defined to use with gene expression tibbles.
  # Gene expression tibbles have the identical data type for all columns except
  # the first column representing rownames.

  # input
  # >>>>>
  #
  # x: A tibble to transpose. The 1st column is considered as rownames.
  # names_x_desc: A character string that describes what names(x) represents.
  #               This character string will be used as the name of the first
  #               column of the trassposed tibble.

  # output
  # >>>>>>
  #
  # A transposed tibble


  # Note: t(x) will set all data to single type

  # If names_x_desc is NULL, it is set as the name of the 1st column of x
  if (is.null(names_x_desc)) {
    names_x_desc <- names(x)[1]
  }

  # Extract rownames of x. This will later used to set column names of t(x)
  rownames_x <- x %>% dplyr::select(1) %>% tibble::deframe()


  tx <- t(x %>% dplyr::select(-1))
  colnames(tx) <- rownames_x

  return(as_tibble(tx, rownames = names_x_desc))
}


# updated get_max_var_annot()
get_max_var_annot <- function(ge, xannot, pset_colname, gene_colname){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # Identifies the maximum Vientiane probe/probeset representing a single gene and
  # its annotation.

  # input
  # >>>>>
  #
  # ge: Gene expression tibble with probe/probeset id as the first column and
  #     column names as sample names (names(ge)).The name of the probe/probeset id column
  #     must be the value of "pset_colname".
  # xannot: Annotation tibble with probe/probeset id and gene id columns named
  #           exactly as the value of "pset_colname" and "gene_colname".
  # pset_colname: Name of column containing probe/probeset ids,
  #                 must present in both ge and xannot.
  # gene_colname: Name of column containing gene ids,
  #                 must present in both ge and xannot.

  # output
  # >>>>>>
  #
  # A tibble representing annotation of maximum variant prode/probeset.


  # Converting all relevant ids to character
  ge <- ge %>%
    dplyr::mutate_at(pset_colname, ~as.character(.x))
  xannot <- xannot %>%
    dplyr::mutate_at(c(pset_colname, gene_colname), ~as.character(.x))


  # print(table((xannot %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe()) %in% (ge %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe())))


  pset <- intersect(ge %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe(),
                    xannot %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe())

  # # ***debugging
  # print(c(length(pset),nrow(xannot)))
  # print(length(pset)/nrow(xannot))


  # Filter and sort @ by
  # Ref: https://dplyr.tidyverse.org/articles/programming.html
  # Ref: https://stackoverflow.com/questions/26497751/pass-a-vector-of-variable-names-to-dplyr::arrange-in-dplyr
  # Both base::as.name(by) and rlang::sym(by) will take the strings as input and
  # return the value encoded by the string as symbols which can be used as
  # column names
  # is.symbol(by) # FALSE
  # is.symbol(as.name(by)) # TRUE
  # is.symbol(sym(by)) # TRUE
  ge <- ge %>%
    # dplyr::mutate_at(pset_colname, ~as.character(.x)) %>%
    dplyr::filter_at(pset_colname, dplyr::all_vars((. %in% pset))) %>%
    dplyr::arrange(!! as.name(pset_colname))

  xannot <- xannot %>%
    # dplyr::mutate_at(pset_colname, ~as.character(.x)) %>%
    dplyr::filter_at(pset_colname, dplyr::all_vars((. %in% pset))) %>%
    dplyr::arrange(!! as.name(pset_colname))

  # Ref !!: https://community.rstudio.com/t/eval-vs-bang-bang-in-functions-using-dplyr/3977/2



  # Pset variance vector
  pset_var <- ge %>%
    t_tibble() %>% # Transposing gene expression tibble (genes on column)
    dplyr::select(-1) %>% # Removing rownames (samplenames)
    purrr::map_dbl(var)

  # Unique gene ids
  ugid <- xannot %>%
    dplyr::select(all_of(gene_colname)) %>%
    dplyr::distinct_at(gene_colname, .keep_all = FALSE) %>%
    tibble::deframe()

  xannot_gene = xannot %>% dplyr::select(gene_colname) %>% tibble::deframe()

  # # ***debugging
  # print(
  #   c(identical(ge %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe(),
  #               xannot %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe()),
  #     identical(names(pset_var),
  #               xannot %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe())) %>%
  #     all()
  # )
  # # Note: names(pset_var) and xannot_gene are in order !!!!!!!!!!!!!!!!!!!!!!!
  # identical(names(pset_var), xannot[, pset_colname] %>% tibble::deframe())
  # TRUE


  # Max var pset per unique gene id
  max_var_pset <-  purrr::map_chr(
    ugid,
    function(gene, pset_var, xannot_gene ){

      # Per gene annot data frame
      # Note: names(pset_var) and xannot_gene are in order !!!!!!!!!!!!!!!!!!!!
      gene_pset_idx  <- which(xannot_gene == gene)
      gene_pset_var <- pset_var[gene_pset_idx]

      # Get max var pset
      names(gene_pset_var)[gene_pset_var == max(gene_pset_var)][1] #if identical maximums select 1st one
    },
    pset_var,
    xannot_gene
  )

  max_var_annot <- tibble(pset = max_var_pset)
  names(max_var_annot) = c(get("pset_colname"))

  max_var_annot <- left_join(max_var_annot,
                             xannot,
                             by = pset_colname)
  return(max_var_annot)

}


head2 <- function(x, n = 5){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # head() will display first 6 rows of a dataframe and all columns irrespective of
  # no.of columns by default. This function will control it by letting the user
  # specify it. Along with a glimpse of a dataframe, this function also prints
  # the dimension. This function is defined to aid in exploratory analysis.
  # Also see the glimpse() function available from tidyverse package

  # input
  # >>>>>
  #
  # x: A tibble or dataframe
  # n: The no.of columns and rows to printout.

  # output
  # >>>>>>
  #
  # A list of subset of the input (head) and its dimension (dim).


  ridx <- min(n, nrow(x))
  cidx <- min(n, ncol(x))
  list(head = x[1:ridx, 1:cidx], dim = dim(x))
}


tiff2 <- function(file = "Rplot.tiff", width = 7.5, height = 7.5){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # To set useful defaults for tiff()

  # input
  # >>>>>
  #
  # file: Character string representing file name to save the figure.
  #        This input will be feed into tiff() function.
  # width: Width of the figure in inches.
  #         This input will be feed into tiff() function.
  # height: Height of the figure in inches.
  #          This input will be feed into tiff() function.

  # output
  # >>>>>>
  #
  # A tiff file device to plot into.

  tiff(filename = file, width = width, height = height, res = 150, units = "in")
}


pdf2 <- function(file = "Rplot.pdf", width = 7.5, height = 7.5, onefile = TRUE){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # To set useful defaults for pdf()

  # input
  # >>>>>
  #
  # file: Character string representing file name to save the figure.
  #        This input will be feed into tiff() function.
  # width: Width of the figure in inches.
  #         This input will be feed into tiff() function.
  # height: Height of the figure in inches.
  #          This input will be feed into tiff() function.
  # onefile: A flg to be feed into pdf() function. TRUE means multiple plots
  #           will be appended to the same pdf file.

  # output
  # >>>>>>
  #
  # A pdf file device to plot into.

  pdf(file = file, width = width, height = height, onefile = onefile)
}


get_module_score <- function(x, module_list, by = "Entrez_Gene"){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # Compute gene module score as a weighted average

  # input
  # >>>>>
  #
  # x: A tibble with genes on column and samples on rows.
  #     1st column is considered as sample names.
  # module_list: List of gene modules. Each gene module must contain at least
  #               two columns; a gene id column with column name exactly as the
  #               value of the "by" input and a direction column with
  #               column name as "Direction". Further, the gene id
  #               in the gene module must match the gene ids of x, ie. column
  #               names of x. The Direction in the gene module represents the
  #               gene's association with the phenotype that the
  #               gene module represents. For instance, for a gene module
  #               representing TP53 mutation srtatus a direction of +1 suggests
  #               up regulated in mutant phenotype and -1 represents
  #               down regulated in mutant phenotype.
  # by: The column name of gene module which contain gene id. The type gene id
  #     in the gene module must match that of column names of x.

  # output
  # >>>>>>
  #
  # A tibble of module score. Column names represnts module names and the 1st
  # column contains sample ids.


  sig_score <- purrr::map_dfr( # map_dfr
    module_list,
    function(sig, xdata, by){

      sig <- sig %>%
        dplyr::filter_at(
          by,
          dplyr::all_vars(. %in% names(xdata))
        )


      xdata <- xdata %>%
        dplyr::select(
          c(names(xdata)[1],
            sig %>% dplyr::select(by) %>% tibble::deframe())
        )


      score <- (xdata %>% dplyr::select(-1) %>% as.matrix()) %*% (sig %>% dplyr::select(Direction) %>% as.matrix())

      # score/nrow(sig) # this old code generates matrix
      as.numeric(score/nrow(sig))
    },
    xdata = x,
    by = by
  )

  sig_score <- dplyr::bind_cols(x %>% dplyr::select(1),
                                sig_score)

}


dataset_gene_scaling <- function(x, dataset_map){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # Performs a per dataset quntile gene rescaling using rescale() function
  # available from genefu R package. The rescaling is performed gene wise
  # (i.e. row wise)with q = 0.05 as the input to rescale function.
  # This methode has been reported to correct dataset effect from pooled dataset.

  # input
  # >>>>>
  #
  # x: A tibble of pooled dataset with genes on rows and samples
  #     on column. The first column contains gene ids.
  #     This function is defined to used with pooled GEO dataset, hence
  #     sample id used was Sample_geo_accession and dataset id used was
  #     Series_matrix_accession.
  # dataset_map: A tibble containing the mapping between sample id and dataset id.
  #               i.e, the mapping between Sample_geo_accession and
  #               Series_matrix_accession.

  # output
  # >>>>>>
  #
  # A tibble of per dataset rescaled gene expression values (gene wise/row wise).


  # make x and dataset_map congruent
  x <- x %>% dplyr::select(1, all_of(dataset_map$Sample_geo_accession))

  # print(identical(names(x)[-1], dataset_map$Sample_geo_accession))

  x_scaled <- purrr::map(
    unique(dataset_map$Series_matrix_accession),
    function(nme, x, dataset_map){
      idx <- which(dataset_map$Series_matrix_accession == nme)
      # print(idx)
      apply(x[,idx+1], 1, function(x){genefu::rescale(x, q = 0.05)}) %>%
        t() %>%
        as_tibble()
    },
    x, dataset_map
  )

  # return
  bind_cols(
    x[, 1],
    bind_cols(x_scaled)
  )

}


global_gene_scaling <- function(x){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # Performs quantile gene rescaling in the pooled dataset using rescale() function
  # available from genefu R package. The dataset id is ignored in this function and
  # the scaling performed on the pooled dataset.The rescaling is performed gene wise
  # (i.e. row wise) with q = 0.05 as the input to rescale function.

  # input
  # >>>>>
  #
  # x: A tibble of pooled dataset with genes on rows and samples
  #     on column. The first column contains gene ids.
  #     This function is defined to used with pooled GEO dataset, hence
  #     sample id used was Sample_geo_accession and dataset id used was
  #     Series_matrix_accession.

  # output
  # >>>>>>
  #
  # A tibble of rescaled gene expression values (gene wise / row wise).


  x_scaled <- apply(x[, -1], 1, function(x){genefu::rescale(x, q = 0.05)}) %>%
    t() %>%
    as_tibble()

  # return
  bind_cols(
    x[, 1],
    x_scaled
  )

}


explore_na <- function(xx){


  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # Function to explore NAs in an expression matrix per gene and per sample

  # input
  # >>>>>
  #
  # xx: expresison matrix tibble; first column ID_REF contains probeset ids

  print("sample")
  idx_sample <- purrr::map_lgl(xx, ~(any(is.na(.x))))
  print("gene")
  idx_gene <- purrr::map_lgl(t(xx) %>% as_tibble(), ~(any(is.na(.x))))

  return(list(na_sample = idx_sample, na_gene = idx_gene))
}


get_patient_summary <- function(xx1){

  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  #
  # Function to generate patient summary


  # input
  # >>>>>
  #
  # xx1: tibble of patient characteristics

  bind_rows(

    # Age
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "Age<=50"),
        (xx1$Age <= 50) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "Age<=50",
              "cat" = "NA",
              n = xx1$Age %>% is.na() %>% sum())
    ),

    # Grade
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "Grade"),
        xx1$Grade %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "Grade",
              "cat" = "NA",
              n = xx1$Grade %>% is.na() %>% sum())
    ),


    # Node
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "Node"),
        (xx1$Node_bin) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "Node",
              "cat" = "NA",
              n = xx1$Node_bin %>% is.na() %>% sum())
    ),

    # Size
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "Size"),
        (xx1$Size_cat) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "Size",
              "cat" = "NA",
              n = xx1$Size_cat %>% is.na() %>% sum())
    ),


    # ER
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "ER"),
        (xx1$Er) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "ER",
              "cat" = "NA",
              n = xx1$Er %>% is.na() %>% sum())
    ),

    # PR
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "PR"),
        (xx1$Pr) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "PR",
              "cat" = "NA",
              n = xx1$Pr %>% is.na() %>% sum())
    ),

    # HR
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "HR"),
        (xx1$Hr) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "HR",
              "cat" = "NA",
              n = xx1$Hr %>% is.na() %>% sum())
    ),

    # HER2
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "HER2"),
        (xx1$Her2) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "HER2",
              "cat" = "NA",
              n = xx1$Her2 %>% is.na() %>% sum())
    ),

    # Subtype IHC
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "Subtype IHC"),
        (xx1$Subtype_ihc) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "Subtype IHC",
              "cat" = "NA",
              n = xx1$Subtype_ihc %>% is.na() %>% sum())
    ),

    # Subtype PAM50
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "Subtype PAM50"),
        (xx1$Subtype_pam50) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "Subtype PAM50",
              "cat" = "NA",
              n = xx1$Subtype_pam50 %>% is.na() %>% sum())
    ),


    # pCR
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "pCR"),
        xx1 %>% group_by(cat = Response) %>% summarise(n = n())
        # The below code will throw error when all Response == NA, if regimen = adj
        # (xx1$Response) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "pCR",
              "cat" = "NA",
              n = xx1$Response %>% is.na() %>% sum())
    ) %>% na.omit(),

    # DFS
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "DFS"),
        (xx1$Event_dfs) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "DFS",
              "cat" = "NA",
              n = xx1$Event_dfs %>% is.na() %>% sum())
    ) %>%
      dplyr::mutate(cat = cat %>%
                      str_replace("0", "no-event") %>%
                      str_replace("1", "event")),

    # Arm_chemo
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "Arm_chemo"),
        (xx1$Arm_chemo) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "Arm_chemo",
              "cat" = "NA",
              n = xx1$Arm_chemo %>% is.na() %>% sum())
    ),

    # Arm_her2
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "Arm_her2"),
        (xx1$Arm_her2) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "Arm_her2",
              "cat" = "NA",
              n = xx1$Arm_her2 %>% is.na() %>% sum())
    ),

    # Arm_hormone
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "Arm_hormone"),
        (xx1$Arm_hormone) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "Arm_hormone",
              "cat" = "NA",
              n = xx1$Arm_hormone %>% is.na() %>% sum())
    ),

    # Arm_other
    bind_rows(
      # non-NAs
      bind_cols(
        tibble(variable = "Arm_other"),
        (xx1$Arm_other) %>% table() %>% as_tibble() %>% dplyr::rename(cat = ".")
      ),
      # NAs
      tibble( variable = "Arm_other",
              "cat" = "NA",
              n = xx1$Arm_other %>% is.na() %>% sum())
    )
  )
}



