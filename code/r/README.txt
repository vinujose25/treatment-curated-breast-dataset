R script structure
==================

When running scripts, keep the order of the scripts in mind, as any artifacts generated in one script are used in other scripts following it. The scripts are not fully automated and require manual intervention while executing.

1. initialize.R
    a. The main script from which other scripts are managed

2. functions.R
    a. Supporting functions for other scripts.

3. s1_geo_search_result_curation.R
    a. Compile the GEO search results.

4. s2_geo_dataset_extraction.R
    a. Extract the sample characteristics and expression matrix from the series matrix.
    b. Format sample characteristics both computationally and manually.

5. s3_sample_characterisitic_curation.R
    a. Manual curation and annotation of sample characteristics names, and re-coding of characteristics values for consistency.

6. s4_expression_data_integration.R
    a. NA filtering, converting to log2 space, maximum probeset collapse, and data integration.

7. s5_clinical_data_integration.R
    a. Selection of relevant clinical-pathological variables, their integration, and the curation of treatment regimens.

8. s6_process_and_validate_integrated_expression_matrix.R
    a. Comparison of methods to correct dataset effect from the pooled gene expression data
      i. Non-corrected (log2)
      ii. Combat (log2 + combat)
      iii. Quantile (log2 + per dataset quantile gene-rescaling)
      iv. Madquantile (log2 + MAD sample-rescaling + per dataset quantile gene-rescaling)
    b. The comparison is made by comparing
      i. RLE plots
      ii. IHC vs. PAM50 molecular subtype agreement (Cohen’s Kappa)
      iii. The ability of  the ER/HER2 gene to predict the respective IHC status (AUC).

9. s7_figures_tables_data.R
    a. Generate journal-specific figures and tables.

