
beta_mat <- NULL
labs_ss_tib <- NULL
for (curDir in mergeDirs_vec) {
  
  # Find Sample Sheet::
  curr_ss_csv <- findFileByPattern(dir=curDir, patter='_AutoSampleSheet.csv.gz$', max=1, recursive=FALSE, verbose=opt$verbose)
  curr_fn_csv <- curr_ss_csv %>% stringr::str_replace('_AutoSampleSheet.csv.gz$','_MergedDataFiles.tib.csv.gz')
  base_dir <- base::dirname(curr_ss_csv)
  stopifnot(file.exists(curr_ss_csv), file.exists(curr_fn_csv))
  
  cat(glue::glue("[{par$prgmTag}]:{TAB} [{outName}] Found sample_csv={curr_ss_csv}.{RET}") )
  cat(glue::glue("[{par$prgmTag}]:{TAB} [{outName}] Found call_table={curr_fn_csv}.{RET}") )
  
  # Load and filter::
  curr_ss_tib <- suppressMessages(suppressWarnings( readr::read_csv(curr_ss_csv) )) %>% 
    dplyr::filter(!!opt$samplePvalName > !!opt$samplePvalPerc) %>% 
    dplyr::filter(!!class_var %in% trainClass_vec)
  
  calls_path_tib <- suppressMessages(suppressWarnings( readr::read_csv(curr_fn_csv) ))
  
  betas_path_tib <- calls_path_tib %>% dplyr::filter(Method %in% c(betaKey))
  pvals_path_tib <- calls_path_tib %>% dplyr::filter(Method %in% c(pvalKey))
  
  stopifnot(base::nrow(betas_path_tib)==1)
  stopifnot(base::nrow(pvals_path_tib)==1)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #            Load Beta/Pval And Merge into Previous Matrix::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  beta_csv <- file.path(base_dir, base::basename(betas_path_tib$Full_Path[1]) )
  pval_csv <- file.path(base_dir, base::basename(pvals_path_tib$Full_Path[1]) )
  
  beta_mat <- loadCallsMatrix(betaCSV=beta_csv, pvalCSV=pval_csv, minPval=pvalMin, mat=beta_mat, 
                              cgn=NULL, ss=curr_ss_tib,
                              verbose=opt$verbose, tc=1, tt=cTracker)
  
  labs_ss_tib <- labs_ss_tib %>% dplyr::bind_rows(curr_ss_tib)
  
  beta_mat %>% dim() %>% print()
  cat(glue::glue("[{par$prgmTag}]:{TAB} [{outName}] Done. {outName}.{RET}{RET}") )
  
  # break
}

# Try sorting Labeled/Merged Sample Sheet by Sample_Class and then re-order matrix by new order::
sort_ss_tib <- labs_tib %>% dplyr::arrange(!!class_var) %>% 
  dplyr::mutate(!!class_var := as.factor(!!class_var),
                !!class_idx := as.integer(!!class_var)-1) %>% 
  dplyr::mutate(!!class_idx := as.integer(!!class_idx) ) %>%
  dplyr::select(Sentrix_Name, !!class_var, !!class_idx) # %>% as.data.frame()

# QC Sanity Check:: Make sure the new ordering works::
if (FALSE) {
  cbind(
    beta_mat %>% colnames(),
    sort_ss_tib %>% dplyr::pull(Sentrix_Name),
    beta_mat[ , dplyr::pull(sort_ss_tib, Sentrix_Name)] %>% colnames()
  ) %>% tibble::as_tibble() %>% dplyr::filter(V2==V3)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Build Raw and Imputed Sorted Matricies::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

beta_raw_matT <- beta_mat[ , dplyr::pull(sort_ss_tib, Sentrix_Name)]

mod_rds <- '/Users/bbarnes/Documents/CustomerFacing/build_models/Sample_Name/BETA-DELTA-Decoder/i-beta_i-poob-0.9/dml-1000-Ivana-145/seed-13/alpha-0.1/glmnet-0.1.model.rds'
mod <- readr::read_rds(mod_rds)




# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             Investigate K-means Clustering between DML/dBL::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD: Lots more to do here!!!
#  - Compare COVIC (C0) probes: performace and rank..
#  - Much more...
#
if (FALSE) {
  kmean_dml_tib <- full_dml_tib %>% dplyr::select(Probe_ID, Pr_t, Rank)
  kmean_dbl_tib <- full_dblMu_tib %>% dplyr::select(Probe_ID, !!dbl_sortMu_key, !!dbl_sdNeg_key, !!dbl_sdPos_key, Rank)
  
  kmean_tib <- dplyr::inner_join(kmean_dml_tib, kmean_dbl_tib, by="Probe_ID", suffix=c("_dml", "_dbl"))
  
  kmean_tib %>% dplyr::filter(!is.na(nSARSCov2_pSARSCov2_CSS_mu)) %>% 
    dplyr::filter(!is.na(nSARSCov2_beta_sd)) %>% 
    dplyr::filter(!is.na(pSARSCov2_beta_sd)) %>% dplyr::arrange(Rank_dbl)
  
  gg <- ggplot2::ggplot(data=kmean_tib, aes(x=Pr_t, y=nSARSCov2_pSARSCov2_CSS_mu)) +
    ggplot2::geom_point() +
    ggplot2::geom_density2d()
  
  #
  # dBL:: Q2 instead of mu::
  #
  kmean_dblQ2_tib <- full_dblMu_tib %>% dplyr::select(Probe_ID, !!dbl_sortQ2_key, !!dbl_sdNeg_key, !!dbl_sdPos_key, Rank)
  kmean_tibQ2 <- dplyr::inner_join(kmean_dml_tib, kmean_dblQ2_tib, by="Probe_ID", suffix=c("_dml", "_dbl"))
  
  gg <- ggplot2::ggplot(data=kmean_tibQ2, aes(x=Pr_t, y=nSARSCov2_pSARSCov2_CSS_q2)) +
    ggplot2::geom_point() +
    ggplot2::geom_density2d()
}


if (FALSE) {
  
  minPval <- 0.9
  cgn <- NULL
  ss <- NULL
  
  class_var <- rlang::sym("Sample_Name")
  class_idx <- rlang::sym("Class_Idx")
  sentrix_name <- rlang::sym("Sentrix_Name")
  
  topDir <- '/Users/bbarnes/Documents/CustomerFacing'
  beta_csv <- file.path(topDir,'merge_builds/Sample_Class/DELTA-8x1-EPIC-Core/EPIC/B4/DELTA-8x1-EPIC-Core_EPIC_B4_i_beta.raw-data.csv.gz')
  poob_csv <- file.path(topDir,'merge_builds/Sample_Class/DELTA-8x1-EPIC-Core/EPIC/B4/DELTA-8x1-EPIC-Core_EPIC_B4_i_poob.raw-data.csv.gz')
  samp_csv <- file.path(topDir,'merge_builds/Sample_Class/DELTA-8x1-EPIC-Core/EPIC/B4/DELTA-8x1-EPIC-Core_EPIC_B4_AutoSampleSheet.csv.gz')
  samp_csv <- '/Users/bbarnes/Documents/CustomerFacing/build_models/Sample_Name/BETA-DELTA-Decoder/i-beta_i-poob-0.9/Sample_Name_BETA-DELTA-Decoder_i-beta_i-poob-0.9.ClasSampleSheet.sorted.csv.gz'
  
  samp_tib <- suppressMessages(suppressWarnings( readr::read_csv(samp_csv) ))
  sort_ss_tib <- samp_tib %>% dplyr::arrange(!!class_var) %>% 
    dplyr::mutate(!!class_var := as.factor(!!class_var),
                  !!class_idx := as.integer(!!class_var)-1) %>% 
    dplyr::mutate(!!class_idx := as.integer(!!class_idx) ) %>%
    dplyr::select(!!sentrix_name, !!class_var, !!class_idx) # %>% as.data.frame()
  
  beta_tib <- suppressMessages(suppressWarnings( readr::read_csv(beta_csv) ))
  if (!is.null(cgn)) beta_tib <- beta_tib %>% dplyr::inner_join(cgn, by="Probe_ID")
  if (!is.null(ss)) beta_tib <- beta_tib %>% dplyr::select('Probe_ID', ss$Sentrix_Name)
  
  pval_tib <- suppressMessages(suppressWarnings( readr::read_csv(poob_csv) ))
  if (!is.null(cgn)) pval_tib <- pval_tib %>% dplyr::inner_join(cgn, by="Probe_ID")
  if (!is.null(ss)) pval_tib <- pval_tib %>% dplyr::select('Probe_ID', ss$Sentrix_Name)
  
  beta_df <- beta_tib %>% tibble::column_to_rownames(var="Probe_ID") %>% as.data.frame()
  pval_df <- pval_tib %>% tibble::column_to_rownames(var="Probe_ID") %>% as.data.frame()
  
  beta_mat <- beta_df %>% as.matrix()
  pval_mat <- pval_df %>% as.matrix()
  
  pval_idx <- which(pval_mat > minPval)
  fail_cnt <- length(pval_idx)
  
  beta_org_vec <- beta_mat[pval_idx]
  pval_org_vec <- pval_mat[pval_idx]
  
  beta_masked_mat <- beta_mat
  beta_masked_mat[pval_idx] <- NA
  masked_cnt <- which(is.na(beta_masked_mat)) %>% length()
  
  nan_tib <- beta_masked_mat %>% getMaskedTib(verbose=10)
  
  beta_impute_mat1 <- beta_masked_mat %>% as.data.frame() %>% makeX_glmnet_imputeNA(na.impute = TRUE)
  beta_impute_mat2 <- impute_matrix_mean(beta_masked_mat)

  #
  # Do imputation by Class::
  #
  mat_samples_vec <- dplyr::select(beta_tib, -Probe_ID) %>% names()
  sort_ss_tib <- sort_ss_tib %>% dplyr::filter(Sentrix_Name %in% mat_samples_vec)
  # split_ss_list <- sort_ss_tib %>% split(.$Sample_Name)
  sort_ss_names <- sort_ss_tib %>% dplyr::distinct(!!class_var) %>% dplyr::pull(!!class_var) %>% as.vector()
  
  beta_masked_mat3 <- NULL
  # for (sName in names(split_ss_list)) {
  for (sName in sort_ss_names) {
    cat(glue::glue("Class={sName}.{RET}"))
    # cur_ss_tib <- split_ss_list[[sName]]
    cur_ss_tib <- sort_ss_tib %>% dplyr::filter(!!class_var == sName)
    
    cur_masked_mat <- beta_masked_mat[ , dplyr::pull(cur_ss_tib, !!sentrix_name) ]
    cur_impute_mat <- impute_matrix_mean(mat=cur_masked_mat)
    cur_impute_mat %>% dim() %>% print()
    
    if (is.null(beta_masked_mat3)) {
      beta_masked_mat3 <- cur_impute_mat
    } else {
      beta_masked_mat3 <- cbind(beta_masked_mat3, cur_impute_mat[match(rownames(beta_masked_mat3), rownames(cur_impute_mat)), ] )
    }
    beta_masked_mat3 %>% dim() %>% print()
  }
  # beta_impute_mat3 <- impute_matrix_mean(beta_masked_mat3)
  
  beta_imp_vec1 <- beta_impute_mat1[pval_idx]
  beta_imp_vec2 <- beta_impute_mat2[pval_idx]
  beta_imp_vec3 <- beta_impute_mat3[pval_idx]
  
  # cpg_avg1_tib <- beta_masked_mat %>% rowMeans(na.rm=TRUE) %>% tibble::enframe(name="Probe_ID","beta_mean")
  
  full_tib <- tibble::tibble(
    Masked_Index=pval_idx,
    Masked_Beta=beta_org_vec,
    Masked_Pval=pval_org_vec,
    Impute_Beta1=beta_imp_vec1,
    Impute_Beta2=beta_imp_vec2,
    Impute_Beta3=beta_imp_vec3
  )
  
  nan_tib %>% dplyr::inner_join(full_tib, by=c("idx"="Masked_Index"))
}

# End of file
