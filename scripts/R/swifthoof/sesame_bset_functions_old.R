


if (FALSE) {
  
  #
  #
  toDo_str= "Infinium I background and HMM-Meth Alignment"
  #
  #
  #  - Rebuild all manifests::
  #     - Include "Full_ID" Probe_ID in Product Manifests
  #     - Reduce all uneeded fields (DESIGN, etc.)
  #  - Signal based ML:: 
  #     - Signal File Output:: Make sure sigs are swapped after inference
  #     - multi-track HMM-meth alignments 
  #
  #  - Add pre-val appending
  #  - 
  
  dplyr::bind_cols(rdat$ssheet_tib, rdat$cur_list$sums_ssh) %>%
    dplyr::select(-starts_with('Iscan_'), starts_with('Iscan_')) %>% 
    dplyr::mutate_if(is.numeric, list(round), 4) %>% getSsheetDataTab(
      minOobPval=min_pval_vec[1], minOobPerc=min_perc_vec[1], 
      minNegPval=min_pval_vec[2], minNegPerc=min_perc_vec[2], verbose = 10)
  
  
  ssheet_tib <- dplyr::bind_cols(rdat$ssheet_tib, rdat$cur_list$sums_ssh) %>%
    dplyr::select(-starts_with('Iscan_'), starts_with('Iscan_')) %>% 
    dplyr::mutate_if(is.numeric, list(round), 4)
  
  ann0_tib <- getSsheetDescTib(
    tib = ssheet_tib, 
    minOobPval = 0.1, minOobPerc = 90, 
    minNegPval = 0.02, minNegPerc = 98, 
    verbose = 10
  )
  
  tmp_ann_tab <- getSsheetDataTab(
    tib = ssheet_tib, 
    minOobPval = 0.1, minOobPerc = 90, 
    minNegPval = 0.02, minNegPerc = 98, 
    verbose = 10)
  
  
  ann0_tib <- getSsheetDescTib(
    tib = ssheet_tib, 
    minOobPval = 0.1, minOobPerc = 90, 
    minNegPval = 0.02, minNegPerc = 98, 
    verbose = 10
  )
  
  ann0_tib %>% gather(Variable, Description) %>% print(n=382)
  rdat$cur_list$sums_ssh %>% gather(Variable, Value) %>% print(n=288)
  
  
  ann_full_tib <- getSsheetDataTab(
    tib=ssheet_tib,
    minOobPval = 0.1, minOobPerc = 90, 
    minNegPval = 0.02, minNegPerc = 98, 
    verbose = 10)
  
  
  #
  # This WORKS!!!!
  #
  report_vec <- c("Requeue","pass_perc","mean")
  
  sum_tab <- rdat$cur_list$sums_dat %>% 
    tidyr::unite(key, Probe_Type,Probe_Design,Metric, sep='_') %>% 
    tidyr::gather(metric, value, -key, -Workflow_key, -Workflow_idx) %>% 
    dplyr::filter(!is.na(value)) %>% 
    dplyr::filter(metric %in% report_vec) %>%
    tidyr::unite(key, key,metric, sep='_') %>% 
    dplyr::mutate(key=paste(key,Workflow_idx, sep='_'))
  
  sum_tab %>% print(n=base::nrow(sum_tab))
  
  sum_tab %>%
    
  
  sum_tab %>%
    # dplyr::group_by(Workflow_key, Workflow_idx) %>%
    tidyr::pivot_wider(id_cols=c(Workflow_key, Workflow_idx), 
                       names_from="key", values_from="value") %>%
    dplyr::select(dplyr::contains("_Requeue_"), 
                  dplyr::contains("_pass_perc_"),
                  dplyr::everything())
  
  #
  # END OF WORKING!!!
  #

  
  #
  # Requeue Flags and Infinium I background calculation
  #  - Split Infinium I/II look at Poob detection p-values
  
  verbose <- opt$verbose
  vt <- 1
  tc <- 0
  tTracker <- NULL
  
  raw_sset <- rdat$raw_sset
  call_dat <- rdat$raw_list$call_dat
  sigs_dat <- rdat$raw_list$sigs_dat
  
  # Split Analytical and Controls::
  #
  raw_call_tib <- call_dat %>%
    expand_ids(verbose=verbose,tc=tc+1,tt=tTracker)
  raw_sigs_tib <- sigs_dat %>%
    expand_ids(verbose=verbose,tc=tc+1,tt=tTracker)
  
  # Order of Operations::
  #  - [i, n, d]
  
  # Requeue
  
  keys <- c('pOOBAH','PnegEcdf')
  mins <- c(opt$minOobPerc, opt$minNegPerc)
  
  raw_sam_tib <- rdat$raw_list$sam_sheet
  req_ssh_tib <- addRequeueFlags(tib=raw_sam_tib, keys=keys, mins=mins,
                                 verbose=verbose,tc=tc+1,tt=tTracker)


  rdat$raw_list$beta$beta_sum
  rdat$raw_list$pval$pOOBAH$pval_sum
  rdat$raw_list$pval$PnegEcdf$pval_sum
  rdat$raw_list$sigs$sigs_sum
  
  
  # Calls Summary:: pvals
  #
  funcTag <- 'testing'
  tabsStr <- '\t'
  by="Probe_ID"
  type="Probe_Type"
  des="Probe_Design"
  
  minNegPval <- 0.1
  percision_pval <- 6
  percision_beta <- 4
  percision_sigs <- 2
  
  verbose <- opt$verbose
  vt <- 1
  tc <- 1
  tt <- NULL
  
  man_tib  <- rdat$sman
  # raw_sset <- rdat$raw_list$sset_dat
  raw_sset <- rdat$raw_sset
  
  #
  # LEFT OFF HERE
  #
  pvalN_tib <- ssetToTib(sset=raw_sset, source='pvals', name='PnegEcdf', 
                         sort=TRUE,
                         percision=6, verbose=4)
  
  pvalN_sum_tib <- ssetTibToSummary(
    tib=pvalN_tib,man=man_tib,
    pval=0.02, perc=98,
    percision=percision_pval,
    save=FALSE, csv=NULL,
    by=by, type=type, des=des,
    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  
  pvalP_tib <- ssetToTib(sset=raw_sset, source='pvals', name='pOOBAH', 
                         sort=TRUE,
                         percision=6, verbose=4)
  
  pvalP_sum_tib <- ssetTibToSummary(
    tib=pvalP_tib,man=man_tib,
    pval=0.1, perc=90,
    percision=percision_pval,
    save=FALSE, csv=NULL,
    by=by, type=type, des=des,
    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  
  betas_tib <- ssetToTib(sset=raw_sset, source='betas', 
                         sort=TRUE,
                         percision=6, verbose=4)
  
  betas_sum_tib <- ssetTibToSummary(
    tib=betas_tib,man=man_tib,
    percision=percision_beta,
    save=FALSE, csv=NULL,
    by=by, type=type, des=des,
    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)

  sigsM_tib <- ssetToTib(sset=raw_sset, source='sigs', man=man_tib,
                         sort=TRUE,
                         percision=6, verbose=4)
  
  sigsM_sum_tib <- ssetTibToSummary(
    tib=sigsM_tib,
    percision=percision_sigs,
    save=FALSE, csv=NULL,
    by=by, type=type, des=des,
    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  
  
  
  
  ssh_cur_tab <- 
    dplyr::bind_rows(
      pvalP_sum_tib,
      pvalN_sum_tib,
      betas_sum_tib,
      sigsM_sum_tib
    ) %>% 
    tidyr::unite(Probe_Group, Probe_Type,Probe_Design,Metric, sep='_') %>% 
    tidyr::gather(Metric, value, -Probe_Group) %>%
    dplyr::filter(!is.na(value)) %>%
    tidyr::unite(Group_Metric, Probe_Group,Metric, sep='_')
  
  ssh_cur_tib <- ssh_cur_tab %>%
    tidyr::spread(Group_Metric, value)
  
  #
  # 
  #
  betas_sum_tib
  pvalP_sum_tib
  pvalN_sum_tib
  sigsM_sum_tib
  
  rdat$raw_list$beta$beta_sum
  rdat$raw_list$pval$pOOBAH$pval_sum
  rdat$raw_list$pval$PnegEcdf$pval_sum
  rdat$raw_list$sigs$sigs_sum
  
  
  
  
  
  
  
  
  #
  # Gather/Summary investigation::
  #
  rdat$raw_list$sigs$sigs_dat %>% 
    dplyr::group_by(Probe_Type, Probe_Design) %>%
    summarise_if(is.numeric, list(min=min, median=median, mean=mean, sd=sd, max=max), na.rm=TRUE) %>%
    dplyr::ungroup()
  
  rdat$raw_list$sigs$sigs_dat %>% 
    dplyr::select(-Probe_ID) %>%
    tidyr::gather(Metric, value, -Probe_Type, -Probe_Design) %>%
    dplyr::group_by(Probe_Type, Probe_Design, Metric) %>%
    summarise_if(is.numeric, list(min=min, median=median, mean=mean, sd=sd, max=max), na.rm=TRUE) %>%
    dplyr::ungroup()
  
  
    rdat$raw_list$pval$pOOBAH$pval_dat %>%
      dplyr::select(-Probe_ID) %>%
      # tidyr::gather(-Probe_Type, -Probe_Design)
      # dplyr::group_by(Probe_Type, Probe_Design) %>%
      tidyr::gather(Metric, value, -Probe_Type, -Probe_Design)
    
  
  

  tmp_sum_tib <- ssetToPvalSheet(sset=raw_sset, man=man_tib, 
                                 key="pOOBAH", min=0.1, perc=90,
                                 percision=6, verbose=4)
  
  cur_pval_ssheet <- 
    sigsSumToSSheet(
      tib=tmp_sum_tib, metric='pass_perc',
      by="Probe_ID", type="Probe_Type", des="Probe_Design", 
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
    purrr::set_names(paste('pOOBAH',names(.), sep='_'))
  
  sigsSumToSSheet(
    tib=tmp_sum_tib, metric='Requeue',
    by="Probe_ID", type="Probe_Type", des="Probe_Design", 
    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
    purrr::set_names(paste('pOOBAH',names(.), sep='_'))
  
  man <- rdat$sman
  call_dat_tib <- rdat$raw_list$call_dat
  
  pval_metric <- 'pass_perc'
  pvals_cols  <- call_dat_tib %>% 
    dplyr::select(!dplyr::ends_with('_beta')) %>% names()
  
  call_pval_ssheet_tib <- NULL
  if (TRUE) {
    
    pval_min_tib <- tibble::tibble(Method=pval_vec, min_pval=min_pval_vec, min_perc=min_perc_vec)
    
    for (pval_idx in c(2:length(pvals_cols)) ) {
      pvals_col <- pvals_cols[pval_idx]
      pvals_key <- pvals_col %>% stringr::str_remove('^[^_]+_')
      
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pval_idx={pval_idx}; pvals_col={pvals_col}...{RET}"))
      
      call_pval_sel_tib <- call_dat_tib %>% 
        dplyr::select(!!by,!!pvals_col)
      
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} call_pval_sel_tib={RET}"))
      print(call_pval_sel_tib)
      
      call_pval_sum_tib <- call_pval_sel_tib %>%
        sigsTibToSummary(man=man, cutoff=minNegPval, 
                         percision=percision_pval,
                         by=by, type=type, des="Probe_Design", 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} call_pval_sum_tib={RET}"))
      print(call_pval_sum_tib)
      
      cur_pval_ssheet <- 
        sigsSumToSSheet(
          tib=call_pval_sum_tib, metric=pval_metric,
          by="Probe_ID", type="Probe_Type", des="Probe_Design", 
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        purrr::set_names(paste(pvals_key,names(.), sep='_'))
      
      call_pval_ssheet_tib <- call_pval_ssheet_tib %>% 
        dplyr::bind_cols(cur_pval_ssheet)
    }
    
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} call_pval_ssheet_tib={RET}"))
    print(call_pval_ssheet_tib)
  }
  
  raw_inf2_tib <- raw_sset@II %>% 
    tibble::as_tibble(rownames = "Probe_ID") %>% 
    dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2)) %>%
    dplyr::filter(Probe_Type=='cg') %>%
    tidyr::separate(Probe_ID, into=c("Seq_ID", "Seq_Pos", "Seq_Str"), 
                    sep='_', remove=FALSE) %>% 
    tidyr::separate(Seq_Str, into=c("Blank","Strand_TB", "Strand_CO", "Design_Type", "Rep_Num"), 
                    sep="", convert=TRUE) %>% dplyr::select(-Blank) %>%
    dplyr::mutate(Probe_Key=paste(Seq_ID,Seq_Pos,Strand_TB,Strand_CO, sep='_'),
                  Probe_Key=as.factor(Probe_Key),
                  Beta=M/(M+U)) %>%
    # dplyr::group_by(Probe_Key) %>% 
    dplyr::arrange(Beta)
  
  
  
  
  
  
  
  
  
  
  rdat$dat_list$sam_sheet$pOOBAH_cg_1_pass_perc_0
  
  rdat$dat_list$sam_sheet %>% addRequeueFlags(minOobPerc = opt$minOobPerc, minNegPerc = opt$minNegPerc)
  
  
  # doc_dir <- '/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/dat/docs'
  doc_dir <- '/Users/bretbarnes/Documents/tools/docs'
  
  # ssh_csv <- file.path(doc_dir, 'AutoSampleSheet_description-table.csv')
  ssh_csv <- '/Users/bretbarnes/Documents/scratch/docker-v.3.46/swifthoof/Customer-Facing/swifthoof/swifthoof_main/202761400007_R01C01_EPIC_B4_AutoSampleSheet.csv.gz'
  
  # ssh_tib <- rdat$ssheet_tib
  ssh_tib <- readr::read_csv(ssh_csv)
  
  sheet_desc_csv <- file.path(doc_dir, 'AutoSampleSheetDescriptionTable.docker.v.3.46.csv.gz')
  sheet_desc_tib <- getSsheetDataTab(tib=ssh_tib, 
                                     minOobPval=opt$minOobPval, minOobPerc=opt$minOobPerc,
                                     minNegPval=opt$minNegPval, minNegPerc=opt$minNegPerc,
                                     verbose = 1) %>%
    dplyr::rename(Example_Value=Value)
  
  sheet_desc_tib %>% print(n=base::nrow(sheet_desc_tib))
  
  readr::write_csv(sheet_desc_tib, sheet_desc_csv)

  
  new_map <- idatToManifestMap(tib = rdat$isig, mans = mans, verbose = 10)  
}



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Old Direct code for calling RAW:: Removed
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  idx <- 0
  cur_workflow <- 'raw'
  cur_sset <- NULL
  
  cur_sset_rds <- paste(out_prefix, paste(cur_workflow,'sset.rds', sep='-'), sep='_')
  
  cur_sset <- mutateSSET_workflow(
    sset=raw_sset, workflow=cur_workflow, pvals=pvals,
    save=opts$save_sset, rds=cur_sset_rds,
    verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
  
  calls_tib <- NULL
  cur_dat_list = ssetToSummary(
    sset=cur_sset, man=top_man_tib, idx=idx, workflow=cur_workflow,
    name=out_name, outDir=opts$outDir,
    
    pvals=pvals, min_pvals=min_pvals, min_percs=min_percs,
    
    write_sset=opts$write_sset, sset_rds=NULL, ret_sset=retData2,
    
    write_beta=opts$write_beta, beta_csv=NULL, ret_beta=retData2,
    write_bsum=opts$write_bsum, bsum_csv=NULL, ret_bsum=retData2,
    
    write_pval=opts$write_pval, pval_csv=NULL, ret_pval=retData2,
    write_psum=opts$write_psum, psum_csv=NULL, ret_psum=retData2,
    
    write_sigs=opts$write_sigs, sigs_csv=NULL, ret_sigs=retData2,
    write_ssum=opts$write_ssum, ssum_csv=NULL, ret_ssum=retData2,
    
    write_call=opts$write_call, call_csv=NULL, ret_call=retData2,
    write_csum=opts$write_csum, csum_csv=NULL, ret_csum=retData2,
    
    percision_sigs=opts$percision_sigs,
    percision_beta=opts$percision_beta,
    percision_pval=opts$percision_pval,
    
    by="Probe_ID", type="Probe_Type", des="Probe_Design",
    fresh=opts$fresh,
    verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
  
  calls_tib = cur_dat_list$call_dat
  sheet_tib = cur_dat_list$sam_sheet
  ssheet_tib <- dplyr::bind_cols(ssheet_tib, sheet_tib)
  
  if (is.null(auto_beta_key)) auto_beta_key <- paste(cur_workflow,'beta', sep=del)
  if (is.null(auto_negs_key)) auto_negs_key <- paste(cur_workflow,'PnegEcdf', sep=del)
  
  if (verbose>=vt) {
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. cur_workflow={cur_workflow}.{RET}{RET}"))
    cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  }
  
  if (retData) {
    ret$raw_list <- cur_dat_list
    # return(ret)
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Add Requeue Flags::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  #
  # TBD:: Make function:: addRequeueFlags()
  #
  if (verbose>=vt) {
    cat(glue::glue("{RET}{RET}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Requeue Flag",
                   "(minOobPerc={opts$minOobPerc},minNegPerc={opts$minNegPerc})...{RET}"))
  }
  
  ssheet_tib <- ssheet_tib %>% dplyr::mutate(
    Requeue_Flag_pOOBAH=dplyr::case_when(
      pOOBAH_cg_1_pass_perc_0 < opts$minOobPerc | pOOBAH_cg_2_pass_perc_0 < opts$minOobPerc ~ TRUE, TRUE ~ FALSE),
    Requeue_Flag_PnegEcdf=dplyr::case_when(
      PnegEcdf_cg_1_pass_perc_0 < opts$minNegPerc | PnegEcdf_cg_2_pass_perc_0 < opts$minNegPerc ~ TRUE, TRUE ~ FALSE)
  )
  
  if (verbose>=vt) {
    cat(glue::glue("[{funcTag}]:{tabsStr} Done. ssheet_tib(all requeue added)={RET}"))
    if (verbose>=vt+4) ssheet_tib %>% print()
    cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}{RET}"))
  }
  
  if (retData) {
    ret$ssheet_tib <- ssheet_tib
    # return(ret)
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Old Functions Removed::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  ssetToPredictions = function(sset, del='_', fresh=FALSE,
                               quality.mask = FALSE, sum.TypeI = FALSE,
                               verbose=0,vt=3,tc=1,tt=NULL) {
    funcTag <- 'ssetToPredictions'
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- system.time({
      
      platform <- NULL
      platform <- sset@platform
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Using platform={platform} for Inference/Prediction calls.{RET}{RET}"))
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                   Sesame Channel Swap Methods:: Extra(sset)
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      igr_cnt <- 0
      irg_cnt <- 0
      igr_key <- 'SwapCntIGR'
      irg_key <- 'SwapCntIRG'
      if (!is.null(sesame::extra(sset)[['IRR']]) && 
          !is.null(sesame::extra(sset)[['IGG']])) {
        igr_cnt <- length(which(!sesame::extra(sset)[['IGG']]))
        irg_cnt <- length(which(!sesame::extra(sset)[['IRR']]))
      }
      ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!igr_key := !!igr_cnt))
      ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!irg_key := !!irg_cnt))
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                       Sesame Inference Methods:: SSET
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if (platform=='EPIC' || platform=='HM450') {
        gct_val <- NULL
        gct_key <- 'GCT'
        gct_val <- safeGCT(sset=sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!gct_key := !!gct_val))
        
        sex_val <- NULL
        sex_key <- 'Sex'
        sex_val <- safeSex(sset=sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!sex_key := !!sex_val))
        
        kar_val <- NULL
        kar_key <- 'SexKaryotype'
        kar_val <- safeSexKaryo(sset=sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!kar_key := !!kar_val))
        
        eth_val <- NULL
        eth_key <- 'Ethnicity'
        eth_val <- safeEthnicity(sset=sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!eth_key := !!eth_val))
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                      Sesame Predictions Methods:: Beta
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if (platform=='EPIC' || platform=='HM450') {
        if (fresh || is.null(sesame::extra(sset)[['betas']]) ) {
          if (verbose>=vt+1)
            cat(glue::glue("[{funcTag}]:{tabsStr} Generating beta calls...{RET}"))
          
          # sesame::extra(sset)[['betas']] <-
          #   sesame::getBetas(sset=sset, mask = quality.mask,sum.TypeI = sum.TypeI)
          sesame::extra(sset)[['betas']] <-
            getBetas2(sset=sset, mask = quality.mask,sum.TypeI = sum.TypeI)
        }
        
        if (platform=='EPIC') {
          skn_val <- NULL
          skn_key <- 'AgeSkinBlood'
          skn_val <- safeSkinAge(beta=sesame::extra(sset)[['betas']], 
                                 verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
          ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!skn_key := !!skn_val))
        }
        
        if (platform=='HM450') {
          phn_val <- NULL
          phn_key <- 'AgePheno'
          phn_val <- safePhenoAge(beta=sesame::extra(sset)[['betas']], 
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
          ret_tib <- ret_tib %>% dplyr::bind_cols(tibble::tibble(!!phn_key := !!phn_val))
        }
      }
      
      ret_cnt <- ret_tib %>% base::nrow()
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
        print(ret_tib)
      }
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                     "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    
    ret_tib
  }
  
  ssetToPvalSheet = function(sset, key, min, perc, man=NULL,
                             percision = -1,
                             by="Probe_ID", type="Probe_Type", des="Probe_Class",
                             verbose=0,vt=3,tc=1,tt=NULL) {
    funcTag <- 'ssetToPvalSheet'
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- system.time({
      
      cur_tib <- sesame::extra(sset)[['pvals']][[key]] %>% 
        tibble::enframe(name='Probe_ID', value=key)
      
      sum_tib <- sigsTibToSummary(tib=cur_tib, man=man, pval=min, perc=perc,
                                  percision=percision,
                                  by=by, type=type, des="Probe_Design", 
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      
      
      ret_tib <- sum_tib
      ret_cnt <- ret_tib %>% base::nrow()
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                     "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    
    ret_tib
  }

  ssetToSigsTib = function(sset, man=NULL, 
                           percision=-1, sort=FALSE, 
                           save=FALSE, csv=NULL,
                           by="Probe_ID", type="Probe_Type", des="Probe_Design", 
                           verbose=0,vt=3,tc=1,tt=NULL) {
    funcTag <- 'ssetToSigsTib'
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- system.time({
      
      by_sym  <- by %>% rlang::sym()
      des_sym <- des %>% rlang::sym()
      
      ret_tib <- dplyr::bind_rows(
        sset@IG   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='IG'),
        sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='OG'),
        
        sset@IR   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='IR'),
        sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='OR'),
        
        sset@II   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='2'),
        
        if (sset@ctl %>% base::nrow() != 0) {
          sset@ctl  %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des_sym:='2') %>%
            dplyr::rename(M=G,U=R) %>% dplyr::select(Probe_ID,!!des_sym,M,U)
        }
      ) %>% dplyr::select(!!by, !!des, dplyr::everything())
      
      if (percision!=-1) ret_tib <- ret_tib %>% 
        dplyr::mutate_if(is.numeric, list(round), percision)
      
      if (!is.null(man)) ret_tib <- dplyr::select(man, !!by, !!type) %>% 
        dplyr::right_join(ret_tib, by=by)
      
      if (sort) ret_tib <- ret_tib %>% dplyr::arrange(!!by_sym)
      
      if (!is.null(save) && save==TRUE && !is.null(csv)) {
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing (percision={percision}) CSV={csv}.{RET}"))
        csv_dir <- base::basename(csv)
        if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
        readr::write_csv(ret_tib, csv)
      }
      
      ret_cnt <- ret_tib %>% base::nrow()
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
        print(ret_tib)
      }
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                     "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    
    ret_tib
  }
  
  ssetToCallTib = function(sset, workflow, fresh=FALSE,
                           quality.mask = FALSE, sum.TypeI = FALSE,
                           percision_beta=-1, percision_pval=-1,
                           save=FALSE, csv=NULL,
                           verbose=0,vt=3,tc=1,tt=NULL) {
    funcTag <- 'ssetToCallTib'
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}.{RET}"))
    if (verbose>=vt+4) print(sset)
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- system.time({
      
      # beta
      #
      name <- paste(workflow,'beta', sep='_')
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Mutating/Settting name={name}...{RET}"))
      
      #
      #
      # TBD:: Ensure we don't have Probe_ID of zero length...
      #
      #
      
      #
      # ssetToBeta (provide return type=tib/dat)
      #
      if (fresh || is.null(sesame::extra(sset)[['betas']]) ) {
        # sesame::extra(sset)[['betas']] <-
        #   sesame::getBetas(sset=sset, mask = quality.mask,sum.TypeI = sum.TypeI)
        sesame::extra(sset)[['betas']] <-
          getBetas2(sset=sset, mask = quality.mask,sum.TypeI = sum.TypeI)
      }
      ret_tib <- tibble::enframe(sesame::extra(sset)[['betas']], 
                                 name='Probe_ID', value=name)
      if (percision_beta!=-1)
        ret_tib <- dplyr::mutate_if(ret_tib, purrr::is_double,round,percision_beta)
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} Current Table(betas)={RET}"))
        ret_tib %>% print()
      }
      
      #
      # Detection p-values::
      #
      pval_names <- sesame::extra(sset)[['pvals']] %>% names()
      for (pval_name in pval_names) {
        out_name <- pval_name
        out_name <- paste(workflow,out_name, sep='_')
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pval_name={pval_name}, out_name={out_name}...{RET}"))
        
        cur_tib <- sesame::extra(sset)[['pvals']][[pval_name]] %>% 
          tibble::enframe(name='Probe_ID', value=out_name)
        
        if (percision_pval!=-1) 
          cur_tib <- dplyr::mutate_if(cur_tib, purrr::is_double,round,percision_pval)
        
        ret_tib <- ret_tib %>% dplyr::left_join(cur_tib,by="Probe_ID")
      }
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} Current Table(pvals)={RET}"))
        ret_tib %>% print()
      }
      
      if (!is.null(save) && save==TRUE && !is.null(csv)) {
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing CSV={csv}.{RET}"))
        csv_dir <- base::basename(csv)
        if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
        readr::write_csv(ret_tib, csv)
      }
      
      ret_cnt <- ret_tib %>% base::nrow()
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
        print(ret_tib)
      }
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                     "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    
    ret_tib
  }
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sex Karyotypes Investigation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  man_tib <- rdat$sman
  raw_sset <- newSset(rdat$prefix, platform = rdat$ssheet_tib$detect_platform[1], 
                      manifest = man_tib, verbose = 20)
  
  cur_sset_rds <- NULL
  cur_sigs_csv <- NULL
  cur_ssum_csv <- NULL
  cur_call_csv <- NULL
  
  idx <- 2
  cur_workflow <- 'ind'
  retData <- FALSE
  out_name <- 'test'
  
  ind_sset <- mutateSSET_workflow(
    sset=raw_sset, workflow=cur_workflow, 
    save=opt$save_sset, rds=cur_sset_rds,
    verbose=10)
  
  ind_lst = ssetToSummary(
    sset=ind_sset, man=man_tib, idx=idx, workflow=cur_workflow,
    name=out_name, outDir=opt$outDir,
    
    write_sset=opt$write_sset, sset_rds=cur_sset_rds, ret_sset=retData,
    write_sigs=opt$write_sigs, sigs_csv=cur_sigs_csv, ret_sigs=retData,
    write_ssum=opt$write_ssum, ssum_csv=cur_ssum_csv, ret_ssum=retData,
    write_call=opt$write_call, call_csv=cur_call_csv, ret_call=retData,
    
    minNegPval=opt$minNegPval,minOobPval=opt$minOobPval,
    percision_sigs=opt$percision_sigs,
    percision_beta=opt$percision_beta,
    percision_pval=opt$percision_pval,
    
    by="Probe_ID", type="Probe_Type", des="Probe_Class",
    fresh=opt$fresh,
    verbose=10
  )
  
  
  safeSexKaryo(sset = raw_sset, verbose = 20)
  
  inferSexKaryotypes2(sset=raw_sset)
  
  sesame::inferSexKaryotypes(sset=raw_sset)
  
  sesame::getSexInfo(sset=raw_sset)
  
}


if (FALSE) {
  ses_dat <- idat2manifest(sigs=idat_list$sig, mans=mans,
                           verbose=verbose,tc=tc+1,tt=tTracker)
  
  # Make sure we have unique names::
  ses_dat$man <- ses_dat$man %>% dplyr::distinct(Probe_ID, .keep_all=TRUE)
  ses_man_tib <- ses_dat$man
  
  ses_add_tib <- ses_dat$add
  platform_key <- ses_dat$platform
  version_key  <- ses_dat$manifest
  
  bead_sum_tib <- ses_add_tib %>% 
    dplyr::filter(Probe_Type=='cg') %>% 
    dplyr::select(Address) %>% dplyr::left_join(idat_list$sig, by="Address") %>% 
    dplyr::summarise(CG_Bead_Count=n(),
                     CG_Bead_Total=sum(Bead_Grn,Bead_Red, na.rm=TRUE), 
                     CG_Bead_AvgRep=round(CG_Bead_Total/CG_Bead_Count/2),
                     .groups='drop')
  if (verbose>=vt+5) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} bead_sum_tib={RET}"))
  if (verbose>=vt+5) bead_sum_tib %>% print()
  
  pool_sum_tib <- ses_man_tib %>% 
    dplyr::filter(Probe_Type=='cg' | Probe_Type=='ch' | Probe_Type=='rs') %>% 
    dplyr::mutate(Probe_Type=stringr::str_to_upper(Probe_Type)) %>%
    dplyr::group_by(Probe_Type) %>%
    dplyr::summarise(Count=n(), .groups='drop') %>% 
    tidyr::spread(Probe_Type,Count) %>% 
    purrr::set_names(paste(names(.),'Manifest_Count',sep='_') ) %>% 
    addBeadPoolToSampleSheet(field='CG_Manifest_Count') %>% dplyr::ungroup()
  if (verbose>=vt+5) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pool_sum_tib={RET}"))
  if (verbose>=vt+5) pool_sum_tib %>% print()
  
  if (retData) {
    ret$sman  <- ses_man_tib
    ret$sadd  <- ses_add_tib
    ret$sbead <- bead_sum_tib
    ret$spool <- pool_sum_tib
  }
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Initializing Auto-Sample-Sheet: [idat/bead/pheno].{RET}"))
  
  chipFormat <- idat_list$ann %>% dplyr::select(Chip_Format) %>% dplyr::pull()
  beadPool   <- pool_sum_tib %>% dplyr::select(Bead_Pool) %>% dplyr::pull()
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} chipFormat={chipFormat}, beadPool={beadPool}.{RET}"))
  
  ssheet_tib <- NULL
  ssheet_tib <- idat_list$ann
  if (!is.null(pool_sum_tib)) ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(pool_sum_tib)
  if (!is.null(bead_sum_tib)) ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(bead_sum_tib)
  ssheet_tib <- ssheet_tib %>% tibble::add_column(minNegPval   = opts$minNegPval)
  ssheet_tib <- ssheet_tib %>% tibble::add_column(minOobPval   = opts$minOobPval)
  ssheet_tib <- ssheet_tib %>% tibble::add_column(minDeltaBeta = opts$minDeltaBeta)
  ssheet_tib <- ssheet_tib %>% tibble::add_column(platformUsed = platform_key)
  ssheet_tib <- ssheet_tib %>% tibble::add_column(platVersUsed = version_key)
}

if (FALSE) {
  
  tar_types <- c('cg','ch','rs')
  tmp_tab <- getManifestBeadStats(dat=tdat$idat$sig, man=tdat$sman, types=tar_types, verbose=10)
  
  
  pool_sum_tib <- tdat$sman %>% 
    dplyr::filter(Probe_Type=='cg' | Probe_Type=='ch' | Probe_Type=='rs') %>% 
    dplyr::mutate(Probe_Type=stringr::str_to_upper(Probe_Type)) %>%
    dplyr::group_by(Probe_Type) %>%
    dplyr::summarise(Count=n(), .groups='drop') %>% 
    tidyr::spread(Probe_Type,Count) %>% 
    purrr::set_names(paste(names(.),'Manifest_Count',sep='_') ) %>% 
    addBeadPoolToSampleSheet(field='CG_Manifest_Count') %>% dplyr::ungroup()
  
  
  tar_types <- c('cg','ch','rs')
  idatAddToSsheet(dat=tdat$idat$sig, add=tdat$tadd, types=tar_types, verbose=10)
  
  
  tar_types <- c('cg','ch','rs')
  getManifestBeadStats(dat=tdat$idat$sig, man=tdat$tman, types=tar_types, verbose=10)
  
  
  tmp_tab <- tdat$tadd %>% 
    dplyr::select(Address,Probe_Type) %>% 
    dplyr::left_join(tdat$idat$sig, by="Address") %>% 
    dplyr::group_by(Probe_Type) %>%
    dplyr::summarise(Bead_Count=n(), 
                     Bead_Total=sum(Bead_Grn,Bead_Red, na.rm=TRUE), 
                     Bead_AvgRep=round(Bead_Total/Bead_Count/2),
                     .groups='drop') %>% 
    tidyr::gather(name, value, -Probe_Type) %>% 
    tidyr::unite(name, Probe_Type,name, sep='_') %>% 
    tidyr::spread(name, value)
  
  
  #
  # Difference between Address Calculations::
  #
  tdat$tadd %>% 
    dplyr::filter(Probe_Type=='cg') %>% 
    dplyr::select(Probe_ID,Design_Type,Probe_Type,Address) %>% 
    dplyr::left_join(tdat$idat$sig, by="Address") %>% 
    dplyr::filter(is.na(Raw_Grn_sig)) %>%
    as.data.frame()
  
  tdat$sadd %>% 
    dplyr::filter(Probe_Type=='cg') %>% 
    dplyr::select(Probe_ID,Design_Type,Probe_Type,Address) %>% 
    dplyr::left_join(tdat$idat$sig, by="Address") %>% 
    dplyr::filter(is.na(Raw_Grn_sig)) %>%
    as.data.frame()
  
  tdat$tadd %>% dplyr::mutate(Add_Len=stringr::str_length(Address)) %>% 
    dplyr::group_by(Add_Len) %>% dplyr::summarise(Len_Count=n(), .groups='drop')
  tdat$sadd %>% dplyr::mutate(Add_Len=stringr::str_length(Address)) %>% 
    dplyr::group_by(Add_Len) %>% dplyr::summarise(Len_Count=n(), .groups='drop')
  
  
  tdat$sadd %>% dplyr::mutate(Add_Len=stringr::str_length(Address)) %>% 
    dplyr::filter(Add_Len==7) %>% dplyr::left_join(tdat$sman, by="Probe_ID")
  

  
  tdat$tadd %>% 
    dplyr::filter(Probe_Type=='cg') %>% 
    dplyr::select(Address) %>% 
    dplyr::left_join(tdat$idat$sig, by="Address") %>% 
    dplyr::summarise(CG_Bead_Count=n(), 
                     CG_Bead_Total=sum(Bead_Grn,Bead_Red, na.rm=TRUE), 
                     CG_Bead_AvgRep=round(CG_Bead_Total/CG_Bead_Count/2),
                     .groups='drop')
  
    
  tdat$sadd %>% 
    dplyr::filter(Probe_Type=='cg') %>% 
    dplyr::select(Address) %>% dplyr::left_join(tdat$idat$sig, by="Address") %>% 
    dplyr::summarise(CG_Bead_Count=n(), 
                     CG_Bead_Total=sum(Bead_Grn,Bead_Red, na.rm=TRUE), 
                     CG_Bead_AvgRep=round(CG_Bead_Total/CG_Bead_Count/2),
                     .groups='drop')
  
  
  
  
  tadd2 <- sesameManToAdd(rdat$tman, verbose = 20)
  
  tadd <- rdat$tman %>% 
    dplyr::select(M,U, Probe_ID, Probe_Type, col) %>% 
    dplyr::mutate(
      Design_Type=dplyr::case_when(
        is.na(M) ~ 'II',
        !is.na(M) & !is.na(col) ~ paste0('I',col),
        TRUE ~ NA_character_)
    ) %>%
    tidyr::gather(MU, Address, -Probe_ID, -col, -Design_Type, -Probe_Type) %>% 
    dplyr::filter(!is.na(Address)) %>%
    dplyr::select(Probe_ID,Address,col,Design_Type,Probe_Type) %>%
    dplyr::arrange(Probe_ID)
  
  
  rdat$tman %>% 
    dplyr::select(M,U, Probe_ID, Probe_Type, col) %>% 
    dplyr::rename(Design_Type=DESIGN, Man_Col=col)
  
  
  
  tadd <- rdat$tman %>% 
    dplyr::select(M,U, Probe_ID, DESIGN, Probe_Type, col) %>% 
    dplyr::rename(Design_Type=DESIGN, Man_Col=col) %>% 
    tidyr::gather(MU, Address, -Probe_ID, -Man_Col, -Design_Type, -Probe_Type) %>% 
    dplyr::filter(!is.na(Address)) %>% 
    dplyr::select(Probe_ID,Address,Man_Col,Design_Type,Probe_Type,MU) %>% 
    dplyr::arrange(Probe_ID)
  
  
  rdat$sadd %>% dplyr::filter(!is.na(Man_Col))
  tadd %>% dplyr::filter(!is.na(Man_Col))
  
  
  
  rdat$tman %>% 
    dplyr::select(M,U, Probe_ID, DESIGN, Probe_Type, col) %>% 
    dplyr::rename(Design_Type=DESIGN, Man_Col=col) %>% 
    tidyr::gather(MU, Address, -Probe_ID, -Man_Col, -Design_Type, -Probe_Type) %>% 
    dplyr::filter(!is.na(Address))
  
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Build Documentation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  doc_dir <- '/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/dat/docs'
  doc_dir <- '/Users/bretbarnes/Documents/tools/docs'
  
  
  sheet_desc_csv <- file.path(doc_dir, 'AutoSampleSheet_description-table.csv')
  
  stmp_tab <- getSsheetDataTab(tib = rdat$ssheet_tib, 
                               minOobPval=opt$minOobPval, minOobPerc=opt$minOobPerc,
                               minNegPval=opt$minNegPval, minNegPerc=opt$minNegPerc,
                               verbose = 1)
  
  
  # Determine Number of Workflows::
  #
  max_idx <- rdat$ssheet_tib %>% 
    dplyr::select(dplyr::starts_with('Method_Idx')) %>% 
    tidyr::gather() %>% dplyr::arrange(-value) %>% 
    head(n=1) %>% dplyr::pull(value)
  
  # Build Description Sheet
  #
  ssheet_desc_tib <- getSsheetCoreAnnoTib()
  for (idx in c(0:max_idx)) {
    ssheet_desc_tib <- ssheet_desc_tib %>% dplyr::bind_cols(getSsheetIndexAnnoTib(idx=idx))
  }
  
  getSsheetDescTib(rdat$ssheet_tib)
  

  # Build Current Sample Sheet Table::
  #
  sheet_data_tab <- dplyr::inner_join(
    rdat$ssheet_tib %>% tidyr::gather(Variable, Example),
    tibble::enframe( unlist( sapply(rdat$ssheet_tib, class) ), name="Variable", value="Data_Type" ),
    by="Variable"
  )
  sheet_data_tab %>% print(n=base::nrow(sheet_data_tab))
  # readr::write_csv(sheet_data_tab,sheet_desc_csv)
  
  sheet_anno_tab
  
  # Join
  sheet_anno_tab <- dplyr::bind_cols(
    getSsheetCoreAnnoTib(),
    getSsheetIndexAnnoTib(idx=0),
    getSsheetIndexAnnoTib(idx=1),
  ) %>% gather(Variable, Description)
  sheet_anno_tab %>% print(n=base::nrow(sheet_anno_tab))
  
  sheet_desc_tab <- sheet_data_tab %>% dplyr::left_join(sheet_anno_tab, by="Variable")
  sheet_desc_tab %>% print(n=base::nrow(sheet_desc_tab))
  
  
  # CG_Bead_Total, CG_Bead_AvgRep
  
  
  
  sheet_anno_tib <- tibble::tibble(
    # Sample Requeue Suggestion::
    #
    Requeue_Flag_pOOBAH = 
      glue::glue("Flag to automatically requeue a faild sample based on percent loci with pOOBAH detection p-value < {opt$minOobPval}. Percent threshold = {opt$minOobPerc}."),
    Requeue_Flag_PnegEcdf = 
      glue::glue("Flag to automatically requeue a faild sample based on percent loci with PnegEcdf detection p-value < {opt$minNegPval}. Percent threshold = {opt$minNegPerc}."),
    
    # Sample Identification::
    #
    Sentrix_Name = "Unique Sentrix Name: Sentrix_Barcode + Sentrix_Poscode.",
    Sentrix_Barcode = "Sentrix Bar Code (AKA Sentrix_ID).",
    Sentrix_Poscode = "Sentrix Position Code (AKA Sentrix_Position).",
    Sentrix_Row = "Sentrix Row on chip.",
    Sentrix_Col = "Sentrix Column on Chip.",
    Chip_Type = "Idat Chip Type.",
    Chip_Format = "Idat Chip Format (e.g. 8x1, 12x1, 24x1, etc.).",
    Bead_Pool = "Automatically identified bead pool based on address overlap (e.g. HM450, EPIC, etc.).",
    
    # Bead and Manifest Loci Summary::
    #
    CG_Manifest_Count = "Total number of sample cg# loci overlapping manifest based on address presence NOT detection p-value.",
    CH_Manifest_Count = "Total number of sample ch# loci overlapping manifest based on address presence NOT detection p-value.",
    RS_Manifest_Count = "Total number of sample rs# loci overlapping manifest based on address presence NOT detection p-value.",
    CG_Bead_Count     = "Total number of sample cg# bead-types found overlapping manifest based on address presence NOT number of bead types.",
    CH_Bead_Count     = "Total number of sample ch# bead-types found overlapping manifest based on address presence NOT number of bead types.",
    RS_Bead_Count     = "Total number of sample rs# bead-types found overlapping manifest based on address presence NOT number of bead types.",
    
    # Parameters Used::
    #
    minNegPval = "Minimum Negative (PnegEcdf) detection p-value threshold used.",
    minOobPval = "Minimum Out-Of-Band (pOOBAH) detection p-value threshold used.",
    minDeltaBeta = "Minimum delta-Beta cutoff used for Sample-Auto-Detection calculations.",
    platformUsed = "Platform name used. Identified during Auto-Manifest-Detection.",
    platVersUsed = "Platform version used. Identified during Auto-Manifest-Detection.",
    
    # Inference and Predictions::
    #
    Method_Key_XXX = glue::glue("Method workflow name used for metrics ending with XXX."),
    Method_Idx_XXX = glue::glue("Method workflow index used for metrics ending with XXX."),
    GCT_XXX = glue::glue("GCT Score for method index XXX. ",
                         "See Bioconductor Package sesame::GCT(sset): ",
                         "Compute GCT score for internal bisulfite conversion control. ",
                         "The function takes a SigSet as input. The higher the GCT score, the more likely the incomplete conversion."),
    
    Sex_XXX = glue::glue("Sex call for method index XXX. ",
                         "See Bioconductor Package ‘sesame’ sesame::inferSex(sset)."),
    SexKaryotype_XXX = glue::glue("Sex Karyotype call for method index XXX. ",
                                  "See Bioconductor Package ‘sesame’ sesame::inferSexKaryotypes(sset)."),
    
    Ethnicity_XXX = glue::glue("Ethnicity call for method index XXX. ",
                               "See Bioconductor Package ‘sesame’ sesame::inferEthnicity(sset).",
                               "This function uses both the built-in rsprobes as well as the type I Color-Channel-Switching probes to infer ethnicity."),
    
    SwapCntIGR_XXX = glue::glue("Number of Grn to Red swapped channels for method index XXX. ",
                                "See Bioconductor Package ‘sesame’ sesame::inferTypeIChannel(sset).",
                                "Infer and reset color channel for Type-I probes instead of using what is specified in manifest."),
    SwapCntIRG_XXX = glue::glue("Number of Red to Grn swapped channels for method index XXX. ",
                                "See Bioconductor Package ‘sesame’ sesame::inferTypeIChannel(sset).",
                                "Infer and reset color channel for Type-I probes instead of using what is specified in manifest."),
    
    AgeSkinBlood_XXX = glue::glue("Horvath Skin and Blood age predictor call for method index XXX. ",
                                  "See Bioconductor Package ‘sesame’ sesame::predictAgeSkinBlood(betas)."),
    
    # Detection P-values:: pOOBAH
    #
    pOOBAH_cg_1_pass_perc_XXX = glue::glue("Percent cg Infinium I loci with pOOBAH detection p-value < {opt$minOobPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
    pOOBAH_cg_2_pass_perc_XXX = glue::glue("Percent cg Infinium II loci with pOOBAH detection p-value < {opt$minOobPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),

    pOOBAH_ch_1_pass_perc_XXX = glue::glue("Percent ch Infinium I loci with pOOBAH detection p-value < {opt$minOobPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
    pOOBAH_ch_2_pass_perc_XXX = glue::glue("Percent ch Infinium II loci with pOOBAH detection p-value < {opt$minOobPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),

    pOOBAH_rs_1_pass_perc_XXX = glue::glue("Percent rs Infinium I loci with pOOBAH detection p-value < {opt$minOobPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
    pOOBAH_rs_2_pass_perc_XXX = glue::glue("Percent rs Infinium II loci with pOOBAH detection p-value < {opt$minOobPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
    
    # Detection P-values:: PnegEcdf
    #
    PnegEcdf_cg_1_pass_perc_XXX = glue::glue("Percent cg Infinium I loci with PnegEcdf detection p-value < {opt$minNegPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
    PnegEcdf_cg_2_pass_perc_XXX = glue::glue("Percent cg Infinium II loci with PnegEcdf detection p-value < {opt$minNegPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
    
    PnegEcdf_ch_1_pass_perc_XXX = glue::glue("Percent ch Infinium I loci with PnegEcdf detection p-value < {opt$minNegPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
    PnegEcdf_ch_2_pass_perc_XXX = glue::glue("Percent ch Infinium II loci with PnegEcdf detection p-value < {opt$minNegPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
    
    PnegEcdf_rs_1_pass_perc_XXX = glue::glue("Percent rs Infinium I loci with PnegEcdf detection p-value < {opt$minNegPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
    PnegEcdf_rs_2_pass_perc_XXX = glue::glue("Percent rs Infinium II loci with PnegEcdf detection p-value < {opt$minNegPval} for method index XXX. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
    

    # Average Beta Values::
    #
    cg_1_mean_XXX = glue::glue("Mean beta-value of cg Infinium I loci for method XXX."),
    cg_2_mean_XXX = glue::glue("Mean beta-value of cg Infinium II loci for method XXX."),
    
    ch_1_mean_XXX = glue::glue("Mean beta-value of ch Infinium I loci for method XXX."),
    ch_2_mean_XXX = glue::glue("Mean beta-value of ch Infinium II loci for method XXX."),
    
    rs_1_mean_XXX = glue::glue("Mean beta-value of rs Infinium I loci for method XXX."),
    rs_2_mean_XXX = glue::glue("Mean beta-value of rs Infinium II loci for method XXX."),
    
    # Average Intensities::
    #
    BISULFITE_CONVERSION_I_2_M_mean_XXX = glue::glue("Mean intensity of BISULFITE_CONVERSION_I Grn loci for method XXX."),
    BISULFITE_CONVERSION_I_2_U_mean_XXX = glue::glue("Mean intensity of BISULFITE_CONVERSION_I Red loci for method XXX."),
    
    BISULFITE_CONVERSION_II_2_M_mean_XXX = glue::glue("Mean intensity of BISULFITE_CONVERSION_II Grn loci for method XXX."),
    BISULFITE_CONVERSION_II_2_U_mean_XXX = glue::glue("Mean intensity of BISULFITE_CONVERSION_II Red loci for method XXX."),
    
    cg_2_M_mean_XXX = glue::glue("Mean intensity of cg Infinium II Grn (methylated) loci for method XXX."),
    cg_2_U_mean_XXX = glue::glue("Mean intensity of cg Infinium II Red (unmethylated) loci for method XXX."),
    
    cg_IG_M_mean_XXX = glue::glue("Mean intensity of cg Infinium I In-Bound Grn (methylated) loci for method XXX."),
    cg_IG_U_mean_XXX = glue::glue("Mean intensity of cg Infinium I In-Bound Grn (unmethylated) loci for method XXX."),
    cg_IR_M_mean_XXX = glue::glue("Mean intensity of cg Infinium I In-Bound Red (methylated) loci for method XXX."),
    cg_IR_U_mean_XXX = glue::glue("Mean intensity of cg Infinium I In-Bound Red (unmethylated) loci for method XXX."),
    
    cg_OG_M_mean_XXX = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Grn (methylated) loci for method XXX."),
    cg_OG_U_mean_XXX = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Grn (unmethylated) loci for method XXX."),
    cg_OR_M_mean_XXX = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Red (methylated) loci for method XXX."),
    cg_OR_U_mean_XXX = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Red (unmethylated) loci for method XXX."),
    
    ch_2_M_mean_XXX = glue::glue("Mean intensity of ch Infinium II Grn (methylated) loci for method XXX."),
    ch_2_U_mean_XXX = glue::glue("Mean intensity of ch Infinium II Red (unmethylated) loci for method XXX."),
    
    EXTENSION_2_M_mean_XXX = glue::glue("Mean intensity of EXTENSION Grn loci for method XXX."),
    EXTENSION_2_U_mean_XXX = glue::glue("Mean intensity of EXTENSION Red loci for method XXX."),
    
    HYBRIDIZATION_2_M_mean_XXX = glue::glue("Mean intensity of HYBRIDIZATION Grn loci for method XXX."),
    HYBRIDIZATION_2_U_mean_XXX = glue::glue("Mean intensity of HYBRIDIZATION Red loci for method XXX."),
    
    NEGATIVE_2_M_mean_XXX = glue::glue("Mean intensity of NEGATIVE Grn loci for method XXX."),
    NEGATIVE_2_U_mean_XXX = glue::glue("Mean intensity of NEGATIVE Red loci for method XXX."),

    NON_POLYMORPHIC_2_M_mean_XXX = glue::glue("Mean intensity of NON_POLYMORPHIC Grn loci for method XXX."),
    NON_POLYMORPHIC_2_U_mean_XXX = glue::glue("Mean intensity of NON_POLYMORPHIC Red loci for method XXX."),
    
    NORM_A_2_M_mean_XXX = glue::glue("Mean intensity of NORM_A Grn loci for method XXX."),
    NORM_A_2_U_mean_XXX = glue::glue("Mean intensity of NORM_A Red loci for method XXX."),
    
    NORM_C_2_M_mean_XXX = glue::glue("Mean intensity of NORM_C Grn loci for method XXX."),
    NORM_C_2_U_mean_XXX = glue::glue("Mean intensity of NORM_C Red loci for method XXX."),
    
    NORM_G_2_M_mean_XXX = glue::glue("Mean intensity of NORM_G Grn loci for method XXX."),
    NORM_G_2_U_mean_XXX = glue::glue("Mean intensity of NORM_G Red loci for method XXX."),
    
    NORM_T_2_M_mean_XXX = glue::glue("Mean intensity of NORM_T Grn loci for method XXX."),
    NORM_T_2_U_mean_XXX = glue::glue("Mean intensity of NORM_T Red loci for method XXX."),
    
    RESTORATION_2_M_mean_XXX = glue::glue("Mean intensity of RESTORATION Grn loci for method XXX."),
    RESTORATION_2_U_mean_XXX = glue::glue("Mean intensity of RESTORATION Red loci for method XXX."),
    
    rs_2_M_mean_XXX = glue::glue("Mean intensity of rs Infinium II Grn (methylated) loci for method XXX."),
    rs_2_U_mean_XXX = glue::glue("Mean intensity of rs Infinium II Red (unmethylated) loci for method XXX."),
    
    rs_IG_M_mean_XXX = glue::glue("Mean intensity of rs Infinium I In-Bound Grn (methylated) loci for method XXX."),
    rs_IG_U_mean_XXX = glue::glue("Mean intensity of rs Infinium I In-Bound Grn (unmethylated) loci for method XXX."),
    rs_IR_M_mean_XXX = glue::glue("Mean intensity of rs Infinium I In-Bound Red (methylated) loci for method XXX."),
    rs_IR_U_mean_XXX = glue::glue("Mean intensity of rs Infinium I In-Bound Red (unmethylated) loci for method XXX."),
    
    rs_OG_M_mean_XXX = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Grn (methylated) loci for method XXX."),
    rs_OG_U_mean_XXX = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Grn (unmethylated) loci for method XXX."),
    rs_OR_M_mean_XXX = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Red (methylated) loci for method XXX."),
    rs_OR_U_mean_XXX = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Red (unmethylated) loci for method XXX."),
    
    SPECIFICITY_I_2_M_mean_XXX = glue::glue("Mean intensity of SPECIFICITY_I Grn loci for method XXX."),
    SPECIFICITY_I_2_U_mean_XXX = glue::glue("Mean intensity of SPECIFICITY_I Red loci for method XXX."),
    
    SPECIFICITY_II_2_M_mean_XXX = glue::glue("Mean intensity of SPECIFICITY_II Grn loci for method XXX."),
    SPECIFICITY_II_2_U_mean_XXX = glue::glue("Mean intensity of SPECIFICITY_II Red loci for method XXX."),
    
    STAINING_2_M_mean_XXX = glue::glue("Mean intensity of STAINING Grn loci for method XXX."),
    STAINING_2_U_mean_XXX = glue::glue("Mean intensity of STAINING Red loci for method XXX."),

    TARGET_REMOVAL_2_M_mean_XXX = glue::glue("Mean intensity of TARGET_REMOVAL Grn loci for method XXX."),
    TARGET_REMOVAL_2_U_mean_XXX = glue::glue("Mean intensity of TARGET_REMOVAL Red loci for method XXX.")
  )
  
  sheet_anno_tib %>% gather(Variable, Description) %>% print(n=base::ncol(sheet_anno_tib))
  
  sheet_data_tab %>% print(n=57)
  
  
    
  
  sheet_data_tab
  sheet_data_tab %>% print(n=base::nrow(sheet_data_tab))
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Bind all Calls/Sample Sheets::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (verbose>=vt+4) {
  cat(glue::glue("{RET}[{funcTag}]:{tabsStr} ret_dat_list={RET}"))
  ret_dat_list %>% print()
  cat(glue::glue("{RET}[{funcTag}]:{tabsStr} names(ret_dat_list)={RET}"))
  names(ret_dat_list) %>% print()
  cat(glue::glue("{RET}[{funcTag}]:{RET}{RET}"))
}

# all_call_tib <- NULL
# all_call_tib <- ses_man_tib %>% dplyr::select(Probe_ID) %>% 
#   dplyr::arrange(Probe_ID)
for (work_name in names(ret_dat_list)) {
  cur_list  <- ret_dat_list[[work_name]]
  cur_data  <- cur_list$call_dat
  cur_sheet <- cur_list$sam_sheet
  
  if (verbose>=vt+5) {
    cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Cur Join; call_dat({work_name})={RET}"))
    cur_data %>% print()
  }
  if (verbose>=vt+5) {
    cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Cur Join; sam_sheet({work_name})={RET}"))
    cur_sheet %>% print()
  }
  
  ssheet_tib <- ssheet_tib %>% 
    dplyr::bind_cols(cur_sheet)
  
  # if (verbose>=vt+5) {
  #   cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}"))
  #   cat(glue::glue("[{funcTag}]:{tabsStr} Cur Join; all_call_tib({work_name})={RET}"))
  #   all_call_tib %>% print()
  # }
  
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; work_name={work_name}{RET}"))
    cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  }
}
# all_call_tib <- all_call_tib %>% dplyr::arrange(Probe_ID)

if (verbose>=vt) 
  cat(glue::glue("{RET}Done Mergeing...{RET}{RET}{RET}{RET}{RET}"))


if (retData) {
  ret$rsum_list <- ret_dat_list
  # ret$all_calls <- all_call_tib
  ret$ssheet_tib <- ssheet_tib
  # return(ret)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Old Test Code::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  if (FALSE) {
    poob_val  <- FALSE
    poob1_val <- ssheet_tib$pOOBAH_cg_1_pass_perc_0[1]
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} poob1_val={poob1_val}.{RET}"))
    
    poob2_val <- ssheet_tib$pOOBAH_cg_2_pass_perc_0[1]
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} poob2_val={poob2_val}.{RET}"))
    
    if (poob1_val<opt$minOobPerc) poob_val <- TRUE
    if (poob2_val<opt$minOobPerc) poob_val <- TRUE
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} poob_val[{poob1_val},{poob2_val}]<{opt$minOobPerc}={poob_val}.{RET}"))
    ssheet_tib <- ssheet_tib %>% tibble::add_column(Requeue_Flag_pOOBAH=poob_val)
    if (verbose>=vt+4) {
      cat(glue::glue("{RET}[{funcTag}]:{tabsStr} ssheet_tib(requeue-poob added)={RET}"))
      ssheet_tib %>% print()
      cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    }
    
    negs_val  <- FALSE
    negs1_val <- ssheet_tib$PnegEcdf_cg_1_pass_perc_0[1]
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} negs1_val={negs1_val}.{RET}"))
    negs2_val <- ssheet_tib$PnegEcdf_cg_2_pass_perc_0[1]
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} negs2_val={negs2_val}.{RET}"))
    
    if (negs1_val<opt$minNegPerc) negs_val <- TRUE
    if (negs2_val<opt$minNegPerc) negs_val <- TRUE
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} negs_val[{negs1_val},{negs2_val}]<{opt$minNegPerc}={negs_val}.{RET}"))
    ssheet_tib <- ssheet_tib %>% tibble::add_column(Requeue_Flag_PnegEcdf=negs_val)
    
    if (verbose>=vt+4) {
      cat(glue::glue("{RET}[{funcTag}]:{tabsStr} ssheet_tib(requeue-negs added)={RET}"))
      ssheet_tib %>% print()
      cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    }
    
    # ssheet_tib <- ssheet_tib %>% dplyr::mutate(
    #   Requeue_Flag_pOOBAH=dplyr::case_when(
    #     poob1_val < opt$minOobPerc | poob2_val < opt$minOobPerc ~ TRUE, TRUE ~ FALSE),
    #   Requeue_Flag_PnegEcdf=dplyr::case_when(
    #     negs1_val < opt$minNegPerc | negs2_val < opt$minNegPerc ~ TRUE, TRUE ~ FALSE)
    # ) %>% dplyr::select(Requeue_Flag_pOOBAH,Requeue_Flag_PnegEcdf, dplyr::everything())
  }
  
  
}
if (FALSE) {
  
  auto_beta_key <- 'ind_beta'
  auto_negs_key <- 'ind_PnegEcdf'
  
  auto_list <- 
    autoDetect_Wrapper(can=rdat$all_calls, ref=auto_sam_tib, man=rdat$sman,
                       
                       minPval=opt$minNegPval, minDelta=opt$minDeltaBeta,
                       dname='Design_Type', pname='Probe_Type', ptype='cg',
                       jval='Probe_ID', field=auto_beta_key, pval=auto_negs_key, suffix='beta', del='_',
                       outDir=opt$outDir, sname=out_name, plotMatrix=opt$plotAuto, writeMatrix=opt$writeAuto,
                       dpi=opt$dpi, format=opt$plotFormat, datIdx=4, non.ref=non_ref,
                       
                       verbose=20)
  
  
  work_data  <- rdat$rsum_list
  work_names <- work_data %>% names()
  
  all_call_tib <- NULL
  for (work_name in work_names) {
    if (is.null(all_call_tib)) {
      all_call_tib <- work_data[[work_name]]$call_dat
    } else {
      all_call_tib <- all_call_tib %>% 
        dplyr::inner_join(work_data[[work_name]]$call_dat, by="Probe_ID")
    }
  }
  
  
  
  
  
  
  
  metric='mean'
  platform <- rdat$platform
  man_tib <- rdat$sman
  
  ww <- 0
  cur_workflow <- 'raw2'
  out_name <- 'test'
  
  opt$load_sset <- FALSE
  opt$save_sset <- FALSE
  
  sset_rds <- file.path(opt$outDir, 'raw.sset.rds')
  raw_sset <- newSset(prefix=chipPrefixes[[prefix]], 
                      platform=platform, manifest=man_tib,
                      load=opt$load_sset, save=opt$save_sset,rds=sset_rds,
                      verbose=4)
  
  
  
  raw_call_tib <- ssetToCallTib(sset=raw_sset, workflow = 'raw', verbose = 10)
  
  cur_sset_rds <- NULL
  cur_sigs_csv <- NULL
  cur_ssum_csv <- NULL
  cur_call_csv <- NULL
  retData <- FALSE
  
  cdat <- 
    ssetToSummary(
      sset=raw_sset, man=man_tib, idx=ww, workflow=cur_workflow,
      name=out_name, outDir=opt$outDir,
      
      write_sset=opt$write_sset, sset_rds=cur_sset_rds, ret_sset=retData,
      write_sigs=opt$write_sigs, sigs_csv=cur_sigs_csv, ret_sigs=retData,
      write_ssum=opt$write_ssum, ssum_csv=cur_ssum_csv, ret_ssum=retData,
      write_call=opt$write_call, call_csv=cur_call_csv, ret_call=TRUE,
      
      minNegPval=opt$minNegPval,minOobPval=opt$minOobPval,
      percision_sigs=opt$percision_sigs,
      percision_beta=opt$percision_beta,
      percision_pval=opt$percision_pval,
      
      by="Probe_ID", type="Probe_Type", des="Probe_Class",
      fresh=opt$fresh,
      verbose=4)
  
  
  
  
  raw_sset <- mutateSesame(sset=raw_sset, method="detectionPnegEcdf", 
                           verbose=4)
  
  # Sigs::Data
  #
  raw_sset_dat_csv <- file.path(opt$outDir, 'raw_sset.dat.csv.gz')
  raw_sigs_dat_tib <- ssetToSigsTib(
    sset=raw_sset, man=man_tib, 
    by="Probe_ID", type="Probe_Type", des="Probe_Class", 
    percision=2, sort=TRUE, save=TRUE, csv=raw_sset_dat_csv, verbose=4)
  
  raw_sigs_dat_tib %>% 
    dplyr::group_by(Probe_Type,Probe_Class) %>% 
    dplyr::summarise(Count=n(), .groups='drop')
  
  # Sigs::Summary
  #
  raw_sset_sum_csv <- file.path(opt$outDir, 'raw_sset.sum.csv.gz')
  raw_sigs_sum_tib <- sigsTibToSummary(
    tib=raw_sigs_dat_tib, 
    by="Probe_ID", type="Probe_Type", des="Probe_Class",
    save=TRUE, csv=raw_sset_sum_csv, verbose=4)
  
  raw_sigs_ssheet_tib <- sigsSumToSSheet(
    tib=raw_sigs_sum_tib, metric=metric,
    by="Probe_ID", type="Probe_Type", des="Probe_Class", 
    verbose=4)
  
  # Calls::
  #
  raw_call_tib <- raw_sset %>% ssetToCallTib(workflow='raw', verbose=4)
  
  # Call::betas
  #
  
  betas_cols <- raw_call_tib %>% 
    dplyr::select(Probe_ID,dplyr::ends_with('_beta')) %>% names()
  
  raw_beta_sum_tib <- raw_call_tib %>% 
    dplyr::select(dplyr::all_of(betas_cols)) %>%
    sigsTibToSummary(man=man_tib, 
                     by="Probe_ID", type="Probe_Type", des="Probe_Design", 
                     verbose=4)
  
  raw_beta_ssheeet_tib <- sigsSumToSSheet(
    tib=raw_beta_sum_tib, metric=metric,
    by="Probe_ID", type="Probe_Type", des="Probe_Design", 
    verbose=4)
  
  # Call::pvals
  #
  pvals_cols <- raw_call_tib %>% 
    dplyr::select(!dplyr::ends_with('_beta')) %>% names()
  
  raw_pval_sum_tib <- raw_call_tib %>% 
    dplyr::select(dplyr::all_of(pvals_cols)) %>%
    sigsTibToSummary(man=man_tib, 
                     by="Probe_ID", type="Probe_Type", des="Probe_Design", 
                     cutoff=0.1,
                     verbose=4)
  
  raw_pval_ssheeet_tib <- sigsSumToSSheet(
    tib=raw_pval_sum_tib, metric='pass_perc',
    by="Probe_ID", type="Probe_Type", des="Probe_Design", 
    verbose=4)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Summary::
  sset_tib <- ssetToTib(rdat$rsum_list$sset_dat)
  sset_tib %>% dplyr::group_by(Probe_Design) %>% dplyr::summarise(Count=n())
  
  sset_sum_tib <- sset_tib %>% sigTibToSummary(man=rdat$sman)
  sset_sum_tib %>% 
    tidyr::unite(Probe_Type, Probe_Type,Probe_Design, sep='_') %>%
    dplyr::select(Probe_Type, ends_with(metric)) %>% tidyr::gather(Metric, Value, -Probe_Type) %>%
    tidyr::unite(Key, Probe_Type,Metric, sep='_') %>%
    tidyr::spread(Key, Value)
  
  
  
  cur_sset <- rdat$cur_sset
  ses_man_tib <- rdat$sman
  ww <- 3
  cur_workflow <- 'ind'
  
  cur_sset_csv <- NULL
  cur_ssum_csv <- NULL
  
  sum_list <- summarySSET_workflow(
    sset=cur_sset, man=ses_man_tib, idx=ww, workflow=cur_workflow,
    writeSset=opt$writeSset,sset_csv=cur_sset_csv,
    writeSsum=opt$writeSsum,ssum_csv=cur_ssum_csv,
    minNegPval=opt$minNegPval,minOobPval=opt$minOobPval,
    percisionBeta=opt$percisionBeta, 
    percisionPval=opt$percisionPval, 
    percisionSigs=opt$percisionSigs,
    verbose=4)
  
  
  rdat$cur_sset@extra$pvals$pOOBAH %>% head()
  rdat$cur_sset@extra$pvals$pOOBAH <- NULL
  
  new_sset <- rdat$raw_sset
  new_sset@extra$pvals$pOOBAH <- NULL
  new_sset@extra$pvals <- NULL
  new_sset <- pOOBAH2(new_sset, force=TRUE)
  new_sset@extra$pvals$pOOBAH %>% head(n=30)
  
  two_sset <- rdat$cur_sset
  two_sset@extra$pvals$pOOBAH <- NULL
  two_sset@extra$pvals <- NULL
  two_sset <- pOOBAH2(two_sset, force=TRUE)
  two_sset@extra$pvals$pOOBAH %>% head(n=30)
  
  rdat$raw_sset@extra$pvals$pOOBAH %>% head(n=30)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Old Auto Sample::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  # Prep all_call_tib to match previous format::
  #  Probe_ID       M        U DESIGN COLOR_CHANNEL col   Probe_Type Probe_Source Next_Base
  auto_call <- ses_man_tib %>% dplyr::rename(Design_Type=DESIGN) %>% 
    dplyr::select(Probe_ID, Design_Type, Probe_Type) %>%
    dplyr::left_join(all_call_tib, by="Probe_ID") # %>% dplyr::select(!!auto_negs_key, !!auto_beta_key)
  
  print(auto_call)
  auto_ref <- ref
  if (grep("_beta", names(ref)) %>% length() == 0) {
    auto_ref <- ref %>% purrr::set_names(paste(names(.),'beta', sep=del) )
    names(auto_ref)[1] <- 'Probe_ID'
  }
  auto_data <- sampleDetect(can=auto_call, ref=auto_ref, minPval=opt$minNegPval, minDelta=opt$minDeltaBeta,
                            dname='Design_Type', pname='Probe_Type', ptype='cg',
                            jval='Probe_ID', field=auto_beta_key, pval=auto_negs_key, suffix='beta', del=del,
                            outDir=opt$outDir, sname=out_name, plotMatrix=opt$plotAuto, writeMatrix=opt$writeAuto,
                            dpi=opt$dpi, format=opt$plotFormat, datIdx=4, non.ref=non_ref,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
  ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(auto_data)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Add Requeue Flags::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
# rdat$ssheet_tib %>% dplyr::select(pOOBAH_cg_1_pass_perc_0, pOOBAH_cg_2_pass_perc_0, PnegEcdf_cg_1_pass_perc_0, PnegEcdf_cg_2_pass_perc_0)

if (FALSE) {
  requ_sum_tib <- tibble::tibble(
    Requeue_Flag_Oob=case_when(ssheet_tib$Poob_Pass_0_Perc <= opt$minOobPerc ~ TRUE, TRUE ~ FALSE),
    Requeue_Flag_Neg=case_when(ssheet_tib$Negs_Pass_0_Perc <= opt$minNegPerc ~ TRUE, TRUE ~ FALSE),
    Requeue_Pass_Perc_Oob=opt$minOobPerc,
    Requeue_Pass_Perc_Neg=opt$minNegPerc)
  
  ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(requ_sum_tib) %>% 
    dplyr::select(Requeue_Flag_Oob, Requeue_Flag_Neg, Requeue_Pass_Perc_Oob, Requeue_Pass_Perc_Neg, everything())
  
  ssheet_ncols <- ssheet_tib %>% base::ncol()
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Binding Sample Sheet (with sample-requeue) ssheet_ncols={ssheet_ncols}.{RET}{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Add Sentrix Names to Calls Files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
if (FALSE && opt$addSentrixID) {
  org_col_key1 <- colnames(all_call_tib)[1]
  all_call_tib <- all_call_tib %>% purrr::set_names(paste(names(.), ssheet_tib$Sentrix_Name[1], sep='_') )
  colnames(all_call_tib)[1] <- org_col_key1
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Add to Full Return Values:: R&D Testing ONLY::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (retData) {
  ret$cur_sset     <- cur_sset
  ret$cur_sset_sum <- sset_sum_tib
  ret$all_call_tib <- all_call_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Sesame Mutation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
initSesameRaw = function(prefix, platform, manifest, load=FALSE, save=FALSE, rds=NULL, 
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'initSesameRaw'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; platform={platform}, prefix={prefix}.{RET}"))
  
  ret_cnt <- 0
  stime <- system.time({
    if (load && !is.null(rds) && file.exists(rds)) {
      if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading RDS={rds}.{RET}"))
      sset <- readr::read_rds(rds)
    } else {
      if (verbose>=vt+5) head(manifest) %>% print()
      sset <- sesame::readIDATpair(prefix, platform=platform, manifest=manifest)
      if (save && !is.null(rds)) {
        if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing RDS={rds}.{RET}"))
        readr::write_rds(sset, rds, compress="gz")
      }
    }
    
    ret_cnt <- sset %>% slotNames() %>% length()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} sset(slots={ret_cnt})={RET}"))
      print(sset)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  sset
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Sesame SSET Sample Sheet Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: rewrite callToSSheet
#

callToSSheet_old1 = function(call, idx, key, pre=NULL, minNegPval, minOobPval, 
                            percision_beta=-1, percision_pval=-1, del='_', 
                            onlyCG=TRUE, id='Probe_ID',
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'callToSSheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; idx={idx}, key={key}, del={del}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    id <- id %>% rlang::sym()
    if (onlyCG) call <- call %>% dplyr::filter(stringr::str_starts(!!id, 'cg'))
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} call={RET}"))
    if (verbose>=vt+4) head(call) %>% print()
    
    beta_tib <- NULL
    beta_key_str <- paste('Beta',idx,'Method', sep=del)
    beta_val_str <- paste('Beta',idx,'Mean', sep=del)
    beta_key_sym <- rlang::sym(beta_key_str)
    beta_val_sym <- rlang::sym(beta_val_str)
    beta_tib <- tibble::tibble(
      !!beta_key_sym := key, 
      !!beta_val_sym := dplyr::select(call, ends_with('_beta')) %>% 
        dplyr::summarise_all(list(mean=mean), na.rm=TRUE) %>% 
        dplyr::pull() %>% round(percision_beta) )
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} beta_key_str={beta_key_str}, beta_val_str={beta_val_str}.{RET}"))
    
    negs_tib <- NULL
    negs_key_str <- paste('Negs_Pass',idx,'Method', sep=del)
    negs_val_str <- paste('Negs_Pass',idx,'Perc', sep=del)
    negs_key_sym <- rlang::sym(negs_key_str)
    negs_val_sym <- rlang::sym(negs_val_str)
    
    negs_tib <- tibble::tibble(
      !!negs_key_sym := key, 
      !!negs_val_sym := dplyr::select(call, ends_with('_negs')) %>%
        dplyr::summarise_all(list(pass_perc=cntPer_lte), min=minNegPval) %>%
        dplyr::pull() %>% round(percision_pval) )
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} negs_key_str={negs_key_str}, negs_val_str={negs_val_str}.{RET}"))
    
    poob_tib <- NULL
    poob_key_str <- paste('Poob_Pass',idx,'Method', sep=del)
    poob_val_str <- paste('Poob_Pass',idx,'Perc', sep=del)
    poob_key_sym <- rlang::sym(poob_key_str)
    poob_val_sym <- rlang::sym(poob_val_str)
    
    poob_tib <- tibble::tibble(!!poob_key_sym := key, !!poob_val_sym := dplyr::select(call, ends_with('_poob')) %>% 
                                 dplyr::summarise_all(list(pass_perc=cntPer_lte), min=minOobPval) %>%
                                 dplyr::pull() %>% round(percision_pval) )
    
    ret_tib <- dplyr::bind_cols(pre, beta_tib, negs_tib, poob_tib)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}


# rdat$rsum_list$call_dat %>% callToSSheet(idx=1, key='bla', minNegPval=0.02, minOobPval=0.1, onlyCG=FALSE)
callToSSheet_old = function(call, idx, key, pre=NULL, minNegPval, minOobPval, 
                            percisionBeta=0, percisionPval=0, del='_', 
                            onlyCG=TRUE, id='Probe_ID',
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'callToSSheet_old'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; idx={idx}, key={key}, del={del}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    id <- id %>% rlang::sym()
    if (onlyCG) call <- call %>% dplyr::filter(stringr::str_starts(!!id, 'cg'))
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} call={RET}"))
    if (verbose>=vt+4) head(call) %>% print()
    
    beta_tib <- NULL
    beta_key_str <- paste('Beta',idx,'Method', sep=del)
    beta_val_str <- paste('Beta',idx,'Mean', sep=del)
    beta_key_sym <- rlang::sym(beta_key_str)
    beta_val_sym <- rlang::sym(beta_val_str)
    beta_tib <- tibble::tibble(
      !!beta_key_sym := key, 
      !!beta_val_sym := dplyr::select(call, ends_with('_beta')) %>% 
        dplyr::summarise_all(list(mean=mean), na.rm=TRUE) %>% 
        dplyr::pull() %>% round(percisionBeta) )
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} beta_key_str={beta_key_str}, beta_val_str={beta_val_str}.{RET}"))
    
    negs_tib <- NULL
    negs_key_str <- paste('Negs_Pass',idx,'Method', sep=del)
    negs_val_str <- paste('Negs_Pass',idx,'Perc', sep=del)
    negs_key_sym <- rlang::sym(negs_key_str)
    negs_val_sym <- rlang::sym(negs_val_str)
    
    negs_tib <- tibble::tibble(
      !!negs_key_sym := key, 
      !!negs_val_sym := dplyr::select(call, ends_with('_negs')) %>%
        dplyr::summarise_all(list(pass_perc=cntPer_lte), min=minNegPval) %>%
        dplyr::pull() %>% round(percisionPval) )
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} negs_key_str={negs_key_str}, negs_val_str={negs_val_str}.{RET}"))
    
    poob_tib <- NULL
    poob_key_str <- paste('Poob_Pass',idx,'Method', sep=del)
    poob_val_str <- paste('Poob_Pass',idx,'Perc', sep=del)
    poob_key_sym <- rlang::sym(poob_key_str)
    poob_val_sym <- rlang::sym(poob_val_str)
    
    poob_tib <- tibble::tibble(!!poob_key_sym := key, !!poob_val_sym := dplyr::select(call, ends_with('_poob')) %>% 
                                 dplyr::summarise_all(list(pass_perc=cntPer_lte), min=minOobPval) %>%
                                 dplyr::pull() %>% round(percisionPval) )
    
    ret_tib <- dplyr::bind_cols(pre, beta_tib, negs_tib, poob_tib)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame SSET To Inference Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssetToInferences_old = function(sset, idx, key, pre=NULL, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToInferences'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    gct_tib <- NULL
    inf_key <- paste('GCT',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('GCT',idx,'Score', sep=del) %>% rlang::sym()
    gct_tib <- tibble::tibble(!!inf_key := key, !!inf_val := safeGCT(sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) )
    
    sex_tib <- NULL
    inf_key <- paste('Sex',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('Sex',idx,'Call', sep=del) %>% rlang::sym()
    sex_tib <- tibble::tibble(!!inf_key := key, !!inf_val := safeSex(sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) )
    
    kar_tib <- NULL
    inf_key <- paste('Karyotype',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('Karyotype',idx,'Call', sep=del) %>% rlang::sym()
    kar_tib <- tibble::tibble(!!inf_key := key, !!inf_val := safeSexKaryo(sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) )
    
    eth_tib <- NULL
    inf_key <- paste('Ethnicity',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('Ethnicity',idx,'Call', sep=del) %>% rlang::sym()
    eth_tib <- tibble::tibble(!!inf_key := key, !!inf_val := safeEthnicity(sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) )
    
    IGR_Swap_Cnt <- 0
    IRG_Swap_Cnt <- 0
    if (!is.null(sesame::extra(sset)[['IRR']]) && 
        !is.null(sesame::extra(sset)[['IGG']])) {
      IGR_Swap_Cnt <- length(which(!sesame::extra(sset)[['IGG']]))
      IRG_Swap_Cnt <- length(which(!sesame::extra(sset)[['IRR']]))
    }
    igr_tib <- NULL
    inf_key <- paste('IGR_Swap',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('IGR_Swap',idx,'Count', sep=del) %>% rlang::sym()
    igr_tib <- tibble::tibble(!!inf_key := key, !!inf_val := IGR_Swap_Cnt)
    
    irg_tib <- NULL
    inf_key <- paste('IRG_Swap',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('IRG_Swap',idx,'Count', sep=del) %>% rlang::sym()
    irg_tib <- tibble::tibble(!!inf_key := key, !!inf_val := IRG_Swap_Cnt)
    
    ret_tib <- dplyr::bind_cols(pre, gct_tib, sex_tib, kar_tib, eth_tib, igr_tib, irg_tib)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame SSET To Inference Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssetToPredictions_old = function(sset, idx, key, pre=NULL, del='_', fresh=FALSE,
                                 quality.mask = FALSE, 
                                 sum.TypeI = FALSE,
                                 verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToPredictions'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}      quality.mask={quality.mask}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}         sum.TypeI={sum.TypeI}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}              sset={RET}"))
  if (verbose>=vt+4) print(sset)
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{RET}{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} sset={RET}"))
    if (verbose>=vt+4) print(sset)
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{RET}{RET}"))
    
    if (fresh || is.null(sesame::extra(sset)[['betas']]) ) {
      # sesame::extra(sset)[['betas']] <-
      #   sesame::getBetas(sset=sset, mask = quality.mask,sum.TypeI = sum.TypeI)
      sesame::extra(sset)[['betas']] <-
        getBetas2(sset=sset, mask = quality.mask,sum.TypeI = sum.TypeI)
    }
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} beta={RET}"))
    if (verbose>=vt+4) head(sesame::extra(sset)[['betas']]) %>% print()
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{RET}{RET}"))
    
    skn_tib <- NULL
    inf_key <- paste('SkinAge',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('SkinAge',idx,'Score', sep=del) %>% rlang::sym()
    skn_tib <- tibble::tibble(
      !!inf_key := key, 
      !!inf_val := safeSkinAge(sesame::extra(sset)[['betas']], 
                               verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) )
    
    phn_tib <- NULL
    inf_key <- paste('PhenoAge',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('PhenoAge',idx,'Score', sep=del) %>% rlang::sym()
    phn_tib <- tibble::tibble(
      !!inf_key := key, 
      !!inf_val := safePhenoAge(sesame::extra(sset)[['betas']], 
                                verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) )
    
    ret_tib <- dplyr::bind_cols(pre, skn_tib, phn_tib)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Sesame SSET to Tibs Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssetToSigsTib_old = function(sset, add, name, del='_', rmAdd=TRUE, 
                             save=FALSE, csv=NULL,
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToSigsTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, name={name}!{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    datTag <- 'sig'
    grnTag <- paste(name,'G', sep=del) %>% rlang::sym()
    redTag <- paste(name,'R', sep=del) %>% rlang::sym()
    inbTag <- paste(name,'Inb_Col', sep=del) %>% rlang::sym()
    swpTag <- paste(name,'Swap', sep=del) %>% rlang::sym()
    
    IG_tib <- sset@IG %>% 
      tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='G')
    
    OR_tib <- sset@oobR %>% 
      tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='G')
    
    IR_tib <- sset@IR %>% 
      tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='R')
    
    OG_tib <- sset@oobG %>% 
      tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='R')
    
    I1 <- dplyr::bind_rows(
      dplyr::inner_join(IG_tib,OR_tib, by=c('Probe_ID', 'Inb_Col', 'Design_Type')),
      dplyr::inner_join(OG_tib,IR_tib, by=c('Probe_ID', 'Inb_Col', 'Design_Type')) ) %>%
      dplyr::mutate(Design_Type=paste0(Design_Type,'I')) %>%
      dplyr::select(Probe_ID,Design_Type,Inb_Col,everything())
    
    I1 <- add %>% dplyr::filter(Design_Type!='II') %>% 
      dplyr::left_join(I1, by=c("Probe_ID", "Design_Type") ) %>%
      dplyr::mutate(Swap=case_when(
        Man_Col!=Inb_Col ~ TRUE,
        Man_Col==Inb_Col ~ FALSE,
        TRUE ~ NA)) %>%
      dplyr::select(-Inb_Col) %>%
      dplyr::select(Probe_ID, Man_Col, Design_Type, Probe_Type, Swap, everything())
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Infinium II::
    I2 <- sset@II %>% 
      tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      dplyr::rename(!!redTag :=U, !!grnTag :=M) %>%
      dplyr::mutate(Inb_Col=NA, Swap=NA, Design_Type='II')
    
    I2 <- add %>% dplyr::filter(Design_Type=='II') %>% 
      dplyr::left_join(I2, by=c("Probe_ID", "Design_Type") ) %>%
      dplyr::select(-Inb_Col) %>%
      dplyr::select(Probe_ID, Man_Col, Design_Type, Probe_Type, Swap, everything())
    
    ret_tib <- dplyr::bind_rows(I1, I2) %>%
      dplyr::rename(!!swpTag := Swap) %>%
      dplyr::arrange(Probe_ID)
    
    if (rmAdd && "Address" %in% names(ret_tib)) 
      ret_tib <- ret_tib %>% dplyr::select(-Address)
    
    if (!is.null(save) && save==TRUE && !is.null(csv)) {
      csv_dir <- base::basename(csv)
      if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
      
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing (percision={percision}) CSV={csv}.{RET}"))
      readr::write_csv(ret_tib, csv)
    }
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

ssetToSigsTib_old = function(sset, man=NULL, by="Probe_ID", 
                             type="Probe_Type", des="Probe_Design", 
                             percision=-1, sort=FALSE, 
                             save=FALSE, csv=NULL,
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToSigsTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    des <- des %>% rlang::sym()
    
    ret_tib <- dplyr::bind_rows(
      sset@IG   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='IG'),
      sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='OG'),
      
      sset@IR   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='IR'),
      sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='OR'),
      
      sset@II   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='II')
    ) %>% dplyr::select(!!by, !!des, everything())
    
    if (percision!=-1) ret_tib <- ret_tib %>% dplyr::mutate_if(is.numeric, list(round), percision)
    if (!is.null(man)) ret_tib <- dplyr::select(man, !!by, !!type) %>% dplyr::right_join(ret_tib, by=by )
    
    by  <- by %>% rlang::sym()
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(!!by)
    
    if (!is.null(save) && save==TRUE && !is.null(csv)) {
      csv_dir <- base::basename(csv)
      if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
      
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing (percision={percision}) CSV={csv}.{RET}"))
      readr::write_csv(ret_tib, csv)
    }
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#            EPIC hg19 Infinium I Inferred SNP Table Construction::
#
#               Should move this scratch to somewhere else::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  # 450k Genome Studio::
  #
  # hm45_man_csv <- '/Users/bretbarnes/Documents/data/manifests/HumanMethylation450_15017482_v.1.2.csv.gz'
  # hm45_man_lst <- loadManifestGenomeStudio(file=hm45_man_csv, verbose=opt$verbose+4)
  # hm45_man_tib <- hm45_man_lst$man %>% dplyr::select(IlmnID,Genome_Build,CHR,MAPINFO,Strand)
  # hm45_ctl_tib <- hm45_man_lst$ctl %>% dplyr::mutate(Address=as.double(Address))
  
  #
  # SNP Table
  #
  epicI_ses_grs <- sesameDataPullVariantAnno_InfiniumI(platform = "EPIC")
  epicR_ses_grs <- sesameDataPullVariantAnno_SNP()
  
  epicI_ses_tib <- epicI_ses_grs %>% 
    as.data.frame() %>% tibble::rownames_to_column(var="IlmnID") %>% tibble::as_tibble()
  
  # Five Examples by hand::
  #  gs_ses_csv <- '/Users/bretbarnes/Documents/tools/notes/inferred-snps.GS-vs-Ses.csv'
  #  gs_ses_tib <- readr::read_csv(gs_ses_csv)
  #
  #
  
  # epic_man_csv <- file.path(par$datDir, 'manifest/base/EPIC-B4.manifest.sesame-base.cpg-sorted.csv.gz')
  epic_man_csv <- '/Users/bretbarnes/Documents/data/manifests/MethylationEPIC_v-1-0_B4.csv.gz'
  epic_man_lst <- loadManifestGenomeStudio(file=epic_man_csv, verbose=opt$verbose+4)
  epic_ctl_tib <- epic_man_lst$ctl %>% dplyr::mutate(Address=as.double(Address))
  epic_man_tib <- epic_man_lst$man %>% 
    dplyr::select(IlmnID,Genome_Build,CHR,MAPINFO,Strand,Next_Base,Infinium_Design_Type) %>%
    dplyr::filter(Infinium_Design_Type=='I') %>% 
    dplyr::mutate(CHR=paste0('chr',CHR)) %>% dplyr::arrange(CHR,MAPINFO)
  
  #
  # Validation Shown Below::
  #
  join_man_tib <- epic_man_tib %>% 
    dplyr::inner_join(epicI_ses_tib, by="IlmnID", suffix=c("_Man", "_Ses")) %>%
    dplyr::mutate(REF_CON=dplyr::case_when(
      strand=='+' & REF=='G' ~ 'A',
      strand=='-' & REF=='C' ~ 'T',
      TRUE ~ REF)
    )
  join_sum_tib <- join_man_tib %>% dplyr::filter(CHR==seqnames) %>% 
    dplyr::mutate(Pos_Dif=MAPINFO-start) %>% 
    dplyr::group_by(Pos_Dif,Strand,strand,Next_Base,REF_CON) %>% 
    dplyr::summarise(Count=n(), .groups='drop')
  
  #
  # Better, but NOT working attempt below::
  #
  if (FALSE) {
    join_man_tib3 <- epic_man_tib %>% 
      dplyr::inner_join(epicI_ses_tib, by="IlmnID", suffix=c("_Man", "_Ses")) %>%
      dplyr::mutate(
        REF_RC=dplyr::case_when(
          strand=='-' ~ revCmp(REF),
          TRUE ~ REF),
        REF_CO=dplyr::case_when(
          strand=='+' & REF_RC=='G' ~ 'A',
          strand=='-' & REF_RC=='C' ~ 'T',
          TRUE ~ REF)
      )
    join_sum_tib3 <- join_man_tib3 %>% 
      dplyr::filter(CHR==seqnames) %>% 
      dplyr::mutate(Pos_Dif=MAPINFO-start) %>% 
      dplyr::group_by(Pos_Dif,Strand,strand,Next_Base,REF_CO) %>% 
      dplyr::summarise(Count=n(), .groups='drop')
  }
  
  #
  # hub19 doesn't have the most up to date dbSNP
  #
  hub_hg19 <- subset(hub, (hub$species == "Homo sapiens") & (hub$genome == "hg19"))
  length(hub_hg19)
  hub_hg19$title[grep("SNP", hub_hg19$title)]
  
  hub_snp137 <- subset(hub_hg19, title=='Common SNPs(137)')
  
  hub_snp137_grs <- hub_snp137[['AH5105']]
  hub_snp137_tib <- hub_snp137_grs %>% as.data.frame() %>% tibble::as_tibble()
  
  #
  # Check overlap of Sesame and dbSNP(137) list
  #
  epicI_ses_tib %>% dplyr::inner_join(hub_snp137_tib, by=c("rs"="name"), suffix=c("_sesI", "_db137"))
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#             Every thing below is not used and can be removed::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Sesame Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #



sigsSumToSSheet_old = function(tib, metric='avg',
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sigsSumToSSheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    pt_tib <- tib %>% split(.$Probe_Type)
    
    tib_names <- pt_tib %>% names()
    core_pts <- c('cg', 'ch', 'rs', 'rp', 'mu')
    
    ss_tib <- NULL
    for (pt in tib_names) {
      if (!c(pt) %in% core_pts) next
      
      nrows = base::nrow(pt_tib[[pt]])
      for (ii in c(1:nrows)) {
        pt_val <- pt_tib[[pt]]$Probe_Type[ii]
        pd_val <- pt_tib[[pt]]$Probe_Design[ii]
        
        key    <- paste(pt_val, pd_val, sep='_')
        valM   <- pt_tib[[pt]]$M_avg[ii];
        valU   <- pt_tib[[pt]]$U_avg[ii];
        
        if (pd_val=='II') {
          keyM <- paste(key,'G',metric, sep='_') %>% rlang::sym()
          keyU <- paste(key,'R',metric, sep='_') %>% rlang::sym()
        } else {
          keyM <- paste(key,'M',metric, sep='_') %>% rlang::sym()
          keyU <- paste(key,'U',metric, sep='_') %>% rlang::sym()
        }
        new_tib <- tibble::tibble(!!keyM := valM, !!keyU := valU) %>%
          dplyr::mutate_if(is.numeric, list(round), 0)
        
        ss_tib <- ss_tib %>% bind_cols(new_tib)
      }
      # cat(glue::glue("1st: pt={pt}.{RET}"))
    }
    # ss_tib %>% as.data.frame() %>% print()
    
    for (pt in tib_names) {
      if (c(pt) %in% core_pts) next
      
      nrows = base::nrow(pt_tib[[pt]])
      for (ii in c(1:nrows)) {
        pt_val <- pt_tib[[pt]]$Probe_Type[ii]
        pd_val <- pt_tib[[pt]]$Probe_Design[ii]
        
        # key    <- paste(pt_val, pd_val, sep='_')
        key    <- paste(pt_val, sep='_')
        valM   <- pt_tib[[pt]]$M_avg[ii];
        valU   <- pt_tib[[pt]]$U_avg[ii];
        
        if (pd_val=='II') {
          keyM <- paste(key,'G',metric, sep='_') %>% rlang::sym()
          keyU <- paste(key,'R',metric, sep='_') %>% rlang::sym()
        } else {
          keyM <- paste(key,'M',metric, sep='_') %>% rlang::sym()
          keyU <- paste(key,'U',metric, sep='_') %>% rlang::sym()
        }
        new_tib <- tibble::tibble(!!keyM := valM, !!keyU := valU) %>%
          dplyr::mutate_if(is.numeric, list(round), 0)
        
        ss_tib <- ss_tib %>% bind_cols(new_tib)
      }
    }
    ret_tib <- ss_tib %>% dplyr::mutate_if(is.numeric, as.integer)
    if (verbose>=vt+4) head(ret_tib) %>% as.data.frame() %>% print()
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

sigsToSummary = function(tib, name, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sigsToSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    cntTag <- 'Count' %>% rlang::sym()
    swapTag <- 'Swap' %>% rlang::sym()
    scntTag <- 'Swap_Count' %>% rlang::sym()
    
    if (! is.null(name)) {
      cntTag  <- paste(name,'Count', sep=del) %>% rlang::sym()
      swapTag <- paste(name,'Swap', sep=del) %>% rlang::sym()
      scntTag <- paste(name,'Swap_Count', sep=del) %>% rlang::sym()
    }
    
    # Remove Address from summary if address exists::
    if ("Address" %in% names(tib)) tib <- tib %>% dplyr::select(-Address) 
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #   Summarize Data::
    cnts <- tib %>% dplyr::group_by(Design_Type, Probe_Type, Man_Col) %>%
      summarise(!!cntTag :=n(), !!scntTag := count(!!swapTag==TRUE), .groups='drop')
    sums <- tib %>% dplyr::group_by(Design_Type, Probe_Type, Man_Col) %>%
      summarise_if(is.numeric, list(min=min, avg=mean, med=median, max=max), 
                   na.rm=TRUE, .groups='drop' )
    
    ret_tib <- cnts %>% 
      dplyr::full_join(sums, by=c("Design_Type","Probe_Type", "Man_Col")) %>%
      dplyr::arrange(Probe_Type)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

getSesameManifest = function(man, sig, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getSesameManifest'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    ret_tib <- man %>% 
      dplyr::right_join(sig, by=c("U"="Address")) %>% 
      dplyr::filter(!is.na(Probe_ID)) %>%
      dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base) %>%
      dplyr::arrange(Probe_ID)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Sesame Swap Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: This should be moved to the same location as similar functions...
#
joinSsetTibInfI = function(tib,
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'joinSsetTibInfI'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_tib <- NULL
  tib_spl <- NULL
  tib_grn <- NULL
  tib_red <- NULL
  
  stime <- system.time({
    tib_spl <- tib %>% dplyr::filter(Probe_Design != 'II') %>% split(.$Probe_Design)
    
    tib_grn <- dplyr::inner_join(tib_spl$IG, tib_spl$OR, by="Probe_ID", suffix=c('_Inb', '_Oob')) %>% 
      dplyr::mutate(Sum_Inb=M_Inb+U_Inb, Sum_Oob=M_Oob+U_Oob, Dif_InOut=Sum_Inb-Sum_Oob)
    
    tib_red <- dplyr::inner_join(tib_spl$IR, tib_spl$OG, by="Probe_ID", suffix=c('_Inb', '_Oob')) %>% 
      dplyr::mutate(Sum_Inb=M_Inb+U_Inb, Sum_Oob=M_Oob+U_Oob, Dif_InOut=Sum_Inb-Sum_Oob)
    
    ret_tib <- dplyr::bind_rows(tib_grn, tib_red) %>% dplyr::arrange(Probe_ID)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

if (FALSE) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Add Swapped Summary Percentages::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (!is.null(opt$skipSwap) && !opt$skipSwap) {
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Building raw_sset_tib...{RET}") )
    
    cur_sset_tib <- NULL
    raw_sset_tib <- NULL
    opt$writeSsetRaw <- FALSE
    raw_sset_tib <- sset2tib(sset=raw_sset, by="Probe_ID", des="Probe_Design",  
                             percision=opt$percisionSigs, sort=FALSE, 
                             save=opt$writeSsetRaw, csv=raw_sset_csv, 
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    if (retData) ret$raw_sset_tib <- raw_sset_tib
  }
  
  if (!is.null(raw_sset_tib) && !is.null(cur_sset_tib)) {
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Inferred Sample Swapped Stats.{RET}"))
    
    swap_sum_tib <- NULL
    swap_sum_tib <- dplyr::inner_join(joinSsetTibInfI(tib=raw_sset_tib), 
                                      joinSsetTibInfI(tib=cur_sset_tib), 
                                      by="Probe_ID", suffix=c("_Raw", "_Cur") ) %>%
      dplyr::mutate(isSwap=dplyr::case_when(
        Probe_Design_Inb_Raw==Probe_Design_Inb_Cur & Probe_Design_Oob_Raw==Probe_Design_Oob_Cur ~ 'Reference',
        Probe_Design_Inb_Raw!=Probe_Design_Inb_Cur & Probe_Design_Oob_Raw!=Probe_Design_Oob_Cur ~ 'Alternate',
        TRUE ~ NA_character_
      )) %>% dplyr::group_by(Probe_Design_Inb_Raw,isSwap) %>% 
      dplyr::summarise(Count=n(), .groups='drop') %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(Total=sum(Count)) %>% 
      tidyr::unite(Type, Probe_Design_Inb_Raw, isSwap, sep='_') %>% 
      dplyr::mutate(Perc=round(100*Count/Total, 3)) %>% 
      dplyr::select(Type, Perc) %>% tidyr::spread(Type, Perc) %>%
      purrr::set_names(paste(names(.),'Perc', sep='_') )
    
    ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(swap_sum_tib)
    
    ssheet_ncols <- ssheet_tib %>% base::ncol()
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Binding Sample Sheet (with channel-swap-stats) ssheet_ncols={ssheet_ncols}.{RET}{RET}"))
  }
}

safeGCT_org = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeGCT_org'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NA
    try_str <- 'pass'
    
    try(val <- sesame::bisConversionControl(sset), silent = TRUE)
    if (is.na(val)) try_str <- 'fail'
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

ssetToPvalTib_old = function(sset, method, name, percision=0, 
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToPvalTib_old'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; method={method}, name={name}.{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  
  # if (verbose>=vt+4) head(sset@pval[method]) %>% print()
  name_sym <- rlang::sym(name)
  
  stime <- system.time({
    if (!is.null(sset@pval) && !is.null(sset@pval[method]) )
      ret_tib <- tibble::tibble(Probe_ID=names(sset@pval[method]), !!name_sym := 1.0)
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. Warning; returning; pval=1.0.{RET}"))
    if (verbose>=vt+4) ret_tib %>% head() %>% print()
    
    ret_tib <- sset@pval[method] %>% tibble::enframe(name='Probe_ID', value=name)
    if (percision!=0) ret_tib <- ret_tib %>% dplyr::mutate_if(purrr::is_double, round, percision)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   OLD Sesame SSET Manipulation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sesameStepAbbreviation = function(x, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sesameStepAbbreviation'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (x=='raw') return('R')
  if (x=='dyeBiasCorrTypeINorm') return('D')
  if (x=='detectionPnegEcdf') return('N')
  if (x=='pOOBAH') return('P')
  if (x=='noob') return('N')
  if (x=='noobsb') return('S')
  if (x=='inferTypeIChannel') return('I')
  stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported Sesame Abbreviation={x}!!!{RET}{RET}") )
  
  return('U')
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     SSET to BeadSET Conversion Methods::
#
#                                 NOT USED!
#                                 NOT USED!
#                                 NOT USED!
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sset2bset = function(sset, addSigs=TRUE, addNegs=TRUE, addPoob=TRUE, addBeta=TRUE,
                     
                     # getBeta Parameters::
                     quality.mask=FALSE, nondetection.mask=FALSE, 
                     mask.use.tcga=FALSE, pval.threshold=1, sum.TypeI=FALSE,
                     as.enframe=TRUE,
                     round_dat=TRUE, round_pval=6, round_beta=4,
                     verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sset2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  bset <- NULL
  stime <- system.time({
    sigs <- dplyr::bind_rows(
      sset@IG   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='IG'),
      sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='OG'),
      
      sset@IR   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='IR'),
      sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='OR'),
      
      sset@II   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col=NA)
    )
    if (round_dat) sigs <- sigs %>% dplyr::mutate_if(is.numeric, list(as.integer))
    
    if (addSigs) { bset <- sigs
    } else { bset <- dplyr::select(sigs, 'Probe_ID', 'col') }
    
    if (addNegs) {
      negs <- sset %>% sesame::detectionPnegEcdf() %>% sesame::pval() %>% tibble::enframe(name='Probe_ID', value='negs')
      if (round_dat && round_pval>0) negs <- negs %>% dplyr::mutate_if(is.numeric, list(round), round_pval)
      bset <- bset %>% dplyr::left_join(negs, by='Probe_ID')
    }
    
    if (addPoob) {
      poob <- sset %>% sesame::pOOBAH() %>% sesame::pval() %>% tibble::enframe(name='Probe_ID', value='poob')
      if (round_dat && round_pval>0) poob <- poob %>% dplyr::mutate_if(is.numeric, list(round), round_pval)
      bset <- bset %>% dplyr::left_join(poob, by='Probe_ID')
    }
    
    if (addBeta) {
      beta <- sesame::getBetas(sset=sset, quality.mask=quality.mask, 
                               nondetection.mask=nondetection.mask, 
                               mask.use.tcga=mask.use.tcga, 
                               pval.threshold=pval.threshold, 
                               sum.TypeI=sum.TypeI)
      if (as.enframe) beta <- beta %>% tibble::enframe(name='Probe_ID', value='beta')
      if (round_dat && round_beta>0) beta <- beta %>% dplyr::mutate_if(is.numeric, list(round), round_beta)
      bset <- bset %>% dplyr::left_join(beta, by='Probe_ID')
    }
    
    # bset <- bset %>% 
    #   dplyr::left_join(negs, by='Probe_ID') %>%
    #   dplyr::left_join(poob, by='Probe_ID') %>%
    #   dplyr::left_join(beta, by='Probe_ID')
  })
  bset_nrows <- bset %>% base::nrow()
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done bset_nrows={bset_nrows}, elapsed={etime}.{RET}{RET}"))
  
  bset
}

ssetBeta2bset = function(sset, bset, nkey, del='_', 
                         quality.mask=FALSE, nondetection.mask=FALSE,
                         mask.use.tcga=FALSE, pval.threshold=1, sum.TypeI=FALSE,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetBeta2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  datTag  <- 'beta'
  betaTag <- paste(nkey,datTag, sep=del) # %>% rlang::sym()
  
  stime <- system.time({
    dat <- NULL
    
    betas <- sesame::getBetas(sset=sset,
                              quality.mask=quality.mask, 
                              nondetection.mask=nondetection.mask, 
                              mask.use.tcga=mask.use.tcga, 
                              pval.threshold=pval.threshold, 
                              sum.TypeI=sum.TypeI) %>%
      tibble::enframe(name='Probe_ID', value=betaTag)
    
    dat <- add2bset(bset=bset, inf1=betas, keyA='Probe_ID', verbose=verbose,vt=1,tc=0,tt=tt)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

ssetPval2bset = function(sset, bset, nkey, pkey, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetPval2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  datTag  <- 'pval'
  pvalTag <- paste(nkey,pkey,datTag, sep=del) # %>% rlang::sym()
  
  stime <- system.time({
    dat <- NULL
    
    pvals <- sset@pval %>% tibble::enframe(name='Probe_ID', value=pvalTag)
    dat <- add2bset(bset=bset, inf1=pvals, keyA='Probe_ID', verbose=verbose,vt=1,tc=0,tt=tt)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

ssetSigs2bset = function(sset, bset, nkey, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetSig2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    dat <- NULL
    
    datTag <- 'sig'
    grnTag <- paste(nkey,'Grn',datTag, sep=del) %>% rlang::sym()
    redTag <- paste(nkey,'Red',datTag, sep=del) %>% rlang::sym()
    
    IG_tib <- sset@IG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    IR_tib <- sset@IR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    OG_tib <- sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    OR_tib <- sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    I1 <- dplyr::bind_rows(dplyr::inner_join(IG_tib,OR_tib, by=c('Probe_ID', 'Design_Type')),
                           dplyr::inner_join(OG_tib,IR_tib, by=c('Probe_ID', 'Design_Type')) )
    
    I2 <- sset@II %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      dplyr::rename(!!redTag :=U, !!grnTag :=M)
    
    dat <- add2bset(bset=bset, inf1=I1, inf2=I2, keyA='Probe_ID',keyB='Design_Type', 
                    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

add2bset = function(bset, inf1, inf2=NULL, keyA=NULL,keyB=NULL,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'add2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    dat <- NULL
    if (is.null(inf2)) inf2 <- inf1
    if (is.null(keyB)) {
      dat[['I1']] <- bset[[1]] %>% dplyr::inner_join(inf1, by=c(keyA))
    } else {
      dat[['I1']] <- bset[[1]] %>% dplyr::inner_join(inf1, by=c(keyA,keyB))
    }
    dat[['I2']] <- bset[[2]] %>% dplyr::inner_join(inf2, by=c(keyA))
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

ssetToBetaTib_old = function(sset, name, 
                             quality.mask = FALSE, nondetection.mask = FALSE, 
                             correct.switch = TRUE, mask.use.tcga = FALSE, 
                             pval.threshold = 1, 
                             pval.method = "pOOBAH", sum.TypeI = FALSE,
                             as.enframe=FALSE, percision=0,
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToBetaTib_old'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}.{RET}"))
  
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}      quality.mask={quality.mask}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} nondetection.mask={nondetection.mask}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}    correct.switch={correct.switch}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}     mask.use.tcga={mask.use.tcga}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}    pval.threshold={pval.threshold}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}         sum.TypeI={sum.TypeI}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}       pval.method={pval.method}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}        as.enframe={as.enframe}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}         percision={percision}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}              sset={RET}"))
  if (verbose>=vt+4) print(sset)
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{RET}{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    ret_tib <- sesame::getBetas(sset=sset,
                                quality.mask = quality.mask, 
                                nondetection.mask = nondetection.mask, 
                                correct.switch = correct.switch, 
                                mask.use.tcga = mask.use.tcga, 
                                pval.threshold = pval.threshold, 
                                pval.method = pval.method, 
                                sum.TypeI = sum.TypeI)
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
    if (verbose>=vt+4) ret_tib %>% head() %>% print()
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{RET}{RET}"))
    
    if (percision!=0) ret_tib <- round(ret_tib, percision)
    if (as.enframe) ret_tib <- ret_tib %>% 
        tibble::enframe(name='Probe_ID', value=name)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}


#
# TBD:: ssetToPrbTib (Probe_ID,beta,pvals...)
#
# ssetToPrbTib(sset=rdat$raw_sset, verbose=10)
ssetToPrbTib_old = function(sset,
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToPrbTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # Validate p-value slots::
    #
    slot_names <- slotNames(sset)
    if (grep("pvals",slot_names) %>% length() == 0) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR; Failed to find pval in sset!!!{RET}{RET}"))
      return(ret_tib)
    }
    
    # Process each p-value::
    #
    pval_names <- sset@pval %>% names
    for (pval_name in pval_names) {
      # list = ret_tib <- sset@pval[pval_name]
      # dobl = ret_tib <- sset@pval[[pval_name]]
      
      ret_tib <- sset@pval[[pval_name]]
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

sset2calls_old = function(sset, workflow, 
                          quality.mask = FALSE, nondetection.mask = FALSE, 
                          correct.switch = TRUE, mask.use.tcga = FALSE, 
                          pval.threshold = 1, 
                          pval.method = "pOOBAH", sum.TypeI = FALSE,
                          
                          as.enframe=FALSE, percisionBeta=0, percisionPval=0,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sset2calls_old'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}.{RET}"))
  if (verbose>=vt+4) print(sset)
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # noob:: beta
    #
    name <- paste(workflow,'beta', sep='_')
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Mutating/Settting name={name}...{RET}"))
    
    #
    # ssetToBeta (provide return type=tib/dat)
    #
    beta <- ssetToBetaTib(sset=sset, name=name, 
                          as.enframe=as.enframe,
                          percision=percisionBeta, 
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    
    
    ret_tib <- tibble::enframe(beta, name='Probe_ID', value=name)
    if (verbose>=vt+4) head(ret_tib) %>% print()
    
    # PnegEcdf:: negs
    #
    name <- paste(workflow,'negs', sep='_')
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Mutating/Setting name={name}...{RET}"))
    
    ssetA <- mutateSesame(sset=sset, method='detectionPnegEcdf', 
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (verbose>=vt+4) print(ssetA)
    
    pvalA <- ssetToPvalTib(sset=ssetA, method='PnegEcdf', name=name, 
                           percision=percisionPval, 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    # if (verbose>=vt+4) head(pval) %>% print()
    ret_tib <- ret_tib %>% dplyr::left_join(pvalA, by="Probe_ID")
    
    # pOOBAH:: poob
    #
    name <- paste(workflow,'poob', sep='_')
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Mutating/Setting name={name}...{RET}"))
    
    ssetB <- mutateSesame(sset=sset, method='pOOBAH', 
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (verbose>=vt+4) print(ssetB)
    
    pvalB <- ssetToPvalTib(sset=sset, method='pOOBAH', name=name, 
                           percision=percisionPval, 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    # if (verbose>=vt+4) head(pval) %>% print()
    ret_tib <- ret_tib %>% dplyr::left_join(pvalB, by="Probe_ID")
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done ret_cnt={ret_cnt}, elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Extracted Sesame Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

inferSexKaryotypes_Copy = function (sset) 
{
  stopifnot(is(sset, "SigSet"))
  sex.info <- getSexInfo_Copy(sset)
  auto.median <- median(sex.info[paste0("chr", seq_len(22))], 
                        na.rm = TRUE)
  XdivAuto <- sex.info["medianX"]/auto.median
  YdivAuto <- sex.info["medianY"]/auto.median
  if (XdivAuto > 1.2) {
    if (sex.info["fracXlinked"] >= 0.5) 
      sexX <- "XaXi"
    else if (sex.info["fracXmeth"] > sex.info["fracXunmeth"]) 
      sexX <- "XiXi"
    else sexX <- "XaXa"
  }
  else {
    if (sex.info["fracXmeth"] > sex.info["fracXunmeth"]) 
      sexX <- "Xi"
    else sexX <- "Xa"
  }
  if ((sexX == "Xi" || sexX == "Xa") && XdivAuto >= 1 && sex.info["fracXlinked"] >= 
      0.5) 
    sexX <- "XaXi"
  if (YdivAuto > 0.3 || sex.info["medianY"] > 2000) 
    sexY <- "Y"
  else sexY <- ""
  karyotype <- paste0(sexX, sexY)
  karyotype
}

inferSex_Copy = function (sset) 
{
  stopifnot(is(sset, "SigSet"))
  sex.info <- getSexInfo_Copy(sset)[seq_len(3)]
  as.character(predict(sesameDataGet("sex.inference"), sex.info))
}

# 
# sesameData::sesameDataGet(paste0(rdat$raw_sset@platform, '.probeInfo'))$chrY.clean
# sex_info <- getSexInfo_Copy(sset=rdat$raw_sset)
# kar_info <- sesame::inferSexKaryotypes(sset=rdat$raw_sset)
getSexInfo_Copy = function (sset) 
{
  if (is(sset, "SigSetList")) 
    return(do.call(cbind, lapply(sset, getSexInfo)))
  stopifnot(is(sset, "SigSet"))
  cleanY <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$chrY.clean
  # cleanY %>% length() %>% print()
  
  xLinked <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$chrX.xlinked
  # xLinked %>% length() %>% print()
  
  probe2chr <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$probe2chr.hg19
  # print(probe2chr)
  
  xLinkedBeta <- sesame::getBetas(sset=sesame::subsetSignal(sset, xLinked), 
                                  quality.mask = FALSE)
  intens <- sesame::totalIntensities(sset)
  probes <- intersect(names(intens), names(probe2chr))
  intens <- intens[probes]
  probe2chr <- probe2chr[probes]
  # print(probe2chr)
  
  # return( sesame::subsetSignal(sset, cleanY) )
  # return( median(sesame::totalIntensities(sesame::subsetSignal(sset, cleanY))) )
  
  c(medianY = median(sesame::totalIntensities(sesame::subsetSignal(sset, cleanY)), na.rm=TRUE), 
    medianX = median(sesame::totalIntensities(sesame::subsetSignal(sset, xLinked)), na.rm=TRUE), fracXlinked = 
      (sum(xLinkedBeta > 0.3 & xLinkedBeta < 0.7, na.rm = TRUE)/sum(!(is.na(xLinkedBeta))) ), 
    fracXmeth = (sum(xLinkedBeta > 0.7, na.rm = TRUE)/sum(!(is.na(xLinkedBeta)))), 
    fracXunmeth = (sum(xLinkedBeta < 0.3, na.rm = TRUE)/sum(!(is.na(xLinkedBeta)))), 
    tapply(intens, probe2chr, median, na.rm=TRUE))
}

getSexInfo_New = function (sset) {
  if (is(sset, "SigSetList")) 
    return(do.call(cbind, lapply(sset, getSexInfo)))
  stopifnot(is(sset, "SigSet"))
  cleanY <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$chrY.clean
  xLinked <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$chrX.xlinked
  probe2chr <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$probe2chr.hg19
  xLinkedBeta <- getBetas(subsetSignal(sset, xLinked), mask = FALSE)
  intens <- totalIntensities(sset)
  probes <- intersect(names(intens), names(probe2chr))
  intens <- intens[probes]
  probe2chr <- probe2chr[probes]
  c(medianY = median(totalIntensities(subsetSignal(sset, cleanY))), 
    medianX = median(totalIntensities(subsetSignal(sset, 
                                                   xLinked))), fracXlinked = (sum(xLinkedBeta > 0.3 & 
                                                                                    xLinkedBeta < 0.7, na.rm = TRUE)/sum(!(is.na(xLinkedBeta)))), 
    fracXmeth = (sum(xLinkedBeta > 0.7, na.rm = TRUE)/sum(!(is.na(xLinkedBeta)))), 
    fracXunmeth = (sum(xLinkedBeta < 0.3, na.rm = TRUE)/sum(!(is.na(xLinkedBeta)))), 
    tapply(intens, probe2chr, median))
}

#
#
#
#
# Original Auto-Annotation Documentation::
#
#
#

if (FALSE) {
  getSsheetDataTab = function(tib,
                              minOobPval,minOobPerc,
                              minNegPval,minNegPerc,
                              verbose=0,vt=3,tc=1,tt=NULL) {
    funcTag <- 'getSsheetDataTab'
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- system.time({
      
      ret_tib <- suppressMessages( suppressWarnings(
        dplyr::inner_join(
          tib %>% tidyr::gather(Variable, Value),
          tibble::enframe( unlist( sapply(tib, class) ), name="Variable", value="Data_Type" ),
          by="Variable"
        )
      ))
      if (verbose>=vt+3) {
        cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
        ret_tib %>% print()
      }
      
      dec_tib <- NULL
      dec_tib <- suppressMessages( suppressWarnings(
        getSsheetDescTib(tib, 
                         minOobPval=minOobPval,minOobPerc=minOobPerc,
                         minNegPval=minNegPval,minNegPerc=minNegPerc,
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
          gather(Variable, Description)
      ))
      if (verbose>=vt+3) {
        cat(glue::glue("[{funcTag}]:{tabsStr} dec_tib={RET}"))
        dec_tib %>% print()
      }
      
      ret_tib <- ret_tib %>% dplyr::left_join(dec_tib, by="Variable")
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
        ret_tib %>% print(n=base::nrow(ret_tib))
      }
      
      ret_cnt <- ret_tib %>% base::nrow()
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
    
    ret_tib
  }
  
  getSsheetDescTib = function(tib,
                              minOobPval,minOobPerc,
                              minNegPval,minNegPerc,
                              verbose=0,vt=3,tc=1,tt=NULL) {
    funcTag <- 'getSsheetDescTib'
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}minOobPval={minOobPval}, minOobPerc={minOobPerc}...{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}minNegPval={minNegPval}, minNegPerc={minNegPerc}...{RET}"))
    }
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- system.time({
      
      # Determine Number of Workflows::
      #
      max_idx <- tib %>% 
        dplyr::select(dplyr::starts_with('Method_Idx')) %>% 
        tidyr::gather() %>% dplyr::arrange(-value) %>% 
        head(n=1) %>% dplyr::pull(value)
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Max Index={max_idx}.{RET}"))
      
      # Build Description Sheet
      #
      ret_tib <- suppressMessages( suppressWarnings( 
        getSsheetCoreAnnoTib(minOobPval,minOobPerc,
                             minNegPval,minNegPerc,
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) ) )
      if (verbose>=vt+3) {
        cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
        ret_tib %>% print(n=base::nrow(ret_tib))
      }
      
      for (idx in c(0:max_idx)) {
        if (verbose>=vt+2)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Current Index={idx}.{RET}"))
        
        ret_tib <- ret_tib %>% 
          dplyr::bind_cols(suppressMessages( suppressWarnings(
            getSsheetIndexAnnoTib(idx=idx, 
                                  minOobPval=minOobPval,
                                  minNegPval=minNegPval,
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)) ) )
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
          ret_tib %>% print(n=base::nrow(ret_tib))
        }
      }
      if (verbose>=vt+3) {
        cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
        ret_tib %>% print(n=base::nrow(ret_tib))
      }
      
      ret_cnt <- ret_tib %>% base::ncol()
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
    
    ret_tib
  }
  
  getSsheetCoreAnnoTib = function(minOobPval,minOobPerc,
                                  minNegPval,minNegPerc,
                                  verbose=0,vt=3,tc=1,tt=NULL) {
    funcTag <- 'getSsheetCoreAnnoTib'
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- system.time({
      
      ret_tib <- tibble::tibble(
        # Sample Requeue Suggestion::
        #
        Requeue_Flag_pOOBAH = 
          glue::glue("Flag to automatically requeue a faild sample based on percent loci with pOOBAH detection p-value < {minOobPerc}. Percent threshold = {minOobPerc}."),
        Requeue_Flag_PnegEcdf = 
          glue::glue("Flag to automatically requeue a faild sample based on percent loci with PnegEcdf detection p-value < {minNegPval}. Percent threshold = {minNegPerc}."),
        
        # Sample Identification::
        #
        Sentrix_Name = "Unique Sentrix Name: Sentrix_Barcode + Sentrix_Poscode.",
        Sentrix_Barcode = "Sentrix Bar Code (AKA Sentrix_ID).",
        Sentrix_Poscode = "Sentrix Position Code (AKA Sentrix_Position).",
        Sentrix_Row = "Sentrix Row on chip.",
        Sentrix_Col = "Sentrix Column on Chip.",
        Chip_Type = "Idat Chip Type.",
        Chip_Format = "Idat Chip Format (e.g. 8x1, 12x1, 24x1, etc.).",
        Bead_Pool = "Automatically identified bead pool based on address overlap (e.g. HM450, EPIC, etc.).",
        
        # Bead and Manifest Loci Summary::
        #
        Bead_AvgRep_cg = "Average number of sample cg# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
        Bead_AvgRep_ch = "Average number of sample ch# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
        Bead_AvgRep_rs = "Average number of sample rs# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
        
        Bead_Count_cg  = "Total number of sample unique cg# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        Bead_Count_ch  = "Total number of sample unique ch# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        Bead_Count_rs  = "Total number of sample unique rs# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        
        Bead_Total_cg  = "Total number of sample cg# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        Bead_Total_ch  = "Total number of sample ch# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        Bead_Total_rs  = "Total number of sample rs# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        
        Loci_Count_cg  = "Total number of sample cg# loci overlapping manifest based on address presence NOT detection p-value.",
        Loci_Count_ch  = "Total number of sample ch# loci overlapping manifest based on address presence NOT detection p-value.",
        Loci_Count_rs  = "Total number of sample rs# loci overlapping manifest based on address presence NOT detection p-value.",
        
        # OLD Metrics from v.3.46:: Bead and Manifest Loci Summary::
        #
        CG_Manifest_Count = "Total number of sample cg# loci overlapping manifest based on address presence NOT detection p-value.",
        CH_Manifest_Count = "Total number of sample ch# loci overlapping manifest based on address presence NOT detection p-value.",
        RS_Manifest_Count = "Total number of sample rs# loci overlapping manifest based on address presence NOT detection p-value.",
        
        CG_Bead_Count     = "Total number of sample unique cg# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        CG_Bead_Total     = "Total number of sample cg# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        CG_Bead_AvgRep    = "Average number of sample cg# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
        
        CH_Bead_Count     = "Total number of sample unique ch# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        CH_Bead_Total     = "Total number of sample ch# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        CH_Bead_AvgRep    = "Average number of sample ch# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
        
        RS_Bead_Count     = "Total number of sample unique rs# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        RS_Bead_Total     = "Total number of sample rs# bead-types found overlapping manifest based on address presence NOT number of bead types.",
        RS_Bead_AvgRep    = "Average number of sample rs# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
        
        # Parameters Used::
        #
        minNegPval = "Minimum Negative (PnegEcdf) detection p-value threshold used.",
        minOobPval = "Minimum Out-Of-Band (pOOBAH) detection p-value threshold used.",
        minDeltaBeta = "Minimum delta-Beta cutoff used for Sample-Auto-Detection calculations.",
        
        # Manifest, Platform and Version Auto Detection Results::
        #
        detect_manifest = "Platform name-version used. Identified during Auto-Manifest-Detection.",
        detect_platform = "Platform name used. Identified during Auto-Manifest-Detection.",
        detect_version  = "Platform version used. Identified during Auto-Manifest-Detection.",
        detect_sample_cnt   = "Number of Sample loci used during Auto-Manifest-Detection.",
        detect_manifest_cnt = "Number of Manifest loci used during Auto-Manifest-Detection.",
        detect_match_cnt    = "Number of Sample/Manifest overlapping loci used during Auto-Manifest-Detection.",
        detect_rc_per       = "Highest reciprical overlap of Sample/Manifest loci identified during Auto-Manifest-Detection.",
        
        # OLD Platform Detection Results::
        #
        platformUsed = "Platform name used. Identified during Auto-Manifest-Detection.",
        platVersUsed = "Platform version used. Identified during Auto-Manifest-Detection.",
        
        # Auto-Sample-Detection Results::
        #
        AutoSample_Total_Cnt = "Loci overlap count with sample and max auto-detect-sample.",
        AutoSample_R2_Key = "Max auto-detect-sample name by R-Squared.",
        AutoSample_R2_Val = "Max auto-detect-sample R-Squared value.",
        AutoSample_dB_Key = "Max auto-detect-sample name by delta-Beta percent",
        AutoSample_dB_Cnt = glue::glue("Count of loci with delta-Beta < {minDeltaBeta} between sample and max auto-detect-sample."),
        AutoSample_dB_Val = glue::glue("Percent of loci with delta-Beta < {minDeltaBeta} between sample and max auto-detect-sample."),
        
        # Iscan Decode and Extration Dates::
        #
        Iscan_Decoding_ScannerID = "Iscan ID at time of decode.",
        Iscan_Decoding_Year = "Time stamp year of chip decode.",
        Iscan_Decoding_Mon  = "Time stamp month of chip decode.",
        Iscan_Decoding_Day  = "Time stamp day of chip decode.",
        Iscan_Decoding_Hour = "Time stamp hour of chip decode.",
        Iscan_Decoding_Min  = "Time stamp minute of chip decode.",
        Iscan_Decoding_Sec  = "Time stamp second of chip decode",
        
        Iscan_Extract_ScannerID = "Iscan ID at time of chip sample intensity extraction.",
        Iscan_Extract_Year = "Time stamp year of chip sample intensity extraction.",
        Iscan_Extract_Mon  = "Time stamp month of chip sample intensity extraction.",
        Iscan_Extract_Day  = "Time stamp day of chip sample intensity extraction.",
        Iscan_Extract_Hour = "Time stamp hour of chip sample intensity extraction.",
        Iscan_Extract_Min  = "Time stamp minute of chip sample intensity extraction.",
        Iscan_Extract_Sec  = "Time stamp second of chip decode"
      )
      
      ret_cnt <- ret_tib %>% base::ncol()
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
    
    ret_tib
  }
  
  getSsheetIndexAnnoTib = function(idx, 
                                   minOobPval,minNegPval,
                                   verbose=0,vt=3,tc=1,tt=NULL) {
    funcTag <- 'getSsheetIndexAnnoTib'
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- system.time({
      
      ret_tib <- tibble::tibble(
        # Inference and Predictions::
        #
        Method_Key = glue::glue("Method workflow name used for metrics ending with  {idx}."),
        Method_Idx = glue::glue("Method workflow index used for metrics ending with  {idx}."),
        GCT = glue::glue("GCT Score for method index  {idx}. ",
                         "See Bioconductor Package sesame::GCT(sset): ",
                         "Compute GCT score for internal bisulfite conversion control. ",
                         "The function takes a SigSet as input. The higher the GCT score, the more likely the incomplete conversion."),
        
        Sex = glue::glue("Sex call for method index  {idx}. ",
                         "See Bioconductor Package ‘sesame’ sesame::inferSex(sset)."),
        SexKaryotype = glue::glue("Sex Karyotype call for method index  {idx}. ",
                                  "See Bioconductor Package ‘sesame’ sesame::inferSexKaryotypes(sset)."),
        
        Ethnicity = glue::glue("Ethnicity call for method index  {idx}. ",
                               "See Bioconductor Package ‘sesame’ sesame::inferEthnicity(sset).",
                               "This function uses both the built-in rsprobes as well as the type I Color-Channel-Switching probes to infer ethnicity."),
        
        SwapCntIGR = glue::glue("Number of Grn to Red swapped channels for method index  {idx}. ",
                                "See Bioconductor Package ‘sesame’ sesame::inferTypeIChannel(sset).",
                                "Infer and reset color channel for Type-I probes instead of using what is specified in manifest."),
        SwapCntIRG = glue::glue("Number of Red to Grn swapped channels for method index  {idx}. ",
                                "See Bioconductor Package ‘sesame’ sesame::inferTypeIChannel(sset).",
                                "Infer and reset color channel for Type-I probes instead of using what is specified in manifest."),
        
        AgeSkinBlood = glue::glue("Horvath Skin and Blood age predictor call for method index  {idx}. ",
                                  "See Bioconductor Package ‘sesame’ sesame::predictAgeSkinBlood(betas)."),
        
        # Detection P-values:: pOOBAH
        #
        pOOBAH_cg_1_pass_perc = glue::glue("Percent cg Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index  {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
        pOOBAH_cg_2_pass_perc = glue::glue("Percent cg Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index  {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
        
        pOOBAH_ch_1_pass_perc = glue::glue("Percent ch Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index  {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
        pOOBAH_ch_2_pass_perc = glue::glue("Percent ch Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index  {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
        
        pOOBAH_rs_1_pass_perc = glue::glue("Percent rs Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index  {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
        pOOBAH_rs_2_pass_perc = glue::glue("Percent rs Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index  {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
        
        # Detection P-values:: PnegEcdf
        #
        PnegEcdf_cg_1_pass_perc = glue::glue("Percent cg Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index  {idx}. ",
                                             "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
        PnegEcdf_cg_2_pass_perc = glue::glue("Percent cg Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index  {idx}. ",
                                             "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
        
        PnegEcdf_ch_1_pass_perc = glue::glue("Percent ch Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index  {idx}. ",
                                             "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
        PnegEcdf_ch_2_pass_perc = glue::glue("Percent ch Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index  {idx}. ",
                                             "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
        
        PnegEcdf_rs_1_pass_perc = glue::glue("Percent rs Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index  {idx}. ",
                                             "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
        PnegEcdf_rs_2_pass_perc = glue::glue("Percent rs Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index  {idx}. ",
                                             "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
        
        # Average Beta Values::
        #
        cg_1_mean = glue::glue("Mean beta-value of cg Infinium I loci for method  {idx}."),
        cg_2_mean = glue::glue("Mean beta-value of cg Infinium II loci for method  {idx}."),
        
        ch_1_mean = glue::glue("Mean beta-value of ch Infinium I loci for method  {idx}."),
        ch_2_mean = glue::glue("Mean beta-value of ch Infinium II loci for method  {idx}."),
        
        rs_1_mean = glue::glue("Mean beta-value of rs Infinium I loci for method  {idx}."),
        rs_2_mean = glue::glue("Mean beta-value of rs Infinium II loci for method  {idx}."),
        
        # Average Intensities::
        #
        BISULFITE_CONVERSION_I_2_M_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_I Grn loci for method  {idx}."),
        BISULFITE_CONVERSION_I_2_U_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_I Red loci for method  {idx}."),
        
        BISULFITE_CONVERSION_II_2_M_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_II Grn loci for method  {idx}."),
        BISULFITE_CONVERSION_II_2_U_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_II Red loci for method  {idx}."),
        
        cg_2_M_mean = glue::glue("Mean intensity of cg Infinium II Grn (methylated) loci for method  {idx}."),
        cg_2_U_mean = glue::glue("Mean intensity of cg Infinium II Red (unmethylated) loci for method  {idx}."),
        
        cg_IG_M_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Grn (methylated) loci for method  {idx}."),
        cg_IG_U_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Grn (unmethylated) loci for method  {idx}."),
        cg_IR_M_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Red (methylated) loci for method  {idx}."),
        cg_IR_U_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Red (unmethylated) loci for method  {idx}."),
        
        cg_OG_M_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Grn (methylated) loci for method  {idx}."),
        cg_OG_U_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Grn (unmethylated) loci for method  {idx}."),
        cg_OR_M_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Red (methylated) loci for method  {idx}."),
        cg_OR_U_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Red (unmethylated) loci for method  {idx}."),
        
        ch_2_M_mean = glue::glue("Mean intensity of ch Infinium II Grn (methylated) loci for method  {idx}."),
        ch_2_U_mean = glue::glue("Mean intensity of ch Infinium II Red (unmethylated) loci for method  {idx}."),
        
        EXTENSION_2_M_mean = glue::glue("Mean intensity of EXTENSION Grn loci for method  {idx}."),
        EXTENSION_2_U_mean = glue::glue("Mean intensity of EXTENSION Red loci for method  {idx}."),
        
        HYBRIDIZATION_2_M_mean = glue::glue("Mean intensity of HYBRIDIZATION Grn loci for method  {idx}."),
        HYBRIDIZATION_2_U_mean = glue::glue("Mean intensity of HYBRIDIZATION Red loci for method  {idx}."),
        
        NEGATIVE_2_M_mean = glue::glue("Mean intensity of NEGATIVE Grn loci for method  {idx}."),
        NEGATIVE_2_U_mean = glue::glue("Mean intensity of NEGATIVE Red loci for method  {idx}."),
        
        NON_POLYMORPHIC_2_M_mean = glue::glue("Mean intensity of NON_POLYMORPHIC Grn loci for method  {idx}."),
        NON_POLYMORPHIC_2_U_mean = glue::glue("Mean intensity of NON_POLYMORPHIC Red loci for method  {idx}."),
        
        NORM_A_2_M_mean = glue::glue("Mean intensity of NORM_A Grn loci for method  {idx}."),
        NORM_A_2_U_mean = glue::glue("Mean intensity of NORM_A Red loci for method  {idx}."),
        
        NORM_C_2_M_mean = glue::glue("Mean intensity of NORM_C Grn loci for method  {idx}."),
        NORM_C_2_U_mean = glue::glue("Mean intensity of NORM_C Red loci for method  {idx}."),
        
        NORM_G_2_M_mean = glue::glue("Mean intensity of NORM_G Grn loci for method  {idx}."),
        NORM_G_2_U_mean = glue::glue("Mean intensity of NORM_G Red loci for method  {idx}."),
        
        NORM_T_2_M_mean = glue::glue("Mean intensity of NORM_T Grn loci for method  {idx}."),
        NORM_T_2_U_mean = glue::glue("Mean intensity of NORM_T Red loci for method  {idx}."),
        
        RESTORATION_2_M_mean = glue::glue("Mean intensity of RESTORATION Grn loci for method  {idx}."),
        RESTORATION_2_U_mean = glue::glue("Mean intensity of RESTORATION Red loci for method  {idx}."),
        
        rs_2_M_mean = glue::glue("Mean intensity of rs Infinium II Grn (methylated) loci for method  {idx}."),
        rs_2_U_mean = glue::glue("Mean intensity of rs Infinium II Red (unmethylated) loci for method  {idx}."),
        
        rs_IG_M_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Grn (methylated) loci for method  {idx}."),
        rs_IG_U_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Grn (unmethylated) loci for method  {idx}."),
        rs_IR_M_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Red (methylated) loci for method  {idx}."),
        rs_IR_U_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Red (unmethylated) loci for method  {idx}."),
        
        rs_OG_M_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Grn (methylated) loci for method  {idx}."),
        rs_OG_U_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Grn (unmethylated) loci for method  {idx}."),
        rs_OR_M_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Red (methylated) loci for method  {idx}."),
        rs_OR_U_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Red (unmethylated) loci for method  {idx}."),
        
        SPECIFICITY_I_2_M_mean = glue::glue("Mean intensity of SPECIFICITY_I Grn loci for method  {idx}."),
        SPECIFICITY_I_2_U_mean = glue::glue("Mean intensity of SPECIFICITY_I Red loci for method  {idx}."),
        
        SPECIFICITY_II_2_M_mean = glue::glue("Mean intensity of SPECIFICITY_II Grn loci for method  {idx}."),
        SPECIFICITY_II_2_U_mean = glue::glue("Mean intensity of SPECIFICITY_II Red loci for method  {idx}."),
        
        STAINING_2_M_mean = glue::glue("Mean intensity of STAINING Grn loci for method  {idx}."),
        STAINING_2_U_mean = glue::glue("Mean intensity of STAINING Red loci for method  {idx}."),
        
        TARGET_REMOVAL_2_M_mean = glue::glue("Mean intensity of TARGET_REMOVAL Grn loci for method  {idx}."),
        TARGET_REMOVAL_2_U_mean = glue::glue("Mean intensity of TARGET_REMOVAL Red loci for method  {idx}.")
      ) %>% purrr::set_names(paste(names(.),idx, sep='_'))
      
      ret_cnt <- ret_tib %>% base::ncol()
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
    
    ret_tib
  }
  
}

# End of file
