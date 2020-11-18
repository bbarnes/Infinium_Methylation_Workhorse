
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("plyr")) )
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Global Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Single Sample Workflow::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sesamizeSingleSample = function(prefix, man, add, ref, opt, workflows, 
                                retData=FALSE, non_ref=FALSE, 
                                del='_', vt=3,tc=0) {
  funcTag <- 'sesamizeSingleSample'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  verbose <- opt$verbose
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} prefix={prefix}.{RET}"))
  
  # TBD:: Validate all options are present in opt!!!
  
  # Retrun Value::
  ret <- NULL
  
  tTracker <- timeTracker$new()
  stime <- system.time({
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Define Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    basecode <- basename(prefix)
    out_name <- paste(basecode, sep=del)
    out_prefix <- file.path(opt$outDir, out_name)
    
    stampBeg_txt <- file.path(opt$outDir, paste(out_name, 'timestamp.beg.txt', sep=del) )
    stampEnd_txt <- file.path(opt$outDir, paste(out_name, 'timestamp.end.txt', sep=del) )
    stampBeg_cmd <- glue::glue('touch {stampBeg_txt}')
    stampEnd_cmd <- glue::glue('touch {stampEnd_txt}')
    system(stampBeg_cmd)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Extract Raw idat::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
    idat_rds  <- paste(out_prefix, 'idat.rds', sep=del)
    idat_list <- prefixToIdat(prefix=prefix, 
                              load=opt$load_idat, save=opt$save_idat, rds=idat_rds,
                              verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    # return(idat_list)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Get Sesame Manifest/Address Tables::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

    #
    # Need to build idat to manifest map workflow::
    #  - map
    #  - idat + man -> new manifest
    
    man_map_tib <- idatToManifestMap(tib=idat_list$sig, mans=mans,
                                     verbose=verbose,tc=tc+1,tt=tTracker)
    
    top_man_tib <- mans[[man_map_tib[1,]$manifest]] %>% 
      dplyr::distinct(Probe_ID, .keep_all=TRUE)
    
    ses_dat <- idat2manifest(sigs=idat_list$sig, mans=mans,
                             verbose=verbose,tc=tc+1,tt=tTracker)
    
    # Make sure we have unique names::
    ses_dat$man <- ses_dat$man %>% dplyr::distinct(Probe_ID, .keep_all=TRUE)
    ses_man_tib <- ses_dat$man
    
    ses_add_tib <- ses_dat$add
    platform_key <- ses_dat$platform
    version_key  <- ses_dat$manifest
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Initialize Data (Structures) Tables::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Build Sample Sheet::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Initializing Auto-Sample-Sheet: [idat/bead/pheno].{RET}"))
    
    requ_sum_tib <- NULL
    call_sum_tib <- NULL
    sset_sum_tib <- NULL
    phen_sum_tib <- NULL
    pred_sum_tib <- NULL
    
    if (retData) {
      ret$idat <- idat_list
      ret$tmap <- man_map_tib
      ret$tman <- top_man_tib
      ret$sman <- ses_man_tib
      ret$sadd <- ses_add_tib
      # return(ret)
    }
    
    bead_sum_tib <- ses_add_tib %>% dplyr::filter(Probe_Type=='cg') %>% 
      dplyr::select(Address) %>% dplyr::left_join(idat_list$sig, by="Address") %>% 
      dplyr::summarise(CG_Bead_Count=n(), 
                       CG_Bead_Total=sum(Bead_Grn,Bead_Red, na.rm=TRUE), 
                       CG_Bead_AvgRep=round(CG_Bead_Total/CG_Bead_Count/2),1,
                       .groups='drop')
    if (verbose>=vt+5) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} bead_sum_tib={RET}"))
    if (verbose>=vt+5) bead_sum_tib %>% print()
    
    pool_sum_tib <- ses_man_tib %>% dplyr::filter(Probe_Type=='cg' | Probe_Type=='ch' | Probe_Type=='rs') %>% 
      dplyr::mutate(Probe_Type=stringr::str_to_upper(Probe_Type)) %>%
      dplyr::group_by(Probe_Type) %>%
      dplyr::summarise(Count=n(), .groups='drop') %>% 
      tidyr::spread(Probe_Type,Count) %>% 
      purrr::set_names(paste(names(.),'Manifest_Count',sep='_') ) %>% 
      addBeadPoolToSampleSheet(field='CG_Manifest_Count') %>% dplyr::ungroup()
    if (verbose>=vt+5) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pool_sum_tib={RET}"))
    if (verbose>=vt+5) pool_sum_tib %>% print()
    
    chipFormat <- idat_list$ann %>% dplyr::select(Chip_Format) %>% dplyr::pull()
    beadPool   <- pool_sum_tib %>% dplyr::select(Bead_Pool) %>% dplyr::pull()
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} chipFormat={chipFormat}, beadPool={beadPool}.{RET}"))
    
    ssheet_tib <- idat_list$ann
    if (!is.null(pool_sum_tib)) ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(pool_sum_tib)
    if (!is.null(bead_sum_tib)) ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(bead_sum_tib)
    ssheet_tib <- ssheet_tib %>% tibble::add_column(minNegPval   = opt$minNegPval)
    ssheet_tib <- ssheet_tib %>% tibble::add_column(minOobPval   = opt$minOobPval)
    ssheet_tib <- ssheet_tib %>% tibble::add_column(minDeltaBeta = opt$minDeltaBeta)
    ssheet_tib <- ssheet_tib %>% tibble::add_column(platformUsed = platform_key)
    ssheet_tib <- ssheet_tib %>% tibble::add_column(platVersUsed = version_key)
    
    if (retData) {
      ret$bead <- bead_sum_tib
      ret$pool <- pool_sum_tib
      ret$prefix     <- prefix
      ret$platform   <- platform_key
      ret$version    <- version_key
      ret$ssheet_tib <- ssheet_tib
      # return(ret)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Output Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Defining outputs...{RET}"))
    
    if (opt$buildSubDir) opt$outDir <- file.path(opt$outDir, chipFormat, beadPool)
    if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
    
    basecode   <- idat_list$ann$Sentrix_Name[1]
    out_name   <- paste(basecode, platform_key, version_key, sep=del)
    out_prefix <- file.path(opt$outDir, out_name)
    
    ssheet_csv <- paste(out_prefix, 'AutoSampleSheet.csv.gz', sep=del)
    # calls_csv  <- paste(out_prefix, 'calls.csv.gz', sep=del)
    calls_csv  <- paste(out_prefix, 'calls.csv', sep=del)
    times_csv  <- paste(out_prefix, 'run-times.csv.gz', sep=del)
    
    raw_sset_csv <- paste(out_prefix, 'raw-sset.csv.gz', sep=del)
    raw_sset_rds <- paste(out_prefix, 'raw-sset.rds', sep=del)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Defined output files: out_prefix={out_prefix}, basecode={basecode}, out_name={out_name}.{RET}"))
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Defined output files: platform_key={platform_key}, version_key={version_key}, outDir={opt$outDir}.{RET}"))
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}  raw_sset_csv={raw_sset_csv}.{RET}") )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}  raw_sset_rds={raw_sset_rds}.{RET}") )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}    ssheet_csv={ssheet_csv}.{RET}") )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}     calls_csv={calls_csv}.{RET}") )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}     times_csv={times_csv}.{RET}{RET}{RET}") )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Initialize RAW SSET::
    #
    #  - Required for initial phenotype calcs before channel inference
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # Current goal:: build out FIXED raw_sset + workflow -> calls file + ssheet update
    #
    
    raw_sset <- NULL
    raw_sset <- newSset(prefix=prefix, 
                        platform=platform_key, manifest=ses_man_tib,
                        load=opt$load_sset,save=opt$save_sset,rds=sset_rds,
                        verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)

    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} raw_sset={RET}") )
    if (verbose>=vt+4) print(raw_sset)
    
    ww <- 0
    cur_workflow <- 'raw'
    cur_sset <- raw_sset
    
    cur_sset_rds <- NULL
    cur_sigs_csv <- NULL
    cur_ssum_csv <- NULL
    cur_call_csv <- NULL
    
    ret_dat_list <- NULL
    ret_dat_list[['raw']] <- ssetToSummary(
      sset=cur_sset, man=ses_man_tib, idx=ww, workflow=cur_workflow,
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
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Built RAW Data.{RET}{RET}"))
    
    if (retData) {
      ret$rsum_list <- ret_dat_list
      # return(ret)
    }
    if (verbose>=vt) 
      cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 SSET to Calls by Order of Operations:: workflows
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting Workflows...{RET}{RET}"))
    
    auto_beta_key <- NULL
    auto_negs_key <- NULL
    workflow_cnt <- length(workflows)
    for (ww in seq(1,workflow_cnt)) {
      cur_workflow <- workflows[ww]
      cur_sset_rds <- paste(out_prefix, paste0(cur_workflow,'-sset.rds'), sep=del)
      
      if (is.null(auto_beta_key)) auto_beta_key <- paste(cur_workflow,'beta', sep=del)
      if (is.null(auto_negs_key)) auto_negs_key <- paste(cur_workflow,'PnegEcdf', sep=del)
      
      cur_sset <- NULL
      if (opt$load_sset && file.exists(cur_sset_rds)) {
        if (verbose>=vt+1) 
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading RDS={cur_sset_rds}.{RET}"))
        cur_sset <- readr::read_rds(cur_sset_rds)
      } else {
        stopifnot(!is.null(raw_sset))
        cur_sset <- mutateSSET_workflow(
          sset=raw_sset, workflow=cur_workflow, 
          save=opt$save_sset, rds=cur_sset_rds,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      }
      stopifnot(!is.null(cur_sset))
      
      cur_sset_rds <- NULL
      cur_sigs_csv <- NULL
      cur_ssum_csv <- NULL
      cur_call_csv <- NULL
      
      ret_dat_list[[cur_workflow]] <- ssetToSummary(
        sset=cur_sset, man=ses_man_tib, idx=ww, workflow=cur_workflow,
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
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_dat_list={RET}"))
        print(ret_dat_list[[cur_workflow]])
      }

      if (verbose>=vt) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. cur_workflow={cur_workflow}.{RET}{RET}"))
        cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
      }
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
    
    all_call_tib <- NULL
    all_call_tib <- ses_man_tib %>% dplyr::select(Probe_ID)
    for (work_name in names(ret_dat_list)) {
      cur_list  <- ret_dat_list[[work_name]]
      cur_data  <- cur_list$call_dat
      cur_sheet <- cur_list$sam_sheet

      if (verbose>=vt+5) {
        cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}"))
        cat(glue::glue("[{funcTag}]:{tabsStr} Cur Join; ret_dat_list[{work_name}]={RET}"))
        cur_list %>% print()
      }
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

      if (is.null(all_call_tib)) {
        all_call_tib <- cur_data
      } else {
        all_call_tib <- all_call_tib %>% 
        #  dplyr::left_join(cur_data, by="Probe_ID")
        #  dplyr::full_join(cur_data, by="Probe_ID")
          dplyr::inner_join(cur_data, by="Probe_ID")
      }
      ssheet_tib <- ssheet_tib %>% 
        dplyr::bind_cols(cur_sheet)
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} Done; work_name={work_name}{RET}"))
        cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
      }
    }
    all_call_tib <- all_call_tib %>% dplyr::arrange(Probe_ID)
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Full Join; all_call_tib={RET}"))
      all_call_tib %>% print()
      cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
      
      cat(glue::glue("[{funcTag}]:{tabsStr} Full Join; ssheet_tib={RET}"))
      ssheet_tib %>% print()
      
      cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    }
          
    
    if (retData) {
      ret$rsum_list <- ret_dat_list
      ret$all_calls <- all_call_tib
      ret$ssheet_tib <- ssheet_tib
      # return(ret)
    }

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Add Requeue Flags::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Requeue Flag",
                     "(minOobPerc={opt$minOobPerc},minNegPerc={opt$minNegPerc})...{RET}"))
    
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
    
    ssheet_tib <- ssheet_tib %>% 
      dplyr::select(Requeue_Flag_pOOBAH,Requeue_Flag_PnegEcdf, dplyr::everything())
    
    if (verbose>=vt+4) {
      cat(glue::glue("{RET}[{funcTag}]:{tabsStr} ssheet_tib(all requeue added)={RET}"))
      ssheet_tib %>% print()
      cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    }
    
    # ssheet_tib <- ssheet_tib %>% dplyr::mutate(
    #   Requeue_Flag_pOOBAH=dplyr::case_when(
    #     poob1_val < opt$minOobPerc | poob2_val < opt$minOobPerc ~ TRUE, TRUE ~ FALSE),
    #   Requeue_Flag_PnegEcdf=dplyr::case_when(
    #     negs1_val < opt$minNegPerc | negs2_val < opt$minNegPerc ~ TRUE, TRUE ~ FALSE)
    # ) %>% dplyr::select(Requeue_Flag_pOOBAH,Requeue_Flag_PnegEcdf, dplyr::everything())
    
    # ssheet_tib <- ssheet_tib %>% dplyr::mutate(
    #   Requeue_Flag_pOOBAH=dplyr::case_when(
    #     pOOBAH_cg_1_pass_perc_0 < opt$minOobPerc | pOOBAH_cg_2_pass_perc_0 < opt$minOobPerc ~ TRUE, TRUE ~ FALSE),
    #   Requeue_Flag_PnegEcdf=dplyr::case_when(
    #     PnegEcdf_cg_1_pass_perc_0 < opt$minNegPerc | PnegEcdf_cg_2_pass_perc_0 < opt$minNegPerc ~ TRUE, TRUE ~ FALSE)
    # ) %>% dplyr::select(Requeue_Flag_pOOBAH,Requeue_Flag_PnegEcdf, dplyr::everything())
    
    if (retData) {
      ret$ssheet_tib <- ssheet_tib
      # return(ret)
    }
    
    #
    # Add formatVCF
    # Sum formatVCF
    #
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Auto-Detect Sample::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (opt$auto_detect) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Detecting Auto-Sample-Sheet: [SampleName, rsquared,deltaBeta].{RET}"))
      
      auto_sheet <- autoDetect_Wrapper(
        can=all_call_tib, ref=ref, man=ses_man_tib,
        minPval=opt$minNegPval, minDelta=opt$minDeltaBeta,
        dname='Design_Type', pname='Probe_Type', ptype='cg',
        jval='Probe_ID', field=auto_beta_key, pval=auto_negs_key, suffix='beta', del=del,
        outDir=opt$outDir, sname=out_name, plotMatrix=opt$plotAuto, writeMatrix=opt$write_auto,
        dpi=opt$dpi, format=opt$plotFormat, datIdx=4, non.ref=non_ref,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      
      ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(auto_sheet)

      ssheet_ncols <- ssheet_tib %>% base::ncol()
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Binding Sample Sheet (with auto-detect-sample) ssheet_ncols={ssheet_ncols}.{RET}{RET}"))
      
      if (retData) {
        ret$ssheet_tib <- ssheet_tib
        # return(ret)
      }
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing Outputs: [call,sigs,signalSummary,sampleSheet,time].{RET}"))
    
    # Arrange the Sample Sheet Columns
    ssheet_tib <- ssheet_tib %>% 
      dplyr::select(-starts_with('Iscan_'), starts_with('Iscan_')) %>% 
      dplyr::mutate_if(is.numeric, list(round), 4)
    
    # Format Time Table
    time_tib <- tTracker$time %>% dplyr::mutate_if(base::is.numeric, round, 4)
    
    t_cmd <- glue::glue("touch {calls_csv}")
    cat(glue::glue("Running; t_cmd={t_cmd}!{RET}{RET}"))
    # system(t_cmd)
    
    all_call_tib %>% print()
    all_call_tib <- all_call_tib %>% 
      dplyr::mutate(p_len=stringr::str_length(Probe_ID)) %>%
      dplyr::filter(p_len>1) %>% dplyr::arrange(-p_len)
    
    cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    cat(glue::glue("Beg; Writing; calls_csv={calls_csv}{RET}"))
    all_call_tib %>% print()
    # readr::write_csv(all_call_tib %>% head(), calls_csv)
    readr::write_csv(all_call_tib, calls_csv)
    cat(glue::glue("End; Writing; calls_csv={calls_csv}{RET}{RET}{RET}"))
    cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    
    readr::write_csv(ssheet_tib, ssheet_csv)
    readr::write_csv(time_tib, times_csv)
    
    cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} calls_csv={calls_csv}={RET}"))
    all_call_tib %>% print()
    cat(glue::glue("[{funcTag}]:{RET}{RET}"))
    cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))

    system(stampEnd_cmd)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (verbose>=vt+3) tTracker %>% print()
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. elasped={etime}.{RET}{RET}"))
  
  if (retData) {
    ret$ssheet_tib <- ssheet_tib
    ret$all_calls  <- all_call_tib
    ret$time       <- tTracker
    return(ret)
  }

  ssheet_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Call (Beta/Pval) Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

maskTibs = function(tib, id, field, pval, minPval, del='_', 
                    rm.suffix=FALSE, only.field=FALSE, no.mask=FALSE,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'maskTibs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; field={field}, pval={pval}, minPval={minPval}.{RET}"))
  if (!no.mask) {
    pval_del <- paste0(del,pval)
    snames <- tib %>% dplyr::select(ends_with(pval_del)) %>% names() %>% stringr::str_remove(pval_del)
    for (sname in snames) {
      fstr <- paste(sname,field, sep=del)
      pstr <- paste(sname,pval, sep=del)
      tib <- maskTib(tib, field=fstr, id=id, pval=pstr, minPval=minPval, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
  }
  if (only.field) tib <- tib %>% dplyr::select(1,ends_with(field))
  if (rm.suffix) tib  <- tib %>% purrr::set_names(stringr::str_remove(names(.), paste0(del,field) ) )
  
  tib
}

maskTib = function(tib, id, field, pval, minPval, rm.pval=FALSE, only.field=FALSE, mask.controls=FALSE,
                   verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'maskTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; id={id}, field={field}, pval={pval}, minPval={minPval}.{RET}"))
  
  field <- field %>% rlang::sym()
  pval  <- pval %>% rlang::sym()
  id    <- id %>% rlang::sym()
  # minPval <- minPval %>% rlang::sym()
  
  stime <- system.time({
    if (mask.controls) {
      tib <- tib %>% dplyr::mutate(!!field := case_when(!!pval > !!minPval ~ NA_real_, TRUE ~ !!field))
    } else {
      tib <- tib %>% dplyr::mutate(!!field := case_when(
        !stringr::str_starts(!!id, 'ctl_') & !!pval > !!minPval ~ NA_real_, 
        TRUE ~ !!field))
    }
    if (rm.pval) tib <- tib %>% dplyr::select(-!!pval)
    if (only.field) tib <- tib %>% dplyr::select(!!field)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

maskCall = function(tib, field, minKey, minVal, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'maskCall'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Masking {field} with {minKey} at {minVal}.{RET}"))
  
  field <- field %>% rlang::sym()  
  minKey  <- minKey %>% rlang::sym()
  # minVal <- minVal %>% rlang::sym()
  
  stime <- system.time({
    tib <- tib %>% dplyr::mutate(!!field := case_when(!!minKey >= !!minVal ~ NA_real_, TRUE ~ !!field))
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

loadCallFile = function(file, selKey, datKey=NULL, minKey=NULL, minVal=NULL, prefix=NULL, 
                        retMin=FALSE, retRaw=FALSE, del='_',
                        verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadCallFile'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) {
    if (!is.null(datKey)) {
      cat(glue::glue("[{funcTag}]:{tabsStr} selKey={selKey}, datKey={datKey} file={file}.{RET}")) 
    } else { 
      cat(glue::glue("[{funcTag}]:{tabsStr} selKey={selKey}, datKey=LOAD_ALL, file={file}.{RET}"))
    }
  }
  
  stime <- system.time({
    selKey <- selKey %>% rlang::sym()
    if (!is.null(datKey)) datKey <- datKey %>% rlang::sym()
    if (!is.null(minKey)) minKey <- minKey %>% rlang::sym()
    
    if (stringr::str_detect(file,'.rds$')) {
      tib <- suppressMessages(suppressWarnings(readr::read_rds(file)) )
    } else {
      tib <- suppressMessages(suppressWarnings(readr::read_csv(file)) )
    }
    if (verbose>=vt+4) tib %>% print()
    
    if (retRaw && !is.null(datKey)) {
      if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} IS_RAW retRaw={retRaw}, datKey={datKey}.{RET}"))
      tib <- tib %>% dplyr::select(!!selKey,ends_with(paste0(del,datKey)) )
      if (verbose>=vt+4) tib %>% print()
    } else {
      if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} NOT_RAW retRaw={retRaw}, datKey={datKey}.{RET}"))
      
      if (!is.null(minVal) && !is.null(datKey)) {
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} USE_PVAL minVal={minVal}, SELECT=[selKey={selKey}, minKey={minKey}, datKey={datKey}].{RET}"))
        tib <- tib %>% dplyr::select(!!selKey,!!minKey, !!datKey)
        if (verbose>=vt+4) tib %>% print()
        
        tot_cnt  <- tib %>% base::nrow()
        pre_cnt  <- tib %>% dplyr::filter(is.na(!!datKey)) %>% base::nrow()
        
        tib <- maskTib(tib=tib, id=selKey, field=datKey, pval=minKey, minPval=minVal, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        # tib <- maskCall(tib=tib, field=datKey, minKey=minKey, minVal=minVal, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        
        pos_cnt  <- tib %>% dplyr::filter(is.na(!!datKey)) %>% base::nrow()
        pos_per  <- round(100*pos_cnt/tot_cnt,3)
        
        if (retMin) tib <- tib %>% dplyr::select(!!selKey,!!datKey)
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Masked ({pos_per}%) Beta Values {pre_cnt} -> {pos_cnt}/{tot_cnt}.{RET}"))
      } else {
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} NO_PVAL.{RET}"))
        if (!is.null(datKey)) tib <- tib %>% dplyr::select(!!selKey,!!datKey)
        if (verbose>=vt+4) tib %>% print()
      }
      if (!is.null(prefix)) tib <- tib %>% purrr::set_names(paste(prefix,names(tib), sep='.') ) %>% dplyr::rename(!!selKey := 1)
    }
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

# TBD:: Add option to include Sentrix_Name
loadCallFiles = function(files, selKey, datKey=NULL, minKey=NULL, minVal=NULL, 
                         prefix=NULL, addRep=TRUE, retMin=FALSE, retRaw=FALSE, 
                         ss=NULL, addSentrix=FALSE, max=NULL, del='_', 
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadCallFiles'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  stopifnot(is.vector(files))
  
  join_tibs <- NULL
  files_cnt <- length(files)
  if (verbose>=vt) {
    if (!is.null(datKey)) cat(glue::glue("[{funcTag}]:{tabsStr} Total number of call files={files_cnt}, selKey={selKey}, datKey={datKey}.{RET}"))
    else cat(glue::glue("[{funcTag}]:{tabsStr} Total number of call files={files_cnt}, selKey={selKey}, datKey=ALL_DATA.{RET}"))
  }
  
  stime <- system.time({
    for (ii in c(1:files_cnt)) {
      if (verbose>=vt+2) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ii={ii}, file={files[ii]}.{RET}"))
      tib <- loadCallFile(file=files[ii], selKey=selKey, datKey=datKey, minKey=minKey, minVal=minVal,
                          retMin=retMin, retRaw=retRaw,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      rep_str <- paste0('Rep',ii)
      if (addSentrix && !is.null(ss) && !is.null(ss$Sentrix_Name[ii])) rep_str <- paste(rep_str,ss$Sentrix_Name[ii], sep=del)
      if (!is.null(prefix)) rep_str <- paste(prefix,rep_str, sep=del)
      if (addRep) tib <- tib %>% purrr::set_names(paste(rep_str,names(tib), sep='.') ) %>% dplyr::rename(!!selKey := 1)
      
      if (is.null(join_tibs)) { join_tibs <- tib
      } else { join_tibs <- join_tibs %>% dplyr::full_join(tib, by=selKey) }
      
      if (!is.null(max) && ii >= max) break
    }
  })
  nrows <- join_tibs %>% base::nrow()
  ncols <- join_tibs %>% base::ncol()
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} nrows={nrows}, ncols={ncols}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  join_tibs
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Auto-detect Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

autoDetect_Wrapper = function(can, ref, man,
                              
                              minPval, minDelta,
                              dname='Design_Type', pname='Probe_Type', ptype='cg',
                              jval='Probe_ID', field='ind_beta', pval='ind_PnegEcdf', 
                              suffix='beta', del='_',
                              outDir=NULL, sname=NULL, plotMatrix=FALSE, writeMatrix=FALSE,
                              dpi=120, format="png", datIdx=4, non.ref=FALSE,
                              
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'autoDetect_Wrapper'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    auto_can <- man %>% dplyr::rename(Design_Type=DESIGN) %>% 
      dplyr::select(Probe_ID, Design_Type, Probe_Type) %>%
      dplyr::left_join(can, by="Probe_ID")
    # print(auto_can)
    
    auto_ref <- ref
    if (grep("_beta", names(ref)) %>% length() == 0) {
      auto_ref <- ref %>% purrr::set_names(paste(names(.),'beta', sep=del) )
      names(auto_ref)[1] <- 'Probe_ID'
    }
    # print(auto_ref)
    
    auto_data <- sampleDetect(
      can=auto_can, ref=auto_ref, minPval=minPval, minDelta=minDelta,
      dname=dname, pname=pname, ptype=ptype,
      jval=jval, field=field, pval=pval, suffix=suffix, del=del,
      outDir=outDir, sname=sname, plotMatrix=plotMatrix, writeMatrix=writeMatrix,
      dpi=dpi, format=format, datIdx=datIdx, non.ref=non.ref,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  # ret_tib$can  <- auto_can
  # ret_tib$ref  <- auto_ref
  # ret_tib$data <- auto_data
  
  auto_data
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Auto Sample Detection Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
deltaMatrix = function(mat, beg=NULL, end=NULL, minDelta,
                       verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'deltaMatrix'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    ncols  <- base::ncol(mat)
    nrows  <- base::nrow(mat)
    cnames <- colnames(mat)
    
    stopifnot(nrows>0)
    
    del_mat <- matrix(0, nrow=ncols ,ncol=ncols)
    colnames(del_mat) <- cnames
    rownames(del_mat) <- cnames
    
    if (is.null(beg)) beg <- 1
    if (is.null(end)) end <- ncols
    
    max_cnt <- 0
    max_idx <- 0
    max_per <- 0
    max_key <- 'Unknown'
    for (ii in c(beg:end)) {
      for (jj in c(beg:end)) {
        if (ii<jj) {
          cur_cnt <- length(which(abs(rowDiffs(mat[,c(ii,jj)])) <= minDelta))
          del_mat[ii,jj] <- cur_cnt
          del_mat[jj,ii] <- cur_cnt
        } else if (ii==jj) {
          del_mat[ii,jj] <- length(mat[,ii])
        }
      }
    }
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  del_mat
}

rsquaredMatrix = function(mat, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'rsquaredMatrix'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  stime <- system.time({
    cor_val <- NULL
    
    cor_val = tryCatch({
      cor(mat, method="pearson", use="pairwise.complete.obs")
    }, warning = function(w) {
      'warning-SBM-Cor'
    }, error = function(e) {
      'error-SBM-Cor'
    }, finally = {
      'cleanup-SBM-Cor'
    })
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  cor_val
}

sampleDetect = function(can, ref, minPval, minDelta, dname, pname, ptype=NULL,
                        jval, pval, field, suffix, del='_',
                        outDir=NULL, sname=NULL, plotMatrix=FALSE, writeMatrix=FALSE,
                        dpi=120, format='png', datIdx=4, roundData=TRUE, non.ref=TRUE,
                        verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'sampleDetect'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting. dname={dname}, pname={pname},",
                                  "field={field}, jval={jval}, pval={pval}.{RET}"))
  
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr} CAN::{RET}") )
    head(can) %>% print()
    cat(glue::glue("[{funcTag}]:{tabsStr} REF::{RET}") )
    head(ref) %>% print()
  }
  
  stime <- system.time({
    sss <- NULL
    can_names <- c(jval,dname,pname, 'Sample_pval','Sample_beta')
    
    # String Values::
    jval_char <- jval
    
    dname <- dname %>% rlang::sym()
    pname <- pname %>% rlang::sym()
    
    # Sym jval to see if it fixes the warning...
    field <- field %>% rlang::sym()
    jval  <- jval %>% rlang::sym()
    pval  <- pval %>% rlang::sym()
    
    if (!is.null(outDir) && !dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
    
    # 1. Format and Filter Data by Probe Type
    if (!is.null(ptype) && !is.null(pname)) {
      can <- can %>% dplyr::filter(!!pname==ptype)
    }
    
    # 2. Format Candidate Data
    can_nrows <- base::nrow(can)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} can_nrows={can_nrows}.{RET}"))
    
    if (verbose>=vt+6) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Pre-join. dname={dname}, pname={pname}, pval={pval}, field={field}, jval={jval}.{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ref=={RET}"))
      head(ref) %>% print()
      
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tib1=={RET}"))
      can %>% dplyr::select(!!jval,!!dname,!!pname,!!pval,!!field) %>% print()
      
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tib2=={RET}"))
      can %>% dplyr::select(!!jval,!!dname,!!pname,!!pval,!!field) %>% purrr::set_names(can_names) %>% print()
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Verbose Message Complete.{RET}{RET}"))
    }
    tib <- can %>% dplyr::select(!!jval,!!dname,!!pname,!!pval,!!field) %>%
      purrr::set_names(can_names) %>%
      dplyr::left_join(ref, by=jval_char)
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joined.{RET}"))
    
    # 3. Filter on P-value
    tib <- maskTibs(tib, id=jval, field=field, pval=pval, minPval=minPval, del=del,
                    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # 4. Build matrix
    mat <- tib %>% dplyr::select(-c(!!jval,!!dname,!!pname)) %>%
      dplyr::select(ends_with(paste0(del,suffix)) ) %>% as.matrix()
    
    # 5. Plot Matrix
    if (plotMatrix && !is.null(outDir) && !is.null(sname))
      gg <- plotPairs(tib=tib, sample=sname, nameA='Sample', nameB='CanonicalReference',
                      field='beta', field_str='beta', detp='pval', outDir=outDir, minPval=minPval, minDelta=minDelta,
                      dpi=dpi, format=format, datIdx=datIdx, non.ref=non.ref,
                      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # 5. Build R-Squared Matrix
    r2m <- rsquaredMatrix(mat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    r2m_tib <- r2m %>% tibble::as_tibble(rownames='Sample')
    
    # 6. Build Deleta Matrix
    dbm <- deltaMatrix(mat, minDelta=minDelta, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    dbm_tib <- dbm %>% tibble::as_tibble(rownames='Sample')
    
    # cat("Building sss r2m/dbm:\n")
    sss <- dplyr::bind_cols(
      r2m_tib %>% head(n=1) %>% dplyr::select(-Sample) %>% 
        tidyr::gather(key='Sample', value='r2') %>% 
        dplyr::mutate(Sample=stringr::str_remove(Sample, paste0('_',suffix)) ) %>%
        dplyr::filter(Sample!='Sample') %>% 
        dplyr::arrange(-r2) %>% head(n=1) %>%
        purrr::set_names('AutoSample_R2_Key', 'AutoSample_R2_Val'),
      
      dbm_tib %>% head(n=1) %>% dplyr::select(-Sample) %>%
        tidyr::gather(key='Sample', value='passCount') %>%
        dplyr::mutate(Sample=stringr::str_remove(Sample, paste0('_',suffix)),
                      failCount=can_nrows-passCount,
                      passPerc=round(100*passCount/can_nrows,1) ) %>%
        dplyr::filter(Sample!='Sample') %>% 
        dplyr::arrange(-passPerc) %>% head(n=1) %>%
        purrr::set_names('AutoSample_dB_Key', 'AutoSample_Total_Cnt', 'AutoSample_dB_Cnt', 'AutoSample_dB_Val')
    ) %>% dplyr::select(AutoSample_Total_Cnt, AutoSample_R2_Key, AutoSample_R2_Val, 
                        AutoSample_dB_Key, AutoSample_dB_Cnt, AutoSample_dB_Val)
    
    # Write matricies::
    if (writeMatrix && !is.null(outDir) && !is.null(sname)) {
      if (is.null(ptype)) ptype <- 'ALL'
      fname <- paste(sname,'Sample','VS','CanonicalReference',ptype,field,
                     paste0('pval-',minPval),paste0('delta-',minDelta), sep='_')
      r2m_csv <- file.path(outDir, paste0(fname,'.rsquaredMatrix.csv.gz'))
      dbm_csv <- file.path(outDir, paste0(fname,'.deltaMatrix.csv.gz'))
      ord_csv <- file.path(outDir, paste0(fname,'.orderedAutoMatrix.csv.gz'))
      
      if (roundData) r2m_tib <- r2m_tib %>% dplyr::mutate_if(base::is.numeric, round, 6)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing RSquared(CSV)={r2m_csv}.{RET}"))
      readr::write_csv(r2m_tib, r2m_csv)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing DeltaBeta(CSV)={dbm_csv}.{RET}"))
      readr::write_csv(dbm_tib, dbm_csv)
      
      # 7. Build r2/delta Table
      ord_tib <- dplyr::inner_join(
        tibble::enframe(r2m[,1], value='r2'), 
        tibble::enframe(dbm[,1], value='delta_PassCount'), by="name") %>%
        dplyr::mutate(delta_PassPerc=round(100*delta_PassCount/can_nrows,1),
                      delta_FailCount=can_nrows-delta_PassCount,
                      name=stringr::str_remove(name, paste0('_',suffix) )
        ) %>%
        dplyr::filter(name!='Sample') %>% dplyr::select(-delta_PassCount)
      
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing OrderAuto(CSV)={ord_csv}.{RET}"))
      readr::write_csv(ord_tib, ord_csv)
    }
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  sss
}

# End of file
