
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

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Single Sample Workflow::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sesamizeSingleSample = function(prefix, man, add, ref, opts, defs=NULL,
                                mask=NULL,
                                pvals=NULL, min_pvals=NULL, min_percs=NULL,
                                workflows=NULL, 
                                
                                btypes=c('cg','ch','rs'),
                                report_vec=c("Requeue","pass_perc","mean"),
                                
                                retData=FALSE, retData2=FALSE, non_ref=FALSE, 
                                del='_', 
                                verbose=0, vt=3,tc=0) {
  funcTag <- 'sesamizeSingleSample'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; prefix={prefix}...{RET}{RET}"))
  
  # TBD:: Validate all options are present in opt!!!
  opt_tib <- dplyr::bind_rows(opts) %>% tidyr::gather("Option", "Value")
  def_tib <- dplyr::bind_rows(defs) %>% tidyr::gather("Option", "Value")
  
  # Retrun Value::
  ret_cnt <- 0
  ret_dat <- NULL
  stime <- system.time({
    
    tTracker <- timeTracker$new()
    if (!dir.exists(opts$outDir)) dir.create(opts$outDir, recursive=TRUE)
    if (opt$fresh) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Cleaning; prefix={prefix}.{RET}{RET}"))
      cleaned_files <- lapply(list.files(opt$outDir, pattern=prefix, full.names=TRUE), unlink)
    }

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Initialize Starting Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    basecode <- basename(prefix)
    out_name <- paste(basecode, sep=del)
    out_prefix <- file.path(opts$outDir, out_name)
    
    stampBeg_txt <- file.path(opts$outDir, paste(out_name, 'timestamp.beg.txt', sep=del) )
    stampEnd_txt <- file.path(opts$outDir, paste(out_name, 'timestamp.end.txt', sep=del) )
    if (!file.exists(stampBeg_txt))
      readr::write_lines(x=date(),file=stampBeg_txt,sep='\n',append=FALSE)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      Write Run Options/Defaults::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    opts_csv <- file.path(opt$outDir, paste(out_name,'program-options.csv', sep=del) )
    opts_tib <- dplyr::bind_rows(opts) %>% tidyr::gather("Option", "Value")
    if (!is.null(defs)) {
      opts_tib <- opts_tib %>% dplyr::left_join(
        dplyr::bind_rows(defs) %>% tidyr::gather("Option", "Value"),
        by="Option") %>% 
        purrr::set_names(c("Option","Value","Default"))
    }
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Writing opts_csv={opts_csv}...{RET}"))
    readr::write_csv(opts_tib, opts_csv)
    
    opt_ssh_tib <- tibble::tibble(
      Min_DeltaBeta = opts$minDeltaBeta,
      Sesame_Version=packageVersion("sesame"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Extract Raw idat::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    idat_sigs_csv <- paste(out_prefix, 'idat.sigs.csv.gz', sep=del)
    idat_info_csv <- paste(out_prefix, 'idat.info.csv.gz', sep=del)
    idat_sigs_tib <- 
      prefixToIdat(prefix=prefix, load=opts$load_idat, save=opts$save_idat, 
                   csv=idat_sigs_csv, ssh=idat_info_csv,
                   verbose=verbose,vt=vt+4,tc=tc+1,tt=tTracker)
    idat_info_tib <- 
      suppressMessages(suppressWarnings( readr::read_csv(idat_info_csv) ))
    
    if (!is.null(opts$runName))
      idat_info_tib <- idat_info_tib %>% dplyr::mutate(Run_Name=opts$runName)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Get Sesame Manifest/Address Tables::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    man_map_tib <- idatToManifestMap(tib=idat_sigs_tib, mans=mans, sortMax=TRUE,
                                     verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    
    top_man_tib <- mans[[man_map_tib[1,]$manifest]] %>% 
      dplyr::distinct(Probe_ID, .keep_all=TRUE)
    
    bead_sum_tib <- 
      getManifestBeadStats(dat=idat_sigs_tib, man=top_man_tib, types=btypes, 
                           verbose=verbose,vt=vt+3,tc=tc+1,tt=tTracker)
    
    #
    # TBD:: Replace or improve addBeadPoolToSampleSheet()
    #
    bead_ssh_tib <- bead_sum_tib %>% 
      tidyr::spread(name, value) %>%
      addBeadPoolToSampleSheet(field="Loci_Count_cg",
                               verbose=verbose,vt=vt+4,tc=tc+1,tt=tTracker)
    
    if (retData) {
      ret_dat$prefix <- prefix
      
      ret_dat$ossh <- opt_ssh_tib
      
      ret_dat$isig <- idat_sigs_tib
      ret_dat$issh <- idat_info_tib
      
      ret_dat$mmap <- man_map_tib
      ret_dat$sman <- top_man_tib
      
      ret_dat$bsum <- bead_sum_tib
      ret_dat$bssh <- bead_ssh_tib
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Build Basic Detection P-value Calling::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    open_sset_dat <- NULL
    open_beta_tib <- NULL
    open_sum1_ssh <- NULL
    open_sum2_ssh <- NULL
    if (man_map_tib$platform[1]=='EPIC') {
      open_sset_dat <- sesame::openSesame(x=prefix, what = 'sigset')
      open_beta_tib <- sesame::getBetas(sset = open_sset_dat, 
                                        mask=FALSE, sum.TypeI=FALSE)

      open_sum1_ssh <- 
        ssetToPassPercSsheet(sset=open_sset_dat, man=NULL,
                             min=min_pvals[1], per=min_percs[1],
                             idx=0, type='cg',
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      
      open_sum2_ssh <- 
        ssetToPassPercSsheet(sset=open_sset_dat, man=top_man_tib,
                             min=min_pvals[1], per=min_percs[1],
                             idx=0, type='cg',
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      
      if (retData) {
        ret_dat$open_sset_dat <- open_sset_dat
        ret_dat$open_beta_tib <- open_beta_tib
        ret_dat$open_sum1_ssh <- open_sum1_ssh
        ret_dat$open_sum2_ssh <- open_sum2_ssh
      }
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Initialize Data (Structures) Tables::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Build Sample Sheet::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Initializing Auto-Sample-Sheet: [idat/bead/pheno].{RET}"))
    
    beadPool     <- bead_ssh_tib  %>% head(n=1) %>% dplyr::select(Bead_Pool)   %>% dplyr::pull()
    chipFormat   <- idat_info_tib %>% head(n=1) %>% dplyr::select(Chip_Format) %>% dplyr::pull()
    version_key  <- man_map_tib$version[1]
    platform_key <- man_map_tib$platform[1]
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} chipFormat={chipFormat}, beadPool={beadPool}.{RET}"))
    
    ssheet_tib <- NULL
    ssheet_tib <- dplyr::bind_cols(
      # Idat Sample Sheet Stats::
      idat_info_tib,
      
      # Options Sample Sheet Stats::
      opt_ssh_tib,
      
      # Bead Sample Sheet Stats::
      bead_ssh_tib,
      
      # Manifest-Detection Sample Sheet Stats::
      man_map_tib %>% purrr::set_names(paste('detect',names(.), sep='_')) %>% head(n=1),
      
      # Basic p-value requeue calling:::
      open_sum1_ssh,
      open_sum2_ssh
    )
    
    if (retData) {
      ret_dat$ssheet_tib <- ssheet_tib
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Output Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Defining outputs...{RET}"))
    
    if (opts$buildSubDir) opts$outDir <- file.path(opts$outDir, chipFormat, beadPool)
    if (!dir.exists(opts$outDir)) dir.create(opts$outDir, recursive=TRUE)
    
    basecode   <- idat_info_tib$Sentrix_Name[1]
    out_name   <- paste(basecode, platform_key, version_key, sep=del)
    out_prefix <- file.path(opts$outDir, out_name)
    
    times_csv    <- paste(out_prefix, 'run-times.csv.gz', sep=del)
    ssheet_csv   <- paste(out_prefix, 'AutoSampleSheet.csv.gz', sep=del)
    dsheet_csv   <- paste(out_prefix, 'AutoSampleSheetDescriptionTable.csv.gz', sep=del)
    requeue_csv  <- paste(out_prefix, 'AutoQC.csv', sep=del)
    new_sset_rds <- paste(out_prefix, 'new.sset.rds', sep=del)
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Defined output files: out_prefix={out_prefix}, basecode={basecode}, out_name={out_name}.{RET}"))
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Defined output files: platform_key={platform_key}, version_key={version_key}, outDir={opts$outDir}.{RET}"))
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}  new_sset_rds={new_sset_rds}.{RET}") )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}    ssheet_csv={ssheet_csv}.{RET}") )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}    dsheet_csv={dsheet_csv}.{RET}") )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}     times_csv={times_csv}.{RET}{RET}{RET}") )
    if (verbose>=vt) cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Initialize RAW SSET::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Building RAW SSET...{RET}"))
    
    new_sset <- NULL
    new_sset <- newSset(prefix=prefix, 
                        platform=platform_key, manifest=top_man_tib,
                        load=opts$load_sset,save=opt$save_sset,rds=new_sset_rds,
                        verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    
    if (retData) {
      ret_dat$new_sset <- new_sset
      # return(ret_dat)
    }

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 SSET to Calls by Order of Operations:: workflows
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt) {
      cat(glue::glue("{tabsStr}{TAB}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting Workflows...{RET}{RET}"))
    }
    
    org_dat_bool <- TRUE
    cur_dat_list <- NULL
    workflow_cnt <- length(workflows)
    for (idx in seq(1,workflow_cnt)) {
      cur_workflow <- workflows[idx]
      if (verbose>=vt) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting idx={idx}, cur_workflow={cur_workflow}.{RET}"))
      }
      
      cur_sset_rds <- 
        paste(out_prefix, paste(cur_workflow,'sset.rds', sep='.'), sep=del)

      cur_sset <- NULL
      if (opts$load_sset && file.exists(cur_sset_rds)) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading RDS={cur_sset_rds}.{RET}"))
        cur_sset <- readr::read_rds(cur_sset_rds)
      } else {
        stopifnot(!is.null(new_sset))
        
        cur_sset <- mutateSSET_workflow(
          sset=new_sset, workflow=cur_workflow, pvals=pvals,
          save=opt$save_sset, rds=cur_sset_rds,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      }
      stopifnot(!is.null(cur_sset))
      
      cur_dat_list = ssetToSummary(
        sset=cur_sset, man=top_man_tib, idx=idx, workflow=cur_workflow,
        name=out_name, outDir=opts$outDir, pre=cur_dat_list, ref=ref,
        mask=mask,
        pvals=pvals, min_pvals=min_pvals, min_percs=min_percs,
        basic=open_beta_tib,
        
        write_sset=opt$save_sset, sset_rds=NULL, ret_sset=retData2,
        
        write_beta=opts$write_beta, beta_csv=NULL, ret_beta=retData2,
        write_bsum=opts$write_bsum, bsum_csv=NULL, ret_bsum=retData2,
        
        write_pval=opts$write_pval, pval_csv=NULL, ret_pval=retData2,
        write_psum=opts$write_psum, psum_csv=NULL, ret_psum=retData2,
        
        write_sigs=opts$write_sigs, sigs_csv=NULL, ret_sigs=retData2,
        write_ssum=opts$write_ssum, ssum_csv=NULL, ret_ssum=retData2,
        
        write_call=opts$write_call, call_csv=NULL, ret_call=retData2,
        write_csum=opts$write_csum, csum_csv=NULL, ret_csum=retData2,
        
        write_snps=opts$write_snps, snps_csv=NULL, ret_snps=retData2,
        write_auto=opts$write_auto, plot_auto=opts$plotAuto,
        makek_pred=opts$make_pred,
        
        percision_sigs=opts$percision_sigs,
        percision_beta=opts$percision_beta,
        percision_pval=opts$percision_pval,
        
        by="Probe_ID", type="Probe_Type", des="Probe_Design",
        minDb=opts$minDeltaBeta, dpi=opts$dpi, plotFormat=opts$plotFormat,
        datIdx=5, non_ref=FALSE, fresh=opts$fresh,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      
      if (retData && org_dat_bool) {
        org_dat_bool <- FALSE
        ret_dat$org_list <- cur_dat_list
      }

      if (verbose>=vt) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. cur_workflow={cur_workflow}.{RET}{RET}"))
        cat(glue::glue("{tabsStr}{TAB}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
      }
      
      # break
    }
    if (verbose>=vt) {
      cat(glue::glue("Done Workflows...{RET}"))
      cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    }
    
    if (retData) {
      ret_dat$cur_list <- cur_dat_list
      # return(ret_dat)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Format Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Arrange the Sample Sheet Columns
    #   cur_dat_list
    ssheet_tib <- dplyr::bind_cols(ssheet_tib, cur_dat_list$sums_ssh) %>%
      dplyr::select(-starts_with('Iscan_'), starts_with('Iscan_')) %>% 
      dplyr::mutate_if(is.numeric, list(round), opts$percision_beta)
    
    # Add Description Sample Sheet Table:: dsheet_tab
    dsheet_tab <- NULL
    dsheet_tab <- getSsheetDataTab(
      tib=ssheet_tib,
      minOobPval=min_pvals[1], minOobPerc=min_percs[1],
      minNegPval=min_pvals[2], minNegPerc=min_percs[2], 
      minDb=opts$minDeltaBeta,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    
    # Format Time Table
    time_tib <- tTracker$time %>% dplyr::mutate_if(base::is.numeric, round, 4)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing Outputs.{RET}"))
    
    if (!is.null(ssheet_tib)) readr::write_csv(ssheet_tib, ssheet_csv)
    if (!is.null(dsheet_tab)) readr::write_csv(dsheet_tab, dsheet_csv)
    
    req_tib <- NULL
    if (is.null(open_sum1_ssh)) {
      req_tib <- requeueFlag(
        tib=cur_dat_list$sums_dat, name=basecode,
        csv=requeue_csv,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    } else {
      req_tib <- requeueFlagOpenSesame(
        tib=ssheet_tib, name=basecode,
        csv=requeue_csv,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    }
    if (!is.null(time_tib))   readr::write_csv(time_tib, times_csv)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Compress Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    # if (verbose>=vt)
    #   cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Compressing Outputs (CSV)...{RET}"))
    #
    # gzip_list <- lapply(list.files(opt$outDir, pattern='.csv$', full.names=TRUE), 
    #                     function(x) { system(glue::glue("gzip -f {x}")) } )
    
    readr::write_lines(x=date(),file=stampEnd_txt,sep='\n',append=FALSE)
    
    ret_cnt <- ssheet_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tTracker)) tTracker$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))

  if (retData) {
    ret_dat$req_tib    <- req_tib
    ret_dat$ssheet_tib <- ssheet_tib
    ret_dat$dsheet_tab <- dsheet_tab
    ret_dat$time       <- tTracker$time
    return(ret_dat)
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
      tib <- maskTib(tib, field=fstr, id=id, pval=pstr, minPval=minPval, 
                     verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
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
        tib <- tib %>% dplyr::select(!!selKey,!!minKey, !!datKey)
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} USE_PVAL minVal={minVal}, ",
                         "SELECT=[selKey={selKey}, minKey={minKey}, datKey={datKey}].{RET}"))
          tib %>% print()
        }
        
        tot_cnt  <- tib %>% base::nrow()
        pre_cnt  <- tib %>% dplyr::filter(is.na(!!datKey)) %>% base::nrow()
        
        tib <- maskTib(tib=tib, id=selKey, field=datKey, pval=minKey, minPval=minVal, 
                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        
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

autoDetect_Wrapper = function(can, ref, man, mask=NULL,
                              
                              minPval, minDelta,
                              dname='Design_Type', pname='Probe_Type', ptype='cg',
                              # jval='Probe_ID', field='ind_beta', pval='ind_PnegEcdf', 
                              jval='Probe_ID', field, pval, suffix, del='_',
                              outDir=NULL, sname=NULL, plotMatrix=FALSE, writeMatrix=FALSE,
                              dpi=120, format="png", datIdx=4, non.ref=FALSE,
                              
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'autoDetect_Wrapper'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    auto_can <- NULL
    auto_ref <- NULL
    mask_can <- NULL
    
    auto_can <- man %>% dplyr::rename(Design_Type=DESIGN) %>% 
      # dplyr::select(Probe_ID, Design_Type, Probe_Type) %>%
      dplyr::select(dplyr::all_of( c(!!jval, !!dname, !!pname) ) ) %>%
      dplyr::left_join(can, by="Probe_ID")
    # print(auto_can)
    
    auto_ref <- ref
    if (grep("_beta", names(ref)) %>% length() == 0) {
      auto_ref <- ref %>% purrr::set_names(paste(names(.),'beta', sep=del) )
      names(auto_ref)[1] <- 'Probe_ID'
    }
    # print(auto_ref)
    
    ret_tib <- sampleDetect(
      can=auto_can, ref=auto_ref, minPval=minPval, minDelta=minDelta,
      dname=dname, pname=pname, ptype=ptype,
      jval=jval, field=field, pval=pval, suffix=suffix, del=del,
      outDir=outDir, sname=sname, plotMatrix=plotMatrix, writeMatrix=writeMatrix,
      dpi=dpi, format=format, datIdx=datIdx, non.ref=non.ref,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)

    if (!is.null(mask)) {
      jval_sym <- rlang::sym(jval)
      mask_can <- auto_can %>% dplyr::filter(! (!!jval_sym %in% mask) )

      mask_tib <- NULL
      mask_tib <- sampleDetect(
        can=mask_can, ref=auto_ref, minPval=minPval, minDelta=minDelta,
        dname=dname, pname=pname, ptype=ptype,
        jval=jval, field=field, pval=pval, suffix=suffix, del=del,
        outDir=outDir, sname=paste(sname,'Mask', sep='_'), 
        plotMatrix=plotMatrix, writeMatrix=writeMatrix,
        dpi=dpi, format=format, datIdx=datIdx, non.ref=non.ref,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        purrr::set_names(paste('Mask',names(.), sep='_'))

      ret_tib <- ret_tib %>% dplyr::bind_cols(mask_tib)
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  # ret_tib$can <- auto_can
  # ret_tib$ref <- auto_ref
  # ret_tib$dat <- ret_tib
  
  ret_tib
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
    jval  <- jval  %>% rlang::sym()
    pval  <- pval  %>% rlang::sym()
    
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
    dbm <- deltaMatrix(mat, minDelta=minDelta,
                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
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
    
    if (TRUE) {
      if (verbose>=vt+6) {
        cat(glue::glue("{RET}{RET}{RET}{RET}"))
        cat(glue::glue("[{funcTag}]:{tabsStr} TIB={RET}"))
        print(tib)
        cat(glue::glue("[{funcTag}]:{tabsStr} done(TIB){RET}"))
        cat(glue::glue("{RET}{RET}{RET}{RET}"))
      }
      
      # 4.1 Build matrix
      #
      tib1 <- tib %>% 
        dplyr::filter(!!dname == "I")
      mat1 <- tib1 %>%
        dplyr::select(-c(!!jval,!!dname,!!pname)) %>%
        dplyr::select(ends_with(paste0(del,suffix)) ) %>% as.matrix()
      if (verbose>=vt+6) mat1 %>% head(n=2) %>% print()
      
      # 5.1 Build R-Squared Matrix
      r2m_1_tib <- 
        rsquaredMatrix(mat1, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
        tibble::as_tibble(rownames='Sample')
      if (verbose>=vt+6) r2m_1_tib %>% head(n=2) %>% print()

      # 6.1 Build Deleta Matrix
      dbm_1_tib <- 
        deltaMatrix(mat1, minDelta=minDelta,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
        tibble::as_tibble(rownames='Sample')
      if (verbose>=vt+6) dbm_1_tib %>% head(n=2) %>% print()
      
      can_1_nrows <- tib1 %>% base::nrow()
      sss_1_tib <- dplyr::bind_cols(
        r2m_1_tib %>% head(n=1) %>% dplyr::select(-Sample) %>% 
          tidyr::gather(key='Sample', value='r2') %>% 
          dplyr::mutate(Sample=stringr::str_remove(Sample, paste0('_',suffix)) ) %>%
          dplyr::filter(Sample!='Sample') %>% 
          dplyr::arrange(-r2) %>% head(n=1) %>%
          purrr::set_names('AutoSample_R2_1_Key', 'AutoSample_R2_1_Val'),
        
        dbm_1_tib %>% head(n=1) %>% dplyr::select(-Sample) %>%
          tidyr::gather(key='Sample', value='passCount') %>%
          dplyr::mutate(Sample=stringr::str_remove(Sample, paste0('_',suffix)),
                        failCount=can_1_nrows-passCount,
                        passPerc=round(100*passCount/can_1_nrows,1) ) %>%
          dplyr::filter(Sample!='Sample') %>% 
          dplyr::arrange(-passPerc) %>% head(n=1) %>%
          purrr::set_names('AutoSample_dB_1_Key', 'AutoSample_Total_1_Cnt', 'AutoSample_dB_1_Cnt', 'AutoSample_dB_1_Val')
      ) %>% dplyr::select(AutoSample_Total_1_Cnt, AutoSample_R2_1_Key, AutoSample_R2_1_Val, 
                          AutoSample_dB_1_Key, AutoSample_dB_1_Cnt, AutoSample_dB_1_Val)
      if (verbose>=vt+6) print(sss_1_tib)

      # 4.2 Build matrix
      #
      tib2 <- tib %>% 
        dplyr::filter(!!dname == "II")
      mat2 <- tib2 %>%
        dplyr::select(-c(!!jval,!!dname,!!pname)) %>%
        dplyr::select(ends_with(paste0(del,suffix)) ) %>% as.matrix()
      if (verbose>=vt+6) mat2 %>% head(n=2) %>% print()
      
      # 5.2 Build R-Squared Matrix
      r2m_2_tib <- 
        rsquaredMatrix(mat2, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
        tibble::as_tibble(rownames='Sample')
      if (verbose>=vt+6) r2m_2_tib %>% head(n=2) %>% print()
      
      # 6.2 Build Deleta Matrix
      dbm_2_tib <- 
        deltaMatrix(mat2, minDelta=minDelta,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
        tibble::as_tibble(rownames='Sample')
      if (verbose>=vt+6) dbm_2_tib %>% head(n=2) %>% print()
      
      can_2_nrows <- tib2 %>% base::nrow()
      sss_2_tib <- dplyr::bind_cols(
        r2m_2_tib %>% head(n=1) %>% dplyr::select(-Sample) %>% 
          tidyr::gather(key='Sample', value='r2') %>% 
          dplyr::mutate(Sample=stringr::str_remove(Sample, paste0('_',suffix)) ) %>%
          dplyr::filter(Sample!='Sample') %>% 
          dplyr::arrange(-r2) %>% head(n=1) %>%
          purrr::set_names('AutoSample_R2_2_Key', 'AutoSample_R2_2_Val'),
        
        dbm_2_tib %>% head(n=1) %>% dplyr::select(-Sample) %>%
          tidyr::gather(key='Sample', value='passCount') %>%
          dplyr::mutate(Sample=stringr::str_remove(Sample, paste0('_',suffix)),
                        failCount=can_2_nrows-passCount,
                        passPerc=round(100*passCount/can_2_nrows,1) ) %>%
          dplyr::filter(Sample!='Sample') %>% 
          dplyr::arrange(-passPerc) %>% head(n=1) %>%
          purrr::set_names('AutoSample_dB_2_Key', 'AutoSample_Total_2_Cnt', 'AutoSample_dB_2_Cnt', 'AutoSample_dB_2_Val')
      ) %>% dplyr::select(AutoSample_Total_2_Cnt, AutoSample_R2_2_Key, AutoSample_R2_2_Val, 
                          AutoSample_dB_2_Key, AutoSample_dB_2_Cnt, AutoSample_dB_2_Val)
      if (verbose>=vt+6) print(sss_2_tib)
      
      sss <- sss %>% dplyr::bind_cols(sss_1_tib) %>% dplyr::bind_cols(sss_2_tib)
    }
    
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
