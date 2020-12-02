
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
sesamizeSingleSample = function(prefix, man, add, ref, opts, defs=NULL,
                                pvals=NULL, min_pvals=NULL, min_percs=NULL,
                                workflows=NULL, 
                                
                                btypes=c('cg','ch','rs'),
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
  ret <- NULL
  
  tTracker <- timeTracker$new()
  stime <- system.time({
    
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
    
    opt_ssh_tib <- tibble::tibble(minDeltaBeta = opts$minDeltaBeta)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Extract Raw idat::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (!dir.exists(opts$outDir)) dir.create(opts$outDir, recursive=TRUE)
    idat_rds  <- paste(out_prefix, 'idat.rds', sep=del)
    idat_list <- prefixToIdat(prefix=prefix, 
                              load=opts$load_idat, save=opts$save_idat, rds=idat_rds,
                              verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    # return(idat_list)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Get Sesame Manifest/Address Tables::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    man_map_tib <- idatToManifestMap(tib=idat_list$sig, mans=mans, sortMax=TRUE,
                                     verbose=verbose,tc=tc+1,tt=tTracker)
    
    top_man_tib <- mans[[man_map_tib[1,]$manifest]] %>% 
      dplyr::distinct(Probe_ID, .keep_all=TRUE)
    
    bead_sum_tib <- 
      getManifestBeadStats(dat=idat_list$sig, man=top_man_tib, types=btypes, 
                           verbose=verbose,tc=tc+1,tt=tTracker)
    
    #
    # TBD:: Replace or improve addBeadPoolToSampleSheet()
    #
    bead_ssh_tib <- bead_sum_tib %>% tidyr::spread(name, value) %>%
      addBeadPoolToSampleSheet(field="Loci_Count_cg",
                               verbose=verbose,tc=tc+1,tt=tTracker)
    
    if (retData) {
      ret$prefix <- prefix
      
      ret$ossh <- opt_ssh_tib
      
      ret$isig <- idat_list$sig
      ret$issh <- idat_list$ann
      
      ret$mmap <- man_map_tib
      ret$sman <- top_man_tib
      
      ret$bsum <- bead_sum_tib
      ret$bssh <- bead_ssh_tib
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
    chipFormat   <- idat_list$ann %>% head(n=1) %>% dplyr::select(Chip_Format) %>% dplyr::pull()
    version_key  <- man_map_tib$version[1]
    platform_key <- man_map_tib$platform[1]
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} chipFormat={chipFormat}, beadPool={beadPool}.{RET}"))
    
    ssheet_tib <- NULL
    ssheet_tib <- dplyr::bind_cols(
      # Idat Sample Sheet Stats::
      idat_list$ann,
      
      # Options Sample Sheet Stats::
      opt_ssh_tib,
      
      # Bead Sample Sheet Stats::
      bead_ssh_tib,
      
      # Manifest-Detection Sample Sheet Stats::
      man_map_tib %>% purrr::set_names(paste('detect',names(.), sep='_')) %>% head(n=1)
    )
    
    if (retData) {
      ret$ssheet_tib <- ssheet_tib
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Output Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Defining outputs...{RET}"))
    
    if (opts$buildSubDir) opts$outDir <- file.path(opts$outDir, chipFormat, beadPool)
    if (!dir.exists(opts$outDir)) dir.create(opts$outDir, recursive=TRUE)
    
    basecode   <- idat_list$ann$Sentrix_Name[1]
    out_name   <- paste(basecode, platform_key, version_key, sep=del)
    out_prefix <- file.path(opts$outDir, out_name)
    
    times_csv    <- paste(out_prefix, 'run-times.csv.gz', sep=del)
    ssheet_csv   <- paste(out_prefix, 'AutoSampleSheet.csv.gz', sep=del)
    dsheet_csv   <- paste(out_prefix, 'AutoSampleSheetDescriptionTable.csv.gz', sep=del)
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
                        load=opts$load_sset,save=opts$save_sset,rds=new_sset_rds,
                        verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    
    if (retData) {
      ret$new_sset <- new_sset
      # return(ret)
    }
    
    ret_dat_list <- NULL
    cur_dat_list <- NULL
    
    auto_beta_key <- NULL
    auto_negs_key <- NULL

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 SSET to Calls by Order of Operations:: workflows
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt) {
      cat(glue::glue("{tabsStr}{TAB}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting Workflows...{RET}{RET}"))
    }
    
    workflow_cnt <- length(workflows)
    for (idx in seq(1,workflow_cnt)) {
      cur_workflow <- workflows[idx]
      if (verbose>=vt) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting idx={idx}, cur_workflow={cur_workflow}.{RET}"))
      }

      cur_sset_rds <- paste(out_prefix, paste(cur_workflow,'sset.rds', sep='.'), sep=del)
      if (is.null(auto_beta_key)) auto_beta_key <- paste(cur_workflow,'beta', sep=del)
      if (is.null(auto_negs_key)) auto_negs_key <- paste(cur_workflow,'PnegEcdf', sep=del)
      
      cur_sset <- NULL
      if (opts$load_sset && file.exists(cur_sset_rds)) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading RDS={cur_sset_rds}.{RET}"))
        cur_sset <- readr::read_rds(cur_sset_rds)
      } else {
        stopifnot(!is.null(new_sset))
        
        cur_sset <- mutateSSET_workflow(
          sset=new_sset, workflow=cur_workflow, pvals=pvals,
          save=opts$save_sset, rds=cur_sset_rds,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      }
      stopifnot(!is.null(cur_sset))
      
      cur_dat_list = ssetToSummary(
        sset=cur_sset, man=top_man_tib, idx=idx, workflow=cur_workflow,
        name=out_name, outDir=opts$outDir, pre=cur_dat_list,
        
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

      if (verbose>=vt) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. cur_workflow={cur_workflow}.{RET}"))
        cat(glue::glue("{tabsStr}{TAB}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
      }
      
      # break
    }
    if (verbose>=vt) {
      cat(glue::glue("Done Workflows...{RET}"))
      cat(glue::glue("# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
    }
    
    if (retData) {
      ret$cur_list <- cur_dat_list
      # return(ret)
    }

    #
    # Add formatVCF
    # Sum formatVCF
    #
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Auto-Detect Sample::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (opts$auto_detect) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Detecting Auto-Sample-Sheet: [SampleName, rsquared,deltaBeta].{RET}"))
      
      auto_sheet <- autoDetect_Wrapper(
        can=cur_dat_list$call_dat, ref=ref, man=top_man_tib,
        minPval=opts$minNegPval, minDelta=opts$minDeltaBeta,
        dname='Design_Type', pname='Probe_Type', ptype='cg',
        jval='Probe_ID', field=auto_beta_key, pval=auto_negs_key, suffix='beta', del=del,
        outDir=opts$outDir, sname=out_name, plotMatrix=opts$plotAuto, writeMatrix=opts$write_auto,
        dpi=opts$dpi, format=opts$plotFormat, datIdx=4, non.ref=non_ref,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      
      ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(auto_sheet)
      
      ssheet_ncols <- ssheet_tib %>% base::ncol()
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Binding Sample Sheet (with auto-detect-sample) ssheet_ncols={ssheet_ncols}.{RET}{RET}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Format Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Arrange the Sample Sheet Columns
    ssheet_tib <- ssheet_tib %>% 
      dplyr::select(-starts_with('Iscan_'), starts_with('Iscan_')) %>% 
      dplyr::mutate_if(is.numeric, list(round), 4) # %>% 
      # dplyr::select(Requeue_Flag_pOOBAH,Requeue_Flag_PnegEcdf, dplyr::everything())
    
    # Add Description Sample Sheet Table:: dsheet_tab
    dsheet_tab <- NULL
      # getSsheetDataTab(tib = ssheet_tib, 
      #                  minOobPval=opts$minOobPval, minOobPerc=opts$minOobPerc,
      #                  minNegPval=opts$minNegPval, minNegPerc=opts$minNegPerc,
      #                  verbose = 1)
    
    # Format Time Table
    time_tib <- tTracker$time %>% dplyr::mutate_if(base::is.numeric, round, 4)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing Outputs.{RET}"))
    
    if (!is.null(ssheet_tib)) readr::write_csv(ssheet_tib, ssheet_csv)
    if (!is.null(dsheet_tab)) readr::write_csv(dsheet_tab, dsheet_csv)
    if (!is.null(time_tib))   readr::write_csv(time_tib, times_csv)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Compress Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Compressing Outputs (CSV)...{RET}"))
    
    # gzip_list <- lapply(list.files(opt$outDir, pattern='.csv$', full.names=TRUE), 
    #                     function(x) { system(glue::glue("gzip -f {x}")) } )
    
    readr::write_lines(x=date(),file=stampEnd_txt,sep='\n',append=FALSE)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (verbose>=vt+3) tTracker %>% print()
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. elasped={etime}.{RET}{RET}"))
  
  if (retData) {
    ret$ssheet_tib <- ssheet_tib
    ret$dsheet_tab <- dsheet_tab
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
