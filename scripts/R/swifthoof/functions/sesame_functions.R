
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

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Sesame Mutation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssetToSummary = function(sset, man, idx, workflow, name, outDir=NULL, pre=NULL,
                         pvals=NULL, min_pvals=NULL, min_percs=NULL,
                         
                         write_sset=FALSE, sset_rds=NULL, ret_sset=FALSE,
                         
                         write_beta=FALSE, beta_csv=NULL, ret_beta=FALSE,
                         write_bsum=FALSE, bsum_csv=NULL, ret_bsum=FALSE,

                         write_pval=FALSE, pval_csv=NULL, ret_pval=FALSE,
                         write_psum=FALSE, psum_csv=NULL, ret_psum=FALSE,
                         
                         write_sigs=FALSE, sigs_csv=NULL, ret_sigs=FALSE,
                         write_ssum=FALSE, ssum_csv=NULL, ret_ssum=FALSE,
                         
                         write_call=FALSE, call_csv=NULL, ret_call=FALSE,
                         write_csum=FALSE, csum_csv=NULL, ret_csum=FALSE,
                         
                         percision_sigs=-1,percision_beta=-1,percision_pval=-1,
                         by="Probe_ID", type="Probe_Type", des="Probe_Design",
                         fresh=FALSE,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}...{RET}"))

  ret_cnt <- 0
  ret_dat <- NULL

  # Validate Inputs::
  #
  pvals_cnt <- 0
  min_pvals_cnt <- 0
  min_percs_cnt <- 0
  if (!is.null(pvals)) pvals_cnt <- length(pvals)
  if (!is.null(min_pvals)) min_pvals_cnt <- length(min_pvals)
  if (!is.null(min_percs)) min_percs_cnt <- length(min_percs)
  
  if (pvals_cnt != min_pvals_cnt || min_pvals_cnt != min_percs_cnt) {
    stop(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR; pval vecs not equal; pvals_cnt={pvals_cnt}, ",
                    "min_pvals_cnt={min_pvals_cnt}, min_percs_cnt={min_percs_cnt}!!!{RET}"))
    return(ret_dat)
  }
  
  sigs_dat_tib <- NULL
  call_dat_tib <- NULL
  sums_dat_tib <- NULL
  
  stime <- system.time({
    
    # Initialize Outputs::
    #
    if (!is.null(name) && !is.null(outDir)) {
      if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
      
      out_name <- paste(name,workflow, sep='_')
      if (is.null(sset_rds)) 
        sset_rds <- file.path(outDir, paste(out_name,'sset.rds', sep='.'))
      
      if (is.null(beta_csv))
        beta_csv <- file.path(outDir, paste(out_name,'beta.dat.csv.gz', sep='.'))
      if (is.null(bsum_csv))
        bsum_csv <- file.path(outDir, paste(out_name,'beta.sum.csv.gz', sep='.'))
      
      if (is.null(sigs_csv))
        sigs_csv <- file.path(outDir, paste(out_name,'sigs.dat.csv.gz', sep='.'))
      if (is.null(ssum_csv))
        ssum_csv <- file.path(outDir, paste(out_name,'sigs.sum.csv.gz', sep='.'))
      
      if (is.null(call_csv)) 
        call_csv <- file.path(outDir, paste(out_name,'call.dat.csv.gz', sep='.'))
      if (is.null(csum_csv)) 
        csum_csv <- file.path(outDir, paste(out_name,'call.sum.csv.gz', sep='.'))
      
      # if (write_sigs) sigs_csv <- clean_file(sigs_csv, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
    
    # Initialize Sample Sheet::
    #
    # sam_sheet <- tibble::tibble(
    #   Method_Key=workflow,
    #   Method_Idx=idx
    # )
    
    # Clear sset to ensure there's no code mistakes
    sset_dat <- sset
    sset <- NULL
    
    call_dat_tib <- man %>% dplyr::select(!!by)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Set/Update/Summarize:: Pvals
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_dat$pval <- NULL
    pval_cut_tib <- tibble::tibble(Method=pvals, min_pval=min_pvals, min_perc=min_percs)
    for (pval_key in pval_cut_tib$Method) {
      
      pval_out_str <- paste('pval',pval_key, sep='-')
      pval_csv <- file.path(outDir, paste(out_name,pval_out_str,'dat.csv.gz', sep='.'))
      psum_csv <- file.path(outDir, paste(out_name,pval_out_str,'sum.csv.gz', sep='.'))

      ret_pval_dat <- NULL
      pval_tib <- pval_cut_tib %>% dplyr::filter(Method==pval_key)
      min_pval <- pval_tib$min_pval[1]
      min_perc <- pval_tib$min_perc[1]
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pval={pval_key}; min={min_pval}, perc={min_perc}.{RET}"))
      
      pval_dat_tib <- ssetToTib(
        sset=sset_dat, source='pvals', name=pval_key,
        percision=percision_pval, sort=FALSE,
        save=write_pval, csv=pval_csv,
        by=by, type=type, des=des,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      if (ret_pval) ret_pval_dat$pval_dat <- pval_dat_tib
      
      pval_sum_tib <- ssetTibToSummary(
        tib=pval_dat_tib,man=man,
        pval=min_pval, perc=min_perc, percision=percision_pval,
        save=write_psum, csv=psum_csv,
        by=by, type=type, des=des,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      if (ret_psum) ret_pval_dat$pval_sum <- pval_sum_tib
      if (ret_pval || ret_psum) ret_dat$pval[[pval_key]] = ret_pval_dat
      
      call_dat_tib <- dplyr::left_join(call_dat_tib,pval_dat_tib, by=by)
      sums_dat_tib <- sums_dat_tib %>% dplyr::bind_rows(pval_sum_tib)
      
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Summarize:: Betas
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_beta_dat <- NULL
    beta_dat_tib <- ssetToTib(
      sset=sset_dat, source='betas',
      percision=percision_beta, sort=TRUE,
      save=write_beta, csv=beta_csv,
      by=by, type=type, des=des,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (ret_beta) ret_beta_dat$beta_dat <- beta_dat_tib
    
    beta_sum_tib <- ssetTibToSummary(
      tib=beta_dat_tib,man=man,
      percision=percision_beta,
      save=write_bsum, csv=bsum_csv,
      by=by, type=type, des=des,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (ret_bsum) ret_beta_dat$beta_sum <- beta_sum_tib
    
    if (ret_beta || ret_bsum) ret_dat$beta <- ret_beta_dat
    
    call_dat_tib <- dplyr::left_join(call_dat_tib,beta_dat_tib, by=by)
    sums_dat_tib <- sums_dat_tib %>% dplyr::bind_rows(beta_sum_tib)

    # Write Updated SSET
    #
    if (write_sset && !is.null(sset_rds) && !file.exists(sset_rds))
      readr::write_rds(sset_dat, sset_rds, compress="gz")
    
    # Write Calls CSV
    #
    if (write_call && !is.null(call_csv))
      readr::write_csv(call_dat_tib, call_csv)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Set/Summarize:: Sigs
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_sigs_dat <- NULL
    sigs_dat_tib <- ssetToTib(
      sset=sset_dat, man=man, 
      source='sigs', name=pval,
      percision=percision_sigs, sort=FALSE, 
      save=write_sigs, csv=sigs_csv, 
      by=by, type=type, des=des,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (ret_sigs) ret_sigs_dat$sigs_dat <- sigs_dat_tib
    
    sigs_sum_tib <- ssetTibToSummary(
      tib=sigs_dat_tib,
      percision=percision_sigs,
      save=write_ssum, csv=ssum_csv, 
      by=by, type=type, des=des,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (ret_ssum) ret_sigs_dat$sigs_sum <- sigs_sum_tib

    if (ret_sigs || ret_ssum) ret_dat$sigs <- ret_sigs_dat
    
    sums_dat_tib <- sums_dat_tib %>% dplyr::bind_rows(sigs_sum_tib)

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Predicted and Inferred Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    pred_dat_tib <- NULL
    pred_dat_tib <- ssetToPredictions(
      sset=sset_dat,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
      purrr::set_names(paste(workflow,names(.), sep='_'))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                Update Calls, Signal and Summary Headers::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    fix_cols <- c(by,type,des)
    sigs_dat_tib <- addColNames(sigs_dat_tib, add=workflow, fix=fix_cols,
                                verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)

    fix_cols <- c(by)
    call_dat_tib <- addColNames(call_dat_tib, add=workflow, fix=fix_cols,
                                verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)

    sums_dat_tib <- sums_dat_tib %>%
      dplyr::mutate(Workflow_key=workflow, Workflow_idx=idx) %>% 
      dplyr::select(Workflow_key, Workflow_idx, dplyr::everything())
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Gather/Merge Results:: Sample Sheet
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (is.null(pre)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Building new data list...{RET}"))
        
      ret_dat$sigs_dat <- sigs_dat_tib
      ret_dat$call_dat <- call_dat_tib
      ret_dat$sums_dat <- sums_dat_tib
      ret_dat$pred_dat <- pred_dat_tib
    } else {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Join with previous data...{RET}"))
        
      ret_dat$sigs_dat <- dplyr::left_join(pre$sigs_dat, sigs_dat_tib, by=c(by,type,des))
      ret_dat$call_dat <- dplyr::left_join(pre$call_dat, call_dat_tib, by=by)
      ret_dat$sums_dat <- dplyr::bind_rows(pre$sums_dat, sums_dat_tib)
      ret_dat$pred_dat <- dplyr::bind_cols(pre$pred_dat, pred_dat_tib)
    }
    
    ret_cnt <- ret_dat %>% names %>% length()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_dat({ret_cnt})={RET}"))
      print(ret_dat)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
                   
  ret_dat
}

mutateSSET_workflow = function(sset, workflow, pvals=NULL,
                               save=FALSE, rds=NULL, 
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutateSSET_workflow'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}.{RET}"))
  
  ret_cnt <- 0
  stime <- system.time({
    
    # Clear extra fields
    sesame::extra(sset)[['pvals']] <- NULL
    sesame::extra(sset)[['betas']] <- NULL
    
    #
    # TBD:: This could be parsed into individual letters and then processed in that order...
    #
    if (workflow=='r' || workflow=='raw') {
      # Do nothing...
    } else if (workflow=='i') {
      sset <- mutateSesame(sset=sset, method = 'inferTypeIChannel',
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='n') {
      sset <- mutateSesame(sset=sset, method = 'noob', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='d') {
      sset <- mutateSesame(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='nd') {
      sset <- mutateSesame(sset=sset, method = 'noob', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSesame(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='dn') {
      sset <- mutateSesame(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSesame(sset=sset, method = 'noob', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='ind') {
      sset <- mutateSesame(sset=sset, method = 'inferTypeIChannel', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSesame(sset=sset, method = 'noob', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSesame(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='ndi') {
      sset <- mutateSesame(sset=sset, method = 'noob', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSesame(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSesame(sset=sset, method = 'inferTypeIChannel', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='din') {
      sset <- mutateSesame(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSesame(sset=sset, method = 'inferTypeIChannel', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSesame(sset=sset, method = 'noob', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='dni') {
      sset <- mutateSesame(sset=sset, method = 'dyeBiasCorrTypeINorm', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSesame(sset=sset, method = 'noob', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      sset <- mutateSesame(sset=sset, method = 'inferTypeIChannel', 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else {
      stop(glue::glue("[{funcTag}]: ERROR: Unsupported workflow={workflow}!{RET}{RET}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Set:: Pvals
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    for (pval in pvals) {
      if (verbose>=vt+2) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Mutating pval={pval}...{RET}"))
      
      sset <- mutateSesame(sset=sset, method=pval,
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Set:: Betas
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt+2) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Mutating betas...{RET}"))
    
    sset <- mutateSesame(sset=sset, method="betas",
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Save Output::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (save && !is.null(rds)) {
      if (verbose>=vt+1) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing RDS={rds}.{RET}"))
      readr::write_rds(sset, rds, compress="gz")
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
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  sset
}

mutateSesame = function(sset, method, 
                        quality.mask = FALSE, nondetection.mask = FALSE, 
                        correct.switch = TRUE, mask.use.tcga = FALSE, 
                        pval.threshold = 1, force=TRUE,
                        pval.method = "pOOBAH", sum.TypeI = FALSE,
                        verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutateSesame'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; Mutate Sesame({method})...{RET}"))
  
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} sset={RET}"))
  if (verbose>=vt+4) print(sset)

  ctl_cnt <- sset@ctl %>% base::nrow()
  if (ctl_cnt==0 && method=='detectionPnegEcdf') {
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Unable to mutate ctl_cnt={ctl_cnt} ",
                     "for method={method}. Returning original sset...{RET}"))
    return(sset)
  }
  
  ret_cnt <- 0
  stime <- system.time({
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} sset(ctl={ctl_cnt})={RET}"))
    if (verbose>=vt+4) print(sset)
    
    if (is.null(method)) stop(glue::glue("{RET}[{funcTag}]: ERROR: Missing method!!!{RET}{RET}"))
    else if (method=='open') 
      sset <- sset %>% sesame::pOOBAH(force=force) %>% sesame::noob() %>% sesame::dyeBiasCorrTypeINorm()
    else if (method=='dyeBiasCorrTypeINorm') 
      sset <- sset %>% sesame::dyeBiasCorrTypeINorm()
    else if (method=='detectionPnegEcdf' || method=='PnegEcdf') 
      sset <- sset %>% sesame::detectionPnegEcdf(force=force)
    else if (method=='pOOBAH') 
      sset <- sset %>% sesame::pOOBAH(force=force)
    else if (method=='noob') 
      sset <- sset %>% sesame::noob()
    else if (method=='noobsb') 
      sset <- sset %>% sesame::noobsb()
    else if (method=='inferTypeIChannel') 
      sset <- sset %>% sesame::inferTypeIChannel(switch_failed=FALSE, verbose=FALSE)
    else if (method=='betas') 
      sesame::extra(sset)[[method]] <- getBetas2(sset=sset, mask=quality.mask,sum.TypeI=sum.TypeI)
    else if (method=='raw') { } # sset <- sset
    else stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported method={method}!!!{RET}{RET}"))
    
    ret_cnt <- sset %>% slotNames() %>% length()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} sset(slots={ret_cnt})={RET}"))
      print(sset)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  sset
}

# End of file
