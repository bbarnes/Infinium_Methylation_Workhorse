
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

ssetToSummary = function(sset, man, idx, workflow, name=NULL, outDir=NULL,
                         write_sset=FALSE, sset_rds=NULL, ret_sset=FALSE,
                         write_sigs=FALSE, sigs_csv=NULL, ret_sigs=FALSE,
                         write_ssum=FALSE, ssum_csv=NULL, ret_ssum=FALSE,
                         write_call=FALSE, call_csv=NULL, ret_call=FALSE,
                         write_csum=FALSE, csum_csv=NULL, ret_csum=FALSE,
                         
                         minOobPval=0.1,minNegPval=0.02,
                         percision_sigs=-1,percision_beta=-1,percision_pval=-1,
                         by="Probe_ID", type="Probe_Type", des="Probe_Class", 
                         fresh=FALSE,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}...{RET}"))
  
  ret_cnt <- 0
  ret_dat <- NULL
  stime <- system.time({
    
    # Initialize Outputs::
    #
    if (!is.null(name) && !is.null(outDir)) {
      if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
      
      out_name <- paste(name,workflow, sep='_')
      if (is.null(sset_rds)) 
        sset_rds <- file.path(outDir, paste(out_name,'sset.rds', sep='.'))
      
      if (is.null(sigs_csv)) 
        sigs_csv <- file.path(outDir, paste(out_name,'sigs.dat.csv.gz', sep='.'))
      if (is.null(ssum_csv)) 
        ssum_csv <- file.path(outDir, paste(out_name,'sigs.sum.csv.gz', sep='.'))
      
      if (is.null(call_csv)) 
        call_csv <- file.path(outDir, paste(out_name,'call.dat.csv.gz', sep='.'))
      if (is.null(csum_csv)) 
        csum_csv <- file.path(outDir, paste(out_name,'call.sum.csv.gz', sep='.'))
    }
    
    # Set/Update Pvals and Betas::
    #
    sset <- mutateSesame(sset=sset, method="betas", 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    sset <- mutateSesame(sset=sset, method="detectionPnegEcdf", 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    sset <- mutateSesame(sset=sset, method="pOOBAH", 
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    if (write_sset && !is.null(sset_rds))
      readr::write_rds(sset, sset_rds, compress="gz")
    
    # Initialize Sample Sheet::
    #
    sam_sheet <- tibble::tibble(
      Method_Key=workflow,
      Method_Idx=idx
    )
    
    # Signal Summary
    #
    sigs_dat_tib <- NULL
    sigs_dat_tib <- ssetToSigsTib(
      sset=sset, man=man,
      percision=percision_sigs, sort=FALSE, 
      save=write_sigs, csv=sigs_csv, 
      by=by, type=type, des="Probe_Class", 
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    sigs_sum_tib <- NULL
    sigs_sum_tib <- sigsTibToSummary(
      tib=sigs_dat_tib,
      save=write_ssum, csv=ssum_csv, 
      by=by, type=type, des="Probe_Class", 
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    sigs_ssheet_tib <- NULL
    sigs_ssheet_tib <- sigsSumToSSheet(
      tib=sigs_sum_tib, metric='mean', 
      by=by, type=type, des="Probe_Class", 
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # Calls Data: pvals/betas
    #
    call_dat_tib <- NULL
    call_dat_tib <- ssetToCallTib(
      sset=sset, workflow=workflow, fresh=fresh,
      save=write_call, csv=call_csv, 
      percision_beta=percision_beta, 
      percision_pval=percision_pval, 
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # Calls Summary:: betas
    #
    beta_metric <- 'mean'
    betas_cols  <- call_dat_tib %>% 
      dplyr::select(!!by,dplyr::ends_with('_beta')) %>% names()
    
    call_beta_sum_tib <- call_dat_tib %>% 
      dplyr::select(dplyr::all_of(betas_cols)) %>%
      sigsTibToSummary(man=man, by=by, type=type, des="Probe_Design", 
                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    call_beta_ssheet_tib <- sigsSumToSSheet(
      tib=call_beta_sum_tib, metric=beta_metric,
      by=by, type=type, des="Probe_Design", 
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # Calls Summary:: pvals
    #
    pval_metric <- 'pass_perc'
    pvals_cols  <- call_dat_tib %>% 
      dplyr::select(!dplyr::ends_with('_beta')) %>% names()
    
    call_pval_ssheet_tib <- NULL
    for (pval_idx in c(2:length(pvals_cols)) ) {
      pvals_col <- pvals_cols[pval_idx]
      pvals_key <- pvals_col %>% stringr::str_remove('^[^_]+_')
      if (verbose>=vt+2)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pval_idx={pval_idx}; pvals_col={pvals_col}...{RET}"))
      
      call_pval_sel_tib <- call_dat_tib %>% 
        dplyr::select(!!by,!!pvals_col)
      if (verbose>=vt+4)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} call_pval_sel_tib={RET}"))
      if (verbose>=vt+4) print(call_pval_sel_tib)
      
      call_pval_sum_tib <- call_pval_sel_tib %>%
        sigsTibToSummary(man=man, by=by, type=type, des="Probe_Design", 
                         cutoff=minNegPval,
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      if (verbose>=vt+4)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} call_pval_sum_tib={RET}"))
      if (verbose>=vt+4) print(call_pval_sum_tib)
      
      cur_pval_ssheet <- 
        sigsSumToSSheet(
          tib=call_pval_sum_tib, metric=pval_metric,
          by="Probe_ID", type="Probe_Type", des="Probe_Design", 
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        purrr::set_names(paste(pvals_key,names(.), sep='_'))
      
      call_pval_ssheet_tib <- call_pval_ssheet_tib %>% 
        dplyr::bind_cols(cur_pval_ssheet)
    }
    if (verbose>=vt+4)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} call_pval_ssheet_tib={RET}"))
    if (verbose>=vt+4) print(call_pval_ssheet_tib)
    
    # Predicted and Inferred Data::
    #
    pred_sum_tib <- NULL
    pred_sum_tib <- ssetToPredictions(
      sset=sset, fresh=fresh,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # Gather all sample sheets::
    #
    sam_sheet <- sam_sheet %>% dplyr::bind_cols(
      pred_sum_tib, call_pval_ssheet_tib, call_beta_ssheet_tib, sigs_ssheet_tib
    ) %>% purrr::set_names(paste(names(.),idx, sep='_'))
    
    # Set Outputs::
    #
    if (ret_sset) ret_dat$sset_dat <- sset
    if (ret_sigs) ret_dat$sigs_dat <- sigs_dat_tib
    if (ret_ssum) ret_dat$sigs_sum <- sigs_sum_tib
    
    if (ret_call) ret_dat$call_dat <- call_dat_tib
    if (ret_csum) ret_dat$beta_sum <- call_beta_sum_tib
    if (ret_csum) ret_dat$pval_sum <- call_pval_sum_tib
    
    if (is.null(ret_dat)) {
      ret_dat <- sam_sheet
      ret_cnt <- ret_dat %>% base::ncol()
    } else {
      ret_dat$sam_sheet <- sam_sheet
      ret_cnt <- ret_dat %>% names %>% length()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_dat
}

mutateSSET_workflow = function(sset, workflow, save=FALSE, rds=NULL, 
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutateSSET_workflow'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}.{RET}"))
  
  stime <- system.time({
    
    #
    # TBD:: This could be parsed into individual letters and then processed in that order...
    #
    if (workflow=='r' || workflow=='raw') {
    } else if (workflow=='i') {
      sset <- sset %>% 
        mutateSesame(method = 'inferTypeIChannel', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt)
    } else if (workflow=='n') {
      sset <- sset %>% 
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='d') {
      sset <- sset %>% 
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt)
    } else if (workflow=='nd') {
      sset <- sset %>% 
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt)
    } else if (workflow=='dn') {
      sset <- sset %>% 
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt)
    } else if (workflow=='ind') {
      sset <- sset %>% 
        mutateSesame(method = 'inferTypeIChannel', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt)
    } else if (workflow=='ndi') {
      sset <- sset %>% 
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'inferTypeIChannel', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt)
    } else if (workflow=='din') {
      sset <- sset %>% 
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'inferTypeIChannel', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='dni') {
      sset <- sset %>% 
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'inferTypeIChannel', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+3,tc=tc+1,tt=tt)
    } else {
      stop(glue::glue("[{funcTag}]: ERROR: Unsupported workflow={workflow}!{RET}{RET}"))
    }
    
    if (save && !is.null(rds)) {
      if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing RDS={rds}.{RET}"))
      readr::write_rds(sset, rds, compress="gz")
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
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
    else if (method=='open') sset <- sset %>% sesame::pOOBAH(force=force) %>% sesame::noob() %>% sesame::dyeBiasCorrTypeINorm()
    else if (method=='dyeBiasCorrTypeINorm') sset <- sset %>% sesame::dyeBiasCorrTypeINorm()
    else if (method=='detectionPnegEcdf') sset <- sset %>% sesame::detectionPnegEcdf(force=force)
    else if (method=='pOOBAH') sset <- sset %>% sesame::pOOBAH(force=force)
    else if (method=='noob') sset <- sset %>% sesame::noob()
    else if (method=='noobsb') sset <- sset %>% sesame::noobsb()
    else if (method=='inferTypeIChannel') sset <- sset %>% sesame::inferTypeIChannel(switch_failed=FALSE, verbose=FALSE)
    else if (method=='betas') sesame::extra(sset)[[method]] <- getBetas2(sset=sset, mask=quality.mask,sum.TypeI=sum.TypeI)
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
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  sset
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Extracted Sesame Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getBetas2 = function (sset, mask = FALSE, sum.TypeI = FALSE) 
{
  if (is(sset, "SigSetList")) {
    return(do.call(cbind, lapply(sset, getBetas2, mask = mask, 
                                 sum.TypeI = sum.TypeI)))
  }
  stopifnot(is(sset, "SigSet"))
  IGs <- IG(sset)
  IRs <- IR(sset)
  if (sum.TypeI) {
    IGs <- IGs + oobR2(sset)
    IRs <- IRs + oobG2(sset)
  }
  else if (!is.null(sset@extra$IGG) && !is.null(sset@extra$IRR)) {
    IGs[!sset@extra$IGG, ] <- sset@oobR[!sset@extra$IGG, ]
    IRs[!sset@extra$IRR, ] <- sset@oobG[!sset@extra$IRR, ]
  }
  betas <- c(pmax(IGs[, "M"], 1)/pmax(IGs[, "M"] + IGs[, "U"], 2), 
             pmax(IRs[, "M"], 1)/pmax(IRs[, "M"] + IRs[, "U"], 2), 
             pmax(II(sset)[, "M"], 1)/pmax(II(sset)[, "M"] + II(sset)[, "U"], 2))
  if (mask) 
    betas[!is.na(match(names(betas), sset@extra$mask))] <- NA
  betas
}

pOOBAH2 = function (sset, force = FALSE) 
{
  stopifnot(is(sset, "SigSet"))
  method <- "pOOBAH"
  if (!force && method %in% names(extra(sset)$pvals)) {
    cat("Retuning original sset\n")
    return(sset)
  }
  funcG <- ecdf(oobG(sset))
  funcR <- ecdf(oobR(sset))
  # funcR <- ecdf(oobR(sset) %>% head())
  
  pIR <- 1 - apply(cbind(funcR(IR(sset)[, "M"]), funcR(IR(sset)[, 
                                                                "U"])), 1, max)
  pIG <- 1 - apply(cbind(funcG(IG(sset)[, "M"]), funcG(IG(sset)[, 
                                                                "U"])), 1, max)
  pII <- 1 - apply(cbind(funcG(II(sset)[, "M"]), funcR(II(sset)[, 
                                                                "U"])), 1, max)
  names(pIR) <- rownames(IR(sset))
  names(pIG) <- rownames(IG(sset))
  names(pII) <- rownames(II(sset))
  if (!("pvals" %in% names(extra(sset)))) {
    cat("Creating fresh list\n")
    extra(sset)[["pvals"]] <- list()
  }
  cat(glue::glue("Setting method={method}{RET}{RET}"))
  pIR %>% head() %>% print()
  IR(sset) %>% head() %>% print()
  # oobR(sset) %>% head() %>% print()
  print(funcR)
  
  extra(sset)[["pvals"]][[method]] <- c(pIR, pIG, pII)
  sset
}


# End of file
