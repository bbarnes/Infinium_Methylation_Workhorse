
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


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#      Predefined SSET to Calls by Order of Operations Workflows::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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

# Predefined SSET to Calls by Order of Operations Workflows::
mutateSSET_workflow = function(sset, workflow, save=FALSE, rds=NULL, 
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutateSSET_workflow'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}.{RET}"))
  
  stime <- system.time({
    
    #
    # TBD:: This could be parsed into individual letters and then processed in that order...
    #
    if (workflow=='r') {
    } else if (workflow=='i') {
      sset <- sset %>% 
        mutateSesame(method = 'inferTypeIChannel', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='n') {
      sset <- sset %>% 
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='d') {
      sset <- sset %>% 
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='nd') {
      sset <- sset %>% 
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='dn') {
      sset <- sset %>% 
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='ind') {
      sset <- sset %>% 
        mutateSesame(method = 'inferTypeIChannel', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='ndi') {
      sset <- sset %>% 
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'inferTypeIChannel', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='din') {
      sset <- sset %>% 
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'inferTypeIChannel', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='dni') {
      sset <- sset %>% 
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'inferTypeIChannel', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'noob', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
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
  
  return(sset)
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
#                    Sesame SSET Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sigsSumToSSheet = function(tib, metric='avg',
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sigsSumToSSheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
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
  
  # ss_tib <- NULL
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
    # cat(glue::glue("2nd: pt={pt}.{RET}"))
  }
  ss_tib <- ss_tib %>% dplyr::mutate_if(is.numeric, as.integer)
  if (verbose>=vt+4) head(ss_tib) %>% as.data.frame() %>% print()
  
  ss_tib
}

sigTibToSSheet = function(sigs, man=NULL, by="Probe_ID", type="Probe_Type", des="Probe_Design", 
                          percision=-1, sort=FALSE, save=FALSE, csv=NULL,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sigTibToSSheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  tib <- NULL
  stime <- system.time({
    des <- des %>% rlang::sym()
    
    # tib <- sigs %>% dplyr::select(-all_of(by)) %>% dplyr::group_by(Probe_Design) %>% summarise_if(is.numeric, list(mu=mean), na.rm=TRUE)
    if (!is.null(man)) {
      type <- type %>% rlang::sym()
      sigs <- dplyr::select(man, !!by, !!type) %>% dplyr::right_join(sigs, by=by ) %>% dplyr::group_by(!!type, !!des) 
      tib <- sigs %>% summarise_if(is.numeric, list(min=min, med=median, avg=mean, sd=sd, max=max), na.rm=TRUE)
    } else {
      stop(glue::glue("\n[{func}]: ERROR: Currently only supported with the manifest!!!{RET}{RET}"))
      tib <- sigs %>% dplyr::group_by(!!des) %>% summarise_if(is.numeric, list(mu=mean), na.rm=TRUE)
    }
    
    if (save && !is.null(csv)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing signal set (percision={percision}) CSV={csv}.{RET}"))
      readr::write_csv(tib, csv)
    }
  })
  nrows <- tib %>% base::nrow()
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done nrows={nrows}, elapsed={etime}.{RET}{RET}"))
  
  tib
}

sset2tib = function(sset, man=NULL, by="Probe_ID", type="Probe_Type", des="Probe_Design", 
                    percision=-1, sort=FALSE, save=FALSE, csv=NULL,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sset2tib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  tib <- NULL
  stime <- system.time({
    
    des <- des %>% rlang::sym()
    
    tib <- dplyr::bind_rows(
      sset@IG   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='IG'),
      sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='OG'),
      
      sset@IR   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='IR'),
      sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='OR'),
      
      sset@II   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='II')
    )
    tib <- tib %>% dplyr::select(!!by, !!des, everything())
    
    if (percision!=-1) tib <- tib %>% dplyr::mutate_if(is.numeric, list(round), percision)
    if (!is.null(man)) tib <- dplyr::select(man, !!by, !!type) %>% dplyr::right_join(tib, by=by )
    
    by  <- by %>% rlang::sym()
    if (sort) tib <- tib %>% dplyr::arrange(!!by)
    
    if (!is.null(save) && save==TRUE && !is.null(csv)) {
      csv_dir <- base::basename(csv)
      if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
      
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing signal set (percision={percision}) CSV={csv}.{RET}"))
      readr::write_csv(tib, csv)
    }
  })
  nrows <- tib %>% base::nrow()
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done nrows={nrows}, elapsed={etime}.{RET}{RET}"))
  
  tib
}

sset2calls = function(sset, workflow, fresh=FALSE,
                      quality.mask = FALSE, sum.TypeI = FALSE,
                      percisionBeta=0, percisionPval=0,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sset2calls'
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
    if (fresh || is.null(sesame::extra(sset)[['betas']]) ) {
      # sesame::extra(sset)[['betas']] <-
      #   sesame::getBetas(sset=sset, mask = quality.mask,sum.TypeI = sum.TypeI)
      sesame::extra(sset)[['betas']] <-
        getBetas2(sset=sset, mask = quality.mask,sum.TypeI = sum.TypeI)
    }
    ret_tib <- tibble::enframe(sesame::extra(sset)[['betas']], 
                               name='Probe_ID', value=name)
    if (percisionBeta!=0) 
      ret_tib <- dplyr::mutate_if(ret_tib, purrr::is_double,round,percisionBeta)
    
    if (verbose>=vt+4) head(ret_tib) %>% print()
    
    #
    # Detection p-values::
    #
    pval_names <- sesame::extra(sset)[['pvals']] %>% names()
    for (pval_name in pval_names) {
      out_name <- pval_name
      if (pval_name=='pOOBAH') out_name <- 'poob'
      if (pval_name=='PnegEcdf') out_name <- 'negs'
      out_name <- paste(workflow,out_name, sep='_')
      if (verbose>=vt+4) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pval_name={pval_name}, out_name={out_name}...{RET}"))
      
       cur_tib <- sesame::extra(sset)[['pvals']][[pval_name]] %>% 
        tibble::enframe(name='Probe_ID', value=out_name)
       
       if (percisionPval!=0) 
         cur_tib <- dplyr::mutate_if(cur_tib, purrr::is_double,round,percisionPval)
       
      ret_tib <- ret_tib %>% dplyr::left_join(cur_tib,by="Probe_ID")
    }
    if (verbose>=vt+4) head(ret_tib) %>% print()
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done ret_cnt={ret_cnt}, elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame SSET To Inference Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: If manifest is provided then seperate by Probe_Type and Probe_Design
callToSSheet = function(call, idx, key, pre=NULL, minNegPval, minOobPval, 
                        percisionBeta=0, percisionPval=0, del='_', onlyCG=TRUE, id='Probe_ID',
                        verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'callToSSheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; idx={idx}, key={key}, del={del}...{RET}"))
  
  tib <- NULL
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
    beta_tib <- tibble::tibble(!!beta_key_sym := key, !!beta_val_sym := dplyr::select(call, ends_with('_beta')) %>% 
                                 dplyr::summarise_all(list(mean=mean), na.rm=TRUE) %>% 
                                 dplyr::pull() %>% round(percisionBeta) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} beta_key_str={beta_key_str}, beta_val_str={beta_val_str}.{RET}"))
    
    negs_tib <- NULL
    negs_key_str <- paste('Negs_Pass',idx,'Method', sep=del)
    negs_val_str <- paste('Negs_Pass',idx,'Perc', sep=del)
    negs_key_sym <- rlang::sym(negs_key_str)
    negs_val_sym <- rlang::sym(negs_val_str)
    
    print(call)
    
    negs_tib <- tibble::tibble(!!negs_key_sym := key, !!negs_val_sym := dplyr::select(call, ends_with('_negs')) %>%
                                 dplyr::summarise_all(list(pass_perc=cntPer_lte), min=minNegPval) %>%
                                 dplyr::pull() %>% round(percisionPval) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} negs_key_str={negs_key_str}, negs_val_str={negs_val_str}.{RET}"))
    
    poob_tib <- NULL
    poob_key_str <- paste('Poob_Pass',idx,'Method', sep=del)
    poob_val_str <- paste('Poob_Pass',idx,'Perc', sep=del)
    poob_key_sym <- rlang::sym(poob_key_str)
    poob_val_sym <- rlang::sym(poob_val_str)
    
    poob_tib <- tibble::tibble(!!poob_key_sym := key, !!poob_val_sym := dplyr::select(call, ends_with('_poob')) %>% 
                                 dplyr::summarise_all(list(pass_perc=cntPer_lte), min=minOobPval) %>%
                                 dplyr::pull() %>% round(percisionPval) )
    
    tib <- dplyr::bind_cols(pre, beta_tib, negs_tib, poob_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

ssetToInferences = function(sset, idx, key, pre=NULL, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToInferences'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  tib <- NULL
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
    
    tib <- dplyr::bind_cols(pre, gct_tib, sex_tib, kar_tib, eth_tib, igr_tib, irg_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

ssetToPredict = function(sset, idx, key, pre=NULL, del='_', fresh=FALSE,
                         quality.mask = FALSE, 
                         sum.TypeI = FALSE,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToPredict'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}      quality.mask={quality.mask}.{RET}"))
  # if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} nondetection.mask={nondetection.mask}.{RET}"))
  # if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}    correct.switch={correct.switch}.{RET}"))
  # if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}     mask.use.tcga={mask.use.tcga}.{RET}"))
  # if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}    pval.threshold={pval.threshold}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}         sum.TypeI={sum.TypeI}.{RET}"))
  # if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}        as.enframe={as.enframe}.{RET}"))
  # if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}         percision={percision}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}              sset={RET}"))
  if (verbose>=vt+4) print(sset)
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{RET}{RET}"))
  
  tib <- NULL
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
    
    tib <- dplyr::bind_cols(pre, skn_tib, phn_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
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

safeGCT = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeGCT'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::bisConversionControl(sset) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

safeEthnicity = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeEthnicity'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::inferEthnicity(sset) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

safeSexKaryo = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSexKaryo'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      # suppressWarnings(inferSexKaryotypes_Copy(sset) )
      suppressWarnings(sesame::inferSexKaryotypes(sset) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

safeSex = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSex'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      # suppressWarnings(inferSex_Copy(sset) )
      suppressWarnings(sesame::inferSex(sset) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame Beta To Predictions Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

safeSkinAge = function(beta, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSkinAge'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::predictAgeSkinBlood(betas=beta) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

safePhenoAge = function(beta, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safePhenoAge'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::predictAgePheno(betas=beta) )
    }, warning = function(w) {
      try_str <- paste('warning',funcTag, sep='-')
      NA
    }, error = function(e) {
      try_str <- paste('error',funcTag, sep='-')
      NA
    }, finally = {
      try_str <- paste('cleanup',funcTag, sep='-')
      NA
    })
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Sesame SSET to Tibs Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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
    
    # ret_tib <- sset@pval[[method]] %>% tibble::enframe(name='Probe_ID', value=name)
    ret_tib <- sset@pval[method] %>% tibble::enframe(name='Probe_ID', value=name)
    if (percision!=0) ret_tib <- ret_tib %>% dplyr::mutate_if(purrr::is_double, round, percision)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

ssetToSigsTib = function(sset, add, name, del='_', rmAdd=TRUE, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToSigsTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, name={name}!{RET}"))
  
  stime <- system.time({
    dat <- NULL
    
    datTag <- 'sig'
    grnTag <- paste(name,'G', sep=del) %>% rlang::sym()
    redTag <- paste(name,'R', sep=del) %>% rlang::sym()
    inbTag <- paste(name,'Inb_Col', sep=del) %>% rlang::sym()
    swpTag <- paste(name,'Swap', sep=del) %>% rlang::sym()
    # grnTag <- 'G' %>% rlang::sym()
    # redTag <- 'R' %>% rlang::sym()
    
    IG_tib <- sset@IG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='G')
    # dplyr::mutate(Design_Type=paste0(Design_Type,'G'))
    
    OR_tib <- sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='G')
    # dplyr::mutate(Design_Type=paste0(Design_Type,'G'))
    
    IR_tib <- sset@IR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='R')
    # dplyr::mutate(Design_Type=paste0(Design_Type,'R'))
    
    OG_tib <- sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='R')
    # dplyr::mutate(Design_Type=paste0(Design_Type,'R'))
    
    I1 <- dplyr::bind_rows(dplyr::inner_join(IG_tib,OR_tib, by=c('Probe_ID', 'Inb_Col', 'Design_Type')),
                           dplyr::inner_join(OG_tib,IR_tib, by=c('Probe_ID', 'Inb_Col', 'Design_Type')) ) %>%
      dplyr::mutate(Design_Type=paste0(Design_Type,'I')) %>%
      dplyr::select(Probe_ID,Design_Type,Inb_Col,everything())
    
    I1 <- add %>% dplyr::filter(Design_Type!='II') %>% 
      dplyr::left_join(I1, by=c("Probe_ID", "Design_Type") ) %>%
      dplyr::mutate(Swap=case_when(Man_Col!=Inb_Col ~ TRUE,
                                   Man_Col==Inb_Col ~ FALSE,
                                   TRUE ~ NA)) %>%
      dplyr::select(-Inb_Col) %>%
      dplyr::select(Probe_ID, Man_Col, Design_Type, Probe_Type, Swap, everything())
    # dplyr::select(Probe_ID, Address, Man_Col, Design_Type, Probe_Type, Inb_Col, Swap, everything())
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Infinium II::
    I2 <- sset@II %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      dplyr::rename(!!redTag :=U, !!grnTag :=M) %>%
      dplyr::mutate(Inb_Col=NA, Swap=NA, Design_Type='II')
    
    I2 <- add %>% dplyr::filter(Design_Type=='II') %>% 
      dplyr::left_join(I2, by=c("Probe_ID", "Design_Type") ) %>%
      dplyr::select(-Inb_Col) %>%
      dplyr::select(Probe_ID, Man_Col, Design_Type, Probe_Type, Swap, everything())
    # dplyr::select(Probe_ID, Address, Man_Col, Design_Type, Probe_Type, Inb_Col, Swap, everything())
    
    dat <- dplyr::bind_rows(I1, I2) %>%
      dplyr::rename(!!swpTag := Swap) %>%
      dplyr::arrange(Probe_ID)
    
    if (rmAdd && "Address" %in% names(dat)) dat <- dat %>% dplyr::select(-Address)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Sesame Tibs Summary Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sigsToSummary = function(tib, name, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sigsToSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    dat <- NULL
    
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
      summarise(!!cntTag :=n(), !!scntTag := count(!!swapTag==TRUE))
    sums <- tib %>% dplyr::group_by(Design_Type, Probe_Type, Man_Col) %>%
      summarise_if(is.numeric, list(min=min, avg=mean, med=median, max=max), na.rm=TRUE )
    
    dat <- cnts %>% dplyr::full_join(sums, by=c("Design_Type","Probe_Type", "Man_Col")) %>%
      dplyr::arrange(Probe_Type)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Sesame Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
getSesameManifest = function(man, sig, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getSesameManifest'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    tib <- man %>% dplyr::right_join(sig, by=c("U"="Address")) %>% 
      dplyr::filter(!is.na(Probe_ID)) %>%
      dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base) %>%
      dplyr::arrange(Probe_ID)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
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
