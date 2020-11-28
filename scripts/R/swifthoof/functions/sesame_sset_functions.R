
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
#                      Sesame SSET Access Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Sesame SSET Prediction Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Sesame Inference Methods:: SSET
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

safeGCT = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeGCT'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({

    try_str <- ''
    ret_val = tryCatch({
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
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))

  ret_val
}

safeEthnicity = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeEthnicity'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({

    try_str <- ''
    ret_val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(inferEthnicity2(sset) )
      # suppressWarnings(sesame::inferEthnicity(sset) )
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
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

safeSexKaryo = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSexKaryo'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({

    try_str <- ''
    ret_val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(inferSexKaryotypes2(sset) )
      # suppressWarnings(sesame::inferSexKaryotypes(sset) )
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
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

safeSex = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSex'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({

    try_str <- ''
    ret_val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(inferSex2(sset) )
      # suppressWarnings(sesame::inferSex(sset) )
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
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame Predictions Methods:: Beta
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

safeSkinAge = function(beta, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSkinAge'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({

    try_str <- ''
    ret_val = tryCatch({
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
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

safePhenoAge = function(beta, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safePhenoAge'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- NA
  stime <- system.time({

    try_str <- ''
    ret_val = tryCatch({
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
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Val={ret_val}; try_str={try_str}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Sesame SSET Sample Sheet Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sigsSumToSSheet = function(tib, metric='mean',
                           by="Probe_ID", type="Probe_Type", des="Probe_Design", 
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sigsSumToSSheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting (metric={metric})...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    by_sym  <- by %>% rlang::sym()
    des_sym  <- des %>% rlang::sym()
    type_sym <- type %>% rlang::sym()

    ret_tib <- tib %>% 
      tidyr::unite(!!type_sym, !!type_sym,!!des_sym, sep='_') %>%
      dplyr::select(!!type_sym, ends_with(!!metric)) %>% 
      tidyr::gather(Metric, Value, -!!type_sym) %>%
      tidyr::unite(Key, !!type_sym,Metric, sep='_') %>%
      tidyr::spread(Key, Value)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame Tibs Summary Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sigsTibToSummary = function(tib, man=NULL,
                           cutoff=NULL, percision=-1,
                           save=FALSE, csv=NULL,
                           by="Probe_ID", type="Probe_Type", des="Probe_Design",
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sigsTibToSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; by={by}, type={type}, des={des}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    by_sym  <- by %>% rlang::sym()
    des_sym  <- des %>% rlang::sym()
    type_sym <- type %>% rlang::sym()
    
    if (!is.null(man)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Will use provided manifest...{RET}"))
      if (verbose>=vt+3) man %>% print()
      
      tib_cols <- names(tib)
      if (!(des %in% tib_cols)) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} right_join3(by,type,des)...{RET}"))
        
        tib <- dplyr::select(man, !!by, !!type, !!des) %>% dplyr::right_join(tib, by=by)
      } else {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} right_join2(by,type)...{RET}"))
        
        tib <- dplyr::select(man, !!by, !!type) %>% dplyr::right_join(tib, by=by)
      }

    } else {
      if (verbose>=vt+1) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} No manifest provided...{RET}"))
    }
    tib_cols <- names(tib)
    if (!(type %in% tib_cols)) 
      tib <- tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(!!by_sym, 1,2))
    if (verbose>=vt+4) print(tib)
    
    # Check for all fields::
    #
    tib_cols <- names(tib)
    if (!(by %in% tib_cols)) {
      stop(glue::glue("{RET}[{funcTag}]:{tabsStr}{TAB} ERROR: Failed to find by={by}!!!{RET}{RET}"))
      return(ret_tib)
    }
    if (!(des %in% tib_cols)) {
      stop(glue::glue("{RET}[{funcTag}]:{tabsStr}{TAB} ERROR: Failed to find des={des}!!!{RET}{RET}"))
      return(ret_tib)
    }
    if (!(type %in% tib_cols)) {
      stop(glue::glue("{RET}[{funcTag}]:{tabsStr}{TAB} ERROR: Failed to find type={type}!!!{RET}{RET}"))
      return(ret_tib)
    }

    if (is.null(cutoff)) {
      ret_tib <- tib %>%
        dplyr::group_by(!!type_sym, !!des_sym) %>%
        summarise_if(is.numeric, list(min=min, median=median, mean=mean, sd=sd, max=max), na.rm=TRUE)
    } else {
      ret_tib <- tib %>% 
        dplyr::group_by(!!type_sym, !!des_sym) %>%
        dplyr::select(-dplyr::all_of(by)) %>%
        dplyr::summarise_all(list(pass_perc=cntPer_lte), min=cutoff)
    }
    ret_tib <- ret_tib %>% dplyr::ungroup()

    if (percision >= 0) 
      ret_tib <- dplyr::mutate_if(ret_tib, is.numeric, list(round), percision)
    
    if (!is.null(save) && save==TRUE && !is.null(csv)) {
      csv_dir <- base::basename(csv)
      if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
      csv <- clean_file(csv, verbose=verbose,vt=vt+3,tc=tc+1,tt=tt)

      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing (percision={percision}) CSV={csv}.{RET}"))
      if (verbose>=vt+4) ret_tib %>% print()
      readr::write_csv(ret_tib, csv)
    }
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Sesame SSET Transform to Tib Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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
      
      ret_tib <- ret_tib %>% dplyr::left_join(cur_tib,by="Probe_ID", copy=TRUE)
      # ret_tib <- ret_tib %>% dplyr::inner_join(cur_tib,by="Probe_ID", copy=TRUE)
    }
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Current Table(pvals)={RET}"))
      ret_tib %>% print()
    }
    
    if (!is.null(save) && save==TRUE && !is.null(csv)) {
      csv_dir <- base::basename(csv)
      if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
      csv <- clean_file(csv, verbose=verbose,vt=vt+3,tc=tc+1,tt=tt)
      
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing CSV={csv}.{RET}"))
      if (verbose>=vt+4) ret_tib %>% print()
      readr::write_csv(ret_tib, csv)
    }
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
      csv_dir <- base::basename(csv)
      if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive=TRUE)
      csv <- clean_file(csv, verbose=verbose,vt=vt+3,tc=tc+1,tt=tt)
      
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing (percision={percision}) CSV={csv}.{RET}"))
      if (verbose>=vt+4) ret_tib %>% print()
      readr::write_csv(ret_tib, csv)
    }
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame SSET Mutation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

mutateSset = function(sset, method, force=TRUE,
                        quality.mask=FALSE, sum.TypeI=FALSE, switch_failed=FALSE,
                        verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutateSset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; Mutate Sesame({method})...{RET}"))
  
  ctl_cnt <- sset@ctl %>% base::nrow()
  if (ctl_cnt==0 && method=='detectionPnegEcdf') {
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Unable to mutate ctl_cnt={ctl_cnt} ",
                     "for method={method}. Returning original sset...{RET}"))
    return(sset)
  }
  
  ret_cnt <- 0
  stime <- system.time({
    if (verbose>=vt+2) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Pre-sset(ctl={ctl_cnt})={RET}"))
      print(sset)
    }
    
    if (is.null(method)) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Missing method!!!{RET}{RET}"))
    } else if (method=='open') {
      sset <- sset %>% sesame::pOOBAH(force=force) %>% 
        sesame::noob() %>% sesame::dyeBiasCorrTypeINorm()
    } else if (method=='dyeBiasCorrTypeINorm') {
      sset <- sset %>% sesame::dyeBiasCorrTypeINorm()
    } else if (method=='detectionPnegEcdf') {
      sset <- sset %>% sesame::detectionPnegEcdf(force=force)
    } else if (method=='pOOBAH') {
      sset <- sset %>% sesame::pOOBAH(force=force)
    } else if (method=='noob') {
      sset <- sset %>% sesame::noob()
    } else if (method=='noobsb') {
      sset <- sset %>% sesame::noobsb()
    } else if (method=='inferTypeIChannel') {
      sset <- sset %>% 
        sesame::inferTypeIChannel(switch_failed=FALSE, verbose=FALSE)
    } else if (method=='betas') {
      sesame::extra(sset)[[method]] <- NULL
      sesame::extra(sset)[[method]] <- 
        getBetas2(sset=sset, mask=quality.mask,sum.TypeI=sum.TypeI)
    } else if (method=='raw') {
      # sset <- sset
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported method={method}!!!{RET}{RET}"))
    }
    
    ret_cnt <- sset %>% slotNames() %>% length()
    if (verbose>=vt+2) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Post-sset(slots={ret_cnt})={RET}"))
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Sesame SSET Initialization Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

newSset = function(prefix, platform, manifest, 
                   load=FALSE, save=FALSE, rds=NULL,
                   verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'newSset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; platform={platform}, prefix={prefix}.{RET}"))
  
  ret_cnt <- 0
  stime <- system.time({
    
    if (load && !is.null(rds) && file.exists(rds)) {
      if (verbose>=vt+1) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading RDS={rds}.{RET}"))
      sset <- readr::read_rds(rds)
    } else {
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} manifest={RET}"))
        manifest %>% print()
      }
      
      sset <- sesame::readIDATpair(prefix, platform=platform, manifest=manifest)
      if (save && !is.null(rds)) {
        if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing RDS={rds}.{RET}"))
        if (verbose>=vt+4) sset %>% print()
        readr::write_rds(sset, rds, compress="gz")
      }
    }
    
    ret_cnt <- sset %>% slotNames() %>% length()
    if (verbose>=vt+5) {
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Extracted Sesame Functions:: Betas
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getBetas2 = function (sset, mask = FALSE, sum.TypeI = FALSE, extra=TRUE) 
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
  else if (extra && !is.null(sset@extra$IGG) && !is.null(sset@extra$IRR)) {
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
    # cat("Retuning original sset\n")
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
    # cat("Creating fresh list\n")
    extra(sset)[["pvals"]] <- list()
  }
  # cat(glue::glue("Setting method={method}{RET}{RET}"))
  pIR %>% head() %>% print()
  IR(sset) %>% head() %>% print()
  # oobR(sset) %>% head() %>% print()
  print(funcR)
  
  extra(sset)[["pvals"]][[method]] <- c(pIR, pIG, pII)
  sset
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Extracted Sesame Functions:: Sex
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

inferSexKaryotypes2 = function (sset) 
{
  stopifnot(is(sset, "SigSet"))
  sex.info <- getSexInfo2(sset)
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

inferSex2 = function (sset) 
{
  stopifnot(is(sset, "SigSet"))
  sex.info <- getSexInfo2(sset)[seq_len(3)]
  as.character(predict(sesameDataGet("sex.inference"), sex.info))
}

getSexInfo2 = function (sset) 
{
  if (is(sset, "SigSetList")) 
    return(do.call(cbind, lapply(sset, getSexInfo)))
  stopifnot(is(sset, "SigSet"))
  cleanY <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$chrY.clean
  # cat("length(cleanY)=\n")
  # cleanY %>% length() %>% print()
  # cleanY %>% head() %>% print()
  
  xLinked <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$chrX.xlinked
  # cat("length(xLinked)=\n")
  # xLinked %>% length() %>% print()
  # xLinked %>% head() %>% print()
  # cat("\n\n\n")
  
  probe2chr <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$probe2chr.hg19
  # cat("length(probe2chr)=\n")
  # probe2chr %>% length() %>% print()
  # probe2chr %>% head() %>% print()
  # cat("\n\n\n")
  
  # xLinkedBeta <- sesame::getBetas(sset=sesame::subsetSignal(sset, xLinked))
  xLinkedBeta <- getBetas2(sset=sesame::subsetSignal(sset, xLinked), extra=FALSE)
  # cat("head(xLinkedBeta)=\n")
  # xLinkedBeta %>% length() %>% print()
  # xLinkedBeta %>% head() %>% print()
  
  intens <- sesame::totalIntensities(sset)
  probes <- intersect(names(intens), names(probe2chr))
  intens <- intens[probes]
  probe2chr <- probe2chr[probes]
  # print(probe2chr)
  # return( sesame::subsetSignal(sset, cleanY) )
  # return( median(sesame::totalIntensities(sesame::subsetSignal(sset, cleanY))) )
  # cat("\n\n\n")
  
  c(medianY = median(sesame::totalIntensities(sesame::subsetSignal(sset, cleanY)), na.rm=TRUE), 
    medianX = median(sesame::totalIntensities(sesame::subsetSignal(sset, xLinked)), na.rm=TRUE), fracXlinked = 
      (sum(xLinkedBeta > 0.3 & xLinkedBeta < 0.7, na.rm = TRUE)/sum(!(is.na(xLinkedBeta))) ), 
    fracXmeth = (sum(xLinkedBeta > 0.7, na.rm = TRUE)/sum(!(is.na(xLinkedBeta)))), 
    fracXunmeth = (sum(xLinkedBeta < 0.3, na.rm = TRUE)/sum(!(is.na(xLinkedBeta)))), 
    tapply(intens, probe2chr, median, na.rm=TRUE))
}

inferEthnicity2 = function (sset) 
{
  if (is(sset, "SigSetList")) 
    return(vapply(sset, inferEthnicity2, character(1)))
  stopifnot(is(sset, "SigSet"))
  stopifnot(sset@platform %in% c("EPIC", "HM450"))
  ethnicity.inference <- sesameDataGet("ethnicity.inference")
  ccsprobes <- ethnicity.inference$ccs.probes
  rsprobes <- ethnicity.inference$rs.probes
  ethnicity.model <- ethnicity.inference$model
  af <- c(getBetas2(sset, mask = FALSE)[rsprobes], 
          getAFTypeIbySumAlleles(subsetSignal(sset,ccsprobes)))
  as.character(predict(ethnicity.model, af))
}

# End of file
