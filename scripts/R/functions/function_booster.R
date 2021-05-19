
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'template_func'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # ret_cnt <- ret_tib %>% base::nrow()
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Common Booster Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

clean_tibble = function(tib,
                        verbose=0,vt=6,tc=1,tt=NULL) {
  funcTag <- 'clean_tibble'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- tib %>%
      select(where(~sum(!is.na(.x)) > 0)) %>%
      utils::type.convert() %>%
      dplyr::mutate(across(where(is.factor), as.character) )
    
    # ret_cnt <- ret_tib %>% base::nrow()
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

makeFieldUnique = function(tib, field, add=NULL,
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'makeFieldUnique'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; field={field}...{RET}"))
  
  field_sym <- rlang::sym(field)
  outkey_sym <- rlang::sym(field)
  if (!is.null(add)) outkey_sym <- rlang::sym(add)
  
  tot_cnt <- tib %>% base::nrow()
  unq_cnt <- tib %>% dplyr::distinct(!!field_sym) %>% base::nrow()
  # unq_cnt <- tib %>% dplyr::distinct(!!field_sym) %>% base::nrow()
  
  if (tot_cnt==unq_cnt) {
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Field({field}) already unique!",
                     "Total={tot_cnt}, Unique={unq_cnt}.{RET}"))
    return(tib)
  }
  stopifnot(tot_cnt>unq_cnt)
  
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Field({field}): ",
                   "Total={tot_cnt}, Unique={unq_cnt}.{RET}"))
  
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # Force ID's to be unique::
    ret_tib <- tib %>% 
      dplyr::group_by(!!field_sym) %>%
      dplyr::mutate(
        !!outkey_sym := paste(!!field_sym,row_number(), sep='')) %>% 
      ungroup()
    
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

get_fileSuffix = function(file,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'get_fileSuffix'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_val <- NULL
  stime <- system.time({
    
    ret_val <- stringr::str_remove(file, "\\.gz$") %>% stringr::str_remove("^.*\\.")
    ret_cnt <- ret_val %>% length()
    # ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Clean Files for docker::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

clean_file = function(file,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'clean_file'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    cmd <- NULL
    if (file.exists(file)) unlink(file)
    cmd <- glue::glue("touch {file}")
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Running cmd='{cmd}'...{RET}"))
    base::system(cmd)
    
    if (!file.exists(file)) {
      stop(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: File doesn't exist; file={file}!!!{RET}"))
    }
    
    ret_cnt <- 1
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  file
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                General Program Initialization Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

program_init = function(name,defs=NULL,
                        opts,opt_reqs=NULL,pars,par_reqs=NULL,
                        libs=TRUE,rcpp=FALSE,
                        verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'program_init'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         Validate Required Options::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  if (!is.null(par_reqs)) {
    pars <- check_list(list=pars, vals=par_reqs, name='params',
                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  }
  if (!is.null(opt_reqs)) {
    opts <- check_list(list=opts, vals=opt_reqs, name='options',
                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  }
  
  stopifnot(!is.null(opts[['outDir']]))
  stopifnot(!is.null(pars[['prgmTag']]))
  stopifnot(!is.null(pars[['exePath']]))
  
  if (!is.null(opt[['verbose']])) verbose <- opts$verbose
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                              Load Libraries::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  if (libs) {
    opts <- load_libraries(opts=opts, pars=pars, rcpp=rcpp,
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  }
  # if (!is.null(defs))
  #   def_tib <- dplyr::bind_rows(defs) %>% tidyr::gather("Params", "Value")
  opt_tib <- dplyr::bind_rows(opts) %>% 
    tidyr::gather("Option", "Value") %>%
    dplyr::distinct()
  
  par_tib <- pars %>%
    unlist(recursive = TRUE) %>%
    dplyr::bind_rows() %>% 
    tidyr::gather("Params", "Value") %>%
    dplyr::distinct()
  if (opts$verbose>=1) par_tib %>% base::print(n=base::nrow(par_tib) )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Build Directories::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # opts$outDir <- file.path(opts$outDir, pars$prgmDir, name)
  opts$outDir <- file.path(opts$outDir, name)
  if (!is.null(opts$runName)) opts$outDir <- file.path(opts$outDir, opts$runName)
  if (!dir.exists(opts$outDir)) dir.create(opts$outDir, recursive=TRUE)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]: Output Directory (TOP)={opts$outDir}...{RET}"))
  
  if (!is.null(opts[['fresh']]) && opts$fresh)
    unlink(list.files(opts$outDir, full.names=TRUE))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         Program Start Time Stamp::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  opts$time_org_txt <- file.path(opts$outDir,paste(name,'time_stamp-org.txt', sep='.'))
  if (!file.exists(opts$time_org_txt) || (!is.null(opts[['fresh']]) && opts$fresh) )
    readr::write_lines(x=date(),file=opts$time_org_txt,sep='\n',append=FALSE)
  
  opts$opt_csv  <- file.path(opts$outDir, paste(pars$prgmTag,'program-options.csv', sep='.') )
  opts$par_csv  <- file.path(opts$outDir, paste(pars$prgmTag,'program-parameters.csv', sep='.') )
  opts$time_csv <- file.path(opts$outDir, paste(pars$prgmTag,'time-tracker.csv.gz', sep='.') )
  
  if (file.exists(opts$opt_csv))  unlink(opts$opt_csv)
  if (file.exists(opts$par_csv))  unlink(opts$par_csv)
  if (file.exists(opts$time_csv)) unlink(opts$time_csv)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                          Program Command Shell::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  pars$cmd_shell <- file.path(opts$outDir,paste(pars$prgmTag,name,'command.sh', sep='.'))
  pars$cmd_shell <- 
    optsToCommand(opts=opt_tib, pre=opts$Rscript,exe=pars$exePath, file=pars$cmd_shell, 
                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  
  ret_cnt <- names(opts) %>% length()
  
  # etime <- stime[3] %>% as.double() %>% round(2)
  etime <- 0
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  opts
}

program_done = function(opts,pars,
                        verbose=0,vt=3,tc=1,tt) {
  funcTag <- 'program_done'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  opts_tib  <- opts %>%
    unlist(recursive=TRUE) %>%
    dplyr::bind_rows() %>% 
    tidyr::gather("Option", "Value")
  
  pars_tib  <- pars %>%
    unlist(recursive=TRUE) %>%
    dplyr::bind_rows() %>% 
    tidyr::gather("Params", "Value")
  time_tib <- tt$time %>% dplyr::mutate_if(is.numeric, list(round), 4)
  
  readr::write_csv(opts_tib, opts$opt_csv)
  readr::write_csv(pars_tib, opts$par_csv)
  readr::write_csv(time_tib, opts$time_csv)
  
  sysTime <- Sys.time()
  cat(glue::glue("{RET}[{pars$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))
  
  return(0)
}

load_libraries = function(opts, pars, rcpp=FALSE,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_libraries'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  
  stopifnot(!is.null(pars$prgmDir))
  stopifnot(!is.null(pars$scrDir))
  
  # list.files("/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/scripts/R", full.names = TRUE) %>% file.path("functions")
  # list.files(pars$scrDir, full.names = TRUE) %>% file.path("functions")
  
  pars$prgm_src_dir <- file.path(pars$scrDir,pars$prgmDir, 'functions')
  if (!dir.exists(pars$prgm_src_dir)) stop(glue::glue("[{pars$prgmTag}]: Program Source={pars$prgm_src_dir} does not exist!{RET}"))
  for (sfile in list.files(path=pars$prgm_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
  cat(glue::glue("[{pars$prgmTag}]: Done. Loading Source Files form Program Source={pars$prgm_src_dir}!{RET}{RET}") )
  
  # Load All other function methods:: Search Method::
  src_files <- base::list.files(pars$scrDir, full.names = TRUE) %>% 
    file.path("functions") %>% base::list.files(pattern='.R$', full.names=TRUE)
  for (sfile in src_files ) {
    base::source(sfile)
    cat(glue::glue("[{pars$prgmTag}]: Done. Loading Source File={sfile}!{RET}") )
  }
  
  if (rcpp) {
    pars$sourceCpp <- file.path(pars$scrDir, 'Rcpp/cpgLociVariation.cpp')
    
    # Don't care about linux now that we use docker::
    # TBD:: Keeping this next line for now, but should be removed if not used
    opts <- setLaunchExe(opts=opts, pars=pars, verbose=opts$verbose, vt=5,tc=0)
    # if (!opts$isLinux) {
    if (file.exists(pars$sourceCpp)) {
      if (!file.exists(pars$sourceCpp)) pars$sourceCpp <- file.path(pars$scrDir, 'Rcpp/cpgLociVariation.cpp')
      if (!file.exists(pars$sourceCpp)) stop(glue::glue("[{pars$prgmTag}]: Source={pars$sourceCpp} does not exist!{RET}"))
      Rcpp::sourceCpp(pars$sourceCpp)
    }
  }
  ret_cnt <- opts %>% names() %>% length()
  
  # etime <- stime[3] %>% as.double() %>% round(2)
  etime <- 0
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  opts
}

check_list = function(list, vals, name="options",
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'check_list'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt  <- 0
  vals_len <- length(vals)
  valid_list <- TRUE
  for (val in vals) {
    if (is.null(list[[val]])) {
      cat(glue::glue("[Usage]: In {name} {val} is NULL [Required]!!!{RET}"))
      valid_list <- FALSE
    }
  }
  if (valid_list) {
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Validated all required options({vals_len}).{RET}"))
  } else {
    stop(glue::glue("[{funcTag}]:{tabsStr} Invalid Options; use --help!!!{RET}{RET}"))
    return(NULL)
  }
  
  ret_cnt <- names(list) %>% length()
  
  # etime <- stime[3] %>% as.double() %>% round(2)
  etime <- 0
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  list
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Basic Tibble Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Check if all files exist and have correct time stamp order
valid_time_stamp = function(files,
                            funcTag='valid_time_stamp',
                            verbose=0,vt=3,tc=1,tt=NULL) {
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_val <- TRUE
  stime <- system.time({
    
    file_cnt <- length(files)
    for (ii in c(1:file_cnt)) {
      if (!file.exists(files[ii])) ret_val <- FALSE
    }
    if (verbose>=vt+4)
      cat(glue::glue("[{funcTag}]:{tabsStr} Files exist={ret_val}.{RET}{RET}"))
    
    if (ret_val) {
      file_cnt <- file_cnt - 1
      for (ii in c(1:file_cnt)) {
        f1 <- files[ii]
        f2 <- files[ii+1]
        
        t1 <- file.mtime(f1)
        t2 <- file.mtime(f2)
        
        if (verbose>=vt+4)
          cat(glue::glue("ii={ii}, t1={t1}; f1={f1}{RET}",
                         "ii={ii}, t2={t2}; f2={f2}{RET}"))
        
        if (t1 > t2) ret_val <- FALSE
        
        if (verbose>=vt+4)
          cat(glue::glue("ii={ii}, Valid Order={ret_val}.{RET}{RET}"))
      }
    }
    ret_cnt <- files %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

# This function is out of date replacing with simpler function above::
#  valid_time_stamp()
check_timeStamps = function(name,outDir,origin=NULL,files,
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'check_timeStamps'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  beg_txt <- NULL
  end_txt <- NULL
  isValid <- TRUE
  stime <- system.time({
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
    beg_txt <- file.path(outDir,paste(name,'time-stamp-beg.txt', sep='.') )
    end_txt <- file.path(outDir,paste(name,'time-stamp-end.txt', sep='.') )
    
    ret_tib <- tibble::tibble(
      valid=isValid,
      beg=beg_txt,
      end=end_txt,
    )
    
    # Check everything to be checked exists::
    if (!is.null(origin) && !file.exists(origin)) isValid <- FALSE
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Origin exist={isValid}.{RET}"))
    
    if (!file.exists(beg_txt)) isValid <- FALSE
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Beg exist={isValid}.{RET}"))
    
    if (!file.exists(end_txt)) isValid <- FALSE
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} End exist={isValid}.{RET}"))
    
    for (file in files) {
      if (!file.exists(file)) isValid <- FALSE
    }
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} All files exist={isValid}.{RET}"))
    
    
    # Check all time stamps::
    if (isValid && !is.null(origin) && !file.exists(origin) &&
        file.mtime(orgin) > file.mtime(beg_txt)) isValid <- FALSE
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Origin stamp check={isValid}.{RET}"))
    
    if (isValid && file.mtime(beg_txt) > file.mtime(end_txt)) isValid <- FALSE
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Beg/End stamp check={isValid}.{RET}"))
    
    if (isValid) {
      for (file in files) {
        if (file.mtime(beg_txt) >= file.mtime(file) ||
            file.mtime(end_txt) <= file.mtime(file)) {
          if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Failed stamp={file}.{RET}"))
          isValid <- FALSE
          break
        }
      }
    }
    
    # Clear everything if validation failed for any reason::
    if (!isValid) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Cleaning all files...{RET}"))
      
      if (file.exists(beg_txt)) unlink(beg_txt)
      if (file.exists(end_txt)) unlink(end_txt)
      for (file in files) {
        if (!file.exists(file)) unlink(file)
      }
      readr::write_lines(x=date(),file=beg_txt,sep='\n',append=FALSE)
    }
    
    # Update return tibble::
    ret_tib$isValid <- isValid
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Basic Tibble Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

reduceSortedTib = function(tib, n=3,
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'reduceSortedTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; n={n}...{RET}"))
  
  tib_cnt <- tib %>% base::nrow()
  if (tib_cnt <= 1) return(tib)
  if (tib_cnt <= n) return(tib)
  
  mid_cnt <- as.integer( tib_cnt/n ) + 1
  sel_vec <- c(1,mid_cnt,tib_cnt) %>% unique()
  tib <- tib[sel_vec, ]
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  tib
}

joinTibbles = function(a, b, by, side="left", verbose=0,vt=3,tc=1) {
  funcTag <- 'joinTibbles'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; platform={platform}.{RET}"))
  
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  
  if (side=='left')  return(dplyr::left_join(a,  b, by=by))
  if (side=='right') return(dplyr::right_join(a, b, by=by))
  if (side=='inner') return(dplyr::inner_join(a, b, by=by))
  if (side=='full')  return(dplyr::full_join(a,  b, by=by))
  stop(glue::glue("[$func]: ERROR: Unsupported join={side}!!!{RET}{RET}"))
  
  # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  NULL
}

addColNames = function(tib, add, fix, prefix=TRUE, del='_',
                       verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'addColNames'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- tib %>% dplyr::select(dplyr::all_of(fix), dplyr::everything())
    new_col <- ret_tib %>% dplyr::select(-dplyr::all_of(fix)) %>% names()
    
    if ( prefix) new_col <- paste(add,new_col, sep=del)
    if (!prefix) new_col <- paste(new_col,add, sep=del)
    
    new_col <- c(fix, new_col)
    ret_tib <- ret_tib %>% purrr::set_names(new_col)
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
#                       Options to Script Commands::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

optsToCommand = function(opts, pre=NULL, exe, rm=NULL, add=NULL, 
                         file=NULL, key="Option", val="Value",
                         cluster=TRUE,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'optsToCommand'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val  <- NULL
  add_str  <- ''
  opt_str  <- ''
  bool_str <- ''
  key <- key %>% as.character() %>% rlang::sym()
  val <- val %>% as.character() %>% rlang::sym()
  
  # Options:: Exclude
  if (!is.null(rm)) opts <- opts %>% dplyr::filter(! (!!key %in% rm) )
  
  # Handel Boolean Variables::
  bool <- NULL
  bool <- opts %>% dplyr::filter(Value=="TRUE")
  opts <- opts %>% dplyr::filter(Value!="TRUE")
  opts <- opts %>% dplyr::filter(Value!="FALSE")
  
  if (!is.null(bool)) {
    bool_str <- bool %>% 
      dplyr::mutate(!!key := stringr::str_c('--',!!key)) %>%
      dplyr::pull(!!key) %>% stringr::str_c(collapse=" ")
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} bool_str='{bool_str}'.{RET}{RET}"))
  }
  
  # Merge Options::
  if (!is.null(opts)) {
    opt_str <- opts %>% dplyr::arrange(!!key) %>%
      tidyr::unite(Param, !!key, !!val, sep='=') %>%
      dplyr::mutate(Param=stringr::str_c('--',Param)) %>% 
      dplyr::pull(Param) %>%
      stringr::str_c(collapse=" ")
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} opt_str='{opt_str}'.{RET}{RET}"))
    
    # Second Removal Attempt of rm fields::
    if (!cluster) {
      opt_str <- opt_str %>% stringr::str_remove('--cluster')
      if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} opt_str(-cluster)='{opt_str}'.{RET}{RET}"))
    }
  }
  
  # Options:: Add
  if (!is.null(add)) {
    add_str <- add %>% tidyr::unite(Param, !!key, !!val, sep='=') %>%
      dplyr::mutate(Param=stringr::str_c('--',Param)) %>% 
      dplyr::pull(Param) %>%
      stringr::str_c(collapse=" ")
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} add_str='{add_str}'.{RET}{RET}"))
  }
  
  # Add Executable and Join Options::
  cmd <- ''
  if (!is.null(pre) && length(pre)!=0) cmd <- pre
  cmd <- stringr::str_c(cmd,exe,opt_str,add_str,bool_str, sep=' ')
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} cmd='{cmd}'.{RET}{RET}"))
  
  ret_val <- cmd
  if (!is.null(file)) {
    dir <- base::dirname(file)
    if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    base::unlink(file)
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing program shell={file}...{RET}"))
    readr::write_lines(cmd, file=file)
    Sys.chmod(file, mode="0777")
    
    ret_val <- file
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  ret_val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          File Searching Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

findFileByPattern = function(dir, pattern, max=1, recursive=FALSE, 
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'findFileByPattern'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; max={max}, pattern={pattern}, dir={dir}.{RET}"))
  
  files <- NULL
  files <- base::list.files(path=dir, pattern=pattern, full.names=TRUE, recursive=recursive)
  
  if (verbose>=vt) {
    cat(glue::glue("[{funcTag}]:{tabsStr} Files={RET}"))
    print(files)
  }
  if (is.null(files) || length(files)==0) {
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed to find ",
                    "pattern='{pattern}' in dir='{dir}'!!!{RET}{RET}"))
    return(NULL)
  }
  
  files_cnt <- length(files)
  if (files_cnt>max)
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Found {files_cnt} > max({max}) with pattern={pattern} in dir={dir}!!!{RET}{RET}"))
  
  files
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Linux Cluster Options Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

setLaunchExe = function(opts, pars, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'setLaunchExe'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  opts$lanExe <- ''
  opts$isLinux <- FALSE
  
  # if (is.null(pars$lixDir))
  if (is.null(pars$lixDir))
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Options lixDir=NULL not allowed!!!{RET}{RET}"))
  
  if (dir.exists(pars$lixDir)) {
    opts$isLinux <- TRUE
    opts$lanExe <- 'qsub -cwd -pe threaded 16 -l excl=true -N'
    if (dir.exists(pars$macDir)) stop(glue::glue("{RET}[{funcTag}]: Linux/Mac directories exist???{RET}{RET}"))
  }
  
  opts
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Summarizing Counting Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
cntPer_lte = function(x, min, prc=3) {
  min <- as.numeric(min)
  round(100*length(x[which(x<=min)])/length(x),prc)
}

cntPer_gte = function(x, min, prc=3) {
  min <- as.numeric(min)
  round(100*length(x[which(x>=min)])/length(x),prc)
}

cntPer_gt = function(x, max, prc=3) {
  max <- as.numeric(max)
  round(100*length(x[which(x> max)])/length(x),prc)
}

cntPer_gt0 = function(x, max, prc=3) {
  max <- as.numeric(max)
  round(length(x[which(x>max)])/length(x),prc)
}

bool2int = function(x) {
  y <- NULL
  for (ii in seq(length(x))) {
    if (is.na(x[ii]) || x[ii]==FALSE) {
      y[ii] <- 0
    } else {
      y[ii] <- 1
    }
  }
  
  y
}

# Guess file delimiter
guess_file_del = function(file, n_max=100,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'guess_file_del'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_val <- NULL
  stime <- system.time({
    
    dat_tib <- readr::read_lines(file, n_max=n_max) %>%
      tibble::as_tibble() %>% 
      dplyr::mutate(Row_Num=dplyr::row_number())
    
    del_tib <- dat_tib %>% 
      dplyr::mutate(
        com_cnt=stringr::str_count(value, COM), 
        tab_cnt=stringr::str_count(value, TAB),
        spc_cnt=stringr::str_count(value, " "))
    
    med_tib <- del_tib %>%
      dplyr::summarise(com=median(com_cnt,na.rm=TRUE),
                       tab=median(tab_cnt,na.rm=TRUE),
                       spc=median(spc_cnt,na.rm=TRUE),
      ) 
    med_key <- med_tib %>% which.max() %>% names()
    
    sum_tib <- del_tib %>%
      dplyr::summarise(com=sum(com_cnt,na.rm=TRUE),
                       tab=sum(tab_cnt,na.rm=TRUE),
                       spc=sum(spc_cnt,na.rm=TRUE),
      )
    sum_key <- sum_tib %>% which.max() %>% names()
    
    if (med_key==sum_key) {
      if (med_key=="com") ret_val=COM
      if (med_key=="tab") ret_val=TAB
      if (med_key=="spc") ret_val=" "
    }
    
    ret_tib <- dplyr::bind_rows(med_tib, sum_tib) %>%
      dplyr::mutate(ret_val=ret_val)
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt,tc)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              String Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

safe_write = function(x, type, file=NULL, 
                      cols=TRUE, append=FALSE, done=FALSE,
                      funcTag="safe_write",
                      verbose=0,vt=4,tc=1,tt=NULL) {
  # If file is provided, usually in the case where the user didn't want it written
  #  simply return invisible and do nothing...
  
  if (is.null(file)) return(base::invisible(NULL))
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting type={type}.{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    build_file_dir(file)
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Writing Data {type}={file}...{RET}"))
    
    if (purrr::is_vector(x)) ret_cnt <- length(x)
    if (purrr::is_list(x)) ret_cnt <- length(x)
    if (base::is.data.frame(x) || tibble::is_tibble(x)) ret_cnt <- base::nrow(x)
    
    if (type=="line") {
      readr::write_lines(x=x, file=file, append=append)
    } else if (type=="csv") {
      readr::write_csv(x=x, file=file, col_names=cols)
    } else if (type=="tsv") {
      readr::write_tsv(x=x, file=file, col_names=cols)
    } else if (type=="rds") {
      readr::write_rds(x, file, compress="gz")
    } else {
      stop(cat(glue::glue("{RET}[{funcTag}]: ERROR: Invalid type={type}{RET}")))
      return(NULL)
    }
    
    if (done) {
      done_file <- paste0(file,".done.txt")
      cmd <- glue::glue("touch {done_file}")
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Running cmd='{cmd}'...{RET}"))
      base::system(cmd)
    }
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_cnt
}

build_file_dir = function(file) {
  dir <- base::dirname(file)
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  invisible(dir)
}

print_tib = function(t, f,  v=0,vt=3,tc=1,l=3, n=NULL,m=NULL) {
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  tib     <- t
  funcTag <- f
  name    <- n
  mssg    <- m
  verbose <- v
  
  ret_cnt <- 0
  ret_tib <- NULL
  
  if (verbose>=vt) {
    if (!is.null(tib) && !is.na(tib)) {
      ret_cnt <- tib %>% base::nrow()
      if (!is.null(mssg)) cat(glue::glue("[{funcTag}]:{tabsStr} {mssg}{RET}"))
      if (!is.null(name)) {
        cat(glue::glue("[{funcTag}]:{tabsStr} tibble({name}) row count=[{ret_cnt}]:{RET}"))
      } else {
        cat(glue::glue("[{funcTag}]:{tabsStr} tibble row count=[{ret_cnt}]:{RET}"))
      }
      print(tib, n=l)
    } else {
      if (!is.null(name)) {
        cat(glue::glue("[{funcTag}]:{tabsStr} tibble({name}) row count=[NULL]:{RET}"))
      } else {
        cat(glue::glue("[{funcTag}]:{tabsStr} tibble row count=[NULL]:{RET}"))
      }
    }
  }
  
  # invisible(ret_tib)
  ret_cnt
}

splitStrToVec = function(x, del=',', unique=TRUE,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'splitStrToVec'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_vec <- NULL
  if (!is.null(x)) 
    ret_vec <- str_split(x, pattern=del, simplify=TRUE) %>% 
    as.vector() %>% unique()
  
  ret_vec
}

spaceToUnderscore = function(x) {
  x %>% stringr::str_replace_all(' ','_')
}

onlyAlphaNumeric = function(x) {
  stringr::str_replace_all(x, "[^[:alnum:]]", "")
}


# End of file
