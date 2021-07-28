
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))

COM  <- ","
TAB  <- "\t"
RET  <- "\n"
RET2 <- "\n\n"
BNG  <- "|"
BRK  <- paste0("# ",
               paste(rep("-----",6),collapse=" "),"|",
               paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func = function(tib,
                         
                         # out_dir, run_tag, re_load=FALSE,pre=c(),del='.',
                         
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='template_func') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   funcTag={funcTag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  # etime   <- 0
  # ret_cnt <- 0
  # ret_tib <- NULL
  # ret_dat <- NULL
  # 
  # if (re_load) {
  #   out_file <- file.path(out_dir,funcTag,paste())
  #   if (valid_time_stamp(pre_tag)) {
  #     
  #   }
  # }
  
  stime <- base::system.time({
    
    # verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Reload for Common Method::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

redata = function(out_dir, 
                  run_tag, 
                  fun_tag,
                  re_load = FALSE,
                  pre_tag = NULL,
                  end_str = 'csv.gz', 
                  sep_chr = '.',
                  verbose=0, vt=3,tc=1,tt=NULL,
                  funcTag='reload_dat') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{fun_tag}]:{tabs}")
  # mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt+8) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+8) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   out_dir={out_dir}.{RET}"))
    cat(glue::glue("{mssg}   run_tag={run_tag}.{RET}"))
    cat(glue::glue("{mssg}   fun_tag={fun_tag}.{RET}"))
    cat(glue::glue("{mssg}   pre_tag={pre_tag}.{RET}"))
    cat(glue::glue("{mssg}   end_str={end_str}.{RET}"))
    cat(glue::glue("{mssg}   sep_chr={sep_chr}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  out_dir <- file.path(out_dir, fun_tag)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_fns <- paste(run_tag, fun_tag, end_str, sep=sep_chr)
  out_csv <- file.path(out_dir, out_fns)

  if (re_load==FALSE) return(out_csv)
  if (!file.exists(out_csv)) return(out_csv)
  pre_tag <- c(pre_tag, out_csv,paste0(out_csv, '.done.txt'))
  
  etime   <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  
  stime <- base::system.time({
    val_ret <- valid_time_stamp(pre_tag, verbose=opt$verbose)
    if (val_ret==TRUE) {
      if (verbose>=vt)
        cat(glue::glue("{mssg} Reloading:: file={out_csv}...{RET}"))
      
      tt$addFile(out_csv)
      ret_tib <- safe_read(out_csv, funcTag=funcTag,
                           verbose=verbose,vt=vt+4,tc=tc+1,tt=tt)
      ret_key <- glue::glue("reloaded-tib({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+8,tc, n=ret_key)
    } else {
      ret_tib <- out_csv
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))

  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Common Booster Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

clean_tibble = function(tib,
                        verbose=0,vt=6,tc=1,tt=NULL,
                        funcTag='clean_tibble') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- tib %>%
      select(where(~sum(!is.na(.x)) > 0)) %>%
      utils::type.convert(as.is=TRUE) %>%
      dplyr::mutate(across(where(is.factor), as.character) )
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

makeFieldUnique = function(tib, field, add=NULL,
                           verbose=0,vt=3,tc=1,tt=NULL,
                           funcTag='makeFieldUnique') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting; field={field}...{RET}"))
  
  field_sym <- rlang::sym(field)
  outkey_sym <- rlang::sym(field)
  if (!is.null(add)) outkey_sym <- rlang::sym(add)
  
  tot_cnt <- tib %>% base::nrow()
  unq_cnt <- tib %>% dplyr::distinct(!!field_sym) %>% base::nrow()
  # unq_cnt <- tib %>% dplyr::distinct(!!field_sym) %>% base::nrow()
  
  if (tot_cnt==unq_cnt) {
    if (verbose>=vt)
      cat(glue::glue("{mssg} Field({field}) already unique!",
                     "Total={tot_cnt}, Unique={unq_cnt}.{RET}"))
    return(tib)
  }
  stopifnot(tot_cnt>unq_cnt)
  
  if (verbose>=vt)
    cat(glue::glue("{mssg} Field({field}): ",
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
      cat(glue::glue("{mssg} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

get_fileSuffix = function(file,
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='get_fileSuffix') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_val <- NULL
  
  ret_val <- stringr::str_remove(file, "\\.gz$") %>% 
    stringr::str_remove("^.*\\.")
  ret_cnt <- ret_val %>% length()
  
  ret_val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Clean Files for docker::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

clean_file = function(file, touch=FALSE,
                      verbose=0,vt=6,tc=1,tt=NULL,
                      funcTag='clean_file') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Cleaning='{file}'{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    cmd <- NULL
    if (file.exists(file)) unlink(file)
    if (touch) {
      cmd <- glue::glue("touch {file}")
      if (verbose>=vt) 
        cat(glue::glue("{mssg} Running cmd='{cmd}'...{RET}"))
      base::system(cmd)
      
      if (!file.exists(file)) {
        stop(glue::glue("{RET}{mssg} ERROR: ",
                        "File doesn't exist; file={file}!!!{RET}"))
        return(ret_tib)
      }
    }
    
    ret_cnt <- 1
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  file
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                General Program Initialization Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

program_init = function(name,defs=NULL,
                        opts,opt_reqs=NULL,pars,par_reqs=NULL,
                        libs=TRUE,rcpp=FALSE,
                        verbose=0,vt=3,tc=1,tt=NULL,
                        funcTag='program_init') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
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
  if (opts$verbose>=1) opt_tib %>% base::print(n=base::nrow(opt_tib) )
  
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
  if (pars$prgmTag==name) cmd_shell_name <- name
  
  pars$cmd_shell <- 
    file.path(opts$outDir,paste(cmd_shell_name,'command.sh', sep='.'))
  pars$cmd_str <- 
    optsToCommand(opts=opt_tib, pre=opts$Rscript,exe=pars$exePath, 
                  file=pars$cmd_shell, 
                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  
  ret_cnt <- names(opts) %>% length()
  
  # etime <- stime[3] %>% as.double() %>% round(2)
  etime <- 0
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  opts
}

program_done = function(opts,pars,
                        verbose=0,vt=3,tc=1,tt,
                        funcTag='program_done') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
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
  cat(glue::glue("{RET}[{pars$prgmTag}]: Finished(time={sysTime}){RET2}"))
  
  return(0)
}

load_libraries = function(opts, pars, rcpp=FALSE,
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='load_libraries') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  
  stopifnot(!is.null(pars$prgmDir))
  stopifnot(!is.null(pars$scrDir))
  
  # list.files("/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/scripts/R", full.names = TRUE) %>% file.path("functions")
  # list.files(pars$scrDir, full.names = TRUE) %>% file.path("functions")
  
  pars$prgm_src_dir <- file.path(pars$scrDir,pars$prgmDir, 'functions')
  if (!dir.exists(pars$prgm_src_dir)) stop(glue::glue("[{pars$prgmTag}]: Program Source={pars$prgm_src_dir} does not exist!{RET}"))
  for (sfile in list.files(path=pars$prgm_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
  cat(glue::glue("[{pars$prgmTag}]: Done. Loading Source Files form Program Source={pars$prgm_src_dir}!{RET2}") )
  
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
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  opts
}

check_list = function(list, vals, name="options",
                      verbose=0,vt=3,tc=1,tt=NULL,
                      funcTag='check_list') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
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
      cat(glue::glue("{mssg} Validated all required options({vals_len}).{RET}"))
  } else {
    stop(glue::glue("{mssg} Invalid Options; use --help!!!{RET2}"))
    return(NULL)
  }
  
  ret_cnt <- names(list) %>% length()
  
  # etime <- stime[3] %>% as.double() %>% round(2)
  etime <- 0
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  list
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Basic Tibble Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Check if all files exist and have correct time stamp order
valid_time_stamp = function(files,
                            funcTag='valid_time_stamp',
                            verbose=0,vt=6,tc=1,tt=NULL) {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_val <- TRUE
  stime <- system.time({
    
    file_cnt <- length(files)
    for (ii in c(1:file_cnt)) {
      if (!file.exists(files[ii])) ret_val <- FALSE
      if (verbose>=vt+4)
        cat(glue::glue("{mssg}{TAB} File({ii})={files[ii]} exist={ret_val}.{RET}"))
      
    }
    if (verbose>=vt+4)
      cat(glue::glue("{mssg} Files exist={ret_val}.{RET2}"))
    
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
          cat(glue::glue("ii={ii}, Valid Order={ret_val}.{RET2}"))
      }
    }
    ret_cnt <- files %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  ret_val
}

# This function is out of date replacing with simpler function above::
#  valid_time_stamp()
check_timeStamps = function(name,outDir,origin=NULL,files,
                            verbose=0,vt=3,tc=1,tt=NULL,
                            funcTag='check_timeStamps') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
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
    if (verbose>=vt+4) cat(glue::glue("{mssg} Origin exist={isValid}.{RET}"))
    
    if (!file.exists(beg_txt)) isValid <- FALSE
    if (verbose>=vt+4) cat(glue::glue("{mssg} Beg exist={isValid}.{RET}"))
    
    if (!file.exists(end_txt)) isValid <- FALSE
    if (verbose>=vt+4) cat(glue::glue("{mssg} End exist={isValid}.{RET}"))
    
    for (file in files) {
      if (!file.exists(file)) isValid <- FALSE
    }
    if (verbose>=vt+4) cat(glue::glue("{mssg} All files exist={isValid}.{RET}"))
    
    
    # Check all time stamps::
    if (isValid && !is.null(origin) && !file.exists(origin) &&
        file.mtime(orgin) > file.mtime(beg_txt)) isValid <- FALSE
    if (verbose>=vt+4) cat(glue::glue("{mssg} Origin stamp check={isValid}.{RET}"))
    
    if (isValid && file.mtime(beg_txt) > file.mtime(end_txt)) isValid <- FALSE
    if (verbose>=vt+4) cat(glue::glue("{mssg} Beg/End stamp check={isValid}.{RET}"))
    
    if (isValid) {
      for (file in files) {
        if (file.mtime(beg_txt) >= file.mtime(file) ||
            file.mtime(end_txt) <= file.mtime(file)) {
          if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabs}{TAB} Failed stamp={file}.{RET}"))
          isValid <- FALSE
          break
        }
      }
    }
    
    # Clear everything if validation failed for any reason::
    if (!isValid) {
      if (verbose>=vt) cat(glue::glue("{mssg} Cleaning all files...{RET}"))
      
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
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Basic Tibble Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

reduceSortedTib = function(tib, n=3,
                           verbose=0,vt=3,tc=1,tt=NULL,
                           funcTag='reduceSortedTib') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting; n={n}...{RET}"))
  
  tib_cnt <- tib %>% base::nrow()
  if (tib_cnt <= 1) return(tib)
  if (tib_cnt <= n) return(tib)
  
  mid_cnt <- as.integer( tib_cnt/n ) + 1
  sel_vec <- c(1,mid_cnt,tib_cnt) %>% unique()
  tib <- tib[sel_vec, ]
  if (verbose>=vt) cat(glue::glue("{mssg} Done.{RET2}"))
  
  tib
}

joinTibbles = function(a, b, by, side="left", 
                       verbose=0,vt=3,tc=1,
                       funcTag='joinTibbles') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  # if (verbose>=vt) cat(glue::glue("{mssg} Starting; platform={platform}.{RET}"))
  
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  
  if (side=='left')  return(dplyr::left_join(a,  b, by=by))
  if (side=='right') return(dplyr::right_join(a, b, by=by))
  if (side=='inner') return(dplyr::inner_join(a, b, by=by))
  if (side=='full')  return(dplyr::full_join(a,  b, by=by))
  stop(glue::glue("[$func]: ERROR: Unsupported join={side}!!!{RET2}"))
  
  # if (verbose>=vt) cat(glue::glue("{mssg} Done.{RET2}"))
  
  NULL
}

addColNames = function(tib, add, fix, prefix=TRUE, del='_',
                       verbose=0,vt=3,tc=1,tt=NULL,
                       funcTag='addColNames') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
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
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Options to Script Commands::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

optsToCommand = function(opts, pre=NULL, exe, rm=NULL, add=NULL, 
                         file=NULL, key="Option", val="Value",
                         cluster=TRUE,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='optsToCommand') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
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
    if (verbose>=vt+4) cat(glue::glue("{mssg} bool_str='{bool_str}'.{RET2}"))
  }
  
  # Merge Options::
  if (!is.null(opts)) {
    opt_str <- opts %>% dplyr::arrange(!!key) %>%
      tidyr::unite(Param, !!key, !!val, sep='=') %>%
      dplyr::mutate(Param=stringr::str_c('--',Param)) %>% 
      dplyr::pull(Param) %>%
      stringr::str_c(collapse=" ")
    if (verbose>=vt+4) cat(glue::glue("{mssg} opt_str='{opt_str}'.{RET2}"))
    
    # Second Removal Attempt of rm fields::
    if (!cluster) {
      opt_str <- opt_str %>% stringr::str_remove('--cluster')
      if (verbose>=vt+4) cat(glue::glue("{mssg} opt_str(-cluster)='{opt_str}'.{RET2}"))
    }
  }
  
  # Options:: Add
  if (!is.null(add)) {
    add_str <- add %>% tidyr::unite(Param, !!key, !!val, sep='=') %>%
      dplyr::mutate(Param=stringr::str_c('--',Param)) %>% 
      dplyr::pull(Param) %>%
      stringr::str_c(collapse=" ")
    if (verbose>=vt+4) cat(glue::glue("{mssg} add_str='{add_str}'.{RET2}"))
  }
  
  # Add Executable and Join Options::
  cmd <- ''
  if (!is.null(pre) && length(pre)!=0) cmd <- pre
  cmd <- stringr::str_c(cmd,exe,opt_str,add_str,bool_str, sep=' ')
  if (verbose>=vt+4) cat(glue::glue("{mssg} cmd='{cmd}'.{RET2}"))
  
  ret_val <- cmd
  if (!is.null(file)) {
    dir <- base::dirname(file)
    if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    base::unlink(file)
    
    if (verbose>=vt) cat(glue::glue("{mssg} Writing program shell={file}...{RET}"))
    readr::write_lines(cmd, file=file)
    Sys.chmod(file, mode="0777")
    
    ret_val <- file
  }
  if (verbose>=vt) cat(glue::glue("{mssg} Done.{RET2}"))
  
  ret_val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          File Searching Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

get_file_list = function(dir, pattern, trim=NULL,
                         max=0, recursive=FALSE, 
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='get_file_list') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    file_vec <- list.files(dir, pattern=pattern, full.names=TRUE, 
                           recursive=recursive)
    name_vec <- file_vec %>% base::basename()
    
    if (!is.null(trim)) {
      for (t in trim) {
        name_vec <- stringr::str_remove(name_vec, t)
      }
    }
    
    ret_tib <- stats::setNames(as.list(file_vec), name_vec)
    ret_cnt <- ret_tib %>% names() %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

get_file_list_old = function(dir, pattern, max=0, recursive=FALSE, 
                             verbose=0,vt=3,tc=1,tt=NULL,
                             funcTag='get_file_list_old') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt)
    cat(glue::glue("{mssg} Starting; max={max}, ",
                   "pattern={pattern}, dir={dir}.{RET}"))
  
  files <- NULL
  files <- base::list.files(path=dir, full.names=TRUE)
  if (verbose>=vt) {
    cat(glue::glue("{mssg} All Files (no recursive)={RET}"))
    print(files)
  }
  
  files <- base::list.files(path=dir, full.names=TRUE, recursive=TRUE)
  if (verbose>=vt) {
    cat(glue::glue("{mssg} All Files (recursive)={RET}"))
    print(files)
  }
  
  files <- base::list.files(path=dir, pattern=pattern, full.names=TRUE, 
                            recursive=recursive)
  if (verbose>=vt) {
    cat(glue::glue("{mssg} Files={RET}"))
    print(files)
  }
  if (is.null(files) || length(files)==0) {
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed to find ",
                    "pattern='{pattern}' in dir='{dir}'!!!{RET2}"))
    return(NULL)
  }
  
  files_cnt <- length(files)
  if (max!=0 && files_cnt>max)
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Found {files_cnt} > max({max}) ",
                    "with pattern={pattern} in dir={dir}!!!{RET2}"))
  
  files
}

findFileByPattern = function(dir, pattern, max=0, recursive=FALSE, 
                             verbose=0,vt=3,tc=1,tt=NULL,
                             funcTag='findFileByPattern') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt)
    cat(glue::glue("{mssg} Starting; max={max}, ",
                   "pattern={pattern}, dir={dir}.{RET}"))
  
  files <- NULL
  files <- base::list.files(path=dir, full.names=TRUE)
  if (verbose>=vt) {
    cat(glue::glue("{mssg} All Files (no recursive)={RET}"))
    print(files)
  }
  
  files <- base::list.files(path=dir, full.names=TRUE, recursive=TRUE)
  if (verbose>=vt) {
    cat(glue::glue("{mssg} All Files (recursive)={RET}"))
    print(files)
  }
  
  files <- base::list.files(path=dir, pattern=pattern, full.names=TRUE, 
                            recursive=recursive)
  if (verbose>=vt) {
    cat(glue::glue("{mssg} Files={RET}"))
    print(files)
  }
  if (is.null(files) || length(files)==0) {
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed to find ",
                    "pattern='{pattern}' in dir='{dir}'!!!{RET2}"))
    return(NULL)
  }
  
  files_cnt <- length(files)
  if (max!=0 && files_cnt>max)
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Found {files_cnt} > max({max}) ",
                    "with pattern={pattern} in dir={dir}!!!{RET2}"))
  
  files
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Linux Cluster Options Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

setLaunchExe = function(opts, pars, 
                        verbose=0,vt=3,tc=1,tt=NULL,
                        funcTag='setLaunchExe') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  opts$lanExe <- ''
  opts$isLinux <- FALSE
  
  # if (is.null(pars$lixDir))
  if (is.null(pars$lixDir))
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Options lixDir=NULL not allowed!!!{RET2}"))
  
  if (dir.exists(pars$lixDir)) {
    opts$isLinux <- TRUE
    opts$lanExe <- 'qsub -cwd -pe threaded 16 -l excl=true -N'
    if (dir.exists(pars$macDir)) stop(glue::glue("{RET}[{funcTag}]: Linux/Mac directories exist???{RET2}"))
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
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='guess_file_del') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) 
    cat(glue::glue("{mssg} Starting...{RET}"))
  
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
    
    ret_cnt <- print_tib(t=ret_tib, f=funcTag, v=verbose,vt=vt+1,tc=tc+1)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  ret_val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              String Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: Add all the other options for reading files and add
#   col_names, col_types options to be passed in. 
safe_read = function(file, type=NULL, clean=TRUE, guess_max=1000,
                     verbose=0,vt=5,tc=1,tt=NULL,
                     funcTag='safe_read') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (is.null(file)) {
    stop(glue::glue("{RET}{mssg} file can't be NULL!!{RET2}"))
    return(NULL)
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # Guess type if not provided
    if (is.null(type)) {
      tmp_file <- file %>% stringr::str_remove(".gz$")
      if (stringr::str_ends(tmp_file,".rds")) {
        type <- "rds"
      } else {
        type <- guess_file_del(file, verbose=verbose,vt=vt+6,tc=tc+1,tt=tt)
      }
    }
    if (verbose>=vt)
      cat(glue::glue("{mssg} Reading Data (sep='{type}') ",
                     "file='{file}'...{RET}"))

    if (type=="line") {
      ret_tib <- suppressMessages(suppressWarnings( 
        readr::read_lines(file=file) ))
    } else if (type=="csv" || type==COM) {
      ret_tib <- suppressMessages(suppressWarnings( 
        readr::read_csv(file=file, guess_max=guess_max) ))
    } else if (type=="tsv" || type==TAB || type==" ") {
      ret_tib <- suppressMessages(suppressWarnings(
        readr::read_tsv(file=file, guess_max=guess_max) ))
    } else if (type=="rds") {
      ret_tib <- readr::read_rds(file)
    } else {
      stop(cat(glue::glue("{RET}[{funcTag}]: ERROR: Invalid type='{type}'!{RET2}")))
      return(NULL)
    }
    if (clean) ret_tib <- clean_tibble(ret_tib)
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg}Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))

  ret_tib
}

safe_write = function(x, type=NULL, file=NULL, 
                      cols=TRUE, append=FALSE, done=FALSE,
                      verbose=0,vt=5,tc=1,tt=NULL,
                      funcTag="safe_write") {
  # If file is provided, usually in the case where the user didn't want 
  #  it written simply return invisible and do nothing...
  
  if (is.null(file)) return(base::invisible(NULL))
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  # if (verbose>=vt) cat(glue::glue("{mssg}Starting type={type}.{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    build_file_dir(file)
    
    # Guess type if not provided
    if (is.null(type)) {
      tmp_file <- file %>% stringr::str_remove(".gz$")
      if (stringr::str_ends(tmp_file,".tsv")) {
        type <- "tsv"
      } else if (stringr::str_ends(tmp_file,".csv")) {
        type <- "csv"
      } else if (stringr::str_ends(tmp_file,".rds")) {
        type <- "rds"
      } else {
        type <- "line"
      }
    }
    if (verbose>=vt)
      cat(glue::glue("{mssg} Writing Data {type}={file}...{RET}"))
    
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
      stop(cat(glue::glue("{mssg} ERROR: Invalid type={type}!{RET}")))
      return(NULL)
    }
    
    if (done) {
      done_file <- paste0(file,".done.txt")
      cmd <- glue::glue("touch {done_file}")
      if (verbose>=vt) 
        cat(glue::glue("{mssg}Running cmd='{cmd}'...{RET}"))
      base::system(cmd)
    }
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg}Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_cnt
}

build_file_dir = function(file) {
  dir <- base::dirname(file)
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  invisible(dir)
}

print_tib = function(t, f="print_tib",  
                     v=0,vt=3,tc=1,l=3, n=NULL,m=NULL) {
  
  tib     <- t
  funcTag=f
  name    <- n
  mssg    <- m
  verbose <- v
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")

  ret_cnt <- 0
  ret_tib <- NULL
  
  if (is.null(tib) || base::nrow(t)==0) return(ret_cnt)
  
  if (is.data.frame(tib)) ret_cnt <- tib %>% base::nrow()
  
  if (verbose>=vt) {
    if (!is.null(tib) && !is.na(tib)) {
      # if (!is.null(mssg)) cat(glue::glue("{mssg} {mssg}{RET}"))
      if (!is.null(name)) {
        cat(glue::glue("{mssg} tibble({name}) row count=[{ret_cnt}]:{RET}"))
      } else {
        cat(glue::glue("{mssg} tibble row count=[{ret_cnt}]:{RET}"))
      }
      if (l==0) l <- tib %>% base::nrow()
      print(tib, n=l)
      cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    } else {
      if (!is.null(name)) {
        cat(glue::glue("{mssg} tibble({name}) row count=[NULL]:{RET}"))
      } else {
        cat(glue::glue("{mssg} tibble row count=[NULL]:{RET}"))
      }
    }
  }
  
  # invisible(ret_tib)
  ret_cnt
}

splitStrToVec = function(x, del=',', unique=TRUE,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='splitStrToVec') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg}Starting...{RET}"))
  
  ret_vec <- NULL
  if (!is.null(x)) 
    ret_vec <- stringr::str_split(x, pattern=del, simplify=TRUE) %>% 
    as.vector() # %>% utils::type.convert() %>% # NEVER DO THIS!!!
  
  if (unique) ret_vec <- ret_vec %>% unique()
  
  ret_vec
}

spaceToUnderscore = function(x) {
  x %>% stringr::str_replace_all(' ','_')
}

onlyAlphaNumeric = function(x) {
  stringr::str_replace_all(x, "[^[:alnum:]]", "")
}

as.num = function(x, na.strings = "NA") {
  stopifnot(is.character(x))
  na = x %in% na.strings
  x[na] = 0
  x = as.numeric(x)
  # x[na] = NA_real_
  x
}

# End of file
