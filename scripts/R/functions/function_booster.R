
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Common Booster Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Options to Script Commands::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# optsToCommand(opts=opt_tib, exe=par$exePath, rm=c("lociBetaKey","lociPvalKey"), verbose=10) %>% print()
#
optsToCommand = function(opts, pre=NULL, exe, rm=NULL, add=NULL, file=NULL, key="Option", val="Value",
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'optsToCommand'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))

  key <- key %>% as.character() %>% rlang::sym()
  val <- val %>% as.character() %>% rlang::sym()
  
  # Options:: Exclude
  if (!is.null(rm)) opts <- opts %>% dplyr::filter(! (!!key %in% rm) )

  # Handel Boolean Variables::
  bool <- opts %>% dplyr::filter(Value=="TRUE")
  opts <- opts %>% dplyr::filter(Value!="TRUE")
  opts <- opts %>% dplyr::filter(Value!="FALSE")
  
  bool_str <- bool %>% 
    dplyr::mutate(!!key := stringr::str_c('--',!!key)) %>%
    dplyr::pull(!!key) %>% stringr::str_c(collapse=" ")
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} bool_str='{bool_str}'.{RET}{RET}"))
  
  # Merge Options::
  opt_str <- opts %>% dplyr::arrange(!!key) %>%
    tidyr::unite(Param, !!key, !!val, sep='=') %>%
    dplyr::mutate(Param=stringr::str_c('--',Param)) %>% 
    dplyr::pull(Param) %>%
    stringr::str_c(collapse=" ")
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} opt_str='{opt_str}'.{RET}{RET}"))
  
  # Second Removal Attempt of rm fields::
  opt_str <- opt_str %>% stringr::str_remove('--cluster')
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} opt_str(-cluster)='{opt_str}'.{RET}{RET}"))
  
  # Options:: Add
  add_str <- add %>% tidyr::unite(Param, !!key, !!val, sep='=') %>%
    dplyr::mutate(Param=stringr::str_c('--',Param)) %>% 
    dplyr::pull(Param) %>%
    stringr::str_c(collapse=" ")
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} add_str='{add_str}'.{RET}{RET}"))
  
  # Add Executable and Join Options::
  cmd <- ''
  if (!is.null(pre) && length(pre)!=0) cmd <- pre
  cmd <- stringr::str_c(cmd,exe,opt_str,add_str,bool_str, sep=' ')
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} cmd='{cmd}'.{RET}{RET}"))

  if (!is.null(file)) {
    dir <- base::dirname(file)
    if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    base::unlink(file)
    
    readr::write_lines(cmd, path=file)
    Sys.chmod(file, mode="0777")
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  file
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          File Searching Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

findFileByPattern = function(dir, pattern, max=1, recursive=FALSE, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'findFileByPattern'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; max={max}, pattern={pattern}, dir={dir}.{RET}"))
  
  files <- NULL
  files <- list.files(path=dir, pattern=pattern, full.names=TRUE, recursive=recursive)
  if (is.null(files) || length(files)==0)
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed to find pattern={pattern} in dir={dir}!!!{RET}{RET}"))
  
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
  
  if (is.null(pars$lixDir))
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Options lixDir=NULL not allowed!!!{RET}{RET}"))
  
  if (dir.exists(pars$lixDir)) {
    opts$isLinux <- TRUE
    opts$lanExe <- 'qsub -cwd -pe threaded 16 -l excl=true -N'
    # if (opts$isLinux) opts$Rscript <- ''
    if (dir.exists(pars$macDir)) stop(glue::glue("{RET}[{funcTag}]: Linux/Mac directories exist???{RET}{RET}"))
  }
  
  opts
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Common Tibble Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Summarizing Counting Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
cntPer_lte = function(x, min, prc=3) {
  round(100*length(x[which(x<=min)])/length(x),prc)
}

cntPer_gt = function(x, max, prc=3) {
  round(100*length(x[which(x>max)])/length(x),prc)
}

cntPer_gt0 = function(x, max, prc=3) {
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              String Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
spaceToUnderscore = function(x) {
  x %>% stringr::str_replace_all(' ','_')
}

onlyAlphaNumeric = function(x) {
  stringr::str_replace_all(x, "[^[:alnum:]]", "")
}


# End of file
