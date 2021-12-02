
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

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
#                           BED/Tabix Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

write_tabix_bed = function(tib, dir, name, 
                           head_char=NULL, s_idx=1, b_idx=2, e_idx=3,
                           header=FALSE,
                           ref=NULL, pipe=NULL,
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'write_tabix_bed'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    
    ret_tib  <- tib
    bed_file <- file.path(dir, paste(name,'sorted.bed', sep='.'))
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing='{bed_file}'...{RET}"))
    if (header) {
      readr::write_tsv(ret_tib, bed_file) # , col_names = FALSE)
    } else {
      readr::write_tsv(ret_tib, bed_file, col_names = FALSE)
    }
    cmd <- glue::glue("bgzip -f {bed_file}")
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Running='{cmd}'...{RET}"))
    system(cmd)
    bed_file <- file.path(dir, paste(name,'sorted.bed.gz', sep='.'))
    
    cmd <- glue::glue("tabix")
    if (!is.null(head_char)) cmd <- glue::glue("{cmd} -c {head_char}")
    cmd <- glue::glue("{cmd} -s {s_idx} -b {b_idx} -e {e_idx} -f {bed_file}")
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Running='{cmd}'...{RET}"))
    system(cmd)
    
    if (!is.null(ref) && file.exists(ref)) {
      ref_name <- ref %>% base::basename() %>%
        stringr::str_remove(".gz$") %>% 
        stringr::str_remove(".bed$") %>%
        stringr::str_remove(".vcf$") %>%
        stringr::str_remove(".csv$") %>%
        stringr::str_remove(".tsv$")
      
      int_file <- file.path(dir, paste(name,ref_name,"intersect.tsv.gz", sep='.'))
      cmd <- glue::glue("tabix -h -R {bed_file} {ref}")
      if (!is.null(pipe)) cmd <- glue::glue("{cmd} | {pipe}")
      cmd <- glue::glue("{cmd} | gzip -c -> {int_file}")
      
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Running='{cmd}'...{RET}"))
      system(cmd)
      
      ret_tib <- load_dbSNP_vcf(int_file, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    } else {
      ret_tib <- bed_file
    }
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}


# End of file
