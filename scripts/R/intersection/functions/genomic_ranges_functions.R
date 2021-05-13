
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("GenomicRanges",quietly=TRUE) ))

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
#                          Genomic Range Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

intersect_GRS = function(can,ref, can_key="unq_can_key", ref_prefix=NULL, # ref_prefix="imp",
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'intersect_GRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    map_tib <- 
      GenomicRanges::findOverlaps(can,ref, ignore.strand=TRUE) %>%
      as.data.frame() %>% tibble::as_tibble()
    map_cnt <- print_tib(map_tib,funcTag, verbose,vt+4,tc, n="map")
    
    can_tib <- can %>% as.data.frame() %>%
      rownames_to_column(var=can_key) %>% 
      tibble::as_tibble() %>%
      dplyr::select(dplyr::all_of(can_key))
    can_cnt <- print_tib(can_tib,funcTag, verbose,vt+4,tc, n="can")
    
    ref_tib <- ref %>% as.data.frame() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(seqnames=as.character(seqnames),
                    strand=as.character(strand)) %>%
      dplyr::rename(chr=seqnames,
                    pos=start,
                    top_srd=strand)
    
    if (!is.null(ref_prefix)) ref_tib <- ref_tib %>% 
      purrr::set_names(paste(ref_prefix,names(.), sep="_"))
    
    ref_cnt <- print_tib(ref_tib,funcTag, verbose,vt+4,tc, n="ref")
    
    # Column Bind Data Sets::
    ret_tib <- dplyr::bind_cols(
      can_tib[map_tib$queryHits, ],
      ref_tib[map_tib$subjectHits,]
    )
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# End of file
