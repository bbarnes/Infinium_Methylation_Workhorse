
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
#                     Genomic Range Methods:: Intersection
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

intersect_GRS = function(ref,can,
                         ref_key=NULL,ref_col=NULL,ref_prefix=NULL,ref_red=TRUE,
                         can_key=NULL,can_col=NULL,can_prefix=NULL, 
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'intersect_GRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    #
    # Perform Intersection::
    #
    map_tib <- 
      GenomicRanges::findOverlaps(can,ref, ignore.strand=TRUE) %>%
      as.data.frame() %>% tibble::as_tibble()
    map_key <- glue::glue("map_tib({funcTag})")
    map_cnt <- print_tib(map_tib,funcTag, verbose,vt+4,tc, n=map_key)
    
    #
    # Return Tibble from Genomic Range: Candidate
    #
    can_tib <- as.data.frame(can) %>% tibble::as_tibble(rownames=can_key)
    
    if (!is.null(can_col) && length(can_col)!=0) can_tib <- can_tib %>%
      dplyr::select(dplyr::all_of(can_col))
    
    if (!is.null(can_prefix)) can_tib <- can_tib %>% 
      purrr::set_names(paste(can_prefix,names(.), sep="_"))
    
    can_key <- glue::glue("can_tib({funcTag})")
    can_cnt <- print_tib(can_tib,funcTag, verbose,vt+4,tc, n=can_key)
    
    #
    # Return Tibble from Genomic Range: Reference
    #
    ref_tib <- as.data.frame(ref) %>% tibble::as_tibble(rownames=ref_key)

    if (ref_red) ref_tib <- ref_tib %>% 
      dplyr::mutate(seqnames=as.character(seqnames),
                    strand=as.character(strand)) %>%
      dplyr::rename(chr=seqnames,
                    pos=start,
                    srd=strand)
    
    if (!is.null(ref_col) && length(ref_col)!=0) ref_tib <- ref_tib %>%
      dplyr::select(dplyr::all_of(ref_col))
    
    if (!is.null(ref_prefix)) ref_tib <- ref_tib %>% 
      purrr::set_names(paste(ref_prefix,names(.), sep="_"))
    
    if (FALSE) {
      ref_tib <- ref %>% as.data.frame() %>%
        tibble::as_tibble(rownames=ref_key) %>%
        dplyr::mutate(seqnames=as.character(seqnames),
                      strand=as.character(strand)) %>%
        dplyr::rename(chr=seqnames,
                      pos=start,
                      top_srd=strand)
      
      if (!is.null(ref_prefix)) ref_tib <- ref_tib %>% 
          purrr::set_names(paste(ref_prefix,names(.), sep="_"))
    }
      
    ref_key <- glue::glue("ref_tib({funcTag})")
    ref_cnt <- print_tib(ref_tib,funcTag, verbose,vt+4,tc, n=ref_key)

    #
    # Bind Candidate and Reference Tibs::
    #
    ret_tib <- dplyr::bind_cols(
      can_tib[map_tib$queryHits, ],
      ref_tib[map_tib$subjectHits,]
    )
    
    ret_key <- glue::glue("ret_FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Intersection Method:: 
#                                 OLD CODE!!!
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

intersectGranges = function(man,ref,
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'intersectGranges'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (! purrr::is_list(ref)) {
      map_tib <- 
        GenomicRanges::findOverlaps(man,ref, ignore.strand=TRUE) %>%
        as.data.frame() %>% tibble::as_tibble()
      
      mani_tib <- man %>% as.data.frame() %>%
        rownames_to_column(var='Seq_ID') %>% tibble::as_tibble() 
      mani_len <- mani_tib %>% base::nrow()
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} mani_tib({mani_len})={RET}"))
        print(mani_tib)
      }
      
      gene_tib <- ref %>% as.data.frame() %>%
        rownames_to_column(var='Name') %>% tibble::as_tibble() # %>% 
      # purrr::set_names(c('Gene','chrom','chromStart','chromEnd','chromLength','chromStrand'))
      gene_len <- gene_tib %>% base::nrow()
      
      #
      # TBD::
      # TBD:: Lame way to fix this; the code commented out above:
      # TBD::
      #
      colnames(gene_tib)[1] <- "Gene"
      colnames(gene_tib)[2] <- "chrom"
      colnames(gene_tib)[3] <- "chromStart"
      colnames(gene_tib)[4] <- "chromEnd"
      colnames(gene_tib)[5] <- "chromLength"
      colnames(gene_tib)[6] <- "chromStrand"
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} gene_tib({gene_len})={RET}"))
        print(gene_tib)
      }
      
      #
      # Last change:: # dplyr::select(Seq_ID),
      #
      cur_tib <- dplyr::bind_cols(
        mani_tib[map_tib$queryHits, ], # %>% dplyr::select(Seq_ID),
        gene_tib[map_tib$subjectHits,] ) # %>% dplyr::mutate(Feature=feat_key)
      cur_len <- cur_tib %>% base::nrow()
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} cur_tib({cur_len})={RET}"))
        print(cur_tib)
      }
      ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
    } else {
      for (feat_key in names(ref)) {
        if (verbose>=vt)
          cat(glue::glue("[{funcTag}]:{tabsStr} GRange Overlap; feature={feat_key}...{RET}"))
        
        map_tib <- 
          GenomicRanges::findOverlaps(man,ref[[feat_key]], ignore.strand=TRUE) %>%
          as.data.frame() %>% tibble::as_tibble()
        map_len <- base::nrow(map_tib)
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; man={RET}"))
          print(man)
          
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; ref={RET}"))
          print(ref[[feat_key]])
          
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; map_tib={map_len}{RET}"))
          print(map_tib)
          
          cat(glue::glue("[{funcTag}]:{tabsStr}{RET}{RET}"))
        }
        
        mani_tib <- man %>% as.data.frame() %>%
          rownames_to_column(var='Seq_ID') %>% tibble::as_tibble() 
        mani_len <- mani_tib %>% base::nrow()
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; mani_tib({mani_len})={RET}"))
          print(mani_tib)
        }
        
        gene_tib <- ref[[feat_key]] %>% as.data.frame() %>% 
          rownames_to_column(var='Name') %>% tibble::as_tibble() %>% 
          purrr::set_names(c('Gene','chrom','chromStart','chromEnd','chromLength','chromStrand'))
        gene_len <- gene_tib %>% base::nrow()
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; gene_tib({gene_len})={RET}"))
          print(gene_tib)
        }
        
        cur_tib <- dplyr::bind_cols(
          # tib[map_tib$queryHits, ] %>% dplyr::select(Seq_ID),
          mani_tib[map_tib$queryHits, ] %>% dplyr::select(Seq_ID),
          gene_tib[map_tib$subjectHits,] ) %>%
          dplyr::mutate(Feature=feat_key)
        cur_len <- cur_tib %>% base::nrow()
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} feature={feat_key}; cur_tib({cur_len})={RET}"))
          print(cur_tib)
        }
        ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
      }
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}



# End of file
