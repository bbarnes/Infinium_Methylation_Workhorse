
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Annotation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Genomic Ranges::
suppressWarnings(suppressPackageStartupMessages( base::require("GenomicRanges",quietly=TRUE) ))

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
#                          Intersection Method::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

intersectGranges = function(man,ref,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'intersectGranges'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    for (feat_key in names(ref)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]: GRange Overlap; feature={feat_key}...{RET}"))
      
      map_tib <- 
        GenomicRanges::findOverlaps(man,ref[[feat_key]], ignore.strand=TRUE) %>%
        as.data.frame() %>% tibble::as_tibble()
      map_len <- base::nrow(map_tib)
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]: Mapped={map_len}{RET}"))
        print(map_tib)
      }
      
      gene_tib <- ref[[feat_key]] %>% as.data.frame() %>% 
        rownames_to_column(var='Name') %>% tibble::as_tibble() %>% 
        purrr::set_names(c('Gene','chrom','chromStart','chromEnd','chromLength','chromStrand'))
      
      cur_tib <- dplyr::bind_cols(
        man_pos_tib[map_tib$queryHits, ] %>% dplyr::select(Seq_ID),
        gene_tib[map_tib$subjectHits,] ) %>%
        dplyr::mutate(Feature=feat_key)
      
      ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
      
      if (verbose>=vt) cat(glue::glue("[{funcTag}]: Done.{RET}{RET}"))
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Standard UCSC Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadUcscGeneGR = function(file,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadUcscGeneGR'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
    colnames(dat_tib)[1] <- stringr::str_remove(colnames(dat_tib)[1], '^#')
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Raw Data={RET}"))
    if (verbose>=vt+4) print(dat_tib)
    
    dat_tib <- dat_tib %>% 
      dplyr::mutate(
        gene_tss=dplyr::case_when(
          strand=='+' ~ txStart,
          strand=='-' ~ txEnd,
          TRUE ~ NA_real_
        ),
        
        tss_200_beg=dplyr::case_when(
          strand=='+' ~ gene_tss-200,
          strand=='-' ~ gene_tss,
          TRUE ~ NA_real_
        ),
        tss_200_end=dplyr::case_when(
          strand=='+' ~ gene_tss,
          strand=='-' ~ gene_tss+200,
          TRUE ~ NA_real_
        ),
        
        tss_1500_beg=dplyr::case_when(
          strand=='+' ~ gene_tss-1500,
          strand=='-' ~ gene_tss+200,
          TRUE ~ NA_real_
        ),
        tss_1500_end=dplyr::case_when(
          strand=='+' ~ gene_tss-200,
          strand=='-' ~ gene_tss+1500,
          TRUE ~ NA_real_
        ),
      )
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
    ret_grs$tss_1500 <- 
      GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand),
              IRanges(start=dat_tib$tss_1500_beg, end=dat_tib$tss_1500_end, names=dat_tib$name) )
    ret_grs$tss_200 <- 
      GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand),
              IRanges(start=dat_tib$tss_200_beg, end=dat_tib$tss_200_end, names=dat_tib$name) )
    ret_grs$tss_body <- 
      GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand),
              IRanges(start=dat_tib$txStart, end=dat_tib$txEnd, names=dat_tib$name) )
    
    ret_cnt <- ret_grs %>% names %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_grs
}



loadUcscCpgsGR = function(file,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadUcscCpgsGR'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) )) %>% 
      dplyr::mutate(name=paste0(chrom,'-',chromStart,'-',chromEnd))
    
    # Build GRange Data Structures:: CpG Islands
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
    ret_grs$NShelf <- GRanges(seqnames=Rle(dat_tib$chrom), 
                              IRanges(start=dat_tib$chromStart-4000, end=dat_tib$chromStart-2000, names=dat_tib$name) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Built North Shelves.{RET}"))
    
    ret_grs$NShore <- GRanges(seqnames=Rle(dat_tib$chrom), 
                              IRanges(start=dat_tib$chromStart-2000, end=dat_tib$chromStart, names=dat_tib$name) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Built North Shores.{RET}"))
    
    ret_grs$Island <- GRanges(seqnames=Rle(dat_tib$chrom), 
                              IRanges(start=dat_tib$chromStart, end=dat_tib$chromEnd, names=dat_tib$name) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Built Islands.{RET}"))
    
    ret_grs$SShore <- GRanges(seqnames=Rle(dat_tib$chrom), 
                              IRanges(start=dat_tib$chromEnd, end=dat_tib$chromEnd+2000, names=dat_tib$name) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Built South Shores.{RET}"))
    
    ret_grs$SShelf <- GRanges(seqnames=Rle(dat_tib$chrom), 
                              IRanges(start=dat_tib$chromEnd+2000, end=dat_tib$chromEnd+4000, names=dat_tib$name) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Built South Shelves.{RET}"))
    
    ret_cnt <- ret_grs %>% names %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_grs
}



# End of file
