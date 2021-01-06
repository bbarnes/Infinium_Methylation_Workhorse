
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
#                      Manifest Intersection Method::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

manifestToAnnotation = function(tib, ann, gen,
                                verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'manifestToAnnotation'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      4.2 Build Genomic Region Set::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Force ID's to be unique::
    tib <- tib %>% dplyr::mutate(
      rep=row_number(),
      Seq_ID=paste(Seq_ID,rep, sep='_')) %>% 
      ungroup()

    man_grs <- GRanges(
      seqnames = Rle(tib$chrom), 
      IRanges(start=tib$chromStart, 
              end=tib$chromEnd, 
              names=tib$Seq_ID) )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                5.1 Intersect Mapped Probes with Annotation::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # NCBI Genes::
    ncbi_mat_tib <- NULL
    ncbi_ann_tsv <- file.path(ann,gen,paste(gen,'ncbi.RefSeqGenes.tsv.gz', sep='.') )
    
    if (!is.null(ncbi_ann_tsv) && file.exists(ncbi_ann_tsv)) {
      if (verbose>=1) cat(glue::glue("[{funcTag}]:{tabsStr} ncbi_ann_tsv={ncbi_ann_tsv}...{RET}"))
      
      ncbi_ann_tib <- suppressMessages(suppressWarnings( readr::read_tsv(ncbi_ann_tsv) ))
      colnames(ncbi_ann_tib)[1] <- stringr::str_remove(colnames(ncbi_ann_tib)[1], '^#')
      
      ncbi_ref_grs <- loadNcbiGeneGR(file=ncbi_ann_tsv, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      ncbi_mat_tib <- intersectGranges(man=man_grs,ref=ncbi_ref_grs,
                                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
        dplyr::left_join(dplyr::select(ncbi_ann_tib,name,name2), by=c("Gene"="name")) %>% 
        dplyr::rename(Transcript=Gene, Gene=name2) %>% 
        dplyr::select(Seq_ID, Gene,Transcript,dplyr::everything()) %>%
        dplyr::mutate(Source="NCBI")
      
      if (verbose>=1) cat(glue::glue("[{funcTag}]:{tabsStr} Done. ncbi.{RET}{RET}"))
    }
    
    # UCSC Genes::
    gene_mat_tib <- NULL
    gene_ann_tsv <- file.path(ann,gen,paste(gen,'ucsc.knownGene.tsv.gz', sep='.') )
    if (!is.null(gene_ann_tsv) && file.exists(gene_ann_tsv)) {
      if (verbose>=1) cat(glue::glue("[{funcTag}]:{tabsStr} gene_ann_tsv={gene_ann_tsv}...{RET}"))
      
      gene_ann_tib <- suppressMessages(suppressWarnings( readr::read_tsv(gene_ann_tsv) )) %>% 
        dplyr::rename(name2=proteinID)
      colnames(gene_ann_tib)[1] <- stringr::str_remove(colnames(gene_ann_tib)[1], '^#')
      
      # proteinID
      gene_ref_grs <- loadUcscGeneGR(file=gene_ann_tsv, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      gene_mat_tib <- intersectGranges(man=man_grs,ref=gene_ref_grs,
                                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        dplyr::left_join(dplyr::select(gene_ann_tib,name,name2), by=c("Gene"="name")) %>% 
        dplyr::rename(Transcript=Gene, Gene=name2) %>% 
        dplyr::select(Seq_ID, Gene,Transcript,dplyr::everything()) %>%
        dplyr::mutate(Source="UCSC")
      
      if (verbose>=1) cat(glue::glue("[{funcTag}]:{tabsStr} Done. ucsc.{RET}{RET}"))
    }
    
    # CpG Islands::
    cpgs_mat_tib <- NULL
    cpgs_ann_tsv <- file.path(ann,gen,paste(gen,'ucsc.CpG-Islands.tsv.gz', sep='.') )
    if (!is.null(cpgs_ann_tsv) && file.exists(cpgs_ann_tsv)) {
      if (verbose>=1) cat(glue::glue("[{funcTag}]:{tabsStr} cpgs_ann_tsv={cpgs_ann_tsv}...{RET}"))
      
      cpgs_ann_tib <- suppressMessages(suppressWarnings( readr::read_tsv(cpgs_ann_tsv) )) %>% dplyr::select(-1)
      cpgs_ref_grs <- loadUcscCpgsGR(file=cpgs_ann_tsv, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      cpgs_mat_tib <- intersectGranges(man=man_grs,ref=cpgs_ref_grs,
                                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        dplyr::mutate(Source="UCSC")
      
      if (verbose>=1) cat(glue::glue("[{funcTag}]:{tabsStr} Done. cpgs.{RET}{RET}"))
    }
    
    ret_tib <- dplyr::bind_rows(ncbi_mat_tib,gene_mat_tib,cpgs_mat_tib) %>% 
      dplyr::arrange(chrom,chromStart,chromEnd) %>%
      dplyr::mutate(Seq_Rep=stringr::str_remove(Seq_ID, '^.*_'),
                    Seq_ID=stringr::str_remove(Seq_ID, '_[0-9]+$'))
    
    if (verbose>=vt) cat(glue::glue("[{par$prgmTag}]: Writing annotation mappings={ann_base_csv}...{RET}"))
    readr::write_csv(ret_tib, ann_base_csv)
    if (verbose>=vt) cat(glue::glue("[{par$prgmTag}]: Done. Writing annotation mappings.{RET}{RET}"))
    
    # man_ana_key <- names(ret_tib)[1]
    # man_ana_beg <- names(ret_tib)[2]
    # man_ana_beg <- names(ret_tib)[2]
    # 
    # man_ana_col <- ret_tib %>% dplyr::select(-1) %>% names()
    # man_all_col <- ret_tib %>% names()
    
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

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
      # Last chage:: # dplyr::select(Seq_ID),
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


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Standard UCSC Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadNcbiGeneGR = function(file,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadUcscGeneGR'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
    colnames(dat_tib)[1] <- stringr::str_remove(colnames(dat_tib)[1], '^#')
    dat_tib <- dat_tib %>% dplyr::distinct(name, .keep_all=TRUE) # %>% dplyr::mutate(transcript=name, name=name2)
    print(dat_tib)
    
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
      GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand), # seqinfo=dat_tib$name2,
              IRanges(start=dat_tib$tss_1500_beg, end=dat_tib$tss_1500_end, names=dat_tib$name) )
    ret_grs$tss_200 <- 
      GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand), # seqinfo=dat_tib$name2,
              IRanges(start=dat_tib$tss_200_beg, end=dat_tib$tss_200_end, names=dat_tib$name) )
    ret_grs$tss_body <- 
      GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand), # seqinfo=dat_tib$name2,
              IRanges(start=dat_tib$txStart, end=dat_tib$txEnd, names=dat_tib$name) )
    
    ret_cnt <- ret_grs %>% names %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_grs
}

loadUcscGeneGR = function(file,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadUcscGeneGR'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
    colnames(dat_tib)[1] <- stringr::str_remove(colnames(dat_tib)[1], '^#')
    dat_tib <- dat_tib %>% dplyr::distinct(name, .keep_all=TRUE)
    print(dat_tib)
    
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
