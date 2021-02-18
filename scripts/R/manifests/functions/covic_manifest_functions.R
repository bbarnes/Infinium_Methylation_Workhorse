# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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
#                            Build Annotations::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadNcbiGene_COVIC = function(file,
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadNcbiGene_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    
    ret_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) )) %>%
      dplyr::mutate(
        strand_FR=stringr::str_replace(strand,'\\+',"F") %>% stringr::str_replace("-","R"), 
        Unique_ID=paste(name,chrom,strand_FR,txStart, sep="_")
      ) %>%
      dplyr::select(Unique_ID,dplyr::everything())
    
    # colnames(ret_tib)[1] <- stringr::str_remove(colnames(ret_tib)[1], '^#')
    # ret_tib <- ret_tib %>% dplyr::distinct(name, .keep_all=TRUE) # %>% dplyr::mutate(transcript=name, name=name2)
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Raw Data={RET}"))
    if (verbose>=vt+4) print(ret_tib)
    
    ret_tib <- ret_tib %>% 
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
    
    ret_tib <- dplyr::bind_rows(
      dplyr::select(ret_tib, Unique_ID, chrom,tss_200_beg,tss_200_end,strand,name,name2) %>% 
        dplyr::mutate(Feature="TSS200",
                      Unique_ID=paste(Unique_ID,Feature, sep="_")) %>% 
        purrr::set_names("Unique_ID","chr","beg","end","srd","tran","gene","feat"),
      
      dplyr::select(ret_tib, Unique_ID, chrom,tss_1500_beg,tss_1500_end,strand,name,name2) %>% 
        dplyr::mutate(Feature="TSS1500",
                      Unique_ID=paste(Unique_ID,Feature, sep="_")) %>% 
        purrr::set_names("Unique_ID","chr","beg","end","srd","tran","gene","feat"),
      
      dplyr::select(ret_tib, Unique_ID, chrom,txStart,txEnd,strand,name,name2) %>% 
        dplyr::mutate(Feature="TSSBody",
                      Unique_ID=paste(Unique_ID,Feature, sep="_")) %>% 
        purrr::set_names("Unique_ID","chr","beg","end","srd","tran","gene","feat")
    )
    
    ret_cnt <- ret_tib %>% names %>% length()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Genomic Ranges Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

intersectGranges_COVIC = function(can,ref,
                                  verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'intersectGranges_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    map_tib <- 
      GenomicRanges::findOverlaps(can,ref, ignore.strand=TRUE) %>%
      as.data.frame() %>% tibble::as_tibble()
    
    can_tib <- can %>% as.data.frame() %>%
      rownames_to_column(var='ID') %>% tibble::as_tibble() %>%
      purrr::set_names(paste(names(.),"CAN", sep="_"))
    can_len <- can_tib %>% base::nrow()
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} can_tib({can_len})={RET}"))
      print(can_tib)
    }
    
    ref_tib <- ref %>% as.data.frame() %>%
      rownames_to_column(var='ID') %>% tibble::as_tibble() %>% 
      purrr::set_names(paste(names(.),"REF", sep="_"))
    ref_len <- ref_tib %>% base::nrow()
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ref_tib({ref_len})={RET}"))
      print(ref_tib)
    }
    
    # Column Bind Datasets::
    ret_tib <- dplyr::bind_cols(
      can_tib[map_tib$queryHits, ],
      ref_tib[map_tib$subjectHits,]
    ) %>% 
      # Add Best Matching Top/Bot & Converted/Opposite Strand::
      dplyr::mutate(
        imp_prbTC=dplyr::case_when(
          prb_des_CAN=="2" ~ stringr::str_remove(prb1U_TC_REF, "[A-Za-z]$"),
          TRUE ~ prb1U_TC_REF
        ),
        imp_prbTO=dplyr::case_when(
          prb_des_CAN=="2" ~ stringr::str_remove(prb1U_TO_REF, "[A-Za-z]$"),
          TRUE ~ prb1U_TO_REF
        ),
        imp_prbBC=dplyr::case_when(
          prb_des_CAN=="2" ~ stringr::str_remove(prb1U_BC_REF, "[A-Za-z]$"),
          TRUE ~ prb1U_BC_REF
        ),
        imp_prbBO=dplyr::case_when(
          prb_des_CAN=="2" ~ stringr::str_remove(prb1U_BO_REF, "[A-Za-z]$"),
          TRUE ~ prb1U_BO_REF
        ),
        imp_mat_srd=dplyr::case_when(
          prb_des_CAN=="U" & prb_aln_50U_CAN==imp_prbTC ~ "TC",
          prb_des_CAN=="U" & prb_aln_50U_CAN==imp_prbTO ~ "TO",
          prb_des_CAN=="U" & prb_aln_50U_CAN==imp_prbBC ~ "BC",
          prb_des_CAN=="U" & prb_aln_50U_CAN==imp_prbBO ~ "BO",
          
          prb_des_CAN=="M" & prb_aln_50M_CAN==imp_prbTC ~ "TC",
          prb_des_CAN=="M" & prb_aln_50M_CAN==imp_prbTO ~ "TO",
          prb_des_CAN=="M" & prb_aln_50M_CAN==imp_prbBC ~ "BC",
          prb_des_CAN=="M" & prb_aln_50M_CAN==imp_prbBO ~ "BO",
          
          prb_des_CAN=="2" & prb_aln_49U_CAN==imp_prbTC ~ "TC",
          prb_des_CAN=="2" & prb_aln_49U_CAN==imp_prbTO ~ "TO",
          prb_des_CAN=="2" & prb_aln_49U_CAN==imp_prbBC ~ "BC",
          prb_des_CAN=="2" & prb_aln_49U_CAN==imp_prbBO ~ "BO",
          TRUE ~ NA_character_)
      )
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
      
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      sum_tib <- ret_tib %>% 
        dplyr::group_by(prb_src_CAN,prb_des_CAN,bsp_din_scr_CAN,
                        bsp_srd_CAN,imp_mat_srd) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      sum_tib %>% print(n=base::nrow(sum_tib))
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
#                       Standard UCSC Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadNcbiGeneGR_COVIC = function(file,
                                verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadNcbiGeneGR_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
    # colnames(dat_tib)[1] <- stringr::str_remove(colnames(dat_tib)[1], '^#')
    dat_tib <- dat_tib %>% dplyr::distinct(name, .keep_all=TRUE)
    # %>% dplyr::mutate(transcript=name, name=name2)
    
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Formatted BSP IO::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bspToGenomicRegion_COVIC = function(bsp, rds=NULL,
                                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'bspToGenomicRegion_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (!is.null(rds)) {
      dir <- base::dirname(rds)
      if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    }
    
    ret_grs <- 
      GRanges(
        seqnames = Rle(paste0("chr",bsp$bsp_chr)),
        strand=Rle(stringr::str_sub( bsp$bsp_srd, 1,1 ) ),
        
        prb_add = bsp$prb_add,
        prb_cgn = bsp$prb_cgn,
        prb_srd = bsp$prb_srd,
        prb_des = bsp$prb_des,
        prb_src = bsp$prb_src,
        
        prb_ord_seq  = bsp$ord_seq,
        prb_aln_49U  = stringr::str_sub(bsp$aln_seq,2),
        prb_aln_50U  = bsp$aln_seq,
        prb_aln_50M  = stringr::str_replace_all(bsp$ord_seq,"G","A"),
        bsp_ref_seq  = bsp$bsp_ref,
        
        bsp_din_scr = bsp$bsp_din_scr,
        bsp_din_ref = bsp$bsp_din_ref,
        bsp_nxb_ref = bsp$bsp_nxb_ref,
        bsp_din_bsc = bsp$bsp_din_bsc,
        bsp_nxb_bsc = bsp$bsp_nxb_bsc,
        
        bsp_tag = bsp$bsp_tag,
        bsp_srd = bsp$bsp_srd,
        bsp_mis   = bsp$bsp_mis,
        bsp_gap   = bsp$bsp_gap,
        
        bsp_hit0=bsp$bsp_hit0,
        bsp_hit1=bsp$bsp_hit1,
        bsp_hit2=bsp$bsp_hit2,
        bsp_hit3=bsp$bsp_hit3,
        bsp_hit4=bsp$bsp_hit4,
        bsp_hit5=bsp$bsp_hit5,
        
        IRanges(start=bsp$bsp_pos,
                end=bsp$bsp_pos+1,
                names=bsp$Unique_ID
        )
      )
    
    if (!is.null(rds)) {
      if (verbose>=1)
        cat(glue::glue("[{par$prgmTag}]: Writing Probe GRS RDS={rds}...{RET}"))
      readr::write_rds(ret_grs, rds, compress="gz")
    }
    ret_cnt <- bsp %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_grs
}

loadBspFormatted_COVIC = function(bsp, src, col=NULL,fields=NULL, sort=TRUE,
                                  verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadBspFormatted_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (is.null(col)) {
      col <- cols(
        # bsp_key = col_character()
        prb_cgn = col_character(),
        prb_srd = col_character(),
        prb_add = col_integer(),
        prb_des = col_character(),
        prb_src = col_character(),
        
        bsp_seq = col_character(),
        bsp_tag = col_character(),
        bsp_chr = col_character(),
        bsp_beg = col_integer(),
        bsp_srd = col_character(),
        bsp_mis = col_integer(),
        
        bsp_ref = col_character(),
        bsp_gap = col_integer(),
        
        # bsp_str = col_character()
        bsp_hit0 = col_integer(),
        bsp_hit1 = col_integer(),
        bsp_hit2 = col_integer(),
        bsp_hit3 = col_integer(),
        bsp_hit4 = col_integer(),
        bsp_hit5 = col_integer()
      )
    }
    
    if (is.null(fields))
      fields=c("prb_cgn","prb_srd","prb_add","prb_des","prb_src")
    
    # Load BSP
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading BSP={bsp}...{RET}"))
    bsp_tib <- NULL
    bsp_tib <- readr::read_tsv(bsp,col_names=names(col$cols),col_types=col)
    
    # Join with source data::
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Joining source data...{RET}"))
    
    ret_tib <- dplyr::inner_join(src,bsp_tib, by=dplyr::all_of(fields))
    ret_tib <- ret_tib %>%
      dplyr::mutate(
        # Generate bisulfite converted Reference Sequences::
        bsp_refU=bscUs(bsp_ref,uc=TRUE),
        bsp_refM=bscMs(bsp_ref,uc=TRUE),
        bsp_refD=bscDs(bsp_ref,uc=TRUE),
        
        # Generate RevComp for Probe Sequence::
        prc_seq=revCmp(aln_seq),
        
        # Confirm Alignment Orientation from Probe Sequence/BSP Seq::
        prb_mat=dplyr::case_when(
          aln_seq==bsp_seq ~ "f",
          prc_seq==bsp_seq ~ "r",
          TRUE ~ NA_character_),
        
        # TBD:: For bsp_ref it should be bisulfite converted to U,M,D
        #  - The actual converted base should then be extracted to 
        #     provide both Ref and Bsc versions. 
        #  - Below only show Ref for now...
        
        # Need both Ref and Converted bsp_refUMD...
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                      Di-Nucleotide & Next Base:: Reference
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        # Extract Potential Next Base from Reference Sequence::
        NB_F1=stringr::str_sub(bsp_ref,-2,-2),
        NB_R1=stringr::str_sub(bsp_ref, 2, 2),
        NB_F2=stringr::str_sub(bsp_ref,-1,-1),
        NB_R2=stringr::str_sub(bsp_ref, 1, 1),
        
        # Extract Potential CG di-nucleotide from Reference Sequence::
        CG_F1=stringr::str_sub(bsp_ref,-4,-3),
        CG_R1=stringr::str_sub(bsp_ref, 3, 4),
        CG_F2=stringr::str_sub(bsp_ref,-3,-2),
        CG_R2=stringr::str_sub(bsp_ref, 2, 3),
        
        # Set Expected target Next Base from Alignment Orientation::
        bsp_nxb_ref=dplyr::case_when(
          prb_des=="M" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ NB_F1,
          prb_des=="M" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ NB_R1,
          
          prb_des=="U" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ NB_F1,
          prb_des=="U" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ NB_R1,
          
          prb_des=="2" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ NB_F2,
          prb_des=="2" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ NB_R2,
          TRUE ~ NA_character_
        ), # %>% stringr::str_to_upper(),
        
        # Set Expected target CG di-nucleotide from Alignment Orientation::
        bsp_din_ref=dplyr::case_when(
          prb_des=="M" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ CG_F1,
          prb_des=="M" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ CG_R1,
          
          prb_des=="U" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ CG_F1,
          prb_des=="U" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ CG_R1,
          
          prb_des=="2" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ CG_F2,
          prb_des=="2" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ CG_R2,
          TRUE ~ NA_character_
        ) %>% stringr::str_to_upper(),
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #             Di-Nucleotide & Next Base:: Bisulfite Converted
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        # Extract Potential Next Base from Reference Sequence::
        NB_F1M=stringr::str_sub(bsp_refM,-2,-2),
        NB_R1M=stringr::str_sub(bsp_refM, 2, 2),
        
        NB_F1U=stringr::str_sub(bsp_refU,-2,-2),
        NB_R1U=stringr::str_sub(bsp_refU, 2, 2),
        
        NB_F2D=stringr::str_sub(bsp_refD,-1,-1),
        NB_R2D=stringr::str_sub(bsp_refD, 1, 1),
        
        # Extract Potential CG di-nucleotide from Reference Sequence::
        CG_F1M=stringr::str_sub(bsp_refM,-4,-3),
        CG_R1M=stringr::str_sub(bsp_refM, 3, 4),
        
        CG_F1U=stringr::str_sub(bsp_refU,-4,-3),
        CG_R1U=stringr::str_sub(bsp_refU, 3, 4),
        
        CG_F2D=stringr::str_sub(bsp_refD,-3,-2),
        CG_R2D=stringr::str_sub(bsp_refD, 2, 3),
        
        # Set Expected target Next Base from Alignment Orientation::
        bsp_nxb_bsc=dplyr::case_when(
          prb_des=="M" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ NB_F1M,
          prb_des=="M" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ NB_R1M,
          
          prb_des=="U" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ NB_F1U,
          prb_des=="U" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ NB_R1U,
          
          prb_des=="2" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ NB_F2D,
          prb_des=="2" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ NB_R2D,
          TRUE ~ NA_character_
        ), # %>% stringr::str_to_upper(),
        
        # Set Expected target CG di-nucleotide from Alignment Orientation::
        bsp_din_bsc=dplyr::case_when(
          prb_des=="M" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ CG_F1M,
          prb_des=="M" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ CG_R1M,
          
          prb_des=="U" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ CG_F1U,
          prb_des=="U" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ CG_R1U,
          
          prb_des=="2" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ CG_F2D,
          prb_des=="2" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ CG_R2D,
          TRUE ~ NA_character_
        ) %>% stringr::str_to_upper(),
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                      Correct CG# Alignment Position::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        # Update Correct Genomic CG# Location based on alignment orientation::
        bsp_pos=dplyr::case_when(
          prb_des=="M" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ bsp_beg +48,
          prb_des=="M" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ bsp_beg + 0,
          
          prb_des=="U" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ bsp_beg +48,
          prb_des=="U" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ bsp_beg + 0,
          
          prb_des=="2" & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ bsp_beg +49,
          prb_des=="2" & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ bsp_beg - 1,
          TRUE ~ NA_real_
        ) %>% as.integer(),
        bsp_din_scr=dplyr::case_when(
          bsp_srd=="--" & bsp_din_ref=="CG" ~ 0,
          bsp_srd=="+-" & bsp_din_ref=="CG" ~ 1,
          bsp_srd=="-+" & stringr::str_starts(bsp_din_ref,"G") ~ 2,
          bsp_srd=="++" & stringr::str_ends(bsp_din_ref,"C") ~ 3,
          TRUE ~ 4) %>% as.integer()
      ) %>% 
      dplyr::group_by(prb_cgn,prb_srd,prb_des,prb_src) %>% 
      dplyr::mutate(Unique_ID=paste(prb_cgn,prb_srd,prb_add,prb_des,prb_src, dplyr::row_number(), sep="_")) %>%
      dplyr::ungroup()
    
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(bsp_chr, bsp_beg)
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} summary level 1={RET}"))
      sum1_tib <- ret_tib %>% 
        dplyr::select(prb_cgn,prb_srd, starts_with("CG"), 
                      bsp_ref,prb_des,bsp_srd,bsp_chr,bsp_beg) %>% 
        dplyr::arrange(prb_cgn) %>% 
        dplyr::group_by(prb_des,bsp_din_ref) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>% 
        dplyr::arrange(-Count)
      sum1_tib %>% print(n=base::nrow(sum1_tib))
      
      cat(glue::glue("[{funcTag}]:{tabsStr} summary level 2={RET}"))
      sum2_tib <- ret_tib %>% 
        dplyr::select(prb_cgn,prb_srd, starts_with("CG"), 
                      bsp_ref,prb_des,bsp_srd,bsp_chr,bsp_beg) %>% 
        dplyr::arrange(prb_cgn) %>% 
        dplyr::group_by(bsp_din_ref) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>% 
        dplyr::arrange(-Count)
      sum2_tib %>% print(n=base::nrow(sum2_tib))
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}{RET}"))
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
#                  Manifest Coordinate Intersection Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

impBuildSummary = function(tib,
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'impBuildSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- tib %>% 
      dplyr::distinct(U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                      chrom,chromStart,chromEnd,CGN_IMP,TB_IMP,CO_IMP, 
                      Probe_Type,Manifest,Unq_Cnt, .keep_all=TRUE) %>%
      dplyr::group_by(Manifest,Probe_Type,chrFlag,difFlag,cgnFlag,srcFlag,
                      Infinium_Design) %>%
      #                srcFlagU,srcFlagM) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} sum_tib({ret_cnt})={RET}"))
      ret_tib %>% print(n=base::nrow(ret_tib))
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# DELETE:: This is old and can probably be removed::
intersectCgnMap_COVIC = function(tsv, man, src=NULL, addBuild=FALSE,
                                 build, datDir, inpDir, runType, ver,
                                 verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'intersectCgnMap_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # TBD:: Make dir part of par$posDir
    #
    cgn_pos_fns <- paste(build,"21092020_improbe-designOutput.cgn-map-seq.cgn-sorted.tsv.gz", sep='-' )
    cgn_pos_pre <- paste(build,runType,ver, sep='-')
    imp_pos_tsv <- file.path(datDir, cgn_pos_fns)
    out_pos_tsv <- file.path(inpDir, paste(cgn_pos_pre,'seq48U_to_cgn.int.improbe-designOutput.cgn-map-seq.cgn-sorted.tsv.gz', sep='.'))
    
    if (!file.exists(out_pos_tsv)) {
      # Build better database::: Combine the two commands below::
      # gzip -dc /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.tsv.gz | cut -f 1,4,5,21-23,45,48 > /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.tsv
      # cat /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.tsv | perl -pe 's/TOP/T/; s/BOT/B/;' | gzip -c - > /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.tsv.gz
      
      int_pos_cmd <- glue::glue("gzip -dc {imp_pos_tsv} | join -t $'\t' -11 -22 - {tsv} | gzip -c - > {out_pos_tsv}")
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]: Running; cmd={int_pos_cmd}...{RET}"))
      system(int_pos_cmd)
    }
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: Loading Coordinates; out_pos_tsv={out_pos_tsv}...{RET}"))
    imp_pos_col <- cols(
      CGN_IMP = col_character(),
      chrom = col_character(),
      chromStart = col_integer(),
      FR_IMP = col_character(),
      TB_IMP = col_character(),
      CO_IMP = col_character(),
      NB_IMP = col_character(),
      Seq50U = col_character(),
      
      Seq48U = col_character(),
      TB_DB2 = col_character(),
      CO_DB2 = col_character(),
      U_DB2 = col_integer(),
      M_DB2 = col_integer()
    )
    imp_pos_tib <- 
      readr::read_tsv(out_pos_tsv, 
                      col_names=names(imp_pos_col$cols), 
                      col_types=imp_pos_col) %>% 
      # dplyr::distinct(CGN_IMP,chrom,chromStart,FR_IMP,TB_IMP,CO_IMP,NB_IMP,Seq50U) %>%
      dplyr::mutate(chrom=stringr::str_remove(chrom,'^chr'),
                    chrom=paste0('chr',chrom),
                    chromEnd=chromStart+1)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #             4.1 Match All Probes To CG# Database Coordinates::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- man %>%
      dplyr::inner_join(imp_pos_tib,
                        by=c("PRB1_U_MAT"="Seq50U"), 
                        suffix=c("_PQC", "_IMP")) %>% 
      dplyr::arrange(chrom,chromStart) %>%
      dplyr::select(chrom,chromStart,chromEnd, 
                    CGN_IMP,FR_IMP,TB_IMP,CO_IMP,NB_IMP,
                    PRB1_U_MAT,
                    dplyr::everything()) %>% 
      dplyr::add_count(U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                       chrom,chromStart,chromEnd,CGN_IMP,TB_IMP,CO_IMP, 
                       name="Unq_Cnt")
    
    if (addBuild) {
      if (build=="GRCh37") {
        ret_tib <- ret_tib %>%
          dplyr::mutate(
            begDif=chromStart-chromStart_GRCh37,
            endDif=chromEnd-chromEnd_GRCh37,
            absDif=base::abs(begDif),
            cgnFlag=dplyr::case_when(
              is.na(CGN_IMP) ~ -2,
              is.na(Seq_ID) ~ -1,
              CGN_IMP==Seq_ID ~ 0,
              TRUE ~ 1
            ),
            chrFlag=dplyr::case_when(
              is.na(chrom) ~ -2,
              is.na(chrom_GRCh37) ~ -1,
              chrom==chrom_GRCh37 ~ 0,
              TRUE ~ 1
            ),
            difFlag=dplyr::case_when(
              is.na(absDif) ~ -1,
              absDif==0 ~ 0,
              TRUE ~ 1
            )
          )
      } else if (build=="GRCh38") {
        ret_tib <- ret_tib %>%
          dplyr::mutate(
            begDif=chromStart-chromStart_GRCh38,
            endDif=chromEnd-chromEnd_GRCh38,
            absDif=base::abs(begDif),
            cgnFlag=dplyr::case_when(
              is.na(CGN_IMP) ~ -2,
              is.na(Seq_ID) ~ -1,
              CGN_IMP==Seq_ID ~ 0,
              TRUE ~ 1
            ),
            chrFlag=dplyr::case_when(
              is.na(chrom) ~ -2,
              is.na(chrom_GRCh38) ~ -1,
              chrom==chrom_GRCh38 ~ 0,
              TRUE ~ 1
            ),
            difFlag=dplyr::case_when(
              is.na(absDif) ~ -1,
              absDif==0 ~ 0,
              TRUE ~ 1
            )
          )
      } else {
        stop(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: Unsupported build type={build}!{RET}{RET}"))
        return(NULL)
      }
      
      if (build=="GRCh38" | build=="GRCh38") {
        ret_tib <- ret_tib %>%
          dplyr::mutate(
            begDif=as.integer(begDif),
            endDif=as.integer(endDif),
            absDif=as.integer(absDif),
            cgnFlag=as.integer(cgnFlag),
            chrFlag=as.integer(chrFlag),
            difFlag=as.integer(difFlag)
          )
      }
      
      if (!is.null(src)) {
        ret_tib <- ret_tib %>%
          dplyr::mutate(
            srcFlagU=dplyr::case_when(
              U %in% src$U ~ 0,
              TRUE ~ 1),
            srcFlagM=dplyr::case_when(
              is.na(M) ~ 0,
              M %in% src$M ~ 0,
              TRUE ~ 1),
            srcFlagU=as.integer(srcFlagU),
            srcFlagM=as.integer(srcFlagM),
            srcFlag=dplyr::case_when(
              srcFlagU==0 & srcFlagM==0 ~ 0,
              TRUE ~ 1
            ),
            srcFlag=as.integer(srcFlag)
          ) %>%
          dplyr::arrange(AlleleA_Probe_Sequence, chrFlag, difFlag, cgnFlag, absDif, 
                         srcFlagM, srcFlagU, srcFlagM)
      } else {
        ret_tib <- ret_tib %>%
          dplyr::arrange(AlleleA_Probe_Sequence, chrFlag, difFlag, cgnFlag, absDif)
      }
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt) {
      sum1_tib <- ret_tib %>% 
        dplyr::group_by(Manifest,Probe_Type) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      cat(glue::glue("[{funcTag}]: sum1_tib({ret_cnt})={RET}"))
      sum1_tib %>% print(n=base::nrow(sum1_tib))
      
      # For hg19 & hg38 this has nrows = 0
      if (addBuild) {
        if (build=="GRCh38" | build=="GRCh38") {
          mis_cnt <- ret_tib %>% dplyr::filter(Unq_Cnt != 1) %>% base::nrow()
          cat(glue::glue("[{funcTag}]: mis_cnt({ret_cnt})={mis_cnt}.{RET}"))
          
          sum2_tib <- ret_tib %>% 
            dplyr::distinct(U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                            chrom,chromStart,chromEnd,CGN_IMP,TB_IMP,CO_IMP, 
                            Probe_Type,Manifest,Unq_Cnt, .keep_all=TRUE) %>%
            dplyr::group_by(Manifest,Probe_Type,chrFlag,difFlag,cgnFlag,Infinium_Design) %>% 
            dplyr::summarise(Count=n(), .groups="drop")
          cat(glue::glue("[{funcTag}]: sum2_tib({ret_cnt})={RET}"))
          sum2_tib %>% print(n=base::nrow(sum2_tib))
        }
      }
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
#                   Standard Manifest Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

extractAddressRow_COVIC = function(tib, 
                                   tarAdd=NULL, remAdd=NULL, remVec=NULL, 
                                   oldNames=NULL, newNames=NULL,
                                   verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'extractAddressRow_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- tib
    
    # Remove NA fields::
    if (!is.null(remAdd)) {
      remAdd_sym <- rlang::sym(remAdd)
      ret_tib <- ret_tib %>% dplyr::filter(is.na(!!remAdd_sym))
    }
    
    # Select NOT NA fields::
    if (!is.null(tarAdd)) {
      tarAdd_sym <- rlang::sym(tarAdd)
      ret_tib <- ret_tib %>% dplyr::filter(!is.na(!!tarAdd_sym))
    }
    
    # Remove fields::
    if (!is.null(remVec)) {
      print(remVec)
      ret_tib <- ret_tib %>%
        dplyr::select(!dplyr::any_of(c(remAdd,remVec)))
    }
    
    # Rename fields::
    if (!is.null(oldNames) && ! is.null(newNames))
      # TBD:: Need to use all_of()
      ret_tib <- ret_tib %>% 
        rename_with(~ newNames[which(dplyr::all_of(oldNames) == .x)], 
                    .cols = dplyr::all_of(oldNames))
    #  rename_with(~ newNames[which(oldNames == .x)], .cols = oldNames)
    
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

manifestToAddress_COVIC = function(tib, 
                                   prbA="AlleleA_Probe_Sequence", 
                                   prbB="AlleleB_Probe_Sequence",
                                   verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'manifestToAddress_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr} prbA={prbA}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} prbB={prbB}.{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    prbA_sym <- rlang::sym(prbA)
    prbB_sym <- rlang::sym(prbB)
    
    tib <- tib %>% dplyr::filter(!is.na(U))
    
    ret_tib <- dplyr::bind_rows(
      
      extractAddressRow_COVIC(
        tib=tib, tarAdd="U",remAdd="M",
        remVec=c("AlleleB_Probe_Sequence"),
        oldNames=c("U","AlleleA_Probe_Sequence"),
        newNames=c("prb_add","ord_seq"),
        verbose=verbose+1,vt=vt+1,tc=tc+1,tt=pTracker) %>%
        dplyr::mutate(Probe_Design="2"),
      
      extractAddressRow_COVIC(
        tib=tib, tarAdd="M",
        remVec=c("M","AlleleB_Probe_Sequence"),
        oldNames=c("U","AlleleA_Probe_Sequence"),
        newNames=c("prb_add","ord_seq"),
        verbose=verbose+1,vt=vt+1,tc=tc+1,tt=pTracker) %>%
        dplyr::mutate(Probe_Design="U"),
      
      extractAddressRow_COVIC(
        tib=tib, tarAdd="M",
        remVec=c("U","AlleleA_Probe_Sequence"),
        oldNames=c("M","AlleleB_Probe_Sequence"),
        newNames=c("prb_add","ord_seq"),
        verbose=verbose+1,vt=vt+1,tc=tc+1,tt=pTracker) %>%
        dplyr::mutate(Probe_Design="M")
      
    ) %>%
      dplyr::mutate(prb_add=as.integer(prb_add),
                    ord_seq=stringr::str_to_upper(ord_seq),
                    aln_seq=ord_seq %>%
                      stringr::str_replace_all("R","A") %>%
                      stringr::str_replace_all("Y","T") ) %>%
      dplyr::rename(
        prb_des=Probe_Design,prb_src=Manifest) %>%
      dplyr::distinct(prb_add,ord_seq, .keep_all=TRUE) %>%
      # dplyr::arrange(prb_add) %>%
      # dplyr::select(Address,Probe_ID,Probe_Type,Probe_Design,Manifest,
      #               ord_seq,aln_seq, dplyr::everything()) %>%
      tidyr::separate(
        Probe_ID, into=c("prb_cgn","prb_srd"), sep="_") %>%
      dplyr::mutate(
        prb_srd=dplyr::case_when(
          is.na(prb_srd) ~ paste0("xxx",prb_des),
          TRUE ~ prb_srd),
        aln_key=paste(prb_cgn,prb_srd,prb_add,prb_des,prb_src, sep="_")
      )
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      ret_tib %>% print()
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
#                   Standard Manifest Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadCoreManifest_COVIC = function(datDir, manDir,
                                  verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadCoreManifest_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; manDir={manDir}...{RET}"))
  
  stopifnot(dir.exists(manDir))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                1.0 Load Sesame hg37/hg38 Manifest:: S4
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ses_hg37_tib <- 
      loadSesameManifest_COVIC(build="GRCh37", tag=TRUE, 
                               verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    ses_hg38_tib <-
      loadSesameManifest_COVIC(build="GRCh38", tag=TRUE,
                               verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_tib <- dplyr::inner_join(
      ses_hg37_tib,ses_hg38_tib, 
      by=c("Probe_ID","Probe_Type","U","M",
           "AlleleA_Probe_Sequence","AlleleB_Probe_Sequence")) %>%
      dplyr::mutate(Manifest="B4") %>%
      dplyr::select(Probe_ID,Probe_Type,U,M,
                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Manifest,
                    dplyr::everything())
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Sesame(hg37/hg38) manifest={RET}"))
      print(ret_tib)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       1.1 Load EPIC Manifest:: B2
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    epic_csv <- file.path(manDir, 'MethylationEPIC_v-1-0_B2.csv.gz')
    if (verbose>=vt+4)
      cat(glue::glue("[{funcTag}]:{tabsStr} EPIC(B2) manifest={epic_csv}{RET}"))
    
    epic_tib <- 
      loadManifestGenomeStudio(file = epic_csv, addSource = TRUE, 
                               normalize = TRUE, retType = "man", 
                               verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} EPIC(B2) manifest={RET}"))
      print(epic_tib)
    }
    
    epic_tib <- epic_tib %>%
      dplyr::rename(Probe_ID=IlmnID, 
                    chrom_GRCh37=Chromosome,
                    chromStart_GRCh37=Coordinate,
                    FR_GRCh37=Strand_FR,
                    nextBase_GRCh37=Next_Base,
                    AlleleA_Probe_Sequence=AlleleA_ProbeSeq,
                    AlleleB_Probe_Sequence=AlleleB_ProbeSeq,
                    U=AddressA_ID, M=AddressB_ID) %>%
      dplyr::mutate(U=as.integer(U),
                    M=as.integer(M),
                    chrom_GRCh37=paste0("chr",chrom_GRCh37),
                    chromEnd_GRCh37=chromStart_GRCh37+1,
                    Manifest="B2") %>%
      dplyr::select(Probe_ID,Probe_Type,U,M,
                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                    Manifest,chrom_GRCh37,chromStart_GRCh37,chromEnd_GRCh37,
                    FR_GRCh37,nextBase_GRCh37)
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} EPIC(B2) manifest(cleaned)={RET}"))
      print(epic_tib)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       1.2 Bind Manifests:: C0/B2
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- ret_tib %>% 
      dplyr::bind_rows(dplyr::anti_join(epic_tib, ret_tib, by="U") ) %>% 
      addSeq48U(field="AlleleA_Probe_Sequence") %>%
      dplyr::mutate(
        Seq_48U=dplyr::case_when( is.na(M) ~ Seq_48U_2,TRUE ~ Seq_48U_1 )
      ) %>%
      dplyr::select(-Seq_48U_1,-Seq_48U_2) %>%
      dplyr::select(Probe_ID,Probe_Type,U,M,
                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Seq_48U,
                    dplyr::everything())
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       1.X Load COVIC Manifest:: C0
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Load Original COVIC (EPIC_C0) manifest and compare to idat_join_tib
    # covic_c0_csv <- file.path(datDir, 'manifest/base/EPIC-C0.manifest.sesame-base.cpg-sorted.csv.gz')
    # covic_c0_tib <- readr::read_csv(covic_c0_csv) %>% 
    #   dplyr::mutate(M=as.integer(M), U=as.integer(U))
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Core Manifests({ret_cnt})={RET}{RET}"))
      print(ret_tib)
    }
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Core Summary({ret_cnt})={RET}{RET}"))
      ret_tib %>% 
        dplyr::group_by(Probe_Type,Manifest) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>% print()
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
#                        Sesame Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadSesameManifest_COVIC = function(build, tag=FALSE, del="_",
                                    verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadSesameManifest_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    genomeBuild_UCSC <- NA_character_
    if (build=='GRCm38') genomeBuild_UCSC <- 'mm10'
    if (build=='GRCh38') genomeBuild_UCSC <- 'hg38'
    if (build=='GRCh37') genomeBuild_UCSC <- 'hg19'
    if (build=='GRCh36') genomeBuild_UCSC <- 'hg18'
    
    if (is.na(genomeBuild_UCSC)) {
      stop(glue::glue("{RET}[{funcTag}]: genomeBuild_UCSC is null!{RET}{RET}"))
      return(ret_tib)
    }
    
    ses_fns <- paste('EPIC',genomeBuild_UCSC,'manifest', sep='.')
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading Sesame Manifest={ses_fns}...{RET}"))
    }
    
    ret_tib <- sesameData::sesameDataGet(ses_fns) %>%
      as.data.frame() %>% 
      tibble::rownames_to_column(var="Probe_ID") %>%
      tibble::as_tibble() %>%
      dplyr::rename(Probe_Type=probeType, 
                    U=address_A, M=address_B,
                    AlleleA_Probe_Sequence=ProbeSeq_A,
                    AlleleB_Probe_Sequence=ProbeSeq_B,
                    chrom=seqnames, chromStart=start, chromEnd=end, 
                    Design_Type=designType) %>% 
      dplyr::mutate(U=as.integer(U), M=as.integer(M),
                    chrom=as.character(chrom),
                    chromStart=as.integer(chromStart),
                    chromEnd=as.integer(chromEnd),
                    FR=dplyr::case_when(
                      strand=='+' ~ "F",
                      strand=='-' ~ "R",
                      TRUE ~ NA_character_
                    ),
                    AlleleB_Probe_Sequence=dplyr::case_when(
                      stringr::str_length(AlleleB_Probe_Sequence)==0 ~ NA_character_, 
                      TRUE ~ AlleleB_Probe_Sequence)
      ) %>%
      dplyr::distinct(U,M, .keep_all=TRUE) %>% 
      dplyr::select(Probe_ID,Probe_Type, 
                    U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, 
                    chrom,chromStart,chromEnd,
                    FR,nextBase, 
                    gene,gene_HGNC,MASK_mapping,MASK_general)
    
    if (tag) {
      ret_tib <- ret_tib %>% purrr::set_names(paste(names(.),build, sep=del))
      
      colnames(ret_tib)[1] = "Probe_ID"
      colnames(ret_tib)[2] = "Probe_Type"
      colnames(ret_tib)[3] = "U"
      colnames(ret_tib)[4] = "M"
      colnames(ret_tib)[5] = "AlleleA_Probe_Sequence"
      colnames(ret_tib)[6] = "AlleleB_Probe_Sequence"
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Core Manifests({ret_cnt})={RET}{RET}"))
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
#                           COVIC Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadAqpWorkflow_COVIC = function(ords, mats, aqps, man=NULL,
                                 verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadAqpWorkflow_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ords_vec <- ords %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
    mats_vec <- mats %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
    aqps_vec <- aqps %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
    
    # Check ORD Data::
    covic_ord_tib <- suppressMessages(suppressWarnings( 
      readr::read_csv(ords_vec[1], skip = 8, guess_max = 50000) )) %>%
      dplyr::rename(Order_ID=Assay_Design_Id) %>%
      dplyr::mutate(AlleleA_Probe_Sequence=stringr::str_to_upper(AlleleA_Probe_Sequence),
                    AlleleB_Probe_Sequence=stringr::str_to_upper(AlleleB_Probe_Sequence),
                    Design_TypeA=dplyr::case_when(!is.na(AlleleB_Probe_Sequence) ~ 'U', TRUE ~ '2' ), 
                    Design_TypeB=dplyr::case_when(!is.na(AlleleB_Probe_Sequence) ~ 'M', TRUE ~ '2' ) ) %>% 
      dplyr::select(-AlleleA_Probe_Id,-AlleleB_Probe_Id,-Normalization_Bin) %>%
      dplyr::arrange(Order_ID) %>%
      dplyr::distinct(AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, .keep_all=TRUE) %>%
      dplyr::mutate(Order_ID=stringr::str_remove(Order_ID, '_[0-9]+$') %>%
                      stringr::str_replace('_','-') %>% 
                      stringr::str_replace('II', '2') %>% 
                      stringr::str_replace('I', '1') %>%
                      stringr::str_remove_all('_') %>%
                      stringr::str_replace('-','_') )
    
    # Build flatten ORD table::
    #   - Non Unique = 10,537; NEW => 10,462
    covic_ord_tab <- dplyr::bind_rows(
      covic_ord_tib %>% dplyr::select(-AlleleB_Probe_Sequence,-Design_TypeB) %>% 
        dplyr::rename(Probe_Seq=AlleleA_Probe_Sequence,Design_Type=Design_TypeA),
      covic_ord_tib %>% dplyr::select(-AlleleA_Probe_Sequence,-Design_TypeA) %>% 
        dplyr::rename(Probe_Seq=AlleleB_Probe_Sequence,Design_Type=Design_TypeB),
    ) %>% dplyr::filter(!is.na(Probe_Seq) & !is.na(Design_Type)) %>% 
      dplyr::arrange(Order_ID) %>%
      dplyr::distinct(Probe_Seq,Design_Type, .keep_all=TRUE)    
    
    # Resolve Duplicates:: INCOMPLETE
    covic_ord_tab %>% 
      dplyr::add_count(Probe_Seq, name="Seq_Count") %>% 
      dplyr::filter(Seq_Count!=1) %>% 
      dplyr::arrange(-Seq_Count,Probe_Seq)
    
    # Check MAT Data::
    covic_mat_tib <- suppressMessages(suppressWarnings( 
      readr::read_tsv(mats_vec[1]) )) %>%
      dplyr::rename(Match_ID=probe_id,
                    Address=address_name,
                    Probe_Seq=sequence) %>%
      dplyr::mutate(Address=stringr::str_remove(Address, '^1') %>% 
                      stringr::str_remove('^0+') %>% as.integer(),
                    Probe_Seq=stringr::str_to_upper(Probe_Seq)) %>%
      dplyr::select(Match_ID,Probe_Seq,Address) %>%
      dplyr::arrange(Address) %>%
      dplyr::distinct(Address, .keep_all=TRUE) %>%
      dplyr::mutate(Match_ID=stringr::str_replace(Match_ID, '_','-') %>%
                      stringr::str_replace('II', '2') %>% 
                      stringr::str_replace('I', '1') %>%
                      stringr::str_remove_all('_') %>%
                      stringr::str_replace('-','_') )
    
    # Check PQC Data::
    #
    covic_aqp_tib <- loadPQC(file = aqps_vec[1], format = 'aqp', verbose = 4) %>% 
      dplyr::mutate(Address=stringr::str_remove(Address, '^1') %>% 
                      stringr::str_remove('^0+') %>% as.integer(),
                    Decode_Status=as.integer(Decode_Status)) %>%
      dplyr::select(Address,Decode_Status) %>%
      dplyr::arrange(Address) %>%
      dplyr::distinct(Address, .keep_all=TRUE)
    
    # Join MAT & AQP::
    covic_mat_aqp_raw_tab <- covic_mat_tib %>% 
      dplyr::inner_join(covic_aqp_tib, by="Address") %>% 
      dplyr::distinct(Probe_Seq,Address, .keep_all=TRUE)
    
    covic_mat_aqp_pas_tab <- covic_mat_aqp_raw_tab %>%
      dplyr::filter(!is.na(Decode_Status) & Decode_Status==0) %>%
      dplyr::distinct(Probe_Seq,Address, .keep_all=TRUE)
    
    # Join all data:: Decode Pass
    covic_ann_aqp_tab <- covic_ord_tab %>% 
      dplyr::inner_join(covic_mat_aqp_pas_tab, by="Probe_Seq") %>% 
      dplyr::select(-Decode_Status)
    
    # Rebuild Pairs:: covic_ord_tib U covic_ann_aqp_tab
    #
    ret_tib <- dplyr::bind_rows(
      covic_ord_tib %>% 
        dplyr::filter(is.na(AlleleB_Probe_Sequence)) %>%
        dplyr::inner_join(
          covic_ann_aqp_tab %>% 
            dplyr::rename(Order_ID_A=Order_ID, Match_ID_A=Match_ID, U=Address), 
          by=c("AlleleA_Probe_Sequence"="Probe_Seq",
               "Design_TypeA"="Design_Type")
        ) %>% dplyr::filter(!is.na(AlleleA_Probe_Sequence) & !is.na(U)),
      
      covic_ord_tib %>% 
        dplyr::filter(!is.na(AlleleB_Probe_Sequence)) %>%
        dplyr::inner_join(
          covic_ann_aqp_tab %>% 
            dplyr::rename(Order_ID_A=Order_ID, Match_ID_A=Match_ID, U=Address), 
          by=c("AlleleA_Probe_Sequence"="Probe_Seq",
               "Design_TypeA"="Design_Type")
        ) %>%
        dplyr::left_join(
          covic_ann_aqp_tab %>% 
            dplyr::rename(Order_ID_B=Order_ID, Match_ID_B=Match_ID, M=Address), 
          by=c("AlleleB_Probe_Sequence"="Probe_Seq",
               "Design_TypeB"="Design_Type")
        ) %>% dplyr::filter(!is.na(AlleleA_Probe_Sequence) & !is.na(U),
                            !is.na(AlleleB_Probe_Sequence) & !is.na(M))
    ) %>% dplyr::distinct(M,U, .keep_all=TRUE) %>% 
      dplyr::select(Order_ID,M,U,
                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                    Design_TypeA,Design_TypeB) %>%
      dplyr::rename(Probe_ID=Order_ID) %>%
      dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2))  %>%
      addSeq48U(field="AlleleA_Probe_Sequence") %>%
      dplyr::mutate(
        Seq_48U=dplyr::case_when( is.na(M) ~ Seq_48U_2,TRUE ~ Seq_48U_1 )
      ) %>%
      dplyr::select(Probe_ID,Probe_Type,U,M,
                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Seq_48U)
    
    if (!is.null(man)) ret_tib <- ret_tib %>% dplyr::mutate(Manifest=man)
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} AQP({ret_cnt})={RET}{RET}"))
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
#                        IDAT Samples Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# loadIdatAddress(prefix="/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R01C01")
loadIdatAddress = function(prefix,
                           verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadIdatAddress'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # idat_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R01C01"
    grn_dat <- illuminaio::readIDAT(file=paste(prefix,"Grn.idat.gz", sep="_"), what="all")
    ret_tib <- grn_dat$Quants %>% as.data.frame() %>% 
      tibble::rownames_to_column(var="Address") %>% 
      dplyr::mutate(Address=as.integer(Address)) %>%
      tibble::as_tibble() %>% 
      dplyr::filter(!is.na(Address)) %>% 
      dplyr::distinct(Address) %>% 
      dplyr::arrange(Address)
    
    if (FALSE) {
      idat1_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R01C01"
      idat1_tib <- prefixToIdat(prefix = idat1_prefix, verbose=verbose)
      
      # idat2_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R02C01"
      # idat2_tib <- prefixToIdat(prefix = idat2_prefix, verbose=verbose)
      # 
      # idat3_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R03C01"
      # idat3_tib <- prefixToIdat(prefix = idat3_prefix, verbose=verbose)
      # 
      # idat4_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R04C01"
      # idat4_tib <- prefixToIdat(prefix = idat4_prefix, verbose=verbose)
      # 
      # idat5_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R05C01"
      # idat5_tib <- prefixToIdat(prefix = idat5_prefix, verbose=verbose)
      # 
      # idat6_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R06C01"
      # idat6_tib <- prefixToIdat(prefix = idat6_prefix, verbose=verbose)
      # 
      # idat7_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R07C01"
      # idat7_tib <- prefixToIdat(prefix = idat7_prefix, verbose=verbose)
      
      idat8_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R08C01"
      idat8_tib <- prefixToIdat(prefix = idat8_prefix, verbose=verbose)
      
      # Combine all tangos::
      idats_tib <- 
        dplyr::select(idat1_tib, Address) %>% dplyr::mutate(Src1=1) %>% dplyr::full_join(
          # dplyr::select(idat2_tib, Address) %>% dplyr::mutate(Src2=2), by="Address" ) %>% dplyr::full_join(
          #   dplyr::select(idat3_tib, Address) %>% dplyr::mutate(Src3=3), by="Address" ) %>% dplyr::full_join(
          #     dplyr::select(idat4_tib, Address) %>% dplyr::mutate(Src4=4), by="Address" ) %>% dplyr::full_join(
          #       dplyr::select(idat5_tib, Address) %>% dplyr::mutate(Src5=5), by="Address" ) %>% dplyr::full_join(
          #         dplyr::select(idat6_tib, Address) %>% dplyr::mutate(Src6=6), by="Address" ) %>% dplyr::full_join(
          #           dplyr::select(idat7_tib, Address) %>% dplyr::mutate(Src7=7), by="Address" ) %>% dplyr::full_join(
          dplyr::select(idat8_tib, Address) %>% dplyr::mutate(Src8=8), by="Address" ) 
      
      unq_idats_tib <- idats_tib %>% dplyr::filter(!is.na(Address)) %>% 
        dplyr::distinct(Address) %>% dplyr::arrange(Address)
      unq_cnt <- unq_idats_tib %>% base::nrow()
      idats_sum_tib <- idats_tib %>% dplyr::summarise(
        Tot_Cnt=n(),
        S1=sum(is.na(Src1)),
        # S2=sum(is.na(Src2)),
        # S3=sum(is.na(Src3)),
        # S4=sum(is.na(Src4)),
        # S5=sum(is.na(Src5)),
        # S6=sum(is.na(Src6)),
        # S7=sum(is.na(Src7)),
        S8=sum(is.na(Src8)),
      ) %>% dplyr::mutate(Unq_Cnt=unq_cnt, Dif_Cnt=Tot_Cnt-Unq_Cnt)
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} idat({ret_cnt})={RET}{RET}"))
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

if (FALSE) {
  
  #
  # NOW JUST Missing the ones below::
  #
  # add_cor_all_tib %>% 
  #   dplyr::anti_join(add_fin_all_tib, by="prb_add") %>% 
  #   dplyr::left_join(add_fin_all_tib, by=c("ord_seq"="prb_ord_seq") )
  #
  # Found earlier::
  #  add_cor_mis_tib %>% dplyr::inner_join(add_imp_all2_tib, by="prb_add") %>% dplyr::distinct(ord_seq, .keep_all=TRUE)
  #
  # No Gene Annotation...
  #  add_cor_mis_tib %>% dplyr::inner_join(add_imp_all_ncbi_tib, by="prb_add") %>% dplyr::distinct(ord_seq, .keep_all=TRUE)
  #
  # Missing Probes Found::
  #  add_imp_all2_tib %>% dplyr::filter(prb_add %in% add_cor_mis_tib$prb_add) %>% dplyr::filter(prb_des != "M") %>% dplyr::distinct(prb_add, .keep_all=TRUE)
  #
  
  
  if (FALSE) {
    
    # Raw::
    add_raw_ses_tib %>% dplyr::group_by(ses_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    add_raw_all_tib %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    
    add_raw_ses_tib %>% dplyr::distinct(ses_add, .keep_all=TRUE) %>% dplyr::group_by(ses_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    add_raw_all_tib %>% dplyr::distinct(prb_add, .keep_all=TRUE) %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    
    # Imp::
    add_imp_ses2_tib %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    add_imp_all2_tib %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    
    add_imp_ses2_tib %>% dplyr::distinct(prb_add, .keep_all=TRUE) %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    add_imp_all2_tib %>% dplyr::distinct(prb_add, .keep_all=TRUE) %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    
    # Should be only B4
    add_imp_ses2_tib %>% dplyr::anti_join(add_imp_all2_tib, by="prb_add") %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    add_imp_ses2_tib %>% dplyr::anti_join(add_imp_all2_tib, by="prb_ord") %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    
    add_imp_ses2_tib %>% dplyr::filter(prb_add %in% add_imp_all2_tib$prb_add) %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    add_imp_ses2_tib %>% dplyr::filter(prb_ord %in% add_imp_all2_tib$prb_ord) %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    
    # Should be only C0
    add_imp_all2_tib %>% dplyr::anti_join(add_imp_ses2_tib, by="prb_add") %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    add_imp_all2_tib %>% dplyr::anti_join(add_imp_ses2_tib, by="prb_ord") %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    
    add_imp_all2_tib %>% dplyr::filter(prb_add %in% add_imp_ses2_tib$prb_add) %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    add_imp_all2_tib %>% dplyr::filter(prb_ord %in% add_imp_ses2_tib$prb_ord) %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    
    # Gene Annotated
    add_imp_ses_ncbi_tib %>% dplyr::distinct(prb_add, .keep_all=TRUE) %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    add_imp_all_ncbi_tib %>% dplyr::distinct(prb_add, .keep_all=TRUE) %>% dplyr::group_by(prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
    
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                Reconcile Sesame/Previous Manifest with 
  #                   New Address Alignment/Annotation::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Orientations Investigation::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Need to add this step back into original source data file...
  # man_imp_int_tib <- man_imp_int_tib %>% 
  #   dplyr::mutate(prb_aln_seq2U=stringr::str_remove(prb_aln_seq, "^[A-Za-z]"), 
  #                 imp_seq2U=stringr::str_remove(imp_seq1U, "[A-Za-z]$"))
  # 
  # man_add_bsp_tib %>% dplyr::select(Unique_ID,prb_add) %>%
  #   dplyr::inner_join(man_imp_int_tib, by=c("Unique_ID"="Seq_ID") ) %>%
  #   dplyr::filter(prb_aln_seq==imp_seq1U | prb_aln_seq2U==imp_seq2U) %>%
  #   dplyr::mutate(CO=stringr::str_remove(prb_srd, "^[A-Z][A-Z]") %>%
  #                   stringr::str_remove("[0-9]+$"),
  #                 FR=stringr::str_remove(prb_srd, "[A-Z][A-Z][0-9]+$")) %>%
  #   dplyr::group_by(CO,FR,bsp_srd,prb_des,bsp_din) %>%
  #   dplyr::summarise(Count=n(), .groups="drop") %>% print(n=1000)
  # 
  # BSMAP Language::
  #
  #  C = [--] CG
  #  C = [+-] CG
  #  O = [-+] G*
  #  O = [++] *C
  #
  
  # Code for Scoring Sesame/Previous Gene Matching::
  #
  # dplyr::mutate(
  #   gene_mat_scr=dplyr::case_when(
  #     is.na(ses_gene) ~ 2,
  #     is.na(ses_HGNC) ~ 3,
  #     is.na(Gene_NCBI_Simple) ~ 4,
  #     stringr::str_detect(ses_gene,Gene_NCBI_Simple) ~ 0,
  #     stringr::str_detect(ses_HGNC,Gene_NCBI_Simple) ~ 1,
  #     TRUE ~ 5
  #   )
  # )
  
}


# End of file
