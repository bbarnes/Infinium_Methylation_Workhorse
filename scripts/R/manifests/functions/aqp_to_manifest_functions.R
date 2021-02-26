
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
#                             BSMAP Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadBspFormatted = function(bsp, src, col=NULL,fields=NULL, sort=TRUE,
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadBspFormatted'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (is.null(col)) {
      col <- cols(
        # bsp_key = col_character()
        prb_cgn = col_character(),
        # prb_srd = col_character(),
        # prb_add = col_integer(),
        prb_des = col_character(),
        # prb_src = col_character(),
        
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
      fields=c("prb_cgn","prb_des")
    #  fields=c("prb_cgn","prb_srd","prb_add","prb_des","prb_src")
    
    # Load BSP
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading BSP={bsp}...{RET}"))
    bsp_tib <- NULL
    bsp_tib <- readr::read_tsv(bsp,col_names=names(col$cols),col_types=col)
    
    print(bsp_tib)
    print(fields)
    
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
        
        #
        # TBD:: Need to understand how bsp_pos can become NA_real !!!
        #
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
      dplyr::group_by(prb_cgn,prb_des) %>% 
      dplyr::mutate(Unique_ID=paste(prb_cgn,prb_des, dplyr::row_number(), sep="_")) %>%
      # dplyr::group_by(prb_cgn,prb_srd,prb_des,prb_src) %>% 
      # dplyr::mutate(Unique_ID=paste(prb_cgn,prb_srd,prb_add,prb_des,prb_src, dplyr::row_number(), sep="_")) %>%
      dplyr::ungroup()
    
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(bsp_chr, bsp_beg)
    
    if (FALSE && verbose>=vt+4) {
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

bspToGenomicRegion = function(bsp, rds=NULL,
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'bspToGenomicRegion'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (!is.null(rds)) {
      dir <- base::dirname(rds)
      if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    }
    
    #
    # TBD:: Temp Fix for bsp_pos == NA !!!
    #
    
    bsp <- bsp %>% dplyr::filter(!is.na(bsp_pos))
    
    ret_grs <- 
      GRanges(
        seqnames = Rle(paste0("chr",bsp$bsp_chr)),
        strand=Rle(stringr::str_sub( bsp$bsp_srd, 1,1 ) ),
        
        ord_id  = bsp$ord_id,
        
        prb_add = bsp$prb_add,
        prb_cgn = bsp$prb_cgn,
        # prb_srd = bsp$prb_srd,
        prb_des = bsp$prb_des,
        # prb_src = bsp$prb_src,
        
        prb_ord_seq  = bsp$ord_seq,
        prb_aln_50U  = bsp$ord_seq %>%
          stringr::str_replace_all("R","A") %>% # A/G
          stringr::str_replace_all("Y","T") %>% # C/T
          
          stringr::str_replace_all("S","C") %>% # G/C
          stringr::str_replace_all("W","A") %>% # A/T
          stringr::str_replace_all("K","T") %>% # G/T
          stringr::str_replace_all("M","A") %>% # A/C
          
          stringr::str_replace_all("B","T") %>% # C/G/T
          stringr::str_replace_all("D","A") %>% # A/G/T
          stringr::str_replace_all("H","A") %>% # A/C/T
          stringr::str_replace_all("V","A") %>% # A/C/G
          
          stringr::str_replace_all("N","A"), # A/C/T/G
        # prb_aln_49U  = stringr::str_sub(bsp$aln_seq,2),
        # prb_aln_50M  = stringr::str_replace_all(bsp$ord_seq,"G","A"),
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Manifest File Methods:: Formatting
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

format_ORD = function(tib,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_ORD'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    tmp_tib <- tib %>%
      dplyr::select(-AlleleA_Probe_Id, -AlleleB_Probe_Id, -Normalization_Bin) %>% 
      dplyr::rename(PrbA=AlleleA_Probe_Sequence,
                    PrbB=AlleleB_Probe_Sequence) %>%
      dplyr::mutate(
        DesA=dplyr::case_when(
          !is.na(PrbA) & !is.na(PrbB) ~ "U",
          !is.na(PrbA) &  is.na(PrbB) ~ "2",
          TRUE ~ NA_character_
        ),
        DesB=dplyr::case_when(
          !is.na(PrbA) & !is.na(PrbB) ~ "M",
          TRUE ~ NA_character_
        )
      )
    
    add_col <- c("ord_id","ord_des","ord_aqp","ord_seq")
    ret_tib <- dplyr::bind_rows(
      tmp_tib %>% dplyr::filter(DesB=="M") %>% 
        dplyr::select(Assay_Design_Id,DesB,ord_aqp,PrbB) %>%
        purrr::set_names(add_col),
      tmp_tib %>% dplyr::filter(DesA=="U") %>% 
        dplyr::select(Assay_Design_Id,DesA,ord_aqp,PrbA) %>%
        purrr::set_names(add_col),
      tmp_tib %>% dplyr::filter(DesA=="2") %>% 
        dplyr::select(Assay_Design_Id,DesA,ord_aqp,PrbA) %>%
        purrr::set_names(add_col)
    ) %>%
      dplyr::mutate(ord_seq=stringr::str_to_upper(ord_seq)) %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) ) %>%
      dplyr::arrange(-ord_aqp) %>%
      dplyr::distinct(ord_seq, .keep_all=TRUE)
    # dplyr::arrange(ord_id)
    
    ret_tib <- ret_tib %>% 
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib=({ret_cnt})={RET}"))
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

format_MAT = function(tib, trim=TRUE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_MAT'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- tib
    if (trim) {
      trim_cnt <- ret_tib %>% 
        dplyr::filter(!stringr::str_starts(Address, '1')) %>% base::nrow()
      if (trim_cnt==0) {
        ret_tib <- ret_tib %>% 
          dplyr::mutate(Address=as.numeric(
            stringr::str_remove(stringr::str_remove(Address, '^1'), '^0+')) )
      } else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: Attempting to trim new tango fomrat, ",
                        "but format doesn't match; trim_cnt={trim_cnt}!!!{RET}{RET}"))
        return(NULL)
      }
    }
    
    ret_tib <- ret_tib %>% 
      dplyr::mutate(tan_seq=stringr::str_sub(Sequence,1,45) %>%
                      stringr::str_to_upper(),
                    mat_seq=stringr::str_sub(Sequence,46)  %>%
                      stringr::str_to_upper() ) %>%
      dplyr::filter(!is.na(Address)) %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) ) %>%
      dplyr::arrange(Address,-mat_aqp,mat_seq) %>% 
      dplyr::rename(mat_add=Address) %>%
      dplyr::select(mat_add,mat_aqp,mat_seq,
                    dplyr::everything()) %>% 
      dplyr::distinct(mat_add, .keep_all=TRUE)
    
    ret_tib <- ret_tib %>% 
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib=({ret_cnt})={RET}"))
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

format_AQP = function(tib, sort=TRUE,filt=TRUE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_AQP'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # Determine if its AQP & PQC::
    #
    ret_tib <- tib %>% 
      dplyr::mutate(Address=as.numeric(
        stringr::str_remove(stringr::str_remove(Address, '^1'), '^0+')) )
    
    if (c("Decode_Status") %in% names(tib)) {
      ret_tib <- ret_tib %>% 
        dplyr::select(Address,Decode_Status,aqp_idx) %>%
        dplyr::rename(aqp_add=Address,
                      aqp_val=Decode_Status) %>%
        dplyr::mutate(aqp_src="AQP")
    } else if (c("Status") %in% names(tib)) {
      ret_tib <- ret_tib %>% 
        dplyr::select(Address,Status,aqp_idx) %>%
        dplyr::rename(aqp_add=Address,
                      aqp_val=Status) %>%
        dplyr::mutate(aqp_src="PQC")
    } else {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Unable to match AQP/PQC!{RET}{RET}"))
      return(NULL)
    }
    pre_cnt <- ret_tib %>% base::nrow()
    
    ret_tib <- ret_tib %>% 
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(-aqp_idx)
    if (filt) ret_tib <- ret_tib %>% dplyr::filter(aqp_val==0) %>%
      dplyr::distinct(aqp_add, .keep_all=TRUE)
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib=({ret_cnt}/{pre_cnt})={RET}"))
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
#                     Manifest File Methods:: Loading
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_manifestBuildFile = function(file, field="Assay_Design_Id", cols=NULL,
                                  n_max=50, guess=1000,
                                  verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_manifestBuildFile'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    skip_cnt <- get_headerSkipCount(
      file=file, field=field, n_max=n_max,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} skip_cnt={skip_cnt}.{RET}"))
    
    file_type <- get_fileSuffix(file)
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} file_type={file_type}.{RET}"))
    
    if (file_type=="csv") {
      ret_tib <- suppressMessages(suppressWarnings(
        readr::read_csv(file, skip=skip_cnt, guess_max=guess, skip_empty_rows=FALSE)))
    } else if (file_type=="tsv" || file_type=="txt") {
      ret_tib <- suppressMessages(suppressWarnings(
        readr::read_tsv(file, skip=skip_cnt, guess_max=guess, skip_empty_rows=FALSE)))
    } else {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Unsupporte fileType={file_type}!!!{RET}{RET}"))
      return(NULL)
    }
    if (!is.null(cols)) {
      cur_len <- length(names(ret_tib))
      col_len <- length(cols)
      if (cur_len==col_len) {
        ret_tib <- ret_tib %>% purrr::set_names(cols)
      } else {
        stop(glue::glue("{RET}[{funcTag}]:ERROR: Unequal Column Lengths {cur_len} != {col_len} !!!{RET}{RET}"))
        return(NULL)
      }
    }
    
    ret_tib <- ret_tib %>% 
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib=({ret_cnt})={RET}"))
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
#                         File Format Predction Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

get_headerSkipCount = function(file, field="Assay_Design_Id", n_max=50,
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'get_headerSkipCount'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; {RET}",
                   "file='{file}'{RET}",
                   "field='{field}'{RET}") )
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- 
      readr::read_lines(file=file, n_max=n_max, skip_empty_rows=FALSE) %>% 
      tibble::as_tibble() %>% dplyr::mutate(value=stringr::str_remove(value, "^\\s+") )
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      print(ret_tib, n=n_max)
    }
    
    ret_cnt <- ret_tib %>%
      dplyr::mutate(
        Row_Num=dplyr::row_number()) %>% 
      dplyr::filter(stringr::str_starts(value,!!field)) %>% 
      # head(n=1) %>% 
      tail(n=1) %>% 
      dplyr::pull(Row_Num) - 1
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_cnt
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Build Improbe Position Genomic Region:: RDS
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Example Build Commands::
#
# cur_pos_db_grs <- write_impTopGenomeGRS(genBuild="GRCm10", verbose=10, tt=pTracker)
# cur_pos_db_grs <- write_impTopGenomeGRS(genBuild="GRCh37", verbose=10, tt=pTracker)
# cur_pos_db_grs <- write_impTopGenomeGRS(genBuild="GRCh38", verbose=10, tt=pTracker)
write_impTopGenomeGRS = function(genBuild,
                              impDir="/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designInput",
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'write_impTopGenomeGRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    imp_pos_col <- NULL
    imp_pos_col <- cols(
      imp_chr = col_integer(),
      imp_pos = col_integer(),
      imp_cgn = col_integer(),
      imp_srd = col_integer(),
      imp_top = col_character()
    )
    print(imp_pos_col)
    
    pos_tsv <- ""
    pos_rds <- ""
    
    pos_tsv <- file.path(impDir, paste0(genBuild,'.cgn-pos-fwd.base.tsv.gz') )
    pos_rds <- file.path(impDir, paste0(genBuild,'.cgn-pos-fwd.base.rds') )
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading TSV={pos_tsv}...{RET}"))
    
    pos_tib <- NULL
    pos_tib <- readr::read_tsv(pos_tsv, 
                               col_names=names(imp_pos_col$cols),
                               col_types=imp_pos_col)
    
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr} pos_tib={RET}"))
      print(pos_tib)
    }
    ret_cnt <- pos_tib %>% base::nrow()
    
    if (FALSE) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Building GRS..{RET}"))
      
      ret_grs <- NULL
      ret_grs <- 
        GRanges(
          seqnames = Rle(paste0("chr",pos_tib$imp_chr)),
          # strand=Rle(stringr::str_sub( bsp$bsp_srd, 1,1 ) ),
          
          imp_cg=pos_tib$imp_cgn,
          imp_fr=pos_tib$imp_fr,
          imp_tb=pos_tib$imp_tb,
          imp_co=pos_tib$imp_co,
          
          IRanges(start=pos_tib$imp_pos,
                  width=2,
                  # end=pos_tib$Coordinate+1,
                  names=paste(pos_tib$imp_cgn,
                              pos_tib$imp_chr,
                              pos_tib$imp_pos,
                              sep="_")
          )
        )
      
      if (verbose>=vt) {
        cat(glue::glue("[{funcTag}]:{tabsStr} ret_grs={RET}"))
        print(ret_grs)
      }
      
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing RDS={pos_rds}...{RET}"))
      readr::write_rds(ret_grs, pos_rds, compress="gz")
      
      ret_cnt <- ret_grs %>% length()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  # ret_grs
  pos_tib
}

# Example Build Commands::
#
# cur_pos_db_grs <- write_impGenomeGRS(genBuild="GRCm10", verbose=10, tt=pTracker)
# cur_pos_db_grs <- write_impGenomeGRS(genBuild="GRCh37", verbose=10, tt=pTracker)
# cur_pos_db_grs <- write_impGenomeGRS(genBuild="GRCh38", verbose=10, tt=pTracker)
write_impGenomeGRS = function(genBuild,
                              impDir="/Users/bretbarnes/Documents/data/improbe/designOutput_21092020",
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'write_impGenomeGRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    imp_pos_col <- NULL
    imp_pos_col <- cols(
      imp_cgn = col_character(),
      imp_chr = col_character(),
      imp_pos = col_integer(),
      
      imp_fr  = col_character(),
      imp_tb  = col_character(),
      imp_co  = col_character()
    )
    
    pos_tsv <- ""
    pos_rds <- ""
    
    pos_tsv <- file.path(impDir, paste0(genBuild,'-21092020_improbe-designOutput.cgn-pos-srd.tsv.gz') )
    pos_rds <- file.path(impDir, paste0(genBuild,'-21092020_improbe-designOutput.cgn-pos-srd.rds') )

    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading TSV={pos_tsv}...{RET}"))
    
    pos_tib <- NULL
    pos_tib <- readr::read_tsv(pos_tsv, 
                               col_names=names(imp_pos_col$cols),
                               col_types=imp_pos_col)
    
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr} pos_tib={RET}"))
      print(pos_tib)
    }
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Building GRS..{RET}"))
    
    
    ret_grs <- NULL
    ret_grs <- 
      GRanges(
        seqnames = Rle(paste0("chr",pos_tib$imp_chr)),
        # strand=Rle(stringr::str_sub( bsp$bsp_srd, 1,1 ) ),
        
        imp_cg=pos_tib$imp_cgn,
        imp_fr=pos_tib$imp_fr,
        imp_tb=pos_tib$imp_tb,
        imp_co=pos_tib$imp_co,
        
        IRanges(start=pos_tib$imp_pos,
                width=2,
                # end=pos_tib$Coordinate+1,
                names=paste(pos_tib$imp_cgn,
                            pos_tib$imp_chr,
                            pos_tib$imp_pos,
                            sep="_")
        )
      )
    
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_grs={RET}"))
      print(ret_grs)
    }

    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Writing RDS={pos_rds}...{RET}"))
    readr::write_rds(ret_grs, pos_rds, compress="gz")
    
    ret_cnt <- ret_grs %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_grs
}

# End of file
