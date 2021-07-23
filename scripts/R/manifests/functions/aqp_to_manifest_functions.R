
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
BNG <- "|"
BRK <- paste0("# ",
              paste(rep("-----",6),collapse=" "),"|",
              paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='template_func') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Function Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   funcTag={funcTag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          AQP to Sesame Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

aqp_to_man = function(tib,
                      verbose=0,vt=3,tc=1,tt=NULL,
                      funcTag='aqp_to_man') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # TBD::
    #
    #  * Investigate BSP alignments::
    #    - Calculate number of hits good/bad
    #    - Calculate extension/color distribution
    #    - Extract BSC Top/Probe-Design
    #
    
    # Workflow::
    #
    #  - Split 2,U,M
    #  - Join U/M
    #  - 
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

aqp_to_sesame1 = function(tib, isMU=FALSE, retData=FALSE,
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='aqp_to_sesame1') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  stime <- base::system.time({
    
    #
    # Infinium I::
    #
    lab_tib <- tib
    
    if (!opt$multi_unique)
      # lab_tib <- lab_tib %>%
      # dplyr::filter(Ord_Din=="cg" & 
      #                 Bsp_Din_Ref_U=="CG" & Bsp_Din_Bsc_U=="TG" &
      #                 Bsp_Din_Ref_M=="CG" & Bsp_Din_Bsc_M=="CG")
      
      lab_tib <- lab_tib %>%
        dplyr::distinct(Ord_Key,Address_U,Address_M,Cgn_Str,Strand_TB_U,Bsp_CO_U,Ord_Prb_U,Ord_Prb_M) %>%
        dplyr::mutate(
          Cgn_Info=paste0(Strand_TB_U,Bsp_CO_U,"1"),
          Probe_ID=paste(Cgn_Str,Cgn_Info, sep="_")
        ) %>%
        dplyr::group_by(Probe_ID) %>%
        dplyr::mutate(Rep_Num=dplyr::row_number(),
                      Probe_ID=paste0(Probe_ID,Rep_Num)) %>%
        dplyr::ungroup() %>% 
        dplyr::select(Ord_Key,Probe_ID,Rep_Num,
                      Address_U,Ord_Prb_U,
                      Address_M,Ord_Prb_M)
    
    cor_tib <- tib %>% 
      dplyr::mutate(
        Name=Cgn_Str,
        U=Address_U,
        AlleleA_ProbeSeq=Ord_Prb_U,
        M=Address_M,
        AlleleB_ProbeSeq=Ord_Prb_M,
        Next_Base=Next_Base_U,
        Color_Channel=dplyr::case_when(
          Next_Base_U=="C" | Next_Base_U=="G" ~ "Grn",
          Next_Base_U=="A" | Next_Base_U=="T" ~ "Red",
          TRUE ~ NA_character_
        ),
        Col=stringr::str_sub(Color_Channel, 1,1),
        Probe_Type=Ord_Din,
        Strand_FR=Bsp_FR_U,
        Strand_TB=Strand_TB_U,
        Strand_CO=Bsp_CO_U,
        Infinium_Design_Type="1",
        CHR=Bsp_Chr,
        MAPINFO=Bsp_Pos,
        Species=opt$Species,
        Genome_Build=opt$genBuild,
        Source_Seq=Probe_Seq_T_U,
        Underlying_CpG_Count=Cpg_Cnt,
        Underlying_CpG_Min_Dist=Cpg_Dis_M,
        Alt_Cgn_Count_U=Alt_Cgn_Cnt_U,
        Alt_Cgn_Count_M=Alt_Cgn_Cnt_M
      ) %>%
      dplyr::select(Ord_Key,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag_U,Bsp_Tag_M)
    
    all_tib <- cor_tib %>%
      dplyr::inner_join(lab_tib,
                        by=c("Ord_Key",
                             "U"="Address_U","M"="Address_M",
                             "AlleleA_ProbeSeq"="Ord_Prb_U",
                             "AlleleB_ProbeSeq"="Ord_Prb_M") ) %>%
      dplyr::select(Probe_ID,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,Rep_Num,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag_U,Bsp_Tag_M) %>%
      clean_tibble()
    
    ret_tib <- dplyr::right_join(cor_tib,lab_tib,
                                 by=c("Ord_Key",
                                      "U"="Address_U","M"="Address_M",
                                      "AlleleA_ProbeSeq"="Ord_Prb_U",
                                      "AlleleB_ProbeSeq"="Ord_Prb_M") )
    if (!opt$multi_unique)
      ret_tib <- ret_tib %>%
      dplyr::mutate(
        CHR=dplyr::case_when(
          Bsp_Tag_U=="UM" && Bsp_Tag_M=="UM" ~ CHR,
          TRUE ~ "0"
        ),
        MAPINFO=dplyr::case_when(
          Bsp_Tag_U=="UM" && Bsp_Tag_M=="UM" ~ MAPINFO,
          TRUE ~ as.integer(0)
        ),
        Strand_FR=dplyr::case_when(
          Bsp_Tag_U=="UM" && Bsp_Tag_M=="UM" ~ Strand_FR,
          TRUE ~ "0"
        )
      )
    
    ret_tib <- ret_tib %>%
      dplyr::select(Probe_ID,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,Rep_Num,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag_U,Bsp_Tag_M) %>%
      dplyr::distinct() %>%
      dplyr::distinct(Probe_ID,Name,
                      U,AlleleA_ProbeSeq,
                      M,AlleleB_ProbeSeq,
                      # Next_Base, 
                      .keep_all = TRUE) %>%
      dplyr::distinct(U,M, .keep_all = TRUE) %>%
      clean_tibble()
    
    # lab_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "rs"))
    # all_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "rs"))
    # ret_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "rs"))
    
    # These should be zero::
    col_cnt1 <- ret_tib %>% dplyr::filter(is.na(Col)) %>% base::nrow()
    col_cnt2 <- ret_tib %>% dplyr::filter(is.na(Color_Channel)) %>% base::nrow()
    max_repn <- max(ret_tib$Rep_Num)
    cat(glue::glue("[{funcTag}]:{tabsStr} col_cnt1={col_cnt1}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} col_cnt2={col_cnt2}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} max_repn={max_repn}{RET}"))
    
    ret_tib %>% dplyr::group_by(Col) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print()
    ret_tib %>% dplyr::group_by(Color_Channel) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print()
    ret_tib %>% dplyr::filter(stringr::str_starts(Probe_ID,"rs")) %>% print()
    
    if (retData) ret_dat$all_tib <- all_tib
    if (retData) ret_dat$ses_tib <- ret_tib
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
}

aqp_to_sesame2 = function(tib, isMU=FALSE, retData=FALSE,
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='aqp_to_sesame2') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  stime <- base::system.time({
    
    # Ord_Din_sym <- rlang::sym(Ord_Din)
    
    lab_tib <- tib
    # if (!isMU)
    #   lab_tib <- lab_tib %>%
    #     dplyr::filter(Bsp_Din_Ref=="CG" & Bsp_Din_Bsc=="YG")
    
    lab_tib <- lab_tib %>%
      dplyr::distinct(Ord_Key,Address,Cgn_Str,Bsp_CO,Ord_Prb, .keep_all = TRUE) %>% 
      dplyr::select(Ord_Key,Address,Cgn_Str,Strand_TB,Bsp_CO,Ord_Prb) %>%
      dplyr::mutate(
        Cgn_Info=paste0(Strand_TB,Bsp_CO,"2"),
        Probe_ID=paste(Cgn_Str,Cgn_Info, sep="_")
      ) %>%
      dplyr::group_by(Probe_ID) %>%
      dplyr::mutate(Rep_Num=dplyr::row_number(),
                    Probe_ID=paste0(Probe_ID,Rep_Num)) %>%
      dplyr::ungroup() %>% 
      dplyr::select(Ord_Key,Probe_ID,Rep_Num,
                    Address,Ord_Prb)
    
    cor_tib <- tib %>% 
      dplyr::mutate(
        Name=Cgn_Str,
        U=Address,
        AlleleA_ProbeSeq=Ord_Prb,
        M="",
        AlleleB_ProbeSeq="",
        Color_Channel="Both",
        Col="",
        Probe_Type=Ord_Din,
        Strand_FR=Bsp_FR,
        Strand_CO=Bsp_CO,
        Infinium_Design_Type="2",
        CHR=Bsp_Chr,
        MAPINFO=Bsp_Pos,
        Species=opt$Species,
        Genome_Build=opt$genBuild,
        Source_Seq=Probe_Seq_T,
        Underlying_CpG_Count=Cpg_Cnt,
        Underlying_CpG_Min_Dist=Cpg_Dis,
        Alt_Cgn_Count_U=Alt_Cgn_Cnt,
        Alt_Cgn_Count_M=0
      ) %>%
      dplyr::select(Ord_Key,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag) %>%
      dplyr::rename(Bsp_Tag_U=Bsp_Tag) %>%
      dplyr::mutate(Bsp_Tag_M="")
    
    all_tib <- cor_tib %>%
      dplyr::inner_join(lab_tib,
                        by=c("Ord_Key","U"="Address",
                             "AlleleA_ProbeSeq"="Ord_Prb") ) %>%
      dplyr::select(Probe_ID,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,Rep_Num,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag_U,Bsp_Tag_M) %>%
      clean_tibble()
    
    ses_tib <- dplyr::right_join(cor_tib,lab_tib, 
                                 by=c("Ord_Key","U"="Address","AlleleA_ProbeSeq"="Ord_Prb"))
    
    if (!isMU)
      ses_tib <- ses_tib %>%
      dplyr::mutate(
        CHR=dplyr::case_when(
          Bsp_Tag_U=="UM" ~ CHR,
          TRUE ~ "0"
        ),
        MAPINFO=dplyr::case_when(
          Bsp_Tag_U=="UM" ~ MAPINFO,
          TRUE ~ as.integer(0)
        ),
        Strand_FR=dplyr::case_when(
          Bsp_Tag_U=="UM" ~ Strand_FR,
          TRUE ~ "0"
        )
      )
    
    ret_tib <- ses_tib %>%
      dplyr::select(Probe_ID,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,Rep_Num,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag_U,Bsp_Tag_M) %>% 
      dplyr::distinct() %>%
      dplyr::distinct(Probe_ID,Name,U,AlleleA_ProbeSeq, .keep_all = TRUE) %>%
      dplyr::distinct(U, .keep_all = TRUE) %>%
      clean_tibble()
    
    if (retData) ret_dat$all_tib <- all_tib
    if (retData) ret_dat$ses_tib <- ret_tib
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Assign Best CGN from:: BSP & SEQ
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

assign_cgn = function(ord, bsp, seq, can, csv=NULL,
                      merge=TRUE, retData=FALSE, join="inner",
                      verbose=0,vt=3,tc=1,tt=NULL,
                      funcTag='assign_cgn') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr} can={can}.{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  stime <- base::system.time({
    
    # Load Canonical CGNs::
    can_tib <- safe_read(file=can, verbose=verbose, vt=vt+1,tc=tc+1,tt=tt) %>% 
      dplyr::select(CGN) %>% 
      dplyr::rename(Cgn=CGN) %>% 
      dplyr::mutate(Can_Cnt=1)
    
    if (retData) ret_dat$can_tib <- can_tib
    
    # Defined Order tib to ord original cgn::
    ord_tib <- ord %>% 
      dplyr::select(Aln_Key,Ord_Cgn) %>%
      dplyr::rename(Cgn=Ord_Cgn) %>% 
      dplyr::mutate(Ord_Cnt=1) %>%
      dplyr::distinct()
    
    if (retData) ret_dat$ord_tib <- ord_tib
    
    # Format BSP::
    bsp_tib <- bsp %>% 
      dplyr::filter(!is.na(Bsp_Cgn)) %>% 
      dplyr::select(Ord_Key, Aln_Key, Ord_Des, Ord_Din, Bsp_Cgn) %>% 
      dplyr::rename(Cgn=Bsp_Cgn) %>%
      dplyr::distinct() %>%
      dplyr::arrange(Aln_Key, Cgn) %>% 
      dplyr::group_by(Ord_Key,Aln_Key,Ord_Des,Ord_Din,Cgn) %>% 
      dplyr::summarise(Bsp_Cnt=n(), .groups = "drop")
    bsp_key <- glue::glue("bsp-tib({funcTag})")
    bsp_cnt <- print_tib(bsp_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    
    if (retData) ret_dat$bsp_tib <- bsp_tib
    
    seq_tib <- seq %>% 
      dplyr::filter(!is.na(Imp_Cgn)) %>% 
      dplyr::select(Address, Ord_Des, Ord_Din, Imp_Cgn) %>% 
      tidyr::unite(Tmp_Key, Ord_Des,Ord_Din, sep="", remove=FALSE) %>%
      tidyr::unite(Aln_Key, Address, Tmp_Key, sep="_", remove=FALSE) %>%
      dplyr::select(-Tmp_Key) %>%
      dplyr::left_join(dplyr::select(ord, Ord_Key,Aln_Key), by="Aln_Key") %>%
      dplyr::select(Ord_Key,Aln_Key,Ord_Des,Ord_Din,Imp_Cgn) %>%
      dplyr::rename(Cgn=Imp_Cgn) %>%
      dplyr::distinct() %>%
      dplyr::arrange(Aln_Key, Cgn) %>% 
      dplyr::group_by(Ord_Key,Aln_Key,Ord_Des,Ord_Din,Cgn) %>% 
      dplyr::summarise(Seq_Cnt=n(), .groups = "drop")
    seq_key <- glue::glue("seq-tib({funcTag})")
    seq_cnt <- print_tib(seq_tib,funcTag, verbose,vt+4,tc, n=seq_key)
    
    if (retData) ret_dat$seq_tib <- seq_tib
    
    # Build and Sort Counts Tables
    cnt_tib <- 
      dplyr::full_join(bsp_tib, seq_tib, 
                       by=c("Ord_Key","Aln_Key","Ord_Des","Ord_Din","Cgn")) %>% 
      dplyr::left_join(can_tib, by="Cgn") %>%
      dplyr::distinct() %>%
      dplyr::left_join(ord_tib, by=c("Aln_Key","Cgn")) %>%
      dplyr::distinct() %>%
      dplyr::mutate(dplyr::across(c(Bsp_Cnt,Seq_Cnt,Can_Cnt,Ord_Cnt), tidyr::replace_na, 0 ),
                    Sum_Cnt=Bsp_Cnt+Seq_Cnt,
                    Max_Cnt=Bsp_Cnt*Seq_Cnt) %>% 
      dplyr::add_count(Aln_Key, name="Cgn_Cnt") %>% 
      dplyr::arrange(-Can_Cnt,-Max_Cnt,-Sum_Cnt,-Ord_Cnt) %>%
      dplyr::mutate(Rank=dplyr::row_number())
    
    if (retData) ret_dat$cnt_tib <- cnt_tib
    
    cnt_list <- cnt_tib %>% split(.$Ord_Des)
    
    # Infinium II::
    inf2_tib <- cnt_list[["2"]] %>%
      dplyr::arrange(Aln_Key, Rank) %>%
      dplyr::distinct(Aln_Key, .keep_all = TRUE)
    
    # Infinium I:: Full Join
    #
    #   TBD:: The joining should really be done by sequence: Ord_Prb
    #
    if (join=="full") {
      inf1_tib <- dplyr::full_join(
        cnt_list[["U"]], cnt_list[["M"]], 
        by=c("Ord_Key","Cgn","Ord_Din"), 
        suffix=c("_U","_M")
      ) %>%
        dplyr::mutate(Rank_Min=pmin(Rank_U,Rank_M)) %>%
        dplyr::arrange(Ord_Key, Rank_Min) %>%
        dplyr::distinct(Ord_Key,Aln_Key_U,Aln_Key_M, .keep_all = TRUE)
    } else if (join=="inner") {
      inf1_tib <- dplyr::inner_join(
        cnt_list[["U"]], cnt_list[["M"]], 
        by=c("Ord_Key","Cgn","Ord_Din"), 
        suffix=c("_U","_M")
      ) %>%
        dplyr::mutate(Rank_Min=pmin(Rank_U,Rank_M)) %>%
        dplyr::arrange(Ord_Key, Rank_Min) %>%
        dplyr::distinct(Ord_Key,Aln_Key_U,Aln_Key_M, .keep_all = TRUE)
    } else {
      stop(glue::glue("[{funcTag}]:{tabsStr} Unsupported join type={join}.{RET}"))
      return(NULL)
    }
    if (retData) ret_dat$inf1_tib <- inf1_tib
    if (retData) ret_dat$inf2_tib <- inf2_tib
    
    ret_tib <- dplyr::bind_rows(
      dplyr::select(inf1_tib, Ord_Key,Aln_Key_U,Cgn,Ord_Des_U,Ord_Din,Can_Cnt_U,Rank_Min) %>% 
        purrr::set_names("Ord_Key","Aln_Key","Cgn","Ord_Des","Ord_Din","Can_Cnt","Rank"),
      
      dplyr::select(inf1_tib, Ord_Key,Aln_Key_M,Cgn,Ord_Des_M,Ord_Din,Can_Cnt_M,Rank_Min) %>% 
        purrr::set_names("Ord_Key","Aln_Key","Cgn","Ord_Des","Ord_Din","Can_Cnt","Rank"),
      
      dplyr::select(inf2_tib, Ord_Key,Aln_Key,Cgn,Ord_Des,Ord_Din,Can_Cnt,Rank)
    ) %>% dplyr::filter(!is.na(Aln_Key)) %>%
      dplyr::distinct()
    
    mul_cnt <- ret_tib %>% dplyr::add_count(Aln_Key,Cgn, name="Multi_Cnt") %>% 
      dplyr::filter(Multi_Cnt != 1) %>% base::nrow()
    mis_cnt <- ret_tib %>% dplyr::filter(is.na(Aln_Key)) %>% base::nrow()
    
    mis_tib <- dplyr::anti_join(ord, ret_tib, by=c("Aln_Key"))
    sig_tib <- dplyr::filter(cnt_tib, Aln_Key %in% mis_tib$Aln_Key) %>%
      dplyr::arrange(Ord_Key,Rank) %>%
      dplyr::distinct(Aln_Key, .keep_all = TRUE)
    sig_cnt <- sig_tib %>% base::nrow()
    
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr}   Miss Count={mis_cnt}.{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}  Multi Count={mul_cnt}.{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr} Single Count={sig_cnt}.{RET}"))
      cat("\n")
    }
    # if (mis_cnt!=0 || mul_cnt!=0 || sig_cnt!=0) {
    #   stop(glue::glue("{RET}[{funcTag}]:{tabsStr} Counts non-zero={mis_cnt},",
    #                   "{mul_cnt},{sig_cnt}!{RET}{RET}"))
    #   return(NULL)
    # }
    
    # Merge all data together::
    #
    ret_cnt <- ret_tib %>% base::nrow()
    ret_tib <- dplyr::bind_rows(
      ret_tib %>% dplyr::mutate(
        Cgn_Tag=dplyr::case_when(
          Ord_Din=="rs" ~ Ord_Din,
          Ord_Din=="ch" ~ Ord_Din,
          TRUE ~ "cg"
        ),
        Cgn_Str=dplyr::case_when(
          Ord_Din=="rs" ~ stringr::str_remove(Ord_Key, "[-_:].*$"),
          Ord_Din=="ch" ~ stringr::str_remove(Ord_Key, "[-_:].*$"),
          TRUE ~ paste0("cg",stringr::str_pad(Cgn,width=8,side="left",pad="0"))
        )),
      mis_tib %>% 
        dplyr::select(Ord_Key, Aln_Key,Ord_Cgn,Ord_Des,Ord_Din) %>% 
        dplyr::rename(Cgn=Ord_Cgn) %>% 
        dplyr::mutate(Can_Cnt=0, 
                      Rank=dplyr::row_number() + ret_cnt,
                      Cgn_Tag="uk",
                      Cgn_Str=paste0(Cgn_Tag,stringr::str_pad(Cgn,width=8,side="left",pad="0"))
        )
    ) %>%
      # TBD:: Capture other CGN's in seperate column:: actual CGN's not Count!!
      dplyr::add_count(Aln_Key, name="Alt_Cgn_Cnt") %>%
      # One Final Clean Up To Ensure Uniqueness::
      dplyr::arrange(Rank) %>% 
      dplyr::distinct(Aln_Key, .keep_all = TRUE)
    
    mul_cnt <- ret_tib %>% 
      dplyr::add_count(Aln_Key,Cgn, name="Multi_Cnt") %>% 
      dplyr::filter(Multi_Cnt != 1) %>% base::nrow()
    
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr}  Multi Count Final={mul_cnt}.{RET}"))
      cat("\n")
    }
    if (mul_cnt!=0) {
      stop(glue::glue("{RET}[{funcTag}]:{tabsStr} Multi-Count Final={mul_cnt} ",
                      "not equal to zero!!!{RET}{RET}"))
      return(NULL)
    }
    
    if (merge) ret_tib <- bsp %>%
      dplyr::left_join(ret_tib, 
                       by=c("Ord_Key","Aln_Key","Ord_Des","Ord_Din"),
                       suffix=c("_bsp","_cgn"))
    
    ret_tib <- ret_tib %>% clean_tibble()
    out_cnt <- safe_write(ret_tib,file=csv, funcTag=funcTag,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    if (retData) ret_dat$ret_tib <- ret_tib
    if (retData) ret_dat$mis_tib <- mis_tib
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                ORD/MAT/AQP/Manifest File Workflows:: Generation
#                         AQP Address Workflow::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

aqp_address_workflow = function(ord, mat, aqp, out=NULL, name=NULL, csv=NULL,
                                prb_key="Ord_Prb", add_key="Address",
                                des_key="Ord_Des", din_key="Ord_Din",
                                verbose=0,vt=3,tc=1,tt=NULL,
                                funcTag='aqp_address_workflow') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                         Process AQP/PQC Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (purrr::is_character(aqp)) {
      aqp_tib <- load_aqp_files(aqp) %>% 
        dplyr::arrange(-Ord_Idx) %>% 
        dplyr::distinct(Address, .keep_all=TRUE)
    } else {
      aqp_tib <- aqp
    }
    if (!tibble::is_tibble(aqp_tib)) {
      stop(glue::glue("{RET}[{funcTag}]: AQP/PQC is NOT a tibble!{RET}{RET}"))
      return(ret_tib)
    }
    aqp_dup_cnt <- aqp_tib %>% 
      dplyr::add_count(Address, name="Add_Cnt") %>% 
      dplyr::filter(Add_Cnt!=1) %>% base::nrow()
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: aqp_dup_cnt={aqp_dup_cnt}.{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Process Match Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (purrr::is_character(mat)) {
      mat_tib <- load_aqp_files(mat) %>% 
        dplyr::arrange(-Ord_Idx) %>% 
        dplyr::distinct(Address, .keep_all=TRUE)
    } else {
      mat_tib <- mat
    }
    if (!tibble::is_tibble(mat_tib)) {
      stop(glue::glue("{RET}[{funcTag}]: Match is NOT a tibble!{RET}{RET}"))
      return(ret_tib)
    }
    
    mat_dup_cnt <- mat_tib %>% 
      dplyr::add_count(Address, name="Add_Cnt") %>% 
      dplyr::filter(Add_Cnt!=1) %>% base::nrow()
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: mat_dup_cnt={mat_dup_cnt}.{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Process Order Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (purrr::is_character(ord)) {
      ord_tib <- load_aqp_files(ord) %>% 
        dplyr::mutate(Ord_Map=Ord_Idx+Ord_Map)
    } else {
      ord_tib <- ord
    }
    if (!tibble::is_tibble(ord_tib)) {
      stop(glue::glue("{RET}[{funcTag}]: Order is NOT a tibble!{RET}{RET}"))
      return(ret_tib)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Build Probes/Address Table::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- aqp_tib %>% 
      dplyr::filter(Decode_Status==0) %>%
      dplyr::select(Address,Ord_Idx) %>%
      dplyr::rename(Aqp_Idx=Ord_Idx) %>%
      dplyr::inner_join(mat_tib, by="Address") %>%
      dplyr::rename(Mat_Idx=Ord_Idx) %>%
      dplyr::inner_join(ord_tib, by=c("Mat_Prb"="Ord_Prb")) %>%
      dplyr::rename(Ord_Prb=Mat_Prb) %>%
      dplyr::mutate(Ord_Cgn=Ord_Key %>% 
                      stringr::str_remove("^[^0-9]+") %>% 
                      stringr::str_remove("[^0-9]+.*$") %>% 
                      stringr::str_remove("^0+") ) %>%
      dplyr::select(Address,Ord_Des,Ord_Din,Ord_Map,Ord_Prb,Ord_Par,
                    Ord_Key,Ord_Col,Ord_Idx,Mat_Idx,Aqp_Idx,Mat_Tan,
                    dplyr::everything())
    
    prb_dup_cnt <- ret_tib %>% 
      dplyr::add_count(Address, name="Add_Cnt") %>% 
      dplyr::filter(Add_Cnt!=1) %>% base::nrow()
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: prb_dup_cnt={prb_dup_cnt}.{RET}"))
    
    par_dup_cnt <- ret_tib %>% 
      dplyr::add_count(Address,Ord_Prb,Ord_Par, name="Dup_Cnt") %>% 
      dplyr::filter(Dup_Cnt!=1) %>% base::nrow()
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: par_dup_cnt={par_dup_cnt}.{RET}"))
    
    # Build Summary::
    aqp_sum <- aqp_tib %>% 
      dplyr::group_by(Ord_Idx, Decode_Status) %>% 
      dplyr::summarise(Decode_Count=n(), .groups="drop")
    print(aqp_sum, n=base::nrow(aqp_sum))
    
    sum_key <- glue::glue("aqp-summary({funcTag})")
    sum_cnt <- print_tib(aqp_sum,funcTag, verbose,vt+4,tc, n=sum_key)
    
    # Write Fasta Output::
    #
    if (!is.null(out)) {
      if (length(name)>0) name <- paste0(name,".")
      dat_csv <- file.path(out, paste0(name,"aqp-pass.address.csv.gz"))
      prb_fas <- file.path(out, paste0(name,"aqp-pass.address.fas.gz"))
      u49_tsv <- file.path(out, paste0(name,"aqp-pass.address-u49.tsv.gz"))
      m49_tsv <- file.path(out, paste0(name,"aqp-pass.address-m49.tsv.gz"))
      sum_csv <- file.path(out, paste0(name,"aqp-decode.summary.csv.gz"))
      
      sum_cnt <- safe_write(x=aqp_sum, file=sum_csv,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      if (!is.null(csv)) dat_csv <- csv
      ret_tib <- add_to_fas(ret_tib,
                            prb_key=prb_key, add_key=add_key,
                            des_key=des_key, din_key=din_key,  
                            prb_fas=prb_fas, dat_csv=dat_csv,
                            u49_tsv=u49_tsv, m49_tsv=m49_tsv,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    } else if (!is.null(csv)) {
      out_cnt <- safe_write(x=ret_tib, file=csv, 
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

load_aqp_files = function(file,
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='load_aqp_files') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr} file={file}.{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_sum <- NULL
  stime <- base::system.time({
    
    if (purrr::is_character(file)) {
      file_vec <- stringr::str_split(file, pattern=",", simplify=TRUE) %>% 
        BiocGenerics::as.vector()
    } else if (purrr::is_vector(file)) {
      file_vec <- file
    } else {
      stop(glue::glue("{RET}[{funcTag}]:{tabsStr} unrecognized variable={file}.{RET}"))
      return(NULL)
    }
    
    ret_tib <- file_vec %>%
      lapply(load_aqp_file,
             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
      dplyr::bind_rows(.id="Ord_Idx") %>%
      clean_tibble()
    
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

load_aqp_file = function(file, idx=NULL,
                         n_max=100, guess=1000,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='load_aqp_file') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  # Saftey Check for Empty Files::
  if (is.null(file) || length(file)==0)
    return(base::invisible(NULL))
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr}  file={file}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} n_max={n_max}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} guess={guess}.{RET}"))
    cat("\n")
  }
  
  val_cols <- list()
  val_cols$ord <- 
    cols(
      Assay_Design_Id        = col_character(),
      AlleleA_Probe_Id       = col_character(),
      AlleleA_Probe_Sequence = col_character(),
      AlleleB_Probe_Id       = col_character(),
      AlleleB_Probe_Sequence = col_character(),
      Normalization_Bin      = col_character()
    )
  
  val_cols$mat1 <- 
    cols(
      Plate    = col_character(),
      Row      = col_character(),
      Col      = col_integer(),
      Address  = col_integer(),
      Mod5     = col_character(),
      Sequence = col_character(),
      Mod3     = col_character(),
      Comments = col_character()
    )
  
  val_cols$mat2 <- 
    cols(
      address_names = col_integer(),
      probe_id      = col_character(),
      sequence      = col_character(),
      type_b        = col_character(),
      address_name  = col_integer(),
      bo_seq        = col_character()
    )
  
  val_cols$mat3 <- 
    cols(
      probe_id      = col_character(),
      sequence      = col_character(),
      type_b        = col_character(),
      address_names = col_integer(),
      bo_seq        = col_character()
    )
  
  val_cols$aqp <- 
    cols(
      Address           = col_integer(),
      Decode_Status     = col_integer(),
      Decode_Error_Code = col_integer(),
      Decode_Score      = col_integer(),
      Func_Status       = col_integer(),
      Func_Error_Code   = col_integer(),
      QC_Action         = col_integer()
    )
  
  val_cols$pqc <- 
    cols(
      Address      = col_integer(),
      Status       = col_integer(),
      Eval_Code    = col_integer(),
      Average_Rep  = col_integer(),
      Expected_Rep = col_integer()
    )
  
  sel_cols <- list()
  sel_cols$ord  <- c("Assay_Design_Id","AlleleA_Probe_Sequence",
                     "AlleleB_Probe_Sequence","Normalization_Bin")
  sel_cols$mat1 <- c("Address","Sequence")
  sel_cols$mat2 <- c("address_names","bo_seq")
  sel_cols$mat3 <- c("address_names","bo_seq")
  sel_cols$aqp  <- c("Address","Decode_Status")
  sel_cols$pqc  <- c("Address","Status")
  
  key_cols <- list()
  key_cols$ord  <- c("Ord_Key","Ord_PrbA","Ord_PrbB","Ord_Norm")
  key_cols$mat1 <- c("Address","Sequence")
  key_cols$mat2 <- c("Address","Sequence")
  key_cols$mat3 <- c("Address","Sequence")
  key_cols$aqp  <- c("Address","Decode_Status")
  key_cols$pqc  <- c("Address","Decode_Status")
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    guess_tib <- guess_aqp_file(file, n_max=n_max,
                                verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    row_num=guess_tib$row_num[1]
    col_num=guess_tib$col_num[1]
    del_key=guess_tib$del_key[1]
    beg_key=guess_tib$beg_key[1]
    
    dat_key <- NULL
    val_col <- NULL
    sel_col <- NULL
    key_col <- NULL
    if (is.null(beg_key) || is.null(col_num)) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Either beg_key OR col_num is ",
                      "NULL!!!{RET}{RET}"))
      return(ret_tib)
    } else if (beg_key==names(val_cols$ord$cols)[1] && 
               col_num==length(val_cols$ord$cols)) {
      dat_key <- "ord"
    } else if (beg_key==names(val_cols$mat1$cols)[1] && 
               col_num==length(val_cols$mat1$cols)) {
      dat_key <- "mat1"
    } else if (beg_key==names(val_cols$mat2$cols)[1] && 
               col_num==length(val_cols$mat2$cols)) {
      dat_key <- "mat2"
    } else if (beg_key==names(val_cols$mat3$cols)[1] && 
               col_num==length(val_cols$mat3$cols)) {
      dat_key <- "mat3"
    } else if (beg_key==names(val_cols$aqp$cols)[1] && 
               col_num==length(val_cols$aqp$cols)) {
      dat_key <- "aqp"
    } else if (beg_key==names(val_cols$pqc$cols)[1] && 
               col_num==length(val_cols$pqc$cols)) {
      dat_key <- "pqc"
    } else {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Failed to match beg_key({beg_key}) AND ",
                      "col_num({col_num}) to known formats!!!{RET}{RET}"))
      return(ret_tib)
    }
    
    # This sets all the proper valid col types, col selection and col renaming::
    val_col <- val_cols[[dat_key]]  # File Format Column Names/Data-Types
    sel_col <- sel_cols[[dat_key]]  # Selected Column Names
    key_col <- key_cols[[dat_key]]  # Selected Column New Names
    
    # Apply Delimiter Type and Column Names/Data-Types::
    if (del_key==COM) {
      ret_tib <- suppressMessages(suppressWarnings(
        readr::read_csv(file, skip=row_num, guess_max=guess, skip_empty_rows=FALSE,
                        col_names=names(val_col$cols), col_types=val_col)))
      
    } else if (del_key==TAB || del_key==" ") {
      ret_tib <- suppressMessages(suppressWarnings(
        readr::read_tsv(file, skip=row_num, guess_max=guess, skip_empty_rows=FALSE,
                        col_names=names(val_col$cols), col_types=val_col)))
      
    } else {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Unsupported del_key={del_key}!!!{RET}{RET}"))
      return(NULL)
    }
    ret1_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n="ret1")
    
    # Apply Selected Columns::
    ret_tib  <- ret_tib %>% dplyr::select(dplyr::all_of(sel_col))
    ret2_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n="ret2")
    
    # Apply Canonical-Renamed Columns::
    ret_tib  <- ret_tib %>% purrr::set_names(key_col)
    ret3_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n="ret3")
    
    # Apply Formatting Rules:: Order Data Type
    if (dat_key=="ord") {
      ret_cnt <- ret_tib %>% base::nrow()
      ord_tib <- ret_tib %>% 
        dplyr::distinct() %>%
        dplyr::mutate(
          Ord_PrbA=stringr::str_to_upper(Ord_PrbA),
          Ord_PrbB=stringr::str_to_upper(Ord_PrbB),
          
          Ord_DesA=dplyr::case_when(
            !is.na(Ord_PrbA) & !is.na(Ord_PrbB) ~ "U",
            !is.na(Ord_PrbA) &  is.na(Ord_PrbB) ~ "2",
            TRUE ~ NA_character_
          ),
          Ord_DesB=dplyr::case_when(
            !is.na(Ord_PrbA) & !is.na(Ord_PrbB) ~ "M",
            TRUE ~ NA_character_
          ),
          Ord_Col=dplyr::case_when(
            Ord_DesA=="2" ~ NA_character_,
            Ord_DesA=="U" & Ord_DesB=="M" & Ord_Norm=="A" ~ "R",
            Ord_DesA=="U" & Ord_DesB=="M" & Ord_Norm=="B" ~ "G",
            TRUE ~ NA_character_
          )
        ) %>% 
        dplyr::distinct(Ord_PrbA,Ord_PrbB, .keep_all=TRUE) %>%
        dplyr::mutate(Ord_Row=dplyr::row_number(),
                      Ord_Map=as.double(Ord_Row/ret_cnt))
      
      ord_cols <- c("Ord_Key","Ord_Des", "Ord_Prb", "Ord_Par", "Ord_Col","Ord_Map")
      ord_colA <- c("Ord_Key","Ord_DesA","Ord_PrbA","Ord_PrbB","Ord_Col","Ord_Map")
      ord_colB <- c("Ord_Key","Ord_DesB","Ord_PrbB","Ord_PrbA","Ord_Col","Ord_Map")
      
      ret_tib <- dplyr::bind_rows(
        ord_tib %>% dplyr::filter(Ord_DesB=="M") %>%
          dplyr::select(dplyr::all_of(ord_colB)) %>%
          purrr::set_names(ord_cols),
        ord_tib %>% dplyr::filter(Ord_DesA=="U") %>%
          dplyr::select(dplyr::all_of(ord_colA)) %>%
          purrr::set_names(ord_cols),
        ord_tib %>% dplyr::filter(Ord_DesA=="2") %>%
          dplyr::select(dplyr::all_of(ord_colA)) %>%
          purrr::set_names(ord_cols)
      ) %>% 
        dplyr::mutate(
          Ord_Din=stringr::str_sub(Ord_Key, 1,2),
          Ord_Key=stringr::str_replace_all(Ord_Key, "_", "-")
        ) %>%
        dplyr::select(
          Ord_Map,Ord_Key,Ord_Des,Ord_Din,Ord_Col,dplyr::everything()
        )
    }
    
    # Apply Formatting Rules:: Match/AQP/PQC Data Types:: Address
    # if (dat_key=="mat1" || dat_key=="mat2" ||
    #     dat_key=="aqp"  || dat_key=="pqc") {
    if (dat_key!="ord") {
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
    
    # Apply Formatting Rules:: Match Files:: Sequence
    if (stringr::str_starts(dat_key,"mat")) {
      ret_tib <- ret_tib %>% 
        dplyr::mutate(Sequence=stringr::str_to_upper(Sequence),
                      Mat_Tan=stringr::str_sub(Sequence,1,45),
                      Mat_Prb=stringr::str_sub(Sequence,46)
        ) %>%
        dplyr::select(-Sequence) %>%
        dplyr::filter(!is.na(Address)) %>% 
        dplyr::distinct(Address,Mat_Prb, .keep_all=TRUE)
    }
    
    # Add Dat_IDX if provided
    if (!is.null(idx))
      ret_tib <- ret_tib %>%
      dplyr::mutate(Dat_IDX=as.integer(idx))
    
    # Standard Column Data Type Clean-Up::
    ret_tib <- ret_tib %>% clean_tibble()
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

guess_aqp_file = function(file, 
                          fields=c("Assay_Design_Id","Plate",
                                   "address_names","Address", "probe_id"),
                          cols=NULL,
                          n_max=100, # guess=100000,
                          verbose=0,vt=6,tc=1,tt=NULL,
                          funcTag='guess_aqp_file') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr}   file={file}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  n_max={n_max}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} fields={RET}"))
    print(fields)
    cat("\n")
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_val <- NULL
  stime <- base::system.time({
    
    file_del_str <- guess_file_del(file, n_max=n_max,
                                   verbose=verbose,vt=vt+4)
    if (is.null(file_del_str)) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: file_del_str=NULL!!!{RET}{RET}"))
      return(ret_tib)
    }
    
    dat_tib <- readr::read_lines(file, n_max=n_max) %>%
      tibble::as_tibble() %>% 
      dplyr::mutate(row_num=dplyr::row_number())
    dat_cnt <- print_tib(dat_tib,funcTag, verbose,vt=vt+6,tc=tc, n="dat_tib")
    
    ret_tib <- NULL
    for (field in fields) {
      if (verbose>=vt+6)
        cat(glue::glue("[{funcTag}]:{tabsStr} field={field}{RET}"))
      
      min_tib <- dat_tib %>% 
        dplyr::filter(stringr::str_starts(value, field)) %>% tail(n=1)
      min_cnt <- base::nrow(min_tib)
      
      # Skip Empty Results::
      if (min_cnt==0) next
      
      # Extract Column Count Step::
      col_num <- min_tib$value[1] %>% 
        stringr::str_split(file_del_str, simplify=TRUE) %>% 
        as.vector() %>% length()
      
      # Combine Step::
      ret_tib <- ret_tib %>% dplyr::bind_rows(
        min_tib %>%
          dplyr::mutate(col_num=col_num,del_key=file_del_str,beg_key=field)
      )
      if (verbose>=vt+6)
        cat(glue::glue("[{funcTag}]:{tabsStr} NEXT...{RET}{RET}"))
    }
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         CGN Mapping Workflow Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cgn_mapping_workflow = function(ref_u49,can_u49,out_u49,
                                ref_m49,can_m49,out_m49,
                                ord=NULL,bed=NULL,org=NULL,out=NULL,
                                idxA=1,idxB=1,reload=FALSE,
                                verbose=0,vt=3,tc=1,tt=NULL,
                                funcTag='cgn_mapping_workflow') {
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (verbose>=vt) {
    cat(glue::glue("[{funcTag}]:{tabsStr} ref_u49={ref_u49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} can_u49={can_u49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} out_u49={out_u49}.{RET}"))
    cat("\n")
    cat(glue::glue("[{funcTag}]:{tabsStr} ref_m49={ref_m49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} can_m49={can_m49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} out_m49={out_m49}.{RET}"))
    cat("\n")
    cat(glue::glue("[{funcTag}]:{tabsStr}    idxA={idxA}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}    idxB={idxB}.{RET}"))
    cat("\n")
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    build_file_dir(out_u49)
    build_file_dir(out_m49)
    
    u49_tib <- intersect_seq(ref=ref_u49,can=can_u49,out=out_u49,
                             idxA=idxA,idxB=idxB, reload=reload,
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    u49_key <- glue::glue("u49-tib({funcTag})")
    u49_cnt <- print_tib(u49_tib,funcTag, verbose,vt+4,tc, n=u49_key)
    
    m49_tib <- intersect_seq(ref=ref_m49,can=can_m49,out=out_m49,
                             idxA=idxA,idxB=idxB, reload=reload,
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    m49_key <- glue::glue("m49-tib({funcTag})")
    m49_cnt <- print_tib(m49_tib,funcTag, verbose,vt+4,tc, n=m49_key)
    
    ret_tib <- join_seq_intersect(u49=u49_tib, m49=m49_tib, 
                                  bed=bed, org=org,
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    ret_key <- glue::glue("ret-tib-0({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    if (!is.null(ord) && file.exists(ord)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Loading order/canonical CSV={ord}.{RET}"))
      
      ord_tib <- suppressMessages(suppressWarnings( readr::read_csv(ord) )) %>%
        purrr::set_names(c("Can_Cgn","Can_Top","Can_Src")) %>%
        dplyr::mutate(Can_Cgn=as.integer(Can_Cgn), 
                      Can_Scr=as.integer(1),
                      Can_Scr=tidyr::replace_na(Can_Scr, 0)) %>%
        clean_tibble()
      ord_cnt <- print_tib(ord_tib,funcTag, verbose,vt+4,tc, n="ord_tib")
      
      ret_tib <- ret_tib %>% 
        dplyr::left_join(ord_tib, by=c("Imp_Cgn"="Can_Cgn")) %>%
        dplyr::mutate(Can_Scr=tidyr::replace_na(Can_Scr, 0)) %>%
        clean_tibble()
      
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Done order/canonical.{RET}{RET}"))
    }
    
    if (!is.null(out)) safe_write(ret_tib,file=out, funcTag=funcTag,
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

intersect_seq = function(ref, can, out, idxA=1, idxB=1, reload=FALSE,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='intersect_seq') {
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr}  ref={ref}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  can={can}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  out={out}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} idxA={idxA}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} idxB={idxB}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} reload={reload}{RET}{RET}"))
  }
  
  int_seq_cols <-
    cols(
      Imp_Seq  = col_character(),
      Imp_Nuc  = col_character(),
      
      Imp_SrdI = col_integer(),
      Imp_Srd3 = col_character(),
      
      Imp_Key  = col_character(),
      Imp_Scr  = col_character(),
      
      Imp_Cnt  = col_integer(),
      Aln_Key  = col_character()
    )
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    if (reload && file.exists(out)) {
      if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Reloading!{RET}"))
    } else {
      clean <- FALSE
      if (stringr::str_ends(can, '.gz')) {
        cmd_str <- glue::glue("gzip -f -k -d {can}")
        if (verbose>=vt)
          cat(glue::glue("[{funcTag}]: Running cmd={cmd_str}...{RET}"))
        
        cmd_ret <- base::system(cmd_str)
        if (cmd_ret!=0) {
          stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed(cmd_ret={cmd_ret}) ",
                          "cmd={cmd_str}!{RET}{RET}"))
          return(ret_tib)
        }
        can <- stringr::str_remove(can, ".gz$")
        clean <- TRUE
      }
      
      cmd_str = glue::glue("gzip -dc {ref} | join -t $'\t' -1{idxA} -2{idxB} - {can} | gzip -c - > {out}")
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]: Running cmd={cmd_str}...{RET}"))
      cmd_ret <- system(cmd_str)
      
      if (clean && !stringr::str_ends(can, '.gz')) {
        cmd_str <- glue::glue("rm {can}")
        if (verbose>=vt)
          cat(glue::glue("[{funcTag}]: Running cmd={cmd_str}...{RET}"))
        
        cmd_ret <- system(cmd_str)
        if (cmd_ret!=0) {
          stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed(cmd_ret={cmd_ret}) ",
                          "cmd={cmd_str}!{RET}{RET}"))
          return(ret_tib)
        }
        can <- paste(can,'gz', sep='.')
      }
    }
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: Loading intersection output={out}...{RET}"))
    
    ret_tib <- suppressMessages(suppressWarnings(
      readr::read_tsv(out, col_names=names(int_seq_cols$cols), 
                      col_types=int_seq_cols) )) # %>% clean_tibble()
    
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

join_seq_intersect = function(u49,m49,bed=NULL,org=NULL,
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'join_seq_intersect'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    imp_col_vec <- c("Address","Ord_Des","Ord_Din",
                     "Imp_Chr","Imp_Pos","Imp_Cgn",
                     "Imp_FR","Imp_TB","Imp_CO","Imp_Nxb",
                     "Bsp_Din_Ref","Bsp_Din_Scr","Aln_Prb")
    
    ret_tib <- 
      dplyr::bind_rows(u49,m49) %>%
      dplyr::select(-Imp_SrdI,-Imp_Scr) %>% 
      dplyr::mutate(Imp_Key=stringr::str_split(Imp_Key, pattern=",") ) %>% 
      tidyr::unnest(Imp_Key)
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret_tib1")
    
    ret_tib <- ret_tib %>%
      tidyr::separate(
        Imp_Key, 
        into=c("Imp_Cgn","Imp_Hit_hg38","Imp_Hit_hg37", "Imp_Hit_hg36", "Imp_Hit_mm10"), 
        sep="_", remove=TRUE)
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret_tib2")
    
    ret_tib <- ret_tib %>%
      tidyr::separate(Aln_Key, into=c("Address", "Tmp_Key"), 
                      sep="_", remove=FALSE) %>%
      tidyr::separate(Tmp_Key, into=c("Ord_Des","Ord_Din"), sep=1) %>%
      dplyr::rename(Aln_Prb=Imp_Seq, Aln_Nuc=Imp_Nuc)
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret_tib3")
    
    ret_tib <- ret_tib %>%
      tidyr::separate(Imp_Srd3, into=c("Imp_TB","Imp_CO", "Imp_Nxb"),
                      sep=c(1,2)) %>%
      dplyr::select(Address, Ord_Des, Ord_Din, Imp_Cgn, 
                    Imp_TB, Imp_CO, Imp_Nxb, Aln_Prb, Aln_Nuc, 
                    dplyr::everything()) %>%
      clean_tibble()
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret_tib4")
    
    if (!is.null(bed)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Adding cgn bed...{RET}"))
      
      ret_tib <- bed %>%
        dplyr::right_join(ret_tib,by=c("Imp_Cgn")) %>%
        dplyr::mutate(
          Imp_FR=dplyr::case_when(
            Imp_Top_Srd=="+" & Imp_TB=="T" ~ "F",
            Imp_Top_Srd=="-" & Imp_TB=="T" ~ "R",
            Imp_Top_Srd=="+" & Imp_TB=="B" ~ "R",
            Imp_Top_Srd=="-" & Imp_TB=="B" ~ "F",
            TRUE ~ NA_character_
          )
        )
    }
    
    if (!is.null(org)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Adding org bed...{RET}"))
      
      ret_tib <- ret_tib %>%
        dplyr::left_join(
          org, 
          by=c("Aln_Prb"="Org_49P",
               "Ord_Des"="Org_Des",
               "Imp_Chr"="Org_Chr",
               "Imp_Pos"="Org_Pos",
               "Imp_FR"="Org_FR",
               "Imp_TB"="Org_TB",
               "Imp_CO"="Org_CO")
        )
    }
    
    ret_tib <- ret_tib  %>%
      dplyr::select(dplyr::any_of(imp_col_vec),dplyr::everything()) %>%
      clean_tibble()
    
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Manifest Mutation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: Pretty Sure this can be removed::
mutate_probe_id = function(tib, 
                           pid="Probe_ID", cgn="Imp_Cgn_Seq",
                           des="Ord_Des",  din="Ord_Din",
                           tb="Imp_TB_Seq", co="Imp_CO_Seq",
                           inf="Infinium_Design",pad=8,
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutate_probe_id'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    pid_sym <- rlang::sym(pid)
    cgn_sym <- rlang::sym(cgn)
    des_sym <- rlang::sym(des)
    din_sym <- rlang::sym(din)
    tb_sym  <- rlang::sym(tb)
    co_sym  <- rlang::sym(co)
    inf_sym <- rlang::sym(inf)
    
    ret_tib <- tib %>% 
      dplyr::mutate(
        !!inf_sym:=dplyr::case_when(
          !!des_sym=="U" | !!des_sym=="M" ~ 1, 
          !!des_sym=="2" ~ 2, 
          TRUE ~ NA_real_) %>% as.integer(), 
        !!pid_sym:=paste(paste0(!!din_sym,stringr::str_pad(!!cgn_sym,width=pad, side="left", pad="0")), 
                         paste0(!!tb_sym,!!co_sym,!!inf_sym), sep="_")
      )
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Address To Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: Pretty Sure this can be removed::
add_comb = function(tibA, tibB, field,
                    join,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'add_comb'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ret_tib <- 
      dplyr::inner_join(
        tibA, tibB, 
        by=dplyr::all_of(join),
        suffix=c("_U","_M")
      ) # %>%
    # dplyr::select(dplyr::starts_with(field)) %>%
    # dplyr::distinct() %>%
    # purrr::set_names(c("U","M"))
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

add_to_man = function(tib, join, runName,
                      des_key="Ord_Des", pid_key="Ord_Key",
                      rep_key=NULL, rep_val=NULL,
                      col_key=NULL, nxb_key=NULL,
                      csv=NULL, validate=TRUE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'add_to_man'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    des_list <- NULL
    des_list <- tib %>% split(.[[des_key]])
    des_cnt  <- des_list %>% names() %>% length()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} des_list[{des_key}]={RET}"))
      print(des_list)
    }
    
    # Build Infinium I::
    ret1_tib <- NULL
    if (!is.null(des_list[["U"]]) && !is.null(des_list[["M"]])) {
      ret1_tib <- dplyr::full_join(
        dplyr::select(des_list[["U"]], -dplyr::all_of(des_key)), 
        dplyr::select(des_list[["M"]], -dplyr::all_of(des_key)), 
        by=dplyr::all_of(join),
        suffix=c("_U","_M")
      ) %>% dplyr::mutate(Infinium_Design=as.integer(1))
      # TBD:: We should allow these "singletons" to pass, but under
      #  a different classification...
      #
      if ("Address_U" %in% names(ret1_tib) && 
          "Address_M" %in% names(ret1_tib)) ret1_tib <- ret1_tib %>% 
          dplyr::filter(!is.na(Address_U) & !is.na(Address_M))
    }
    ret1_cnt <- print_tib(ret1_tib,funcTag, verbose,vt+4,tc, n="InfI")
    
    # Build Infinium II::
    ret2_tib <- NULL
    if (!is.null(des_list[["2"]])) {
      ret2_tib <- dplyr::bind_cols(
        dplyr::select(des_list[["2"]],  dplyr::all_of(join)),
        dplyr::select(des_list[["2"]], -dplyr::all_of(join)) %>% 
          purrr::set_names(paste(names(.),"U", sep="_"))
      ) %>% dplyr::mutate(Infinium_Design=as.integer(2))
      ret2_cnt <- print_tib(ret2_tib,funcTag, verbose,vt+4,tc, n="InfII")
    }
    
    # Bind Infinium I/II into single manifest::
    ret_tib <- dplyr::bind_rows(ret1_tib, ret2_tib) %>%
      dplyr::mutate(
        # Infinium_Design=dplyr::case_when(
        #   is.na(Address_M) ~ 2,
        #   !is.na(Address_M) ~ 1,
        #   TRUE ~ NA_real_
        # ) %>% as.integer(),
        Infinium_Design_Type=dplyr::case_when(
          Infinium_Design==1 ~ 'I',
          Infinium_Design==2 ~ 'II',
          TRUE ~ NA_character_)
      )
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="InfI+II.0")
    
    if (!is.null(col_key)) {
      col_sym <- rlang::sym(col_key)
      ret_tib <- ret_tib %>%
        dplyr::mutate(
          col=dplyr::case_when(
            Infinium_Design==1 ~ !!col_sym,
            TRUE ~ NA_character_),
          Color_Channel=dplyr::case_when(
            Infinium_Design==2 ~ 'Both',
            col=="R" ~ 'Red',
            col=='G' ~ 'Grn',
            TRUE ~ NA_character_),
          Next_Base=dplyr::case_when(
            col=="R" ~ "A",
            col=="G" ~ "C",
            TRUE ~ NA_character_
          ),
          Probe_Source=runName
        )
    } else if (!is.null(nxb_key)) {
      nxb_sym <- rlang::sym(nxb_key)
      
      ret_tib <- ret_tib %>%
        dplyr::mutate(
          Color_Channel=dplyr::case_when(
            Infinium_Design==2 ~ 'Both',
            Infinium_Design==1 & (!!nxb_sym=='A' | !!nxb_sym=='T' |
                                    !!nxb_sym=='a' | !!nxb_sym=='t') ~ 'Red',
            Infinium_Design==1 & (!!nxb_sym=='C' | !!nxb_sym=='G' |
                                    !!nxb_sym=='c' | !!nxb_sym=='g') ~ 'Grn',
            TRUE ~ NA_character_),
          col=dplyr::case_when(
            Infinium_Design==2 ~ NA_character_,
            Infinium_Design==1 ~ stringr::str_sub(Color_Channel, 1,1),
            TRUE ~ NA_character_),
          Next_Base=dplyr::case_when(
            col=="R" ~ !!nxb_sym,
            col=="G" ~ !!nxb_sym,
            TRUE ~ NA_character_
          ),
          # Next_Base=!!nxb_sym,
          Probe_Source=runName
        )
    } else {
      # Do nothing for now, but this should not be a standard use case
      return(NULL)
    }
    
    if (!is.null(pid_key)) {
      pid_key_sym = rlang::sym(pid_key)
      ret_tib <- ret_tib %>% dplyr::arrange(pid_key)
    }
    
    if (!is.null(pid_key) && !is.null(rep_key) && !is.null(rep_val)) {
      rep_key_sym = rlang::sym(rep_key)
      rep_val_sym = rlang::sym(rep_val)
      
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Adding Probe Replicate({rep_key}/{rep_val})...{RET}"))
      ret_tib <- ret_tib %>% dplyr::group_by(!!rep_key_sym) %>%
        dplyr::mutate(!!rep_val_sym := dplyr::row_number(),
                      !!pid_key_sym := paste0(!!rep_key_sym,!!rep_val_sym)) %>%
        dplyr::ungroup()
    }
    
    if (!is.null(csv)) {
      safe_write(ret_tib,"csv",csv,
                 funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="InfI+II")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             FASTA File Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_to_fas = function(tib, prb_key="Prb_Seq", 
                      add_key="Address", des_key="Prb_Des", din_key="prb_type",
                      prb_fas=NULL,dat_csv=NULL,
                      u49_tsv=NULL,m49_tsv=NULL,
                      del="_",
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'add_to_fas'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; key={prb_key}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    prb_sym <- rlang::sym(prb_key)
    des_sym <- rlang::sym(des_key)
    din_sym <- rlang::sym(din_key)
    aln_vec <- c(add_key,des_key)
    
    # aln_vec <- c(des_key,din_key)
    #
    # To make all Aln_Key's <=15 characters for improbe we'll 
    #   unite the Aln_Key in two steps to remove one underscore ("_")
    #
    # Old Definition: aln_vec <- c(add_key,des_key,din_key)
    #
    
    # Build Alignment Keys/Seqs
    ret_tib <- tib %>% 
      tidyr::unite(Tmp_Key, dplyr::all_of(aln_vec),sep=del,remove=FALSE) %>%
      tidyr::unite(Aln_Key, Tmp_Key,!!din_sym,sep="",remove=FALSE) %>%
      dplyr::select(-Tmp_Key) %>%
      dplyr::mutate(
        Aln_Prb=deMs(!!prb_sym, uc=TRUE),
        Aln_Rev=revCmp(Aln_Prb),
        Aln_P49=dplyr::case_when(
          !!des_sym == '2' ~ stringr::str_sub(Aln_Prb, 2),
          !!des_sym == 'U' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          !!des_sym == 'M' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::distinct(Aln_Key,Aln_Prb, .keep_all=TRUE) %>%
      clean_tibble()
    
    # Generate Remaining Data::
    if (!is.null(prb_fas))
      fas_vec <- ret_tib %>%
      dplyr::mutate(fas_line=paste0(">",Aln_Key,"\n",Aln_Prb) ) %>%
      dplyr::pull(fas_line)
    
    u49_tib <- NULL
    if (!is.null(u49_tsv))
      u49_tib <- ret_tib %>% 
      dplyr::filter(!!des_sym == '2' | !!des_sym=='U') %>%
      dplyr::filter(!is.na(Aln_P49)) %>%
      dplyr::select(Aln_P49,Aln_Key) %>%
      dplyr::arrange(Aln_P49)
    
    m49_tib <- NULL
    if (!is.null(m49_tsv))
      m49_tib <- ret_tib %>% 
      dplyr::filter(!!des_sym=='M') %>%
      dplyr::filter(!is.na(Aln_P49)) %>%
      dplyr::select(Aln_P49,Aln_Key) %>%
      dplyr::arrange(Aln_P49)
    
    # Safe Write Outputs::
    safe_write(ret_tib,"csv",dat_csv,  
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    safe_write(fas_vec,"line",prb_fas, 
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    safe_write(u49_tib,"tsv",u49_tsv,cols=FALSE, 
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    safe_write(m49_tib,"tsv",m49_tsv,cols=FALSE, 
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Genome FASTA Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_genome = function(file, chr_key="Chromosome", nrec=0, ret_map=TRUE,
                       verbose=0,vt=3,tc=1,tt=NULL,
                       funcTag='load_genome') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("[{funcTag}]:{tabsStr} Function Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   file={file}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   nrec={nrec}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  stime <- base::system.time({
    
    if (nrec==0) {
      seqs <- Biostrings::readDNAStringSet(filepath=file, format="fasta")
    } else {
      seqs <-  Biostrings::readDNAStringSet(filepath=file, format="fasta", 
                                            nrec=nrec)
    }
    ret_cnt <- seqs %>% length()
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Loaded {ret_cnt} chromosomes.{RET}"))
    if (verbose>=vt+2) print(seqs)
    
    if (ret_map) {
      chr_sym <- rlang::sym(chr_key)
      
      ret_tib <- seqs %>% 
        names() %>% 
        stringr::str_remove(" .*$") %>% 
        stringr::str_remove("^chr") %>%
        tibble::tibble() %>% 
        purrr::set_names(chr_key) %>% 
        dplyr::mutate(Idx=dplyr::row_number(),
                      !!chr_sym := paste0("chr",!!chr_sym) )

      ret_dat$seqs <- seqs
      ret_dat$maps <- ret_tib
      
      ret_key <- glue::glue("ret-FIN({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  if (ret_map) return(ret_dat)
  
  ret_tib
}

parse_genomic_seqs = function(tib, seq,
                              srd_str="F", 
                              pos_key="Coordinate",
                              chr_str="chr0", 
                              ups_len=60,
                              seq_len=122, 
                              verbose=0,vt=3,tc=1,tt=NULL,
                              funcTag='parse_genomic_seqs') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Field Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   srd_str={srd_str}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   pos_key={pos_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   chr_str={chr_str}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   ups_len={ups_len}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   seq_len={seq_len}{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    pos_sym  <- rlang::sym(pos_key)
    
    cur_key <- glue::glue("cur_tib_1=chr_str={chr_str}")
    cur_cnt <- print_tib(tib,funcTag, verbose,vt+4,tc, n=cur_key)
    
    # Format Coordinates::
    #   Subtract two for cg upstream 122mer formation and add two+two as well
    begs <- NULL
    ends <- NULL
    begs <- tib %>% dplyr::pull(!!pos_sym) - ups_len - 2
    ends <- begs + seq_len - 1 + 4
    
    if (verbose>=vt+6) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} begs={RET}"))
      begs %>% head() %>% print()
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ends={RET}"))
      ends %>% head() %>% print()
    }
    
    # Sub-string ref sequences::
    seq_vec <- NULL
    seq_vec <- stringr::str_sub( as.character(seq), begs, ends ) %>%
      stringr::str_to_upper()
    ret_cnt <- seq_vec %>% length()
    
    # *** NOTE:: DO NOT use this if loading BSC genomes ***
    #
    #         [ they're already strand specific ]
    #
    if (srd_str=="R") seq_vec <- seq_vec %>% cmpl()
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Parsed Sequences({ret_cnt}):{RET}"))
      seq_vec %>% head(n=2) %>% print()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  seq_vec
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Common Conversion Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fas_to_seq = function(tib, 
                      
                      nrec=0, 
                      gen_bld="na", 
                      gen_ref_fas, 
                      bsc_ref_fas=NULL,
                      gen_snp_fas=NULL, 
                      bsc_snp_fas=NULL,
                      
                      imp_tsv=NULL,
                      seq_csv=NULL,
                      
                      build=c(""),
                      
                      ids_key="Aln_Key_Unq",
                      din_key="Ord_Din",
                      
                      ext_seq="Ext_Forward_Seq",
                      iup_seq="Iupac_Forward_Sequence",
                      imp_seq="Forward_Sequence",
                      
                      srd_str="F",
                      pos_key="Coordinate",
                      chr_key="Chromosome",

                      srsplit=FALSE,
                      srd_key=NULL,
                      cosplit=FALSE,
                      cos_key=NULL,
                      
                      ref_col="Ref",
                      alt_col="Alt",
                      iup_col="Iupac",

                      ups_len=60, 
                      seq_len=122, 
                      iupac=NULL,
                      del="_",
                      
                      reload=FALSE,
                      retData=FALSE,
                      
                      parallel=FALSE,
                      r_improbe=FALSE,
                      s_improbe=FALSE,
                      add_flanks=FALSE, 
                      add_matseq=TRUE, 
                      
                      verbose=0,vt=4,tc=1,tt=NULL,
                      funcTag='fas_to_seq') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("[{funcTag}]:{tabsStr} Genome Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}          nrec={nrec}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       gen_bld={gen_bld}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   gen_ref_fas={gen_ref_fas}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   bsc_ref_fas={bsc_ref_fas}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   gen_snp_fas={gen_snp_fas}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   bsc_snp_fas={bsc_snp_fas}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Output File Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       imp_tsv={imp_tsv}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       seq_csv={seq_csv}{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Field Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       ids_key={ids_key}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       din_key={din_key}{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       ext_seq={ext_seq}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       iup_seq={iup_seq}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       imp_seq={imp_seq}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       srd_str={srd_str}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       pos_key={pos_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       chr_key={chr_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       srsplit={srsplit}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       srd_key={srd_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       cosplit={cosplit}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       cos_key={cos_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       ref_col={ref_col}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       alt_col={alt_col}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       iup_col={iup_col}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}           del={del}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}         iupac={iupac}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       ups_len={ups_len}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       seq_len={seq_len}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}      reload={reload}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     retData={retData}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}    parallel={parallel}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   r_improbe={r_improbe}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   s_improbe={s_improbe}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  add_flanks={add_flanks}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  add_matseq={add_matseq}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  etime   <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # Reload if all data is present::
  #
  if (reload &&
      !purrr::is_null(imp_tsv) && file.exists(imp_tsv) &&
      !purrr::is_null(seq_csv) && file.exists(seq_csv)) {
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Reloading seq_csv{seq_csv}...{RET}"))
    
    stime <- base::system.time({
      ret_tib <- safe_read(seq_csv)
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    return(ret_tib)
  }
  
  stime <- base::system.time({
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Building fresh...{RET}"))
    
    # Define symbolic variables::
    #
    ids_sym  <- rlang::sym(ids_key)
    ext_sym  <- rlang::sym(ext_seq)
    iup_sym  <- rlang::sym(iup_seq)
    imp_sym  <- rlang::sym(imp_seq)
    
    pos_sym  <- rlang::sym(pos_key)
    chr_sym  <- rlang::sym(chr_key)

    ref_col_sym  <- rlang::sym(ref_col)
    alt_col_sym  <- rlang::sym(alt_col)
    iup_col_sym  <- rlang::sym(iup_col)
    
    # Load Genome::
    #
    dna_dat <- load_genome(file=gen_ref_fas, nrec=nrec,
                           chr_key=chr_key, ret_map=TRUE,
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # Split Data by Chromosome::
    #
    chr_list <- tib %>% 
      dplyr::arrange(!!chr_sym,!!pos_sym) %>%
      split(.[[chr_key]])
    
    chr_maps <- dna_dat$maps %>%
      split(.[[chr_key]])
    
    # Process each chromosome::
    #  TBD:: Add parallel computing over chromosomes
    #
    chr_vec_1 <- names(chr_list)
    chr_vec_2 <- names(chr_maps)
    chr_names <- chr_vec_1[chr_vec_1 %in% chr_vec_2]
    
    if (parallel) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Extacting sequence templates ",
                       "from genome (Parallel)...{RET}"))
      
      ret_tib <- foreach (chr_str=chr_names, .combine=rbind) %dopar% {
        chr_idx <- chr_maps[[chr_str]] %>% head(n=1) %>% pull(Idx) %>% as.integer()
        cur_tib <- s_improbe_template_workflow(
          tib=chr_list[[chr_str]], seq=dna_dat$seqs[[chr_idx]],
          srd_str=srd_str, pos_key=pos_key, chr_key=chr_key, chr_str=chr_str,
          ext_seq=ext_seq, iup_seq=iup_seq, imp_seq=imp_seq,
          ups_len=ups_len, seq_len=seq_len, iupac=iupac, del=del,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      }

    } else {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Extacting sequence templates ",
                       "from genome (Linear)...{RET}"))

      for (chr_str in chr_names) {
        chr_idx <- chr_maps[[chr_str]] %>% head(n=1) %>% pull(Idx) %>% as.integer()
        cur_tib <- s_improbe_template_workflow(
          tib=chr_list[[chr_str]], seq=dna_dat$seqs[[chr_idx]], 
          srd_str=srd_str, pos_key=pos_key, chr_key=chr_key, chr_str=chr_str,
          ext_seq=ext_seq, iup_seq=iup_seq, imp_seq=imp_seq,
          ups_len=ups_len, seq_len=seq_len, iupac=iupac, del=del,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                          Merge Probe Designs::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
        cur_key <- glue::glue("cur-tib({funcTag}-{chr_str})")
        cur_cnt <- print_tib(cur_tib %>% dplyr::select(!!ids_sym,!!imp_sym),
                             funcTag, verbose,vt=vt+4,tc=tc+1, n=cur_key)
        
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Substring ",
                         "chr_str={chr_str}.{RET}{RET}"))
      }
    }
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Extacting sequence templates ",
                   "from genome.{RET}{RET}"))
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose+100,vt=0,tc, n=ret_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #      Forward Template Sequence Generation:: Tri-fecta s-improbe
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (add_flanks) {
      
      tri_tib <- NULL
      tri_tib <- s_improbe_trifecta(tib=ret_tib,
                                    
                                    tar_din=tar_din, 
                                    ids_key=ids_key, 
                                    din_key=din_key, 
                                    
                                    pos_key=pos_key, 
                                    chr_str=chr_str, 
                                    
                                    ext_seq=ext_seq, 
                                    iup_seq=iup_seq, 
                                    imp_seq=imp_seq, 
                                    
                                    verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
      
      ret_tib <- dplyr::bind_rows(ret_tib, tri_tib)
    }
    ret_tib <- ret_tib %>% dplyr::arrange(!!ids_sym)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Probe Design:: s_improbe()
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (s_improbe && length(build)>0) {
      
      ret_tib <- s_improbe(tib=ret_tib, 
                           build=build, 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Probe Design:: r_improbe()
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (r_improbe) {
      
      retData   <- TRUE
      r_imp_tib <- r_improbe(tib=ret_tib,
                             
                             sr_str='FR', 
                             co_str='CO',
                             
                             ids_key=ids_key,
                             seq_key=iup_seq,
                             prb_key=din_key,
                             
                             srsplit=srsplit,
                             srd_key=srd_key,
                             cosplit=cosplit,
                             cos_key=cos_key,
                             
                             prb_len=ups_len,
                             seq_len=seq_len,
                             
                             parallel=parallel, 
                             add_matseq=add_matseq,
                             
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      ret_dat$r_imp <- r_imp_tib
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Write improbe Design Input File:: TSV
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (!is.null(imp_tsv)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing improbe input TSV={imp_tsv}...{RET}"))
      
      imp_col <- c("Seq_ID","Sequence","Genome_Build",
                   "Chromosome","Coordinate","CpG_Island")
      imp_tib <- ret_tib %>% 
        dplyr::mutate(Genome_Build=!!gen_bld, 
                      CpG_Island="FALSE") %>%
        dplyr::select(dplyr::all_of(
          c(!!ids_key,  !!imp_seq, "Genome_Build", 
            !!chr_key, !!pos_key, "CpG_Island")) ) %>%
        purrr::set_names(imp_col)
      
      out_cnt <- safe_write(x=imp_tib, file=imp_tsv, funcTag=funcTag, 
                            verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Calculate Data Summary:: CSV
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Calculating Summary...{RET}"))
    
    sum_tib <- ret_tib %>% 
      dplyr::group_by(up61,dn61,Des_Din) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    sum_key <- glue::glue("ret-summary({funcTag})")
    sum_cnt <- print_tib(sum_tib,funcTag, verbose,vt+4,tc, n=sum_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Write Full Data Set:: CSV
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (!is.null(seq_csv)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing seq table CSV={seq_csv}...{RET}"))
      
      sum_csv <- seq_csv %>% 
        stringr::str_remove(".gz$") %>%
        stringr::str_remove(".[tcsv]+$") %>%
        paste(".summary.csv.gz")
      
      seq_cnt <- safe_write(x=ret_tib, file=seq_csv, funcTag=funcTag, 
                            verbose=verbose,vt=vt,tc=tc,append=FALSE)
      sum_cnt <- safe_write(x=sum_tib, file=sum_csv, funcTag=funcTag, 
                            verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    
    if (retData) ret_dat$s_imp <- ret_tib
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
}

# End of file
