
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='template_func') {
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # ret_cnt <- ret_tib %>% base::nrow()
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
    "{tabsStr}# ----- ----- ----- ----- |----- ----- ----- ----- #{RET}{RET}"))

  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Assign Best CGN from:: BSP & SEQ
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

assign_cgn = function(add, bsp, seq, can, csv=NULL,
                      merge=TRUE, retData=FALSE,
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
    
    # Defined Order tib to add original cgn::
    ord_tib <- add %>% 
      dplyr::select(Aln_Key,Ord_Cgn) %>%
      dplyr::rename(Cgn=Ord_Cgn) %>% 
      dplyr::mutate(Ord_Cnt=1) %>%
      dplyr::distinct()
    
    # Format BSP::
    bsp_tib <- bsp %>% 
      dplyr::filter(!is.na(Bsp_Cgn)) %>% 
      dplyr::select(Ord_Key,Aln_Key, Ord_Des, Ord_Din, Bsp_Cgn) %>% 
      dplyr::rename(Cgn=Bsp_Cgn) %>%
      dplyr::distinct() %>%
      dplyr::arrange(Aln_Key, Cgn) %>% 
      dplyr::group_by(Ord_Key,Aln_Key,Ord_Des,Ord_Din,Cgn) %>% 
      dplyr::summarise(Bsp_Cnt=n(), .groups = "drop")
    bsp_key <- glue::glue("bsp-tib({funcTag})")
    bsp_cnt <- print_tib(bsp_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    
    seq_tib <- seq %>% 
      dplyr::filter(!is.na(Imp_Cgn)) %>% 
      dplyr::select(Address, Ord_Des, Ord_Din, Imp_Cgn) %>% 
      tidyr::unite(Aln_Key, Address, Ord_Des, Ord_Din, sep="_", remove=FALSE) %>%
      dplyr::left_join(dplyr::select(aqp_add_tib, Ord_Key,Aln_Key), by="Aln_Key") %>%
      dplyr::select(Ord_Key,Aln_Key,Ord_Des,Ord_Din,Imp_Cgn) %>%
      dplyr::rename(Cgn=Imp_Cgn) %>%
      dplyr::distinct() %>%
      dplyr::arrange(Aln_Key, Cgn) %>% 
      dplyr::group_by(Ord_Key,Aln_Key,Ord_Des,Ord_Din,Cgn) %>% 
      dplyr::summarise(Seq_Cnt=n(), .groups = "drop")
    seq_key <- glue::glue("seq-tib({funcTag})")
    seq_cnt <- print_tib(seq_tib,funcTag, verbose,vt+4,tc, n=seq_key)
    
    # Build and Sort Counts Tables
    cnt_tib <- 
      dplyr::full_join(bsp_tib, seq_tib, by=c("Ord_Key","Aln_Key","Ord_Des","Ord_Din","Cgn")) %>% 
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
    
    cnt_list <- cnt_tib %>% split(.$Ord_Des)
    
    # Infinium II::
    inf2_tib <- cnt_list[["2"]] %>%
      dplyr::arrange(Aln_Key, Rank) %>%
      dplyr::distinct(Aln_Key, .keep_all = TRUE)
    
    # Infinium I:: Full Join
    #
    #   TBD:: The joining should really be done by sequence: Ord_Prb
    #
    inf1_tib <- dplyr::full_join(
      cnt_list[["U"]], cnt_list[["M"]], 
      by=c("Ord_Key","Cgn","Ord_Din"), 
      suffix=c("_U","_M")
    ) %>%
      dplyr::mutate(Rank_Min=pmin(Rank_U,Rank_M)) %>%
      dplyr::arrange(Ord_Key, Rank_Min) %>%
      dplyr::distinct(Ord_Key,Aln_Key_U,Aln_Key_M, .keep_all = TRUE)
    
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
    
    print(ret_tib)
    
    mul_cnt <- ret_tib %>% dplyr::add_count(Aln_Key,Cgn, name="Multi_Cnt") %>% 
      dplyr::filter(Multi_Cnt != 1) %>% base::nrow()
    mis_cnt <- ret_tib %>% dplyr::filter(is.na(Aln_Key)) %>% base::nrow()

    mis_tib <- dplyr::anti_join(aqp_add_tib, ret_tib, by=c("Aln_Key"))
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
    
    safe_write(ret_tib,file=csv, funcTag=funcTag,
               verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    if (retData) ret_dat$ret_tib <- ret_tib
    if (retData) ret_dat$mis_tib <- mis_tib
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  if (retData) return(ret_dat)

  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                ORD/MAT/AQP/Manifest File Workflows:: Generation
#                         AQP Address Workflow::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

aqp_address_workflow = function(ord, 
                                mat=NULL, aqp=NULL,
                                bpn=NULL, aqn=NULL,
                                add_csv=NULL, man_csv=NULL,
                                add_fas=NULL,
                                u49_tsv=NULL, m49_tsv=NULL,
                                runName=NA, retData=FALSE,
                                verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'aqp_address_workflow'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr} ord={ord}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} mat={mat}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} aqp={aqp}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} bpn={bpn}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} aqn={aqn}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} runName={runName}{RET}{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  stime <- base::system.time({
    
    ord_tib <- NULL
    mat_tib <- NULL
    aqp_tib <- NULL
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 Validate and Load Data:: Order Files
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    ord_vec <- stringr::str_split(ord, pattern=",", simplify=TRUE) %>% 
      BiocGenerics::as.vector()
    ord_len <- length(ord_vec)
    bpn_vec <- c(1:ord_len)
    aqn_vec <- c(1:ord_len)
    ord_tib <- load_aqp_files(ord_vec, verbose=verbose, vt=vt+1,tc=tc+1, tt=tt)
    if (!is.null(bpn)) bpn_vec <- stringr::str_split(bpn, pattern=",", simplify=TRUE) %>% 
      BiocGenerics::as.vector()
    
    if (retData) ret_dat$ord <- ord_tib

    if (!is.null(mat)) {
      if (is.null(aqp)) {
        stop(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: ",
                        "Must provide aqp if provided match files!",
                        "aqp = NULL! Exiting...{RET}"))
        return(ret_tib)
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                 Validate and Load Data:: Match Files
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      mat_vec <- stringr::str_split(mat, pattern=",", simplify=TRUE) %>% 
        BiocGenerics::as.vector()
      mat_len <- length(mat_vec)
      # Suspending this rule for now...
      #
      # if (ord_len != mat_len) {
      #   stop(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: ",
      #                   "order length != match length! ",
      #                   "({ord_len} != {mat_len})! Exiting...{RET}"))
      #   return(ret_tib)
      # }
      mat_tib <- load_aqp_files(mat_vec, verbose=verbose, vt=vt+1,tc=tc+1, tt=tt)
      
      if (retData) ret_dat$mat <- mat_tib
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                 Validate and Load Data:: AQP/PQC Files
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      aqp_vec <- stringr::str_split(aqp, pattern=",", simplify=TRUE) %>% 
        BiocGenerics::as.vector()
      aqp_len <- length(aqp_vec)
      aqn_vec <- c(1:aqp_len)
      if (!is.null(aqn)) aqn_vec <- 
        stringr::str_split(aqn, pattern=",", simplify=TRUE) %>% 
        BiocGenerics::as.vector()
      
      if (aqp_len>1 && aqp_len != mat_len) {
        stop(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: ",
                        "AQP>1 must match match length! ",
                        "({aqp_len} != {mat_len})! Exiting...{RET}"))
        return(ret_tib)
      }
      aqp_tib <- load_aqp_files(aqp_vec, verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
      
      if (retData) ret_dat$aqp <- aqp_tib
    }
    if (verbose>=vt+4)
      cat(glue::glue("{RET}[{funcTag}]:{tabsStr} Done. Loading data.{RET}{RET}"))
    ord_key <- glue::glue("ret-ord({funcTag})")
    ord_cnt <- print_tib(ord_tib,funcTag, verbose,vt+4,tc, n=ord_key)

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Build Bead Pool/AQP/PQC Info Data Structure::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    exp_tib <- 
      tibble::tibble(Bpn_Idx=bpn_vec, Aqp_Idx=aqn_vec) %>%
      dplyr::mutate(Ord_Idx=dplyr::row_number()) %>% 
      dplyr::select(Ord_Idx, dplyr::everything()) %>% 
      clean_tibble()
    exp_key <- glue::glue("ret-exp({funcTag})")
    exp_cnt <- print_tib(exp_tib,funcTag, verbose,vt+4,tc, n=exp_key)
    if (retData) ret_dat$exp <- exp_tib
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                                Join Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    ret_tib <- ord_tib %>% 
      dplyr::arrange(-Ord_Idx,Ord_Key) %>%
      dplyr::left_join(exp_tib, by="Ord_Idx")
    ret_key <- glue::glue("ret-ord({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)

    if (!is.null(mat_tib) && !is.null(aqp_tib)) {
      if (verbose>=vt+4)
        cat(glue::glue("[{funcTag}]:{tabsStr} Joining mat/aqp...{RET}"))
      
      mat_tib <- mat_tib %>% dplyr::left_join(exp_tib, by="Ord_Idx")
      aqp_tib <- aqp_tib %>% dplyr::left_join(exp_tib, by="Ord_Idx")
      
      if (mat_len==aqp_len) {
        if (verbose>=vt+2)
          cat(glue::glue("[{funcTag}]:{tabsStr} AQP MATCHING...{RET}"))
        
        ret_tib <- ret_tib %>%
          dplyr::full_join(mat_tib, by=c("Ord_Prb"="Mat_Prb",
                                         "Ord_Idx","Bpn_Idx","Aqp_Idx")) %>%
          dplyr::left_join(aqp_tib, by=c("Address",
                                         "Ord_Idx","Bpn_Idx","Aqp_Idx") )
      } else {
        if (verbose>=vt+2)
          cat(glue::glue("[{funcTag}]:{tabsStr} PQC MATCHING...{RET}"))
        
        aqp_tib <- aqp_tib %>% 
          dplyr::select(-dplyr::all_of( c("Ord_Idx","Bpn_Idx","Aqp_Idx") ))
        
        ret_tib <- ret_tib %>%
          dplyr::full_join(mat_tib, by=c("Ord_Prb"="Mat_Prb",
                                         "Ord_Idx","Bpn_Idx","Aqp_Idx")) %>%
          dplyr::left_join(aqp_tib, by=c("Address") )
      }
      ret_key <- glue::glue("pre-pass-decode({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)

      ret_tib <- ret_tib %>%
        dplyr::filter(!is.na(Decode_Status)) %>%
        dplyr::filter(Decode_Status==0) %>%
        dplyr::distinct(Address, .keep_all=TRUE)
      
      ret_key <- glue::glue("post-pass-decode({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    }
    #
    # Add Artifical Address, Decode_Status
    #
    print(ret_tib)
    if (!"Address" %in% names(ret_tib)) ret_tib <- ret_tib %>% 
      dplyr::mutate(Address=dplyr::row_number())
    if (!"Decode_Status" %in% names(ret_tib)) ret_tib <- ret_tib %>% 
      dplyr::mutate(Decode_Status=as.integer(0))
    
    # Add Order Predicted CGN::
    ret_tib <- ret_tib %>%
      dplyr::mutate(Ord_Cgn=Ord_Key %>% 
                      stringr::str_remove("^[a-zA-Z]+") %>% 
                      stringr::str_remove("-.*$") %>% 
                      stringr::str_remove("_.*$") %>% 
                      as.integer())

    ret_tib <- ret_tib %>%
      dplyr::add_count(Ord_Prb, name="Ord_Prb_Rep") %>%
      dplyr::add_count(Ord_Prb,Ord_Par, name="Ord_Par_Rep") %>%
      dplyr::select(Ord_Key,Ord_Des,Ord_Din,Ord_Col,Ord_Idx,Ord_Prb,
                    Ord_Prb_Rep,Ord_Par_Rep,
                    dplyr::everything())
    ret_key <- glue::glue("ret-tib3({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)

    
    # Write Summary::
    if (verbose>=vt+4 &&
        !is.null(mat_tib) && !is.null(aqp_tib)) {
      ret_sum <- 
        ret_tib %>% 
        dplyr::group_by(Ord_Idx,Bpn_Idx,Aqp_Idx,
                        Ord_Des,Ord_Col,Ord_Din,
                        Decode_Status) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      ret_sum %>% print(n=base::nrow(ret_sum))
      # sum_cnt <- print_tib(ret_sum,funcTag, verbose,vt+4,tc, n="summary")
    }
    
    # Write Manifest Output::
    if (!is.null(man_csv)) {
      man_vec <- c("Ord_Key","Ord_Din","Ord_Col",
                   dplyr::all_of(base::names(exp_tib) ) )
      print(man_vec)
      print(ret_tib)
      
      man_tib <- ret_tib %>%
        add_to_man(join=man_vec,
                   runName=runName,
                   des_key="Ord_Des", pid="Ord_Key",
                   col_key="Ord_Col",
                   csv=man_csv,
                   validate=TRUE,
                   verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      if (retData) ret_dat$man <- man_tib
    }
    
    # Write Outputs::
    if (!is.null(add_fas)) {
      # Set Look Up Keys::
      if ("Address" %in% names(ret_tib)) {
        add_key <- "Address"
      } else if ("Ord_Key" %in% names(ret_tib)) {
        add_key <- "Ord_Key"
      } else if ("Ord_Prb_U" %in% names(ret_tib)) {
        add_key <- "Ord_Prb_U"
      } else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR Failed to find valid add_key!!!{RET}{RET}"))
        return(NULL)
      }
      
      print(ret_tib)
      # Write Fasta Output::
      ret_tib <- ret_tib %>%
        add_to_fas(
          prb_key="Ord_Prb", add_key=add_key, 
          des_key="Ord_Des", din_key="Ord_Din",
          prb_fas=add_fas, dat_csv=add_csv,
          u49_tsv=u49_tsv, m49_tsv=m49_tsv,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
    }
    if (!is.null(add_csv)) 
      safe_write(ret_tib,file=add_csv, funcTag=funcTag,
                 verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    if (retData) ret_dat$add <- ret_tib
    
    # Overall Summary::
    # if (verbose>=vt+4) {
    #   aqp_add_sum <- aqp_add_tib %>% 
    #     dplyr::group_by(dplyr::any_of(c("Aqp_Idx","Ord_Des","Ord_Din"))) %>%
    #     dplyr::summarise(Count=n(), .groups="drop")
    #   aqp_add_sum %>% print(n=base::nrow(aqp_add_sum))
    # }
    
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  if (retData) return(ret_dat)
  
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
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
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
      ord_tib <- ret_tib %>% 
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
        )
      
      ord_cols <- c("Ord_Key","Ord_Des", "Ord_Prb", "Ord_Par", "Ord_Col")
      ord_colA <- c("Ord_Key","Ord_DesA","Ord_PrbA","Ord_PrbB","Ord_Col")
      ord_colB <- c("Ord_Key","Ord_DesB","Ord_PrbB","Ord_PrbA","Ord_Col")
      
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
          Ord_Key,Ord_Des,Ord_Din,Ord_Col,dplyr::everything()
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
    # if (dat_key=="mat1" || dat_key=="mat2") {
    if (stringr::str_starts(dat_key,"mat")) {
      ret_tib <- ret_tib %>% 
        dplyr::mutate(Sequence=stringr::str_to_upper(Sequence),
                      Mat_Tan=stringr::str_sub(Sequence,1,45),
                      Mat_Prb=stringr::str_sub(Sequence,46)
        ) %>%
        dplyr::select(-Sequence) %>%
        dplyr::filter(!is.na(Address)) # Not sure if we need this...
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
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
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
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
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
    
    if (!is.null(out)) safe_write(ret_tib,file=out, funcTag=par$prgmTag,
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)

    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
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
      aln_key  = col_character()
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
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
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
      tidyr::separate(
        aln_key, 
        into=c("Address", "Ord_Des", "Ord_Din"), 
        sep="_") %>%
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
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Manifest Mutation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Address To Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
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
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
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
    aln_vec <- c(add_key,des_key,din_key)
    
    # Build Alignment Keys/Seqs
    ret_tib <- tib %>% # tail(n=200) %>%
      tidyr::unite(Aln_Key, dplyr::all_of(aln_vec),sep=del,remove=FALSE) %>%
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
    # utils::type.convert() %>% 
    # dplyr::mutate(across(where(is.factor), as.character) )
    
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
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Common Conversion Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

seq_to_prbs = function(tib, seq,
                       cgn="Seq_ID",chr="Chromosome", pos="Coordinate",
                       srd="F",
                       ups_len=60, seq_len=122,
                       iupac=NULL, add_flank=FALSE,
                       del="_", 
                       verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'seq_to_prbs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    cgn_sym <- rlang::sym(cgn)
    chr_sym <- rlang::sym(chr)
    pos_sym <- rlang::sym(pos)
    
    cur_begs <- NULL
    cur_ends <- NULL
    cur_begs <- tib %>% dplyr::pull(!!pos_sym) - ups_len - 2
    cur_ends <- cur_begs + seq_len - 1 + 4
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_begs={RET}"))
      cur_begs %>% head() %>% print()
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_ends={RET}"))
      cur_ends %>% head() %>% print()
    }
    
    # Sub-string ref sequences::
    #
    ref_seqs <- NULL
    ref_seqs <- 
      stringr::str_sub( as.character(dna_seqs[[chr_idx]]), cur_begs, cur_ends) %>%
      stringr::str_to_upper()
    
    if (srd=="R") ref_seqs <- ref_seqs %>% cmpl()
    
    # Build individual parts of the template sequence::
    # up01.up02...up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn02.dn01
    #                             [  C   G  ]
    #
    ref_up01 <- stringr::str_sub(ref_seqs, 1,1)
    ref_up02 <- stringr::str_sub(ref_seqs, 2,2)
    ref_up58 <- stringr::str_sub(ref_seqs, 3,ups_len-2+2)
    
    ref_up59 <- stringr::str_sub(ref_seqs, ups_len-2+3,ups_len-2+3)
    ref_up60 <- stringr::str_sub(ref_seqs, ups_len-2+4,ups_len-2+4)
    ref_up61 <- stringr::str_sub(ref_seqs, ups_len-2+5,ups_len-2+5)
    
    iupac_vec <- ref_up61
    if (!is.null(iupac)) iupac_vec <- cur_tib %>% dplyr::pull(!!iupac)
    
    ref_dn61 <- stringr::str_sub(ref_seqs, ups_len-2+6,ups_len-2+6)
    ref_dn60 <- stringr::str_sub(ref_seqs, ups_len-2+7,ups_len-2+7)
    ref_dn59 <- stringr::str_sub(ref_seqs, ups_len-2+8,ups_len-2+8)
    
    ref_dn58 <- stringr::str_sub(ref_seqs, ups_len-2+9,ups_len-2+ups_len-2+8)
    ref_dn02 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+8,ups_len-2+ups_len-2+8)
    ref_dn01 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+9,ups_len-2+ups_len-2+9)
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ref_seqs({chr_str})={RET}"))
      ref_seqs %>% head() %>% print()
    }
    
    # Add all fields to current tibble::
    #
    cur_tib <- cur_tib %>%
      dplyr::mutate(
        !!ext_sym := ref_seqs %>% addBrac(),
        
        up01=ref_up01,
        up02=ref_up02,
        up58=ref_up58,
        
        up59=ref_up59,
        up60=ref_up60,
        up61=ref_up61,
        
        iupac=iupac_vec,
        
        dn61=ref_dn61,
        dn60=ref_dn60,
        dn59=ref_dn59,
        
        dn58=ref_dn58,
        dn02=ref_dn02,
        dn01=ref_dn01,
        
        # Now we can assemble to optimal template::
        #
        !!des_sym := paste0(up58,
                            up59,up60,
                            iupac,dn61,
                            dn60,dn59,
                            dn58) %>%
          addBrac(),
        !!imp_sym := paste0(up58,
                            up59,up60,
                            "CG",
                            dn60,dn59,
                            dn58) %>%
          addBrac(),
        
        # TBD:: Pretty sure this right, but really needs more testing...
        #
        Din_Str=stringr::str_to_lower(paste0(iupac,dn61)),
        Des_Din=dplyr::case_when(
          up61=='C' & dn61=='G' ~ "cg",
          up61=='C' & dn61!='G' ~ "ch",
          up61!='C' ~ "rs",
          TRUE ~ NA_character_
        )
      )
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_tib({chr_str})={RET}"))
      cur_tib %>% print()
    }
    
    #  
    # TBD:: We can also break out tri-fecta sites if probe type == rs
    #
    ups_tib <- NULL
    dns_tib <- NULL
    
    if (add_flank) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Upstream CG Flank Seuqnces({chr_str})...{RET}"))
      
      ups_tib <- cur_tib %>% 
        dplyr::filter(up59=="C" & up60=="G") %>%
        dplyr::mutate(
          !!des_sym := paste0(up01,up02,up58,up59,up60,iupac,dn61,dn60,dn59,dn58) %>%
            stringr::str_sub(1,seq_len) %>%
            addBrac(),
          Din_Str=stringr::str_to_lower(paste0(up59,up60)),
          Des_Din="cg",
          Name=paste(!!name_sym,"CG_UP", sep=del),
          Coordinate=as.integer(Coordinate-2)
        )
      
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Downstream CG Flank Seuqnces({chr_str})...{RET}"))
      
      dns_tib <- cur_tib %>% 
        dplyr::filter(dn61=="C" & dn60=="G") %>%
        dplyr::mutate(
          !!des_sym := paste0(up58,up59,up60,iupac,dn61,dn60,dn59,dn58,dn02) %>%
            stringr::str_sub(2) %>%
            addBrac(),
          Din_Str=stringr::str_to_lower(paste0(dn61,dn60)),
          Des_Din="cg",
          Name=paste(!!name_sym,"CG_DN", sep=del),
          Coordinate=as.integer(Coordinate+2)
        )
    }
    
    # Add New Tri-fecta Probes::
    #  - NOTE:: These will have their names changed later after alignment...
    #
    if (verbose>=vt+1)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joining and sorting...{RET}"))
    
    cur_tib <- dplyr::bind_rows(cur_tib, ups_tib, dns_tib) %>%
      dplyr::arrange(!!name_sym)
    
    ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
    
    if (verbose>=vt+1)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Substring chr_str={chr_str}.{RET}{RET}"))
    
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

bed_to_prbs = function(tib, fas,
                       file=NULL,
                       din="Probe_Type",gen="na",
                       cgn="Seq_ID",chr="Chromosome", pos="Coordinate",
                       fwd_seq="Fwd_Temp_Seq",
                       iup_seq="Iup_Temp_Seq",
                       imp_seq="Imp_Temp_Seq",
                       ups_len=60, seq_len=122, nrec=0,
                       iupac=NULL, add_flank=FALSE,
                       del="_", 
                       verbose=0,vt=3,tc=1,tt=NULL,
                       funcTag='bed_to_prbs') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    cgn_sym <- rlang::sym(cgn)
    chr_sym <- rlang::sym(chr)
    pos_sym <- rlang::sym(pos)
    
    # Get List of BSC Genomes::
    # fas_dir <- fas %>% base::dirname()
    fas_pre <- fas %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$")
    if (is.null(gen)) gen <- base::basename(fas_pre)
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Prefix({gen})={fas_pre}.{RET}"))
    
    # NOTE:: There isn't an opposite genome, build that from the converted::
    #  This is why we don't need 'srd_cos <- c("C","O")' below::
    #
    srd_cos <- c("C")
    srd_frs <- c("F","R")
    srd_mus <- c("U","M","D")
    
    for (fr in srd_frs) {
      # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} fr={fr}...{RET}"))
      
      for (co in srd_cos) {
        # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} co={co}...{RET}"))
        
        for (mu in srd_mus) {
          if (verbose>=vt) 
            cat(glue::glue("[{funcTag}]:{tabsStr} fr={fr}, co={co}, mu={mu}...{RET}"))
          gen_fas <- paste0(fas_pre,".", fr,co,mu, ".fa.gz")
          
          if (!file.exists(gen_fas)) {
            if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Skipping; gen_fas=({fr},{co},{mu})",
                                            "={gen_fas} Does NOT Exist!{RET}"))
            next
          } else {
            if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading; gen_fas({fr},{co},{mu})",
                                            "={gen_fas}...{RET}"))
          }
          
          cur_tib <- fas_to_seq(
            tib=tib, fas=gen_fas,
            file=file,
            name=cgn,din=din,gen=gen,
            chr1=chr,pos=pos,
            srd=fr,
            fwd_seq=fwd_seq,iup_seq=iup_seq,imp_seq=imp_seq,
            
            # See Above:: Will Pass these variables in later...
            # fwd_seq="Fwd_Temp_Seq",des_seq="Iup_Temp_Seq",imp_seq="Imp_Temp_Seq",
            
            iupac=iupac,
            nrec=nrec,
            add_flank=FALSE,
            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt
          ) %>%
            dplyr::mutate(
              Srd_FR=fr,
              Srd_CO=co,
              Prb_Des=mu
            )
          cur_cnt <- print_tib(cur_tib,funcTag, verbose,vt+4,tc, n="cur-pre-join")
          cgn_vec <- cur_tib %>% dplyr::pull(!!cgn_sym)

          # TBD:: Add the current FR/CO/TB fields to cur_tib
          # TBD:: Add the opposite designs
          # TBD:: Remove previous information from the input tibble
          
          prb_tib <- NULL
          if (fr=="F") {
            prb_tib <- tibble::tibble(
              !!cgn_sym := cgn_vec,
              
              # Forward Converted Designs:: Infinium I/II
              Prb_SeqC1 = paste0(
                cur_tib$up61,
                cur_tib$dn61,
                cur_tib$dn60,
                cur_tib$dn59,
                cur_tib$dn58 %>% stringr::str_sub(1,46)
              ) %>% revCmp(),
              Prb_SeqC2 = paste0(
                # cur_tib$up61,
                cur_tib$dn61,
                cur_tib$dn60,
                cur_tib$dn59,
                cur_tib$dn58 %>% stringr::str_sub(1,47)
              ) %>% revCmp(),
              
              # Forward Opposite Designs:: Infinium I/II
              #   Just need to back up by one base::
              Prb_SeqO1 = paste0(
                cur_tib$up58 %>% stringr::str_sub(12,58), # 13 -> 12
                cur_tib$up59,
                cur_tib$up60,
                cur_tib$up61
                # cur_tib$dn61 # Remove dn61
              ), # No complement
              Prb_SeqO2 = paste0(
                cur_tib$up58 %>% stringr::str_sub(11,58), # 12 -> 11
                cur_tib$up59,
                cur_tib$up60
                # cur_tib$up61 # Remove dn61
                # cur_tib$dn61,
                # cur_tib$dn60
              ) # No Complement
              
            )
          } else if (fr=="R") {
            #
            # NOTE: From the ch/rs perspective the reverse strand 
            #   coordinates frame should be moved back by 1
            #
            # if ()
            prb_tib <- tibble::tibble(
              # Seq_ID = cur_tib$Seq_ID,
              # !!cgn_sym := cur_tib[[!!cgn_sym]],
              !!cgn_sym := cgn_vec,
              
              # Reverse Converted Designs:: Infinium I/II
              Prb_SeqC1 = paste0(
                cur_tib$up58 %>% stringr::str_sub(13,58),
                cur_tib$up59,
                cur_tib$up60,
                cur_tib$up61,
                cur_tib$dn61
              ) %>% cmpl(),
              Prb_SeqC2 = paste0(
                cur_tib$up58 %>% stringr::str_sub(12,58),
                cur_tib$up59,
                cur_tib$up60,
                cur_tib$up61
                # cur_tib$dn61,
                # cur_tib$dn60
              ) %>% cmpl(),
              
              # Reverse Opposite Designs:: Infinium I/II
              #   Just need to move forward by 1::
              Prb_SeqO1 = paste0(
                # cur_tib$up61, # Remove up61
                cur_tib$dn61,
                cur_tib$dn60,
                cur_tib$dn59,
                cur_tib$dn58 %>% stringr::str_sub(1,47) # 46 -> 47
              ), # %>% revCmp(), # No revCmp
              Prb_SeqO2 = paste0(
                # cur_tib$up61,
                # cur_tib$dn61, # Remove dn61
                cur_tib$dn60,
                cur_tib$dn59,
                cur_tib$dn58 %>% stringr::str_sub(1,48) # 47 -> 48
              ) # %>% revCmp(), # No revCmp
              
            )
          } else {
            stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported fr={fr}...{RET}"))
            return(NULL)
          }
          prb_cnt <- print_tib(prb_tib,funcTag, verbose,vt+4,tc, n="prb")

          cur_tib <- cur_tib %>% dplyr::left_join(prb_tib, by=cgn)
          ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
        }
      }
    }
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

fas_to_seq = function(tib, fas, 
                      file=NULL,
                      name,din,gen="na",
                      chr1="Chromosome", pos="Coordinate",
                      chr2="Chrom_Char", srd="F",
                      fwd_seq="Fwd_Temp_Seq",
                      iup_seq="Iup_Temp_Seq",
                      imp_seq="Imp_Temp_Seq",
                      ups_len=60, seq_len=122, nrec=0,
                      iupac=NULL,
                      ref_col="Ref",alt_col="Alt",iup_col="Iupac",
                      add_flank=FALSE,del="_", 
                      verbose=0,vt=4,tc=1,tt=NULL,
                      funcTag='fas_to_seq') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    if (nrec==0) {
      dna_seqs <- 
        Biostrings::readDNAStringSet(filepath=fas, format="fasta")
    } else {
      dna_seqs <- 
        Biostrings::readDNAStringSet(filepath=fas, format="fasta", nrec=nrec)
    }
    
    dna_maps <- dna_seqs %>% 
      names() %>% 
      stringr::str_remove(" .*$") %>% 
      stringr::str_remove("^chr") %>%
      tibble::tibble() %>% 
      purrr::set_names("Chrom_Char") %>% 
      dplyr::mutate(Idx=dplyr::row_number(),
                    Chrom_Char=paste0("chr",Chrom_Char)
      )
    
    pos_sym  <- rlang::sym(pos)
    name_sym <- rlang::sym(name)
    chr1_sym <- rlang::sym(chr1)
    chr2_sym <- rlang::sym(chr2)
    fwd_sym  <- rlang::sym(fwd_seq)
    iup_sym  <- rlang::sym(iup_seq)
    imp_sym  <- rlang::sym(imp_seq)
    din_sym  <- rlang::sym(din)
    
    ref_col_sym  <- rlang::sym(ref_col)
    alt_col_sym  <- rlang::sym(alt_col)
    iup_col_sym  <- rlang::sym(iup_col)
    
    chr_list <- tib %>% 
      dplyr::arrange(!!chr1_sym,!!pos_sym) %>%
      split(.[[chr1]])
    
    chr_maps <- dna_maps %>%
      split(.[[chr2]])
    
    chr_names <- chr_list %>% names()
    for (chr_str in chr_names) {
      cur_tib <- chr_list[[chr_str]]
      
      if (is.null(chr_maps[[chr_str]])) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Failed to find chr_str={chr_str} ",
                         "in map; Skipping...{RET}"))
        next
      }
      
      chr_idx <- chr_maps[[chr_str]] %>% 
        dplyr::filter(!!chr2_sym==chr_str) %>% 
        head(n=1) %>% pull(Idx) %>% as.integer()
      
      cur_key <- glue::glue("cur_tib_1=chr(str/idx)=({chr_str}/{chr_idx})")
      cur_cnt <- print_tib(cur_tib,funcTag, verbose,vt+4,tc, n=cur_key)

      # Subtract two for cg upstream 122mer formation and add two+two as well
      #
      cur_begs <- NULL
      cur_ends <- NULL
      cur_begs <- cur_tib %>% dplyr::pull(!!pos_sym) - ups_len - 2
      cur_ends <- cur_begs + seq_len - 1 + 4
      
      if (verbose>=vt+6) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_begs={RET}"))
        cur_begs %>% head() %>% print()
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_ends={RET}"))
        cur_ends %>% head() %>% print()
      }
      
      # Sub-string ref sequences::
      #
      ref_seqs <- NULL
      ref_seqs <- 
        stringr::str_sub( as.character(dna_seqs[[chr_idx]]), cur_begs, cur_ends) %>%
        stringr::str_to_upper()
      
      if (srd=="R") ref_seqs <- ref_seqs %>% cmpl()
      
      # Build individual parts of the template sequence::
      # up01.up02...up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn02.dn01
      #                             [  C   G  ]
      #
      ref_up01 <- stringr::str_sub(ref_seqs, 1,1)
      ref_up02 <- stringr::str_sub(ref_seqs, 2,2)
      ref_up58 <- stringr::str_sub(ref_seqs, 3,ups_len-2+2)
      
      ref_up59 <- stringr::str_sub(ref_seqs, ups_len-2+3,ups_len-2+3)
      ref_up60 <- stringr::str_sub(ref_seqs, ups_len-2+4,ups_len-2+4)
      ref_up61 <- stringr::str_sub(ref_seqs, ups_len-2+5,ups_len-2+5)
      
      iupac_vec <- ref_up61
      if (!is.null(iupac)) iupac_vec <- cur_tib %>% dplyr::pull(!!iupac)
      
      ref_dn61 <- stringr::str_sub(ref_seqs, ups_len-2+6,ups_len-2+6)
      ref_dn60 <- stringr::str_sub(ref_seqs, ups_len-2+7,ups_len-2+7)
      ref_dn59 <- stringr::str_sub(ref_seqs, ups_len-2+8,ups_len-2+8)
      
      ref_dn58 <- stringr::str_sub(ref_seqs, ups_len-2+9,ups_len-2+ups_len-2+8)
      ref_dn02 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+8,ups_len-2+ups_len-2+8)
      ref_dn01 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+9,ups_len-2+ups_len-2+9)
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ref_seqs({chr_str})={RET}"))
        ref_seqs %>% head(n=2) %>% print()
      }
      
      # Add all fields to current tibble::
      #
      cur_tib <- cur_tib %>%
        dplyr::mutate(
          !!fwd_sym := ref_seqs %>% addBrac(),
          
          up01=ref_up01,
          up02=ref_up02,
          up58=ref_up58,
          
          up59=ref_up59,
          up60=ref_up60,
          up61=ref_up61,
          
          iupac=iupac_vec,
          
          dn61=ref_dn61,
          dn60=ref_dn60,
          dn59=ref_dn59,
          
          dn58=ref_dn58,
          dn02=ref_dn02,
          dn01=ref_dn01,
          
          # Now we can assemble to optimal template::
          #
          !!iup_sym := paste0(up58,
                              up59,up60,
                              iupac,dn61,
                              dn60,dn59,
                              dn58) %>%
            addBrac(),
          !!imp_sym := paste0(up58,
                              up59,up60,
                              "CG",
                              dn60,dn59,
                              dn58) %>%
            addBrac(),
          
          # TBD:: Pretty sure this right, but really needs more testing...
          #
          Din_Str=stringr::str_to_lower(paste0(iupac,dn61)),
          Des_Din=dplyr::case_when(
            up61=='C' & dn61=='G' ~ "cg",
            up61=='C' & dn61!='G' ~ "ch",
            up61!='C' ~ "rs",
            TRUE ~ NA_character_
          )
        )
      cur_key <- glue::glue("cur_tib_2=chr(str/idx)=({chr_str}/{chr_idx})")
      cur_cnt <- print_tib(cur_tib,funcTag, verbose,vt+4,tc, n=cur_key)
      
      #  
      # TBD:: We can also break out tri-fecta sites if probe type == rs
      #
      ups_tib <- NULL
      dns_tib <- NULL
      
      if (add_flank) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Upstream CG Flank Seuqnces({chr_str})...{RET}"))
        
        ups_tib <- cur_tib %>% 
          dplyr::filter(!!din_sym == "rs") %>%
          dplyr::filter(up59=="C" & up60=="G") %>%
          dplyr::mutate(
            !!din_sym := "cg",
            !!fwd_sym := paste0(up01,up02,up58,up59,up60,up61,dn61,dn60,dn59,dn58) %>%
              stringr::str_sub(1,seq_len) %>%
              addBrac(),
            !!iup_sym := paste0(up01,up02,up58,up59,up60,iupac,dn61,dn60,dn59,dn58) %>%
              stringr::str_sub(1,seq_len) %>%
              addBrac(),
            !!imp_sym := paste0(up01,up02,up58,"CG",up61,dn61,dn60,dn59,dn58) %>%
              stringr::str_sub(1,seq_len) %>%
              addBrac(),
            !!ref_col_sym := "C",
            !!alt_col_sym := "T",
            !!iup_col_sym := "Y",
            Din_Str=stringr::str_to_lower(paste0(up59,up60)),
            Des_Din="cg",
            !!name_sym:=paste(!!name_sym,"CG-UP", sep='-'),
            Coordinate=as.integer(Coordinate-2)
          )
        
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Downstream CG Flank Seuqnces({chr_str})...{RET}"))
        
        dns_tib <- cur_tib %>% 
          dplyr::filter(!!din_sym == "rs") %>%
          dplyr::filter(dn61=="C" & dn60=="G") %>%
          dplyr::mutate(
            !!din_sym := "cg",
            !!fwd_sym := paste0(up58,up59,up60,up61,dn61,dn60,dn59,dn58,dn02) %>%
              stringr::str_sub(2) %>%
              addBrac(),
            !!iup_sym := paste0(up58,up59,up60,iupac,dn61,dn60,dn59,dn58,dn02) %>%
              stringr::str_sub(2) %>%
              addBrac(),
            !!imp_sym := paste0(up58,up59,up60,up61,"CG",dn59,dn58,dn02) %>%
              stringr::str_sub(2) %>%
              addBrac(),
            !!ref_col_sym := "C",
            !!alt_col_sym := "T",
            !!iup_col_sym := "Y",
            Din_Str=stringr::str_to_lower(paste0(dn61,dn60)),
            Des_Din="cg",
            !!name_sym:=paste(!!name_sym,"CG-DN", sep='-'),
            Coordinate=as.integer(Coordinate+2)
          )
      }
      
      # Add New Tri-fecta Probes::
      #  - NOTE:: These will have their names changed later after alignment...
      #
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joining and sorting...{RET}"))
      
      cur_tib <- dplyr::bind_rows(cur_tib, ups_tib, dns_tib) %>%
        dplyr::arrange(!!name_sym)
      
      ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
      
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Substring chr_str={chr_str}.{RET}{RET}"))
    }
    
    if (!is.null(file)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing improbe input TSV={file}...{RET}"))
      
      imp_col <- c("Seq_ID","Sequence","Genome_Build","Chromosome","Coordinate","CpG_Island")
      imp_tib <- ret_tib %>% 
        dplyr::mutate(Genome_Build=!!gen, CpG_Island="FALSE") %>%
        dplyr::select(dplyr::all_of(c(!!name, !!imp_seq, "Genome_Build", !!chr1, !!pos, "CpG_Island")) ) %>%
        purrr::set_names(imp_col)
      
      safe_write(x=imp_tib,type="tsv",file=file,funcTag=funcTag, 
                 verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} DIN Summary={RET}"))
      
      sum_tib <- ret_tib %>% 
        dplyr::group_by(up61,dn61,Des_Din) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      sum_tib %>% print(n=base::nrow(sum_tib))
    }
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret_fin")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# End of file
