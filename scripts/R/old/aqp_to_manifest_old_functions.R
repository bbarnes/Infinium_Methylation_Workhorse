
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
#                             BSMAP Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsp_to_grs = function(bsp, rds=NULL,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'bsp_to_grs'
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
    # TBD:: Temp Fix for Bsp_Pos == NA !!!
    #
    
    bsp <- bsp %>% dplyr::filter(!is.na(Bsp_Pos))
    
    ret_grs <- 
      GRanges(
        seqnames = Rle(paste0("chr",bsp$Bsp_Chr)),
        strand=Rle(stringr::str_sub( bsp$Bsp_Srd, 1,1 ) ),
        
        # Ord_Id  = bsp$Ord_Id,
        
        # Bsp_Add = bsp$Bsp_Add,
        # Bsp_Add = bsp$Address,
        
        # Bsp_Srd = bsp$Bsp_Srd,
        # Bsp_Src = bsp$Bsp_Src,
        
        Prb_Cgn = bsp$Prb_Cgn,
        Prb_Des = bsp$Prb_Des,
        
        Prb_Ord_Seq  = bsp$Prb_Seq,
        Prb_Aln_50U  = bsp$Prb_Seq %>%
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
        
        Bsp_Ref_Seq = bsp$Bsp_Ref,
        Bsp_Din_Scr = bsp$Bsp_Din_Scr,
        Bsp_Din_Ref = bsp$Bsp_Din_Ref,
        Bsp_Nxb_Ref = bsp$Bsp_Nxb_Ref,
        Bsp_Din_Bsc = bsp$Bsp_Din_Bsc,
        Bsp_Nxb_Bsc = bsp$Bsp_Nxb_Bsc,
        
        Bsp_Tag = bsp$Bsp_Tag,
        Bsp_Srd = bsp$Bsp_Srd,
        Bsp_Mis = bsp$Bsp_Mis,
        Bsp_Gap = bsp$Bsp_Gap,
        
        Bsp_Hit0 = bsp$Bsp_Hit0,
        Bsp_Hit1 = bsp$Bsp_Hit1,
        Bsp_Hit2 = bsp$Bsp_Hit2,
        Bsp_Hit3 = bsp$Bsp_Hit3,
        Bsp_Hit4 = bsp$Bsp_Hit4,
        Bsp_Hit5 = bsp$Bsp_Hit5,
        
        IRanges(start=bsp$Bsp_Pos,
                end=bsp$Bsp_Pos+1,
                names=bsp$Unique_ID
        )
      )
    
    if (!is.null(rds)) {
      if (verbose>=1)
        cat(glue::glue("[{funcTag}]: Writing Probe GRS RDS={rds}...{RET}"))
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

format_ORD = function(tib, idx_key="IDX", prb_key="Prb_Seq", uniq=TRUE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_ORD'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    idx_sym <- rlang::sym(idx_key)
    prb_sym <- rlang::sym(prb_key)
    
    add_col  <- c(idx_key,"Seq_ID","Prb_Des","Prb_Seq","Prb_Par","Prb_Col")
    sel_colA <- c(idx_key,"Assay_Design_Id","DesA","PrbA","PrbB","PrbCol")
    sel_colB <- c(idx_key,"Assay_Design_Id","DesB","PrbB","PrbA","PrbCol")
    
    tmp_tib <- tib %>%
      dplyr::select(-AlleleA_Probe_Id, -AlleleB_Probe_Id) %>%
      dplyr::rename(PrbA=AlleleA_Probe_Sequence,
                    PrbB=AlleleB_Probe_Sequence) %>%
      dplyr::mutate(
        PrbA=stringr::str_to_upper(PrbA),
        PrbB=stringr::str_to_upper(PrbB),
        DesA=dplyr::case_when(
          !is.na(PrbA) & !is.na(PrbB) ~ "U",
          !is.na(PrbA) &  is.na(PrbB) ~ "2",
          TRUE ~ NA_character_
        ),
        DesB=dplyr::case_when(
          !is.na(PrbA) & !is.na(PrbB) ~ "M",
          TRUE ~ NA_character_
        ),
        PrbCol=dplyr::case_when(
          DesA=="2" ~ NA_character_,
          DesA=="U" & DesB=="M" & Normalization_Bin=="A" ~ "R",
          DesA=="U" & DesB=="M" & Normalization_Bin=="B" ~ "G",
          TRUE ~ NA_character_
        )
      )
    tmp_cnt <- print_tib(tmp_tib,funcTag,verbose,vt+6,tc, n="tmp")
    
    ret_tib <- dplyr::bind_rows(
      tmp_tib %>% dplyr::filter(DesB=="M") %>% 
        dplyr::select(dplyr::all_of(sel_colB)) %>%
        purrr::set_names(add_col),
      tmp_tib %>% dplyr::filter(DesA=="U") %>% 
        dplyr::select(dplyr::all_of(sel_colA)) %>%
        purrr::set_names(add_col),
      tmp_tib %>% dplyr::filter(DesA=="2") %>% 
        dplyr::select(dplyr::all_of(sel_colA)) %>%
        purrr::set_names(add_col)
    ) %>%
      dplyr::arrange(plyr::desc(!!idx_sym))
    
    ret1_cnt <- print_tib(ret_tib,funcTag,verbose,vt+6,tc, n="ret1")
    if (uniq) ret_tib <- ret_tib %>% dplyr::distinct(!!prb_sym, .keep_all=TRUE)
    unq1_cnt <- print_tib(ret_tib,funcTag,verbose,vt+6,tc, n="retUnq")
    
    ret_tib <- ret_tib %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    ret_cnt <- print_tib(ret_tib,funcTag,verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

format_MAT = function(tib, idx_key="IDX", uniq=TRUE, trim=TRUE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_MAT'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    idx_sym <- rlang::sym(idx_key)
    
    ret_tib <- tib
    if (c("address_names") %in% names(tib)) {
      
      if (trim) {
        trim_cnt <- ret_tib %>% 
          dplyr::filter(!stringr::str_starts(address_names, '1')) %>% base::nrow()
        if (trim_cnt==0) {
          ret_tib <- ret_tib %>% 
            dplyr::mutate(address_names=as.numeric(
              stringr::str_remove(stringr::str_remove(address_names, '^1'), '^0+')) )
        } else {
          stop(glue::glue("{RET}[{funcTag}]: ERROR: Attempting to trim new tango fomrat, ",
                          "but format doesn't match; trim_cnt={trim_cnt}!!!{RET}{RET}"))
          return(NULL)
        }
      }
      
      ret_tib <- ret_tib %>% 
        dplyr::rename(Address=address_names,
                      Sequence=bo_seq)
      
    } else if (c("Address") %in% names(tib)) {
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
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed to match 'addres_names' or 'Address'!!!{RET}{RET}"))
      return(NULL)
    }
    
    ret_tib <- ret_tib %>% 
      dplyr::mutate(Sequence=stringr::str_to_upper(Sequence),
                    tan_seq=stringr::str_sub(Sequence,1,45),
                    Prb_Seq=stringr::str_sub(Sequence,46)
      ) %>%
      dplyr::filter(!is.na(Address)) %>%
      dplyr::select(Address,!!idx_sym,Prb_Seq, dplyr::everything()) %>% 
      dplyr::arrange(plyr::desc(!!idx_sym))
    
    
    if (uniq) ret_tib <- ret_tib %>% dplyr::distinct(Address, .keep_all=TRUE)
    
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

format_AQP = function(tib, idx_key="IDX", uniq=TRUE,filt=TRUE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_AQP'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    idx_sym <- rlang::sym(idx_key)
    
    # Determine if its AQP & PQC::
    #
    ret_tib <- tib %>% 
      dplyr::mutate(Address=as.numeric(
        stringr::str_remove(stringr::str_remove(Address, '^1'), '^0+')) )
    
    if (c("Decode_Status") %in% names(tib)) {
      ret_tib <- ret_tib %>% 
        dplyr::select(Address,Decode_Status,!!idx_sym) %>%
        dplyr::rename(aqp_val=Decode_Status) %>%
        dplyr::mutate(aqp_src="AQP")
    } else if (c("Status") %in% names(tib)) {
      ret_tib <- ret_tib %>% 
        dplyr::select(Address,Status,!!idx_sym) %>%
        dplyr::rename(aqp_val=Status) %>%
        dplyr::mutate(aqp_src="PQC")
    } else {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Unable to match AQP/PQC!{RET}{RET}"))
      return(NULL)
    } %>% dplyr::arrange(plyr::desc(!!idx_sym))
    
    pre_cnt <- ret_tib %>% base::nrow()
    ret_tib <- ret_tib %>% 
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    if (filt) ret_tib <- ret_tib %>% dplyr::filter(aqp_val==0)
    if (uniq) ret_tib <- ret_tib %>% dplyr::distinct(Address, .keep_all=TRUE)
    
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
#                      OLD Manifest File Methods:: Loading
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_manifestBuildFile = function(file, field="Assay_Design_Id", cols=NULL,
                                  # aqp=NULL, bpn=NULL,
                                  n_max=50, guess=100000,
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
    } else if (file_type=="tsv" || file_type=="txt" || file_type=="match") {
      ret_tib <- suppressMessages(suppressWarnings(
        readr::read_tsv(file, skip=skip_cnt, guess_max=guess, skip_empty_rows=FALSE)))
    } else {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Unsupported fileType={file_type}!!!{RET}{RET}"))
      return(NULL)
    }
    ret1_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n="ret1")
    
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
      tibble::as_tibble() %>% 
      dplyr::mutate(value=stringr::str_remove(value, "^\\s+") )
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="raw")
    
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
          # strand=Rle(stringr::str_sub( bsp$Bsp_Srd, 1,1 ) ),
          
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
        # strand=Rle(stringr::str_sub( bsp$Bsp_Srd, 1,1 ) ),
        
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
