
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Alignment Methods::
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
#                         BSP Alignment Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_bsmap = function(exe, fas, gen, bsp, 
                     add=NULL, bed=NULL, org=NULL,
                     
                     join_key="Aln_Key",join_type="inner",
                     des_key="Ord_Des",din_key="Ord_Din",
                     sort=TRUE,
                     
                     opt=NULL, lan=NULL, run=FALSE,
                     verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'run_bsmap'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; bsp={bsp}.{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  aln_ssh <- NULL
  ret_tsv <- NULL
  sys_ret <- "Not_Run"
  stime <- system.time({
    if (is.null(opt) || length(opt)==0)
      opt="-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
    
    build_file_dir(bsp)
    aln_ssh <- paste(bsp, 'run.sh', sep='.')
    ret_tsv <- paste(bsp, 'tsv.gz', sep='.')
    
    aln_cmd <- glue::glue("{exe} -a {fas} -d {gen} {opt} -o {bsp}{RET}",
                          "cat {bsp} | ",
                          # "cut -f 1,2,4-11 | ",
                          "gzip -c - > {ret_tsv}{RET}",
                          "rm -f {bsp}{RET}")
    
    safe_write(aln_cmd,"line",aln_ssh,funcTag=funcTag,
               verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    Sys.chmod(paths=aln_ssh, mode="0777")
    
    if (run) {
      exe_method="Running"
      if (!is.null(lan)) {
        exe_method="Launching"
        aln_ssh <- paste(lan,aln_ssh, sep=' ')
      }
      
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} {exe_method}; CMD={aln_cmd}...{RET}"))
      sys_ret <- base::system(aln_ssh)
      if (sys_ret!=0) {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: sys_ret({sys_ret}) != 0!!!{RET}{RET}"))
        return(NULL)
      }
      
      if (!is.null(add)) {
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                      Join Address and Alignment Data::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        ret_tib <-
          join_bsmap(
            add=add, bed=bed, org=org, file=ret_tsv,
            join_key=join_key,join_type="inner",
            prb_des_key=des_key,prb_din_key=din_key,
            sort=TRUE,
            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        
        ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
      }
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; sys={sys_ret}; elapsed={etime}.{RET}{RET}"))
  
  if (!is.null(ret_tib)) return(ret_tib)
  if (run) return(ret_tsv)
  aln_ssh
}

load_bsmap = function(bsp, sort=FALSE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_bsmap'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    bsp_cols <- cols(
      Bsp_Key  = col_character(),
      Bsp_Seq  = col_character(),
      Bsp_Qual = col_character(),
      Bsp_Tag  = col_character(),
      
      Bsp_Chr  = col_character(),
      Bsp_Beg  = col_integer(),
      Bsp_Srd  = col_character(),
      Bsp_Mis  = col_integer(),
      
      Bsp_Ref  = col_character(),
      Bsp_Gap  = col_integer(),
      
      Bsp_Str  = col_character()
    )
    
    # Load BSP
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading BSP={bsp}...{RET}"))
    
    ret_tib <- 
      readr::read_tsv(bsp,col_names=names(bsp_cols$cols),col_types=bsp_cols) %>%
      dplyr::select(-Bsp_Qual) %>%
      dplyr::mutate(Bsp_Chr=paste0('chr',stringr::str_remove(Bsp_Chr,'^chr')))
    
    # Sort by genomic position::
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(Bsp_Chr, Bsp_Beg)
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

join_bsmap = function(add, bsp=NULL, bed=NULL, org=NULL, file=NULL, 
                      join_key, join_type="inner", 
                      prb_des_key="Ord_Des", prb_din_key="Ord_Din",
                      fields=NULL, sort=FALSE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'join_bsmap'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}",
                   "[{funcTag}]:{tabsStr}  join_key={join_key}{RET}",
                   "[{funcTag}]:{tabsStr} join_type={join_type}{RET}")
    )
  
  if (is.null(bsp) && is.null(file)) {
    stop(glue::glue("[{funcTag}]:{tabsStr} ERROR: bsp AND file are null!!!{RET}{RET}"))
    return(NULL)
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    imp_col_vec <- c("Address","Ord_Des","Ord_Din",
                     "Imp_Chr","Imp_Pos","Imp_Cgn",
                     "Imp_FR","Imp_TB","Imp_CO","Imp_Nxb",
                     "Bsp_Din_Ref","Bsp_Din_Scr","Aln_Prb")
    
    # Load from file if provided::
    if (!is.null(file))
      bsp <- load_bsmap(bsp=file, sort=sort,
                        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    if (is.null(bsp)) {
      stop(glue::glue("[{funcTag}]:{tabsStr} ERROR: bsp is null!!!{RET}{RET}"))
      return(ret_tib)
    }
    print_tib(add,funcTag, verbose,vt+4,tc, n="address_tib")
    
    prb_des_sym  <- rlang::sym(prb_des_key)
    prb_din_sym  <- rlang::sym(prb_din_key)
    bsp_join_key <- names(bsp)[1]
    join_key_sym <- rlang::sym(join_key)
    bsp_join_sym <- rlang::sym( bsp_join_key )
    
    # Rename Bsp_Key to join key and join::
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Joining source data...{RET}",
                     "[{funcTag}]:{tabsStr} Join Keys::{RET}",
                     "[{funcTag}]:{tabsStr} add_join_key={join_key}{RET}",
                     "[{funcTag}]:{tabsStr} bsp_join_key={bsp_join_key}{RET}{RET}" )
      )
    
    if (join_type=="inner") {
      ret_tib <- bsp %>%
        dplyr::rename(!!join_key_sym := !!bsp_join_sym) %>%
        dplyr::inner_join(add, by=join_key) # %>% head(n=100)
    } else if (join_type=="left") {
      ret_tib <- bsp %>%
        dplyr::rename(!!join_key_sym := !!bsp_join_sym) %>%
        dplyr::left_join(add, by=join_key) # %>% head(n=100)
    } else if (join_type=="right") {
      ret_tib <- bsp %>%
        dplyr::rename(!!join_key_sym := !!bsp_join_sym) %>%
        dplyr::right_join(add, by=join_key) # %>% head(n=100)
    } else if (join_type=="full") {
      ret_tib <- bsp %>%
        dplyr::rename(!!join_key_sym := !!bsp_join_sym) %>%
        dplyr::full_join(add, by=join_key) # %>% head(n=100)
    } else {
      stop(glue::glue("[{funcTag}]:{tabsStr} ERROR: Unsupported join_type={join_type}!!!{RET}{RET}"))
      return(ret_tib)
    }
    
    #
    # Calculate new fields from join::
    #
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Calculating new fields from joined data...{RET}"))
    
    ret_tib <- ret_tib %>%
      dplyr::mutate(
        # Generate bisulfite converted Reference Sequences::
        Bsp_RefU=bscUs(Bsp_Ref,uc=TRUE),
        Bsp_RefM=bscMs(Bsp_Ref,uc=TRUE),
        Bsp_RefD=bscDs(Bsp_Ref,uc=TRUE),
        
        # Add Alignment Orientation wrt Probe Sequence/BSP Seq::
        #  This is used to determine NxB, CpG matching, etc.
        Bsp_Prb_Dir=dplyr::case_when(
          Aln_Prb==Bsp_Seq ~ "f",
          Aln_Rev==Bsp_Seq ~ "r",
          TRUE ~ NA_character_),
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                      Di-Nucleotide & Next Base:: Reference
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        # Extract Potential Next Base from Reference Sequence::
        NB_F1=stringr::str_sub(Bsp_Ref,-2,-2),
        NB_R1=stringr::str_sub(Bsp_Ref, 2, 2),
        NB_F2=stringr::str_sub(Bsp_Ref,-1,-1),
        NB_R2=stringr::str_sub(Bsp_Ref, 1, 1),
        
        # Extract Potential CG di-nucleotide from Reference Sequence::
        CG_F1=stringr::str_sub(Bsp_Ref,-4,-3),
        CG_R1=stringr::str_sub(Bsp_Ref, 3, 4),
        CG_F2=stringr::str_sub(Bsp_Ref,-3,-2),
        CG_R2=stringr::str_sub(Bsp_Ref, 2, 3),
        
        # Set Expected target Next Base from Alignment Orientation::
        Bsp_Nxb_Ref=dplyr::case_when(
          !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F1,
          !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R1,
          
          !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F1,
          !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R1,
          
          !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F2,
          !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R2,
          TRUE ~ NA_character_
        ), # %>% stringr::str_to_upper(),
        
        # Set Expected target CG di-nucleotide from Alignment Orientation::
        Bsp_Din_Ref=dplyr::case_when(
          !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F1,
          !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R1,
          
          !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F1,
          !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R1,
          
          !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F2,
          !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R2,
          TRUE ~ NA_character_
        ) %>% stringr::str_to_upper(),
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #             Di-Nucleotide & Next Base:: Bisulfite Converted
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        # Extract Potential Next Base from Reference Sequence::
        NB_F1M=stringr::str_sub(Bsp_RefM,-2,-2),
        NB_R1M=stringr::str_sub(Bsp_RefM, 2, 2),
        
        NB_F1U=stringr::str_sub(Bsp_RefU,-2,-2),
        NB_R1U=stringr::str_sub(Bsp_RefU, 2, 2),
        
        NB_F2D=stringr::str_sub(Bsp_RefD,-1,-1),
        NB_R2D=stringr::str_sub(Bsp_RefD, 1, 1),
        
        # Extract Potential CG di-nucleotide from Reference Sequence::
        CG_F1M=stringr::str_sub(Bsp_RefM,-4,-3),
        CG_R1M=stringr::str_sub(Bsp_RefM, 3, 4),
        
        CG_F1U=stringr::str_sub(Bsp_RefU,-4,-3),
        CG_R1U=stringr::str_sub(Bsp_RefU, 3, 4),
        
        CG_F2D=stringr::str_sub(Bsp_RefD,-3,-2),
        CG_R2D=stringr::str_sub(Bsp_RefD, 2, 3),
        
        # Set Expected target Next Base from Alignment Orientation::
        Bsp_Nxb_Bsc=dplyr::case_when(
          !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F1M,
          !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R1M,
          
          !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F1U,
          !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R1U,
          
          !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F2D,
          !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R2D,
          TRUE ~ NA_character_
        ), # %>% stringr::str_to_upper(),
        
        # Set Expected target CG di-nucleotide from Alignment Orientation::
        Bsp_Din_Bsc=dplyr::case_when(
          !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F1M,
          !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R1M,
          
          !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F1U,
          !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R1U,
          
          !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F2D,
          !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R2D,
          TRUE ~ NA_character_
        ) %>% stringr::str_to_upper(),
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                      Correct CG# Alignment Position::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        #
        # Update Correct Genomic CG# Location based on alignment orientation::
        #
        # TBD:: Need to correct for non-CpG sites::
        #       - rs:: [+-] => +0
        #              [--] => +1
        #
        Bsp_Pos=dplyr::case_when(
          # prb_din_sym=="cg":: Original Code for only CpG! TO BE DELETED...
          #
          # !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          # !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          # 
          # !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          # !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          # 
          # !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          # !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          #
          # prb_din_sym=="cg"::
          #
          !!prb_din_sym=="cg" & !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          !!prb_din_sym=="cg" & !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="cg" & !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          !!prb_din_sym=="cg" & !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="cg" & !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="cg" & !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          #
          # prb_din_sym=="ch":: for now assume ch==rs
          #
          !!prb_din_sym=="ch" & !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="ch" & !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="ch" & !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="ch" & !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="ch" & !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          !!prb_din_sym=="ch" & !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          #
          # prb_din_sym=="rs"::
          #
          # !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          # !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          # 
          # !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          # !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          # 
          # !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          # !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          #
          # prb_din_sym=="rs":: --/+-
          #
          !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="--") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="+-") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="--") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="+-") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="--") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="+-") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          # 
          # prb_din_sym=="rs":: --/+-
          #
          !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 0,
          
          TRUE ~ NA_real_
        ) %>% as.integer(),
        Bsp_Din_Scr=dplyr::case_when(
          Bsp_Srd=="--" & Bsp_Din_Ref=="CG" ~ 0,
          Bsp_Srd=="+-" & Bsp_Din_Ref=="CG" ~ 1,
          Bsp_Srd=="-+" & stringr::str_starts(Bsp_Din_Ref,"G") ~ 2,
          Bsp_Srd=="++" & stringr::str_ends(Bsp_Din_Ref,"C") ~ 3,
          TRUE ~ 9) %>% as.integer()
      ) %>% 
      # Create Unique Aln Key for multiple hits::
      #
      dplyr::group_by(!!join_key_sym) %>%
      dplyr::mutate(Aln_Key_Unq=paste(!!join_key_sym, dplyr::row_number(), sep="_")) %>%
      dplyr::ungroup() %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    # Sort by genomic position::
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(Bsp_Chr, Bsp_Beg)
    
    if (!is.null(bed)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Adding cgn bed...{RET}"))
      
      ret_tib <- bed %>% 
        dplyr::right_join(ret_tib, by=c("Imp_Chr"="Bsp_Chr", "Imp_Pos"="Bsp_Pos")) %>%
        dplyr::mutate(
          Imp_FR=dplyr::case_when(
            Bsp_Srd=="+-" ~ "R",
            Bsp_Srd=="--" ~ "F",
            Bsp_Srd=="++" ~ "R",
            Bsp_Srd=="-+" ~ "F",
            TRUE ~ NA_character_
          ),
          Imp_CO=dplyr::case_when(
            Bsp_Srd=="+-" ~ "C",
            Bsp_Srd=="--" ~ "C",
            Bsp_Srd=="++" ~ "O",
            Bsp_Srd=="-+" ~ "O",
            TRUE ~ NA_character_
          ),
          Imp_TB=dplyr::case_when(
            Imp_Top_Srd=="+" & Imp_FR=="F" ~ "T",
            TRUE ~ "B"
          ),
          Imp_Nxb=stringr::str_to_upper(Bsp_Nxb_Ref)
        )
    }
    
    if (!is.null(org)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Adding org bed...{RET}"))
      
      ret_tib <- ret_tib %>%
        dplyr::left_join(
          org, 
          by=c("Aln_P49"="Org_49P",
               "Ord_Des"="Org_Des",
               "Imp_Chr"="Org_Chr",
               "Imp_Pos"="Org_Pos",
               "Imp_FR"="Org_FR",
               "Imp_TB"="Org_TB",
               "Imp_CO"="Org_CO")
        )
    }
    
    ret_tib <- ret_tib  %>%
      dplyr::mutate(Cgn_Id=stringr::str_remove(Aln_Key_Unq, "-.*$")) %>%
      dplyr::select(dplyr::any_of(imp_col_vec),
                    dplyr::everything()) %>%
      clean_tibble()

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
