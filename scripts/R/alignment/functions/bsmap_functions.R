
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
#                         Load CGN/Position Database::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_cgn_dB = function(file,
                       verbose=0,vt=3,tc=1,tt=NULL,
                       funcTag='load_cgn_dB') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr} file={file}.{RET}"))
  }
  
  cgd_bed_cols <-
    cols(
      Cgd_Chr = col_character(),
      Cgd_Pos = col_integer(),
      Cgd_Cgn = col_integer(),
      Cgd_Top = col_character()
    )
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ret_tib <- suppressMessages(suppressWarnings(
      readr::read_tsv(run$cgn_bed_tsv, 
                      col_names=names(cgd_bed_cols$cols),
                      col_types=cgd_bed_cols) )) %>% 
      dplyr::mutate(Cgd_Chr=paste0("chr",Cgd_Chr))
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; ",
                   "Return Count={ret_cnt}; ",
                   "elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ",
                   "----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         BSP Alignment Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsp_mapping_workflow = 
  function(ref, can, ord, seq, out, exe, cgn=NULL, csv=NULL, canonical=NULL,
           sort=TRUE, light=FALSE, reload=FALSE, full=FALSE, merge=TRUE,
           join_key="Aln_Key", join_type="inner",
           des_key="Ord_Des", din_key="Ord_Din",
           opt="-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R",
           verbose=0,vt=3,tc=1,tt=NULL, funcTag='bsp_mapping_workflow') {
    
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}   ref={ref}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}   can={can}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}   out={out}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}   exe={exe}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}   csv={csv}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}   opt={opt}{RET}{RET}"))
      
      cat(glue::glue("[{funcTag}]:{tabsStr}   sort={sort}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}  light={light}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr} reload={reload}{RET}{RET}"))
    }
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- base::system.time({
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                      Probe Alignment:: BSMAP
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ret_tib <-run_bsmap(ref=ref,
                          can=can,
                          out=out,
                          exe=exe,
                          sort=sort,
                          light=light, 
                          reload=reload,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                      Join Address and Alignment Data::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ret_tib <- join_bsmap(bsp=ret_tib,
                            ord=ord,
                            cgn=cgn,
                            join_key=join_key,
                            join_type=join_type,
                            des_key=des_key,
                            din_key=din_key,
                            sort=sort,
                            full=full,
                            csv=NULL,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                         Consolidate/Assign CGNs::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ret_tib <- assign_cgn(ord=ord,
                            bsp=ret_tib,
                            seq=seq, 
                            can=canonical,
                            csv=csv,
                            merge=TRUE,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                           Alignment Summary::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      bsp_hit_csv <- NULL
      bsp_hit_sum <- ret_tib %>% 
        dplyr::group_by(Address) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>% 
        dplyr::group_by(Count) %>% 
        dplyr::summarise(His_Count=n(), .groups="drop")
      print_tib(bsp_hit_sum, funcTag, verbose,vt, l=20, n="bsp_hit_sum")

      # Top Ranked Offfenders::
      top_hit_csv <- NULL
      top_hit_sum <- ret_tib %>% 
        dplyr::group_by(Address) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>%
        dplyr::filter(Count!=1) %>%
        dplyr::arrange(-Count)
      # print(top_hit_sum, n=base::nrow(top_hit_sum))
      
      if (!is.null(csv)) {
        sum_out_dir <- base::dirname(csv)
        bsp_hit_csv <- file.path(sum_out_dir, "bsp-alignment-histogram.csv.gz")
        top_hit_csv <- file.path(sum_out_dir, "bsp-multi-alignment-loci.csv.gz")
      }
      
      sum_cnt <- safe_write(x=bsp_hit_sum, file=bsp_hit_csv, funcTag=funcTag,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      top_cnt <- safe_write(x=top_hit_sum, file=top_hit_csv, funcTag=funcTag,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
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

run_bsmap = function(ref, can, out, exe, 
                     sort=FALSE, 
                     light=FALSE,
                     reload=FALSE,
                     opt="-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R",
                     verbose=0,vt=3,tc=1,tt=NULL,
                     funcTag='run_bsmap') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr}   ref={ref}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   can={can}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   out={out}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   exe={exe}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   exe={exe}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   opt={opt}{RET}{RET}"))
    
    cat(glue::glue("[{funcTag}]:{tabsStr}   sort={sort}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  light={light}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} reload={reload}{RET}{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    if (!dir.exists(out)) dir.create(out, recursive=TRUE)
    
    out_bsp <- file.path(out, "probe.bsp")
    bsp_ssh <- file.path(out, "probe.bsp.run.sh")
    bsp_tsv <- file.path(out, "probe.bsp.tsv.gz")
    
    if (!reload || !file.exists(bsp_tsv)) {
      add_cmd <- ""
      if (light) add_cmd <- "cut -f 1,2,4-11 | "
      
      bsp_cmd <- glue::glue("{exe} -a {can} -d {ref} {opt} -o {out_bsp}{RET}",
                            "cat {out_bsp} | {add_cmd} gzip -c -> {bsp_tsv}{RET}",
                            "rm -f {out_bsp}{RET}")
      
      out_cnt <- safe_write(bsp_cmd,"line",bsp_ssh,funcTag=funcTag,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      Sys.chmod(paths=bsp_ssh, mode="0777")
      
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} CMD={bsp_cmd}...{RET}"))
      sys_ret <- 0
      sys_ret <- base::system(bsp_ssh)
      
      if (sys_ret!=0) {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: sys_ret({sys_ret}) != 0!!!{RET}{RET}"))
        return(ret_tib)
      }
    }
    
    ret_tib <- load_bsmap(file=bsp_tsv, sort=sort, light=light,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
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

load_bsmap = function(file,
                      sort=FALSE,
                      light=FALSE,
                      verbose=0,vt=3,tc=1,tt=NULL,
                      funcTag='load_bsmap') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr}  file={file}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  sort={sort}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} light={light}{RET}{RET}"))
  }
  
  stopifnot(file.exists(file))
  
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
    
    low_cols <- cols(
      Bsp_Key  = col_character(),
      Bsp_Seq  = col_character(),
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
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading BSP={file}...{RET}"))
    
    if (light) {
      ret_tib <- 
        readr::read_tsv(file,col_names=names(low_cols$cols),col_types=low_cols) %>%
        dplyr::mutate(Bsp_Chr=paste0('chr',stringr::str_remove(Bsp_Chr,'^chr')))
      
    } else {
      ret_tib <- 
        readr::read_tsv(file,col_names=names(bsp_cols$cols),col_types=bsp_cols) %>%
        dplyr::select(-Bsp_Qual) %>%
        dplyr::mutate(Bsp_Chr=paste0('chr',stringr::str_remove(Bsp_Chr,'^chr')))
    }
    
    # Sort by genomic position::
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(Bsp_Chr, Bsp_Beg)
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

join_bsmap = function(bsp, ord, cgn, 
                      join_key,
                      join_type="inner", 
                      des_key="Ord_Des",
                      din_key="Ord_Din",
                      sort=TRUE,
                      full=FALSE,
                      csv=NULL,
                      verbose=0,vt=3,tc=1,tt=NULL,
                      funcTag='join_bsmap') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr}   des_key={des_key}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   din_key={din_key}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  join_key={join_key}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} join_type={join_type}{RET}{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # Excessive fields to be removed generally
    rem_sel_vec <- 
      c("CG_F1","CG_R1","CG_F2","CG_R2","CG_F1M",
        "CG_R1M","CG_F1U","CG_R1U","CG_F2D","CG_R2D",
        "NB_F1","NB_R1","NB_F2","NB_R2","NB_F1M",
        "NB_R1M","NB_F1U","NB_R1U","NB_F2D","NB_R2D",
        "Bsp_RefU","Bsp_RefM","Bsp_RefD")
    
    # Order of fields to displayed
    imp_col_vec <- c("Address","Ord_Des","Ord_Din",
                     "Bsp_Chr","Bsp_Pos","Bsp_Cgn",
                     "Bsp_FR","Bsp_TB","Bsp_CO","Bsp_Nxb",
                     "Ord_Prb",
                     "Aln_Key_Unq",
                     "Bsp_Nxb_Ref","Bsp_Nxb_Bsc",
                     "Bsp_Din_Ref","Bsp_Din_Bsc",
                     "Bsp_Prb_Dir","Bsp_Din_Scr",
                     "Aln_Prb")
    
    prb_des_sym  <- rlang::sym(des_key)
    prb_din_sym  <- rlang::sym(din_key)
    
    ord_join_key <- join_key
    ord_join_sym <- rlang::sym(ord_join_key)
    
    bsp_join_key <- names(bsp)[1]
    bsp_join_sym <- rlang::sym( bsp_join_key )
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}    join_type={join_type}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}      des_key={des_key}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}      din_key={din_key}{RET}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr} bsp_join_key={bsp_join_key}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr} ord_join_key={ord_join_key}{RET}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Join BSP Alignment with Probe Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    bsp_key <- glue::glue("bsp-tib({funcTag})")
    bsp_cnt <- print_tib(bsp,funcTag, verbose,vt+4,tc, n=bsp_key)
    
    ord_key <- glue::glue("ord-tib({funcTag})")
    ord_cnt <- print_tib(ord,funcTag, verbose,vt+4,tc, n=ord_key)
    
    if (join_type=="inner") {
      ret_tib <- bsp %>%
        dplyr::rename(!!ord_join_sym := !!bsp_join_sym) %>%
        dplyr::inner_join(ord, by=ord_join_key) # %>% head(n=100)
    } else if (join_type=="left") {
      ret_tib <- bsp %>%
        dplyr::rename(!!ord_join_sym := !!bsp_join_sym) %>%
        dplyr::left_join(ord, by=ord_join_key) # %>% head(n=100)
    } else if (join_type=="right") {
      ret_tib <- bsp %>%
        dplyr::rename(!!ord_join_sym := !!bsp_join_sym) %>%
        dplyr::right_join(ord, by=ord_join_key) # %>% head(n=100)
    } else if (join_type=="full") {
      ret_tib <- bsp %>%
        dplyr::rename(!!ord_join_sym := !!bsp_join_sym) %>%
        dplyr::full_join(ord, by=ord_join_key) # %>% head(n=100)
    } else {
      stop(glue::glue("[{funcTag}]:{tabsStr} ERROR: Unsupported join_type={join_type}!!!{RET}{RET}"))
      return(ret_tib)
    }
    ret_key <- glue::glue("bsp/ord-join({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      Calculate New Fields from Join::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Calculating new fields from joined data...{RET}"))
    
    ret_tib <- ret_tib %>%
      dplyr::mutate(
        
        # Add more common FR/CO values::
        Bsp_FR=dplyr::case_when(
          Bsp_Srd=="+-" ~ "R",
          Bsp_Srd=="--" ~ "F",
          Bsp_Srd=="++" ~ "R",
          Bsp_Srd=="-+" ~ "F",
          TRUE ~ NA_character_
        ),
        Bsp_CO=dplyr::case_when(
          Bsp_Srd=="+-" ~ "C",
          Bsp_Srd=="--" ~ "C",
          Bsp_Srd=="++" ~ "O",
          Bsp_Srd=="-+" ~ "O",
          TRUE ~ NA_character_
        ),
        
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
        # TBD:: Need to correct for non-CpG sites:: I think this is done???
        #       - rs:: [+-] => +0
        #              [--] => +1
        #
        Bsp_Pos=dplyr::case_when(
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
          # prb_din_sym=="mu"::
          #
          !!prb_din_sym=="mu" & !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          !!prb_din_sym=="mu" & !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="mu" & !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          !!prb_din_sym=="mu" & !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="mu" & !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="mu" & !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          #
          # prb_din_sym=="bc"::
          #
          !!prb_din_sym=="bc" & !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          !!prb_din_sym=="bc" & !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="bc" & !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          !!prb_din_sym=="bc" & !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="bc" & !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="bc" & !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
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
          TRUE ~ 9) %>% as.integer(),
        Bsp_Nxb=stringr::str_to_upper(Bsp_Nxb_Ref)
      ) %>% 
      # Create Unique Aln Key for multiple hits::
      #
      dplyr::group_by(!!ord_join_sym) %>%
      dplyr::mutate(
        Bsp_Rank=dplyr::row_number(),
        Aln_Key_Unq=paste(!!ord_join_sym, Bsp_Rank, sep="_")) %>%
      dplyr::ungroup() %>%
      clean_tibble()
    
    if (!full) ret_tib <- ret_tib %>% 
      dplyr::select(- dplyr::all_of(rem_sel_vec))
    
    # Sort by genomic position::
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(Bsp_Chr, Bsp_Beg)
    ann_key <- glue::glue("ret-annotated({funcTag})")
    ann_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ann_key)
    
    ret_key <- glue::glue("calculated-fields({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Done. Calculating new fields.{RET}{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      Annotate CGN by Alignment Position::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (tibble::is_tibble(cgn)) {
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} Usung tibble={RET}"))
        print(cgn)
      }
    } else if (file.exists(cgn) && !dir.exists(cgn)) {
      if (verbose>=vt+4)
        cat(glue::glue("[{funcTag}]:{tabsStr} Loading from file={cgn}.{RET}"))
      cgn <- load_cgn_dB(file=cgn, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
    } else {
      cat(glue::glue("{RET}[{funcTag}]:{tabsStr} Warning: Unknown cgn type!{RET}"))
      cgn <- NULL
    }
    
    if (!is.null(cgn)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Adding cgn by alignment position...{RET}"))
      
      cgn <- cgn %>% purrr::set_names(c("Bsp_Chr","Bsp_Pos","Bsp_Cgn","Cgd_Top"))
      
      mat_tib <- dplyr::bind_rows(
        dplyr::inner_join(
          ret_tib %>% dplyr::select(Aln_Key_Unq,Bsp_Chr,Bsp_Pos),
          cgn %>% dplyr::mutate(Cgd_Nuc="up"),
          by=c("Bsp_Chr","Bsp_Pos"),
          suffix=c("_bsp","_cgn")),
        dplyr::inner_join(
          ret_tib %>% dplyr::select(Aln_Key_Unq,Bsp_Chr,Bsp_Pos),
          cgn %>% dplyr::mutate(Bsp_Pos=Bsp_Pos+1, Cgd_Nuc="dn"),
          by=c("Bsp_Chr","Bsp_Pos"),
          suffix=c("_bsp","_cgn"))
      ) %>% dplyr::distinct(Aln_Key_Unq, .keep_all = TRUE)
      
      ret_tib <- ret_tib %>% 
        dplyr::left_join(mat_tib, by=c("Aln_Key_Unq","Bsp_Chr","Bsp_Pos")) %>%
        dplyr::mutate(
          Bsp_TB=dplyr::case_when(
            Cgd_Top=="+" & Bsp_FR=="F" ~ "T",
            Cgd_Top=="-" & Bsp_FR=="R" ~ "T",
            
            Cgd_Top=="+" & Bsp_FR=="R" ~ "B",
            Cgd_Top=="-" & Bsp_FR=="F" ~ "B",
            TRUE ~ NA_character_
          )
        )
      
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Done. Adding cgn by position.{RET}{RET}"))
    }
    
    ret_tib <- ret_tib %>% 
      dplyr::select(dplyr::any_of(imp_col_vec),
                    dplyr::everything()) %>% clean_tibble()
    
    out_cnt <- safe_write(x=ret_tib, file=csv, funcTag=funcTag,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
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
#                          BSMAP Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# This isn't ready yet, but it is another way to validate CG#'s via coordinate
#   alignment with: 
#   BED File: /Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-input/GRCh37.cgn.bed.gz
#
if (FALSE) {
  
  bsmap_to_bed = function(tib,bed,
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='bsmap_to_bed') {
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- base::system.time({
      
      
      
      safe_write(bed_tib,"tsv",bed,funcTag=funcTag,
                 vt=vt+1,tc=tc+1,tt=tt)
      
      
      # ret_cnt <- ret_tib %>% base::nrow()
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    ret_tib
  }
}



# End of file
