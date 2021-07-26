
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Alignment Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

COM  <- ","
TAB  <- "\t"
RET  <- "\n"
RET2 <- "\n\n"
BNG  <- "|"
BRK  <- paste0("# ",
               paste(rep("-----",6),collapse=" "),"|",
               paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='template_func') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   funcTag={funcTag}.{RET}"))
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
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         BSP Alignment Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsp_mapping_workflow = 
  function(ref_fas, 
           ref_tib = NULL,
           can_fas = NULL,
           can_tib = NULL,
           seq_tib = NULL,
           
           cgn_src = NULL,
           canonical = NULL,
           
           ids_key = "Prb_Key",
           unq_key = "Unq_Key",
           prb_key = "Prb_Seq",
           des_key = "Ord_Des",
           din_key = "Ord_Din",
           
           join_key="Prb_Key",  # Previously Aln_Key
           join_type="inner",
           
           sort   = TRUE, 
           full   = FALSE, 
           merge  = TRUE,
           light  = FALSE, 
           reload = FALSE, 
           
           bsp_exe,
           bsp_opt = "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R",
           
           out_csv=NULL, out_dir, run_tag, 
           re_load=FALSE, pre_tag=NULL,
           end_str='csv.gz', sep_chr='.',
           verbose=0,vt=3,tc=1,tt=NULL, funcTag='bsp_mapping_workflow') {
    
    tabs <- paste0(rep(TAB, tc), collapse='')
    mssg <- glue::glue("[{funcTag}]:{tabs}")
    
    out_csv <- redata(out_dir, run_tag, funcTag, re_load, 
                      pre_tag, end_str=end_str, sep=sep_chr,
                      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (tibble::is_tibble(out_csv)) return(out_csv)
    if (is.null(out_csv)) {
      stop(glue::glue("{RET}{mssg} ERROR: out_csv is NULL!{RET2}"))
      return(out_csv)
    }
    out_dir <- base::dirname(out_csv)
    
    bsp_hit_csv <- out_csv %>% 
      stringr::str_remove(paste0(sep_chr,end_str,"$") ) %>%
      paste("bsp-alignment-histogram.csv.gz", sep=sep_chr)
    
    top_hit_csv <- out_csv %>% 
      stringr::str_remove(paste0(sep_chr,end_str,"$") ) %>%
      paste("bsp-multi-alignment-loci.csv.gz", sep=sep_chr)
    
    if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
    if (verbose>=vt+4) {
      cat(glue::glue("{RET}"))
      cat(glue::glue("{mssg} File IO Parameters::{RET}"))
      cat(glue::glue("{mssg}          ref_fas={ref_fas}.{RET}"))
      cat(glue::glue("{mssg}          can_fas={can_fas}.{RET}"))
      if (purrr::is_character(cgn_src))
        cat(glue::glue("{mssg}        cgn_src={cgn_src}.{RET}"))
      cat(glue::glue("{mssg}          out_dir={out_dir}.{RET}"))
      cat(glue::glue("{mssg}          out_csv={out_csv}.{RET}"))
      cat(glue::glue("{mssg}      bsp_hit_csv={bsp_hit_csv}.{RET}"))
      cat(glue::glue("{mssg}      top_hit_csv={top_hit_csv}.{RET}"))
      cat(glue::glue("{mssg}    canonical={canonical}.{RET}"))
      cat(glue::glue("{RET}"))
      
      cat(glue::glue("{mssg} BSMAP Parameters::{RET}"))
      cat(glue::glue("{mssg}          bsp_exe={bsp_exe}.{RET}"))
      cat(glue::glue("{mssg}          bsp_opt={bsp_opt}.{RET}"))
      cat(glue::glue("{RET}"))
      
      cat(glue::glue("{mssg} Field Parameters::{RET}"))
      cat(glue::glue("{mssg}      ids_key={ids_key}.{RET}"))
      cat(glue::glue("{mssg}      unq_key={unq_key}.{RET}"))
      cat(glue::glue("{mssg}      prb_key={prb_key}.{RET}"))
      cat(glue::glue("{mssg}      des_key={des_key}.{RET}"))
      cat(glue::glue("{mssg}      din_key={din_key}.{RET}"))
      cat(glue::glue("{mssg}     join_key={join_key}.{RET}"))
      cat(glue::glue("{mssg}    join_type={join_type}.{RET}"))
      cat(glue::glue("{RET}"))
      
      cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
      cat(glue::glue("{mssg}         sort={sort}.{RET}"))
      cat(glue::glue("{mssg}         full={full}.{RET}"))
      cat(glue::glue("{mssg}        merge={merge}.{RET}"))
      cat(glue::glue("{mssg}        light={light}.{RET}"))
      cat(glue::glue("{mssg}       reload={reload}.{RET}"))
      cat(glue::glue("{RET}"))
    }
    
    ret_cnt <- 0
    ret_tib <- NULL
    
    if (!is.null(ref_fas) && !is.null(ref_tib)) {
      cat(glue::glue("{RET}{mssg} Both ref_fas AND ref_tib can NOT be NULL!{RETs}"))
      return(ret_tib)
    }
    
    if (!is.null(can_fas) && !is.null(can_tib)) {
      cat(glue::glue("{RET}{mssg} Both can_fas AND can_tib can NOT be NULL!{RETs}"))
      return(ret_tib)
    }
    
    stime <- base::system.time({
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                    Pre-processing:: Write FASTA if Needed
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ids_sym <- rlang::sym(ids_key)
      prb_sym <- rlang::sym(prb_key)
      
      if (is.null(ref_fas) || !file.exists(ref_fas)) {
        ref_fas <- ref_tib %>%
          dplyr::filter( Alphabet=="dna" & Strand_BSC=="N" & Strand_FR=="F") %>% 
          head(n=1) %>% dplyr::pull(Path)
        
        if (is.null(ref_fas) || !file.exists(ref_fas)) {
          cat(glue::glue("{RET}{mssg} ERROR: Failed to find ref_fas from ",
                         "user input or reference genome tibble!{RETs}"))
          print(ref_tib)
          return(ret_tib)
        }
      }
      
      # Ensure we have Aln_Prb/Aln_Rev::
      can_tib <- can_tib %>%
        dplyr::mutate(Aln_Prb=deMs(!!prb_sym, uc=TRUE),
                      Aln_Rev=revCmp(Aln_Prb) )
      
      if (is.null(can_fas) || !file.exists(can_fas)) {
        fas_vec <- can_tib %>%
          dplyr::mutate( fas_line=paste0(">",!!ids_sym,"\n",Aln_Prb) ) %>%
          dplyr::pull(fas_line)
        
        can_fas <- out_csv %>% 
          stringr::str_remove(paste0(sep_chr,end_str,"$") ) %>%
          paste("fa.gz", sep=sep_chr)
        
        safe_write(x=fas_vec, type="line", file=can_fas, funcTag=funcTag,
                   verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
        
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                      Probe Alignment:: BSMAP
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ret_tib <-run_bsmap(ref_fas = ref_fas,
                          can_fas = can_fas,
                          out_dir=out_dir,
                          bsp_exe=bsp_exe,
                          sort=sort,
                          light=light, 
                          reload=reload,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                      Join Address and Alignment Data::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ret_tib <- join_bsmap(bsp_tib = ret_tib,
                            can_tib = can_tib,
                            bsp_csv = NULL,
                            cgn_src = cgn_src,
                            
                            unq_key = unq_key,
                            des_key = des_key,
                            din_key = din_key,
                            
                            join_key  = join_key,
                            join_type = join_type,
                            
                            sort = sort,
                            full = full,
                            
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                         Consolidate/Assign CGNs::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      #
      # TBD:: This function needs serious re-writing!!!!
      #
      ret_tib <- assign_cgn(ord_tib = can_tib,
                            bsp_tib = ret_tib,
                            seq_tib = seq_tib, 
                            can_csv = canonical,
                            
                            ids_key = join_key,
                            bsp_csv = out_csv,
                            merge   = merge,
                            
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                           Alignment Summary::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      bsp_hit_sum <- ret_tib %>% 
        dplyr::group_by(Address) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>% 
        dplyr::group_by(Count) %>% 
        dplyr::summarise(His_Count=n(), .groups="drop")
      hit_key <- glue::glue("bsp_hit_sum({funcTag})")
      hit_cnt <- print_tib(bsp_hit_sum, funcTag, verbose,vt+4,tc, n=hit_key)
      
      # Top Ranked Offfenders::
      top_hit_sum <- ret_tib %>% 
        dplyr::group_by(Address) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>%
        dplyr::filter(Count!=1) %>%
        dplyr::arrange(-Count)
      # print(top_hit_sum, n=base::nrow(top_hit_sum))
      
      sum_cnt <- safe_write(x=bsp_hit_sum, file=bsp_hit_csv, funcTag=funcTag,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      top_cnt <- safe_write(x=top_hit_sum, file=top_hit_csv, funcTag=funcTag,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      out_cnt <- safe_write(x=ret_tib, file=out_csv, funcTag=funcTag, done=TRUE,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      tt$addFile(out_csv)
      
      ret_key <- glue::glue("ret-FIN({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
      
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{mssg}{BRK}{RET2}"))
    
    ret_tib
  }

run_bsmap = function(ref_fas, 
                     can_fas, 
                     out_dir, 
                     
                     sort   = FALSE, 
                     light  = FALSE,
                     reload = FALSE,
                     
                     bsp_exe,
                     bsp_opt = "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R",
                     
                     verbose=0,vt=3,tc=1,tt=NULL,
                     funcTag='run_bsmap') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} improbe Parameters::{RET}"))
    cat(glue::glue("{mssg}   ref_fas={ref_fas}.{RET}"))
    cat(glue::glue("{mssg}   can_fas={can_fas}.{RET}"))
    cat(glue::glue("{mssg}   out_dir={out_dir}.{RET}"))
    cat(glue::glue("{mssg}   bsp_exe={bsp_exe}.{RET}"))
    cat(glue::glue("{mssg}   bsp_opt={bsp_opt}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("{mssg} Run Time Parameters::{RET}"))
    cat(glue::glue("{mssg}      sort={sort}.{RET}"))
    cat(glue::glue("{mssg}     light={light}.{RET}"))
    cat(glue::glue("{mssg}    reload={reload}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
    
    out_bsp <- file.path(out_dir, "probe.bsp")
    bsp_ssh <- file.path(out_dir, "probe.bsp.run.sh")
    bsp_tsv <- file.path(out_dir, "probe.bsp.tsv.gz")
    
    if (!reload || !file.exists(bsp_tsv)) {
      add_cmd <- ""
      if (light) add_cmd <- "cut -f 1,2,4-11 | "
      
      bsp_cmd <- glue::glue(
        "{bsp_exe} -a {can_fas} -d {ref_fas} {bsp_opt} -o {out_bsp}{RET}",
        "cat {out_bsp} | {add_cmd} gzip -c -> {bsp_tsv}{RET}",
        "rm -f {out_bsp}{RET}")
      
      out_cnt <- safe_write(bsp_cmd, "line", bsp_ssh, funcTag = funcTag,
                            verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
      Sys.chmod(paths=bsp_ssh, mode="0777")
      
      if (verbose>=vt)
        cat(glue::glue("{mssg} CMD={bsp_cmd}...{RET}"))
      
      sys_ret <- 0
      sys_ret <- base::system(bsp_ssh)
      
      if (sys_ret!=0) {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: sys_ret({sys_ret}) != 0!!!{RET2}"))
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
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  ret_tib
}

load_bsmap = function(file,
                      sort=FALSE,
                      light=FALSE,
                      
                      verbose=0,vt=3,tc=1,tt=NULL,
                      funcTag='load_bsmap') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}    file={file}{RET}"))
    cat(glue::glue("{mssg}    sort={sort}{RET}"))
    cat(glue::glue("{mssg}   light={light}{RET}"))
    cat(glue::glue("{RET}"))
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
      cat(glue::glue("{mssg} Loading BSP={file}...{RET}"))
    
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
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  ret_tib
}

join_bsmap = function(bsp_tib,
                      can_tib, 
                      bsp_csv = NULL,
                      cgn_src = NULL,
                      
                      unq_key = "Unq_Key",
                      des_key = "Ord_Des",
                      din_key = "Ord_Din",
                      
                      join_key,
                      join_type = "inner", 
                      
                      sort = TRUE,
                      full = FALSE,
                      
                      verbose=0,vt=3,tc=1,tt=NULL,
                      funcTag='join_bsmap') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+4) {
    
    cat(glue::glue("{mssg} Field Parameters::{RET}"))
    cat(glue::glue("{mssg}     des_key={des_key}{RET}"))
    cat(glue::glue("{mssg}     din_key={din_key}{RET}"))
    cat(glue::glue("{mssg}    join_key={join_key}{RET}"))
    cat(glue::glue("{mssg}   join_type={join_type}{RET2}"))
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
    
    aln_unq_sym  <- rlang::sym(unq_key)
    
    prb_des_sym  <- rlang::sym(des_key)
    prb_din_sym  <- rlang::sym(din_key)
    
    ord_join_key <- join_key
    ord_join_sym <- rlang::sym(ord_join_key)
    
    bsp_join_key <- names(bsp_tib)[1]
    bsp_join_sym <- rlang::sym( bsp_join_key )
    
    if (verbose>=vt+4) {
      cat(glue::glue("{mssg}    join_type={join_type}{RET}"))
      cat(glue::glue("{mssg}      des_key={des_key}{RET}"))
      cat(glue::glue("{mssg}      din_key={din_key}{RET2}"))
      cat(glue::glue("{mssg} bsp_join_key={bsp_join_key}{RET}"))
      cat(glue::glue("{mssg} ord_join_key={ord_join_key}{RET}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Join BSP Alignment with Probe Data::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    bsp_key <- glue::glue("bsp-tib({funcTag})")
    bsp_cnt <- print_tib(bsp_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    
    ord_key <- glue::glue("ord-tib({funcTag})")
    ord_cnt <- print_tib(can_tib,funcTag, verbose,vt+4,tc, n=ord_key)
    
    if (join_type=="inner") {
      ret_tib <- bsp_tib %>%
        dplyr::rename(!!ord_join_sym := !!bsp_join_sym) %>%
        dplyr::inner_join(can_tib, by=ord_join_key) # %>% head(n=100)
    } else if (join_type=="left") {
      ret_tib <- bsp_tib %>%
        dplyr::rename(!!ord_join_sym := !!bsp_join_sym) %>%
        dplyr::left_join(can_tib, by=ord_join_key) # %>% head(n=100)
    } else if (join_type=="right") {
      ret_tib <- bsp_tib %>%
        dplyr::rename(!!ord_join_sym := !!bsp_join_sym) %>%
        dplyr::right_join(can_tib, by=ord_join_key) # %>% head(n=100)
    } else if (join_type=="full") {
      ret_tib <- bsp_tib %>%
        dplyr::rename(!!ord_join_sym := !!bsp_join_sym) %>%
        dplyr::full_join(can_tib, by=ord_join_key) # %>% head(n=100)
    } else {
      stop(glue::glue("{mssg} ERROR: Unsupported join_type={join_type}!!!{RET2}"))
      return(ret_tib)
    }
    ret_key <- glue::glue("bsp/ord-join({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      Calculate New Fields from Join::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt)
      cat(glue::glue("{mssg} Calculating new fields from joined data...{RET}"))
    
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
        !!aln_unq_sym := paste(!!ord_join_sym, Bsp_Rank, sep="_")
        # Aln_Key_Unq=paste(!!ord_join_sym, Bsp_Rank, sep="_")
      ) %>%
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
      cat(glue::glue("{mssg} Done. Calculating new fields.{RET2}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      Annotate CGN by Alignment Position::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (!is.null(cgn_src) && tibble::is_tibble(cgn_src)) {
      if (verbose>=vt+4) {
        cat(glue::glue("{mssg} Usung tibble={RET}"))
        print(cgn_src)
      }
    } else if (!is.null(cgn_src) && 
               file.exists(cgn_src) && 
               !dir.exists(cgn_src)) {
      if (verbose>=vt+4)
        cat(glue::glue("{mssg} Loading from file={cgn_src}.{RET}"))
      cgn_src <- load_cgn_dB(file=cgn_src, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
    } else {
      cat(glue::glue("{RET}{mssg} Warning: Unknown cgn_src type!{RET}"))
      cgn_src <- NULL
    }
    
    if (!is.null(cgn_src)) {
      if (verbose>=vt)
        cat(glue::glue("{mssg} Adding cgn by alignment position...{RET}"))
      
      cgn_src <- cgn_src %>% purrr::set_names(c("Bsp_Chr","Bsp_Pos","Bsp_Cgn","Cgd_Top"))
      
      mat_tib <- dplyr::bind_rows(
        dplyr::inner_join(
          # ret_tib %>% dplyr::select(Aln_Key_Unq,Bsp_Chr,Bsp_Pos),
          ret_tib %>% dplyr::select(!!aln_unq_sym,Bsp_Chr,Bsp_Pos),
          cgn_src %>% dplyr::mutate(Cgd_Nuc="up"),
          by=c("Bsp_Chr","Bsp_Pos"),
          suffix=c("_bsp","_cgn")),
        dplyr::inner_join(
          # ret_tib %>% dplyr::select(Aln_Key_Unq,Bsp_Chr,Bsp_Pos),
          ret_tib %>% dplyr::select(!!aln_unq_sym,Bsp_Chr,Bsp_Pos),
          cgn_src %>% dplyr::mutate(Bsp_Pos=Bsp_Pos+1, Cgd_Nuc="dn"),
          by=c("Bsp_Chr","Bsp_Pos"),
          suffix=c("_bsp","_cgn"))
      ) %>% 
        dplyr::distinct(!!aln_unq_sym, .keep_all = TRUE)
      # dplyr::distinct(Aln_Key_Unq, .keep_all = TRUE)
      
      ret_tib <- ret_tib %>% 
        # dplyr::left_join(mat_tib, by=c("Aln_Key_Unq","Bsp_Chr","Bsp_Pos")) %>%
        dplyr::left_join(mat_tib, by=c(unq_key,"Bsp_Chr","Bsp_Pos")) %>%
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
        cat(glue::glue("{mssg} Done. Adding cgn by position.{RET2}"))
    }
    
    ret_tib <- ret_tib %>% 
      dplyr::select(dplyr::any_of(imp_col_vec),
                    dplyr::everything()) %>% clean_tibble()
    
    out_cnt <- safe_write(x=ret_tib, file=bsp_csv, funcTag=funcTag,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Assign Best CGN from:: BSP & SEQ
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: This function needs to be cleaned up! Logic probably needs to be fixed
#   as well. 
#
assign_cgn = function(ord_tib,
                      bsp_tib,
                      seq_tib,
                      
                      can_csv,
                      can_tib = NULL,
                      bsp_csv = NULL,
                      
                      ids_key = ids_key,  # Use to be Aln_Key
                      
                      join    = "inner",
                      merge   = TRUE, 
                      retData = FALSE, 
                      
                      verbose=0,vt=3,tc=1,tt=NULL,
                      funcTag='assign_cgn') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}    can_csv={can_csv}.{RET}"))
    cat(glue::glue("{mssg}    bsp_csv={bsp_csv}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}      ids_key={ids_key}.{RET}"))
    cat(glue::glue("{mssg}      join={join}.{RET}"))
    cat(glue::glue("{mssg}     merge={merge}.{RET}"))
    cat(glue::glue("{mssg}   retData={retData}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ids_sym <- rlang::sym(run$ids_key)
    
    #
    # Load Canonical CGNs::
    #
    if (is.null(can_tib)) {
      can_tib <- safe_read(file=run$canonical_csv, verbose=opt$verbose) %>% 
        dplyr::select(CGN) %>% 
        dplyr::rename(Cgn=CGN) %>% 
        dplyr::mutate(Can_Cnt=1)
    }
    
    bsp_key <- glue::glue("can_tib({funcTag})")
    bsp_cnt <- print_tib(can_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    if (verbose>=vt)
      cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    if (retData) ret_dat$can_tib <- can_tib
    
    # Defined Order tib to ord original cgn::
    ord_tib <- aqp_ord_tib %>% 
      dplyr::select(Ord_Map,Prb_Key,Ord_Cgn) %>%
      dplyr::rename(Cgn=Ord_Cgn) %>% 
      dplyr::mutate(Ord_Cnt=1) %>%
      dplyr::distinct() %>%
      clean_tibble()
    
    bsp_key <- glue::glue("ord_tib({funcTag})")
    bsp_cnt <- print_tib(ord_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    if (verbose>=vt)
      cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    if (retData) ret_dat$ord_tib <- ord_tib
    
    #
    # Format BSP::
    #
    bsp_tib <- aqp_bsp_tib %>% 
      dplyr::filter(!is.na(Bsp_Cgn)) %>% 
      dplyr::select(Ord_Map, Prb_Key, Ord_Des, Ord_Din, Bsp_Cgn) %>% 
      dplyr::rename(Cgn=Bsp_Cgn) %>%
      dplyr::distinct() %>%
      dplyr::arrange(Prb_Key, Cgn) %>% 
      dplyr::group_by(Ord_Map,Prb_Key,Ord_Des,Ord_Din,Cgn) %>% 
      dplyr::summarise(Bsp_Cnt=n(), .groups = "drop") %>%
      clean_tibble()
    
    bsp_key <- glue::glue("bsp_tib({funcTag})")
    bsp_cnt <- print_tib(bsp_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    if (verbose>=vt)
      cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    #
    # Format Seq::
    #
    seq_tib <- aqp_seq_tib %>% 
      dplyr::filter(!is.na(Imp_Cgn)) %>% 
      dplyr::left_join(dplyr::distinct(aqp_bsp_tib,Address, Ord_Map), by=c("Address")) %>%
      dplyr::select(Address, Ord_Des, Ord_Din, Ord_Map, Imp_Cgn) %>% 
      tidyr::unite(Tmp_Key, Ord_Des,Ord_Din, sep="", remove=FALSE) %>%
      tidyr::unite(Prb_Key, Address, Tmp_Key, sep="_", remove=FALSE) %>%
      dplyr::select(-Tmp_Key) %>%
      dplyr::select(Ord_Map,Prb_Key,Ord_Des,Ord_Din,Imp_Cgn) %>%
      dplyr::rename(Cgn=Imp_Cgn) %>%
      dplyr::distinct() %>%
      dplyr::arrange(Prb_Key, Cgn) %>% 
      dplyr::group_by(Ord_Map,Prb_Key,Ord_Des,Ord_Din,Cgn) %>% 
      dplyr::summarise(Seq_Cnt=n(), .groups = "drop") %>%
      clean_tibble()
    
    bsp_key <- glue::glue("seq_tib({funcTag})")
    bsp_cnt <- print_tib(seq_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    if (verbose>=vt)
      cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    #
    # Critical Step:: Build and Sort Counts Tables
    #  NOTE:: Should probably save this...
    #
    cnt_tib <- 
      dplyr::full_join(
        bsp_tib, seq_tib, by=c("Ord_Map",run$ids_key,"Ord_Des","Ord_Din","Cgn")
      ) %>%
      dplyr::left_join(can_tib, by="Cgn") %>%
      dplyr::distinct() %>%
      dplyr::left_join(ord_tib, by=c(run$ids_key,"Cgn","Ord_Map")
      ) %>%
      dplyr::distinct() %>%
      dplyr::mutate(dplyr::across(c(Bsp_Cnt,Seq_Cnt,Can_Cnt,Ord_Cnt), tidyr::replace_na, 0 ),
                    Sum_Cnt=Bsp_Cnt+Seq_Cnt,
                    Max_Cnt=Bsp_Cnt*Seq_Cnt) %>% 
      dplyr::add_count(!!ids_sym, name="Cgn_Cnt") %>% 
      dplyr::arrange(-Can_Cnt,-Max_Cnt,-Sum_Cnt,-Ord_Cnt) %>%
      dplyr::mutate(Rank=dplyr::row_number()) %>%
      clean_tibble()
    
    bsp_key <- glue::glue("cnt_tib({funcTag})")
    bsp_cnt <- print_tib(cnt_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    tib_cnt <- cnt_tib %>% base::nrow()
    if (verbose>=vt)
      cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    cnt_list <- cnt_tib %>% split(.$Ord_Des)
    #
    # Infinium II::
    #
    inf2_tib <- cnt_list[["2"]] %>%
      dplyr::arrange(!!ids_sym, Rank) %>%
      dplyr::distinct(!!ids_sym, .keep_all = TRUE)
    
    bsp_key <- glue::glue("inf2_tib({funcTag})")
    bsp_cnt <- print_tib(inf2_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    if (verbose>=vt)
      cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    # Infinium I:: Full Join
    #
    #   TBD:: The joining should really be done by sequence: Ord_Prb
    #
    run$join <- "inner"
    if (run$join=="full") {
      inf1_tib <- dplyr::full_join(
        cnt_list[["U"]], 
        cnt_list[["M"]], 
        by=c("Ord_Map","Cgn","Ord_Din"), 
        suffix=c("_U","_M")
      ) %>%
        dplyr::mutate(Rank_Min=pmin(Rank_U,Rank_M)) %>%
        dplyr::arrange(Ord_Map, Rank_Min) %>%
        dplyr::distinct(Ord_Map,Prb_Key_U,Prb_Key_M, .keep_all = TRUE)
    } else if (run$join=="inner") {
      inf1_tib <- dplyr::inner_join(
        cnt_list[["U"]], 
        cnt_list[["M"]], 
        by=c("Ord_Map","Cgn","Ord_Din"), 
        suffix=c("_U","_M")
      ) %>%
        dplyr::mutate(Rank_Min=pmin(Rank_U,Rank_M)) %>%
        dplyr::arrange(Ord_Map, Rank_Min) %>%
        dplyr::distinct(Ord_Map,Prb_Key_U,Prb_Key_M, .keep_all = TRUE)
    } else {
      stop(glue::glue("{mssg} Unsupported join type={join}.{RET}"))
      return(NULL)
    }
    bsp_key <- glue::glue("inf1_tib({funcTag})")
    bsp_cnt <- print_tib(inf1_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    if (retData) ret_dat$inf1_tib <- inf1_tib
    if (retData) ret_dat$inf2_tib <- inf2_tib
    
    ret_tib <- dplyr::bind_rows(
      dplyr::select(inf1_tib, Ord_Map,Prb_Key_U,Cgn,Ord_Des_U,Ord_Din,Can_Cnt_U,Rank_Min) %>% 
        purrr::set_names("Ord_Map",run$ids_key,"Cgn","Ord_Des","Ord_Din","Can_Cnt","Rank"),
      
      dplyr::select(inf1_tib, Ord_Map,Prb_Key_M,Cgn,Ord_Des_M,Ord_Din,Can_Cnt_M,Rank_Min) %>% 
        purrr::set_names("Ord_Map",run$ids_key,"Cgn","Ord_Des","Ord_Din","Can_Cnt","Rank"),
      
      dplyr::select(inf2_tib, Ord_Map,!!ids_sym,Cgn,Ord_Des,Ord_Din,Can_Cnt,Rank)
    ) %>% dplyr::filter(!is.na(!!ids_sym)) %>%
      dplyr::distinct()
    bsp_key <- glue::glue("ret_tib({funcTag})")
    bsp_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    if (verbose>=vt)
      cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    
    mul_cnt <- ret_tib %>% dplyr::add_count(!!ids_sym,Cgn, name="Multi_Cnt") %>% 
      dplyr::filter(Multi_Cnt != 1) %>% base::nrow()
    mis_cnt <- ret_tib %>% dplyr::filter(is.na(!!ids_sym)) %>% base::nrow()
    
    mis_tib <- dplyr::anti_join(ord_tib, ret_tib, by=c(run$ids_key)) %>%
      dplyr::left_join(aqp_ord_tib %>% dplyr::select(Prb_Key,Ord_Des,Ord_Din), by=c("Prb_Key"))
    sig_tib <- dplyr::filter(cnt_tib, run$ids_key %in% mis_tib[[run$ids_key]]) %>%
      dplyr::arrange(Ord_Map,Rank) %>%
      dplyr::distinct(!!ids_sym, .keep_all = TRUE)
    sig_cnt <- sig_tib %>% base::nrow()
    
    bsp_key <- glue::glue("ret_tib({funcTag})")
    bsp_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    if (verbose>=vt)
      cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    if (verbose>=vt) {
      cat(glue::glue("{mssg}   Miss Count={mis_cnt}.{RET}"))
      cat(glue::glue("{mssg}  Multi Count={mul_cnt}.{RET}"))
      cat(glue::glue("{mssg} Single Count={sig_cnt}.{RET2}"))
    }
    
    # if (mis_cnt!=0 || mul_cnt!=0 || sig_cnt!=0) {
    #   stop(glue::glue("{RET}{mssg} Counts non-zero={mis_cnt},",
    #                   "{mul_cnt},{sig_cnt}!{RET2}"))
    #   return(NULL)
    # }
    
    # Merge all data together::
    #
    ret_cnt <- ret_tib %>% base::nrow()
    ret_tib <- dplyr::bind_rows(
      #
      # Add Formatted cg#'s
      #
      ret_tib %>% 
        dplyr::mutate(
          
          Cgn_Tag=dplyr::case_when(
            Ord_Din=="rs" ~ Ord_Din,
            Ord_Din=="ch" ~ Ord_Din,
            TRUE ~ "cg"
          ),
          Cgn_Str=dplyr::case_when(
            Ord_Din=="rs" ~ stringr::str_remove(Ord_Map, "[-_:].*$"),
            Ord_Din=="ch" ~ stringr::str_remove(Ord_Map, "[-_:].*$"),
            TRUE ~ paste0("cg",stringr::str_pad(Cgn,width=8,side="left",pad="0"))
          )
        ),
      mis_tib %>% 
        dplyr::select(Ord_Map, !!ids_sym,Cgn,Ord_Des,Ord_Din) %>% 
        
        # Ord_Cgn -> Cgn naming is done above
        # dplyr::select(Ord_Map, !!ids_sym,Ord_Cgn,Ord_Des,Ord_Din) %>% 
        # dplyr::rename(Cgn=Ord_Cgn) %>% 
        dplyr::mutate(Can_Cnt=0, 
                      Rank=dplyr::row_number() + ret_cnt,
                      Cgn_Tag="uk",
                      Cgn_Str=paste0(Cgn_Tag,stringr::str_pad(Cgn,width=8,side="left",pad="0")
                      )
        )
    ) %>%
      # TBD:: Capture other CGN's in seperate column:: actual CGN's not Count!!
      dplyr::add_count(!!ids_sym, name="Alt_Cgn_Cnt") %>%
      
      # One Final Clean Up To Ensure Uniqueness::
      dplyr::arrange(Rank) %>% 
      dplyr::distinct(!!ids_sym, .keep_all = TRUE)
    
    mul_cnt <- ret_tib %>% 
      dplyr::add_count(!!ids_sym,Cgn, name="Multi_Cnt") %>% 
      dplyr::filter(Multi_Cnt != 1) %>% base::nrow()
    
    if (verbose>=vt)
      cat(glue::glue("{mssg}  Multi Count Final={mul_cnt}.{RET2}"))
    
    if (mul_cnt!=0) {
      stop(glue::glue("{RET}{mssg} Multi-Count Final={mul_cnt} ",
                      "not equal to zero!!!{RET2}"))
      return(NULL)
    }
    
    if (merge) ret_tib <- bsp_tib %>%
      dplyr::left_join(ret_tib, 
                       by=c("Ord_Map",ids_key,"Ord_Des","Ord_Din"),
                       suffix=c("_bsp","_cgn"))
    
    ret_tib <- ret_tib %>% clean_tibble()
    out_cnt <- safe_write(ret_tib,file=bsp_csv, funcTag=funcTag,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Load CGN/Position Database::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_cgn_dB = function(file,
                       verbose=0,vt=3,tc=1,tt=NULL,
                       funcTag='load_cgn_dB') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   file={file}.{RET}"))
    cat(glue::glue("{RET}"))
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
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# End of file
