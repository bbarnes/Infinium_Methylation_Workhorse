
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
#                         CGN Mapping Workflow Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cgn_mapping_workflow = function(tib, dir,
                                pattern_u,
                                pattern_m,
                                
                                prb_key = "Prb_Seq",
                                add_key = "Address", 
                                des_key = "Prb_Des", 
                                din_key = "Prb_Din",
                                ids_key = "Prb_Key",
                                aln_key = "Aln_P49",

                                out    = NULL,
                                prefix = NULL,
                                suffix = NULL,
                                
                                idxA=1, idxB=1,
                                reload   = FALSE,
                                parallel = FALSE,
                                
                                del="_",
                                
                                out_csv=NULL, out_dir, run_tag, 
                                re_load=FALSE, pre_tag=NULL,
                                end_str='csv.gz', sep_chr='.',
                                verbose=0,vt=3,tc=1,tt=NULL,
                                funcTag='cgn_mapping_workflow') {
  
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
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}         dir={dir}.{RET}"))
    cat(glue::glue("{mssg}   pattern_u={pattern_u}.{RET}"))
    cat(glue::glue("{mssg}   pattern_m={pattern_m}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}     prb_key={prb_key}.{RET}"))
    cat(glue::glue("{mssg}     add_key={add_key}.{RET}"))
    cat(glue::glue("{mssg}     des_key={des_key}.{RET}"))
    cat(glue::glue("{mssg}     din_key={din_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}     out_dir={out_dir}.{RET}"))
    cat(glue::glue("{mssg}     out_csv={out_csv}.{RET}"))
    cat(glue::glue("{mssg}      prefix={prefix}.{RET}"))
    cat(glue::glue("{mssg}      suffix={suffix}.{RET}"))
    cat(glue::glue("{mssg}        idxA={idxA}.{RET}"))
    cat(glue::glue("{mssg}        idxB={idxB}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}      reload={reload}.{RET}"))
    cat(glue::glue("{mssg}    parallel={parallel}.{RET}"))
    cat(glue::glue("{mssg}         del={del}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ids_sym <- rlang::sym(ids_key)
    prb_sym <- rlang::sym(prb_key)
    des_sym <- rlang::sym(des_key)
    din_sym <- rlang::sym(din_key)
    aln_sym <- rlang::sym(aln_key)
    
    if (!dir.exists(out)) dir.create(out, recursive = TRUE)

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 Build Candidate Sub-String Probes:: U49/M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # NOTE: Build U49/M49 Sub-sequences and split by 5' two nucelotide prefix  
    #  into sorted output files. Splitting by prefix makes the join later much 
    #  faster...
    #
    
    ret_tib <- tib %>%
      dplyr::mutate(
        Aln_Prb = deMs(!!prb_sym, uc=TRUE),
        !!aln_sym := dplyr::case_when(
          !!des_sym == '2' ~ stringr::str_sub(Aln_Prb, 2),
          !!des_sym == 'U' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          !!des_sym == 'M' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::distinct(!!ids_sym,Aln_Prb, .keep_all=TRUE) %>%
      clean_tibble()
    
    ret_key <- glue::glue("add-P49({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Intersecting by Probe Sequence:: U49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    u49_tib <- 
      intersect_seq_workflow(tib = ret_tib, dir = dir, pattern = pattern_u,
                             
                             tar_bsc = "U49", tar_des = c("2","U"), 
                             
                             prb_key = prb_key, 
                             des_key = des_key,
                             din_key = din_key,
                             ids_key = ids_key,
                             aln_key = aln_key,

                             idxA = idxA, idxB = idxB, 
                             out = out_dir, prefix = prefix, suffix = suffix, 
                             
                             reload   = reload, parallel = parallel, 
                             
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Intersecting by Probe Sequence:: M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    m49_tib <- 
      intersect_seq_workflow(tib = ret_tib, dir = dir, pattern = pattern_m,
                             
                             tar_bsc = "M49", tar_des = c("M"), 
                             
                             prb_key = prb_key, 
                             des_key = des_key,
                             din_key = din_key,
                             ids_key = ids_key,
                             aln_key = aln_key,
                             
                             idxA = idxA, idxB = idxB, 
                             out = out_dir, prefix = prefix, suffix = suffix, 
                             
                             reload   = reload, parallel = parallel, 
                             
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Intersecting by Probe Sequence:: M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- join_seq_intersect(u49 = u49_tib, 
                                  m49 = m49_tib, 
                                  ids_key = ids_key,
                                  verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    # Defined Above::
    #
    # out_csv <- 
    #   file.path(out, paste(prefix,suffix,"intersect.tsv.gz", sep='.'))
    
    out_cnt <- safe_write(ret_tib, file=out_csv, funcTag=funcTag, done=TRUE,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    tt$addFile(out_csv)
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET}{RET}"))
  
  ret_tib
}

intersect_seq_workflow = function(tib, dir, pattern, tar_bsc, tar_des,

                                  prb_key = "Prb_Seq", 
                                  des_key = "Prb_Des", 
                                  din_key = "Prb_Din",
                                  ids_key = "Prb_Key",
                                  aln_key = "Aln_P49",

                                  out, prefix, suffix,
                                  
                                  idxA=1,idxB=1,
                                  
                                  reload=FALSE, parallel=FALSE,
                                  
                                  verbose=0,vt=3,tc=1,tt=NULL,
                                  funcTag='intersect_seq_workflow') {
  
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
    
    prb_sym <- rlang::sym(prb_key)
    des_sym <- rlang::sym(des_key)
    din_sym <- rlang::sym(din_key)
    aln_sym <- rlang::sym(aln_key)

    tar_des_str <- paste(tar_des, collapse="-")
    
    if (!dir.exists(out)) dir.create(out, recursive = TRUE)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #      Find All Pre-Built Reference Prefix-Partition Files:: U49/M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    cgn_tsvs <- get_file_list(dir = dir, 
                              pattern = pattern, 
                              trim    = pattern,
                              verbose = opt$verbose)
    
    tsv_cnt <- cgn_tsvs %>% names() %>% length()
    pre_len <- cgn_tsvs %>% names() %>% stringr::str_length() %>% max()
    
    if (verbose>=vt+6) {
      cat(glue::glue("{mssg} Found {tsv_cnt} {tar_bsc}-files, ",
                     "prefix length={pre_len}, tar_des={tar_des_str}.{RET}"))
      cgn_tsvs %>% head(n=3) %>% print()
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Intersecting by Probe Sequence:: U49/M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    prb_tibs <- tib %>% 
      dplyr::filter(!!des_sym %in% tar_des) %>%
      dplyr::filter(!is.na(!!aln_sym)) %>%
      dplyr::select(dplyr::all_of(c(aln_key,ids_key)) ) %>%
      dplyr::mutate(
        Pre_Nuc=stringr::str_sub(!!aln_sym, 1, pre_len)
      ) %>%
      dplyr::arrange(!!aln_sym) %>%
      split(f=.$Pre_Nuc, drop = TRUE)

    if (verbose>=vt+6) {
      nuc_cnt <- prb_tibs %>% names() %>% length()
      cat(glue::glue("{mssg} Found {nuc_cnt} {tar_bsc}-tibs, ",
                     "prefix length={pre_len}, tar_des={tar_des_str}.{RET}"))
      ret_key <- glue::glue("{tar_bsc}/{tar_des_str}-{tar_des}-pre1-int-tib({funcTag})")
      ret_cnt <- print_tib(prb_tibs[[2]],funcTag, verbose,vt+6,tc, n=ret_key)
      ret_key <- glue::glue("{tar_bsc}/{tar_des_str}-{tar_des}-pre2-int-tib({funcTag})")
      ret_cnt <- print_tib(prb_tibs[[2]],funcTag, verbose,vt+6,tc, n=ret_key)
    }
    
    ret_tib <- NULL
    if (parallel) {
      if (verbose>=vt) 
        cat(glue::glue("{mssg}{TAB} Intersecting probe sequences ",
                       "{tar_bsc}-{tar_des_str} (Parallel)...{RET}"))
      
      ret_tib <- foreach (pre_nuc=names(prb_tibs), .combine=rbind) %dopar% {
        cur_tib <- prb_tibs[[pre_nuc]] %>% dplyr::select(-Pre_Nuc)
        intersect_seq_strand(
          can = cur_tib, ref = cgn_tsvs[[pre_nuc]],
          bsc_str = tar_bsc, pre_nuc = pre_nuc, idxA = idxA, idxB = idxB,
          out = out, prefix = prefix, suffix = suffix, reload = reload,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      }
    } else {
      if (verbose>=vt) 
        cat(glue::glue("{mssg}{TAB} Intersecting probe sequences ",
                       "{tar_bsc}-{tar_des_str} (Linear)...{RET}"))
      
      for (pre_nuc in names(prb_tibs)) {
        cur_tib <- prb_tibs[[pre_nuc]] %>% dplyr::select(-Pre_Nuc)
        
        prb_tib <- intersect_seq_strand(
          can = cur_tib, ref = cgn_tsvs[[pre_nuc]],
          bsc_str = tar_bsc, pre_nuc = pre_nuc, idxA = idxA, idxB = idxB,
          out = out, prefix = prefix, suffix = suffix, reload = reload,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        ret_tib <- dplyr::bind_rows(ret_tib, prb_tib)
        
        if (verbose>=vt) 
          cat(glue::glue("{mssg}{TAB}{TAB} Done. Intersecting probe ",
                         "sequences {tar_bsc}/{tar_des_str} nuc={pre_nuc}{RET2}"))
      }
    }
    if (verbose>=vt) 
      cat(glue::glue("{mssg}{TAB} Done. Intersecting probe sequences ",
                     "{tar_bsc}/{tar_des_str}!{RET2}"))
    
    ret_key <- glue::glue("{tar_bsc}/{tar_des_str}-{tar_des}-int-tib({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET}{RET}"))
  
  ret_tib
}

intersect_seq_strand = function(can, ref,
                                bsc_str, pre_nuc, idxA, idxB,
                                out, prefix, suffix, 
                                reload = FALSE,
                                verbose=0,vt=3,tc=1,tt=NULL,
                                funcTag='intersect_seq_strand') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")

  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}       ref={ref}.{RET}"))
    cat(glue::glue("{mssg}   bsc_str={bsc_str}.{RET}"))
    cat(glue::glue("{mssg}   pre_nuc={pre_nuc}.{RET}"))
    cat(glue::glue("{mssg}      idxA={idxA}.{RET}"))
    cat(glue::glue("{mssg}      idxB={idxB}.{RET}"))
    cat(glue::glue("{mssg}       out={out}.{RET}"))
    cat(glue::glue("{mssg}    prefix={prefix}.{RET}"))
    cat(glue::glue("{mssg}    suffix={suffix}.{RET}"))
    cat(glue::glue("{mssg}    reload={reload}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ref_tsv <- ref
    can_tsv <- file.path(
      out, paste(prefix, suffix,pre_nuc, bsc_str,"tsv.gz", sep='.') )
    out_tsv <- file.path(
      out, paste(prefix, suffix,pre_nuc, bsc_str, "intersect.tsv.gz", sep='.') )
    
    safe_write(x = can, file = can_tsv, cols = FALSE, 
               funcTag = funcTag, verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    ret_tib <- intersect_seq(ref = ref_tsv,
                             can = can_tsv,
                             out = out_tsv,

                             idxA = idxA,
                             idxB = idxB,
                             reload = reload,

                             verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("{pre_nuc}-{bsc_str}-int-tib({funcTag})")
    ret_cnt <- print_tib(ret_tib, funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET}{RET}"))
  
  ret_tib
}

intersect_seq = function(ref, can, out, 
                         idxA=1, idxB=1, 
                         reload=FALSE,
                         
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='intersect_seq') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("{mssg}  ref={ref}{RET}"))
    cat(glue::glue("{mssg}  can={can}{RET}"))
    cat(glue::glue("{mssg}  out={out}{RET}"))
    cat(glue::glue("{mssg} idxA={idxA}{RET}"))
    cat(glue::glue("{mssg} idxB={idxB}{RET}"))
    cat(glue::glue("{mssg} reload={reload}{RET}{RET}"))
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
      Prb_Key  = col_character()
    )
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    if (reload && file.exists(out)) {
      if (verbose>=vt+1) cat(glue::glue("{mssg} Reloading!{RET}"))
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
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET}{RET}"))
  
  ret_tib
}

join_seq_intersect = function(u49, m49,
                              ids_key = "Prb_Key",
                              
                              verbose=0,vt=3,tc=1,tt=NULL,
                              funcTag='join_seq_intersect') {

  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
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
    ret_key <- glue::glue("ret_tib1({funcTag})")
    ret_cnt <- print_tib(ret_tib, funcTag, verbose,vt+6,tc, n=ret_key)
    
    ret_tib <- ret_tib %>%
      tidyr::separate(Imp_Key, 
                      into=c("Imp_Cgn",
                             "Imp_Hit_hg38", "Imp_Hit_hg37", 
                             "Imp_Hit_hg36", "Imp_Hit_mm10"), 
                      sep="_", remove=TRUE)
    ret_key <- glue::glue("ret_tib2({funcTag})")
    ret_cnt <- print_tib(ret_tib, funcTag, verbose,vt+6,tc, n=ret_key)
    
    ret_tib <- ret_tib %>%
      tidyr::separate(ids_key, into=c("Address", "Tmp_Key"), 
                      sep="_", remove=FALSE) %>%
      tidyr::separate(Tmp_Key, into=c("Ord_Des","Ord_Din"), sep=1) %>%
      dplyr::rename(Aln_Prb=Imp_Seq, Aln_Nuc=Imp_Nuc)
    ret_key <- glue::glue("ret_tib3({funcTag})")
    ret_cnt <- print_tib(ret_tib, funcTag, verbose,vt+6,tc, n=ret_key)
    
    ret_tib <- ret_tib %>%
      tidyr::separate(Imp_Srd3, into=c("Imp_TB","Imp_CO", "Imp_Nxb"),
                      sep=c(1,2)) %>%
      dplyr::select(Address, Ord_Des, Ord_Din, Imp_Cgn, 
                    Imp_TB, Imp_CO, Imp_Nxb, Aln_Prb, Aln_Nuc, 
                    dplyr::everything()) %>%
      clean_tibble()
    ret_key <- glue::glue("ret_tib4({funcTag})")
    ret_cnt <- print_tib(ret_tib, funcTag, verbose,vt+6,tc, n=ret_key)
    
    ret_tib <- ret_tib  %>%
      dplyr::select(dplyr::any_of(imp_col_vec), dplyr::everything()) %>%
      clean_tibble()
    
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET}{RET}"))
  
  ret_tib
}

# End of file
