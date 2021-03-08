
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
#                       Address To Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

addressToManifest = function(tib, des="Ord_Des", ord_id="Ord_Key",
                             join,
                             csv=NULL, validate=TRUE,
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'addressToManifest'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (!is.null(csv)) {
      outDir <- base::dirname(csv)
      if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
    }
    
    des_list <- NULL
    des_list <- tib %>% split(.[[des]])
    des_cnt  <- des_list %>% names() %>% length()
    
    # Build Infinium I::
    ret1_tib <- NULL
    if (!is.null(des_list[["U"]]) && !is.null(des_list[["M"]])) {
      ret1_tib <- dplyr::full_join(
        dplyr::select(des_list[["U"]], -dplyr::all_of(des)), 
        dplyr::select(des_list[["M"]], -dplyr::all_of(des)), 
        by=dplyr::all_of(join),
        suffix=c("_U","_M")
      )
    }
    ret1_cnt <- print_tib(ret1_tib,funcTag, verbose,vt+4,tc, n="InfI")
    
    # Build Infinium II::
    ret2_tib <- NULL
    ret2_tib <- dplyr::bind_cols(
      dplyr::select(des_list[["2"]],  dplyr::all_of(join)),
      dplyr::select(des_list[["2"]], -dplyr::all_of(join)) %>% 
        purrr::set_names(paste(names(.),"U", sep="_"))
    )
    ret2_cnt <- print_tib(ret2_tib,funcTag, verbose,vt+4,tc, n="InfII")
    
    # Bind Infinium I/II into single manifest::
    ret_tib <- dplyr::bind_rows(ret1_tib, ret2_tib) %>%
      dplyr::arrange(ord_id)
    
    safe_write(ret_tib,"csv",csv,
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
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

tib_to_fasta = function(tib, prb_key="prb_seq", 
                        add_key="Address", des_key="prb_des", type_key="prb_type",
                        prb_fas=NULL,dat_csv=NULL,u49_tsv=NULL,u50_tsv=NULL,
                        # dir=NULL, runName=NULL, 
                        del="_",
                        verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'tib_to_fasta'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; key={prb_key}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    prb_sym <- rlang::sym(prb_key)
    des_sym <- rlang::sym(des_key)
    aln_vec <- c(add_key,des_key,type_key)
    
    # Build Alignment Keys/Seqs
    ret_tib <- tib %>% # tail(n=200) %>%
      tidyr::unite(aln_key, dplyr::all_of(aln_vec),sep=del,remove=FALSE) %>%
      dplyr::mutate(
        aln_seq=deMs(!!prb_sym, uc=TRUE),
        aln_rev=revCmp(aln_seq),
        aln_u49=stringr::str_sub(aln_seq, 2)
      ) %>%
      dplyr::distinct(aln_key,aln_seq, .keep_all=TRUE) %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor), as.character) )
    
    # Generate Remaining Data::
    if (!is.null(prb_fas))
      fas_vec <- ret_tib %>%
      dplyr::mutate(fas_line=paste0(">",aln_key,"\n",aln_seq) ) %>%
      dplyr::pull(fas_line)
    
    if (!is.null(u49_tsv))
      u49_tib <- ret_tib %>% 
      dplyr::filter(!!des_sym == '2') %>%
      dplyr::select(aln_key,aln_u49) %>%
      dplyr::arrange(aln_u49)
    
    if (!is.null(u50_tsv))
      u50_tib <- ret_tib %>% 
      dplyr::filter(!!des_sym == 'U') %>%
      dplyr::select(aln_key,aln_seq) %>%
      dplyr::arrange(aln_seq)
    
    # Safe Write Outputs::
    safe_write(ret_tib,"csv",dat_csv,  
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    safe_write(fas_vec,"line",prb_fas, 
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    safe_write(u49_tib,"tsv",u49_tsv,cols=FALSE, 
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    safe_write(u50_tib,"tsv",u50_tsv,cols=FALSE, 
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
#                        Bowtie Alignment Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_bsmap = function(exe, fas, gen, bsp, 
                     opt=NULL, lan=NULL, run=FALSE,
                     verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'run_bsmap'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; bsp={bsp}.{RET}"))
  
  aln_bsp <- bsp
  
  aln_ssh <- NULL
  ret_tsv <- NULL
  sys_ret <- "Not_Run"
  stime <- system.time({
    if (is.null(opt) || length(opt)==0)
      opt="-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
    
    build_file_dir(aln_bsp)
    aln_ssh <- paste(aln_bsp, 'run.sh', sep='.')
    ret_tsv <- paste(aln_bsp, 'tsv.gz', sep='.')
    
    aln_cmd <- glue::glue("{exe} -a {fas} -d {gen} {opt} -o {aln_bsp}{RET}",
                          "cat {aln_bsp} | ",
                          # "cut -f 1,2,4-11 | ",
                          "gzip -c - > {ret_tsv}{RET}",
                          "rm -f {aln_bsp}{RET}")
    
    safe_write(aln_cmd,"line",aln_ssh,
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    Sys.chmod(paths=aln_ssh, mode="0777")
    
    if (run) {
      fas_name <- fas %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
      gen_name <- gen %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Will Execute bsmap: {fas_name} vs. {gen_name}...{RET}"))
      
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
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; sys={sys_ret}; elapsed={etime}.{RET}{RET}"))
  
  if (run) return(ret_tsv)
  aln_ssh
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             BSMAP Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_bsmap = function(bsp, sort=FALSE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_bsmap'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    bsp_cols <- cols(
      bsp_key  = col_character(),
      bsp_seq  = col_character(),
      bsp_qual = col_character(),
      bsp_tag  = col_character(),
      
      bsp_chr  = col_character(),
      bsp_beg  = col_integer(),
      bsp_srd  = col_character(),
      bsp_mis  = col_integer(),
      
      bsp_ref  = col_character(),
      bsp_gap  = col_integer(),
      
      bsp_str  = col_character()
    )
    
    # Load BSP
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading BSP={bsp}...{RET}"))
    
    ret_tib <- 
      readr::read_tsv(bsp,col_names=names(bsp_cols$cols),col_types=bsp_cols) %>%
      dplyr::select(-bsp_qual) %>%
      dplyr::mutate(bsp_chr=paste0('chr',stringr::str_remove(bsp_chr,'^chr')))
    
    # Sort by genomic position::
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(bsp_chr, bsp_beg)
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

join_bsmap = function(add, bsp=NULL, file=NULL, 
                      join_key, join_type="inner", prb_des="Ord_Des",
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
    
    # Load from file if provided::
    if (!is.null(file))
      bsp <- load_bsmap(bsp=file, sort=sort,
                        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    if (is.null(bsp)) {
      stop(glue::glue("[{funcTag}]:{tabsStr} ERROR: bsp is null!!!{RET}{RET}"))
      return(ret_tib)
    }
    print_tib(add,funcTag, verbose,vt+4,tc, n="address_tib")

    prb_des_sym  <- rlang::sym(prb_des)
    bsp_join_key <- names(bsp)[1]
    join_key_sym <- rlang::sym(join_key)
    bsp_join_sym <- rlang::sym( bsp_join_key )
    
    # Rename bsp_key to join key and join::
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
    print(ret_tib)
    
    #
    # Calculate new fields from join::
    #
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Calculating new fields from joined data...{RET}"))
    
    ret_tib <- ret_tib %>%
      dplyr::mutate(
        # Generate bisulfite converted Reference Sequences::
        bsp_refU=bscUs(bsp_ref,uc=TRUE),
        bsp_refM=bscMs(bsp_ref,uc=TRUE),
        bsp_refD=bscDs(bsp_ref,uc=TRUE),
        
        # Add Alignment Orientation wrt Probe Sequence/BSP Seq::
        #  This is used to determine NxB, CpG matching, etc.
        bsp_prb_dir=dplyr::case_when(
          aln_seq==bsp_seq ~ "f",
          aln_rev==bsp_seq ~ "r",
          TRUE ~ NA_character_),
        
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
          !!prb_des_sym=="M" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ NB_F1,
          !!prb_des_sym=="M" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ NB_R1,
          
          !!prb_des_sym=="U" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ NB_F1,
          !!prb_des_sym=="U" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ NB_R1,
          
          !!prb_des_sym=="2" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ NB_F2,
          !!prb_des_sym=="2" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ NB_R2,
          TRUE ~ NA_character_
        ), # %>% stringr::str_to_upper(),
        
        # Set Expected target CG di-nucleotide from Alignment Orientation::
        bsp_din_ref=dplyr::case_when(
          !!prb_des_sym=="M" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ CG_F1,
          !!prb_des_sym=="M" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ CG_R1,
          
          !!prb_des_sym=="U" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ CG_F1,
          !!prb_des_sym=="U" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ CG_R1,
          
          !!prb_des_sym=="2" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ CG_F2,
          !!prb_des_sym=="2" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ CG_R2,
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
          !!prb_des_sym=="M" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ NB_F1M,
          !!prb_des_sym=="M" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ NB_R1M,
          
          !!prb_des_sym=="U" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ NB_F1U,
          !!prb_des_sym=="U" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ NB_R1U,
          
          !!prb_des_sym=="2" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ NB_F2D,
          !!prb_des_sym=="2" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ NB_R2D,
          TRUE ~ NA_character_
        ), # %>% stringr::str_to_upper(),
        
        # Set Expected target CG di-nucleotide from Alignment Orientation::
        bsp_din_bsc=dplyr::case_when(
          !!prb_des_sym=="M" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ CG_F1M,
          !!prb_des_sym=="M" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ CG_R1M,
          
          !!prb_des_sym=="U" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ CG_F1U,
          !!prb_des_sym=="U" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ CG_R1U,
          
          !!prb_des_sym=="2" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ CG_F2D,
          !!prb_des_sym=="2" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ CG_R2D,
          TRUE ~ NA_character_
        ) %>% stringr::str_to_upper(),
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                      Correct CG# Alignment Position::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        #
        # TBD:: Need to understand how bsp_pos can become NA_real !!!
        #
        # Update Correct Genomic CG# Location based on alignment orientation::
        #
        bsp_pos=dplyr::case_when(
          !!prb_des_sym=="M" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ bsp_beg +48,
          !!prb_des_sym=="M" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ bsp_beg + 0,
          
          !!prb_des_sym=="U" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ bsp_beg +48,
          !!prb_des_sym=="U" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ bsp_beg + 0,
          
          !!prb_des_sym=="2" & (bsp_srd=="--" | bsp_srd=="++") & bsp_prb_dir=="f" ~ bsp_beg +49,
          !!prb_des_sym=="2" & (bsp_srd=="+-" | bsp_srd=="-+") & bsp_prb_dir=="r" ~ bsp_beg - 1,
          TRUE ~ NA_real_
        ) %>% as.integer(),
        bsp_din_scr=dplyr::case_when(
          bsp_srd=="--" & bsp_din_ref=="CG" ~ 0,
          bsp_srd=="+-" & bsp_din_ref=="CG" ~ 1,
          bsp_srd=="-+" & stringr::str_starts(bsp_din_ref,"G") ~ 2,
          bsp_srd=="++" & stringr::str_ends(bsp_din_ref,"C") ~ 3,
          TRUE ~ 4) %>% as.integer()
      ) %>% 
      # Create Unique Aln Key for multiple hits::
      #
      dplyr::group_by(!!join_key_sym) %>%
      dplyr::mutate(unq_aln_key=paste(!!join_key_sym, dplyr::row_number(), sep="_")) %>%
      dplyr::ungroup() %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    # Sort by genomic position::
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(bsp_chr, bsp_beg)
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

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
    # TBD:: Temp Fix for bsp_pos == NA !!!
    #
    
    bsp <- bsp %>% dplyr::filter(!is.na(bsp_pos))
    
    ret_grs <- 
      GRanges(
        seqnames = Rle(paste0("chr",bsp$bsp_chr)),
        strand=Rle(stringr::str_sub( bsp$bsp_srd, 1,1 ) ),
        
        # ord_id  = bsp$ord_id,
        
        # prb_add = bsp$prb_add,
        # prb_add = bsp$Address,
        
        prb_cgn = bsp$prb_cgn,
        # prb_srd = bsp$prb_srd,
        prb_des = bsp$prb_des,
        # prb_src = bsp$prb_src,
        
        prb_ord_seq  = bsp$prb_seq,
        prb_aln_50U  = bsp$prb_seq %>%
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
        
        # ord_id  = bsp$ord_id,
        
        # prb_add = bsp$prb_add,
        # prb_add = bsp$Address,
        
        prb_cgn = bsp$prb_cgn,
        # prb_srd = bsp$prb_srd,
        prb_des = bsp$prb_des,
        # prb_src = bsp$prb_src,
        
        prb_ord_seq  = bsp$prb_seq,
        prb_aln_50U  = bsp$prb_seq %>%
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

format_ORD = function(tib, idx_key="IDX", uniq=TRUE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_ORD'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    idx_sym <- rlang::sym(idx_key)
    add_col  <- c(idx_key,"ord_id","prb_des","prb_seq","prb_par","prb_col")
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
    if (uniq) ret_tib <- ret_tib %>% dplyr::distinct(prb_seq, .keep_all=TRUE)
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
                      Sequence=bo_seq) # %>%
      # dplyr::mutate(Sequence=stringr::str_to_upper(Sequence),
      #               tan_seq=stringr::str_sub(Sequence,1,45),
      #               prb_seq=stringr::str_sub(Sequence,46)
      # ) %>%
      # dplyr::filter(!is.na(Address)) %>%
      # dplyr::select(Address,!!idx_sym,prb_seq, dplyr::everything()) %>% 
      # dplyr::arrange(plyr::desc(!!idx_sym))
      
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
      
      # ret_tib <- ret_tib %>% 
      #   dplyr::mutate(Sequence=stringr::str_to_upper(Sequence),
      #                 tan_seq=stringr::str_sub(Sequence,1,45),
      #                 prb_seq=stringr::str_sub(Sequence,46)
      #   ) %>%
      #   dplyr::filter(!is.na(Address)) %>%
      #   dplyr::select(Address,!!idx_sym,prb_seq, dplyr::everything()) %>% 
      #   dplyr::arrange(plyr::desc(!!idx_sym))
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed to match 'addres_names' or 'Address'!!!{RET}{RET}"))
      return(NULL)
    }
    
    ret_tib <- ret_tib %>% 
      dplyr::mutate(Sequence=stringr::str_to_upper(Sequence),
                    tan_seq=stringr::str_sub(Sequence,1,45),
                    prb_seq=stringr::str_sub(Sequence,46)
      ) %>%
      dplyr::filter(!is.na(Address)) %>%
      dplyr::select(Address,!!idx_sym,prb_seq, dplyr::everything()) %>% 
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
#                     Manifest File Methods:: Loading
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

guess_man_file = function(file, 
                          fields=c("Assay_Design_Id","Plate","address_names","Address"), 
                          cols=NULL,
                          n_max=100, # guess=100000,
                          verbose=0,vt=6,tc=1,tt=NULL) {
  funcTag <- 'guess_man_file'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_val <- NULL
  stime <- system.time({
    
    file_del_str <- guess_file_del(file, n_max=n_max, verbose=verbose,vt=vt+4)
    if (is.null(file_del_str)) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: file_del_str=NULL!!!{RET}{RET}"))
      return(NULL)
    }

    dat_tib <- readr::read_lines(file, n_max=n_max) %>%
      tibble::as_tibble() %>% 
      dplyr::mutate(row_num=dplyr::row_number())

    ret_tib <- NULL
    for (field in fields) {
      if (verbose>=vt+6)
        cat(glue::glue("[{funcTag}]:{tabsStr} field={field}{RET}"))
      
      min_tib <- dat_tib %>% dplyr::filter(stringr::str_starts(value, field)) %>% tail(n=1)
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

load_man_file = function(file, idx=NULL,
                         n_max=100, guess=1000,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_man_file'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  # Saftey Check for Empty Files::
  if (is.null(file) || length(file)==0)
    return(base::invisible(NULL))
  
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
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
  sel_cols$ord  <- c("Assay_Design_Id","AlleleA_Probe_Sequence","AlleleB_Probe_Sequence","Normalization_Bin")
  sel_cols$mat1 <- c("Address","Sequence")
  sel_cols$mat2 <- c("address_names","bo_seq")
  # sel_cols$mat2 <- c("address_names","sequence","probe_id")
  sel_cols$aqp  <- c("Address","Decode_Status")
  sel_cols$pqc  <- c("Address","Status")
  
  key_cols <- list()
  key_cols$ord  <- c("Ord_Key","Ord_PrbA","Ord_PrbB","Ord_Norm")
  key_cols$mat1 <- c("Address","Sequence")
  key_cols$mat2 <- c("Address","Sequence")
  key_cols$aqp  <- c("Address","Decode_Status")
  key_cols$pqc  <- c("Address","Decode_Status")

  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    guess_tib <- guess_man_file(file, n_max=n_max,
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
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Both beg_key AND col_num are NULL!!!{RET}{RET}"))
      return(ret_tib)
    } else if (beg_key==names(val_cols$ord$cols)[1] && col_num==length(val_cols$ord$cols)) {
      dat_key <- "ord"
      
      # val_col <- val_cols$ord
      # sel_col <- sel_cols$ord
      # key_col <- key_cols$ord
    } else if (beg_key==names(val_cols$mat1$cols)[1] && col_num==length(val_cols$mat1$cols)) {
      dat_key <- "mat1"
      
      # val_col <- val_cols$mat1
      # sel_col <- sel_cols$mat1
      # key_col <- key_cols$mat1
    } else if (beg_key==names(val_cols$mat2$cols)[1] && col_num==length(val_cols$mat2$cols)) {
      dat_key <- "mat2"
      
      # val_col <- val_cols$mat2
      # sel_col <- sel_cols$mat2
      # key_col <- key_cols$mat2
    } else if (beg_key==names(val_cols$aqp$cols)[1] && col_num==length(val_cols$aqp$cols)) {
      dat_key <- "aqp"
      
      # val_col <- val_cols$aqp
      # sel_col <- sel_cols$aqp
      # key_col <- key_cols$aqp
    } else if (beg_key==names(val_cols$pqc$cols)[1] && col_num==length(val_cols$pqc$cols)) {
      dat_key <- "pqc"
      
      # val_col <- val_cols$pqc
      # sel_col <- sel_cols$pqc
      # key_col <- key_cols$pqc
    } else {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Failed to match beg_key({beg_key}) AND ",
                      "col_num({col_num}) to known formats!!!{RET}{RET}"))
      return(ret_tib)
    }
    
    # This sets all the proper valid col types, col selection and col renaming::
    #   TBD:: Remove all the commented stuff above...
    #
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
        dplyr::mutate(Ord_Din=stringr::str_sub(Ord_Key, 1,2))
    }
    
    # Apply Formatting Rules:: Match/AQP/PQC Data Types:: Address
    if (dat_key=="mat1" || dat_key=="mat2" ||
        dat_key=="aqp"  || dat_key=="pqc") {
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
    if (dat_key=="mat1" || dat_key=="mat2") {
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
    # if (verbose>=vt+4) {
    #   cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
    #   print(ret_tib, n=n_max)
    # }
    
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
