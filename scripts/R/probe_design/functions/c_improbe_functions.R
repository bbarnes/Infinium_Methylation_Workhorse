
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   improbe (Infinium Methylation) Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

suppressWarnings(suppressPackageStartupMessages(require("Biostrings")) )

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
#
#                Infinium Methylation Probe Design Methods::
#             c-improbe = c++ improbe (traditional) via docker image
#                    Includes Thermodynamic Calculations
#                   Only Designs Infinium I U/M Probes 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

c_improbe_workflow = function(tib,
                              verbose=0,vt=3,tc=1,tt=NULL,
                              funcTag='c_improbe_workflow') {
  
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
#                          Docker improbe Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_improbe_docker = function(file, out=NULL, name="unk", image, shell,
                              suffix='improbe-design',
                              reload=FALSE,
                              verbose=1,vt=3,tc=1,tt=NULL,
                              funcTag = 'run_improbe_docker') {
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}      out={out}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     file={file}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     name={name}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}    image={image}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}    shell={shell}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   reload={reload}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  out_dir <- out
  if (is.null(out)) out_dir <- base::dirname(file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
  
  ret_log <- file.path(out_dir, paste(name,suffix,'log', sep='.'))
  ret_tsv <- file.path(out_dir, paste(name,suffix,'tsv.gz', sep='.'))
  base_file <- base::basename( file )
  
  if (reload && 
      file.exists(file) &&
      file.exists(ret_log) &&
      file.exists(ret_tsv) &&
      file.mtime(file) <= file.mtime(ret_log) &&
      file.mtime(file) <= file.mtime(ret_tsv)) {
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Reloading...{RET}"))
    
    ret_cnt <- 0
    etime   <- 0
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    return(ret_cnt)
  }
  
  ret_cnt <- 1
  stime <- base::system.time({
    
    if (file.exists(ret_log)) unlink(ret_log)
    if (file.exists(ret_tsv)) unlink(ret_tsv)
    
    base::system(glue::glue("touch {ret_log}"))
    base::system(glue::glue("touch {ret_tsv}"))
    
    imp_doc_cmd <- glue::glue("docker run -i --rm ",
                              "-v {out_dir}:/input -v {out_dir}:/output -w /work ",
                              "{image} {shell} {base_file} {name}")
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: Running improbe cmd='{imp_doc_cmd}'...{RET}"))
    ret_cnt <- base::system(imp_doc_cmd)
    
    if (ret_cnt != 0) {
      cat(glue::glue("{RET}[{funcTag}]: ERROR: cmd return={ret_cnt} cmd='{imp_doc_cmd}'!{RET}{RET}"))
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_cnt
}

load_improbe_design = function(
  file, out=NULL, join=NULL, 
  join_new=c("Aln_Key_Unq","Bsp_Chr","Bsp_Pos","Bsp_FR","Bsp_CO"),
  join_old=c("Seq_ID","Chromosome","Coordinate","Strand_FR","Strand_CO"),
  level=0, add_inf=TRUE, verbose=0,vt=3,tc=1,tt=NULL, 
  funcTag='load_improbe_design') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; level={level}...{RET}"))
  
  # improbe output original fields::
  #
  org_cols <- NULL
  new_cols <- NULL
  if (level>0) {
    # Level 1::
    ids_org_cols <- c("Seq_ID", "Genome_Build", "Chromosome", "Coordinate")
    seq_org_cols <- c("Forward_Sequence", "Top_Sequence")
    prb_org_cols <- c("Methyl_Probe_Covered_Top_Sequence", "UnMethyl_Probe_Sequence", "Methyl_Probe_Sequence")
    srd_org_cols <- c("Methyl_Allele_FR_Strand", "Methyl_Allele_TB_Strand", "Methyl_Allele_CO_Strand", "Methyl_Next_Base")
    
    # Level 2:: Final Scores::
    scrF_org_cols <- c("UnMethyl_Final_Score", "Methyl_Final_Score")
    
    # Level 3:: Both Metric/Score::
    rawB_org_cols <- c("Methyl_Underlying_CpG_Count", "Methyl_Underlying_CpG_Min_Dist")
    scrB_org_cols <- c("Methyl_Underlying_CpG_Score", "Methyl_Next_Base_Score")
    
    # Level 4:: Metrics
    rawM_org_cols <- c("Methyl_Tm",    "Methyl_GC_Percent", "Methyl_13mer_Count", "Methyl_Address_Count",
                       "Methyl_Self_Complementarity",       "Methyl_Mono_Run",    "Methyl_Ectopic_Count")
    rawU_org_cols <- paste0("Un",rawM_org_cols)
    
    # Level 5:: Scores
    scrM_org_cols <- c("Methyl_Tm_Score","Methyl_GC_Score",   "Methyl_13mer_Score",  "Methyl_Address_Score",
                       "Methyl_Self_Complementarity_Score", "Methyl_Mono_Run_Score", "Methyl_Ectopic_Score")
    scrU_org_cols <- paste0("Un",scrM_org_cols)
    
    # improbe output new fields::
    #
    ids_new_cols <- c("Seq_ID", "Genome_Build", "Chromosome", "Coordinate")
    seq_new_cols <- c("Forward_Sequence", "Top_Sequence")
    prb_new_cols <- c("Probe_Seq_T", "Probe_Seq_U", "Probe_Seq_M")
    srd_new_cols <- c("Strand_FR", "Strand_TB", "Strand_CO", "Next_Base")
    
    # Level 2:: Final Scores
    scrF_new_cols <- c("Scr_U", "Scr_M")
    
    # Level 3: Both Metric/Score::
    rawB_new_cols <- c("Cpg_Cnt", "Cpg_Dis")
    scrB_new_cols <- c("Cpg_Scr", "Nxb_Scr")
    
    # Level 4:: New Metrics::
    met_new_cols  <- c("Tm", "GC", "Mer13", "Address", "SelfCmpl", "Mono", "Ectopic")
    rawU_new_cols <- paste(met_new_cols, "Raw_U", sep="_")
    rawM_new_cols <- paste(met_new_cols, "Raw_M", sep="_")
    
    # Level 5:: New Metrics::
    scrU_new_cols <- paste(met_new_cols, "Scr_U", sep="_")
    scrM_new_cols <- paste(met_new_cols, "Scr_M", sep="_")
    
    #
    # Old/New Field Names to Select/Return::
    #
    if (level>=1) {
      org_cols <- c(org_cols, ids_org_cols,seq_org_cols,prb_org_cols,srd_org_cols)
      new_cols <- c(new_cols, ids_new_cols,seq_new_cols,prb_new_cols,srd_new_cols)
    }
    if (level>=2) {
      org_cols <- c(org_cols, scrF_org_cols)
      new_cols <- c(new_cols, scrF_new_cols)
    }
    if (level>=3) {
      org_cols <- c(org_cols, rawB_org_cols,scrB_org_cols)
      new_cols <- c(new_cols, rawB_new_cols,scrB_new_cols)
    }
    if (level>=4) {
      org_cols <- c(org_cols, rawU_org_cols,rawM_org_cols)
      new_cols <- c(new_cols, rawU_new_cols,rawM_new_cols)
    }
    if (level>=5) {
      org_cols <- c(org_cols, scrU_org_cols,scrM_org_cols)
      new_cols <- c(new_cols, scrU_new_cols,scrM_new_cols)
    }
    
    if (verbose>=vt+6) {
      org_len <- org_cols %>% length()
      cat(glue::glue("[{funcTag}]:{tabsStr} org_cols(org_len)={RET}"))
      print(org_cols)
      
      new_len <- new_cols %>% length()
      cat(glue::glue("[{funcTag}]:{tabsStr} new_cols(new_len)={RET}"))
      print(new_cols)
    }
  } else {
    cat(glue::glue("[{funcTag}]:{tabsStr} Returning full data (level={level})...{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading designs={file}...{RET}"))
    
    ret_tib <- 
      suppressMessages(suppressWarnings( readr::read_tsv(file) ))
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="Raw_Improbe_Data")
    
    # Now subset and rename based on input level::
    #
    if (!is.null(org_cols) && !is.null(new_cols)) {
      ret_tib <- ret_tib %>% 
        dplyr::select(dplyr::all_of(org_cols)) %>%
        dplyr::mutate(Methyl_Allele_TB_Strand=stringr::str_sub(Methyl_Allele_TB_Strand, 1,1)) %>%
        purrr::set_names(new_cols)
    }
    ret_tib <- ret_tib %>% clean_tibble()
    
    if (add_inf && level>=3) {
      #
      # Adding some extra required fields for matching later::
      #
      ret_tib <- ret_tib %>%
        dplyr::rename(Strand_Ref_FR=Strand_FR) %>%
        dplyr::mutate(
          Scr_Min=pmin(Scr_U,Scr_M),
          Inf_Type=dplyr::case_when(
            Scr_Min<0.2                  ~ 0,
            Scr_Min<0.3 & Strand_CO=="C" ~ 0,
            
            Scr_Min< 0.3 & Strand_CO=="O" & Cpg_Cnt==0 ~ 1,
            Scr_Min< 0.4 & Strand_CO=="O" & Cpg_Cnt==1 ~ 1,
            Scr_Min< 0.5 & Strand_CO=="O" & Cpg_Cnt==2 ~ 1,
            Scr_Min< 0.6 & Strand_CO=="O" & Cpg_Cnt==3 ~ 1,
            
            Scr_Min< 0.4 & Strand_CO=="C" & Cpg_Cnt==0 ~ 1,
            Scr_Min< 0.5 & Strand_CO=="C" & Cpg_Cnt==1 ~ 1,
            Scr_Min< 0.6 & Strand_CO=="C" & Cpg_Cnt==2 ~ 1,
            Scr_Min< 0.7 & Strand_CO=="C" & Cpg_Cnt==3 ~ 1,
            
            Scr_Min>=0.7 ~ 1,
            
            TRUE ~ 2
          ),
          Imp_U49=dplyr::case_when(
            Strand_CO=="C" ~ stringr::str_sub(Probe_Seq_U,1,49),
            Strand_CO=="O" ~ stringr::str_sub(Probe_Seq_U,2,50),
            TRUE ~ NA_character_
          ),
          Imp_M49=dplyr::case_when(
            Strand_CO=="C" ~ stringr::str_sub(Probe_Seq_M,1,49),
            Strand_CO=="O" ~ stringr::str_sub(Probe_Seq_M,2,50),
            TRUE ~ NA_character_
          ),
          Strand_FR=dplyr::case_when(
            Strand_Ref_FR=="F" ~ "R",
            Strand_Ref_FR=="R" ~ "F",
            TRUE ~ NA_character_
          )
        )
    }
    
    # Merge data back with BSP results::
    if (!is.null(join)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Joining data with improbe data...{RET}"))
      
      ret_tib <- dplyr::left_join(
        join, dplyr::rename_with(ret_tib, ~ join_new, dplyr::all_of(join_old) ),
        by=join_new, suffix=c("_bsp","_imp")
      )
    }
    
    out_cnt <- safe_write(x=ret_tib,file=out,funcTag=funcTag,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("Clean_Improbe_Data({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# End of file
