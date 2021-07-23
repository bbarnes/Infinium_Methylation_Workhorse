
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   improbe (Infinium Methylation) Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

suppressWarnings(suppressPackageStartupMessages(require("Biostrings")) )

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
#                           Comparison Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

compare_probes = function(ref, can, ref_keys, can_keys,
                          by=c("Seq_ID","Strand_FR","Strand_CO",
                               "Ord_Des","Ord_Din"),
                          grp=c("Ord_Des","Ord_Din","Strand_FR","Strand_CO"),
                          retData=FALSE, pivot=FALSE,
                          verbose=0,vt=4,tc=1,tt=NULL,
                          funcTag='compare_probes') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  ret_tab <- NULL
  stime <- base::system.time({
    
    # Inner Join::
    inn_tib <- dplyr::inner_join(ref,can, by=by,suffix=c("_ref","_can"))
    ret_key <- glue::glue("inn-tib({funcTag})")
    ret_cnt <- print_tib(inn_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    for (ii in c(1:length(ref_keys))) {
      fieldA <- ref_keys[ii]
      fieldB <- can_keys[ii]
      
      # Fix any names that were duplicated::
      if (!fieldA %in% names(inn_tib))
        fieldA <- paste0(fieldA,"_ref")
      
      if (!fieldB %in% names(inn_tib))
        fieldB <- paste0(fieldB,"_can")
      
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} ii={ii}, A={fieldA}, B={fieldB}{RET}"))
      
      sel_vec <- c(by, fieldA, fieldB)
      new_vec <- c(by, "Probe_A", "Probe_B", 
                   "Man_MisMatch", "DI_NUC_AB", "Man_TarMatch")
      
      ret_tib <- inn_tib %>% 
        dplyr::select(dplyr::all_of(sel_vec)) %>%
        purrr::set_names(sel_vec) %>%
        cmpInfII_MisMatch(fieldA=fieldA,fieldB=fieldB, 
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
        purrr::set_names(new_vec) %>%
        dplyr::bind_rows(ret_tib)
      
      # Exact Match::
      #
      # inn_tib %>% 
      #   dplyr::select(dplyr::all_of(sel_vec)) %>%
      #   purrr::set_names(sel_vec) %>%
      #   cmpInfII(fieldA=fieldA,fieldB=fieldB, 
      #            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% print()
      
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr}{BRK}{RET}{RET}"))
    }
    
    if (pivot) ret_tab <- ret_tib %>% tidyr::pivot_longer(
      cols=c("Probe_A","Probe_B"),names_to="Prb_Source", values_to="Probe_Seq")
    
    if (retData) ret_dat$cmp <- ret_tib
    if (retData) ret_dat$tab <- ret_tab
    
    sum_key_vec <- c("Man_MisMatch", "Man_TarMatch", "DI_NUC_AB")
    for (key in sum_key_vec) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Summarizing {key}...{RET}"))
      
      grp_vec <- c(key,grp)
      sum_tib <- ret_tib %>% 
        dplyr::group_by(dplyr::across(dplyr::all_of(grp_vec)) ) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      
      if (retData) ret_dat[[key]] <- sum_tib
      
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Done.Summarizing {key}.{RET}{RET}"))
    }
    
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
#                           IMP Designs Workflow::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imp_designs_workflow = function(
  tib, max=0,
  
  # improbe File Parameters::
  imp_inp_tsv, imp_des_tsv=NULL, imp_fin_tsv=NULL, imp_seq_csv=NULL,
  
  # Genomes Parameters::
  gen_bld="UNK",
  gen_ref_fas, 
  bsc_ref_fas=NULL,
  gen_snp_fas=NULL, 
  bsc_snp_fas=NULL,
  
  # Field Parameters::
  ids_key="Aln_Key_Unq", 
  din_key="Ord_Din", 
  pos_key="Bsp_Pos", 
  chr_key="Bsp_Chr", 
  
  srsplit=FALSE,
  srd_key=NULL,
  cosplit=FALSE,
  cos_key=NULL,
  
  # Docker Parameters::
  run_name, 
  doc_image,
  doc_shell="run_improbe.sh",
  
  join_new=c("Aln_Key_Unq","Bsp_Chr","Bsp_Pos","Bsp_FR","Bsp_CO"),
  join_old=c("Seq_ID","Chromosome","Coordinate","Strand_FR","Strand_CO"),
  
  reload=FALSE,
  retData=FALSE, 
  parallel=FALSE,
  r_improbe=FALSE,
  s_improbe=FALSE,
  
  verbose=0,vt=3,tc=1,tt=NULL, funcTag='imp_designs_workflow') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  # Set file names if they were NULL
  if (is.null(imp_seq_csv)) imp_seq_csv <- imp_inp_tsv %>% 
      stringr::str_remove(".gz$") %>%
      stringr::str_remove(".[ctsv]+$") %>%
      stringr::str_remove("improbe-inputs$") %>%
      paste0("improbe-sequence.csv.gz")
  
  if (is.null(imp_des_tsv)) imp_des_tsv <- imp_inp_tsv %>% 
      stringr::str_remove(".gz$") %>%
      stringr::str_remove(".[ctsv]+$") %>%
      stringr::str_remove("improbe-inputs$") %>%
      paste0("improbe-designOutput.tsv.gz")
  
  if (is.null(imp_fin_tsv)) imp_fin_tsv <- imp_des_tsv %>% 
      stringr::str_remove(".gz$") %>%
      stringr::str_remove(".[ctsv]+$") %>%
      paste0(".clean.tsv.gz")
  
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} improbe Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   imp_inp_tsv={imp_inp_tsv}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   imp_des_tsv={imp_des_tsv}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   imp_fin_tsv={imp_fin_tsv}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Genome Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       gen_bld={gen_bld}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   gen_ref_fas={gen_ref_fas}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   bsc_ref_fas={bsc_ref_fas}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   gen_snp_fas={gen_snp_fas}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   bsc_snp_fas={bsc_snp_fas}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Field Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       ids_key={ids_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       din_key={din_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       pos_key={pos_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       chr_key={chr_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       srsplit={srsplit}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       srd_key={srd_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       cosplit={cosplit}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       cos_key={cos_key}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Docker Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}      run_name={run_name}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     doc_shell={doc_shell}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     doc_image={doc_image}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Time Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}           max={max}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}        reload={reload}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       retData={retData}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}      parallel={parallel}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     r_improbe={r_improbe}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     s_improbe={s_improbe}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  stime <- base::system.time({
    
    if (max>0) tib <- tib %>% head(n=max)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                1.0 Load Ref Genomes and run r-improbe::
    #
    #                (r-improbe = improbe re-implemented in R)
    #
    #                          Write improbe input::
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- fas_to_seq(tib=tib,
                          gen_bld=gen_bld,
                          
                          gen_ref_fas=gen_ref_fas,
                          bsc_ref_fas=bsc_ref_fas,
                          gen_snp_fas=gen_snp_fas,
                          bsc_snp_fas=bsc_snp_fas,
                          
                          imp_tsv=imp_inp_tsv, 
                          seq_csv=imp_seq_csv,
                          
                          ids_key=ids_key, 
                          din_key=din_key,
                          pos_key=pos_key,
                          chr_key=chr_key,
                          
                          srsplit=srsplit,
                          srd_key=srd_key,
                          cosplit=cosplit,
                          cos_key=cos_key,
                          
                          reload=reload,
                          # retData=retData,
                          retData=TRUE,
                          parallel=parallel,
                          r_improbe=r_improbe,
                          s_improbe=s_improbe,
                          
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    if (FALSE) {
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #
      #                 1.1 Run improbe designs:: via docker
      #
      #                     (improbe = original c++ version)
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ret_val <- run_improbe_docker(file=imp_inp_tsv, 
                                    name=run_name, 
                                    image=doc_image, 
                                    shell=doc_shell,
                                    reload=reload,
                                    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                     1.2 Load improbe designs:: filter
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ret_tib <- load_improbe_design(file=imp_des_tsv, out=imp_fin_tsv, 
                                     join=ret_tib,
                                     join_new=join_new,join_old=join_old,
                                     level=3, add_inf=TRUE, 
                                     verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      if (retData) ret_dat$fwd <- fwd_tib
      if (retData) ret_dat$imp <- ret_tib
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #
      #                1.3 Load BSC Genomes and run s-improbe::
      #
      #                    (s-improbe = sub-string improbe)
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # build=c("Prb1C","Prb2C","Prb1O","Prb2O"),
      
      
      ret_key <- glue::glue("ret-FIN({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
      
    }
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Infinium Methylation Probe Design Methods::
#               s-improbe = improbe by BSC genome sub-string
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

s_improbe = function(tib, build=c("Prb1C","Prb2C","Prb1O","Prb2O"),
                     verbose=0,vt=3,tc=1,tt=NULL,
                     funcTag='s_improbe') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   build={build}{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    #
    # Build individual parts of the template sequence:: 
    #
    #                                            iupac
    #     up01.up02...up11.up12.up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn12.dn11.dn02.dn01
    #                                       Nxb [  C   G  ] Nxb
    #
    # Probe Design Formulas::
    #
    # Inf1C                               Nxb60* up61------------------dn58
    # Inf2C                                     ext61* dn61-------------------dn12
    #
    # Inf1O                up12------------------up61 Nxb61*
    # Inf2O           up11------------------up60 ext61*
    #
    #
    #  TBD:: Add reverse/complements::
    #
    # Build requested sub-string probes::
    #
    ret_tib <- tib
    if ("Prb1C" %in% build)
      ret_tib <- ret_tib %>% dplyr::mutate(
        Prb1C=paste0(iupac,dn61,dn60,dn59,dn58),
        Nxb1C=paste0(up60),
        Prb1C_Len=stringr::str_length(Prb1C) )
    
    if ("Prb2C" %in% build)
      ret_tib <- ret_tib %>% dplyr::mutate(
        Prb2C=paste0(dn61,dn60,dn59,dn58,dn12),
        Nxb2C=paste0(iupac),
        Prb2C_Len=stringr::str_length(Prb2C))
    
    if ("Prb1O" %in% build)
      ret_tib <- ret_tib %>% dplyr::mutate(
        Prb1O=paste0(up12,up58,up59,up60,iupac),
        Nxb1O=paste0(dn61),
        Prb1O_Len=stringr::str_length(Prb1O) )
    
    if ("Prb2O" %in% build)
      ret_tib <- ret_tib %>% dplyr::mutate(
        Prb2O=paste0(up11,up12,up58,up59,up60),
        Nxb2O=paste0(iupac),
        Prb2O_Len=stringr::str_length(Prb2O))
    
    # NOTE:: I think you can ignore FR strand since the template sequence
    #   should already be 5' -> 3' for the strand of interest
    #
    # if ("Prb2_RC" %in% build)
    #   ret_tib <- ret_tib %>% dplyr::mutate(
    #     Prb2_RC=paste0(up12,up58,up59,up60,up61), 
    #     Prb2_RC_Len=stringr::str_length(Prb2_RC))
    # 
    # if ("Prb1_RC" %in% build)
    #   ret_tib <- ret_tib %>% dplyr::mutate(
    #     Prb1_RC=paste0(up58,up59,up60,up61,iupac), 
    #     Prb1_RC_Len=stringr::str_length(Prb1_RC))
    
    # NOTE::This is just for graphical sanity checks and should be removed
    #   once its validated...
    #
    ret_tib <- ret_tib %>% dplyr::mutate(
      Temp=paste(up01,up02,up11,up12,up58,up59,up60,up61,
                 dn61,dn60,dn59,dn58,dn12,dn11,dn02,dn01, sep=""),
      Temp_Len=stringr::str_length(Temp),
      Temp=addBrac(Temp),
      Temp=paste(up01,up02,up11,up12,up58,up59,up60,up61,
                 dn61,dn60,dn59,dn58,dn12,dn11,dn02,dn01, sep=" ") )    
    
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

# s_improbe_template_workflow(tibs=chr_list[[chr_str]], maps=chr_maps[[chr_str]], seqs=dna_dat$seqs, chr_str=chr_str, chr_key=chr_key)
#
s_improbe_template_workflow = function(tib, 
                                       seq,
                                       
                                       srd_str="F",
                                       pos_key="Coordinate",
                                       chr_key="Chromosome",
                                       chr_str,
                                       
                                       ext_seq="Ext_Forward_Seq",
                                       iup_seq="Iupac_Forward_Sequence",
                                       imp_seq="Forward_Sequence",
                                       
                                       ups_len=60, 
                                       seq_len=122, 
                                       iupac=NULL,
                                       del="_",
                                       
                                       verbose=0,vt=3,tc=1,tt=NULL,
                                       funcTag='s_improbe_template_workflow') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Field Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   srd_str={srd_str}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   pos_key={pos_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   chr_key={chr_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   chr_str={chr_str}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   ext_seq={ext_seq}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   iup_seq={iup_seq}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   imp_seq={imp_seq}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     iupac={iupac}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   ups_len={ups_len}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   seq_len={seq_len}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # map_key <- glue::glue("map-tib({funcTag})")
    # map_cnt <- print_tib(map,funcTag, verbose,vt=0,tc, n=map_key)
    
    
    chr_sym <- rlang::sym(chr_key)
    cur_tib <- tib
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 Parse Forward Template Sequence from Genome::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ref_seqs <- parse_genomic_seqs(tib=cur_tib,
                                   seq=as.character(seq),
                                   
                                   srd_str=srd_str,
                                   pos_key=pos_key,
                                   chr_str=chr_str,
                                   
                                   ups_len=ups_len, 
                                   seq_len=seq_len,
                                   verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                Forward Template Sequence Generation:: s-improbe
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- s_improbe_template(tib=cur_tib,
                                  seqs=ref_seqs,
                                  
                                  chr_str=chr_str, 
                                  
                                  ext_seq=ext_seq,
                                  iup_seq=iup_seq,
                                  imp_seq=imp_seq,
                                  
                                  iupac=iupac, 
                                  ups_len=ups_len,
                                  seq_len=seq_len,
                                  verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
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

# s_improbe_parse(tib=cur_tib, seqs=ref_seqs, chr_str=chr_str, iupac=iupac, ext_seq=ext_seq, iup_seq=iup_seq, imp_seq=imp_seq)
s_improbe_template = function(tib, seqs, 
                              chr_str="chr0", 
                              
                              ext_seq="Ext_Forward_Seq",
                              iup_seq="Iupac_Forward_Sequence",
                              imp_seq="Forward_Sequence",
                              
                              iupac=NULL,
                              ups_len=60,
                              seq_len=122,
                              verbose=0,vt=3,tc=1,tt=NULL,
                              funcTag='s_improbe_template') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Field Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   chr_str={chr_str}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   ext_seq={ext_seq}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   iup_seq={iup_seq}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   imp_seq={imp_seq}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       iupac={iupac}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     ups_len={ups_len}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     seq_len={seq_len}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ext_sym  <- rlang::sym(ext_seq)
    iup_sym  <- rlang::sym(iup_seq)
    imp_sym  <- rlang::sym(imp_seq)
    
    #
    # Build individual parts of the template sequence:: 
    #
    #                                            iupac
    #     up01.up02...up11.up12.up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn12.dn11.dn02.dn01
    #                                       Nxb [  C   G  ] Nxb
    #
    
    ref_up01 <- stringr::str_sub(seqs, 1,1)
    ref_up02 <- stringr::str_sub(seqs, 2,2)
    
    ref_up11 <- stringr::str_sub(seqs, 3,13)
    ref_up12 <- stringr::str_sub(seqs, 14,14)
    
    ref_up58 <- stringr::str_sub(seqs, 15,ups_len-2+2)
    
    ref_up59 <- stringr::str_sub(seqs, ups_len-2+3,ups_len-2+3)
    ref_up60 <- stringr::str_sub(seqs, ups_len-2+4,ups_len-2+4)
    ref_up61 <- stringr::str_sub(seqs, ups_len-2+5,ups_len-2+5)
    
    iupac_vec <- ref_up61
    if (!is.null(iupac)) iupac_vec <- tib %>% dplyr::pull(!!iupac)
    
    ref_dn61 <- stringr::str_sub(seqs, ups_len-2+6,ups_len-2+6)
    ref_dn60 <- stringr::str_sub(seqs, ups_len-2+7,ups_len-2+7)
    ref_dn59 <- stringr::str_sub(seqs, ups_len-2+8,ups_len-2+8)
    
    ref_dn58 <- stringr::str_sub(seqs, ups_len-2+9,ups_len-2+ups_len-2+8-12)
    
    ref_dn12 <- stringr::str_sub(seqs, ups_len-2+ups_len-2+8-12+1,ups_len-2+ups_len-2+8-12+1)
    ref_dn11 <- stringr::str_sub(seqs, ups_len-2+ups_len-2+8-12+2,ups_len-2+ups_len-2+8)
    
    ref_dn02 <- stringr::str_sub(seqs, ups_len-2+ups_len-2+8,ups_len-2+ups_len-2+8)
    ref_dn01 <- stringr::str_sub(seqs, ups_len-2+ups_len-2+9,ups_len-2+ups_len-2+9)
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} seqs({chr_str})={RET}"))
      seqs %>% head(n=2) %>% print()
    }
    
    # Add all fields to current tibble::
    #
    ret_tib <- tib %>%
      dplyr::mutate(
        !!ext_sym := seqs %>% addBrac(),
        
        up01=ref_up01,
        up02=ref_up02,
        
        up11=ref_up11,
        up12=ref_up12,
        
        up58=ref_up58,
        
        up59=ref_up59,
        up60=ref_up60,
        up61=ref_up61,
        
        iupac=iupac_vec,
        
        dn61=ref_dn61,
        dn60=ref_dn60,
        dn59=ref_dn59,
        
        dn58=ref_dn58,
        
        dn12=ref_dn12,
        dn11=ref_dn11,
        
        dn02=ref_dn02,
        dn01=ref_dn01,
        
        # Now we can assemble to optimal template::
        #
        !!iup_sym := paste0(up11,up12, up58, up59,up60, iupac,dn61, dn60,dn59, dn58, dn12,up11) %>% addBrac(),
        !!imp_sym := paste0(up11,up12, up58, up59,up60,     "CG",   dn60,dn59, dn58, dn12,up11) %>% addBrac(),
        
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
    ret_key <- glue::glue("ret_tib_2=chr_str={chr_str}")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

s_improbe_trifecta = function(tib, 
                              tar_din="rs",
                              ids_key="Aln_Key_Unq",
                              din_key="Ord_Din",
                              
                              pos_key="Coordinate",
                              chr_str="chr0", 
                              
                              ext_seq="Ext_Forward_Seq",
                              iup_seq="Iupac_Forward_Sequence",
                              imp_seq="Forward_Sequence",
                              
                              verbose=0,vt=3,tc=1,tt=NULL,
                              funcTag='s_improbe_trifecta') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Field Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   tar_din={tar_din}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   ids_key={ids_key}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   din_key={din_key}{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   pos_key={pos_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   chr_str={chr_str}{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   ext_seq={ext_seq}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   iup_seq={iup_seq}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   imp_seq={imp_seq}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ids_sym  <- rlang::sym(ids_key)
    din_sym  <- rlang::sym(din_key)
    ext_sym  <- rlang::sym(ext_seq)
    iup_sym  <- rlang::sym(iup_seq)
    imp_sym  <- rlang::sym(imp_seq)
    
    ref_col_sym  <- rlang::sym(ref_col)
    alt_col_sym  <- rlang::sym(alt_col)
    iup_col_sym  <- rlang::sym(iup_col)
    
    pos_sym  <- rlang::sym(pos_key)
    
    if (verbose>=vt+1)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Upstream CG ",
                     "Flank Seuqnces({chr_str})...{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 Check for Trifecta Probes:: Upstream of SNP
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ups_tib <- tib %>% 
      dplyr::filter(!!din_sym==tar_din) %>%
      dplyr::filter(up59=="C" & up60=="G") %>%
      dplyr::mutate(
        !!din_sym := "cg",
        !!ext_sym := paste0(up01,up02,up11,up12,up58,up59,up60,  up61,dn61, dn60,dn59,dn58,dn12,dn11) %>% stringr::str_sub(1,seq_len) %>% addBrac(),
        !!iup_sym := paste0(up01,up02,up11,up12,up58,up59,up60, iupac,dn61, dn60,dn59,dn58,dn12,dn11) %>% stringr::str_sub(1,seq_len) %>% addBrac(),
        !!imp_sym := paste0(up01,up02,up11,up12,up58,               "CG",   up61,dn61,dn60,dn59,dn58,dn12,dn11) %>% stringr::str_sub(1,seq_len) %>% addBrac(),
        !!ref_col_sym := "C",
        !!alt_col_sym := "T",
        !!iup_col_sym := "Y",
        Din_Str=stringr::str_to_lower(paste0(up59,up60)),
        Des_Din="cg",
        !!ids_sym:=paste(!!ids_sym,"CG-UP", sep='-'),
        !!pos_sym:=as.integer(!!pos_sym - 2)    # Previously:: Coordinate=as.integer(Coordinate+2)
      )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #               Check for Trifecta Probes:: Downstream of SNP
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt+1)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Downstream CG ",
                     "Flank Seuqnces({chr_str})...{RET}"))
    
    dns_tib <- tib %>% 
      dplyr::filter(!!din_sym==tar_din) %>%
      dplyr::filter(dn61=="C" & dn60=="G") %>%
      dplyr::mutate(
        !!din_sym := "cg",
        !!ext_sym := paste0(up11,up12,up58,up59,up60,up61,  dn61,dn60, dn59,dn58,dn12,dn11,dn02) %>% stringr::str_sub(2) %>% addBrac(),
        !!iup_sym := paste0(up11,up12,up58,up59,up60,iupac, dn61,dn60, dn59,dn58,dn12,dn11,dn02) %>% stringr::str_sub(2) %>% addBrac(),
        !!imp_sym := paste0(up11,up12,up58,up59,up60,up61,     "CG",   dn59,dn58,dn12,dn11,dn02) %>% stringr::str_sub(2) %>% addBrac(),
        !!ref_col_sym := "C",
        !!alt_col_sym := "T",
        !!iup_col_sym := "Y",
        Din_Str=stringr::str_to_lower(paste0(dn61,dn60)),
        Des_Din="cg",
        !!ids_sym:=paste(!!ids_sym,"CG-DN", sep='-'),
        !!pos_sym:=as.integer(!!pos_sym - 2)    # Previously:: Coordinate=as.integer(Coordinate+2)
      )
    
    ret_tib <- dplyr::bind_rows(ups_tib,dns_tib)
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
#                           r-improbe re-implemented
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# NOTE:: Function has been renamed desSeq_to_prbs() -> r_improbe()
#
r_improbe = function(tib,
                     
                     sr_str='FR',
                     co_str='CO',
                     
                     ids_key, 
                     seq_key, 
                     prb_key,
                     
                     srsplit=FALSE,
                     srd_key=NULL,
                     cosplit=FALSE,
                     cos_key=NULL,
                     
                     prb_len=60, 
                     seq_len=122, 
                     
                     parallel=FALSE, 
                     add_matseq=TRUE,
                     
                     del='_',
                     max=0,
                     verbose=0,vt=3,tc=1,tt=NULL, 
                     funcTag='r_improbe') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}      sr_str={sr_str}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}      co_str={co_str}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     ids_key={ids_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     prb_key={prb_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     seq_key={seq_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     srsplit={srsplit}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     srd_key={srd_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     cosplit={cosplit}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     cos_key={cos_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     prb_len={prb_len}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     seq_len={seq_len}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}    parallel={parallel}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  add_matseq={add_matseq}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  #
  # TBD:: Add target strands and extract strands from input tib!
  #   NOTE:: This can be done by spliting the input tib and then 
  #    calling the function recursively using an input parameter::
  #  --srsplit=[TRUE/FALSE]
  #  --tar_srd=[NULL, FC, RC, FO, RO, TC, BC, TO, BO ]
  #    if (srsplit==TRUE) split and call with tar_srd
  #    if (srsplit==FALSE) run function with tar_srd if present...
  #
  
  ret_cnt <- 0
  ret_tib <- NULL
  if ( ( srsplit && !is.null(srd_key) ) || 
       ( cosplit && !is.null(cos_key) ) ) {
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Will split by strand...{RET}"))
    
    if (srsplit &&!is.null(srd_key)) {
      tib_list <- split(tib, f=dplyr::pull(tib, srd_key))
      
      for (srd in names(tib_list)) {
        if (verbose>=vt) 
          cat(glue::glue("{RET}{BRK}{RET}[{funcTag}]:{tabsStr} ",
                         "Starting::split={tc}/{srd}{RET}{RET}{RET}"))
        if (verbose>=vt+4)
          tib_list[[srd]] %>% dplyr::select(
            dplyr::all_of(c(!!ids_key,!!seq_key,!!prb_key,
                            !!srd_key,!!cos_key))) %>% print()
        
        cur_tib <- NULL
        cur_tib <- r_improbe(tib=tib_list[[srd]], 
                             sr_str=sr_str, 
                             co_str=co_str,
                             
                             ids_key=ids_key,
                             seq_key=seq_key,
                             prb_key=prb_key,
                             
                             srsplit=FALSE,
                             srd_key=srd_key,
                             cosplit=cosplit,
                             cos_key=cos_key, 
                             
                             prb_len=prb_len, 
                             seq_len=seq_len, 
                             
                             parallel=parallel, 
                             add_matseq=add_matseq,
                             
                             del=del, 
                             max=max, 
                             verbose=verbose,vt=vt,tc=tc+1,tt=tt)
        ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
        
        if (verbose>=vt) 
          cat(glue::glue("{RET}[{funcTag}]:{tabsStr} DONE::split={srd}{RET}{BRK}{RET}"))
      }
      
    } else if (cosplit && !is.null(cos_key)) {
      tib_list <- split(tib, f=dplyr::pull(tib, cos_key))
      
      for (srd in names(tib_list)) {
        if (verbose>=vt) 
          cat(glue::glue("{RET}{BRK}{RET}[{funcTag}]:{tabsStr} ",
                         "Starting::split={tc}/{srd}{RET}{RET}{RET}"))
        if (verbose>=vt+4)
          tib_list[[srd]] %>% dplyr::select(
            dplyr::all_of(c(!!ids_key,!!seq_key,!!prb_key,
                            !!srd_key,!!cos_key))) %>% print()
        
        cur_tib <- NULL
        cur_tib <- r_improbe(tib=tib_list[[srd]], 
                             sr_str=sr_str, 
                             co_str=co_str,
                             
                             ids_key=ids_key, 
                             seq_key=seq_key, 
                             prb_key=prb_key,
                             
                             srsplit=srsplit,
                             srd_key=srd_key,
                             cosplit=FALSE,
                             cos_key=cos_key,
                             
                             prb_len=prb_len, 
                             seq_len=seq_len, 
                             
                             parallel=parallel, 
                             add_matseq=add_matseq,
                             
                             del=del, 
                             max=max,
                             
                             verbose=verbose,vt=vt,tc=tc+1,tt=tt)
        ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
        
        if (verbose>=vt) 
          cat(glue::glue("{RET}[{funcTag}]:{tabsStr} DONE::split={srd}{RET}{BRK}{RET}"))
      }
    }
  } else {
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Running r-improbe...{RET}{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- system.time({
      
      valid_sr_vec <- c("FR","TB")
      if (!sr_str %in% valid_sr_vec) {
        stop(glue::glue("{RET}[{funcTag}]:{tabsStr} Invalid SR string={sr_str}!{RET}{RET}"))
        return(ret_tib)
      }
      
      # Ambiguous Source Design Sequence Strand
      sr_sym <- rlang::sym(sr_str)
      co_sym <- rlang::sym(co_str)
      sr_vec <- stringr::str_split(sr_str, '', simplify=TRUE) %>% as.vector()
      co_vec <- stringr::str_split(co_str, '', simplify=TRUE) %>% as.vector()
      
      # Extract all SR/CO's from data tib::
      tar_srd_vec <- sr_vec
      if (!is.null(srd_key))
        tar_srd_vec <- tib %>% dplyr::pull(!!srd_key) %>% unique() %>% as.vector()
      
      tar_cos_vec <- co_vec
      if (!is.null(cos_key))
        tar_cos_vec <- tib %>% dplyr::pull(!!cos_key) %>% unique() %>% as.vector()
      
      tar_srd_vec <- expand.grid(tar_srd_vec, tar_cos_vec) %>% 
        tibble::as_tibble() %>% 
        tidyr::unite(srd, Var1,Var2, sep='', remove=TRUE)
      
      ids_sym <- rlang::sym(ids_key)
      seq_sym <- rlang::sym(seq_key)
      prb_sym <- rlang::sym(prb_key)
      
      if (purrr::is_character(tib))
        tib <- safe_read(tib, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      if (max>0) {
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{tabsStr} Will subset input to max={max}.{RET}"))
        tib <- tib %>% head(n=max)
      }
      
      src_man_tib <- tib %>% 
        dplyr::select(!!ids_sym, !!prb_sym, !!seq_sym) %>% 
        dplyr::mutate(Seq_ID:=!!ids_sym, PRB_DES:=!!prb_sym)
      #  dplyr::mutate(Seq_ID:=!!ids_sym, PRB_DES:=!!prb_sym)
      
      ret_key <- glue::glue("src_man_tib-1({funcTag})")
      ret_cnt <- print_tib(src_man_tib,funcTag, verbose,vt+6,tc, n=ret_key)
      
      # Ensure we have 122 mer format 60[NN]60
      src_man_tib <- validate_templates(tib=src_man_tib,
                                        seq_key=seq_key,
                                        prb_len=prb_len,
                                        seq_len=seq_len,
                                        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      if (FALSE) {
        if (verbose>=vt)
          cat(glue::glue("[{funcTag}]:{tabsStr} Validating design sequneces...{RET}"))
        
        src_man_tib <- src_man_tib %>%
          dplyr::mutate(
            !!seq_sym := stringr::str_replace(!!seq_sym, '\\[','_') %>% 
              stringr::str_replace('\\]','_')) %>%
          tidyr::separate(!!seq_sym, into=c("PreSeqN", "MidSeqN", "PosSeqN"), sep='_') %>%
          dplyr::mutate(
            PreSeqN=stringr::str_sub(PreSeqN,   -60),
            PosSeqN=stringr::str_sub(PosSeqN, 1, 60),
            PreSeqN=stringr::str_pad(string=PreSeqN, width=60, side='left', pad='N'),
            PosSeqN=stringr::str_pad(string=PosSeqN, width=60, side='right', pad='N'),
            DesNucA=stringr::str_sub(MidSeqN, 1,1), DesNucB=stringr::str_sub(MidSeqN, 2,2),
            !!seq_sym :=paste0(PreSeqN,'[',MidSeqN,']',PosSeqN) )
        
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{tabsStr} Done. Validating design sequneces.{RET}{RET}"))
        
        ret_key <- glue::glue("src_man_tib-2({funcTag})")
        ret_cnt <- print_tib(src_man_tib,funcTag, verbose,vt+6,tc, n=ret_key)
        
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #
      #          Build Ref/BSC Template Sequences for Target Strands::
      #                              bsc_templates()
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Building forward & reverse ",
                       "reference template sequneces...{RET}"))
      
      bsc_tibs <- NULL
      bsc_tibs <- c(
        bsc_templates(tib=src_man_tib, 
                      srd_str="F",
                      srd_vec=tar_srd_vec, 
                      seq_key=seq_key,
                      
                      srd_bol=TRUE,
                      srd_key=sr_str,  # Previously "SR"
                      cos_bol=TRUE,
                      cos_key=co_str,  # Previously "CO"
                      
                      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
        
        bsc_templates(tib=src_man_tib, 
                      srd_str="R",
                      srd_vec=tar_srd_vec, 
                      seq_key=seq_key, 
                      
                      srd_bol=FALSE,
                      srd_key=sr_str,  # Previously "SR"
                      cos_bol=TRUE,
                      cos_key=co_str,  # Previously "CO"
                      
                      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) )
      
      srd_names <- names(bsc_tibs)
      srd_count <- srd_names %>% length()
      
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Strand Names({srd_names})={RET}"))
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Done. Building forward & reverse ",
                       "reference template sequneces(strands={srd_count})!{RET}{RET}"))
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                    Build all Probes on Each Strand::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if (parallel) {
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Building probes for each ",
                         "strand (Parallel)...{RET}"))
        
        ret_tib <- foreach (srd=srd_names, .combine=rbind) %dopar% {
          lapply(split(bsc_tibs[[srd]], dplyr::pull(bsc_tibs[[srd]],prb_key)), 
                 des_all_prbs, srd_key=sr_str, cos_key=co_str, prb_key=prb_key,
                 verbose=verbose, vt=vt+4,tc=tc+1) %>% dplyr::bind_rows()
        }
      } else {
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Building probes for each ",
                         "strand (Linear)...{RET}"))
        
        for (srd in srd_names) {
          ret_tib <- ret_tib %>% dplyr::bind_rows(
            lapply(split(bsc_tibs[[srd]], dplyr::pull(bsc_tibs[[srd]],prb_key)), 
                   des_all_prbs, srd_key=sr_str, cos_key=co_str, prb_key=prb_key,
                   verbose=verbose, vt=vt+4,tc=tc+1) %>% 
              dplyr::bind_rows() )
        }
      }
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Building probes for ",
                       "each strand.{RET}{RET}"))
      ret_key <- glue::glue("prbs-on-all-srds:ret-tib({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n=ret_key)
      
      # Update Keys::
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Updating keys and strands.{RET}"))
      
      sr_out_str <- paste("Strand",sr_str, sep=del)
      co_out_str <- paste("Strand",co_str, sep=del)
      sr_out_sym <- rlang::sym(sr_out_str)
      co_out_sym <- rlang::sym(co_out_str)
      
      ret_tib <- ret_tib %>% 
        dplyr::mutate(
          !!sr_out_sym:=case_when(!!sr_sym ~ sr_vec[1], !(!!sr_sym) ~ sr_vec[2], TRUE ~ NA_character_),
          !!co_out_sym:=case_when(!!co_sym ~ co_vec[1], !(!!co_sym) ~ co_vec[2], TRUE ~ NA_character_),
          Seq_ID_Unq=paste(!!ids_sym,paste0(!!sr_out_sym,!!co_out_sym), sep=del)
          # Seq_ID_Unq=paste(Seq_ID,paste0(!!sr_out_sym,!!co_out_sym), sep=del)
        ) %>% dplyr::arrange(Seq_ID_Unq)
      
      ret_key <- glue::glue("Updated-keys/srds:ret-tib({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n=ret_key)
      
      # Add match probe sequences
      if (add_matseq) {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding match seqqunes.{RET}"))
        ret_tib <- ret_tib %>% 
          dplyr::mutate(PRB1_U_MAT=stringr::str_to_upper(PRB1_U),
                        PRB1_M_MAT=stringr::str_to_upper(PRB1_M),
                        PRB2_D_MAT=stringr::str_to_upper(PRB2_D) )
      }
      
      ret_key <- glue::glue("ret-FIN({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n=ret_key)
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    return(ret_tib)
  }
  ret_key <- glue::glue("ret-FIN({funcTag})")
  ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n=ret_key)
  
  etime <- 0
  # etime <- stime[3] %>% as.double() %>% round(2)
  # if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

validate_templates = function(tib,
                              seq_key,
                              prb_len=60, 
                              seq_len=122, 
                              del="_",
                              pad="N",
                              verbose=0,vt=3,tc=1,tt=NULL,
                              funcTag='validate_templates') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   seq_key={seq_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   prb_len={prb_len}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   seq_len={seq_len}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    seq_sym <- rlang::sym(seq_key)
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Validating design sequneces...{RET}"))
    
    ret_tib <- tib %>%
      dplyr::mutate(
        !!seq_sym := stringr::str_replace(!!seq_sym, '\\[','_') %>% 
          stringr::str_replace('\\]','_')) %>%
      tidyr::separate(!!seq_sym, 
                      into=c("PreSeqN", "MidSeqN", "PosSeqN"), sep=del) %>%
      dplyr::mutate(
        PreSeqN=stringr::str_sub(PreSeqN,   -prb_len),
        PosSeqN=stringr::str_sub(PosSeqN, 1, prb_len),
        PreSeqN=stringr::str_pad(string=PreSeqN, 
                                 width=prb_len, side='left', pad=pad),
        PosSeqN=stringr::str_pad(string=PosSeqN, 
                                 width=prb_len, side='right', pad=pad),
        DesNucA=stringr::str_sub(MidSeqN, 1,1), 
        DesNucB=stringr::str_sub(MidSeqN, 2,2),
        !!seq_sym :=paste0(PreSeqN,'[',MidSeqN,']',PosSeqN) )
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Done. Validating design ",
                     "sequneces.{RET}{RET}"))
    
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

bsc_templates = function(tib,
                         srd_str,
                         srd_vec,
                         seq_key,
                         
                         srd_bol,
                         srd_key="SR",
                         cos_bol,
                         cos_key="CO",
                         
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='bsc_templates') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   srd_str={srd_str}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   srd_vec={srd_vec}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   seq_key={seq_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   srd_bol={srd_bol}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   srd_key={srd_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   cos_bol={cos_bol}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   cos_key={cos_key}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  stime <- base::system.time({
    
    srdC <- paste0(srd_str,"C")
    srdO <- paste0(srd_str,"O")
    
    seq_sym <- rlang::sym(seq_key)
    srd_sym <- rlang::sym(srd_key)
    cos_sym <- rlang::sym(cos_key)
    
    des_seq_tib <- NULL
    if (srdC %in% srd_vec || srdO %in% srd_vec)
      des_seq_tib <- tib %>% dplyr::mutate(
        !!srd_sym:=srd_bol,
        !!cos_sym:=cos_bol, 
        DesSeqN=shearBrac(!!seq_sym) )
    ret_key <- glue::glue("des_seq_tib({funcTag})")
    ret_cnt <- print_tib(des_seq_tib,funcTag, verbose,vt+6,tc, n=ret_key)
    
    # BSC-Forward-Converted::
    #
    if (!is.null(des_seq_tib)) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Building bisulfite converted ",
                       "forward strand template sequence(s)...{RET}"))
      
      ret_dat[[srdC]] <- des_seq_tib %>% 
        dplyr::mutate(
          DesBscU = bscUs(DesSeqN),
          DesBscM = bscMs(DesSeqN),
          DesBscD = bscDs(DesSeqN) )
      
      if (verbose>=vt+1) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Done. Building bisulfite ",
                       "template sequences ({srdC}).{RET}"))
      ret_key <- glue::glue("bsc_tibs-{srdC}({funcTag})")
      ret_cnt <- print_tib(ret_dat[[srdC]],funcTag, verbose,vt+6,tc, n=ret_key)
      
      # BSC-Forward-Opposite::
      #
      if (srdO %in% srd_vec) {
        ret_dat[[srdO]] <- ret_dat[[srdC]] %>% 
          dplyr::mutate(
            !!srd_sym:=srd_bol,
            !!cos_sym:=!cos_bol,
            DesBscU=revCmp(DesBscU),
            DesBscM=revCmp(DesBscM),
            DesBscD=revCmp(DesBscD) )
        
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr} Done. Building bisulfite ",
                         "template sequences ({srdO}).{RET}"))
        ret_key <- glue::glue("bsc_tibs-{srdO}({funcTag})")
        ret_cnt <- print_tib(ret_dat[[srdO]],funcTag, verbose,vt+6,tc,n=ret_key)
      }
      if (!srdC %in% srd_vec) ret_dat[[srdC]] <- NULL
      if (verbose>=vt) cat(glue::glue("{RET}"))
    }
    
    ret_key <- glue::glue("des_seq_tib-FIN({funcTag})")
    ret_cnt <- print_tib(des_seq_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_dat
}

# Function to design all probes in a single call
des_all_prbs = function(tib,
                        srd_key="SR",
                        cos_key="CO",
                        prb_key,
                        verbose=0,vt=5,tc=1,tt=NULL, funcTag='des_all_prbs') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   srd_key={srd_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   cos_key={cos_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   prb_key={prb_key}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_tib <- NULL
  stime <- system.time({
    
    srd_sym <- rlang::sym(srd_key)
    cos_sym <- rlang::sym(cos_key)
    prb_sym <- rlang::sym(prb_key)
    
    fr <- tib %>% dplyr::distinct(!!srd_sym) %>% base::as.logical()
    co <- tib %>% dplyr::distinct(!!cos_sym) %>% base::as.logical()
    pr <- tib %>% dplyr::distinct(!!prb_sym) %>% base::as.character()
    
    ret_tib <- dplyr::bind_rows(
      tib %>% 
        des_to_prbs(fwd=fr, con=co, pr=pr, mu='N', des_seq='DesSeqN', 
                    verbose=verbose,vt=vt+1,tc=tc+1) %>%
        des_to_prbs(fwd=fr, con=co, pr=pr, mu='U', des_seq='DesBscU', 
                    verbose=verbose,vt=vt+1,tc=tc+1) %>%
        des_to_prbs(fwd=fr, con=co, pr=pr, mu='M', des_seq='DesBscM', 
                    verbose=verbose,vt=vt+1,tc=tc+1) %>%
        des_to_prbs(fwd=fr, con=co, pr=pr, mu='D', des_seq='DesBscD', 
                    verbose=verbose,vt=vt+1,tc=tc+1) )
    
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

# Function to design probes from single design strand and orientations::
#   TBD:: Previous default was 'QC_CPN=TRUE' Not sure if that is needed...
#
des_to_prbs = function(tib, 
                       fwd, 
                       con, 
                       pr, 
                       mu, 
                       des_seq='DesSeqN', 
                       len=48, 
                       del='_',
                       QC_CPN=FALSE,
                       verbose=0,vt=5,tc=1,tt=NULL, funcTag='des_to_prbs') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       pr={pr}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       mu={mu}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}      fwd={fwd}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}      con={con}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  des_seq={des_seq}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    stopifnot(is.logical(fwd))
    stopifnot(is.logical(con))
    
    if (!is.logical(fwd)) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: fwd={fwd} Must be ",
                      "logical!{RET}{RET}"))
      return(ret_tib)
    }
    if (!is.logical(con)) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: con={con} Must be ",
                      "logical!{RET}{RET}"))
      return(ret_tib)
    }
    
    if (mu!='N' && mu!='U' && mu!='M' && mu!='D') {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: mu={mu} Only ",
                      "Supported=[N,U,M,D]!{RET}{RET}"))
      return(ret_tib)
    }
    if (pr!='cg' && pr!='ch' && pr!='rs' && pr!='rp' && pr!='mu' && pr!='bc') {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: pr={pr} Only ",
                      "Supported=[cg,ch,rp,rs]!{RET}{RET}"))
      return(ret_tib)
    }
    
    des_sym <- rlang::sym(des_seq)
    if (pr=='rs') {
      if      ( fwd &&  con) nxb_pos <- 60
      else if (!fwd &&  con) nxb_pos <- 61
      else if ( fwd && !con) nxb_pos <- 61
      else if (!fwd && !con) nxb_pos <- 60
      else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: unsupported combination ",
                        "fwd={fwd}, con={con}!{RET}{RET}"))
        return(ret_tib)
      }
    } else if (pr=='ch') {
      # Originally this was identical to rs format above, but for forward 
      #  sequences needs to be shifted upstream for converted and downstream 
      #  for opposite::
      if      ( fwd &&  con) nxb_pos <- 61 # Previously = 60
      else if (!fwd &&  con) nxb_pos <- 61
      else if ( fwd && !con) nxb_pos <- 60 # Previously = 61
      else if (!fwd && !con) nxb_pos <- 60
      else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: unsupported combination ",
                        "fwd={fwd}, con={con}!{RET}{RET}"))
      }
      
      # NEW:: CpH Code::
      #
      if      ( fwd &&  con) nxb_pos <- 60 # Previously = 60
      else if (!fwd &&  con) nxb_pos <- 61
      else if ( fwd && !con) nxb_pos <- 61 # Previously = 61
      else if (!fwd && !con) nxb_pos <- 60
      else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: unsupported combination ",
                        "fwd={fwd}, con={con}!{RET}{RET}"))
        return(ret_tib)
      }
      
    } else if (pr=='cg' || pr=='rp' || pr=='mu' || pr=='bc' || 
               stringr::str_starts(pr,'ct')) {
      # $prb_F_C_I  = revCmp(substr($des_F_C, 60, 50));
      # $prb_F_O_I  = revCmp(substr($des_F_O, 61, 50));
      # $prb_R_C_I  = revCmp(substr($des_R_C, 60, 50));
      # $prb_R_O_I  = revCmp(substr($des_R_O, 61, 50));
      # 
      # $prb_F_C_II  = revCmp(substr($des_F_C, 61, 50));
      # $prb_F_O_II  = revCmp(substr($des_F_O, 62, 50));
      # $prb_R_C_II  = revCmp(substr($des_R_C, 61, 50));
      # $prb_R_O_II  = revCmp(substr($des_R_O, 62, 50));
      nxb_pos <- 60
      if (!con) nxb_pos <- 61
      
    } else {
      stop(glue::glue("[{funcTag}]: ERROR: Probe_Type={pr} is currently NOT ",
                      "supported!{RET}{RET}"))
      return(ret_tib)
    }
    cpg_pos <- nxb_pos + 1
    sec_pos <- cpg_pos + 1
    bod_pos <- sec_pos + 1
    end_pos <- bod_pos + len
    
    # Special consideration is needed for U/M strands at the query site. 
    #  For CN (i.e. cg or ch) this is actually done naturally in U/M conversion
    #  However, for  non-CN probes (i.e. rs) this needs to be forced to U/M
    #
    # This is handled by the TAR (Target/Query Nucleotide). This should only 
    #  change for U/M (QMAP_U/QMAP_M) for D its just itself.
    #
    ret_tib <- tib %>% dplyr::mutate(
      NXB=stringr::str_sub(!!des_sym, nxb_pos, nxb_pos),
      CPN=stringr::str_sub(!!des_sym, cpg_pos, cpg_pos),
      TAR=qmaps(CPN, mu=mu),
      SEC=stringr::str_sub(!!des_sym, sec_pos, sec_pos),
      BOD=stringr::str_sub(!!des_sym, bod_pos, end_pos-1),
      END=stringr::str_sub(!!des_sym, end_pos, end_pos)
    )
    
    #  QC TEST:: for CpN (cg or ch) verify that the probes are equal. Well call
    #   this PRB0 (CGN) and PRB1 (TAR). After testing remove PRB0
    #
    if (QC_CPN && (pr=='cg')) {
      ret_tib <- ret_tib %>%
        tidyr::unite(PRB0, CPN,SEC,BOD, sep='', remove=FALSE) %>%
        tidyr::unite(PRB1, TAR,SEC,BOD, sep='', remove=FALSE) %>%
        tidyr::unite(PRB2, SEC,BOD,END, sep='', remove=FALSE) %>%
        dplyr::mutate(PRB0=revCmp(PRB0), PRB1=revCmp(PRB1), PRB2=revCmp(PRB2))
      
      qc_tib <- ret_tib %>% filter(PRB0!=PRB1)
      qc_len <- qc_tib %>% base::nrow()
      if (qc_len != 0) {
        qc_tib %>% dplyr::select(1,PRB0,PRB1) %>% print()
        stop(glue::glue("{RET}[{funcTag}]: ERROR: pr={pr}, qc_len={qc_len} ",
                        "!= 0!{RET}{RET}"))
        return(NULL)
      }
    } else {
      ret_tib <- ret_tib %>%
        tidyr::unite(PRB1, TAR,SEC,BOD, sep='', remove=FALSE) %>%
        tidyr::unite(PRB2, SEC,BOD,END, sep='', remove=FALSE) %>%
        dplyr::mutate(PRB1=revCmp(PRB1), PRB2=revCmp(PRB2))
    }
    
    # Add suffix to sequences for merging later
    ret_tib <- ret_tib %>%
      dplyr::select(PRB1,PRB2, NXB,CPN,TAR,BOD,END, everything()) %>%
      dplyr::rename(!!paste('PRB1',mu, sep=del):=PRB1,
                    !!paste('PRB2',mu, sep=del):=PRB2,
                    !!paste('NXB', mu, sep=del):=NXB,
                    !!paste('CPN', mu, sep=del):=CPN,
                    !!paste('TAR', mu, sep=del):=TAR,
                    !!paste('SEC', mu, sep=del):=SEC,
                    !!paste('BOD', mu, sep=del):=BOD,
                    !!paste('END', mu, sep=del):=END)
    
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

des2prbsNOTES = function(srd, desSeq,
                         verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'des2prbs'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  # bscNewU <- bscNewU(desSeq)
  # bscNewM <- bscNewM(desSeq)
  # bscNewD <- bscNewD(desSeq)
  
  # my @desSetU = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewU);
  # my @desSetM = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewM);
  # my @desSetD = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewD);
  
  #    my @desSetU=unpack("A11"."A48"."AAAA"."A47"."A12",bscNewU($desSeq,$retUC));                                                                                                             
  #    my @desSetM=unpack("A11"."A48"."AAAA"."A47"."A12",bscNewM($desSeq,$retUC));                                                                                                             
  
  # my ($desNxbU, $desCpgU, $desSeqU, $desEndU);
  # my ($desNxbM, $desCpgM, $desSeqM, $desEndM);
  # my ($desNxbD, $desCpgD, $desSeqD, $desEndD);
  
  #         (  $desNxbU,    $desCpgU,               $desSeqU,                $desEndU) =                                                                                                       
  #   return( $desSetU[2], $desSetU[3], $desSetU[4].$desSetU[5].$desSetU[6], $desSetU[7],                                                                                                      
  #           $desSetM[2], $desSetM[3], $desSetM[4].$desSetM[5].$desSetM[6], $desSetM[7]) if ($desCO eq $C);                                                                                   
  #          ( $desNxbM,    $desCpgM,               $desSeqM,                $desEndM) =                                                                                                       
  
  #                  ( $desNxbU,            $desCpgU,                    $desSeqU,                  $desEndU) =                                                                                
  #   return( revCmpl($desSetU[4]), revCmpl($desSetU[3]), revCmpl($desSetU[1].$desSetU[2]), revCmpl($desSetU[0]),                                                                              
  #           revCmpl($desSetM[4]), revCmpl($desSetM[3]), revCmpl($desSetM[1].$desSetM[2]), revCmpl($desSetM[0])) if ($desCO eq $O);                                                           
  #                  ( $desNxbM,            $desCpgM,                    $desSeqM,                  $desEndM) =                                                                                
  
  
  # $$prbRef[$srd][$iU]=[ $desSetU[2], $desSetU[3], $desSetU[4].$desSetU[5].$desSetU[6], $desSetU[7] ];
  # $$prbRef[$srd][$iM]=[ $desSetM[2], $desSetM[3], $desSetM[4].$desSetM[5].$desSetM[6], $desSetM[7] ];
  # $$prbRef[$srd][$iD]=[ $desSetD[2], $desSetD[3], $desSetD[4].$desSetD[5].$desSetD[6], $desSetD[7] ];
  # 
  # $srd++;
  # $$prbRef[$srd][$iU]=[ revCmpl($desSetU[4]), revCmpl($desSetU[3]), revCmpl($desSetU[1].$desSetU[2]), revCmpl($desSetU[0]) ];
  # $$prbRef[$srd][$iM]=[ revCmpl($desSetM[4]), revCmpl($desSetM[3]), revCmpl($desSetM[1].$desSetM[2]), revCmpl($desSetM[0]) ];
  # $$prbRef[$srd][$iD]=[ revCmpl($desSetD[4]), revCmpl($desSetD[3]), revCmpl($desSetD[1].$desSetD[2]), revCmpl($desSetD[0]) ];
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  NULL
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Run Both improbe and r-improbe Simultaneously!
#
#  NOTE:: This code hasn't been tested much and has more or less been 
#   re-implemented in seperate workflows, but keeping it around for now...
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  improbe_design_all = function(tib, ptype, outDir, gen, image, shell,
                                seqKey="IUPAC_Sequence", sr_str="FR",
                                reduce_imp=TRUE, parse_din=FALSE,
                                parallel=TRUE, retData=FALSE,
                                sidx=2, plen=50,del='_',
                                verbose=0,vt=3,tc=1,tt=NULL,
                                funcTag='improbe_design_all') {
    
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Starting; ptype={ptype}...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    ret_dat <- NULL
    stime <- system.time({
      idx1 <- sidx
      len1 <- plen - 1
      idx2 <- sidx + 1
      len2 <- plen
      
      if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
      
      imp_fwd_name <- paste(ptype,gen,"improbe_fwd-seq.tsv.gz", sep=del)
      imp_fwd_path <- file.path(outDir, imp_fwd_name)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                                 improbe::
      #                                docker/c++
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      imp_fwd_tib <- tib %>% 
        dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island)
      readr::write_tsv(imp_fwd_tib, imp_fwd_path)
      
      imp_des_tib <- run_improbe_docker(
        dir=outDir, file=imp_fwd_name, 
        name=gen, image=image, shell=shell, 
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        dplyr::mutate(
          Seq_48U_1=stringr::str_sub(UnMethyl_Probe_Sequence, idx1,len1) %>% 
            stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% 
            stringr::str_replace_all('Y', 'T'),
          
          Seq_48U_2=stringr::str_sub(UnMethyl_Probe_Sequence, idx2,len2) %>% 
            stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% 
            stringr::str_replace_all('Y', 'T'),
          
          Chromosome=as.character(Chromosome),
          Strand_SR=Methyl_Allele_FR_Strand,
          Strand_TB=stringr::str_sub(Methyl_Allele_TB_Strand,1,1),
          Strand_CO=Methyl_Allele_CO_Strand,
          Probe_Score_Min=pmin(Methyl_Final_Score,UnMethyl_Final_Score),
          Underlying_CpG_Count=as.integer(Methyl_Underlying_CpG_Count),
          Underlying_CpG_Min_Dist=as.integer(Methyl_Underlying_CpG_Min_Dist)
        )
      
      if (parse_din) imp_des_tib <- imp_des_tib %>%
        tidyr::separate(Seq_ID, into=c("Seq_ID","DiNuc"), sep=del)
      
      #
      # TBD:: Add better reduced returned value...
      #
      if (reduce_imp) {
        imp_des_tib <- imp_des_tib %>%
          dplyr::select(
            Seq_ID,Forward_Sequence,Genome_Build,Chromosome,Coordinate,Top_Sequence,
            Strand_SR,Strand_TB,Strand_CO,Probe_Score_Min,
            Underlying_CpG_Count,Underlying_CpG_Min_Dist,
            Methyl_Probe_Sequence,UnMethyl_Probe_Sequence,Seq_48U_1,Seq_48U_2)
      }
      imp_key <- glue::glue("imp_des_tib({funcTag})")
      imp_cnt <- print_tib(imp_des_tib,funcTag, verbose,vt+4,tc, n=imp_key)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                               r-improbe::
      #                           de-novo IUPAC designs
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      iup_des_tib <- r_improbe(
        tib=tib,
        ids_key="Seq_ID",seqKey=seqKey,prbKey="Probe_Type",
        sr_str=sr_str, parallel=parallel,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        dplyr::mutate(
          
          Seq_48U_1=stringr::str_sub(PRB1_U_MAT, idx1,len1) %>% 
            stringr::str_to_upper() %>% 
            stringr::str_replace_all('R', 'A') %>% 
            stringr::str_replace_all('Y', 'T'),
          
          Seq_48U_2=stringr::str_sub(PRB1_U_MAT, idx2,len2) %>% 
            stringr::str_to_upper() %>% 
            stringr::str_replace_all('R', 'A') %>% 
            stringr::str_replace_all('Y', 'T'),
        )
      
      if (parse_din) iup_des_tib <- iup_des_tib %>%
        tidyr::separate(Seq_ID, into=c("Seq_ID","DiNuc"), sep=del)
      
      if (!reduce_imp) imp_des_tib <- imp_des_tib %>% 
        dplyr::rename(Probe_Type_IMP=Probe_Type)
      
      ret_tib <- dplyr::inner_join(
        imp_des_tib, iup_des_tib,
        by=c("Seq_ID", "Strand_SR", "Strand_CO",
             "Seq_48U_1", "Seq_48U_2"),
        suffix=c("_IMP", "_IUP")
      )
      
      if (retData) ret_dat$imp <- imp_des_tib
      if (retData) ret_dat$iup <- iup_des_tib
      
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
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#          Infinium Methylation Probe toString/printing Methodss::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Replace Seq_ID and PRB_DES with symbolic links
#  Seq_ID::  seq_sym -> rlang::sym(seq_key)
#  PRB_DES:: des_sym -> rlang::sym(des_key)
#
print_prbs = function(tib, 
                      pr='cg', 
                      org=NULL, 
                      outDir, 
                      plotName, 
                      max=0,
                      verbose=0,vt=5,tc=1,tt=NULL, funcTag='print_prbs') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:{tabsStr}       pr={pr}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   outDir={outDir}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} plotName={plotName}.{RET}{RET}"))
  }
  
  prb_mat_tibs <- tib %>% dplyr::filter(PRB_DES==pr) %>% dplyr::distinct()
  prb_mat_cnt  <- prb_mat_tibs %>% base::nrow()
  
  plot_ord_tib <- prb_mat_tibs %>% dplyr::distinct(Seq_ID, .keep_all=TRUE)
  plot_ord_cnt <- plot_ord_tib %>% base::nrow()
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} plot_ord_cnt={plot_ord_cnt}, ",
                   "prb_mat_cnt={prb_mat_cnt}.{RET}"))
  
  if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
  
  passFile <- file.path(outDir, paste(plotName,pr,'pass.txt', sep='_') )
  failFile <- file.path(outDir, paste(plotName,pr,'fail.txt', sep='_') )
  unlink(passFile)
  unlink(failFile)
  
  if (!is.null(org)) org <- org %>% dplyr::distinct()
  
  tibs <- NULL
  for (ii in seq(1,plot_ord_cnt)) {
    tag_tib <- plot_ord_tib[ii,]
    cur_SeqID <- tag_tib$Seq_ID[1]
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ii={ii}; cur_SeqID={cur_SeqID}.{RET}"))
    
    cur_org <- NULL
    org_str <- ''
    if (!is.null(org)) {
      cur_org <- org %>% dplyr::filter(Seq_ID==cur_SeqID) %>% head(n=1)
      
      if (FALSE) {
        for (jj in c(1:length(cur_org))) {
          org_str[jj] <- stringr::str_c(cur_org[jj], collapse=',')
          # org_str[jj] <- paste(cur_org[jj], collapse=', ')
          break
        }
      }
      org_str <- cur_org %>% paste(collapse=', ')
    }
    
    cur_tib <- prb_mat_tibs %>% dplyr::filter(Seq_ID==cur_SeqID)
    tibs <- cur_tib %>% srds_to_brac()
    
    # Build Printable Strings for reach strand::
    strs <- NULL
    strs$F <- prbsToStr(tibs$F, pr=tag_tib$PRB_DES[1], verbose=verbose, vt=6)
    strs$R <- prbsToStr(tibs$R, pr=tag_tib$PRB_DES[1], verbose=verbose, vt=6)
    
    # Determine which file to write to::
    #   - Failed if both F/R strands have equal PRB1_U==PRB1_M
    #   - Passed otherwise
    fwd_cnt <- tibs$F %>% dplyr::filter(PRB1_U!=PRB1_M) %>% base::nrow()
    rev_cnt <- tibs$R %>% dplyr::filter(PRB1_U!=PRB1_M) %>% base::nrow()
    
    space_str <- paste(stringr::str_pad("# ", width=150, side="right", pad="-"), "\n", sep='')
    seq_id_str <- ''
    for (jj in c(1:length(org_str))) {
      
      # seq_id_str <- paste0("Probe_Type=",pr,"; Seq_ID: ", 
      #                      paste(unique(cur_tib$Seq_ID), collapse='\t'),
      #                      "; Original: ",org_str[jj],"\n")
      
      seq_id_str <- 
        glue::glue("Probe_Type={pr}: Seq_ID: ",
                   paste(unique(cur_tib$Seq_ID), collapse='\t'),
                   "; Original: {org_str[jj]}{RET}")
    }
    out_lines <- NULL
    out_lines[1] <- space_str
    out_lines[2] <- seq_id_str
    out_lines[3] <- stringr::str_c(strs$F)
    out_lines[4] <- stringr::str_c(strs$R)
    
    if (verbose>=vt) cat(out_lines)
    
    if (fwd_cnt>0 || rev_cnt>0) {
      readr::write_lines(out_lines, passFile, append=TRUE)
    } else {
      readr::write_lines(out_lines, failFile, append=TRUE)
    }
    
    # This was done for known comparisons::
    #
    if (FALSE) {
      if (FALSE) { # && cur_tib$Infinium_Design_Type[1]=='II') {
        cur_tib %>% dplyr::select(FR,CO, PRB2_D_IUP, PRB2_D_IMP,Man_MisMatch,Man_TarMatch,Bad_Design) %>% print()
      } else {
        prb1U_iup <- cur_tib %>% dplyr::select(FR,CO, PRB1_U_IUP, PRB1_U_IMP,Man_MisMatch,Man_TarMatch,Bad_Design)
        prb1M_iup <- cur_tib %>% dplyr::select(FR,CO, PRB1_M_IUP, PRB1_M_IMP,Man_MisMatch,Man_TarMatch,Bad_Design)
        # readr::write_file(str_c(prb1U_iup), passFile, append=TRUE)
        # readr::write_file(str_c(prb1M_iup), passFile, append=TRUE)
      }
    }
    
    if (max!=0 && ii>=max) break
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  tibs
}

srds_to_brac = function(tib, 
                        beg1=1, end1=60, mid1=61,
                        beg2=63,end2=122,mid2=62) {
  # TBD:: Calculate all data points based on sequence length
  
  tib %>% 
    dplyr::mutate(StrandFR=case_when(FR ~ 'F', !FR ~ 'R', TRUE ~ NA_character_),
                  StrandCO=case_when(CO ~ 'C', !CO ~ 'O', TRUE ~ NA_character_),
                  DesSeqN=paste0(stringr::str_sub(DesSeqN,beg1,end1),{BNG},
                                 stringr::str_sub(DesSeqN,mid1,mid2),{BNG},
                                 stringr::str_sub(DesSeqN,beg2,end2)),
                  
                  DesBscU=paste0(stringr::str_sub(DesBscU,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscU,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscU,beg2,end2)),
                  
                  DesBscM=paste0(stringr::str_sub(DesBscM,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscM,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscM,beg2,end2)),
                  
                  DesBscD=paste0(stringr::str_sub(DesBscD,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscD,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscD,beg2,end2))) %>%
    dplyr::arrange(StrandFR, StrandCO) %>% split(.$StrandFR)
  
}

prbsToStr = function(tib, pr, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'prbsToStr'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; pr={pr}.{RET}"))
  
  # fr1Key <- tib %>% dplyr::select(!!frKey) %>% head(n=1) %>% pull()
  # fr2Key <- tib$StrandFR[2]
  # co1Key <- tib$StrandCO[1]
  # co2Key <- tib$StrandCO[2]
  
  fr1Key <- tib$StrandFR[1]
  fr2Key <- tib$StrandFR[2]
  co1Key <- tib$StrandCO[1]
  co2Key <- tib$StrandCO[2]
  
  if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} fr1Key={fr1Key}, fr2Key={fr2Key}, co1Key={co1Key}, co2Key={co2Key}.{RET}"))
  
  mud <- list('U'='U','M'='M','D'='D')
  
  desSeq <- 'DesSeqN'
  bscKey <- lapply(mud, function(x) { paste('DesBsc',x,sep='')} )
  nxbKey <- lapply(mud, function(x) { paste('NXB',x,sep='_')} )
  cpnKey <- lapply(mud, function(x) { paste('CPN',x,sep='_')} )
  tarKey <- lapply(mud, function(x) { paste('TAR',x,sep='_')} )
  secKey <- lapply(mud, function(x) { paste('SEC',x,sep='_')} )
  bodKey <- lapply(mud, function(x) { paste('BOD',x,sep='_')} )
  endKey <- lapply(mud, function(x) { paste('END',x,sep='_')} )
  
  # Dertermine if the probes on this strand were designable::
  pass_cnt <- tib %>% dplyr::filter(PRB1_U!=PRB1_M) %>% base::nrow()
  pass_str <- 'PASS: '
  if (pass_cnt==0) pass_str <- 'FAIL: '
  
  # TBD:: Note on the Opposite Strand we should reverse all framents, but currently fragLen==1 are left alone for effiecntcy...
  # Sketch Output::
  #
  # F_C_N    DesSeqN[CG]DesSeqN
  #
  #                  D2 22222
  #                N M1 1111
  #                N U1 1111
  # F_C_U    DesBscU[tG]DesBscU
  # F_O_U    DesBscU[Ca]DesBscU
  #             1111 1U N
  #             1111 1M N
  #            22222 2D
  #
  if (pr=='rs'||pr=='ch') {
    if (fr1Key=='F' && fr2Key=='F') {
      bufC <- 0
      bufO <- 0
      str <- glue::glue(
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 61+bufC), collapse=''),"{tib[[tarKey$D]][1]}{tib[[secKey$D]][1]}{BNG}{tib[[bodKey$D]][1]}{tib[[endKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 59+bufC), collapse=''),"{tib[[nxbKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[secKey$M]][1]}{BNG}{tib[[bodKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 59+bufC), collapse=''),"{tib[[nxbKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[secKey$U]][1]}{BNG}{tib[[bodKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_U {TAB}{tib[[bscKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_M {TAB}{tib[[bscKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{tib[[bscKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{tib[[desSeq]][1]}{RET}",
        # "FwdSeq{TAB}{tib$Forward_Sequence[1]}{RET}",
        "{RET}",
        "{pass_str}{fr1Key}_{co2Key}_N {TAB}{cmpl(tib[[desSeq]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_M {TAB}{Biostrings::reverse(tib[[bscKey$M]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_U {TAB}{Biostrings::reverse(tib[[bscKey$U]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 11-bufO), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][2])}{tib[[secKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[nxbKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 11-bufO), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][2])}{tib[[secKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[nxbKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 10-bufO), collapse=''),"{tib[[endKey$D]][2]}{Biostrings::reverse(tib[[bodKey$D]][2])}{tib[[secKey$D]][2]}{BNG}{tib[[tarKey$D]][2]}{RET}",
        "{RET}")
    } else if (fr1Key=='R' && fr2Key=='R') {
      bufO <- 0
      bufC <- 0
      str <- glue::glue(
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 61+bufO), collapse=''),"{tib[[tarKey$D]][2]}{tib[[secKey$D]][2]}{BNG}{tib[[bodKey$D]][2]}{tib[[endKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 59+bufO), collapse=''),"{tib[[nxbKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[secKey$M]][2]}{BNG}{tib[[bodKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 59+bufO), collapse=''),"{tib[[nxbKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[secKey$U]][2]}{BNG}{tib[[bodKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{tib[[bscKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_N {TAB}{revCmp(tib[[desSeq]][2])}{RET}",
        # "{fr2Key}_{co2Key}_N {TAB}{tib[[desSeq]][2]}{RET}",
        "{RET}",
        # "{fr1Key}_{co1Key}_N {TAB}{cmpl(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{Biostrings::reverse(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 11-bufC), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][1])}{tib[[secKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[nxbKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 11-bufC), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][1])}{tib[[secKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[nxbKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 10-bufC), collapse=''),"{tib[[endKey$D]][1]}{Biostrings::reverse(tib[[bodKey$D]][1])}{tib[[secKey$D]][1]}{BNG}{tib[[tarKey$D]][1]}{RET}",
        "{RET}")
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: fr1Key={fr1Key}, fr2Key={fr2Key}, Allowed Values=[F,R]!{RET}{RET}"))
    }
  } else {
    if (fr1Key=='F' && fr2Key=='F') {
      buf <- 0
      str <- glue::glue(
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 61+buf), collapse=''),"{tib[[tarKey$D]][1]}{tib[[secKey$D]][1]}{BNG}{tib[[bodKey$D]][1]}{tib[[endKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 59+buf), collapse=''),"{tib[[nxbKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[secKey$M]][1]}{BNG}{tib[[bodKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 59+buf), collapse=''),"{tib[[nxbKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[secKey$U]][1]}{BNG}{tib[[bodKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{tib[[bscKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{tib[[desSeq]][1]}{RET}",
        "{RET}",
        "{pass_str}{fr1Key}_{co2Key}_N {TAB}{cmpl(tib[[desSeq]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 11-buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][2])}{tib[[secKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[nxbKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 11-buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][2])}{tib[[secKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[nxbKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 10-buf), collapse=''),"{tib[[endKey$D]][2]}{Biostrings::reverse(tib[[bodKey$D]][2])}{tib[[secKey$D]][2]}{BNG}{tib[[tarKey$D]][2]}{RET}",
        "{RET}")
    } else if (fr1Key=='R' && fr2Key=='R') {
      buf <- 0
      buf <- 1
      str <- glue::glue(
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 61+buf), collapse=''),"{tib[[tarKey$D]][2]}{BNG}{tib[[secKey$D]][2]}{tib[[bodKey$D]][2]}{tib[[endKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 60+buf), collapse=''),"{tib[[nxbKey$M]][2]}{tib[[tarKey$M]][2]}{BNG}{tib[[secKey$M]][2]}{tib[[bodKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 60+buf), collapse=''),"{tib[[nxbKey$U]][2]}{tib[[tarKey$U]][2]}{BNG}{tib[[secKey$U]][2]}{tib[[bodKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{tib[[bscKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_N {TAB}{revCmp(tib[[desSeq]][2])}{RET}",
        # "{fr2Key}_{co2Key}_N {TAB}{tib[[desSeq]][2]}{RET}",
        "{RET}",
        # "{fr1Key}_{co1Key}_N {TAB}{cmpl(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{Biostrings::reverse(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 11+buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][1])}{BNG}{tib[[secKey$U]][1]}{tib[[tarKey$U]][1]}{BNG}{tib[[nxbKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 11+buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][1])}{BNG}{tib[[secKey$M]][1]}{tib[[tarKey$M]][1]}{BNG}{tib[[nxbKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 10+buf), collapse=''),"{tib[[endKey$D]][1]}{Biostrings::reverse(tib[[bodKey$D]][1])}{BNG}{tib[[secKey$D]][1]}{tib[[tarKey$D]][1]}{RET}",
        "{RET}")
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: fr1Key={fr1Key}, fr2Key={fr2Key}, Allowed Values=[F,R]!{RET}{RET}"))
    }
  }
  if (verbose>=vt) cat(str)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  str
}

prbsToStrMUD = function(tib, mu='U', verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'prbsToStr'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  fr1Key <- tib$StrandFR[1]
  fr2Key <- tib$StrandFR[2]
  
  co1Key <- tib$StrandCO[1]
  co2Key <- tib$StrandCO[2]
  
  bscKey <- paste('DesBsc',mu, sep='')
  
  nxbKey <- paste('NXB',mu, sep='_')
  cpnKey <- paste('CPN',mu, sep='_')
  tarKey <- paste('TAR',mu, sep='_')
  secKey <- paste('SEC',mu, sep='_')
  bodKey <- paste('BOD',mu, sep='_')
  endKey <- paste('END',mu, sep='_')
  
  str <- glue::glue(
    paste0(rep(" ", 59+10), collapse=''),"{tib[[cpnKey]][1]}{tib[[secKey]][1]}{BNG}{tib[[bodKey]][1]}{tib[[endKey]][1]}{RET}",
    paste0(rep(" ", 59+8),  collapse=''),"{tib[[nxbKey]][1]}{BNG}{tib[[cpnKey]][1]}{tib[[secKey]][1]}{BNG}{tib[[bodKey]][1]}{RET}",
    "{fr1Key}_{co1Key}_{mu}{TAB}{tib[[bscKey]][1]}{RET}",
    "{fr2Key}_{co2Key}_{mu}{TAB}{Biostrings::reverse(tib[[bscKey]][2])}{RET}",
    paste0(rep(" ", 9+10), collapse=''),"{Biostrings::reverse(tib[[bodKey]][2])}{tib[[secKey]][2]}{BNG}{tib[[cpnKey]][2]}{tib[[nxbKey]][2]}{RET}",
    paste0(rep(" ", 9+9),  collapse=''),"{tib[[endKey]][2]}{Biostrings::reverse(tib[[bodKey]][2])}{tib[[secKey]][2]}{BNG}{tib[[cpnKey]][2]}{RET}",
    "{RET}")
  
  if (verbose>=vt) cat(str)
  
  str
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Genomic Substring Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

parse_design_seq = function(pos, chr, off=60, len=122,
                            verbose=0,vt=3,tc=1,tt=NULL,
                            funcTag='parse_design_seq') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, pos={pos}, len={len}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    chr_len <- chr %>% length()
    cat(glue::glue("[{funcTag}]:{tabsStr} chr_len={chr_len}.{RET}"))
    
    ret_tib <- Biostrings::subseq(chr, start=pos, width=len) %>% 
      as.character() # %>% addBrac()
    
    ret_cnt <- ret_tib %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      MisMatch Probe Comparison Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cmpInfIMU_MisMatch = function(tib, fieldA, fieldB, mu, del='_',
                              verbose=0,vt=5,tc=1,tt=NULL,
                              funcTag='cmpInfIMU_MisMatch') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]:     mu={mu}{RET}"))
    cat(glue::glue("[{funcTag}]: fieldA={fieldA}.{RET}"))
    cat(glue::glue("[{funcTag}]: fieldB={fieldB}.{RET}{RET}"))
  }
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(
      BOD_NumMM=mapply(
        adist,
        stringr::str_sub(!!fieldA,1,stringr::str_length(!!fieldA)-1),
        stringr::str_sub(!!fieldB,1,stringr::str_length(!!fieldB)-1),
        ignore.case=TRUE ),
      
      DI_NUC_AB=paste0(
        stringr::str_to_upper(
          stringr::str_sub(!!fieldA,
                           stringr::str_length(!!fieldA),
                           stringr::str_length(!!fieldA)) ),
        stringr::str_to_upper(
          stringr::str_sub(!!fieldB,
                           stringr::str_length(!!fieldB),
                           stringr::str_length(!!fieldB)) )
      ),
      TAR_EQU=cmpIUPACs(DI_NUC_AB)
    ) %>%
    dplyr::rename(!!paste('BOD_NumMM',mu, sep=del):=BOD_NumMM,
                  !!paste('TAR_EQU',  mu, sep=del):=TAR_EQU)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: Done.{RET}{RET}"))
  
  tib
}

cmpInfI_MisMatch = function(tib, fieldAU, fieldBU, fieldAM, fieldBM, del='_',
                            verbose=0,vt=5,tc=1,tt=NULL,
                            funcTag='cmpInfI_MisMatch') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]: fieldAU={fieldAU}.{RET}"))
    cat(glue::glue("[{funcTag}]: fieldBU={fieldBU}.{RET}"))
    cat(glue::glue("[{funcTag}]: fieldAM={fieldAM}.{RET}"))
    cat(glue::glue("[{funcTag}]: fieldBM={fieldBM}.{RET}"))
  }
  
  tib <- tib %>% 
    cmpInfIMU_MisMatch(fieldAU, fieldBU, mu='U', del=del,
                       verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>% 
    cmpInfIMU_MisMatch(fieldAM, fieldBM, mu='M', del=del,
                       verbose=verbose, vt=vt+1,tc=tc+1)
  
  tib <- tib %>% dplyr::mutate(
    Man_MisMatch=(BOD_NumMM_U+BOD_NumMM_M)/2, #, na.rm=TRUE),
    Man_TarMatch=case_when(TAR_EQU_U & TAR_EQU_M ~ TRUE, TRUE ~ FALSE) ) %>%
    dplyr::select(-c(BOD_NumMM_U,BOD_NumMM_M,TAR_EQU_U,TAR_EQU_M))
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: Done.{RET}{RET}"))
  
  tib
}

cmpInfII_MisMatch = function(tib, fieldA, fieldB, mu='D', del='_',
                             verbose=0,vt=4,tc=1,tt=NULL,
                             funcTag='cmpInfII_MisMatch') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("[{funcTag}]: fieldA={fieldA}.{RET}"))
    cat(glue::glue("[{funcTag}]: fieldB={fieldB}.{RET}{RET}"))
  }
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    cmpInfIMU_MisMatch(fieldA, fieldB, mu=mu, del=del,
                       verbose=verbose, vt=vt+1,tc=tc+1) %>%
    dplyr::rename(
      Man_MisMatch=BOD_NumMM_D,
      Man_TarMatch=TAR_EQU_D)
  # dplyr::select(-c(BOD_NumMM_U,BOD_NumMM_M,TAR_EQU_U,TAR_EQU_M))
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: Done.{RET}{RET}"))
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Exact Probe Comparison Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cmpInfIMU= function(tib, fieldA, fieldB, mu, del='_',
                    verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'cmpInfMU'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldA={fieldA}, fieldB={fieldB} mu={mu}{RET}"))
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(
      SUB_SEQ_A=stringr::str_to_upper(stringr::str_sub(!!fieldA,1,stringr::str_length(!!fieldA)-1)),
      SUB_SEQ_B=stringr::str_to_upper(stringr::str_sub(!!fieldB,1,stringr::str_length(!!fieldB)-1)),
      DI_NUC_AB=paste0(
        stringr::str_to_upper(stringr::str_sub(!!fieldA,stringr::str_length(!!fieldA),stringr::str_length(!!fieldA)) ),
        stringr::str_to_upper(stringr::str_sub(!!fieldB,stringr::str_length(!!fieldB),stringr::str_length(!!fieldB)) )
      ),
      
      BOD_EQU=case_when(SUB_SEQ_A==SUB_SEQ_B ~ TRUE, TRUE ~ FALSE),
      TAR_EQU=cmpIUPACs(DI_NUC_AB),
      Inf1_Match=case_when(BOD_EQU & BOD_EQU==TAR_EQU ~ TRUE, TRUE ~ FALSE)
    ) %>%
    dplyr::select(-c(SUB_SEQ_A,SUB_SEQ_B,BOD_EQU,DI_NUC_AB,TAR_EQU)) %>%
    dplyr::rename(!!paste('Inf1_Match',mu, sep=del):=Inf1_Match)
  
  if (verbose>=vt+1) print(tib)
  
  tib
}

cmpInfI = function(tib, fieldAU, fieldBU, fieldAM, fieldBM, del='_',
                   verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'cmpInfI'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldAU={fieldAU}, fieldBU={fieldBU}{RET}"))
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldAM={fieldAM}, fieldBM={fieldBM}{RET}"))
  
  tib <- tib %>% cmpInfIMU(fieldAU, fieldBU, mu='U', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>% cmpInfIMU(fieldAM, fieldBM, mu='M', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>%
    dplyr::mutate(Man_Match=case_when(Inf1_Match_U & Inf1_Match_M ~ TRUE, TRUE ~ FALSE) ) %>%
    dplyr::select(-c(Inf1_Match_U,Inf1_Match_M))
  
  tib
}

cmpInfII = function(tib, fieldA, fieldB, mu='D', del='_',
                    verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'cmpInfI'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(Man_Match=stringr::str_to_upper(!!fieldA)==stringr::str_to_upper(!!fieldB) )
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Illumina Strand Methods:: TOP/BOT
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

setTopBot_tib = function(tib, seqKey, srdKey, topKey=NULL, max=0,
                         verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'setTopBot_tib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; select field={seqKey}, new field={srdKey}...{RET}"))
  
  seqKey <- rlang::sym(seqKey)
  srdKey <- rlang::sym(srdKey)
  if (!is.null(topKey)) topKey <- rlang::sym(topKey)
  
  # TBD:: Should remove any seqKey that is na or not the correct length...
  
  ret_tib <- NULL
  stime <- system.time({
    if (max!=0) tib <- tib %>% head(n=max)
    
    bit_tib <- tib %>% 
      dplyr::mutate(
        pre_seq = !!seqKey %>% stringr::str_remove('\\[.*$') %>% 
          Biostrings::reverse() %>%
          stringr::str_to_upper() %>% 
          stringr::str_replace_all('A', '1') %>% 
          stringr::str_replace_all('T', '1') %>% 
          stringr::str_replace_all('C', '2') %>% 
          stringr::str_replace_all('G', '2'),
        
        pos_seq = !!seqKey %>% stringr::str_remove('^.*\\]') %>%
          stringr::str_to_upper() %>% 
          stringr::str_replace_all('A', '1') %>% 
          stringr::str_replace_all('T', '1') %>% 
          stringr::str_replace_all('C', '2') %>% 
          stringr::str_replace_all('G', '2')
      ) %>% dplyr::select(pre_seq, pos_seq)
    
    pre_bit <- bit_tib$pre_seq %>% stringr::str_split('', simplify = TRUE)
    pos_bit <- bit_tib$pos_seq %>% stringr::str_split('', simplify = TRUE)
    
    pre_mat <- pre_bit %>% as.data.frame() %>% data.matrix()
    pos_mat <- pos_bit %>% as.data.frame() %>% data.matrix()
    dif_mat <- pre_mat-pos_mat
    dif_mat <- dif_mat %>% cbind(rep(2, dim(dif_mat)[1]))
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} dif_mat={RET}"))
    if (verbose>=vt+4) print(dif_mat)
    
    # This is just for a sanity check...
    #  dif_mat[1,1:60] = 0
    
    bit_vec <- apply(dif_mat,1, function(x) head(x[x!=0],1)) + 1
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} bit_vec={RET}"))
    if (verbose>=vt+4) print(bit_vec)
    
    # Convert numeric values back to characters (T=TOP,B=BOT,U=Unknown)
    tb_vec <- bit_vec %>% as.character() %>% 
      stringr::str_replace('0', 'T') %>% 
      stringr::str_replace('2', 'B') %>%
      stringr::str_replace('3', 'U')
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tb_vec={RET}"))
    if (verbose>=vt+4) print(tb_vec)
    
    ret_tib <- tib %>% dplyr::mutate(!!srdKey := tb_vec)
    
    if (!is.null(topKey)) {
      revCmp_vec <- ret_tib %>% 
        dplyr::pull(!!seqKey) %>% 
        shearBrac() %>% 
        revCmp() %>% 
        addBrac()
      ret_tib <- ret_tib %>% dplyr::mutate(TMP_REVCMP_SEQ=revCmp_vec)
      
      ret_tib <- ret_tib %>% dplyr::mutate(
        !!topKey:=dplyr::case_when(
          !!srdKey=='T' ~ !!seqKey,
          !!srdKey=='B' ~ TMP_REVCMP_SEQ,
          TRUE ~ NA_character_)
      ) %>% dplyr::select(-TMP_REVCMP_SEQ)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# Single variable method that is OUT-OF-DATE; Should use the method above...
isTopBot_single = function(seq, verbose=0, vt=4) {
  seq <- seq %>% stringr::str_to_upper() %>% 
    stringr::str_replace_all('T', 'A') %>% 
    stringr::str_replace_all('G', 'C')
  
  seq_vec <- stringr::str_split(seq,'\\[', simplify=TRUE) %>% as.vector()
  pre_seq <- seq_vec[1] %>% Biostrings::reverse()
  seq_vec <- stringr::str_split(seq,'\\]', simplify=TRUE) %>% as.vector()
  pos_seq <- seq_vec[2]
  
  pre_vec <- pre_seq %>% stringr::str_split('', simplify=TRUE) %>% as.vector()
  pos_vec <- pos_seq %>% stringr::str_split('', simplify=TRUE) %>% as.vector()
  
  if (verbose>=vt) {
    print(seq)
    print(pre_seq)
    print(pos_seq)
  }
  
  min_len <- min(stringr::str_length(pre_seq),stringr::str_length(pos_seq))
  for (ii in c(1:min_len)) {
    if (pre_vec[ii]==pos_vec[ii]) next
    if (pre_vec[ii]!=pos_vec[ii]) {
      if (pre_vec[ii]=='A') return('T')
      return('B')
    }
  }
  return('U')
}

isTopBots = function(x, verbose=0, vt=4) {
  x <- lapply(x,isTopBot_single)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Basic Bisulfite Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

MAPDi = function(x) {
  if (length(MAP_DI[[x]])==0) return(NA)
  MAP_DI[[x]]
}
mapDIs = function(x) {
  x <- lapply(x, MAPDi) %>% BiocGenerics::unlist()
}

# De-Methylate Probe Sequence for 4-base aligners
#
deM = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  
  # aln_seq=prb_seq %>%
  #   stringr::str_replace_all("R","A") %>% # A/G
  #   stringr::str_replace_all("Y","T") %>% # C/T
  #   
  #   stringr::str_replace_all("S","C") %>% # G/C
  #   stringr::str_replace_all("W","A") %>% # A/T
  #   stringr::str_replace_all("K","T") %>% # G/T
  #   stringr::str_replace_all("M","A") %>% # A/C
  #   
  #   stringr::str_replace_all("B","T") %>% # C/G/T
  #   stringr::str_replace_all("D","A") %>% # A/G/T
  #   stringr::str_replace_all("H","A") %>% # A/C/T
  #   stringr::str_replace_all("V","A") %>% # A/C/G
  #   
  #   stringr::str_replace_all("N","A"), # A/C/T/G
  
  if (uc) x <- tr(x, 'RYSWKMBDHV', 'ATCATATAAA')
  else    x <- tr(x, 'RYSWKMBDHV', 'atcatataaa')
  x
}
deMs = function(x, uc=FALSE) { deM(x, uc) }

bscU = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  if (uc) x <- tr(x, 'CYSMBHV', 'TTKWKWD')
  else    x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  x
}
bscUs = function(x, uc=FALSE) { bscU(x, uc) }

MAPM = function(x) {
  if (length(MAP_M[[x]])==0) return(x)
  MAP_M[[x]]
}
bscM = function(x) { stringr::str_replace_all(x, '([CYSMBHV][GRSKBDV])', MAPM) }
bscMs = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  x <- lapply(x, bscM) %>% BiocGenerics::unlist()
  x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  if (uc) x <- stringr::str_to_upper(x)
  x
}

MAPD = function(x) {
  if (length(MAP_D[[x]])==0) return(x)
  MAP_D[[x]]
}
bscD = function(x) { stringr::str_replace_all(x, '([CYSMBHV][GRSKBDV])', MAPD) }
bscDs = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  x <- lapply(x, bscD) %>% BiocGenerics::unlist()
  x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  if (uc) x <- stringr::str_to_upper(x)
  x
}

QMAP = function(x, mu) {
  if (mu=='U') {
    return(QMAP_U[[x]])
  } else if (mu=='M') {
    return(QMAP_M[[x]])
  }
  x
}
qmaps = function(x, mu) {
  x <- lapply(x, QMAP, mu) %>% BiocGenerics::unlist()
}

cmpIUPAC = function(x) {
  if (base::is.element(x, names(IUPAC_EQ) )) return(IUPAC_EQ[[x]])
  # if (is.null(IUPAC_EQ[[x]])) return(FALSE)
  # if (length(IUPAC_EQ[[x]])==0) return(FALSE)
  FALSE
}

cmpIUPACs = function(x) {
  x <- lapply(x, cmpIUPAC) %>% BiocGenerics::unlist()
}

shearBrac = function(x) {
  x %>% stringr::str_remove('\\[') %>% stringr::str_remove('\\]')
}

addBrac = function(seq) {
  seq_len <- stringr::str_length(seq)
  mid_len <- 2
  snp_idx <- as.integer(seq_len / 2)
  
  pre_beg <- 1
  pre_end <- snp_idx-1
  
  mid_beg <- snp_idx
  mid_end <- mid_beg+mid_len-1
  
  pos_beg <- mid_end+1
  pos_end <- seq_len
  
  pre_seq <- seq %>% stringr::str_sub(pre_beg,pre_end)
  mid_seq <- seq %>% stringr::str_sub(mid_beg,mid_end)
  pos_seq <- seq %>% stringr::str_sub(pos_beg,pos_end)
  
  new_seq <- paste0(pre_seq,'[',mid_seq,']',pos_seq)
  
  new_seq
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Basic Reverse/Complement Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
revCmp = function(x) {
  Biostrings::reverse(x) %>% cmpl()
}

revSeq = function(x) {
  Biostrings::reverse(x)
}

cmpl = function(x) {
  tr(x, 'ACTGRYSWKMBDHVactgryswkmbdhv[]','TGACYRSWMKVHDBtgacyrswmkvhdb][')
  # x <- tr(x, 'ACTGRYSWKMBDHV','TGACYRSWMKVHDB')
  # x <- tr(x, 'actgryswkmbdhv','tgacyrswmkvhdb')
  # x
}

tr = function(x, old, new) {
  Biostrings::chartr(old, new, x)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Improbe Header Mapping Structure::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

INIT_IMP_HEADER = function() {
  IMP_TYPES = cols(
    Seq_ID = col_character(),
    Forward_Sequence = col_character(),
    Genome_Build = col_character(),
    Chromosome = col_character(),
    Coordinate = col_double(),
    Design_State = col_character(),
    Seq_Length = col_double(),
    Forward_CpG_Coord = col_double(),
    TB_Strand = col_character(),
    Top_Sequence = col_character(),
    Top_CpG_Coord = col_double(),
    Probe_Type = col_character(),
    Probeset_ID = col_character(),
    Probeset_Score = col_double(),
    Methyl_Probe_ID = col_character(),
    Methyl_Probe_Sequence = col_character(),
    Methyl_Probe_Length = col_double(),
    Methyl_Start_Coord = col_double(),
    Methyl_End_Coord = col_double(),
    Methyl_Probe_Covered_Top_Sequence = col_character(),
    Methyl_Allele_FR_Strand = col_character(),
    Methyl_Allele_TB_Strand = col_character(),
    Methyl_Allele_CO_Strand = col_character(),
    Methyl_Allele_Type = col_character(),
    Methyl_Final_Score = col_double(),
    Methyl_Tm = col_double(),
    Methyl_Tm_Score = col_double(),
    Methyl_GC_Percent = col_double(),
    Methyl_GC_Score = col_double(),
    Methyl_13mer_Count = col_double(),
    Methyl_13mer_Score = col_double(),
    Methyl_Address_Count = col_double(),
    Methyl_Address_Score = col_double(),
    Methyl_Self_Complementarity = col_double(),
    Methyl_Self_Complementarity_Score = col_double(),
    Methyl_Mono_Run = col_double(),
    Methyl_Mono_Run_Score = col_double(),
    Methyl_Ectopic_Count = col_double(),
    Methyl_Ectopic_Score = col_double(),
    Methyl_Underlying_CpG_Count = col_double(),
    Methyl_Underlying_CpG_Min_Dist = col_double(),
    Methyl_Underlying_CpG_Score = col_double(),
    Methyl_In_CpG_Island_Relaxed = col_logical(),
    Methyl_CpG_Island_Score = col_double(),
    Methyl_Next_Base = col_character(),
    Methyl_Next_Base_Score = col_double(),
    UnMethyl_Probe_ID = col_character(),
    UnMethyl_Probe_Sequence = col_character(),
    UnMethyl_Probe_Length = col_double(),
    UnMethyl_Start_Coord = col_double(),
    UnMethyl_End_Coord = col_double(),
    UnMethyl_Probe_Covered_Top_Sequence = col_character(),
    UnMethyl_Allele_FR_Strand = col_character(),
    UnMethyl_Allele_TB_Strand = col_character(),
    UnMethyl_Allele_CO_Strand = col_character(),
    UnMethyl_Allele_Type = col_character(),
    UnMethyl_Final_Score = col_double(),
    UnMethyl_Tm = col_double(),
    UnMethyl_Tm_Score = col_double(),
    UnMethyl_GC_Percent = col_double(),
    UnMethyl_GC_Score = col_double(),
    UnMethyl_13mer_Count = col_double(),
    UnMethyl_13mer_Score = col_double(),
    UnMethyl_Address_Count = col_double(),
    UnMethyl_Address_Score = col_double(),
    UnMethyl_Self_Complementarity = col_double(),
    UnMethyl_Self_Complementarity_Score = col_double(),
    UnMethyl_Mono_Run = col_double(),
    UnMethyl_Mono_Run_Score = col_double(),
    UnMethyl_Ectopic_Count = col_double(),
    UnMethyl_Ectopic_Score = col_double(),
    UnMethyl_Underlying_CpG_Count = col_double(),
    UnMethyl_Underlying_CpG_Min_Dist = col_double(),
    UnMethyl_Underlying_CpG_Score = col_double(),
    UnMethyl_In_CpG_Island_Relaxed = col_logical(),
    UnMethyl_CpG_Island_Score = col_double(),
    UnMethyl_Next_Base = col_character(),
    UnMethyl_Next_Base_Score = col_double()
  )
  
  IMP_TYPES
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Nucelotide Mapping Structures::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# mapA = function(x) {
#   if (length(MAP_A[[x]])==0) return(x)
#   MAP_A[[x]]
# }
# INIT_MAP_A = function() {
#   MAP_A <- NULL
#   MAP_A[['AA']] <- 'aa'
#   MAP_A[['Aa']] <- 'aa'
#   MAP_A[['aA']] <- 'aa'
#   MAP_A[['aa']] <- 'aa'
#   
#   MAP_A
# }
# MAP_A <- INIT_MAP_A()

# Di-Nuc to IUPAC::
INIT_MAP_DI = function() {
  MAP <- NULL
  
  # Mono-Nuc Map::
  MAP[['A']] <- 'A'
  MAP[['C']] <- 'C'
  MAP[['G']] <- 'G'
  MAP[['T']] <- 'T'
  
  # Di-Nuc Map::
  MAP[['AC']] <- 'M'
  MAP[['AG']] <- 'R'
  MAP[['AT']] <- 'W'
  MAP[['AA']] <- 'A'
  
  MAP[['CA']] <- 'M'
  MAP[['CT']] <- 'Y'
  MAP[['CG']] <- 'S'
  MAP[['CC']] <- 'C'
  
  MAP[['GA']] <- 'R'
  MAP[['GT']] <- 'K'
  MAP[['GC']] <- 'S'
  MAP[['GG']] <- 'G'
  
  MAP[['TC']] <- 'Y'
  MAP[['TG']] <- 'K'
  MAP[['TA']] <- 'W'
  MAP[['TT']] <- 'T'
  
  #
  # Tri-Nuc Map:: A
  #
  MAP[['AAC']] <- 'M'
  MAP[['AAG']] <- 'R'
  MAP[['AAT']] <- 'W'
  MAP[['AAA']] <- 'A'
  
  MAP[['ACA']] <- 'M'
  MAP[['ACT']] <- 'H' # NEW-3
  MAP[['ACG']] <- 'V' # NEW-3
  MAP[['ACC']] <- 'M' # NEW-2
  
  MAP[['AGA']] <- 'R'
  MAP[['AGT']] <- 'D' # NEW-3
  MAP[['AGC']] <- 'V' # NEW-3
  MAP[['AGG']] <- 'R' # NEW-2
  
  MAP[['ATC']] <- 'H' # NEW-3
  MAP[['ATG']] <- 'D' # NEW-3
  MAP[['ATA']] <- 'W'
  MAP[['ATT']] <- 'W' # NEW-2
  
  #
  # Tri-Nuc Map:: C
  #
  MAP[['CAC']] <- 'M'
  MAP[['CAG']] <- 'V' # NEW-3
  MAP[['CAT']] <- 'H' # NEW-3
  MAP[['CAA']] <- 'M' # NEW-2
  
  MAP[['CCA']] <- 'M'
  MAP[['CCT']] <- 'Y'
  MAP[['CCG']] <- 'S'
  MAP[['CCC']] <- 'C'
  
  MAP[['CGA']] <- 'V' # NEW-3
  MAP[['CGT']] <- 'B' # NEW-3
  MAP[['CGC']] <- 'S'
  MAP[['CGG']] <- 'S' # NEW-2
  
  MAP[['CTC']] <- 'Y'
  MAP[['CTG']] <- 'B' # NEW-3
  MAP[['CTA']] <- 'H' # NEW-3
  MAP[['CTT']] <- 'Y' # NEW-2
  
  #
  # Tri-Nuc Map:: G
  #
  MAP[['GAC']] <- 'V' # NEW-3
  MAP[['GAG']] <- 'R'
  MAP[['GAT']] <- 'D' # NEW-3
  MAP[['GAA']] <- 'R' # NEW-2
  
  MAP[['GCA']] <- 'V' # NEW-3
  MAP[['GCT']] <- 'B' # NEW-3
  MAP[['GCG']] <- 'S'
  MAP[['GCC']] <- 'S' # NEW-2
  
  MAP[['GGA']] <- 'R'
  MAP[['GGT']] <- 'K'
  MAP[['GGC']] <- 'S'
  MAP[['GGG']] <- 'G'
  
  MAP[['GTC']] <- 'B' # NEW-3
  MAP[['GTG']] <- 'K'
  MAP[['GTA']] <- 'D' # NEW-3
  MAP[['GTT']] <- 'K' # NEW-2
  
  #
  # Tri-Nuc Map:: T
  #
  MAP[['TAC']] <- 'H' # NEW-3
  MAP[['TAG']] <- 'D' # NEW-3
  MAP[['TAT']] <- 'W'
  MAP[['TAA']] <- 'W' # NEW-2
  
  MAP[['TCA']] <- 'H' # NEW-3
  MAP[['TCT']] <- 'Y'
  MAP[['TCG']] <- 'B' # NEW-3
  MAP[['TCC']] <- 'Y' # NEW-2
  
  MAP[['TGA']] <- 'D' # NEW-3
  MAP[['TGT']] <- 'K'
  MAP[['TGC']] <- 'B' # NEW-3
  MAP[['TGG']] <- 'K' # NEW-2
  
  MAP[['TTC']] <- 'Y'
  MAP[['TTG']] <- 'K'
  MAP[['TTA']] <- 'W'
  MAP[['TTT']] <- 'T'
  
  MAP
}
MAP_DI <- INIT_MAP_DI()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Initialize M/D Maps for Bisulfite Conversion::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Data generated with:: /Users/bbarnes/Documents/Projects/scripts/mapMD.pl
INIT_MAP_M = function() {
  MAP <- NULL
  MAP[['CG']] <- 'cG'
  MAP[['CR']] <- 'cR'
  MAP[['CS']] <- 'cS'
  MAP[['CK']] <- 'cK'
  MAP[['CB']] <- 'cB'
  MAP[['CD']] <- 'cD'
  MAP[['CV']] <- 'cV'
  MAP[['YG']] <- 'yG'
  MAP[['YR']] <- 'yR'
  MAP[['YS']] <- 'yS'
  MAP[['YK']] <- 'yK'
  MAP[['YB']] <- 'yB'
  MAP[['YD']] <- 'yD'
  MAP[['YV']] <- 'yV'
  MAP[['SG']] <- 'sG'
  MAP[['SR']] <- 'sR'
  MAP[['SS']] <- 'sS'
  MAP[['SK']] <- 'sK'
  MAP[['SB']] <- 'sB'
  MAP[['SD']] <- 'sD'
  MAP[['SV']] <- 'sV'
  MAP[['MG']] <- 'mG'
  MAP[['MR']] <- 'mR'
  MAP[['MS']] <- 'mS'
  MAP[['MK']] <- 'mK'
  MAP[['MB']] <- 'mB'
  MAP[['MD']] <- 'mD'
  MAP[['MV']] <- 'mV'
  MAP[['BG']] <- 'bG'
  MAP[['BR']] <- 'bR'
  MAP[['BS']] <- 'bS'
  MAP[['BK']] <- 'bK'
  MAP[['BB']] <- 'bB'
  MAP[['BD']] <- 'bD'
  MAP[['BV']] <- 'bV'
  MAP[['HG']] <- 'hG'
  MAP[['HR']] <- 'hR'
  MAP[['HS']] <- 'hS'
  MAP[['HK']] <- 'hK'
  MAP[['HB']] <- 'hB'
  MAP[['HD']] <- 'hD'
  MAP[['HV']] <- 'hV'
  MAP[['VG']] <- 'vG'
  MAP[['VR']] <- 'vR'
  MAP[['VS']] <- 'vS'
  MAP[['VK']] <- 'vK'
  MAP[['VB']] <- 'vB'
  MAP[['VD']] <- 'vD'
  MAP[['VV']] <- 'vV'
  
  MAP
}
MAP_M <- INIT_MAP_M()

INIT_MAP_D = function() {
  MAP <- NULL
  
  MAP[['CG']] <- 'yG'
  MAP[['CR']] <- 'yR'
  MAP[['CS']] <- 'yS'
  MAP[['CK']] <- 'yK'
  MAP[['CB']] <- 'yB'
  MAP[['CD']] <- 'yD'
  MAP[['CV']] <- 'yV'
  MAP[['YG']] <- 'yG'
  MAP[['YR']] <- 'yR'
  MAP[['YS']] <- 'yS'
  MAP[['YK']] <- 'yK'
  MAP[['YB']] <- 'yB'
  MAP[['YD']] <- 'yD'
  MAP[['YV']] <- 'yV'
  MAP[['SG']] <- 'bG'
  MAP[['SR']] <- 'bR'
  MAP[['SS']] <- 'bS'
  MAP[['SK']] <- 'bK'
  MAP[['SB']] <- 'bB'
  MAP[['SD']] <- 'bD'
  MAP[['SV']] <- 'bV'
  MAP[['MG']] <- 'hG'
  MAP[['MR']] <- 'hR'
  MAP[['MS']] <- 'hS'
  MAP[['MK']] <- 'hK'
  MAP[['MB']] <- 'hB'
  MAP[['MD']] <- 'hD'
  MAP[['MV']] <- 'hV'
  MAP[['BG']] <- 'bG'
  MAP[['BR']] <- 'bR'
  MAP[['BS']] <- 'bS'
  MAP[['BK']] <- 'bK'
  MAP[['BB']] <- 'bB'
  MAP[['BD']] <- 'bD'
  MAP[['BV']] <- 'bV'
  MAP[['HG']] <- 'hG'
  MAP[['HR']] <- 'hR'
  MAP[['HS']] <- 'hS'
  MAP[['HK']] <- 'hK'
  MAP[['HB']] <- 'hB'
  MAP[['HD']] <- 'hD'
  MAP[['HV']] <- 'hV'
  MAP[['VG']] <- 'nG'
  MAP[['VR']] <- 'nR'
  MAP[['VS']] <- 'nS'
  MAP[['VK']] <- 'nK'
  MAP[['VB']] <- 'nB'
  MAP[['VD']] <- 'nD'
  MAP[['VV']] <- 'nV'
  
  MAP
}
MAP_D <- INIT_MAP_D()



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Initialize Query Maps::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# These are for mapping Query positions to their possible outcomes based on
#  the assumption of unmethylated and methylated
#
# OLD_INIT_QMAP_* are previous versions that did not account for Infinium I
#  probes with some degenerate bases like W={A,T} that would still work
#

# TBD:: Address the tri-nucelotide version!

# Unmethylated will always use the letter earliest in the alphabet
#  e.g. 
#   A over T, 
#   A over C, 
#   A over G, 
#   C over G, 
#   C over T, 
#   G over T,
#
INIT_QMAP_U = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'
  MAP[['C']] <- 'C'  # '-'
  MAP[['G']] <- 'G'  # '-'
  MAP[['T']] <- 'T'
  
  MAP[['R']] <- 'A'
  MAP[['Y']] <- 'T'
  MAP[['S']] <- 'C'  # C over G 
  MAP[['W']] <- 'A'  # A over T
  # MAP[['S']] <- 'S'  # '-' Old Version
  # MAP[['W']] <- 'W'  # '-' Old Version
  
  MAP[['K']] <- 'T'
  MAP[['M']] <- 'A'
  MAP[['B']] <- 'T'
  MAP[['D']] <- 'W'
  MAP[['H']] <- 'W'
  MAP[['V']] <- 'A'
  
  MAP[['N']] <- 'W'
  
  # Lower Case::
  MAP[['a']] <- 'a'
  MAP[['c']] <- 'c'  # '-'
  MAP[['g']] <- 'g'  # '-'
  MAP[['t']] <- 't'
  
  MAP[['r']] <- 'a'
  MAP[['y']] <- 't'
  MAP[['s']] <- 'c'  # c over g
  MAP[['w']] <- 'a'  # a over t
  # MAP[['s']] <- 's'  # '-' Old Version
  # MAP[['w']] <- 'w'  # '-' Old Version
  
  MAP[['k']] <- 't'
  MAP[['m']] <- 'a'
  MAP[['b']] <- 't'
  MAP[['d']] <- 'w'
  MAP[['h']] <- 'w'
  MAP[['v']] <- 'a'
  
  MAP[['n']] <- 'w'
  
  MAP
}
QMAP_U <- INIT_QMAP_U()

# Methylated will always use the letter latest in the alphabet
#  e.g. 
#   T over A, 
#   T over C, 
#   T over G,
#   G over A, 
#   G over C, 
#   C over A, 
#
INIT_QMAP_M = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'  # '-'
  MAP[['C']] <- 'C'
  MAP[['G']] <- 'G'
  MAP[['T']] <- 'T'  # '-'
  
  MAP[['R']] <- 'G'
  MAP[['Y']] <- 'C'
  MAP[['S']] <- 'G'  # G over C
  MAP[['W']] <- 'T'  # T over A
  # MAP[['S']] <- 'S'  # '-' Old Version
  # MAP[['W']] <- 'W'  # '-' Old Version
  
  MAP[['K']] <- 'G'
  MAP[['M']] <- 'C'
  MAP[['B']] <- 'S'
  MAP[['D']] <- 'G'
  MAP[['H']] <- 'C'
  MAP[['V']] <- 'S'
  
  MAP[['N']] <- 'S'
  
  # Lower Case::
  MAP[['a']] <- 'a'  # '-'
  MAP[['c']] <- 'c'
  MAP[['g']] <- 'g'
  MAP[['t']] <- 't'  # '-'
  
  MAP[['r']] <- 'g'
  MAP[['y']] <- 'c'
  MAP[['s']] <- 'g'  # g over c
  MAP[['w']] <- 't'  # t over a
  # MAP[['s']] <- 's'  # '-' Old Version
  # MAP[['w']] <- 'w'  # '-' Old Version
  
  MAP[['k']] <- 'g'
  MAP[['m']] <- 'c'
  MAP[['b']] <- 's'
  MAP[['d']] <- 'g'
  MAP[['h']] <- 'c'
  MAP[['v']] <- 's'
  
  MAP[['n']] <- 's'
  
  MAP
}
QMAP_M <- INIT_QMAP_M()

# OLD Versions to be used now::
#
INIT_QMAP_U = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'
  MAP[['C']] <- 'C'  # '-'
  MAP[['G']] <- 'G'  # '-'
  MAP[['T']] <- 'T'
  
  MAP[['R']] <- 'A'
  MAP[['Y']] <- 'T'
  MAP[['S']] <- 'S'  # '-'
  MAP[['W']] <- 'W'  # '-'
  
  MAP[['K']] <- 'T'
  MAP[['M']] <- 'A'
  MAP[['B']] <- 'T'
  MAP[['D']] <- 'W'
  MAP[['H']] <- 'W'
  MAP[['V']] <- 'A'
  
  MAP[['N']] <- 'W'
  
  # Lower Case::
  MAP[['a']] <- 'a'
  MAP[['c']] <- 'c'  # '-'
  MAP[['g']] <- 'g'  # '-'
  MAP[['t']] <- 't'
  
  MAP[['r']] <- 'a'
  MAP[['y']] <- 't'
  MAP[['s']] <- 's'  # '-'
  MAP[['w']] <- 'w'  # '-'
  
  MAP[['k']] <- 't'
  MAP[['m']] <- 'a'
  MAP[['b']] <- 't'
  MAP[['d']] <- 'w'
  MAP[['h']] <- 'w'
  MAP[['v']] <- 'a'
  
  MAP[['n']] <- 'w'
  
  MAP
}
# OLD_QMAP_U <- OLD_INIT_QMAP_U()

OLD_INIT_QMAP_M = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'  # '-'
  MAP[['C']] <- 'C'
  MAP[['G']] <- 'G'
  MAP[['T']] <- 'T'  # '-'
  
  MAP[['R']] <- 'G'
  MAP[['Y']] <- 'C'
  MAP[['S']] <- 'S'  # '-'
  MAP[['W']] <- 'W'  # '-'
  
  MAP[['K']] <- 'G'
  MAP[['M']] <- 'C'
  MAP[['B']] <- 'S'
  MAP[['D']] <- 'G'
  MAP[['H']] <- 'C'
  MAP[['V']] <- 'S'
  
  MAP[['N']] <- 'S'
  
  # Lower Case::
  MAP[['a']] <- 'a'  # '-'
  MAP[['c']] <- 'c'
  MAP[['g']] <- 'g'
  MAP[['t']] <- 't'  # '-'
  
  MAP[['r']] <- 'g'
  MAP[['y']] <- 'c'
  MAP[['s']] <- 's'  # '-'
  MAP[['w']] <- 'w'  # '-'
  
  MAP[['k']] <- 'g'
  MAP[['m']] <- 'c'
  MAP[['b']] <- 's'
  MAP[['d']] <- 'g'
  MAP[['h']] <- 'c'
  MAP[['v']] <- 's'
  
  MAP[['n']] <- 's'
  
  MAP
}
# QMAP_M <- OLD_INIT_QMAP_M()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Initialize IUPAC Equivelency Tables::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# This is really for convience...
#
INIT_IUPAC_EQ = function() {
  MAP <- NULL
  
  MAP[['AA']] <- TRUE
  MAP[['AD']] <- TRUE
  MAP[['AH']] <- TRUE
  MAP[['AM']] <- TRUE
  MAP[['AN']] <- TRUE
  MAP[['AR']] <- TRUE
  MAP[['AV']] <- TRUE
  MAP[['AW']] <- TRUE
  MAP[['BB']] <- TRUE
  MAP[['BC']] <- TRUE
  MAP[['BG']] <- TRUE
  MAP[['BK']] <- TRUE
  MAP[['BN']] <- TRUE
  MAP[['BS']] <- TRUE
  MAP[['BT']] <- TRUE
  MAP[['BY']] <- TRUE
  MAP[['CB']] <- TRUE
  MAP[['CC']] <- TRUE
  MAP[['CH']] <- TRUE
  MAP[['CM']] <- TRUE
  MAP[['CN']] <- TRUE
  MAP[['CS']] <- TRUE
  MAP[['CV']] <- TRUE
  MAP[['CY']] <- TRUE
  MAP[['DA']] <- TRUE
  MAP[['DD']] <- TRUE
  MAP[['DG']] <- TRUE
  MAP[['DK']] <- TRUE
  MAP[['DN']] <- TRUE
  MAP[['DR']] <- TRUE
  MAP[['DT']] <- TRUE
  MAP[['DW']] <- TRUE
  MAP[['GB']] <- TRUE
  MAP[['GD']] <- TRUE
  MAP[['GG']] <- TRUE
  MAP[['GK']] <- TRUE
  MAP[['GN']] <- TRUE
  MAP[['GR']] <- TRUE
  MAP[['GS']] <- TRUE
  MAP[['GV']] <- TRUE
  MAP[['HA']] <- TRUE
  MAP[['HC']] <- TRUE
  MAP[['HH']] <- TRUE
  MAP[['HM']] <- TRUE
  MAP[['HN']] <- TRUE
  MAP[['HT']] <- TRUE
  MAP[['HW']] <- TRUE
  MAP[['HY']] <- TRUE
  MAP[['KB']] <- TRUE
  MAP[['KD']] <- TRUE
  MAP[['KG']] <- TRUE
  MAP[['KK']] <- TRUE
  MAP[['KN']] <- TRUE
  MAP[['KT']] <- TRUE
  MAP[['MA']] <- TRUE
  MAP[['MC']] <- TRUE
  MAP[['MH']] <- TRUE
  MAP[['MM']] <- TRUE
  MAP[['MN']] <- TRUE
  MAP[['MV']] <- TRUE
  MAP[['NA']] <- TRUE
  MAP[['NB']] <- TRUE
  MAP[['NC']] <- TRUE
  MAP[['ND']] <- TRUE
  MAP[['NG']] <- TRUE
  MAP[['NH']] <- TRUE
  MAP[['NK']] <- TRUE
  MAP[['NM']] <- TRUE
  MAP[['NN']] <- TRUE
  MAP[['NR']] <- TRUE
  MAP[['NS']] <- TRUE
  MAP[['NT']] <- TRUE
  MAP[['NV']] <- TRUE
  MAP[['NW']] <- TRUE
  MAP[['NY']] <- TRUE
  MAP[['RA']] <- TRUE
  MAP[['RD']] <- TRUE
  MAP[['RG']] <- TRUE
  MAP[['RN']] <- TRUE
  MAP[['RR']] <- TRUE
  MAP[['RV']] <- TRUE
  MAP[['SB']] <- TRUE
  MAP[['SC']] <- TRUE
  MAP[['SG']] <- TRUE
  MAP[['SN']] <- TRUE
  MAP[['SS']] <- TRUE
  MAP[['SV']] <- TRUE
  MAP[['TB']] <- TRUE
  MAP[['TD']] <- TRUE
  MAP[['TH']] <- TRUE
  MAP[['TK']] <- TRUE
  MAP[['TN']] <- TRUE
  MAP[['TT']] <- TRUE
  MAP[['TW']] <- TRUE
  MAP[['TY']] <- TRUE
  MAP[['VA']] <- TRUE
  MAP[['VC']] <- TRUE
  MAP[['VG']] <- TRUE
  MAP[['VM']] <- TRUE
  MAP[['VN']] <- TRUE
  MAP[['VR']] <- TRUE
  MAP[['VS']] <- TRUE
  MAP[['VV']] <- TRUE
  MAP[['WA']] <- TRUE
  MAP[['WD']] <- TRUE
  MAP[['WH']] <- TRUE
  MAP[['WN']] <- TRUE
  MAP[['WT']] <- TRUE
  MAP[['WW']] <- TRUE
  MAP[['YB']] <- TRUE
  MAP[['YC']] <- TRUE
  MAP[['YH']] <- TRUE
  MAP[['YN']] <- TRUE
  MAP[['YT']] <- TRUE
  MAP[['YY']] <- TRUE
  
  MAP
}
IUPAC_EQ <- INIT_IUPAC_EQ()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Temp Code to be deleted::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  ret_tib <- readr::read_tsv(file, # guess_max=1000000)
                             col_types = cols(
                               Seq_ID = col_character(),
                               Forward_Sequence = col_character(),
                               Genome_Build = col_character(),
                               Chromosome = col_character(),
                               Coordinate = col_double(),
                               Design_State = col_character(),
                               Seq_Length = col_double(),
                               Forward_CpG_Coord = col_double(),
                               TB_Strand = col_character(),
                               Top_Sequence = col_character(),
                               Top_CpG_Coord = col_double(),
                               Probe_Type = col_character(),
                               Probeset_ID = col_character(),
                               Probeset_Score = col_double(),
                               Methyl_Probe_ID = col_character(),
                               Methyl_Probe_Sequence = col_character(),
                               Methyl_Probe_Length = col_double(),
                               Methyl_Start_Coord = col_double(),
                               Methyl_End_Coord = col_double(),
                               Methyl_Probe_Covered_Top_Sequence = col_character(),
                               Methyl_Allele_FR_Strand = col_character(),
                               Methyl_Allele_TB_Strand = col_character(),
                               Methyl_Allele_CO_Strand = col_character(),
                               Methyl_Allele_Type = col_character(),
                               Methyl_Final_Score = col_double(),
                               Methyl_Tm = col_double(),
                               Methyl_Tm_Score = col_double(),
                               Methyl_GC_Percent = col_double(),
                               Methyl_GC_Score = col_double(),
                               Methyl_13mer_Count = col_double(),
                               Methyl_13mer_Score = col_double(),
                               Methyl_Address_Count = col_double(),
                               Methyl_Address_Score = col_double(),
                               Methyl_Self_Complementarity = col_double(),
                               Methyl_Self_Complementarity_Score = col_double(),
                               Methyl_Mono_Run = col_double(),
                               Methyl_Mono_Run_Score = col_double(),
                               Methyl_Ectopic_Count = col_double(),
                               Methyl_Ectopic_Score = col_double(),
                               Methyl_Underlying_CpG_Count = col_double(),
                               Methyl_Underlying_CpG_Min_Dist = col_double(),
                               Methyl_Underlying_CpG_Score = col_double(),
                               Methyl_In_CpG_Island_Relaxed = col_logical(),
                               Methyl_CpG_Island_Score = col_double(),
                               Methyl_Next_Base = col_character(),
                               Methyl_Next_Base_Score = col_double(),
                               UnMethyl_Probe_ID = col_character(),
                               UnMethyl_Probe_Sequence = col_character(),
                               UnMethyl_Probe_Length = col_double(),
                               UnMethyl_Start_Coord = col_double(),
                               UnMethyl_End_Coord = col_double(),
                               UnMethyl_Probe_Covered_Top_Sequence = col_character(),
                               UnMethyl_Allele_FR_Strand = col_character(),
                               UnMethyl_Allele_TB_Strand = col_character(),
                               UnMethyl_Allele_CO_Strand = col_character(),
                               UnMethyl_Allele_Type = col_character(),
                               UnMethyl_Final_Score = col_double(),
                               UnMethyl_Tm = col_double(),
                               UnMethyl_Tm_Score = col_double(),
                               UnMethyl_GC_Percent = col_double(),
                               UnMethyl_GC_Score = col_double(),
                               UnMethyl_13mer_Count = col_double(),
                               UnMethyl_13mer_Score = col_double(),
                               UnMethyl_Address_Count = col_double(),
                               UnMethyl_Address_Score = col_double(),
                               UnMethyl_Self_Complementarity = col_double(),
                               UnMethyl_Self_Complementarity_Score = col_double(),
                               UnMethyl_Mono_Run = col_double(),
                               UnMethyl_Mono_Run_Score = col_double(),
                               UnMethyl_Ectopic_Count = col_double(),
                               UnMethyl_Ectopic_Score = col_double(),
                               UnMethyl_Underlying_CpG_Count = col_double(),
                               UnMethyl_Underlying_CpG_Min_Dist = col_double(),
                               UnMethyl_Underlying_CpG_Score = col_double(),
                               UnMethyl_In_CpG_Island_Relaxed = col_logical(),
                               UnMethyl_CpG_Island_Score = col_double(),
                               UnMethyl_Next_Base = col_character(),
                               UnMethyl_Next_Base_Score = col_double()
                             )
  )
}


# End of file
