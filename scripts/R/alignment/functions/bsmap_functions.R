
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    2.0 Align All Probe Sequence:: BSMAP
#
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
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    2.0 Align All Probe Sequence:: BSMAP
#                           Main Workflow Driver
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsp_mapping_workflow = 
  function(ref_fas = NULL, 
           ref_tib = NULL,
           
           can_fas = NULL,
           can_tib,
           
           cgn_src = NULL,

           ids_key = "Prb_Key",
           unq_key = "Unq_Key",
           prb_key = "Prb_Seq",
           des_key = "Ord_Des",
           din_key = "Ord_Din",
           
           join_key="Prb_Key",
           join_type="inner",
           
           sort   = TRUE, 
           full   = FALSE, 
           merge  = TRUE,
           light  = FALSE, 
           reload = FALSE,
           retData = FALSE,
           
           bsp_dir = NULL,
           bsp_exe,
           bsp_opt = "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R",
           
           out_csv = NULL,
           out_dir,
           run_tag, 
           re_load = FALSE,
           pre_tag = NULL,
           end_str = 'csv.gz',
           sep_chr = '.',
           out_col = c("Prb_Key","Address","Ord_Des",
                       "Ord_Din","Ord_Map","Ord_Prb"),

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
    
    can_fas <- out_csv %>%
      stringr::str_remove(paste0(sep_chr,end_str,"$") ) %>%
      paste("bsp-probes.fa.gz", sep=sep_chr)
    
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
    ret_dat <- NULL
    
    if (is.null(ref_fas) && is.null(ref_tib)) {
      cat(glue::glue("{RET}{mssg} Both ref_fas AND ref_tib can NOT be NULL!{RET2}"))
      return(ret_tib)
    }
    
    if (is.null(can_fas) && is.null(can_tib)) {
      cat(glue::glue("{RET}{mssg} Both can_fas AND can_tib can NOT be NULL!{RET2}"))
      return(ret_tib)
    }
    
    stime <- base::system.time({
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                        Pre-processing BSMAP Inputs:: 
      #                       Reference/Candidate Fasta Files
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ref_tib <- format_reference_bsmap(fas = ref_fas,
                                        tib = ref_tib,
                                        
                                        out_dir = out_dir,
                                        bsc_key = "N", 
                                        frs_key = "F", 
                                        alphabet_key = "dna",
                                        
                                        verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
      
      can_tib <- format_candidate_bsmap(fas = can_fas, 
                                        tib = can_tib,
                                        
                                        ids_key = ids_key,
                                        prb_key = prb_key,
                                        aln_key = "Aln_Prb", 
                                        rev_key = "Aln_Rev",
                                        
                                        verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
        
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #         Align Candidate Probes Against all Reference Fastas:: BSMAP
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ref_cnt <- ref_tib %>% base::nrow()
      if (verbose>=vt+2)
        cat(glue::glue("{mssg} Will align candidate fasta={can_fas} against ",
                       "{ref_cnt} reference fasta files...{RET}"))
        
      for (ref_idx in c(1:ref_cnt)) {
        cur_ref_fas <- ref_tib$Path[ref_idx]
        cur_out_dir <- ref_tib$Out_Dir[ref_idx]
        
        if (verbose>=vt+2) {
          cat(glue::glue("{mssg}{TAB} Current Ref Fas = {cur_ref_fas}.{RET}"))
          cat(glue::glue("{mssg}{TAB} Current Out Dir = {cur_out_dir}.{RET}"))
        }
        
        cur_tib <- NULL
        cur_tib <- run_bsmap(ref_fas = cur_ref_fas,
                             can_fas = can_fas,
                             out_dir = cur_out_dir,
                             bsp_exe = bsp_exe,
                             bsp_dir = bsp_dir,
                             bsp_opt = bsp_opt,
                             sort = sort,
                             light = light, 
                             reload = reload,
                             verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
        
        ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
      }
      if (retData) ret_dat$bsp_raw <- ret_tib
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                      Join Address and Alignment Data::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      ret_tib <- join_bsmap(bsp_tib = ret_tib,
                            can_tib = can_tib,
                            
                            out_csv = NULL,
                            cgn_src = cgn_src,
                            
                            unq_key = unq_key,
                            des_key = des_key,
                            din_key = din_key,
                            
                            join_key  = join_key,
                            join_type = join_type,
                            
                            sort = sort,
                            full = full,
                            
                            verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
      
      if (retData) ret_dat$bsp_join <- ret_tib

      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                           Alignment Summary::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # Multiple Hit Summary::
      bsp_hit_sum <- ret_tib %>% 
        dplyr::group_by(Address) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>%
        dplyr::group_by(Count) %>% 
        dplyr::summarise(His_Count=n(), .groups="drop")
      hit_key <- glue::glue("bsp_hit_sum({funcTag})")
      hit_cnt <- print_tib(bsp_hit_sum, funcTag, verbose,vt=vt+4,tc=tc+1, n=hit_key)
      
      # Top Ranked Offenders::
      top_hit_sum <- ret_tib %>% 
        dplyr::group_by(Address) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>%
        dplyr::filter(Count!=1) %>%
        dplyr::arrange(-Count)
      # print(top_hit_sum, n=base::nrow(top_hit_sum))
      
      # Format Final Output::
      ret_tib <- ret_tib %>% 
        dplyr::select(dplyr::all_of(out_col), dplyr::everything())
      
      sum_cnt <- safe_write(x=bsp_hit_sum, file=bsp_hit_csv, funcTag=funcTag,
                            verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
      top_cnt <- safe_write(x=top_hit_sum, file=top_hit_csv, funcTag=funcTag,
                            verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
      out_cnt <- safe_write(x=ret_tib, file=out_csv, funcTag=funcTag, done=TRUE,
                            verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
      tt$addFile(out_csv)

      ret_key <- glue::glue("ret-FIN({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
      
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{mssg}{BRK}{RET2}"))
    
    if (retData) return(ret_dat)
    
    ret_tib
  }


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Pre-processing BSMAP Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

format_reference_bsmap = function(fas = NULL,
                                  tib = NULL,
                                  
                                  out_dir,
                                  bsc_key = "N",
                                  frs_key = "F",
                                  alphabet_key = "dna",
                                  
                                  verbose=0, vt=3,tc=1,tt=NULL,
                                  funcTag='format_reference_bsmap') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}            fas = {fas}.{RET}"))
    cat(glue::glue("{mssg}        bsc_key = {bsc_key}.{RET}"))
    cat(glue::glue("{mssg}        frs_key = {frs_key}.{RET}"))
    cat(glue::glue("{mssg}   alphabet_key = {alphabet_key}.{RET}"))
    # if (verbose>=vt+4) {
    #   cat(glue::glue("{RET}"))
    #   cat(glue::glue("{mssg}   tib = {tib}.{RET}"))
    # }
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    if (!is.null(fas) && file.exists(fas)) {
      ret_tib <- tibble::tibble(
        Path=fas,
        Genome_Base_Name = base::basename(Path) %>% 
          stringr::str_remove(".fa.gz"),
        Genome_Version = "Unknown",
        Molecule_Type  = "Whole_Genome",
        Out_Dir = file.path(
          out_dir,
          Genome_Version,
          Molecule_Type,
          Genome_Base_Name) 
      )
    }
    
    if (!is.null(tib)) {
      ret_tib <- tib %>%
        dplyr::filter( 
          Genome_Alphabet==alphabet_key & 
            Genome_Strand_BSC == bsc_key & 
            Genome_Strand_FR  == frs_key) %>%
        dplyr::mutate(
          Out_Dir=file.path(
            out_dir,
            Genome_Version,
            Molecule_Type,
            Genome_Base_Name)
        ) %>% dplyr::bind_rows(ret_tib)
    }
    
    # Validate that each fasta path exists::
    for (cur_fas in ret_tib$Path) {
      if (is.null(cur_fas) || !file.exists(cur_fas)) {
        stop(glue::glue("{RET}{mssg} ERROR: Failed to find fasta={cur_fas}! ",
                        "user input or reference genome tibble!{RET2}"))
        print(ret_tib)
        return(NULL)
      }
    }
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

format_candidate_bsmap = function(fas,
                                  tib,
                                  
                                  ids_key,
                                  prb_key,
                                  aln_key,
                                  rev_key,
                                  
                                  verbose=0, vt=3,tc=1,tt=NULL,
                                  funcTag='format_candidate_bsmap') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}       fas = {fas}.{RET}"))
    cat(glue::glue("{mssg}   ids_key = {ids_key}.{RET}"))
    cat(glue::glue("{mssg}   prb_key = {prb_key}.{RET}"))
    cat(glue::glue("{mssg}   aln_key = {aln_key}.{RET}"))
    cat(glue::glue("{mssg}   rev_key = {rev_key}.{RET}"))
    # if (verbose>=vt+4) {
    #   cat(glue::glue("{RET}"))
    #   cat(glue::glue("{mssg}   tib = {tib}.{RET}"))
    # }
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ids_sym <- rlang::sym(ids_key)
    prb_sym <- rlang::sym(prb_key)
    aln_sym <- rlang::sym(aln_key)
    rev_sym <- rlang::sym(rev_key)
    
    # Ensure we have Aln_Prb/Aln_Rev::
    ret_tib <- tib %>%
      dplyr::mutate(!!aln_sym := deMs(!!prb_sym, uc=TRUE),
                    !!rev_sym := revCmp(!!aln_sym) )
    
    fas_vec <- ret_tib %>%
      dplyr::mutate( fas_line=paste0(">",!!ids_sym,"\n", !!aln_sym) ) %>%
      dplyr::pull(fas_line)
    
    safe_write(x=fas_vec, type="line", file=fas,
               funcTag=funcTag, done = TRUE,
               verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Execute BSMAP Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_bsmap = function(ref_fas, 
                     can_fas, 
                     out_dir, 
                     
                     sort   = FALSE, 
                     light  = FALSE,
                     reload = FALSE,
                     
                     bsp_exe,
                     bsp_dir = NULL,
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
      
      if (!is.null(bsp_dir)) {
        bsp_exe_test <- file.path( bsp_dir, base::basename(bsp_exe) )
        if (file.exists(bsp_exe_test)) bsp_exe_test
      }
      
      if (!file.exists(bsp_exe)) {
        cat(glue::glue("{mssg} Warning: Unable to locate bsp_exe={bsp_exe}. ",
                       "Will try docker version: '/repo/BSMAPz/bsmapz'.{RET2}"))
        bsp_exe <- '/repo/BSMAPz/bsmapz'
      }
      
      if (!file.exists(bsp_exe)) {
        stop(glue::glue("{RET}{mssg} ERROR: bsp_exe='{bsp_exe}' does ",
                        "NOT exist!{RET2}"))
        return(NULL)
      }
      
      bsp_cmd <- glue::glue(
        "{bsp_exe} -a {can_fas} -d {ref_fas} {bsp_opt} -o {out_bsp}{RET}",
        "cat {out_bsp} | {add_cmd} gzip -c -> {bsp_tsv}{RET}")
        # "rm -f {out_bsp}{RET}")
      
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
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
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
    
    if (is.null(ret_tib) || length(ret_tib)==0) {
      stop(glue::glue("{RET}{mssg} ERROR: bsmap results are NULL!{RET2}"))
      return(NULL)
    }
    
    # Sort by genomic position::
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(Bsp_Chr, Bsp_Beg)
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n="ret")
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
                      
                      out_csv = NULL,
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
    bsp_cnt <- print_tib(bsp_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=bsp_key)
    
    ord_key <- glue::glue("ord-tib({funcTag})")
    ord_cnt <- print_tib(can_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ord_key)
    
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
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
    
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
    ann_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ann_key)
    
    ret_key <- glue::glue("calculated-fields({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
    
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
    
    out_cnt <- safe_write(x=ret_tib, file=out_csv, funcTag=funcTag,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  ret_tib
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
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# End of file
