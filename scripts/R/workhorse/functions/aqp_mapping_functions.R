
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   1.0 AQP Address Manifest Workflow: 
#                           Order/Match/AQP/PQC
#
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

# NOTE:: spec() will return you with the cols() needed to load a file 
#  rapidly...

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
#                   1.0 AQP Address Manifest Workflow: 
#                           Order/Match/AQP/PQC
#                           Main Workflow Driver
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

aqp_mapping_workflow = function(ord_dat, 
                                mat_dat, 
                                aqp_dat,

                                prb_key = "Ord_Prb",
                                add_key = "Address",
                                des_key = "Ord_Des",
                                din_key = "Ord_Din",
                                ids_key = "Prb_Key",

                                del = "_",
                                
                                out_csv = NULL,
                                out_dir,
                                run_tag, 
                                re_load = FALSE,
                                pre_tag = NULL,
                                end_str = 'csv.gz',
                                sep_chr = '.',
                                out_col = c("Prb_Key","Address","Ord_Des",
                                            "Ord_Din","Ord_Map","Ord_Prb"),

                                verbose=0, vt=3,tc=1,tt=NULL,
                                funcTag='aqp_mapping_workflow') {
  
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
  
  sum_csv <- out_csv %>% 
    stringr::str_remove(paste0(sep_chr,end_str,"$") ) %>%
    paste("summary.csv.gz", sep=sep_chr)
  
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}       del={del}.{RET}"))
    cat(glue::glue("{mssg}   out_dir={out_dir}.{RET}"))
    cat(glue::glue("{mssg}   out_csv={out_csv}.{RET}"))
    cat(glue::glue("{mssg}   sum_csv={sum_csv}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}   prb_key={prb_key}.{RET}"))
    cat(glue::glue("{mssg}   add_key={add_key}.{RET}"))
    cat(glue::glue("{mssg}   des_key={des_key}.{RET}"))
    cat(glue::glue("{mssg}   din_key={din_key}.{RET}"))
    cat(glue::glue("{mssg}   ids_key={ids_key}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  etime   <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL

  stime <- base::system.time({
    
    prb_sym <- rlang::sym(prb_key)
    din_sym <- rlang::sym(din_key)
    ids_sym <- rlang::sym(ids_key)
    aln_vec <- c(add_key,des_key)

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                         Process AQP/PQC Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (purrr::is_character(aqp_dat)) {
      aqp_tib <- load_aqp_files(aqp_dat, 
                                verbose=verbose, vt=vt+1,tc=tc+1,tt=tt) %>% 
        dplyr::arrange(-Ord_Idx) %>% 
        dplyr::distinct(Address, .keep_all=TRUE)
    } else {
      aqp_tib <- aqp_dat
    }
    if (!tibble::is_tibble(aqp_tib)) {
      stop(glue::glue("{RET}[{funcTag}]: AQP/PQC is NOT a tibble!{RET}{RET}"))
      return(ret_tib)
    }
    aqp_dup_cnt <- aqp_tib %>% 
      dplyr::add_count(Address, name="Add_Cnt") %>% 
      dplyr::filter(Add_Cnt!=1) %>% base::nrow()
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: aqp_dup_cnt={aqp_dup_cnt}.{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Process Match Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (purrr::is_character(mat_dat)) {
      mat_tib <- load_aqp_files(mat_dat,
                                verbose=verbose, vt=vt+1,tc=tc+1,tt=tt) %>% 
        dplyr::arrange(-Ord_Idx) %>% 
        dplyr::distinct(Address, .keep_all=TRUE)
    } else {
      mat_tib <- mat_dat
    }
    if (!tibble::is_tibble(mat_tib)) {
      stop(glue::glue("{RET}[{funcTag}]: Match is NOT a tibble!{RET}{RET}"))
      return(ret_tib)
    }
    
    mat_dup_cnt <- mat_tib %>% 
      dplyr::add_count(Address, name="Add_Cnt") %>% 
      dplyr::filter(Add_Cnt!=1) %>% base::nrow()
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: mat_dup_cnt={mat_dup_cnt}.{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Process Order Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (purrr::is_character(ord_dat)) {
      ord_tib <- load_aqp_files(ord_dat,
                                verbose=verbose, vt=vt+1,tc=tc+1,tt=tt) %>% 
        dplyr::mutate(Ord_Map=Ord_Idx+Ord_Map)
    } else {
      ord_tib <- ord_dat
    }
    if (!tibble::is_tibble(ord_tib)) {
      stop(glue::glue("{RET}[{funcTag}]: Order is NOT a tibble!{RET}{RET}"))
      return(ret_tib)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Build Probes/Address Table::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- aqp_tib %>% 
      dplyr::filter(Decode_Status==0) %>%
      dplyr::select(Address,Ord_Idx) %>%
      dplyr::rename(Aqp_Idx=Ord_Idx) %>%
      dplyr::inner_join(mat_tib, by=add_key) %>%
      dplyr::rename(Mat_Idx=Ord_Idx) %>%
      dplyr::rename_with(~ c(prb_key), dplyr::all_of(c("Mat_Prb")) ) %>%
      dplyr::inner_join(ord_tib, by=c(prb_key)) %>%
      tidyr::unite(Tmp_Key, dplyr::all_of(aln_vec), sep=del, remove=FALSE) %>%
      tidyr::unite(!!ids_sym, Tmp_Key,!!din_sym, sep="", remove=FALSE) %>%
      dplyr::select(-Tmp_Key) %>%
      dplyr::mutate(Ord_Cgn=Ord_Key %>% 
                      stringr::str_remove("^[^0-9]+") %>% 
                      stringr::str_remove("[^0-9]+.*$") %>% 
                      stringr::str_remove("^0+") ) %>%
      dplyr::select(Address,Ord_Des,Ord_Din,Ord_Map,Ord_Prb,Ord_Par,
                    Ord_Key,Ord_Col,Ord_Idx,Mat_Idx,Aqp_Idx,Mat_Tan,
                    dplyr::everything())

    prb_dup_cnt <- ret_tib %>% 
      dplyr::add_count(Address, name="Add_Cnt") %>% 
      dplyr::filter(Add_Cnt!=1) %>% base::nrow()
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: prb_dup_cnt={prb_dup_cnt}.{RET}"))
    
    par_dup_cnt <- ret_tib %>% 
      dplyr::add_count(Address,Ord_Prb,Ord_Par, name="Dup_Cnt") %>% 
      dplyr::filter(Dup_Cnt!=1) %>% base::nrow()
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: par_dup_cnt={par_dup_cnt}.{RET}"))
    
    # Build Summary::
    aqp_sum <- aqp_tib %>% 
      dplyr::group_by(Ord_Idx, Decode_Status) %>% 
      dplyr::summarise(Decode_Count=n(), .groups="drop")
    print(aqp_sum, n=base::nrow(aqp_sum))
    
    sum_key <- glue::glue("aqp-summary({funcTag})")
    sum_cnt <- print_tib(aqp_sum,funcTag, verbose,vt=vt+4,tc=tc+1, n=sum_key)
    
    # Format Final Output::
    ret_tib <- ret_tib %>% 
      dplyr::select(dplyr::all_of(out_col), dplyr::everything())
    
    # Write Probe and Summary Output::
    sum_cnt <- safe_write(x=aqp_sum, file=sum_csv, done = TRUE,
                          verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    out_cnt <- safe_write(x=ret_tib, file=out_csv, done = TRUE,
                          verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    tt$addFile(out_csv)
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose, vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Input Parameter Validation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

valid_dir_files = function(dir, csv, type,

                           verbose=0,vt=6,tc=1,tt=NULL,
                           funcTag='valid_dir_files') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
  etime   <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  ret_vec <- NULL
  
  dir_vec <- splitStrToVec(dir)
  csv_vec <- splitStrToVec(csv)
  dir_len <- length(dir_vec)
  csv_len <- length(csv_vec)
  
  if (type == "order" && csv_len==0 && dir_len==1) {
    if (verbose>=vt) cat(glue::glue("{mssg} Search {type} directory = ",
                                    "{dir_vec} for '.csv.gz' files.{RET}"))
    csv_vec <- list.files(dir_vec, pattern = ".csv.gz$")
    csv_len <- length(csv_vec)
  }
  
  if (dir_len == 0) {
    stop(glue::glue(
      "{RET}{mssg} ERROR: Must provide at least one directory for {type}. ",
      "dir_len={dir_len}; csv_len={csv_len}.{RET}"))
    return(ret_tib)
  }
  
  if (csv_len == 0) {
    stop(glue::glue(
      "{RET}{mssg} ERROR: Must provide at least one file name for {type}. ",
      "dir_len={dir_len}; csv_len={csv_len}.{RET}"))
    return(ret_tib)
  }
  
  if (dir_len != 1 && dir_len != csv_len) {
    stop(glue::glue(
      "{RET}{mssg} ERROR: Must provide either a single {type} directory that ",
      "contains all {type} files, or provide an individual {type} directory ",
      "for each {type} file. dir_len={dir_len}; csv_len={csv_len}.{RET}"))
    return(ret_tib)
  }
  ret_vec <- paste(dir_vec, csv_vec, sep="/")
  for (file in ret_vec) {
    if (!file.exists(file)) {
      cat(glue::glue("{RET}{mssg} ERROR: Failed to find {type} ",
                     "file={file}! Please check your inputs.{RET2}"))
      return(ret_tib)
    }
  }
  ret_cnt <- ret_vec %>% length()
  
  ret_str <- paste(ret_vec, sep=",")
  
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET}{RET}"))

  ret_str
}

valid_aqp_inputs = function(ord_dir, ord_csv,
                            mat_dir, mat_tsv,
                            aqp_dir, aqp_tsv,
                            
                            verbose=0,vt=3,tc=1,tt=NULL,
                            funcTag='valid_aqp_inputs') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} File IO Parameters::{RET}"))
    cat(glue::glue("{mssg}   ord_dir={ord_dir}.{RET}"))
    cat(glue::glue("{mssg}   ord_csv={ord_csv}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}   mat_dir={mat_dir}.{RET}"))
    cat(glue::glue("{mssg}   mat_tsv={mat_tsv}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg}   aqp_dir={aqp_dir}.{RET}"))
    cat(glue::glue("{mssg}   aqp_tsv={aqp_tsv}.{RET}"))
    cat(glue::glue("{RET}"))
  }

  etime   <- 0
  ret_cnt <- 0
  ret_dat <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         Validate Order Files::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  type <- "order"
  ord_str <- 
    valid_dir_files(ord_dir, ord_csv, type, verbose=verbose, vt=vt+1,tc=tc+1)
  
  if (is.null(ord_str)) {
    stop(glue::glue("{RET}{mssg} ERROR: Failed type={type}!{RET2}"))
    return(ret_dat)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         Validate Match Files::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  type <- "match"
  mat_str <- 
    valid_dir_files(mat_dir, mat_tsv, type, verbose=verbose, vt=vt+1,tc=tc+1)
  
  if (is.null(mat_str)) {
    stop(glue::glue("{RET}{mssg} ERROR: Failed type={type}!{RET2}"))
    return(ret_dat)
  }

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Validate AQP/PQC Files::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  type <- "AQP/PQC"
  aqp_str <- 
    valid_dir_files(aqp_dir, aqp_tsv, type, verbose=verbose, vt=vt+1,tc=tc+1)

  if (is.null(aqp_str)) {
    stop(glue::glue("{RET}{mssg} ERROR: Failed type={type}!{RET2}"))
    return(ret_dat)
  }
  
  ret_dat$ord_str <- ord_str
  ret_dat$mat_str <- mat_str
  ret_dat$aqp_str <- aqp_str
  
  ret_cnt <- ret_dat %>% length()
  
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET}{RET}"))

  return(ret_dat)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Order, Match, AQP/PQC File IO Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_aqp_files = function(file,
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='load_aqp_files') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   file={file}.{RET}"))
    cat(glue::glue("{RET}"))
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
    } else {
      stop(glue::glue("{RET}{mssg} unrecognized variable={file}.{RET}"))
      return(NULL)
    }
    
    ret_tib <- file_vec %>%
      lapply(load_aqp_file,
             verbose=verbose, vt=vt+1,tc=tc+1,tt=tt) %>% 
      dplyr::bind_rows(.id = "Ord_Idx") %>%
      clean_tibble()
    
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET}{RET}"))
  
  ret_tib
}

load_aqp_file = function(file, 
                         idx=NULL,
                         n_max=100, 
                         guess=1000,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='load_aqp_file') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  # Saftey Check for Empty Files::
  if (is.null(file) || length(file)==0)
    return(base::invisible(NULL))
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     idx = {idx}.{RET}"))
    cat(glue::glue("{mssg}    file = {file}.{RET}"))
    cat(glue::glue("{mssg}   n_max = {n_max}.{RET}"))
    cat(glue::glue("{mssg}   guess = {guess}.{RET}"))
    cat(glue::glue("{RET}"))
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
  
  val_cols$mat4 <- 
    cols(
      address_name = col_integer(),
      probe_id      = col_character(),
      sequence      = col_character(),
      type_b        = col_character(),
      address_name2 = col_integer(),
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
  sel_cols$mat4 <- c("address_name","bo_seq")
  sel_cols$aqp  <- c("Address","Decode_Status")
  sel_cols$pqc  <- c("Address","Status")
  
  key_cols <- list()
  key_cols$ord  <- c("Ord_Key","Ord_PrbA","Ord_PrbB","Ord_Norm")
  key_cols$mat1 <- c("Address","Sequence")
  key_cols$mat2 <- c("Address","Sequence")
  key_cols$mat3 <- c("Address","Sequence")
  key_cols$mat4 <- c("Address","Sequence")
  key_cols$aqp  <- c("Address","Decode_Status")
  key_cols$pqc  <- c("Address","Decode_Status")
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    guess_tib <- guess_aqp_file( file, n_max=n_max,
                                 verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
    
    row_num = guess_tib$row_num[1]
    col_num = guess_tib$col_num[1]
    del_key = guess_tib$del_key[1]
    beg_key = guess_tib$beg_key[1]
    
    dat_key <- NULL
    val_col <- NULL
    sel_col <- NULL
    key_col <- NULL
    if (is.null(beg_key) || is.null(col_num)) {
      cat(glue::glue("{RET}[{funcTag}]: ERROR: Either beg_key OR col_num is ",
                      "NULL!{RET}"))
      cat(glue::glue("{RET}[{funcTag}]: ERROR: file='{file}'{RET}"))
      cat(glue::glue("[{funcTag}]: ERROR: guess_tib={RET}"))
      print(guess_tib)
      stop(glue::glue("{RET2}[{funcTag}]: ERROR: Either beg_key OR col_num is ",
                      "NULL!{RET}{TAB}file='{file}'!{RET2}"))
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
    } else if (beg_key==names(val_cols$mat4$cols)[1] && 
               col_num==length(val_cols$mat4$cols)) {
      dat_key <- "mat4"
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
    ret1_key <- glue::glue("ret-tib1({funcTag})")
    ret1_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+6,tc=tc+1, n=ret1_key)
    
    # More goofy format verification::
    if (dat_key=="mat2" || dat_key=="mat4") {
      if (dat_key=="mat2") {
        add_mat_tib <- ret_tib %>% dplyr::filter(address_names != address_name)
        add_mat_cnt <- add_mat_tib %>% base::nrow()
      }
      if (dat_key=="mat4") {
        add_mat_tib <- ret_tib %>% dplyr::filter(address_name != address_name2)
        add_mat_cnt <- add_mat_tib %>% base::nrow()
      }
      if (add_mat_cnt != 0) {
        cat(glue::glue("{RET}[{funcTag}]: ERROR: Goofy match format ",
                       "({dat_key}) has non-zero ({add_mat_cnt}) address ",
                       "miss-matching! add_mat_tib={RET}"))
        print(add_mat_tib)
        stop(glue::glue("{RET}[{funcTag}]: ERROR: Goofy match format ",
                        "({dat_key}) has non-zero ({add_mat_cnt}) address ",
                        "miss-matching!{RET2}"))
        return(NULL)
      }
      if (verbose>=vt+6)
        cat(glue::glue("{RET}[{funcTag}]: Success: Goofy match format ",
                       "({dat_key}) has zero ({add_mat_cnt}) address ",
                       "miss-matching!{RET2}"))
    }

    # Apply Selected Columns::
    ret_tib  <- ret_tib %>% dplyr::select(dplyr::all_of(sel_col))
    ret2_key <- glue::glue("ret-tib2({funcTag})")
    ret2_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+6,tc=tc+1, n=ret2_key)
    
    # Apply Canonical-Renamed Columns::
    ret_tib  <- ret_tib %>% purrr::set_names(key_col)
    ret3_key <- glue::glue("ret-tib2({funcTag})")
    ret3_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+6,tc=tc+1, n=ret3_key)
    
    # Apply Formatting Rules:: Order Data Type
    if (dat_key=="ord") {
      ret_cnt <- ret_tib %>% base::nrow()
      ord_tib <- ret_tib %>% 
        dplyr::distinct() %>%
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
        ) %>% 
        dplyr::distinct(Ord_PrbA,Ord_PrbB, .keep_all=TRUE) %>%
        dplyr::mutate(Ord_Row=dplyr::row_number(),
                      Ord_Map=as.double(Ord_Row/ret_cnt))
      
      ord_cols <- c("Ord_Key","Ord_Des", "Ord_Prb", "Ord_Par", "Ord_Col","Ord_Map")
      ord_colA <- c("Ord_Key","Ord_DesA","Ord_PrbA","Ord_PrbB","Ord_Col","Ord_Map")
      ord_colB <- c("Ord_Key","Ord_DesB","Ord_PrbB","Ord_PrbA","Ord_Col","Ord_Map")
      
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
          Ord_Map,Ord_Key,Ord_Des,Ord_Din,Ord_Col,dplyr::everything()
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
                        "but format doesn't match; trim_cnt={trim_cnt}!!!{RET2}"))
        return(NULL)
      }
    }
    
    # Apply Formatting Rules:: Match Files:: Sequence
    if (stringr::str_starts(dat_key,"mat")) {
      ret_tib <- ret_tib %>% 
        dplyr::mutate(Sequence=stringr::str_to_upper(Sequence),
                      Mat_Tan=stringr::str_sub(Sequence,1,45),
                      Mat_Prb=stringr::str_sub(Sequence,46)
        ) %>%
        dplyr::select(-Sequence) %>%
        dplyr::filter(!is.na(Address)) %>% 
        dplyr::distinct(Address,Mat_Prb, .keep_all=TRUE)
      
      # Verify the goofy naming schemes::
      add_mat_cnt <- 0
      add_mat_tib <- NULL
    }
    
    # Add Dat_IDX if provided
    if (!is.null(idx))
      ret_tib <- ret_tib %>%
      dplyr::mutate(Dat_IDX=as.integer(idx))
    
    # Standard Column Data Type Clean-Up::
    ret_tib <- ret_tib %>% clean_tibble()
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET}{RET}"))
  
  ret_tib
}

guess_aqp_file = function(file, 
                          fields=c("Assay_Design_Id","Plate",
                                   "address_names","address_name",
                                   "Address", "probe_id"),
                          cols=NULL,
                          n_max=100,
                          verbose=0,vt=4,tc=1,tt=NULL,
                          funcTag='guess_aqp_file') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+4) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}     file={file}.{RET}"))
    cat(glue::glue("{mssg}   fields={fields}.{RET}"))
    cat(glue::glue("{mssg}     cols={cols}.{RET}"))
    cat(glue::glue("{mssg}    n_max={n_max}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_val <- NULL
  stime <- base::system.time({
    
    file_del_str <- guess_file_del(file, n_max=n_max,
                                   verbose=verbose, vt=vt+4,tc=tc+1)
    if (is.null(file_del_str)) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: file_del_str=NULL!!!{RET}{RET}"))
      return(ret_tib)
    }
    if (verbose>=vt)
      cat(glue::glue("{mssg} file_del_strfile={file_del_str}.{RET}"))

    dat_tib <- readr::read_lines(file, n_max=n_max) %>%
      tibble::as_tibble() %>% 
      dplyr::mutate(row_num=dplyr::row_number())
    dat_cnt <- print_tib(dat_tib,funcTag, verbose,vt=vt+6,tc=tc, n="dat_tib")
    
    ret_tib <- NULL
    for (field in fields) {
      if (verbose>=vt+6)
        cat(glue::glue("{mssg} field={field}{RET}"))
      
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
        cat(glue::glue("{mssg} NEXT...{RET}{RET}"))
    }
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Manifest Mutation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: Pretty Sure this can be removed::
mutate_probe_id = function(tib, 
                           
                           pid="Probe_ID", 
                           cgn="Imp_Cgn_Seq",
                           des="Ord_Des",  
                           din="Ord_Din",
                           tb="Imp_TB_Seq", 
                           co="Imp_CO_Seq",
                           inf="Infinium_Design",
                           
                           pad=8,
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutate_probe_id'
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
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Noob-Mask Manifest Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# VERY HACKY QUICK FIX FOR NOOB CONTROLS
# TBD:: Make this clean!!!
#
# noob_mask(noob_csv = opt$noob, ctl_csv = opt$ses_ctl_csv, verbose=opt$verbose, tt=pTracker)
noob_mask = function(noob_csv,
                     ctl_csv = NULL,
                     
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='noob_mask') {
  
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
    
    if (!is.null(ctl_csv) && file.exists(ctl_csv))
      ctl_tib <- safe_read(ctl_csv, verbose=verbose)
    
    # pqc_tib  <- load_aqp_files(opt$aqps, verbose=opt$verbose)
    noob_tib <- safe_read(noob_csv, verbose=verbose) %>%
      dplyr::mutate(
        Probe_Design=dplyr::case_when(
          DESIGN=="I"  ~ 1,
          DESIGN=="II" ~ 2,
          TRUE ~ NA_real_
        )
      )
    
    # pqc_pas_tib <- dplyr::filter(pqc_tib,Decode_Status==0)
    # pqc_mis_tib <- dplyr::filter(pqc_tib,Decode_Status!=0)
    
    # noob_mis_tib <- dplyr::bind_rows(
    #   noob_tib %>% dplyr::filter(!is.na(U)) %>% dplyr::filter(U %in% pqc_mis_tib$Address),
    #   noob_tib %>% dplyr::filter(!is.na(M)) %>% dplyr::filter(M %in% pqc_mis_tib$Address)
    # ) %>% dplyr::distinct(Probe_ID, .keep_all = TRUE)
    
    # noob_pas_tib <- noob_tib %>% dplyr::filter(!Probe_ID %in% noob_mis_tib$Probe_ID)
    
    # noob_man_tib <- dplyr::bind_rows(
    #   noob_pas_tib %>% dplyr::filter(!is.na(U)) %>% dplyr::filter(U %in% pqc_pas_tib$Address),
    #   noob_pas_tib %>% dplyr::filter(!is.na(M)) %>% dplyr::filter(M %in% pqc_pas_tib$Address)
    # ) %>% dplyr::distinct(Probe_ID, .keep_all = TRUE)
    
    # noob_man_sum <- noob_man_tib %>% 
    #   dplyr::group_by(SubSet,DESIGN) %>% 
    #   dplyr::summarise(Count=n(), .groups="drop")
    
    # For everything::
    # noob_man_tib <- noob_tib
    
    noob_sel_tib <- noob_tib %>%
      dplyr::filter(SubSet) %>%
      dplyr::select(-Probe_ID,-SubSet) %>%
      dplyr::rename(Probe_ID=Hash_ID) %>%
      dplyr::select(Probe_ID, dplyr::everything()) %>%
      clean_tibble()
    
    # All Noob & GS Controls::
    #
    ret_tib <- dplyr::bind_rows(noob_sel_tib,ctl_tib)

    # verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
# TBD:: These first three manifest related functions need to be cleaned up
#   and re-written and moved to manfiest functions()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          AQP to Sesame Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

aqp_to_man = function(tib,
                      verbose=0,vt=3,tc=1,tt=NULL,
                      funcTag='aqp_to_man') {
  
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
    
    # TBD::
    #
    #  * Investigate BSP alignments::
    #    - Calculate number of hits good/bad
    #    - Calculate extension/color distribution
    #    - Extract BSC Top/Probe-Design
    #
    
    # Workflow::
    #
    #  - Split 2,U,M
    #  - Join U/M
    #  - 
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET}{RET}"))
  
  ret_tib
}

aqp_to_sesame1 = function(tib, isMU=FALSE, retData=FALSE,
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='aqp_to_sesame1') {
  
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
  ret_dat <- NULL
  stime <- base::system.time({
    
    #
    # Infinium I::
    #
    lab_tib <- tib
    
    if (!opt$multi_unique)
      # lab_tib <- lab_tib %>%
      # dplyr::filter(Ord_Din=="cg" & 
      #                 Bsp_Din_Ref_U=="CG" & Bsp_Din_Bsc_U=="TG" &
      #                 Bsp_Din_Ref_M=="CG" & Bsp_Din_Bsc_M=="CG")
      
      lab_tib <- lab_tib %>%
        dplyr::distinct(Ord_Key,Address_U,Address_M,Cgn_Str,Strand_TB_U,Bsp_CO_U,Ord_Prb_U,Ord_Prb_M) %>%
        dplyr::mutate(
          Cgn_Info=paste0(Strand_TB_U,Bsp_CO_U,"1"),
          Probe_ID=paste(Cgn_Str,Cgn_Info, sep="_")
        ) %>%
        dplyr::group_by(Probe_ID) %>%
        dplyr::mutate(Rep_Num=dplyr::row_number(),
                      Probe_ID=paste0(Probe_ID,Rep_Num)) %>%
        dplyr::ungroup() %>% 
        dplyr::select(Ord_Key,Probe_ID,Rep_Num,
                      Address_U,Ord_Prb_U,
                      Address_M,Ord_Prb_M)
    
    cor_tib <- tib %>% 
      dplyr::mutate(
        Name=Cgn_Str,
        U=Address_U,
        AlleleA_ProbeSeq=Ord_Prb_U,
        M=Address_M,
        AlleleB_ProbeSeq=Ord_Prb_M,
        Next_Base=Next_Base_U,
        Color_Channel=dplyr::case_when(
          Next_Base_U=="C" | Next_Base_U=="G" ~ "Grn",
          Next_Base_U=="A" | Next_Base_U=="T" ~ "Red",
          TRUE ~ NA_character_
        ),
        Col=stringr::str_sub(Color_Channel, 1,1),
        Probe_Type=Ord_Din,
        Strand_FR=Bsp_FR_U,
        Strand_TB=Strand_TB_U,
        Strand_CO=Bsp_CO_U,
        Infinium_Design_Type="1",
        CHR=Bsp_Chr,
        MAPINFO=Bsp_Pos,
        Species=opt$Species,
        Genome_Build=opt$genBuild,
        Source_Seq=Probe_Seq_T_U,
        Underlying_CpG_Count=Cpg_Cnt,
        Underlying_CpG_Min_Dist=Cpg_Dis_M,
        Alt_Cgn_Count_U=Alt_Cgn_Cnt_U,
        Alt_Cgn_Count_M=Alt_Cgn_Cnt_M
      ) %>%
      dplyr::select(Ord_Key,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag_U,Bsp_Tag_M)
    
    all_tib <- cor_tib %>%
      dplyr::inner_join(lab_tib,
                        by=c("Ord_Key",
                             "U"="Address_U","M"="Address_M",
                             "AlleleA_ProbeSeq"="Ord_Prb_U",
                             "AlleleB_ProbeSeq"="Ord_Prb_M") ) %>%
      dplyr::select(Probe_ID,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,Rep_Num,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag_U,Bsp_Tag_M) %>%
      clean_tibble()
    
    ret_tib <- dplyr::right_join(cor_tib,lab_tib,
                                 by=c("Ord_Key",
                                      "U"="Address_U","M"="Address_M",
                                      "AlleleA_ProbeSeq"="Ord_Prb_U",
                                      "AlleleB_ProbeSeq"="Ord_Prb_M") )
    if (!opt$multi_unique)
      ret_tib <- ret_tib %>%
      dplyr::mutate(
        CHR=dplyr::case_when(
          Bsp_Tag_U=="UM" && Bsp_Tag_M=="UM" ~ CHR,
          TRUE ~ "0"
        ),
        MAPINFO=dplyr::case_when(
          Bsp_Tag_U=="UM" && Bsp_Tag_M=="UM" ~ MAPINFO,
          TRUE ~ as.integer(0)
        ),
        Strand_FR=dplyr::case_when(
          Bsp_Tag_U=="UM" && Bsp_Tag_M=="UM" ~ Strand_FR,
          TRUE ~ "0"
        )
      )
    
    ret_tib <- ret_tib %>%
      dplyr::select(Probe_ID,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,Rep_Num,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag_U,Bsp_Tag_M) %>%
      dplyr::distinct() %>%
      dplyr::distinct(Probe_ID,Name,
                      U,AlleleA_ProbeSeq,
                      M,AlleleB_ProbeSeq,
                      # Next_Base, 
                      .keep_all = TRUE) %>%
      dplyr::distinct(U,M, .keep_all = TRUE) %>%
      clean_tibble()
    
    # lab_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "rs"))
    # all_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "rs"))
    # ret_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "rs"))
    
    # These should be zero::
    col_cnt1 <- ret_tib %>% dplyr::filter(is.na(Col)) %>% base::nrow()
    col_cnt2 <- ret_tib %>% dplyr::filter(is.na(Color_Channel)) %>% base::nrow()
    max_repn <- max(ret_tib$Rep_Num)
    cat(glue::glue("{mssg} col_cnt1={col_cnt1}{RET}"))
    cat(glue::glue("{mssg} col_cnt2={col_cnt2}{RET}"))
    cat(glue::glue("{mssg} max_repn={max_repn}{RET}"))
    
    ret_tib %>% dplyr::group_by(Col) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print()
    ret_tib %>% dplyr::group_by(Color_Channel) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print()
    ret_tib %>% dplyr::filter(stringr::str_starts(Probe_ID,"rs")) %>% print()
    
    if (retData) ret_dat$all_tib <- all_tib
    if (retData) ret_dat$ses_tib <- ret_tib
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET}{RET}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
}

aqp_to_sesame2 = function(tib, isMU=FALSE, retData=FALSE,
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='aqp_to_sesame2') {
  
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
  ret_dat <- NULL
  stime <- base::system.time({
    
    # Ord_Din_sym <- rlang::sym(Ord_Din)
    
    lab_tib <- tib
    # if (!isMU)
    #   lab_tib <- lab_tib %>%
    #     dplyr::filter(Bsp_Din_Ref=="CG" & Bsp_Din_Bsc=="YG")
    
    lab_tib <- lab_tib %>%
      dplyr::distinct(Ord_Key,Address,Cgn_Str,Bsp_CO,Ord_Prb, .keep_all = TRUE) %>% 
      dplyr::select(Ord_Key,Address,Cgn_Str,Strand_TB,Bsp_CO,Ord_Prb) %>%
      dplyr::mutate(
        Cgn_Info=paste0(Strand_TB,Bsp_CO,"2"),
        Probe_ID=paste(Cgn_Str,Cgn_Info, sep="_")
      ) %>%
      dplyr::group_by(Probe_ID) %>%
      dplyr::mutate(Rep_Num=dplyr::row_number(),
                    Probe_ID=paste0(Probe_ID,Rep_Num)) %>%
      dplyr::ungroup() %>% 
      dplyr::select(Ord_Key,Probe_ID,Rep_Num,
                    Address,Ord_Prb)
    
    cor_tib <- tib %>% 
      dplyr::mutate(
        Name=Cgn_Str,
        U=Address,
        AlleleA_ProbeSeq=Ord_Prb,
        M="",
        AlleleB_ProbeSeq="",
        Color_Channel="Both",
        Col="",
        Probe_Type=Ord_Din,
        Strand_FR=Bsp_FR,
        Strand_CO=Bsp_CO,
        Infinium_Design_Type="2",
        CHR=Bsp_Chr,
        MAPINFO=Bsp_Pos,
        Species=opt$Species,
        Genome_Build=opt$genBuild,
        Source_Seq=Probe_Seq_T,
        Underlying_CpG_Count=Cpg_Cnt,
        Underlying_CpG_Min_Dist=Cpg_Dis,
        Alt_Cgn_Count_U=Alt_Cgn_Cnt,
        Alt_Cgn_Count_M=0
      ) %>%
      dplyr::select(Ord_Key,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag) %>%
      dplyr::rename(Bsp_Tag_U=Bsp_Tag) %>%
      dplyr::mutate(Bsp_Tag_M="")
    
    all_tib <- cor_tib %>%
      dplyr::inner_join(lab_tib,
                        by=c("Ord_Key","U"="Address",
                             "AlleleA_ProbeSeq"="Ord_Prb") ) %>%
      dplyr::select(Probe_ID,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,Rep_Num,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag_U,Bsp_Tag_M) %>%
      clean_tibble()
    
    ses_tib <- dplyr::right_join(cor_tib,lab_tib, 
                                 by=c("Ord_Key","U"="Address","AlleleA_ProbeSeq"="Ord_Prb"))
    
    if (!isMU)
      ses_tib <- ses_tib %>%
      dplyr::mutate(
        CHR=dplyr::case_when(
          Bsp_Tag_U=="UM" ~ CHR,
          TRUE ~ "0"
        ),
        MAPINFO=dplyr::case_when(
          Bsp_Tag_U=="UM" ~ MAPINFO,
          TRUE ~ as.integer(0)
        ),
        Strand_FR=dplyr::case_when(
          Bsp_Tag_U=="UM" ~ Strand_FR,
          TRUE ~ "0"
        )
      )
    
    ret_tib <- ses_tib %>%
      dplyr::select(Probe_ID,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                    Next_Base,Color_Channel,Col,Probe_Type,
                    Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,Rep_Num,
                    CHR,MAPINFO,Species,Genome_Build,
                    Source_Seq,Forward_Sequence,Top_Sequence,
                    Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                    Alt_Cgn_Count_U,Alt_Cgn_Count_M, Bsp_Tag_U,Bsp_Tag_M) %>% 
      dplyr::distinct() %>%
      dplyr::distinct(Probe_ID,Name,U,AlleleA_ProbeSeq, .keep_all = TRUE) %>%
      dplyr::distinct(U, .keep_all = TRUE) %>%
      clean_tibble()
    
    if (retData) ret_dat$all_tib <- all_tib
    if (retData) ret_dat$ses_tib <- ret_tib
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET}{RET}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#              These functions should probably be deprecated
#                        and moved to function graveyard!
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Order/Match/AQP/PQC Expected Columns::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# NEW SOLUTION FOR FILE HEADERS::
#  - Store all file headers in a single rData (rds) file. They're 
#    minature, but required for files without headers and they 
#    massively speed up reading of very large data files since
#    variable types can be expected!
#
# TBD:: This is copied in the code, so it should be removed from here
#   or moved to a config file... Maybe a single r-data structure (RDS)
#   for manifest parameter defaults is the way to go...
#

# par_cols <- list()
# par_cols$ord <- 
#   cols(
#     Assay_Design_Id        = col_character(),
#     AlleleA_Probe_Id       = col_character(),
#     AlleleA_Probe_Sequence = col_character(),
#     AlleleB_Probe_Id       = col_character(),
#     AlleleB_Probe_Sequence = col_character(),
#     Normalization_Bin      = col_character()
#   )
# 
# par_cols$mat <- 
#   cols(
#     Plate    = col_character(),
#     Row      = col_character(),
#     Col      = col_integer(),
#     Address  = col_integer(),
#     Mod5     = col_character(),
#     Sequence = col_character(),
#     Mod3     = col_character(),
#     Comments = col_character()
#   )
# 
# par_cols$ma2 <- 
#   cols(
#     address_names = col_integer(),
#     probe_id      = col_character(),
#     sequence      = col_character(),
#     type_b        = col_character(),
#     address_name  = col_integer(),
#     bo_seq        = col_character()
#   )
# 
# par_cols$aqp <- 
#   cols(
#     Address           = col_integer(),
#     Decode_Status     = col_integer(),
#     Decode_Error_Code = col_integer(),
#     Decode_Score      = col_integer(),
#     Func_Status       = col_integer(),
#     Func_Error_Code   = col_integer(),
#     QC_Action         = col_integer()
#   )
# 
# par_cols$pqc <- 
#   cols(
#     Address      = col_integer(),
#     Status       = col_integer(),
#     Eval_Code    = col_integer(),
#     Average_Rep  = col_integer(),
#     Expected_Rep = col_integer()
#   )
# 
# par$ord_col <- par_cols$ord$cols %>% names()
# 
# par$mat_col <- par_cols$mat$cols %>% names()
# par$ma2_col <- par_cols$ma2$cols %>% names()
# 
# par$aqp_col <- par_cols$aqp$cols %>% names()
# par$pqc_col <- par_cols$pqc$cols %>% names()
#

# End of file
