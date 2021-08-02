
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
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Manifest Comparison Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

manifest_column_agreement = function(a, b, source_a, source_b,
                                     
                                     join_cols = c("Probe_ID","Prb_Des","Prb_Din"),
                                     
                                     verbose=0,vt=3,tc=1,tt=NULL,
                                     funcTag='manifest_column_agreement') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}    source_a={source_a}.{RET}"))
    cat(glue::glue("{mssg}    source_b={source_b}.{RET}"))
    cat(glue::glue("{mssg}   join_cols={join_cols}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    int_cols <- 
      intersect( names(a), names(b) )
    
    dif_cols <-
      setdiff( names(a), names(b) )
    join_cols <- intersect(join_cols, int_cols)
    
    if (verbose>=vt) {
      cat(glue::glue("{mssg} SetDiff Columns::{RET}"))
      cat(glue::glue("{mssg}   dif_cols={dif_cols}.{RET}"))
      cat(glue::glue("{RET}"))
      cat(glue::glue("{mssg} Intersect Columns::{RET}"))
      cat(glue::glue("{mssg}   dif_cols={dif_cols}.{RET}"))
      cat(glue::glue("{RET}"))
    }
    
    int_tibs <- dplyr::inner_join(
      a %>% dplyr::select(dplyr::all_of(int_cols)),
      b %>% dplyr::select(dplyr::all_of(int_cols)),
      by=join_cols,
      suffix=c("_a","_b")
    )
    int_key <- glue::glue("int-tib({funcTag})")
    int_cnt <- print_tib(int_tibs,funcTag, verbose,vt=vt+4,tc=tc+1, n=int_key)

    ret_tib <- NULL
    for (col_key in int_cols) {
      if (col_key %in% join_cols) next
      
      if (verbose>=vt+1)
        cat(glue::glue("{mssg} Comparing column = {col_key}...{RET}"))
      
      colA_key <- paste0(col_key,"_a")
      colB_key <- paste0(col_key,"_b")
      
      colA_sym <- rlang::sym(colA_key)
      colB_sym <- rlang::sym(colB_key)
      
      cur_tib <- int_tibs %>% 
        dplyr::select(dplyr::all_of( c(colA_key,colB_key) ) ) %>%
        dplyr::summarise(
          Mat_Cnt=sum(!!colA_sym == !!colB_sym, na.rm = TRUE),
          Mis_Cnt=sum(!!colA_sym != !!colB_sym, na.rm = TRUE),
          Total=sum(!is.na(!!colA_sym) & !is.na(!!colB_sym) ),
          Mat_Per=round(100*Mat_Cnt/Total, 2),
          Mis_Per=round(100*Mis_Cnt/Total, 2)
        ) %>% 
        dplyr::mutate(Column = col_key,
                      Source_A = source_a,
                      Source_B = source_b)
      
      cur_key <- glue::glue("cur-key({funcTag})")
      cur_cnt <- print_tib(cur_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=cur_key)

      ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Sesame Manifest IO Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_sesame_repo_address = function(
  name,
  add_decoy = FALSE,
  add_masks = FALSE,
  normalize = TRUE,
  old_cols = c("Probe_ID", "seqnames", "start", "end", "strand", "designType", 
               "channel", "nextBase", "nextBaseRef", "probeType", "gene",
               "gene_HGNC"),
  new_cols = c("Probe_ID", "Chromosome", "Coordinate", "CoordinateG", 
               "Strand_FR", "Infinium_Design_Type", "Color", "Prb_Nxb", 
               "Prb_Nxb_Ref", "Prb_Din", "gene", "gene_HGNC"),
  
  verbose=0,vt=3,tc=1,tt=NULL,
  funcTag='load_sesame_repo_manifest') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}        name={name}.{RET}"))
    cat(glue::glue("{mssg}   add_decoy={add_decoy}.{RET}"))
    cat(glue::glue("{mssg}   add_masks={add_masks}.{RET}"))
    cat(glue::glue("{RET}"))
    if (verbose>=vt+4) {
      cat(glue::glue("{mssg} Old Sesame Columns::{RET}"))
      cat(glue::glue("{mssg}    old_cols={old_cols}.{RET}"))
      cat(glue::glue("{RET}"))
      cat(glue::glue("{mssg} New Sesame Columns::{RET}"))
      cat(glue::glue("{mssg}    new_cols={new_cols}.{RET}"))
      cat(glue::glue("{RET}"))
    }
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # Cache Sesame Manifest Data::
    data_key <- name %>% stringr::str_remove("\\..*$")
    sesameData::sesameDataCache(data_key)
    
    man_tib <- sesameData::sesameDataGet( name ) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var="Probe_ID") %>% 
      tibble::as_tibble()
    man_key <- glue::glue("man-tib({funcTag})")
    man_cnt <- print_tib(man_tib,funcTag, verbose, vt=vt+4,tc=tc+1, n=man_key)
    
    mask_cols <- NULL
    if (add_masks) mask_cols <- man_tib %>% 
      dplyr::select( dplyr::starts_with("MASK_") ) %>% names()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Extract/Format Probe A::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    prb_key <- "A"
    old_prbA_cols <- man_tib %>% 
      dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>% 
      dplyr::select(-dplyr::starts_with("wDecoy_")) %>% 
      names()
    new_prbA_cols <- old_prbA_cols %>% 
      stringr::str_remove(paste0("_",prb_key,"$")) %>% 
      stringr::str_to_title() %>%
      paste("Prb",., sep="_") %>%
      stringr::str_replace("Prb_Address", "Address") %>%
      stringr::str_replace("Prb_Probeseq", "Prb_Seq") %>%
      stringr::str_replace("Prb_Chrm", "Prb_Chr")
    
    old_decoyA_cols <- NULL
    new_decoyA_cols <- NULL
    if (add_decoy) {
      old_decoyA_cols <- man_tib %>%
        dplyr::select(dplyr::starts_with("wDecoy_")) %>%
        dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>%
        names()
      
      new_decoyA_cols <- old_decoyA_cols %>%
        stringr::str_remove(paste0("_",prb_key,"$")) %>% 
        stringr::str_remove("^wDecoy_") %>% 
        stringr::str_to_title() %>%
        stringr::str_replace("Chrm", "Chr") %>%
        paste("Prb_Decoy",., sep="_")
    }
    
    oldA_cols <- c(old_cols, old_prbA_cols, old_decoyA_cols, mask_cols)
    newA_cols <- c(new_cols, new_prbA_cols, new_decoyA_cols, mask_cols)
    
    oldA_cnt <- oldA_cols %>% length()
    newA_cnt <- newA_cols %>% length()
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} Old Sesame Columns:: Probe A({oldA_cnt}).{RET}"))
    if (verbose>=vt+4)
      cat(glue::glue("{mssg}    oldA_cols={oldA_cols}.{RET2}"))
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} New Sesame Columns:: Probe A({newA_cnt}).{RET}"))
    if (verbose>=vt+4)
      cat(glue::glue("{mssg}    newA_cols={newA_cols}.{RET2}"))
    
    addA_tib <- NULL
    addA_tib <- man_tib %>% 
      dplyr::select(dplyr::all_of(oldA_cols)) %>%
      purrr::set_names(newA_cols) %>% 
      dplyr::mutate(Prb_Des=dplyr::case_when(
        Infinium_Design_Type=="I"  ~ "U",
        Infinium_Design_Type=="II" ~ "2",
        TRUE ~ NA_character_)
      )
    
    addA_cnt <- addA_tib %>% base::nrow()
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} Probe Set A: I/U Count={addA_cnt}.{RET2}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Extract/Format Probe B::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    prb_key <- "B"
    old_prbB_cols <- man_tib %>% 
      dplyr::filter(designType=="I") %>%
      dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>% 
      dplyr::select(-dplyr::starts_with("wDecoy_")) %>% 
      names()
    new_prbB_cols <- old_prbB_cols %>% 
      stringr::str_remove(paste0("_",prb_key,"$")) %>% 
      stringr::str_to_title() %>%
      paste("Prb",., sep="_") %>%
      stringr::str_replace("Prb_Address", "Address") %>%
      stringr::str_replace("Prb_Probeseq", "Prb_Seq") %>%
      stringr::str_replace("Prb_Chrm", "Prb_Chr")
    
    old_decoyB_cols <- NULL
    new_decoyB_cols <- NULL
    if (add_decoy) {
      old_decoyB_cols <- man_tib %>%
        dplyr::filter(designType=="I") %>%
        dplyr::select(dplyr::starts_with("wDecoy_")) %>%
        dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>%
        names()
      
      new_decoyB_cols <- old_decoyB_cols %>%
        stringr::str_remove(paste0("_",prb_key,"$")) %>% 
        stringr::str_remove("^wDecoy_") %>% 
        stringr::str_to_title() %>%
        stringr::str_replace("Chrm", "Chr") %>%
        paste("Prb_Decoy",., sep="_")
    }
    
    oldB_cols <- c(old_cols, old_prbB_cols, old_decoyB_cols, mask_cols)
    newB_cols <- c(new_cols, new_prbB_cols, new_decoyB_cols, mask_cols)
    
    oldB_cnt <- oldB_cols %>% length()
    newB_cnt <- newB_cols %>% length()
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} Old Sesame Columns:: Probe B({oldB_cnt}){RET}"))
    if (verbose>=vt+4)
      cat(glue::glue("{mssg}    oldB_cols={oldB_cols}.{RET2}"))
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} New Sesame Columns:: Probe B({newB_cnt}){RET}"))
    if (verbose>=vt+4)
      cat(glue::glue("{mssg}    newB_cols={newB_cols}.{RET2}"))
    
    addB_tib <- NULL
    addB_tib <- man_tib %>% 
      dplyr::filter(designType=="I") %>%
      dplyr::select(dplyr::all_of(oldB_cols)) %>%
      purrr::set_names(newB_cols) %>% 
      dplyr::mutate(Prb_Des="M")
    
    addB_cnt <- addB_tib %>% base::nrow()
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} Probe Set B: M Count={addB_cnt}.{RET2}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Stack Probe Designs:: A/B
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- dplyr::bind_rows(addA_tib,addB_tib) %>%
      dplyr::mutate(
        Strand_FR=dplyr::case_when(
          Strand_FR=="+" ~ "F",
          Strand_FR=="-" ~ "R",
          TRUE ~ NA_character_),
        Prb_Nxb=stringr::str_remove_all(Prb_Nxb,"[^a-zA-Z]") %>% mapDIs(),
        Prb_Inf=dplyr::case_when(
          Infinium_Design_Type=="I"  ~ 1,
          Infinium_Design_Type=="II" ~ 2,
          TRUE ~ NA_real_
        ) %>% as.integer(),
        Color=dplyr::case_when(
          Color=="Both" ~ NA_character_,
          TRUE ~ Color)
      )  %>% dplyr::arrange(Probe_ID)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Summary and Sanity Checks::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # This better be zero::
    dup_cnt <- ret_tib %>% 
      dplyr::arrange(Address) %>% 
      dplyr::add_count(Address, name="Add_Cnt") %>% 
      dplyr::filter(Add_Cnt!=1) %>% base::nrow()
    if (dup_cnt != 0) {
      stop(glue::glue("{mssg} ERROR: Duplicates Found = {dup_cnt}!{RET2}"))
      return(NULL)
    }
    if (verbose>=vt+1)
      cat(glue::glue("{mssg} No duplicates detected = {dup_cnt}.{RET}"))
    
    if (verbose>=vt+2) {
      cat(glue::glue("{mssg} Probe Design Coverage={RET}"))
      sum_tib <- ret_tib %>% 
        dplyr::group_by(Prb_Des) %>% 
        dplyr::summarise(Count=n(), .groups = "drop")
      sum_tib %>% print(n=base::nrow(sum_tib))
      cat(glue::glue("{RET}"))
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Manifest I/O::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_genome_studio_address = function(
  file,
  load_clean     = TRUE,
  load_controls  = FALSE,
  write_clean    = TRUE,
  overwrite      = FALSE,
  add_annotation = FALSE,
  ret_data       = FALSE,
  old_cols = c("IlmnID", "Name", 
               "AddressA_ID", "AlleleA_ProbeSeq", 
               "AddressB_ID", "AlleleB_ProbeSeq", "Next_Base",
               "Color_Channel", "CHR", "MAPINFO", "Strand"),
  new_cols = c("Probe_ID", "Prb_Cgn", 
               "Address_A", "Prb_Seq_A", "Address_B", "Prb_Seq_B",
               "Prb_Nxb", "Color", "Chromosome", "Coordinate", "Strand_FR"),
  non_ann_cols = c("Infinium_Design_Type", "Forward_Sequence",
                   "Genome_Build", "SourceSeq"),
  
  verbose=0, vt=3,tc=1,tt=NULL,
  funcTag='load_genome_studio_address') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}          file={file}.{RET}"))
    cat(glue::glue("{mssg}    load_clean={load_clean}.{RET}"))
    cat(glue::glue("{mssg}   write_clean={write_clean}.{RET}"))
    cat(glue::glue("{mssg}     overwrite={overwrite}.{RET}"))
    cat(glue::glue("{mssg}      ret_data={ret_data}.{RET}"))
    cat(glue::glue("{RET}"))
    if (verbose>=vt+4) {
      old_cnt <- old_cols %>% length()
      new_cnt <- new_cols %>% length()
      cat(glue::glue("{mssg} Old Genome Studio Columns({old_cnt})::{RET}"))
      cat(glue::glue("{mssg}    old_cols={old_cols}.{RET}"))
      cat(glue::glue("{RET}"))
      cat(glue::glue("{mssg} New Genome Studio Columns({new_cnt})::{RET}"))
      cat(glue::glue("{mssg}    new_cols={new_cols}.{RET}"))
      cat(glue::glue("{RET}"))
    }
  }  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    man_tib <-
      load_genome_studio_manifest(file = file,
                                  load_clean    = load_clean,
                                  load_controls = load_controls,
                                  write_clean   = write_clean,
                                  overwrite     = overwrite,
                                  ret_data      = ret_data,
                                  verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    man_cnt <- man_tib %>% base::nrow()
    if (verbose>=vt+1)
      cat(glue::glue("{mssg} Loaded Genome Studio manifest({man_cnt})!{RET2}"))
    
    if (is.null(man_tib) || !tibble::is_tibble(man_tib)) {
      stop(glue::glue("{RET}{mssg} ERROR: Failed to load Genome Studio ",
                      "manifest!{RET2}"))
      return(ret_tib)
    }
    man_key <- glue::glue("genome-studio-manifest-1({funcTag})")
    man_cnt <- print_tib(man_tib,funcTag, verbose, vt=vt+4,tc=tc+1, n=man_key)

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Swap Old/New Column Names::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    man_tib <- dplyr::rename_with(man_tib, ~ new_cols, dplyr::all_of(old_cols) )
    
    man_key <- glue::glue("genome-studio-manifest-renamed-2({funcTag})")
    man_cnt <- print_tib(man_tib,funcTag, verbose, vt=vt+4,tc=tc+1, n=man_key)
    
    if (!add_annotation) {
      man_tib <- man_tib %>% 
        dplyr::select(dplyr::any_of( c(new_cols,non_ann_cols) ))
    }
    
    man_key <- glue::glue("genome-studio-manifest-annotation-3({funcTag})")
    man_cnt <- print_tib(man_tib,funcTag, verbose, vt=vt+4,tc=tc+1, n=man_key)

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Extract/Format Probe A::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    prb_key <- "A"
    old_prbA_cols <- man_tib %>% 
      dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>% 
      names()
    new_prbA_cols <- old_prbA_cols %>%
      stringr::str_remove(paste0("_",prb_key,"$"))
    
    oldA_cnt <- old_prbA_cols %>% length()
    newA_cnt <- new_prbA_cols %>% length()
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} Old Sesame Columns:: Probe A({oldA_cnt}){RET}"))
    if (verbose>=vt+4)
      cat(glue::glue("{mssg}    old_prbA_cols={old_prbA_cols}.{RET2}"))
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} New Sesame Columns:: Probe A({newA_cnt}){RET}"))
    if (verbose>=vt+4)
      cat(glue::glue("{mssg}    new_prbA_cols={new_prbA_cols}.{RET2}"))
    
    prb_key <- "B"
    addA_tib <- NULL
    addA_tib <- man_tib %>% 
      dplyr::select(-dplyr::ends_with(paste0("_",prb_key))) %>% 
      dplyr::rename_with( ~ new_prbA_cols, dplyr::all_of(old_prbA_cols) ) %>%
      dplyr::mutate(Prb_Des=dplyr::case_when(
        Infinium_Design_Type=="I"  ~ "U",
        Infinium_Design_Type=="II" ~ "2",
        TRUE ~ NA_character_)
      )
    
    addA_cnt <- addA_tib %>% base::nrow()
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} Probe Set A: I/U Count={addA_cnt}.{RET2}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Extract/Format Probe B::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    prb_key <- "B"
    old_prbB_cols <- man_tib %>% 
      dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>% 
      names()
    new_prbB_cols <- old_prbB_cols %>%
      stringr::str_remove(paste0("_",prb_key,"$"))
    
    oldB_cnt <- old_prbB_cols %>% length()
    newB_cnt <- new_prbB_cols %>% length()
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} Old Sesame Columns:: Probe B({oldB_cnt}){RET}"))
    if (verbose>=vt+4)
      cat(glue::glue("{mssg}    old_prbB_cols={old_prbB_cols}.{RET2}"))
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} New Sesame Columns:: Probe B({newB_cnt}){RET}"))
    if (verbose>=vt+4)
      cat(glue::glue("{mssg}    new_prbB_cols={new_prbB_cols}.{RET2}"))
    
    prb_key <- "A"
    addB_tib <- NULL
    addB_tib <- man_tib %>% 
      dplyr::filter(Infinium_Design_Type=="I") %>%
      dplyr::select(-dplyr::ends_with(paste0("_",prb_key))) %>% 
      dplyr::rename_with( ~ new_prbB_cols, dplyr::all_of(old_prbB_cols) ) %>% 
      dplyr::mutate(Prb_Des="M")
    
    addB_cnt <- addB_tib %>% base::nrow()
    if (verbose>=vt+2)
      cat(glue::glue("{mssg} Probe Set B: M Count={addB_cnt}.{RET2}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Stack Probe Designs:: A/B
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- dplyr::bind_rows(addA_tib,addB_tib) %>%
      dplyr::mutate(
        Prb_Inf=dplyr::case_when(
          Infinium_Design_Type=="I"  ~ 1,
          Infinium_Design_Type=="II" ~ 2,
          TRUE ~ NA_real_
        ) %>% as.integer(),
        Prb_Din = Probe_ID %>%
          stringr::str_sub(1,2),
        Prb_Cgn = Prb_Cgn %>%
          stringr::str_remove("^ch.[0-9A-Z]+.") %>%
          stringr::str_remove("[A-Z]+$") %>%
          stringr::str_remove("^[^0-9]+") %>%
          stringr::str_remove("^0+") %>%
          as.integer(),
        Chromosome=Chromosome %>%
          stringr::str_remove("^chr") %>%
          paste0("chr",.)
      ) %>% dplyr::arrange(Probe_ID)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Summary and Sanity Checks::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # This better be zero::
    dup_cnt <- ret_tib %>% 
      dplyr::arrange(Address) %>% 
      dplyr::add_count(Address, name="Add_Cnt") %>% 
      dplyr::filter(Add_Cnt!=1) %>% base::nrow()
    if (dup_cnt != 0) {
      stop(glue::glue("{mssg} ERROR: Duplicates Found = {dup_cnt}!{RET2}"))
      return(NULL)
    }
    if (verbose>=vt+1)
      cat(glue::glue("{mssg} No duplicates detected = {dup_cnt}.{RET}"))

    if (verbose>=vt+2) {
      cat(glue::glue("{mssg} Probe Design Coverage={RET}"))
      sum_tib <- ret_tib %>% 
        dplyr::group_by(Prb_Des) %>% 
        dplyr::summarise(Count=n(), .groups = "drop")
      sum_tib %>% print(n=base::nrow(sum_tib))
      cat(glue::glue("{RET}"))
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

#
# epic_gst_dat <- 
#   load_genome_studio_manifest(file = opt$genome_manifest_csv,
#                               
#                               load_clean    = TRUE,
#                               load_controls = TRUE,
#                               write_clean   = TRUE,
#                               overwrite     = TRUE,
#                               ret_data      = TRUE,
#                               
#                               verbose = opt$verbose, tt = pTracker)
#

load_genome_studio_manifest = function(file,
                                       
                                       load_clean    = TRUE,
                                       load_controls = FALSE,
                                       write_clean   = TRUE,
                                       overwrite     = FALSE,
                                       ret_data      = FALSE,
                                       
                                       verbose=0, vt=3,tc=1,tt=NULL,
                                       funcTag='load_genome_studio_manifest') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}          file={file}.{RET}"))
    cat(glue::glue("{mssg}    load_clean={load_clean}.{RET}"))
    cat(glue::glue("{mssg}   write_clean={write_clean}.{RET}"))
    cat(glue::glue("{mssg}     overwrite={overwrite}.{RET}"))
    cat(glue::glue("{mssg}      ret_data={ret_data}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  probes_tib   <- NULL
  controls_tib <- NULL
  
  stime <- base::system.time({
    
    # Set Column Repair Prefix::
    col_prefix <- "Val_"
    if (stringr::str_detect(file,"EPIC")) col_prefix <- "HM"
    
    gs_control_names <- 
      c("Address","Control_Group","GS_Color","Control_Type")
    
    analytical_suffix_csv <- ".analytical.csv.gz"
    analytical_suffix_rds <- ".analytical.rds"
    control_suffix_csv    <- ".controls.csv.gz"
    
    # if (!R.utils::isGzipped(file, ".gz"))
    #   file <- R.utils::gzip(file) %>% as.character()
    if (!stringr::str_ends(file,".gz")) {
      system(glue::glue("gzip {file}"))
      file <- paste0(file,".gz")
    }
    
    if (verbose>=vt+3) {
      cat(glue::glue("{mssg} Setting suffix variables::{RET}"))
      cat(glue::glue("{mssg}       control_suffix_csv=",
                     "{control_suffix_csv}.{RET}"))
      cat(glue::glue("{mssg}    analytical_suffix_csv=",
                     "{analytical_suffix_csv}.{RET}"))
      cat(glue::glue("{mssg}    analytical_suffix_rds=",
                     "{analytical_suffix_rds}.{RET}"))
      cat(glue::glue("{RET}"))
    }
    
    # If the input file is the analytical file load it!
    if (stringr::str_ends(file, analytical_suffix_csv)) {
      clean_prb_csv <- file
      clean_col_rds <- clean_basename %>%
        paste(analytical_suffix_csv)
      clean_prb_done <- paste(clean_prb_csv,'done.txt', sep='.')
      
      if (verbose>=vt+2)
        cat(glue::glue("{mssg} Found clean manifest! Will load clean ",
                       "manifest={clean_prb_csv}...{RET}"))
      
      if (file.exists(clean_col_rds) &&
          file.mtime(clean_prb_csv) <= file.mtime(clean_col_rds)) {
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Using column header={clean_col_rds}...{RET}"))
        
        probes_type_cols <- readr::read_rds(clean_col_rds)
        probes_tib  <- readr::read_csv(clean_prb_csv,
                                       # col_names = names(probes_type_cols$cols),
                                       col_types = probes_type_cols)
      } else {
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Will build column header...{RET}"))
        
        probes_tib <- safe_read( clean_prb_csv, clean = TRUE, 
                                 verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
        
        # Get Column Types String::
        col_types_str <- sapply(probes_tib, typeof) %>% 
          tibble::as_tibble(rownames = "Key") %>% 
          dplyr::mutate(value=stringr::str_sub(value, 1,1)) %>% 
          dplyr::pull(value) %>% paste0(collapse = '')
        
        # Build cols() object for future loading::
        probes_type_cols <- spec(readr::read_csv(
          readr::format_csv(probes_tib), col_types = col_types_str) )
        
        if (write_clean) {
          if (verbose>=vt+2) 
            cat(glue::glue("{mssg} Writing clean manifest cols RDS...{RET}"))
          
          safe_write( probes_type_cols, file = clean_col_rds, done = TRUE,
                      verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
        }
      }
      
      if (load_controls) {
        controls_csv <- file %>% 
          stringr::str_replace(analytical_suffix_csv, control_suffix_csv)
        
        if (!file.exists(controls_csv)) {
          cat(glue::glue("{mssg} Warning: No controls file present. Try ",
                         "rebuilding genome studio manifest from scratch.{RET}"))
        } else {
          if (verbose>=vt+2)
            cat(glue::glue("{mssg} Loading controls={controls_csv}...{RET}"))
          
          controls_tib <- 
            suppressMessages(suppressWarnings( readr::read_csv( controls_csv )))
        }
      }
      
    } else {
      
      # Define clean formatted files::
      clean_basename <- file %>% 
        stringr::str_remove(".gz$") %>%
        stringr::str_remove(".csv$")
      
      clean_prb_csv <- clean_basename %>%
        paste0(analytical_suffix_csv)
      clean_col_rds <- clean_basename %>%
        paste0(analytical_suffix_rds)
      clean_prb_done <- paste(clean_prb_csv,'done.txt', sep='.')
      controls_csv   <- clean_basename %>%
        paste0(control_suffix_csv)
      
      if (verbose>=vt+3) {
        cat(glue::glue("{mssg} Setting output variables::{RET}"))
        cat(glue::glue("{mssg}    clean_basename={clean_basename}.{RET}"))
        cat(glue::glue("{mssg}      controls_csv={controls_csv}.{RET}"))
        cat(glue::glue("{mssg}     clean_col_rds={clean_col_rds}.{RET}"))
        cat(glue::glue("{mssg}     clean_prb_csv={clean_prb_csv}.{RET}"))
        cat(glue::glue("{mssg}    clean_prb_done={clean_prb_done}.{RET}"))
        cat(glue::glue("{RET}"))
      }
      
      if (load_clean && file.exists(clean_prb_csv) &&
          file.mtime(file) <= file.mtime(clean_prb_csv) &&
          file.mtime(clean_prb_csv) <= file.mtime(clean_prb_done) ) {
        
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Loading clean manifest={clean_prb_csv}...{RET}"))
        
        if (file.exists(clean_col_rds) &&
            file.mtime(clean_prb_csv) <= file.mtime(clean_col_rds)) {
          if (verbose>=vt+2)
            cat(glue::glue("{mssg} Using column header={clean_col_rds}...{RET}"))
          
          probes_type_cols <- readr::read_rds(clean_col_rds)
          
          probes_tib  <- readr::read_csv(clean_prb_csv,
                                         # col_names = names(probes_type_cols$cols),
                                         col_types = probes_type_cols)
        } else {
          if (verbose>=vt+2) 
            cat(glue::glue("{mssg} Will build column header...{RET}"))
          
          probes_tib <- safe_read( clean_prb_csv, clean = TRUE, 
                                   verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
          
          # Get Column Types String::
          col_types_str <- sapply(probes_tib, typeof) %>% 
            tibble::as_tibble(rownames = "Key") %>% 
            dplyr::mutate(value=stringr::str_sub(value, 1,1)) %>% 
            dplyr::pull(value) %>% paste0(collapse = '')
          
          # Build cols() object for future loading::
          probes_type_cols <- spec(readr::read_csv(
            readr::format_csv(probes_tib), col_types = col_types_str) )
          
          if (write_clean) {
            if (verbose>=vt+2) 
              cat(glue::glue("{mssg} Writing clean manifest cols RDS...{RET}"))
            
            safe_write( probes_type_cols, file = clean_col_rds, done = TRUE,
                        verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
          }
        }
        
        # Load Controls if requested
        if (load_controls) {
          if (!file.exists(controls_csv))
            cat(glue::glue("{mssg} Warning: No controls file present. Try ",
                           "rebuilding genome studio manifest from scratch.{RET}"))
          else
            controls_tib <- 
              suppressMessages(suppressWarnings( readr::read_csv( controls_csv )))
        }
        
      } else {
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Will load raw manifest={file}...{RET}"))
        
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Will build column header...{RET}"))
        
        lines_vec <- readr::read_lines(file)
        if (ret_data) ret_dat$lines <- lines_vec
        
        # Get Analytical and Control Start & End Indexes::
        beg_idx <- head(
          which( stringr::str_starts(head(lines_vec, n=20), "IlmnID" ) ), n=1)
        con_idx <- tail(
          which( stringr::str_starts(lines_vec, "\\[Controls\\]" ) ), n=1)
        end_idx <- lines_vec %>% length()
        
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Setting Indexes=({beg_idx}, {con_idx}, ",
                         "{end_idx}).{RET}"))
        
        # Get Analytical Column Names::
        probes_name_cols <- 
          lines_vec[beg_idx] %>% 
          stringr::str_split(pattern = ",", simplify = TRUE) %>% 
          as.vector() %>% 
          stringr::str_replace("^([0-9])", paste0(col_prefix, "\\$1") ) %>%
          stringr::str_remove_all("\\\\")
        
        probes_name_cnt <- probes_name_cols %>% length()
        if (verbose>=vt+3) {
          name1_str <- probes_name_cols[1]
          name2_str <- probes_name_cols[probes_name_cnt]
          out_str <- glue::glue("({name1_str} ... {name2_str})")
          cat(glue::glue("{mssg} Extracted Probe Names Vec({probes_name_cnt})",
                         " = {out_str}.{RET2}"))
        }
        
        # Build Analytical Probes Tibble
        probes_tib <- 
          lines_vec[c( (beg_idx+1):(con_idx-1) ) ] %>% # head() %>%
          tibble::as_tibble() %>% 
          tidyr::separate(value, into=c(probes_name_cols), sep=',') %>%
          clean_tibble()
        ret_key <- glue::glue("probes-tib({funcTag})")
        ret_cnt <- print_tib(probes_tib,funcTag, verbose,vt=vt+4,tc=tc+1, 
                             n=ret_key)
        
        # Get Column Types String::
        col_types_str <- sapply(probes_tib, typeof) %>% 
          tibble::as_tibble(rownames = "Key") %>% 
          dplyr::mutate(value=stringr::str_sub(value, 1,1)) %>% 
          dplyr::pull(value) %>% paste0(collapse = '')
        
        # Build cols() object for future loading::
        probes_type_cols <- spec(readr::read_csv(
          readr::format_csv(probes_tib), col_types = col_types_str) )
        if (verbose>=vt+4) {
          cat(glue::glue("{mssg} probes_type_cols={RET}"))
          probes_type_cols %>% print()
        }
        
        # Build Controls Tibble::
        controls_tib <- 
          lines_vec[c( (con_idx+1):(end_idx) ) ] %>%
          stringr::str_remove(",+$") %>%
          stringr::str_remove(",AVG$") %>%
          tibble::as_tibble() %>% 
          tidyr::separate(value, into=c(gs_control_names), sep=',') %>% 
          clean_tibble()
        ret_key <- glue::glue("controls-tib({funcTag})")
        ret_cnt <- print_tib(controls_tib,funcTag, verbose,vt=vt+4,tc=tc+1,
                             n=ret_key)
        
        if (ret_data) ret_dat$probes   <- probes_tib
        if (ret_data) ret_dat$controls <- controls_tib
        
        if (ret_data) ret_dat$probes_key_str <- probes_name_cols
        if (ret_data) ret_dat$probes_col_str <- col_types_str
        if (ret_data) ret_dat$probes_cols    <- probes_type_cols
        
        # Write clean files if requested::
        if (write_clean) {
          if (verbose>=vt+2) 
            cat(glue::glue("{mssg} Writing clean manifest files...{RET}"))
          
          out_tag <- glue::glue("write-controls-CSV({funcTag})")
          safe_write( controls_tib, file = controls_csv, done = TRUE, 
                      funcTag = out_tag, verbose=verbose, vt=vt+3,tc=tc+1,tt=tt)
          
          out_tag <- glue::glue("write-probes-CSV({funcTag})")
          safe_write( probes_tib, file = clean_prb_csv, done = TRUE,
                      funcTag = out_tag, verbose=verbose, vt=vt+3,tc=tc+1,tt=tt)
          
          out_tag <- glue::glue("write-probes-col-RDS({funcTag})")
          safe_write( probes_type_cols, file = clean_col_rds, done = TRUE,
                      funcTag = out_tag, verbose=verbose, vt=vt+3,tc=tc+1,tt=tt)
        }
      }
    }
    
    probes_tib <- probes_tib %>% dplyr::arrange(IlmnID)
    
    if (load_controls && is.null(controls_tib)) {
      ret_dat$probes   <- probes_tib
      ret_dat$controls <- controls_tib
      ret_data <- TRUE
    }
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(probes_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  if (ret_data) return(ret_dat)
  
  probes_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Noob Probe_ID Masking::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# changing default mod=100000000 to mod=999999999
# changing default prefix="nb" to prefix="cg9999"
# NOTE: 5 9's is optimal because improbe can only handle 15 characters
#  as input 2 character (cg) + 8 digits + 2 character ([TB][CO]) = 12
# However, we can reduce characters to a single leter 1+8 = 9 
#  Using 5 would put us at 14 leaving one charcter to expand 8 digits to 9
#
# None of that really matters since we'll never use these probes in improbe
#  The issue is running into your growing cg# (8 to 9 to 10) space. 
#  We'll make 5 9's
# Wait why don't we just make cgX
noob_mask = function(x, seed=21L, mod=100000000, prefix="cgBK",
                     funcTag='noob_mask') {
  
  # Clean input::
  x <- x %>% stringr:: str_remove_all("[^0-9]+") %>% as.integer()
  
  hash <- digest::digest2int(as.character(x), seed) %% mod
  
  # Waiting to see if this ever fails
  m_len <- stringr::str_length(format(mod, scientific = FALSE))
  h_len <- stringr::str_length(hash)
  
  if (h_len>=m_len) {
    stop(glue::glue("{RET}[{funcTag}]: ERROR: ",
                    "h_len({h_len}) >= m_len({m_len})!!!{RET}",
                    "hash={hash}{RET}",
                    "mod={mod}{RET}{RET}"))
    return(NULL)
  }
  
  # Now remake cg style number::
  hash <- paste0(prefix,
                 stringr::str_pad(hash, width=m_len-1,side="left", pad="0"))
  
  hash
}

noob_mask_manifest = function(tib,
                              key = "Probe_ID", 
                              out = NULL, 
                              prefix = "cg",
                              
                              verbose=0, vt=3,tc=1,tt=NULL,
                              funcTag='noob_mask_manifest') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}      key={key}.{RET}"))
    cat(glue::glue("{mssg}      out={out}.{RET}"))
    cat(glue::glue("{mssg}   prefix={prefix}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (is.null(out)) out <- key
    key_sym <- rlang::sym(key)
    out_sym <- rlang::sym(out)
    
    noob_vec <- tib %>% 
      dplyr::pull(key_sym) %>% 
      lapply( noob_mask) %>% unlist()
    
    ret_tib <- tib %>% dplyr::mutate(!!out_sym := noob_vec)
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# End of file
