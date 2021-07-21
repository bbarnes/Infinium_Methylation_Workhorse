
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   improbe (Infinium Methylation) Methods::
#
#                           [FUNCTION GRAVEYARD]
#
#  The functions have been replaced. This file should be deleted, but is 
#           just a resting place for these functions for now :)
#
#              These functions should probably be deprecated
#
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
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
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
#                        Modern improbe IO Methods:: 
#                          Not so Modern Any More
#
#              These functions should probably be deprecated
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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
#                        Modern improbe IO Methods:: 
#                          Not so Modern Any More
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

addSeq48U = function(tib, field,
                     sidx=2, plen=50,
                     verbose=0,vt=3,tc=1,tt=NULL,
                     funcTag='addSeq48U') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    idx1 <- sidx
    len1 <- plen - 1
    idx2 <- sidx + 1
    len2 <- plen
    
    field_sym <- rlang::sym(field)
    
    ret_tib <- tib %>% dplyr::mutate(
      Seq_48U_1=stringr::str_sub(!!field_sym, idx1,len1) %>% 
        stringr::str_to_upper() %>% 
        stringr::str_replace_all('R', 'A') %>% 
        stringr::str_replace_all('Y', 'T'),
      
      Seq_48U_2=stringr::str_sub(!!field_sym, idx2,len2) %>% 
        stringr::str_to_upper() %>% 
        stringr::str_replace_all('R', 'A') %>% 
        stringr::str_replace_all('Y', 'T'),
    )
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

isTibIMP = function(tib,
                    verbose=0,vt=3,tc=1,tt=NULL,
                    funcTag='isTibIMP') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- FALSE
  stime <- system.time({
    
    imp_col <- INIT_IMP_HEADER()
    imp_key <- imp_col$cols %>% names()
    sam_col <- tib %>% names()
    
    if (base::all.equal(imp_key,sam_col)) ret_val <- TRUE
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_val
}

loadIMP = function(file, max=Inf, format=NULL, del='\t',
                   verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadIMP'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    imp_col <- INIT_IMP_HEADER()
    imp_key <- imp_col$cols %>% names()
    sam_col <- utils::read.table(file, header=TRUE, sep=del, nrows=1) %>% names()
    
    if (max==Inf) {
      if (base::all.equal(imp_key,sam_col)) {
        ret_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file, col_types = imp_col) ))
      } else {
        ret_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
      }
    } else {
      if (!is.numeric(max)) {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: max must be an integer; max={max}!!!{RET}{RET}"))
        return(ret_tib)
      }
      max <- as.integer(max)
      
      if (base::all.equal(imp_key,sam_col)) {
        ret_tib <- suppressMessages(suppressWarnings( 
          utils::read.table(file, header=TRUE, sep=del, nrows=max, colClasses=imp_col) )) %>% 
          tibble::as_tibble() %>% dplyr::mutate(Chromosome=as.character(Chromosome))
      } else {
        ret_tib <- suppressMessages(suppressWarnings( 
          utils::read.table(file, header=TRUE, sep=del, nrows=max) )) %>% tibble::as_tibble()
      }
    }
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

tibToFas = function(tib, key,seq, prefix,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'tibToFas'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; key={key}, seq={seq}, prefix={prefix}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    dir <- base::dirname(prefix)
    if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    
    file_fas <- paste0(prefix, ".fa.gz")
    file_csv <- paste0(prefix, ".csv.gz")
    
    key_sym <- rlang::sym(key)
    seq_sym <- rlang::sym(seq)
    
    ret_tib <- tib %>% dplyr::mutate(
      # fas_line=paste0(">",!!key, "\n",!!seq)
      fas_line=paste0(">",!!key_sym, "\n",!!seq_sym)
    )  %>%
      dplyr::filter(!is.na(fas_line)) %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    fas_vec <- ret_tib %>%
      # dplyr::distinct(fas_line) %>%
      dplyr::pull(fas_line)
    
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Writing Data CSV={file_csv}...{RET}"))
    readr::write_csv(ret_tib, file_csv)
    
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Writing Data FAS={file_fas}...{RET}"))
    readr::write_lines(x=fas_vec, file=file_fas)
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      ret_tib %>% print()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

ordToFas = function(tib, dir, name, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'ordToFas'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}, dir={dir}.{RET}"))
  
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  fas_file <- file.path(dir, paste0(name,'.cgn.fa.gz'))
  
  fas_tib <- tib %>% dplyr::mutate(
    line=dplyr::case_when(
      Normalization_Bin=='C' ~ paste0('>',AlleleA_Probe_Id,'\n',AlleleA_Probe_Sequence),
      Normalization_Bin!='C' ~ paste0('>',AlleleA_Probe_Id,'\n',AlleleA_Probe_Sequence,'\n',
                                      '>',AlleleB_Probe_Id,'\n',AlleleB_Probe_Sequence),
      TRUE ~ NA_character_ ) ) %>%
    dplyr::filter(!is.na(line)) %>% dplyr::pull(line)
  
  readr::write_lines(x=fas_tib, file=fas_file)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  fas_file
}

writeBedFas = function(tib, dir, name, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'writeBedFas'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}, dir={dir}.{RET}"))
  
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  fas_file <- file.path(dir, paste0(name,'.fa.gz'))
  bed_file <- file.path(dir, paste0(name,'.sorted.bed.gz'))
  
  tib <- tib %>% dplyr::arrange(-MinScore) %>% dplyr::distinct(PRB1_U, .keep_all=TRUE) %>%
    dplyr::arrange(Chromosome,Coordinate) %>% dplyr::rename(Chrom=Chromosome, Beg=Coordinate) %>%
    dplyr::mutate(End=Beg+stringr::str_length(PRB1_U), 
                  Full_ID=paste(Seq_ID,FR_Imp,TB_Imp,CO_Imp,Beg, sep='_') )
  
  fas_tib <- tib %>% dplyr::select(Full_ID, PRB1_U) %>% 
    dplyr::mutate(line=paste0('>',Full_ID,'\n',PRB1_U)) %>% dplyr::pull(line)
  bed_tib <- tib %>% dplyr::select(Chrom, Beg, End, Full_ID, MinScore, FR_Imp, TB_Imp, CO_Imp, everything())
  
  readr::write_lines(x=fas_tib, file=fas_file)
  readr::write_tsv(bed_tib, bed_file)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  bed_tib
}

writeImprobeInput = function(tib, name, dir, run=FALSE, 
                             exe=NULL, impTango=NULL, imp13mer=NULL, 
                             tbVar='BOTH', coVar='BOTH',
                             verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'writeImprobeInput'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}, dir={dir}.{RET}"))
  
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  
  imp_run_sh  <- file.path(dir, paste0(name,'.improbe-input.sh'))
  imp_inp_tsv <- file.path(dir, paste0(name,'.improbe-input.tsv'))
  imp_out_tsv <- file.path(dir, paste0(name,'.improbe-output.tsv'))
  imp_log_txt <- file.path(dir, paste0(name,'.improbe-output.log'))
  imp_fin_txt <- file.path(dir, paste0(name,'.improbe-output.fin.txt'))
  imp_cmd_str <- ''
  
  readr::write_tsv(tib, imp_inp_tsv)
  
  if (!is.null(exe) && ! is.null(impTango) && !is.null(imp13mer)) {
    # stopifnot(file.exists(exe))
    # stopifnot(file.exists(impTango))
    # stopifnot(file.exists(imp13mer))
    
    # ${CMD} -oASPE -tBOTH -cBoth -n${fn_13mer} -a${fn_tango} -V ${fn_in} >$fn_out 2>${fn_log}
    imp_cmd_str <- glue::glue("{exe} -oASPE -t{tbVar} -c{coVar} -n{imp13mer} -a{impTango} -V ",
                              "{imp_inp_tsv} >{imp_out_tsv} 2>{imp_log_txt}{RET}",
                              "touch {imp_fin_txt}{RET}")
    
    readr::write_file(imp_cmd_str, imp_run_sh)
    Sys.chmod(imp_run_sh, mode='0777')
    
    if (run) system(imp_run_sh, wait = TRUE)
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  imp_out_tsv
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Load/Format Input File Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
loadManifestRS = function(file, swap, revAllele=FALSE, 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(RS)={file}.{RET}"))
  
  snp_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(file) ))
  snp_swp_tib <- suppressMessages(suppressWarnings(readr::read_tsv(swap) ))
  
  # NOTE: Rediculous Swapping Needed for Infinium I M/U -> U/M from old manufacturing mistake::
  #  - This method should not be used...
  if (revAllele) {
    snp_man_tib <- snp_man_tib %>% 
      dplyr::mutate(
        TMP_ADD_A=AddressA_ID,
        TMP_ADD_B=AddressB_ID,
        TMP_SEQ_A=AlleleA_ProbeSeq,
        TMP_SEQ_B=AlleleB_ProbeSeq,
        IS_SWAP=case_when(
          Infinium_Design_Type=='I' & stringr::str_ends(AlleleA_ProbeSeq,'C') ~ TRUE,
          TRUE ~ FALSE
        ),
        AddressA_ID=case_when( IS_SWAP ~ TMP_ADD_B, TRUE ~ AddressA_ID),
        AddressB_ID=case_when( IS_SWAP ~ TMP_ADD_A, TRUE ~ AddressB_ID),
        AlleleA_ProbeSeq=case_when( IS_SWAP ~ TMP_SEQ_B, TRUE ~ AlleleA_ProbeSeq ),
        AlleleB_ProbeSeq=case_when( IS_SWAP ~ TMP_SEQ_A, TRUE ~ AlleleB_ProbeSeq )
      )
  }
  
  snp_man_tib <-  snp_man_tib %>% 
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$') ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq, 
                  Infinium_Design_Type, Next_Base, Color_Channel, CGN)
  
  snp_swp_tib <- snp_swp_tib %>% tidyr::separate(Seq_ID, into=c('IlmnID', 'diNUC'), sep='_') %>% 
    dplyr::rename(Forward_Sequence=Sequence, CHR=Chromosome, MAPINFO=Coordinate) %>% 
    dplyr::mutate(Forward_Sequence=stringr::str_replace(Forward_Sequence,'\\[CG\\]', paste0('[',diNUC,']') ) ) %>% 
    dplyr::select(-CpG_Island)
  
  snp_man_tib <- snp_man_tib %>% dplyr::inner_join(snp_swp_tib, by="IlmnID") %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  snp_man_tib
}                          

loadManifestCG = function(file, pr='cg', 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestCG'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CG)={file}.{RET}"))
  
  cpg_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(file, skip=7) )) %>%
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[([A-Za-z][A-Za-z])].*$', "\\$1"),'^\\\\') ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cpg_man_tib
}                          

loadManifestCH = function(file, pr='ch', ry='R', 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestCH'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CH)={file}.{RET}"))
  
  cph_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(cph_man_csv ))) %>% 
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(Strand=stringr::str_sub(IlmnID, -1), FR=Strand, CO='O',
                  TP=case_when( is.na(AlleleB_ProbeSeq) ~ 'II', TRUE ~ 'I' )) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=paste0(ry,stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[[A-Za-z]([A-Za-z])].*$', "\\$1"),'^\\\\')),
                  Forward_Sequence=stringr::str_replace(Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', paste0('[',diNUC,']')) ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cph_man_tib
}

loadManifestYH = function(file, pr='ch', 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestYH'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CH)={file}.{RET}"))
  
  cph_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(cph_man_csv ))) %>% 
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(Strand=stringr::str_sub(IlmnID, -1), FR=Strand, CO='O',
                  TP=case_when( is.na(AlleleB_ProbeSeq) ~ 'II', TRUE ~ 'I' )) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=paste0('Y',stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[[A-Za-z]([A-Za-z])].*$', "\\$1"),'^\\\\')),
                  Forward_Sequence=stringr::str_replace(Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', paste0('[',diNUC,']')) ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cph_man_tib
}                          

#
# TBD:: This function should probably be deprecated::
#
loadImprobeDesign = function(file=NULL, src_des_tib=NULL, 
                             verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadImprobeDesign'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading improbe={file}.{RET}"))
  
  if (!is.null(file)) src_des_tib <- suppressMessages(suppressWarnings(readr::read_tsv(file)))
  stopifnot(!is.null(src_des_tib))
  
  src_des_tib <- src_des_tib %>%
    dplyr::rename(PRB1_U=UnMethyl_Probe_Sequence,
                  PRB1_M=Methyl_Probe_Sequence,
                  NXB_U=UnMethyl_Next_Base,
                  NXB_M=Methyl_Next_Base) %>%
    dplyr::mutate(Probe_ID=paste(Seq_ID, Methyl_Allele_FR_Strand, stringr::str_sub(Methyl_Allele_TB_Strand,1,1), 
                                 Methyl_Allele_CO_Strand, sep='_'),
                  FR=case_when(Methyl_Allele_FR_Strand=='F'   ~ TRUE, Methyl_Allele_FR_Strand=='R'   ~ FALSE, TRUE ~ NA),
                  TB=case_when(Methyl_Allele_TB_Strand=='TOP' ~ TRUE, Methyl_Allele_TB_Strand=='BOT' ~ FALSE, TRUE ~ NA),
                  CO=case_when(Methyl_Allele_CO_Strand=='C'   ~ TRUE, Methyl_Allele_CO_Strand=='O'   ~ FALSE, TRUE ~ NA),
                  diNUC='CG',
                  NXB_IMP=case_when(NXB_U==NXB_M ~ NXB_U,
                                    TRUE ~ NA_character_),
                  COL_U=case_when(NXB_U=='A'|NXB_U=='T' ~ FALSE, # Red == FALSE
                                  NXB_U=='C'|NXB_U=='G' ~ TRUE,  # Grn == TRUE
                                  TRUE ~ NA),
                  COL_M=case_when(NXB_M=='A'|NXB_M=='T' ~ FALSE, # Red == FALSE
                                  NXB_M=='C'|NXB_M=='G' ~ TRUE,  # Grn == TRUE
                                  TRUE ~ NA),
                  COL_IMP=case_when(COL_U==COL_M ~ COL_U,
                                    TRUE ~ NA)
                  
    ) %>%
    
    # Design Score Parameters::
    dplyr::rename(
      PRB_SCR_U=UnMethyl_Final_Score,
      PRB_SCR_M=Methyl_Final_Score,
      PRB_SCR_S=Probeset_Score,
      
      TM_RAW_M=Methyl_Tm,
      TM_SCR_M=Methyl_Tm_Score,
      GC_RAW_M=Methyl_GC_Percent,
      GC_SCR_M=Methyl_GC_Score,
      KM_RAW_M=Methyl_13mer_Count,
      KM_SCR_M=Methyl_13mer_Score,
      AD_RAW_M=Methyl_Address_Count,
      AD_SCR_M=Methyl_Address_Score,
      CM_RAW_M=Methyl_Self_Complementarity,
      CM_SCR_M=Methyl_Self_Complementarity_Score,
      MO_RAW_M=Methyl_Mono_Run,
      MO_SCR_M=Methyl_Mono_Run_Score,
      EP_RAW_M=Methyl_Ectopic_Count,
      EP_SCR_M=Methyl_Ectopic_Score,
      CG_RAW_M=Methyl_Underlying_CpG_Count,
      MD_RAW_M=Methyl_Underlying_CpG_Min_Dist,
      CG_SCR_M=Methyl_Underlying_CpG_Score,
      NB_SCR_M=Methyl_Next_Base_Score,
      
      TM_RAW_U=UnMethyl_Tm,
      TM_SCR_U=UnMethyl_Tm_Score,
      GC_RAW_U=UnMethyl_GC_Percent,
      GC_SCR_U=UnMethyl_GC_Score,
      KM_RAW_U=UnMethyl_13mer_Count,
      KM_SCR_U=UnMethyl_13mer_Score,
      AD_RAW_U=UnMethyl_Address_Count,
      AD_SCR_U=UnMethyl_Address_Score,
      CM_RAW_U=UnMethyl_Self_Complementarity,
      CM_SCR_U=UnMethyl_Self_Complementarity_Score,
      MO_RAW_U=UnMethyl_Mono_Run,
      MO_SCR_U=UnMethyl_Mono_Run_Score,
      EP_RAW_U=UnMethyl_Ectopic_Count,
      EP_SCR_U=UnMethyl_Ectopic_Score,
      CG_RAW_U=UnMethyl_Underlying_CpG_Count,
      MD_RAW_U=UnMethyl_Underlying_CpG_Min_Dist,
      CG_SCR_U=UnMethyl_Underlying_CpG_Score,
      NB_SCR_U=UnMethyl_Next_Base_Score
    ) %>%
    dplyr::select(Probe_ID, Seq_ID, FR, TB, CO, diNUC, NXB_IMP, COL_IMP, PRB1_U, PRB1_M, 
                  Forward_Sequence, 
                  Genome_Build, Chromosome, Coordinate,
                  dplyr::contains("_RAW_"), dplyr::contains("_SCR_"))
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  src_des_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Format Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

addDesignSeqCG = function(tib, seq, add, din=NULL,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'addDesignSeqCG'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    add_sym <- add %>% rlang::sym()
    seq_sym <- seq %>% rlang::sym()
    
    ret_tib <- tib %>% 
      dplyr::mutate(!!add:=stringr::str_replace(!!seq_sym, "\\[[a-zA-Z][a-zA-Z]\\]", "[CG]"))
    
    if (!is.null(din))
      ret_tib <- ret_tib %>% dplyr::mutate(
        !!din := !!seq_sym %>%
          stringr::str_remove("\\].*$") %>%
          stringr::str_remove("^.*\\[") )
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

replaceDesignSeqCG = function(tib, seq, add, nuc,
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'replaceDesignSeqCG'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    add_sym <- add %>% rlang::sym()
    seq_sym <- seq %>% rlang::sym()
    nuc_sym <- nuc %>% rlang::sym()
    
    ret_tib <- tib %>% 
      dplyr::mutate(PreSeq=stringr::str_remove(!!seq_sym, "\\[.*$"),
                    PosSeq=stringr::str_remove(!!seq_sym, "^.*\\]"),
                    !!add:=paste0(PreSeq,"[",!!nuc_sym,"]",PosSeq) ) %>%
      dplyr::select(-PreSeq, -PosSeq)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          BSMAP Alignment Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsmapProbeAlign = function(exe, fas, gen, dir, 
                           opt="-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R", 
                           lan=NULL, run=FALSE,
                           verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'bsmapProbeAlign'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} dir={dir}.{RET}"))
  
  stime <- system.time({
    if (is.null(opt) || length(opt)==0)
      opt="-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
    
    fas_name <- fas %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
    gen_name <- gen %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
    
    # Create Shell and Alignment directories::
    sh_dir <- file.path(dir, 'shells')
    al_dir <- file.path(dir, 'align')
    
    if (!dir.exists(sh_dir)) dir.create(sh_dir, recursive=TRUE)
    if (!dir.exists(al_dir)) dir.create(al_dir, recursive=TRUE)
    
    aln_ssh <- file.path(sh_dir, paste0('run_bsp-',fas_name,'-',gen_name,'.sh') )
    aln_bsp <- file.path(al_dir, paste0(fas_name,'-',gen_name,'.bsmap.bsp') )
    ret_tsv <- file.path(al_dir, paste0(fas_name,'-',gen_name,'.bsmap.formatted.tsv.gz') )
    aln_cmd <- paste(exe, '-a',fas, '-d',gen, opt,'-o',aln_bsp, sep=' ')
    aln_cmd <- glue::glue("{aln_cmd}{RET}",
                          "gzip {aln_bsp}{RET}",
                          "gzip -dc {aln_bsp}.gz | cut -f 1,2,4-11 | perl -pe 's/_/\t/gi; s/:/\t/gi' | ",
                          "gzip -c - > {ret_tsv}{RET}")
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{TAB} Launching; CMD={aln_cmd}...{RET}"))
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{TAB} Launching bsmap alignments: {fas_name} vs. {gen_name}...{RET}"))
    readr::write_lines(x=aln_cmd, file=aln_ssh, append=FALSE)
    Sys.chmod(paths=aln_ssh, mode="0777")
    
    if (!is.null(lan)) aln_ssh <- paste(lan,aln_ssh, sep=' ')
    if (run) base::system(aln_ssh)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  if (run) return(ret_tsv)
  aln_ssh
}

# End of file
