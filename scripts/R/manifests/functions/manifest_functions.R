
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

#
# TBD: 
#  1. Move functions from sample_sheet_functions.R to manifest_functions.R
#  2. genomeStudioToSesameManifest()
#     - Sesame base CSV file
#     - Sesame full CSV file with full names, bscU seqs
#     - Fasta file with full names, addresses, bscU seqs
#  3. bspToFullManifest()
#     - Write BED file
#     - Write Full manifest
#     - Rename multi-unique probes
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Specialized Mouse Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# fixOrderProbeIDs(full_man_tib, verbose=opt$verbose)
fixOrderProbeIDs = function(tib, field="Probe_Type",
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'manToBeadSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  # Add Infinium Type Definition:: Also, fail any probe that breaks this method
  tib <- tib %>% dplyr::mutate(
    M=dplyr::case_when(M=='NA' ~ NA_real_, TRUE ~ M),
    Infinium_Design=dplyr::case_when(
      !is.na(AlleleA_Probe_Sequence) &  is.na(AlleleB_Probe_Sequence) & !is.na(U) &  is.na(M) ~ 1,
      !is.na(AlleleA_Probe_Sequence) & !is.na(AlleleB_Probe_Sequence) & !is.na(U) & !is.na(M) ~ 2,
      TRUE ~ NA_real_
    )
  )
  
  # Add Probe Rep_Num, Old_Probe_ID, and Probe Lengths
  tib <- tib %>%
    dplyr::filter(!is.na(Infinium_Design)) %>%
    dplyr::add_count(Probe_ID, name="Rep_Num") %>%
    dplyr::mutate(AlleleA_Probe_Length=stringr::str_length(AlleleA_Probe_Sequence),
                  AlleleB_Probe_Length=stringr::str_length(AlleleB_Probe_Sequence),
                  Old_Probe_ID=paste(Probe_ID,paste0('r',Rep_Num), sep='_'))

  # Probe_Type Group_Count
  # <chr>            <int>
  # 2 cg              278216  cg36638866_F_T_C_II
  # 3 ch                3763  ch110207106_F_B_C_II
  # 7 rp                4518  rp10563_R_B_C_II
  #
  # 4 mu                6307  mu31334197_F_T_C_s5_h2_II
  #
  # 8 rs                1669  rs249350499_YC_F_R_B_C_I
  #
  # Skipping for now::
  # 1 BS                1041  BSC_CCND3-CT-R-40170_F_B_O_II
  # 6 NO                 815  NON_CCND3.87768_GT_2_R_B_O_II
  # 5 ne                 395  neg_cg42568514_R_B_C
  tib_types <- tib %>% split(.$Probe_Type)
  
  ret_tib <- NULL
  for (type in names(tib_types)) {
    cur_tib <- tib_types[[type]]
    if (type=='cg' || type=='ch' || type=='rp') {
      cur_tib <- cur_tib %>% 
        dplyr::mutate(
          Probe_ID=stringr::str_replace(Probe_ID, '_([FR])_([CO]_[I]+)$',  '_\\R1_N_\\$2')
          # Probe_ID=stringr::str_replace(Probe_ID, '_F_C_II$', '_F_N_C_II'),
          # Probe_ID=stringr::str_replace(Probe_ID, '_R_C_I$',  '_R_N_C_I'),
          # Probe_ID=stringr::str_replace(Probe_ID, '_R_C_II$', '_R_N_C_II')
        ) %>%
        tidyr::separate(Probe_ID, into=c('Seq_ID','FR','TB','CO','PD'), sep='_', remove=FALSE)

    } else if (type=='mu') {
      cur_tib <- cur_tib %>% dplyr::mutate(Probe_ID=stringr::str_remove(Probe_ID, '_[0-9]+M$')) %>%
        tidyr::separate(Probe_ID, into=c('Seq_ID','FR','TB','CO','DS','HS', 'PD'), sep='_', remove=FALSE) %>%
        dplyr::mutate(DS=stringr::str_remove(DS, 's'), HS=stringr::str_remove(HS, 'h'),
                      PD=dplyr::case_when(is.na(PD) ~ DS, TRUE ~ PD),
                      DS=dplyr::case_when(DS=='I' | DS=='II' ~ NA_character_, TRUE ~ DS)
        )
    } else if (type=='rs') {
      cur_tib <- cur_tib %>% dplyr::filter(Probe_Type=='rs') %>% dplyr::mutate(
        Probe_ID=stringr::str_replace(Probe_ID, '^rs([A-Z0-9]+)_([0-9]+)_', 'rs-\\$1-\\$2_') %>% stringr::str_remove_all('\\\\'),
        Probe_ID=stringr::str_replace(Probe_ID, '_([FR])_([CO]_[I]+)$',  '_NN_N_\\$1_N_\\$2') %>% stringr::str_remove_all('\\\\')
      ) %>% 
        tidyr::separate(Probe_ID, into=c('Seq_ID','Di','FN','FR','TB','CO','PD'), sep='_', remove=FALSE)
    } else if (type=='BS') {
      cur_tib <- cur_tib %>% dplyr::mutate(
        Di=stringr::str_replace(Probe_ID, '^.*-([ACTG][ACTG])-.*$', '\\$1') %>% stringr::str_remove_all('\\\\')
      ) %>% 
        dplyr::mutate(Probe_ID=stringr::str_replace(Probe_ID, '^BSC_', 'bs-')) %>% 
        tidyr::separate(Probe_ID, into=c('Seq_ID','FR','TB','CO','PD'), sep='_', remove=FALSE) 
    } else if (type=='NO') {
      cur_tib <- cur_tib %>% dplyr::mutate(
        Di=stringr::str_replace(Probe_ID, '^.*-([ACTG][ACTG])-.*$', '\\$1') %>% stringr::str_remove_all('\\\\')
      ) %>% 
        dplyr::mutate(Probe_ID=stringr::str_replace(Probe_ID, '^NON_', 'no-') %>% 
                        stringr::str_replace('\\.', '-') %>% stringr::str_replace('_','-') %>% stringr::str_replace('_','-')) %>% 
        tidyr::separate(Probe_ID, into=c('Seq_ID','FR','TB','CO','PD'), sep='_', remove=FALSE) 
    } else if (type=='ne') {
      cur_tib <- cur_tib %>% 
        dplyr::mutate(Probe_ID=paste(stringr::str_replace(Probe_ID, '^neg_', 'ne-'),'_I') ) %>% 
        tidyr::separate(Probe_ID, into=c('Seq_ID','FR','TB','CO','PD'), sep='_', remove=FALSE) 
    } else {
      cat(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: Unsupported type={type}!{RET}{RET}"))
    }    
    ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
    
  }
  
  # Add Probe Sequences for matching::
  ret_tib <- ret_tib %>% dplyr::mutate(
    Mat_Prb=dplyr::case_when(
      Infinium_Design=='I'  ~ stringr::str_sub(AlleleA_Probe_Sequence, 2,50),
      Infinium_Design=='II' ~ stringr::str_sub(AlleleA_Probe_Sequence, 3,50),
      TRUE ~ NA_character_
    ) %>% stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% stringr::str_replace_all('Y', 'T'),
    
    Mat_PrbA=dplyr::case_when(
      Infinium_Design=='I'  ~ stringr::str_sub(AlleleA_Probe_Sequence, 2,49),
      Infinium_Design=='II' ~ stringr::str_sub(AlleleA_Probe_Sequence, 3,50),
      TRUE ~ NA_character_
    ) %>% stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% stringr::str_replace_all('Y', 'T')
  ) %>% dplyr::select(Seq_ID, FR,TB,CO,PD,Infinium_Design,Mat_PrbA,Mat_Prb, everything())
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Manifest Stats Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

manToBeadSummary = function(man, field="Infinium_Design",
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'manToBeadSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} field={field}.{RET}"))
  
  field <- rlang::sym(field)
  inf1_cnt <- man %>% dplyr::group_by(!!field) %>% dplyr::summarise(Count=n()) %>% 
    dplyr::filter(!!field=='I') %>% dplyr::summarise(Sum_Count=sum(Count)) %>% dplyr::pull(Sum_Count) %>% as.integer()
  inf2_cnt <- man %>% dplyr::group_by(!!field) %>% dplyr::summarise(Count=n()) %>% 
    dplyr::filter(!!field=='II') %>% dplyr::summarise(Sum_Count=sum(Count)) %>% dplyr::pull(Sum_Count) %>% as.integer()
  
  bead1_cnt <- inf1_cnt * 2
  bead2_cnt <- inf2_cnt
  beads_cnt <- bead1_cnt + bead2_cnt
  beads_per <- bead1_cnt / beads_cnt
  
  sum_tib <- tibble::tibble(Inf1=inf1_cnt, Inf2=inf2_cnt, 
                            Bead1=bead1_cnt, Bead2=bead2_cnt, Beads=beads_cnt,
                            BeadPerc=beads_per)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done...{RET}"))
  
  sum_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        AQP Manifest Generation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

decodeToManifest = function(ord, mat, aqp=NULL, pqc=NULL, platform, version,
                            ordSkip=8, matSkip=40, aqpSkip=7,
                            guessCnt=50000,
                            full=FALSE, original=FALSE, cleanAdds=TRUE,
                            verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'decodeToManifest'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  if (is.null(aqp) && is.null(pqc))
    stop(glue::glue("{RET}[{funcTag}]: Both aqp and pqc are NULL!!!{RET}{RET}"))

  manifest_str <- paste(platform,version, sep='-')
  match_cnames <- c("Plate", "Row", "Col", "Address", "Mod5", "Full_Seq", "Mod3", "Comment")
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading Order(skip={ordSkip})={ord}.{RET}"))
  ord_tib <- suppressMessages(suppressWarnings( readr::read_csv(ord, skip=ordSkip, guess_max=guessCnt)  )) %>%
    dplyr::mutate(AlleleB_Probe_Id=as.character(AlleleB_Probe_Id), 
                  AlleleB_Probe_Sequence=as.character(AlleleB_Probe_Sequence))
  if (verbose>=vt) print(ord_tib)
  
  if (original) {
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading Match(org)={mat}.{RET}"))
    
    mat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(mat) )) %>%
      dplyr::rename(Probe_ID=probe_id, Full_Seq=bo_seq, Address=address_name, Ord_Seq=sequence) %>% 
      dplyr::mutate(Tango_Seq=stringr::str_sub(Full_Seq,1,45),
                    Probe_Seq=stringr::str_sub(Full_Seq,46))
    
  } else {
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading Match(new; skip={matSkip})={mat}.{RET}"))
    
    mat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(mat, skip=matSkip) )) %>% 
      purrr::set_names(match_cnames) %>% 
      dplyr::mutate(Tango_Seq=stringr::str_sub(Full_Seq,1,45), 
                    Probe_Seq=stringr::str_sub(Full_Seq,46))
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} mat_tib::{RET}"))
  if (verbose>=vt) print(mat_tib)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading AQP(skip={aqpSkip})={aqp}.{RET}"))
  if (!is.null(aqp)) {
    aqp_tib <- suppressMessages(suppressWarnings( readr::read_tsv(aqp, skip=aqpSkip) )) %>%
      dplyr::select(Address, Decode_Status)
  } else if (!is.null(pqc)) {
    aqp_tib  <- suppressMessages(suppressWarnings( readr::read_tsv(pqc,skip=aqpSkip) )) %>% 
      purrr::set_names(stringr::str_replace_all(names(.),' ','_')) %>%
      dplyr::rename(Decode_Status=Status) %>%
      dplyr::select(Address, Decode_Status)
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} aqp_tib::{RET}"))
  if (verbose>=vt) print(aqp_tib)
  
  if (full) {
    ret <- NULL
    ret$aqp <- aqp_tib
    ret$mat <- mat_tib
    ret$ord <- ord_tib
    ret$man <- man_tib
    
    return(ret)
  }
  
  mat_tib <- mat_tib %>% dplyr::inner_join(aqp_tib, by="Address")
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} join(mat_tib, aqp_tib)::{RET}"))
  if (verbose>=vt) print(mat_tib)
  
  man_tib <- ord_tib %>% dplyr::rename(Probe_ID=Assay_Design_Id) %>% dplyr::select(-AlleleA_Probe_Id, -AlleleB_Probe_Id) %>%
    dplyr::left_join(dplyr::select(mat_tib, Address, Tango_Seq, Probe_Seq, Decode_Status), by=c("AlleleA_Probe_Sequence"="Probe_Seq")) %>%
    dplyr::rename(Address_A=Address, Address_A_Seq=Tango_Seq, QC_A_Action=Decode_Status) %>%
    dplyr::left_join(dplyr::select(mat_tib, Address, Tango_Seq, Probe_Seq, Decode_Status), by=c("AlleleB_Probe_Sequence"="Probe_Seq")) %>%
    dplyr::rename(Address_B=Address, Address_B_Seq=Tango_Seq, QC_B_Action=Decode_Status)
  
  # TBD:: Should redo the logic to remove -999 too
  #  - probably easier to remove !=0 or NA...
  man_tib <- man_tib %>% dplyr::filter(QC_A_Action!=-1) %>% dplyr::filter(is.na(QC_B_Action) | QC_B_Action!=-1) %>% 
    dplyr::rename(U=Address_A, M=Address_B) %>% 
    dplyr::mutate(
      DESIGN=case_when(Normalization_Bin=='C' ~ 'II',
                       Normalization_Bin=='A' ~ 'I',
                       Normalization_Bin=='B' ~ 'I',
                       # Normalization_Bin=='A' ~ 'IR',
                       # Normalization_Bin=='B' ~ 'IG',
                       TRUE ~ NA_character_),
      COLOR_CHANNEL=case_when(Normalization_Bin=='A' ~ 'Red',
                              Normalization_Bin=='B' ~ 'Grn',
                              TRUE ~ NA_character_),
      col=stringr::str_sub(COLOR_CHANNEL, 1,1),
      Probe_Type=stringr::str_sub(Probe_ID, 1,2),
      Probe_Source=manifest_str,
      Next_Base=case_when(Normalization_Bin=='A' ~ 'T',
                          Normalization_Bin=='B' ~ 'C',
                          TRUE ~ NA_character_)
    ) %>% 
    dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, everything())
  
  if (cleanAdds) {
    man_tib <- man_tib %>%
      dplyr::mutate(M=as.numeric(stringr::str_remove(stringr::str_remove(M, '^1'), '^0+')),
                    U=as.numeric(stringr::str_remove(stringr::str_remove(U, '^1'), '^0+')) )
  }
  
  if (full) {
    ret <- NULL
    ret$aqp <- aqp_tib
    ret$mat <- mat_tib
    ret$ord <- ord_tib
    ret$man <- man_tib
    
    return(ret)
  }
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  
  man_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Manifest I/O::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

idat2manifest = function(sigs, mans, verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'idat2manifest'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret <- NULL
  min_man_tib <- NULL
  min_add_tib <- NULL
  top_platform <- NULL
  top_manifest <- NULL
  
  stime <- system.time({
    min_rec_per <- 0
    sigs_add_cnt <- sigs %>% dplyr::distinct(Address) %>% base::nrow()
    for (cur_man_key in names(mans)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_man_key={cur_man_key}.{RET}"))
      
      cur_key_tib <- stringr::str_split(cur_man_key, pattern='-', simplify=TRUE) %>% as.data.frame() %>% tibble::as_tibble() %>% 
        purrr::set_names(c('platform','manifest')) %>% dplyr::mutate(platform=as.character(platform), manifest=as.character(manifest))
      
      cur_platform <- cur_key_tib %>% dplyr::pull(platform) %>% head(n=1) 
      cur_manifest <- cur_key_tib %>% dplyr::pull(manifest) %>% head(n=1) 
      
      # Probe_ID, Address, Man_Col, Design_Type, Probe_Type
      cur_add_tib <- dplyr::bind_rows(
        dplyr::select(mans[[cur_man_key]], Probe_ID, U, col, DESIGN, Probe_Type) %>% 
          purrr::set_names('Probe_ID','Address','Man_Col','Design_Type','Probe_Type') %>% 
          dplyr::filter(!is.na(Address)) %>%
          dplyr::mutate(Design_Type=dplyr::case_when(!is.na(Man_Col) ~ paste0(Design_Type, Man_Col), TRUE ~ Design_Type)),
        dplyr::select(mans[[cur_man_key]], Probe_ID, M, col, DESIGN, Probe_Type) %>% 
          purrr::set_names('Probe_ID','Address','Man_Col','Design_Type','Probe_Type') %>% 
          dplyr::filter(!is.na(Address)) %>%
          dplyr::mutate(Design_Type=dplyr::case_when(!is.na(Man_Col) ~ paste0(Design_Type, Man_Col), TRUE ~ Design_Type))
      ) %>% 
        dplyr::mutate(Address=stringr::str_remove(Address,'^1') %>% stringr::str_remove('^0+') %>% as.numeric() %>% as.integer() ) %>%
        dplyr::arrange(Probe_ID)
      
      man_add_cnt <- cur_add_tib %>% dplyr::distinct(Address) %>% base::nrow()
      add_mat_cnt <- sigs %>% dplyr::inner_join(cur_add_tib, by="Address") %>% base::nrow()
      add_rec_per <- round(100*add_mat_cnt / base::min(sigs_add_cnt,man_add_cnt), 3)
      
      if (add_rec_per > min_rec_per) {
        min_rec_per <- add_rec_per
        min_man_tib <- mans[[cur_man_key]]
        min_add_tib <- cur_add_tib
        
        top_platform <- cur_platform
        top_manifest <- cur_manifest
      }
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} platform={cur_platform}, manifest={cur_manifest}, ",
                                      "add_mat_cnt={add_mat_cnt}, add_rec_per={add_rec_per}.{RET}"))
    }
    ret$man <- min_man_tib
    ret$add <- min_add_tib
    ret$platform <- top_platform
    ret$manifest <- top_manifest
  })
  if (verbose>vt+4) print(ret)
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; top_platform={top_platform}, top_manifest={top_manifest}, elapsed={etime}.{RET}{RET}"))
  
  ret
}

getManifestList = function(path=NULL, platform=NULL, manifest=NULL, dir=NULL, 
                           verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'getManifestList'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))

  if ( is.null(path) && (is.null(platform) || is.null(manifest) || is.null(dir) ) ) {
    stop(glue::glue("[{funcTag}]: ERROR: Path, platform, manifest, dir can't all be null!{RET}{RET}"))
    return(NULL)
  }
  
  stime <- system.time({
    paths <- NULL
    if (!is.null(path) && file.exists(path)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Match 1; path={path}.{RET}"))
      paths <- list.files(dirname(path), pattern=basename(path), full.names=TRUE)
      
    } else if (!is.null(platform) && !is.null(manifest) && !is.null(dir) ) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Match 2.{RET}"))
      fileName <- paste0(opt$platform,'-',opt$manifest,'.manifest.sesame-base.cpg-sorted.csv.gz')
      filePath <- file.path(dir, fileName)
      stopifnot(file.exists(filePath))
      paths <- list.files(dir, pattern=fileName, full.names=TRUE)
      
    } else if ( !is.null(path) && !is.null(dir) && path=='auto' ) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Match 3.{RET}"))
      paths <- list.files( file.path(dir), pattern=paste0('.manifest.sesame-base.cpg-sorted.csv.gz'), full.names=TRUE )
      
    } else {
      stop(glue::glue("[{funcTag}]: ERROR: Path, platform, manifest, dir FAILED CHECK!{RET}{RET}"))
      return(NULL)
    }

    stopifnot(!is.null(paths))
    stopifnot(length(paths)>0)
    
    man_tibs <- NULL
    for (ii in c(1:length(paths))) {
      cur_key_tib <- paths[ii] %>% basename() %>% stringr::str_remove('\\.manifest.*$') %>% 
        stringr::str_split(pattern='-', simplify=TRUE) %>% as.data.frame() %>% 
        purrr::set_names(c('platform', 'manifest'))
      
      cur_platform <- cur_key_tib$platform
      cur_manifest <- cur_key_tib$manifest
      cur_key <- paste(cur_platform, cur_manifest, sep='-')
      
      man_tibs[[cur_key]] <- suppressMessages(suppressWarnings( readr::read_csv(paths[ii]) )) %>%
        dplyr::mutate(M=as.integer(M), U=as.integer(U))
    }
  })
  if (verbose>vt+4) print(man_tibs)
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  man_tibs
}

loadManifestGenomeStudio = function(file, addSource=FALSE, normalize=FALSE, retType=NULL, max=0,
                                    verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadManifestGenomeStudio'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading file={file}.{RET}"))
  
  ctl_cols <- c('Address', 'Control_Group', 'Control_Color', 'Control_Name')
  
  ret_dat <- NULL
  stime <- system.time({
    tib <- NULL
    man_tib <- NULL
    ctl_tib <- NULL

    src <- base::basename(file) %>% stringr::str_remove('\\.gz$') %>% stringr::str_remove('\\.csv$')
    cmd <- glue::glue("gzip -dc {file} | head -n 50 | grep -n '^IlmnID' ")
    cnt <- suppressMessages(suppressWarnings( system(as.character(cmd), intern=TRUE, ignore.stderr=TRUE) )) %>% 
      stringr::str_remove('\\r$') %>% stringr::str_remove(':.*$') %>% as.integer()
    
    cnt <- cnt - 1
    if (max>0) {
      tib <- suppressMessages(suppressWarnings(readr::read_csv(file, skip=cnt, n_max=max)))
    } else {
      tib <- suppressMessages(suppressWarnings(readr::read_csv(file, skip=cnt)))
    }
    
    tib_len <- base::nrow(tib)
    ctl_idx <- which(tib$IlmnID=='[Controls]') %>% head(n=1) %>% as.integer()
    man_tib <- tib %>% head(n=ctl_idx-1)
    if (ctl_idx>0)
      ctl_tib <- tib %>% tail(tib_len-ctl_idx) %>% dplyr::select(1:4) %>% purrr::set_names(ctl_cols)

    man_tib <- man_tib %>% dplyr::mutate(
      Probe_Type=stringr::str_sub(IlmnID, 1,2),
      Strand_CO=case_when(
        Probe_Type=='cg' ~ 'C',
        Probe_Type=='rs' ~ 'C',
        Probe_Type=='ch' ~ 'O',
        TRUE ~ NA_character_
      )
    )
    man_tib <- man_tib %>% dplyr::mutate(AddressA_ID=stringr::str_remove(AddressA_ID, '^0+') %>% as.double(),
                                         AddressB_ID=stringr::str_remove(AddressB_ID, '^0+') %>% as.double() )
    
    if (normalize) {
      # Check for field that are unique to known manifests::
      if (grep("TopGenomicSeq", names(man_tib)) %>% length() == 1) {
        man_tib <- man_tib %>% 
          dplyr::rename(Top_Sequence=TopGenomicSeq,
                        Genome_Build=GenomeBuild,Chromosome=Chr, Coordinate=MapInfo) %>%
          dplyr::mutate(IlmnStrand=stringr::str_sub(IlmnStrand, 1,1),
                        SourceStrand=stringr::str_sub(SourceStrand, 1,1)) %>%
          dplyr::mutate(Infinium_Design='I') %>% 
          dplyr::select('IlmnID', 'AddressA_ID', 'AlleleA_ProbeSeq', 'AddressB_ID', 'AlleleB_ProbeSeq', 
                        'Top_Sequence', 'SourceSeq', 
                        'Probe_Type','Infinium_Design','Next_Base', 'Color_Channel',
                        'Genome_Build', 'Chromosome', 'Coordinate', 
                        'Strand_CO', 'IlmnStrand', 'SourceStrand')
        
      } else if (grep("Forward_Sequence", names(man_tib)) %>% length() == 1) {
        man_tib <- man_tib %>% 
          dplyr::rename(Strand_FR=Strand,Infinium_Design=Infinium_Design_Type,
                        Chromosome=CHR, Coordinate=MAPINFO, 
                        Strand_FR=Strand) %>%
          dplyr::select('IlmnID', 'AddressA_ID', 'AlleleA_ProbeSeq', 'AddressB_ID', 'AlleleB_ProbeSeq', 
                        'Forward_Sequence', 'SourceSeq', 
                        'Probe_Type','Infinium_Design','Next_Base', 'Color_Channel',
                        'Genome_Build', 'Chromosome', 'Coordinate', 
                        'Strand_FR', 'Strand_CO')
        
      } else {
        cat(glue::glue("[{funcTag}]:{tabsStr} ERROR; Unknown Genome Studio Manifest Type!!!{RET}{RET}"))
        return(NULL)
      }
    }
    
    if (addSource) {
      man_tib <- man_tib %>% dplyr::mutate(Man_Source=!!src)
      ctl_tib <- ctl_tib %>% dplyr::mutate(Man_Source=!!src)
    }
    
    if (is.null(retType)) {
      ret_dat$man <- man_tib
      ret_dat$ctl <- ctl_tib
    } else if (!is.null(retType) && retType=='man') {
      ret_dat <- man_tib
    } else if (!is.null(retType) && retType=='ctl') {
      ret_dat <- ctl_tib
    } else {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Warning: Unsupported retType request={retType}. Defaulting to list.{RET}"))
      ret_dat$man <- man_tib
      ret_dat$ctl <- ctl_tib
    }
    })
  if (verbose>vt+4) print(tib)

  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  ret_dat
}

# TBD:: Should be renamed loadManifestSource -> loadManifestSesame
loadManifestSource = function(file,addSource=FALSE, verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadManifestSource'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading file={file}.{RET}"))
  
  stime <- system.time({
    tib <- NULL
    if ( stringr::str_ends(file, '.tsv') || stringr::str_ends(file, '.tsv.gz') ) {
      tib <- suppressMessages(suppressWarnings(readr::read_tsv(file)))
      source <- base::basename(file) %>% stringr::str_remove('\\.gz$') %>% stringr::str_remove('\\.tsv$')
    } else if ( stringr::str_ends(file, '.csv') || stringr::str_ends(file, '.csv.gz') ) {
      tib <- suppressMessages(suppressWarnings(readr::read_csv(file)))
      source <- base::basename(file) %>% stringr::str_remove('\\.gz$') %>% stringr::str_remove('\\.csv$')
    } else if ( stringr::str_ends(file, '.rds')) {
      tib <- suppressMessages(suppressWarnings(readr::read_rds(file)))
      source <- base::basename(file) %>% stringr::str_remove('\\.rds$')
    } else {
      stop("{RET}[{funcTag}]: ERROR: Unsupported manifest format suffix (only csv/csv.gz or rds): file={file}!!!{RET}{RET}")
    }
    
    # Fix Genecode Fields if its Genecode...
    if ( length( grep('genesUniq', names(tib) ) ) > 0 && 
         length( grep('distToTSS', names(tib) ) ) > 0)
      tib <- tib %>% dplyr::mutate(genesUniq=as.character(genesUniq), distToTSS=as.integer(distToTSS) )
    
    if (addSource) {
      tib <- tib %>% dplyr::mutate(Man_Source=!!source)
    }
  })
  if (verbose>vt+4) print(tib)

  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  tib
}

loadAddressSource = function(file, man, fresh=FALSE, save=TRUE, split=FALSE,
                             verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadAddressSource'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  stime <- system.time({
    tibs <- NULL
    if (!fresh && file.exists(file)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading file={file}.{RET}"))
      ## TBD:: CHeck Type!!!
      if (stringr::str_ends(file,'.rds')) {
        tibs <- suppressMessages(suppressWarnings(readr::read_rds(file) ))
      } else {
        tibs <- suppressMessages(suppressWarnings(readr::read_csv(file) ))
      }
    } else {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building file={file}.{RET}"))
      
      if (split) {
        tibs[['I1']] <- man %>% dplyr::filter(DESIGN=='I') %>% 
          dplyr::select(Probe_ID,Probe_Type,U,M) %>% 
          tidyr::gather(Design_Type,Address, -Probe_Type, -Probe_ID) %>% 
          dplyr::arrange(Address)
        
        tibs[['I2']] <- man %>% dplyr::filter(DESIGN=='II') %>% 
          dplyr::rename(Address=U) %>% dplyr::select(Probe_ID,Probe_Type,Address) %>%
          dplyr::arrange(Address)
      } else {
        tibs <- dplyr::bind_rows(
          man %>% dplyr::filter(DESIGN=='I') %>% 
            dplyr::select(Probe_ID,Probe_Type,U,M, col) %>% 
            tidyr::gather(Design_Type,Address, -Probe_Type, -Probe_ID, -col) %>%
            dplyr::rename(Man_Col=col) %>% 
            dplyr::mutate(Design_Type=paste0(Design_Type,'I')) %>%
            dplyr::select(Probe_ID,Address,Man_Col,Design_Type,Probe_Type),
          
          man %>% dplyr::filter(DESIGN=='II') %>% 
            dplyr::rename(Address=U, Man_Col=col) %>%
            dplyr::mutate(Design_Type='II') %>%
            dplyr::select(Probe_ID,Address,Man_Col,Design_Type,Probe_Type)
        ) %>% dplyr::arrange(Probe_ID)
      }
      
      if (save) {
        if (stringr::str_ends(file,'.rds')) {
          if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving file(RDS)={file}.{RET}"))
          readr::write_rds(tibs, file, compress='gz')
        } else {
          if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving file(CSV)={file}.{RET}"))
          readr::write_csv(tibs, file)
        }
      }
    }
  })
  if (verbose>=vt+4) print(tibs)
  
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  tibs
}

# End of file
