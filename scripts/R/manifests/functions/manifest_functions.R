
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
#                          Standard Function Method::
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
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        AQP/PQC Workflow Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

reducedProbeMapping = function(man,
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'reducedProbeMapping'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

decodeAqpPqcWrapper = function(ord_vec, mat_vec, aqp_vec=NULL, pqc_vec=NULL,
                               platform, version, matFormat="new", sidx=2, plen=50,
                               name=NULL, outDir=NULL, full=FALSE, trim=TRUE,
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'decodeAqpPqcWrapper'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # Load PQC::
    pqc_man_tib <- NULL
    pqc_fix_tib <- NULL
    if (!is.null(pqc_vec) && length(pqc_vec)!=0) {
      pqc_man_tib <- decodeToManifestWrapper(
        ords=ord_vec, mats=mat_vec, pqcs=pqc_vec, aqps=aqp_vec, 
        platform=platform, version=version,
        matFormat=matFormat,
        pqcFormat='pqc',
        full=par$retData, trim=TRUE,
        verbose=verbose,tc=0,tt=pTracker)
      
      # Add Probe Sequences for matching::
      pqc_fix_tib <- 
        fixOrderProbeIDs(pqc_man_tib, verbose=verbose,tc=0,tt=pTracker) %>%
        addReducedMapProbe(sidx=sidx,plen=plen, verbose=verbose,tc=0,tt=pTracker) %>%
        dplyr::select(Seq_ID, FR,TB,CO,PD,Infinium_Design,Mat_PrbA, everything())
    }
    
    # Load AQP::
    aqp_man_tib <- NULL
    aqp_fix_tib <- NULL
    if (!is.null(aqp_vec) && length(aqp_vec)!=0) {
      aqp_man_tib <- decodeToManifestWrapper(
        ords=ord_vec, mats=mat_vec, pqcs=pqc_vec, aqps=aqp_vec, 
        platform=platform, version=version, 
        matFormat=matFormat,
        pqcFormat='aqp',
        full=par$retData, trim=TRUE,
        verbose=verbose,tc=0,tt=pTracker)
      
      # Add Probe Sequences for matching::
      aqp_fix_tib <- 
        fixOrderProbeIDs(aqp_man_tib, verbose=verbose,tc=0,tt=pTracker) %>%
        addReducedMapProbe(sidx=sidx,plen=plen, verbose=verbose,tc=0,tt=pTracker) %>%
        dplyr::select(Seq_ID, FR,TB,CO,PD,Infinium_Design,Mat_PrbA, everything())
    }
    
    # QC Sanity Checks for AQP/PQC if present::
    qc_man_tib  <- NULL
    if (!is.null(pqc_man_tib) && !is.null(aqp_man_tib)) {
      qc_man_tib <- AQP_PQC_QC(aqp=aqp_man_tib, pqc=pqc_man_tib, 
                               aqp_vec=aqp_vec, pqc_vec=pqc_vec,
                               aqp_fix=aqp_fix_tib, pqc_fix=pqc_fix_tib,
                               name=name,outDir=outDir,
                               verbose=verbose,tc=tc+1,tt=pTracker)
    }
    
    # Use the PQC manifest if not use AQP::
    ret_tib <- NULL
    if (length(pqc_vec)!=0) {
      ret_tib <- pqc_fix_tib
    } else if (length(aqp_vec)!=0) {
      ret_tib <- aqp_fix_tib
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Niether PQC or AQP tibble exists!!!}{RET}{RET}"))
      return(NULL)
    }
    ret_tib <- ret_tib %>% dplyr::arrange(Mat_PrbA) %>% dplyr::mutate(CO=stringr::str_remove_all(CO,' '))
    ret_cnt <- ret_tib %>% base::nrow()
    
    if (!is.null(outDir) && !is.null(name)) {
      if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
      ret_csv <- file.path(outDir, paste(name,'manifest.raw.tsv.gz', sep='_') )
      
      cat(glue::glue("[{funcTag}]: Writing Raw Manifest(cnt={ret_cnt})={ret_csv}...{RET}") )
      readr::write_csv(ret_tib,ret_csv)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Quality Control Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

AQP_PQC_QC = function(aqp, pqc, aqp_vec, pqc_vec,aqp_fix, pqc_fix,
                      name=NULL, outDir=NULL,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'AQP_PQC_QC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    aqp_unq_man_tib <- aqp %>% dplyr::anti_join(pqc, by=c("M","U"))
    pqc_unq_man_tib <- pqc %>% dplyr::anti_join(aqp, by=c("M","U"))
    
    aqp_list <- rev(aqp_vec)
    names(aqp_list) <- base::basename(aqp_vec) %>% stringr::str_remove('.gz$') %>% stringr::str_remove('.txt$')
    pqc_name <- base::basename(pqc_vec[1]) %>% stringr::str_remove('.gz$') %>% stringr::str_remove('.txt$')
    
    aqp_unq_tib <- lapply(aqp_list, loadPQC, format='aqp', trim=TRUE) %>% dplyr::bind_rows(.id="AQP_Name") %>%
      dplyr::select(Address,Decode_Status,AQP_Name) %>% dplyr::distinct(Address, .keep_all=TRUE) %>% dplyr::arrange(Address)
    pqc_unq_tib <- loadPQC(pqc_vec[1], format='pqc', trim=TRUE) %>% dplyr::mutate(AQP_Name=pqc_name) %>% 
      dplyr::distinct(Address, .keep_all=TRUE) %>% dplyr::arrange(Address)
    
    # Conclusion:: All missing calls from PQC to AQP are 0 at the most recent AQP stage!!!
    #  - This is what we want to see
    pqc_unq_sum_tib <- dplyr::bind_rows(
      pqc_unq_man_tib %>% dplyr::rename(Decode_Status_PQC=Decode_Status_A) %>%
        dplyr::inner_join(aqp_unq_tib, by=c("U"="Address") ) %>% 
        dplyr::select(Decode_Status_PQC,Decode_Status,AQP) %>% dplyr::group_by_all() %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>%
        dplyr::mutate(Allele="U"),
      pqc_unq_man_tib %>% dplyr::rename(Decode_Status_PQC=Decode_Status_A) %>%
        dplyr::inner_join(aqp_unq_tib, by=c("M"="Address") ) %>% 
        dplyr::select(Decode_Status_PQC,Decode_Status,AQP) %>% dplyr::group_by_all() %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>%
        dplyr::mutate(Allele="M") )
    
    pqc_unq_mis_cnt <- pqc_unq_sum_tib %>% dplyr::filter(Decode_Status_PQC != 0 | Decode_Status != 0) %>% base::nrow()
    if (pqc_unq_mis_cnt!=0) {
      print(pqc_unq_sum_tib)
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed QC Sanity Check; pqc_unq_mis_cnt={pqc_unq_mis_cnt}"))
      return(NULL)
    }
    
    # Conclusion:: Need to investigate more; take away is that the AQP results will have failures
    #  NOTE: its surprising how other tangos come up often... Implying Biziprobe has a ranked order of probes to use
    aqp_unq_sum_tib <- aqp_unq_man_tib %>% 
      dplyr::rename(Decode_Status_AQP=Decode_Status_A) %>%
      dplyr::inner_join(pqc_unq_tib, by=c("U"="Address") ) %>% 
      dplyr::select(Decode_Status_AQP,Decode_Status,AQP) %>% dplyr::group_by_all() %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    
    aqp_unq_mis_cnt <- aqp_unq_sum_tib %>% dplyr::filter(Decode_Status_AQP != 0 & Decode_Status != 0 & Decode_Status != 1 ) %>% base::nrow()
    if (aqp_unq_mis_cnt!=0) {
      print(pqc_unq_sum_tib)
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed QC Sanity Check; aqp_unq_mis_cnt={aqp_unq_mis_cnt}"))
      return(NULL)
    }
    if (verbose>=vt) cat(glue::glue("[{funcTag}]: Passed all AQP/PQC Sanity Validation; ",
                                    "pqc_unq_mis_cnt={pqc_unq_mis_cnt}, aqp_unq_mis_cnt={aqp_unq_mis_cnt}.{RET}"))
    
    fix_man_sum_tib <- dplyr::full_join(
      dplyr::group_by(pqc_fix,Probe_Type) %>% dplyr::summarise(Count=n(), .groups="drop"),
      dplyr::group_by(aqp_fix,Probe_Type) %>% dplyr::summarise(Count=n(), .groups="drop"),
      by="Probe_Type", suffix=c("_PQC","_AQP")) %>%
      dplyr::mutate(Count_Dif=Count_AQP-Count_PQC)
    fix_man_mis_cnt <- fix_man_sum_tib %>% dplyr::filter(Count_Dif<=0) %>% base::nrow()
    
    if (fix_man_mis_cnt>0) {
      fix_man_sum_tib %>% print()
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed QC Sanity Check; PQC>AQP; fix_man_mis_cnt={fix_man_mis_cnt}"))
      return(NULL)
    }
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]: Passed all QC Sanity Check; PQC>AQP Validation; fix_man_mis_cnt={fix_man_mis_cnt}.{RET}"))
    
    ret_tib <- tibble::tibble(aqp_unq_mis_cnt=aqp_unq_mis_cnt,
                              pqc_unq_mis_cnt=pqc_unq_mis_cnt,
                              fix_man_mis_cnt=fix_man_mis_cnt )
    if (!is.null(outDir) && !is.null(name)) {
      if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
      ret_csv <- file.path(outDir,paste(name,'AQP-PQC.quality-control.summary.csv.gz',sep='_') )
      readr::write_csv(ret_tib,ret_csv)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Wrote QC Summary={ret_csv}.{RET}"))
    }
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Specialized Mouse Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fixOrderProbeIDs = function(tib, field="Probe_Type", sidx=2,plen=50,
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'fixOrderProbeIDs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # Add Infinium Type Definition:: Also, fail any probe that breaks this method
    tib <- tib %>% dplyr::mutate(
      M=dplyr::case_when(M=='NA' ~ NA_real_, TRUE ~ M),
      Infinium_Design=dplyr::case_when(
        !is.na(AlleleA_Probe_Sequence) & !is.na(AlleleB_Probe_Sequence) & !is.na(U) & !is.na(M) ~ 1,
        !is.na(AlleleA_Probe_Sequence) &  is.na(AlleleB_Probe_Sequence) & !is.na(U) &  is.na(M) ~ 2,
        TRUE ~ NA_real_
      )
    )
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} tib(with Inf Des)=.{RET}"))
    if (verbose>=vt+4) print(tib)
    
    # Add Probe Rep_Num, Old_Probe_ID, and Probe Lengths
    tib <- tib %>%
      dplyr::filter(!is.na(Infinium_Design)) %>%
      dplyr::add_count(Probe_ID, name="Rep_Max") %>%
      dplyr::group_by(Probe_ID) %>% dplyr::mutate(Rep_Num=dplyr::row_number()) %>% dplyr::ungroup() %>%
      dplyr::mutate(AlleleA_Probe_Length=stringr::str_length(AlleleA_Probe_Sequence),
                    AlleleB_Probe_Length=stringr::str_length(AlleleB_Probe_Sequence),
                    Old_Probe_ID=paste(Probe_ID,paste0('r',Rep_Num), sep='_'))
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} tib(with Rep, Old-ID, Len)=.{RET}"))
    if (verbose>=vt+4) print(tib)
    
    tib <- tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2))
    
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
    
    for (type in names(tib_types)) {
      cur_tib <- tib_types[[type]]
      
      if (type=='cg' || type=='ch' || type=='rp') {
        cur_tib <- cur_tib %>% 
          dplyr::mutate(
            Probe_ID=stringr::str_replace(Probe_ID, '_([FR])_([CO]_[I]+)$',  '_\\$1_N_\\$2') %>% stringr::str_remove_all('\\\\')
          ) %>%
          tidyr::separate(Probe_ID, into=c('Seq_ID','FR','TB','CO','PD'), sep='_', remove=FALSE)
        
      } else if (type=='mu') {
        # Clean up for some mm10 mu early designs...
        cur_tib <- dplyr::bind_rows(
          dplyr::filter(cur_tib, !stringr::str_detect(Probe_ID, '_s[0-9]')) %>%
            dplyr::mutate(Probe_ID=stringr::str_replace(Probe_ID, '_(I.*)$', '_s0_h0_\\$1') %>% stringr::str_remove('\\\\')),
          dplyr::filter(cur_tib, stringr::str_detect(Probe_ID, '_s[0-9]')) )
        
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
      } else if (type=='BS' || type=='bs') {
        cur_tib <- cur_tib %>% dplyr::mutate(
          Di=stringr::str_replace(Probe_ID, '^.*-([ACTG][ACTG])-.*$', '\\$1') %>% stringr::str_remove_all('\\\\')
        ) %>% 
          dplyr::mutate(Probe_ID=stringr::str_replace(Probe_ID, '^BSC_', 'bs-')) %>% 
          tidyr::separate(Probe_ID, into=c('Seq_ID','FR','TB','CO','PD'), sep='_', remove=FALSE) 
      } else if (type=='NO' || type=='no') {
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
    # ret_tib <- addReducedMapProbe(tib=ret_tib,sidx=sidx, plen=plen) %>%
    #   dplyr::select(Seq_ID, FR,TB,CO,PD,Infinium_Design,Mat_PrbA,Mat_Prb, everything())
    
    # Add Probe Sequences for matching::
    # ret_tib <- ret_tib %>% dplyr::mutate(
    #   Mat_Prb=dplyr::case_when(
    #     Infinium_Design==1 ~ stringr::str_sub(AlleleA_Probe_Sequence, 2,50),
    #     Infinium_Design==2 ~ stringr::str_sub(AlleleA_Probe_Sequence, 3,50),
    #     TRUE ~ NA_character_
    #   ) %>% stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% stringr::str_replace_all('Y', 'T'),
    #   
    #   Mat_PrbA=dplyr::case_when(
    #     Infinium_Design==1 ~ stringr::str_sub(AlleleA_Probe_Sequence, 2,49),
    #     Infinium_Design==2 ~ stringr::str_sub(AlleleA_Probe_Sequence, 3,50),
    #     TRUE ~ NA_character_
    #   ) %>% stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% stringr::str_replace_all('Y', 'T')
    # ) %>% dplyr::select(Seq_ID, FR,TB,CO,PD,Infinium_Design,Mat_PrbA,Mat_Prb, everything())
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tib({ret_cnt})::{RET}"))
    if (verbose>=vt+4) print(ret_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

addReducedMapProbe = function(tib, sidx=2, plen=50,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'template_func'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; Assumed Probe Length={plen}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    idx1 <- sidx
    len1 <- plen - 1
    idx2 <- sidx + 1
    len2 <- plen
    
    ret_tib <- tib %>% dplyr::mutate(
      Mat_PrbA=dplyr::case_when(
        Infinium_Design==1 ~ stringr::str_sub(AlleleA_Probe_Sequence, idx1,len1),
        Infinium_Design==2 ~ stringr::str_sub(AlleleA_Probe_Sequence, idx2,len2),
        TRUE ~ NA_character_
      ) %>% stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% stringr::str_replace_all('Y', 'T')
      # Mat_Prb=dplyr::case_when(
      #   Infinium_Design==1 ~ stringr::str_sub(AlleleA_Probe_Sequence, 2,plen),
      #   Infinium_Design==2 ~ stringr::str_sub(AlleleA_Probe_Sequence, 3,plen),
      #   TRUE ~ NA_character_
      # ) %>% stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% stringr::str_replace_all('Y', 'T'),
    )
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tib({ret_cnt})::{RET}"))
    if (verbose>=vt+4) print(ret_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
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
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    field <- rlang::sym(field)
    inf1_cnt <- man %>% dplyr::group_by(!!field) %>% dplyr::summarise(Count=n()) %>% 
      dplyr::filter(!!field=='I') %>% dplyr::summarise(Sum_Count=sum(Count)) %>% dplyr::pull(Sum_Count) %>% as.integer()
    inf2_cnt <- man %>% dplyr::group_by(!!field) %>% dplyr::summarise(Count=n()) %>% 
      dplyr::filter(!!field=='II') %>% dplyr::summarise(Sum_Count=sum(Count)) %>% dplyr::pull(Sum_Count) %>% as.integer()
    
    bead1_cnt <- inf1_cnt * 2
    bead2_cnt <- inf2_cnt
    beads_cnt <- bead1_cnt + bead2_cnt
    beads_per <- bead1_cnt / beads_cnt
    
    ret_tib <- tibble::tibble(Inf1=inf1_cnt, Inf2=inf2_cnt, 
                              Bead1=bead1_cnt, Bead2=bead2_cnt, Beads=beads_cnt,
                              BeadPerc=beads_per)
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tib({ret_cnt})::{RET}"))
    if (verbose>=vt+4) print(ret_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Order, Match, PQC/AQP File IO Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadORD = function(file, format='old', skip=8, guess=50000, trim=FALSE,
                   verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadORD'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Staring; format={format}; skip={skip}, file={file}.{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    if (format=='old') {
      ret_tib <- suppressMessages(suppressWarnings( readr::read_csv(file, skip=skip, guess_max=guess)  )) %>%
        dplyr::mutate(AlleleB_Probe_Id=as.character(AlleleB_Probe_Id), 
                      AlleleB_Probe_Sequence=as.character(AlleleB_Probe_Sequence))
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported format={format}!!!{RET}{RET}"))
      return(NULL)
    }
    ret_tib <- ret_tib %>% dplyr::filter(!is.na(Assay_Design_Id))
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tib({ret_cnt})::{RET}"))
    if (verbose>=vt+4) print(ret_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

loadMAT = function(file, format='new', skip=40, guess=1000, trim=FALSE,
                   verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadMAT'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Staring; format={format}; skip={skip}, file={file}.{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  match_cnames <- c("Plate", "Row", "Col", "Address", "Mod5", "Full_Seq", "Mod3", "Comment")
  
  stime <- system.time({
    if (format=='old') {
      ret_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file, skip=skip, guess_max=guess) )) %>%
        dplyr::rename(Probe_ID=probe_id, Full_Seq=bo_seq, Address=address_name, Ord_Seq=sequence) %>% 
        dplyr::mutate(Tango_Seq=stringr::str_sub(Full_Seq,1,45),
                      Probe_Seq=stringr::str_sub(Full_Seq,46))
    } else if (format=='new') {
      ret_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file, skip=skip, guess_max=guess) )) %>% 
        purrr::set_names(match_cnames) %>% 
        dplyr::mutate(Tango_Seq=stringr::str_sub(Full_Seq,1,45), 
                      Probe_Seq=stringr::str_sub(Full_Seq,46))
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported format={format}!!!{RET}{RET}"))
      return(NULL)
    }
    ret_tib <- ret_tib %>% dplyr::filter(!is.na(Address))
    if (trim) {
      trim_cnt <- ret_tib %>% dplyr::filter(!stringr::str_starts(Address, '1')) %>% base::nrow()
      if (trim_cnt==0) {
        ret_tib <- ret_tib %>% dplyr::mutate(Address=as.numeric(stringr::str_remove(stringr::str_remove(Address, '^1'), '^0+')) )
      } else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: Attempting to trim new tango fomrat, ",
                        "but format doesn't match; trim_cnt={trim_cnt}!!!{RET}{RET}"))
        return(NULL)
      }
    }
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tib({ret_cnt})::{RET}"))
    if (verbose>=vt+4) print(ret_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

loadPQC = function(file, format='pqc', skip=7, guess=1000, trim=FALSE,
                   verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadPQC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Staring; format={format}; skip={skip}, file={file}.{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  
  stime <- system.time({
    if (format=='aqp') {
      ret_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file, skip=skip, guess_max=guess) )) %>%
        dplyr::select(Address, Decode_Status)
    } else if (format=='pqc') {
      ret_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file, skip=skip, guess_max=guess) )) %>% 
        purrr::set_names(stringr::str_replace_all(names(.),' ','_')) %>%
        dplyr::rename(Decode_Status=Status) %>%
        dplyr::select(Address, Decode_Status)
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported format={format}!!!{RET}{RET}"))
      return(NULL)
    }
    ret_tib <- ret_tib %>% dplyr::filter(!is.na(Address)) %>% dplyr::filter(!is.na(Decode_Status)) %>%
      dplyr::ungroup()
    if (trim) {
      trim_cnt <- ret_tib %>% dplyr::filter(!stringr::str_starts(Address, '1')) %>% base::nrow()
      if (trim_cnt==0) {
        ret_tib <- ret_tib %>% dplyr::mutate(Address=as.numeric(stringr::str_remove(stringr::str_remove(Address, '^1'), '^0+')) )
      } else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: Attempting to trim new tango fomrat, ",
                        "but format doesn't match; trim_cnt={trim_cnt}!!!{RET}{RET}"))
        return(NULL)
      }
    }
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tib({ret_cnt})::{RET}"))
    if (verbose>=vt+4) print(ret_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   AQP Manifest Generation Wrappers::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

decodeToManifest = function(ord, mat, pqc, platform, version, round=NULL,
                            ordFormat='old', ordSkip=8,  ordGuess=50000,
                            matFormat='new', matSkip=40, matGuess=1000,
                            pqcFormat='pqc', pqcSkip=7,  pqcGuess=1000,
                            full=FALSE, trim=TRUE,
                            verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'decodeToManifest'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  ret <- NULL
  ret_tib <- NULL
  stime <- system.time({
    
    manifest_str <- paste(platform,version, sep='-')
    
    ord_tib <- NULL
    mat_tib <- NULL
    pqc_tib <- NULL
    ord_tib <- loadORD(file=ord, format=ordFormat, skip=ordSkip, guess=ordGuess, trim=trim,
                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    mat_tib <- loadMAT(file=mat, format=matFormat, skip=matSkip, guess=matGuess, trim=trim,
                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    pqc_tib <- loadPQC(file=pqc, format=pqcFormat, skip=pqcSkip, guess=pqcGuess, trim=trim,
                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    if (full) {
      ret$ord <- ord_tib
      ret$mat <- mat_tib
      ret$aqp <- pqc_tib
    }
    
    # Reformat Order File::
    #  - Upper case sequences
    ord_tib <- dplyr::rename(ord_tib, Probe_ID=Assay_Design_Id) %>% 
      dplyr::select(-AlleleA_Probe_Id, -AlleleB_Probe_Id) %>%
      dplyr::rename(AlleleA_Probe_Sequence_MIX=AlleleA_Probe_Sequence,
                    AlleleB_Probe_Sequence_MIX=AlleleB_Probe_Sequence ) %>%
      dplyr::mutate(AlleleA_Probe_Sequence=stringr::str_to_upper(AlleleA_Probe_Sequence_MIX),
                    AlleleB_Probe_Sequence=stringr::str_to_upper(AlleleB_Probe_Sequence_MIX) )
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Order Tibble={RET}"))
    if (verbose>=vt+4) print(ord_tib)
    
    # Join Match with PQC Results into Experimental tibble::
    #
    exp_tib <- dplyr::inner_join(mat_tib,pqc_tib, by="Address") %>%
      dplyr::rename(Probe_Seq_Mix=Probe_Seq) %>%
      dplyr::mutate(Probe_Seq=stringr::str_to_upper(Probe_Seq_Mix)) %>%
      dplyr::select(Address, Tango_Seq, Probe_Seq, Decode_Status)
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Experiment Tibble={RET}"))
    if (verbose>=vt+4) print(exp_tib)
    
    # Build Manifest Tibble::
    #
    ret_tib <- NULL
    ret_tib <- ord_tib %>%
      dplyr::left_join(exp_tib, by=c("AlleleA_Probe_Sequence"="Probe_Seq")) %>% 
      dplyr::rename(U=Address, Decode_Status_A=Decode_Status, Address_Seq_A=Tango_Seq) %>%
      dplyr::left_join(exp_tib, by=c("AlleleB_Probe_Sequence"="Probe_Seq")) %>% 
      dplyr::rename(M=Address, Decode_Status_B=Decode_Status, Address_Seq_B=Tango_Seq) %>%
      dplyr::filter(!is.na(Decode_Status_A)) %>%
      dplyr::filter(Decode_Status_A!=-1) %>% 
      dplyr::filter(is.na(Decode_Status_B) | Decode_Status_B!=-1) %>%
      dplyr::select(Probe_ID, M, U, everything() )
    if (!is.null(round)) ret_tib <- ret_tib %>% dplyr::mutate(BP=round, AQP=round)
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Manifest Tibble={RET}"))
    if (verbose>=vt+4) print(ret_tib)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

decodeToManifestWrapper = function(ords, mats, pqcs=NULL, aqps=NULL, platform, version,
                                   ordFormat='old', ordSkip=8,  ordGuess=50000,
                                   matFormat='new', matSkip=40, matGuess=1000,
                                   pqcFormat='pqc', pqcSkip=7,  pqcGuess=1000,
                                   full=FALSE, trim=TRUE,
                                   verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'decodeToManifestWrapper'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  if (is.null(pqcs) && is.null(aqps)) {
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Both AQPs and PQCs can't be null!!!{RET}{RET}"))
    return(ret_tib)
  }
  
  stime <- system.time({
    
    ord_cnt <- length(ords)
    for (idx in c(1:ord_cnt)) {
      pqc_file <- NULL
      if (pqcFormat=='pqc') {
        pqc_file <- pqcs[1]
      } else if (pqcFormat=='aqp') {
        pqc_file <- aqps[idx]
      } else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: Unrecognized pqc format={pqcFormat}!!!{RET}{RET}"))
        return(NULL)
      }
      if (is.null(pqc_file)) {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: PQC File Logic Failed useAQP={useAQP}!!!{RET}{RET}"))
        return(NULL)
      }
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Round={idx}, pqc={pqc_file}.{RET}"))
      
      ret_tib <- ret_tib %>% dplyr::bind_rows(
        decodeToManifest(ord=ords[idx], mat=mats[idx], pqc=pqc_file,
                         platform=platform, version=version, round=idx,
                         
                         ordFormat=ordFormat, ordSkip=ordSkip, ordGuess=ordGuess,
                         matFormat=matFormat, matSkip=matSkip, matGuess=matGuess,
                         pqcFormat=pqcFormat, pqcSkip=pqcSkip, pqcGuess=pqcGuess,
                         
                         full=full, trim=trim,
                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      )
    }
    pre_cnt <- ret_tib %>% base::nrow()
    ret_tib <- ret_tib %>% dplyr::arrange(-AQP) %>% dplyr::group_by(U) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup()
    ret_cnt <- ret_tib %>% base::nrow()
    del_cnt <- pre_cnt - ret_cnt
    cat(glue::glue("[{funcTag}]: Full AQP History={pre_cnt}, AQP Ordered Unique={ret_cnt}, delta={del_cnt}.{RET}") )
    
    # QC Sanity Checks below::
    #
    tan_mat_cnt <- ret_tib %>% dplyr::filter(M==U) %>% base::nrow()
    if (tan_mat_cnt!=0) {
      ret_tib %>% dplyr::filter(M==U) %>% print()
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Non-zero value for M=U tangos={tan_mat_cnt}!!!{RET}{RET}"))
      return(NULL)
    }
    non_unq_cnt <- ret_tib %>% dplyr::add_count(U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, name="Unique_Tan_Seq_Count") %>% 
      dplyr::filter(Unique_Tan_Seq_Count!=1) %>% base::nrow()
    if (non_unq_cnt!=0) {
      ret_tib %>% dplyr::add_count(M,U,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, name="Unique_Tan_Seq_Count") %>% 
        dplyr::filter(Unique_Tan_Seq_Count!=1) %>% print()
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Non-zero value for unique(U,M,PrbU,PrbM) tangos={non_unq_cnt}!!!{RET}{RET}"))
      return(NULL)
    }
    if (verbose>=vt) cat(glue::glue("[{funcTag}]: Valid QC Sanity; tan_mat_cnt={tan_mat_cnt}, non_unq_cnt={non_unq_cnt}.{RET}") )
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
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
      fileName <- paste0(platform,'-',manifest,'.manifest.sesame-base.cpg-sorted.csv.gz')
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
