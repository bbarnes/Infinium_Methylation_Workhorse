

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Sample Sheet Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Sample Sheet I/O Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadAutoSampleSheets = function(dir, platform=NULL, manifest=NULL, suffix='AutoSampleSheet.csv.gz', 
                                addSampleName=FALSE, addPathsCall=FALSE, addPathsSigs=FALSE,
                                flagDetectPval=FALSE, flagSampleDetect=FALSE, flagRefMatch=FALSE,
                                pvalDetectMinKey=NULL, pvalDetectMinVal=0,
                                
                                dbMin=90, r2Min=0.9,
                                dbKey='AutoSample_dB_Key', r2Key='AutoSample_R2_Key',
                                dbVal='AutoSample_dB_Val', r2Val='AutoSample_R2_Val',
                                verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadAutoSampleSheets'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, dir={dir}.{RET}"))
  
  stime <- system.time({
    dat <- NULL
    
    if (!is.null(pvalDetectMinKey)) pvalDetectMinKey <- pvalDetectMinKey %>% rlang::sym()
    # pvalDetectMinVal <- pvalDetectMinVal %>% rlang::sym()
    
    pattern <- NULL
    if (!is.null(platform)) pattern <- paste(pattern,platform, sep='_')
    if (!is.null(manifest)) pattern <- paste(pattern,manifest, sep='_')
    pattern <- paste(pattern,suffix, sep='_')
    
    auto_ss_list <- list.files(dir, pattern=pattern, recursive=TRUE, full.names=TRUE)
    auto_ss_llen <- auto_ss_list %>% length()
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, SampleSheetCount={auto_ss_llen}, pattern={pattern} dir={dir}.{RET}"))
    stopifnot(auto_ss_llen>0)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Load Samples::
    auto_ss_tibs <- suppressMessages(suppressWarnings(lapply(auto_ss_list, readr::read_csv) )) %>% 
      dplyr::bind_rows()
    
    auto_ss_tlen <- base::nrow(auto_ss_tibs)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(bPool)={auto_ss_tlen}{RET}"))
    # print(auto_ss_tibs)
    
    if (addSampleName || flagDetectPval || flagSampleDetect || flagRefMatch) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Adding flags...{RET}"))
      if (verbose>=vt+4) print(auto_ss_tibs)
      
      dbKey <- dbKey %>% rlang::sym()
      dbVal <- dbVal %>% rlang::sym()
      r2Key <- r2Key %>% rlang::sym()
      r2Val <- r2Val %>% rlang::sym()
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #  Add General Sample Name Field::
      if (addSampleName) {
        auto_ss_tibs <- auto_ss_tibs %>% dplyr::mutate(Auto_Sample_Name=!!dbKey) %>%
          dplyr::select(Auto_Sample_Name, dplyr::everything())
      }
      if (verbose>=vt+4) print(auto_ss_tibs)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #  Flag Probe Detected (pval)::
      if (flagDetectPval && !is.null(pvalDetectMinKey) && !is.null(pvalDetectMinVal)) {
        fail_tag <- paste0("Failed<",pvalDetectMinVal)
        pass_tag <- paste0("Passed>=",pvalDetectMinVal)
        
        auto_ss_tibs <- auto_ss_tibs %>% 
          dplyr::mutate(
            detectPvalCase=case_when(
              !!pvalDetectMinKey < !!pvalDetectMinVal ~ fail_tag,
              TRUE ~ pass_tag),
            detectPvalFlag=case_when(
              !!pvalDetectMinKey < !!pvalDetectMinVal ~ FALSE,
              TRUE ~ TRUE)
          )
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #  Flag Auto-Detected Samples::
      if (flagSampleDetect) {
        fail_r2  <- paste0('Failed_r2<',r2Min)
        fail_db  <- paste0('Failed_db<',dbMin)
        fail_mt  <- paste0('Failed_r21-db')
        pass_tag <- 'Passed'
        
        auto_ss_tibs <- auto_ss_tibs %>% dplyr::mutate(detectedSample=case_when(
          !!r2Val < r2Min ~ fail_r2,
          !!dbVal < dbMin ~ fail_db,
          !!dbKey != !!r2Key ~ fail_mt,
          TRUE ~ pass_tag)
        )
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #  Flag Auto Detections with poor matching (pval)::
      if (flagRefMatch) {
        auto_ss_tibs <- auto_ss_tibs %>% dplyr::filter() %>% dplyr::filter(!!dbVal >= dbMin) %>% dplyr::filter(!!r2Val >= r2Min)
        auto_ss_flen <- base::nrow(auto_ss_tibs)
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(filt)={auto_ss_flen}{RET}"))
        
        # Remove Unidentifiable
        # if (rmOdd) auto_ss_tibs <- auto_ss_tibs %>% dplyr::filter(Bead_Pool!='Odd')
      }
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Add Paths::
    if (addPathsCall) {
      auto_ss_tibs <- auto_ss_tibs %>%
        addPathsToSampleSheet(dir=dir, platform=platform, manifest=manifest, del='_',
                              field='Calls_Path', suffix='calls.csv.gz$', verbose=verbose)
      
      auto_ss_tlen <- base::nrow(auto_ss_tibs)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(paths)={auto_ss_tlen}{RET}"))
    }
    if (addPathsSigs) {
      auto_ss_tibs <- auto_ss_tibs %>%
        addPathsToSampleSheet(dir=dir, platform=platform, manifest=manifest, del='_',
                              field='Sigs_Path', suffix='sigs.csv.gz$', verbose=verbose)
      
      auto_ss_tlen <- base::nrow(auto_ss_tibs)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(paths)={auto_ss_tlen}{RET}"))
    }
    
    dat <- auto_ss_tibs
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Sample Sheet Manipulation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getUniqueFields = function(ss, keys, verbose=0,vt=1,tc=1) {
  funcTag <- 'getUniqueFields'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  # Remove Variables with a single value
  uniq_exp_tib <- ss %>% dplyr::select(keys) %>% 
    # dplyr::summarise_each(list(n_distinct) ) %>% 
    dplyr::summarise_all(list(n_distinct) ) %>% 
    tidyr::gather() %>% dplyr::filter(value>1) %>%
    dplyr::arrange(value) %>% 
    tidyr::spread(key, value)
  
  uniq_exp_keys <- uniq_exp_tib %>% names()
  
  # Special Treatment for Chip_Format (put in front)
  field <- 'Chip_Format'  
  isPresent <- grep(field, uniq_exp_keys) %>% length() > 0
  if (isPresent)
    uniq_exp_tib <- uniq_exp_tib %>% dplyr::select(field, everything())
  
  # Special Treatment for Bead_Pool (put in last)
  field <- 'Bead_Pool'  
  isPresent <- grep(field, uniq_exp_keys) %>% length() > 0
  if (isPresent)
    uniq_exp_tib <- uniq_exp_tib %>% dplyr::select(-field, everything())
  
  uniq_exp_keys <- uniq_exp_tib %>% names()
  
  uniq_exp_keys
}

addPathsToSampleSheet = function(ss, dir, platform=NULL, manifest=NULL, field, suffix, del='.',
                                 verbose=0,vt=4,tc=1) {
  funcTag <- 'addPathsToSampleSheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; suffix={suffix}.{RET}"))
  
  pattern <- NULL
  if (!is.null(platform)) pattern <- paste(pattern,platform, sep=del)
  if (!is.null(manifest)) pattern <- paste(pattern,manifest, sep=del)
  pattern <- paste(pattern,suffix, sep=del)
  
  file_list <- NULL
  file_list <- list.files(dir, pattern=pattern, recursive=TRUE, full.names=TRUE)
  file_cnts <- length(file_list)
  if (file_cnts==0) file_list <- NULL
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Found; pattern={pattern}, file_cnts={file_cnts}.{RET}"))
  
  ss_path_tibs <- tibble::tibble(Sentrix_Name=basename(file_list) %>% 
                                   stringr::str_replace('^([^_]+_[^_]+)_.*$', '\\$1') %>%
                                   stringr::str_remove_all('\\\\'),
                                 !!field := file_list)
  ss_path_cnts <- ss_path_tibs %>% base::nrow()
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Foubd; ss_path_cnts={ss_path_cnts}.{RET}"))
  if (verbose>=vt+3) print(ss_path_tibs)
  if (verbose>=vt+3) print(ss)
  
  tib <- ss_path_tibs %>%
    dplyr::distinct(Sentrix_Name, .keep_all=TRUE) %>% 
    dplyr::inner_join(ss, by="Sentrix_Name")

  tib_len <- tib %>% base::nrow()
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; tib_len={tib_len}.{RET}{RET}"))

  tib
}

addBeadPoolToSampleSheet = function(ss, field, verbose=0, vt=1, tc=1) {
  funcTag <- 'addBeadPoolToSampleSheet'
  
  field <- field %>% rlang::sym()
  ss <- ss %>% dplyr::mutate(
    Bead_Pool=case_when(!!field >862926+1000 ~ 'EPIC_Plus',
                        !!field >862926-1000 ~ 'EPIC',
                        !!field >612329-1000 ~ 'BP1234',
                        !!field >448696-1000 ~ 'BP123',
                        !!field >148977-1000 ~ 'BP2',
                        !!field >103113-1000 ~ 'Odd',
                        !!field<=103113-1000 ~ 'UNK',
                        TRUE ~ NA_character_) ) %>%
    dplyr::select(Bead_Pool, everything())
  
  ss
}

listChipIds = function(dir, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'listChipIds'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, dir={dir}.{RET}"))
  
  tib <- NULL
  dirs <- list.dirs(dir, recursive = FALSE)
  for (d in dirs) {
    exp_name <- base::basename(d)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Experiment={exp_name}, dir={d}.{RET}"))
    
    chip_paths <- list.dirs(d, recursive = FALSE)
    chip_names <- chip_paths %>% base::basename()
    
    cur_tib <- tibble::tibble(Build_Paths=chip_paths, 
                              Sentrix_Barcode=chip_names,
                              Experiment_Name=exp_name) %>%
      dplyr::filter(Sentrix_Barcode!='shells') %>%
      dplyr::mutate(Sentrix_Barcode=as.double(Sentrix_Barcode))
    
    # cur_tib %>% print()
    tib <- tib %>% dplyr::bind_rows(cur_tib)
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  tib <- tib %>% dplyr::select(Sentrix_Barcode, Experiment_Name, everything())
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Unique Format Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

formatSS_COVID = function(csv, skip=0, addID,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'formatSS1'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (!file.exists(csv))
    stop(glue::glue("{RET}[{funcTag}]: Failed to find humman annotation sample sheet={csv}!!!{RET}{RET}"))
  
  ss_tib <- suppressMessages(suppressWarnings( readr::read_csv(csv, skip=skip) )) %>% 
    purrr::set_names(stringr::str_replace_all(names(.), ' ', '_'))
  
  if (addID) ss_tib <- ss_tib %>% dplyr::mutate(Sample_ID=Sample_Name)
  if (verbose>=vt+2) print(ss_tib)
  
  ss_tib <- ss_tib %>%
    dplyr::rename(COVID_Status=COVID_status) %>%
  dplyr::mutate(
    Sentrix_Name=paste(Sentrix_ID,Sentrix_Position, sep='_'),
    Sample_Class=dplyr::case_when(
      Sample_Group=='Control'       ~ 'nSARSCov2',
      Sample_Group=='Negative'      ~ 'nSARSCov2',
      
      Sample_Group=='Case'          ~ 'pSARSCov2',
      Sample_Group=='Positive'      ~ 'pSARSCov2',
      
      Sample_Group=='Hyper_Control' ~ 'T99UCD',
      Sample_Group=='Hypo_Control'  ~ 'T00UCD',
      TRUE ~ NA_character_
    ),
    Source_Sample_ID=dplyr::case_when(
      Sample_Group=='Hyper_Control' ~ 'hih',
      Sample_Group=='Hypo_Control' ~ 'low',
      TRUE ~ Sample_ID
    ),
    Source_Sample_Name=Sample_Name,
    Sample_Name=stringr::str_sub(paste0('S',stringr::str_pad(Sample_Name, 3, side="left", pad='0') ), 1, 4)
  ) %>% 
    dplyr::mutate_if(is.character, spaceToUnderscore ) %>% 
    dplyr::select(Sentrix_Name,Sample_Class,COVID_Status,Source_Sample_ID,Source_Sample_Name,Tissue_Source,
                  Sample_Well,Sample_Plate,Sample_Group,Species)
  
  ss_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Merge Calls Files Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

mergeCallsFromSS = function(ss, max=0, outName, outDir, 
                            verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'mergeCallsFromSS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))

  files_tib <- tibble::tibble()
  file_list <- base::list()
  data_list <- base::list()
  stime <- system.time({

    ss_cnt <- ss %>% base::nrow()
    for (ssIdx in c(1:ss_cnt)) {
      sentrix_name <- ss$Sentrix_Name[ssIdx]
      call_path <- ss$Calls_Path[ssIdx]
      ss_perc = round(100*ssIdx/ss_cnt, 3)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading Calls; idx={ssIdx}/{ss_cnt}={ss_perc}, sentrix_name={sentrix_name}; path={call_path}.{RET}"))

      call_tib  <- suppressMessages(suppressWarnings( readr::read_csv(call_path) ))
      call_cols <- call_tib %>% dplyr::select(-1) %>% names()
      call_lens <- length(call_cols)
      
      for (cIdx in (1:call_lens)) {
        cur_col <- call_cols[cIdx]
        cur_tib <- dplyr::select(call_tib, 1,cur_col)
        colnames(cur_tib)[2] <- sentrix_name
        if (is.null(data_list[[cur_col]])) {
          data_list[[cur_col]] <- cur_tib
        } else {
          data_list[[cur_col]] <- dplyr::full_join(data_list[[cur_col]],cur_tib, by="Probe_ID")
        }
      }
      
      if (max>0 && ssIdx>=max) break
    }
    
    for (name in names(data_list)) {
      out_name <- paste(outName,name, sep='_')
      file_list[[name]] <- file.path(outDir, paste(out_name,'raw-data.csv.gz', sep='.') )
      cat(glue::glue("[{funcTag}]: Writing Calls File; name={name}, out_name={out_name}; CSV={file_list[[name]]}...{RET}"))
      
      readr::write_csv(data_list[[name]], file_list[[name]])
    }
    
    # Build return tibble of mreged/split calls files::
    files_tib <- file_list %>% dplyr::bind_rows() %>% tidyr::gather(Method, Full_Path)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  files_tib
}

loadCallsMatrix = function(betaCSV, pvalCSV, minPval=NULL, mat=NULL, cgn=NULL, ss=NULL, 
                           idKey="Probe_ID", addPval=FALSE, betaName='beta', pvalName='pval', del='.',
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadCallsMatrix'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  if (!file.exists(betaCSV)) {
    cat(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: File does not exist; betaCSV={betaCSV}.{RET}"))
    return(NULL)
  }
  if (!file.exists(pvalCSV)) {
    cat(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: File does not exist; pvalCSV={pvalCSV}.{RET}"))
    return(NULL)
  }
  stopifnot(file.exists(betaCSV))
  stopifnot(file.exists(pvalCSV))
  
  if (!is.null(minPval)) minPval <- as.double(minPval)
  
  stime <- system.time({
    
    idKey_sym <- rlang::sym(idKey)
    
    beta_tib <- suppressMessages(suppressWarnings( readr::read_csv(betaCSV) ))
    if (!is.null(cgn)) beta_tib <- beta_tib %>% dplyr::inner_join(cgn, by=idKey)
    if (!is.null(ss)) beta_tib <- beta_tib %>% dplyr::select(!!idKey_sym, ss$Sentrix_Name)
    if (verbose>=vt+4) print(beta_tib)
    
    pval_tib <- suppressMessages(suppressWarnings( readr::read_csv(pvalCSV) ))
    if (!is.null(cgn)) pval_tib <- pval_tib %>% dplyr::inner_join(cgn, by=idKey)
    if (!is.null(ss)) pval_tib <- pval_tib %>% dplyr::select(!!idKey_sym, ss$Sentrix_Name)
    if (verbose>=vt+4) print(pval_tib)
    
    if (addPval) {
      if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Will combine pvals with beta values...{RET}"))
      
      beta_col <- paste(names(beta_tib),betaName, sep=del)
      beta_col[1] <- idKey
      beta_tib <- beta_tib %>% purrr::set_names(beta_col)

      pval_col <- paste(names(pval_tib),pvalName, sep=del)
      pval_col[1] <- idKey
      pval_tib <- pval_tib %>% purrr::set_names(pval_col)
    }
    
    # Convert to matrix::
    if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Converting to matricies...{RET}"))
    beta_mat <- beta_tib %>% tibble::column_to_rownames(var=idKey) %>% as.data.frame() %>% as.matrix()
    pval_mat <- pval_tib %>% tibble::column_to_rownames(var=idKey) %>% as.data.frame() %>% as.matrix()
    
    pval_idx <- c()
    if (!is.null(minPval)) pval_idx <- which(pval_mat > minPval)
    fail_cnt <- length(pval_idx)
    
    if (fail_cnt>0) beta_mat[pval_idx] <- NA
    nana_cnt <- which(is.na(beta_mat)) %>% length()
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} minPval={minPval}, failCnt={fail_cnt}, naCnt={nana_cnt}.{RET}"))
    
    if (addPval) {
      beta_mat <- cbind(beta_mat, pval_mat[ match(rownames(beta_mat), rownames(pval_mat)), ])
      beta_col <- colnames(beta_mat) %>% sort()
      beta_mat <- beta_mat[,beta_col]
    }
    
    if (is.null(mat)) {
      mat <- beta_mat
    } else {
      mat <- cbind(mat, beta_mat[ match(rownames(mat), rownames(beta_mat)), ] )
    }
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  mat
}

getCallsMatrixFiles = function(betaKey,pvalKey,pvalMin, dirs, cgn=NULL, classes=NULL,
                               class_var, class_idx, pval_name, pval_perc,
                               beta_rds, pval_rds=NULL, ss_csv, mask_csv,
                               clean=FALSE, addPval=FALSE, 
                               sentrix_name="Sentrix_Name", idKey="Probe_ID", 
                               betaName='beta', pvalName='pval', del='.', exp_name='Experiment_Key',
                               sam_suffix="_AutoSampleSheet.csv.gz$",dat_suffix="_MergedDataFiles.tib.csv.gz",
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getCallsMatrixFiles'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))

  ret <- NULL
  beg_txt <- NULL
  end_txt <- NULL
  
  beta_mat <- NULL
  labs_tib <- NULL
  file_tib <- tibble::tibble( Type=character(), File_Name=character(), Dir=character() )
  
  if (opt$verbose>=vt && clean) cat(glue::glue("[{funcTag}]:{tabsStr} Will clean all matrix files...{RET}"))
  
  if (addPval && is.null(pval_rds)) {
    cat(glue::glue("[{funcTag}]:{tabsStr} WARNING: Calling addPval without defined pval_rds output!!!{RET}"))
    addPval <- FALSE
  }
  
  sort_ss_tibs <- NULL
  masked_idx_tibs <- NULL
  impt_mats <- NULL
  pval_mats <- NULL
  
  stime <- system.time({
    
    beg_txt <- paste(beta_rds,'begTime.txt', sep='.')
    end_txt <- paste(beta_rds,'endTime.txt', sep='.')
    nan_csv <- paste(beta_rds,'rm-cgn.class.csv.gz', sep='.')

    if (!clean && 
        file.exists(beg_txt) && file.exists(end_txt) &&
        file.exists(mask_csv) && file.exists(ss_csv) && file.exists(beta_rds) && 
        file.mtime(beg_txt)  < file.mtime(ss_csv) &&
        file.mtime(ss_csv)   < file.mtime(mask_csv) &&
        file.mtime(mask_csv) < file.mtime(beta_rds) &&
        file.mtime(beta_rds) < file.mtime(end_txt) ) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} All Beta Matrix files exist in correct order.{RET}"))
      
    } else {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Will build fresh version of matrix files...{RET}"))
      
      if (file.exists(beg_txt) && file.exists(end_txt) &&
          file.exists(mask_csv) && file.exists(ss_csv) && file.exists(beta_rds)) {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} All files do not exist...{RET}"))
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}  beg_txt={beg_txt}.{RET}"))
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} mask_csv={mask_csv}.{RET}"))
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}   ss_csv={ss_csv}.{RET}"))
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} beta_rds={beta_rds}.{RET}"))
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} end_txt={end_txt}.{RET}"))

        if (verbose>=vt+4) {
          if (file.mtime(beg_txt)  < file.mtime(ss_csv))   cat(glue::glue("[{funcTag}]:{tabsStr} beg_txt < ss_csv{RET}"))
          if (file.mtime(ss_csv)   < file.mtime(mask_csv)) cat(glue::glue("[{funcTag}]:{tabsStr} ss_csv < mask_csv.{RET}"))
          if (file.mtime(mask_csv) < file.mtime(beta_rds)) cat(glue::glue("[{funcTag}]:{tabsStr} mask_csv < beta_rds.{RET}"))
          if (file.mtime(beta_rds) < file.mtime(end_txt))  cat(glue::glue("[{funcTag}]:{tabsStr} beta_rds < end_txt{RET}"))
        }
      }
      
      if (file.exists(beg_txt)) {
        unlink(beg_txt)
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Cleaned  beg_txt={beg_txt}.{RET}"))
      }
      if (file.exists(nan_csv)) {
        unlink(nan_csv)
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Cleaned  nan_csv={nan_csv}.{RET}"))
      }
      if (file.exists(mask_csv)) {
        unlink(mask_csv)
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Cleaned mask_csv={mask_csv}.{RET}"))
      }
      if (file.exists(ss_csv)) {
        unlink(ss_csv)
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Cleaned   ss_csv={ss_csv}.{RET}"))
      }
      if (file.exists(beta_rds)) {
        unlink(beta_rds)
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Cleaned beta_rds={beta_rds}.{RET}"))
      }
      if (file.exists(pval_rds)) {
        unlink(pval_rds)
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Cleaned pval_rds={pval_rds}.{RET}"))
      }
      if (file.exists(end_txt)) {
        unlink(end_txt)
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Cleaned  end_txt={end_txt}.{RET}"))
      }
      if (verbose>=vt+4) cat(glue::glue("{RET}{RET}{RET}"))
      
      if (base::typeof(exp_name)=="character") exp_name <- rlang::sym(exp_name)
      if (base::typeof(class_var)=="character") class_var <- rlang::sym(class_var)
      if (base::typeof(class_idx)=="character") class_idx <- rlang::sym(class_idx)
      if (base::typeof(sentrix_name)=="character") sentrix_name <- rlang::sym(sentrix_name)
      stopifnot(base::typeof(exp_name)=="symbol")
      stopifnot(base::typeof(class_var)=="symbol")
      stopifnot(base::typeof(class_idx)=="symbol")
      stopifnot(base::typeof(sentrix_name)=="symbol")
      # pval_name <- rlang::sym(pval_name)
      # pval_perc <- rlang::sym(pval_perc)
      
      trainClass_vec <- NULL
      if (!is.null(classes)) trainClass_vec <- classes %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
      if (!is.null(classes)) {
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{TAB} classes={classes}; trainClass_vec={RET}") )
          print(trainClass_vec)
          cat(glue::glue("{RET}{RET}"))
        }
      }
      
      for (curDir in dirs) {
        
        # Find Sample Sheet::
        cur_ss_csv <- findFileByPattern(dir=curDir, patter=sam_suffix, max=1, recursive=FALSE, verbose=verbose)
        cur_fn_csv <- cur_ss_csv %>% stringr::str_replace(sam_suffix,dat_suffix)
        
        base_dir <- base::dirname(cur_ss_csv)
        base_key <- base::basename(base_dir)
        stopifnot(file.exists(cur_ss_csv), file.exists(cur_fn_csv))
        
        # Load and filter::
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Loading Sample CSV; cur_ss_csv={cur_ss_csv}...{RET}") )
        cur_ss_tib <- suppressMessages(suppressWarnings( readr::read_csv(cur_ss_csv) )) %>% dplyr::mutate(!!exp_name:=base_key)
        if (!is.null(pval_name) && !is.null(pval_perc)) cur_ss_tib <- cur_ss_tib %>% dplyr::filter(!!pval_name > !!pval_perc)
        if (verbose>=vt+4) print(cur_ss_tib)

        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Filtering Sample CSV; class_var={class_var}...{RET}") )
        if (!is.null(classes) && length(trainClass_vec)>0)
          cur_ss_tib <- cur_ss_tib %>% dplyr::filter(!!class_var %in% trainClass_vec)
        
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{TAB} Filtered Sample Sheet; cur_ss_tib=...{RET}") )
        if (verbose>=vt+4) print(cur_ss_tib)
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Loading call_table={cur_fn_csv}...{RET}") )
        calls_path_tib <- suppressMessages(suppressWarnings( readr::read_csv(cur_fn_csv) ))
        if (verbose>=vt+4) print(calls_path_tib)
        
        betas_path_tib <- calls_path_tib %>% dplyr::filter(Method %in% c(betaKey))
        pvals_path_tib <- calls_path_tib %>% dplyr::filter(Method %in% c(pvalKey))
        
        stopifnot(base::nrow(betas_path_tib)==1)
        stopifnot(base::nrow(pvals_path_tib)==1)
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #            Load Beta/Pval And Merge into Previous Matrix::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Load Beta/Pval And Merge into Previous Matrix...{RET}") )
        
        beta_csv <- file.path(base_dir, base::basename(betas_path_tib$Full_Path[1]) )
        pval_csv <- file.path(base_dir, base::basename(pvals_path_tib$Full_Path[1]) )
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} beta_csv={beta_csv}.{RET}") )
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pval_csv={pval_csv}.{RET}") )
        beta_mat <- loadCallsMatrix(betaCSV=beta_csv, pvalCSV=pval_csv, minPval=pvalMin, mat=beta_mat, 
                                    cgn=cgn, ss=cur_ss_tib,
                                    idKey=idKey, addPval=addPval, betaName=betaName, pvalName=pvalName, del='.',
                                    verbose=verbose, vt=vt,tc=1, tt=tt)
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} beta_mat=...{RET}") )
        if (verbose>=vt+4) beta_mat %>% head(n=1) %>% print()
        
        labs_tib <- labs_tib %>% dplyr::bind_rows(cur_ss_tib)
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Finished loading; outName={outName}.{RET}") )
        if (verbose>=vt+4) print(labs_tib)
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                      Sort Sample Sheet by class_var::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Sorting Sample Sheet by class_var={class_var}.{RET}") )
        sort_ss_tib <- labs_tib %>% dplyr::arrange(!!class_var) %>% 
          dplyr::mutate(!!class_var := as.factor(!!class_var),
                        !!class_idx := as.integer(!!class_var)-1) %>% 
          dplyr::mutate(!!class_idx := as.integer(!!class_idx) ) %>% 
          dplyr::select(!!sentrix_name,!!class_var,!!class_idx,!!exp_name, everything()) 
        # %>% dplyr::select(!!sentrix_name, !!class_var, !!class_idx)

        if (verbose>=vt+4) print(sort_ss_tib)
        
        # Bit of tricks going on below; order of operations matters::
        pval_mat <- NULL
        beta_col <- dplyr::pull(sort_ss_tib, !!sentrix_name)
        if (addPval) {
          pval_col <- paste(beta_col,pvalName,sep=del)
          beta_col <- paste(beta_col,betaName,sep=del)
          pval_mat <- beta_mat[, pval_col];
        }
        mask_mat <- beta_mat[ , beta_col ]
        
        # Strip any suffix 
        if (addPval) {
          beta_col <- mask_mat %>% colnames() %>% stringr::str_remove(paste0(del,betaName))
          colnames(mask_mat) <- beta_col
          
          pval_col <- pval_mat %>% colnames() %>% stringr::str_remove(paste0(del,pvalName))
          colnames(pval_mat) <- pval_col
        }
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Masked Matrix={RET}") )
        if (verbose>=vt+4) beta_mat %>% head(n=1) %>% print()
        if (verbose>=vt+4) print(beta_col)
        if (verbose>=vt+4) mask_mat %>% head(n=1) %>% print()
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                   Build Raw and Imputed Sorted Matricies::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building Raw and Imputed Sorted Matricies; sentrix_name={sentrix_name}...{RET}") )
        
        #
        # Imputation needs to be done on a class basis::
        #
        sort_ss_names <- sort_ss_tib %>% dplyr::distinct(!!class_var) %>% dplyr::pull(!!class_var) %>% as.vector()
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} sort_ss_names={RET}") )
        if (verbose>=vt+4) print(sort_ss_names)
        
        impt_mat <- NULL
        for (sName in sort_ss_names) {
          if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Imputing Class={sName}...{RET}") )
          
          cur_ss_tib <- sort_ss_tib %>% dplyr::filter(!!class_var == sName)
          if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_ss_tib::{RET}") )
          if (verbose>=vt+4) print(cur_ss_tib)
          
          cur_sort_vec <- dplyr::pull(cur_ss_tib, !!sentrix_name)
          if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} sentrix_name={sentrix_name}; cur_sort_vec={RET}") )

          cur_mask_mat <- mask_mat[ , cur_sort_vec ]
          if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_mask_mat::{RET}") )
          if (verbose>=vt+4) cur_mask_mat %>% head(n=1) %>% print()
          
          cur_impt_mat <- impute_matrix_mean(mat=cur_mask_mat)
          if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_impt_mat(dim-1)::{RET}") )
          if (verbose>=vt+4) cur_impt_mat %>% dim() %>% print()

          if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_impt_mat::{RET}") )
          if (verbose>=vt+4) cur_impt_mat %>% head(n=3) %>% print()
          
          if (is.null(impt_mat)) {
            impt_mat <- cur_impt_mat
          } else {
            impt_mat <- cbind(impt_mat, cur_impt_mat[match(rownames(impt_mat), rownames(cur_impt_mat)), ] )
          }
          if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_impt_mat(dim-2)::{RET}") )
          if (verbose>=vt+4) impt_mat %>% dim() %>% print()
        }
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Finished imputation across all classes.{RET}{RET}") )
        
        # Need to remove any probes that could not be imputed!!
        #  NOTE:: Found these in cancer samples training where the probe is likely deleted and never found!
        #
        rm_row_cnt <- 0
        rm_cgn_cnt <- 0
        rm_impute_tib <- NULL
        
        rm_cgn_tib <- getMaskedTib(impt_mat,verbose=verbose,vt=vt,tc=tc,tt=tt)
        if (verbose>=vt+4) rm_cgn_tib %>% head(n=3) %>% print()
        rm_tot_cnt <- rm_cgn_tib %>% base::nrow()
        if (rm_tot_cnt>0) {
          rm_row_cnt <- rm_cgn_tib %>% dplyr::distinct(row) %>% base::nrow()
          rm_cgn_cnt <- rm_cgn_tib %>% dplyr::distinct(Probe_ID) %>% base::nrow()
        }
        if (verbose>=vt+4) 
          cat(glue::glue("[{funcTag}]:{tabsStr} Recalculated Masked(NA) Probes; rm_tot_cnt={rm_tot_cnt}, ",
                         "rm_row_cnt={rm_row_cnt}, rm_cgn_cnt={rm_cgn_cnt}.{RET}") )
        
        if (rm_cgn_cnt!=0) {
          if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Forced to remove rows={rm_row_cnt}, cgns={rm_cgn_cnt}, total={rm_tot_cnt} ",
                                          "from training data due to too many NA's to impute.{RET}") )
          rm_row_vec    <- rm_cgn_tib %>% dplyr::distinct(row) %>% dplyr::pull(row) %>% as.vector() %>% unique()
          # rm_row_vec %>% head() %>% print()
          
          if (FALSE) {
            rm_impute_tib <- impt_mat[  rm_row_vec, ] %>% as.data.frame() %>% 
              tibble::rownames_to_column(var="Probe_ID") %>% tibble::tibble()
            # rm_impute_tib %>% head() %>% print()
          }
          
          impt_mat <- impt_mat[ -rm_row_vec, ]
          if (verbose>=vt+4) impt_mat %>% head() %>% print()
          
          mask_mat <- mask_mat[ -rm_row_vec, ]
          if (verbose>=vt+4)  mask_mat %>% head() %>% print()
        } else {
          if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Not Forced to remove any probes.{RET}") )
          if (verbose>=vt+4) impt_mat %>% head(n=1) %>% print()
          if (verbose>=vt+4) impt_mat %>% dim() %>% print()
        }
        
        # Calculate Summary Stats for QC::
        #
        masked_idx <- which(is.na(mask_mat))
        impute_idx <- which(is.na(impt_mat))
        masked_cnt <- length(masked_idx)
        impute_cnt <- length(impute_idx)
        total_cnt  <- length(mask_mat)
        masked_per <- round(100*masked_cnt/total_cnt, 3)
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Imputation Results: masked_cnt={masked_cnt}, ",
                                        "post-impute_masked-cnt={impute_cnt}, masked-perc={masked_per}.{RET}") )
        if (impute_cnt!=0)
          stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed Imputation: impute_cnt={impute_cnt}, masked_cnt={masked_cnt}{RET}{RET}"))
        
        # masked_idx_tib <- tibble::tibble(Masked_Index=masked_idx)
        masked_idx_tib <- mask_mat %>% getMaskedTib(verbose=verbose,vt=vt,tc=tc,tt=tt)
        impute_idx_tib <- tibble::tibble(Impute_Index=impute_idx)
        masked_idx_cnt <- masked_idx_tib %>% base::nrow()
        impute_idx_cnt <- impute_idx_tib %>% base::nrow()
        
        stopifnot(masked_cnt==masked_idx_cnt)
        stopifnot(impute_cnt==impute_idx_cnt)
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #               Write Outputs:: Sorted Matrix AND Sample Sheet
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # sort_ss_tibs    <- dplyr::bind_rows(sort_ss_tibs,sort_ss_tib)
      # masked_idx_tibs <- dplyr::bind_rows(masked_idx_tibs, masked_idx_tib)
      # impt_mats <- dplyr::bind_rows(impt_mats, impt_mat)
      # pval_mats <- dplyr::bind_rows(pval_mats, pval_mat)
      
      system(glue::glue("touch {beg_txt}"))

      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing Sorted Sample Sheet={ss_csv}.{RET}") )
      readr::write_csv(sort_ss_tib,ss_csv)
      Sys.sleep(1)
      
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing Beta Masked Vector={mask_csv}.{RET}") )
      readr::write_csv(masked_idx_tib,mask_csv)
      Sys.sleep(1)
      
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing Beta Imputed Matrix={beta_rds}.{RET}") )
      readr::write_rds(impt_mat,beta_rds, compress="gz")
      Sys.sleep(1)
      
      if (addPval) {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing Pval Matrix={beta_rds}.{RET}") )
        readr::write_rds(pval_mat,pval_rds, compress="gz")
        Sys.sleep(1)
      }
      system(glue::glue("touch {end_txt}"))
      
    }
    file_tib <- file_tib %>% tibble::add_row( Type="Mask", File_Name=base::basename(mask_csv), Dir=base::dirname(mask_csv) )
    file_tib <- file_tib %>% tibble::add_row( Type="SampleSheet", File_Name=base::basename(ss_csv), Dir=base::dirname(ss_csv) )
    file_tib <- file_tib %>% tibble::add_row( Type="Beta", File_Name=base::basename(beta_rds), Dir=base::dirname(beta_rds) )
    if (addPval) file_tib <- file_tib %>% tibble::add_row( Type="Pval", File_Name=base::basename(pval_rds), Dir=base::dirname(pval_rds) )
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  file_tib
}


# End of file
