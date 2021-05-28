
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Sample Sheet Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

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
#                        Sample Sheet I/O Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadAutoSampleSheets = function(dir, platform=NULL, manifest=NULL, workflow=NULL,
                                suffix='AutoSampleSheet.csv.gz',
                                addSampleName=FALSE, addPathsCall=FALSE, addPathsSset=FALSE,
                                flagDetectPval=FALSE, flagSampleDetect=FALSE, flagRefMatch=FALSE,
                                pvalDetectMinKey=NULL, pvalDetectMinVal=0,
                                
                                dbMin=90, r2Min=0.9,clean_gta=FALSE,
                                dbKey='AutoSample_dB_Key', r2Key='AutoSample_R2_Key',
                                dbVal='AutoSample_dB_Val', r2Val='AutoSample_R2_Val',
                                verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadAutoSampleSheets'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, dir={dir}.{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  pattern <- NULL
  if (!is.null(platform)) pattern <- paste(pattern,platform, sep='_')
  if (!is.null(manifest)) pattern <- paste(pattern,manifest, sep='_')
  pattern <- paste(pattern,suffix, sep='_')
  
  auto_ss_list <- list.files(dir, pattern=pattern, recursive=TRUE, full.names=TRUE)
  auto_ss_cnts <- auto_ss_list %>% length()
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting, SampleSheetCount={auto_ss_cnts}, pattern={pattern} dir={dir}.{RET}"))
  
  if (is.null(auto_ss_cnts) || auto_ss_cnts==0) {
    if (verbose>=vt+1)
      cat(glue::glue("{RET}[{funcTag}]:{tabsStr} Warning: Failed to find auto sample sheets!{RET}{RET}"))
    return(ret_tib)
  }

  pvalDetectMinKey <- NULL
  stime <- base::system.time({

    if (!is.null(pvalDetectMinKey)) pvalDetectMinKey <- pvalDetectMinKey %>% rlang::sym()
    # pvalDetectMinVal <- pvalDetectMinVal %>% rlang::sym()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Load Samples::
    ret_tib <- suppressMessages(suppressWarnings(lapply(auto_ss_list, readr::read_csv) )) %>% 
      dplyr::bind_rows()
    
    if (clean_gta) {
      ret_tib <- ret_tib %>% dplyr::mutate(
        Sentrix_Name=stringr::str_remove(Sentrix_Name,'_[0-9]+$'),
        Sentrix_Poscode=stringr::str_remove(Sentrix_Poscode,'_[0-9]+$')
      )
    }
    
    auto_ss_tlen <- base::nrow(ret_tib)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(bPool)={auto_ss_tlen}{RET}"))
    # print(ret_tib)
    
    if (addSampleName || flagDetectPval || flagSampleDetect || flagRefMatch) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Adding flags...{RET}"))
      if (verbose>=vt+4) print(ret_tib)
      
      dbKey <- dbKey %>% rlang::sym()
      dbVal <- dbVal %>% rlang::sym()
      r2Key <- r2Key %>% rlang::sym()
      r2Val <- r2Val %>% rlang::sym()
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #  Add General Sample Name Field::
      if (addSampleName) {
        ret_tib <- ret_tib %>% dplyr::mutate(Auto_Sample_Name=!!dbKey) %>%
          dplyr::select(Auto_Sample_Name, dplyr::everything())
      }
      if (verbose>=vt+4) print(ret_tib)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #  Flag Probe Detected (pval)::
      if (flagDetectPval && !is.null(pvalDetectMinKey) && !is.null(pvalDetectMinVal)) {
        fail_tag <- paste0("Failed<",pvalDetectMinVal)
        pass_tag <- paste0("Passed>=",pvalDetectMinVal)
        
        ret_tib <- ret_tib %>% 
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
        
        ret_tib <- ret_tib %>% dplyr::mutate(detectedSample=case_when(
          !!r2Val < r2Min ~ fail_r2,
          !!dbVal < dbMin ~ fail_db,
          !!dbKey != !!r2Key ~ fail_mt,
          TRUE ~ pass_tag)
        )
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #  Flag Auto Detection with poor matching (pval)::
      if (flagRefMatch) {
        ret_tib <- ret_tib %>% 
          dplyr::filter() %>% 
          dplyr::filter(!!dbVal >= dbMin) %>% 
          dplyr::filter(!!r2Val >= r2Min)
        auto_ss_flen <- base::nrow(ret_tib)
        
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(filt)={auto_ss_flen}{RET}"))
        
        # Remove Unidentifiable
        # if (rmOdd) ret_tib <- ret_tib %>% dplyr::filter(Bead_Pool!='Odd')
      }
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Add Paths::
    if (addPathsCall) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Adding Calls Paths...{RET}"))
      
      ret_tib <- 
        addPathsToSampleSheet(ss=ret_tib, dir=dir, 
                              platform=platform, manifest=manifest, workflow=workflow,
                              del='_',field='Calls_Path', suffix='call.dat.csv.gz$',
                              verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      auto_ss_tlen <- base::nrow(ret_tib)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(calls)={auto_ss_tlen}.{RET}{RET}"))
    }
    
    if (addPathsSset) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Adding SSETs Paths...{RET}"))
      
      ret_tib <- 
        addPathsToSampleSheet(ss=ret_tib, dir=dir, 
                              platform=platform, manifest=manifest, workflow=workflow,
                              del='_',field='Ssets_Path', suffix='sigs.dat.csv.gz$', 
                              verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      auto_ss_tlen <- base::nrow(ret_tib)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(ssets)={auto_ss_tlen}.{RET}{RET}"))
    }
    
    # dat <- ret_tib
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  ret_tib
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

addPathsToSampleSheet = function(ss, dir, 
                                 platform=NULL, manifest=NULL, workflow=NULL,
                                 field, suffix, del='.',
                                 verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'addPathsToSampleSheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; suffix={suffix}.{RET}"))
  
  pattern <- NULL
  if (!is.null(platform)) pattern <- paste(pattern,platform, sep=del)
  if (!is.null(manifest)) pattern <- paste(pattern,manifest, sep=del)
  if (!is.null(workflow)) pattern <- paste(pattern,workflow, sep=del)
  pattern <- paste(pattern,suffix, sep='.')
  
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
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Found; ss_path_cnts={ss_path_cnts}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} ss_path_tibs={RET}"))
  if (verbose>=vt+4) print(ss_path_tibs)
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} ss={RET}"))
  if (verbose>=vt+4) print(ss)
  
  tib <- ss_path_tibs %>%
    dplyr::distinct(Sentrix_Name, .keep_all=TRUE) %>% 
    dplyr::inner_join(ss, by="Sentrix_Name")
  
  tib_len <- tib %>% base::nrow()
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; tib_len={tib_len}.{RET}{RET}"))
  
  tib
}

addBeadPoolToSampleSheet = function(ss, field, 
                                    verbose=0,vt=4,tc=1,tt=NULL) {
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
                            chipName='Sentrix_Name', pathName='Calls_Path', 
                            joinNameA="Probe_ID", joinNameB="Probe_Design",
                            verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'mergeCallsFromSS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; chipName={chipName}, pathName={pathName}, joinName={joinNameA}.{RET}"))
  
  file_list <- base::list()
  data_list <- base::list()
  
  ret_cnt <- 0
  ret_tib <- NULL
  
  stime <- base::system.time({
    
    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
    
    ss_cnt <- ss %>% base::nrow()
    for (ssIdx in c(1:ss_cnt)) {
      sentrix_name <- ss[[chipName]][ssIdx]
      call_path    <- ss[[pathName]][ssIdx]
      
      ss_perc = round(100*ssIdx/ss_cnt, 3)
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading {pathName}; ",
                       "idx={ssIdx}/{ss_cnt}={ss_perc}, ",
                       "sentrix_name={sentrix_name}; path={call_path}.{RET}"))
      
      call_tib  <- suppressMessages(suppressWarnings( readr::read_csv(call_path) ))
      
      if (pathName=='Calls_Path') {
        call_cols <- call_tib %>% dplyr::select(-c(!!joinNameA)) %>% names()
        call_lens <- length(call_cols)
      } else {
        call_cols <- call_tib %>% dplyr::select(-c(!!joinNameA,!!joinNameB)) %>% names()
        call_lens <- length(call_cols)
      }
      # print(call_cols)
      # print(call_tib)
      
      for (cIdx in (1:call_lens)) {
        cur_col <- call_cols[cIdx]
        # cur_tib <- dplyr::select(call_tib, 1,cur_col)
        
        if (pathName=='Calls_Path') {
          cur_tib <- dplyr::select(call_tib, !!joinNameA, dplyr::all_of(cur_col) )
          colnames(cur_tib)[2] <- sentrix_name
          if (is.null(data_list[[cur_col]])) {
            data_list[[cur_col]] <- cur_tib
          } else {
            data_list[[cur_col]] <- dplyr::full_join(data_list[[cur_col]],cur_tib, by=joinNameA)
          }
        } else {
          cur_tib <- dplyr::select(call_tib, !!joinNameA,!!joinNameB, dplyr::all_of(cur_col) )
          colnames(cur_tib)[3] <- sentrix_name
          
          if (is.null(data_list[[cur_col]])) {
            data_list[[cur_col]] <- cur_tib
          } else {
            data_list[[cur_col]] <- dplyr::full_join(data_list[[cur_col]],cur_tib, by=c(joinNameA,joinNameB) )
          }
        }
      }
      
      if (!is.null(max) && max>0 && ssIdx>=max) break
    }
    
    for (name in names(data_list)) {
      out_name <- paste(outName,name, sep='_')
      file_list[[name]] <- file.path(outDir, paste(out_name,'dat.csv.gz', sep='.') )
      
      safe_write(data_list[[name]],"csv",file_list[[name]],
                 funcTag=funcTag, verbose=verbose,tt=tt)
      
      # if (verbose>=vt)
      #   cat(glue::glue("[{funcTag}]: Writing {pathName} File; name={name}, out_name={out_name}; CSV={file_list[[name]]}...{RET}"))
      # readr::write_csv(data_list[[name]], file_list[[name]])
    }
    
    # Build return tibble of mreged/split calls files::
    ret_tib <- file_list %>% 
      dplyr::bind_rows() %>% 
      tidyr::gather(Method, Full_Path)
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
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
  
  stime <- base::system.time({
    
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
      if (verbose>=vt+1) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Will combine pvals with beta values...{RET}"))
      
      beta_col <- paste(names(beta_tib),betaName, sep=del)
      beta_col[1] <- idKey
      beta_tib <- beta_tib %>% purrr::set_names(beta_col)
      
      pval_col <- paste(names(pval_tib),pvalName, sep=del)
      pval_col[1] <- idKey
      pval_tib <- pval_tib %>% purrr::set_names(pval_col)
    }
    
    # Convert to matrix::
    if (verbose>=vt+1)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Converting to matricies...{RET}"))
    beta_mat <- beta_tib %>% tibble::column_to_rownames(var=idKey) %>% as.data.frame() %>% as.matrix()
    pval_mat <- pval_tib %>% tibble::column_to_rownames(var=idKey) %>% as.data.frame() %>% as.matrix()
    
    pval_idx <- c()
    if (!is.null(minPval)) pval_idx <- which(pval_mat > minPval)
    fail_cnt <- length(pval_idx)
    
    if (fail_cnt>0) beta_mat[pval_idx] <- NA
    nana_cnt <- which(is.na(beta_mat)) %>% length()
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} minPval={minPval}, failCnt={fail_cnt}, naCnt={nana_cnt}.{RET}"))
    
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
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} dim(mat)={RET}"))
      mat %>% dim() %>% print()
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
  
  if (verbose>=vt && clean) cat(glue::glue("[{funcTag}]:{tabsStr} Will clean all matrix files...{RET}"))
  
  if (addPval && is.null(pval_rds)) {
    cat(glue::glue("[{funcTag}]:{tabsStr} WARNING: Calling addPval without defined pval_rds output!!!{RET}"))
    addPval <- FALSE
  }
  
  sort_ss_tibs <- NULL
  masked_idx_tibs <- NULL
  impt_mats <- NULL
  pval_mats <- NULL
  
  stime <- base::system.time({
    
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
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{TAB} Loading Sample CSV; cur_ss_csv={cur_ss_csv}...{RET}") )
        
        cur_ss_tib <- suppressMessages(suppressWarnings( readr::read_csv(cur_ss_csv) )) %>% 
          dplyr::mutate(!!exp_name:=base_key)
        if (!is.null(pval_name) && !is.null(pval_perc)) cur_ss_tib <- cur_ss_tib %>% dplyr::filter(!!pval_name > !!pval_perc)
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{TAB} cur_ss_tib={RET}") )
          print(cur_ss_tib)
        }
        
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{TAB} Filtering Sample CSV; class_var={class_var}...{RET}") )
        if (!is.null(classes) && length(trainClass_vec)>0)
          cur_ss_tib <- cur_ss_tib %>% dplyr::filter(!!class_var %in% trainClass_vec)
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{TAB} Filtered Sample Sheet; cur_ss_tib={RET}") )
          print(cur_ss_tib)
        }
        
        if (verbose>=vt)
          cat(glue::glue("[{funcTag}]:{TAB} Loading call_table={cur_fn_csv}...{RET}") )
        calls_path_tib <- suppressMessages(suppressWarnings( readr::read_csv(cur_fn_csv) ))
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{TAB} calls_path_tib={RET}") )
          print(calls_path_tib)
        }
        
        betas_path_tib <- calls_path_tib %>% dplyr::filter(Method %in% c(betaKey))
        pvals_path_tib <- calls_path_tib %>% dplyr::filter(Method %in% c(pvalKey))
        
        stopifnot(base::nrow(betas_path_tib)==1)
        stopifnot(base::nrow(pvals_path_tib)==1)
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #            Load Beta/Pval And Merge into Previous Matrix::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        if (verbose>=vt)
          cat(glue::glue("[{funcTag}]:{tabsStr} Load Beta/Pval And Merge into Previous Matrix...{RET}") )
        
        beta_csv <- file.path(base_dir, base::basename(betas_path_tib$Full_Path[1]) )
        pval_csv <- file.path(base_dir, base::basename(pvals_path_tib$Full_Path[1]) )
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} beta_csv={beta_csv}.{RET}") )
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} pval_csv={pval_csv}.{RET}") )
        beta_mat <- loadCallsMatrix(betaCSV=beta_csv, pvalCSV=pval_csv, minPval=pvalMin, mat=beta_mat, 
                                    cgn=cgn, ss=cur_ss_tib,
                                    idKey=idKey, addPval=addPval, betaName=betaName, pvalName=pvalName, del='.',
                                    verbose=verbose, vt=vt,tc=1, tt=tt)
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} beta_mat={RET}") )
          beta_mat %>% head(n=1) %>% print()
        }
        
        labs_tib <- labs_tib %>% dplyr::bind_rows(cur_ss_tib)
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Finished loading; outName={outName}.{RET}") )
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} labs_tib={RET}") )
          print(labs_tib)
        }
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                      Sort Sample Sheet by class_var::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{tabsStr} Sorting Sample Sheet by class_var={class_var}.{RET}") )
        
        sort_ss_tib <- labs_tib %>% dplyr::arrange(!!class_var) %>% 
          dplyr::mutate(!!class_var := as.factor(!!class_var),
                        !!class_idx := as.integer(!!class_var)-1) %>% 
          dplyr::mutate(!!class_idx := as.integer(!!class_idx) ) %>% 
          dplyr::select(!!sentrix_name,!!class_var,!!class_idx,!!exp_name, everything()) 
        # %>% dplyr::select(!!sentrix_name, !!class_var, !!class_idx)
        
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} sort_ss_tib={RET}") )
          print(sort_ss_tib)
        }
        
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
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} Masked Matrix beta_mat={RET}") )
          beta_mat %>% head(n=1) %>% print()
          cat(glue::glue("[{funcTag}]:{tabsStr} Masked Matrix beta_col={RET}") )
          print(beta_col)
          cat(glue::glue("[{funcTag}]:{tabsStr} Masked Matrix mask_mat={RET}") )
          mask_mat %>% head(n=1) %>% print()
        }
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                   Build Raw and Imputed Sorted Matricies::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building Raw and Imputed Sorted Matricies; sentrix_name={sentrix_name}...{RET}") )
        
        #
        # Imputation needs to be done on a class basis::
        #
        sort_ss_names <- sort_ss_tib %>% dplyr::distinct(!!class_var) %>% dplyr::pull(!!class_var) %>% as.vector()
        if (verbose>=vt+4) {
          cat(glue::glue("[{funcTag}]:{tabsStr} sort_ss_names={RET}") )
          print(sort_ss_names)
        }
        
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Sample Sheet Documentation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getSsheetDataTab = function(tib,
                            minOobPval,minOobPerc,
                            minNegPval,minNegPerc,
                            minDb,
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getSsheetDataTab'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ret_tib <- suppressMessages( suppressWarnings(
      dplyr::inner_join(
        tib %>% tidyr::gather(Variable, Value),
        tibble::enframe( unlist( sapply(tib, class) ), name="Variable", value="Data_Type" ),
        by="Variable"
      ) ) %>%
        dplyr::mutate(
          Variable=stringr::str_replace_all(
            stringr::str_squish((stringr::str_replace_all(Variable, regex("\\W+"), " ")) ), " ", "_") )
    )

    if (verbose>=vt+3) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
      ret_tib %>% print()
    }
    
    dec_tib <- NULL
    dec_tib <- suppressMessages( suppressWarnings(
      getSsheetDescTib(tib, 
                       minOobPval=minOobPval,minOobPerc=minOobPerc,
                       minNegPval=minNegPval,minNegPerc=minNegPerc,
                       minDb=minDb,
                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
        gather(Variable, Description)
    ))
    if (verbose>=vt+3) {
      cat(glue::glue("[{funcTag}]:{tabsStr} dec_tib={RET}"))
      dec_tib %>% print(n=base::nrow(dec_tib))
    }
    
    ret_tib <- ret_tib %>% dplyr::left_join(dec_tib, by="Variable")
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
      ret_tib %>% print(n=base::nrow(ret_tib))
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

getSsheetDescTib = function(tib,
                            minOobPval,minOobPerc,
                            minNegPval,minNegPerc,
                            minDb,
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getSsheetDescTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) {
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}minOobPval={minOobPval}, minOobPerc={minOobPerc}...{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}minNegPval={minNegPval}, minNegPerc={minNegPerc}...{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # Determine Number of Workflows::
    #
    max_idx <- tib %>% 
      dplyr::select(dplyr::starts_with('Method_Idx')) %>% 
      tidyr::gather() %>% 
      dplyr::arrange(-value) %>% 
      head(n=1) %>% dplyr::pull(value)
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Max Index={max_idx}.{RET}"))
    
    # Build Description Sheet
    #
    ret_tib <- suppressMessages( suppressWarnings( 
      getSsheetCoreAnnoTib(minOobPval,minOobPerc,
                           minNegPval,minNegPerc,
                           minDb=minDb,
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) ) )
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
      ret_tib %>% print(n=base::nrow(ret_tib))
    }
    
    for (idx in c(0:max_idx)) {
      if (verbose>=vt+2)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Current Index={idx}.{RET}"))
      
      ret_tib <- ret_tib %>% 
        dplyr::bind_cols(suppressMessages( suppressWarnings(
          getSsheetIndexAnnoTib(idx=idx, 
                                minOobPval=minOobPval,minOobPerc=minOobPerc,
                                minNegPval=minNegPval,minNegPerc=minNegPerc,
                                minDb=minDb,
                                verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)) ) )
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
        ret_tib %>% print(n=base::nrow(ret_tib))
      }
    }
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
      ret_tib %>% print(n=base::nrow(ret_tib))
    }
    
    ret_cnt <- ret_tib %>% base::ncol()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

getSsheetCoreAnnoTib = function(minOobPval,minOobPerc,
                                minNegPval,minNegPerc,
                                minDb,
                                verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getSsheetCoreAnnoTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    idx_zero <- 0
    ret_tib <- tibble::tibble(
      # Sample Requeue Suggestion::
      #
      Requeue_Flag_pOOBAH = 
        glue::glue("Flag to automatically requeue a faild sample based on percent loci with pOOBAH detection p-value < {minOobPval}. Percent threshold = {minOobPerc}."),
      Requeue_Flag_PnegEcdf = 
        glue::glue("Flag to automatically requeue a faild sample based on percent loci with PnegEcdf detection p-value < {minNegPval}. Percent threshold = {minNegPerc}."),
      
      # Sample Identification::
      #
      Sentrix_Name    = glue::glue("Unique Sentrix Name: Sentrix_Barcode + Sentrix_Poscode."),
      Sentrix_Barcode = glue::glue("Sentrix Bar Code (AKA Sentrix_ID)."),
      Sentrix_Poscode = glue::glue("Sentrix Position Code (AKA Sentrix_Position)."),
      Sentrix_Row     = glue::glue("Sentrix Row on chip."),
      Sentrix_Col     = glue::glue("Sentrix Column on Chip."),
      Chip_Type       = glue::glue("Idat Chip Type."),
      Chip_Format     = glue::glue("Idat Chip Format (e.g. 8x1/12x1/24x1/etc.)."),
      Bead_Pool       = glue::glue("Automatically identified bead pool based on address overlap (e.g. HM450/EPIC/etc.)."),
      Run_Name        = glue::glue("Run Name provided by user."),
      Sesame_Version  = glue::glue("Sesame Version."),
      
      # Bead and Manifest Loci Summary::
      #
      Bead_AvgRep_cg = "Average number of sample cg# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
      Bead_AvgRep_ch = "Average number of sample ch# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
      Bead_AvgRep_rs = "Average number of sample rs# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
      
      Bead_Count_cg  = "Total number of sample unique cg# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      Bead_Count_ch  = "Total number of sample unique ch# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      Bead_Count_rs  = "Total number of sample unique rs# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      
      Bead_Total_cg  = "Total number of sample cg# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      Bead_Total_ch  = "Total number of sample ch# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      Bead_Total_rs  = "Total number of sample rs# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      
      Loci_Count_cg  = "Total number of sample cg# loci overlapping manifest based on address presence NOT detection p-value.",
      Loci_Count_ch  = "Total number of sample ch# loci overlapping manifest based on address presence NOT detection p-value.",
      Loci_Count_rs  = "Total number of sample rs# loci overlapping manifest based on address presence NOT detection p-value.",
      
      # OLD Metrics from v.3.46:: Bead and Manifest Loci Summary::
      #
      CG_Manifest_Count = "Total number of sample cg# loci overlapping manifest based on address presence NOT detection p-value.",
      CH_Manifest_Count = "Total number of sample ch# loci overlapping manifest based on address presence NOT detection p-value.",
      RS_Manifest_Count = "Total number of sample rs# loci overlapping manifest based on address presence NOT detection p-value.",
      
      CG_Bead_Count     = "Total number of sample unique cg# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      CG_Bead_Total     = "Total number of sample cg# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      CG_Bead_AvgRep    = "Average number of sample cg# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
      
      CH_Bead_Count     = "Total number of sample unique ch# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      CH_Bead_Total     = "Total number of sample ch# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      CH_Bead_AvgRep    = "Average number of sample ch# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
      
      RS_Bead_Count     = "Total number of sample unique rs# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      RS_Bead_Total     = "Total number of sample rs# bead-types found overlapping manifest based on address presence NOT number of bead types.",
      RS_Bead_AvgRep    = "Average number of sample rs# bead-replicates found overlapping manifest based on address presence NOT number of bead types.",
      
      # Parameters Used::
      #
      # minNegPval = "Minimum Negative (PnegEcdf) detection p-value threshold used.",
      # minOobPval = "Minimum Out-Of-Band (pOOBAH) detection p-value threshold used.",
      minDeltaBeta  = "Minimum delta-Beta cutoff used for Sample-Auto-Detection calculations.",
      Min_DeltaBeta = "Minimum delta-Beta cutoff used for Sample-Auto-Detection calculations.",
      
      # Manifest, Platform and Version Auto Detection Results::
      #
      detect_manifest = "Platform name-version used. Identified during Auto-Manifest-Detection.",
      detect_platform = "Platform name used. Identified during Auto-Manifest-Detection.",
      detect_version  = "Platform version used. Identified during Auto-Manifest-Detection.",
      detect_sample_cnt   = "Number of Sample loci used during Auto-Manifest-Detection.",
      detect_manifest_cnt = "Number of Manifest loci used during Auto-Manifest-Detection.",
      detect_match_cnt    = "Number of Sample/Manifest overlapping loci used during Auto-Manifest-Detection.",
      detect_rc_per       = "Highest reciprical overlap of Sample/Manifest loci identified during Auto-Manifest-Detection.",
      
      # Open Sesame Basic Calling:: idx=0
      #
      cg_total_cnt_basic_0 = glue::glue("Only EPIC; Total cg loci using EPIC openSesame() for method index {idx_zero}. ",
                                          "See Bioconductor Package ‘sesame’ openSesame(sset)."),
      cg_pass_count_basic_0 = glue::glue("Only EPIC; Cout of recalculated cg loci with pOOBAH detection p-value <= {minOobPval} for method index {idx_zero}. ",
                                         "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset) and openSesame(sset)."),
      cg_pass_cnt_basic_0  = glue::glue("Only EPIC; Cout of cg loci with pOOBAH detection p-value <= {minOobPval} for method index {idx_zero}. ",
                                        "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset) and openSesame(sset)."),
      cg_pass_perc_basic_0 = glue::glue("Only EPIC; Percent of cg loci with pOOBAH detection p-value <= {minOobPval} for method index {idx_zero}. ",
                                        "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset) and openSesame(sset)."),
      
      cg_Failed_QC_cnt_basic_0 = glue::glue("Only EPIC; Flag count to automatically requeue a faild sample based on cg_pass_perc_basic."),
      cg_Failed_QC_basic_0     = glue::glue("Only EPIC; Flag to automatically requeue a faild sample based on cg_pass_perc_basic."),
      
      cg_1_total_cnt_basic_0 = glue::glue("Only EPIC; Total Infinium I cg loci using EPIC openSesame() for method index {idx_zero}. ",
                                          "See Bioconductor Package ‘sesame’ openSesame(sset)."),
      cg_2_total_cnt_basic_0 = glue::glue("Only EPIC; Total Infinium II cg loci using EPIC openSesame() for method index {idx_zero}. ",
                                          "See Bioconductor Package ‘sesame’ openSesame(sset)."),
      
      cg_1_pass_cnt_basic_0  = glue::glue("Only EPIC; Cout of Infinium I cg loci with pOOBAH detection p-value <= {minOobPval} for method index {idx_zero}. ",
                                          "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset) and openSesame(sset)."),
      cg_2_pass_cnt_basic_0  = glue::glue("Only EPIC; Cout of Infinium II cg loci with pOOBAH detection p-value <= {minOobPval} for method index {idx_zero}. ",
                                          "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset) and openSesame(sset)."),
      
      cg_1_pass_perc_basic_0 = glue::glue("Only EPIC; Percent of Infinium I cg loci with pOOBAH detection p-value <= {minOobPval} for method index {idx_zero}. ",
                                          "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset) and openSesame(sset)."),
      cg_2_pass_perc_basic_0 = glue::glue("Only EPIC; Percent of Infinium II cg loci with pOOBAH detection p-value <= {minOobPval} for method index {idx_zero}. ",
                                          "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset) and openSesame(sset)."),
      
      # OLD Platform Detection Results::
      #
      platformUsed = "Platform name used. Identified during Auto-Manifest-Detection.",
      platVersUsed = "Platform version used. Identified during Auto-Manifest-Detection.",
      
      # Auto-Sample-Detection Results::
      #
      AutoSample_Total_Cnt = "Loci overlap count with sample and max auto-detect-sample.",
      AutoSample_R2_Key = "Max auto-detect-sample name by R-Squared.",
      AutoSample_R2_Val = "Max auto-detect-sample R-Squared value.",
      AutoSample_dB_Key = "Max auto-detect-sample name by delta-Beta percent",
      AutoSample_dB_Cnt = glue::glue("Count of loci with delta-Beta < {minDb} between sample and max auto-detect-sample."),
      AutoSample_dB_Val = glue::glue("Percent of loci with delta-Beta < {minDb} between sample and max auto-detect-sample."),
      
      # Iscan Decode and Extration Dates::
      #
      Iscan_Decoding_ScannerID = "Iscan ID at time of decode.",
      Iscan_Decoding_Year = "Time stamp year of chip decode.",
      Iscan_Decoding_Mon  = "Time stamp month of chip decode.",
      Iscan_Decoding_Day  = "Time stamp day of chip decode.",
      Iscan_Decoding_Hour = "Time stamp hour of chip decode.",
      Iscan_Decoding_Min  = "Time stamp minute of chip decode.",
      Iscan_Decoding_Sec  = "Time stamp second of chip decode",
      
      Iscan_Extract_ScannerID = "Iscan ID at time of chip sample intensity extraction.",
      Iscan_Extract_Year = "Time stamp year of chip sample intensity extraction.",
      Iscan_Extract_Mon  = "Time stamp month of chip sample intensity extraction.",
      Iscan_Extract_Day  = "Time stamp day of chip sample intensity extraction.",
      Iscan_Extract_Hour = "Time stamp hour of chip sample intensity extraction.",
      Iscan_Extract_Min  = "Time stamp minute of chip sample intensity extraction.",
      Iscan_Extract_Sec  = "Time stamp second of chip decode"
    )
    
    ret_cnt <- ret_tib %>% base::ncol()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

getSsheetIndexAnnoTib = function(idx, 
                                 minOobPval,minOobPerc,
                                 minNegPval,minNegPerc,
                                 minDb,
                                 verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getSsheetIndexAnnoTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    ret_tib <- tibble::tibble(
      # Inference and Predictions::
      #
      Method_Key = glue::glue("Method workflow name used for metrics ending with {idx}."),
      Method_Idx = glue::glue("Method workflow index used for metrics ending with {idx}."),
      
      # Auto-Sample-Detection Results:: All
      #
      AutoSample_Total_Cnt = 
        glue::glue("Loci overlap count with sample and max auto-detect-sample ",
                   "for method index {idx}."),
      AutoSample_R2_Key = 
        glue::glue("Max auto-detect-sample name by R-Squared ",
                   "for method index {idx}."),
      AutoSample_R2_Val = 
        glue::glue("Max auto-detect-sample R-Squared value ",
                   "for method index {idx}."),
      AutoSample_dB_Key = 
        glue::glue("Max auto-detect-sample name by delta-Beta percent ",
                   "for method index {idx}."),
      AutoSample_dB_Cnt = 
        glue::glue("Count of loci with delta-Beta < {minDb} between sample ",
                   "and max auto-detect-sample for method index {idx}."),
      AutoSample_dB_Val = 
        glue::glue("Percent of loci with delta-Beta < {minDb} between sample ",
                   "and max auto-detect-sample for method index {idx}."),
      
      # Auto-Sample-Detection Results:: Infinium I
      #
      AutoSample_Total_1_Cnt = 
        glue::glue("Infinium I loci overlap count with sample and max auto-detect-sample ",
                   "for method index {idx}."),
      AutoSample_R2_1_Key = 
        glue::glue("Infinium I max auto-detect-sample name by R-Squared ",
                   "for method index {idx}."),
      AutoSample_R2_1_Val = 
        glue::glue("Infinium I max auto-detect-sample R-Squared value ",
                   "for method index {idx}."),
      AutoSample_dB_1_Key = 
        glue::glue("Infinium I max auto-detect-sample name by delta-Beta percent ",
                   "for method index {idx}."),
      AutoSample_dB_1_Cnt = 
        glue::glue("Infinium I count of loci with delta-Beta < {minDb} between sample ",
                   "and max auto-detect-sample for method index {idx}."),
      AutoSample_dB_1_Val = 
        glue::glue("Infinium I percent of loci with delta-Beta < {minDb} between sample ",
                   "and max auto-detect-sample for method index {idx}."),
      
      # Auto-Sample-Detection Results:: Infinium II
      #
      AutoSample_Total_2_Cnt = 
        glue::glue("Infinium II loci overlap count with sample and max auto-detect-sample ",
                   "for method index {idx}."),
      AutoSample_R2_2_Key = 
        glue::glue("Infinium II max auto-detect-sample name by R-Squared ",
                   "for method index {idx}."),
      AutoSample_R2_2_Val = 
        glue::glue("Infinium II max auto-detect-sample R-Squared value ",
                   "for method index {idx}."),
      AutoSample_dB_2_Key = 
        glue::glue("Infinium II max auto-detect-sample name by delta-Beta percent ",
                   "for method index {idx}."),
      AutoSample_dB_2_Cnt = 
        glue::glue("Infinium II count of loci with delta-Beta < {minDb} between sample ",
                   "and max auto-detect-sample for method index {idx}."),
      AutoSample_dB_2_Val = 
        glue::glue("Infinium II percent of loci with delta-Beta < {minDb} between sample ",
                   "and max auto-detect-sample for method index {idx}."),
      
      # Predicted and Inferred Calls::
      #
      GCT = glue::glue("GCT Score for method index {idx}. ",
                       "See Bioconductor Package sesame::GCT(sset): ",
                       "Compute GCT score for internal bisulfite conversion control. ",
                       "The function takes a SigSet as input. The higher the GCT score, the more likely the incomplete conversion."),
      
      Sex = glue::glue("Sex call for method index {idx}. ",
                       "See Bioconductor Package ‘sesame’ sesame::inferSex(sset)."),
      SexKaryotype = glue::glue("Sex Karyotype call for method index {idx}. ",
                                "See Bioconductor Package ‘sesame’ sesame::inferSexKaryotypes(sset)."),
      
      Ethnicity = glue::glue("Ethnicity call for method index {idx}. ",
                             "See Bioconductor Package ‘sesame’ sesame::inferEthnicity(sset).",
                             "This function uses both the built-in rsprobes as well as the type I Color-Channel-Switching probes to infer ethnicity."),
      
      SwapCntIGR = glue::glue("Number of Grn to Red swapped channels for method index {idx}. ",
                              "See Bioconductor Package ‘sesame’ sesame::inferTypeIChannel(sset).",
                              "Infer and reset color channel for Type-I probes instead of using what is specified in manifest."),
      SwapCntIRG = glue::glue("Number of Red to Grn swapped channels for method index {idx}. ",
                              "See Bioconductor Package ‘sesame’ sesame::inferTypeIChannel(sset).",
                              "Infer and reset color channel for Type-I probes instead of using what is specified in manifest."),
      
      AgeSkinBlood = glue::glue("Horvath Skin and Blood age predictor call for method index {idx}. ",
                                "See Bioconductor Package ‘sesame’ sesame::predictAgeSkinBlood(betas)."),
      
      # Detection P-values:: pOOBAH:: NEW Both
      #
      cg_calls_pass_perc = 
        glue::glue("Percent cg loci with pOOBAH detection p-value < {minOobPval} for method index {idx}, ",
                   "calculated directly from the calls file. ",
                   "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      cg_calls_pass_count = 
        glue::glue("Number cg loci with pOOBAH detection p-value < {minOobPval} for method index {idx}, ",
                   "calculated directly from the calls file. ",
                   "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      cg_calls_total_count =
        glue::glue("Total number cg loci calculated directly from the calls file."),
      cg_calls_min_cutoff = 
        glue::glue("Minium detection p-value for method index {idx} from direct calls file."),
      cg_calls_metric =
        glue::glue("Minium detection p-value metric for method index {idx} from direct calls file."),
      
      cg_pvals_PnegEcdf_pass_perc =
        glue::glue("Percent cg loci with PnegEcdf detection p-value < {minOobPval} for method index {idx}, ",
                   "calculated directly from the calls file. ",
                   "See Bioconductor Package ‘sesame’ sesame::PnegEcdf(sset)."),
      cg_pvals_pOOBAH_pass_perc =
        glue::glue("Percent cg loci with pOOBAH detection p-value < {minOobPval} for method index {idx}, ",
                   "calculated directly from the calls file. ",
                   "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      
      ch_pvals_PnegEcdf_pass_perc =
        glue::glue("Percent ch loci with PnegEcdf detection p-value < {minOobPval} for method index {idx}, ",
                   "calculated directly from the calls file. ",
                   "See Bioconductor Package ‘sesame’ sesame::PnegEcdf(sset)."),
      ch_pvals_pOOBAH_pass_perc =
        glue::glue("Percent ch loci with pOOBAH detection p-value < {minOobPval} for method index {idx}, ",
                   "calculated directly from the calls file. ",
                   "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      
      rs_pvals_PnegEcdf_pass_perc =
        glue::glue("Percent rs loci with PnegEcdf detection p-value < {minOobPval} for method index {idx}, ",
                   "calculated directly from the calls file. ",
                   "See Bioconductor Package ‘sesame’ sesame::PnegEcdf(sset)."),
      rs_pvals_pOOBAH_pass_perc =
        glue::glue("Percent rs loci with pOOBAH detection p-value < {minOobPval} for method index {idx}, ",
                   "calculated directly from the calls file. ",
                   "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),

      # Detection P-values:: pOOBAH:: NEW
      #
      cg_1_pvals_pOOBAH_pass_perc = glue::glue("Percent cg Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                               "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      cg_2_pvals_pOOBAH_pass_perc = glue::glue("Percent cg Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                               "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      
      ch_1_pvals_pOOBAH_pass_perc = glue::glue("Percent ch Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                               "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      ch_2_pvals_pOOBAH_pass_perc = glue::glue("Percent ch Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                               "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      
      rs_1_pvals_pOOBAH_pass_perc = glue::glue("Percent rs Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                               "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      rs_2_pvals_pOOBAH_pass_perc = glue::glue("Percent rs Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                               "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      
      # Detection P-values:: pOOBAH:: NEW:: mean
      #
      cg_1_pvals_pOOBAH_mean = glue::glue("Percent cg Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                          "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      cg_2_pvals_pOOBAH_mean = glue::glue("Percent cg Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                          "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      
      ch_1_pvals_pOOBAH_mean = glue::glue("Percent ch Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                          "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      ch_2_pvals_pOOBAH_mean = glue::glue("Percent ch Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                          "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      
      rs_1_pvals_pOOBAH_mean = glue::glue("Percent rs Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                          "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      rs_2_pvals_pOOBAH_mean = glue::glue("Percent rs Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                          "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      
      # Detection P-values:: pOOBAH:: NEW:: Requeue
      #
      cg_1_pvals_pOOBAH_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with pOOBAH detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for cg Infinium I loci."),
      cg_2_pvals_pOOBAH_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with pOOBAH detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for cg Infinium II loci."),
      
      ch_1_pvals_pOOBAH_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with pOOBAH detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for ch Infinium I loci."),
      ch_2_pvals_pOOBAH_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with pOOBAH detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for ch Infinium II loci."),
      
      rs_1_pvals_pOOBAH_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with pOOBAH detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for rs Infinium I loci."),
      rs_2_pvals_pOOBAH_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with pOOBAH detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for rs Infinium II loci."),
      
      # Detection P-values:: pOOBAH:: OLD
      #
      pOOBAH_cg_1_pass_perc = glue::glue("Percent cg Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                         "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      pOOBAH_cg_2_pass_perc = glue::glue("Percent cg Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                         "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      
      pOOBAH_ch_1_pass_perc = glue::glue("Percent ch Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                         "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      pOOBAH_ch_2_pass_perc = glue::glue("Percent ch Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                         "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      
      pOOBAH_rs_1_pass_perc = glue::glue("Percent rs Infinium I loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                         "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      pOOBAH_rs_2_pass_perc = glue::glue("Percent rs Infinium II loci with pOOBAH detection p-value < {minOobPval} for method index {idx}. ",
                                         "See Bioconductor Package ‘sesame’ sesame::pOOBAH(sset)."),
      
      # Detection P-values:: PnegEcdf:: NEW
      #
      cg_1_pvals_PnegEcdf_pass_perc = glue::glue("Percent cg Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                                 "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      cg_2_pvals_PnegEcdf_pass_perc = glue::glue("Percent cg Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                                 "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      
      ch_1_pvals_PnegEcdf_pass_perc = glue::glue("Percent ch Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                                 "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      ch_2_pvals_PnegEcdf_pass_perc = glue::glue("Percent ch Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                                 "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      
      rs_1_pvals_PnegEcdf_pass_perc = glue::glue("Percent rs Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                                 "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      rs_2_pvals_PnegEcdf_pass_perc = glue::glue("Percent rs Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                                 "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      
      # Detection P-values:: PnegEcdf:: NEW:: mean
      #
      cg_1_pvals_PnegEcdf_mean = glue::glue("Percent cg Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                            "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      cg_2_pvals_PnegEcdf_mean = glue::glue("Percent cg Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                            "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      
      ch_1_pvals_PnegEcdf_mean = glue::glue("Percent ch Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                            "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      ch_2_pvals_PnegEcdf_mean = glue::glue("Percent ch Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                            "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      
      rs_1_pvals_PnegEcdf_mean = glue::glue("Percent rs Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                            "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      rs_2_pvals_PnegEcdf_mean = glue::glue("Percent rs Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                            "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      
      # Detection P-values:: PnegEcdf:: NEW:: Requeue
      #
      cg_1_pvals_PnegEcdf_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with PnegEcdf detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for cg Infinium I loci."),
      cg_2_pvals_PnegEcdf_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with PnegEcdf detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for cg Infinium II loci."),
      
      ch_1_pvals_PnegEcdf_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with PnegEcdf detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for ch Infinium I loci."),
      ch_2_pvals_PnegEcdf_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with PnegEcdf detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for ch Infinium II loci."),
      
      rs_1_pvals_PnegEcdf_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with PnegEcdf detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for rs Infinium I loci."),
      rs_2_pvals_PnegEcdf_Requeue = 
        glue::glue("Flag to automatically requeue a faild sample based on ",
                   "percent loci with PnegEcdf detection p-value < {minOobPerc}. ",
                   "Percent threshold = {minOobPerc} for rs Infinium II loci."),      
      # Detection P-values:: PnegEcdf:: OLD
      #
      PnegEcdf_cg_1_pass_perc = glue::glue("Percent cg Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      PnegEcdf_cg_2_pass_perc = glue::glue("Percent cg Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      
      PnegEcdf_ch_1_pass_perc = glue::glue("Percent ch Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      PnegEcdf_ch_2_pass_perc = glue::glue("Percent ch Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      
      PnegEcdf_rs_1_pass_perc = glue::glue("Percent rs Infinium I loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      PnegEcdf_rs_2_pass_perc = glue::glue("Percent rs Infinium II loci with PnegEcdf detection p-value < {minNegPval} for method index {idx}. ",
                                           "See Bioconductor Package ‘sesame’ sesame::detectionPnegEcdf(sset)."),
      
      # Average Beta Values::
      #
      cg_1_betas_mean = glue::glue("Mean beta-value of cg Infinium I loci for method {idx}."),
      cg_2_betas_mean = glue::glue("Mean beta-value of cg Infinium II loci for method {idx}."),
      cg_1_mean = glue::glue("Mean beta-value of cg Infinium I loci for method {idx}."),
      cg_2_mean = glue::glue("Mean beta-value of cg Infinium II loci for method {idx}."),
      
      ch_1_betas_mean = glue::glue("Mean beta-value of ch Infinium I loci for method {idx}."),
      ch_2_betas_mean = glue::glue("Mean beta-value of ch Infinium II loci for method {idx}."),
      ch_1_mean = glue::glue("Mean beta-value of ch Infinium I loci for method {idx}."),
      ch_2_mean = glue::glue("Mean beta-value of ch Infinium II loci for method {idx}."),
      
      rs_1_betas_mean = glue::glue("Mean beta-value of rs Infinium I loci for method {idx}."),
      rs_2_betas_mean = glue::glue("Mean beta-value of rs Infinium II loci for method {idx}."),
      rs_1_mean = glue::glue("Mean beta-value of rs Infinium I loci for method {idx}."),
      rs_2_mean = glue::glue("Mean beta-value of rs Infinium II loci for method {idx}."),
      
      # Average Intensities:: NEW
      #
      BISULFITE_CONVERSION_I_2_sig_M_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_I Grn loci for method {idx}."),
      BISULFITE_CONVERSION_I_2_sig_U_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_I Red loci for method {idx}."),
      
      BISULFITE_CONVERSION_II_2_sig_M_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_II Grn loci for method {idx}."),
      BISULFITE_CONVERSION_II_2_sig_U_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_II Red loci for method {idx}."),
      
      cg_2_sig_M_mean = glue::glue("Mean intensity of cg Infinium II Grn (methylated) loci for method {idx}."),
      cg_2_sig_U_mean = glue::glue("Mean intensity of cg Infinium II Red (unmethylated) loci for method {idx}."),
      
      cg_IG_sig_M_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Grn (methylated) loci for method {idx}."),
      cg_IG_sig_U_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Grn (unmethylated) loci for method {idx}."),
      cg_IR_sig_M_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Red (methylated) loci for method {idx}."),
      cg_IR_sig_U_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Red (unmethylated) loci for method {idx}."),
      
      cg_OG_sig_M_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Grn (methylated) loci for method {idx}."),
      cg_OG_sig_U_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Grn (unmethylated) loci for method {idx}."),
      cg_OR_sig_M_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Red (methylated) loci for method {idx}."),
      cg_OR_sig_U_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Red (unmethylated) loci for method {idx}."),
      
      ch_2_sig_M_mean = glue::glue("Mean intensity of ch Infinium II Grn (methylated) loci for method {idx}."),
      ch_2_sig_U_mean = glue::glue("Mean intensity of ch Infinium II Red (unmethylated) loci for method {idx}."),
      
      EXTENSION_2_sig_M_mean = glue::glue("Mean intensity of EXTENSION Grn loci for method {idx}."),
      EXTENSION_2_sig_U_mean = glue::glue("Mean intensity of EXTENSION Red loci for method {idx}."),
      
      HYBRIDIZATION_2_sig_M_mean = glue::glue("Mean intensity of HYBRIDIZATION Grn loci for method {idx}."),
      HYBRIDIZATION_2_sig_U_mean = glue::glue("Mean intensity of HYBRIDIZATION Red loci for method {idx}."),
      
      NEGATIVE_2_sig_M_mean = glue::glue("Mean intensity of NEGATIVE Grn loci for method {idx}."),
      NEGATIVE_2_sig_U_mean = glue::glue("Mean intensity of NEGATIVE Red loci for method {idx}."),
      
      NON_POLYMORPHIC_2_sig_M_mean = glue::glue("Mean intensity of NON_POLYMORPHIC Grn loci for method {idx}."),
      NON_POLYMORPHIC_2_sig_U_mean = glue::glue("Mean intensity of NON_POLYMORPHIC Red loci for method {idx}."),
      
      NORM_A_2_sig_M_mean = glue::glue("Mean intensity of NORM_A Grn loci for method {idx}."),
      NORM_A_2_sig_U_mean = glue::glue("Mean intensity of NORM_A Red loci for method {idx}."),
      
      NORM_C_2_sig_M_mean = glue::glue("Mean intensity of NORM_C Grn loci for method {idx}."),
      NORM_C_2_sig_U_mean = glue::glue("Mean intensity of NORM_C Red loci for method {idx}."),
      
      NORM_G_2_sig_M_mean = glue::glue("Mean intensity of NORM_G Grn loci for method {idx}."),
      NORM_G_2_sig_U_mean = glue::glue("Mean intensity of NORM_G Red loci for method {idx}."),
      
      NORM_T_2_sig_M_mean = glue::glue("Mean intensity of NORM_T Grn loci for method {idx}."),
      NORM_T_2_sig_U_mean = glue::glue("Mean intensity of NORM_T Red loci for method {idx}."),
      
      RESTORATION_2_sig_M_mean = glue::glue("Mean intensity of RESTORATION Grn loci for method {idx}."),
      RESTORATION_2_sig_U_mean = glue::glue("Mean intensity of RESTORATION Red loci for method {idx}."),
      
      rs_2_sig_M_mean = glue::glue("Mean intensity of rs Infinium II Grn (methylated) loci for method {idx}."),
      rs_2_sig_U_mean = glue::glue("Mean intensity of rs Infinium II Red (unmethylated) loci for method {idx}."),
      
      rs_IG_sig_M_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Grn (methylated) loci for method {idx}."),
      rs_IG_sig_U_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Grn (unmethylated) loci for method {idx}."),
      rs_IR_sig_M_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Red (methylated) loci for method {idx}."),
      rs_IR_sig_U_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Red (unmethylated) loci for method {idx}."),
      
      rs_OG_sig_M_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Grn (methylated) loci for method {idx}."),
      rs_OG_sig_U_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Grn (unmethylated) loci for method {idx}."),
      rs_OR_sig_M_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Red (methylated) loci for method {idx}."),
      rs_OR_sig_U_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Red (unmethylated) loci for method {idx}."),
      
      SPECIFICITY_I_2_sig_M_mean = glue::glue("Mean intensity of SPECIFICITY_I Grn loci for method {idx}."),
      SPECIFICITY_I_2_sig_U_mean = glue::glue("Mean intensity of SPECIFICITY_I Red loci for method {idx}."),
      
      SPECIFICITY_II_2_sig_M_mean = glue::glue("Mean intensity of SPECIFICITY_II Grn loci for method {idx}."),
      SPECIFICITY_II_2_sig_U_mean = glue::glue("Mean intensity of SPECIFICITY_II Red loci for method {idx}."),
      
      STAINING_2_sig_M_mean = glue::glue("Mean intensity of STAINING Grn loci for method {idx}."),
      STAINING_2_sig_U_mean = glue::glue("Mean intensity of STAINING Red loci for method {idx}."),
      
      TARGET_REMOVAL_2_sig_M_mean = glue::glue("Mean intensity of TARGET_REMOVAL Grn loci for method {idx}."),
      TARGET_REMOVAL_2_sig_U_mean = glue::glue("Mean intensity of TARGET_REMOVAL Red loci for method {idx}."),
      
      
      # Average Intensities:: OLD
      #
      BISULFITE_CONVERSION_I_2_M_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_I Grn loci for method {idx}."),
      BISULFITE_CONVERSION_I_2_U_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_I Red loci for method {idx}."),
      
      BISULFITE_CONVERSION_II_2_M_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_II Grn loci for method {idx}."),
      BISULFITE_CONVERSION_II_2_U_mean = glue::glue("Mean intensity of BISULFITE_CONVERSION_II Red loci for method {idx}."),
      
      cg_2_M_mean = glue::glue("Mean intensity of cg Infinium II Grn (methylated) loci for method {idx}."),
      cg_2_U_mean = glue::glue("Mean intensity of cg Infinium II Red (unmethylated) loci for method {idx}."),
      
      cg_IG_M_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Grn (methylated) loci for method {idx}."),
      cg_IG_U_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Grn (unmethylated) loci for method {idx}."),
      cg_IR_M_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Red (methylated) loci for method {idx}."),
      cg_IR_U_mean = glue::glue("Mean intensity of cg Infinium I In-Bound Red (unmethylated) loci for method {idx}."),
      
      cg_OG_M_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Grn (methylated) loci for method {idx}."),
      cg_OG_U_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Grn (unmethylated) loci for method {idx}."),
      cg_OR_M_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Red (methylated) loci for method {idx}."),
      cg_OR_U_mean = glue::glue("Mean intensity of cg Infinium I Out-Of-Bound Red (unmethylated) loci for method {idx}."),
      
      ch_2_M_mean = glue::glue("Mean intensity of ch Infinium II Grn (methylated) loci for method {idx}."),
      ch_2_U_mean = glue::glue("Mean intensity of ch Infinium II Red (unmethylated) loci for method {idx}."),
      
      EXTENSION_2_M_mean = glue::glue("Mean intensity of EXTENSION Grn loci for method {idx}."),
      EXTENSION_2_U_mean = glue::glue("Mean intensity of EXTENSION Red loci for method {idx}."),
      
      HYBRIDIZATION_2_M_mean = glue::glue("Mean intensity of HYBRIDIZATION Grn loci for method {idx}."),
      HYBRIDIZATION_2_U_mean = glue::glue("Mean intensity of HYBRIDIZATION Red loci for method {idx}."),
      
      NEGATIVE_2_M_mean = glue::glue("Mean intensity of NEGATIVE Grn loci for method {idx}."),
      NEGATIVE_2_U_mean = glue::glue("Mean intensity of NEGATIVE Red loci for method {idx}."),
      
      NON_POLYMORPHIC_2_M_mean = glue::glue("Mean intensity of NON_POLYMORPHIC Grn loci for method {idx}."),
      NON_POLYMORPHIC_2_U_mean = glue::glue("Mean intensity of NON_POLYMORPHIC Red loci for method {idx}."),
      
      NORM_A_2_M_mean = glue::glue("Mean intensity of NORM_A Grn loci for method {idx}."),
      NORM_A_2_U_mean = glue::glue("Mean intensity of NORM_A Red loci for method {idx}."),
      
      NORM_C_2_M_mean = glue::glue("Mean intensity of NORM_C Grn loci for method {idx}."),
      NORM_C_2_U_mean = glue::glue("Mean intensity of NORM_C Red loci for method {idx}."),
      
      NORM_G_2_M_mean = glue::glue("Mean intensity of NORM_G Grn loci for method {idx}."),
      NORM_G_2_U_mean = glue::glue("Mean intensity of NORM_G Red loci for method {idx}."),
      
      NORM_T_2_M_mean = glue::glue("Mean intensity of NORM_T Grn loci for method {idx}."),
      NORM_T_2_U_mean = glue::glue("Mean intensity of NORM_T Red loci for method {idx}."),
      
      RESTORATION_2_M_mean = glue::glue("Mean intensity of RESTORATION Grn loci for method {idx}."),
      RESTORATION_2_U_mean = glue::glue("Mean intensity of RESTORATION Red loci for method {idx}."),
      
      rs_2_M_mean = glue::glue("Mean intensity of rs Infinium II Grn (methylated) loci for method {idx}."),
      rs_2_U_mean = glue::glue("Mean intensity of rs Infinium II Red (unmethylated) loci for method {idx}."),
      
      rs_IG_M_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Grn (methylated) loci for method {idx}."),
      rs_IG_U_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Grn (unmethylated) loci for method {idx}."),
      rs_IR_M_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Red (methylated) loci for method {idx}."),
      rs_IR_U_mean = glue::glue("Mean intensity of rs Infinium I In-Bound Red (unmethylated) loci for method {idx}."),
      
      rs_OG_M_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Grn (methylated) loci for method {idx}."),
      rs_OG_U_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Grn (unmethylated) loci for method {idx}."),
      rs_OR_M_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Red (methylated) loci for method {idx}."),
      rs_OR_U_mean = glue::glue("Mean intensity of rs Infinium I Out-Of-Bound Red (unmethylated) loci for method {idx}."),
      
      SPECIFICITY_I_2_M_mean = glue::glue("Mean intensity of SPECIFICITY_I Grn loci for method {idx}."),
      SPECIFICITY_I_2_U_mean = glue::glue("Mean intensity of SPECIFICITY_I Red loci for method {idx}."),
      
      SPECIFICITY_II_2_M_mean = glue::glue("Mean intensity of SPECIFICITY_II Grn loci for method {idx}."),
      SPECIFICITY_II_2_U_mean = glue::glue("Mean intensity of SPECIFICITY_II Red loci for method {idx}."),
      
      STAINING_2_M_mean = glue::glue("Mean intensity of STAINING Grn loci for method {idx}."),
      STAINING_2_U_mean = glue::glue("Mean intensity of STAINING Red loci for method {idx}."),
      
      TARGET_REMOVAL_2_M_mean = glue::glue("Mean intensity of TARGET_REMOVAL Grn loci for method {idx}."),
      TARGET_REMOVAL_2_U_mean = glue::glue("Mean intensity of TARGET_REMOVAL Red loci for method {idx}."),
      
      # Comparison to Open Sesame by R2/deltaBeta for each method
      cg_basic_r2 = glue::glue("R-squared of core open-sesame loci vs method {idx}."),
      cg_basic_dB = glue::glue("DeltaBeta of core open-sesame loci vs method {idx}."),
      
      unknown = "Unknown"
    ) %>% purrr::set_names(paste(names(.),idx, sep='_'))
    
    ret_cnt <- ret_tib %>% base::ncol()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}


# End of file
