
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Machine Learning Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )
suppressWarnings(suppressPackageStartupMessages(require("grid")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

# XGBoost::
# suppressWarnings(suppressPackageStartupMessages(require("xgboost")) )

# glmnet:: Elastic Search::
suppressWarnings(suppressPackageStartupMessages(require("glmnet")) )

# randomForest::
suppressWarnings(suppressPackageStartupMessages(require("randomForest")) )
suppressWarnings(suppressPackageStartupMessages(require("ROCR")) )

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Masking and Imputation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getMaskedTib = function(mat, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getMaskedTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  mis_tib <- tibble::tibble(Probe_ID=character(), idx=integer())
  stime <- system.time({
    
    mis_len <- 0
    mis_vec <- NULL
    mis_dat <- which(is.na(mat) | !is.finite(mat), arr.ind=TRUE)
    mis_vec <- which(is.na(mat) | !is.finite(mat) )
    # cut_tib <- mis_dat %>% as.data.frame() %>% tibble::as_tibble(rownames="Probe_ID") %>% tibble::add_column(idx=mis_vec)
    # mis_tib <- dplyr::bind_rows(mis_tib,cut_tib)
    
    if (!is.null(mis_vec) && length(mis_vec)>0){
      mis_tib <- mis_dat %>% as.data.frame() %>% tibble::as_tibble(rownames="Probe_ID") %>% tibble::add_column(idx=mis_vec)
      mis_len <- mis_tib %>% base::nrow()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; mis_len={mis_len}; elapsed={etime}.{RET}{RET}"))
  
  mis_tib
}

impute_matrix_mean <- function(mat, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'impute_matrix_mean'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  stime <- system.time({
    missing <- which(is.na(mat) | !is.finite(mat), arr.ind=TRUE)
    mis_cnt <- length(missing)
    if (length(missing)!=0){
      for (j in 1:nrow(missing)){
        mean <- mean(mat[missing[j,1],][is.finite(mat[missing[j,1],])], na.rm=TRUE)
        mat[missing[j,1],missing[j,2]] <- mean
        # cat(glue::glue("j={j}/{mis_cnt}, mean={mean}{RET}"))
      }
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))

  mat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Machine Learning IO::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadFromFileTib = function(tib, type, dir=NULL,
                           verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'loadFromFileTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  data <- NULL
  stime <- system.time({
    sel_tib <- tib %>% dplyr::filter(Type==type) %>% head(n=1)
    
    file <- file.path(sel_tib$Dir[1], sel_tib$File_Name[1])
    if (!is.null(dir)) file <- file.path(dir, sel_tib$File_Name[1])
    if (!file.exists(file))
      stop(glue::glue("{RET}[{funcTag}]: ERROR: File does NOT exist; file={file}!{RET}{RET}"))
    
    if (stringr::str_ends(file, '.rds')) {
      data <- readr::read_rds(file)
    } else if (stringr::str_ends(file, '.csv') || stringr::str_ends(file, '.csv.gz')) {
      data <- suppressMessages(suppressWarnings( readr::read_csv(file) ))
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported file type={file}!{RET}{RET}"))
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  data
}

saveModel = function(name, dir, type='model', 
                     mod=NULL, cgns=NULL, ss=NULL, pars=NULL, call=NULL, sums=NULL, ms=NULL,
                     method=NULL, alpha=NULL, lambda=NULL,
                     verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'saveModel'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  file_tib <- tibble::tibble( Type=character(), File_Name=character(), Dir=character() )
  stime <- system.time({
    mod_rds <- paste(name, paste(type, sep='-'), 'rds', sep='.')
    sam_csv <- paste(name, paste(type,'SampleSheet',sep='-'), 'csv.gz', sep='.')
    cgn_csv <- paste(name, paste(type,'features',sep='-'), 'csv.gz', sep='.')
    par_csv <- paste(name, paste(type,'params',sep='-'), 'csv.gz', sep='.')
    mss_csv <- paste(name, paste(type,'MethodSheet',sep='-'), 'csv.gz', sep='.')
    fns_csv <- paste(name, paste(type,'files',sep='-'), 'csv.gz', sep='.')

    call_csv <- paste(name, paste(type,'calls',sep='-'), 'csv.gz', sep='.')
    sums_csv <- paste(name, paste(type,'calls-summary',sep='-'), 'csv.gz', sep='.')
    
    if (!is.null(mod)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving {type} Structs RDS={mod_rds}...{RET}") )
      readr::write_rds(mod, file.path(dir,mod_rds), compress="gz" )
      file_tib <- file_tib %>% tibble::add_row( Type="Model",File_Name=mod_rds, Dir=dir )
    }
    
    if (!is.null(cgns)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving {type} CG Features={cgn_csv}...{RET}") )
      readr::write_csv(cgns, file.path(dir,cgn_csv) )
      file_tib <- file_tib %>% tibble::add_row( Type="Features",File_Name=cgn_csv, Dir=dir )
    }
    
    if (!is.null(ss)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving {type} SampleSheet={sam_csv}...{RET}") )
      readr::write_csv(ss, file.path(dir,sam_csv) )
      file_tib <- file_tib %>% tibble::add_row( Type="SampleSheet",File_Name=sam_csv, Dir=dir )
    }
    
    if (!is.null(pars)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving {type} Parameters={par_csv}...{RET}") )
      readr::write_csv(pars, file.path(dir,par_csv) )
      file_tib <- file_tib %>% tibble::add_row( Type="Params",File_Name=par_csv, Dir=dir )
    }
    
    if (!is.null(call)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving {type} Calls={call_csv}...{RET}") )
      readr::write_csv(call, file.path(dir,call_csv) )
      file_tib <- file_tib %>% tibble::add_row( Type="Calls",File_Name=call_csv, Dir=dir )
    }
    
    if (!is.null(sums)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving {type} Summary={sums_csv}...{RET}") )
      readr::write_csv(sums, file.path(dir,sums_csv) )
      file_tib <- file_tib %>% tibble::add_row( Type="Summary",File_Name=sums_csv, Dir=dir )
    }

    if (!is.null(ms)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving {type} Method Sheet={mss_csv}...{RET}") )
      readr::write_csv(ms, file.path(dir,mss_csv) )
      file_tib <- file_tib %>% tibble::add_row( Type="MethodSheet",File_Name=mss_csv, Dir=dir )
    }
    
    if (!is.null(method)) file_tib <- file_tib %>% tibble::add_column(method=method)
    if (!is.null(alpha))  file_tib <- file_tib %>% tibble::add_column(alpha=alpha)
    if (!is.null(lambda)) file_tib <- file_tib %>% tibble::add_column(lambda=lambda)
    
    # Write Files tib::
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving Model Files Table={fns_csv}...{RET}") )
    readr::write_csv(file_tib, file.path(dir,fns_csv) )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  file_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Cross Train Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# grp_ss_tib <- addCrossTrainGroups(ss=sampleSheet_tib, class=class_var, verbose=10)
addCrossTrainGroups = function(ss, class, min=18,
                                verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'addCrossTrainGroups'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; min={min}.{RET}"))

  # Determine the Sample Partition (2 or 3 based on minimum Sample Class Count)
  class_sum_ss <- ss %>% 
    dplyr::group_by(!!class) %>% 
    dplyr::summarise(Class_Count=n()) %>% 
    dplyr::ungroup() %>%
    dplyr::arrange(Class_Count) %>%
    purrr::set_names(c('Class_Name', 'Class_Count'))
  
  min_class_cnt <- class_sum_ss %>% head(n=1) %>% dplyr::pull(Class_Count) %>% as.integer()
  
  # Add Group Partition to Sample Sheet
  mod_idx <- NULL
  ss_tib <- NULL
  if (min_class_cnt<min) {
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Will do two-fold cross valudation min_class_cnt={min_class_cnt} < {min}.{RET}"))
    mod_idx <- 2

  } else {
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Will do three-fold cross valudation min_class_cnt={min_class_cnt} >= {min}.{RET}"))
    mod_idx <- 3
  }
  ss_tib <- ss %>% dplyr::mutate(Row=row_number(), Cross_Group=Row %% !!mod_idx)
  
  return(ss_tib)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Simplified glmnet->RandomForest Training::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

glmnetRforestWrapper = function(beta, ss, cgns, pars,
                                class_idx, seed, outName, dir,
                                alpha_min=1, alpha_max=11, alpha_inc=0.1, cgnMax=15000, cgnMin=6,
                                
                                parallel=FALSE,crossTrain=FALSE,class_var=NULL,cross_perc_min=NULL,
                                verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'glmnetRforestWrapper'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; seed={seed}.{RET}"))

  file_tib  <- tibble::tibble()
  stime <- system.time({
    cgns_cnt  <- cgns %>% base::nrow()
    alpha_vec <- c(alpha_min, alpha_max, seq(alpha_min, alpha_max, alpha_inc) ) %>% unique() %>% sort()

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Build Cross Validation Sets::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    tests_samp_tib <- list()
    train_samp_tib <- list()
    if (crossTrain && !is.null(class_var)) {
      grp_ss_tib <- addCrossTrainGroups(ss=ss, class=class_var, 
                                        verbose=verbose,vt=4,tc=tc+1,tt=tt)
      groups <- grp_ss_tib %>% dplyr::distinct(Cross_Group) %>% dplyr::pull(Cross_Group) %>% as.vector()
      for (group in groups) {
        group_key <- as.character(group)
        
        tests_ss_tib = grp_ss_tib %>% dplyr::filter(Cross_Group==group)
        train_ss_tib = grp_ss_tib %>% dplyr::filter(Cross_Group!=group)
        
        tests_samp_tib[[group_key]] = tests_ss_tib
        train_samp_tib[[group_key]] = train_ss_tib
        
        # train_beta_mat[[group_key]] = beta[ train_ss_tib$Sentrix_Name, ]
        # tests_beta_mat[[group_key]] = beta[ tests_ss_tib$Sentrix_Name, ]
      }
    } else {
      # return(FALSE)
    }

    verbose <- verbose + 20
    parallel <- FALSE
    if (parallel && cgns_cnt <= cgnMax) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Parallel Processing; cgns_cnt={cgns_cnt}...{RET}"))
      
      file_tib <- foreach (alpha_val=alpha_vec, .combine = rbind) %dopar% {
        alpha_dbl <- as.double( alpha_val )
        alpha_dir <- file.path(dir, paste('alpha',alpha_dbl, sep='-') )
        if (!dir.exists(alpha_dir)) dir.create(alpha_dir, recursive=TRUE)
        unlink(list.files(alpha_dir, full.names=TRUE))
        
        cgn_min_tib <- cgns
        cgn_1se_tib <- cgns
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                        Train Glmnet:: Elastic Net
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        if (cgns_cnt <= cgnMax) {
          files_glm_tib <- trainGlmnet(alpha=alpha_dbl, data=beta, ss=ss, cgns=cgn_1se_tib, class_idx=class_idx, seed=seed,
                                       dir=alpha_dir,retFiles=TRUE, pars=pars,
                                       crossTrain=TRUE,tests_ss=tests_samp_tib,train_ss=train_samp_tib,
                                       cross_perc_min=cross_perc_min,
                                       verbose=verbose,vt=3,tc=tc+1,tt=tt)
          file_tib <- file_tib %>% dplyr::bind_rows(files_glm_tib)
          
          glm_cgn_tib <- loadFromFileTib(tib=files_glm_tib, type="Features", verbose=verbose,vt=4,tc=tc+1,tt=tt)
          cgn_min_tib <- glm_cgn_tib %>% dplyr::filter(lambda.min!=0)
          cgn_1se_tib <- glm_cgn_tib %>% dplyr::filter(lambda.1se!=0)
        }
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                        Train RandomForest:: lambda.min
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        cgn_min_cnt <- cgn_min_tib %>% base::nrow()
        if (cgn_min_cnt > cgnMin) {
          files_min_tib <- trainRandomForest(alpha=alpha_dbl, data=beta, ss=ss, cgns=cgn_min_tib, class_idx, seed=seed,
                                             dir=alpha_dir,retFiles=TRUE,lambda="lambda-min", pars=pars,
                                             crossTrain=TRUE,tests_ss=tests_samp_tib,train_ss=train_samp_tib,
                                             cross_perc_min=cross_perc_min,
                                             verbose=verbose,vt=3,tc=tc+1,tt=tt)
          file_tib <- file_tib %>% dplyr::bind_rows(files_min_tib)
        }
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                        Train RandomForest:: lambda.1se
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        cgn_1se_cnt <- cgn_1se_tib %>% base::nrow()
        if (cgn_1se_cnt > cgnMin) {
          files_1se_tib <- trainRandomForest(alpha=alpha_dbl, data=beta, ss=ss, cgns=cgn_1se_tib, class_idx, seed=seed,
                                             dir=alpha_dir,retFiles=TRUE,lambda="lambda-1se", pars=pars,
                                             crossTrain=TRUE,tests_ss=tests_samp_tib,train_ss=train_samp_tib,
                                             cross_perc_min=cross_perc_min,
                                             verbose=verbose,vt=3,tc=tc+1,tt=tt)
          file_tib <- file_tib %>% dplyr::bind_rows(files_1se_tib)
        }
        
        if (cgns_cnt > cgnMax) break
        # if (alpha_idx>=2) break
      }
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. Parallel Processing.{RET}{RET}"))
      
    } else {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Linear Processing; cgns_cnt={cgns_cnt}...{RET}"))
      
      for (alpha_val in alpha_vec) {
        alpha_dbl <- as.double( alpha_val )
        alpha_dir <- file.path(dir, paste('alpha',alpha_dbl, sep='-') )
        if (!dir.exists(alpha_dir)) dir.create(alpha_dir, recursive=TRUE)
        unlink(list.files(alpha_dir, full.names=TRUE))
        
        cgn_min_tib <- cgns
        cgn_1se_tib <- cgns
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                        Train Glmnet:: Elastic Net
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        if (cgns_cnt <= cgnMax) {
          files_glm_tib <- trainGlmnet(alpha=alpha_dbl, data=beta, ss=ss, cgns=cgn_1se_tib, class_idx=class_idx, seed=seed,
                                       dir=alpha_dir,retFiles=TRUE, pars=pars,
                                       crossTrain=TRUE,tests_ss=tests_samp_tib,train_ss=train_samp_tib,
                                       cross_perc_min=cross_perc_min,
                                       verbose=verbose,vt=3,tc=tc+1,tt=tt)
          file_tib <- file_tib %>% dplyr::bind_rows(files_glm_tib)

          glm_cgn_tib <- loadFromFileTib(tib=files_glm_tib, type="Features", verbose=verbose,vt=4,tc=tc+1,tt=tt)
          cgn_min_tib <- glm_cgn_tib %>% dplyr::filter(lambda.min!=0)
          cgn_1se_tib <- glm_cgn_tib %>% dplyr::filter(lambda.1se!=0)
        }
        # return(files_glm_tib)
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                        Train RandomForest:: lambda.min
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        cgn_min_cnt <- cgn_min_tib %>% base::nrow()
        if (cgn_min_cnt > cgnMin) {
          files_min_tib <- trainRandomForest(alpha=alpha_dbl, data=beta, ss=ss, cgns=cgn_min_tib, class_idx, seed=seed,
                                             dir=alpha_dir,retFiles=TRUE,lambda="lambda-min", pars=pars,
                                             crossTrain=TRUE,tests_ss=tests_samp_tib,train_ss=train_samp_tib,
                                             cross_perc_min=cross_perc_min,
                                             verbose=verbose,vt=3,tc=tc+1,tt=tt)
          file_tib <- file_tib %>% dplyr::bind_rows(files_min_tib)
        }
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                        Train RandomForest:: lambda.1se
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        cgn_1se_cnt <- cgn_1se_tib %>% base::nrow()
        if (cgn_1se_cnt > cgnMin) {
          files_1se_tib <- trainRandomForest(alpha=alpha_dbl, data=beta, ss=ss, cgns=cgn_1se_tib, class_idx, seed=seed,
                                             dir=alpha_dir,retFiles=TRUE,lambda="lambda-1se", pars=pars,
                                             crossTrain=TRUE,tests_ss=tests_samp_tib,train_ss=train_samp_tib,
                                             cross_perc_min=cross_perc_min,
                                             verbose=verbose,vt=3,tc=tc+1,tt=tt)
          file_tib <- file_tib %>% dplyr::bind_rows(files_1se_tib)
        }
        
        if (cgns_cnt > cgnMax) break
        # if (alpha_idx>=2) break
      }
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. Linear Processing.{RET}{RET}"))
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  file_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Cross Training Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


crossTrain = function(name, alpha, data, ss_train, ss_tests, 
                      class_idx, seed=NULL,
                      type.measure="class", family="binomial", parallel=TRUE,
                      
                      dir=NULL,lambda=NULL,
                      modName="glmnet",lambda_vec=c("lambda.min","lambda.1se"),
                      
                      verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'crossTrain'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))

  cross_min_perc <- NULL
  cross_call_tib <- NULL
  cross_sums_tib <- NULL
  
  stime <- system.time({
    
    # Set seed for reproducibility
    if (!is.null(seed)) set.seed(seed)
    
    # TBD: Loop over lambda's properly::
    if (is.null(lambda)) lambda <- lambda_vec[2]
    
    grps <- ss_train %>% names()
    for (grp in grps) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}{TAB} Cross Group={grp}.{RET}"))
      
      cur_group <- paste(name,"-Part",grp, sep='')
      
      train_sam_tib = ss_train[[grp]]
      train_lab_idx = train_sam_tib %>% dplyr::pull(!!class_idx) %>% as.vector()
      train_dat_mat = data[ train_sam_tib$Sentrix_Name, ]
      
      tests_sam_tib = ss_tests[[grp]]
      tests_lab_idx = tests_sam_tib %>% dplyr::pull(!!class_idx) %>% as.vector()
      tests_dat_mat = data[ tests_sam_tib$Sentrix_Name, ]
      
      # trainGlmnet(alpha=alpha, data=train_dat_mat, labs=train_sam_tib, cgns=cgns, testData=tests_dat_mat, testLabs=tests_sam_tib, seed=seed, type.measure=type.measure, family=family, parallel=parallel, outName=name, dir=dir, lambda=lambda)
      
      if (modName=='glmnet') {
        cur_mod = glmnet::cv.glmnet(x=train_dat_mat, y=train_lab_idx,
                                    type.measure=type.measure,alpha=alpha,family=family,parallel=parallel)
        
        cur_pred = predGlmnet(mod=cur_mod, data=tests_dat_mat, labs=tests_lab_idx, 
                              name=cur_group, lambda="lambda.1se", type=type.measure,
                              verbose=verbose,vt=4,tc=tc+1,tt=tt) %>% dplyr::mutate(CrossGroup=cur_group)
        
        cur_call_tib <- predToCalls(pred=cur_pred, labs=tests_lab_idx, pred_lab="Pred_Class",
                                    verbose=verbose,vt=4,tc=tc+1,tt=tt)
      } else if (modName=='rforest') {

        data_df  <- cbind(train_dat_mat, train_lab_idx) %>% as.data.frame()
        df_ncols <- data_df %>% base::ncol()
        colnames(data_df)[df_ncols] <- 'Species'
        data_df$Species <- factor(data_df$Species)
        
        cur_mod = randomForest(Species ~ ., data=data_df, ntree=df_ncols-1, importance=TRUE)
        
        cur_pred = predRandomForest(mod=cur_mod, data=tests_dat_mat, labs=tests_lab_idx, 
                              name=cur_group, # lambda="lambda.1se", type=type.measure,
                              verbose=verbose,vt=4,tc=tc+1,tt=tt) %>% dplyr::mutate(CrossGroup=cur_group)
        
        cur_call_tib <- predToCalls(pred=cur_pred, labs=tests_lab_idx, pred_lab="Pred_Class",
                                    verbose=verbose,vt=4,tc=tc+1,tt=tt)
        
      } else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported modName={modName}!!!{RET}{RET}"))
      }
      
      cur_sum_tib <- callToSumTib(call=cur_call_tib, name_lab=cur_group, true_lab="True_Class",call_lab="Call",
                                  verbose=verbose,vt=4,tc=tc+1,tt=tt)
      
      cross_call_tib <- cross_call_tib %>% dplyr::bind_rows(cur_call_tib)
      cross_sums_tib <- cross_sums_tib %>% dplyr::bind_rows(cur_sum_tib)
    }
    
    # Build Cross Valdation Sample Sheet::
    cross_samp_tib <- dplyr::bind_cols(
      cross_sums_tib %>% tidyr::unite(Key, CrossGroup,True_Class,Call, sep='_') %>% 
        dplyr::select(Key, Group_Perc) %>% tidyr::spread(Key, Group_Perc) %>% purrr::set_names(paste(names(.),'Perc', sep='_')),
      cross_sums_tib %>% tidyr::unite(Key, CrossGroup,True_Class,Call, sep='_') %>% 
        dplyr::select(Key, Total_Count) %>% tidyr::spread(Key, Total_Count) %>% purrr::set_names(paste(names(.),'Total_Count', sep='_'))
    )
    
    # Write all Cross Valudation::
    cross_file_tib <- tibble::tibble()
    if (!is.null(dir)) {
      cross_file_tib <- saveModel(name=name, dir=dir, type='crossValidation',
                                  ms=cross_samp_tib, call=cross_call_tib, sums=cross_sums_tib,
                                  method=modName, alpha=alpha,
                                  verbose=verbose,vt=4,tc=tc+1,tt=tt)
    }
    
    # cross_min_perc <- cross_sums_tib %>% dplyr::filter(Call=='TP') %>% 
    #   dplyr::summarise(Min_Perc=min(Group_Perc)) %>% dplyr::pull(Min_Perc) %>% as.double()
    # if (is.null(cross_min_perc)) cross_min_perc <- 0.0
    # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Cross Validation Min={cross_min_perc}.{RET}{RET}"))

    cross_file_cnt <- base::nrow(cross_file_tib)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Cross Validation File Count={cross_file_cnt}.{RET}{RET}"))
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  cross_file_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Model Training Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

trainGlmnet = function(alpha, data, ss, cgns=NULL, class_idx, seed=NULL,
                       type.measure="class", family="binomial", parallel=TRUE,
                       
                       outName=NULL, dir=NULL,retFiles=FALSE,lambda=NULL,pars=NULL,
                       modName="glmnet",lambda_vec=c("lambda.min","lambda.1se"),
                       
                       crossTrain=FALSE,tests_ss=NULL,train_ss=NULL,cross_perc_min=NULL,
                       verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'trainGlmnet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  stime <- system.time({
    
    ## Set seed for reproducibility
    if (!is.null(seed)) set.seed(seed)
    if (!is.null(cgns)) data <- data[,cgns$Probe_ID]
    
    labs_idx_vec <- ss %>% dplyr::pull(!!class_idx) %>% as.vector()
    labs_idx_cnt <- labs_idx_vec %>% unique() %>% length()
    if (labs_idx_cnt>2) family <- "multinomial"
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} labs_idx_cnt={labs_idx_cnt}, family={family}.{RET}"))

    if (is.null(outName) ) {
      outName <- paste(modName, sep='-')
      if (!is.null(alpha)) outName <- paste(outName,alpha, sep='-')
      if (!is.null(lambda)) outName <- paste(outName,lambda, sep='-')
    }
    
    #
    # TBD:: 
    #    1. [Done]: Return files_tib below and not cross_min_perc
    #    2. [Done]: Re-calculate cross_min_perc from files() in files_tib
    #
    cross_file_tib <- NULL
    if (crossTrain) {
      cross_file_tib <- crossTrain(
        name=outName, alpha=alpha, data=data, ss_train=train_ss, ss_tests=tests_ss, 
        class_idx=class_idx, seed=seed,
        type.measure=type.measure, family=family, parallel=TRUE,
        dir=dir,lambda=lambda_vec[2],
        modName=modName,lambda_vec=lambda_vec,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      cross_min_perc <- NULL
      cross_min_perc <- 0.0
      cross_min_perc <- loadFromFileTib(tib=cross_file_tib, type="Summary", verbose=verbose,vt=4,tc=tc+1,tt=tt) %>%
        dplyr::filter(Call=='TP') %>%
        dplyr::summarise(Min_Perc=min(Group_Perc)) %>% dplyr::pull(Min_Perc) %>% as.double()
      if (is.null(cross_min_perc)) cross_min_perc <- 0.0
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Cross Validation Min={cross_min_perc} > cross_perc_min???{RET}{RET}"))
    }
    
    # Build full model::
    mod <- glmnet::cv.glmnet(x=data,y=labs_idx_vec,type.measure=type.measure,alpha=alpha,family=family,parallel=parallel)
    
    # Need to loop over all dgCMatrix coef for labs_idx_cnt > 2
    cgns <-NULL
    if (labs_idx_cnt>2) {
      if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Class Size={labs_idx_cnt}>2 Multi.{RET}") )
      cgns <- glmnetModToCoefN(mod, verbose=verbose, vt=vt+1, tc=tc+1)
    } else {
      if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Class Size={labs_idx_cnt} Pair.{RET}") )
      cgns <- glmnetModToCoef2(mod, verbose=verbose, vt=vt+1, tc=tc+1)
    }

    # Update Parameters::
    pars <- pars %>% tibble::add_row( Option="featureSize", Value=as.character( base::ncol(data) ), Type="Original" )
    for (lam_key in lambda_vec) {
      lam_key <- rlang::sym(lam_key)
      cnt_str <- cgns %>% dplyr::filter(!!lam_key != 0) %>% base::nrow() %>% as.character()
      pars <- pars %>% tibble::add_row( Option="featureSize", Value=cnt_str, Type=as.character(lam_key) )
    }
    pars <- pars %>% tibble::add_row( Option="model", Value=modName, Type="train" )
    if (!is.null(alpha)) pars <- pars %>% tibble::add_row( Option="alpha", Value=as.character(alpha), Type="train" )

    # Save Model::
    if (!is.null(dir)) {
      if (!is.null(cross_min_perc) && cross_min_perc < cross_perc_min) {
        if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Failed Cross Fold Validation ",
                                          "cross_min_perc={cross_min_perc} < {cross_perc_min}.{RET}") )
        mod <- NULL
      }
      
      file_tib <- saveModel(name=outName, dir=dir, type='model',
                            mod=mod, cgns=cgns, ss=ss, pars=pars,
                            method=modName, alpha=alpha,
                            verbose=verbose,vt=4,tc=tc+1,tt=tt) %>%
        dplyr::bind_rows(cross_file_tib)
    } else {
      retFiles=FALSE
    }
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  if (retFiles) return(file_tib)
  
  mod
}

trainRandomForest = function(alpha=NULL, data, ss, cgns, class_idx, seed=NULL, 
                             outName=NULL, dir=NULL,retFiles=FALSE,lambda=NULL, pars=NULL,
                             modName="rforest",
                             crossTrain=FALSE,tests_ss=NULL,train_ss=NULL,cross_perc_min=NULL,
                             verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'trainRandomForest'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  stime <- system.time({
    
    ## Set seed for reproducibility
    if (!is.null(seed)) set.seed(seed)
    if (!is.null(cgns)) data <- data[,cgns$Probe_ID]
    
    labs_idx_vec <- ss %>% dplyr::pull(!!class_idx) %>% as.vector()
    labs_idx_cnt <- labs_idx_vec %>% unique() %>% length()

    if (is.null(outName) ) {
      outName <- paste(modName, sep='-')
      if (!is.null(alpha)) outName <- paste(outName,alpha, sep='-')
      if (!is.null(lambda)) outName <- paste(outName,lambda, sep='-')
    }
    
    cross_min_perc <- NULL
    if (crossTrain) {
      # cross_min_perc <- crossTrain(
      #   name=outName, alpha=alpha, data=data, ss_train=train_ss, ss_tests=tests_ss, 
      #   class_idx=class_idx, seed=seed,
      #   dir=dir, modName=modName, 
      #   verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      cross_file_tib <- crossTrain(
        name=outName, alpha=alpha, data=data, ss_train=train_ss, ss_tests=tests_ss, 
        class_idx=class_idx, seed=seed,
        dir=dir, modName=modName, 
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
        # name=outName, alpha=alpha, data=data, ss_train=train_ss, ss_tests=tests_ss,
        # class_idx=class_idx, seed=seed,
        # type.measure=type.measure, family=family, parallel=TRUE,
        # dir=dir,lambda=lambda_vec[2],
        # modName=modName,lambda_vec=lambda_vec,
        # verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      cross_min_perc <- NULL
      cross_min_perc <- 0.0
      cross_min_perc <- loadFromFileTib(tib=cross_file_tib, type="Summary", verbose=verbose,vt=4,tc=tc+1,tt=tt) %>%
        dplyr::filter(Call=='TP') %>%
        dplyr::summarise(Min_Perc=min(Group_Perc)) %>% dplyr::pull(Min_Perc) %>% as.double()
      if (is.null(cross_min_perc)) cross_min_perc <- 0.0
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Cross Validation Min={cross_min_perc} > cross_perc_min???{RET}{RET}"))
    }

    # Convert Training::
    data_df  <- cbind(data, labs_idx_vec) %>% as.data.frame()
    df_ncols <- data_df %>% base::ncol()
    colnames(data_df)[df_ncols] <- 'Species'
    data_df$Species <- factor(data_df$Species)
    
    #  The first parameter specifies our formula: Species ~ . (we want to predict Species using each of the remaining columns of data).
    #  ntree defines the number of trees to be generated. It is typical to test a range of values for this parameter (i.e. 100,200,300,400,500) and choose the one that minimises the OOB estimate of error rate.
    #    - First attempt: df_ncols-1, might want to use default
    #
    #  mtry is the number of features used in the construction of each tree. These features are selected at random, which is where the “random” in “random forests” comes from. The default value for this parameter, when performing classification, is sqrt(number of features).
    #    - First attempt: using default
    #  importance enables the algorithm to calculate variable importance.
    #
    mod = randomForest(Species ~ ., data=data_df, ntree=df_ncols-1, importance=TRUE)
    
    # Example of reporting confusion matrix...
    #
    # mod$confusion %>% tibble::as_tibble(rownames='Class') %>% 
    #   purrr::set_names(c('Class_Reported', 'Class_0', 'Class_1', 'Class_Error')) %>% 
    #   tidyr::gather(Class_Predicted, 'Class_Count', -Class_Reported, -Class_Error)
    
    # Variables by Importance plot::
    #  TBD:: This needs to be captured and plotted in a png/pdf
    #
    # var_imp_plot <- varImpPlot(mod)
    
    # Update Parameters::
    pars <- pars %>% tibble::add_row( Option="featureSize", Value=as.character( base::ncol(data) ), Type="Original" )
    pars <- pars %>% tibble::add_row( Option="model", Value=modName, Type="train" )
    if (!is.null(alpha)) pars <- pars %>% tibble::add_row( Option="alpha", Value=as.character(alpha), Type="train" )
    if (!is.null(lambda)) pars <- pars %>% tibble::add_row( Option="lambda", Value=lambda, Type="train" )

    # Save Model::
    if (!is.null(dir)) {
      if (!is.null(cross_min_perc) && cross_min_perc < cross_perc_min) {
        if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Failed Cross Fold Validation ",
                                          "cross_min_perc={cross_min_perc} < {cross_perc_min}.{RET}") )
        mod <- NULL
      }
      
      file_tib <- saveModel(name=outName, dir=dir, type='model',
                            mod=mod, cgns=cgns, ss=ss, pars=pars,
                            method=modName, alpha=alpha,
                            verbose=verbose,vt=4,tc=tc+1,tt=tt) %>%
        dplyr::bind_rows(cross_file_tib)
    } else {
      retFiles=FALSE
    }
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  if (retFiles) return(file_tib)
  
  mod
}

trainXGBoost = function(data, labs, tData, tLabs, num_class, seed=NULL,
                        subsample_perc=0.75, objective_func="multi:softprob", eval_met="mlogloss",
                        booster_method="gbtree", eta=0.001, max_depth=5, gamma=3, colsample_bytree=1,
                        nrounds=100000, nthreads=4, early_stopping_rounds=20,
                        verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'trainXGBoost'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  stime <- system.time({
    
    ## Set seed for reproducibility
    if (!is.null(seed)) set.seed(seed)
    
    # Transform the two data sets into xgb.Matrix
    xgb_train = xgboost::xgb.DMatrix(data=data,  label=labs)
    xgb_test  = xgboost::xgb.DMatrix(data=tData, label=tLabs)
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Built Matricies.{RET}"))
    if (verbose>=vt+4) print(xgb_train)
    
    # num_class = boolLab %>% unique() %>% length()
    xgb_params = list(
      booster=booster_method,
      eta=eta,
      max_depth=max_depth,
      gamma=gamma,
      subsample=subsample_perc,
      colsample_bytree=colsample_bytree,
      objective=as.character(objective_func),
      eval_metric=as.character(eval_met),
      num_class=num_class
    )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Set Params; num_class={num_class}, nthreads={nthreads}, nrounds={nrounds}.{RET}"))
    if (verbose>=vt+4) print(xgb_params)
    
    # Train the XGBoost classifer
    xgb_verbose <- 0
    if (verbose>=vt+6) xgb_verbose = 1
    xgb_fit = xgboost::xgb.train(
      params=xgb_params,
      data=xgb_train,
      nrounds=nrounds,
      nthreads=nthreads,
      early_stopping_rounds=early_stopping_rounds,
      watchlist=list(val1=xgb_train, val2=xgb_test),
      verbose=xgb_verbose
    )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Trained Classifier.{RET}"))
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  xgb_fit
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Model Predicting Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

predGlmnet = function(mod, data, labs, name, lambda="lambda.1se", type='class',
                      verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'predGlmnet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  tib <- NULL
  stime <- system.time({
    # Different return values for glmnet: “link”, “response”, “coefficients”, “class”, “nonzero”
    
    labs_idx_cnt <- unique(labs) %>% length()
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} labs_idx_cnt={labs_idx_cnt}.{RET}"))

    if (labs_idx_cnt>2) {
      tib <- predict(object=mod, newx=data) %>%
        tibble::as_tibble(rownames='ample') %>%
        purrr::set_names(paste('S',stringr::str_remove(names(.), '.[0-9]+$'), sep='')) %>%
        dplyr::rename(Sentrix_Name=Sample) %>%
        dplyr::mutate(Pred_Class=as.integer( as.vector(predict(object=mod, newx=data, type="class")) ),
                      Pred_Class=as.factor(Pred_Class),
                      Model_Name=name)
    } else {
      tib <- predict(object=mod, newx=data, s=lambda, type=type) %>%
        tibble::as_tibble(rownames='Sentrix_Name') %>%
        purrr::set_names(c('Sentrix_Name','Pred_Class')) %>%
        dplyr::mutate(Pred_Class=as.factor(as.integer(Pred_Class)),
                      Model_Name=paste(lambda,name, sep='-'))
    }
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  tib
}

predRandomForest = function(mod, data, labs, name,
                            verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'predRandomForest'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  tib <- NULL
  stime <- system.time({
    
    labs_idx_cnt <- unique(labs) %>% length()
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} labs_idx_cnt={labs_idx_cnt}.{RET}"))
    
    pred_prob_mat  <- predict(mod, newdata=data, reshape=T, type="prob")
    pred_class_tib <- predict(mod, newdata=data, reshape=T, type="class") %>% 
      tibble::enframe(name="Sentrix_Name", value="Pred_Class_Raw") %>% 
      dplyr::mutate(Pred_Class_Raw=as.character(Pred_Class_Raw))

    if (labs_idx_cnt>2) {
      tib <- pred_prob_mat %>% as.data.frame() %>%
        tibble::as_tibble(rownames='ample') %>% 
        purrr::set_names(paste('S',stringr::str_remove(names(.), '.[0-9]+$'), sep='')) %>%
        dplyr::rename(Sentrix_Name=Sample) %>%
        dplyr::left_join(pred_class_tib, by="Sentrix_Name") %>%
        dplyr::rename(Pred_Class=Pred_Class_Raw) %>%
        dplyr::mutate(Pred_Class=as.factor( as.integer( Pred_Class ) ),
                      Model_Name=name)
    } else {
      pred_prob_tib <- pred_prob_mat %>% 
        as.data.frame() %>%
        tibble::as_tibble(rownames = "Sentrix_Name") %>% 
        purrr::set_names(c("Sentrix_Name", "Pred_Prob_0", "Pred_Prob_1") ) %>% 
        dplyr::mutate(Pred_Class=apply(pred_prob_mat,1,function(x) base::colnames(pred_prob_mat)[base::which.max(x)]),
                      Pred_Class=as.factor(as.integer(Pred_Class)))
      tib <- pred_prob_tib %>% dplyr::left_join(pred_class_tib, by="Sentrix_Name") %>% dplyr::mutate(Model_Name=name)
    }
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  tib
}

predXGBoost = function(mod, data, labs, name, 
                       verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'predXGBoost'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  tib <- NULL
  stime <- system.time({
    
    cur_pred <- predict(mod, newdata=data, reshape=T)
    colnames(cur_pred) <- labs
    
    tib <- tibble::tibble(
      Sentrix_Name=data %>% rownames(),
      Pred_Class = apply(cur_pred,1,function(x) base::colnames(cur_pred)[base::which.max(x)])
    ) %>% dplyr::mutate(Model_Name=name)
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Model Helper Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

glmnetModToCoef2 = function(m, lambda_vec=c("lambda.min","lambda.1se"), 
                            verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'glmnetModToCoef2'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  coef_tib <- NULL
  for (lambda_key in lambda_vec) {
    coef_db <- coef(m, s=lambda_key)[,1] %>% as.data.frame()
    cur_tib <- coef_db %>% tail(n=base::nrow(coef_db)-1) %>% 
      tibble::as_tibble(rownames="Probe_ID") %>% 
      purrr::set_names(c("Probe_ID", lambda_key))
    
    if (is.null(coef_tib)) {
      coef_tib <- cur_tib
    } else {
      coef_tib <- dplyr::full_join(coef_tib,cur_tib, by="Probe_ID")
    }
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  coef_tib
}

glmnetModToCoefN = function(m, lambda_vec=c("lambda.min","lambda.1se"), 
                            verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'glmnetModToCoefN'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  coef_tib <- tibble::tibble(Probe_ID=rownames(coef(m)[[1]])[-1])
  for (lambda_key in lambda_vec) {
    if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} lambda={lambda_key}.{RET}"))
    
    coef_list <- coef(m, s=lambda_key)
    coef_keys <- coef_list %>% names()
    
    cgns <- NULL
    for (coef_key in coef_keys) {
      cgns <- c( cgns, rownames(coef_list[[coef_key]])[ which(coef_list[[coef_key]] != 0) ][-1] )
    }
    
    lambda_key <- rlang::sym( lambda_key )
    cur_tib <- tibble::tibble(Probe_ID=unique(cgns), !!lambda_key := 1 ) %>% dplyr::arrange(Probe_ID)
    
    coef_tib <- dplyr::full_join(coef_tib,cur_tib, by="Probe_ID")
    # if (is.null(coef_tib)) { coef_tib <- cur_tib
    # } else { coef_tib <- dplyr::full_join(coef_tib,cur_tib, by="Probe_ID") }
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  coef_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Prediction Helper Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

callToSumTib = function(call, name_lab, true_lab="True_Class",call_lab="Call",
                        verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'callToSumTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  true_str <- true_lab %>% as.character()
  name_str <- name_lab %>% as.character()
  true_lab <- true_lab %>% as.character() %>% rlang::sym()
  name_lab <- name_lab %>% as.character() %>% rlang::sym()
  call_lab <- call_lab %>% as.character() %>% rlang::sym()
  
  sum_tib <- dplyr::left_join(
    dplyr::group_by(call, !!true_lab,!!call_lab) %>% dplyr::summarise(Group_Count=n()),
    dplyr::group_by(call, !!true_lab) %>% dplyr::summarise(Total_Count=n()),
    by=true_str) %>%
    dplyr::mutate(Group_Perc=round(100*Group_Count/Total_Count, 3)) %>%
    dplyr::ungroup()
  
  sum_tib <- sum_tib %>%
    dplyr::mutate(CrossGroup=!!name_str) %>%
    dplyr::select(CrossGroup, everything())
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  sum_tib
}

predToCalls = function(pred, labs, pred_lab="Pred_Class",
                       verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'predToCalls'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  pred_lab <- pred_lab %>% as.character() %>% rlang::sym()
  
  call_tib <- pred %>% dplyr::mutate(
    True_Class=labs,
    !!pred_lab := as.double(!!pred_lab)-1,
    Call=dplyr::case_when(
      !!pred_lab==True_Class ~ 'TP',
      !!pred_lab!=True_Class ~ 'FP',
      TRUE ~ NA_character_)
  )
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  call_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Rebranded Copied glmnet Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

makeX_glmnet_imputeNA=function(train, test=NULL, na.impute=FALSE,
                               verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'makeX_glmnet_imputeNA'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  df=train

  istest=!is.null(test)
  ntr=base::nrow(train)
  if(istest){
    nte=nrow(test)
    df=rbind(df,test)
  }
  x=prepareX_glmnet(df,sparse=FALSE)
  # if(na.impute){
  if (TRUE) {
    xbar=colMeans(x[seq(ntr),],na.rm=TRUE)
    x=na_replace_glmnet(x,xbar)
  }
  if(istest){
    xt=x[seq(ntr+1,ntr+nte),]
    x=x[seq(ntr),]
  }
  # if(na.impute)attr(x,"means")=xbar
  if(TRUE) attr(x,"means")=xbar
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  if(istest)list(x=x,xtest=xt) else x
}

prepareX_glmnet = function(df,sparse=FALSE,
                           verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'prepareX_glmnet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if(!inherits(df,"data.frame"))stop("first argument must be of class `data.frame`")
  whichfac=sapply(df,inherits,"factor")
  oldna=options()$na.action
  cna=as.character(substitute(na.action))
  options(na.action=na.pass)
  on.exit(options(na.action=oldna))
  if(any(whichfac))
    ctr=lapply(df[,whichfac], contrasts,contrast=FALSE)
  else ctr=NULL
  
  if(sparse){
    # m=sparse.model.matrix(~.-1,data=df,contrasts.arg=ctr,...)
    m=sparse.model.matrix(~.-1,data=df,contrasts.arg=ctr)
    m=na_sparse_fix_glmnet(m,names(df)) # sparse.model.matrix is faulty
  } else {
    # m = model.matrix(~.-1,data=df,contrasts.arg=ctr,...)
    m = model.matrix(~.-1,data=df,contrasts.arg=ctr)
  }
  if(any(whichfac))attr(m,"contrasts")=NULL
  attr(m,"assign")=NULL
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  m
}

na_sparse_fix_glmnet = function(x,dfnames) {
  a=attributes(x)
  ac=a$contrasts
  as=a$assign
  if(is.null(ac))return(x)
  acn=names(ac)
  whichn=match(acn,dfnames)
  for(i in whichn){
    xi=x[,as==i]
    rowtot=rowSums(xi)
    if(sum(rowtot)<length(rowtot)){# Nas present
      x[rowtot==0,as==i]=NA
    }
  }
  x
}

na_replace_glmnet=function(x,m=rowSums(x,na.rm=TRUE)){
  if(base::inherits(x,"sparseMatrix")){
    x=as(x,"CsparseMatrix")
    ccount=diff(x@p)
    cindex=rep(1:length(ccount),ccount)
    nas=is.na(x@x)
    if(any(nas))x@x[nas]=m[cindex[nas]]
  } else{
    d=dim(x)
    cols=rep(1:d[2],rep(d[1],d[2]))
    nas=is.na(x)
    if(any(nas))x[nas]=m[cols[nas]]
  }
  x
}

# End of file
