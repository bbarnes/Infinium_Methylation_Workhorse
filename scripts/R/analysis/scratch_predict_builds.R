

if (FALSE) {
  
  
  totCnt <- 875431
  newCnt <- 7903
  
  newCnt_100k <- 763
  newCnt_10k  <- 124
  newCnt_1k   <- 20
  
  newCnt/totCnt
  
  newCnt_1k/1000
  newCnt_10k/10000
  newCnt_100k/100000
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #               Latest:: July-13th: Full Builds from Cluster
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  tests_key <- 'COVIC-Set6-03062020'
  tests_key <- 'COVIC-Set5-10062020'
  pattern   <- paste0(tests_key,'.*.csv.gz$')
  
  # Sub set::
  # c0_csv <- list.files('/Users/bbarnes/Documents/CustomerFacing/results/C0', pattern=pattern, full.names=TRUE)
  # b4_csv <- list.files('/Users/bbarnes/Documents/CustomerFacing/results/B4', pattern=pattern, full.names=TRUE)
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #               Latest:: July-13th: Full Builds from Cluster:: C0
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  outDir <- '/Users/bbarnes/Documents/CustomerFacing/results'
  c0_sum_csv <- file.path(outDir, 'summary-c0.csv')
  b4_sum_csv <- file.path(outDir, 'summary-b4.csv')
  
  c0_csv <- list.files('/Users/bbarnes/Documents/CustomerFacing/results/latest/C0', pattern=pattern, full.names=TRUE, recursive=TRUE)
  c0_tib <- lapply(c0_csv, function(x) { readr::read_csv(x) %>% dplyr::mutate(User_lociPvalMin=as.double(User_lociPvalMin)) } ) %>% 
    dplyr::bind_rows() # %>% dplyr::distinct()
  
  c0_tib %>% base::nrow()
  c0_tib %>% dplyr::distinct() %>% base::nrow()
  
  c0_sum_tib <- c0_tib %>% dplyr::group_by(train_alpha, train_model, train_lambda,
                             User_lociBetaKey, User_lociPvalMin, # User_lociPvalKey, 
                             Original_featureSize) %>% #, User_runName, User_samplePvalName,Feature_featureSizeDml ) %>%
    purrr::set_names(stringr::str_remove(names(.),'User_') ) %>% 
    purrr::set_names(stringr::str_remove(names(.),'train_') ) %>%
    dplyr::rename(fSize=Original_featureSize, dmlSize=Feature_featureSizeDml, beta=lociBetaKey, pval=lociPvalMin ) %>% #, pvalKey=lociPvalKey) %>%
    dplyr::rename(
      NegFP=Known_5_tests_0_FP_Pred_Accuracy,
      NegTP=Known_5_tests_0_TP_Pred_Accuracy,
      PosTP=Known_5_tests_1_TP_Pred_Accuracy,
      PosFP=Known_5_tests_1_FP_Pred_Accuracy) %>%
    dplyr::summarise(dmlSize=as.integer(mean(dmlSize, na.rm=TRUE) ),
                     # NegTP.min=min(NegTP,na.rm=TRUE),NegTP.avg=mean(NegTP,na.rm=TRUE),NegTP.med=median(NegTP,na.rm=TRUE),NegTP.max=max(NegTP,na.rm=TRUE),
                     # PosTP.min=min(PosTP,na.rm=TRUE),PosTP.avg=mean(PosTP,na.rm=TRUE),PosTP.med=median(PosTP,na.rm=TRUE),PosTP.max=max(PosTP,na.rm=TRUE),
                     NegTP.min=min(NegTP,na.rm=TRUE),NegTP.avg=mean(NegTP,na.rm=TRUE),NegTP.max=max(NegTP,na.rm=TRUE),NegTP.sd=sd(NegTP,na.rm=TRUE),
                     PosTP.min=min(PosTP,na.rm=TRUE),PosTP.avg=mean(PosTP,na.rm=TRUE),PosTP.max=max(PosTP,na.rm=TRUE),PosTP.sd=sd(PosTP,na.rm=TRUE),
                     Count=n()) %>% dplyr::mutate(AvgSum=(NegTP.avg+PosTP.avg)/2 )
  
  c0_sum_tib %>% dplyr::filter(alpha==0.3 & model=="rforest" & lambda=="lambda-min" & beta=='ind_beta') %>% 
    dplyr::select(model,alpha,beta,pval,fSize,dmlSize,NegTP.avg,PosTP.avg,NegTP.sd,PosTP.sd,everything()) %>% # dplyr::arrange(model,alpha,pval)
    dplyr::filter(!is.na(PosTP.sd)) %>% as.data.frame()
  
  c0_sum_tib %>% dplyr::arrange(-AvgSum) %>% dplyr::filter(AvgSum>80)

  c0_out_tib <- c0_sum_tib %>% dplyr::select(model,alpha,beta,pval,fSize,dmlSize,NegTP.avg,PosTP.avg,NegTP.sd,PosTP.sd,everything()) %>% dplyr::arrange(model,alpha,pval)
  
  readr::write_csv(c0_out_tib,c0_sum_csv)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #               Latest:: July-13th: Full Builds from Cluster:: B4
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  b4_csv <- list.files('/Users/bbarnes/Documents/CustomerFacing/results/latest/B4', pattern=pattern, full.names=TRUE, recursive=TRUE)
  b4_tib <- lapply(b4_csv, function(x) { readr::read_csv(x) %>% dplyr::mutate(User_lociPvalMin=as.double(User_lociPvalMin)) } ) %>% 
    dplyr::bind_rows() # %>% dplyr::distinct()
  
  b4_tib %>% base::nrow()
  b4_tib %>% dplyr::distinct() %>% base::nrow()
  
  b4_sum_tib <- b4_tib %>% dplyr::group_by(train_alpha, train_model, train_lambda,
                                           User_lociBetaKey, User_lociPvalMin, # User_lociPvalKey, 
                                           Original_featureSize) %>% #, User_runName, User_samplePvalName,Feature_featureSizeDml ) %>%
    purrr::set_names(stringr::str_remove(names(.),'User_') ) %>% 
    purrr::set_names(stringr::str_remove(names(.),'train_') ) %>%
    dplyr::rename(fSzie=Original_featureSize, dmlSize=Feature_featureSizeDml, beta=lociBetaKey, pval=lociPvalMin ) %>% #, pvalKey=lociPvalKey) %>%
    dplyr::rename(
      NegFP=Known_5_tests_0_FP_Pred_Accuracy,
      NegTP=Known_5_tests_0_TP_Pred_Accuracy,
      PosTP=Known_5_tests_1_TP_Pred_Accuracy,
      PosFP=Known_5_tests_1_FP_Pred_Accuracy) %>%
    dplyr::summarise(dmlSize=as.integer(mean(dmlSize, na.rm=TRUE) ),
                     # NegTP.min=min(NegTP,na.rm=TRUE),NegTP.avg=mean(NegTP,na.rm=TRUE),NegTP.med=median(NegTP,na.rm=TRUE),NegTP.max=max(NegTP,na.rm=TRUE),
                     # PosTP.min=min(PosTP,na.rm=TRUE),PosTP.avg=mean(PosTP,na.rm=TRUE),PosTP.med=median(PosTP,na.rm=TRUE),PosTP.max=max(PosTP,na.rm=TRUE),
                     NegTP.min=min(NegTP,na.rm=TRUE),NegTP.avg=mean(NegTP,na.rm=TRUE),NegTP.max=max(NegTP,na.rm=TRUE),
                     PosTP.min=min(PosTP,na.rm=TRUE),PosTP.avg=mean(PosTP,na.rm=TRUE),PosTP.max=max(PosTP,na.rm=TRUE),
                     Count=n()) %>% dplyr::mutate(AvgSum=(NegTP.avg+PosTP.avg)/2 ) %>% dplyr::arrange(-AvgSum) %>% dplyr::filter(AvgSum>80)
  
  readr::write_csv(b4_sum_tib,b4_sum_csv)
  
  
  
  # For Set 6 ONLY::
  b4_tib %>% dplyr::group_by(train_alpha, train_model, train_lambda,
                             User_lociBetaKey, User_lociPvalMin, # User_lociPvalKey, 
                             Original_featureSize) %>% #, User_runName, User_samplePvalName,Feature_featureSizeDml ) %>%
    purrr::set_names(stringr::str_remove(names(.),'User_') ) %>% 
    purrr::set_names(stringr::str_remove(names(.),'train_') ) %>%
    dplyr::rename(fSzie=Original_featureSize, dmlSize=Feature_featureSizeDml, beta=lociBetaKey, pval=lociPvalMin ) %>% #, pvalKey=lociPvalKey) %>%
    dplyr::rename(
      NegFP=Known_NA_tests_0_FP_Pred_Accuracy,
      NegTP=Known_NA_tests_0_TP_Pred_Accuracy) %>%
    # PosTP=Known_NA_tests_1_TP_Pred_Accuracy,
    # PosFP=Known_NA_tests_1_FP_Pred_Accuracy) %>%
    # dplyr::rename(NegFP=Known_5_tests_0_FP_Pred_Accuracy,NegTP=Known_5_tests_0_TP_Pred_Accuracy,
    #               PosTP=Known_5_tests_1_TP_Pred_Accuracy,PosFP=Known_5_tests_1_FP_Pred_Accuracy) %>%
    dplyr::summarise(dmlSize=as.integer(mean(dmlSize, na.rm=TRUE) ),
                     # NegTP.min=min(NegTP,na.rm=TRUE),NegTP.avg=mean(NegTP,na.rm=TRUE),NegTP.med=median(NegTP,na.rm=TRUE),NegTP.max=max(NegTP,na.rm=TRUE),
                     # PosTP.min=min(PosTP,na.rm=TRUE),PosTP.avg=mean(PosTP,na.rm=TRUE),PosTP.med=median(PosTP,na.rm=TRUE),PosTP.max=max(PosTP,na.rm=TRUE),
                     NegTP.min=min(NegTP,na.rm=TRUE),NegTP.avg=mean(NegTP,na.rm=TRUE),NegTP.max=max(NegTP,na.rm=TRUE),
                     # PosTP.min=min(PosTP,na.rm=TRUE),PosTP.avg=mean(PosTP,na.rm=TRUE),PosTP.max=max(PosTP,na.rm=TRUE),
                     Count=n()) %>% dplyr::mutate(AvgSum=(NegTP.avg) ) %>% dplyr::arrange(-AvgSum) %>% 
    dplyr::filter(AvgSum>80)
  
  
  
  b4_tib %>% dplyr::filter(train_alpha==0.1 & train_model=='rforest' & train_lambda=='lambda-min' & 
                             User_lociBetaKey=='ind_beta', User_lociPvalMin==1) %>% as.data.frame()
  
  b4_tib %>% dplyr::filter(train_alpha==0.7 & train_model=='rforest' & train_lambda=='lambda-min' & 
                             User_lociBetaKey=='ind_beta', User_lociPvalMin==1) %>% as.data.frame()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #               Latest:: July-12th: Full Builds from Cluster
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  all_line_csv <- list.files('/Users/bbarnes/Documents/CustomerFacing/results/july/full-Ivana', pattern=".csv.gz$", full.names=TRUE)
  all_line_csv <- list.files('/Users/bbarnes/Documents/CustomerFacing/results/july/full', pattern=".csv.gz$", full.names=TRUE)
  
  all_line_tib <- lapply(all_line_csv, function(x) { readr::read_csv(x) %>% dplyr::mutate(User_lociPvalMin=as.double(User_lociPvalMin)) } ) %>% 
    dplyr::bind_rows()
  
  start_str <- 'Known_5_tests'
  all_line_tib %>% dplyr::select(1:17, starts_with(start_str)) %>% 
    dplyr::select(1:17,ends_with("Accuracy")) %>% 
    dplyr::select(1:17,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    # purrr::set_names( new_cols ) %>% 
    mutate_all(~replace(., is.na(.), 0)) %>%
    # dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    # dplyr::arrange(-TP_Sum,-PosTP,-NegTP) %>%
    dplyr::select(c(2:8,16:20)) %>%
    dplyr::rename(NegFP=Known_5_tests_0_FP_Pred_Accuracy,NegTP=Known_5_tests_0_TP_Pred_Accuracy,
                  PosTP=Known_5_tests_1_TP_Pred_Accuracy,PosFP=Known_5_tests_1_FP_Pred_Accuracy) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::select(NegTP,NegFP, PosTP,PosFP, everything()) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP) %>% dplyr::select(-Known_5_tests_0_FP_Pred_Mode_Count) %>% print(n=528)
  
  
  
  
  new_cols <- c(all_line_tib %>% dplyr::select(1:17) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  all_line_tib %>% dplyr::select(1:17, starts_with(start_str)) %>% 
    dplyr::select(1:17,ends_with("Accuracy")) %>% 
    dplyr::select(1:17,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    # purrr::set_names( new_cols ) %>% 
    mutate_all(~replace(., is.na(.), 0)) %>%
    # dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    # dplyr::arrange(-TP_Sum,-PosTP,-NegTP) %>%
    dplyr::select(c(4,6:11,16:20))
  # dplyr::select(c(5,7,8,11,18:24))
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Latest:: July-12th
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # 17 intro columns::
  #  all_line_tib %>% dplyr::select(Feature_featureNamePre:User_version)
  
  all_line_csv <- '/Users/bbarnes/Documents/CustomerFacing/predict_builds/Sample_Class/COVIC-Set1-15052020.SamplePrediction.sheet.csv.gz'
  all_line_tib <- read_csv(all_line_csv)
  
  all_line_dir <- '/Users/bbarnes/Documents/CustomerFacing/predict_builds/Sample_Class/COVIC-Set1-15052020'
  
  all_line_csv <- list.files('/Users/bbarnes/Documents/CustomerFacing/predict_builds/Sample_Class/COVIC-Set1-15052020', pattern='SamplePrediction.sheet.csv.gz$', recursive=TRUE, full.names=TRUE)
  # all_line_tib <- lapply(all_line_csv, readr::read_csv) %>% dplyr::bind_rows()
  all_line_tib <- lapply(all_line_csv, function(x) { readr::read_csv(x) %>% dplyr::mutate(User_lociPvalMin=as.double(User_lociPvalMin)) } ) %>% 
    dplyr::bind_rows()
  
  start_str <- 'Known_5_tests'
  new_cols <- c(all_line_tib %>% dplyr::select(1:17) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  all_line_tib %>% dplyr::select(1:17, starts_with(start_str)) %>% 
    dplyr::select(1:17,ends_with("Accuracy")) %>% 
    dplyr::select(1:17,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP) %>%
    dplyr::select(c(4,6:11,18:24))
    # dplyr::select(c(5,7,8,11,18:24))
  

  start_str <- 'Known_1_train'
  new_cols <- c(all_line_tib %>% dplyr::select(1:17) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  all_line_tib %>% dplyr::select(1:17, starts_with(start_str)) %>% 
    dplyr::select(1:17,ends_with("Accuracy")) %>% 
    dplyr::select(1:17,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
    # dplyr::select(ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% 
    mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP) %>%
    dplyr::select(c(5,7,8,11,18:24))
  
  
  
  
  
  
  
  
  
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Latest:: July-9th
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  resDir <- file.path(par$topDir, 'results/june')
  # dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set1-15052020.Test-1-5.SamplePrediction.sheet.csv.gz") )
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set1-15052020.Test-2.SamplePrediction.sheet.csv.gz") )
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set1-15052020.Test-3.SamplePrediction.sheet.csv.gz") )
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set1-15052020.Test-4.SamplePrediction.sheet.csv.gz") )
  
  # dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set5-10062020.Test-1-5.SamplePrediction.sheet.csv.gz") )
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set5-10062020.Test-2.SamplePrediction.sheet.csv.gz") )
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set5-10062020.Test-3.SamplePrediction.sheet.csv.gz") )
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set5-10062020.Test-4.SamplePrediction.sheet.csv.gz") )
  
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Case by Case::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set1-15052020.Test-1-5.SamplePrediction.sheet.csv.gz") )
  
  start_str <- 'Known_5_tests'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  start_str <- 'Known_1_train'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Case by Case::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set1-15052020.Test-2.SamplePrediction.sheet.csv.gz") )
  
  start_str <- 'Known_2_tests'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  start_str <- 'Known_1_train'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Case by Case::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set1-15052020.Test-3.SamplePrediction.sheet.csv.gz") )
  
  start_str <- 'Known_3_tests'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  start_str <- 'Known_1_train'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Case by Case::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set1-15052020.Test-4.SamplePrediction.sheet.csv.gz") )
  
  start_str <- 'Known_4_tests'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  start_str <- 'Known_1_train'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  
  
  
  
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Case by Case::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set5-10062020.Test-1-5.SamplePrediction.sheet.csv.gz") )
  
  start_str <- 'Known_1_tests'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  start_str <- 'Known_5_train'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Case by Case::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set5-10062020.Test-2.SamplePrediction.sheet.csv.gz") )
  
  start_str <- 'Known_2_tests'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  start_str <- 'Known_5_train'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Case by Case::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set5-10062020.Test-3.SamplePrediction.sheet.csv.gz") )
  
  start_str <- 'Known_3_tests'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  start_str <- 'Known_5_train'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Case by Case::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set5-10062020.Test-4.SamplePrediction.sheet.csv.gz") )
  
  start_str <- 'Known_4_tests'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  start_str <- 'Known_4_train'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP","PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Case by Case::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  dat_tib <- readr::read_csv(file.path(resDir, "COVIC-Set5-10062020.Test-1-5.SamplePrediction.sheet.csv.gz") )
  
  start_str <- 'Known_1_tests'
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  start_str <- 'Known_5_train'
  new_cols <- c(dat_tib %>% dplyr::select(1:7) %>% names(), "NegTP","NegFP", "PosTP") #,"PosFP")
  dat_tib %>% dplyr::select(1:7, starts_with(start_str)) %>% 
    dplyr::select(1:7,ends_with("Accuracy")) %>% 
    dplyr::select(1:7,ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names( new_cols ) %>% mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate(PosFP=0,NegSum=NegTP+NegFP, PosSum=PosTP+PosFP, TP_Sum=(NegTP+PosTP)/2 ) %>%
    dplyr::arrange(-TP_Sum,-PosTP,-NegTP)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  dat_tib %>% dplyr::select(starts_with("Known_5_tests")) %>% dplyr::select(ends_with("Accuracy")) %>% 
    dplyr::select(ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    dplyr::arrange(Known_5_tests_1_FP_Pred_Accuracy) %>% 
    dplyr::mutate(S1=Known_5_tests_0_TP_Pred_Accuracy+Known_5_tests_0_FP_Pred_Accuracy, 
                  S2=Known_5_tests_1_FP_Pred_Accuracy+Known_5_tests_1_TP_Pred_Accuracy) %>% as.data.frame()
  
  
  dat_tib %>% dplyr::select(starts_with("Known_1_train")) %>% 
    dplyr::select(ends_with("Accuracy")) %>% 
    dplyr::select(ends_with("0_TP_Pred_Accuracy"), ends_with("0_FP_Pred_Accuracy"),
                  ends_with("1_TP_Pred_Accuracy"), ends_with("1_FP_Pred_Accuracy")) %>%
    purrr::set_names(c("NegTP","NegFP", "PosTP","PosFP") )
  dplyr::arrange(1,3) %>% 
    dplyr::mutate(S1=.[1]+.[2])
  
  dplyr::arrange(Known_5_tests_1_FP_Pred_Accuracy) %>% 
    dplyr::mutate(S1=Known_5_tests_0_TP_Pred_Accuracy+Known_5_tests_0_FP_Pred_Accuracy, 
                  S2=Known_5_tests_1_FP_Pred_Accuracy+Known_5_tests_1_TP_Pred_Accuracy) %>% as.data.frame()
  
  
  
  ss_csv <- '/Users/bbarnes/Documents/CustomerFacing/results/june/COVIC-Set1-15052020.Test-1-5.SamplePrediction.sheet.csv.gz'
  ss_tib <- readr::read_csv(ss_csv)
  
  ss_tib %>% dplyr::arrange(Known_1_train_0_TP_Pred_Accuracy) %>% as.data.frame()
  ss_tib %>% dplyr::arrange(Known_1_tests_1_FP_Pred_Accuracy) %>% as.data.frame()
  
  
  ss_tib %>% dplyr::select(ends_with("Accuracy")) %>% names()
  
  ss_tib %>% dplyr::select(starts_with("Known_1_train")) %>% dplyr::select(ends_with("Accuracy")) %>% names()
  
  ss_tib %>% dplyr::select(starts_with("Known_5_tests")) %>% dplyr::select(ends_with("Accuracy"))
  
  
  
  ss_tib %>% dplyr::select(starts_with("Known_5_tests")) %>% dplyr::select(ends_with("Accuracy")) %>% dplyr::arrange(Known_5_tests_1_FP_Pred_Accuracy) %>% dplyr::mutate(S1=Known_5_tests_0_TP_Pred_Accuracy+Known_5_tests_0_FP_Pred_Accuracy, S2=Known_5_tests_1_FP_Pred_Accuracy+Known_5_tests_1_TP_Pred_Accuracy) %>% as.data.frame()
  

  
  
  
  
  
  
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Really-Really Old Testing::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  pred_files <- list.files(cur_opt_dir, pattern='SamplePrediction.sheet.csv.gz$', recursive=TRUE, full.names=TRUE)
  
  pred_dir <- file.path(par$topDir, 'predict_builds/COVIC-Set5-10062020/i-beta_i-poob-1')
  pred_dir <- file.path(par$topDir, 'results/fromCluster/predict_builds/COVIC-Set1-15052020/i-beta_i-poob-1')
  
  pred_files <- list.files(pred_dir, pattern='SamplePrediction.sheet.csv.gz$', recursive=TRUE, full.names=TRUE)
  
  all_pred_tib <- lapply(pred_files, readr::read_csv ) %>% dplyr::bind_rows()
  
  all_pred_tib %>% dplyr::arrange(Known_tests_0_FP_Pred_Accuracy) %>% as.data.frame()
  all_pred_tib %>% dplyr::arrange(Known_tests_1_FP_Pred_Accuracy) %>% as.data.frame()
  
  #
  # Karyotype_1_call Predictions::
  #
  pred_files <- list.files('/Users/bbarnes/Documents/CustomerFacing/predict_builds/COVIC-Set1-5', pattern='.SamplePrediction.sheet.csv.gz$', recursive=TRUE, full.names=TRUE)
  all_pred_tib <- lapply(pred_files, readr::read_csv ) %>% dplyr::bind_rows()
  
  all_pred_tib %>% dplyr::arrange(Known_tests_0_TP_Pred_Accuracy) %>% as.data.frame()
  all_pred_tib %>% dplyr::arrange(Known_tests_1_FP_Pred_Accuracy) %>% as.data.frame()
  
  #
  # From Cluster::
  #
  cov1_csv <- '/Users/bbarnes/Documents/CustomerFacing/results/july/COVIC-Set1-15052020.tib.csv.gz'
  cov1_tib <- readr::read_csv(cov1_csv)
  
  cov1_tib %>% dplyr::arrange(Known_tests_0_FP_Pred_Accuracy) %>% as.data.frame()
  cov1_tib %>% dplyr::arrange(Known_tests_1_FP_Pred_Accuracy) %>% as.data.frame()
  
  
  cov5_csv <- '/Users/bbarnes/Documents/CustomerFacing/results/july/COVIC-Set5-10062020.tib.csv.gz'
  cov5_tib <- readr::read_csv(cov1_csv)
  
  
  cov5_tib %>% dplyr::arrange(Known_tests_0_FP_Pred_Accuracy) %>% as.data.frame()
  cov5_tib %>% dplyr::arrange(Known_tests_1_FP_Pred_Accuracy) %>% as.data.frame()
  
  
  
  # Full Auto Sample Sheet Summary::
  
  ss_pre_tib <- readr::read_csv('/Users/bbarnes/Documents/CustomerFacing/transfer_AutoSampleSheet.csv.gz')
  ss_all_tib <- readr::read_csv('/Users/bbarnes/Documents/CustomerFacing/transfer.07072020_AutoSampleSheet.csv.gz')
  
}


# End of file
