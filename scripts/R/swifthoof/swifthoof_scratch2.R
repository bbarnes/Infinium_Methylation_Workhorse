
if (FALSE) {
  rdat1 <- rdat
  
  rdat$ssheet_tib %>% dplyr::select(
    cg_calls_mins_count_1,cg_calls_pass_count_1,cg_calls_pass_perc_1,cg_calls_total_count_1
  )
  
  # rdat$ssheet_tib$cg_pvals_pOOBAH_pass_perc_2
  # rdat$ssheet_tib$cg_pvals_pOOBAH_pass_perc_1
  # 
  # rdat$ssheet_tib$cg_2_pass_perc_basic_0
  # rdat$ssheet_tib$cg_2_pass_cnt_basic_0
  # 
  # rdat$ssheet_tib$cg_pass_perc_basic_0
  # rdat$ssheet_tib$cg_1_pass_perc_basic_0
  # 
  # rdat$cur_list$sums_dat
  # 
  # 
  rdat$ssheet_tib %>% 
    dplyr::select(Sentrix_Name,
                  cg_pass_perc_basic_0,
                  cg_1_pass_cnt_basic_0,cg_1_pass_perc_basic_0,
                  cg_2_pass_cnt_basic_0,cg_2_pass_perc_basic_0,
                  cg_pvals_pOOBAH_pass_perc_1,cg_pvals_pOOBAH_pass_perc_2
    )
  
  
  # cg_calls_mins_count_1 cg_calls_pass_count_1 cg_calls_pass_perc_1 cg_calls_total_count_1
  # <dbl>                 <dbl>                <dbl>                  <dbl>
  #   1                799276                838830                 97.2                 862927

  # > rdat$ssheet_tib$cg_pass_perc_basic_0
  # [1] 93.029
  
  #
  # callToPassPerc
  
  call_dirs <- '/Users/bretbarnes/Documents/data/VA_MVP/transfer/data.v.4.3'
  call_files <- list.files(call_dirs, pattern='raw.call.dat.csv.gz$', full.names=TRUE)
  
  for (cidx in c(1:length(call_files))) {
    
  }
}





if (FALSE) {
  rdat2 <- rdat
  
  cluster_ssh_csv <- '/Users/bretbarnes/Documents/data/VA_MVP/transfer/docker.v.1.9.AutoSampleSheet.csv.gz'
  cluster_ssh_tib <- readr::read_csv(cluster_ssh_csv)
  
  va43_ss_dir <- '/Users/bretbarnes/Documents/data/VA_MVP/transfer/data'
  va43_ss_csv <- file.path(va43_ss_dir,'docker.v.4.3.AutoSampleSheet.csv.gz')
  va43_ss_tib <- readr::read_csv(va43_ss_csv)
  
  sentrix_tar <- '203962710025_R08C01'
  sentrix_tar <- '203962710025_R01C01'
  sentrix_tar <- '203962710025_R02C01'
  var43_call_csv <- file.path(va43_ss_dir, paste(sentrix_tar, 'EPIC_B4_raw.call.dat.csv.gz', sep='_'))
  var43_call_ssh <- callToPassPerc(file =var43_call_csv, 
                                   key="pvals_pOOBAH", name="v.4.6", min=cur_min_pval)

  dir19_mlk_dir <- '/Users/bretbarnes/Documents/scratch/mlk2/swifthoof_main/CNTL-Samples_VendA_10092020'
  dir19_dsh_csv <- file.path(dir19_mlk_dir, paste(sentrix_tar,'EPIC_B4_AutoSampleSheetDescriptionTable.csv.gz', sep='_'))
  dir19_dsh_tib <- suppressMessages(suppressWarnings( readr::read_csv(dir19_dsh_csv) ))
  dir19_ssh_csv <- file.path(dir19_mlk_dir, paste(sentrix_tar,'EPIC_B4_AutoSampleSheet.csv.gz', sep='_'))
  dir19_ssh_tib <- suppressMessages(suppressWarnings( readr::read_csv(dir19_ssh_csv) ))
  
  dir19_ssh_tib %>% dplyr::select( Sentrix_Name, cg_pvals_pOOBAH_pass_perc_1, cg_pvals_pOOBAH_pass_perc_2 ) %>% print()
  var43_call_ssh %>% print()
  
  # This needs to be true:: from rdat$ssheet_tib
  #  cg_calls_pass_perc_1 == cg_pvals_pOOBAH_pass_perc_1
  #
  
  rdat$cur_list$sums_dat %>% 
    dplyr::select( dplyr::all_of(c("Workflow_key", "Workflow_idx", "Probe_Type", "full_pass_perc") ) ) %>% 
    dplyr::filter(!is.na(full_pass_perc))
  
  del='_'
  rdat$cur_list$sums_dat %>%
    dplyr::filter(Workflow_key=='raw') %>%
    tidyr::unite(key, Probe_Type,Metric, sep=del) %>% 
    tidyr::gather(metric, value, -key, -Workflow_key, -Workflow_idx, -Probe_Design) %>% 
    dplyr::filter(!is.na(value)) %>% 
    dplyr::filter(metric %in% c('full_pass_perc')) %>%
    dplyr::mutate(metric='pass_perc') %>%
    tidyr::unite(key, key,metric, sep=del) %>% 
    dplyr::select(-Workflow_key, -Workflow_idx, -Probe_Design) %>%
    dplyr::distinct() %>%
    tidyr::spread(key, value) %>%
    dplyr::select(dplyr::contains('_pass_perc_'), 
                  dplyr::everything())
  
  # TBD::
  # Check 1:
  #  - Compare beta and pvals vs. extra slots in ssetToSummary()
  #
  # Check 2:
  #  - Functionize callToPassPerc()
  #  - Auto Check Detection P-values against calls file
  #
  # Check 3:
  #  - compare methods auto detect...
  #
  callToPassPerc = function(tib, key, name, min=0.1, perc=90) {
    key_sys <- rlang::sym(key)
    tot_cnt <- tib %>% dplyr::filter(stringr::str_starts(Probe_ID, 'cg')) %>% base::nrow()
    pas_cnt <- tib %>% dplyr::filter(stringr::str_starts(Probe_ID, 'cg')) %>% 
      dplyr::filter(!!key_sys <= min) %>% base::nrow()
    pas_per  <- round(100*pas_cnt/tot_cnt, 3)
    cat(glue::glue("per={pas_per}, pas={pas_cnt}, tot={tot_cnt}{RET}"))
    
    ret_tib <- tibble(
      name=name,
      pas_per=pas_per,
      pas_cnt=pas_cnt,
      tot_cnt=tot_cnt,
      cut_off=min
    )
    
    ret_tib
  }
  
  cur_min_pval <- min_pval_vec[1]
  cur_min_perc <- min_perc_vec[1]
  
  va_full_ss_csv <- '/Users/bretbarnes/Documents/data/VA-CNTL-AB-mergedAutoSampleSheet_v.4.3.csv.gz'
  va_full_ss_tib <- readr::read_csv(va_full_ss_csv)
  
  va_full_ss_tib %>% dplyr::select(Sentrix_Name,
                                   Poob_Pass_1_Perc,
                                   Poob_Pass_2_Perc,
                                   cg_1_pvals_pOOBAH_pass_perc,
                                   cg_2_pvals_pOOBAH_pass_perc) %>% 
    dplyr::filter(Sentrix_Name==rdat$ssheet_tib$Sentrix_Name)
  
  # Prev Calls::
  # prev_ssh1_csv <- '/Users/bretbarnes/Documents/scratch/swifthoof_main.jan6/CNTL-Samples_VendA_10092020/203962710025_R08C01_EPIC_B4_AutoSampleSheet.csv.gz'
  # prev_ssh1_tib <- readr::read_csv(prev_ssh1_csv)
  
  # Prev Calls::
  pr46_call_dir <- '/Users/bretbarnes/Documents/scratch/swifthoof_main.jan6/CNTL-Samples_VendA_10092020'
  pr46_call_csv <- file.path(pr46_call_dir, paste(rdat$ssheet_tib$Sentrix_Name, 'EPIC_B4_raw.call.dat.csv.gz', sep='_'))
  pr46_call_tib <- readr::read_csv(pr46_call_csv)
  pr46_sum_tib  <- callToPassPerc(tib=pr46_call_tib, key="pvals_pOOBAH", name="v.4.6", min=cur_min_pval)
  
  pr43_call_dir <- '/Users/bretbarnes/Documents/tools/bk/docker-repo/Infinium_Methylation_Workhorse.v.4.3/scratch/swifthoof/swifthoof_main/_v.4.3'
  pr43_ssh1_csv <- file.path(pr43_call_dir, paste(rdat$ssheet_tib$Sentrix_Name, 'EPIC_B4_AutoSampleSheet.csv.gz', sep='_'))
  pr43_call_csv <- file.path(pr43_call_dir, paste(rdat$ssheet_tib$Sentrix_Name, 'EPIC_B4_raw.call.dat.csv.gz', sep='_'))
  pr43_ssh1_tib <- readr::read_csv(pr43_ssh1_csv)
  pr43_call_tib <- readr::read_csv(pr43_call_csv)
  pr43_sum_tib  <- callToPassPerc(tib=pr43_call_tib, key="pvals_pOOBAH", name="v.4.3", min=cur_min_pval)
  
  # Cur Calls::
  curr_call_csv <- file.path(opt$outDir, paste(rdat$ssheet_tib$Sentrix_Name, 'EPIC_B4_raw.call.dat.csv.gz', sep='_'))
  curr_call_tib <- readr::read_csv(curr_call_csv)
  curr_sum_tib  <- callToPassPerc(tib=curr_call_tib, key="pvals_pOOBAH", name="v.1.8", min=cur_min_pval)
  
  # Open Calls::
  open_call_tib <- rdat$open_sset_dat@pval %>% tibble::enframe(name="Probe_ID", value="pvals_pOOBAH")
  open_sum_tib  <- callToPassPerc(tib=open_call_tib, key="pvals_pOOBAH", name="v.1.8", min=cur_min_pval)
  
  dplyr::bind_rows(pr46_sum_tib,
                   pr43_sum_tib,
                   curr_sum_tib,
                   open_sum_tib) %>% print()
  
  
  open_ssh1_tib <- ssetToPassPercSsheet(sset = rdat$open_sset_dat, min=cur_min_pval, per=cur_min_perc,
                                        idx=0, type='cg', verbose=1)
  
  open_ssh1_tib %>% print()
  rdat$open_sum1_ssh %>% print()
  rdat$open_sum2_ssh %>% print()
  
  callToPassPerc(file = pr46_call_csv, key = "pvals_pOOBAH", name = "v.4.6", min = 0.1, type = "cg", idx = 0, verbose = opt$verbose)
  callToPassPerc(tib = pr46_call_tib, key = "pvals_pOOBAH", name = "v.4.6", min = 0.1, type = "cg", idx = 0, verbose = opt$verbose)
  
  #
  # Compare Beta Values::
  #
  sample <- rdat$ssheet_tib$AutoSample_R2_Key_2
  sample <- rdat$ssheet_tib$AutoSample_dB_Key_2
  
  auto_beta_tib <- auto_sam_tib %>% dplyr::select(dplyr::all_of(c("Probe_ID",sample))) %>% purrr::set_names(c("Probe_ID", "betas"))
  
  r2_2tibs = function(a, b) {
    a <- a %>% dplyr::select(Probe_ID, dplyr::contains("betas")) %>% dplyr::select(1,2) %>% purrr::set_names(c("Probe_ID", "betas"))
    b <- b %>% dplyr::select(Probe_ID, dplyr::contains("betas")) %>% dplyr::select(1,2) %>% purrr::set_names(c("Probe_ID", "betas"))
    
    # print(a)
    # print(b)
    
    ab_tib <- dplyr::inner_join(a,b, by="Probe_ID", suffix=c("_ref","_can")) %>%
      dplyr::filter(stringr::str_starts(Probe_ID,'cg')) %>% 
      # print(ab_tib)
      
      ab_mat <- ab_tib %>%
      tibble::column_to_rownames(var="Probe_ID") %>% 
      as.matrix() 
    # ab_mat %>% head() %>% print()
    
    r2_tib <- cor(ab_mat, method="pearson", use="pairwise.complete.obs") %>% as_tibble()
    # print(r2_tib)
    
    r2_val <- r2_tib %>% head(n=1) %>% dplyr::pull(2)
    
    r2_val
  }
  pr46_r2_val <- r2_2tibs(auto_beta_tib, pr46_call_tib)
  pr43_r2_val <- r2_2tibs(auto_beta_tib, pr43_call_tib)
  curr_r2_val <- r2_2tibs(auto_beta_tib, curr_call_tib)
  
  rdat$ssheet_tib$AutoSample_dB_Val_1
  rdat$ssheet_tib$AutoSample_R2_Val_1        
  
  dB_2tibs = function(a,b) {
    a <- a %>% dplyr::select(Probe_ID, dplyr::contains("betas")) %>% dplyr::select(1,2) %>% purrr::set_names(c("Probe_ID", "betas"))
    b <- b %>% dplyr::select(Probe_ID, dplyr::contains("betas")) %>% dplyr::select(1,2) %>% purrr::set_names(c("Probe_ID", "betas"))
    
    dB_val <- 
      dplyr::inner_join(a,b, by="Probe_ID", suffix=c("_ref","_can")) %>%
      dplyr::filter(stringr::str_starts(Probe_ID,'cg')) %>% 
      dplyr::mutate(dB=base::abs(betas_ref-betas_can)) %>% 
      dplyr::summarise(pass_perc=cntPer_lte(dB,opt$minDeltaBeta)) %>%
      head(n=1) %>% dplyr::pull(1)
    
    dB_val
  }
  pr46_dB_val <- dB_2tibs(auto_beta_tib, pr46_call_tib)
  pr43_dB_val <- dB_2tibs(auto_beta_tib, pr43_call_tib)
  curr_dB_val <- dB_2tibs(auto_beta_tib, curr_call_tib)
  
  
  
  # Get their difference::
  mis1_cnt <- cur_call_tib %>% dplyr::anti_join(open_sset_tib, by="Probe_ID") %>% dplyr::filter(stringr::str_starts(Probe_ID, 'cg')) %>% base::nrow()
  mis2_cnt <- open_sset_tib %>% dplyr::anti_join(cur_call_tib, by="Probe_ID") %>% dplyr::filter(stringr::str_starts(Probe_ID, 'cg')) %>% base::nrow()
  
  prev_mis1_cnt <- cur_call_tib %>% dplyr::anti_join(prev_call_tib, by="Probe_ID") %>% dplyr::filter(stringr::str_starts(Probe_ID, 'cg')) %>% base::nrow()
  prev_mis2_cnt <- prev_call_tib %>% dplyr::anti_join(cur_call_tib, by="Probe_ID") %>% dplyr::filter(stringr::str_starts(Probe_ID, 'cg')) %>% base::nrow()
  
  
  prev_join_call_tib <- prev_call_tib %>% 
    dplyr::inner_join(cur_call_tib, by="Probe_ID", suffix=c("_Old","_New")) %>% 
    dplyr::mutate(Diff_Poob=pvals_pOOBAH_Old-pvals_pOOBAH_New)
  
  curr_join_call_tib <- cur_call_tib %>% dplyr::inner_join(open_sset_tib, by="Probe_ID", suffix=c("_Old","_New")) %>% 
    dplyr::mutate(Diff_Poob=pvals_pOOBAH-Pval)
  
  prev_join_call_tib %>% dplyr::filter(Diff_Poob>0.01) %>% base::nrow()
  curr_join_call_tib %>% dplyr::filter(Diff_Poob>0.01) %>% base::nrow()
  
  
  
  cur_min_pval <- 0.2
  cur_min_perc <- 90
  new2_tib <- ssetToPassPercSsheet(sset = rdat$open_sset_dat, min=cur_min_pval, per=cur_min_perc,
                                   idx=0, type='cg', verbose=10)
  
  cur_min_pval <- 0.05
  cur_min_perc <- 90
  new05_tib <- ssetToPassPercSsheet(sset = rdat$open_sset_dat, min=cur_min_pval, per=cur_min_perc,
                                    idx=0, type='cg', verbose=10)
  
  
  dplyr::inner_join(
    tidyr::gather(rdat$ssheet_tib, Variable, Value),
    tibble::enframe( unlist( sapply(rdat$ssheet_tib, class) ), name="Variable", value="Data_Type" ),
    by="Variable"
  )
  
  ano_ssh1_tib <- getSsheetCoreAnnoTib(minOobPval = 0.1, minOobPerc = 90, minNegPval = 0.02, minNegPerc = 98, minDb = 0.2, verbose = 10)
  ano_ssh2_tib <- getSsheetIndexAnnoTib(idx = 1, minOobPval = 0.1, minOobPerc = 90, minNegPval = 0.02, minNegPerc = 98, minDb = 0.2, verbose = 10)
  
  dtmp_tib <- getSsheetDataTab(tib = rdat$ssheet_tib, 
                               minOobPval = 0.1,  minOobPerc = 90, 
                               minNegPval = 0.02, minNegPerc = 98, minDb = 0.2, 
                               verbose = 10)
  print(dtmp_tib, n=base::nrow(dtmp_tib))
  
  
  rtib_tib <- getSsheetDataTab(tib = rdat$ssheet_tib, 
                               minOobPval = 0.1,  minOobPerc = 90, 
                               minNegPval = 0.02, minNegPerc = 98, minDb = 0.2, 
                               verbose = 10)
  
  
  sesame::formatVCF(sset = rdat$open_sset_dat, vcf = file.path(opt$outDir,'tmp.vcf'))
  
  # Sentrix_Name
  # Failed_QC
  # Min_Pass_Perc
  rdat$ssheet_tib$Sentrix_Name
  rdat$ssheet_tib$cg_Failed_QC_basic_0
  rdat$ssheet_tib$cg_pass_perc_basic_0
  
  open_req_tib <- requeueFlagOpenSesame(tib=rdat$ssheet_tib, name=rdat$ssheet_tib$Sentrix_Name,
                                        verbose = 10)
  
  #
  #
  # TBD:: Two things::
  #   - Capture r2_val/db_val described below
  #   - Report numerator and denominator in pval perc passing calculations
  #   - Add run_name to Auto Sample Sheet
  #
  
  #
  # rdat$ssheet_tib %>% dplyr::select(dplyr::contains("basic"))
  #
  
  dat_dir <- '/Users/bretbarnes/Documents/scratch/mlk/swifthoof_main/CNTL-Samples_VendA_10092020'
  dat_lst <- list.files(dat_dir, pattern='_AutoSampleSheet.csv.gz$', full.names=TRUE)
  ssh_tib <- lapply(dat_lst, readr::read_csv) %>% dplyr::bind_rows()
  
  ssh_tib %>% 
    dplyr::select(Sentrix_Name,AutoSample_dB_1_Key_1,dplyr::contains('pass_perc')) %>% 
    dplyr::select(Sentrix_Name,AutoSample_dB_1_Key_1,dplyr::starts_with('cg')) %>%
    dplyr::filter(!stringr::str_starts(AutoSample_dB_1_Key_1,'T')) %>%
    dplyr::arrange(AutoSample_dB_1_Key_1,cg_pass_perc_basic_0) %>%
    dplyr::select(!dplyr::contains('_PnegEcdf_')) %>%
    as.data.frame()
  
  #
  # Capture both of these values to demonstrate the value of normalization via Sesame::
  #
  r2_val <- rdat$open_beta_tib %>% 
    tibble::enframe(name="Probe_ID", value="betas") %>%
    dplyr::inner_join(rdat$org_list$beta$beta_dat, 
                      by="Probe_ID", suffix=c("_ref","_can")) %>%
    tibble::column_to_rownames(var="Probe_ID") %>% 
    as.matrix() %>% cor() %>% as_tibble() %>% 
    head(n=1) %>% dplyr::pull(2)
  
  dB_val <- rdat$open_beta_tib %>% 
    tibble::enframe(name="Probe_ID", value="betas") %>%
    dplyr::inner_join(rdat$org_list$beta$beta_dat, 
                      by="Probe_ID", suffix=c("_ref","_can")) %>%
    dplyr::filter(stringr::str_starts(Probe_ID,'cg')) %>% 
    dplyr::mutate(dB=base::abs(betas_ref-betas_can)) %>% 
    dplyr::summarise(pass_perc=cntPer_lte(dB,opt$minDeltaBeta)) %>%
    head(n=1) %>% dplyr::pull(1)
  
  
  
  
  
  dB_tib <- 
    rdat$open_beta_tib %>% 
    tibble::enframe(name="Probe_ID", value="betas") %>%
    dplyr::inner_join(rdat$org_list$beta$beta_dat, 
                      by="Probe_ID", suffix=c("_ref","_can")) %>%
    dplyr::filter(stringr::str_starts(Probe_ID,'cg')) %>% 
    dplyr::mutate(dB=betas_ref-betas_can)
  
  # dplyr::mutate(dB=base::abs(betas_ref-betas_can))
  
  dB_tib %>% dplyr::summarise(pass_perc=cntPer_lte(dB,opt$minDeltaBeta), pcount=count(dB<0.2),tot_cnt=n(),per2=pcount/tot_cnt)
  
  
  
  # Generate open sset
  #  - Get open_call_tib
  #
  # Get manifest = null stats + requeue
  # Get manifest = EPIC stats
  #
  
  if (stringr::str_starts(rdat$ssheet_tib$detect_manifest, 'EPIC') ) {
    open_sset_dat <- sesame::openSesame(x=rdat$prefix, what = 'sigset')
    
    open_sum2_tib <- ssetToPassPercSsheet(sset=open_sset_dat, 
                                          man=rdat$sman, min=0.1, per=90,
                                          verbose=10)
    
    open_sum1_tib <- ssetToPassPercSsheet(sset=open_sset_dat,
                                          min=0.1, per=90,
                                          verbose=10)
    
    print(open_sum1_tib)
    print(open_sum2_tib)
    
  }
  
  open_dat %>% dplyr::mutate(Workflow_idx=0) %>% requeueFlag(idx=0)
  
  # Open Sesame Comparison::
  open_sset <- sesame::openSesame(x=rdat$prefix, what = 'sigset')
  open_pval_tib <- open_sset@pval %>% tibble::enframe(name="Probe_ID", value="raw_pvals_pOOBAH")
  
  
  # Add the full filter without Inf I/II spliting...
  open_pval_tib %>% 
    dplyr::filter(stringr::str_starts(Probe_ID,'cg')) %>% 
    dplyr::summarise(pass_perc=cntPer_lte(raw_pvals_pOOBAH, min=0.1), 
                     total_cnt=n(), 
                     pass_cnt=count(raw_pvals_pOOBAH<0.1), 
                     pass_perc2=round(pass_cnt/total_cnt, 3)) %>%
    purrr::set_names(paste('cg',names(.),'basic_0', sep='_'))
  
  # TBD:: Manifest loading function needs to retain original IDs
  #
  open_pval_tib %>% dplyr::left_join(rdat$sman, by="Probe_ID") %>% 
    dplyr::group_by(Probe_Type,Probe_Design) %>% 
    dplyr::summarise(pass_perc=cntPer_lte(raw_pvals_pOOBAH, min=0.1), 
                     Total_Count=n(), Pass_Count=count(raw_pvals_pOOBAH<0.1), 
                     Pass_Perc2=round(Pass_Count/Total_Count, 3))
  
  tmp_outDir <- '/Users/bretbarnes/Documents/tmp'
  open_sum_dat <- ssetToSummary(sset = open_sset, man = rdat$sman, idx=2, 
                                workflow='open', name='open', outDir=tmp_outDir,
                                min_percs = c(98), min_pvals = c(0.1), pvals = c('pOOBAH'),
                                makek_pred=FALSE, fresh=TRUE,
                                verbose=10)
  
  
  
  val_rdat <- rdat
  
  # v.4.6 values::
  # val_rdat$ssheet_tib$cg_1_pvals_pOOBAH_pass_perc_1 # 87.903
  # val_rdat$ssheet_tib$cg_2_pvals_pOOBAH_pass_perc_1 # 83.219
  
  val_rdat$ssheet_tib$cg_1_pvals_pOOBAH_pass_perc_1
  val_rdat$ssheet_tib$cg_2_pvals_pOOBAH_pass_perc_1
  
  cur_sentrix_name <- rdat$ssheet_tib$Sentrix_Name
  
  old_dir <- '/Users/bretbarnes/Documents/tools/bk/docker-repo/Infinium_Methylation_Workhorse.v.4.6/scratch/swifthoof/swifthoof_main'
  old_dir <- "/Users/bretbarnes/Documents/scratch/RStudio/swifthoof_main/CNTL-Samples_VendA_10092020"
  call_old_csv <- file.path(old_dir, paste(cur_sentrix_name, 'EPIC_B4_raw.call.dat.csv.gz', sep="_"))
  
  call_old_tib <- readr::read_csv(call_old_csv) %>%
    dplyr::left_join(rdat$sman, by="Probe_ID") %>%
    dplyr::arrange(Probe_ID)
  
  # Calculate Pass Percent With/Without NA's
  call_old_tib %>% 
    # dplyr::filter(!is.na(betas)) %>% 
    dplyr::group_by(Probe_Type,Probe_Design) %>%
    dplyr::summarise(Total_Count=n(), 
                     Pass_Perc_Poob=cntPer_lte(pvals_pOOBAH,0.1),
                     Pass_Perc_Negs=cntPer_lte(pvals_PnegEcdf,0.05),
                     .groups="drop")
  
  rdat$new_sset %>% ssetToTib(source="pvals", name="pOOBAH") %>% 
    dplyr::left_join(rdat$sman, by="Probe_ID") %>% 
    dplyr::group_by(Probe_Type,Probe_Design) %>% 
    dplyr::summarise(Total_Count=n(), Pass_Cnt=count(pvals_pOOBAH<=0.1), Pass_Per=Pass_Cnt/Total_Count, Pass_Perc=cntPer_lte(pvals_pOOBAH,0.1), .groups="drop")
  
  
  rdat$cur_list$call_dat %>% dplyr::select(1,2) %>% ssetTibToSummary()
  
  
  
  # Basic Calculation::
  cg_inf1_cnt <- call_old_tib %>% dplyr::filter(Probe_Type=='cg' & Probe_Design==1) %>% base::nrow()
  cg_inf2_cnt <- call_old_tib %>% dplyr::filter(Probe_Type=='cg' & Probe_Design==2) %>% base::nrow()
  
  cg_pas1_cnt <- call_old_tib %>% dplyr::filter(Probe_Type=='cg' & Probe_Design==1 & pvals_pOOBAH <= 0.1) %>% base::nrow()
  cg_pas2_cnt <- call_old_tib %>% dplyr::filter(Probe_Type=='cg' & Probe_Design==2 & pvals_pOOBAH <= 0.1) %>% base::nrow()
  
  cg_pas1_cnt / cg_inf1_cnt
  cg_pas2_cnt / cg_inf2_cnt
  
  #
  # Old analysis
  #
  rdat1 <- rdat
  dB1_tib <- dplyr::inner_join( 
    purrr::set_names(rdat1$org_list$call_dat, c("Probe_ID", "can_poob", "can_negs", "can_beta")),
    dplyr::select(auto_sam_tib, Probe_ID, rdat1$ssheet_tib$AutoSample_dB_Key_1) %>% 
      purrr::set_names(c("Probe_ID","ref_beta")), 
    by="Probe_ID") %>% 
    dplyr::mutate(
      del_beta=can_beta-ref_beta,
      abs_beta=base::abs(can_beta-ref_beta)) %>%
    dplyr::inner_join(rdat$sman, by="Probe_ID") # %>% dplyr::mutate(DESIGN=as.factor(DESIGN))
  
  ind1_tib <- dplyr::inner_join(
    rdat1$cur_list$call_dat %>% dplyr::select(Probe_ID, dplyr::starts_with('ind_')) %>%
      purrr::set_names(c("Probe_ID", "can_poob", "can_negs", "can_beta")),
    
    dplyr::select(auto_sam_tib, Probe_ID, rdat1$ssheet_tib$AutoSample_dB_Key_1) %>% 
      purrr::set_names(c("Probe_ID","ref_beta")), 
    by="Probe_ID") %>% 
    dplyr::mutate(
      del_beta=can_beta-ref_beta,
      abs_beta=base::abs(can_beta-ref_beta)) %>%
    dplyr::inner_join(rdat$sman, by="Probe_ID")
  
  dB1_tib %>%
    ggplot2::ggplot(aes(x=can_poob, y=del_beta, color = DESIGN)) + 
    # ggplot2::geom_point() + 
    # ggplot2::geom_abline() +
    ggplot2::geom_density2d()
  
  ind1_tib %>%
    dplyr::filter(Probe_Type=='cg') %>%
    ggplot2::ggplot(aes(x=can_poob, y=del_beta, color = DESIGN)) + 
    ggplot2::geom_point() +
    ggplot2::geom_abline() +
    # ggplot2::geom_density2d() +
    ggplot2::facet_grid(rows = "DESIGN")
  
  
  
  #
  # Plot histogram of delta-beta split by detection p-value
  #
  dB1_tib %>% 
    # dplyr::filter(can_poob<0.05) %>%
    ggplot2::ggplot(group = DESIGN) + 
    ggplot2::geom_density(aes(x=del_beta))
  
  ggplot2::ggplot(data=dB1_tib) + 
    ggplot2::geom_density(aes(x=del_beta))
  
  
  
  
  dB7_tib <- dplyr::inner_join( 
    purrr::set_names(rdat7$org_list$call_dat, c("Probe_ID", "can_poob", "can_negs", "can_beta")),
    dplyr::select(auto_sam_tib, Probe_ID, rdat7$ssheet_tib$AutoSample_dB_Key_1) %>% 
      purrr::set_names(c("Probe_ID","ref_beta")), 
    by="Probe_ID") %>% 
    dplyr::mutate(del_beta=base::abs(can_beta-ref_beta))
  
  ggplot2::ggplot(data=dB7_tib, aes(x=can_poob, y=del_beta)) + 
    # ggplot2::geom_point() + 
    ggplot2::geom_density2d() +
    ggplot2::geom_abline()
}

# tim_csv <- file.path(opt$outDir, paste(par$prgmTag,'time-tracker.csv.gz', sep='.') )
# tim_tib <- pTracker$time %>% dplyr::mutate_if(is.numeric, list(round), 4)
# readr::write_csv(tim_tib, tim_csv)

# rdat$cur_list$call_dat %>% dplyr::filter(stringr::str_starts(Probe_ID, 'ch'))

# rdat$cur_list$call_dat %>% dplyr::filter(stringr::str_starts(Probe_ID, 'cg')) %>% dplyr::filter(raw_pvals_pOOBAH<=0.1) %>% base::nrow()
# rdat$cur_list$call_dat %>% dplyr::filter(stringr::str_starts(Probe_ID, 'cg')) %>% base::nrow()



# End of file
