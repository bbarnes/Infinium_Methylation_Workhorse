

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Test Methods for COVID Direct Detection::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  covid_man_csv <- '/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/dat/manifest/base/COVID-C1.manifest.sesame-base.cpg-sorted.csv.gz'
  covid_man_tib <- readr::read_csv(covid_man_csv)
  
  top_man_tib <- rdat$sman
  raw_sset <- rdat$raw_sset
  
  retData <- TRUE
  verbose=10
  tc=1
  vt=1
  tTracker=NULL
  opts = opt
  
  out_name="test"
  
  idx <- 0
  cur_workflow <- 'raw'
  cur_sset <- raw_sset
  
  cur_sset_rds <- NULL
  cur_sigs_csv <- NULL
  cur_ssum_csv <- NULL
  cur_call_csv <- NULL
  
  ret_dat_list <- NULL
  cur_dat_list <- NULL
  
  auto_beta_key <- NULL
  auto_negs_key <- NULL
  
  calls_tib <- NULL
  cur_dat_list = ssetToSummary(
    sset=cur_sset, man=top_man_tib, idx=idx, workflow=cur_workflow,
    name=out_name, outDir=opts$outDir,
    
    write_sset=opts$write_sset, sset_rds=cur_sset_rds, ret_sset=retData,
    write_sigs=opts$write_sigs, sigs_csv=cur_sigs_csv, ret_sigs=retData,
    write_ssum=opts$write_ssum, ssum_csv=cur_ssum_csv, ret_ssum=retData,
    write_call=opts$write_call, call_csv=cur_call_csv, ret_call=retData,
    
    minNegPval=opts$minNegPval,minOobPval=opts$minOobPval,
    percision_sigs=opts$percision_sigs,
    percision_beta=opts$percision_beta,
    percision_pval=opts$percision_pval,
    
    by="Probe_ID", type="Probe_Type", des="Probe_Class",
    fresh=opts$fresh,
    verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
  
  
  
  
  
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #            Look at Infinium II Intensity Distribution by Beta::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  
  raw_inf2_tib <- raw_sset@II %>% tibble::as_tibble(rownames = "Probe_ID") %>% 
    dplyr::filter(stringr::str_starts(Probe_ID, 'cg')) %>%
    tidyr::separate(Probe_ID, into=c("Seq_ID", "Seq_Pos", "Seq_Str"), 
                    sep='_', remove=FALSE) %>% 
    tidyr::separate(Seq_Str, into=c("Blank","Strand_TB", "Strand_CO", "Design_Type", "Rep_Num"), 
                    sep="", convert=TRUE) %>% dplyr::select(-Blank) %>%
    dplyr::mutate(Probe_Key=paste(Seq_ID,Seq_Pos,Strand_TB,Strand_CO, sep='_'),
                  Probe_Key=as.factor(Probe_Key),
                  Beta=M/(M+U)) %>%
    # dplyr::group_by(Probe_Key) %>% 
    dplyr::arrange(Beta)
  
  
  ggplot2::ggplot(raw_inf2_tib, aes(x=Probe_Key, y=Beta)) + 
    ggplot2::geom_boxplot()
  
}

# End of file
