
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Probe Order Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Probe Order Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

order2stats = function(tib, verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'order2stats'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  prb_seqs_tib <- NULL
  prb_seqs_tib <- tib %>% dplyr::select(AlleleA_Probe_Sequence, AlleleB_Probe_Sequence) %>% 
    tidyr::gather(Probe_Src, Probe_Seq) %>% dplyr::filter(!is.na(Probe_Seq)) %>% 
    dplyr::mutate(Probe_Seq=stringr::str_to_upper(Probe_Seq))
  
  prb_sum_tib <- tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(Assay_Design_Id, 1,2)) %>% 
    dplyr::group_by(Probe_Type) %>% dplyr::summarise(Type_Count=n()) %>% tidyr::spread(Probe_Type, Type_Count) %>%
    purrr::set_names(paste(names(.),"Count", sep='_'))
  print(prb_sum_tib)
  
  tot_ids_cnt <- tib %>% base::nrow()
  unq_ids_cnt <- tib %>% dplyr::distinct(Assay_Design_Id) %>% base::nrow()
  
  tot_seq_cnt <- prb_seqs_tib %>% base::nrow()
  unq_seq_cnt <- prb_seqs_tib %>% dplyr::distinct(Probe_Seq) %>% base::nrow()
  
  inf1_cnt <- tib %>% dplyr::filter(!is.na(AlleleB_Probe_Id) & stringr::str_length(AlleleB_Probe_Id)!=0) %>% base::nrow()
  inf2_cnt <- tib %>% dplyr::filter( is.na(AlleleB_Probe_Id) | stringr::str_length(AlleleB_Probe_Id)==0) %>% base::nrow()
  inf_tot_cnt <- (2*inf1_cnt) + inf2_cnt
  
  stats <- tibble::tibble(tot_ids_cnt := tot_ids_cnt,
                          unq_ids_cnt := unq_ids_cnt,
                          inf1_cnt    := inf1_cnt,
                          int2_cnt    := inf2_cnt,
                          inf_tot_cnt := inf_tot_cnt,
                          tot_seq_cnt := tot_seq_cnt,
                          unq_seq_cnt := unq_seq_cnt) %>%
    dplyr::bind_cols(prb_sum_tib)
  
  if (verbose>=vt) print(stats)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} DONE...{RET}{RET}"))
  
  
  stats
}

# End of file
