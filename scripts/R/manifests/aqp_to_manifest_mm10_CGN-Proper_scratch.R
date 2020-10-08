
if (FALSE) {
  #
  # Most (Basic) Approach::
  #
  
  # Most Common::
  most_des_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCm10-21092020_improbe-designOutput.LEGX-core.tsv.gz'
  most_des_tib <- loadIMP(file=most_des_tsv,verbose=opt$verbose,vt=4,tc=0,tt=pTracker)
  
  most_mat_tib <- most_des_tib %>%
    dplyr::mutate(Mat_CGN=stringr::str_remove(Seq_ID,'^cg') %>% stringr::str_remove('^0+') %>% as.double(),
                  Mat_TB=stringr::str_sub(Methyl_Allele_TB_Strand, 1,1),
                  Mat_CO=Methyl_Allele_CO_Strand) %>% 
    dplyr::select(Mat_CGN,Mat_TB,Mat_CO, everything())
  
  # Most Unique::
  top_unq_tib <- raw_man_tib %>% dplyr::distinct(Mat_PrbA,M,U,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence)
  
  # Join Data::
  most_int_tib <- most_mat_tib %>% dplyr::inner_join(raw_s48_int_tib, by=c("Mat_CGN","Mat_TB","Mat_CO"), suffix=c("_IMP","_DES"))
  
  most_int_tib %>% dplyr::filter(!is.na(Seq_ID_IMP)) %>% dplyr::filter(!is.na(Seq_ID_DES))
  most_int_tib %>% dplyr::filter(!is.na(Seq_ID_IMP)) %>% dplyr::filter(!is.na(Seq_ID_DES)) %>% dplyr::filter(Seq_ID_IMP == Seq_ID_DES) %>% dplyr::select(Seq_ID_IMP,Seq_ID_DES)
  most_int_tib %>% dplyr::filter(!is.na(Seq_ID_IMP)) %>% dplyr::filter(!is.na(Seq_ID_DES)) %>% dplyr::filter(Seq_ID_IMP != Seq_ID_DES) %>% dplyr::select(Seq_ID_IMP,Seq_ID_DES)
  
  most_int_rev_tib <- raw_s48_int_tib %>% dplyr::inner_join(most_mat_tib, by=c("Mat_CGN","Mat_TB","Mat_CO"), suffix=c("_DES","_IMP"))
  most_prb_inp_tib <- most_int_rev_tib %>% dplyr::select(Seq_ID_IMP,Top_Sequence) %>% dplyr::distinct() %>% 
    dplyr::rename(Seq_ID=Seq_ID_IMP,Forward_Sequence=Top_Sequence) %>% dplyr::arrange(Seq_ID) %>% dplyr::mutate(Probe_Type='cg')
  
  #
  # BUILD PROBES::
  #
  opt$load_probes  <- TRUE
  opt$load_probes  <- FALSE
  most_prb_des_csv <- file.path(opt$outDir, paste(opt$runName,'most-all-strands-current-probes.csv.gz', sep='.') )
  if (opt$load_probes && file.exists(most_prb_des_csv)) {
    most_prb_des_tib <- readr::read_csv(most_prb_des_csv)
  } else {
    most_prb_des_tib <- tib2prbs(tib=most_prb_inp_tib, idsKey="Seq_ID", prbKey="Probe_Type", 
                                 seqKey="Forward_Sequence", verbose=opt$verbose+10)
    
    readr::write_csv(most_prb_des_tib,most_prb_des_csv)
  }
  raw_man_tib %>% dplyr::group_by(Probe_Type,Infinium_Design) %>% dplyr::summarise(Count=n())
  
  #
  # Check if upper case makes a difference::
  #
  most_prb_des_MAT_tib <- most_prb_des_tib %>% 
    dplyr::mutate(PRB1_U_MAT=stringr::str_to_upper(PRB1_U),
                  PRB1_M_MAT=stringr::str_to_upper(PRB1_M),
                  PRB2_D_MAT=stringr::str_to_upper(PRB2_D),
                  TB_Str=dplyr::case_when(
                    FR_Str=='F' ~ 'C', FR_Str=='R' ~ 'O', TRUE ~ NA_character_
                  )
    )
  
  mat_inf1_tib <-  dplyr::inner_join(raw_man_tib, most_prb_des_MAT_tib, 
                                     by=c("AlleleA_Probe_Sequence"="PRB1_U_MAT","AlleleB_Probe_Sequence"="PRB1_M_MAT"),
                                     suffix=c("_RAW","_DES"))
  mat_inf1_tib %>% dplyr::group_by(Probe_Type_RAW,Infinium_Design) %>% 
    dplyr::summarise(Count=n(), .groups='drop') %>% print()
  
  mat_inf2_tib <- dplyr::inner_join(raw_man_tib,most_prb_des_MAT_tib, 
                                    by=c("AlleleA_Probe_Sequence"="PRB2_D_MAT"), 
                                    suffix=c("_RAW","_DES"))
  mat_inf2_tib %>% dplyr::group_by(Probe_Type_RAW,Infinium_Design) %>% 
    dplyr::summarise(Count=n(), .groups='drop') %>% print()
  
  raw_man_tib %>% dplyr::filter(Probe_Type=='cg') %>% dplyr::group_by(Probe_Type,Infinium_Design) %>% 
    dplyr::summarise(Count=n(), .groups='drop') %>% print()
  
  fin_inf1_cnt <- mat_inf1_tib %>% base::nrow()
  fin_inf2_cnt <- mat_inf2_tib %>% base::nrow()
  raw_man_cnt  <- raw_man_tib  %>% base::nrow()
  
  #
  # Determine non-matched probes
  #
  mis_man_tib <- raw_man_tib %>% dplyr::anti_join(
    dplyr::bind_rows(mat_inf1_tib,mat_inf2_tib),
    by=c("M","U") )
  mis_sum_tib <- mis_man_tib %>% dplyr::group_by(Probe_Type) %>% 
    dplyr::summarise(PT_Miss_Count=n(), .groups='drop')
  
  mat_man_tib <- raw_man_tib %>% dplyr::inner_join(
    dplyr::bind_rows(mat_inf1_tib,mat_inf2_tib),
    by=c("M","U") )
  mat_sum_tib <- mat_man_tib %>% dplyr::group_by(Probe_Type) %>% 
    dplyr::summarise(PT_Miss_Count=n(), .groups='drop')
  
  #
  # The two above should allow clear separation 
  #
  
  
  
  #
  # CURRENTLY HERE::
  #
  
  dplyr::filter(raw_man_tib, Probe_Type=='cg' & Infinium_Design==1) %>% 
    dplyr::anti_join(dplyr::filter(mat_inf1_tib, Probe_Type_DES=='cg' & Infinium_Design==1), 
                     by=c("AlleleA_Probe_Sequence"="PRB1_U_MAT","AlleleB_Probe_Sequence"="PRB1_M_MAT"))
  
  man1_tib <- dplyr::filter(raw_man_tib, Probe_Type=='cg' & Infinium_Design==1)
  fin1_tib <- dplyr::filter(mat_inf1_tib, Probe_Type_DES=='cg' & Infinium_Design==1)
  
  dplyr::anti_join(man1_tib, fin1_tib, by=c("AlleleA_Probe_Sequence","AlleleB_Probe_Sequence") )
  
  
  man2_tib <- dplyr::filter(raw_man_tib, Probe_Type=='cg' & Infinium_Design==2)
  fin2_tib <- dplyr::filter(mat_inf2_tib, Probe_Type_DES=='cg' & Infinium_Design==2)
  
  
  dplyr::anti_join(man2_tib, fin2_tib, by=c("AlleleA_Probe_Sequence","PRB2_D_MAT") )
  dplyr::anti_join(man2_tib, fin2_tib, by=c("AlleleA_Probe_Sequence","PRB2_D_MAT") )
  
  
  
  #
  #
  # RUNNING ABOVE....
  #
  #
  
  
  
  
  
  
  
  
  
  
  # Full Join::
  most_int_all_tib <- raw_s48_int_tib %>% dplyr::full_join(most_mat_tib, by=c("Mat_CGN","Mat_TB","Mat_CO"), suffix=c("_DES","_IMP"))
  
  
  
  
  most_des_tib %>% dplyr::mutate()
  
  raw_s48_int_tib %>% dplyr::select(Mat_PrbA,Mat_CGN,Mat_TB,Mat_CO) %>% dplyr::distinct()
  raw_s48_int_tib %>% dplyr::select(Mat_PrbA) %>% dplyr::distinct()
  raw_s48_int_tib %>% dplyr::select(Mat_PrbA,Mat_CGN) %>% dplyr::distinct()
  
  # Split Groups::
  #  - Match  by s48
  #  - Failed by s48
  mat_s48_man_tib <- raw_man_tib %>% dplyr::inner_join(dplyr::select(raw_s48_int_tib, M,U), by=c("M","U") )
  mtt_s48_man_tib <- raw_man_tib %>% dplyr::inner_join(raw_s48_int_tib, by=c("M","U") )
  mis_s48_man_tib <- raw_man_tib %>% dplyr::anti_join(raw_s48_int_tib, by=c("M","U") )
  
  raw_s48_int_tib %>% base::nrow()
  mat_s48_man_tib %>% base::nrow()
  mtt_s48_man_tib %>% base::nrow()
  mis_s48_man_tib %>% base::nrow()
  
  # Build Summaries::
  #
  mat_s48_sum_tib <- mat_s48_man_tib %>% dplyr::group_by(Probe_Type,AQP) %>% dplyr::summarise(Count=n(), .groups='drop')
  mis_s48_sum_tib <- mis_s48_man_tib %>% dplyr::group_by(Probe_Type,AQP) %>% dplyr::summarise(Count=n(), .groups='drop')
  raw_man_tib %>% dplyr::group_by(Probe_Type,AQP) %>% dplyr::summarise(Count=n(), .groups='drop')
  
  
  # Build intersect M,U keys::
  #
  mus_s48_int_tib <- raw_s48_int_tib %>% dplyr::distinct(M,U)
  
  mus_s48_int_tib <- raw_s48_int_tib %>% dplyr::distinct(Mat_CGN,Mat_TB,Mat_CO,M,U,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence)
  
  
  
  # Build basic stats
  #
  raw_man_tib %>% dplyr::group_by(Probe_Type,AQP) %>% dplyr::summarise(Count=n(), .groups='drop') %>% print()
  raw_s48_int_tib %>% dplyr::group_by(Probe_Type,AQP) %>% dplyr::summarise(Count=n(), .groups='drop') %>% print()
  
  # Split Groups::
  #  - Match  by s48
  #  - Failed by s48
  raw_man_tib %>% dplyr::inner_join(mus_s48_int_tib, by=c("M","U") )
  raw_man_tib %>% dplyr::inner_join(dplyr::distinct(raw_s48_int_tib, M,U), by=c("M","U") )
  raw_man_tib %>% dplyr::inner_join(dplyr::select(raw_s48_int_tib, M,U), by=c("M","U") )
  
  
  
  
  # We need to pull up design data for these cg's::
  raw_s48_int_tib %>% dplyr::distinct(Mat_CGN,Mat_TB,Mat_CO) %>% dplyr::arrange(Mat_CGN)
  raw_s48_int_tib %>% dplyr::distinct(Mat_CGN,Mat_TB,Mat_CO,Seq_ID)
  raw_s48_int_tib %>% dplyr::distinct(Mat_CGN,Mat_TB,Mat_CO,Seq_ID,Probe_Type) %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(PTCount=n())
  
  pull_seqIds_csv <- file.path(opt$outDir, 'pull-seqIds.csv.gz')
  pull_seqIds_tib <- dplyr::distinct(raw_s48_int_tib, Seq_ID)
  readr::write_csv(pull_seqIds_tib, pull_seqIds_csv)
  
  raw_s48_int_tib %>% dplyr::filter(Mat_CGN!=Seq_CGN) %>% # dplyr::select(Mat_CGN,Seq_CGN,Seq_ID,Mat_TB,Mat_CO,Probe_ID,Probe_Type) %>%
    dplyr::group_by(Probe_Type,AQP) %>% dplyr::summarise(PTCount=n())
  raw_s48_int_tib %>% dplyr::filter(Mat_CGN==Seq_CGN) %>% # dplyr::select(Mat_CGN,Seq_CGN,Seq_ID,Mat_TB,Mat_CO,Probe_ID,Probe_Type) %>%
    dplyr::group_by(Probe_Type,AQP) %>% dplyr::summarise(PTCount=n())
  
  #
  # Left off here1 ::
  #
  
  #
  # Can we load the design database??? NO!!!
  #
  
  mm10_des_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCm10-21092020_improbe-designOutput.tsv.gz'
  # Direct subsample command::
  #  ./scripts/subset/getSubset.simple.pl -header -t /Users/bretbarnes/Documents/scratch/manifests/mm10-LEGX-S3/pull-seqIds.csv.gz -d /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCm10-21092020_improbe-designOutput.tsv.gz | gzip -c -> /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCm10-21092020_improbe-designOutput.LEGX-core.tsv.gz
  #
  subset_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCm10-21092020_improbe-designOutput.LEGX-core.tsv.gz'
  subset_cmd <- glue::glue("{par$topDir}/scripts/subset/getSubset.simple.pl -header -t {pull_seqIds_tib} -d {mm10_des_tsv} | gzip -c -> {subset_tsv}")
  
  # mm10_imp_des_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCm10-21092020_improbe-designOutput.LEGX-core.tsv.gz'
  # mm10_imp_des_tib <- readr::read_tsv(mm10_des_tsv)
  
  mm10_s48_sel_tib <- dplyr::distinct(raw_s48_int_tib, Seq_ID,Mat_CGN,Mat_TB,Mat_CO) %>% dplyr::arrange(Mat_CGN)
  mm10_tar_des_tib <- dplyr::mutate(mm10_imp_des_tib, Mat_CGN=stringr::str_remove(Seq_ID, 'cg') %>% stringr::str_remove('^0+') %>% as.double(),
                                    Mat_TB=stringr::str_sub(Methyl_Allele_TB_Strand,1,1),
                                    Mat_CO=Methyl_Allele_CO_Strand) %>% dplyr::select(Seq_ID,Mat_CGN, Mat_TB, Mat_CO) %>% dplyr::arrange(Mat_CGN)
  
  
  
  #
  # Left off here2 :: Need to ensure joining...
  #
  
  mm10_s48_cgn_tib <- dplyr::inner_join(
    dplyr::distinct(raw_s48_int_tib, Mat_CGN,Mat_TB,Mat_CO, .keep_all=TRUE) %>% dplyr::arrange(Mat_CGN),
    dplyr::mutate(mm10_des_tib, Mat_CGN=stringr::str_remove(Seq_ID, 'cg') %>% stringr::str_remove('^0+') %>% as.double(),
                  Mat_TB=stringr::str_sub(Methyl_Allele_TB_Strand,1,1),
                  Mat_CO=Methyl_Allele_CO_Strand),
    by=c("Mat_CGN","Mat_TB","Mat_CO"))
  
  
  
  
  mm10_s48_cgn_tib %>% dplyr::group_by(Probe_Type.x) %>% dplyr::summarise(PT_Count=n())
  
  
  mm10_s48_des_tib <- dplyr::inner_join(
    dplyr::distinct(raw_s48_int_tib, Mat_CGN,Mat_TB,Mat_CO, .keep_all=TRUE) %>% dplyr::arrange(Mat_CGN),
    dplyr::mutate(mm10_des_tib, Mat_CGN=stringr::str_remove(Seq_ID, 'cg') %>% stringr::str_remove('^0+') %>% as.double(),
                  Mat_TB=stringr::str_sub(Methyl_Allele_TB_Strand,1,1),
                  Mat_CO=Methyl_Allele_CO_Strand),
    by=c("Mat_CGN"="Mat_CGN", "Mat_TB","Mat_CO"))
  
  mm10_s48_des_tib %>% dplyr::group_by(Probe_Type.x) %>% dplyr::summarise(PT_Count=n())
  raw_s48_int_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(PT_Count=n())
  
  # Some Stat summaries::
  #
  raw_s48_int_tib %>% dplyr::filter(Mat_CGN==Seq_CGN) %>% dplyr::group_by(AQP) %>% dplyr::summarise(AQP_Count=n())
  raw_s48_int_tib %>% dplyr::filter(Mat_CGN!=Seq_CGN) %>% dplyr::group_by(AQP) %>% dplyr::summarise(AQP_Count=n())
  
  #
  # Split by matching and mis-matching 
  #
  mm10_s48_mat_tib <- raw_s48_int_tib %>% dplyr::distinct(Seq_ID,U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, .keep_all=TRUE) %>% dplyr::filter(!is.na(Seq_ID))
  mm10_s48_mis_tib <- raw_man_tib %>% dplyr::anti_join(mm10_s48_mat_tib, by="Seq_ID")
  #                   raw_man_tib %>% dplyr::anti_join(raw_s48_int_tib, by="Mat_PrbA")
  
  
  
  
  mm10_s48_mat_tib %>% dplyr::filter(Mat_CGN!=Seq_CGN) %>% dplyr::group_by(Probe_Type,AQP) %>% dplyr::summarise(PTCount=n())
  
  mm10_s48_mat_tib %>% dplyr::group_by(Probe_Type,AQP) %>% dplyr::summarise(PTCount=n())
  mm10_s48_mis_tib %>% dplyr::group_by(Probe_Type,AQP) %>% dplyr::summarise(PTCount=n())
  
  
  
  
  mm10_s48_mat_cnt <- mm10_s48_mat_tib %>% base::nrow()
  mm10_s48_mis_cnt <- mm10_s48_mis_tib %>% base::nrow()
  mm10_s48_sum_cnt <- mm10_s48_mat_cnt + mm10_s48_mis_cnt
  
  cat(glue::glue("[{par$prgmTag}]: Match={raw_mat_cnt}={mm10_s48_sum_cnt}; mm10_s48_mat_cnt={mm10_s48_mat_cnt}, mm10_s48_mis_cnt={mm10_s48_mis_cnt}.{RET}"))
  
  mm10_s48_mat_tib %>% group_by(Probe_Type) %>% dplyr::summarise(Count=n())
  mm10_s48_mis_tib %>% group_by(Probe_Type) %>% dplyr::summarise(Count=n())
  
  # mm10_s48_mat_tib %>% dplyr::distinct(Seq_ID, .keep_all=TRUE) %>% dplyr::filter(!stringr::str_starts(Probe_ID, 'cg')) %>% dplyr::select(Seq_ID,Probe_ID)
  
  
  # Verification::
  mm10_s48_mat_tib %>% dplyr::mutate(Seq_ID2=stringr::str_remove(Probe_ID, '_.*$') ) %>%
    dplyr::filter(Seq_ID!=Seq_ID2)
  
  mm10_s48_mat_tib %>% dplyr::mutate(Seq_ID2=stringr::str_remove(Probe_ID, '_.*$') ) %>%
    dplyr::filter(Seq_ID==Seq_ID2) %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(PCount=n())
  
  # The confusing part::
  mm10_s48_mat_tib %>% dplyr::filter(!stringr::str_starts(Seq_ID, 'cg'))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if (FALSE) {
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     mm10 CGN Genome Count:: Global Data
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # mm10_cgn_cnt_col <- c('Genomic_CGN_Count', 'Seq_ID')
    # mm10_cgn_cnt_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-counts.tsv.gz'
    # mm10_cgn_cnt_tib <- readr::read_tsv(mm10_cgn_cnt_tsv, col_names=mm10_cgn_cnt_col)
    # 
    # raw_s48_int_tib %>% raw_s48_int_tib %>%
    #   dplyr::left_join(mm10_cgn_cnt_tib, by=c("Mat_CGN"="Seq_ID"))
    
    
    # Adjust mu and rp with with cg prefix::
    #
    if (FALSE) {
      raw_s48_int_tib <- mm10_s48_org_tib %>% dplyr::mutate(
        Seq_ID=dplyr::case_when(
          Probe_Type=='mu' ~ stringr::str_replace(Seq_ID, '^mu', 'cg'),
          TRUE ~ Seq_ID) ) %>% 
        dplyr::mutate(Mat_Seq_ID=stringr::str_remove(Seq_ID,'[a-zA-Z][a-zA-A]') %>% stringr::str_remove('^0+') %>% as.integer()) # %>%
      # dplyr::filter(Mat_CGN != Mat_Seq_ID)
      # dplyr::select(Mat_CGN, Mat_Seq_ID)
    } else {
      raw_s48_int_tib <- mm10_s48_org_tib %>%
        dplyr::mutate(Mat_Seq_ID=stringr::str_remove(Seq_ID,'[a-zA-Z][a-zA-A]') %>% stringr::str_remove('^0+') %>% as.integer())
    }
    
    # raw_s48_int_tib %>% 
    #   dplyr::mutate(Mat_Seq_ID=stringr::str_remove(Seq_ID,'cg') %>% stringr::str_remove('^0+') %>% as.integer()) %>% 
    #   dplyr::select(Mat_CGN, Mat_Seq_ID)
    
    # Split by matching and mismatch Seq_ID
    #
    mm10_s48_mat_tib <- raw_s48_int_tib %>% dplyr::filter(Mat_CGN==Mat_Seq_ID) %>% dplyr::arrange(Seq_ID)
    mm10_s48_mis_tib <- raw_s48_int_tib %>% dplyr::filter(Mat_CGN!=Mat_Seq_ID) %>% dplyr::arrange(Seq_ID)
    #  AND make a list of missing probes from original data::
    mm10_s48_off_tib <- raw_man_tib %>% dplyr::anti_join(raw_s48_int_tib, by="Mat_PrbA")
    
    # QC Counts Matching::
    raw_s48_int_tib %>% base::nrow()
    mm10_s48_mat_tib %>% base::nrow()
    mm10_s48_mis_tib %>% base::nrow()
    mm10_s48_off_tib %>% base::nrow()
    
    # QC: Check type distributions::
    #  mat = cg
    #  mis = rp,mu
    mm10_s48_mat_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Grp_Count=n())
    mm10_s48_mis_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Grp_Count=n())
    mm10_s48_off_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Grp_Count=n())
    
    # QC:: This should be equal to zero!!!
    #
    mm10_s48_mat_tib %>% dplyr::filter(Mat_CGN!=Seq_ID) %>% base::nrow()
    
    #
    # MU Investigation:: This just verifies all the MU probes are found in their cg analogous matching codes::
    #
    if (FALSE) {
      raw_s48_int_tib %>% dplyr::filter(Probe_Type=='mu') %>% base::nrow()
      mm10_s48_mat_tib %>% dplyr::filter(Probe_Type=='mu') %>% base::nrow()
      mm10_s48_mis_tib %>% dplyr::filter(Probe_Type=='mu') %>% base::nrow()
      
      mm10_s48_mis_tib %>% dplyr::filter(Probe_Type=='mu') %>% dplyr::filter(Mat_CGN %in% raw_s48_int_tib$Mat_CGN) %>% base::nrow()
      raw_s48_int_tib %>% dplyr::filter(Probe_Type=='mu') %>% dplyr::filter(Mat_CGN %in% mm10_s48_mis_tib$Mat_CGN) %>% base::nrow()
      raw_s48_int_tib %>% dplyr::filter(Probe_Type=='mu') %>% dplyr::anti_join(mm10_s48_mis_tib, by="Mat_CGN") %>% base::nrow()
      
      # True Unique::
      raw_s48_int_tib %>% dplyr::filter(Probe_Type=='mu') %>% dplyr::distinct(Mat_CGN) %>% base::nrow()
      mm10_s48_mat_tib %>% dplyr::filter(Probe_Type=='mu') %>% dplyr::distinct(Mat_CGN) %>% base::nrow()
      
      raw_s48_int_tib %>% dplyr::filter(Probe_Type=='mu') %>% dplyr::distinct(Mat_CGN) %>% base::nrow()
      raw_s48_int_tib %>% dplyr::filter(Probe_Type=='mu') %>% dplyr::distinct(Mat_CGN,Infinium_Design,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence) %>% base::nrow()
      raw_s48_int_tib %>% dplyr::filter(Probe_Type=='mu') %>% dplyr::distinct(Mat_CGN,Infinium_Design,M,U,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence) %>% base::nrow()
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Manifest Improbe Matching:: CGN ONLY
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Write improbe lookup file::
    mm10_s48_int_mat_cgn_tsv <- file.path(opt$outDir, paste(opt$runName,'s48_int_mat_cgn_tsv.gz', sep='.') )
    mm10_des_int_mat_cgn_tsv <- file.path(opt$outDir, paste(opt$runName,'des_int_mat_cgn_tsv.gz', sep='.') )
    mm10_s48_int_mat_cgn_tib <- mm10_s48_mat_tib %>% dplyr::distinct(Mat_CGN) %>% dplyr::arrange(Mat_CGN) %>% dplyr::rename(Seq_ID=Mat_CGN)
    
    opt$subExe <- '/Users/bbarnes/Documents/Projects/scripts/subset/getSubset.simple.pl'
    opt$impDesOutTsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.tsv.gz'
    
    run_sub_cmd <- TRUE
    run_sub_cmd <- FALSE
    if (run_sub_cmd) {
      readr::write_tsv(mm10_s48_int_mat_cgn_tib, mm10_s48_int_mat_cgn_tsv)
      sub_cmd <- glue::glue("{opt$subExe} -header -t {mm10_s48_int_mat_cgn_tsv} -d {opt$impDesOutTsv} | gzip -c - > {mm10_des_int_mat_cgn_tsv}")
      system(sub_cmd)
    }
    mm10_s48_int_mat_imp_tib <- read_tsv(mm10_des_int_mat_cgn_tsv) %>% 
      dplyr::mutate(Imp_FR=Methyl_Allele_FR_Strand,
                    Imp_TB=stringr::str_sub(Methyl_Allele_TB_Strand, 1,1),
                    Imp_CO=Methyl_Allele_CO_Strand)
    
    #
    # Join on Seq_ID for MAT only
    #
    mm10_s48_int_all_cgnTop_tib <- dplyr::distinct(mm10_s48_int_mat_imp_tib, Seq_ID,Top_Sequence) %>% 
      dplyr::add_count(Seq_ID, name="CGN_Count") %>% dplyr::add_count(Top_Sequence, name="Top_Count")
    
    mm10_s48_int_mat_cgnTop_tib <- dplyr::distinct(mm10_s48_int_mat_imp_tib, Seq_ID,Top_Sequence) %>% 
      dplyr::inner_join(mm10_s48_mat_tib, by="Seq_ID") %>% dplyr::distinct(Seq_ID, Top_Sequence) %>% 
      dplyr::add_count(Seq_ID, name="CGN_Count") %>% dplyr::add_count(Top_Sequence, name="Top_Count")
    
    mm10_s48_int_ant_cgnTop_tib <- dplyr::distinct(mm10_s48_int_mat_imp_tib, Seq_ID,Top_Sequence) %>% 
      dplyr::anti_join(mm10_s48_mat_tib, by="Seq_ID") %>% dplyr::distinct(Seq_ID, Top_Sequence) %>%
      dplyr::add_count(Seq_ID, name="CGN_Count") %>% dplyr::add_count(Top_Sequence, name="Top_Count")
    
    mm10_s48_int_mat_imp_tib %>% base::nrow()
    mm10_s48_int_all_cgnTop_tib %>% base::nrow()
    mm10_s48_int_mat_cgnTop_tib %>% base::nrow()
    mm10_s48_int_ant_cgnTop_tib %>% base::nrow()
    
    mm10_s48_int_all_cgnTop_tib %>% dplyr::filter(CGN_Count!=1) %>% dplyr::filter(Top_Count!=1)
    mm10_s48_int_mat_cgnTop_tib %>% dplyr::filter(CGN_Count!=1) %>% dplyr::filter(Top_Count!=1)
    mm10_s48_int_ant_cgnTop_tib %>% dplyr::filter(CGN_Count!=1) %>% dplyr::filter(Top_Count!=1)
    
    #
    # Write the mouse forced assigned CGN->TopSeqs
    #
    mm10_cgnToTop_csv <- file.path(opt$outDir, paste(opt$runName,'forced-defined-cgnToTopSeq.csv.gz', sep='.') )
    readr::write_csv(mm10_s48_int_all_cgnTop_tib,mm10_cgnToTop_csv)
    
    #
    # Need to reconcile Horvath vs. mm10 mismatches::
    #
    #  Conclusion:: Need to use the Horvath Assignments!!!
    #
    err_map_csv <- '/Users/bbarnes/Documents/Projects/methylation/scratch/manifest_scratch/error-cgnToTopSeq.csv.gz'
    err_map_tib <- readr::read_csv(err_map_csv)
    
    mm10_s48_int_mat_imp_tib %>% dplyr::filter(Seq_ID %in% err_map_tib$Seq_ID_mm10) %>% dplyr::pull(Probe_Type)
    mm10_s48_mat_tib %>% dplyr::filter(Seq_ID %in% err_map_tib$Seq_ID_mm10) %>% dplyr::pull(Probe_Type)
    
    #
    #
    # Left off here::
    #
    #
    
    
    
    
    
    
    #
    # Join improbe data to designs:: MAT Only
    #
    mm10_man_mat_imp_tib <- mm10_s48_mat_tib %>% 
      dplyr::left_join(mm10_s48_int_mat_imp_tib, by=c("Seq_ID", "Probe_Type","Mat_TB"="Imp_TB", "Mat_CO"="Imp_CO")) %>% 
      dplyr::mutate(Inter_Type='T') %>% 
      tidyr::unite(Probe_ID_Suffix, Mat_TB,Mat_CO,Infinium_Design,Inter_Type,Rep_Num, sep='', remove=FALSE) %>% 
      tidyr::unite(Probe_ID, Seq_ID,Probe_ID_Suffix, sep='_')
    
    mm10_man_mat_imp_nunq_tib <- mm10_man_mat_imp_tib %>% dplyr::add_count(Probe_ID, name="Probe_Count") %>% dplyr::filter(Probe_Count!=1)
    
    
    
    mm10_man_mat_ses_tib <- mm10_man_mat_imp_tib %>%
      dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, AlleleA_Probe_Sequence,AlleleB_Probe_Sequence)
    
    # QC Counts Matching::
    mm10_man_mat_imp_tib %>% base::nrow()
    mm10_s48_mat_tib %>% base::nrow()
    mm10_man_mat_imp_tib %>% dplyr::distinct(Probe_ID) %>% base::nrow()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Manifest Improbe Matching:: Off ONLY
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Split ch/mu/rp analysis
    mm10_s48_off_types <- mm10_s48_off_tib %>% split(.$Probe_Type)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Manifest Improbe Matching:: rs ONLY
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Design Source::
    #  rs = /Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/Redesign/data/SNP/selected_SNP_probes.bed
    
    #
    # TBD:: Add back previous SNPs!!!!!       dplyr::filter(TB ===== 'N' & !is.na(TB) )         Missing most SNPs now...
    #
    mm10_man_newRS_ses_tib <- mm10_s48_off_types[['rs']] %>% dplyr::filter(TB != 'N' & !is.na(TB) ) %>%
      dplyr::mutate(Inter_Type='T') %>% 
      tidyr::unite(Probe_ID_Suffix, TB,CO,Infinium_Design,Inter_Type,Rep_Num, sep='', remove=FALSE) %>% 
      tidyr::unite(Probe_ID, Seq_ID,Probe_ID_Suffix, sep='_') %>%
      dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, AlleleA_Probe_Sequence,AlleleB_Probe_Sequence)
    
    mm10_man_newRS_ses_tib %>% base::nrow()
    mm10_man_newRS_ses_tib %>% dplyr::distinct(Probe_ID) %>% base::nrow()
    
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Sesame Manifest::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # current = B2
    outName <- paste(opt$platform,opt$version, sep='-')
    outFile <- paste(outName,'manifest.sesame-base.cpg-sorted.csv.gz', sep='.')
    mm10_ses_man_cgn_out_csv <- file.path(opt$outDir, outFile)
    mm10_ses_man_cgn_git_csv <- file.path(par$datDir,'manifest/base',outFile)
    mm10_ses_man_cgn_tib <- dplyr::bind_rows(mm10_man_mat_ses_tib,
                                             mm10_man_newRS_ses_tib,
                                             ses_unq_ctl_tib) %>% dplyr::arrange(Probe_ID)
    readr::write_csv(mm10_ses_man_cgn_tib,mm10_ses_man_cgn_out_csv)
    readr::write_csv(mm10_ses_man_cgn_tib,mm10_ses_man_cgn_git_csv)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Manifest Improbe Matching:: Non ONLY
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Split ch/mu/rp analysis
    mm10_s48_mis_types <- mm10_s48_mis_tib %>% split(.$Probe_Type)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Manifest Improbe Matching:: mu ONLY
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # Conclusion:: mu probes can be converted to cg and added back to mat_tib with extra cg's added as extra field...
    #
    clean_mu_all_tib <- mm10_s48_mis_types[['mu']] %>% dplyr::mutate(Seq_ID=stringr::str_replace(Seq_ID, '^mu', 'cg'))
    clean_mu_mat_tib <- clean_mu_all_tib %>% dplyr::filter(Mat_CGN==Seq_ID) %>% dplyr::add_count(Seq_ID,Mat_TB,Mat_CO,Infinium_Design, name="ID_Count")
    clean_mu_mis_tib <- clean_mu_all_tib %>% dplyr::filter(! Seq_ID %in% clean_mu_mat_tib$Seq_ID)
    
    clean_mu_mat_tib %>% dplyr::filter(ID_Count!=1) %>% as.data.frame()
    
    # Write improbe lookup file::
    mm10_s48_int_mis_cgn_tsv <- file.path(opt$outDir, 'mm10_s48_int_mis_cgn_tsv.gz')
    mm10_s48_int_mis_cgn_tib <- clean_mu_all_tib %>% dplyr::distinct(Mat_CGN) %>% dplyr::arrange(Mat_CGN) %>%
      dplyr::rename(Seq_ID=Mat_CGN)
    # readr::write_tsv(mm10_s48_int_mis_cgn_tib, mm10_s48_int_mis_cgn_tsv)
    
    
    
    
    
    
    
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #               Quick Look at Auto Sample Sheets from CG Only::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # mm10_hum_ss_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/Laird-IDs_to_SampleNames_basic.csv'
    # mm10_hum_ss_tib <- readr::read_csv(mm10_hum_ss_csv) %>% purrr::set_names(c('Laird_ID', 'Sentrix_Pos', 'Cell_Type', 'Sample_Class'))
    
    # Possible mismatch: 20364 = 20382
    RepAC_ID <- 20364
    RepAC_ID <- 20382
    
    sam_map_tib <- tibble::tibble(
      Laird_ID   =c( RepAC_ID,  20384,  20385,  21026,   20010,   20015,   20012),
      Sample_Name=c('RepAC','RepS3','RepM1','RepSA', 'T00DZ', 'T50DZ', 'T99DZ')
    )
    
    tit_map_tib <- tibble::tibble(
      Laird_ID   =c( 20010,   20015,   20012),
      Sample_Name=c('T00DZ', 'T50DZ', 'T99DZ')
    )
    
    rep_map_tib <- tibble::tibble(
      Laird_ID   =c( RepAC_ID,  20384,  20385,  21026),
      Sample_Name=c('RepAC','RepS3','RepM1','RepSA')
    )
    
    # Load ILS Sample Sheet::
    #
    ils_map_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/ILS-Mouse_Methylation_samplesheet.csv'
    ils_map_tib <- readr::read_csv(ils_map_csv) %>%
      dplyr::mutate(Laird_ID=as.integer(stringr::str_remove(Sample_ID, '_.*$')),
                    Sentrix_ID=SentrixBarcode_A,
                    Sentrix_Pos=SentrixPosition_A,
                    Sentrix_Name=paste(Sentrix_ID,Sentrix_Pos, sep='_')
      ) %>% 
      dplyr::select(Laird_ID,Sentrix_Name,Sentrix_ID,Sentrix_Pos) %>% dplyr::mutate(Lab='ILS')
    
    # Load VAI Sample Sheet::
    #
    vai_map_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/VAI_SampleSheetPlate1-2.no-header.csv'
    vai_map_tib <- readr::read_csv(vai_map_csv) %>% 
      dplyr::rename(Laird_Str=Sample_Name, Sentrix_Pos=Sentrix_Position) %>% 
      dplyr::mutate(Laird_ID=stringr::str_remove(Laird_Str, '_.*$') %>% as.integer(),
                    Sentrix_Name=paste(Sentrix_ID,Sentrix_Pos, sep='_')) %>% 
      dplyr::select(Laird_ID,Sentrix_Name,Sentrix_ID,Sentrix_Pos) %>% dplyr::mutate(Lab='VAI')
    
    #
    # Build beta full Sample Sheet::
    #
    
    # Full Analytical Sample Sheet::
    beta_ss_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/analytical/ILS-VAI.analytical_SampleSheet.csv.gz'
    beta_ss_tib <- dplyr::bind_rows(
      dplyr::inner_join(ils_map_tib, sam_map_tib, by="Laird_ID"),
      dplyr::inner_join(vai_map_tib, sam_map_tib, by="Laird_ID") ) %>%
      dplyr::select(Sentrix_Name,Sample_Name, everything())
    readr::write_csv(beta_ss_tib,beta_ss_csv)
    
    # Titration Analytical Sample Sheet::
    beta_tit_ss_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/analytical/Titration.ILS-VAI.analytical_SampleSheet.csv.gz'
    beta_tit_ss_tib <- beta_ss_tib %>% dplyr::filter(Laird_ID %in% tit_map_tib$Laird_ID)
    readr::write_csv(beta_tit_ss_tib,beta_tit_ss_csv)
    
    # Replicate Analytical Sample Sheet::
    beta_rep_ss_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/analytical/Replicates.ILS-VAI.analytical_SampleSheet.csv.gz'
    beta_rep_ss_tib <- beta_ss_tib %>% dplyr::filter(Laird_ID %in% rep_map_tib$Laird_ID)
    readr::write_csv(beta_rep_ss_tib,beta_rep_ss_csv)
    
    # cp /Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/analytical/* /Users/bbarnes/Documents/Projects/methylation/tools/Infinium_Methylation_Workhorse/dat/sampleSheets/mm10/
    
    
    #
    # Done with Sample Sheets!!!
    #
    
    #
    #
    # Need to write Titration and Replicates Sample Sheets seperately...
    #
    #
    
    #
    # Actual Idats and plots::
    #
    ils_dat_dir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/scratch/swifthoof_main/ILMN_mm10_betaTest_17082020'
    ils_ss_tib  <- loadAutoSampleSheets(dir=ils_dat_dir, verbose=opt$verbose) %>% dplyr::mutate(Lab='ILS')
    
    van_dat_dir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/scratch/swifthoof_main/VanAndel_mm10_betaTest_31082020'
    van_ss_tib  <- loadAutoSampleSheets(dir=van_dat_dir, verbose=opt$verbose) %>% dplyr::mutate(Lab='Van')
    
    fox_dat_dir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/scratch/swifthoof_main/MURMETVEP_mm10_betaTest_06082020'
    fox_ss_tib  <- loadAutoSampleSheets(dir=fox_dat_dir, verbose=opt$verbose) %>% dplyr::mutate(Lab='FOX')
    
    mm10_ss_tib <- dplyr::bind_rows(ils_ss_tib, van_ss_tib, fox_ss_tib) %>% dplyr::group_by(Lab)
    ggplot2::ggplot(data=mm10_ss_tib, aes(x=Beta_2_Mean, y=Poob_Pass_0_Perc, color=Lab)) + ggplot2::geom_point() # + ylim(80,100)
    ggplot2::ggplot(data=mm10_ss_tib, aes(x=Beta_2_Mean, y=Poob_Pass_0_Perc, color=Lab)) + ggplot2::geom_point() + ylim(75,100)
    
    
    #
    # Evidence For Missing Replicate Group::
    #
    ils_ss_tib %>% dplyr::inner_join(beta_ss_tib, by="Sentrix_Name") %>% dplyr::group_by(Sample_Name) %>% dplyr::summarise(SGroup_Count=n())
    
    #
    # Current Merged Data::
    #
    lifeDir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics'
    
    mm10_tit_ss_merge_csv <- file.path(lifeDir, 'scratch/merge_builds/LEGX/S1/Sample_Name/mm10-ILS-VAI.Titration/mm10-ILS-VAI.Titration_LEGX_S1_AutoSampleSheet.csv.gz')
    mm10_tit_ss_merge_tib <- readr::read_csv(mm10_tit_ss_merge_csv)
    
    mm10_rep_ss_merge_csv <- file.path(lifeDir, 'scratch/merge_builds/LEGX/S1/Sample_Name/mm10-ILS-VAI.Replicates/mm10-ILS-VAI.Replicates_LEGX_S1_AutoSampleSheet.csv.gz')
    mm10_rep_ss_merge_tib <- readr::read_csv(mm10_rep_ss_merge_csv)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Source and Design Search::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (FALSE) {
      # Objectives::
      #
      #  1. Gather all mouse idats tangos from a few recent idats
      #
      #  2. Gather all public product improbe GGNs and designs 
      #     - Add Horvath Manifests
      #     - Add NZT Manifest
      #     - Add Excalibur Manifest
      #     - Add Custom Manifests (Genknowme, TruDX, Univ Chicago)
      #  3. Gather all mouse methylation improbe CGNs and designs
      #
      #  4. Comapre Previous Products with Mouse for CGN -> Top Integrity
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                   Identify Address Not Found in IDATs::
      #
      #  1. Gather all mouse idats tangos from a few recent idats
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      par$compare_idats <- FALSE
      if (par$compare_idats) {
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                         Load Tangos From IDATS::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        idat_vec <- NULL
        if (!is.null(opt$idats)) idat_vec <- opt$idats %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
        
        max_idat_load <- 3
        prefix_list <- sesame::searchIDATprefixes(idat_vec[1], recursive=TRUE, use.basename=TRUE)
        idat_add_tib <- NULL
        for (pIdx in c(1:length(prefix_list))) {
          idat_cur_dat <- prefixToIdat(prefix=prefix_list[[pIdx]], verbose=opt$verbose)
          idat_add_tib <- idat_add_tib %>% dplyr::bind_rows(idat_cur_dat$sig %>% dplyr::select(Address))
          
          if (pIdx >= max_idat_load) break
        }
        adds_cnt_tib <- idat_add_tib %>% dplyr::group_by(Address) %>% dplyr::summarise(Count=n())
        idat_add_vec <- adds_cnt_tib %>% dplyr::pull(Address)
        
        # Sanity Check
        adds_mis_cnt <- adds_cnt_tib %>% dplyr::filter(Count != max_idat_load) %>% base::nrow()
        stopifnot(adds_mis_cnt==0)
        
        add_mis_man_tib <- dplyr::bind_rows(
          raw_man_tib %>% dplyr::filter(! U %in% idat_add_vec),
          raw_man_tib %>% dplyr::filter(!is.na(M)) %>% dplyr::filter(! M %in% idat_add_vec)
        )
        add_mis_ctl_tib <- dplyr::bind_rows(
          ses_ctl_tib %>% dplyr::filter(! U %in% idat_add_vec),
          ses_ctl_tib %>% dplyr::filter(!is.na(M)) %>% dplyr::filter(! M %in% idat_add_vec)
        )
        
        # List of Probe IDs with missing tangos in idats from Wanding
        add_mis_cgn_wand_tsv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests/tango-issues/probe_with_tango_issue.tsv.gz'
        add_mis_cgn_wand_tib <- readr::read_tsv(add_mis_cgn_wand_tsv)
        
        #
        # CONCLUSION: I don't see any missing addresses: Only checked ILMN so far...
        #  TBD:: Check IDATs from other BETA sites
        #
        
        cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Identify Address Not Found in IDATs!{RET}{RET}") )
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #           Preprocessing:: Improbe Previous Product Designs
      #
      #  2. Gather all public product improbe GGNs and designs 
      #     - Add Horvath Manifests
      #     - Add NZT Manifest
      #     - Add Excalibur Manifest
      #     - Add Custom Manifests (Genknowme, TruDX, Univ Chicago)
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      imp_prd_tsv <- '/Users/bbarnes/Documents/Projects/methylation/EWAS/improbe/EPIC-reorder.improbe-design.tsv.gz'
      imp_prd_tib <- suppressMessages(suppressWarnings( readr::read_tsv(imp_prd_tsv) )) %>% 
        dplyr::mutate(Min_Final_Score=pmin(Methyl_Final_Score,UnMethyl_Final_Score),
                      Methyl_Allele_TB_Strand=stringr::str_sub(Methyl_Allele_TB_Strand, 1,1))
      
      #
      # Investigate Integrety of CGN -> Top_Sequence in Previous Products::
      #
      
      # Build Summary:: i.e. 27k vs. 450k/EPIC
      imp_prd_tib %>% dplyr::group_by(Genome_Build) %>% dplyr::summarise(Count=n()) %>% print()
      
      # Make sure CGN -> Top is unique::
      prd_uniq_cgn_top_tib <- imp_prd_tib %>% dplyr::select(Seq_ID, Top_Sequence) %>% dplyr::distinct()
      
      # Check Top_Sequences all have a Count == 1
      prd_top_cnt_tib <- prd_uniq_cgn_top_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Count=n())
      
      # There should be 7 and they should be from the 27k:: i.e. Genome_Build == 36
      prd_top_mis_tib <- prd_top_cnt_tib %>% dplyr::filter(Count != 1)
      prd_bad_cgn_tib <- imp_prd_tib %>% dplyr::filter(Top_Sequence %in% prd_top_mis_tib$Top_Sequence) %>% 
        dplyr::select(Seq_ID, Genome_Build) %>% dplyr::distinct()
      prb_bad_cgn_cnt <- prd_bad_cgn_tib %>% base::nrow()
      
      # Check CGNs all have a Count == 1
      prd_cgn_cnt_tib <- prd_uniq_cgn_top_tib %>% dplyr::group_by(Seq_ID) %>% dplyr::summarise(Count=n())
      
      # This should be zero::
      prd_cgn_mis_tib <- prd_cgn_cnt_tib %>% dplyr::filter(Count != 1)
      
      cat(glue::glue("[{par$prgmTag}]: Done. Investigating Previous Product CGN -> Top_Sequence Designs; prb_bad_cgn_cnt={prb_bad_cgn_cnt}!{RET}{RET}") )
      
      #
      # Investigate uniqueness of probes::
      #
      imp_unq_tib <- imp_prd_tib %>% dplyr::distinct(Seq_ID, Methyl_Allele_CO_Strand, Methyl_Allele_TB_Strand, UnMethyl_Probe_Sequence)
      imp_mis_tib <- imp_unq_tib %>% dplyr::group_by(UnMethyl_Probe_Sequence) %>% dplyr::summarise(Seq_Count=n())
      imp_ann_tib <- imp_unq_tib %>% dplyr::left_join(imp_mis_tib, by="UnMethyl_Probe_Sequence")
      
      # imp_ann_tib %>% dplyr::filter(Seq_Count!=1, Methyl_Allele_CO_Strand=='C')
      imp_bad_prb_tib <- imp_ann_tib %>% dplyr::filter(Seq_Count!=1, Methyl_Allele_CO_Strand=='C') %>% dplyr::group_by(Seq_ID) %>% 
        dplyr::summarise(CGN_Count=n()) %>% dplyr::filter(CGN_Count!=1)
      prd_bad_prb_tib <- imp_prd_tib %>% dplyr::filter(Seq_ID %in% imp_bad_prb_tib$Seq_ID) %>% dplyr::distinct(Seq_ID, Genome_Build)
      prd_bad_prb_cnt <- prd_bad_prb_tib %>% base::nrow()
      
      cat(glue::glue("[{par$prgmTag}]: Done. Investigating Previous Product Non-Unique Prb_Seq_U; prd_bad_prb_cnt={prd_bad_prb_cnt}!{RET}{RET}") )
      
      #
      # Join all Failed Probes due to CGN->Top or UnMethyl_Probe_Sequeunce Integrity Issues::
      #
      prd_fin_bad_tib <- dplyr::bind_rows( prd_bad_cgn_tib, prd_bad_prb_tib ) %>% dplyr::distinct()
      prd_fin_bad_cnt <- prd_fin_bad_tib %>% base::nrow()
      
      cat(glue::glue("[{par$prgmTag}]: Done. Investigating Previous Product Integrity Failures; prd_fin_bad_cnt={prd_fin_bad_cnt}!{RET}{RET}") )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                Preprocessing:: Improbe Mouse Designs
      #
      #  3. Gather all mouse methylation improbe CGNs and designs
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # mm10_imp_AL_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top-seqU.tsv.gz'
      # mm10_imp_AL_tib <- readr::read_tsv(mm10_imp_AL_tsv)
      
      # Ensure we have the correct matching database between our designs, our design database (v1) and the full database (AL)
      mm10_imp_v1_tsv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/data/dropbox/merged_with_raw_ordered_withHeader.tsv.gz'
      mm10_imp_v1_tib <- suppressMessages(suppressWarnings( readr::read_tsv(mm10_imp_v1_tsv) )) %>% 
        dplyr::mutate(
          Seq_ID = stringr::str_remove_all(Seq_ID, ' '),
          Mat_Prb1=stringr::str_sub(UnMethyl_Probe_Sequence, 2,50), 
          Mat_Prb2=stringr::str_sub(UnMethyl_Probe_Sequence, 2,49) )
      
      # Non Traditional Probes::
      #  mm10_imp_v1_tib %>% dplyr::filter(Customer_Probe_IU!=UnMethyl_Probe_Sequence) 
      
      #
      # Make sure both V1 have CGN->TOP Integrity::
      #
      mm10_uniq_cgn_top_tib <- mm10_imp_v1_tib %>% dplyr::distinct(Seq_ID,Top_Sequence)
      mm10_uniq_top_cnt_tib <- mm10_uniq_cgn_top_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Top_Count=n())
      mm10_uniq_top_mis_cnt <- mm10_uniq_top_cnt_tib %>% dplyr::filter(Top_Count != 1) %>% base::nrow()
      
      cat(glue::glue("[{par$prgmTag}]: Done. Investigating Mouse Product Integrity Failures; ",
                     "mm10_uniq_top_mis_cnt={mm10_uniq_top_mis_cnt}!{RET}{RET}") )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                Preprocessing:: Improbe Mouse Designs
      #
      #  4. Comapre Previous Products with Mouse for CGN -> Top Integrity
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      prd_mm10_uniq_cgn_top_tib <- dplyr::bind_rows( prd_uniq_cgn_top_tib, mm10_uniq_cgn_top_tib ) %>% dplyr::distinct()
      prd_mm10_uniq_top_cnt_tib <- prd_mm10_uniq_cgn_top_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Top_Count=n())
      
      # Validation:: Only loss of integrity comes from 27k::
      prd_mm10_uniq_top_min_tib <- prd_mm10_uniq_top_cnt_tib %>% dplyr::filter(Top_Count != 1) %>% 
        dplyr::inner_join(imp_prd_tib %>% dplyr::select(Seq_ID,Top_Sequence,Genome_Build), by="Top_Sequence") %>% dplyr::distinct(Seq_ID, Genome_Build)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                Preprocessing:: Update mouse CGN names::
      #
      #  5. Use improbe mouse designs to generate new proper CGNs
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      #
      # LAST POINT BEFORE RESTART: TBD:
      #  - Validate designs for all probes on raw_man_tib
      #
      par$tbd_step <- FALSE
      if (par$tbd_step) {
        mm10_man_v1_matPrbCgn_tib <- dplyr::bind_rows(
          raw_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb1", "Seq_ID")),
          raw_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb2", "Seq_ID")) )
        
        mm10_man_v1_matPrb_tib <- dplyr::bind_rows(
          raw_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb1")),
          raw_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb2")) )
        
        # Can't Anit-join until we remove the non Seq_ID matches::
        # raw_man_tib %>% dplyr::anti_join(mm10_man_v1_matPrb_tib, by='')
        mm10_man_v1_matPrb_tib 
        
        
        # Probes with non-match Seq_ID's
        mm10_man_v1_matPrb_tib %>% dplyr::filter(Seq_ID.x != Seq_ID.y) %>% dplyr::select(Seq_ID.x, Seq_ID.y, Top_Sequence, Mat_Prb)
        
        mm10_man_v1_matPrb_tib %>% dplyr::filter(Seq_ID.x != Seq_ID.y) %>% dplyr::select(Seq_ID.x, Seq_ID.y) %>% as.data.frame()
        
      }
      
      #
      # NEXT STEP:: 
      #
      #  - Check missing probes from basic sort above
      #
      #  - Perl PrbSeqU_48 and sort by it from full design file.
      #  - Intersect by Mat_Prb
      #  - Checking missing probes
      #  - Identify RP/MU probes and add their values
      #
      mm10_des_AL_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top-seqU-to-sequ48.tsv.gz'
      mm10_des_AL_tib <- readr::read_tsv(mm10_des_AL_tsv)
      
      
      #
      # Full Unique:: 
      #
      par$full_cross_verification <- FALSE
      if (par$full_cross_verification) {
        #
        # Takes too much memory locally::
        #  mm10_imp_AL_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top-seqU.tsv.gz'
        #  mm10_imp_AL_tib <- readr::read_tsv(mm10_imp_AL_tsv)
        #
        # Next step: cross validate all CGN/TOP against all manifest designs...
        #  Unique Count = 20239851
        #
        mm10_imp_AL_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top.uniq.tsv'
        mm10_imp_AL_rds <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top.uniq.rds'
        # mm10_imp_AL_tib <- readr::read_tsv(mm10_imp_AL_tsv)
        # readr::write_rds(mm10_imp_AL_tib, mm10_imp_AL_rds,  compress="gz")
        mm10_imp_AL_tib <- readr::read_rds(mm10_imp_AL_rds)
        
        prd_mm10_full_cgn_top_tib <- dplyr::bind_rows( prd_uniq_cgn_top_tib, mm10_uniq_cgn_top_tib, mm10_imp_AL_tib ) %>% dplyr::distinct()
        prd_mm10_full_top_cnt_tib <- prd_mm10_full_cgn_top_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Top_Count=n())
        prd_mm10_full_top_cnt_csv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.prd_mm10_full_top_cnt.csv.gz'
        readr::write_csv(prd_mm10_full_top_cnt_tib, prd_mm10_full_top_cnt_csv)
        
        # Same ones as 27k::
        prd_mm10_full_top_cnt_tib %>% dplyr::filter(Top_Count != 1)
        
        # Conclusion: Cross Validation Validated!
      }
      
    }
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if (FALSE) {
    
    #
    # Cannot trust below becuase coordinate system is off for columns:: Also, the above seems to work...
    #
    if (FALSE) {
      
      #
      # TBD:: Repeat above with AL improbe mouse designs for integrity check...
      #  4. Comapre Previous Products with Mouse for CGN -> Top Integrity
      #
      
      #
      # Complete mouse improbe design::
      #
      mm10_imp_org_csv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/complicated-header.csv'
      mm10_imp_org_col <- suppressMessages(suppressWarnings( readr::read_csv(mm10_imp_org_csv) )) %>% names() %>% as.vector()
      # mm10_imp_org_col %>% length() = 92
      mm10_imp_ext_col <- c('Add_Forward_Sequence','Add_Genome_Build','Add_Chromosome','Add_Position','Add_Top_Sequence',
                            'Add_Probe_Seq_M','Add_FR','Add_TB','Add_CO','Add_Score_M','Add_CpgCnt','Add_CpgPos','Add_Probe_Seq_U',
                            'Add_Score_U')
      mm10_imp_AL_col <- c(mm10_imp_org_col,mm10_imp_ext_col)
      # mm10_imp_AL_col %>% length() = 106
      
      # mm10_imp_AL_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests/joined.designed.dat.head.csv'
      mm10_imp_AL_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests/joined.designed.dat.csv.gz'
      mm10_imp_AL_tib <- suppressMessages(suppressWarnings( readr::read_csv(mm10_imp_AL_csv, col_names=mm10_imp_AL_col) )) %>% 
        dplyr::mutate(
          Mat_Prb1=stringr::str_sub(Add_Probe_Seq_U, 2,50), 
          Mat_Prb2=stringr::str_sub(Add_Probe_Seq_U, 2,49) )
    }
    
    # Join manifest probes and improbe based on Allele_A probe
    if (FALSE) {
      # v1 improbe::
      mm10_man_v1_mat_tib <- dplyr::bind_rows(
        raw_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb1")),
        raw_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb2")) )
      
      # AL improbe::
      mm10_man_AL_mat_tib <- dplyr::bind_rows(
        raw_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb1")),
        raw_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb2")) )
      
      
      #
      # Check improbe design uniqueness by CGN->TOP Integrity::
      #
      mm10_imp_v1_tib %>% base::nrow() # 285604
      mm10_imp_AL_tib %>% base::nrow() # 1108632
      
      mm10_unq_v1_cgTop_tib <- mm10_imp_v1_tib %>% dplyr::distinct(Seq_ID, Top_Sequence)
      mm10_unq_AL_cgTop_tib <- mm10_imp_AL_tib %>% dplyr::distinct(Probe_ID, Top_Sequence) %>% dplyr::rename(Seq_ID=Probe_ID)
      
      # Failed in translation: Probe_ID=='cg00554227'
      
      #
      # By the numbers::
      #
      raw_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb1")) %>% base::nrow()
      raw_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb1")) %>% base::nrow()
      
      raw_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb2")) %>% base::nrow()
      raw_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb2")) %>% base::nrow()
      
      
      #
      # Probe IDs that don't match::
      #
      raw_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb1")) %>% dplyr::filter(Seq_ID.x != Seq_ID.y) %>% as.data.frame()
      raw_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb2")) %>% dplyr::filter(Seq_ID.x != Seq_ID.y) %>% as.data.frame()
      
      
    }
    
    #
    # Left off checking the other mouse designs here::
    #
    
    # improbe_AL_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.core-v1.tsv.gz'
    # imp_des_AL_tib <- readr::read_tsv(improbe_AL_tsv)
    # improbe_AL_cols <- c('Probe_ID', 'Chrom1', 'Pos1', 'Pos2', 'Chip_Str1', 
    #                      'Customer_Score', 'Customer_Underlying_CpG_Count', 'Customer_Pos', 'Customer_Design_Char', 'Customer_Dist1', 'Customer_Dist2',
    #                      'Customer_Long_ID', 'Customer_ID_CO', 'Customer_Prb_A', 'Customer_Prb_B', 'Customer_Design_Type'
    # )
    
    mm10_imp_org_csv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/complicated-header.csv'
    mm10_imp_org_col <- suppressMessages(suppressWarnings( readr::read_csv(mm10_imp_org_csv) )) %>% names() %>% as.vector()
    # mm10_imp_org_col %>% length() = 92
    mm10_imp_ext_col <- c('Add_Forward_Sequence','Add_Genome_Build','Add_Chromosome','Add_Position','Add_Top_Sequence',
                          'Add_Probe_Seq_M','Add_FR','Add_TB','Add_CO','Add_Score_M','Add_CpgCnt','Add_CpgPos','Add_Probe_Seq_U',
                          'Add_Score_U')
    
    mm10_imp_AL_col <- c(mm10_imp_org_col,mm10_imp_ext_col)
    # mm10_imp_AL_col %>% length() = 106
    
    # mm10_imp_AL_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests/joined.designed.dat.head.csv'
    mm10_imp_AL_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests/joined.designed.dat.csv.gz'
    mm10_imp_AL_tib <- suppressMessages(suppressWarnings( readr::read_csv(mm10_imp_AL_csv, col_names=mm10_imp_AL_col) )) %>% 
      dplyr::mutate(
        Mat_Prb1=stringr::str_sub(UnMethyl_Probe_Sequence, 2,50), 
        Mat_Prb2=stringr::str_sub(UnMethyl_Probe_Sequence, 2,49) )
    
    # Join manifest probes and improbe based on Allele_A probe
    # raw_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb1"))
    # raw_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb2"))
    # raw_man_tib %>% dplyr::full_join(mm10_imp_AL_tib, by=c("Seq_ID"="Probe_ID") ) %>% 
    #   dplyr::select(Seq_ID, PD, Address_A_Seq, Mat_Prb1, Mat_Prb2)
    
    
    
    #
    # Make sure both V1/AL have CGN->TOP Integrity::
    #
    unq_imp_v1_tib <- imp_des_v1_tib %>% dplyr::distinct(Seq_ID,Top_Sequence)
    top_cnt_v1_tib <- unq_imp_v1_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Top_Count=n())
    
    unq_imp_AL_tib <- improbe_AL_tib %>% dplyr::distinct(Probe_ID,Top_Sequence)
    top_cnt_AL_tib <- unq_imp_AL_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Top_Count=n())
    
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                Preprocessing:: Improbe Mouse Designs
    #
    #  4. Comapre Previous Products with Mouse for CGN -> Top Integrity
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #              Dirty and Quick Selection of Optimal Content::
    #
    #                           EXAMPLE CODE BELOW::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (FALSE) {
      man_cn_des_all_tib <- dplyr::left_join(man_sel_cn_tib,
                                             imp_des_tib %>% dplyr::filter(Methyl_Allele_CO_Strand=='C') %>% dplyr::select(-Probe_Type),
                                             by=c("IlmnID"="Seq_ID"), suffix=c("_Man", "_Des")
      ) %>% dplyr::mutate(
        Mat_Prb1_Man=stringr::str_sub(AlleleA_ProbeSeq, 2,50) %>% stringr::str_replace_all('R', 'A'), 
        Mat_Prb2_Man=stringr::str_sub(AlleleA_ProbeSeq, 3,50) %>% stringr::str_replace_all('R', 'A'), 
        Mat_Prb1_Des=stringr::str_sub(UnMethyl_Probe_Sequence, 2,50), 
        Mat_Prb2_Des=stringr::str_sub(UnMethyl_Probe_Sequence, 2,49) ) %>% dplyr::filter(! IlmnID %in% man_27k_tib$IlmnID)
      
      mat_des1_tib <- man_cn_des_all_tib %>% dplyr::filter(Mat_Prb1_Man==Mat_Prb1_Des & Forward_Sequence_Man == Forward_Sequence_Des)
      mat_des2_tib <- man_cn_des_all_tib %>% dplyr::filter(Mat_Prb2_Man==Mat_Prb2_Des & Forward_Sequence_Man == Forward_Sequence_Des)
      
      man_cn_des_mat_tib <- dplyr::bind_rows( mat_des1_tib,mat_des2_tib ) %>% 
        dplyr::arrange(IlmnID) %>% dplyr::rename(Forward_Sequence=Forward_Sequence_Des) %>% 
        dplyr::select(-Forward_Sequence_Man)
      
      man_cn_des_mat_tib %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise(Count=n()) %>% print()
    }
    
    
  }
  
  
  
  
  
  
  
  
  
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Simplified Output Manifest::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  out_man_tib <- raw_man_tib %>% 
    dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, 
                  AlleleA_Probe_Sequence,AlleleB_Probe_Sequence) %>%
    dplyr::mutate(Version=opt$tar_version)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Extract Probe CGN Positions for cg & mu::
  #                              NOT USED YET!
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (FALSE) {
    cgn_key_tib <- out_man_tib %>% dplyr::filter(Probe_Type=='mu' | Probe_Type=='cg') %>% 
      dplyr::mutate(CGN=stringr::str_remove(Probe_ID,'_.*$')) %>% 
      dplyr::mutate(CGN=stringr::str_replace(CGN,'^mu','cg')) %>% 
      dplyr::select(CGN,Probe_ID,Probe_Type) %>% dplyr::arrange(CGN)
    
    cgn_key_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Old Control Probe Loading::
  #                              NOT USED YET!
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (FALSE) {
    opt$manifestPath <- file.path(opt$manDir, 'HM450-B2.manifest.sesame-base.cpg-sorted.csv.gz')
    if (is.null(opt$manifestPath))
      opt$manifestPath <- file.path(opt$manDir, paste0(opt$platform,'-',opt$manifest,'.manifest.sesame-base.cpg-sorted.csv.gz') )
    
    full_org_man_tib <- loadManifestSource(opt$manifestPath, verbose=opt$verbose,tt=NULL)
    org_man_tib <- loadManifestSource(opt$manifestPath, verbose=opt$verbose,tt=NULL) %>%
      dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base) %>% 
      dplyr::mutate(Version=opt$manifest)
    org_ctl_tib <- org_man_tib %>% dplyr::filter(Probe_Type!='cg', Probe_Type!='ch', Probe_Type!='rs') %>% dplyr::arrange(Probe_ID)
    
    if (opt$addControls) out_man_tib <- dplyr::bind_rows(out_man_tib, org_ctl_tib) %>% dplyr::arrange(Probe_ID, DESIGN)
    if (opt$addManifest) out_man_tib <- dplyr::bind_rows(out_man_tib, org_man_tib) %>% dplyr::arrange(Probe_ID, DESIGN)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Sesame Manifest::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (FALSE) {
    ses_neg_tib <- ses_ctl_tib %>% dplyr::filter(Probe_Type=='NEGATIVE')
    # ses_man_tib <- dplyr::bind_rows(out_man_tib) %>% 
    # ses_man_tib <- dplyr::bind_rows(out_man_tib,ses_ctl_tib) %>% 
    ses_man_tib <- dplyr::bind_rows(out_man_tib,ses_neg_tib) %>% 
      # dplyr::select(-AlleleA_Probe_Sequence,-AlleleB_Probe_Sequence) %>%
      dplyr::arrange(Probe_ID)
    readr::write_csv(ses_man_tib, ses_man_csv)
    
    # Local Mac Copy Command For Testing::
    #  cp /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-cp/mm10_LEGX_cp.manifest.sesame-base.cpg-sorted.csv.gz tools/Infinium_Methylation_Workhorse/dat/manifest/base/LEGX-B0..manifest.sesame-base.cpg-sorted.csv.gz
    #
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Genome Studio Manifest::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  loci_cnt <- out_man_tib %>% base::nrow()
  gs_head_tib <- tibble::tibble(Key=character(), Value=character())
  gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="Illumina, Inc.", Value="")
  gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="[Heading]", Value="")
  gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="Descriptor File Name", Value=paste0('Infinium_Methylation_',opt$runName,'.csv.gz'))
  gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="Assay Format", Value="Infinium 2")
  gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="Date Manufactured", Value=as.character(Sys.Date()))
  gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="Loci Count", Value=as.character(loci_cnt))
  gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="[Assay]", Value="")
  
  #  IlmnID,Name,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,
  #   Infinium_Design_Type,Next_Base,Color_Channel,Forward_Sequence,
  #   Genome_Build,CHR,MAPINFO,SourceSeq,Strand
  #
  #   Missing: Source_Seq
  #
  src_seq <- as.character(paste(rep('N',50), collapse='') )
  gs_body_tib <- out_man_tib %>% # head(n=1000) %>%
    dplyr::rename(IlmnID=Probe_ID,
                  AddressA_ID=U,AlleleA_ProbeSeq=AlleleA_Probe_Sequence,
                  AddressB_ID=M,AlleleB_ProbeSeq=AlleleB_Probe_Sequence,
                  Infinium_Design_Type=DESIGN,Color_Channel=col) %>%
    dplyr::mutate(Name=IlmnID,Forward_Sequence='N',Genome_Build=opt$genomeBuild,
                  CHR='chr1',MAPINFO=1,Strand='F',Source_Seq=!!src_seq) %>%
    dplyr::select(IlmnID,Name,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,
                  Infinium_Design_Type,Next_Base,Color_Channel,Forward_Sequence,
                  Genome_Build,CHR,MAPINFO,Source_Seq,Strand,Source_Seq,Probe_Type,Probe_Source) %>%
    dplyr::arrange(IlmnID)
  ctl_head_df <- tibble::tibble(Head="[Controls]")
  
  # Silly addition of commas to the end of the controls section...
  #
  org_ctl_tib <- org_ctl_tib %>% dplyr::mutate(Probe_ID=paste0(Probe_ID,",,,,,,"))
  
  # Silly incosistent padding of Tango Addresses..
  gs_body_tib <- gs_body_tib %>% dplyr::mutate(
    AddressA_ID=stringr::str_pad(string=AddressA_ID, pad='0', side="left", width=10), 
    AddressB_ID=stringr::str_pad(string=AddressB_ID, pad='0', side="left", width=10),
    Color_Channel=dplyr::case_when(
      Color_Channel=='R' ~ 'Red',
      Color_Channel=='G' ~ 'Grn',
      TRUE ~ ''
    )
  )
  
  # Write Genome Studio CSV::
  write_delim(x=gs_head_tib, path=gss_man_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=FALSE)
  write_delim(x=gs_body_tib, path=gss_man_csv, delim=',', col_names=TRUE,  quote_escape=FALSE, na='', append=TRUE)
  write_delim(x=ctl_head_df, path=gss_man_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)
  write_delim(x=org_ctl_tib, path=gss_man_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)
  
  # Follow up silly command to remove quotes::
  #  gzip -dc /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-B4/mm10_LEGX_B4.manifest.GenomeStudio.cpg-sorted.csv.gz | perl -pe 's/\"//gi' | gzip -c - > /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-B4/mm10_LEGX_B4.manifest.GenomeStudio.cpg-sorted.clean.csv.gz
  
  # gss_man_zip <- paste(gss_man_csv,'gz', sep='.')
  # cmd <- glue::glue("cat {gss_man_csv} | perl -pe 's/\"//gi' | gzip -c - > {gss_man_zip}")
  # system(cmd)
  
  # Follow up silly command to remove quotes::
  #  cat /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-B5/mm10_LEGX_B5.manifest.GenomeStudio.cpg-sorted.csv | perl -pe 's/\"//gi' | gzip -c - > /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-B5/mm10_LEGX_B5.manifest.GenomeStudio.cpg-sorted.clean.csv.gz
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Write mm10 Manifest::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (opt$make_addresss) {
    out_add_tib <- NULL
    out_add_tib <- dplyr::bind_rows(
      out_man_tib %>% dplyr::filter(!is.na(M)) %>% dplyr::select(Probe_ID, M, col, DESIGN, Probe_Type) %>% dplyr::rename(Address=M),
      out_man_tib %>% dplyr::filter(!is.na(U)) %>% dplyr::select(Probe_ID, U, col, DESIGN, Probe_Type) %>% dplyr::rename(Address=U)
    ) %>% dplyr::rename(Man_Col=col, Design_Type=DESIGN)
    
    readr::write_csv(out_add_tib, out_add_csv)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Load Annotation/Alignment Results::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (FALSE) {
    ann_cnames <- c('Probe_ID', 'Probe_ID_A', 'Probe_Seq_A', 'Probe_ID_B', 'Probe_Seq_B', 
                    'Normalization_Bin', 'Probe_Chrom', 'Probe_Beg', 'Probe_End', 'Strand_FR', 'Probe_Cnt',
                    'Reg_Type', 'Transcript', 'Gene_Name', 'Reg_Beg', 'Reg_End', 'Reg_Idx', 
                    'Gene_Chrom', 'Gene_Beg', 'Gene_End', 'Gene_Score')
    ann1_tsv <- file.path(opt$alnDir, 'COVID_EPIC_Round1.COVID-Genes.intersection.tsv')
    ann1_tib <- readr::read_tsv(file=ann1_tsv, col_names=ann_cnames)
    
    gene_ann_tib <- ann1_tib %>% dplyr::rename(Aln_Probe_ID=Probe_ID) %>% 
      dplyr::mutate(Probe_ID=stringr::str_remove(Probe_ID_A, '_A$')) %>% 
      tidyr::separate(Probe_ID, into=c("CG","FR","TB","CO","DT","Hits"), sep='_', remove=FALSE) %>%
      dplyr::select(Probe_ID, CG, FR, TB, CO, DT, Hits, everything())
    
    hits_ann_tib <- gene_ann_tib %>% dplyr::select(Probe_ID:DT, Probe_Chrom, Probe_Beg, Probe_End, Strand_FR) %>% dplyr::distinct() %>% 
      dplyr::add_count(Probe_ID, name='Aln_Count') %>% dplyr::arrange(-Aln_Count) %>% dplyr::distinct(Probe_ID, Aln_Count)
    
    gene_man_tib <- raw_man_tib %>% dplyr::inner_join(gene_ann_tib, by="Probe_ID", suffix=c("_Man", "_Ann"))
    misA_man_tib <- raw_man_tib %>% dplyr::anti_join(gene_ann_tib, by="Probe_ID", suffix=c("_Man", "_Ann"))
    misB_man_tib <- gene_ann_tib %>% dplyr::anti_join(raw_man_tib, by="Probe_ID", suffix=c("_Man", "_Ann"))
    
    # raw_man_tib %>% dplyr::filter(stringr::str_starts( 'cg00065388' ) )
    
    
    #
    # Quick look up for Fasial for chr18 HLA region coverage "ESTIMATE"
    #
    #  hg19: chr5 140164938 140864474
    #  mm10: chr18 36930065 37815275
    #
    #
    # mm_tib <- readr::read_tsv('/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/to-order/designW.Final_Dec-10-2019.cpg-only/Mus_musculus.annotation.genomic.sorted-cgn.chr-pos.tsv', col_names = c("Probe_ID",'Pos'))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Genome Studio Controls Swap::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    color_tib <- suppressMessages(suppressWarnings( readr::read_csv(file.path(par$datDir, 'params/GenomeStudioColors.csv')) ))
    color_vec <- color_tib %>% dplyr::pull(Color) %>% as.vector()
    color_len <- length(color_vec)
    
    man_gs_csv  <- file.path(par$topDir, 'data/manifests/tmp/mm10_LEGX_B13.manifest.GenomeStudio.cpg-sorted.clean.csv.gz')
    man_head_df <- readr::read_lines(man_gs_csv, n_max = 6) %>% as.data.frame()
    
    man_gs_csv  <- file.path(par$topDir, 'tmp/LifeEpigentics/manifests/mm10-LEGX-B5/mm10_LEGX_B5.manifest.GenomeStudio.cpg-sorted.clean.csv.gz')
    man_gs_list <- loadManifestGenomeStudio(file = man_gs_csv, addSource = TRUE, normalize = FALSE, verbose = 20)
    ctl_head_df <- tibble::tibble(Head="[Controls]") %>% tibble::add_column(BL1='',BL2='',BL3='',BL4='') %>% as.data.frame()
    
    man_gs_list$ctl %>% dplyr::group_by(Control_Group) %>% dplyr::summarise(Count=n())
    man_gs_list$man %>% dplyr::group_by(Probe_Type, Infinium_Design_Type) %>% dplyr::summarise(Type_Count=n())
    
    man_ct_tib <- man_gs_list$man %>% dplyr::filter(Probe_Type=='BS' | Probe_Type=='NO')
    man_gs_tib <- man_gs_list$man %>% dplyr::filter(Probe_Type!='BS') %>% dplyr::filter(Probe_Type!='NO')
    ctl_gs_tib <- man_gs_list$ctl
    
    ctl_cols <- ctl_gs_tib %>% dplyr::select(1:4) %>% names()
    
    
    
    
    
    
    ctl_gs_tib <- dplyr::bind_rows(
      ctl_fin_tib <- fin_core_split[['Control']] %>% dplyr::arrange(Probe_ID) %>%
        dplyr::distinct(M,U, .keep_all=TRUE) %>% dplyr::select(U,Control_Group,Probe_ID) %>%
        dplyr::filter(!is.na(U)),
      ctl_fin_tib <- fin_core_split[['Control']] %>% dplyr::arrange(Probe_ID) %>%
        dplyr::distinct(M,U, .keep_all=TRUE) %>% dplyr::select(U,Control_Group,Probe_ID) %>%
        dplyr::filter(!is.na(M))
    ) dplyr::arrange(Control_Group)
    
    
    # Naming investigation::
    ctl_gs_tib %>% dplyr::filter(Control_Group=="BISULFITE CONVERSION I") %>% dplyr::arrange(Control_Name) %>% dplyr::select(-Control_Color, -Man_Source) %>% as.data.frame()  
    ctl_gs_tib %>% dplyr::filter(Control_Group=="BISULFITE CONVERSION II") %>% dplyr::arrange(Control_Name) %>% dplyr::select(-Control_Color, -Man_Source) %>% as.data.frame()
    ctl_gs_tib %>% dplyr::filter(Control_Group=="NON-POLYMORPHIC") %>% dplyr::arrange(Control_Name) %>% dplyr::select(-Control_Color, -Man_Source) %>% as.data.frame()
    
    # Make new controls::
    new_ct_list <- man_ct_tib %>% dplyr::arrange(Probe_Type) %>%
      dplyr::mutate(
        Control_Group=dplyr::case_when(
          Probe_Type=='BS' & Infinium_Design_Type=='I'  ~ "BISULFITE CONVERSION I",
          Probe_Type=='BS' & Infinium_Design_Type=='II' ~ "BISULFITE CONVERSION II",
          Probe_Type=='NO' ~ 'NON-POLYMORPHIC',
          TRUE ~ NA_character_
        ),
        DiNuc=dplyr::case_when(
          Probe_Type=='BS' ~ stringr::str_replace(IlmnID, '^.*-([ACTG][ACTG])-.*$', '\\$1') %>% stringr::str_remove_all('\\\\'),
          Probe_Type=='NO' ~ stringr::str_replace(IlmnID, '^.*_([ACTG][ACTG])_.*$', '\\$1') %>% stringr::str_remove_all('\\\\'),
          TRUE ~ NA_character_
        ),
        N1=stringr::str_sub(DiNuc, 1,1),
        N2=stringr::str_sub(DiNuc, 2,2),
        
        Control_ID=paste(IlmnID, stringr::str_sub(Color_Channel, 1,1), sep='_')
      ) %>% split(.$Infinium_Design_Type)
    # dplyr::select(IlmnID,DiNuc,N1,N2,Control_Group,Control_ID,Color_Channel,Probe_Type,Infinium_Design_Type)
    
    # Notes:: 
    #  BISULFITE CONVERSION I; All M's get 'SkyBlue' the C's get a unique color
    #     N2.B == A(G) => C
    #     N2.B == T(R) -> U
    #     N1.A == C    -> M
    #
    bs1_both_tib <- new_ct_list$I %>% dplyr::mutate(
      Row_Idx=dplyr::row_number() + 100,
      Row_Str=Row_Idx,
      Col_Idx=Row_Idx %% (color_len-100),
      
      Control_Name_A=paste0('BS Conversion I M',Row_Str),
      Control_Name_B=dplyr::case_when(
        Color_Channel=='Grn' ~ paste0('BS Conversion I-C',Row_Str),
        Color_Channel=='Red' ~ paste0('BS Conversion I-U',Row_Str),
        TRUE ~ NA_character_
      ),
      Control_ColorA='Red',
      Control_ColorB=color_vec[Col_Idx+1]
    )
    
    # %>% dplyr::select(Control_Name_A,Control_Name_B,
    #                     IlmnID,DiNuc,N1,N2,Control_Group,Row,Control_ID,
    #                     Color_Channel,Probe_Type,Infinium_Design_Type) %>% dplyr::select(1,2) %>% as.data.frame()
    
    fin_bs1_tib <- dplyr::bind_rows(
      bs1_both_tib %>% dplyr::select(AddressA_ID,Control_Group,Control_ColorA,Control_Name_A) %>% purrr::set_names(ctl_cols),
      bs1_both_tib %>% dplyr::select(AddressB_ID,Control_Group,Control_ColorB,Control_Name_B) %>% purrr::set_names(ctl_cols)
    ) %>% dplyr::arrange(Control_Name)
    
    fin_bs2_tib <- new_ct_list$II %>% 
      dplyr::filter(Control_Group=='BISULFITE CONVERSION II') %>%
      dplyr::mutate(Row_Idx=row_number()+100,
                    Row_Str=Row_Idx,
                    Control_Color='Blue',
                    Control_Name=paste('BS Conversion II',Row_Str, sep='-')
      ) %>% 
      dplyr::select(AddressA_ID,Control_Group,Control_Color,Control_Name) %>% 
      purrr::set_names(ctl_cols) %>%
      dplyr::arrange(Control_Name)
    
    fin_non_tib <- new_ct_list$II %>% 
      dplyr::filter(Control_Group=='NON-POLYMORPHIC') %>%
      dplyr::mutate(Row_Idx=row_number()+100,
                    Row_Str=Row_Idx,
                    Control_Color='Blue',
                    Control_Name=paste0('NP (',N1,') ',Row_Str)
      ) %>% 
      dplyr::select(AddressA_ID,Control_Group,Control_Color,Control_Name) %>% 
      purrr::set_names(ctl_cols) %>%
      dplyr::arrange(Control_Name)
    
    fin_hum_tib <- ctl_gs_tib %>% dplyr::select(1:4) %>% 
      purrr::set_names(ctl_cols) %>% dplyr::arrange(Control_Name) %>%
      dplyr::mutate(Address=as.double(Address))
    
    fin_ctl_tib <- dplyr::bind_rows(fin_bs1_tib,fin_bs2_tib,fin_non_tib,fin_hum_tib)
    
    #
    # Investigat Duplicates::
    #
    add_cnt_tib <- dplyr::bind_rows(
      dplyr::select(man_gs_tib,AddressA_ID) %>% dplyr::rename(Address=AddressA_ID), 
      dplyr::select(man_gs_tib,AddressB_ID) %>% dplyr::rename(Address=AddressB_ID),
      dplyr::select(fin_ctl_tib,Address)
    ) %>% dplyr::filter(!is.na(Address)) %>%
      dplyr::group_by(Address) %>% dplyr::summarise(Add_Dup_Count=n()) %>% dplyr::arrange(-Add_Dup_Count)
    
    add_tib_tib <- add_cnt_tib %>% dplyr::filter(Add_Dup_Count>1)
    
    man_fin_tib <- man_gs_tib %>% dplyr::filter(!AddressA_ID %in% add_tib_tib$Address) %>% dplyr::filter(!AddressB_ID %in% add_tib_tib$Address)
    ctl_fin_tib <- fin_ctl_tib %>% dplyr::filter(! Address %in% add_tib_tib$Address) %>% 
      tibble::add_column(BL1='',BL2='',BL3='',BL4='',BL5='',BL6='')
    # tibble::add_column(BL1='',BL2='',BL3='',BL4='',BL5='',BL6='',BL7='',BL8='')
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Write Control Swapped Manifest::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    gs_man_ver  <- 8
    gs_man_ver  <- 9
    gs_swap_dir <- file.path('/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests', paste0('mm10-LEGX-B',gs_man_ver))
    gs_swap_csv <- file.path(gs_swap_dir, paste0('mm10_LEGX_B',gs_man_ver,'.manifest.GenomeStudio.cpg-sorted.clean.csv') )
    gz_swap_csv <- file.path(gs_swap_dir, paste0('mm10_LEGX_B',gs_man_ver,'.manifest.GenomeStudio.cpg-sorted.clean.csv.gz') )
    
    if (!dir.exists(gs_swap_dir)) dir.create(gs_swap_dir, recursive=TRUE)
    
    man_head_replace_str <- paste0('_B',gs_man_ver,'.csv.gz')
    man_head_df <- man_head_df %>% as_tibble() %>% purrr::set_names(c("Head_Var")) %>% 
      dplyr::mutate(Head_Var=as.character(Head_Var),
                    Head_Var=stringr::str_replace(Head_Var, '_B[0-9]+.csv.gz$', man_head_replace_str) ) %>%
      as.data.frame()
    
    # Write Genome Studio CSV::
    write_delim(x=man_head_df, path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=FALSE)
    write_delim(x=man_fin_tib, path=gs_swap_csv, delim=',', col_names=TRUE,  quote_escape=FALSE, na='', append=TRUE)
    write_delim(x=ctl_head_df, path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)
    write_delim(x=ctl_fin_tib, path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)
    
    cmd <- glue::glue("cat {gs_swap_csv} | perl -pe 's/\"//gi' | gzip -c - > {gz_swap_csv}")
    # system(cmd)
  }
}




# End of file
