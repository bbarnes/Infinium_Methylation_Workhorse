
#
# Validation NOTES::
#
#
# Validate Against Order File::
#
old_man_csv <- "/Users/bretbarnes/Documents/data/manifests/methylation/Chicago-Ober-Custom.original/Chicago-S38.manifest.sesame-base.cpg-sorted.csv"
old_man_tib <- safe_read(old_man_csv)

ses_man_tib %>% dplyr::filter(!is.na(U) & ! U %in% old_man_tib$U) %>% as.data.frame()
ses_man_tib %>% dplyr::filter(!is.na(M) & ! M %in% old_man_tib$M) %>% as.data.frame()

ses_mis_tib <- dplyr::bind_rows(
  old_man_tib %>% dplyr::filter(!is.na(U) & ! U %in% ses_man_tib$U),
  old_man_tib %>% dplyr::filter(!is.na(M) & ! M %in% ses_man_tib$M)
) %>% dplyr::distinct() %>% 
  dplyr::mutate(
    Infinium_Design_Type=dplyr::case_when(
      Infinium_Design_Type==1 ~ "I",
      Infinium_Design_Type==2 ~ "II",
      TRUE ~ NA_character_)
  ) %>% clean_tibble()


#
# Validation that all the unpaired probes have failing partners!!!
#
if (FALSE) {
  # These seem to be singletons::
  #
  misM_tib <- aqp_add_tib %>% 
    dplyr::filter(Ord_Des=="M") %>% 
    dplyr::filter(! Address %in% ses_man_tib$M)
  misU_tib <- aqp_add_tib %>% 
    dplyr::filter(Ord_Des!="M") %>% 
    dplyr::filter(! Address %in% ses_man_tib$U)
  
  # No Overlap::
  #
  misM_tib %>% dplyr::inner_join(misU_tib, by="Ord_Cgn")
  
  mat_vec  <- splitStrToVec(opt$mats)
  mats_tib <- dplyr::bind_rows(load_aqp_files(mat_vec[1]),load_aqp_files(mat_vec[2]))
  
  mat_mis1_tib <- dplyr::bind_rows(
    mats_tib %>% dplyr::filter(Mat_Prb %in% misM_tib$Ord_Prb),
    mats_tib %>% dplyr::filter(Mat_Prb %in% misU_tib$Ord_Prb)
  )
  pqc_tib %>% dplyr::filter(Address %in% mat_mis1_tib$Address) %>% print(n=10000)
  
  # The Pairs all failed!!!
  mat_mis2_tib <- dplyr::bind_rows(
    mats_tib %>% dplyr::filter(Mat_Prb %in% misM_tib$Ord_Par),
    mats_tib %>% dplyr::filter(Mat_Prb %in% misU_tib$Ord_Par)
  )
  pqc_tib %>% dplyr::filter(Address %in% mat_mis2_tib$Address) %>% print(n=10000)
  
  pass_mis_csv <- file.path(par$topDir, "data/CustomContent/McMaster/manifest/passing-but-missing.10072021.txt")
  pass_mis_tib <- readr::read_csv(pass_mis_csv)
  
  mis1_tib <- aqp_add_tib %>% 
    dplyr::mutate(Ord_Key2=stringr::str_remove(Ord_Key,"[-_:].*$")) %>% 
    dplyr::filter(Ord_Key2 %in% pass_mis_tib$Probe_ID) %>% 
    dplyr::distinct(Ord_Key, .keep_all = TRUE)
  
  ord_tib <- load_aqp_file( opt$ords )
  mis2_tib <- ord_tib %>% 
    dplyr::filter(Ord_Key %in% pass_mis_tib$Probe_ID) %>% 
    dplyr::distinct(Ord_Key, .keep_all = TRUE)
  
  mis2_tib %>% dplyr::filter(Ord_Key %in% mis1_tib$Ord_Key)
  mis1_tib %>% dplyr::filter(Ord_Key %in% mis2_tib$Ord_Key)
  
  aqp_add_tib %>% dplyr::filter(Ord_Prb %in% mis2_tib$Ord_Prb)
  
  mis2_tib %>% dplyr::inner_join(mats_tib, by=c("Ord_Prb"="Mat_Prb"))
  
  load_aqp_files(mat_vec)
  
}


#
# Mapping to original coordinates with Genomic Regions::
#
imp_trim_tib <- NULL
if (FALSE) {
  
  pos_cols <- 
    cols(
      Cgn    = col_integer(),
      Chr    = col_character(),
      Pos    = col_integer(),
      PrbA   = col_character(),
      ExtA   = col_character(),
      Srd_FR = col_character(),
      Srd_TB = col_character(),
      Srd_CO = col_character(),
      Nxb    = col_character(),
      PrbB   = col_character(),
      ExtB   = col_character()
    )
  
  # ord_pos_tib <- safe_read(par$ord_pos_csv, verbose=opt$verbose)
  ord_pos_tib <- readr::read_tsv(par$ord_pos_csv, 
                                 col_names = names(pos_cols$cols), 
                                 col_types=pos_cols) %>%
    dplyr::mutate(Chr=paste0("chr",Chr)) %>% 
    dplyr::distinct(Chr,Pos, .keep_all = TRUE) %>%
    dplyr::mutate(Rank=dplyr::row_number(),
                  Unq_Key=paste(Chr,Pos,Rank,sep="_"))
  
  buf_len <- 50
  ord_grs =
    GenomicRanges::GRanges(
      seqnames=Rle(ord_pos_tib$Chr),
      # strand=Rle(ret_tib$srd),
      
      PrbA=ord_pos_tib$PrbA,
      
      IRanges(start=ord_pos_tib$Pos-buf_len,
              end=ord_pos_tib$Pos+buf_len,
              names=ord_pos_tib$Unq_Key)
    )
  
  imp_des_list <- aqp_imp_tib %>% split(.$Ord_Des)
  
  imp_grs =
    GenomicRanges::GRanges(
      seqnames=Rle(aqp_imp_tib$Bsp_Chr),
      # strand=Rle(ret_tib$srd),
      
      Cgn_Str=aqp_imp_tib$Cgn_Str,
      Aln_Key=aqp_imp_tib$Aln_Key,
      Ord_Key=aqp_imp_tib$Ord_Key,
      Cpg_Pos=aqp_imp_tib$Bsp_Pos,
      
      IRanges(start=aqp_imp_tib$Bsp_Pos,
              width = 2,
              names=aqp_imp_tib$Aln_Key_Unq)
    )
  
  imp_int_tib <- intersect_GRS(ord_grs,imp_grs, 
                               can_prefix="Imp", 
                               can_key="Unq_Key", 
                               verbose=opt$verbose)
  
  imp_trim_tib <- aqp_imp_tib %>% dplyr::filter(Aln_Key %in% imp_int_tib$Imp_Aln_Key)
  
  imp_des_list <- NULL
  if (!is.null(imp_trim_tib)) {
    imp_des_list <- imp_trim_tib %>% split(.$Ord_Des)
  } else {
    imp_des_list <- aqp_imp_tib %>% split(.$Ord_Des)
  }
}



#
# Current thoughts are true joining:: Two ways::
#   A: by coordinate and orientation
#      - TBD: Should add coordinate cgn look up!!
#      - TBD: Use the extension distribution some how
#   C: by cgn in case alignment failed
#
# This should all be tested with a small multi-unique set like NZT!
#
# Both bsp_cgn_tib and aqp_bsp_tib need to be joined seperately
#  with aqp_add_tib
#
if (FALSE) {
  #   A: by coordinate and orientation
  aqp_bsp_tibA <- aqp_bsp_tib %>% 
    dplyr::select(Address,Bsp_Chr,Bsp_Pos,Bsp_Din_Bsc,Bsp_Nxb_Bsc,
                  Bsp_Tag,Bsp_Srd,Ord_Key) %>%
    dplyr::arrange(Bsp_Chr,Bsp_Pos)
  
  aqp_bsp_tibB <- aqp_bsp_tib %>% 
    dplyr::select(Address,Bsp_Chr,Bsp_Pos,Bsp_Din_Bsc,Bsp_Nxb_Bsc,
                  Bsp_Tag,Bsp_Srd,Ord_Key) %>%
    dplyr::arrange(Bsp_Chr,Bsp_Pos)
  
  aqp_bsp_inn <- aqp_bsp_tibA %>% 
    dplyr::left_join(aqp_bsp_tibB, by=c("Bsp_Chr","Bsp_Pos"), suffix=c("_A", "_B")) %>% 
    dplyr::distinct() %>%
    dplyr::filter(Address_A!=Address_B)
  
  
  aqp_bsp_pas <- dplyr::bind_rows(
    aqp_bsp_inn %>% dplyr::filter(Ord_Key_A==Ord_Key_B),
    aqp_bsp_inn %>% dplyr::filter(Ord_Key_A!=Ord_Key_B)
  )
  aqp_bsp_bad <- dplyr::anti_join(
    aqp_bsp_inn,
    aqp_bsp_pas, 
    by=c("Address_A","Address_B")
  )
  
  # Reload the original order file::
  ord_tib <- load_aqp_file(opt$ords)
  mat_tib <- load_aqp_file(opt$mats)
  aqp_tib <- load_aqp_file(opt$aqps)
  aqp_tib %>% dplyr::group_by(Decode_Status) %>% dplyr::summarise(Count=n(), .groups="drop")
  
  bsp_tsv <- file.path(par$topDir, "scratch/workhorse_main_dev_latest/NZT-NZT-C0-GRCh37/aln/NZT-NZT-C0-GRCh37.aqp-pass.address.bsp.tsv.gz")
  bsp_tib <- load_bsmap(bsp_tsv, verbose = opt$verbose)
  
  bed_cols <- 
    cols(
      Imp_Chr = col_character(),
      Imp_Pos = col_integer(),
      imp_End = col_integer(),
      Imp_Cgn = col_integer(),
      Imp_Bld = col_character(),
      Imp_Srd = col_character()
    )
  
  imp_bed <- file.path(par$topDir, "data/improbe/scratch/cgnDB/dbSNP_Core4/design-input/GRCh37.cgn.bed.gz")
  imp_tib <- readr::read_tsv(imp_bed, col_names=names(bed_cols), col_types=bed_cols)
}






# End of file
