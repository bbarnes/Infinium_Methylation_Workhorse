
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           workhorse_function()::
#
#                           [FUNCTION GRAVEYARD]
#
#  The functions have been replaced. This file should be deleted, but is 
#           just a resting place for these functions for now :)
#
#              These functions should probably be deprecated
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Manifest Mutation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


if (!full) ret_tib <- ret_tib %>% 
    dplyr::select(- dplyr::all_of(rem_sel_vec))

# Sort by genomic position:: bsp_chr_sym
#  if (sort) ret_tib <- ret_tib %>% dplyr::arrange(Bsp_Chr, Bsp_Beg)
if (sort) ret_tib <- ret_tib %>% dplyr::arrange(!!bsp_chr_sym, !!bsp_beg_sym)

ann_key <- glue::glue("ret-annotated({funcTag})")
ann_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ann_key)

ret_key <- glue::glue("calculated-fields({funcTag})")
ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)

if (verbose>=vt)
  cat(glue::glue("{mssg} Done. Calculating new fields.{RET2}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Annotate CGN by Alignment Position::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

map_tib <- 
  load_cgn_map_tsv(file=cgn_src, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
ret_key <- glue::glue("map_tib({funcTag})")
ret_cnt <- print_tib(map_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)

if (!is.null(map_tib)) {
  if (verbose>=vt)
    cat(glue::glue("{mssg} Adding cgn by alignment position...{RET}"))
  
  # map_tib <- map_tib %>% 
  #   purrr::set_names(c(chr_key, pos_key, cgn_int_key, top_srd_key))
  # ret_key <- glue::glue("map_tib; updated-fields({funcTag})")
  # ret_cnt <- print_tib(map_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  
  sel_tib <-                 # Address,    # Chromosome   # Coordinate
    dplyr::select(ret_tib, !!ord_join_sym, !!bsp_chr_sym, !!bsp_pos_sym )
  sel_col <- sel_tib %>% names()
  sel_tib <- sel_tib %>% purrr::set_names("Id", "Chr", "Pos")
  
  # Chromosome, Coordinate, Cgn_Int, Top_Srd
  map_col <- map_tib %>% names()
  map_tib <- purrr::set_names("Chr","Pos","Cgn","Top") %>%
    dplyr::mutate(CG_Nuc="up")
  
  int_col <- c(sel_col[1],sel_col[2],sel_col[3],map_col[3],map_col[4],"Query_CG")
  ups_tib <- dplyr::inner_join(sel_tib,map_tib, by=c(Chr,Pos)) %>% 
    purrr::set_names(int_col)
  print(ups_tib)
  
  dns_tib <- dplyr::inner_join(
    sel_tib, map_tib %>% dplyr::mutate(Pos=Pos+1, CG_Nuc="dn"),
    by=c(Chr,Pos)) %>% purrr::set_names(int_col)
  print(dns_tib)
  
  # ups_tib <- dplyr::inner_join(
  #   dplyr::select(ret_tib, !!ord_join_sym, !!bsp_chr_sym, !!bsp_pos_sym ),
  #   dplyr::mutate(map_tib, !!cg_nuc_sym := "up"),
  #   by=c(chr_key, pos_key),
  #   suffix=c("_bsp","_map")
  # )
  
  dns_tib <- dplyr::inner_join(
    dplyr::select(ret_tib, !!ord_join_sym, !!bsp_chr_sym, !!bsp_pos_sym ),
    dplyr::mutate(map_tib, 
                  !!bsp_pos_sym := !!pos_key + 1,
                  # !!bsp_pos_sym:=!!bsp_pos_sym+1,
                  !!cg_nuc_sym := "dn"),
    by=c(chr_key, pos_key),
    suffix=c("_bsp","_map")
  )
  print(dns_tib)
  

  
  #
  #
  # IMPORTANT:: 
  #
  #
  ret_tib <- ret_tib %>% 
    # dplyr::left_join(mat_tib, by=c("Aln_Key_Unq","!!bsp_chr_sym","Bsp_Pos")) %>%
    dplyr::left_join(mat_tib, by=c(unq_key,"!!bsp_chr_sym","Bsp_Pos")) %>%
    dplyr::mutate(
      Bsp_TB=dplyr::case_when(
        Cgd_Top=="+" & Bsp_FR=="F" ~ "T",
        Cgd_Top=="-" & Bsp_FR=="R" ~ "T",
        
        Cgd_Top=="+" & Bsp_FR=="R" ~ "B",
        Cgd_Top=="-" & Bsp_FR=="F" ~ "B",
        TRUE ~ NA_character_
      )
    )
  
  if (verbose>=vt)
    cat(glue::glue("{mssg} Done. Adding cgn by position.{RET2}"))
}









# TBD:: Pretty Sure this can be removed::
mutate_probe_id = function(tib, 
                           
                           pid="Probe_ID", 
                           cgn="Imp_Cgn_Seq",
                           des="Ord_Des",  
                           din="Ord_Din",
                           tb="Imp_TB_Seq", 
                           co="Imp_CO_Seq",
                           inf="Infinium_Design",
                           
                           pad=8,
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutate_probe_id'
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   funcTag={funcTag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    pid_sym <- rlang::sym(pid)
    cgn_sym <- rlang::sym(cgn)
    des_sym <- rlang::sym(des)
    din_sym <- rlang::sym(din)
    tb_sym  <- rlang::sym(tb)
    co_sym  <- rlang::sym(co)
    inf_sym <- rlang::sym(inf)
    
    ret_tib <- tib %>% 
      dplyr::mutate(
        !!inf_sym:=dplyr::case_when(
          !!des_sym=="U" | !!des_sym=="M" ~ 1, 
          !!des_sym=="2" ~ 2, 
          TRUE ~ NA_real_) %>% as.integer(), 
        !!pid_sym:=paste(paste0(!!din_sym,stringr::str_pad(!!cgn_sym,width=pad, side="left", pad="0")), 
                         paste0(!!tb_sym,!!co_sym,!!inf_sym), sep="_")
      )
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET}{RET}"))
  
  ret_tib
}



# TBD:: Write load_genome_studio_address()
# TBD:: Test with annotation!!!

old_cols <- c("IlmnID", "Name", 
              "AddressA_ID", "AlleleA_ProbeSeq", "AddressB_ID", "AlleleB_ProbeSeq",
              "Color_Channel", "CHR", "MAPINFO", "Strand")
new_cols <- c("Probe_ID", "Prb_Cgn", 
              "Address_A", "Prb_Seq_A", "Address_B", "Prb_Seq_B",
              "Color", "Chromosome", "Coordinate", "Strand_FR")

# genome_manifest_vec  <- splitStrToVec(opt$genome_manifest_csv)
genome_manifest_list <- get_file_list(files=opt$genome_manifest_csv,
                                      trim = c(".csv.gz"), 
                                      alpha_numeric = TRUE, del = COM)

genome_manifest_dat <- lapply(genome_manifest_list, load_genome_studio_manifest,
                              load_clean    = TRUE,
                              load_controls = TRUE,
                              write_clean   = TRUE,
                              overwrite     = TRUE,
                              ret_data      = FALSE,
                              verbose = opt$verbose, tt = pTracker)


genome_manifest_vec  <- splitStrToVec(opt$genome_manifest_csv)

genome_address_dat <- 
  load_genome_studio_address(file = genome_manifest_vec[1],
                             load_clean     = TRUE,
                             load_controls  = TRUE,
                             write_clean    = TRUE,
                             overwrite      = TRUE, 
                             add_annotation = TRUE,
                             ret_data       = FALSE,
                             verbose = opt$verbose+10, tt = pTracker)

join_cols <- c("Probe_ID","Prb_Des","Prb_Din")
int_cols <- 
  intersect( names(genome_manifest_dat[[1]]), names(genome_manifest_dat[[2]]) )

int_tibs <- dplyr::inner_join(
  genome_manifest_dat[[1]] %>% dplyr::select(dplyr::all_of(int_cols)),
  genome_manifest_dat[[2]] %>% dplyr::select(dplyr::all_of(int_cols)),
  by=join_cols,
  suffix=c("_a","_b")
)

cmp_tib <- NULL
for (col_key in int_cols) {
  if (col_key %in% join_cols) next
  cat(glue::glue("col_key={col_key}.{RET}"))
  
  colA_key <- paste0(col_key,"_a")
  colB_key <- paste0(col_key,"_b")
  
  colA_sym <- rlang::sym(colA_key)
  colB_sym <- rlang::sym(colB_key)
  
  cur_tib <- int_tibs %>% 
    dplyr::select(dplyr::all_of( c(colA_key,colB_key) ) ) %>%
    dplyr::summarise(
      Mat_Cnt=sum(!!colA_sym == !!colB_sym, na.rm = TRUE),
      Mis_Cnt=sum(!!colA_sym != !!colB_sym, na.rm = TRUE),
      # Mat_Cnt=sum(!!colA_key == !!colB_key),
      # Mis_Cnt=sum(!!colA_key != !!colB_key), 
      Total=sum(!is.na(!!colA_sym) & !is.na(!!colB_sym) ),
      Mat_Per=round(100*Mat_Cnt/Total, 2),
      Mis_Per=round(100*Mis_Cnt/Total, 2)
    ) %>% dplyr::mutate(Column=col_key)
  
  cmp_tib %>% print()
  
  cmp_tib <- cmp_tib %>% dplyr::bind_rows(cur_tib)
}

#
# Get Sesame Data::
#
sesame_data_keys  <- sesameData::sesameDataList()

sesameDataCache("EPIC")
sesameDataCache("HM450")
sesameDataCacheAll()

sesame_manifest_vec  <- splitStrToVec(opt$sesame_manifest_dat)
sesmae_manifest_list <- lapply(sesame_manifest_vec, sesameData::sesameDataGet)

hm450_probeInfo_key <- "HM450.probeInfo"
hm450_hg19_man_key  <- "HM450.hg19.manifest"
hm450_probeInfo_dat <- sesameData::sesameDataGet( hm450_probeInfo_key )
hm450_hg19_man_grs  <- sesameData::sesameDataGet( hm450_hg19_man_key )
hm450_hg19_man_tib  <- hm450_hg19_man_grs %>% as.data.frame() %>% 
  tibble::rownames_to_column(var="Probe_ID") %>% tibble::as_tibble()

epic_probeInfo_key <- "EPIC.probeInfo"
epic_hg19_man_key  <- "EPIC.hg19.manifest"
epic_probeInfo_dat <- sesameData::sesameDataGet( epic_probeInfo_key )
epic_hg19_man_grs  <- sesameData::sesameDataGet( epic_hg19_man_key )
epic_hg19_man_tib  <- epic_hg19_man_grs %>% as.data.frame() %>% 
  tibble::rownames_to_column(var="Probe_ID") %>% tibble::as_tibble()

gen_info_hg19_key   <- "genomeInfo.hg19"
gen_info_hg19_dat   <- sesameData::sesameDataGet( gen_info_hg19_key )

sesame_address_list <- get_file_list(files=opt$sesame_manifest_dat,
                                     trim = c(".csv.gz"), del = COM)


hm450_ses_add_tib <- load_sesame_repo_address(hm450_hg19_man_key, 
                                              verbose=opt$verbose+10, tt=pTracker)

epic_ses_add_tib <- load_sesame_repo_address(epic_hg19_man_key, 
                                             verbose=opt$verbose, tt=pTracker)

epic_hg19_gst_man_csv <- file.path(par$topDir, "data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B4-Beadpool_ID.csv.gz")
# epic_hg19_gst_out_csv <- file.path(par$topDir, "data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B4-Beadpool_ID.analytical-only.csv.gz")

epic_hg19_gst_man_lines <- readr::read_lines(epic_hg19_gst_man_csv)

beg_idx <- which( stringr::str_starts(head(epic_hg19_gst_man_lines, n=20), "IlmnID" ) ) %>% head(n=1)
con_idx <- which( stringr::str_starts(epic_hg19_gst_man_lines, "\\[Controls\\]" ) )
end_idx <- epic_hg19_gst_man_lines %>% length()


if (FALSE) {
  
  epic_ses_add_tib %>% 
    dplyr::arrange(Address) %>% 
    dplyr::add_count(Address, name="Add_Cnt") %>% 
    dplyr::filter(Add_Cnt!=1) %>% 
    head() %>% as.data.frame()
  
  epic_ses_add_tib %>% 
    dplyr::group_by(Prb_Des) %>% 
    dplyr::summarise(Count=n(), .groups = "drop")
  
  #
  # Consolidate 450k and EPIC Sesame data::
  #
  # epic_hg19_man_tib %>% dplyr::select(Probe_ID, seqnames, start, end, strand, address_A, ProbeSeq_A, address_B, ProbeSeq_B, )
  
  old_names <- c("Probe_ID", "seqnames", "start", "end", "strand", "designType", "channel",
                 "nextBase", "nextBaseRef", "probeType", "gene", "gene_HGNC")
  new_names <- c("Probe_ID", "Chromosome", "Coordinate", "CoordinateG", "Strand_FR", "Infinium_Design_Type", "Color",
                 "Prb_Nxb", "Prb_Nxb_Ref", "Prb_Din", "gene", "gene_HGNC")
  
  prb_key <- "A"
  oldA_names <- 
    epic_hg19_man_tib %>% 
    dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>% 
    dplyr::select(-dplyr::starts_with("wDecoy_")) %>% 
    names()
  newA_names <- oldA_names %>% 
    stringr::str_remove(paste0("_",prb_key,"$")) %>% 
    stringr::str_to_title() %>%
    paste("Prb",., sep="_") %>%
    stringr::str_replace("Prb_Address", "Address") %>%
    stringr::str_replace("Prb_Probeseq", "Prb_Seq") %>%
    stringr::str_replace("Prb_Chrm", "Prb_Chr")
  
  epic_hg19_man_tib %>% 
    dplyr::select(dplyr::starts_with("wDecoy_")) %>%
    dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>%
    names() %>%
    stringr::str_remove(paste0("_",prb_key,"$")) %>% 
    stringr::str_remove("^wDecoy_") %>% 
    stringr::str_to_title() %>%
    stringr::str_replace("Chrm", "Chr") %>%
    paste("Prb_Decoy",., sep="_")
  
  ses_epic_man_tib <- epic_hg19_man_tib %>%
    dplyr::select(Probe_ID, seqnames, start, end, strand, designType, channel,
                  nextBase, nextBaseRef, probeType, gene, gene_HGNC,
                  dplyr::ends_with("_A"), dplyr::starts_with("MASK_") ) %>%
    dplyr::select(-dplyr::starts_with("wDecoy_")) %>%
    dplyr::rename(Chromosome=seqnames, Coordinate=start, CoordinateG=end,
                  Strand_FR=strand, Infinium_Design_Type=designType,
                  Color=channel, Prb_Nxb=nextBase, Prb_Nxb_Ref=nextBaseRef, 
                  Prb_Din=probeType, 
                  
                  Address=address_A, Prb_Seq=ProbeSeq_A, 
                  Prb_Chr=chrm_A, Prb_Beg=beg_A, Prb_Flag=flag_A, 
                  Prb_MapQ=mapQ_A, Prb_Cigar=cigar_A, Prb_NM=NM_A) %>%
    
    dplyr::mutate(Strand_FR=dplyr::case_when(
      Strand_FR=="+" ~ "F",
      Strand_FR=="-" ~ "R",
      TRUE ~ NA_character_),
      Prb_Nxb=stringr::str_remove_all(Prb_Nxb,"[^a-zA-Z]") %>% mapDIs(),
      Prb_Inf=dplyr::case_when(
        Infinium_Design_Type=="I"  ~ 1,
        Infinium_Design_Type=="II" ~ 2,
        TRUE ~ NA_real_
      ) %>% as.integer()
    )
  
  ses_epic_man_tib %>% 
    dplyr::group_by(Prb_Inf, Prb_Nxb, Prb_Nxb_Ref, Color) %>% 
    dplyr::summarise(Count=n(), .groups = "drop")
  
  #
  # Load Genome Studio Manifests::
  #
  # epic_hg19_gst_man_csv <- file.path(par$topDir, "data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B4-Beadpool_ID.csv.gz")
  # epic_hg19_test_man_csv <- "/Users/bretbarnes/Documents/tmp/MethylationEPIC_v-1-0_B4-Beadpool_ID.csv.gz"
  
  
}

if (FALSE) {
  
  col_prefix <- "HM"
  probe_name_cols <- 
    epic_hg19_gst_man_lines[beg_idx] %>% 
    stringr::str_split(pattern = ",", simplify = TRUE) %>% 
    as.vector() %>% 
    stringr::str_replace("^([0-9])", paste0(col_prefix, "\\$1") ) %>%
    stringr::str_remove_all("\\\\")
  
  probes_tib <- 
    epic_hg19_gst_man_lines[c( (beg_idx+1):(con_idx-1) ) ] %>% head() %>%
    tibble::as_tibble() %>% 
    tidyr::separate(value, into=c(probe_name_cols), sep=',') %>% 
    clean_tibble()
  
  col_types_str <- sapply(probes_tib, typeof) %>% 
    tibble::as_tibble(rownames = "Key") %>% 
    dplyr::mutate(value=stringr::str_sub(value, 1,1)) %>% 
    dplyr::pull(value) %>% paste0(collapse = '')
  
  probes_type_cols <- spec(readr::read_csv(
    readr::format_csv(probes_tib), col_types = col_types_str))
  
  gs_control_names <- c("Address","Control_Group","GS_Color","Control_Type")
  controls_tib <- 
    epic_hg19_gst_man_lines[c( (con_idx+1):(end_idx) ) ] %>%
    stringr::str_remove(",+$") %>%
    stringr::str_remove(",AVG$") %>%
    tibble::as_tibble() %>% 
    tidyr::separate(value, into=c(gs_control_names), sep=',') %>% 
    clean_tibble()
}

if (FALSE) {
  epic_hg19_gst_man_tib <- readr::read_csv(epic_hg19_gst_out_csv)
  epic_hg19_gst_cln_tib <- epic_hg19_gst_man_tib %>% clean_tibble()
  
  key1_tib <- sapply(epic_hg19_gst_man_tib, typeof) %>% tibble::as_tibble(rownames = "Key")
  key2_tib <- sapply(epic_hg19_gst_cln_tib, typeof) %>% tibble::as_tibble(rownames = "Key")
  epic_hg19_gst_spec <- spec(epic_hg19_gst_man_tib)
  
  key2_tib %>% dplyr::inner_join(key1_tib, by="Key", suffix=c("_new", "_old")) %>% dplyr::filter(value_new != value_old)
  
  dot_prefix <- "HM"
  epic_hg19_gst_man_tib %>% tibble::as_tibble(.name_repair = "universal") %>% purrr::set_names(stringr::str_replace(names(.), "^..", dot_prefix)) %>% names()
  
  epic_hg19_gst_spec$cols$AddressA_ID <- col_integer()
  epic_hg19_gst_spec$cols$AddressB_ID <- col_integer()
  epic_hg19_gst_spec$cols$Genome_Build <- col_integer()
  epic_hg19_gst_spec$cols$DNase_Hypersensitivity_Evidence_Count <- col_integer()
  epic_hg19_gst_spec$cols$TFBS_Evidence_Count <- col_integer()
  epic_hg19_gst_spec$cols$Coordinate_36 <- col_integer()
  
  cols(
    IlmnID = col_character(),
    Name = col_character(),
    AddressA_ID = col_integer(),
    AlleleA_ProbeSeq = col_character(),
    AddressB_ID = col_integer(),
    AlleleB_ProbeSeq = col_character(),
    Infinium_Design_Type = col_character(),
    Next_Base = col_character(),
    Color_Channel = col_character(),
    Forward_Sequence = col_character(),
    Genome_Build = col_integer(),
    CHR = col_character(),
    MAPINFO = col_integer(),
    SourceSeq = col_character(),
    Strand = col_character(),
    UCSC_RefGene_Name = col_character(),
    UCSC_RefGene_Accession = col_character(),
    UCSC_RefGene_Group = col_character(),
    UCSC_CpG_Islands_Name = col_character(),
    Relation_to_UCSC_CpG_Island = col_character(),
    Phantom4_Enhancers = col_character(),
    Phantom5_Enhancers = col_character(),
    DMR = col_character(),
    HM450_Enhancer = col_logical(),
    HMM_Island = col_character(),
    Regulatory_Feature_Name = col_character(),
    Regulatory_Feature_Group = col_character(),
    GencodeBasicV12_NAME = col_character(),
    GencodeBasicV12_Accession = col_character(),
    GencodeBasicV12_Group = col_character(),
    GencodeCompV12_NAME = col_character(),
    GencodeCompV12_Accession = col_character(),
    GencodeCompV12_Group = col_character(),
    DNase_Hypersensitivity_NAME = col_character(),
    DNase_Hypersensitivity_Evidence_Count = col_integer(),
    OpenChromatin_NAME = col_character(),
    OpenChromatin_Evidence_Count = col_integer(),
    TFBS_NAME = col_character(),
    TFBS_Evidence_Count = col_integer(),
    Methyl27_Loci = col_logical(),
    Methyl450_Loci = col_logical(),
    Chromosome_36 = col_character(),
    Coordinate_36 = col_integer(),
    SNP_ID = col_character(),
    SNP_DISTANCE = col_character(),
    SNP_MinorAlleleFrequency = col_character(),
    Random_Loci = col_logical(),
    Beadpool_ID = col_character()
  )
}

# End of file
