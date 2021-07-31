
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


epic_hg19_gst_man_csv <- file.path(par$topDir, "data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B4-Beadpool_ID.csv.gz")
# epic_hg19_gst_out_csv <- file.path(par$topDir, "data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B4-Beadpool_ID.analytical-only.csv.gz")

epic_hg19_gst_man_lines <- readr::read_lines(epic_hg19_gst_man_csv)

beg_idx <- which( stringr::str_starts(head(epic_hg19_gst_man_lines, n=20), "IlmnID" ) ) %>% head(n=1)
con_idx <- which( stringr::str_starts(epic_hg19_gst_man_lines, "\\[Controls\\]" ) )
end_idx <- epic_hg19_gst_man_lines %>% length()


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
