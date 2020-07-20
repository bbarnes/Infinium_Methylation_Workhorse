
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Basic Controls Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("illuminaio") ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Controls I/O::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

addControlType = function(tib,
                          verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'addControlType'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  tib <- tib %>% dplyr::mutate(
    Control_Type=case_when(
      stringr::str_starts(Probe_ID, 'neg')  ~ 'NEGATIVE',
      stringr::str_starts(Probe_ID, 'Norm') ~ 'Normalization',
      
      stringr::str_starts(Probe_ID, 'BS_Conversion_II') ~ 'BsConversion_II',
      stringr::str_starts(Probe_ID, 'BS_Conversion_I_') ~ 'BsConversion_I',
      
      stringr::str_starts(Probe_ID, 'BS_conversion_NoBiasHIII_ASPE') ~ 'BsConversion_NoBiasHIII_ASPE',
      
      stringr::str_starts(Probe_ID, 'Non_Specific_II_') ~ 'NonSpecific_II',
      stringr::str_starts(Probe_ID, 'Non_Specific_I_')  ~ 'NonSpecific_I',
      
      stringr::str_starts(Probe_ID, 'A_Hairpin1')   ~ 'Hairpin_A_1',
      stringr::str_starts(Probe_ID, 'C_Hairpin1')   ~ 'Hairpin_C_1',
      stringr::str_starts(Probe_ID, 'G_Hairpin1')   ~ 'Hairpin_G_1',
      stringr::str_starts(Probe_ID, 'T_Hairpin1')   ~ 'Hairpin_T_1',
      
      stringr::str_starts(Probe_ID, 'A_Hairpin2')   ~ 'Hairpin_A_2',
      stringr::str_starts(Probe_ID, 'C_Hairpin2')   ~ 'Hairpin_C_2',
      stringr::str_starts(Probe_ID, 'G_Hairpin2')   ~ 'Hairpin_G_2',
      stringr::str_starts(Probe_ID, 'T_Hairpin2')   ~ 'Hairpin_T_2',
      
      stringr::str_starts(Probe_ID, 'GT_mismatch')  ~ 'Mismatch_GT',
      
      stringr::str_starts(Probe_ID, '2_HIGH_MM')      ~ 'HIGH_MM_2',
      stringr::str_starts(Probe_ID, '3_HIGH_MM')      ~ 'HIGH_MM_3',
      stringr::str_starts(Probe_ID, '18_MEDIUM_MM')   ~ 'MEDIUM_MM_18',
      stringr::str_starts(Probe_ID, '39_MEDIUM_MM')   ~ 'MEDIUM_MM_39',
      stringr::str_starts(Probe_ID, '60_LOW_MM')      ~ 'LOW_MM_68',
      stringr::str_starts(Probe_ID, '74_LOW_MM')      ~ 'LOW_MM_74',
      
      stringr::str_starts(Probe_ID, '74_YEAST_3MM')   ~ 'YEAST_3MM_74',
      stringr::str_starts(Probe_ID, '90_YEAST_3MM')   ~ 'YEAST_3MM_90',
      
      stringr::str_starts(Probe_ID, 'nonPolyA_ATG2')  ~ 'NonPoly_A_ATG2',
      stringr::str_starts(Probe_ID, 'nonPolyA_HK1')   ~ 'NonPoly_A_HK1',
      stringr::str_starts(Probe_ID, 'nonPolyA_PPIH')  ~ 'NonPoly_A_PPIH',
      
      stringr::str_starts(Probe_ID, 'nonPolyC_HK2')   ~ 'NonPoly_C_HK2',
      stringr::str_starts(Probe_ID, 'nonPolyC_PPID')  ~ 'NonPoly_C_PPID',
      stringr::str_starts(Probe_ID, 'nonPolyC_PPIE')  ~ 'NonPoly_C_PPIE',
      
      stringr::str_starts(Probe_ID, 'nonPolyG_HK2')   ~ 'NonPoly_G_HK2',
      stringr::str_starts(Probe_ID, 'nonPolyG_HK3')   ~ 'NonPoly_G_HK3',
      stringr::str_starts(Probe_ID, 'nonPolyG_PPIE')  ~ 'NonPoly_G_PPIE',
      stringr::str_starts(Probe_ID, 'nonPolyG_PPIH')  ~ 'NonPoly_G_PPIH',
      stringr::str_starts(Probe_ID, 'nonPolyG_PPIG')  ~ 'NonPoly_G_PPIG',
      
      stringr::str_starts(Probe_ID, 'nonPolyT_ALDOB') ~ 'NonPoly_T_ALDOB',
      stringr::str_starts(Probe_ID, 'nonPolyT_HK2')   ~ 'NonPoly_T_HK2',
      stringr::str_starts(Probe_ID, 'nonPolyT_PPIH')  ~ 'NonPoly_T_PPIH',
      
      stringr::str_starts(Probe_ID, 'CA8')  ~ 'CA8',
      stringr::str_starts(Probe_ID, 'CA7')  ~ 'CA7',
      stringr::str_starts(Probe_ID, 'CA6')  ~ 'CA6',
      stringr::str_starts(Probe_ID, 'CA5')  ~ 'CA5',
      stringr::str_starts(Probe_ID, 'CA1')  ~ 'CA1',
      
      stringr::str_starts(Probe_ID, 'CC8')  ~ 'CC8',
      stringr::str_starts(Probe_ID, 'CC7')  ~ 'CC7',
      stringr::str_starts(Probe_ID, 'CC5')  ~ 'CC5',
      stringr::str_starts(Probe_ID, 'CC4')  ~ 'CC4',
      stringr::str_starts(Probe_ID, 'CC3')  ~ 'CC3',
      stringr::str_starts(Probe_ID, 'CC1')  ~ 'CC1',
      
      stringr::str_starts(Probe_ID, 'CG12') ~ 'CG12',
      stringr::str_starts(Probe_ID, 'CG10') ~ 'CG10',
      stringr::str_starts(Probe_ID, 'CG8')  ~ 'CG8',
      stringr::str_starts(Probe_ID, 'CG2')  ~ 'CG2',
      stringr::str_starts(Probe_ID, 'CG1')  ~ 'CG1',
      
      stringr::str_starts(Probe_ID, 'CT12') ~ 'CT12',
      stringr::str_starts(Probe_ID, 'CT10') ~ 'CT10',
      stringr::str_starts(Probe_ID, 'CT9')  ~ 'CT9',
      stringr::str_starts(Probe_ID, 'CT5')  ~ 'CT5',
      stringr::str_starts(Probe_ID, 'CT2')  ~ 'CT2',
      
      stringr::str_starts(Probe_ID, 'TM183T')  ~ 'TM183T',
      stringr::str_starts(Probe_ID, 'TM182T')  ~ 'TM182T',
      stringr::str_starts(Probe_ID, 'TM169T')  ~ 'TM169T',
      stringr::str_starts(Probe_ID, 'TM167T')  ~ 'TM167T',
      stringr::str_starts(Probe_ID, 'TM150T')  ~ 'TM150T',
      stringr::str_starts(Probe_ID, 'TM148T')  ~ 'TM148T',
      
      stringr::str_starts(Probe_ID, 'rs') ~ 'SNP',
      
      TRUE ~ NA_character_)
  )
  
  if (verbose>=vt) {
    mat_cnt <- tib %>% dplyr::group_by(Control_Type) %>% dplyr::summarise(Probe_Count=n()) %>% base::nrow()
    unk_cnt <- tib %>% dplyr::filter(is.na(Control_Type)) %>% dplyr::select(Probe_ID) %>% base::nrow()
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Match Count={mat_cnt}, Unknown Type={unk_cnt}.{RET}"))
    tib %>% dplyr::group_by(Control_Type) %>% dplyr::summarise(Probe_Count=n()) %>% print()
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}"))

  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       BSP Map Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

combineBspBed = function(bsp, bed, src,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'combineBspBed'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))

  stime <- system.time({
    bsp_cols <- c('Probe_ID', 'BSP_PrbSeq', 'BSP_Qual', 'BSP_Tag', 'BSP_Chrom', 'BSP_Pos', 'BSP_SRD', 'BSP_GapCnt', 'BSP_RefSeq', 'BSP_TopMisCnt', 'BSP_AllMisStr')
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading bsp={bsp}...{RET}"))
    bsp_tib <- suppressMessages(suppressWarnings( readr::read_tsv(bsp, col_names=bsp_cols) )) %>% dplyr::select(-BSP_Qual)
    #  dplyr::mutate(Design_Method="BSC") %>% 
    bsp_cnt <- bsp_tib %>% base::nrow()
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} bsp_cnt={bsp_cnt}...{RET}{RET}"))
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading bed={bed}...{RET}"))
    bed_tib <- suppressMessages(suppressWarnings( readr::read_tsv(bed) ))
    
    tib <- bsp_tib %>% dplyr::filter(BSP_Tag=='UM') %>% dplyr::inner_join(bed_tib, by=c("Probe_ID"="Full_ID") ) %>% 
      dplyr::mutate(Design_Method=src,
                    Probe_New_Seq_A=paste0(PRB1_U_SS,'G'),
                    Probe_New_Seq_B=paste0(PRB1_M_SS,'A'),
                    Probe_New_Seq_D=PRB2_D_UC)
    bed_cnt <- bed_tib %>% base::nrow()
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} bed_cnt={bed_cnt}...{RET}{RET}"))

  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  tib
}

bsp2prbsControls = function(tib, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'addControlType'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (verbose>=vt) tib %>% dplyr::group_by(Control_Type) %>% dplyr::summarise(Count=n()) %>% as.data.frame() %>% print()
  
  prbs <- tib  %>%
    dplyr::mutate(
      BSP_RefBscU=case_when(
        # Control_Type=='BsConversion_II' & BSP_SRD=='--' ~ BSP_RefSeq,
        # Control_Type=='BsConversion_II' & BSP_SRD=='-+' ~ revCmp(BSP_RefSeq),
        # Control_Type=='BsConversion_II' & BSP_SRD=='+-' ~ revCmp(BSP_RefSeq),
        # Control_Type=='BsConversion_II' & BSP_SRD=='++' ~ BSP_RefSeq,
        
        stringr::str_starts(Control_Type, 'NonPoly_G') & BSP_SRD=='--' ~ BSP_RefSeq,          # NOT CHECKED
        stringr::str_starts(Control_Type, 'NonPoly_G') & BSP_SRD=='-+' ~ bscUs(revCmp(BSP_RefSeq)),
        stringr::str_starts(Control_Type, 'NonPoly_G') & BSP_SRD=='+-' ~ revCmp(BSP_RefSeq),  # NOT CHECKED
        stringr::str_starts(Control_Type, 'NonPoly_G') & BSP_SRD=='++' ~ bscUs(BSP_RefSeq),
        
        # TBD:: NORM_A/G Seems to be working now::
        stringr::str_starts(Control_Type, 'NORM_G') & BSP_SRD=='--' ~ BSP_RefSeq,          # Can Probably Remove
        stringr::str_starts(Control_Type, 'NORM_A') & BSP_SRD=='-+' ~ bscUs(revCmp(BSP_RefSeq)),
        stringr::str_starts(Control_Type, 'NORM_G') & BSP_SRD=='-+' ~ bscUs(revCmp(BSP_RefSeq)),
        stringr::str_starts(Control_Type, 'NORM_G') & BSP_SRD=='+-' ~ revCmp(BSP_RefSeq),  # Can Probably Remove
        stringr::str_starts(Control_Type, 'NORM_A') & BSP_SRD=='++' ~ bscUs(BSP_RefSeq),
        stringr::str_starts(Control_Type, 'NORM_G') & BSP_SRD=='++' ~ bscUs(BSP_RefSeq),
        
        BSP_SRD=='--' ~ revCmp(bscUs(revCmp(BSP_RefSeq))), # 728
        BSP_SRD=='-+' ~ revCmp(BSP_RefSeq),
        BSP_SRD=='+-' ~ revCmp(bscUs(BSP_RefSeq)), # 993
        BSP_SRD=='++' ~ BSP_RefSeq,  # Working So Far
        TRUE ~ NA_character_
      ),
      
      BSP_PrePosB= 2,
      BSP_NxbPosB= 3+stringr::str_length(BSP_PrbSeq),
      BSP_NxbPosE= 3+stringr::str_length(BSP_PrbSeq)+1,
      BSP_RefPreU=BSP_RefBscU %>% stringr::str_to_upper() %>% stringr::str_sub(BSP_PrePosB,BSP_PrePosB),
      BSP_RefNxbU=BSP_RefBscU %>% stringr::str_to_upper() %>% stringr::str_sub(BSP_NxbPosB,BSP_NxbPosB),
      BSP_RefPrbU=BSP_RefBscU %>% stringr::str_to_upper() %>% stringr::str_sub(3) %>% stringr::str_sub(1,stringr::str_length(BSP_PrbSeq)),
      BSP_CmpPrbU=BSP_RefPrbU %>% stringr::str_to_upper() %>% stringr::str_sub(1,stringr::str_length(BSP_RefPrbU)-1)
    )

  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}"))
  
  prbs
}

loadControlProbes = function(name, dir, opt,
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadControlProbes'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))

  # name <- '01152015_DarkMatterControls.unique.probe.match'
  ctl_tsv <- file.path(dir, name)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading ctl_tsv={ctl_tsv}...{RET}"))
  ctl_tib <- suppressMessages(suppressWarnings( readr::read_tsv(ctl_tsv) )) %>% dplyr::distinct()
  
  ctl_cnt <- ctl_tib %>% base::nrow()
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loaded Control Sequences(n={ctl_cnt}).{RET}"))
  
  seq_tib <- ctl_tib %>% dplyr::rename(Probe_ID=probe_id, Probe_Seq=sequence, Address=address_name) %>% 
    dplyr::mutate(Control_ID=stringr::str_remove(Probe_ID, '_[AB]$'),
                  AB_Str=stringr::str_replace(Probe_ID, '^.*_([AB])$', '\\$1') %>% stringr::str_remove('\\\\'),
                  Last_Base_Src=stringr::str_sub(Probe_Seq, -1),
                  Address=stringr::str_remove(Address, '^1') %>% stringr::str_remove('^0+') %>% as.numeric(),
                  Address_Src=Address,
                  Probe_Seq_S=Probe_Seq %>% stringr::str_to_upper() %>% stringr::str_sub(1,stringr::str_length(Probe_Seq)-1)
    ) %>% addControlType() %>% 
    # dplyr::select(Control_Type, Control_ID, Probe_ID, AB_Str, Last_Base_Src, Address, Address_Src, Probe_Seq, Probe_Seq_S) %>%
    dplyr::arrange(Control_Type, Control_ID)
  
  # Format from stack to paired probes
  #
  prb_tib <- dplyr::full_join(seq_tib %>% dplyr::filter(AB_Str=='A'),
                                  seq_tib %>% dplyr::filter(AB_Str=='B'),
                                  by=c("Control_ID", "Control_Type"), suffix=c('_A', '_B')) %>%
    dplyr::mutate(Design_Type=case_when(!is.na(AB_Str_B) ~ 'I',
                                        is.na(AB_Str_B) ~ 'II',
                                        TRUE ~ NA_character_),
                  Last_Base_A=stringr::str_sub(Probe_Seq_A, -1),
                  Last_Base_B=stringr::str_sub(Probe_Seq_B, -1) )
  
  # Sanity Checks::
  #  Check that all Infinium I probes are equal except the last base::
  #  Check the purposely swapped base distributions::
  #
  if (verbose>=vt) {
    prb_tib %>% dplyr::filter(Design_Type=='I') %>% dplyr::filter(Probe_Seq_S_A!=Probe_Seq_S_B) %>% print()
    prb_tib %>%
      dplyr::filter(Design_Type=='I') %>%
      dplyr::group_by(AB_Str_A, Last_Base_Src_A, AB_Str_B, Last_Base_Src_B) %>% dplyr::summarise(Count=n()) %>% print()
    prb_tib %>% dplyr::filter(Design_Type=='I') %>%
      dplyr::group_by(AB_Str_A, Last_Base_A, AB_Str_B, Last_Base_B) %>% dplyr::summarise(Count=n()) %>% print()
  }
  
  # Last Base Mapping::
  last_base_map <- prb_tib %>% dplyr::group_by(Control_Type, AB_Str_A, AB_Str_B, Last_Base_A, Last_Base_B) %>% 
    dplyr::summarise(Count=n()) %>% arrange(AB_Str_A, AB_Str_B)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Write Original Manifest Fasta::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  opt$write_original_fasta <- FALSE
  if (opt$write_original_fasta) {
    ord_fas <- file.path(opt$outDir, paste(name,'bscU.fa.gz', sep='.') )
    fas_tib <- prb_tib %>% 
      dplyr::select(Probe_ID_A, Address_A, Probe_Seq_A) %>% 
      dplyr::mutate(line=paste0('>',Address_A,'_',Probe_ID_A,'\n',Probe_Seq_A)) %>% dplyr::pull(line)
    if (opt$writeOutput) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing ord_fas={ord_fas}...{RET}"))
      readr::write_lines(x=fas_tib, path=ord_fas)
    }
    
    # Print all Data::
    ord_fas <- file.path(opt$outDir, paste(name,'fa.gz', sep='.') )
    fas_tib <- c(prb_tib %>% dplyr::select(Probe_ID_A, Address_A, Probe_Seq_A) %>% dplyr::filter(!is.na(Probe_ID_A)) %>%
                       dplyr::mutate(line=paste0('>',Address_A,'_',Probe_ID_A,'\n',Probe_Seq_A)) %>% dplyr::pull(line),
                     prb_tib %>% dplyr::select(Probe_ID_B, Address_B, Probe_Seq_B) %>% dplyr::filter(!is.na(Probe_ID_B)) %>%
                       dplyr::mutate(line=paste0('>',Address_B,'_',Probe_ID_B,'\n',Probe_Seq_B)) %>% dplyr::pull(line) )
    if (opt$writeOutput) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing ord_fas={ord_fas}...{RET}"))
      readr::write_lines(x=fas_tib, path=ord_fas)
    }
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}"))

  prb_tib
}

formatBSP = function(bsp) {
  bsp_bed_cols <- c('BSP_Chrom','BSP_Beg','BSP_End','Probe_ID','BSP_Tag','BSP_Srd','BSP_PrbSeq','BSP_RefSeq','BSP_GapCnt','BSP_MisCnt','BSP_MisStr',
                    'GEN_Chrom','GEN_Beg','GEN_End','GEN_Tran','GEN_Srd','GEN_Class')
  
  tib <- readr::read_tsv(bsp, col_names=bsp_bed_cols) %>% 
    tidyr::separate(BSP_MisStr, into=c('H0','H1','H2','H3','H4','H5'), sep=':', convert=TRUE) %>%
    dplyr::filter(H0>0) %>%
    tidyr::separate(Probe_ID, into=c('CGN', 'FR', 'TB', 'CO', 'ScrU', 'ScrM', 'ImpSumCnt', 'ImpCgnCnt'), sep='_') %>% 
    tidyr::unite(Probe_ID, CGN, FR, TB, CO, remove=FALSE, sep='_') %>%
    tidyr::separate(GEN_Tran, into=c('Gene', 'Transcript'), sep=':') %>%
    dplyr::mutate(MinScore=as.double(pmin(ScrU,ScrM)))
  
  tib
}


# End of file
