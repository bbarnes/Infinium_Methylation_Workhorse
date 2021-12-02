
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   improbe (Infinium Methylation) Methods::
#
#                           [FUNCTION GRAVEYARD]
#
#  The functions have been replaced. This file should be deleted, but is 
#           just a resting place for these functions for now :)
#
#              These functions should probably be deprecated
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

suppressWarnings(suppressPackageStartupMessages(require("Biostrings")) )

COM <- ","
TAB <- "\t"
RET <- "\n"
BNG <- "|"
BRK <- paste0("# ",
              paste(rep("-----",6),collapse=" "),"|",
              paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='template_func') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Function Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   funcTag={funcTag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    String Diff Method from the Internet::
#
#  NOTE:: Maybe useful, so we'll just keep it here for now. 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

list.string.diff<-function(a="ATTCGA-",b="attTGTT",exclude=c("-","?"),ignore.case=TRUE)
{
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  if(ignore.case==TRUE)
  {
    a<-toupper(a)
    b<-toupper(b)
  }
  seq.a<-unlist(strsplit(a,split=""))
  seq.b<-unlist(strsplit(b,split=""))
  diff.d<-rbind(seq.a,seq.b)
  only.diff<-diff.d[,diff.d[1,]!=diff.d[2,]]
  pos<-which(diff.d[1,]!=diff.d[2,])
  only.diff<-rbind(pos,only.diff)
  for(ex.loop in 1:length(exclude))
  {
    only.diff<-only.diff[,!(only.diff["seq.a",]==exclude[ex.loop]|only.diff["seq.b",]==exclude[ex.loop])]
  }
  return(only.diff)
}



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#              These functions should probably be deprecated
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Run Time Version Options:: 
  #                       Platform, Genome Build, etc
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

  opt$run_name     <- NULL
  opt$platform     <- NULL
  opt$version      <- NULL
  opt$genome_build <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Run Time User Input Directories:: 
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opt$out_dir <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Run Time User Input Files:: 
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opt$ords <- NULL
  opt$mats <- NULL
  opt$aqps <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Run Time User Input Executable(s):: 
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opt$Rscript   <- NULL
  opt$bsmap_opt <- NULL
  opt$bsmap_exe <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Pre-defined Static Data Directories:: 
  #            improbe, Annotation, Genomic, Manifest, Validation Idats
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opt$imp_dir  <- NULL
  opt$ann_dir  <- NULL
  opt$gen_dir  <- NULL
  opt$man_dir  <- NULL
  opt$idat_dir <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                  Pre-defined Static External File Options:: 
  #                   Manifest, Controls, Design Coordinates
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opt$load_sesame_manfiest  <- FALSE
  opt$sesame_manifest_csv   <- NULL
  opt$genome_manifest_csv   <- NULL
  opt$sesame_controls_csv   <- NULL
  opt$genome_controlst_csv  <- NULL
  opt$noob_controls_csv     <- NULL
  opt$source_coordinate_csv <- NULL
  opt$canonical_cgn_csv     <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Run Time File Options:: Time Stamps
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  opt$time_org_txt <- NULL

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                          Run Time Mode Options::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opt$single    <- FALSE
  opt$parallel  <- FALSE
  opt$cluster   <- FALSE
  opt$trackTime <- FALSE
  opt$fresh     <- FALSE
  opt$reload    <- FALSE
  
  opt$verbose   <- 3
  

  #
  #
  # Old Options Definitions::
  #
  #
  
  
  # Executable::
  opt$Rscript <- NULL
  
  # BSMAP Parameters::
  opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
  opt$bsmap_opt <- "\"-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R\""
  opt$bsmap_exe <- "/Users/bretbarnes/Documents/tools/programs/BSMAPz/bsmapz"
  
  # Run Parameters::
  opt$run_name    <- NULL
  
  # Null Place Holders::
  # opt$cpg_top_tsv <- NULL
  # opt$cpg_pos_tsv <- NULL
  # opt$cph_pos_tsv <- NULL
  # opt$snp_pos_tsv <- NULL
  # opt$ord_des_csv <- NULL
  
  # Directories::
  opt$out_dir  <- NULL
  opt$imp_dir  <- NULL
  opt$ann_dir  <- NULL
  opt$gen_dir  <- NULL
  opt$man_dir  <- NULL
  
  # Manufacturing Files:: Required
  opt$ords <- NULL
  opt$mats <- NULL
  opt$aqps <- NULL
  
  # Boolean Flag to load pre-defined Sesame manifest from Sesame source
  opt$load_sesame_manfiest <- FALSE
  
  # Pre-defined manifest(s) to be re-built and/or added to new manifest::
  opt$sesame_manifest_csv <- NULL
  opt$genome_manifest_csv <- NULL
  
  # Pre-defined manifest control(s) to be added to new manifest::
  opt$sesame_controls_csv <- NULL
  opt$genome_controls_csv <- NULL
  
  # Pre-defined noob-masked control(s) to be added to new manifest::
  opt$noob_controls_csv <- NULL
  
  # Validation existing idats directory to confirm Addresses against::
  opt$idat_dir <- NULL
  
  # Original source design file used for canonical position selection.::
  opt$source_coordinate_csv <- NULL
  
  # Validation existing idats directory to confirm Addresses against::
  opt$canonical_cgn_csv <- NULL
  
  # Platform/Method Options::
  opt$genome_build <- NULL
  opt$platform <- NULL
  opt$version  <- NULL
  
  # Run-time Files
  opt$time_org_txt <- NULL
  
  # Process Parallel/Cluster Parameters::
  opt$single   <- FALSE
  opt$parallel <- TRUE
  opt$cluster  <- FALSE
  
  # Run-time Options::
  opt$trackTime    <- FALSE
  opt$fresh  <- FALSE
  opt$reload <- FALSE
  
  # verbose Options::
  opt$verbose <- 3
  
  
  
  
  
  
  
  if (verbose>=vt) {
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    cat(glue::glue("{mssg} names(inf1_tib)={RET}"))
    inf1_names <- names(inf1_tib)
    print(inf1_names)
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} inf1_names ! in selU_cols={RET}"))
    inf1_names[!inf1_names %in% selU_cols] %>% print()
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} selU_cols ! in inf1_names={RET}"))
    selU_cols[!selU_cols %in% inf1_names] %>% print()
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} selM_cols ! in inf1_names={RET}"))
    selU_cols[!selM_cols %in% inf1_names] %>% print()
    cat(glue::glue("{mssg}{RET}"))
    
    selU_len <- selU_cols %>% length()
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    cat(glue::glue("{mssg} selU_cols({selU_len})={RET}"))
    print(selU_cols)
    # dplyr::select(inf1_tib, dplyr::all_of(selU_cols) ) %>% 
    #   purrr::set_names("Ord_Map",ids_key,"Cgn","Ord_Des",
    #                    "Ord_Din","Can_Cnt","Rank") %>% head(n=2) %>% print()
    
    selM_len <- selM_cols %>% length()
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    cat(glue::glue("{mssg} selM_cols({selM_len})={RET}"))
    print(selM_cols)
    # dplyr::select(inf1_tib, dplyr::all_of(selM_cols) ) %>% 
    #   purrr::set_names("Ord_Map",ids_key,"Cgn","Ord_Des",
    #                    "Ord_Din","Can_Cnt","Rank") %>% head(n=2) %>% print()
    
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    cat(glue::glue("{mssg} names(inf1_tib)={RET}"))
    inf2_names <- names(inf2_tib)
    print(inf2_names)
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg} inf2_names ! in sel2_cols={RET}"))
    inf2_names[!inf2_names %in% sel2_cols] %>% print()
    cat(glue::glue("{mssg} sel2_cols ! in inf2_names={RET}"))
    sel2_cols[!sel2_cols %in% inf2_names] %>% print()
    cat(glue::glue("{mssg}{RET}"))
    
    sel2_len <- sel2_cols %>% length()
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    cat(glue::glue("{mssg} sel2_cols({sel2_len})={RET}"))
    print(sel2_cols)
    # dplyr::select(inf2_tib, dplyr::all_of(sel2_cols) ) %>% head(n=2) %>% print()
    
    name_len <- name_cols %>% length()
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    cat(glue::glue("{mssg} name_cols({name_len})={RET}"))
    print(name_cols)
    cat(glue::glue("{mssg}{RET}"))
    
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
  }
  
  # assign cgn scratch::
  #
  if (FALSE) {
    
    #
    #
    #
    #
    #   LEFT OFF HERE !!!!!!
    #
    #
    #
    
    
    
    can_tib_bk <- can_tib
    
    run$merge <- FALSE
    cgn_out_csv <- file.path(opt$outDir, "assigned-cgns.csv.gz")
    ret_dat <- assign_cgn_old(ord_tib = aqp_ord_tib,
                              bsp_tib = aqp_bsp_tib,
                              
                              seq_tib = aqp_seq_tib, 
                              can_csv = run$canonical,
                              can_tib = can_tib_bk,
                              
                              ids_key = run$ids_key,
                              
                              bsp_csv = cgn_out_csv,
                              merge = run$merge,
                              
                              retData = TRUE,
                              
                              verbose=opt$verbose+100)
    
  }
  

  
  
  ord_des_tib <- aqp_bsp_tib %>% 
    dplyr::select(run$unq_key, run$ids_key, run$add_key, 
                  run$des_key, run$din_key, 
                  run$srd_key, run$cos_key,
                  run$pos_key, run$chr_key, 
                  run$prb_key, "Bsp_Str", "Bsp_Tag", "Bsp_Prb_Dir")
  
  
  dplyr::inner_join(ord_des_tib, r_improbe_fwd_tib, by=c("Prb_Key_Unq", "Ord_Din")) %>% 
    dplyr::filter((Ord_Prb==Prb_1U & Ord_Des=="U") | (Ord_Prb==Prb_1M & Ord_Des=="M") | (Ord_Prb==Prb_2D & Ord_Des=="2")) %>% 
    dplyr::distinct(Prb_Key_Unq, .keep_all = TRUE) %>% dplyr::group_by(Strand_FR, Strand_CO) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  dplyr::inner_join(ord_des_tib, r_improbe_fwd_tib, by=c("Prb_Key_Unq"="Seq_ID", "Ord_Din")) %>%
    dplyr::filter((Ord_Prb==PRB1_U_MAT & Ord_Des=="U") | (Ord_Prb==PRB1_M_MAT & Ord_Des=="M") | (Ord_Prb==PRB2_D_MAT & Ord_Des=="2")) %>% 
    dplyr::distinct(Prb_Key_Unq, .keep_all = TRUE) %>% 
    dplyr::group_by(Strand_SR, Strand_CO) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  
  
  
  
  
  # dplyr::inner_join(ord_des_tib, r_improbe_fwd_tib, by=c("Prb_Key_Unq"="Seq_ID", "Ord_Din")) %>% 
  dplyr::inner_join(ord_des_tib, r_improbe_fwd_tib, by=c("Prb_Key_Unq"="Seq_ID", "Ord_Din")) %>%
    dplyr::mutate(
      Ord_Match=dplyr::case_when(
        Ord_Des=="2" & Ord_Prb == PRB2_D_MAT ~ "D0",
        Ord_Des=="2" & Ord_Prb != PRB2_D_MAT ~ "D1",
        
        Ord_Des=="U" & Ord_Prb == PRB1_U_MAT ~ "U0",
        Ord_Des=="U" & Ord_Prb != PRB1_U_MAT ~ "U1",
        
        Ord_Des=="M" & Ord_Prb == PRB1_M_MAT ~ "M0",
        Ord_Des=="M" & Ord_Prb != PRB1_M_MAT ~ "M1",
        
        TRUE ~ NA_character_
      )
    ) %>% dplyr::group_by(Ord_Match) %>% dplyr::summarise(Count=n(), .groups = "drop")
  
  
  
  
  dplyr::inner_join(
    dplyr::select(c_imp_tib, Seq_ID,Strand_Ref_FR,Strand_TB,Strand_CO, Probe_Seq_U,Probe_Seq_M, Top_Sequence),
    dplyr::select(s_imp_tib_FCM_tb, Raw_TB,Bsp_FR, Top_Key, starts_with("Prb")),
    by=c("Seq_ID"="Prb_Key_Unq")
  ) %>% dplyr::filter(Probe_Seq_M == Prb1C) %>%
    dplyr::select(Top_Sequence, Top_Key)  
  
  
  s_imp_tib_FCN_tb <- readr::read_csv(file_tab %>% dplyr::filter(Alphabet=="dna" & improbe_type=="s" & Strand_CO=="C", Strand_FR=="F" & Strand_BSC=="N") %>% dplyr::pull(Path)) %>%
    setTopBot_tib(seqKey="Forward_Sequence", srdKey="Raw_TB", topKey = "Top_Key") 
  
  dplyr::inner_join(
    dplyr::select(c_imp_tib, Seq_ID,Strand_Ref_FR,Strand_TB,Strand_CO, Probe_Seq_U,Probe_Seq_M, Top_Sequence),
    dplyr::select(s_imp_tib_FCN_tb, Raw_TB,Bsp_FR, Top_Key, starts_with("Prb")),
    by=c("Seq_ID"="Prb_Key_Unq")
  ) %>% dplyr::filter(Top_Sequence==Top_Key)  
  
  
  
  s_all_tib <- s_imp_tib %>% dplyr::select(starts_with("Prb")) %>% 
    dplyr::mutate(Prb1Ccm=cmpl(Prb1C),
                  Prb2Ccm=cmpl(Prb1C),
                  Prb1Ocm=cmpl(Prb1C),
                  PrbwOcm=cmpl(Prb1C),
                  
                  Prb1Crc=revCmp(Prb1C),
                  Prb2Crc=revCmp(Prb1C),
                  Prb1Orc=revCmp(Prb1C),
                  PrbwOrc=revCmp(Prb1C),
                  
                  Prb1Crv=reverse(Prb1C),
                  Prb2Crv=reverse(Prb1C),
                  Prb1Orv=reverse(Prb1C),
                  PrbwOrv=reverse(Prb1C) )
  
  cs_inn <- c_imp_tib %>% dplyr::inner_join(s_all_tib, by=c("Seq_ID"="Prb_Key_Unq"))
  
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1C)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1C)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1C)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1C)
  
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Ccm)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Ccm)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Ccm)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Ccm)
  
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Crc)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Crc)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Crc)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Crc)
  
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Crv)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Crv)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Crv)
  cs_inn %>% dplyr::filter(Probe_Seq_M==Prb1Crv)
  
  
  
  c_imp_tib %>% dplyr::inner_join(
    r_imp_tib, by=c("Seq_ID"="Prb_Key_Unq", 
                    "Forward_Sequence", "Probe_Seq_M"="Prb_1M",
                    "Strand_CO") ) %>%
    dplyr::select(Strand_FR.x, Strand_FR.y)
  
  
  c_imp_tib %>% dplyr::inner_join(
    r_imp_tib, by=c("Seq_ID"="Prb_Key_Unq", "Probe_Seq_M"="Prb_1M"), suffix=c("_c","_r"))
  
  s_imp_tib %>% dplyr::inner_join(
    r_imp_tib, by=c("Prb_Key_Unq")) %>% dplyr::filter("Prb1C"=="Prb_1M")
  
  s_imp_tib %>% dplyr::inner_join(
    r_imp_tib, by=c("Prb_Key_Unq")) %>% dplyr::filter("Prb1O"=="Prb_1M")
  
  s_imp_tib %>% 
    dplyr::inner_join(
      r_imp_tib, by=c("Prb_Key_Unq")) %>% dplyr::select(Prb_Key_Unq, Prb1C, Prb1O, Prb_1M, Prb_2D,
                                                        Forward_Sequence, 
                                                        Iupac_Forward_Sequence.x,
                                                        Iupac_Forward_Sequence.y,
                                                        Tmp_Pad)
  
  
  
  
  # NOTE:: Need to run r_improbe on Forward_Sequence
  
  
  all_dat <- NULL
  for (imp in names(file_tab_list)) {
    imp_tib <- NULL
    for (ii in c(1:base::nrow(file_tab_list[[imp]]))) {
      cat("ii=",ii,"\n")
      cur_tib <- safe_read(file_tab[ii, ]$Path) %>% dplyr::mutate(
        Key=file_tab[ii,]$Key, improbe_type=file_tab[ii,]$improbe_type, 
        Alpha=file_tab[ii,]$Alphabet, Bsc_Tag=file_tab[ii,]$Bsc_Tag, 
        FR=file_tab[ii,]$Strand_FR, CO=file_tab[ii,]$Strand_CO, 
        BSC=file_tab[ii,]$Strand_BSC)
      
      imp_tib <- dplyr::bind_rows(imp_tib, cur_tib)
      
      print(cur_tib)
    }
    all_dat[[imp]] <- imp_tib
  }
  
  dplyr::inner_join(all_dat$c,all_dat$r, 
                    by=c("Probe_Seq_M"="Prb_1M"), 
                    suffix=c("_c","_r") ) %>% dplyr::select(starts_with("Strand"))
  dplyr::group_by(Strand_Ref_FR_c, Strand_TB_c, Strand_CO_c, Strand_FR_c,  Strand_CO_r, Strand_FR_r,
                  Ord_Din)
  
  
  
  
  
  if (FALSE) {
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                         OLD COMPLEX VERSION::
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    
    aqp_prb_tib <- prb_designs_workflow_complex(
      tib = aqp_bsp_tib %>% 
        dplyr::select(run$unq_key, run$ids_key, run$add_key, 
                      run$des_key, run$din_key, 
                      run$srd_key, run$cos_key,
                      run$pos_key, run$chr_key), 
      max = 0,
      
      out_dir  = run$seq_dir,
      run_name = opt$runName,
      
      imp_level  = 3,
      
      # imp_inp_tsv = run$imp_inp_tsv,
      # imp_des_tsv = run$imp_des_tsv,
      # imp_fin_tsv = run$imp_fin_tsv,
      # imp_seq_csv = run$imp_seq_csv,
      
      # Genomes Parameters::
      gen_bld  = opt$genBuild,
      gen_nrec = run$gen_nrec,
      gen_key  = run$gen_key,
      gen_tib  = all_gen_tib,
      
      # Field Parameters:: general
      # ids_key = run$ids_key, 
      ids_key = run$unq_key, 
      din_key = run$din_key, 
      pos_key = run$pos_key,
      chr_key = run$chr_key,
      
      # Field Parameters:: s-improbe
      ext_seq = run$ext_seq,
      iup_seq = run$iup_seq,
      imp_seq = run$imp_seq,
      
      # Field Parameters:: r-improbe
      srsplit = run$srsplit,
      srd_key = run$srd_key,
      srd_str = run$srd_str,
      
      cosplit = run$cosplit,
      cos_key = run$cos_key,
      cos_str = run$cos_str,
      
      # Docker Parameters::
      doc_image = run$doc_image,
      doc_shell = run$doc_shell,
      
      join     = FALSE,
      join_new = c("Aln_Key_Unq","Bsp_Chr","Bsp_Pos","Bsp_FR","Bsp_CO"),
      join_old = c("Seq_ID","Chromosome","Coordinate","Strand_FR","Strand_CO"),
      
      subset   = TRUE,
      sub_cols = NULL,
      
      # reload=opt$reload,
      reload  = TRUE,
      retData = TRUE,
      # retData=FALSE,
      
      parallel=opt$parallel,
      # parallel=FALSE,
      
      r_improbe = TRUE,
      s_improbe = FALSE,
      c_improbe = TRUE,
      
      add_flanks = TRUE,
      add_matseq = TRUE,
      
      verbose=opt$verbose, tt=pTracker
    )
  }
  
  
  
  
  
  
  
  
  if (FALSE) {
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                         4.0 improbe fwd design::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    stamp_vec <- c(stamp_vec,
                   run$imp_inp_tsv,
                   run$imp_des_tsv,
                   run$imp_fin_tsv)
    
    if (opt$fresh || !valid_time_stamp(stamp_vec)) {
      
      # TBD:: Add a safty check for Aln_Key_Unq > 15 characters!!!
      # TBD:: Add converted genomes and substring probe checking!!!
      # TBD:: Add r-improbe designer
      #
      
      #    aqp_imp_list <- imp_designs_workflow(
      aqp_imp_tib_pre <- aqp_imp_tib
      aqp_imp_tib <- imp_designs_workflow(
        tib=aqp_cgn_tib, # max=20,
        
        imp_inp_tsv=run$imp_inp_tsv,
        imp_des_tsv=run$imp_des_tsv,
        imp_fin_tsv=run$imp_fin_tsv,
        imp_seq_csv=run$imp_seq_csv,
        
        # Genomes::
        gen_bld=opt$genBuild,
        
        gen_ref_fas=run$gen_ref_fas,
        bsc_ref_fas=run$bsc_ref_fas,
        gen_snp_fas=run$gen_snp_fas,
        bsc_snp_fas=run$bsc_snp_fas,
        
        # Field Parameters::
        ids_key="Aln_Key_Unq",
        din_key="Ord_Din",
        pos_key="Bsp_Pos",
        chr_key="Bsp_Chr",
        
        srsplit=TRUE,
        srd_key="Bsp_FR",
        cosplit=TRUE,
        cos_key="Bsp_CO",
        
        # Docker Parameters::
        run_name=opt$runName,
        doc_image=image_str,
        doc_shell=image_ssh,
        
        join_new=c("Aln_Key_Unq","Bsp_Chr","Bsp_Pos","Bsp_FR","Bsp_CO"),
        join_old=c("Seq_ID","Chromosome","Coordinate","Strand_FR","Strand_CO"),
        
        subset=TRUE,
        sub_cols=NULL,
        
        # reload=opt$reload,
        reload=FALSE,
        retData=TRUE,
        # retData=FALSE,
        
        parallel=opt$parallel,
        # parallel=FALSE,
        r_improbe=TRUE,
        s_improbe=FALSE,
        
        add_flanks=TRUE,
        add_matseq=TRUE,
        
        verbose=opt$verbose, tt=pTracker)
      
    } else {
      
      aqp_imp_tib <- safe_read(
        run$imp_fin_tsv, funcTag="aqp-imp", clean=TRUE,
        verbose=opt$verbose,tt=pTracker)
      
      if (FALSE) {
        imp_fmt_tib <- aqp_imp_tib %>% purrr::set_names(
          c("Imp_Cgn","Imp_Bld","Imp_Chr","Imp_Pos","Imp_Fwd","Imp_Top",
            "Imp_PrbT","Imp_PrbU","Imp_PrbM",
            "Imp_FR","Imp_TB","Imp_CO","Imp_Nxt",
            "Imp_ScrU","Imp_ScrM",
            "Imp_CpgCnt","Imp_CpgDis","Imp_Scr","Imp_NxbScr",
            "Imp_ScrMin","Imp_Inf") )
      }
    }
    
    # Qucik QC::
    #  aqp_imp_tib %>% dplyr::filter(Inf_Type != 0) %>% dplyr::distinct(Seq_ID)
    #  aqp_imp_tib %>% dplyr::group_by(Ord_Des,Ord_Din) %>% dplyr::summarise(Count=n(), .groups = "drop")
    #
  }
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #
  #                              END OF ROUND 1::
  #
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  
  # Load Genome::
  #
  chr_key <- "Bsp_Chr"
  gen_ref_fas <- run$gen_ref_fas
  nrec <- 1
  
  
  dna_dat <- load_genome(file=gen_ref_fas, nrec=nrec,
                         chr_key=chr_key, ret_map=TRUE,
                         verbose=opt$verbose)
  
  prb_seqs <- min_bsp_tib %>% 
    dplyr::mutate(Prb_Beg=dplyr::case_when(
      Bsp_FR=="F" & Bsp_CO=="O" & Ord_Des=="U" ~ Bsp_Pos + 0.0,
      Bsp_FR=="F" & Bsp_CO=="O" & Ord_Des=="M" ~ Bsp_Pos + 0.0,
      Bsp_FR=="F" & Bsp_CO=="O" & Ord_Des=="2" ~ Bsp_Pos + 1.0,
      TRUE ~ NA_real_ ) ) %>%
    dplyr::filter(Bsp_FR=="F") %>% 
    dplyr::filter(Bsp_CO=="O") %>% 
    parse_genomic_seqs(seq = as.character(dna_dat$seqs[1]), 
                       srd_str = "F",
                       pos_key = "Prb_Beg",
                       ups_len = 0,
                       seq_len = 50,
                       pad_len = 0,
                       chr_str = "chr1") %>% revCmp() 
  
  dplyr::mutate(Prb_Seq=prb_seqs)
  
  #
  # - fas_to_seq() reference should be moved into bsp workflow
  # - setTopBot_tib() on reference seqs 
  # - validate TB calls against improbe
  #
  
  #
  # - Select Probe Order Design from s_improbe
  # - Select U/M to compare against improbe
  #
  
  #
  #  N => r_improbe()
  # !N => s_improbe()
  #
  
  sel_cols <- c("Aln_Key_Unq", 
                "Ord_Des","Ord_Din", 
                "Bsp_Chr","Bsp_Pos","Bsp_Cgn","Bsp_FR","Bsp_CO","Ord_Prb")
  
  all_prb_tib <- NULL
  min_bsp_tib <- aqp_cgn_tib %>% dplyr::select(dplyr::all_of(sel_cols))
  all_gen_list <- all_gen_tib %>% split(f=all_gen_tib$Genome_Key)
  #for (ii in c(1:base::nrow(all_gen_tib))) {
  for (gen_key in names(all_gen_list)) {
    # gen_key  <- all_gen_tib[[]]$Genome_Key
    gen_fas  <- all_gen_list[[gen_key]]$Path
    gen_FR   <- all_gen_list[[gen_key]]$Strand_FR
    gen_CO   <- all_gen_list[[gen_key]]$Strand_CO
    gen_BSC  <- all_gen_list[[gen_key]]$Strand_BSC
    gen_prep <- all_gen_list[[gen_key]]$Alphabet
    
    if (gen_prep=="N") r_improbe_val <- TRUE
    if (gen_prep=="N") s_improbe_val <- FALSE
    if (gen_prep!="N") r_improbe_val <- FALSE
    if (gen_prep!="N") s_improbe_val <- TRUE
    
    cat(glue::glue("[{par$prgmTag}]: Genome({gen_key})={gen_fas}...{RET}"))
    
    gen_key <- "RCM_dna"
    # all_prb_tib[[gen_key]] <- fas_to_seq(tib=min_bsp_tib,
    RCM_dna_revCmp_tib <- fas_to_seq(tib=min_bsp_tib, srd_str = "R",
                                     gen_bld=paste(opt$genBuild,gen_key, sep="_"),
                                     
                                     gen_ref_fas=gen_fas,
                                     # imp_tsv=imp_inp_tsv,
                                     # seq_csv=imp_seq_csv,
                                     
                                     build=c("Prb1C","Prb2C","Prb1O","Prb2O"),
                                     
                                     ids_key="Aln_Key_Unq",
                                     din_key="Ord_Din",
                                     tar_din="rs",
                                     pos_key="Bsp_Pos",
                                     chr_key="Bsp_Chr",
                                     
                                     # subset=subset,
                                     # sub_cols=sub_cols,
                                     
                                     reload=FALSE,
                                     retData=FALSE,
                                     parallel=opt$parallel,
                                     r_improbe=r_improbe_val,
                                     s_improbe=s_improbe_val,
                                     add_flanks=TRUE,
                                     
                                     verbose=opt$verbose+1,tt=pTracker)
    
    cat(glue::glue("[{par$prgmTag}]: Genome({gen_key})={gen_fas}...{RET}{RET}"))
  }
  
  
  
  
  pairwiseAlignment(FCN_r_improbe_tib$Prb_1M %>% revCmp(), FCN_r_improbe_tib$DesBscM %>% stringr::str_to_upper(), type="local-global")
  
  all_prb_tib$RCM_dna %>% dplyr::filter(Ord_Des=="M") %>% dplyr::select(Aln_Key_Unq, Bsp_FR, Bsp_CO, Ord_Prb, Ord_Prb, Prb1C, Temp, Iupac_Forward_Sequence) %>% head() %>% as.data.frame()
  RCM_dna_revCmp_tib %>% dplyr::filter(Ord_Des=="M") %>% dplyr::select(Aln_Key_Unq, Bsp_FR, Bsp_CO, Ord_Prb, Ord_Prb, Prb1C, Temp, Iupac_Forward_Sequence) %>% head() %>% as.data.frame()
  
  FCN_r_improbe_tib <- r_improbe(tib=aqp_imp_tib$s_ref, ids_key="Aln_Key_Unq", seq_key="Forward_Sequence", din_key="Ord_Din", srsplit = TRUE, srd_key = "Bsp_FR", cosplit = TRUE, cos_key = "Bsp_CO", parallel = TRUE, add_matseq = TRUE, verbose = opt$verbose, tt=pTracker)
  
  print_prbs(tib = FCN_r_improbe_tib, 
             tar_des = "cg", 
             ids_key = "Aln_Key_Unq", 
             prb_key = "Ord_Prb",
             des_key = "Ord_Des", 
             din_key = "Ord_Din",
             outDir = file.path(opt$outDir),
             plotName = "FCN_dna",
             verbose = opt$verbose+10)
  
  FCN_r_improbe_tib %>% head(n=1) %>%
    dplyr::rename(
      StrandFR=Strand_FR,
      StrandCO=Strand_CO,
    ) %>%
    prbsToStr(pr="cg", verbose = opt$verbose+5)
  
  all_prb_tib$FCM_iupac %>% 
    cmpInfII_MisMatch(fieldA="Ord_Prb", 
                      fieldB="Prb1C", 
                      verbose=opt$verbose) %>% 
    dplyr::group_by(Ord_Des, Ord_Din, Bsp_FR, Bsp_CO, Man_MisMatch) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  
  
  
  
  
  
  
  
  aqp_top_tib <- setTopBot_tib(aqp_imp_tib$s_ref, 
                               seqKey="Forward_Sequence", 
                               srdKey="Bsp_TB2", 
                               verbose=opt$verbose,tt=pTracker)
  
  aqp_top_tib <- setTopBot_tib(aqp_imp_tib$s_ref, 
                               seqKey="Forward_Sequence", 
                               srdKey="Bsp_TB2", 
                               verbose=opt$verbose,tt=pTracker)
  
  imp_top_tib <- setTopBot_tib(aqp_imp_tib$i_imp, 
                               seqKey="Forward_Sequence_imp", 
                               srdKey="Bsp_TB2",
                               verbose=opt$verbose,tt=pTracker)
  
  all_top_tib <- setTopBot_tib(aqp_imp_tib$i_imp,
                               seqKey="Top_Sequence", 
                               srdKey="Bsp_TB2",
                               verbose=opt$verbose,tt=pTracker)
  
  #
  # Design Steps::
  #  gen_ref_fas
  #    - Don't build probes (--build_prbs=FALSE)
  #
  #  r-improbe(Forward_Sequence=Imp_Temp_Seq)
  #    - Compare Prb1_[U/M] against improbe [U/M]
  #    - Compare by Ord_Des (2/U/M) against Ord_Prb
  #
  #  s-improbe (substring improbe) [BSC U/M/D]
  #    - Compare Prb1_[U/M] against improbe [U/M]
  #    - Compare by Ord_Des (2/U/M) against Ord_Prb
  #
  # SNP Check::
  #    - Pos:Iupac (Include Next Base: pos=0)
  #
  # Update Probe_ID
  #    - rs/ch database
  #    - mu = multiple zero mismatch hits
  #    - ma = multiple non-zero mismatch hits
  #    - um = un-paired Infinium I probes
  #
  # Extension/Color Summary::
  #    - Extension Summary (Cpg, Nxb)
  #    - Color Summary (Red, Grn)
  #
  # Clean-Up Steps::
  #  - Remove/Rename Temp_Seq
  #  - Remove *_Len
  #  - Only build probes when nescessary
  #
  
  tmp_join_vec <- c("Aln_Key_Unq", "Address","Ord_Des", "Ord_Din")
  
  tmp_join_tib <- dplyr::inner_join(
    aqp_imp_list$fwd, aqp_imp_list$ret, 
    by=tmp_join_vec,
    suffix=c("_fwd","_imp") )
  
  
  
  
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       1.0 Write improbe input::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  #
  # TBD:: 
  #  - Compare designs from r-improbe to substr-BSC-genomes!!!
  #
  
  nrec <- 2
  
  fwd_seq_tsv <- file.path(par$topDir, "tmp/imp.fwd-seq.tsv.gz")
  snp_seq_tsv <- file.path(par$topDir, "tmp/imp.snp-seq.tsv.gz")
  fwd_des_tsv <- file.path(par$topDir, "tmp/test.improbe-designOutput.tsv.gz")
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   0.0 Get Order Probes for Comparison::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ord_imp_tib <- aqp_imp_tib %>% 
    dplyr::mutate(Seq_ID=Aln_Key_Unq,
                  Len=stringr::str_length(Seq_ID))
  # Sanity Check length(Seq_ID) <= 15::
  ord_imp_tib %>% dplyr::select(Seq_ID, Len) %>% dplyr::arrange(-Len)
  
  ord_prb_tib <- ord_imp_tib %>% 
    dplyr::select(Seq_ID,Ord_Des,Ord_Din, 
                  Bsp_FR,Strand_Ref_FR,Bsp_CO, Ord_Prb) %>%
    dplyr::rename(Strand_FR=Strand_Ref_FR,
                  Strand_CO=Bsp_CO)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       1.0 Get Templates from Genome::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  fwd_seq_tib <- ord_imp_tib %>% 
    dplyr::arrange(Bsp_Chr,Bsp_Pos) %>%
    fas_to_seq(fas=run$gen_ref_fas, file=fwd_seq_tsv, 
               name="Seq_ID", din="Ord_Din", 
               gen=opt$genBuild,
               chr1="Bsp_Chr", pos="Bsp_Pos",
               nrec=nrec, verbose=opt$verbose)
  
  # snp_seq_org <- snp_seq_tib
  # snp_seq_org2 <- snp_seq_tib
  
  snp_seq_tib <- ord_imp_tib %>% 
    dplyr::arrange(Bsp_Chr,Bsp_Pos) %>%
    fas_to_seq(fas=run$gen_snp_fas, file=snp_seq_tsv, 
               name="Seq_ID", din="Ord_Din", 
               gen=opt$genBuild,
               chr1="Bsp_Chr", pos="Bsp_Pos",
               nrec=nrec, verbose=opt$verbose)
  
  # Compare probes against designs::
  #
  fwd_seq_tib %>% dplyr::filter(Bsp_FR=="F", Bsp_CO=="C") %>% 
    dplyr::select(Aln_Key_Unq, Probe_Seq_T,Prb1_FC)
  snp_seq_tib %>% dplyr::filter(Bsp_FR=="F", Bsp_CO=="C") %>% 
    dplyr::select(Aln_Key_Unq, Probe_Seq_T,Prb1_FC)
  
  #
  # Get Converted Genomes::
  #
  
  cur_gen_dir <- file.path(opt$genDir, opt$genBuild,"Sequence/WholeGenomeFasta")
  ref_gen_pat <- paste0(opt$genBuild,".genome.[FR]C[MUD].fa.gz$")
  bsc_gen_pat <- paste0(opt$genBuild,".genome.dbSNP-151.iupac.[FR]C[MUD].fa.gz$")
  
  ref_bsc_files <- list.files(cur_gen_dir, pattern=ref_gen_pat, full.names=TRUE)
  snp_bsc_files <- list.files(cur_gen_dir, pattern=bsc_gen_pat, full.names=TRUE)
  
  #
  # Validate Comparison of difference (i.e. contains SNPs)
  #
  if (FALSE) {
    fwdSnp_seq_tib <- dplyr::inner_join(fwd_seq_tib,snp_seq_tib, 
                                        by="Seq_ID",
                                        suffix=c("_fwd","_snp"))
    
    fwdSnp_seq_tib %>% 
      dplyr::filter(Fwd_Temp_Seq_fwd != Fwd_Temp_Seq_snp) %>% 
      dplyr::select(Fwd_Temp_Seq_fwd,Fwd_Temp_Seq_snp)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         2.0 improbe/docker::fwd
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ret_val <- 
    run_improbe_docker(
      file=fwd_seq_tsv, 
      name="test", image=image_str, shell=image_ssh,
      verbose=opt$verbose)
  
  fwd_ides_tib <-
    load_improbe_design(
      file=fwd_des_tsv, join=NULL, out=NULL,
      level=3, add_inf=TRUE,
      verbose=opt$verbose)
  fwd_iprb_tib <- fwd_ides_tib %>% 
    dplyr::select(Seq_ID,Strand_Ref_FR,Strand_TB,
                  Strand_CO,Probe_Seq_U,Probe_Seq_M) %>% 
    dplyr::rename(Strand_FR=Strand_Ref_FR) %>%
    dplyr::left_join(dplyr::select(ord_imp_tib, Seq_ID,Ord_Des,Ord_Din),
                     by="Seq_ID", suffix=c("_fwd","_ord")) %>%
    dplyr::select(Seq_ID,Ord_Des,Ord_Din, dplyr::everything())
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           3.0 r-improbe::fwd
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  fwd_rdes_tib <- desSeq_to_prbs(
    tib=fwd_seq_tib, 
    ids_key="Seq_ID", seq_key="Imp_Temp_Seq", prb_key="Ord_Din", 
    strsSR="FR", strsCO="CO", 
    parallel=TRUE, verbose=opt$verbose )
  fwd_rprb_tib <- fwd_rdes_tib %>% 
    dplyr::select(Seq_ID, Strand_SR,Strand_CO, 
                  PRB1_U_MAT,PRB1_M_MAT,PRB2_D_MAT) %>%
    dplyr::rename(Strand_FR=Strand_SR) %>%
    dplyr::left_join(dplyr::select(ord_imp_tib, Seq_ID,Ord_Des,Ord_Din),
                     by="Seq_ID", suffix=c("_fwd","_ord")) %>%
    dplyr::select(Seq_ID,Ord_Des,Ord_Din, dplyr::everything())
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           4.0 r-improbe::snp
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  snp_rdes_tib <- desSeq_to_prbs(
    tib=snp_seq_tib,
    ids_key="Seq_ID", seq_key="Iup_Temp_Seq", prb_key="Ord_Din", 
    strsSR="FR", strsCO="CO", 
    parallel=TRUE, verbose=opt$verbose )
  snp_rprb_tib <- snp_rdes_tib %>% 
    dplyr::select(Seq_ID, Strand_SR,Strand_CO, 
                  PRB1_U_MAT,PRB1_M_MAT,PRB2_D_MAT) %>%
    dplyr::rename(Strand_FR=Strand_SR) %>%
    dplyr::left_join(dplyr::select(ord_imp_tib, Seq_ID,Ord_Des,Ord_Din),
                     by="Seq_ID", suffix=c("_fwd","_ord")) %>%
    dplyr::select(Seq_ID,Ord_Des,Ord_Din, dplyr::everything())
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           4.0 Compare Results::
  #
  #  - fwd_iprb_tib vs. [fwd_rprb_tib, snp_rprb_tib]
  #  - ord_prb_tib  vs. [fwd_iprb_tib, fwd_rprb_tib, snp_rprb_tib]
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  #
  # TBD:: Functionize Intersect/Comparison Methods::
  #
  
  bys_vec  <- c("Seq_ID","Ord_Des","Ord_Din","Strand_FR","Strand_CO")
  grp_vec  <- c("Ord_Des","Ord_Din","Strand_FR","Strand_CO")
  
  ref_keys <- c("Probe_Seq_U","Probe_Seq_M")
  can_keys <- c("PRB1_U_MAT","PRB1_M_MAT")
  
  ir_fwd_cmp <- compare_probes(fwd_iprb_tib,fwd_rprb_tib, 
                               ref_keys=ref_keys, can_keys=can_keys,
                               by=bys_vec, grp=grp_vec,
                               retData=TRUE,
                               verbose=opt$verbose)
  
  ref_keys <- c("PRB1_U_MAT","PRB1_M_MAT","PRB2_D_MAT")
  can_keys <- c("PRB1_U_MAT","PRB1_M_MAT","PRB2_D_MAT")
  
  ii_fwdSnp_cmp <- compare_probes(fwd_rprb_tib,snp_rprb_tib,
                                  ref_keys=ref_keys, can_keys=can_keys,
                                  by=bys_vec, grp=grp_vec,
                                  retData=TRUE, pivot=TRUE,
                                  verbose=opt$verbose)
  
  ref_keys <- c("Ord_Prb")
  can_keys <- c("Probe_Seq")
  
  snp_rprb_tab <- dplyr::bind_rows(
    dplyr::filter(snp_rprb_tib,Ord_Des=="U") %>% dplyr::select(-PRB1_M_MAT,-PRB2_D_MAT) %>% dplyr::rename(Probe_Seq=PRB1_U_MAT),
    dplyr::filter(snp_rprb_tib,Ord_Des=="M") %>% dplyr::select(-PRB1_U_MAT,-PRB2_D_MAT) %>% dplyr::rename(Probe_Seq=PRB1_M_MAT),
    dplyr::filter(snp_rprb_tib,Ord_Des=="2") %>% dplyr::select(-PRB1_U_MAT,-PRB1_M_MAT) %>% dplyr::rename(Probe_Seq=PRB2_D_MAT)
  )
  
  oi_ordSnp_cmp <- compare_probes(ord_prb_tib,snp_rprb_tab,
                                  ref_keys=ref_keys, can_keys=can_keys,
                                  by=bys_vec, grp=grp_vec,
                                  retData=TRUE, pivot=TRUE,
                                  verbose=opt$verbose)
  
  
  
  snp_rprb_tib_lc <- snp_rprb_tib %>% dplyr::mutate(PRB1_U_MAT=stringr::str_to_lower(PRB1_U_MAT),PRB1_M_MAT=stringr::str_to_lower(PRB1_M_MAT),PRB2_D_MAT=stringr::str_to_lower(PRB2_D_MAT))
  
  mapply(list.string.diff, fwd_rprb_tib$PRB1_U_MAT,snp_rprb_tib$PRB1_U_MAT)
  
  mapply(list.string.diff, fwd_rprb_tib$PRB1_U_MAT,snp_rprb_tib$PRB1_U_MAT)
  
  # Pivot for human viewing::
  ii_fwdSnp_cmp$cmp %>% 
    # dplyr::filter(Ord_Des=="M" & Man_MisMatch==3) %>% head() %>% 
    tidyr::pivot_longer(cols=c("Probe_A","Probe_B"), 
                        names_to="Prb_Source", values_to="Probe_Seq") %>% 
    as.data.frame()
  
  
  
  
  
  # r-improbe:: improbe_design_all()
  #
  
  #
  # TBD::
  #
  #  Manifest Generation::
  #    - Code Clean Up
  #       - [done] bsp_mapping_workflow()
  #
  #    - Calculate extension/color distribution
  #    - Extract BSC Top/Probe-Design
  #
  #    - Join by position
  #    - Add masked-controls from source
  #
  #  Annotation::
  #    - Incorporate
  #    - Annotation Summary
  #    - Implement ftp for missing files
  #
  #  Cluster::
  #    - Transfer all files to cluster
  #
  
  
  
  
  assign_cgn_old = function(ord_tib,
                            bsp_tib,
                            seq_tib,
                            
                            can_csv,
                            can_tib = NULL,
                            bsp_csv = NULL,
                            
                            ids_key = run$ids_key,  # Use to be Aln_Key
                            
                            join    = "inner",
                            merge   = TRUE, 
                            retData = FALSE, 
                            
                            verbose=0,vt=3,tc=1,tt=NULL,
                            funcTag='assign_cgn_old') {
    
    tabs <- paste0(rep(TAB, tc), collapse='')
    mssg <- glue::glue("[{funcTag}]:{tabs}")
    
    if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
    if (verbose>=vt+4) {
      cat(glue::glue("{RET}"))
      cat(glue::glue("{mssg} File IO Parameters::{RET}"))
      cat(glue::glue("{mssg}    can_csv={can_csv}.{RET}"))
      cat(glue::glue("{mssg}    bsp_csv={bsp_csv}.{RET}"))
      cat(glue::glue("{RET}"))
      
      cat(glue::glue("{mssg} Function Parameters::{RET}"))
      cat(glue::glue("{mssg}      ids_key={ids_key}.{RET}"))
      cat(glue::glue("{mssg}      join={join}.{RET}"))
      cat(glue::glue("{mssg}     merge={merge}.{RET}"))
      cat(glue::glue("{mssg}   retData={retData}.{RET}"))
      cat(glue::glue("{RET}"))
    }
    
    ret_cnt <- 0
    ret_tib <- NULL
    ret_dat <- NULL
    # stime <- base::system.time({
    
    ids_sym <- rlang::sym(ids_key)
    
    # Load Canonical CGNs::
    if (is.null(can_tib)) {
      can_tib <- 
        safe_read(file=can_csv, 
                  verbose=verbose, vt=vt+1,tc=tc+1,tt=tt) %>% 
        dplyr::select(CGN) %>% 
        dplyr::rename(Cgn=CGN) %>% 
        dplyr::mutate(Can_Cnt=1) %>%
        clean_tibble()
    }
    ret_key <- glue::glue("can_tib({funcTag})")
    ret_cnt <- print_tib(can_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    if (retData) ret_dat$can_tib <- can_tib
    
    # Defined Order tib to ord original cgn::
    ord_tib <- ord_tib %>% 
      dplyr::select(!!ids_sym,Ord_Cgn,Ord_Map) %>%
      dplyr::rename(Cgn=Ord_Cgn) %>% 
      dplyr::mutate(Ord_Cnt=1) %>%
      dplyr::distinct() %>%
      clean_tibble()
    
    ret_key <- glue::glue("ord_tib({funcTag})")
    ret_cnt <- print_tib(ord_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    if (retData) ret_dat$ord_tib <- ord_tib
    
    # Format BSP::
    bsp_tib <- bsp_tib %>% 
      dplyr::filter(!is.na(Bsp_Cgn)) %>% 
      dplyr::select(Address,Ord_Map, !!ids_sym, Ord_Des, Ord_Din, Bsp_Cgn) %>% 
      dplyr::rename(Cgn=Bsp_Cgn) %>%
      dplyr::distinct() %>%
      dplyr::arrange(!!ids_sym, Cgn) %>% 
      dplyr::group_by(Address,Ord_Map,!!ids_sym,Ord_Des,Ord_Din,Cgn) %>% 
      dplyr::summarise(Bsp_Cnt=n(), .groups = "drop") %>%
      clean_tibble()
    
    bsp_key <- glue::glue("bsp_tib({funcTag})")
    bsp_cnt <- print_tib(bsp_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    if (retData) ret_dat$bsp_tib <- bsp_tib
    
    # return(ret_dat)
    
    seq_tib <- seq_tib %>% 
      dplyr::filter(!is.na(Imp_Cgn)) %>% 
      dplyr::left_join(dplyr::distinct(bsp_tib,Address, Ord_Map), by=c("Address")) %>%
      dplyr::select(Address, Ord_Des, Ord_Din, Imp_Cgn) %>% 
      tidyr::unite(Tmp_Key, Ord_Des,Ord_Din, sep="", remove=FALSE) %>%
      tidyr::unite(!!ids_sym, Address, Tmp_Key, sep="_", remove=FALSE) %>%
      dplyr::select(-Tmp_Key) %>%
      # dplyr::left_join(dplyr::select(ord_tib, Ord_Map,!!ids_sym), by=ids_key) %>%
      dplyr::select(Ord_Map,!!ids_sym,Ord_Des,Ord_Din,Imp_Cgn) %>%
      dplyr::rename(Cgn=Imp_Cgn) %>%
      dplyr::distinct() %>%
      dplyr::arrange(!!ids_sym, Cgn) %>% 
      dplyr::group_by(Ord_Map,!!ids_sym,Ord_Des,Ord_Din,Cgn) %>% 
      dplyr::summarise(Seq_Cnt=n(), .groups = "drop") %>%
      clean_tibble()
    
    seq_key <- glue::glue("seq-tib({funcTag})")
    seq_cnt <- print_tib(seq_tib,funcTag, verbose,vt+4,tc, n=seq_key)
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    if (retData) ret_dat$seq_tib <- seq_tib
    
    # Build and Sort Counts Tables
    cnt_tib <- 
      dplyr::full_join(bsp_tib, seq_tib, 
                       by=c("Ord_Map",ids_key,"Ord_Des","Ord_Din","Cgn")) %>% 
      dplyr::left_join(can_tib, by="Cgn") %>%
      dplyr::distinct() %>%
      dplyr::left_join(ord_tib, by=c(ids_key,"Cgn")) %>%
      dplyr::distinct() %>%
      dplyr::mutate(dplyr::across(c(Bsp_Cnt,Seq_Cnt,Can_Cnt,Ord_Cnt), tidyr::replace_na, 0 ),
                    Sum_Cnt=Bsp_Cnt+Seq_Cnt,
                    Max_Cnt=Bsp_Cnt*Seq_Cnt) %>% 
      dplyr::add_count(!!ids_sym, name="Cgn_Cnt") %>% 
      dplyr::arrange(-Can_Cnt,-Max_Cnt,-Sum_Cnt,-Ord_Cnt) %>%
      dplyr::mutate(Rank=dplyr::row_number()) %>%
      clean_tibble()
    
    seq_key <- glue::glue("cnt_tib({funcTag})")
    seq_cnt <- print_tib(cnt_tib,funcTag, verbose,vt+4,tc, n=seq_key)
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    
    if (retData) ret_dat$cnt_tib <- cnt_tib
    
    cnt_list <- cnt_tib %>% split(.$Ord_Des)
    
    cnt_tib_cnt <- cnt_tib %>% base::nrow()
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    cat(glue::glue("{mssg}    cnt_tib(n={cnt_tib_cnt})::{RET}"))
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    cat(glue::glue("{mssg}{RET}"))
    cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
    
    # Infinium II::
    inf2_tib <- cnt_list[["2"]] %>%
      dplyr::arrange(!!ids_sym, Rank) %>%
      dplyr::distinct(!!ids_sym, .keep_all = TRUE)
    
    # Infinium I:: Full Join
    #
    #   TBD:: The joining should really be done by sequence: Ord_Prb
    #
    if (join=="full") {
      inf1_tib <- dplyr::full_join(
        cnt_list[["U"]], cnt_list[["M"]], 
        by=c("Ord_Map","Cgn","Ord_Din"), 
        suffix=c("_U","_M")
      ) %>%
        dplyr::mutate(Rank_Min=pmin(Rank_U,Rank_M)) %>%
        dplyr::arrange(Ord_Map, Rank_Min) %>%
        dplyr::distinct(Ord_Map,Aln_Key_U,Aln_Key_M, .keep_all = TRUE)
    } else if (join=="inner") {
      inf1_tib <- dplyr::inner_join(
        cnt_list[["U"]], cnt_list[["M"]], 
        by=c("Ord_Map","Cgn","Ord_Din"), 
        suffix=c("_U","_M")
      ) %>%
        dplyr::mutate(Rank_Min=pmin(Rank_U,Rank_M)) %>%
        dplyr::arrange(Ord_Map, Rank_Min) %>%
        dplyr::distinct(Ord_Map,Aln_Key_U,Aln_Key_M, .keep_all = TRUE)
    } else {
      stop(glue::glue("{mssg} Unsupported join type={join}.{RET}"))
      return(NULL)
    }
    if (retData) ret_dat$inf1_tib <- inf1_tib
    if (retData) ret_dat$inf2_tib <- inf2_tib
    
    ret_tib <- dplyr::bind_rows(
      dplyr::select(inf1_tib, Ord_Map,Aln_Key_U,Cgn,Ord_Des_U,Ord_Din,Can_Cnt_U,Rank_Min) %>% 
        purrr::set_names("Ord_Map",ids_key,"Cgn","Ord_Des","Ord_Din","Can_Cnt","Rank"),
      
      dplyr::select(inf1_tib, Ord_Map,Aln_Key_M,Cgn,Ord_Des_M,Ord_Din,Can_Cnt_M,Rank_Min) %>% 
        purrr::set_names("Ord_Map",ids_key,"Cgn","Ord_Des","Ord_Din","Can_Cnt","Rank"),
      
      dplyr::select(inf2_tib, Ord_Map,!!ids_sym,Cgn,Ord_Des,Ord_Din,Can_Cnt,Rank)
    ) %>% dplyr::filter(!is.na(!!ids_sym)) %>%
      dplyr::distinct()
    
    mul_cnt <- ret_tib %>% dplyr::add_count(!!ids_sym,Cgn, name="Multi_Cnt") %>% 
      dplyr::filter(Multi_Cnt != 1) %>% base::nrow()
    mis_cnt <- ret_tib %>% dplyr::filter(is.na(!!ids_sym)) %>% base::nrow()
    
    mis_tib <- dplyr::anti_join(ord_tib, ret_tib, by=c(ids_key))
    sig_tib <- dplyr::filter(cnt_tib, ids_key %in% mis_tib[[ids_key]]) %>%
      dplyr::arrange(Ord_Map,Rank) %>%
      dplyr::distinct(!!ids_sym, .keep_all = TRUE)
    sig_cnt <- sig_tib %>% base::nrow()
    
    if (verbose>=vt) {
      cat(glue::glue("{mssg}   Miss Count={mis_cnt}.{RET}"))
      cat(glue::glue("{mssg}  Multi Count={mul_cnt}.{RET}"))
      cat(glue::glue("{mssg} Single Count={sig_cnt}.{RET}"))
      cat("\n")
    }
    # if (mis_cnt!=0 || mul_cnt!=0 || sig_cnt!=0) {
    #   stop(glue::glue("{RET}{mssg} Counts non-zero={mis_cnt},",
    #                   "{mul_cnt},{sig_cnt}!{RET2}"))
    #   return(NULL)
    # }
    
    # Merge all data together::
    #
    ret_cnt <- ret_tib %>% base::nrow()
    ret_tib <- dplyr::bind_rows(
      ret_tib %>% dplyr::mutate(
        Cgn_Tag=dplyr::case_when(
          Ord_Din=="rs" ~ Ord_Din,
          Ord_Din=="ch" ~ Ord_Din,
          TRUE ~ "cg"
        ),
        Cgn_Str=dplyr::case_when(
          Ord_Din=="rs" ~ stringr::str_remove(Ord_Map, "[-_:].*$"),
          Ord_Din=="ch" ~ stringr::str_remove(Ord_Map, "[-_:].*$"),
          TRUE ~ paste0("cg",stringr::str_pad(Cgn,width=8,side="left",pad="0"))
        )),
      mis_tib %>% 
        dplyr::select(Ord_Map, !!ids_sym,Ord_Cgn,Ord_Des,Ord_Din) %>% 
        dplyr::rename(Cgn=Ord_Cgn) %>% 
        dplyr::mutate(Can_Cnt=0, 
                      Rank=dplyr::row_number() + ret_cnt,
                      Cgn_Tag="uk",
                      Cgn_Str=paste0(Cgn_Tag,stringr::str_pad(Cgn,width=8,side="left",pad="0"))
        )
    ) %>%
      # TBD:: Capture other CGN's in seperate column:: actual CGN's not Count!!
      dplyr::add_count(!!ids_sym, name="Alt_Cgn_Cnt") %>%
      # One Final Clean Up To Ensure Uniqueness::
      dplyr::arrange(Rank) %>% 
      dplyr::distinct(!!ids_sym, .keep_all = TRUE)
    
    mul_cnt <- ret_tib %>% 
      dplyr::add_count(!!ids_sym,Cgn, name="Multi_Cnt") %>% 
      dplyr::filter(Multi_Cnt != 1) %>% base::nrow()
    
    if (verbose>=vt) {
      cat(glue::glue("{mssg}  Multi Count Final={mul_cnt}.{RET}"))
      cat("\n")
    }
    if (mul_cnt!=0) {
      stop(glue::glue("{RET}{mssg} Multi-Count Final={mul_cnt} ",
                      "not equal to zero!!!{RET2}"))
      return(NULL)
    }
    
    if (merge) ret_tib <- bsp_tib %>%
      dplyr::left_join(ret_tib, 
                       by=c("Ord_Map",ids_key,"Ord_Des","Ord_Din"),
                       suffix=c("_bsp","_cgn"))
    
    ret_tib <- ret_tib %>% clean_tibble()
    out_cnt <- safe_write(ret_tib,file=bsp_csv, funcTag=funcTag,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    if (retData) ret_dat$ret_tib <- ret_tib
    if (retData) ret_dat$mis_tib <- mis_tib
    # })
    etime <- 0
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{mssg}{BRK}{RET2}"))
    
    if (retData) return(ret_dat)
    
    ret_tib
  }
  
  
}


if (FALSE) {
  

  
  run$merge <- FALSE
  cgn_out_csv <- file.path(opt$outDir, "assigned-cgns.csv.gz")
  
  can_tib_bk <- can_tib
  
  ids_sym <- rlang::sym(run$ids_key)
  
  can_tib <- NULL
  ord_tib <- aqp_ord_tib
  bsp_tib = aqp_bsp_tib
  
  seq_tib = aqp_seq_tib
  can_csv = run$canonical_csv
  can_tib = can_tib_bk
  
  ids_key = run$ids_key
  
  bsp_csv = cgn_out_csv
  merge = run$merge
  
  retData <- FALSE
  
  
  tt <- NULL
  opt$verbose <- 5
  verbose <- opt$verbose
  tc <- 1
  vt <- 0
  funcTag <- "Testing"
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  # Load Canonical CGNs::
  if (is.null(can_tib)) {
    can_tib <- safe_read(file=run$canonical_csv, verbose=opt$verbose) %>% 
      dplyr::select(CGN) %>% 
      dplyr::rename(Cgn=CGN) %>% 
      dplyr::mutate(Can_Cnt=1)
  }
  # can_tib = can_tib_bk
  
  
  bsp_key <- glue::glue("can_tib({funcTag})")
  bsp_cnt <- print_tib(can_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
  cat(glue::glue("{mssg}{RET}"))
  cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
  
  if (retData) ret_dat$can_tib <- can_tib
  
  # Defined Order tib to ord original cgn::
  ord_tib <- aqp_ord_tib %>% 
    dplyr::select(Ord_Map,Prb_Key,Ord_Cgn) %>%
    dplyr::rename(Cgn=Ord_Cgn) %>% 
    dplyr::mutate(Ord_Cnt=1) %>%
    dplyr::distinct() %>%
    clean_tibble()
  
  bsp_key <- glue::glue("ord_tib({funcTag})")
  bsp_cnt <- print_tib(ord_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
  cat(glue::glue("{mssg}{RET}"))
  cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
  
  # if (retData) ret_dat$ord_tib <- ord_tib
  
  # Format BSP::
  bsp_tib <- aqp_bsp_tib %>% 
    dplyr::filter(!is.na(Bsp_Cgn)) %>% 
    dplyr::select(Ord_Map, Prb_Key, Ord_Des, Ord_Din, Bsp_Cgn) %>% 
    dplyr::rename(Cgn=Bsp_Cgn) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Prb_Key, Cgn) %>% 
    dplyr::group_by(Ord_Map,Prb_Key,Ord_Des,Ord_Din,Cgn) %>% 
    dplyr::summarise(Bsp_Cnt=n(), .groups = "drop") %>%
    clean_tibble()
  
  bsp_key <- glue::glue("bsp_tib({funcTag})")
  bsp_cnt <- print_tib(bsp_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
  cat(glue::glue("{mssg}{RET}"))
  cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
  
  
  # Format Seq::
  seq_tib <- aqp_seq_tib %>% 
    dplyr::filter(!is.na(Imp_Cgn)) %>% 
    dplyr::left_join(dplyr::distinct(aqp_bsp_tib,Address, Ord_Map), by=c("Address")) %>%
    dplyr::select(Address, Ord_Des, Ord_Din, Ord_Map, Imp_Cgn) %>% 
    tidyr::unite(Tmp_Key, Ord_Des,Ord_Din, sep="", remove=FALSE) %>%
    tidyr::unite(Prb_Key, Address, Tmp_Key, sep="_", remove=FALSE) %>%
    dplyr::select(-Tmp_Key) %>%
    # dplyr::left_join(dplyr::select(ord_tib, Ord_Map,Prb_Key), by=run$ids_key) %>%
    dplyr::select(Ord_Map,Prb_Key,Ord_Des,Ord_Din,Imp_Cgn) %>%
    dplyr::rename(Cgn=Imp_Cgn) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Prb_Key, Cgn) %>% 
    dplyr::group_by(Ord_Map,Prb_Key,Ord_Des,Ord_Din,Cgn) %>% 
    dplyr::summarise(Seq_Cnt=n(), .groups = "drop") %>%
    clean_tibble()
  
  bsp_key <- glue::glue("seq_tib({funcTag})")
  bsp_cnt <- print_tib(seq_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
  cat(glue::glue("{mssg}{RET}"))
  cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
  
  # Build and Sort Counts Tables
  cnt_tib <- 
    dplyr::full_join(
      bsp_tib, seq_tib, by=c("Ord_Map",run$ids_key,"Ord_Des","Ord_Din","Cgn")
    ) %>%
    dplyr::left_join(can_tib, by="Cgn") %>%
    dplyr::distinct() %>%
    dplyr::left_join(ord_tib, by=c(run$ids_key,"Cgn","Ord_Map")
    ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(dplyr::across(c(Bsp_Cnt,Seq_Cnt,Can_Cnt,Ord_Cnt), tidyr::replace_na, 0 ),
                  Sum_Cnt=Bsp_Cnt+Seq_Cnt,
                  Max_Cnt=Bsp_Cnt*Seq_Cnt) %>% 
    dplyr::add_count(!!ids_sym, name="Cgn_Cnt") %>% 
    dplyr::arrange(-Can_Cnt,-Max_Cnt,-Sum_Cnt,-Ord_Cnt) %>%
    dplyr::mutate(Rank=dplyr::row_number()) %>%
    clean_tibble()
  
  bsp_key <- glue::glue("cnt_tib({funcTag})")
  bsp_cnt <- print_tib(cnt_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
  cat(glue::glue("{mssg}{RET}"))
  cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
  
  
  cnt_list <- cnt_tib %>% split(.$Ord_Des)
  
  #
  # Infinium II::
  #
  inf2_tib <- cnt_list[["2"]] %>%
    dplyr::arrange(!!ids_sym, Rank) %>%
    dplyr::distinct(!!ids_sym, .keep_all = TRUE)
  
  bsp_key <- glue::glue("inf2_tib({funcTag})")
  bsp_cnt <- print_tib(inf2_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
  cat(glue::glue("{mssg}{RET}"))
  cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
  
  # Infinium I:: Full Join
  #
  #   TBD:: The joining should really be done by sequence: Ord_Prb
  #
  run$join <- "inner"
  if (run$join=="full") {
    inf1_tib <- dplyr::full_join(
      cnt_list[["U"]], 
      cnt_list[["M"]], 
      by=c("Ord_Map","Cgn","Ord_Din"), 
      suffix=c("_U","_M")
    ) %>%
      dplyr::mutate(Rank_Min=pmin(Rank_U,Rank_M)) %>%
      dplyr::arrange(Ord_Map, Rank_Min) %>%
      dplyr::distinct(Ord_Map,Prb_Key_U,Prb_Key_M, .keep_all = TRUE)
  } else if (run$join=="inner") {
    inf1_tib <- dplyr::inner_join(
      cnt_list[["U"]], 
      cnt_list[["M"]], 
      by=c("Ord_Map","Cgn","Ord_Din"), 
      suffix=c("_U","_M")
    ) %>%
      dplyr::mutate(Rank_Min=pmin(Rank_U,Rank_M)) %>%
      dplyr::arrange(Ord_Map, Rank_Min) %>%
      dplyr::distinct(Ord_Map,Prb_Key_U,Prb_Key_M, .keep_all = TRUE)
  } else {
    stop(glue::glue("{mssg} Unsupported join type={join}.{RET}"))
    return(NULL)
  }
  bsp_key <- glue::glue("inf1_tib({funcTag})")
  bsp_cnt <- print_tib(inf1_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
  cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
  
  # retData <- FALSE
  # if (retData) ret_dat$inf1_tib <- inf1_tib
  # if (retData) ret_dat$inf2_tib <- inf2_tib
  
  ret_tib <- dplyr::bind_rows(
    dplyr::select(inf1_tib, Ord_Map,Prb_Key_U,Cgn,Ord_Des_U,Ord_Din,Can_Cnt_U,Rank_Min) %>% 
      purrr::set_names("Ord_Map",run$ids_key,"Cgn","Ord_Des","Ord_Din","Can_Cnt","Rank"),
    
    dplyr::select(inf1_tib, Ord_Map,Prb_Key_M,Cgn,Ord_Des_M,Ord_Din,Can_Cnt_M,Rank_Min) %>% 
      purrr::set_names("Ord_Map",run$ids_key,"Cgn","Ord_Des","Ord_Din","Can_Cnt","Rank"),
    
    dplyr::select(inf2_tib, Ord_Map,!!ids_sym,Cgn,Ord_Des,Ord_Din,Can_Cnt,Rank)
  ) %>% dplyr::filter(!is.na(!!ids_sym)) %>%
    dplyr::distinct()
  bsp_key <- glue::glue("ret_tib({funcTag})")
  bsp_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
  cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
  
  
  mul_cnt <- ret_tib %>% dplyr::add_count(!!ids_sym,Cgn, name="Multi_Cnt") %>% 
    dplyr::filter(Multi_Cnt != 1) %>% base::nrow()
  mis_cnt <- ret_tib %>% dplyr::filter(is.na(!!ids_sym)) %>% base::nrow()
  
  # bsp_key <- glue::glue("ret_tib({funcTag})")
  # bsp_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
  
  mis_tib <- dplyr::anti_join(ord_tib, ret_tib, by=c(run$ids_key)) %>%
    dplyr::left_join(aqp_ord_tib %>% dplyr::select(Prb_Key,Ord_Des,Ord_Din), by=c("Prb_Key"))
  sig_tib <- dplyr::filter(cnt_tib, run$ids_key %in% mis_tib[[run$ids_key]]) %>%
    dplyr::arrange(Ord_Map,Rank) %>%
    dplyr::distinct(!!ids_sym, .keep_all = TRUE)
  sig_cnt <- sig_tib %>% base::nrow()
  
  bsp_key <- glue::glue("ret_tib({funcTag})")
  bsp_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=bsp_key)
  cat(glue::glue("{mssg}{RET}"))
  cat(glue::glue("{mssg}{RET}{mssg}{BRK}{RET2}"))
  
  if (verbose>=vt) {
    cat(glue::glue("{mssg}   Miss Count={mis_cnt}.{RET}"))
    cat(glue::glue("{mssg}  Multi Count={mul_cnt}.{RET}"))
    cat(glue::glue("{mssg} Single Count={sig_cnt}.{RET}"))
    cat("\n")
  }
  # if (mis_cnt!=0 || mul_cnt!=0 || sig_cnt!=0) {
  #   stop(glue::glue("{RET}{mssg} Counts non-zero={mis_cnt},",
  #                   "{mul_cnt},{sig_cnt}!{RET2}"))
  #   return(NULL)
  # }
  
  # Merge all data together::
  #
  ret_cnt <- ret_tib %>% base::nrow()
  ret_tib <- dplyr::bind_rows(
    #
    # Add Formatted cg#'s
    #
    ret_tib %>% 
      dplyr::mutate(
        
        Cgn_Tag=dplyr::case_when(
          Ord_Din=="rs" ~ Ord_Din,
          Ord_Din=="ch" ~ Ord_Din,
          TRUE ~ "cg"
        ),
        Cgn_Str=dplyr::case_when(
          Ord_Din=="rs" ~ stringr::str_remove(Ord_Map, "[-_:].*$"),
          Ord_Din=="ch" ~ stringr::str_remove(Ord_Map, "[-_:].*$"),
          TRUE ~ paste0("cg",stringr::str_pad(Cgn,width=8,side="left",pad="0"))
        )
      ),
    mis_tib %>% 
      dplyr::select(Ord_Map, !!ids_sym,Cgn,Ord_Des,Ord_Din) %>% 
      
      # Ord_Cgn -> Cgn naming is done above
      # dplyr::select(Ord_Map, !!ids_sym,Ord_Cgn,Ord_Des,Ord_Din) %>% 
      # dplyr::rename(Cgn=Ord_Cgn) %>% 
      dplyr::mutate(Can_Cnt=0, 
                    Rank=dplyr::row_number() + ret_cnt,
                    Cgn_Tag="uk",
                    Cgn_Str=paste0(Cgn_Tag,stringr::str_pad(Cgn,width=8,side="left",pad="0")
                    )
      )
  ) %>%
    # TBD:: Capture other CGN's in seperate column:: actual CGN's not Count!!
    dplyr::add_count(!!ids_sym, name="Alt_Cgn_Cnt") %>%
    # One Final Clean Up To Ensure Uniqueness::
    dplyr::arrange(Rank) %>% 
    dplyr::distinct(!!ids_sym, .keep_all = TRUE)
  
  mul_cnt <- ret_tib %>% 
    dplyr::add_count(!!ids_sym,Cgn, name="Multi_Cnt") %>% 
    dplyr::filter(Multi_Cnt != 1) %>% base::nrow()
  
  if (verbose>=vt) {
    cat(glue::glue("{mssg}  Multi Count Final={mul_cnt}.{RET}"))
    cat("\n")
  }
  if (mul_cnt!=0) {
    stop(glue::glue("{RET}{mssg} Multi-Count Final={mul_cnt} ",
                    "not equal to zero!!!{RET2}"))
    return(NULL)
  }
  
  if (merge) ret_tib <- bsp_tib %>%
    dplyr::left_join(ret_tib, 
                     by=c("Ord_Map",ids_key,"Ord_Des","Ord_Din"),
                     suffix=c("_bsp","_cgn"))
  
  ret_tib <- ret_tib %>% clean_tibble()
  out_cnt <- safe_write(ret_tib,file=bsp_csv, funcTag=funcTag,
                        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
  
  ret_key <- glue::glue("ret-FIN({funcTag})")
  ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  
  if (retData) ret_dat$ret_tib <- ret_tib
  if (retData) ret_dat$mis_tib <- mis_tib
  # })
  etime <- 0
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{mssg}{BRK}{RET2}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
  
  
  
  
}
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             FASTA File Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ord_to_fas = function(tib, 
                      prb_key = "Prb_Seq", 
                      add_key = "Address", 
                      des_key = "Prb_Des", 
                      din_key = "Prb_Din",
                      
                      prb_csv = NULL,
                      prb_fas = NULL,
                      
                      out    = NULL,
                      prefix = NULL,
                      suffix = NULL,
                      pre_len = 0,
                      
                      # u49_tsv=NULL,
                      # m49_tsv=NULL,
                      del="_",
                      verbose=0, vt=3,tc=1,tt=NULL, funcTag='ord_to_fas') {
  
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
    
    prb_sym <- rlang::sym(prb_key)
    des_sym <- rlang::sym(des_key)
    din_sym <- rlang::sym(din_key)
    aln_vec <- c(add_key,des_key)
    
    # aln_vec <- c(des_key,din_key)
    #
    # To make all Aln_Key's <=15 characters for improbe we'll 
    #   unite the Aln_Key in two steps to remove one underscore ("_")
    #
    # Old Definition: aln_vec <- c(add_key,des_key,din_key)
    #
    
    # Build Alignment Keys/Seqs and Write Data::
    #
    ret_tib <- tib %>% 
      tidyr::unite(Tmp_Key, dplyr::all_of(aln_vec),sep=del, remove=FALSE) %>%
      tidyr::unite(Aln_Key, Tmp_Key,!!din_sym,sep="", remove=FALSE) %>%
      dplyr::select(-Tmp_Key) %>%
      dplyr::mutate(
        Aln_Prb=deMs(!!prb_sym, uc=TRUE),
        Aln_Rev=revCmp(Aln_Prb),
        Aln_P49=dplyr::case_when(
          !!des_sym == '2' ~ stringr::str_sub(Aln_Prb, 2),
          !!des_sym == 'U' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          !!des_sym == 'M' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::distinct(Aln_Key,Aln_Prb, .keep_all=TRUE) %>%
      clean_tibble()
    
    safe_write(x=ret_tib, file=dat_csv, funcTag=funcTag,
               verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    # Build/Write FASTA File::
    #
    if (!is.null(prb_fas))
      fas_vec <- ret_tib %>%
      dplyr::mutate(fas_line=paste0(">",Aln_Key,"\n",Aln_Prb) ) %>%
      dplyr::pull(fas_line)
    
    safe_write(x=fas_vec, type="line", file=prb_fas, funcTag=funcTag,
               verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    
    # Build U49/M49 Sub-sequences and split by 5' two nucelotide prefix into 
    #  sorted output files. Splitting by prefix makes the join later much 
    #  faster...
    #
    
    u49_tibs <- NULL
    if (!is.null(u49_tsv))
      u49_tib <- ret_tib %>% 
      dplyr::filter(!!des_sym == '2' | !!des_sym=='U') %>%
      dplyr::filter(!is.na(Aln_P49)) %>%
      dplyr::select(Aln_P49,Aln_Key) %>%
      dplyr::arrange(Aln_P49)
    
    m49_tib <- NULL
    if (!is.null(m49_tsv))
      m49_tib <- ret_tib %>% 
      dplyr::filter(!!des_sym=='M') %>%
      dplyr::filter(!is.na(Aln_P49)) %>%
      dplyr::select(Aln_P49,Aln_Key) %>%
      dplyr::arrange(Aln_P49)
    
    if (pre_len>0) {
      u49_tibs <- u49_tib %>% 
        dplyr::mutate(Pre_Nuc=stringr::str_sub(Aln_P49, 1,pre_len)) %>%
        split(f=.$Pre_Nuc)
      
      for (pre_nuc in names(u49_tibs)) {
        prb_csv <- file.path(out, paste())
        
        prb_tsv <- file.path(out, paste(prefix,suffix,pre_nuc,"U49.tsv.gz", sep='.'))
        
        
      }
    }
    
    
    
    safe_write(u49_tib,"tsv",u49_tsv,cols=FALSE, 
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    safe_write(m49_tib,"tsv",m49_tsv,cols=FALSE, 
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

if (FALSE) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Address To Manifest Methods::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # TBD:: Pretty Sure this can be removed::
  add_comb = function(tibA, tibB, field,
                      join,
                      verbose=0,vt=3,tc=1,tt=NULL) {
    funcTag <- 'add_comb'
    tabs <- paste0(rep(TAB, tc), collapse='')   mssg <- glue::glue("[{funcTag}]:{tabs}")
    if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- base::system.time({
      
      ret_tib <- 
        dplyr::inner_join(
          tibA, tibB, 
          by=dplyr::all_of(join),
          suffix=c("_U","_M")
        ) # %>%
      # dplyr::select(dplyr::starts_with(field)) %>%
      # dplyr::distinct() %>%
      # purrr::set_names(c("U","M"))
      
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    ret_tib
  }
  
  add_to_man = function(tib, join, runName,
                        des_key="Ord_Des", pid_key="Ord_Key",
                        rep_key=NULL, rep_val=NULL,
                        col_key=NULL, nxb_key=NULL,
                        csv=NULL, validate=TRUE,
                        verbose=0,vt=3,tc=1,tt=NULL,
                        funcTag='add_to_man') {
    
    tabs <- paste0(rep(TAB, tc), collapse='')
    mssg <- glue::glue("[{funcTag}]:{tabs}")
    
    if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- base::system.time({
      
      des_list <- NULL
      des_list <- tib %>% split(.[[des_key]])
      des_cnt  <- des_list %>% names() %>% length()
      if (verbose>=vt+4) {
        cat(glue::glue("{mssg} des_list[{des_key}]={RET}"))
        print(des_list)
      }
      
      # Build Infinium I::
      ret1_tib <- NULL
      if (!is.null(des_list[["U"]]) && !is.null(des_list[["M"]])) {
        ret1_tib <- dplyr::full_join(
          dplyr::select(des_list[["U"]], -dplyr::all_of(des_key)), 
          dplyr::select(des_list[["M"]], -dplyr::all_of(des_key)), 
          by=dplyr::all_of(join),
          suffix=c("_U","_M")
        ) %>% dplyr::mutate(Infinium_Design=as.integer(1))
        # TBD:: We should allow these "singletons" to pass, but under
        #  a different classification...
        #
        if ("Address_U" %in% names(ret1_tib) && 
            "Address_M" %in% names(ret1_tib)) ret1_tib <- ret1_tib %>% 
            dplyr::filter(!is.na(Address_U) & !is.na(Address_M))
      }
      ret1_cnt <- print_tib(ret1_tib,funcTag, verbose,vt+4,tc, n="InfI")
      
      # Build Infinium II::
      ret2_tib <- NULL
      if (!is.null(des_list[["2"]])) {
        ret2_tib <- dplyr::bind_cols(
          dplyr::select(des_list[["2"]],  dplyr::all_of(join)),
          dplyr::select(des_list[["2"]], -dplyr::all_of(join)) %>% 
            purrr::set_names(paste(names(.),"U", sep="_"))
        ) %>% dplyr::mutate(Infinium_Design=as.integer(2))
        ret2_cnt <- print_tib(ret2_tib,funcTag, verbose,vt+4,tc, n="InfII")
      }
      
      # Bind Infinium I/II into single manifest::
      ret_tib <- dplyr::bind_rows(ret1_tib, ret2_tib) %>%
        dplyr::mutate(
          # Infinium_Design=dplyr::case_when(
          #   is.na(Address_M) ~ 2,
          #   !is.na(Address_M) ~ 1,
          #   TRUE ~ NA_real_
          # ) %>% as.integer(),
          Infinium_Design_Type=dplyr::case_when(
            Infinium_Design==1 ~ 'I',
            Infinium_Design==2 ~ 'II',
            TRUE ~ NA_character_)
        )
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="InfI+II.0")
      
      if (!is.null(col_key)) {
        col_sym <- rlang::sym(col_key)
        ret_tib <- ret_tib %>%
          dplyr::mutate(
            col=dplyr::case_when(
              Infinium_Design==1 ~ !!col_sym,
              TRUE ~ NA_character_),
            Color_Channel=dplyr::case_when(
              Infinium_Design==2 ~ 'Both',
              col=="R" ~ 'Red',
              col=='G' ~ 'Grn',
              TRUE ~ NA_character_),
            Next_Base=dplyr::case_when(
              col=="R" ~ "A",
              col=="G" ~ "C",
              TRUE ~ NA_character_
            ),
            Probe_Source=runName
          )
      } else if (!is.null(nxb_key)) {
        nxb_sym <- rlang::sym(nxb_key)
        
        ret_tib <- ret_tib %>%
          dplyr::mutate(
            Color_Channel=dplyr::case_when(
              Infinium_Design==2 ~ 'Both',
              Infinium_Design==1 & (!!nxb_sym=='A' | !!nxb_sym=='T' |
                                      !!nxb_sym=='a' | !!nxb_sym=='t') ~ 'Red',
              Infinium_Design==1 & (!!nxb_sym=='C' | !!nxb_sym=='G' |
                                      !!nxb_sym=='c' | !!nxb_sym=='g') ~ 'Grn',
              TRUE ~ NA_character_),
            col=dplyr::case_when(
              Infinium_Design==2 ~ NA_character_,
              Infinium_Design==1 ~ stringr::str_sub(Color_Channel, 1,1),
              TRUE ~ NA_character_),
            Next_Base=dplyr::case_when(
              col=="R" ~ !!nxb_sym,
              col=="G" ~ !!nxb_sym,
              TRUE ~ NA_character_
            ),
            # Next_Base=!!nxb_sym,
            Probe_Source=runName
          )
      } else {
        # Do nothing for now, but this should not be a standard use case
        return(NULL)
      }
      
      if (!is.null(pid_key)) {
        pid_key_sym = rlang::sym(pid_key)
        ret_tib <- ret_tib %>% dplyr::arrange(pid_key)
      }
      
      if (!is.null(pid_key) && !is.null(rep_key) && !is.null(rep_val)) {
        rep_key_sym = rlang::sym(rep_key)
        rep_val_sym = rlang::sym(rep_val)
        
        if (verbose>=vt) 
          cat(glue::glue("{mssg} Adding Probe Replicate({rep_key}/{rep_val})...{RET}"))
        ret_tib <- ret_tib %>% dplyr::group_by(!!rep_key_sym) %>%
          dplyr::mutate(!!rep_val_sym := dplyr::row_number(),
                        !!pid_key_sym := paste0(!!rep_key_sym,!!rep_val_sym)) %>%
          dplyr::ungroup()
      }
      
      if (!is.null(csv)) {
        safe_write(ret_tib,"csv",csv,
                   funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      }
      
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="InfI+II")
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    ret_tib
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          BSMAP Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# This isn't ready yet, but it is another way to validate CG#'s via coordinate
#   alignment with: 
#   BED File: /Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-input/GRCh37.cgn.bed.gz
#
if (FALSE) {
  
  bsmap_to_bed = function(tib,bed,
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='bsmap_to_bed') {
    
    tabs <- paste0(rep(TAB, tc), collapse='')
    mssg <- glue::glue("[{funcTag}]:{tabs}")
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- base::system.time({
      
      
      
      safe_write(bed_tib,"tsv",bed,funcTag=funcTag,
                 vt=vt+1,tc=tc+1,tt=tt)
      
      
      # ret_cnt <- ret_tib %>% base::nrow()
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    ret_tib
  }
}

# cgn_mapping_workflow(dir=run$imp_prb_dir, pattern_u="-probe_U49_cgn-table.tsv.gz", pattern_m="-probe_M49_cgn-table.tsv.gz")
cgn_mapping_workflow = function(tib, dir,
                                pattern_u,
                                pattern_m,
                                
                                prb_key = "Prb_Seq", 
                                add_key = "Address", 
                                des_key = "Prb_Des", 
                                din_key = "Prb_Din",
                                unq_key = "Unq_Key",
                                
                                out    = NULL,
                                prefix = NULL,
                                suffix = NULL,
                                
                                idxA=1,idxB=1,
                                reload=FALSE,
                                parallel=FALSE,
                                
                                verbose=0, vt=3,tc=1,tt=NULL,
                                funcTag='cgn_mapping_workflow') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (verbose>=vt) {
    cat(glue::glue("[{funcTag}]:{tabsStr} ref_u49={ref_u49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} can_u49={can_u49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} out_u49={out_u49}.{RET}"))
    cat("\n")
    cat(glue::glue("[{funcTag}]:{tabsStr} ref_m49={ref_m49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} can_m49={can_m49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} out_m49={out_m49}.{RET}"))
    cat("\n")
    cat(glue::glue("[{funcTag}]:{tabsStr}    idxA={idxA}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}    idxB={idxB}.{RET}"))
    cat("\n")
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    prb_sym <- rlang::sym(prb_key)
    des_sym <- rlang::sym(des_key)
    din_sym <- rlang::sym(din_key)
    unq_sym <- rlang::sym(unq_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #        Find All Pre-Built Reference Prefix-Partition Files:: U49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    cgn_U49_tsvs <- get_file_list(dir = dir, 
                                  pattern = pattern_u, 
                                  trim    = pattern_u,
                                  verbose = opt$verbose)
    u49_tsv_cnt <- cgn_U49_tsvs %>% names() %>% length()
    u49_pre_len <- cgn_U49_tsvs %>% names() %>% stringr::str_length() %>% max()
    
    if (verbose>=vt+2) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Found {u49_tsv_cnt} U49-files, ",
                     "prefix length={u49_pre_len}.{RET}"))
      cgn_U49_tsvs %>% head(n=3) %>% print()
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #        Find All Pre-Built Reference Prefix-Partition Files:: M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    cgn_M49_tsvs <- get_file_list(dir = dir,
                                  pattern = pattern_m, 
                                  trim    = pattern_m,
                                  verbose = opt$verbose)
    m49_tsv_cnt <- cgn_M49_tsvs %>% names() %>% length()
    m49_pre_len <- cgn_M49_tsvs %>% names() %>% stringr::str_length() %>% max()
    
    if (verbose>=vt+2) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Found {m49_tsv_cnt} M49-files, ",
                     "prefix length={m49_pre_len}.{RET}"))
      cgn_M49_tsvs %>% head(n=3) %>% print()
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #               Sanity Check:: Prefix Lengths Must be Equal!!!
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (u49_pre_len != m49_pre_len) {
      stop(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: prefix lengths for ",
                      "U49 and M49 are not equal: ",
                      "{u49_pre_len} != {m49_pre_len}!{RET2}"))
      return(ret_tib)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 Build Candidate Sub-String Probes:: U49 & M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # NOTE: Build U49/M49 Sub-sequences and split by 5' two nucelotide prefix  
    #  into sorted output files. Splitting by prefix makes the join later much 
    #  faster...
    #
    
    ret_tib <- tib %>%
      dplyr::mutate(
        Aln_Prb=deMs(!!prb_sym, uc=TRUE),
        # Aln_Rev=revCmp(Aln_Prb),
        Aln_P49=dplyr::case_when(
          !!des_sym == '2' ~ stringr::str_sub(Aln_Prb, 2),
          !!des_sym == 'U' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          !!des_sym == 'M' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          TRUE ~ NA_character_
        ),
        Pre_Nuc=stringr::str_sub(Aln_P49, 1,u49_pre_len)
      ) %>%
      dplyr::distinct(Aln_Key,Aln_Prb, .keep_all=TRUE) %>%
      clean_tibble()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Intersecting by Probe Sequence:: U49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    u49_tibs <- ret_tib %>% 
      dplyr::filter(!!des_sym == '2' | !!des_sym=='U') %>%
      dplyr::filter(!is.na(Aln_P49)) %>%
      dplyr::select(Aln_P49,Aln_Key) %>%
      dplyr::arrange(Aln_P49) %>%
      split(f=.$Pre_Nuc)
    
    u49_tib <- NULL
    if (parallel) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Intersecting probe sequences ",
                       "U49 (Parallel)...{RET}"))
      
      u49_tib <- foreach (pre_nuc=u49_tibs, .combine=rbind) %dopar% {
        intersect_seq_strand(
          can = u49_tibs[[pre_nuc]], ref = cgn_U49_tsvs[[pre_nuc]],
          bsc_str = "U49", pre_nuc = pre_nuc, idxA = idxA, idxB = idxB,
          out = out, prefix = prefix, suffix = suffix, reload = reload,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      }
    } else {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Intersecting probe sequences ",
                       "U49 (Linear)...{RET}"))
      
      for (pre_nuc in names(u49_tibs)) {
        cur_u49_tib <- intersect_seq_strand(
          can = u49_tibs[[pre_nuc]], ref = cgn_U49_tsvs[[pre_nuc]],
          bsc_str = "U49", pre_nuc = pre_nuc, idxA = idxA, idxB = idxB,
          out = out, prefix = prefix, suffix = suffix, reload = reload,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        u49_tib <- dplyr::bind_rows(u49_tib, cur_u49_tib)
        
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}{TAB} Done. Intersecting probe ",
                         "sequences U49 nuc={pre_nuc}{RET2}"))
      }
    }
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Intersecting probe sequences ",
                     "U49!{RET2}"))
    
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Intersecting by Probe Sequence:: M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    m49_tibs <- ret_tib %>% 
      dplyr::filter(!!des_sym=='M') %>%
      dplyr::filter(!is.na(Aln_P49)) %>%
      dplyr::select(Aln_P49,Aln_Key) %>%
      dplyr::arrange(Aln_P49) %>%
      split(f=.$Pre_Nuc)
    
    m49_tib <- NULL
    if (parallel) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Intersecting probe sequences ",
                       "M49 (Parallel)...{RET}"))
      
      m49_tib <- foreach (pre_nuc=m49_tibs, .combine=rbind) %dopar% {
        intersect_seq_strand(
          can = m49_tibs[[pre_nuc]], ref = cgn_M49_tsvs[[pre_nuc]],
          bsc_str = "M49", pre_nuc = pre_nuc, idxA = idxA, idxB = idxB,
          out = out, prefix = prefix, suffix = suffix, reload = reload,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      }
    } else {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Intersecting probe sequences ",
                       "M49 (Linear)...{RET}"))
      
      for (pre_nuc in names(m49_tibs)) {
        cur_m49_tib <- 
          intersect_seq_strand(
            can = m49_tibs[[pre_nuc]], ref = cgn_M49_tsvs[[pre_nuc]],
            bsc_str = "M49", pre_nuc = pre_nuc, idxA = idxA, idxB = idxB,
            out = out, prefix = prefix, suffix = suffix, reload = reload,
            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        m49_tib <- dplyr::bind_rows(m49_tib, cur_m49_tib)
        
        if (verbose>=vt) 
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}{TAB} Done. Intersecting probe ",
                         "sequences M49 nuc={pre_nuc}{RET2}"))
      }
    }
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Intersecting probe sequences ",
                     "M49!{RET2}"))
    
    
    # CODE IS A FUNCTION ABOVE NOW:: DELET THIS::
    #
    # ref_tsv <- cgn_U49_tsvs[[pre_nuc]]
    # can_tsv <- 
    #   file.path(out, paste(prefix,suffix,pre_nuc,"U49.tsv.gz", sep='.'))
    # out_tsv <- 
    #   file.path(out, paste(prefix,suffix,pre_nuc,"U49.intersect.tsv.gz", sep='.'))
    # 
    # safe_write(x = u49_tibs[[pre_nuc]], file = can_tsv, cols = FALSE, 
    #            funcTag = funcTag, verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    # 
    # cur_u49_tib <- intersect_seq(ref=ref_tsv,
    #                              can=can_tsv,
    #                              out=out_tsv,
    #                              
    #                              idxA=idxA,
    #                              idxB=idxB,
    #                              reload=reload,
    #                              
    #                              verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    # 
    # u49_key <- glue::glue("u49-tib({funcTag})")
    # u49_cnt <- print_tib(cur_u49_tib, funcTag, verbose,vt+4,tc, n=u49_key)
    #
    # ref_tsv <- cgn_M49_tsvs[[pre_nuc]]
    # can_tsv <- 
    #   file.path(out, paste(prefix,suffix,pre_nuc,"M49.tsv.gz", sep='.'))
    # out_tsv <- 
    #   file.path(out, paste(prefix,suffix,pre_nuc,"M49.intersect.tsv.gz", sep='.'))
    # 
    # safe_write(x = m49_tibs[[pre_nuc]], file = can_tsv, cols = FALSE, 
    #            funcTag = funcTag, verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    # 
    # cur_m49_tib <- intersect_seq(ref=ref_tsv,
    #                          can=can_tsv,
    #                          out=out_tsv,
    #                          
    #                          idxA=idxA,
    #                          idxB=idxB,
    #                          reload=reload,
    #                          
    #                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    # 
    # m49_key <- glue::glue("m49-tib({funcTag})")
    # m49_cnt <- print_tib(cur_m49_tib, funcTag, verbose,vt+4,tc, n=m49_key)
    
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 Join Probe Sequence Intersection:: U49/M49
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # build_file_dir(out_u49)
    # build_file_dir(out_m49)
    # 
    # u49_tib <- intersect_seq(ref=ref_u49,can=can_u49,out=out_u49,
    #                          idxA=idxA,idxB=idxB, reload=reload,
    #                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    # u49_key <- glue::glue("u49-tib({funcTag})")
    # u49_cnt <- print_tib(u49_tib,funcTag, verbose,vt+4,tc, n=u49_key)
    # 
    # m49_tib <- intersect_seq(ref=ref_m49,can=can_m49,out=out_m49,
    #                          idxA=idxA,idxB=idxB, reload=reload,
    #                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    # m49_key <- glue::glue("m49-tib({funcTag})")
    # m49_cnt <- print_tib(m49_tib,funcTag, verbose,vt+4,tc, n=m49_key)
    
    out_csv <- NULL
    ret_tib <- join_seq_intersect(u49 = u49_tib, 
                                  m49 = m49_tib, 
                                  bed = bed, org = org,
                                  verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    out_csv <- 
      file.path(out, paste(prefix,suffix,"intersect.tsv.gz", sep='.'))
    
    out_cnt <- safe_write(ret_tib, file=out_csv, funcTag=funcTag,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

cgn_mapping_workflow_single = function(ref_u49,can_u49,out_u49,
                                       ref_m49,can_m49,out_m49,
                                       ord=NULL,bed=NULL,org=NULL,out=NULL,
                                       idxA=1,idxB=1,reload=FALSE,
                                       verbose=0,vt=3,tc=1,tt=NULL,
                                       funcTag='cgn_mapping_workflow_single') {
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (verbose>=vt) {
    cat(glue::glue("[{funcTag}]:{tabsStr} ref_u49={ref_u49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} can_u49={can_u49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} out_u49={out_u49}.{RET}"))
    cat("\n")
    cat(glue::glue("[{funcTag}]:{tabsStr} ref_m49={ref_m49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} can_m49={can_m49}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} out_m49={out_m49}.{RET}"))
    cat("\n")
    cat(glue::glue("[{funcTag}]:{tabsStr}    idxA={idxA}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}    idxB={idxB}.{RET}"))
    cat("\n")
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    build_file_dir(out_u49)
    build_file_dir(out_m49)
    
    u49_tib <- intersect_seq(ref=ref_u49,can=can_u49,out=out_u49,
                             idxA=idxA,idxB=idxB, reload=reload,
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    u49_key <- glue::glue("u49-tib({funcTag})")
    u49_cnt <- print_tib(u49_tib,funcTag, verbose,vt+4,tc, n=u49_key)
    
    m49_tib <- intersect_seq(ref=ref_m49,can=can_m49,out=out_m49,
                             idxA=idxA,idxB=idxB, reload=reload,
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    m49_key <- glue::glue("m49-tib({funcTag})")
    m49_cnt <- print_tib(m49_tib,funcTag, verbose,vt+4,tc, n=m49_key)
    
    ret_tib <- join_seq_intersect(u49=u49_tib, m49=m49_tib, 
                                  bed=bed, org=org,
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    ret_key <- glue::glue("ret-tib-0({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    if (!is.null(ord) && file.exists(ord)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Loading order/canonical CSV={ord}.{RET}"))
      
      ord_tib <- suppressMessages(suppressWarnings( readr::read_csv(ord) )) %>%
        purrr::set_names(c("Can_Cgn","Can_Top","Can_Src")) %>%
        dplyr::mutate(Can_Cgn=as.integer(Can_Cgn), 
                      Can_Scr=as.integer(1),
                      Can_Scr=tidyr::replace_na(Can_Scr, 0)) %>%
        clean_tibble()
      ord_cnt <- print_tib(ord_tib,funcTag, verbose,vt+4,tc, n="ord_tib")
      
      ret_tib <- ret_tib %>% 
        dplyr::left_join(ord_tib, by=c("Imp_Cgn"="Can_Cgn")) %>%
        dplyr::mutate(Can_Scr=tidyr::replace_na(Can_Scr, 0)) %>%
        clean_tibble()
      
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Done order/canonical.{RET}{RET}"))
    }
    
    if (!is.null(out)) safe_write(ret_tib,file=out, funcTag=funcTag,
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Common Conversion Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fas_to_seq = function(tib, 
                      
                      nrec=0, 
                      gen_bld="na", 
                      gen_ref_fas, 
                      bsc_ref_fas=NULL,
                      gen_snp_fas=NULL, 
                      bsc_snp_fas=NULL,
                      
                      imp_tsv=NULL,
                      seq_csv=NULL,
                      
                      build=c(""),
                      
                      s_dat_key,
                      r_dat_key,
                      
                      ids_key="Aln_Key_Unq",
                      din_key="Ord_Din",
                      tar_din="rs",
                      
                      ext_seq="Ext_Forward_Seq",
                      iup_seq="Iupac_Forward_Sequence",
                      imp_seq="Forward_Sequence",
                      
                      srd_str="F",
                      pos_key="Coordinate",
                      chr_key="Chromosome",
                      
                      srsplit=FALSE,
                      srd_key=NULL,
                      cosplit=FALSE,
                      cos_key=NULL,
                      
                      ref_col="Ref",
                      alt_col="Alt",
                      iup_col="Iupac",
                      
                      ups_len=60,
                      seq_len=122,
                      iupac=NULL,
                      del="_",
                      
                      subset=FALSE,
                      sub_cols=NULL,
                      
                      reload=FALSE,
                      retData=FALSE,
                      
                      parallel=FALSE,
                      r_improbe=FALSE,
                      s_improbe=FALSE,
                      add_flanks=FALSE, 
                      add_matseq=TRUE, 
                      
                      verbose=0,vt=4,tc=1,tt=NULL,
                      funcTag='fas_to_seq') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Genome Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}          nrec={nrec}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       gen_bld={gen_bld}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   gen_ref_fas={gen_ref_fas}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   bsc_ref_fas={bsc_ref_fas}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   gen_snp_fas={gen_snp_fas}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   bsc_snp_fas={bsc_snp_fas}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Output File Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       imp_tsv={imp_tsv}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       seq_csv={seq_csv}{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Field Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       ids_key={ids_key}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       din_key={din_key}{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       tar_din={tar_din}{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       ext_seq={ext_seq}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       iup_seq={iup_seq}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       imp_seq={imp_seq}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       srd_str={srd_str}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       pos_key={pos_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       chr_key={chr_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       srsplit={srsplit}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       srd_key={srd_key}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       cosplit={cosplit}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       cos_key={cos_key}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       ref_col={ref_col}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       alt_col={alt_col}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       iup_col={iup_col}.{RET}"))
    cat(glue::glue("{RET}"))
    
    cat(glue::glue("[{funcTag}]:{tabsStr} Run Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}           del={del}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}         iupac={iupac}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       ups_len={ups_len}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}       seq_len={seq_len}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}      subset={subset}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}    sub_cols={sub_cols}.{RET}"))
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}      reload={reload}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}     retData={retData}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}    parallel={parallel}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   r_improbe={r_improbe}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   s_improbe={s_improbe}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  add_flanks={add_flanks}.{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}  add_matseq={add_matseq}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  etime   <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  # Reload if all data is present::
  #
  if (reload &&
      !purrr::is_null(imp_tsv) && file.exists(imp_tsv) &&
      !purrr::is_null(seq_csv) && file.exists(seq_csv)) {
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Reloading seq_csv{seq_csv}...{RET}"))
    
    stime <- base::system.time({
      ret_tib <- safe_read(seq_csv)
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    return(ret_tib)
  }
  
  stime <- base::system.time({
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Building fresh...{RET}"))
    
    # Define symbolic variables::
    #
    ids_sym  <- rlang::sym(ids_key)
    ext_sym  <- rlang::sym(ext_seq)
    iup_sym  <- rlang::sym(iup_seq)
    imp_sym  <- rlang::sym(imp_seq)
    
    pos_sym  <- rlang::sym(pos_key)
    chr_sym  <- rlang::sym(chr_key)
    
    ref_col_sym  <- rlang::sym(ref_col)
    alt_col_sym  <- rlang::sym(alt_col)
    iup_col_sym  <- rlang::sym(iup_col)
    
    # Load Genome::
    #
    dna_dat <- load_genome(file=gen_ref_fas, nrec=nrec,
                           chr_key=chr_key, ret_map=TRUE,
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # Split Data by Chromosome::
    #
    chr_list <- tib %>% 
      dplyr::arrange(!!chr_sym,!!pos_sym) %>%
      split(.[[chr_key]])
    
    chr_maps <- dna_dat$maps %>%
      split(.[[chr_key]])
    
    # Process each chromosome::
    #  TBD:: Add parallel computing over chromosomes
    #
    chr_vec_1 <- names(chr_list)
    chr_vec_2 <- names(chr_maps)
    chr_names <- chr_vec_1[chr_vec_1 %in% chr_vec_2]
    
    if (parallel) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Extacting sequence templates ",
                       "from genome (Parallel)...{RET}"))
      
      ret_tib <- foreach (chr_str=chr_names, .combine=rbind) %dopar% {
        chr_idx <- chr_maps[[chr_str]] %>% head(n=1) %>% pull(Idx) %>% as.integer()
        cur_tib <- s_improbe_template_workflow(
          tib=chr_list[[chr_str]], seq=dna_dat$seqs[[chr_idx]],
          srd_str=srd_str, pos_key=pos_key, chr_key=chr_key, chr_str=chr_str,
          ext_seq=ext_seq, iup_seq=iup_seq, imp_seq=imp_seq,
          ups_len=ups_len, seq_len=seq_len, iupac=iupac, del=del,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      }
      
    } else {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Extacting sequence templates ",
                       "from genome (Linear)...{RET}"))
      
      for (chr_str in chr_names) {
        chr_idx <- chr_maps[[chr_str]] %>% head(n=1) %>% pull(Idx) %>% as.integer()
        cur_tib <- s_improbe_template_workflow(
          tib=chr_list[[chr_str]], seq=dna_dat$seqs[[chr_idx]], 
          srd_str=srd_str, pos_key=pos_key, chr_key=chr_key, chr_str=chr_str,
          ext_seq=ext_seq, iup_seq=iup_seq, imp_seq=imp_seq,
          ups_len=ups_len, seq_len=seq_len, iupac=iupac, del=del,
          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                          Merge Probe Designs::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
        cur_key <- glue::glue("cur-tib({funcTag}-{chr_str})")
        cur_cnt <- print_tib(cur_tib %>% dplyr::select(!!ids_sym,!!imp_sym),
                             funcTag, verbose,vt=vt+4,tc=tc+1, n=cur_key)
        
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Substring ",
                         "chr_str={chr_str}.{RET}{RET}"))
      }
    }
    cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Extacting sequence templates ",
                   "from genome.{RET}{RET}"))
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc, n=ret_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #      Forward Template Sequence Generation:: Tri-fecta s-improbe
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (add_flanks) {
      
      tri_tib <- NULL
      tri_tib <- s_improbe_trifecta(tib=ret_tib,
                                    
                                    tar_din=tar_din, 
                                    ids_key=ids_key, 
                                    din_key=din_key,
                                    
                                    pos_key=pos_key,
                                    # chr_str=chr_str, 
                                    
                                    ext_seq=ext_seq, 
                                    iup_seq=iup_seq, 
                                    imp_seq=imp_seq, 
                                    
                                    ref_col=ref_col,
                                    alt_col=alt_col,
                                    iup_col=iup_col,
                                    
                                    verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
      
      ret_tib <- dplyr::bind_rows(ret_tib, tri_tib)
    }
    ret_tib <- ret_tib %>% dplyr::arrange(!!ids_sym)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Probe Design:: s_improbe()
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (s_improbe && length(build)>0) {
      
      ret_tib <- s_improbe(tib=ret_tib, 
                           build=build, 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Probe Design:: r_improbe()
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (r_improbe) {
      
      r_imp_tib <- r_improbe(tib=ret_tib,
                             
                             sr_str='FR', 
                             co_str='CO',
                             
                             ids_key=ids_key,
                             seq_key=iup_seq,
                             din_key=din_key,
                             
                             srsplit=srsplit,
                             srd_key=srd_key,
                             cosplit=cosplit,
                             cos_key=cos_key,
                             
                             prb_len=ups_len,
                             seq_len=seq_len,
                             
                             subset=subset,
                             sub_cols=sub_cols,
                             
                             parallel=parallel, 
                             add_matseq=add_matseq,
                             
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      # Join s_imp_tib (ret_tib) and r_imp_tib::
      # by_cols <- c(ids_key,din_key,iup_seq, srd_key,cos_key)
      # print(by_cols)
      # 
      # dplyr::rename(r_imp_tib, srd_key="Strand_FR", cos_key="Strand_CO") %>% 
      #   head() %>% as.data.frame() %>% print()
      # 
      # ret_tib %>% head() %>% as.data.frame() %>%  print()
      # 
      # ret_tib <- ret_tib %>% dplyr::inner_join(
      #   dplyr::rename(r_imp_tib, srd_key="Strand_FR", cos_key="Strand_CO"),
      #   by=by_cols, suffix=c("_simp","_rimp"))
      # 
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Write improbe Design Input File:: TSV
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (!is.null(imp_tsv)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing improbe input TSV={imp_tsv}...{RET}"))
      
      imp_col <- c("Seq_ID","Sequence","Genome_Build",
                   "Chromosome","Coordinate","CpG_Island")
      imp_tib <- ret_tib %>% 
        dplyr::mutate(Genome_Build=!!gen_bld, 
                      CpG_Island="FALSE") %>%
        dplyr::select(dplyr::all_of(
          c(!!ids_key,  !!imp_seq, "Genome_Build", 
            !!chr_key, !!pos_key, "CpG_Island")) ) %>%
        purrr::set_names(imp_col)
      
      out_cnt <- safe_write(x=imp_tib, file=imp_tsv, funcTag=funcTag, 
                            verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Calculate Data Summary:: CSV
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Calculating Summary...{RET}"))
    
    sum_tib <- ret_tib %>% 
      dplyr::group_by(up61,dn61,Des_Din) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    sum_key <- glue::glue("ret-summary({funcTag})")
    sum_cnt <- print_tib(sum_tib,funcTag, verbose,vt+4,tc, n=sum_key)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Write Full Data Set:: CSV
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (!is.null(seq_csv)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing seq table CSV={seq_csv}...{RET}"))
      
      sum_csv <- seq_csv %>% 
        stringr::str_remove(".gz$") %>%
        stringr::str_remove(".[tcsv]+$") %>%
        paste(".summary.csv.gz")
      
      seq_cnt <- safe_write(x=ret_tib, file=seq_csv, funcTag=funcTag, 
                            verbose=verbose,vt=vt,tc=tc,append=FALSE)
      sum_cnt <- safe_write(x=sum_tib, file=sum_csv, funcTag=funcTag, 
                            verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    
    if (retData) ret_dat[[s_dat_key]] <- ret_tib
    if (retData) ret_dat[[r_dat_key]] <- r_imp_tib
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  if (retData) return(ret_dat)
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Run Both improbe and r-improbe Simultaneously!
#
#  NOTE:: This code hasn't been tested much and has more or less been 
#   re-implemented in seperate workflows, but keeping it around for now...
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  improbe_design_all = function(tib, ptype, outDir, gen, image, shell,
                                seqKey="IUPAC_Sequence", sr_str="FR",
                                reduce_imp=TRUE, parse_din=FALSE,
                                parallel=TRUE, retData=FALSE,
                                sidx=2, plen=50,del='_',
                                verbose=0,vt=3,tc=1,tt=NULL,
                                funcTag='improbe_design_all') {
    
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Starting; ptype={ptype}...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    ret_dat <- NULL
    stime <- system.time({
      idx1 <- sidx
      len1 <- plen - 1
      idx2 <- sidx + 1
      len2 <- plen
      
      if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
      
      imp_fwd_name <- paste(ptype,gen,"improbe_fwd-seq.tsv.gz", sep=del)
      imp_fwd_path <- file.path(outDir, imp_fwd_name)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                                 improbe::
      #                                docker/c++
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      imp_fwd_tib <- tib %>% 
        dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island)
      readr::write_tsv(imp_fwd_tib, imp_fwd_path)
      
      imp_des_tib <- run_improbe_docker(
        dir=outDir, file=imp_fwd_name, 
        name=gen, image=image, shell=shell, 
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        dplyr::mutate(
          Seq_48U_1=stringr::str_sub(UnMethyl_Probe_Sequence, idx1,len1) %>% 
            stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% 
            stringr::str_replace_all('Y', 'T'),
          
          Seq_48U_2=stringr::str_sub(UnMethyl_Probe_Sequence, idx2,len2) %>% 
            stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% 
            stringr::str_replace_all('Y', 'T'),
          
          Chromosome=as.character(Chromosome),
          Strand_SR=Methyl_Allele_FR_Strand,
          Strand_TB=stringr::str_sub(Methyl_Allele_TB_Strand,1,1),
          Strand_CO=Methyl_Allele_CO_Strand,
          Probe_Score_Min=pmin(Methyl_Final_Score,UnMethyl_Final_Score),
          Underlying_CpG_Count=as.integer(Methyl_Underlying_CpG_Count),
          Underlying_CpG_Min_Dist=as.integer(Methyl_Underlying_CpG_Min_Dist)
        )
      
      if (parse_din) imp_des_tib <- imp_des_tib %>%
        tidyr::separate(Seq_ID, into=c("Seq_ID","DiNuc"), sep=del)
      
      #
      # TBD:: Add better reduced returned value...
      #
      if (reduce_imp) {
        imp_des_tib <- imp_des_tib %>%
          dplyr::select(
            Seq_ID,Forward_Sequence,Genome_Build,Chromosome,Coordinate,Top_Sequence,
            Strand_SR,Strand_TB,Strand_CO,Probe_Score_Min,
            Underlying_CpG_Count,Underlying_CpG_Min_Dist,
            Methyl_Probe_Sequence,UnMethyl_Probe_Sequence,Seq_48U_1,Seq_48U_2)
      }
      imp_key <- glue::glue("imp_des_tib({funcTag})")
      imp_cnt <- print_tib(imp_des_tib,funcTag, verbose,vt+4,tc, n=imp_key)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                               r-improbe::
      #                           de-novo IUPAC designs
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      iup_des_tib <- r_improbe(
        tib=tib,
        ids_key="Seq_ID",seqKey=seqKey,prbKey="Probe_Type",
        sr_str=sr_str, parallel=parallel,
        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
        dplyr::mutate(
          
          Seq_48U_1=stringr::str_sub(PRB1_U_MAT, idx1,len1) %>% 
            stringr::str_to_upper() %>% 
            stringr::str_replace_all('R', 'A') %>% 
            stringr::str_replace_all('Y', 'T'),
          
          Seq_48U_2=stringr::str_sub(PRB1_U_MAT, idx2,len2) %>% 
            stringr::str_to_upper() %>% 
            stringr::str_replace_all('R', 'A') %>% 
            stringr::str_replace_all('Y', 'T'),
        )
      
      if (parse_din) iup_des_tib <- iup_des_tib %>%
        tidyr::separate(Seq_ID, into=c("Seq_ID","DiNuc"), sep=del)
      
      if (!reduce_imp) imp_des_tib <- imp_des_tib %>% 
        dplyr::rename(Probe_Type_IMP=Probe_Type)
      
      ret_tib <- dplyr::inner_join(
        imp_des_tib, iup_des_tib,
        by=c("Seq_ID", "Strand_SR", "Strand_CO",
             "Seq_48U_1", "Seq_48U_2"),
        suffix=c("_IMP", "_IUP")
      )
      
      if (retData) ret_dat$imp <- imp_des_tib
      if (retData) ret_dat$iup <- iup_des_tib
      
      ret_key <- glue::glue("ret-FIN({funcTag})")
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    if (retData) return(ret_dat)
    
    ret_tib
  }
  
}

if (FALSE) {
  if (FALSE) {
    cur_tib <- chr_list[[chr_str]]
    chr_idx <- chr_maps[[chr_str]] %>% 
      dplyr::filter(!!chr_sym==chr_str) %>% 
      head(n=1) %>% pull(Idx) %>% as.integer()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 Parse Forward Template Sequence from Genome::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ref_seqs <- parse_genomic_seqs(tib=cur_tib,
                                   seq=as.character(dna_dat$seqs[[chr_idx]]), 
                                   
                                   srd_str=srd_str,
                                   pos_key=pos_key,
                                   chr_str=chr_str,
                                   chr_idx=chr_idx,
                                   
                                   ups_len=ups_len, 
                                   seq_len=seq_len,
                                   verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                Forward Template Sequence Generation:: s-improbe
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    cur_tib <- s_improbe_template(tib=cur_tib,
                                  seqs=ref_seqs,
                                  
                                  chr_str=chr_str, 
                                  chr_idx=chr_idx,
                                  
                                  ext_seq=ext_seq,
                                  iup_seq=iup_seq,
                                  imp_seq=imp_seq,
                                  
                                  iupac=iupac, 
                                  ups_len=ups_len,
                                  seq_len=seq_len,
                                  verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
  }
  
}
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Modern improbe IO Methods:: 
#                          Not so Modern Any More
#
#              These functions should probably be deprecated
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Common Conversion Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# NOTE:: Pretty sure these function is old and not needed. Its probably not
#   needed anymore and can be moved to the graveyard. Seems useful though...
#
if (FALSE) {
  
  seq_to_prbs = function(tib, seq, ids,
                         cgn="Seq_ID",chr="Chromosome", pos="Coordinate",
                         srd="F",
                         ups_len=60, seq_len=122,
                         iupac=NULL, add_flank=FALSE,
                         del="_", 
                         verbose=0,vt=3,tc=1,tt=NULL) {
    funcTag <- 'seq_to_prbs'
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    
    stop(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: This function isn't ready!!!{RET}{RET}"))
    return(ret_tib)
    
    stime <- base::system.time({
      
      cgn_sym <- rlang::sym(cgn)
      chr_sym <- rlang::sym(chr)
      pos_sym <- rlang::sym(pos)
      ids_sym <- rlang::sym(ids)
      
      cur_begs <- NULL
      cur_ends <- NULL
      cur_begs <- tib %>% dplyr::pull(!!pos_sym) - ups_len - 2
      cur_ends <- cur_begs + seq_len - 1 + 4
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_begs={RET}"))
        cur_begs %>% head() %>% print()
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_ends={RET}"))
        cur_ends %>% head() %>% print()
      }
      
      # Sub-string ref sequences::
      #
      ref_seqs <- NULL
      ref_seqs <- 
        stringr::str_sub( as.character(dna_seqs[[chr_idx]]), cur_begs, cur_ends) %>%
        stringr::str_to_upper()
      
      if (srd=="R") ref_seqs <- ref_seqs %>% cmpl()
      
      # Build individual parts of the template sequence::
      # up01.up02...up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn02.dn01
      #                             [  C   G  ]
      #
      ref_up01 <- stringr::str_sub(ref_seqs, 1,1)
      ref_up02 <- stringr::str_sub(ref_seqs, 2,2)
      ref_up58 <- stringr::str_sub(ref_seqs, 3,ups_len-2+2)
      
      ref_up59 <- stringr::str_sub(ref_seqs, ups_len-2+3,ups_len-2+3)
      ref_up60 <- stringr::str_sub(ref_seqs, ups_len-2+4,ups_len-2+4)
      ref_up61 <- stringr::str_sub(ref_seqs, ups_len-2+5,ups_len-2+5)
      
      iupac_vec <- ref_up61
      if (!is.null(iupac)) iupac_vec <- cur_tib %>% dplyr::pull(!!iupac)
      
      ref_dn61 <- stringr::str_sub(ref_seqs, ups_len-2+6,ups_len-2+6)
      ref_dn60 <- stringr::str_sub(ref_seqs, ups_len-2+7,ups_len-2+7)
      ref_dn59 <- stringr::str_sub(ref_seqs, ups_len-2+8,ups_len-2+8)
      
      ref_dn58 <- stringr::str_sub(ref_seqs, ups_len-2+9,ups_len-2+ups_len-2+8)
      ref_dn02 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+8,ups_len-2+ups_len-2+8)
      ref_dn01 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+9,ups_len-2+ups_len-2+9)
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ref_seqs({chr_str})={RET}"))
        ref_seqs %>% head() %>% print()
      }
      
      # Add all fields to current tibble::
      #
      cur_tib <- cur_tib %>%
        dplyr::mutate(
          !!ext_sym := ref_seqs %>% addBrac(),
          
          up01=ref_up01,
          up02=ref_up02,
          up58=ref_up58,
          
          up59=ref_up59,
          up60=ref_up60,
          up61=ref_up61,
          
          iupac=iupac_vec,
          
          dn61=ref_dn61,
          dn60=ref_dn60,
          dn59=ref_dn59,
          
          dn58=ref_dn58,
          dn02=ref_dn02,
          dn01=ref_dn01,
          
          # Now we can assemble to optimal template::
          #
          !!des_sym := paste0(up58,
                              up59,up60,
                              iupac,dn61,
                              dn60,dn59,
                              dn58) %>%
            addBrac(),
          !!imp_sym := paste0(up58,
                              up59,up60,
                              "CG",
                              dn60,dn59,
                              dn58) %>%
            addBrac(),
          
          # TBD:: Pretty sure this right, but really needs more testing...
          #
          Din_Str=stringr::str_to_lower(paste0(iupac,dn61)),
          Des_Din=dplyr::case_when(
            up61=='C' & dn61=='G' ~ "cg",
            up61=='C' & dn61!='G' ~ "ch",
            up61!='C' ~ "rs",
            TRUE ~ NA_character_
          )
        )
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_tib({chr_str})={RET}"))
        cur_tib %>% print()
      }
      
      #  
      # TBD:: We can also break out tri-fecta sites if probe type == rs
      #
      ups_tib <- NULL
      dns_tib <- NULL
      
      if (add_flank) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Upstream CG Flank Seuqnces({chr_str})...{RET}"))
        
        ups_tib <- cur_tib %>% 
          dplyr::filter(up59=="C" & up60=="G") %>%
          dplyr::mutate(
            !!des_sym := paste0(up01,up02,up58,up59,up60,iupac,dn61,dn60,dn59,dn58) %>%
              stringr::str_sub(1,seq_len) %>%
              addBrac(),
            Din_Str=stringr::str_to_lower(paste0(up59,up60)),
            Des_Din="cg",
            Name=paste(!!ids_sym,"CG_UP", sep=del),
            Coordinate=as.integer(Coordinate-2)
          )
        
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Downstream CG Flank Seuqnces({chr_str})...{RET}"))
        
        dns_tib <- cur_tib %>% 
          dplyr::filter(dn61=="C" & dn60=="G") %>%
          dplyr::mutate(
            !!des_sym := paste0(up58,up59,up60,iupac,dn61,dn60,dn59,dn58,dn02) %>%
              stringr::str_sub(2) %>%
              addBrac(),
            Din_Str=stringr::str_to_lower(paste0(dn61,dn60)),
            Des_Din="cg",
            Name=paste(!!ids_sym,"CG_DN", sep=del),
            Coordinate=as.integer(Coordinate+2)
          )
      }
      
      # Add New Tri-fecta Probes::
      #  - NOTE:: These will have their names changed later after alignment...
      #
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joining and sorting...{RET}"))
      
      cur_tib <- dplyr::bind_rows(cur_tib, ups_tib, dns_tib) %>%
        dplyr::arrange(!!ids_sym)
      
      ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
      
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Substring chr_str={chr_str}.{RET}{RET}"))
      
      
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    ret_tib
  }
  
  bed_to_prbs = function(tib, fas,
                         file=NULL,
                         din="Probe_Type",gen="na",
                         cgn="Seq_ID",chr="Chromosome", pos="Coordinate",
                         ext_seq="Fwd_Temp_Seq",
                         iup_seq="Iup_Temp_Seq",
                         imp_seq="Imp_Temp_Seq",
                         ups_len=60, seq_len=122, nrec=0,
                         iupac=NULL, add_flank=FALSE,
                         del="_", 
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='bed_to_prbs') {
    
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- base::system.time({
      
      cgn_sym <- rlang::sym(cgn)
      chr_sym <- rlang::sym(chr)
      pos_sym <- rlang::sym(pos)
      
      # Get List of BSC Genomes::
      # fas_dir <- fas %>% base::dirname()
      fas_pre <- fas %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$")
      if (is.null(gen)) gen <- base::basename(fas_pre)
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Prefix({gen})={fas_pre}.{RET}"))
      
      # NOTE:: There isn't an opposite genome, build that from the converted::
      #  This is why we don't need 'srd_cos <- c("C","O")' below::
      #
      srd_cos <- c("C")
      srd_frs <- c("F","R")
      srd_mus <- c("U","M","D")
      
      for (fr in srd_frs) {
        # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} fr={fr}...{RET}"))
        
        for (co in srd_cos) {
          # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} co={co}...{RET}"))
          
          for (mu in srd_mus) {
            if (verbose>=vt) 
              cat(glue::glue("[{funcTag}]:{tabsStr} fr={fr}, co={co}, mu={mu}...{RET}"))
            gen_fas <- paste0(fas_pre,".", fr,co,mu, ".fa.gz")
            
            if (!file.exists(gen_fas)) {
              if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Skipping; gen_fas=({fr},{co},{mu})",
                                              "={gen_fas} Does NOT Exist!{RET}"))
              next
            } else {
              if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading; gen_fas({fr},{co},{mu})",
                                              "={gen_fas}...{RET}"))
            }
            
            cur_tib <- fas_to_seq(tib=tib, fas=gen_fas, file=file,
                                  name=cgn, din=din, gen=gen,
                                  chr1=chr, pos=pos, srd=fr,
                                  ext_seq=ext_seq,iup_seq=iup_seq,imp_seq=imp_seq,
                                  
                                  # See Above:: Will Pass these variables in later...
                                  # ext_seq="Fwd_Temp_Seq",des_seq="Iup_Temp_Seq",
                                  # imp_seq="Imp_Temp_Seq",
                                  
                                  iupac=iupac, mnrec=nrec, add_flank=FALSE,
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt
            ) %>%
              dplyr::mutate(Srd_FR=fr, Srd_CO=co, Prb_Des=mu)
            
            cur_cnt <- print_tib(cur_tib,funcTag, verbose,vt+4,tc, n="cur-pre-join")
            cgn_vec <- cur_tib %>% dplyr::pull(!!cgn_sym)
            
            # TBD:: Add the current FR/CO/TB fields to cur_tib
            # TBD:: Add the opposite designs
            # TBD:: Remove previous information from the input tibble
            
            prb_tib <- NULL
            if (fr=="F") {
              prb_tib <- tibble::tibble(
                !!cgn_sym := cgn_vec,
                
                # Forward Converted Designs:: Infinium I/II
                Prb_SeqC1 = paste0(
                  cur_tib$up61,
                  cur_tib$dn61,
                  cur_tib$dn60,
                  cur_tib$dn59,
                  cur_tib$dn58 %>% stringr::str_sub(1,46)
                ) %>% revCmp(),
                Prb_SeqC2 = paste0(
                  # cur_tib$up61,
                  cur_tib$dn61,
                  cur_tib$dn60,
                  cur_tib$dn59,
                  cur_tib$dn58 %>% stringr::str_sub(1,47)
                ) %>% revCmp(),
                
                # Forward Opposite Designs:: Infinium I/II
                #   Just need to back up by one base::
                Prb_SeqO1 = paste0(
                  cur_tib$up58 %>% stringr::str_sub(12,58), # 13 -> 12
                  cur_tib$up59,
                  cur_tib$up60,
                  cur_tib$up61
                  # cur_tib$dn61 # Remove dn61
                ), # No complement
                Prb_SeqO2 = paste0(
                  cur_tib$up58 %>% stringr::str_sub(11,58), # 12 -> 11
                  cur_tib$up59,
                  cur_tib$up60
                  # cur_tib$up61 # Remove dn61
                  # cur_tib$dn61,
                  # cur_tib$dn60
                ) # No Complement
                
              )
            } else if (fr=="R") {
              #
              # NOTE: From the ch/rs perspective the reverse strand 
              #   coordinates frame should be moved back by 1
              #
              # if ()
              prb_tib <- tibble::tibble(
                # Seq_ID = cur_tib$Seq_ID,
                # !!cgn_sym := cur_tib[[!!cgn_sym]],
                !!cgn_sym := cgn_vec,
                
                # Reverse Converted Designs:: Infinium I/II
                Prb_SeqC1 = paste0(
                  cur_tib$up58 %>% stringr::str_sub(13,58),
                  cur_tib$up59,
                  cur_tib$up60,
                  cur_tib$up61,
                  cur_tib$dn61
                ) %>% cmpl(),
                Prb_SeqC2 = paste0(
                  cur_tib$up58 %>% stringr::str_sub(12,58),
                  cur_tib$up59,
                  cur_tib$up60,
                  cur_tib$up61
                  # cur_tib$dn61,
                  # cur_tib$dn60
                ) %>% cmpl(),
                
                # Reverse Opposite Designs:: Infinium I/II
                #   Just need to move forward by 1::
                Prb_SeqO1 = paste0(
                  # cur_tib$up61, # Remove up61
                  cur_tib$dn61,
                  cur_tib$dn60,
                  cur_tib$dn59,
                  cur_tib$dn58 %>% stringr::str_sub(1,47) # 46 -> 47
                ), # %>% revCmp(), # No revCmp
                Prb_SeqO2 = paste0(
                  # cur_tib$up61,
                  # cur_tib$dn61, # Remove dn61
                  cur_tib$dn60,
                  cur_tib$dn59,
                  cur_tib$dn58 %>% stringr::str_sub(1,48) # 47 -> 48
                ) # %>% revCmp(), # No revCmp
                
              )
            } else {
              stop(glue::glue("{RET}[{funcTag}]: Unsupported fr={fr}...{RET}"))
              return(NULL)
            }
            prb_cnt <- print_tib(prb_tib,funcTag, verbose,vt+4,tc, n="prb")
            
            cur_tib <- cur_tib %>% dplyr::left_join(prb_tib, by=cgn)
            ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
          }
        }
      }
      
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    ret_tib
  }
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Modern improbe IO Methods:: 
#                          Not so Modern Any More
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

addSeq48U = function(tib, field,
                     sidx=2, plen=50,
                     verbose=0,vt=3,tc=1,tt=NULL,
                     funcTag='addSeq48U') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    idx1 <- sidx
    len1 <- plen - 1
    idx2 <- sidx + 1
    len2 <- plen
    
    field_sym <- rlang::sym(field)
    
    ret_tib <- tib %>% dplyr::mutate(
      Seq_48U_1=stringr::str_sub(!!field_sym, idx1,len1) %>% 
        stringr::str_to_upper() %>% 
        stringr::str_replace_all('R', 'A') %>% 
        stringr::str_replace_all('Y', 'T'),
      
      Seq_48U_2=stringr::str_sub(!!field_sym, idx2,len2) %>% 
        stringr::str_to_upper() %>% 
        stringr::str_replace_all('R', 'A') %>% 
        stringr::str_replace_all('Y', 'T'),
    )
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

isTibIMP = function(tib,
                    verbose=0,vt=3,tc=1,tt=NULL,
                    funcTag='isTibIMP') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- FALSE
  stime <- system.time({
    
    imp_col <- INIT_IMP_HEADER()
    imp_key <- imp_col$cols %>% names()
    sam_col <- tib %>% names()
    
    if (base::all.equal(imp_key,sam_col)) ret_val <- TRUE
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_val
}

loadIMP = function(file, max=Inf, format=NULL, del='\t',
                   verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadIMP'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    imp_col <- INIT_IMP_HEADER()
    imp_key <- imp_col$cols %>% names()
    sam_col <- utils::read.table(file, header=TRUE, sep=del, nrows=1) %>% names()
    
    if (max==Inf) {
      if (base::all.equal(imp_key,sam_col)) {
        ret_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file, col_types = imp_col) ))
      } else {
        ret_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
      }
    } else {
      if (!is.numeric(max)) {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: max must be an integer; max={max}!!!{RET}{RET}"))
        return(ret_tib)
      }
      max <- as.integer(max)
      
      if (base::all.equal(imp_key,sam_col)) {
        ret_tib <- suppressMessages(suppressWarnings( 
          utils::read.table(file, header=TRUE, sep=del, nrows=max, colClasses=imp_col) )) %>% 
          tibble::as_tibble() %>% dplyr::mutate(Chromosome=as.character(Chromosome))
      } else {
        ret_tib <- suppressMessages(suppressWarnings( 
          utils::read.table(file, header=TRUE, sep=del, nrows=max) )) %>% tibble::as_tibble()
      }
    }
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

tibToFas = function(tib, key,seq, prefix,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'tibToFas'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; key={key}, seq={seq}, prefix={prefix}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    dir <- base::dirname(prefix)
    if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    
    file_fas <- paste0(prefix, ".fa.gz")
    file_csv <- paste0(prefix, ".csv.gz")
    
    key_sym <- rlang::sym(key)
    seq_sym <- rlang::sym(seq)
    
    ret_tib <- tib %>% dplyr::mutate(
      # fas_line=paste0(">",!!key, "\n",!!seq)
      fas_line=paste0(">",!!key_sym, "\n",!!seq_sym)
    )  %>%
      dplyr::filter(!is.na(fas_line)) %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    fas_vec <- ret_tib %>%
      # dplyr::distinct(fas_line) %>%
      dplyr::pull(fas_line)
    
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Writing Data CSV={file_csv}...{RET}"))
    readr::write_csv(ret_tib, file_csv)
    
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Writing Data FAS={file_fas}...{RET}"))
    readr::write_lines(x=fas_vec, file=file_fas)
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      ret_tib %>% print()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

ordToFas = function(tib, dir, name, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'ordToFas'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}, dir={dir}.{RET}"))
  
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  fas_file <- file.path(dir, paste0(name,'.cgn.fa.gz'))
  
  fas_tib <- tib %>% dplyr::mutate(
    line=dplyr::case_when(
      Normalization_Bin=='C' ~ paste0('>',AlleleA_Probe_Id,'\n',AlleleA_Probe_Sequence),
      Normalization_Bin!='C' ~ paste0('>',AlleleA_Probe_Id,'\n',AlleleA_Probe_Sequence,'\n',
                                      '>',AlleleB_Probe_Id,'\n',AlleleB_Probe_Sequence),
      TRUE ~ NA_character_ ) ) %>%
    dplyr::filter(!is.na(line)) %>% dplyr::pull(line)
  
  readr::write_lines(x=fas_tib, file=fas_file)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  fas_file
}

writeBedFas = function(tib, dir, name, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'writeBedFas'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}, dir={dir}.{RET}"))
  
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  fas_file <- file.path(dir, paste0(name,'.fa.gz'))
  bed_file <- file.path(dir, paste0(name,'.sorted.bed.gz'))
  
  tib <- tib %>% dplyr::arrange(-MinScore) %>% dplyr::distinct(PRB1_U, .keep_all=TRUE) %>%
    dplyr::arrange(Chromosome,Coordinate) %>% dplyr::rename(Chrom=Chromosome, Beg=Coordinate) %>%
    dplyr::mutate(End=Beg+stringr::str_length(PRB1_U), 
                  Full_ID=paste(Seq_ID,FR_Imp,TB_Imp,CO_Imp,Beg, sep='_') )
  
  fas_tib <- tib %>% dplyr::select(Full_ID, PRB1_U) %>% 
    dplyr::mutate(line=paste0('>',Full_ID,'\n',PRB1_U)) %>% dplyr::pull(line)
  bed_tib <- tib %>% dplyr::select(Chrom, Beg, End, Full_ID, MinScore, FR_Imp, TB_Imp, CO_Imp, everything())
  
  readr::write_lines(x=fas_tib, file=fas_file)
  readr::write_tsv(bed_tib, bed_file)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  bed_tib
}

writeImprobeInput = function(tib, name, dir, run=FALSE, 
                             exe=NULL, impTango=NULL, imp13mer=NULL, 
                             tbVar='BOTH', coVar='BOTH',
                             verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'writeImprobeInput'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}, dir={dir}.{RET}"))
  
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  
  imp_run_sh  <- file.path(dir, paste0(name,'.improbe-input.sh'))
  imp_inp_tsv <- file.path(dir, paste0(name,'.improbe-input.tsv'))
  imp_out_tsv <- file.path(dir, paste0(name,'.improbe-output.tsv'))
  imp_log_txt <- file.path(dir, paste0(name,'.improbe-output.log'))
  imp_fin_txt <- file.path(dir, paste0(name,'.improbe-output.fin.txt'))
  imp_cmd_str <- ''
  
  readr::write_tsv(tib, imp_inp_tsv)
  
  if (!is.null(exe) && ! is.null(impTango) && !is.null(imp13mer)) {
    # stopifnot(file.exists(exe))
    # stopifnot(file.exists(impTango))
    # stopifnot(file.exists(imp13mer))
    
    # ${CMD} -oASPE -tBOTH -cBoth -n${fn_13mer} -a${fn_tango} -V ${fn_in} >$fn_out 2>${fn_log}
    imp_cmd_str <- glue::glue("{exe} -oASPE -t{tbVar} -c{coVar} -n{imp13mer} -a{impTango} -V ",
                              "{imp_inp_tsv} >{imp_out_tsv} 2>{imp_log_txt}{RET}",
                              "touch {imp_fin_txt}{RET}")
    
    readr::write_file(imp_cmd_str, imp_run_sh)
    Sys.chmod(imp_run_sh, mode='0777')
    
    if (run) system(imp_run_sh, wait = TRUE)
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  imp_out_tsv
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Load/Format Input File Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
loadManifestRS = function(file, swap, revAllele=FALSE, 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(RS)={file}.{RET}"))
  
  snp_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(file) ))
  snp_swp_tib <- suppressMessages(suppressWarnings(readr::read_tsv(swap) ))
  
  # NOTE: Rediculous Swapping Needed for Infinium I M/U -> U/M from old manufacturing mistake::
  #  - This method should not be used...
  if (revAllele) {
    snp_man_tib <- snp_man_tib %>% 
      dplyr::mutate(
        TMP_ADD_A=AddressA_ID,
        TMP_ADD_B=AddressB_ID,
        TMP_SEQ_A=AlleleA_ProbeSeq,
        TMP_SEQ_B=AlleleB_ProbeSeq,
        IS_SWAP=case_when(
          Infinium_Design_Type=='I' & stringr::str_ends(AlleleA_ProbeSeq,'C') ~ TRUE,
          TRUE ~ FALSE
        ),
        AddressA_ID=case_when( IS_SWAP ~ TMP_ADD_B, TRUE ~ AddressA_ID),
        AddressB_ID=case_when( IS_SWAP ~ TMP_ADD_A, TRUE ~ AddressB_ID),
        AlleleA_ProbeSeq=case_when( IS_SWAP ~ TMP_SEQ_B, TRUE ~ AlleleA_ProbeSeq ),
        AlleleB_ProbeSeq=case_when( IS_SWAP ~ TMP_SEQ_A, TRUE ~ AlleleB_ProbeSeq )
      )
  }
  
  snp_man_tib <-  snp_man_tib %>% 
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$') ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq, 
                  Infinium_Design_Type, Next_Base, Color_Channel, CGN)
  
  snp_swp_tib <- snp_swp_tib %>% tidyr::separate(Seq_ID, into=c('IlmnID', 'diNUC'), sep='_') %>% 
    dplyr::rename(Forward_Sequence=Sequence, CHR=Chromosome, MAPINFO=Coordinate) %>% 
    dplyr::mutate(Forward_Sequence=stringr::str_replace(Forward_Sequence,'\\[CG\\]', paste0('[',diNUC,']') ) ) %>% 
    dplyr::select(-CpG_Island)
  
  snp_man_tib <- snp_man_tib %>% dplyr::inner_join(snp_swp_tib, by="IlmnID") %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  snp_man_tib
}                          

loadManifestCG = function(file, pr='cg', 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestCG'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CG)={file}.{RET}"))
  
  cpg_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(file, skip=7) )) %>%
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[([A-Za-z][A-Za-z])].*$', "\\$1"),'^\\\\') ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cpg_man_tib
}                          

loadManifestCH = function(file, pr='ch', ry='R', 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestCH'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CH)={file}.{RET}"))
  
  cph_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(cph_man_csv ))) %>% 
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(Strand=stringr::str_sub(IlmnID, -1), FR=Strand, CO='O',
                  TP=case_when( is.na(AlleleB_ProbeSeq) ~ 'II', TRUE ~ 'I' )) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=paste0(ry,stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[[A-Za-z]([A-Za-z])].*$', "\\$1"),'^\\\\')),
                  Forward_Sequence=stringr::str_replace(Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', paste0('[',diNUC,']')) ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cph_man_tib
}

loadManifestYH = function(file, pr='ch', 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestYH'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CH)={file}.{RET}"))
  
  cph_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(cph_man_csv ))) %>% 
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(Strand=stringr::str_sub(IlmnID, -1), FR=Strand, CO='O',
                  TP=case_when( is.na(AlleleB_ProbeSeq) ~ 'II', TRUE ~ 'I' )) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=paste0('Y',stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[[A-Za-z]([A-Za-z])].*$', "\\$1"),'^\\\\')),
                  Forward_Sequence=stringr::str_replace(Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', paste0('[',diNUC,']')) ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cph_man_tib
}                          

#
# TBD:: This function should probably be deprecated::
#
loadImprobeDesign = function(file=NULL, src_des_tib=NULL, 
                             verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadImprobeDesign'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading improbe={file}.{RET}"))
  
  if (!is.null(file)) src_des_tib <- suppressMessages(suppressWarnings(readr::read_tsv(file)))
  stopifnot(!is.null(src_des_tib))
  
  src_des_tib <- src_des_tib %>%
    dplyr::rename(PRB1_U=UnMethyl_Probe_Sequence,
                  PRB1_M=Methyl_Probe_Sequence,
                  NXB_U=UnMethyl_Next_Base,
                  NXB_M=Methyl_Next_Base) %>%
    dplyr::mutate(Probe_ID=paste(Seq_ID, Methyl_Allele_FR_Strand, stringr::str_sub(Methyl_Allele_TB_Strand,1,1), 
                                 Methyl_Allele_CO_Strand, sep='_'),
                  FR=case_when(Methyl_Allele_FR_Strand=='F'   ~ TRUE, Methyl_Allele_FR_Strand=='R'   ~ FALSE, TRUE ~ NA),
                  TB=case_when(Methyl_Allele_TB_Strand=='TOP' ~ TRUE, Methyl_Allele_TB_Strand=='BOT' ~ FALSE, TRUE ~ NA),
                  CO=case_when(Methyl_Allele_CO_Strand=='C'   ~ TRUE, Methyl_Allele_CO_Strand=='O'   ~ FALSE, TRUE ~ NA),
                  diNUC='CG',
                  NXB_IMP=case_when(NXB_U==NXB_M ~ NXB_U,
                                    TRUE ~ NA_character_),
                  COL_U=case_when(NXB_U=='A'|NXB_U=='T' ~ FALSE, # Red == FALSE
                                  NXB_U=='C'|NXB_U=='G' ~ TRUE,  # Grn == TRUE
                                  TRUE ~ NA),
                  COL_M=case_when(NXB_M=='A'|NXB_M=='T' ~ FALSE, # Red == FALSE
                                  NXB_M=='C'|NXB_M=='G' ~ TRUE,  # Grn == TRUE
                                  TRUE ~ NA),
                  COL_IMP=case_when(COL_U==COL_M ~ COL_U,
                                    TRUE ~ NA)
                  
    ) %>%
    
    # Design Score Parameters::
    dplyr::rename(
      PRB_SCR_U=UnMethyl_Final_Score,
      PRB_SCR_M=Methyl_Final_Score,
      PRB_SCR_S=Probeset_Score,
      
      TM_RAW_M=Methyl_Tm,
      TM_SCR_M=Methyl_Tm_Score,
      GC_RAW_M=Methyl_GC_Percent,
      GC_SCR_M=Methyl_GC_Score,
      KM_RAW_M=Methyl_13mer_Count,
      KM_SCR_M=Methyl_13mer_Score,
      AD_RAW_M=Methyl_Address_Count,
      AD_SCR_M=Methyl_Address_Score,
      CM_RAW_M=Methyl_Self_Complementarity,
      CM_SCR_M=Methyl_Self_Complementarity_Score,
      MO_RAW_M=Methyl_Mono_Run,
      MO_SCR_M=Methyl_Mono_Run_Score,
      EP_RAW_M=Methyl_Ectopic_Count,
      EP_SCR_M=Methyl_Ectopic_Score,
      CG_RAW_M=Methyl_Underlying_CpG_Count,
      MD_RAW_M=Methyl_Underlying_CpG_Min_Dist,
      CG_SCR_M=Methyl_Underlying_CpG_Score,
      NB_SCR_M=Methyl_Next_Base_Score,
      
      TM_RAW_U=UnMethyl_Tm,
      TM_SCR_U=UnMethyl_Tm_Score,
      GC_RAW_U=UnMethyl_GC_Percent,
      GC_SCR_U=UnMethyl_GC_Score,
      KM_RAW_U=UnMethyl_13mer_Count,
      KM_SCR_U=UnMethyl_13mer_Score,
      AD_RAW_U=UnMethyl_Address_Count,
      AD_SCR_U=UnMethyl_Address_Score,
      CM_RAW_U=UnMethyl_Self_Complementarity,
      CM_SCR_U=UnMethyl_Self_Complementarity_Score,
      MO_RAW_U=UnMethyl_Mono_Run,
      MO_SCR_U=UnMethyl_Mono_Run_Score,
      EP_RAW_U=UnMethyl_Ectopic_Count,
      EP_SCR_U=UnMethyl_Ectopic_Score,
      CG_RAW_U=UnMethyl_Underlying_CpG_Count,
      MD_RAW_U=UnMethyl_Underlying_CpG_Min_Dist,
      CG_SCR_U=UnMethyl_Underlying_CpG_Score,
      NB_SCR_U=UnMethyl_Next_Base_Score
    ) %>%
    dplyr::select(Probe_ID, Seq_ID, FR, TB, CO, diNUC, NXB_IMP, COL_IMP, PRB1_U, PRB1_M, 
                  Forward_Sequence, 
                  Genome_Build, Chromosome, Coordinate,
                  dplyr::contains("_RAW_"), dplyr::contains("_SCR_"))
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  src_des_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Format Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

addDesignSeqCG = function(tib, seq, add, din=NULL,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'addDesignSeqCG'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    add_sym <- add %>% rlang::sym()
    seq_sym <- seq %>% rlang::sym()
    
    ret_tib <- tib %>% 
      dplyr::mutate(!!add:=stringr::str_replace(!!seq_sym, "\\[[a-zA-Z][a-zA-Z]\\]", "[CG]"))
    
    if (!is.null(din))
      ret_tib <- ret_tib %>% dplyr::mutate(
        !!din := !!seq_sym %>%
          stringr::str_remove("\\].*$") %>%
          stringr::str_remove("^.*\\[") )
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

replaceDesignSeqCG = function(tib, seq, add, nuc,
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'replaceDesignSeqCG'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    add_sym <- add %>% rlang::sym()
    seq_sym <- seq %>% rlang::sym()
    nuc_sym <- nuc %>% rlang::sym()
    
    ret_tib <- tib %>% 
      dplyr::mutate(PreSeq=stringr::str_remove(!!seq_sym, "\\[.*$"),
                    PosSeq=stringr::str_remove(!!seq_sym, "^.*\\]"),
                    !!add:=paste0(PreSeq,"[",!!nuc_sym,"]",PosSeq) ) %>%
      dplyr::select(-PreSeq, -PosSeq)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          BSMAP Alignment Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsmapProbeAlign = function(exe, fas, gen, dir, 
                           opt="-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R", 
                           lan=NULL, run=FALSE,
                           verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'bsmapProbeAlign'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} dir={dir}.{RET}"))
  
  stime <- system.time({
    if (is.null(opt) || length(opt)==0)
      opt="-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
    
    fas_name <- fas %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
    gen_name <- gen %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
    
    # Create Shell and Alignment directories::
    sh_dir <- file.path(dir, 'shells')
    al_dir <- file.path(dir, 'align')
    
    if (!dir.exists(sh_dir)) dir.create(sh_dir, recursive=TRUE)
    if (!dir.exists(al_dir)) dir.create(al_dir, recursive=TRUE)
    
    aln_ssh <- file.path(sh_dir, paste0('run_bsp-',fas_name,'-',gen_name,'.sh') )
    aln_bsp <- file.path(al_dir, paste0(fas_name,'-',gen_name,'.bsmap.bsp') )
    ret_tsv <- file.path(al_dir, paste0(fas_name,'-',gen_name,'.bsmap.formatted.tsv.gz') )
    aln_cmd <- paste(exe, '-a',fas, '-d',gen, opt,'-o',aln_bsp, sep=' ')
    aln_cmd <- glue::glue("{aln_cmd}{RET}",
                          "gzip {aln_bsp}{RET}",
                          "gzip -dc {aln_bsp}.gz | cut -f 1,2,4-11 | perl -pe 's/_/\t/gi; s/:/\t/gi' | ",
                          "gzip -c - > {ret_tsv}{RET}")
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{TAB} Launching; CMD={aln_cmd}...{RET}"))
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{TAB} Launching bsmap alignments: {fas_name} vs. {gen_name}...{RET}"))
    readr::write_lines(x=aln_cmd, file=aln_ssh, append=FALSE)
    Sys.chmod(paths=aln_ssh, mode="0777")
    
    if (!is.null(lan)) aln_ssh <- paste(lan,aln_ssh, sep=' ')
    if (run) base::system(aln_ssh)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  if (run) return(ret_tsv)
  aln_ssh
}

# End of file
