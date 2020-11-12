

if (FALSE) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                              Source Packages::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # rm(list=ls(all=TRUE))
  
  # Load sesame:: This causes issues with "ExperimentHub Caching causes a warning"
  # suppressWarnings(suppressPackageStartupMessages( base::require("sesame") ))
  # suppressWarnings(suppressPackageStartupMessages( base::require("dbplyr") ))
  
  suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))
  
  suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
  suppressWarnings(suppressPackageStartupMessages( base::require("plyr")) )
  suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
  suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))
  
  suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
  suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))
  
  # Parallel Computing Packages
  suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #            EPIC hg19 Infinium I Inferred SNP Table Construction::
  #
  #               Should move this scratch to somewhere else::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # 450k Genome Studio::
  #
  # hm45_man_csv <- '/Users/bretbarnes/Documents/data/manifests/HumanMethylation450_15017482_v.1.2.csv.gz'
  # hm45_man_lst <- loadManifestGenomeStudio(file=hm45_man_csv, verbose=opt$verbose+4)
  # hm45_man_tib <- hm45_man_lst$man %>% dplyr::select(IlmnID,Genome_Build,CHR,MAPINFO,Strand)
  # hm45_ctl_tib <- hm45_man_lst$ctl %>% dplyr::mutate(Address=as.double(Address))
  
  #
  # SNP Table
  #
  epicI_ses_grs <- sesameDataPullVariantAnno_InfiniumI(platform = "EPIC")
  epicR_ses_grs <- sesameDataPullVariantAnno_SNP()
  
  epicI_ses_tib <- epicI_ses_grs %>% 
    as.data.frame() %>% tibble::rownames_to_column(var="IlmnID") %>% tibble::as_tibble()
  
  # Five Examples by hand::
  #  gs_ses_csv <- '/Users/bretbarnes/Documents/tools/notes/inferred-snps.GS-vs-Ses.csv'
  #  gs_ses_tib <- readr::read_csv(gs_ses_csv)
  #
  #
  
  # epic_man_csv <- file.path(par$datDir, 'manifest/base/EPIC-B4.manifest.sesame-base.cpg-sorted.csv.gz')
  epic_man_csv <- '/Users/bretbarnes/Documents/data/manifests/MethylationEPIC_v-1-0_B4.csv.gz'
  epic_man_lst <- loadManifestGenomeStudio(file=epic_man_csv, verbose=opt$verbose+4)
  epic_ctl_tib <- epic_man_lst$ctl %>% dplyr::mutate(Address=as.double(Address))
  epic_man_tib <- epic_man_lst$man %>% 
    dplyr::select(IlmnID,Genome_Build,CHR,MAPINFO,Strand,Next_Base,Infinium_Design_Type) %>%
    dplyr::filter(Infinium_Design_Type=='I') %>% 
    dplyr::mutate(CHR=paste0('chr',CHR)) %>% dplyr::arrange(CHR,MAPINFO)
  
  #
  # Validation Shown Below::
  #
  join_man_tib <- epic_man_tib %>% 
    dplyr::inner_join(epicI_ses_tib, by="IlmnID", suffix=c("_Man", "_Ses")) %>%
    dplyr::mutate(REF_CON=dplyr::case_when(
      strand=='+' & REF=='G' ~ 'A',
      strand=='-' & REF=='C' ~ 'T',
      TRUE ~ REF)
    )
  join_sum_tib <- join_man_tib %>% dplyr::filter(CHR==seqnames) %>% 
    dplyr::mutate(Pos_Dif=MAPINFO-start) %>% 
    dplyr::group_by(Pos_Dif,Strand,strand,Next_Base,REF_CON) %>% 
    dplyr::summarise(Count=n(), .groups='drop')
  
  #
  # Better, but NOT working attempt below::
  #
  if (FALSE) {
    join_man_tib3 <- epic_man_tib %>% 
      dplyr::inner_join(epicI_ses_tib, by="IlmnID", suffix=c("_Man", "_Ses")) %>%
      dplyr::mutate(
        REF_RC=dplyr::case_when(
          strand=='-' ~ revCmp(REF),
          TRUE ~ REF),
        REF_CO=dplyr::case_when(
          strand=='+' & REF_RC=='G' ~ 'A',
          strand=='-' & REF_RC=='C' ~ 'T',
          TRUE ~ REF)
      )
    join_sum_tib3 <- join_man_tib3 %>% 
      dplyr::filter(CHR==seqnames) %>% 
      dplyr::mutate(Pos_Dif=MAPINFO-start) %>% 
      dplyr::group_by(Pos_Dif,Strand,strand,Next_Base,REF_CO) %>% 
      dplyr::summarise(Count=n(), .groups='drop')
  }
  
  #
  # hub19 doesn't have the most up to date dbSNP
  #
  hub_hg19 <- subset(hub, (hub$species == "Homo sapiens") & (hub$genome == "hg19"))
  length(hub_hg19)
  hub_hg19$title[grep("SNP", hub_hg19$title)]
  
  hub_snp137 <- subset(hub_hg19, title=='Common SNPs(137)')
  
  hub_snp137_grs <- hub_snp137[['AH5105']]
  hub_snp137_tib <- hub_snp137_grs %>% as.data.frame() %>% tibble::as_tibble()
  
  #
  # Check overlap of Sesame and dbSNP(137) list
  #
  epicI_ses_tib %>% dplyr::inner_join(hub_snp137_tib, by=c("rs"="name"), suffix=c("_sesI", "_db137"))
  
  #
  # Waiting to download full hg19 dbSNP 2018-04-18 vcf.gz (likely to be put in trash)
  #
  
  # genom_snp_csv <- '/Users/bretbarnes/Documents/data/annotation/GRCh37/All_20180423.cut1-7.vcf.gz'
  # genom_snp_tib <- readr::read_csv(file=genom_snp_csv, comment="##")
  
  db151_hs37_col <- cols(CHR = col_character(),
                         POS = col_double(),
                         SNP = col_character(),
                         REF = col_character(),
                         ALT = col_character())
  
  # db151_hs37_tsv <- '/Users/bretbarnes/Documents/data/annotation/GRCh37/All_20180423.snp-5.tsv.gz'
  db151_hs37_tsv <- '/Users/bretbarnes/Documents/data/annotation/GRCh37/All_20180423.snp-5.chr1Stsv.gz'
  db151_hs37_tib <- readr::read_tsv(db151_hs37_tsv, col_names=names(db151_hs37_col$cols), col_types=db151_hs37_col)
}

# End of file
