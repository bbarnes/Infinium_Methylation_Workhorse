
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Big Science::
#
# Sequencing -> Discovery -> Array Dominace (Can't compete)
# COVID: Sample Prep -> Gary; do you know about Direct Detection
#   - Existing array
#   - GC Content difference with infectious cells
#   - Normalization (Noob) applied to both arrays and sequencing
#   - SPIT F-ING FIRE!!!
#
# VACINE PASSPORT REPORT CARD!!!
#
# Should improve after new cgnDB is built...
#
# use Digest::SHA qw(sha1_hex);
# my $var = 123;
# my $sha1_hash = sha1_hex($var);
# print $sha1_hash;
#

rm(list=ls(all=TRUE))

# Genomic Ranges::
suppressWarnings(suppressPackageStartupMessages( base::require("GenomicRanges",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("plyr")) )
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringi") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("data.table") ))
suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Global Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Define Default Params and Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
opt <- NULL
par <- NULL

# Illumina based directories::
par$date    <- Sys.Date() %>% as.character()
par$runMode <- ''
par$maxTest <- NULL
par$macDir1 <- '/Users/bbarnes/Documents/Projects/methylation'
par$macDir2 <- '/Users/bretbarnes/Documents'
par$lixDir1 <- '/illumina/scratch/darkmatter'
par$lixDir  <- '/illumina/scratch/darkmatter'

# Program Parameters::
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'workhorse'
par$prgmTag <- 'workhorse_main_dev_seq'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))


# Executables::
opt$Rscript <- NULL

# BSMAP Parameters::
opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
opt$bsmap_exe <- "/Users/bretbarnes/Documents/programs/BSMAPz/bsmapz"

# Run Parameters::
opt$runName    <- NULL
opt$Species    <- NULL

# Null Place Holders::
opt$cpg_top_tsv <- NULL
opt$cpg_pos_tsv <- NULL
opt$cph_pos_tsv <- NULL
opt$snp_pos_tsv <- NULL
opt$org_des_tsv <- NULL

# Directories::
opt$outDir  <- NULL
opt$impDir  <- NULL
opt$annDir  <- NULL
opt$genDir  <- NULL
opt$manDir  <- NULL

# Manufacturing Files:: Required
opt$ords <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib
opt$mats <- NULL
opt$aqps <- NULL

opt$mans <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib
opt$vcfs <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib
opt$beds <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib
opt$snps <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib

# Manufacturing Info:: Required
opt$bpns <- NULL
opt$aqpn <- NULL

# Pre-defined files (Controls & IDAT Validation):: Optional
opt$ctls <- NULL
opt$idat <- NULL

# Platform/Method Options::
opt$genBuild <- NULL
opt$platform <- NULL
opt$version  <- NULL

# Process Parallel/Cluster Parameters::
opt$single   <- FALSE
opt$parallel <- TRUE
opt$cluster  <- FALSE

# Run-time Files
opt$opt_csv  <- NULL
opt$par_csv  <- NULL
opt$time_csv <- NULL
opt$time_org_txt <- NULL

# Run-time Options::
opt$fresh   <- FALSE

# verbose Options::
opt$verbose <- 3

#
# Data Structures to be pre-defined
#
org_des_tib <- NULL
add_org_tib <- NULL
cgn_bed_tib <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Order/Match/AQP/PQC Expected Columns::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_cols <- list()
par_cols$ord <- 
  cols(
    Assay_Design_Id        = col_character(),
    AlleleA_Probe_Id       = col_character(),
    AlleleA_Probe_Sequence = col_character(),
    AlleleB_Probe_Id       = col_character(),
    AlleleB_Probe_Sequence = col_character(),
    Normalization_Bin      = col_character()
  )

par_cols$mat <- 
  cols(
    Plate    = col_character(),
    Row      = col_character(),
    Col      = col_integer(),
    Address  = col_integer(),
    Mod5     = col_character(),
    Sequence = col_character(),
    Mod3     = col_character(),
    Comments = col_character()
  )

par_cols$ma2 <- 
  cols(
    address_names = col_integer(),
    probe_id      = col_character(),
    sequence      = col_character(),
    type_b        = col_character(),
    address_name  = col_integer(),
    bo_seq        = col_character()
  )

par_cols$aqp <- 
  cols(
    Address           = col_integer(),
    Decode_Status     = col_integer(),
    Decode_Error_Code = col_integer(),
    Decode_Score      = col_integer(),
    Func_Status       = col_integer(),
    Func_Error_Code   = col_integer(),
    QC_Action         = col_integer()
  )

par_cols$pqc <- 
  cols(
    Address      = col_integer(),
    Status       = col_integer(),
    Eval_Code    = col_integer(),
    Average_Rep  = col_integer(),
    Expected_Rep = col_integer()
  )

par$ord_col <- par_cols$ord$cols %>% names()

par$mat_col <- par_cols$mat$cols %>% names()
par$ma2_col <- par_cols$ma2$cols %>% names()

par$aqp_col <- par_cols$aqp$cols %>% names()
par$pqc_col <- par_cols$pqc$cols %>% names()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Parse Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
args.dat <- commandArgs(trailingOnly = FALSE)
if (args.dat[1]=='RStudio') {
  
  if (dir.exists(par$macDir1)) par$topDir <- par$macDir1
  if (dir.exists(par$macDir2)) par$topDir <- par$macDir2
  if (dir.exists(par$lixDir1)) par$topDir <- par$lixDir1
  if (dir.exists(par$lixDir))  par$topDir <- par$lixDir
  
  if (!dir.exists(par$topDir)) dir.create(par$topDir, recursive=TRUE)
  
  par$runMode    <- args.dat[1]
  cat(glue::glue("[{par$prgmTag}]: Local args.dat[1]={args.dat[1]}.{RET}"))
  cat(glue::glue("[{par$prgmTag}]: Local      runMode={par$runMode}.{RET}"))
  cat(glue::glue("[{par$prgmTag}]: Local       topDir={par$topDir}.{RET}"))
  
  # Default Parameters for local Mac::
  par$srcDir     <- file.path(par$topDir, 'tools', par$codeDir)
  par$scrDir     <- file.path(par$srcDir, 'scripts')
  par$exePath    <- file.path(par$scrDir, 'R', par$prgmDir, paste0(par$prgmTag,'.R'))
  
  par$prgmTag <- base::sub('\\.R$', '', base::basename(par$exePath))
  par$locPath <- base::dirname(par$exePath)
  par$scrDir  <- base::dirname(base::normalizePath(par$locPath) )
  par$srcDir  <- base::dirname(base::normalizePath(par$scrDir) )
  par$datDir  <- file.path(base::dirname(base::normalizePath(par$srcDir)), 'dat')
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda2-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda3-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/conda_4.6.8/bin/Rscript'
  
  #
  # End of local parameter definitions::
  #
  
  opt$outDir  <- file.path(par$topDir, 'scratch')
  opt$impDir  <- file.path(par$topDir, 'data/improbe')
  opt$annDir  <- file.path(par$topDir, 'data/annotation')
  opt$manDir  <- file.path(par$topDir, 'data/manifests')
  opt$genDir  <- file.path(par$topDir, 'data/iGenomes/Homo_sapiens/NCBI')
  # opt$idatDir <- file.path(par$topDir, 'data/idats')
  
  opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
  opt$bsmap_exe <- "/Users/bretbarnes/Documents/programs/BSMAPz/bsmapz"
  
  # Pre-defined local options runTypes::
  #
  par$local_runType <- NULL
  par$local_runType <- 'NZT'
  par$local_runType <- 'COVIC'
  par$local_runType <- 'GSA'
  par$local_runType <- 'HM450'
  par$local_runType <- 'TruDx'
  par$local_runType <- 'EWAS'
  par$local_runType <- 'Chicago'
  par$local_runType <- 'GRCm10'
  
  if (par$local_runType=='EWAS') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'EPIC'
    opt$version  <- 'A2'
    opt$Species  <- "Human"
    
    opt$idat   <- NULL
    opt$ords   <- NULL
    # opt$ords   <- paste(
    #   "/Users/bretbarnes/Documents/data/CustomContent/EWAS/orders/round1/EPIC-reorder.partition-*.order.csv.gz",
    #   sep=','
    # )
    tsv_files <- list.files(file.path(par$topDir, "data/CustomContent/EPIC_v2/11102020"), pattern=".tsv.gz$", full.names=TRUE)
    csv_files <- c(
      list.files(file.path(par$topDir, "data/CustomContent/EWAS/orders/round1"), pattern=".order.csv.gz$", full.names=TRUE),
      list.files(file.path(par$topDir, "data/CustomContent/EWAS/orders/round2"), pattern=".order.csv.gz$", full.names=TRUE)
    )
    
    new_dir <- file.path(par$topDir, "data/CustomContent/EPIC_v2/11102020/csv")
    if (!dir.exists(new_dir)) dir.create(new_dir, recursive=TRUE)
    
    for (tsv in tsv_files) {
      cur_nam <- base::basename(tsv) %>% 
        stringr::str_remove('.gz') %>% 
        stringr::str_remove('.tsv') %>%
        stringr::str_remove("^bp_")
      
      cur_nam <- paste0("9",cur_nam)
      
      out_nam <- paste("EPIC-reorder.partition",cur_nam, sep='-')
      out_csv <- file.path(new_dir, paste(out_nam,'order.csv.gz', sep='.'))
      cur_tib <- suppressMessages(suppressWarnings( readr::read_tsv(tsv, guess_max = 50000) ))
      
      cat(glue::glue("[Pre-process]: name={cur_nam}; inp_csv={tsv}{RET}"))
      
      readr::write_csv(cur_tib, out_csv)
      
      cat(glue::glue("[Pre-process]: name={cur_nam}; out_csv={out_csv}{RET}{RET}"))
      # print(cur_tib)
    }
    new_files <- list.files(new_dir, pattern=".csv.gz$", full.names=TRUE)
    
    opt$ords <- c(csv_files,new_files)
    par$ordn <- c(csv_files,new_files) %>% base::basename() %>% 
      stringr::str_remove("^EPIC-reorder.partition-") %>%
      stringr::str_remove(".order.csv.gz$")
    
    opt$mats   <- NULL
    opt$aqps   <- NULL
    
    opt$bpns   <- NULL
    opt$aqpn   <- NULL
    
    opt$mans   <- NULL
    # opt$mans   <- paste(
    #   file.path(opt$manDir, "genotyping/GSA-24v2-0_A1.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz"),
    # file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz"),
    # sep=',')
    
    opt$vcfs   <- NULL
    # opt$vcfs   <- paste(
    # "/Users/bretbarnes/Documents/data/annotation/dbSNP/dbSNP-151/GRCh37/common_all_20180423.snps.csv.gz",
    # "/Users/bretbarnes/Documents/data/annotation/dbSNP/dbSNP-151/GRCh37/All_20180423.vcf.gz",
    #  sep=',')
    opt$beds   <- NULL
    opt$snps   <- NULL
    # opt$snps   <- paste(
    #   "/Users/bretbarnes/Documents/data/CustomContent/TruDx/SNPs/TruDx_target_SNPs.csv",
    #   sep=',')
    
    opt$org_des_tsv <- NULL
    
  } else if (par$local_runType=='TruDx') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'GSA'
    opt$version  <- 'A4'
    opt$Species  <- "Human"
    
    opt$idat   <- NULL
    opt$ords   <- NULL
    opt$mats   <- NULL
    opt$aqps   <- NULL
    
    opt$bpns   <- NULL
    opt$aqpn   <- NULL
    
    opt$mans   <- NULL
    # opt$mans   <- paste(
    #   file.path(opt$manDir, "genotyping/GSA-24v2-0_A1.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz"),
    # file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz"),
    # sep=',')
    
    opt$vcfs   <- paste(
      # "/Users/bretbarnes/Documents/data/annotation/dbSNP/dbSNP-151/GRCh37/common_all_20180423.snps.csv.gz",
      "/Users/bretbarnes/Documents/data/annotation/dbSNP/dbSNP-151/GRCh37/All_20180423.vcf.gz",
      sep=',')
    opt$beds   <- NULL
    opt$snps   <- paste(
      "/Users/bretbarnes/Documents/data/CustomContent/TruDx/SNPs/TruDx_target_SNPs.csv",
      sep=',')
    
    opt$org_des_tsv <- NULL
    
  } else if (par$local_runType=='HM450') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'GSA'
    opt$version  <- 'A2'
    opt$Species  <- "Human"
    
    opt$idat   <- NULL
    opt$ords   <- NULL
    opt$mats   <- NULL
    opt$aqps   <- NULL
    
    opt$bpns   <- NULL
    opt$aqpn   <- NULL
    
    opt$mans   <- NULL
    opt$mans   <- paste(
      #   file.path(opt$manDir, "genotyping/GSA-24v2-0_A1.csv.gz"),
      #   file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz"),
      file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz"),
      #   file.path(opt$manDir, "methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz"),
      sep=',')
    
    opt$vcfs   <- NULL
    opt$beds   <- NULL
    
    opt$org_des_tsv <- NULL
  } else if (par$local_runType=='GSA') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'GSA'
    opt$version  <- 'A2'
    opt$Species  <- "Human"
    
    opt$idat   <- NULL
    opt$ords   <- NULL
    opt$mats   <- NULL
    opt$aqps   <- NULL
    
    opt$bpns   <- NULL
    opt$aqpn   <- NULL
    
    # opt$mans   <- paste(
    #   file.path(opt$manDir, "genotyping/GSA-24v2-0_A1.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz"),
    #   sep=',')
    
    opt$mans   <- NULL
    opt$vcfs   <- NULL
    opt$beds   <- NULL
    
    opt$org_des_tsv <- NULL
    
  } else if (par$local_runType=='Chicago') {
    opt$genBuild <- 'GRCh38'
    opt$genBuild <- 'GRCh37'
    
    opt$platform <- 'EPIC'
    opt$version  <- 'A1'
    opt$version  <- 'A2'
    opt$version  <- 'A3'
    opt$version  <- 'A4'
    opt$version  <- 'A5'
    opt$version  <- 'A6'
    opt$version  <- 'A7'
    opt$version  <- 'A8'
    opt$version  <- 'A9'
    
    opt$Species  <- "Human"
    
    opt$idat   <- NULL
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/UnivChicago/latest')
    
    opt$ords <- paste(
      file.path(par$aqpDir, 'UofChicago-A_A_Array-CpG-order-FINAL.csv'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20504790_probes.match.tsv'),
      sep=',')
    
    opt$aqps <- paste(
      file.path(par$aqpDir, '329922X374054_A_ProductQC.txt'),
      sep=',')
    
    opt$bpns <- paste(1, sep=",")
    opt$aqpn <- paste(1, sep=",")
    
    opt$org_des_tsv <- file.path(par$topDir, "data/CustomContent/UnivChicago/improbe_input/CpGs_UnivChicago_alldesigns_55860sites.cgn-pos-srd-prbs.tsv.gz")
    
  } else if (par$local_runType=='COVIC') {
    opt$genBuild <- 'GRCh36'
    opt$genBuild <- 'GRCh38'
    opt$genBuild <- 'GRCh37'
    
    opt$platform    <- 'COVIC'
    opt$version     <- 'C0'
    opt$version     <- 'B1'
    opt$version     <- 'C2'
    opt$version     <- 'C11'
    opt$version     <- 'C12'
    
    opt$Species <- "Human"
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    opt$pre_man_csv <- '/Users/bretbarnes/Documents/data/manifests/MethylationEPIC_v-1-0_B2.csv.gz'
    
    opt$idat  <- NULL
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/COVID-19_HLA/AQP/COVIC-Host-Immune-Detection')
    opt$ords <- paste(
      file.path(par$aqpDir, 'COVID_EPIC_Round1.03172020.unique.order_AP.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20447043_probes.match.gz'),
      sep=',')
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'BS0032581-AQP.txt.gz'),
      sep=',')
    
    opt$bpns <- paste(1, sep=",")
    opt$aqpn <- paste(1, sep=",")
    
  } else if (par$local_runType=='GRCm10') {
    opt$platform <- 'LEGX'
    opt$version  <- 'C0'
    opt$version  <- 'C30'
    opt$version  <- 'C31'
    opt$version  <- 'C32'
    opt$version  <- 'C0'
    opt$version  <- 'P0'
    opt$version  <- 'P1'
    
    opt$genBuild <- 'GRCm38'
    opt$genBuild <- 'GRCm10'
    
    opt$Species <- "Mouse"
    
    opt$genDir  <- file.path(par$topDir, 'data/iGenomes/Mus_musculus/NCBI')
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
    opt$cph_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.chn-sorted.tsv.gz') )
    opt$snp_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.snp-sorted.tsv.gz') )
    opt$org_des_tsv <- file.path(par$topDir, "data/CustomContent/LifeEpigentics/data/dropbox/merged_with_raw_ordered.cgn-pos-srd-prbs.tsv.gz")
    
    # opt$snp_des_csv <- file.path(opt$impDir, 'cph-snp-designs/LEGX_SpikeIn_Reorder-SNP-Only.designs.csv.gz')
    # opt$cph_des_csv <- file.path(opt$impDir, 'cph-snp-designs/LEGX_SpikeIn_Reorder-CpH-Only.designs.csv.gz')
    
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/LifeEpigentics/AQP')
    opt$ords <- paste(
      file.path(par$aqpDir, 'orders/Mus_musculus.order_BP1.csv.gz'),
      file.path(par$aqpDir, 'orders/Mus_musculus.order_BP2.csv.gz'),
      file.path(par$aqpDir, 'orders/mm10_LEGX_nonCpG_probes.Jan16-2020.order.csv.gz'),
      file.path(par$aqpDir, 'orders/LEGX_SpikeIn_Reorder-All-06052020.order.withHeader.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, 'BP1/20420178_AQP1_LifeEpigen_BP1.txt.gz'),
      file.path(par$aqpDir, 'BP2/20420260_AQP1_LifeEpigen_BP2.txt.gz'),
      file.path(par$aqpDir, 'BP3/20420260_AQP2_LifeEpigen_BP2.txt.gz'),
      file.path(par$aqpDir, 'BP4/20455357_AQP1_LifeEpigen_BP4.txt.gz'),
      sep=',')
    
    if (opt$version == 'P0') {
      opt$aqps <- paste(
        file.path(par$aqpDir, 'PQC/20042400_A_ProductQC.txt.gz'),
        sep=',')
      opt$aqpn <- paste(1, sep=",")
      
    } else {
      
      opt$aqps <- paste(
        file.path(par$aqpDir, 'AQP_Copy/BS0032527-AQP.txt.gz'),
        file.path(par$aqpDir, 'AQP_Copy/BS0032533-AQP.txt.gz'),
        file.path(par$aqpDir, 'AQP_Copy/BS0032545-AQP.txt.gz'),
        file.path(par$aqpDir, 'AQP_Copy/BS0032636-AQP.txt.gz'),
        sep=',')
      opt$aqpn <- paste(1,1,2,1, sep=",")
    }
    
    opt$bpns <- paste(1,2,2,3, sep=",")
    
  } else if (par$local_runType=='NZT') {
    opt$genBuild <- 'GRCh36'
    opt$genBuild <- 'GRCh38'
    opt$genBuild <- 'GRCh37'
    opt$platform    <- 'NZT'
    opt$version     <- 'N0'
    
    opt$Species <- "Human"
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/NZT/decode')
    opt$ords <- paste(
      file.path(par$aqpDir, 'selected.order1.csv.gz'),
      file.path(par$aqpDir, 'selected.order2.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20297484_probes.match1.tsv.gz'),
      file.path(par$aqpDir, '20297484_probes.match2.tsv.gz'),
      sep=',')
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'BS0031918-AQP1.txt.gz'),
      file.path(par$aqpDir, 'BS0032272-AQP2.txt.gz'),
      sep=',')
    
    opt$bpns <- paste(1,2, sep=",")
    opt$aqpn <- paste(1,2, sep=",")
    
    opt$idat <- NULL
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$parallel <- TRUE
  opt$runName <- paste(par$local_runType,opt$platform,opt$version,opt$genBuild, sep='-')
  
  # opt$fresh <- TRUE
  opt$fresh <- FALSE
  
  opt$verbose <- 10
  opt$verbose <- 3
  
} else {
  par$runMode    <- 'CommandLine'
  par$exePath <- base::substring(args.dat[grep("--file=", args.dat)], 8)
  
  par$prgmTag <- base::sub('\\.R$', '', base::basename(par$exePath))
  par$locPath <- base::dirname(par$exePath)
  par$scrDir  <- base::dirname(base::normalizePath(par$locPath) )
  par$srcDir  <- base::dirname(base::normalizePath(par$scrDir) )
  par$datDir  <- file.path(base::dirname(base::normalizePath(par$srcDir)), 'dat')
  
  args.dat <- commandArgs(trailingOnly = TRUE)
  option_list = list(
    # Executables::
    make_option(c("--Rscript"), type="character", default=opt$Rscript, 
                help="Rscript path [default= %default]", metavar="character"),
    
    make_option(c("--bsmap_opt"), type="character", default=opt$bsmap_opt, 
                help="BSMAP Options [default= %default]", metavar="character"),
    make_option(c("--bsmap_exe"), type="character", default=opt$bsmap_exe, 
                help="BSMAP Executable path [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--Species"), type="character", default=opt$Species, 
                help="Target Species Name [default= %default]", metavar="character"),
    make_option(c("--pre_man_csv"), type="character", default=opt$pre_man_csv, 
                help="Previously defined manifest [default= %default]", metavar="character"),
    
    make_option(c("--cpg_pos_tsv"), type="character", default=opt$cpg_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--cpg_top_tsv"), type="character", default=opt$cpg_top_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--cph_pos_tsv"), type="character", default=opt$cph_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--snp_pos_tsv"), type="character", default=opt$snp_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--org_des_tsv"), type="character", default=opt$org_des_tsv, 
                help="Original design file used for canonical position selection. [default= %default]", metavar="character"),
    
    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    
    make_option(c("--impDir"), type="character", default=opt$impDir, 
                help="improbe data directory [default= %default]", metavar="character"),
    make_option(c("--annDir"), type="character", default=opt$iannDir, 
                help="Annotation data directory [default= %default]", metavar="character"),
    make_option(c("--genDir"), type="character", default=opt$genDir, 
                help="Genomic data directory [default= %default]", metavar="character"),
    make_option(c("--manDir"), type="character", default=opt$manDir, 
                help="Pre-built Manifest data directory [default= %default]", metavar="character"),
    
    # Manufacturing Files:: Required
    make_option(c("--ords"), type="character", default=opt$ords, 
                help="Order file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mats"), type="character", default=opt$mats, 
                help="Match (format 1 or 2) file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqps"), type="character", default=opt$aqps, 
                help="AQP file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mans"), type="character", default=opt$mans,
                help="Manifest CSV file(s) (comma seperated) [default= %default]", metavar="character"),
    
    # Not fully supported yet::
    make_option(c("--vcfs"), type="character", default=opt$vcfs, 
                help="Target Design VCF file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--beds"), type="character", default=opt$beds, 
                help="Target Design Coordinate BED file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--snps"), type="character", default=opt$snps, 
                help="Target Design Coordinate SNP file(s) (comma seperated) [default= %default]", metavar="character"),
    
    # Manufacturing Info:: Required
    make_option(c("--bpns"), type="character", default=opt$bpns, 
                help="Bead Pool Numbers (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqpn"), type="character", default=opt$aqpn, 
                help="AQP Numbers (comma seperated) [default= %default]", metavar="character"),
    
    # Pre-defined files (Controls & IDAT Validation):: Optional
    make_option(c("--ctls"), type="character", default=opt$ctls, 
                help="Pre-Defined Infinium Methylation Controls (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--idat"), type="character", default=opt$idat, 
                help="idat directories (comma seperated) [default= %default]", metavar="character"),
    
    # Platform/Method Options::
    make_option(c("--genBuild"), type="character", default=opt$genBuild, 
                help="Genome Build (e.g. GRCh36, GRCh37, GRCh38, GRCm38) [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform (e.g. HM450, EPIC, LEGX, NZT, COVIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest Version (e.g. B0,B1,B2,B3,B4,C0) [default= %default]", metavar="character"),
    
    
    # Process Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
    # Run-time Files::
    make_option(c("--opt_csv"), type="character", default=opt$opt_csv, 
                help="Unused variable opt_csv [default= %default]", metavar="character"),
    make_option(c("--par_csv"), type="character", default=opt$par_csv, 
                help="Unused variable par_csv [default= %default]", metavar="character"),
    make_option(c("--time_csv"), type="character", default=opt$time_csv, 
                help="Unused variable time_csv [default= %default]", metavar="character"),
    make_option(c("--time_org_txt"), type="character", default=opt$time_org_txt, 
                help="Unused variable time_org_txt [default= %default]", metavar="character"),
    
    # Run=time Options::
    make_option(c("--fresh"), action="store_true", default=opt$fresh, 
                help="Boolean variable to run a fresh build [default= %default]", metavar="boolean"),
    
    # verbose::
    make_option(c("-v", "--verbose"), type="integer", default=opt$verbose, 
                help="0-5 (5 is very verbose) [default= %default]", metavar="integer")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Program Initialization::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_reqs <- c('runMode','prgmTag','scrDir','datDir','exePath')
opt_reqs <- c('outDir','impDir','Species',
              # 'ords',
              'genBuild','platform','version','bsmap_exe', # 'bsmap_opt',
              'Rscript','verbose')

# par$gen_src_dir <- file.path(par$scrDir, 'functions')
# if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
# for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
# cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )

par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir))
  stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))

for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', 
                         full.names=TRUE, recursive=TRUE)) base::source(sfile)
if (opt$verbose>=0)
  cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form ",
                 "General Source={par$gen_src_dir}!{RET}{RET}") )

opt <- program_init(name=par$prgmTag,
                    opts=opt, opt_reqs=opt_reqs, 
                    pars=par, par_reqs=par_reqs,
                    libs=TRUE,rcpp=FALSE,
                    verbose=opt$verbose,vt=3,tc=0,tt=NULL)

par_tib <- par %>%
  unlist(recursive = TRUE) %>%
  dplyr::bind_rows() %>% 
  tidyr::gather("Params", "Value")
opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")

stamp_vec <- NULL
stamp_vec <- c(stamp_vec,opt$time_org_txt)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Parse List Options
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

image_key <- "bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse_Centos"
image_ver <- "v.1.15"
image_ver <- "v.1.25"
image_ssh <- "run_improbe.sh"
image_str <- glue::glue("{image_key}.{image_ver}")

# TBD:: Make this more general Manifest Control Defaults::
#
if (is.null(opt$ctls)) {
  opt$ctls <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Pre-processing:: Run Time:: Output Directories
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run <- NULL

run$manDir <- file.path(opt$outDir, 'man')
if (!dir.exists(run$manDir)) dir.create(run$manDir, recursive=TRUE)

run$desDir <- file.path(opt$outDir, 'des')
if (!dir.exists(run$desDir)) dir.create(run$desDir, recursive=TRUE)

run$addDir <- file.path(opt$outDir, 'add')
if (!dir.exists(run$addDir)) dir.create(run$addDir, recursive=TRUE)

run$intDir <- file.path(opt$outDir, 'int')
if (!dir.exists(run$intDir)) dir.create(run$intDir, recursive=TRUE)

run$fasDir <- file.path(opt$outDir, 'fas')
if (!dir.exists(run$fasDir)) dir.create(run$fasDir, recursive=TRUE)

run$alnDir <- file.path(opt$outDir, 'aln')
if (!dir.exists(run$alnDir)) dir.create(run$alnDir, recursive=TRUE)

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Building Output Directories.{RET}{RET}"))


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Pre-processing:: Run Time:: Intermediate Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Define Run Time:: Ref Alignment Genome
run$gen_ref_fas <- 
  file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
            paste0(opt$genBuild,".genome.fa.gz"))

# Define Run Time:: SNP IUPAC Genome
run$gen_snp_fas <- 
  file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
            paste0(opt$genBuild,".dbSNP151-genome.fa.gz"))

# Define Pre-built improbe files
run$imp_prb_dir <- file.path(opt$impDir, "scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49")
run$imp_u49_tsv <- file.path(run$imp_prb_dir, paste("probe_U49_cgn-table.csv.gz", sep="-") )
run$imp_m49_tsv <- file.path(run$imp_prb_dir, paste("probe_M49_cgn-table.csv.gz", sep="-") )
stopifnot(dir.exists(run$imp_prb_dir))
stopifnot(file.exists(run$imp_u49_tsv))
stopifnot(file.exists(run$imp_m49_tsv))

# Defined Run Time:: Intermediate Files
run$aqp_man_csv  <- file.path(run$manDir, paste(opt$runName,"aqp-pass.manifest-sesame.csv.gz", sep="."))

run$aqp_add_csv  <- file.path(run$addDir, paste(opt$runName,"aqp-pass.address.csv.gz", sep="."))
run$aqp_prb_fas  <- file.path(run$fasDir, paste(opt$runName,"aqp-pass.address.fas.gz",  sep='.'))
run$aqp_u49_tsv  <- file.path(run$intDir, paste(opt$runName,"aqp-pass.address-u49.tsv.gz", sep='.') )
run$aqp_m49_tsv  <- file.path(run$intDir, paste(opt$runName,"aqp-pass.address-m49.tsv.gz", sep='.') )

run$aqp_prb_bsp  <- file.path(run$alnDir, paste(opt$runName,"aqp-pass.address.bsp",  sep='.') )
run$aqp_bsp_tsv  <- file.path(run$alnDir, paste(opt$runName,"aqp-pass.address-bsp.tsv.gz",  sep='.') )

run$int_u49_tsv <- file.path(run$intDir, paste(opt$runName, "aqp-u49.intersect.tsv.gz", sep='.') )
run$int_m49_tsv <- file.path(run$intDir, paste(opt$runName, "aqp-m49.intersect.tsv.gz", sep='.') )
run$int_seq_tsv <- file.path(run$intDir, paste(opt$runName, "aqp-seq.intersect.tsv.gz", sep='.') )

run$imp_inp_tsv  <- file.path(run$desDir, paste(opt$runName, 'improbe-inputs.tsv.gz', sep='.') )
run$imp_des_tsv  <- file.path(run$desDir, paste(opt$runName, 'improbe-design.tsv.gz', sep='.') )

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Defining Run Time Files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  2.1 Functional Manifest Generation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par$retData <- FALSE

#
# TBD:: Critical:: Pass in Ord_Key or Address as main key!!!
#   - Testing right now
#   - Issue is passing in Ord_Key can create too long (20 char) for improbe
#   - Issue is poor uniqueness of keys!!!
#

stamp_vec <- c(stamp_vec,run$aqp_add_csv)
if (opt$fresh || !valid_time_stamp(stamp_vec)) {
  aqp_add_tib <- 
    aqp_address_workflow(
      ord=opt$ords,
      mat=opt$mats, aqp=opt$aqps, 
      bpn=opt$bpns, aqn=opt$aqns,
      add_csv=run$aqp_add_csv, 
      man_csv=run$aqp_man_csv,
      add_fas=run$aqp_prb_fas,
      u49_tsv=run$aqp_u49_tsv, m49_tsv=run$aqp_m49_tsv,
      runName=opt$runName, retData=par$retData,
      verbose=opt$verbose, tt=pTracker)
} else {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading aqp_add={run$aqp_add_csv}...{RET}"))
  aqp_add_tib <- suppressMessages(suppressWarnings( 
    readr::read_csv(run$aqp_add_csv, guess_max=100000) )) %>%
    clean_tibble()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    3.2 Align All Probe Sequence:: BSMAP
#                  3.3 Join Address and Alignment Data:: BSMAP
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

stamp_vec <- c(stamp_vec, run$aqp_bsp_tsv)
if (opt$fresh || !valid_time_stamp(stamp_vec)) {
  
  aqp_bsp_tib <- 
    run_bsmap(
      exe=opt$bsmap_exe, 
      fas=run$aqp_prb_fas, 
      gen=run$gen_ref_fas,
      bsp=run$aqp_prb_bsp,
      
      add=aqp_add_tib,bed=cgn_bed_tib,org=add_org_tib,
      
      join_key="Aln_Key",join_type="inner",
      des_key="Ord_Des", din_key="Ord_Din",
      sort=TRUE,

      opt=NULL, lan=NULL, run=TRUE,
      verbose=opt$verbose,tt=pTracker)
  
  safe_write(aqp_bsp_tib,"csv",run$aqp_bsp_tsv, funcTag=par$prgmTag,
             verbose=opt$verbose, tt=pTracker)
  
} else {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading aqp_bsp_tsv={run$aqp_bsp_tsv}...{RET}"))
  
  aqp_bsp_tib <- suppressMessages(suppressWarnings(
    readr::read_csv(run$aqp_bsp_tsv, guess_max=100000) )) %>%
    clean_tibble()
}

if (FALSE) {
  man_bsp_tib <- aqp_bsp_tib %>% 
    dplyr::filter(Bsp_Tag=="UM") %>%
    add_to_man(join = "Ord_Key", runName = opt$runName, 
               des_key = "Ord_Des", pid_key = "Aln_Key", 
               # rep_key = "Ord_Par_Rep", rep_val = "Ord_Par_Rep", 
               # col_key = "Ord_Col", 
               # nxb_key = "Bsp_Nxb_Ref", 
               verbose = opt$verbose, tt=pTracker)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#         3.4 Intersect Sequences Address and improbe:: U49/M49
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

seq_imp_tib <- NULL
stamp_vec <- c(stamp_vec, 
               run$int_u49_tsv,
               run$int_m49_tsv,
               run$int_seq_tsv)
if (opt$fresh || !valid_time_stamp(stamp_vec)) {
  
  int_u49_tib2 <- 
    intersect_seq(ref=run$imp_u49_tsv,
                  can=run$aqp_u49_tsv,
                  out=run$int_u49_tsv,
                  idxA=1, idxB=1, 
                  verbose=opt$verbose,tt=pTracker)
  
  int_m49_tib2 <- 
    intersect_seq(ref=run$imp_m49_tsv,
                  can=run$aqp_m49_tsv, 
                  out=run$int_m49_tsv,
                  idxA=1, idxB=1, 
                  verbose=opt$verbose,tt=pTracker)
  
  seq_imp_tib <- 
    join_seq_intersect(u49=int_u49_tib2, m49=int_m49_tib2, 
                       bed=cgn_bed_tib, org=add_org_tib,
                       verbose=opt$verbose, tt=pTracker)
  
  safe_write(seq_imp_tib,"tsv",run$int_seq_tsv, funcTag=par$prgmTag,
             verbose=opt$verbose)
  
} else {
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading imp_seq_tsv={run$int_seq_tsv}...{RET}"))
  seq_imp_tib <- 
    suppressMessages(suppressWarnings( readr::read_tsv(run$int_seq_tsv) )) %>%
    clean_file()

}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 3.4.1 Bind BSMAP & Seq-Match into Table::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Row Binding all_imp_tab=(seq_imp_tib/bsp_imp_tib)...{RET}"))

com_des_cols <- c("Address","Ord_Des","Ord_Din")
com_gen_cols <- c("Imp_Chr","Imp_Pos","Imp_FR","Imp_TB","Imp_CO")
com_scr_cols <- c("Can_Mis_Scr","Mat_Src_Scr","Org_Mis_Scr","Bsp_Din_Scr")
com_cgn_cols <- c("Imp_Cgn","Imp_Nxb","Imp_Din") # ,"Bsp_Din_Ref")

imp_col_vec <- c(com_des_cols,com_gen_cols,com_scr_cols,com_cgn_cols)

seq_cgn_tib <- seq_imp_tib %>%
  dplyr::filter(is.na(Imp_Chr))

seq_dat_tib <- seq_imp_tib %>%
  dplyr::filter(!is.na(Imp_Chr)) %>%
  dplyr::mutate(Imp_Din="CG", Bsp_Din_Scr=as.integer(8)) %>%
  dplyr::mutate(Mat_Src_Scr=as.integer(0), Mat_Src_Key="Seq") %>%
  dplyr::select(dplyr::any_of(imp_col_vec),"Mat_Src_Key")

bsp_dat_tib <- bsp_imp_tib %>%
  dplyr::rename(Imp_Din=Bsp_Din_Ref) %>%
  dplyr::mutate(Mat_Src_Scr=as.integer(1), Mat_Src_Key="Bsp") %>%
  dplyr::select(dplyr::any_of(imp_col_vec),"Mat_Src_Key")

imp_col_vec2 <- c(com_des_cols,com_gen_cols,com_scr_cols,com_cgn_cols,"Mat_Src_Key")

all_imp_tab <- NULL
all_imp_tab <- dplyr::bind_rows(
  dplyr::select( seq_dat_tib, dplyr::any_of(imp_col_vec2) ),
  dplyr::select( bsp_dat_tib, dplyr::any_of(imp_col_vec2) ),
)  %>% 
  dplyr::mutate(
    Can_Mis_Scr=dplyr::case_when(
      is.na(Can_Mis_Scr) ~ 9,
      TRUE ~ as.double(Can_Mis_Scr)
    ) %>% as.integer(),
    Org_Mis_Scr=dplyr::case_when(
      is.na(Org_Mis_Scr) ~ 9,
      TRUE ~ as.double(Org_Mis_Scr)
    ) %>% as.integer(),
    Bsp_Din_Scr=dplyr::case_when(
      is.na(Bsp_Din_Scr) ~ 9,
      TRUE ~ as.double(Bsp_Din_Scr)
    ) %>% as.integer()
  ) %>% 
  dplyr::arrange( dplyr::across( dplyr::all_of(c(com_scr_cols)) ) ) %>%
  dplyr::distinct()

all_imp_tab %>% 
  dplyr::arrange(Address,Ord_Des,Ord_Din,Imp_Chr,Imp_Pos,Imp_FR,Imp_TB,Imp_CO) %>%
  dplyr::add_count(Address,Ord_Des,Ord_Din,Imp_Chr,Imp_Pos,Imp_FR,Imp_TB,Imp_CO, name="Group_Count") %>% 
  dplyr::filter(Group_Count>1) %>%
  print(n=1000)

all_imp_tab %>% 
  dplyr::add_count(Address,Ord_Des,Ord_Din,Imp_Chr,Imp_Pos,Imp_FR,Imp_TB,Imp_CO, name="Group_Count") %>% 
  dplyr::group_by(Group_Count,Mat_Src_Key,Bsp_Din_Scr,Ord_Des) %>% 
  dplyr::summarise(Count=n(), .groups="drop") %>% 
  print(n=1000)

# Max Best Extension::
#   - This requires the template sequence and probe designs
#     It can be inferred for BSP alignments, but requires full templates and 
#     designs for Seq data. 
#   - Question; is there ever a case where a probe matches non-seq matching
#     more than seq-matching. Yes; poorly designed probes. This is why we need
#     all the template sequences followed by designs.
#
# Max Canonical CGN
#
# Max Source (Seq,Bsp)
# Max U over M
#

# TBD::
#  - Validate against Wanding's manifest:: cgn/pos
#  - Intersect with Chrom HMM/UCSC
#  - Pick Max Wanding Annotation
#  - Compare Internal Annotation vs. Wanding Annotation
#




if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Row Binding all_imp_tab=(seq_imp_tib/bsp_imp_tib).{RET}{RET}"))



if (FALSE) {
  #
  # Spare scratch for TruResponse
  #
  
  tru_tib <- readr::read_tsv("/Users/bretbarnes/Documents/data/CustomContent/TruReponse/sampleSheets/saliva-barcodes.txt")
  
  tru_tib %>% 
    dplyr::group_by(Sample_Name,Sample_Long,Sample_Day,Sample_Hour,Sample_Prep) %>% 
    dplyr::summarise(Count=n(), .groups="drop") %>% 
    dplyr::arrange(-Count,Sample_Name) %>% 
    print(n=1000)
  
  tru_tib %>% 
    dplyr::group_by(Sample_Name,Sample_Long,Sample_Day,Sample_Prep) %>% 
    dplyr::summarise(Count=n(), .groups="drop") %>% 
    dplyr::arrange(-Count,Sample_Name) %>% 
    print(n=1000)
  
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       4.0 Write improbe input::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
fwd_des_tib <- aqp_bsp_tib %>% # head(n=10) %>%
  fas_to_seq(fas = run$gen_ref_fas, file = run$imp_inp_tsv, name = "Aln_Key", 
             din = "Ord_Din", gen = opt$genBuild, 
             chr1 = "Bsp_Chr", pos = "Bsp_Pos",
             # nrec = 1,
             verbose=opt$verbose+10, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     5.0 Run improbe designs:: docker
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imp_des_tib <- 
  improbe_docker(
    dir=run$desDir, file=run$imp_inp_tsv, 
    name=opt$runName, image=image_str, shell=image_ssh,
    level=3, add_inf=TRUE,
    verbose=opt$verbose, tt=pTracker)

# %>%
#   dplyr::mutate(Methyl_Allele_TB_Strand=stringr::str_sub(Methyl_Allele_TB_Strand,1,1))
#

# Qucik QC::
#  imp_des_tib %>% dplyr::filter(Inf_Type != 0) %>% dplyr::distinct(Seq_ID)


# All srd probe extraction::
#  - Seesm to be working, will validate later...
if (FALSE) {

  all_srd_prb_tib <- aqp_bsp_tib %>% head(n=10) %>%
    bed_to_prbs(fas=run$gen_ref_fas, 
                din="Ord_Din", cgn="Aln_Key",
                chr="Bsp_Chr", pos="Bsp_Pos",
                nrec = 1,
                verbose=opt$verbose+10, tt=pTracker)

}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
