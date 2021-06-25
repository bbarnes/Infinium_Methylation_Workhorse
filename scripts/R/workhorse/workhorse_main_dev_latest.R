
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
BNG <- "|"

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
par$prgmTag <- 'workhorse_main_dev_latest'
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
opt$ord_des_csv <- NULL

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

par$buildManifest <- FALSE
par$run_improbe   <- FALSE
par$load_ann <- FALSE

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
  par$local_runType <- 'GRCm10'
  par$local_runType <- 'EWAS'
  par$local_runType <- 'Chicago'
  
  if (par$local_runType=='EWAS') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'EWAS'
    opt$version  <- 'A1'
    opt$version  <- 'A2'
    opt$Species  <- "Human"
    
    opt$idat   <- NULL
    opt$ords   <- NULL
    # opt$ords   <- paste(
    #   "/Users/bretbarnes/Documents/data/CustomContent/EWAS/orders/round1/EPIC-reorder.partition-*.order.csv.gz",
    #   sep=','
    # )
    new_files <- NULL
    tsv_files <- list.files(file.path(par$topDir, "data/CustomContent/EPIC_v2/11102020"), pattern=".tsv.gz$", full.names=TRUE)
    csv_files <- c(
      list.files(file.path(par$topDir, "data/CustomContent/EWAS/orders/round1"), pattern=".order.csv.gz$", full.names=TRUE),
      list.files(file.path(par$topDir, "data/CustomContent/EWAS/orders/round2"), pattern=".order.csv.gz$", full.names=TRUE)
    )
    
    if (FALSE) {
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
    }
    
    # opt$ords <- c(csv_files,new_files)
    opt$ords <- c(csv_files)
    par$ordn <- opt$ords %>% base::basename() %>% 
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
    
    # opt$ord_des_csv <- file.path(par$datDir, "manifest/cgnDB/canonical-assignment.cgn-top-grp.csv.gz")
    opt$ord_des_csv <- NULL
    
    par$buildManifest <- TRUE
    par$run_improbe   <- TRUE

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
    
    opt$ord_des_csv <- NULL
    
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
    
    opt$ord_des_csv <- NULL
    
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
    
    opt$ord_des_csv <- NULL
    
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
    opt$version  <- 'A10'
    opt$version  <- 'A11'
    
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
    
    opt$ord_des_csv <- file.path(par$topDir, "data/CustomContent/UnivChicago/improbe_input/CpGs_UnivChicago_alldesigns_55860sites.cgn-pos-srd-prbs.tsv.gz")
    
    par$buildManifest <- TRUE
    par$run_improbe   <- TRUE

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
    opt$version  <- 'P2'
    
    opt$genBuild <- 'GRCm38'
    opt$genBuild <- 'GRCm10'
    
    opt$Species <- "Mouse"
    
    opt$genDir  <- file.path(par$topDir, 'data/iGenomes/Mus_musculus/NCBI')
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
    opt$cph_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.chn-sorted.tsv.gz') )
    opt$snp_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.snp-sorted.tsv.gz') )
    opt$ord_des_csv <- file.path(par$topDir, "data/CustomContent/LifeEpigentics/data/dropbox/merged_with_raw_ordered.cgn-pos-srd-prbs.tsv.gz")
    
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
    
    par$load_ann <- TRUE
    
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
  opt$fresh   <- FALSE
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
    make_option(c("--ord_des_csv"), type="character", default=opt$ord_des_csv, 
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

pTracker <- timeTracker$new()
# pTracker <- timeTracker$new(verbose=opt$verbose)

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

run$annDir <- file.path(opt$outDir, 'ann')
if (!dir.exists(run$annDir)) dir.create(run$annDir, recursive=TRUE)

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

if (is.null(opt$ord_des_csv))
  opt$ord_des_csv <- file.path(par$datDir, "manifest/cgnDB/canonical-assignment.cgn-top-grp.csv.gz")

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
run$imp_des_tsv  <- file.path(run$desDir, paste(opt$runName, 'improbe-designOutput.tsv.gz', sep='.') )
run$cln_des_csv  <- file.path(run$desDir, paste(opt$runName, 'improbe-designOutput.clean.tsv.gz', sep='.') )

# TBD:: This needs to be changed in the docker image to replace the lines above!
# run$imp_des_tsv  <- file.path(run$desDir, paste(opt$runName, 'improbe-design.tsv.gz', sep='.') )
# run$cln_des_csv  <- file.path(run$desDir, paste(opt$runName, 'improbe-design.clean.csv.gz', sep='.') )

run$ann_int_csv  <- file.path(run$annDir, paste(opt$runName, 'cpg-pass.annotation.csv.gz', sep='.') )

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Defining Run Time Files.{RET}{RET}"))

if (par$buildManifest) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                  2.1 Functional Manifest Generation::
  #                         AQP Address Workflow()
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  par$retData <- TRUE
  par$retData <- FALSE
  # opt$verbose <- 10
  
  aqp_add_tib <- NULL
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
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]: Loading aqp_add={run$aqp_add_csv}...{RET}"))
    aqp_add_tib <- suppressMessages(suppressWarnings( 
      readr::read_csv(run$aqp_add_csv, guess_max=100000) )) %>%
      clean_tibble()
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #         3.4 Intersect Sequences Address and improbe:: U49/M49
  #                         CGN Mapping Workflow()
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  seq_cgn_tib <- NULL
  stamp_vec <- c(stamp_vec, 
                 run$int_u49_tsv,
                 run$int_m49_tsv,
                 run$int_seq_tsv)
  
  if (opt$fresh || !valid_time_stamp(stamp_vec)) {
    
    if (par$local_runType=="Chicago") opt$ord_des_csv <- NULL
    
    seq_cgn_tib <- cgn_mapping_workflow(
      ref_u49=run$imp_u49_tsv,can_u49=run$aqp_u49_tsv,out_u49=run$int_u49_tsv,
      ref_m49=run$imp_m49_tsv,can_m49=run$aqp_m49_tsv,out_m49=run$int_m49_tsv,
      idxA=1, idxB=1,
      ord=opt$ord_des_csv,
      verbose=opt$verbose,tt=pTracker
    )

    # TBD:: The writing function should be moved to cgn_mapping_workflow()
    safe_write(seq_cgn_tib,"tsv",run$int_seq_tsv, funcTag=par$prgmTag,
               verbose=opt$verbose)
  } else {
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]: Loading imp_seq_tsv={run$int_seq_tsv}...{RET}"))
    seq_cgn_tib <- 
      suppressMessages(suppressWarnings( readr::read_tsv(run$int_seq_tsv) )) %>%
      clean_tibble()
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    3.2 Align All Probe Sequence:: BSMAP
  #                  3.3 Join Address and Alignment Data:: BSMAP
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  aqp_bsp_tib <- NULL
  stamp_vec <- c(stamp_vec, run$aqp_bsp_tsv)
  if (opt$fresh || !valid_time_stamp(stamp_vec)) {
    
    aqp_bsp_tib <- 
      run_bsmap(
        exe=opt$bsmap_exe,fas=run$aqp_prb_fas, 
        gen=run$gen_ref_fas,bsp=run$aqp_prb_bsp,
        add=aqp_add_tib,bed=cgn_bed_tib,org=add_org_tib,
        join_key="Aln_Key",join_type="inner",
        des_key="Ord_Des", din_key="Ord_Din",
        sort=TRUE,opt=NULL,lan=NULL,run=TRUE,
        verbose=opt$verbose,tt=pTracker)
    
    # TBD:: The writing function should be moved to run_bsmap()
    safe_write(aqp_bsp_tib,"csv",run$aqp_bsp_tsv, funcTag=par$prgmTag,
               verbose=opt$verbose, tt=pTracker)
  } else {
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]: Loading aqp_bsp_tsv={run$aqp_bsp_tsv}...{RET}"))
    
    aqp_bsp_tib <- suppressMessages(suppressWarnings(
      readr::read_csv(run$aqp_bsp_tsv, guess_max=100000) )) %>%
      clean_tibble()
  }
  
  #
  # TBD:: 
  #  - CGN Look up by coordinates from BSP (this will take a bit more time to code)
  #    BED File: /Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-input/GRCh37.cgn.bed.gz
  #

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       4.0 improbe fwd design::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  imp_des_tib <- NULL
  stamp_vec <- c(stamp_vec, 
                 run$imp_inp_tsv,
                 run$imp_des_tsv,
                 run$cln_des_csv)
  
  if (opt$fresh || !valid_time_stamp(stamp_vec)) {
    
    imp_des_tib <- imp_designs_workflow(
      tib=aqp_bsp_tib,fas=run$gen_ref_fas,
      imp=run$imp_inp_tsv,gen=opt$genBuild,
      des=run$imp_des_tsv,out=run$cln_des_csv,
      
      name=opt$runName,image=image_str,shell=image_ssh,
      verbose=opt$verbose, tt=pTracker)

  } else {
    
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]: Loading cln_des_csv={run$cln_des_csv}...{RET}"))
    imp_des_tib <- 
      suppressMessages(suppressWarnings( readr::read_csv(run$cln_des_csv) )) %>%
      clean_tibble()
  }
  
  # Qucik QC::
  #  imp_des_tib %>% dplyr::filter(Inf_Type != 0) %>% dplyr::distinct(Seq_ID)
  
}

# NOTE: Now we have four core tables::
#  Address:   aqp_add_tib  => aqp_address_workflow()
#  CG Number: seq_cgn_tib  => cgn_mapping_workflow()
#  Alignment: aqp_bsp_tib  => run_bspmap()
#  Design:    imp_des_tib  => imp_designs_workflow()

# aqp_add_tib %>% print()
# seq_cgn_tib %>% print()
# aqp_bsp_tib %>% print()
# imp_des_tib %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 3.4.1 aqp_add_tib with seq/unq_cgn_tib::
#                           Add Canonical Order::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Add cgn verification by coordinate matching???
#
add_cgn_inn <- NULL
add_cgn_inn <- dplyr::inner_join(
  aqp_add_tib,seq_cgn_tib,
  by=c("Ord_Des","Ord_Din","Address","Aln_P49"="Aln_Prb")
)

add_cgn_ant <- dplyr::anti_join(
  aqp_add_tib,seq_cgn_tib,
  by=c("Ord_Des","Ord_Din","Address","Aln_P49"="Aln_Prb")
)
add_cgn_ant_cnt <- add_cgn_ant %>% base::nrow()
if (opt$verbose>0) 
  cat(glue::glue("[{par$prgmTag}]: add_cgn_ant_cnt={add_cgn_ant_cnt}. ",
                 "Should be zero.{RET}"))

# TBD: Add check for Canonical data (Can_Scr)
#
# Determine Optimal Sorting::
#
if (!is.null(add_cgn_inn)) {
  add_cgn_inn1 <- add_cgn_inn %>%
    dplyr::arrange(Can_Scr,-Imp_Hit_hg37) %>%
    dplyr::distinct(Address,Ord_Des,Ord_Din,Ord_Prb, .keep_all=TRUE) %>% 
    dplyr::mutate(
      Mat_Cgn=dplyr::case_when(
        Imp_Cgn==Ord_Cgn ~ 0,
        TRUE ~ 1) %>% as.integer()
    )
  add_cgn_sum1 <- add_cgn_inn1 %>%
    # dplyr::group_by(Mat_Cgn,Aqp_Idx,Ord_Des,Ord_Din) %>%
    dplyr::group_by(Mat_Cgn) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  add_cgn_sum1 %>% print(n=base::nrow(add_cgn_sum1))
  rm(add_cgn_inn1)
  rm(add_cgn_sum1)
  
  add_cgn_inn2 <- add_cgn_inn %>%
    dplyr::arrange(Can_Scr,Imp_Hit_hg37) %>%
    dplyr::distinct(Address,Ord_Des,Ord_Din,Ord_Prb, .keep_all=TRUE) %>% 
    dplyr::mutate(
      Mat_Cgn=dplyr::case_when(
        Imp_Cgn==Ord_Cgn ~ 0,
        TRUE ~ 1) %>% as.integer()
    )
  add_cgn_sum2 <- add_cgn_inn2 %>%
    # dplyr::group_by(Mat_Cgn,Aqp_Idx,Ord_Des,Ord_Din) %>%
    dplyr::group_by(Mat_Cgn) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  add_cgn_sum2 %>% print(n=base::nrow(add_cgn_sum2))
  rm(add_cgn_inn2)
  rm(add_cgn_sum2)
  
  add_cgn_inn3 <- add_cgn_inn %>%
    dplyr::arrange(-Can_Scr,Imp_Hit_hg37) %>%
    dplyr::distinct(Address,Ord_Des,Ord_Din,Ord_Prb, .keep_all=TRUE) %>% 
    dplyr::mutate(
      Mat_Cgn=dplyr::case_when(
        Imp_Cgn==Ord_Cgn ~ 0,
        TRUE ~ 1) %>% as.integer()
    )
  add_cgn_sum3 <- add_cgn_inn3 %>%
    # dplyr::group_by(Mat_Cgn,Aqp_Idx,Ord_Des,Ord_Din) %>%
    dplyr::group_by(Mat_Cgn) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  add_cgn_sum3 %>% print(n=base::nrow(add_cgn_sum3))
  rm(add_cgn_inn3)
  rm(add_cgn_sum3)

  add_cgn_inn4 <- add_cgn_inn %>%
    dplyr::arrange(-Can_Scr,-Imp_Hit_hg37) %>%
    dplyr::distinct(Address,Ord_Des,Ord_Din,Ord_Prb, .keep_all=TRUE) %>% 
    dplyr::mutate(
      Mat_Cgn=dplyr::case_when(
        Imp_Cgn==Ord_Cgn ~ 0,
        TRUE ~ 1) %>% as.integer()
    )
  add_cgn_sum4 <- add_cgn_inn4 %>%
    # dplyr::group_by(Mat_Cgn,Aqp_Idx,Ord_Des,Ord_Din) %>%
    dplyr::group_by(Mat_Cgn) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  add_cgn_sum4 %>% print(n=base::nrow(add_cgn_sum4))

  #
  # CONCLUSION::
  #   BEST=   dplyr::arrange(-Can_Scr,-Imp_Hit_hg37)
  # Select Best:: currently = add_cgn_inn4
  #
  add_cgn_inn <- add_cgn_inn4
  rm(add_cgn_inn4)
  rm(add_cgn_sum4)
}

if (FALSE) {
  # A tibble: 1,141,866 x 89
  add_cgn_inn %>% dplyr::inner_join(aqp_bsp_tib, by=c("Aln_Prb"))
  
  # A tibble: 1,033,006 x 88
  add_cgn_inn %>% dplyr::inner_join(aqp_bsp_tib, by=c("Address","Aln_Prb"))
  
  # A tibble: 1,033,006 x 86
  add_cgn_inn %>% dplyr::inner_join(aqp_bsp_tib, by=c("Address","Ord_Des","Ord_Din","Aln_Prb"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   3.4.2 aqp_add_tib with aqp_bsp_tib::
#                        Add Sesame pos comparison::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# BSP Success summary:: Not technically needed, should be moved to run_bsmap()
#
# aqp_bsp_sum <- aqp_bsp_tib %>% 
#   dplyr::group_by(Address,Ord_Des,Ord_Din,Ord_Prb) %>%
#   dplyr::summarise(Bsp_Din_Scr_Min=min(Bsp_Din_Scr, na.rm=TRUE),
#                    Bsp_Din_Scr_Avg=mean(Bsp_Din_Scr, na.rm=TRUE),
#                    Bsp_Din_Scr_Med=median(Bsp_Din_Scr, na.rm=TRUE),
#                    Bsp_Din_Scr_Max=max(Bsp_Din_Scr, na.rm=TRUE),
#                    Bsp_Din_Scr_Cnt=n(),
#                    .groups="drop"
#   )

#
# Full Join of add, cgn, imp, bsp::
# NOTE: For Chicago this blows up real quick::
#
# A tibble: 66,880,588 x 47
add_cgn_imp_bsp_inn <- 
  dplyr::inner_join(add_cgn_inn,imp_des_tib,
                    by=c("Aln_Key"="Seq_ID",
                         "Imp_TB"="Strand_TB",
                         "Imp_CO"="Strand_CO")
  ) %>% 
  #
  # TBD:: aqp_bsp_tib is already joined with aqp_add_tib earlier in the bsp
  #  generation. We should remove some of this redundancy...
  #
  dplyr::left_join(
    aqp_bsp_tib %>% dplyr::select(Address,Ord_Des,Ord_Din,
                                  Aln_Key,Ord_Prb, dplyr::starts_with("Bsp_")),
    by=c("Address","Ord_Des","Ord_Din","Aln_Key","Ord_Prb",
         "Chromosome"="Bsp_Chr","Coordinate"="Bsp_Pos")
  )

# These three numbers below demonstrate uniqueness across the board::
#
um_cnt1 <- add_cgn_imp_bsp_inn %>%
  dplyr::filter(Bsp_Tag=="UM") %>% base::nrow()

um_cnt2 <- add_cgn_imp_bsp_inn %>%
  dplyr::filter(Bsp_Tag=="UM") %>%
  dplyr::distinct(Address,Ord_Des,Ord_Din,Ord_Prb,Chromosome,Coordinate) %>% 
  base::nrow()

um_cnt3 <- add_cgn_imp_bsp_inn %>%
  dplyr::filter(Bsp_Tag=="UM") %>%
  dplyr::distinct(Address,Ord_Des,Ord_Din,Ord_Prb) %>% base::nrow()

#
# Here are the Multi Unique Probes::
#
ma_cnt1 <- add_cgn_imp_bsp_inn %>%
  dplyr::filter(Bsp_Tag!="UM") %>%
  dplyr::distinct(Address,Ord_Des,Ord_Din,Ord_Prb) %>% base::nrow()

cat(glue::glue("[{par$prgmTag}]: um_cnt={um_cnt1}={um_cnt2}={um_cnt3}, ",
               "ma_cnt={ma_cnt1}.{RET}{RET}"))

# add_cgn_imp_bsp_inn %>% dplyr::select(Ord_Key,Bsp_Din_Ref,Top_Sequence)
# add_cgn_imp_bsp_inn %>% dplyr::select(Ord_Key,Bsp_Din_Ref,Bsp_Din_Bsc,Bsp_Din_Scr,Top_Sequence) %>% dplyr::filter(Bsp_Din_Ref=="CG") %>% as.data.frame()

# Good summary check::
#
add_cgn_imp_bsp_sum <- add_cgn_imp_bsp_inn %>% 
  dplyr::arrange(Bsp_Din_Scr) %>%
  dplyr::distinct(Address,Ord_Des,Ord_Din,Ord_Prb, .keep_all=TRUE) %>%
  # dplyr::group_by(Bsp_Din_Scr,Ord_Des,Ord_Din,Aqp_Idx) %>%
  dplyr::group_by(Bsp_Din_Scr) %>%
  dplyr::summarise(Count=n(), .groups="drop")
add_cgn_imp_bsp_sum %>% print(n=base::nrow(add_cgn_imp_bsp_sum))

#
# Now Filter/Summarize::
#
# add_cgn_imp_bsp_inn %>%
#   dplyr::group_by(Address,Ord_Des,Ord_Din,Ord_Prb) %>%
#   dplyr::summarise(Bsp_Din_Scr_Min=min(Bsp_Din_Scr, na.rm=TRUE),
#                    Bsp_Din_Scr_Avg=mean(Bsp_Din_Scr, na.rm=TRUE),
#                    Bsp_Din_Scr_Med=median(Bsp_Din_Scr, na.rm=TRUE),
#                    Bsp_Din_Scr_Max=max(Bsp_Din_Scr, na.rm=TRUE),
#                    Bsp_Din_Scr_Cnt=n(),
#                    .groups="drop"
#   )

#
# Simple Manifest Generation::
#  TBD:: This needs to be improved ALOT!!!
#
add_cgn_imp_bsp_man <- NULL
add_cgn_imp_bsp_man <- 
  add_cgn_imp_bsp_inn %>% 
  dplyr::filter(Bsp_Tag=="UM") %>%
  add_to_man(join=c("Ord_Key","Ord_Din","Ord_Col"),
             runName=opt$runName,
             des_key="Ord_Des", pid="Ord_Key",
             col_key="Ord_Col",
             # csv=man_csv,
             validate=TRUE,
             verbose=10) %>% 
  dplyr::group_by(Ord_Prb_U) %>%
  dplyr::mutate(
    Rank=dplyr::row_number()
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Strand_FR_U=dplyr::case_when(
      Strand_FR_U=="F" ~ "+", 
      Strand_FR_U=="R" ~ "-", 
      TRUE ~ NA_character_),
    IlmnID=paste0(
      Imp_Cgn_U,"_",Imp_TB_U,Imp_CO_U,Infinium_Design,Rank)
      # Imp_Cgn_U,"_",Imp_TB_U,Imp_CO_U,Infinium_Design,Ord_Prb_Rep_U)
  ) %>% 
  dplyr::distinct(IlmnID, .keep_all=TRUE) %>%
  dplyr::arrange(Chromosome_U,Coordinate_U)

man_pos_grs <- 
  GenomicRanges::GRanges(
    seqnames=Rle(add_cgn_imp_bsp_man$Chromosome_U), 
    strand=Rle(add_cgn_imp_bsp_man$Strand_FR_U),
    Probe_Type=add_cgn_imp_bsp_man$Inf_Type_U,
    
    IRanges(start=add_cgn_imp_bsp_man$Coordinate_U, 
            end=add_cgn_imp_bsp_man$Coordinate_U+1, 
            names=add_cgn_imp_bsp_man$IlmnID)
  )

# Now we need to review scores::
#
# add_cgn_imp_bsp_man$Bsp_Din_Scr_M
# add_cgn_imp_bsp_man$Bsp_Din_Scr_U
#
# add_cgn_imp_bsp_man$Cpg_Scr_M
# add_cgn_imp_bsp_man$Cpg_Scr_U
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       6.0 Annotation Conformation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
if (opt$genBuild=="GRCh37" || opt$genBuild=="GRCh38") {
  hmm_ann_drn <- paste("GRCh36",opt$genBuild, sep='-')
  hmm_ann_dir <- file.path(par$topDir, "data/annotation/liftOver/chrom_hmm/wgEncodeBroadHmm/ucsc_liftover_main",hmm_ann_drn)
  hmm_ann_fns <- list.files(hmm_ann_dir, pattern=".map.bed.gz$", full.names=TRUE, recursive=FALSE)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       6.3 Load Annotation:: NCBI/UCSC
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  core_anno_dir <- file.path(par$topDir, "data/annotation", opt$genBuild)
  ncib_gene_tsv <- file.path(core_anno_dir, paste(opt$genBuild,"ncbi.RefSeqGenes.tsv.gz", sep='.'))
  ucsc_gene_tsv <- file.path(core_anno_dir, paste(opt$genBuild,"ucsc.knownGene.tsv.gz", sep='.'))
  ucsc_cpgs_tsv <- file.path(core_anno_dir, paste(opt$genBuild,"ucsc.CpG-Islands.tsv.gz", sep='.'))
  
  ncbi_gene_tib <- load_ncbi_gene(file=ncib_gene_tsv, verbose=opt$verbose, tt=pTracker)
  ucsc_gene_tib <- load_ucsc_gene(file=ucsc_gene_tsv, verbose=opt$verbose, tt=pTracker)
  ucsc_cpgs_tib <- load_ucsc_cpgs(file=ucsc_cpgs_tsv, verbose=opt$verbose, tt=pTracker)
  
  ncbi_gene_grs <- load_ncbi_gene(file=ncib_gene_tsv, grs=TRUE, verbose=opt$verbose, tt=pTracker)
  ucsc_gene_grs <- load_ucsc_gene(file=ucsc_gene_tsv, grs=TRUE, verbose=opt$verbose, tt=pTracker)
  ucsc_cpgs_grs <- load_ucsc_cpgs(file=ucsc_cpgs_tsv, grs=TRUE, verbose=opt$verbose, tt=pTracker)
  
  ann_int_list <- NULL
  # NCBI Gene Comparison::
  #
  ref_pre_str <- "NCBI_Gene"
  ann_key_str <- ref_pre_str
  ref_pre_str <- NULL
  # ncbi_gene_int_tib <- 
  ann_int_list[[ann_key_str]] <-
    intersect_GRS(can=man_pos_grs, ref=ncbi_gene_grs, 
                  can_key="IlmnID", ref_prefix=ref_pre_str, 
                  verbose=opt$verbose, tt=pTracker)
  
  if (FALSE) {
    ncbi_gene_crs_tib <- 
      ann_int_list[[ann_key_str]] %>% 
      dplyr::inner_join(ses_ann_tab, by="IlmnID")
    ncbi_gene_crs_sum <- ncbi_gene_crs_tib %>% 
      dplyr::group_by(class,Value_Str) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% 
      dplyr::arrange(-Count)
    ncbi_gene_crs_sum %>% print(n=base::nrow(ncbi_gene_crs_sum))
  }
  
  # UCSC Gene Comparison::
  #
  ref_pre_str <- "UCSC_Gene"
  ann_key_str <- ref_pre_str
  ref_pre_str <- NULL
  # ucsc_gene_int_tib <- 
  ann_int_list[[ann_key_str]] <-
    intersect_GRS(can=man_pos_grs, ref=ucsc_gene_grs, 
                  can_key="IlmnID", ref_prefix=ref_pre_str, 
                  verbose=opt$verbose, tt=pTracker)
  
  if (FALSE) {
    ucsc_gene_crs_tib <- 
      ann_int_list[[ann_key_str]] %>% 
      dplyr::inner_join(ses_ann_tab, by="IlmnID")
    ucsc_gene_crs_sum <- ucsc_gene_crs_tib %>% 
      dplyr::group_by(class,Value_Str) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% 
      dplyr::arrange(-Count)
    ucsc_gene_crs_sum %>% print(n=base::nrow(ucsc_gene_crs_sum))
  }
  
  # UCSC Islands Comparison::
  #
  ref_pre_str <- "UCSC_Islands"
  ann_key_str <- ref_pre_str
  ref_pre_str <- NULL
  # ucsc_cpgs_int_tib <- 
  ann_int_list[[ann_key_str]] <-
    intersect_GRS(can=man_pos_grs, ref=ucsc_cpgs_grs, 
                  can_key="IlmnID", ref_prefix=ref_pre_str, 
                  verbose=opt$verbose, tt=pTracker)
  
  if (FALSE) {
    ucsc_cpgs_crs_tib <- 
      ann_int_list[[ann_key_str]] %>% 
      dplyr::inner_join(ses_ann_tab, by="IlmnID")
    
    ucsc_cpgs_crs_sum <- ucsc_cpgs_crs_tib %>% 
      dplyr::group_by(class,Value_Str) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% 
      dplyr::arrange(-Count)
    ucsc_cpgs_crs_sum %>% print(n=base::nrow(ucsc_cpgs_crs_sum))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       6.4 Load Annotation:: Chrom HMM
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  hmm_cols <- 
    cols(
      chr    = col_character(),
      beg    = col_integer(),
      end    = col_integer(),
      class  = col_character(),
      val1   = col_integer(),
      val2   = col_character(),
      val3   = col_integer(),
      val4   = col_integer(),
      val5   = col_character()
    )
  
  hmm_sufix <- ".map.bed.gz"
  hmm_names <- hmm_ann_fns %>% base::basename() %>% 
    stringr::str_remove(hmm_sufix) %>% 
    stringr::str_remove("-.*$") %>% 
    stringr::str_remove("^wgEncodeBroadHmm") %>% 
    stringr::str_remove("HMM$")
  
  hmm_list <- NULL
  for (ii in c(1:length(hmm_ann_fns))) {
    name <- hmm_names[ii]
    hmm_list[[name]] <- hmm_ann_fns[ii]
  }
  
  # lapply(hmm_ann_fns, suppressMessages(suppressWarnings(readr::read_tsv)))
  hmm_dat_list <- lapply(hmm_list, 
                         readr::read_tsv, 
                         col_names=names(hmm_cols$cols), 
                         col_types=hmm_cols)
  
  for (samp in base::names(hmm_dat_list)) {
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]: Sample={samp}...{RET}"))
    
    cur_key <- paste(samp,"hmm", sep="_")
    cur_grp <- paste(cur_key,"class", sep="_")
    cur_grp_sym <- rlang::sym(cur_grp)
    cur_name <- paste(opt)
    
    cur_tib <- hmm_dat_list[[samp]] %>%
      dplyr::arrange(chr,beg) %>% 
      dplyr::group_by(class) %>% 
      dplyr::mutate(class=stringr::str_replace_all(class, '_','-'), 
                    Class_Rank=dplyr::row_number(), 
                    Uniq_Id=paste(class,Class_Rank, sep="_")) %>%
      dplyr::ungroup()
    
    cur_grs <- 
      GenomicRanges::GRanges(
        seqnames=Rle(cur_tib$chr),
        # strand=Rle(ses_pos_tib$srd),
        # Probe_Type=ses_pos_tib$Probe_Type,
        name=paste(cur_tib$chr,cur_tib$beg,cur_tib$end, sep='-'),
        name2=NA_character_,
        class=cur_tib$class,
        source="Chrom-HMM",
        tissue=samp,
        rank=cur_tib$Class_Rank,
        
        IRanges(start=cur_tib$beg,
                end=cur_tib$end,
                names=cur_tib$Uniq_Id)
      )
    
    if (opt$verbose>0) {
      cat(glue::glue("[{par$prgmTag}]: Sample={samp}; cur_grs={RET}"))
      print(cur_grs)
    }
    
    # Comparison::
    #
    ref_pre_str <- cur_key
    ann_key_str <- cur_key
    ref_pre_str <- NULL
    
    # cur_int_tib <- 
    ann_int_list[[ann_key_str]] <-
      intersect_GRS(can=man_pos_grs, ref=cur_grs,
                    can_key="IlmnID", ref_prefix=ref_pre_str, 
                    verbose=opt$verbose, tt=pTracker)
    
    if (FALSE) {
      
      cur_crs_tib <- 
        ann_int_list[[ann_key_str]] %>% 
        dplyr::inner_join(ses_ann_tab, by="IlmnID")
      cur_crs_sum <- cur_crs_tib %>% 
        # dplyr::group_by(!!cur_grp_sym,Value_Str) %>% 
        dplyr::group_by(class,Value_Str) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>% 
        dplyr::arrange(-Count)
      cur_crs_sum %>% print(n=base::nrow(cur_crs_sum))
      
    }
    
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]: Done. Sample={samp}.{RET}{RET}"))
    
    # break
  }
}











#
#  Validation of Manifest/BSP Coordinates via Sesame
#
par$validateSesame <- FALSE
if (par$validateSesame) {
  
  #
  # Load Sesame equivlent and compare coordinates and names by sequence
  #
  par$sesBuild <- NULL
  if (opt$genBuild=="GRCh37") par$sesBuild <- "hg19"
  if (opt$genBuild=="GRCh38") par$sesBuild <- "hg18"
  ses_man_grs <- sesameData::sesameDataGet(paste("EPIC",par$sesBuild,"manifest", sep='.'))
  ses_man_tib <- ses_man_grs %>% as.data.frame() %>% 
    tibble::rownames_to_column(var="Ses_Cgn") %>%
    tibble::as_tibble()
  
  #
  # Compare aqp_bsp_tib vs. ses_man_tib
  #
  ses_pos_col <- c("Ses_Cgn", "seqnames", "start", "end", "probeBeg", "probeEnd", "strand", 
                   "designType", "probeType")
  ses_pos_colA <- c(ses_pos_col,"ProbeSeq_A")
  ses_pos_colB <- c(ses_pos_col,"ProbeSeq_B")
  
  ses_pos_tab <- dplyr::bind_rows(
    ses_man_tib %>% 
      dplyr::filter(designType=="I") %>%
      dplyr::select(dplyr::all_of(ses_pos_colA)) %>% 
      dplyr::rename(Ord_Prb=ProbeSeq_A) %>%
      dplyr::mutate(Ord_Prb=stringr::str_to_upper(Ord_Prb),
                    Ord_Des="U") %>%
      tibble::as_tibble(),
    ses_man_tib %>% 
      dplyr::filter(designType=="I") %>%
      dplyr::select(dplyr::all_of(ses_pos_colB)) %>% 
      dplyr::rename(Ord_Prb=ProbeSeq_B) %>%
      dplyr::mutate(Ord_Prb=stringr::str_to_upper(Ord_Prb),
                    Ord_Des="M") %>%
      tibble::as_tibble(),
    ses_man_tib %>% 
      dplyr::filter(designType=="II") %>%
      dplyr::select(dplyr::all_of(ses_pos_colA)) %>% 
      dplyr::rename(Ord_Prb=ProbeSeq_A) %>%
      dplyr::mutate(Ord_Prb=stringr::str_to_upper(Ord_Prb),
                    Ord_Des="2") %>%
      tibble::as_tibble()
  ) %>% 
    dplyr::rename(Ord_Din=probeType) %>%
    clean_tibble() %>%
    dplyr::select(Ses_Cgn,Ord_Des,Ord_Din,Ord_Prb, dplyr::everything())
  
  #
  # TBD:: We should use the manifest from above rather than raw aqp_bsp_tib
  #   results...
  #
  aqp_bsp_tib %>% 
    dplyr::inner_join(ses_pos_tab, by=c("Ord_Prb","Ord_Des","Ord_Din")) %>%
    dplyr::mutate(Pos_Dif=Bsp_Pos-start) %>% 
    dplyr::filter(Pos_Dif!=0) %>% 
    dplyr::filter(Bsp_Tag!="UM") %>%
    dplyr::group_by(Ord_Des,Ord_Din,
                    # Aqp_Idx,
                    Bsp_Nxb_Bsc,
                    Bsp_Din_Bsc,
                    # CG_F1,CG_R1,
                    Bsp_Tag) %>% 
    dplyr::summarise(Count=n(), .groups="drop") %>%
    print(n=1000)
  
  #
  # Matching by manifest::
  #
  ses_man_inn_tib <- add_cgn_imp_bsp_man %>% 
    dplyr::mutate(Ord_Prb=stringr::str_to_upper(Ord_Prb_U),
                  Ord_Des=Ord_Des_U,
                  Bsp_Chr=Chromosome_U,
                  Bsp_Pos=Coordinate_U,
                  Bsp_Tag=Bsp_Tag_U,
                  Aqp_Idx=Aqp_Idx_U,
                  Bsp_Din_Scr=Bsp_Din_Scr_U) %>%
    dplyr::select(Ord_Key,Ord_Des,Ord_Din,Ord_Prb,Bsp_Chr,Bsp_Pos,Bsp_Tag,
                  Aqp_Idx,Bsp_Din_Scr) %>%
    dplyr::inner_join(ses_pos_tab, by=c("Ord_Prb","Ord_Des","Ord_Din"))
  
  # Match Cases::
  #
  ses_man_mat_tib <- ses_man_inn_tib %>% 
    dplyr::mutate(Pos_Dif=Bsp_Pos-start) %>% 
    dplyr::filter(Pos_Dif==0) %>% 
    dplyr::select(Bsp_Din_Scr,Ord_Key,Ses_Cgn,Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx,
                  Ord_Prb,Pos_Dif,Bsp_Chr,seqnames)
  ses_man_mat_sum <- ses_man_mat_tib %>% 
    # dplyr::group_by(Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx) %>% 
    dplyr::group_by(Bsp_Tag,Bsp_Din_Scr,Ord_Des,Ord_Din,Aqp_Idx) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  ses_man_mat_sum %>% print(n=base::nrow(ses_man_mat_sum))
  
  # Miss Cases::
  #
  ses_man_mis_tib <- ses_man_inn_tib %>% 
    dplyr::mutate(Pos_Dif=Bsp_Pos-start) %>% 
    dplyr::filter(Pos_Dif!=0) %>% 
    dplyr::select(Bsp_Din_Scr,Ord_Key,Ses_Cgn,Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx,
                  Ord_Prb,Pos_Dif,Bsp_Chr,seqnames)
  ses_man_mis_sum <- ses_man_mis_tib %>% 
    # dplyr::group_by(Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx) %>% 
    dplyr::group_by(Bsp_Tag,Bsp_Din_Scr,Ord_Des,Ord_Din,Aqp_Idx) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  ses_man_mis_sum %>% print(n=base::nrow(ses_man_mis_sum))
  
  
  #
  # Probe Order Sequence match to check position:: BSP vs. Sesame::
  #
  ses_bsp_inn_tib <- aqp_bsp_tib %>% 
    dplyr::select(Ord_Key,Ord_Des,Ord_Din,Ord_Prb,Bsp_Chr,Bsp_Pos,Bsp_Tag,
                  Aqp_Idx,Bsp_Din_Scr) %>%
    dplyr::inner_join(ses_pos_tab, by=c("Ord_Prb","Ord_Des","Ord_Din"))
  
  # Match Cases::
  #
  ses_bsp_mat_tib <- ses_bsp_inn_tib %>% 
    dplyr::mutate(Pos_Dif=Bsp_Pos-start) %>% 
    dplyr::filter(Pos_Dif==0) %>% 
    dplyr::select(Bsp_Din_Scr,Ord_Key,Ses_Cgn,Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx,
                  Ord_Prb,Pos_Dif,Bsp_Chr,seqnames)
  ses_bsp_mat_sum <- ses_bsp_mat_tib %>% 
    # dplyr::group_by(Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx) %>% 
    dplyr::group_by(Bsp_Tag,Bsp_Din_Scr,Ord_Des,Ord_Din,Aqp_Idx) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  ses_bsp_mat_sum %>% print(n=base::nrow(ses_bsp_mat_sum))

  # Miss Cases::
  #
  ses_bsp_mis_tib <- ses_bsp_inn_tib %>% 
    dplyr::mutate(Pos_Dif=Bsp_Pos-start) %>% 
    dplyr::filter(Pos_Dif!=0) %>% 
    dplyr::select(Bsp_Din_Scr,Ord_Key,Ses_Cgn,Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx,
                  Ord_Prb,Pos_Dif,Bsp_Chr,seqnames)
  ses_bsp_mis_sum <- ses_bsp_mis_tib %>% 
    # dplyr::group_by(Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx) %>% 
    dplyr::group_by(Bsp_Tag,Bsp_Din_Scr,Ord_Des,Ord_Din,Aqp_Idx) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  ses_bsp_mis_sum %>% print(n=base::nrow(ses_bsp_mis_sum))
  
  #
  # BSP:: Missing Probes???
  #
  aqp_bsp_tib %>% 
    dplyr::select(Ord_Key,Ord_Des,Ord_Din,Ord_Prb,Bsp_Chr,Bsp_Pos,Aqp_Idx) %>%
    dplyr::anti_join(ses_pos_tab, by=c("Ord_Prb","Ord_Des","Ord_Din")) %>% 
    dplyr::add_count(Ord_Prb, name="Aln_Cnt") %>%
    dplyr::group_by(Aqp_Idx,Aln_Cnt,Ord_Des,Ord_Din) %>%
    dplyr::summarise(Count=n(), .groups="drop") %>% print(n=10000)
  
  #
  # Compare seq_cgn_tib vs. ses_man_tib
  #
  aqp_bsp_tib %>% dplyr::filter(Ord_Key %in% seq_cgn_tib$Address) %>%
    dplyr::select(Bsp_Seq,Aln_Prb)

  aqp_bsp_tib %>% dplyr::filter(!Bsp_Seq %in% seq_cgn_tib$Aln_Prb)
  aqp_bsp_tib %>% dplyr::filter(!Aln_Prb %in% seq_cgn_tib$Aln_Prb)
  
  aqp_bsp_tib %>% dplyr::select(Ord_Des:Aln_Key) %>% dplyr::arrange(Aln_Prb)
  seq_cgn_tib %>% dplyr::select(Address:Aln_Prb)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 3.4.3 Bind BSMAP & Seq-Match into Table::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


#
# Compare seq_cgn_tib vs. imp_des_tib for TB/CO strand names
#

#
# Compare aqp_ann_tib vs. ses_man_tib
#   NOTE: Done below doesn't exist yet...
#


#
# improbe matching testing code below::
#
if (FALSE) {
  
  if (FALSE) {
    
    imp_ext_tib <- imp_des_tib %>% 
      # dplyr::select(Seq_ID,Probe_Seq_U,Probe_Seq_M) %>% 
      dplyr::mutate(Aln_U49=stringr::str_sub(Probe_Seq_U, 1,49), 
                    Aln_M49=stringr::str_sub(Probe_Seq_M, 1,49))

    aqp_cgn_vec <- aqp_add_tib %>% dplyr::mutate(Ord_Cgn=stringr::str_remove(Ord_Key, "-.*$")) %>% dplyr::pull(Ord_Cgn) %>% unique()

    #
    # Matching only "Mat_Prb" Plus
    #
    aqp_add_des_tib6 <- dplyr::bind_rows(
      
      # UC::
      aqp_add_tib %>% 
        dplyr::filter(Ord_Des=="U") %>% 
        dplyr::mutate(Mat_Prb=Aln_Prb) %>%
        # head(n=2) %>% 
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=Probe_Seq_U), 
                          by=c("Mat_Prb")),
      
      # UO::
      aqp_add_tib %>% 
        dplyr::filter(Ord_Des=="U") %>% 
        dplyr::mutate(Mat_Prb=stringr::str_sub(Aln_Prb, 2)) %>%
        # head(n=2) %>% 
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=stringr::str_sub(Aln_U49,1)),
                          by=c("Mat_Prb")),
      
      # MC::
      aqp_add_tib %>% 
        dplyr::filter(Ord_Des=="M") %>% 
        dplyr::mutate(Mat_Prb=Aln_Prb) %>%
        # head(n=2) %>% 
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=Probe_Seq_M), 
                          by=c("Mat_Prb")),
      
      # MO::
      aqp_add_tib %>% 
        dplyr::filter(Ord_Des=="M") %>% 
        dplyr::mutate(Mat_Prb=stringr::str_sub(Aln_Prb, 2)) %>%
        # head(n=2) %>% 
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=stringr::str_sub(Aln_U49,1)),
                          by=c("Mat_Prb")),
      
      # 2C::
      aqp_add_tib %>% 
        dplyr::filter(Ord_Des=="2") %>% 
        dplyr::mutate(Mat_Prb=Aln_P49) %>%
        # head(n=2) %>% 
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=Aln_U49), 
                          by=c("Mat_Prb")),
      
      # 2O::
      aqp_add_tib %>%
        dplyr::filter(Ord_Des=="2") %>%
        dplyr::mutate(Mat_Prb=stringr::str_sub(Aln_Prb, 1,49)) %>%
        # head(n=2) %>%
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=Aln_U49),
                          by=c("Mat_Prb"))
      
    )
    aqp_add_tib$Ord_Prb %>% unique() %>% length()
    aqp_add_des_tib6$Ord_Prb %>% unique() %>% length()
    aqp_add_tib %>% dplyr::filter( Ord_Prb %in% aqp_add_des_tib6$Ord_Prb) %>% dplyr::distinct(Ord_Prb) %>% base::nrow()
    aqp_add_tib %>% dplyr::filter(!Ord_Prb %in% aqp_add_des_tib6$Ord_Prb) %>% dplyr::distinct(Ord_Prb) %>% base::nrow()
    
    aqp_add_tib %>% dplyr::filter(!Ord_Prb %in% aqp_add_des_tib6$Ord_Prb) %>%
      dplyr::mutate(Ord_Srd=stringr::str_remove(Ord_Key, "^.*-")) %>%
      dplyr::group_by(Ord_Idx,Ord_Srd,Ord_Des) %>%
      dplyr::summarise(Count=n(), .groups="drop") %>% print(n=1000)
    
    aqp_add_tib %>% dplyr::filter(Ord_Prb %in% aqp_add_des_tib6$Ord_Prb) %>%
      dplyr::mutate(Ord_Srd=stringr::str_remove(Ord_Key, "^.*-")) %>%
      dplyr::group_by(Ord_Idx,Ord_Srd,Ord_Des) %>%
      dplyr::summarise(Count=n(), .groups="drop") %>% print(n=1000)
    
    # aqp_add_des_tib6 %>% dplyr::mutate(Ord_Cgn=stringr::str_remove(Ord_Key, "-.*$")) %>% dplyr::pull(Ord_Cgn) %>% unique()
    aqp_add_des_tib6 %>% dplyr::mutate(Ord_Cgn=stringr::str_remove(Ord_Key, "-.*$")) %>%
      dplyr::filter(!Ord_Cgn %in% aqp_cgn_vec)
    
    #
    #
    # CONCLUSION: WE GET ALL THE CGN's!!!
    #
    #
  }
}

#
# All srd probe extraction::
#  - Seems to be working, will validate later...
#  - Not really needed at this moment, but useful...
#
if (FALSE) {
  
  all_srd_prb_tib <- aqp_bsp_tib %>% head(n=10) %>%
    bed_to_prbs(fas=run$gen_ref_fas, 
                din="Ord_Din", cgn="Aln_Key",
                chr="Bsp_Chr", pos="Bsp_Pos",
                nrec = 1,
                verbose=opt$verbose+10, tt=pTracker)
  
}

if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Binding all annotation...{RET}"))

ann_int_tib <- ann_int_list %>% 
  dplyr::bind_rows() %>%
  dplyr::arrange(IlmnID)

ann_int_sum <- ann_int_tib %>% 
  dplyr::group_by(source,class,tissue) %>% 
  dplyr::summarise(Count=n(), .groups="drop")
ann_int_sum %>% print(n=base::nrow(ann_int_sum))

if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Writing all annotation(CSV)={run$ann_int_csv}...{RET}"))
readr::write_csv(ann_int_tib, run$ann_int_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
