
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
par$prgmTag <- 'workhorse_main_db151'
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
opt$ords <- NULL
opt$mats <- NULL
opt$aqps <- NULL
opt$pqcs <- NULL
opt$mans <- NULL
opt$vcfs <- NULL

opt$beds <- NULL

# Manufacturing Info:: Required
opt$bpns <- NULL
opt$aqpn <- NULL
opt$pqcn <- NULL

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
  par$local_runType <- 'GRCm10'
  par$local_runType <- 'Chicago'
  par$local_runType <- 'GSA'
  par$local_runType <- 'dbSNP'
  
  if (par$local_runType=='dbSNP') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'GSA'
    opt$version  <- 'A1'
    opt$Species  <- "Human"
    
    opt$idat   <- NULL
    opt$ords   <- NULL
    opt$mats   <- NULL
    opt$aqps   <- NULL
    opt$pqcs   <- NULL
    
    opt$bpns   <- NULL
    opt$aqpn   <- NULL
    opt$pqcn   <- NULL
    
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
    
  } else if (par$local_runType=='GSA') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'GSA'
    opt$version  <- 'A1'
    opt$version  <- 'A2'
    opt$Species  <- "Human"
    
    opt$idat   <- NULL
    opt$ords   <- NULL
    opt$mats   <- NULL
    opt$aqps   <- NULL
    opt$pqcs   <- NULL
    
    opt$bpns   <- NULL
    opt$aqpn   <- NULL
    opt$pqcn   <- NULL
    
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
    
    opt$Species  <- "Human"
    
    opt$idat   <- NULL
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/UnivChicago/latest')
    
    opt$ords <- paste(
      file.path(par$aqpDir, 'UofChicago-A_A_Array-CpG-order-FINAL.csv'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20504790_probes.match.tsv'),
      sep=',')
    
    opt$aqps <- NULL
    
    opt$pqcs <- paste(
      file.path(par$aqpDir, '329922X374054_A_ProductQC.txt'),
      sep=',')
    
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
    
    opt$pqcs <- NULL
    
    opt$bpns <- paste(1, sep=",")
    opt$aqpn <- paste(1, sep=",")
    opt$pqcn <- NULL
    
  } else if (par$local_runType=='GRCm10') {
    opt$platform <- 'LEGX'
    opt$version  <- 'C0'
    opt$version  <- 'C30'
    opt$version  <- 'C31'
    opt$version  <- 'C32'
    
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
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'AQP_Copy/BS0032527-AQP.txt.gz'),
      file.path(par$aqpDir, 'AQP_Copy/BS0032533-AQP.txt.gz'),
      file.path(par$aqpDir, 'AQP_Copy/BS0032545-AQP.txt.gz'),
      file.path(par$aqpDir, 'AQP_Copy/BS0032636-AQP.txt.gz'),
      sep=',')
    
    opt$pqcs <- paste(
      file.path(par$aqpDir, 'PQC/20042400_A_ProductQC.txt.gz'),
      sep=',')
    
    opt$bpns <- paste(1,2,2,3, sep=",")
    opt$aqpn <- paste(1,1,2,1, sep=",")
    opt$pqcn <- paste(1, sep=",")
    
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
    opt$pqcn <- NULL
    
    opt$pqcs <- NULL
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
    make_option(c("--pqcs"), type="character", default=opt$pqcs, 
                help="PQC file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mans"), type="character", default=opt$mans,
                help="Manifest CSV file(s) (comma seperated) [default= %default]", metavar="character"),
    
    # Not fully supported yet::
    make_option(c("--vcfs"), type="character", default=opt$vcfs, 
                help="Target Design VCF file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--beds"), type="character", default=opt$beds, 
                help="Target Design Coordinate BED file(s) (comma seperated) [default= %default]", metavar="character"),
    
    # Manufacturing Info:: Required
    make_option(c("--bpns"), type="character", default=opt$bpns, 
                help="Bead Pool Numbers (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqpn"), type="character", default=opt$aqpn, 
                help="AQP Numbers (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--pqcn"), type="character", default=opt$pqcn, 
                help="PQC Numbers; should be 1 (comma seperated) [default= %default]", metavar="character"),
    
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
image_ssh <- "run_improbe.sh"
image_str <- glue::glue("{image_key}.{image_ver}")

# TBD:: Make this more general Manifest Control Defaults::
#
if (is.null(opt$ctls)) {
  opt$ctls <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
}

# TBD:: Add method validation for this step
#
if (is.null(opt$ords) && is.null(opt$mans) && is.null(opt$beds)) {
  stop(glue::glue("{RET}[{par$prgmTag}]: Must provide order, or manifest or coordinate files!!!.{RET}{RET}"))
}

if (!is.null(opt$ords)) {
  ords_vec <- NULL
  mats_vec <- NULL
  aqps_vec <- NULL
  pqcs_vec <- NULL
  if (!is.null(opt$ords)) ords_vec <- opt$ords %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  if (!is.null(opt$mats)) mats_vec <- opt$mats %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  if (!is.null(opt$aqps)) aqps_vec <- opt$aqps %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  if (!is.null(opt$pqcs)) pqcs_vec <- opt$pqcs %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  
  ords_len <- length(ords_vec)
  mats_len <- length(mats_vec)
  stopifnot(ords_len>0)
  stopifnot(mats_len>0)
  stopifnot(ords_len==mats_len)
  
  bpns_vec <- NULL
  aqpn_vec <- NULL
  pqcn_vec <- NULL
  if (!is.null(opt$bpns)) bpns_vec <- opt$bpns %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  if (!is.null(opt$aqpn)) aqpn_vec <- opt$aqpn %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  if (!is.null(opt$pqcn)) pqcn_vec <- opt$pqcn %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  
  # TBD:: Should throw an error if pqcn_vec > 1...
  #  Maybe not, I guess there could be more than on PQC...
  
  if (is.null(bpns_vec)) bpns_vec <- c(1:ords_len)
  if (is.null(aqpn_vec)) aqpn_vec <- c(1:ords_len)
  # if (is.null(pqcn_vec)) pqcn_vec <- c(1:ords_len)
  
  man_info_tib <- 
    tibble::tibble(Dat_BPN=bpns_vec, Dat_AQP=aqpn_vec) %>%
    dplyr::mutate(Dat_IDX=dplyr::row_number()) %>% 
    dplyr::select(Dat_IDX, dplyr::everything()) %>% 
    utils::type.convert()
  
  ctls_vec <- NULL
  idat_vec <- NULL
  if (!is.null(opt$ctls)) ctls_vec <- opt$ctls %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  if (!is.null(opt$idat)) idat_vec <- opt$idat %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
}

mans_vec <- NULL
if (!is.null(opt$mans)) {
  if (!is.null(opt$mans)) mans_vec <- opt$mans %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
}

if (!is.null(opt$vcfs)) {
  
}
if (!is.null(opt$beds)) {
  
}

cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

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

# Defined Run Time:: Intermediate Files
run$add_pas_csv <- file.path(run$addDir, paste(opt$runName,"address-pass.csv.gz", sep="."))
run$man_fun_csv <- file.path(run$manDir, paste(opt$runName,"functional-sesame.manifest.csv.gz", sep="."))
run$add_fin_csv <- file.path(run$manDir, paste(opt$runName,"final-address-aligned.csv.gz", sep="."))

run$man_mas_csv <- file.path(run$manDir, paste(opt$runName,"final-master.manifest.csv.gz", sep="."))
run$man_ses_csv <- file.path(run$manDir, paste(opt$runName,"final-sesame.manifest.csv.gz", sep="."))
run$man_gsm_csv <- file.path(run$manDir, paste(opt$runName,"final-GenomeStudio.manifest.csv.gz", sep="."))

run$add_prb_fas <- file.path(run$fasDir, paste(opt$runName, "aln-seq.fa.gz",  sep='.') )
run$add_dat_csv <- file.path(run$addDir, paste(opt$runName, "add_dat.csv.gz", sep='.') )

run$add_u49_tsv <- file.path(run$intDir, paste(opt$runName, "map-u49.tsv", sep='.') )
run$add_m49_tsv <- file.path(run$intDir, paste(opt$runName, "map-u50.tsv", sep='.') )

run$int_u49_tsv <- file.path(run$intDir, paste(opt$runName, "int-u49.tsv.gz", sep='.') )
run$int_m49_tsv <- file.path(run$intDir, paste(opt$runName, "int-m49.tsv.gz", sep='.') )
run$int_seq_tsv <- file.path(run$intDir, paste(opt$runName, "int-seq-imp.tsv.gz", sep='.') )

run$add_prb_bsp  <- file.path(run$alnDir, paste(opt$runName, "bsp",  sep='.') )
run$add_prb_bspz <- paste(run$add_prb_bsp, 'tsv.gz', sep='.')

run$imp_inp_tsv  <- file.path(run$desDir, paste(opt$runName, 'improbe-inputs.tsv.gz', sep='.') )
run$imp_des_tsv  <- file.path(run$desDir, paste(opt$runName, 'improbe-design.tsv.gz', sep='.') )

run$add_pas_bsp_csv <- file.path(run$alnDir, paste(opt$runName, "add_pas_bsp.csv.gz",  sep='.') )
run$add_pas_grs_rds <- file.path(run$alnDir, paste(opt$runName, "add_pas_grs.rds",  sep='.') )

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Defining Run Time Files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    0.0 Pre-processing:: Three Groups
#
#  0. Validation idats
#  1. Original Order File Genomic Positions
#  2. Canonical CGN Preference
#     - Need to pull canonical sources to see what's old and new!!!
#  3. Target Genome BED File
#
#  4. Load Reference FASTA Genome(s)
#     - Ref (ACTG)
#     - SNP (IUPAC)
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #




# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                0.4 Extract Forward Sequence from Genome::
#                        Reference/IUPAC Reference
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#            0.4.1 Pre-processing:: Run Time:: Ref Alignment Genome
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run$gen_ref_fas <- 
  file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
            paste0(opt$genBuild,".genome.fa.gz"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          VALIDATION CODE BELOW::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  #
  # Load Validation/Testing Manifest::
  #
  tar_chr <- "chr1"
  
  # Load Manifest for comparison::
  #
  opt$mans   <- paste(
    # file.path(opt$manDir, "genotyping/GSA-24v2-0_A1.csv.gz"),
    # file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz"),
    file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz"),
    # file.path(opt$manDir, "methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz"),
    sep=',')
  
  mans_vec <- NULL
  if (!is.null(opt$mans)) mans_vec <- opt$mans %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  
  mans_dat_tib <- mans_vec %>%
    lapply(loadManifestGenomeStudio,
           addSource=TRUE, normalize=TRUE, retType="man", 
           verbose=opt$verbose+10,tt=pTracker) %>%
    dplyr::bind_rows()
  
  # Extract Only Chr1 hits:: top 10
  #
  cur_man_tib <- mans_dat_tib %>% 
    dplyr::filter(Chromosome==tar_chr) %>%
    dplyr::rename(Seq_ID=IlmnID) %>% head(n=10)
  
  #
  # Validate BSC Conversions::
  #
  run$ref_FCU_fas <- run$gen_ref_fas %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$") %>% paste("FCU.fa.gz", sep=".")
  run$ref_FCM_fas <- run$gen_ref_fas %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$") %>% paste("FCM.fa.gz", sep=".")
  run$ref_FCD_fas <- run$gen_ref_fas %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$") %>% paste("FCD.fa.gz", sep=".")

  run$ref_RCU_fas <- run$gen_ref_fas %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$") %>% paste("RCU.fa.gz", sep=".")
  run$ref_RCM_fas <- run$gen_ref_fas %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$") %>% paste("RCM.fa.gz", sep=".")
  run$ref_RCD_fas <- run$gen_ref_fas %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$") %>% paste("RCD.fa.gz", sep=".")
  
  ref_FCU_seq <- Biostrings::readDNAStringSet(filepath = run$ref_FCU_fas, format = "fasta", nrec = 2)
  ref_FCM_seq <- Biostrings::readDNAStringSet(filepath = run$ref_FCM_fas, format = "fasta", nrec = 2)
  ref_FCD_seq <- Biostrings::readDNAStringSet(filepath = run$ref_FCD_fas, format = "fasta", nrec = 2)

  ref_RCU_seq <- Biostrings::readDNAStringSet(filepath = run$ref_RCU_fas, format = "fasta", nrec = 2)
  ref_RCM_seq <- Biostrings::readDNAStringSet(filepath = run$ref_RCM_fas, format = "fasta", nrec = 2)
  ref_RCD_seq <- Biostrings::readDNAStringSet(filepath = run$ref_RCD_fas, format = "fasta", nrec = 2)
  
  # ref_FCU_seq[19] %>% stringr::str_sub(5236016 - 60, 5236016 + 61) %>% addBrac()
  
  
  # Extract top 10 Forward_Sequences::
  #     off_set <- 60
  off_set <- 60

  ref_fwd_tib <- tibble::tibble(
    Seq_ID  = cur_man_tib$Seq_ID,
    Forward_Sequence = stringr::str_sub(
      ref_fwd_seq[1], 
      cur_man_tib$Coordinate - off_set, 
      cur_man_tib$Coordinate + off_set + 1),
    PRB1=stringr::str_sub(Forward_Sequence, off_set+1, off_set+1+off_set-11) %>% revCmp(),
    PRB2=stringr::str_sub(Forward_Sequence, off_set+2, off_set+2+off_set-11) %>% revCmp(),
    # Probe_Type = cur_man_tib$Probe_Type,
    Strand_FR=cur_man_tib$Strand_FR,
    Inf_Design=cur_man_tib$Infinium_Design
  ) %>%
    dplyr::mutate(
      Forward_Sequence=Forward_Sequence %>% addBrac(),
    )
  
  cur_FCU_tib <- tibble::tibble(
    Seq_ID  = cur_man_tib$Seq_ID,
    Forward_Sequence = stringr::str_sub(
      ref_FCU_seq[1], 
      cur_man_tib$Coordinate - off_set, 
      cur_man_tib$Coordinate + off_set + 1),
    PRB1=stringr::str_sub(Forward_Sequence, off_set+1, off_set+1+off_set-11) %>% revCmp(),
    PRB2=stringr::str_sub(Forward_Sequence, off_set+2, off_set+2+off_set-11) %>% revCmp(),
    # Probe_Type = cur_man_tib$Probe_Type,
    Strand_FR=cur_man_tib$Strand_FR,
    Inf_Design=cur_man_tib$Infinium_Design
  ) %>%
    dplyr::mutate(
      Forward_Sequence=Forward_Sequence %>% addBrac(),
    )
  
  cur_FCM_tib <- tibble::tibble(
    Seq_ID  = cur_man_tib$Seq_ID,
    Forward_Sequence = stringr::str_sub(
      ref_FCM_seq[1], 
      cur_man_tib$Coordinate - off_set, 
      cur_man_tib$Coordinate + off_set + 1),
    PRB1=stringr::str_sub(Forward_Sequence, off_set+1, off_set+1+off_set-11) %>% revCmp(),
    PRB2=stringr::str_sub(Forward_Sequence, off_set+2, off_set+2+off_set-11) %>% revCmp(),
    # Probe_Type = cur_man_tib$Probe_Type,
    Strand_FR=cur_man_tib$Strand_FR,
    Inf_Design=cur_man_tib$Infinium_Design
  ) %>%
    dplyr::mutate(
      Forward_Sequence=Forward_Sequence %>% addBrac(),
    )
  
  cur_FCD_tib <- tibble::tibble(
    Seq_ID  = cur_man_tib$Seq_ID,
    Forward_Sequence = stringr::str_sub(
      ref_FCD_seq[1], 
      cur_man_tib$Coordinate - off_set, 
      cur_man_tib$Coordinate + off_set + 1),
    PRB1=stringr::str_sub(Forward_Sequence, off_set+1, off_set+1+off_set-11) %>% revCmp(),
    PRB2=stringr::str_sub(Forward_Sequence, off_set+2, off_set+2+off_set-11) %>% revCmp(),
    # Probe_Type = cur_man_tib$Probe_Type,
    Strand_FR=cur_man_tib$Strand_FR,
    Inf_Design=cur_man_tib$Infinium_Design
  ) %>%
    dplyr::mutate(
      Forward_Sequence=Forward_Sequence %>% addBrac(),
    )
  
  print(ref_fwd_tib)
  print(cur_FCU_tib)
  print(cur_FCM_tib)
  print(cur_FCD_tib)
  
}

#
# Validation for SNPs on the 450k::
#
if (FALSE) {
  
  if (FALSE) {
    epic_snp_ids <- readr::read_csv("/Users/bretbarnes/Documents/data/manifests/raw/manifests/methylation/rs-repair/450k.rs.txt", col_names = c("Seq_ID"))
    epic_snp_tib <- vcf_tib %>% dplyr::filter(Seq_ID %in% epic_snp_ids$Seq_ID)
    
    off_set <- 60
    tmp_seq <- stringr::str_sub(new_seq, vcf_chrs[[1]]$Coordinate - off_set, vcf_chrs[[1]]$Coordinate + off_set + 1) %>% addBrac()
    tmp_seq %>% print()
    tmp_seq %>% stringr::str_length() %>% print()
    
    # 110359535
    snp_tib <- tibble::tibble(
      Seq_ID  = vcf_chrs[[1]]$Seq_ID,
      Forward_Sequence = stringr::str_sub(
        new_seq, 
        vcf_chrs[[1]]$Coordinate - off_set, 
        vcf_chrs[[1]]$Coordinate + off_set + 1) %>% addBrac(),
      Probe_Type = "rs"
    )
    
    snp_des_tib <- desSeq_to_prbs(
      tib = snp_tib, 
      ids_key="Seq_ID", seq_key="Forward_Sequence", prb_key="Probe_Type", 
      strsSR="FR", strsCO="CO", 
      addMatSeq=TRUE, parallel=TRUE, verbose=opt$verbose, tt=pTracker)
  }
  
  opt$mans   <- paste(
    # file.path(opt$manDir, "genotyping/GSA-24v2-0_A1.csv.gz"),
    # file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz"),
    file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz"),
    # file.path(opt$manDir, "methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz"),
    sep=',')
  
  mans_vec <- NULL
  if (!is.null(opt$mans)) mans_vec <- opt$mans %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  
  mans_dat_tib <- mans_vec %>%
    lapply(loadManifestGenomeStudio,
           addSource=TRUE, normalize=TRUE, retType="man", 
           verbose=opt$verbose+10,tt=pTracker) %>%
    dplyr::bind_rows()
  
  epic_rs_tib <- mans_dat_tib %>% 
    dplyr::filter(Probe_Type=="rs") %>%
    dplyr::rename(Seq_ID=IlmnID)
  
  mat_snp_tib <- dplyr::bind_rows(
    epic_rs_tib %>% dplyr::inner_join(snp_des_tib, by=c("Seq_ID","Probe_Type")) %>%
      dplyr::filter(AlleleA_ProbeSeq==PRB1_M_MAT),
    epic_rs_tib %>% dplyr::inner_join(snp_des_tib, by=c("Seq_ID","Probe_Type")) %>%
      dplyr::filter(AlleleA_ProbeSeq==PRB2_D_MAT)
  )
  
}

if (FALSE) {
  #
  # LIFE GOT BETTER:: VCF Reading Scratch Space::
  #
  opt$dbSNP   <- "dbSNP-151"
  run$snp_vcf <- file.path(opt$annDir, "dbSNP", opt$dbSNP, opt$genBuild) %>% list.files(pattern = ".vcf.gz$", full.names=TRUE) %>% head(n=1)
  run$snp_csv <- run$snp_vcf %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".vcf$") %>% paste0(".snps.csv.gz")
  run$snp_fas <- run$gen_ref_fas %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$") %>% paste(opt$dbSNP,"iupac.fa.gz", sep=".")
  
  vcf_tib <- load_dbSNP_vcf(vcf=run$snp_vcf, file=run$snp_csv,
                            fresh=opt$fresh,
                            verbose=opt$verbose+10, tt=pTracker)
  
  if (FALSE) {
    
    # Summary by chromosome coverage::
    vcf_sum <- vcf_tib %>% 
      dplyr::group_by(Chromosome) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    vcf_sum %>% print(n=base::nrow(vcf_sum))
    
    # Summary by SNP Allele coverage::
    vcf_sum <- vcf_tib %>% 
      dplyr::mutate(Probe_Type=stringr::str_sub(Seq_ID, 1,2)) %>% 
      dplyr::group_by(AlleleA_Str,AlleleB_Str2,AlleleC_Iup) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    vcf_sum %>% print(n=base::nrow(vcf_sum))
    
  }
  
  #
  # TBD:: Turn this into a function...
  # TBD:: Add perl MUD call..
  # Build dbSNP/IUPAC Genome
  #
  vcf_chrs <- vcf_tib %>% split(.$Chromosome)
  chr_maps <- ref_fwd_tab %>% split(.[["Chrom_Char"]])
  
  snp_seqs <- NULL
  snp_seqs <- Biostrings::DNAStringSet()
  chr_seqs <- NULL
  
  # for (chr_str in names(vcf_chrs)) {
  for (chr_str in names(chr_maps)) {
    
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Mutating chromosome={chr_str}...{RET}"))
    
    chr_idx <- chr_maps[[chr_str]] %>% pull(Idx) %>% as.integer()
    
    # This should never happen...
    if (is.null(ref_fwd_seq[chr_idx])) {
      if (opt$verbose>=1)
        cat(glue::glue("[{par$prgmTag}]: No Genomic Sequence data for chromosome={chr_str}. Skipping...{RET}"))
      next
    }
    chr_seqs <- c(chr_seqs, chr_str)
    
    snp_seq <- NULL
    if (is.null(vcf_chrs[[chr_str]])) {
      snp_seq <- ref_fwd_seq[[chr_idx]]
    } else {
      val_len <- length(vcf_chrs[[chr_str]]$AlleleC_Iup)
      pos_len <- length(vcf_chrs[[chr_str]]$Coordinate)
      
      snp_seq <- mutate_chrom_seq(seq=ref_fwd_seq[[chr_idx]], 
                                  pos=vcf_chrs[[chr_str]]$Coordinate, 
                                  val=vcf_chrs[[chr_str]]$AlleleC_Iup, 
                                  verbose=opt$verbose, tt=pTracker)
    }
    
    snp_seqs <- snp_seqs %>%
      BiocGenerics::append(Biostrings::DNAStringSet(snp_seq) )
    
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Done. Mutating chromosome={chr_str} ",
                     "total positions/values={pos_len}/{val_len}.{RET}{RET}"))
    
    # if (chr_idx>=3) break
    # if (chr_str=="chr11") break
  }
  
  # Add names back, but need to track names above to make sure nothing was skipped!!!
  #
  names(snp_seqs) <- chr_seqs
  
  # Sort names before printing::
  #
  Biostrings::writeXStringSet(
    x=snp_seqs[ref_fwd_tab$Chrom_Char[ref_fwd_tab$Chrom_Char %in% names(snp_seqs)]], 
    filepath=run$snp_fas, 
    append=FALSE, compress=TRUE)
  
}

#
#
# STOPPPED HERE for preperation. All of this should exist now.
#
#

















# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     0.0 Pre-processing:: IDAT Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_idat_tib <- NULL
if (!is.null(opt$idat)) {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading idats from dir={opt$idat}.{RET}"))
  
  grn_idat_dat <- 
    loadIdat(prefix=file.path(opt$idat, paste(base::basename(opt$idat),"R01C01", sep="_")), 
             col="Grn", gzip=TRUE, verbose=opt$verbose, tt=pTracker )
  
  add_idat_tib <- 
    grn_idat_dat$Quants %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="Address") %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(Address=as.integer(Address)) %>%
    dplyr::distinct(Address)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#         0.1 Pre-processing:: Original Order File Genomic Positions
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

org_des_tib <- NULL
add_org_tib <- NULL
if (!is.null(opt$org_des_tsv) && file.exists(opt$org_des_tsv)) {
  
  org_des_col <-
    cols(
      Org_Cgn      = col_integer(),
      Org_Chr      = col_character(),
      Org_Pos      = col_integer(),
      Org_M49P     = col_character(),
      Org_M49N     = col_character(),
      
      Org_FR       = col_character(),
      Org_TB       = col_character(),
      Org_CO       = col_character(),
      
      Org_NXB      = col_character(),
      Org_U49P     = col_character(),
      Org_U49N     = col_character()
    )
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading org_des_tsv={opt$org_des_tsv}...{RET}"))
  
  org_des_tib <- suppressMessages(suppressWarnings(
    readr::read_tsv(opt$org_des_tsv, 
                    col_names=names(org_des_col$cols),
                    col_types=org_des_col) )) %>%
    dplyr::mutate(Org_Chr=paste0("chr",stringr::str_remove(Org_Chr, "^chr")) ) #,Org_Srd3=paste0(Org_TB,Org_CO,Org_NXB)
  
  add_org_tib <- dplyr::bind_rows(
    dplyr::select(org_des_tib, Org_Cgn,Org_Chr,Org_Pos,Org_FR,Org_TB,Org_CO,Org_NXB,Org_U49P,Org_U49N) %>%
      dplyr::rename(Org_49P=Org_U49P, Org_49N=Org_U49N) %>%
      dplyr::mutate(Org_49P=Org_49P %>% stringr::str_sub(1,49),
                    Org_49N=NA_character_,
                    Org_Des="2"),
    dplyr::select(org_des_tib, Org_Cgn,Org_Chr,Org_Pos,Org_FR,Org_TB,Org_CO,Org_NXB,Org_U49P,Org_U49N) %>%
      dplyr::rename(Org_49P=Org_U49P, Org_49N=Org_U49N) %>%
      dplyr::mutate(Org_Des="U"),
    dplyr::select(org_des_tib, Org_Cgn,Org_Chr,Org_Pos,Org_FR,Org_TB,Org_CO,Org_NXB,Org_M49P,Org_M49N) %>%
      dplyr::rename(Org_49P=Org_M49P, Org_49N=Org_M49N) %>%
      dplyr::mutate(Org_Des="M")
  ) %>% 
    dplyr::mutate(Org_Mis_Scr=0, Org_Mis_Scr=as.integer(Org_Mis_Scr)) %>%
    dplyr::arrange(Org_Cgn,Org_Des) %>% 
    dplyr::select(Org_Cgn,Org_Des, dplyr::everything())
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Done.{RET}{RET}"))  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                0.2 Pre-processing:: Canonical CGN Preference
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run$can_cgn_top_csv <- file.path(par$datDir, "manifest/cgnDB/canonical-assignment.cgn-top-grp.csv.gz")

can_cgn_top_tib <- NULL
if (file.exists(run$can_cgn_top_csv)) {
  can_cgn_top_tib <- 
    suppressMessages(suppressWarnings( readr::read_csv(run$can_cgn_top_csv) )) %>% 
    dplyr::mutate(CGN=as.integer(CGN)) %>%
    dplyr::rename(Can_Cgn=CGN, Can_Top=TOP, Can_Src=SRC) %>%
    clean_tibble()
}
  
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                0.3 Pre-processing:: Target Genome BED File
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cgn_bed_dir <- "/Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-input"
cgn_bed_tsv <- file.path(cgn_bed_dir, paste(opt$genBuild, "cgn.bed.gz", sep="."))

cgn_bed_tib <- NULL
if (file.exists(cgn_bed_tsv)) {
  cgn_bed_col <-
    cols(
      Imp_Chr      = col_character(),
      Imp_Beg      = col_integer(),
      Imp_End      = col_integer(),
      Imp_Cgn      = col_integer(),
      Imp_Src      = col_character(),
      Imp_Top_Srd  = col_character()
    )
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading cgn_bed_tsv={cgn_bed_tsv}...{RET}"))
  
  cgn_bed_tib <- suppressMessages(suppressWarnings( 
    readr::read_tsv(cgn_bed_tsv, 
                    col_names=names(cgn_bed_col$cols),
                    col_types=cgn_bed_col) )) %>% 
    dplyr::select(-Imp_Beg,-Imp_Src) %>% 
    dplyr::rename(Imp_Pos=Imp_End) %>%
    dplyr::left_join(can_cgn_top_tib, by=c("Imp_Cgn"="Can_Cgn")) %>%
    dplyr::mutate(
      Can_Mis_Scr=dplyr::case_when(
        !is.na(Can_Src) ~ 0,
        TRUE ~ 1) %>% as.integer()
    )
  
  cgn_bed_grs <- 
    GRanges(
      seqnames = cgn_bed_tib$Imp_Chr,
      strand   = Rle(cgn_bed_tib$Imp_Top_Srd),
      cgn      = cgn_bed_tib$Imp_Cgn,
      
      # top_seq  = cgn_bed_tib$Can_Top,
      # top_srd   = cgn_bed_tib$Imp_Top_Srd,
      # top_src  = cgn_bed_tib$Can_Src,
      
      IRanges(start=cgn_bed_tib$Imp_Pos, width = 1, 
              names=paste(cgn_bed_tib$Imp_Cgn, 
                          cgn_bed_tib$Imp_Chr, 
                          cgn_bed_tib$Imp_Pos, sep='_')
      )
    )
  
  # Remove can_cgn_top_tib from memory, or make a function...
  #
  # rm(can_cgn_top_tib)
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                X.X Load Candidate Target Regions & build GRS
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# 1. EWAS**
# 2. dbSNP-151
# 3. Manifest Validation
#

opt$dbSNP   <- "dbSNP-151"
run$snp_vcf <- file.path(opt$annDir, "dbSNP", opt$dbSNP, opt$genBuild) %>% list.files(pattern = ".vcf.gz$", full.names=TRUE) %>% head(n=1)
run$snp_csv <- run$snp_vcf %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".vcf$") %>% paste0(".snps.csv.gz")

vcf_tib <- load_dbSNP_vcf(vcf=run$snp_vcf, file=run$snp_csv,
                          fresh=opt$fresh,
                          verbose=opt$verbose+10, tt=pTracker)

if (FALSE) {
  
  new_cpg_dir <- '/Users/bretbarnes/Documents/data/CustomContent/EWAS/CORSIVs_and_others/Not_Already_Designed_Content'
  new_wat_csv <- file.path(new_cpg_dir, 'Waterland_CORSIVS/Waterland_CORSIV_Genomic_Coordinates_S3-Table_1.csv.gz')
  
  new_reg_tib <- suppressMessages(suppressWarnings(readr::read_csv(new_wat_csv) )) %>% 
    dplyr::select(CorSIV_ID,CpG.Count,UCSC.Coordinates,USCS_Coordinates_CoRSIV) %>%
    dplyr::rename(CpG_Count=CpG.Count) %>%
    tidyr::separate(UCSC.Coordinates, into=c("Chr","range"), sep=":") %>%
    tidyr::separate(range, into=c("Beg","End"), sep="-") %>%
    clean_tibble() %>%
    dplyr::group_by(CorSIV_ID) %>%
    dplyr::mutate(Rep_Idx=row_number()) %>%
    dplyr::ungroup() %>% 
    tidyr::unite(Reg_ID, CorSIV_ID,Rep_Idx, sep="-", remove=FALSE)
  
  new_reg_grs <- 
    GRanges(
      seqnames=Rle(new_reg_tib$Chr), # strand=Rle(new_reg_tib$strand),
      CpG_Count=new_reg_tib$CpG_Count, 
      USCS_Coordinates_CoRSIV=new_reg_tib$USCS_Coordinates_CoRSIV,
      IRanges(start=new_reg_tib$Beg, 
              end=new_reg_tib$End, 
              names=new_reg_tib$Reg_ID)
    )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                 X.X Intersect Target Regions with known CpGs
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # A tibble: 34,974 x 7
  new_imp_int_tib <- 
    intersect_GRS(new_reg_grs, cgn_bed_grs, can_key="Reg_ID", 
                  verbose=opt$verbose, tt=pTracker) %>%
    dplyr::distinct(chr,pos,end, .keep_all=TRUE) %>%
    dplyr::rename(Chromosome=chr, Coordinate=pos)
  
  # Region Summary Coverage::
  #  - TBD:: Add missing regsions, can be done later...
  #
  new_imp_int_reg_sum <- new_imp_int_tib %>% 
    dplyr::group_by(Reg_ID) %>%
    dplyr::summarise(Reg_Cnt=n(), .groups="drop")
  # new_imp_int_reg_sum %>% print(n=base::nrow(new_imp_int_reg_sum))
  
  new_imp_int_cnt_sum <- new_imp_int_reg_sum %>%
    dplyr::group_by(Reg_Cnt) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  new_imp_int_cnt_sum %>% print(n=base::nrow(new_imp_int_cnt_sum))
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              X.X Extract Target Forward Design Sequences 122-mers
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imp_dna_tib <- 
  fas_to_seq(
    tib = vcf_tib, # %>% head(n=10000), 
    fas = run$gen_ref_fas,
    file=run$imp_inp_tsv,
    name="Seq_ID",gen=opt$genBuild,
    chr1="Chromosome",pos="Coordinate",chr2="Chrom_Char",
    ext_seq="Fwd_Seq",des_seq="Des_Seq",imp_seq="Sequence",
    iupac = "AlleleC_Iup",
    nrec=2,
    add_flank=FALSE,
    verbose = opt$verbose, tt=pTracker)

imp_des_tib <- 
  improbe_docker(
    dir=run$desDir,
    file=run$imp_inp_tsv,
    name=opt$runName, 
    image=image_str, shell=image_ssh,level=3,
    verbose=opt$verbose, tt=pTracker)

# [done]: dna_to_template() 
# [done]: improbe_docker()
# filter???
# bsmap
#


#
# Validation of BSC Chromosomes via coordinates::
#

validate_bsc_strands <- FALSE
if (validate_bsc_strands) {

  imp_des_fr <- imp_des_tib %>% split(.$Strand_FR)
  
  ref_fas_list <- list()
  ref_fas_prbs <- list()
  
  for (fr in names(imp_des_fr)) {
    file_str_cu <- paste0(fr,"CU")
    file_str_cm <- paste0(fr,"CM")
    srd_str_vec <- c(file_str_cu,file_str_cm)
    
    ref_fas_list[[file_str_cu]] <- run$gen_ref_fas %>% 
      stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$") %>% 
      paste(file_str_cu,"fa.gz", sep=".")
    
    ref_fas_list[[file_str_cm]] <- run$gen_ref_fas %>% 
      stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$") %>% 
      paste(file_str_cm,"fa.gz", sep=".")
    
    for (srd in srd_str_vec) {
      if (opt$verbose>=1)
        cat(glue::glue("[{par$prgmTag}]: Extracting fr={fr}, srd={srd}...{RET}"))
      
      ref_fas_prbs[[srd]] <-
        fas_to_seq(
          tib = new_imp_int_tib, 
          fas = ref_fas_list[[srd]],
          name="cgn",gen=opt$genBuild,
          chr1="Chromosome",pos="Coordinate",chr2="Chrom_Char",
          srd=fr,
          ext_seq="Fwd_Seq",des_seq="Des_Seq",imp_seq="Sequence",
          # iupac = "QI_T",
          # nrec=2,
          add_flank=FALSE,
          verbose = opt$verbose, tt=pTracker
        )
      
      prb_tib <- NULL
      if (fr=="F") {
        prb_tib <- tibble::tibble(
          Seq_ID = ref_fas_prbs[[srd]]$cgn,
          Reg_ID = ref_fas_prbs[[srd]]$Reg_ID,
          Strand_FR = fr,
          Strand_CO = "C",
          
          Chromosome=ref_fas_prbs[[srd]]$Chromosome,
          Coordinate=ref_fas_prbs[[srd]]$Coordinate,
          
          Probe_Seq = paste0(
            ref_fas_prbs[[srd]]$up61,
            ref_fas_prbs[[srd]]$dn61,
            ref_fas_prbs[[srd]]$dn60,
            ref_fas_prbs[[srd]]$dn59,
            ref_fas_prbs[[srd]]$dn58 %>% stringr::str_sub(1,46)
          ) %>% revCmp()
        )
      } else if (fr=="R") {
        prb_tib <- tibble::tibble(
          Seq_ID = ref_fas_prbs[[srd]]$cgn,
          Reg_ID = ref_fas_prbs[[srd]]$Reg_ID,
          Strand_FR = fr,
          Strand_CO = "C",
          
          Chromosome=ref_fas_prbs[[srd]]$Chromosome,
          Coordinate=ref_fas_prbs[[srd]]$Coordinate,
          
          Probe_Seq = paste0(
            ref_fas_prbs[[srd]]$up58 %>% stringr::str_sub(13,58),
            ref_fas_prbs[[srd]]$up59,
            ref_fas_prbs[[srd]]$up60,
            ref_fas_prbs[[srd]]$up61,
            ref_fas_prbs[[srd]]$dn61
          ) %>% cmpl()
        )
      } else {
        stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported fr={fr}...{RET}"))
        return(NULL)
      }
      
      mat_tib <- NULL
      if (stringr::str_ends(srd,"U")) {
        mat_tib <- imp_des_tib %>% 
          dplyr::inner_join(
            prb_tib, 
            by=c("Seq_ID","Strand_FR","Strand_CO",
                 "Chromosome","Coordinate",
                 "Probe_Seq_U"="Probe_Seq")
          )
      } else {
        mat_tib <- imp_des_tib %>% 
          dplyr::inner_join(
            prb_tib, 
            by=c("Seq_ID","Strand_FR","Strand_CO",
                 "Chromosome","Coordinate",
                 "Probe_Seq_M"="Probe_Seq")
          )
      }
      
      if (opt$verbose>=6) {
        mat_tib %>% dplyr::select(Seq_ID, Strand_FR, Strand_CO,
                                  Probe_Seq_M) %>% print()
      }
      
      if (opt$verbose>=1)
        cat(glue::glue("[{par$prgmTag}]: Done. Extracting fr={fr}, srd={srd}.{RET}{RET}"))
      
      # break
    }
    
    # break
  }
  
}

#
# TBD: 
#   - Check Top Sequence from imp_des_tib vs. can_cgn_top_tib
#   - [Done] Substring Infinium I probes from U/M genomes
#




#
# Temp Fix for EWAS::
#
opt$min_scr <- 0.2
des_ord_tib <- imp_des_tib %>% 
  dplyr::mutate(Ord_Cgn=Seq_ID, 
                Seq_ID=paste(Seq_ID,Strand_TB,Strand_CO, sep="-"),
                Scr_Min=pmin(Scr_U,Scr_M)
  ) %>%
  dplyr::filter(Scr_Min>=opt$min_scr) %>%
  dplyr::distinct(Seq_ID,Probe_Seq_U,Probe_Seq_M, .keep_all=TRUE)

add_pas_tib <- dplyr::bind_rows(
  des_ord_tib %>% 
    dplyr::select(Seq_ID, Probe_Seq_U) %>% 
    dplyr::mutate(Ord_Des="U", Ord_Din="cg") %>%
    dplyr::rename(Ord_Prb=Probe_Seq_U),
  des_ord_tib %>% 
    dplyr::select(Seq_ID, Probe_Seq_M) %>% 
    dplyr::mutate(Ord_Des="M", Ord_Din="cg") %>%
    dplyr::rename(Ord_Prb=Probe_Seq_M)
) %>%
  dplyr::mutate(Dat_IDX=1)

#
# TBD:: Add Infinium II Probes???
#

if (TRUE) {

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #                    3.0 Alignment/Sequence Based CGN Mapping::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                  3.1 Format Fasta AND U49/U50 Probe Seqs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  stamp_vec <- c(stamp_vec,
                 run$add_dat_csv,run$add_prb_fas,
                 run$add_u49_tsv,run$add_m49_tsv)
  
  add_pas_fas_tib <- NULL
  if (opt$fresh || !valid_time_stamp(stamp_vec)) {
    
    add_pas_fas_tib <-
      add_to_fas(
        tib=add_pas_tib, 
        prb_key="Ord_Prb",
        add_key="Seq_ID", des_key="Ord_Des", type_key="Ord_Din",
        prb_fas=run$add_prb_fas, dat_csv=run$add_dat_csv,
        u49_tsv=run$add_u49_tsv, m49_tsv=run$add_m49_tsv,
        verbose=opt$verbose, tt=pTracker)
    
  } else {
    
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Loading add_dat_csv={run$add_dat_csv}...{RET}"))
    
    add_pas_fas_tib <- suppressMessages(suppressWarnings(
      readr::read_csv(run$add_dat_csv, guess_max=100000) )) %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
  }
  
  add_pas_fas_sum <- add_pas_fas_tib %>% 
    dplyr::group_by(Dat_IDX, Ord_Din, Ord_Des) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  add_pas_fas_sum %>% print(n=base::nrow(add_pas_fas_sum))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    3.2 Align All Probe Sequence:: BSMAP
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  stamp_vec <- c(stamp_vec, run$add_prb_bspz)
  if (opt$fresh || !valid_time_stamp(stamp_vec)) {
    
    run$add_prb_bspz <- 
      run_bsmap(
        exe=opt$bsmap_exe, 
        fas=run$add_prb_fas, 
        gen=run$gen_ref_fas,
        bsp=run$add_prb_bsp,
        opt=NULL, lan=NULL, run=TRUE,
        verbose=opt$verbose,tt=pTracker)
    
  } else {
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Using existing file add_prb_bspz={run$add_prb_bspz}...{RET}"))
  }
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  3.3 Join Address and Alignment Data:: BSMAP
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsp_imp_tib <- NULL
stamp_vec <- c(stamp_vec, run$add_pas_bsp_csv)
if (opt$fresh || !valid_time_stamp(stamp_vec)) {
  
  bsp_imp_tib <-
    join_bsmap(
      add=add_pas_fas_tib,
      bed=cgn_bed_tib, org=add_org_tib,
      file=run$add_prb_bspz,
      join_key="Aln_Key", 
      prb_des_key="Ord_Des", prb_din_key="Ord_Din",
      join_type="inner",
      sort=TRUE,
      verbose=opt$verbose,tt=pTracker)
  
  safe_write(bsp_imp_tib,"csv",run$add_pas_bsp_csv, funcTag=par$prgmTag, 
             verbose=opt$verbose)
  
} else {
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading add_pas_bsp_csv={run$add_pas_bsp_csv}...{RET}"))
  
  bsp_imp_tib <- suppressMessages(suppressWarnings(
    readr::read_csv(run$add_pas_bsp_csv, guess_max=100000) )) %>%
    utils::type.convert() %>% 
    dplyr::mutate(across(where(is.factor),  as.character) )
  
}


#
# Understand the CGN naming differences::
#
if (FALSE) {

  bsp_imp_tib %>% 
    dplyr::mutate(Aln_Cgn=stringr::str_remove(Aln_Key, "-.*$") %>% as.integer()) %>%
    dplyr::inner_join(can_cgn_top_tib, by=c("Imp_Cgn"="Can_Cgn"), suffix=c("_BSP","_CAN")) %>%
    dplyr::select(Imp_Cgn, Aln_Cgn, Can_Top_BSP, Can_Top_CAN)
  
  bsp_imp_tib %>% 
    dplyr::mutate(Aln_Cgn=stringr::str_remove(Aln_Key, "-.*$") %>% as.integer()) %>%
    dplyr::inner_join(can_cgn_top_tib, by=c("Aln_Cgn"="Can_Cgn"), suffix=c("_BSP","_CAN")) %>%
    dplyr::select(Imp_Cgn, Aln_Cgn, Can_Top_BSP, Can_Top_CAN)
  
}

# Same Counts below::
#   add_pas_fas_tib %>% dplyr::distinct(Ord_Prb)
#   bsp_imp_tib %>% dplyr::distinct(Aln_Prb)

bsp_sep_tib <- bsp_imp_tib %>%
  tidyr::separate(Aln_Key, into=c("Aln_Cgn", "Aln_TB", "Aln_CO"), sep='-') %>%
  tidyr::separate(Aln_CO, into=c("Aln_CO", "Aln_Des", "Aln_Din"), sep="_") %>%
  clean_tibble()

bsp_sep_des <- bsp_sep_tib %>% split(.$Ord_Des)

cur_fin_tib <- dplyr::inner_join(
  bsp_sep_des$U,bsp_sep_des$M, 
  by=c("Imp_Cgn",
       "Aln_Cgn",
       "Imp_Chr","Imp_Pos",
       "Ord_Din"),
  suffix=c("_U","_M")
)

mid_csv <- file.path(run$manDir, "mid-set.csv.gz")
readr::write_csv(cur_fin_tib, mid_csv)

# LEFT OFF HERE::

cur_fin_tib %>% dplyr::inner_join(new_imp_int_tib, by=c("Imp_Chr"="Chromosome", "Imp_Pos"="Coordinate")) %>% dplyr::distinct(Imp_Cgn,Imp_TB_U,Imp_CO_U)
cur_fin_tib %>% dplyr::inner_join(new_imp_int_tib, by=c("Imp_Chr"="Chromosome", "Imp_Pos"="Coordinate")) %>% dplyr::distinct(Aln_Cgn)

bsp_tmp_tib <- bsp_imp_tib %>% 
  dplyr::distinct(Aln_Prb, .keep_all = TRUE) %>% 
  dplyr::select(Imp_Cgn,Imp_TB,Imp_CO,Aln_Key,Bsp_Tag,Aln_Prb) %>% 
  tidyr::separate(Aln_Key, into=c("Aln_Cgn", "Aln_TB", "Aln_CO"), sep='-') %>%
  tidyr::separate(Aln_CO, into=c("Aln_CO", "Aln_Des", "Aln_Din"), sep="_") %>%
  clean_tibble()

bsp_tmp_tib %>% 
  dplyr::filter(Imp_Cgn==Aln_Cgn) %>% 
  dplyr::distinct(Imp_Cgn,Aln_Des,Aln_Din, .keep_all=TRUE) %>% 
  dplyr::group_by(Aln_Des,Aln_Din,Bsp_Tag) %>% 
  dplyr::summarise(Count=n(), .groups="drop")

#
# 1. Split
# 2. Join (cg, pos)
# 3. Unique (pos)
bsp_imp_tib









#
#
# LEFT OFF HERE: April 11 2021
#
#


if (!is.null(opt$ords)) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #               1.0 Data Collection:: Ord/Mat/AQP/PQC/Manifest
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  stamp_vec <- c(stamp_vec,run$add_pas_csv)
  if (opt$fresh || !valid_time_stamp(stamp_vec)) {
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     1.1.0 Data Collection:: Order Files
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # TBD:: Pass in all pre-defined cols (val, sel, key)
    #
    
    # This method could be replaced with 
    #  purrr::imap(vec, function(x,i) { load_man_file(x,i, ...) })
    #  - Current method is simple and works, but need to check for
    #    AQP/PQC null vecs up front...
    #
    ords_dat_tib <- ords_vec %>%
      lapply(load_man_file,
             verbose=opt$verbose,tt=pTracker) %>% 
      dplyr::bind_rows(.id="Dat_IDX") %>%
      utils::type.convert() %>%
      dplyr::mutate(across(where(is.factor),  as.character) ) %>%
      dplyr::left_join(man_info_tib, by="Dat_IDX") %>%
      dplyr::add_count(Ord_Prb, name="Ord_Prb_Rep") %>%
      dplyr::add_count(Ord_Prb,Ord_Par, name="Ord_Par_Rep")
    
    ords_sum_tib <- ords_dat_tib %>% 
      dplyr::group_by(Dat_IDX,Dat_BPN,Dat_AQP, 
                      Ord_Din,Ord_Des,Ord_Col,
                      Ord_Prb_Rep,Ord_Par_Rep) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    ords_sum_tib %>% print(n=base::nrow(ords_sum_tib))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    1.2.0 Data Collection:: Match Files
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    mats_dat_tib <- mats_vec %>%
      lapply(load_man_file,
             verbose=opt$verbose,tt=pTracker) %>% 
      dplyr::bind_rows(.id="Dat_IDX") %>%
      utils::type.convert() %>%
      dplyr::mutate(across(where(is.factor),  as.character) ) %>%
      dplyr::left_join(man_info_tib, by="Dat_IDX") %>% 
      dplyr::add_count(Address, name="Mat_Add_Rep")
    
    mats_sum_tib <- mats_dat_tib %>%
      dplyr::group_by(Dat_IDX,Dat_BPN,Dat_AQP,Mat_Add_Rep) %>%
      dplyr::summarise(Count=n(), .groups="drop")
    mats_sum_tib %>% print(n=base::nrow(mats_sum_tib))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     1.3.1 Data Collection:: AQP Files
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    aqps_dat_tib <- NULL
    if (!is.null(aqpn_vec) && length(aqps_vec)!=0) {
      
      aqps_dat_tib <- aqps_vec %>%
        lapply(load_man_file,
               verbose=opt$verbose,tt=pTracker) %>% 
        dplyr::bind_rows(.id="Dat_IDX") %>%
        utils::type.convert() %>%
        dplyr::mutate(across(where(is.factor),  as.character) ) %>%
        dplyr::left_join(man_info_tib, by="Dat_IDX")
      
      aqps_sum_tib <- aqps_dat_tib %>%
        dplyr::group_by(Dat_IDX,Dat_BPN,Dat_AQP,Decode_Status) %>%
        dplyr::summarise(Count=n(), .groups="drop")
      aqps_sum_tib %>% print(n=base::nrow(aqps_sum_tib))
      
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     1.3.2 Data Collection:: PQC Files
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # NOTE:: There's only one PQC, so no need to join with man_info_tib
    #
    pqcs_dat_tib <- NULL
    pqcs_dat_tib <- pqcs_vec %>%
      lapply(load_man_file,
             verbose=opt$verbose,tt=pTracker) %>% 
      dplyr::bind_rows(.id="Dat_PQC") %>%
      utils::type.convert() %>%
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    pqcs_sum_tib <- pqcs_dat_tib %>% 
      dplyr::group_by(Dat_PQC,Decode_Status) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    pqcs_sum_tib %>% print(n=base::nrow(pqcs_sum_tib))
    
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                    2.0 Data Merging:: Ord/Mat/AQP/PQC
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    add_pas_tib <- NULL
    if (!is.null(pqcs_dat_tib)) {
      
      add_pqc_dat_tib <- 
        ords_dat_tib %>% 
        dplyr::full_join(mats_dat_tib, by=c("Ord_Prb"="Mat_Prb","Dat_IDX","Dat_BPN","Dat_AQP")) %>%
        dplyr::left_join(pqcs_dat_tib, by="Address")
      
      add_pqc_sum_tib <- 
        add_pqc_dat_tib %>% 
        dplyr::group_by(Dat_IDX,Dat_BPN,Dat_AQP,Dat_PQC,
                        Ord_Des,Ord_Col,Ord_Din,Decode_Status) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      add_pqc_sum_tib %>% print(n=base::nrow(add_pqc_sum_tib))
      
      add_pas_tib <- add_pqc_dat_tib
      # add_pas_tib <- add_pqc_dat_tib %>% 
      #   dplyr::filter(Decode_Status==0) %>% 
      #   dplyr::arrange(-Dat_IDX,Ord_Key) %>%
      #   dplyr::distinct(Address, .keep_all=TRUE) %>%
      #   dplyr::select(Ord_Key,Ord_Des,Ord_Prb, dplyr::everything())
      
    } else if (!is.null(aqps_dat_tib)) {
      
      add_aqp_dat_tib <- 
        ords_dat_tib %>% 
        dplyr::full_join(mats_dat_tib, by=c("Ord_Prb"="Mat_Prb","Dat_IDX","Dat_BPN","Dat_AQP")) %>%
        dplyr::left_join(aqps_dat_tib, by=c("Address","Dat_IDX","Dat_BPN","Dat_AQP") )
      
      add_aqp_sum_tib <- 
        add_aqp_dat_tib %>% 
        dplyr::group_by(Dat_IDX,Dat_BPN,Dat_AQP,
                        Ord_Des,Ord_Col,Ord_Din,Decode_Status) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      add_aqp_sum_tib %>% print(n=base::nrow(add_aqp_sum_tib))
      
      add_pas_tib <- add_aqp_dat_tib
      # add_pas_tib <- add_aqp_dat_tib %>% 
      #   dplyr::filter(Decode_Status==0) %>% 
      #   dplyr::arrange(-Dat_IDX,Ord_Key) %>% 
      #   dplyr::distinct(Address, .keep_all=TRUE) %>%
      #   dplyr::select(Ord_Key,Ord_Des,Ord_Prb, dplyr::everything())
      
    } else {
      stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Both AQP and PQC are length zero!!!{RET}{RET}"))
    }
    
    add_pas_tib <- add_pas_tib %>% 
      dplyr::filter(Decode_Status==0) %>% 
      dplyr::arrange(-Dat_IDX,Ord_Key) %>%
      dplyr::distinct(Address, .keep_all=TRUE) %>%
      dplyr::select(Ord_Key,Ord_Des,Ord_Prb, dplyr::everything())
    
    add_pas_sum_tib <- 
      add_pas_tib %>% 
      dplyr::group_by(Dat_IDX,Dat_BPN,Dat_AQP,Dat_PQC,
                      Ord_Des,Ord_Col,Ord_Din,Decode_Status) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    add_pas_sum_tib %>% print(n=base::nrow(add_pas_sum_tib))
    
    safe_write(add_pas_tib,"csv",run$add_pas_csv, funcTag=par$prgmTag, 
               verbose=opt$verbose)
    
  } else {
    
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Loading add_pas_csv={run$add_pas_csv}...{RET}"))
    add_pas_tib <- suppressMessages(suppressWarnings( 
      readr::read_csv(run$add_pas_csv, guess_max=100000) )) %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
  }
}

if (!is.null(opt$mans)) {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     1.4.0 Data Collection:: Manifest Files
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  mans_dat_tib <- NULL
  if (!is.null(mans_vec) && length(mans_vec)!=0) {

    mans_dat_tib <- mans_vec %>%
      lapply(loadManifestGenomeStudio,
             addSource=TRUE, normalize=TRUE, retType="man", 
             verbose=opt$verbose+10,tt=pTracker)
    
    # Left off here fill in option; Add Souce...
    
    # TBD:: Add summary to manifest loading...
    # mans_sum_tib <- mans_dat_tib %>%
    #   dplyr::group_by(Dat_IDX,Dat_BPN,Dat_AQP,Mat_Add_Rep) %>%
    #   dplyr::summarise(Count=n(), .groups="drop")
    # mans_sum_tib %>% print(n=base::nrow(mans_sum_tib))
    
  }
  
}






# AQP/Manifest -> Address -> manifest
# Coordinates  -> Designs -> manifest/order

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  2.1 Functional Manifest Generation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# add_comb(add_pas_list[['U']], add_pas_list[['M']],
#          field="Address",
#          join=c("Ord_Key","Ord_Din","Ord_Col","Decode_Status"),
#          verbose=opt$verbose, tt=pTracker)

man_fun_vec <- c("Ord_Key","Ord_Din","Ord_Col",
                 dplyr::all_of(base::names(man_info_tib) ) )
stamp_vec <- c(stamp_vec,run$man_fun_csv)

if (opt$fresh || !valid_time_stamp(stamp_vec)) {
  man_fun_tib <- add_pas_tib %>%
    add_to_man(join=man_fun_vec,
               runName=opt$runName,
               des_key="Ord_Des", pid="Ord_Key",
               col_key="Ord_Col",
               csv=run$man_fun_csv, 
               validate=TRUE, 
               verbose=opt$verbose, tt=pTracker)
} else {
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading man_fun_csv={opt$man_fun_csv}...{RET}"))
  man_fun_tib <- suppressMessages(suppressWarnings( 
    readr::read_csv(run$man_fun_csv, guess_max=100000) )) %>%
    utils::type.convert() %>% 
    dplyr::mutate(across(where(is.factor),  as.character) )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    3.0 Alignment/Sequence Based CGN Mapping::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  3.1 Format Fasta AND U49/U50 Probe Seqs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

stamp_vec <- c(stamp_vec,
               run$add_dat_csv,run$add_prb_fas,
               run$add_u49_tsv,run$add_m49_tsv)

add_pas_fas_tib <- NULL
if (opt$fresh || !valid_time_stamp(stamp_vec)) {
  
  add_pas_fas_tib <-
    add_to_fas(
      tib=add_pas_tib, 
      prb_key="Ord_Prb",
      add_key="Address", des_key="Ord_Des", type_key="Ord_Din",
      prb_fas=run$add_prb_fas, dat_csv=run$add_dat_csv,
      u49_tsv=run$add_u49_tsv, m49_tsv=run$add_m49_tsv,
      verbose=opt$verbose, tt=pTracker)
  
} else {
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading add_dat_csv={run$add_dat_csv}...{RET}"))
  
  add_pas_fas_tib <- suppressMessages(suppressWarnings(
    readr::read_csv(run$add_dat_csv, guess_max=100000) )) %>%
    utils::type.convert() %>% 
    dplyr::mutate(across(where(is.factor),  as.character) )
  
}

add_pas_fas_sum <- add_pas_fas_tib %>% 
  dplyr::group_by(Dat_IDX, Ord_Din, Ord_Des) %>% 
  dplyr::summarise(Count=n(), .groups="drop")
add_pas_fas_sum %>% print(n=base::nrow(add_pas_fas_sum))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    3.2 Align All Probe Sequence:: BSMAP
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

stamp_vec <- c(stamp_vec, run$add_prb_bspz)
if (opt$fresh || !valid_time_stamp(stamp_vec)) {
  
  run$add_prb_bspz <- 
    run_bsmap(
      exe=opt$bsmap_exe, 
      fas=run$add_prb_fas, 
      gen=run$gen_ref_fas,
      bsp=run$add_prb_bsp,
      opt=NULL, lan=NULL, run=TRUE,
      verbose=opt$verbose,tt=pTracker)
  
} else {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Using existing file add_prb_bspz={run$add_prb_bspz}...{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  3.3 Join Address and Alignment Data:: BSMAP
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsp_imp_tib <- NULL
stamp_vec <- c(stamp_vec, run$add_pas_bsp_csv)
if (opt$fresh || !valid_time_stamp(stamp_vec)) {
  
  bsp_imp_tib <-
    join_bsmap(
      add=add_pas_fas_tib,
      bed=cgn_bed_tib, org=add_org_tib,
      file=run$add_prb_bspz,
      join_key="Aln_Key", 
      prb_des_key="Ord_Des", prb_din_key="Ord_Din",
      join_type="inner",
      sort=TRUE,
      verbose=opt$verbose,tt=pTracker)
  
  safe_write(bsp_imp_tib,"csv",run$add_pas_bsp_csv, funcTag=par$prgmTag, 
             verbose=opt$verbose)
  
} else {
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading add_pas_bsp_csv={run$add_pas_bsp_csv}...{RET}"))
  
  bsp_imp_tib <- suppressMessages(suppressWarnings(
    readr::read_csv(run$add_pas_bsp_csv, guess_max=100000) )) %>%
    utils::type.convert() %>% 
    dplyr::mutate(across(where(is.factor),  as.character) )
  
}

if (FALSE) {
  # Create basic Genomic Regions RDS as well::
  #  NOTE:: Currently not used, so skipping for now...
  #
  add_pas_bsp_grs <- NULL
  stamp_vec <- c(stamp_vec, run$add_pas_grs_rds)
  if (opt$fresh || !valid_time_stamp(stamp_vec)) {
    add_pas_bsp_grs <-
      GRanges(
        seqnames = Rle(bsp_imp_tib$Bsp_Chr),
        strand   = Rle(stringr::str_sub(bsp_imp_tib$Bsp_Srd, 1,1)),
        
        IRanges(start = bsp_imp_tib$Bsp_Pos,
                width = 1,
                names=bsp_imp_tib$Aln_Key_Unq)
      )
    safe_write(add_pas_bsp_grs,"rds",run$add_pas_grs_rds, funcTag=par$prgmTag,
               verbose=opt$verbose)
    
  } else {
    # NOTE:: Currently not used, so skipping for now...
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Loading add_pas_grs_rds={run$add_pas_grs_rds}...{RET}"))
    add_pas_bsp_grs <- readr::read_rds(run$add_pas_grs_rds)
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#         3.4 Intersect Sequences Address and improbe:: U49/M49
#
#                               NEW VERSION::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imp_prb_dir <- file.path(opt$impDir, "scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49")
imp_u49_tsv <- file.path(imp_prb_dir, paste("probe_U49_cgn-table.csv.gz", sep="-") )
imp_m49_tsv <- file.path(imp_prb_dir, paste("probe_M49_cgn-table.csv.gz", sep="-") )

seq_imp_tib <- NULL
stamp_vec <- c(stamp_vec, 
               run$int_u49_tsv,
               run$int_m49_tsv,
               run$int_seq_tsv)
if (opt$fresh || !valid_time_stamp(stamp_vec)) {
  
  int_u49_tib2 <- 
    intersect_seq(ref=imp_u49_tsv,can=run$add_u49_tsv, out=run$int_u49_tsv,
                  idxA=1, idxB=1, verbose=opt$verbose,tt=pTracker)
  
  int_m49_tib2 <- 
    intersect_seq(ref=imp_m49_tsv,can=run$add_m49_tsv, out=run$int_m49_tsv,
                  idxA=1, idxB=1, verbose=opt$verbose,tt=pTracker)
  
  seq_imp_tib <- 
    join_seq_intersect(u49=int_u49_tib2, m49=int_m49_tib2, 
                       bed=cgn_bed_tib, org=add_org_tib,
                       verbose=opt$verbose, tt=pTracker)
  
  safe_write(seq_imp_tib,"tsv",run$int_seq_tsv, funcTag=par$prgmTag,
             verbose=opt$verbose)
  
} else {
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading imp_seq_tsv={run$int_seq_tsv}...{RET}"))
  seq_imp_tib <- suppressMessages(suppressWarnings( readr::read_tsv(run$int_seq_tsv) )) %>%
    utils::type.convert() %>% 
    dplyr::mutate(across(where(is.factor),  as.character) )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#              3.5 Join improbe Seq and bsmap Pos Data:: U49/M49/BSP
#
#   CGN Identification::
#    - seq_imp_tib = improbe CGN's matched by sequence
#    - Bsp_Pos_tib = genomic alignment coordinates (BSMAP)
#   Further Processing::
#    - bsp_imp_tib = annotated Bsp_Pos_tib
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Need to normalize matching method fields!!!
#
#   - BSP:: Add,Des,Din, [Bsp_Chr,Bsp_Pos,Bsp_Srd]  [aln_seq(50)]
#   - SEQ:: Add,Des,Din, [Imp_Cgn,Imp_Srd3],        [Imp_Seq=aln_seq49,Imp_Nuc]

# Remove::
#   -Mat_Src_Bsp, -Mat_Src_Seq, -Bsp_Srd

# What we need is a match score:: 
#   - A Null value is the highest score
#   - Perfect match is zero
#  Srd_Mis_Scr = count(Imp_FR/Imp_TB/Imp_CO/Imp_Top_Srd_Bsp)
#  Cgn_Mis_Scr = Imp_Cgn
#  Nxb_Mis_Scr = Imp_Nxb
#  Can_Mis_Scr = Can_Mis_Scr ???

#  Bsp_Mis_Scr = is.na(Imp_Cgn_Seq)
#  Seq_Mis_Scr = is.na(Imp_Cgn_Bsp)


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
  dplyr::filter(Group_Count>1)

all_imp_tab %>% 
  dplyr::add_count(Address,Ord_Des,Ord_Din,Imp_Chr,Imp_Pos,Imp_FR,Imp_TB,Imp_CO, name="Group_Count") %>% 
  dplyr::group_by(Group_Count,Mat_Src_Key,Bsp_Din_Scr,Ord_Des) %>% 
  dplyr::summarise(Count=n(), .groups="drop")

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

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Row Binding all_imp_tab=(seq_imp_tib/bsp_imp_tib).{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 3.4.1 Join BSMAP & Seq-Match into Tibble::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Row Joining all_imp_tib=(seq_imp_tib/bsp_imp_tib)...{RET}"))

all_imp_tib <- NULL
all_imp_tib <- all_imp_tab %>%
  tidyr::pivot_wider(id_cols = dplyr::all_of(c(com_des_cols,com_gen_cols)), 
                     names_from = "Mat_Src_Key", 
                     values_from = dplyr::all_of(c(com_scr_cols,com_cgn_cols)) ) %>%
  # Now Add Scores for pair matching::
  dplyr::mutate(
    Cgn_Mis_Scr=dplyr::case_when(
      is.na(Imp_Cgn_Seq) ~ 3,
      is.na(Imp_Cgn_Bsp) ~ 2,
      Imp_Cgn_Bsp != Imp_Cgn_Seq ~ 1,
      Imp_Cgn_Bsp == Imp_Cgn_Seq ~ 0,
      TRUE ~ 9
    ) %>% as.integer(),
    
    Nxb_Mis_Scr=dplyr::case_when(
      Ord_Des=="2" ~ 0,
      is.na(Imp_Nxb_Seq) ~ 3,
      is.na(Imp_Nxb_Bsp) ~ 2,
      Imp_Nxb_Bsp != Imp_Nxb_Seq ~ 1,
      Imp_Nxb_Bsp == Imp_Nxb_Seq ~ 0,
      TRUE ~ 9
    ) %>% as.integer(),
    
    Org_Mis_Scr=dplyr::case_when(
      is.na(Org_Mis_Scr_Seq) ~ 3,
      is.na(Org_Mis_Scr_Bsp) ~ 2,
      !is.na(Org_Mis_Scr_Seq) & !is.na(Org_Mis_Scr_Bsp) ~ 0,
      TRUE ~ 9) %>% as.integer(),
    
    Imp_Cgn_Max=dplyr::case_when(
      !is.na(Imp_Cgn_Seq) & Can_Mis_Scr_Seq==0 ~ Imp_Cgn_Seq,
      !is.na(Imp_Cgn_Bsp) & Can_Mis_Scr_Bsp==0 ~ Imp_Cgn_Bsp,
      !is.na(Imp_Cgn_Seq) ~ Imp_Cgn_Seq,
      !is.na(Imp_Cgn_Bsp) ~ Imp_Cgn_Bsp,
      TRUE ~ NA_integer_
    ) %>% as.integer(),
    
    Infinium_Design=dplyr::case_when(
      Ord_Des=="U" ~ 1,
      Ord_Des=="M" ~ 1,
      Ord_Des=="2" ~ 2,
      TRUE ~ NA_real_
    ) %>% as.integer()
    
  ) %>%
  mutate_probe_id(
    pid="Reference_Cgn_ID",cgn="Imp_Cgn_Max",
    des="Ord_Des",din="Ord_Din",
    tb="Imp_TB",co="Imp_CO",
    inf="Infinium_Design",
    verbose=opt$verbose, tt=pTracker) %>%
  dplyr::select(Address,Ord_Des,Ord_Din,
                Reference_Cgn_ID,dplyr::everything())

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Row Joining all_imp_tib=(seq_imp_tib/bsp_imp_tib).{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 3.4.1 Build Summary Statistic Tables::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Building Stats Tables...{RET}"))

# Genomic Hits Table::
bsp_aln_cnt_tib <- NULL
bsp_aln_cnt_tib <- bsp_imp_tib %>% 
  dplyr::distinct(Address,Ord_Des,Ord_Din, 
                  Imp_Chr,Imp_Pos, Imp_FR,Imp_TB,Imp_CO) %>%
  dplyr::select(Address,Ord_Des,Ord_Din) %>%
  dplyr::group_by_all() %>% 
  dplyr::summarise(Aln_Hits=n(), .groups="drop")

# CGN Hits Table::
bsp_cgn_cnt_tib <- NULL
bsp_cgn_cnt_tib <- all_imp_tab %>% 
  dplyr::distinct(Address,Ord_Des,Ord_Din,Imp_Cgn) %>%
  dplyr::select(Address,Ord_Des,Ord_Din) %>%
  dplyr::group_by_all() %>% 
  dplyr::summarise(Cgn_Hits=n(), .groups="drop")

# Add all Din Counts::
bsp_din_cnt_tib <- NULL
bsp_din_cnt_tib <- bsp_imp_tib %>% 
  dplyr::distinct(Address,Ord_Des,Ord_Din, 
                  Imp_Chr,Imp_Pos, Imp_FR,Imp_TB,Imp_CO, 
                  Bsp_Din_Ref) %>%
  dplyr::select(Address,Ord_Des,Ord_Din,Bsp_Din_Ref) %>%
  dplyr::group_by_all() %>% 
  dplyr::summarise(Count=n(), .groups="drop") %>% 
  tidyr::pivot_wider(
    id_cols = c("Address","Ord_Des","Ord_Din"), 
    names_from="Bsp_Din_Ref", 
    names_prefix="Din_Cnt_",
    values_from="Count", 
    values_fill=0) # %>% dplyr::arrange(-TG)

# Add all Din Counts::
bsp_nxb_cnt_tib <- NULL
bsp_nxb_cnt_tib <- bsp_imp_tib %>% 
  dplyr::distinct(Address,Ord_Des,Ord_Din, 
                  Imp_Chr,Imp_Pos, Imp_FR,Imp_TB,Imp_CO, 
                  Imp_Nxb) %>%
  dplyr::select(Address,Ord_Des,Ord_Din,Imp_Nxb) %>%
  dplyr::group_by_all() %>% 
  dplyr::summarise(Count=n(), .groups="drop") %>% 
  tidyr::pivot_wider(
    id_cols = c("Address","Ord_Des","Ord_Din"), 
    names_from="Imp_Nxb", 
    names_prefix="Nxb_Cnt_",
    values_from="Count", 
    values_fill=0)

bsp_sum_tab <- NULL
bsp_sum_tab <- bsp_cgn_cnt_tib %>% 
  dplyr::inner_join(bsp_aln_cnt_tib, by=c("Address","Ord_Des","Ord_Din")) %>% 
  dplyr::inner_join(bsp_din_cnt_tib, by=c("Address","Ord_Des","Ord_Din")) %>% 
  dplyr::inner_join(bsp_nxb_cnt_tib, by=c("Address","Ord_Des","Ord_Din"))

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Building Stats Tables.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              3.4.1 Select "Top" Canonical (Unique) Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Selecting Top Canonical Hits...{RET}"))

#
# TBD:: Pick best cgn, tally stats...
#   - move top after pivot???
#   - add matching type to top-score???
#

# Top Hits for Manifest Skeleton::
top_scr_col <- names(all_imp_tib %>% dplyr::select(Cgn_Mis_Scr, dplyr::contains("_Scr")))
top_imp_tib <- NULL
top_imp_tib <- all_imp_tib %>% 
  dplyr::arrange( dplyr::across( dplyr::all_of(top_scr_col) )) %>%
  dplyr::distinct( 
    dplyr::across( dplyr::all_of(com_des_cols) ), .keep_all=TRUE ) %>%
  mutate_probe_id(
    pid="Canonical_Cgn_ID",cgn="Imp_Cgn_Max",
    des="Ord_Des",din="Ord_Din",
    tb="Imp_TB",co="Imp_CO",
    inf="Infinium_Design",
    verbose=opt$verbose, tt=pTracker) %>%
  dplyr::select(Address,Ord_Des,Ord_Din,
                Canonical_Cgn_ID,Reference_Cgn_ID,
                dplyr::everything()) %>%
  dplyr::right_join(add_pas_fas_tib, 
                    by=com_des_cols) %>%
  dplyr::mutate(Species=opt$Species,
                Genome_Build=opt$genBuild)

#
# TBD:: Validate that all original passing probes were matched...
#

# Top Hits Summary::
top_imp_sum <- NULL
top_imp_sum <- top_imp_tib %>% 
  dplyr::group_by( dplyr::across( dplyr::all_of( c("Ord_Des","Ord_Din", top_scr_col) )) ) %>%
  dplyr::summarise(Count=n(), .groups="drop") %>% 
  dplyr::arrange(-Count)
top_imp_sum %>% print(n=base::nrow(top_imp_sum))

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Selecting Top Canonical Hits.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                3.4.1 Merge Back Canonical (Unique) Cgn ID::
#                         Also merge summary stats
#                           Write full alignments
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fin_imp_tib <- NULL
fin_imp_tib <- top_imp_tib %>% 
  dplyr::select(Address,Ord_Des,Ord_Din,Canonical_Cgn_ID) %>% 
  dplyr::right_join(all_imp_tib, by=c("Address","Ord_Des","Ord_Din")) %>%
  dplyr::left_join(bsp_sum_tab, by=c("Address","Ord_Des","Ord_Din") ) %>%
  dplyr::mutate(Seq_ID=paste(Address,Ord_Des,Ord_Din, sep="_"))

fin_imp_cnt <- safe_write(fin_imp_tib,"tsv",run$add_fin_csv, funcTag=par$prgmTag,
                          verbose=opt$verbose)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#           3.4.2 Extract improbe template sequence from Genome::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run$use_top <- TRUE

if (run$use_top) {
  imp_pos_tib <- top_imp_tib %>%
    dplyr::mutate(Seq_ID=paste(Address,Ord_Des,Ord_Din, sep="_")) %>%
    dplyr::select(Seq_ID, Imp_Chr,Imp_Pos) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Imp_Chr,Imp_Pos)
  
} else {
  imp_pos_tib <- fin_imp_tib %>% 
    dplyr::select(Seq_ID, Imp_Chr,Imp_Pos) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Imp_Chr,Imp_Pos)
}

imp_dna_tib <- 
  dna_to_template(tib=imp_pos_tib, dna=ref_fwd_seq, map=ref_fwd_tab, 
                  file=run$imp_inp_tsv,
                  name="Seq_ID",gen=opt$genBuild,
                  chr1="Imp_Chr",pos="Imp_Pos",chr2="Chrom_Char",
                  ext_seq="Fwd_Seq",des_seq="Des_Seq",imp_seq="Sequence",
                  # iupac = "QI_T",
                  add_flank=FALSE,
                  verbose = opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                3.4.3 Run improbe on genome template 122mers::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: Return Strand_[FR/CO/TB] in improbe_docker...
#
imp_des_tib <- 
  improbe_docker(dir = run$desDir, file = run$imp_inp_tsv, 
                 name = opt$runName, 
                 image = image_str, shell = image_ssh,
                 verbose = opt$verbose, tt=pTracker) %>%
  dplyr::mutate(Methyl_Allele_TB_Strand=stringr::str_sub(Methyl_Allele_TB_Strand,1,1))

# TBD:: Need to design all possible designs, but with iupac R code so we can match
#  probes via Infinium I/II, etc...

iup_des_tib <- imp_des_tib %>% 
  dplyr::select(Seq_ID,Forward_Sequence) %>% 
  dplyr::mutate(Probe_Type="cg") %>%
  desSeq_to_prbs(ids_key = "Seq_ID", seq_key = "Forward_Sequence", prb_key = "Probe_Type",
                 strsSR = "FR",strsCO = "CO", 
                 addMatSeq = TRUE, parallel = TRUE, # max = 10,
                 verbose = opt$verbose, tt = pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#        3.4.4 Join original designs with template 122mers designs::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD: Join top designs with iup_des_tib and then with matching imp_des_tib
#

imp_des_tib %>% dplyr::inner_join(iup_des_tib, by=c("Methyl_Probe_Sequence"="PRB1_M_MAT"))

tmp_des_tib <- dplyr::inner_join(
  imp_des_tib, iup_des_tib, 
  by=c("Seq_ID",
       "UnMethyl_Probe_Sequence"="PRB1_U_MAT",
       "Methyl_Probe_Sequence"="PRB1_M_MAT"
  )
)


# Left off here...
#
#
# Now join with top hits...
#

tmp_des_tib %>% 
  dplyr::inner_join(dplyr::mutate(top_imp_tib, (Seq_ID=paste(Address,Ord_Des,Ord_Din, sep="_")), by="Seq_ID"))




# NOTE: (Seperate Idea...) Instead of defininng regions by general genomic (GC% for islands) defined them by sample epigenics

top_imp_tib %>% 
  dplyr::mutate(Seq_ID=paste(Address,Ord_Des,Ord_Din, sep="_")) %>%
  dplyr::left_join(imp_des_tib, by="Seq_ID")

# Try joining with TB/CO as well...

top_imp_tib %>% 
  dplyr::mutate(Seq_ID=paste(Address,Ord_Des,Ord_Din, sep="_")) %>%
  dplyr::inner_join(imp_des_tib, by="Seq_ID") %>%
  dplyr::filter(Imp_TB==Methyl_Allele_TB_Strand &
                  Imp_CO==Methyl_Allele_CO_Strand) %>%
  dplyr::filter(Ord_Prb==UnMethyl_Probe_Sequence)
head() %>% as.data.frame()


#  dplyr::select(Seq_ID,)





























#
# - top_imp_
#
# - [done]: Merge Canonical IDs and Summary Tables::
#
# - Address/CGN Manifest Pairing
#   - Check for any missing probes
#   - Add Infinium Type Integer, Proper Probe ID, 
#     Original Function Sesame Probe Name (cg########),
#     Original Design ID (Horvath, etc.)
#
# - [done] Add Source Design Coordinate Checks
# - Add Multiple Genome Support
# - Add Controls
# - Add Manifest Input
# - Add Core/HMM Annotations
# - Copy all new data source files to central location
#   Transfer to cluster
#
# - Performance Summary::
#    [done] Din Table
#    - Nxb.I, Cpg.I/II Pivot Table
#
# - Alignment Summary::
#    [done] Hit Table
#    [done] Cgn Table
#
# - SNP Summary::
#   - Degenerate Designs
#   - Next Base Ratio
#
# - Extract SNP/Ref 122mer/Probe
# - Redesign R from 122mer
#
# - Write Master Manifest
# - Write Sesame Manifest
# - Write Genome Studio Manifest
#
# - Add Comments
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#       3.6.0 Rebuild Manifest From Address Alignments:: Top/Master
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_join_pid_vec <- c("Canonical_Cgn_ID","Ord_Din")

#
# Need to ensure that there is a mate for Infinum I pairs::
#

man_top_master_tib <- NULL
man_top_master_tib <-
  add_to_man(top_imp_tib, 
             join=man_join_pid_vec,
             runName=opt$runName,
             des_key="Ord_Des", 
             pid_key="Probe_ID",
             rep_key="Canonical_Cgn_ID",
             rep_val="Replicate_Count",
             nxb_key="Imp_Nxb_Bsp_U",
             csv=run$man_mas_csv,
             validate=TRUE,
             verbose=opt$verbose, tt=pTracker)

man_top_master_tib %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#       3.6.1 Rebuild Manifest From Address Alignments:: Top/Sesame
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_top_sesame_tib <- man_top_master_tib %>% 
  dplyr::rename(
    # Add Name...
    U=Address_U,
    M=Address_M,
    AlleleA_ProbeSeq=Ord_Prb_U,
    AlleleB_ProbeSeq=Ord_Prb_M,
    Probe_Type=Ord_Din,
    Strand_FR=Imp_FR_U,
    Strand_TB=Imp_TB_U,
    Strand_CO=Imp_CO_U,
    CHR=Imp_Chr_U,
    MAPINFO=Imp_Pos_U,
    Species=Species_U,
    Genome_Build=Genome_Build_U,
    Rep_Num=Replicate_Count
  ) %>%
  dplyr::arrange(Probe_ID) %>%
  dplyr::select(
    Probe_ID,
    U,AlleleA_ProbeSeq,
    M,AlleleB_ProbeSeq,
    Next_Base,
    Color_Channel,
    col,
    Probe_Type,
    Strand_FR,Strand_TB,Strand_CO,
    Infinium_Design,
    Infinium_Design_Type,
    CHR,MAPINFO,Species,Genome_Build,
    
    Rep_Num
  )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    3.7 Extract 122mer & SNP Probes::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  # Define Run Time:: Ref Alignment Genome
  ref_fwd_seq <- NULL
  run$gen_ref_fas <- 
    file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
              paste0(opt$genBuild,".genome.fa.gz"))
  ref_fwd_seq <- 
    Biostrings::readDNAStringSet(filepath = run$gen_ref_fas, format = "fasta") # , nrec = 2)
  
  
  # Define Run Time:: SNP IUPAC Genome
  snp_fwd_seq <- NULL
  run$gen_snp_fas <- 
    file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
              paste0(opt$genBuild,".dbSNP151-genome.fa.gz"))
  if (!is.null(run$gen_snp_fas) && file.exists(run$gen_snp_fas)) {
    snp_fwd_seq <- 
      Biostrings::readDNAStringSet(filepath = run$gen_snp_fas, format = "fasta") # , nrec = 2)
  }
  
  ref_fwd_tab <- ref_fwd_seq %>% names() %>% 
    stringr::str_remove(" .*$") %>% 
    stringr::str_remove("^chr") %>%
    tibble::tibble() %>% 
    purrr::set_names("Chrom_Char") %>% 
    dplyr::mutate(Idx=dplyr::row_number(),
                  Chrom_Char=paste0("chr",Chrom_Char) )
  print(ref_fwd_tab)
  
  snp_fwd_tab <- snp_fwd_seq %>% names() %>% 
    stringr::str_remove(" .*$") %>% 
    stringr::str_remove("^chr") %>%
    tibble::tibble() %>% 
    purrr::set_names("Chrom_Char") %>% 
    dplyr::mutate(Idx=dplyr::row_number(),
                  Chrom_Char=paste0("chr",Chrom_Char) )
  print(snp_fwd_tab)
  
  #
  # Should Probably Join Everything First...
  #
  # tar_add_list <- all_imp_tib %>% 
  #   dplyr::filter(!is.na(Imp_Chr) & !is.na(Imp_Pos)) %>%
  #   split(.$Imp_Chr)
  
  # TBD:: Full List for later...
  #
  # tar_add_list <- 
  #   all_imp_tab %>% 
  #   dplyr::distinct(Address,Ord_Des,Ord_Din, Imp_Chr,Imp_Pos) %>% 
  #   dplyr::arrange(Imp_Chr,Imp_Pos) %>%
  #   split(.$Imp_Chr)
  
  tar_add_list <- 
    man_top_sesame_tib %>% 
    dplyr::distinct(Probe_ID, U, M, 
                    Probe_Type,Infinium_Design,
                    AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                    Strand_FR,Strand_TB,Strand_CO,
                    CHR,MAPINFO,Genome_Build) %>% 
    dplyr::rename(Imp_Chr=CHR,Imp_Pos=MAPINFO) %>%
    dplyr::arrange(Imp_Chr,Imp_Pos) %>%
    split(.$Imp_Chr)
  
  fwd_des_tib <- NULL
  fet_chroms <- names(tar_add_list)
  for (chr_str in fet_chroms) {
    cat(glue::glue("[{par$prgmTag}]: chr={chr_str}.{RET}"))
    
    if (is.null(tar_add_list[[chr_str]])) next
    
    chr_idx <- snp_fwd_tab %>% 
      dplyr::filter(Chrom_Char==chr_str) %>% 
      head(n=1) %>% pull(Idx) %>% as.integer()
    print(chr_idx)
    
    # bsp_begs <- head(tar_add_list[[chr_str]]$imp_pos) - 60
    bsp_begs <- tar_add_list[[chr_str]]$Imp_Pos - 60
    bsp_ends <- bsp_begs + 122 - 1
    # print(bsp_begs)
    
    ref_seqs <- stringr::str_sub( as.character(ref_fwd_seq[[chr_idx]]), bsp_begs, bsp_ends) %>% addBrac()
    snp_seqs <- stringr::str_sub( as.character(snp_fwd_seq[[chr_idx]]), bsp_begs, bsp_ends) %>% addBrac()
    # print(snp_seqs)
    
    cur_des_tib <- 
      tibble::tibble(
        Des_Key=tar_add_list[[chr_str]]$Probe_ID,
        Des_AddU=tar_add_list[[chr_str]]$U,
        Des_AddM=tar_add_list[[chr_str]]$M,
        # Des_Key=paste(tar_add_list[[chr_str]]$Address,
        #               tar_add_list[[chr_str]]$Ord_Des,
        #               tar_add_list[[chr_str]]$Ord_Din, sep="_"),
        # Des_Bld=opt$genBuild,
        Des_Bld=tar_add_list[[chr_str]]$Genome_Build,
        Des_Chr=chr_str,
        Des_Pos=tar_add_list[[chr_str]]$Imp_Pos,
        # Des_Din=tar_add_list[[chr_str]]$Ord_Din,
        Des_Inf=tar_add_list[[chr_str]]$Infinium_Design,
        Des_Din=tar_add_list[[chr_str]]$Probe_Type,
        Des_FR=tar_add_list[[chr_str]]$Strand_FR,
        Des_TB=tar_add_list[[chr_str]]$Strand_TB,
        Des_CO=tar_add_list[[chr_str]]$Strand_CO,
        
        Des_PrbA=tar_add_list[[chr_str]]$AlleleA_ProbeSeq,
        Des_PrbB=tar_add_list[[chr_str]]$AlleleB_ProbeSeq,
        
        Des_Ref_Seq=ref_seqs,
        Des_Snp_Seq=snp_seqs
      ) %>% 
      dplyr::mutate(
        Des_Chr=as.character(Des_Chr),
        Des_Pos=as.integer(Des_Pos)
      )
    
    # TBD:: Should Split based on Des_Din (cg, ch, rs, etc...)
    #
    des_ref_tib <- 
      desSeq_to_prbs(cur_des_tib, 
                     ids_key = "Des_Key", seq_key = "Des_Ref_Seq", prb_key = "Des_Din",
                     strsSR = "FR",strsCO = "CO", 
                     addMatSeq = TRUE, parallel = TRUE, # max = 10,
                     verbose = opt$verbose, tt = pTracker)
    
    # des_snp_tib <- 
    #   desSeq_to_prbs(cur_des_tib, 
    #                  idsKey = "Des_Key", seqKey = "Des_Snp_Seq", prbKey = "Des_Din", # prbKey = "cg",
    #                  strsSR = "FR",strsCO = "CO", 
    #                  addMatSeq = TRUE, parallel = TRUE, max = 10,
    #                  verbose = opt$verbose, tt = pTracker)
    
    cur_join_tib <- cur_des_tib %>% dplyr::inner_join(des_ref_tib, by=c("Des_Key") ) %>% 
      dplyr::mutate(Match=dplyr::case_when(
        Des_Inf==1 & Des_PrbA==PRB1_U_MAT & Des_PrbB==PRB1_M_MAT ~ 0,
        Des_Inf==2 & Des_PrbA==PRB2_D_MAT ~ 0,
        TRUE ~ 1)
      ) %>% dplyr::filter(Match==0) # %>% as.data.frame()
    
    if (FALSE) {
      #
      # Run improbe full design for thermodynamic scores::
      #
      opt$desDir <- file.path(opt$outDir, "des")
      if (!dir.exists(opt$desDir)) dir.create(opt$desDir, recursive = TRUE)
      ref_des_tsv <- file.path(opt$desDir, paste(opt$genBuild,chr_str,"improbe-ref-design-input.tsv", sep='-'))
      
      ref_des_tib <- tibble::tibble(
        Seq_ID=cur_des_tib$Des_Key,
        Sequence=cur_des_tib$Des_Ref_Seq,
        Genome_Build=cur_des_tib$Des_Bld,
        Chromosome=cur_des_tib$Des_Chr,
        Coordinate=cur_des_tib$Des_Pos,
        CpG_Island="FALSE"
      )
      readr::write_tsv(ref_des_tib,ref_des_tsv)
      
      # Not sure why this isn't running...
      ref_imp_des_tib <- 
        improbe_docker(dir = opt$desDir, file = ref_des_tsv, 
                       name = paste(opt$genBuild,chr_str, sep="-"), 
                       image = image_str, shell = image_ssh,
                       verbose = opt$verbose, tt=pTracker)
    }
    
    
    fwd_des_tib <- bind_rows(fwd_des_tib, cur_join_tib)
    # print(fwd_des_tib)
    
    cat(glue::glue("[{par$prgmTag}]: Done.{RET}{RET}"))
    
    # break
  }
  # fwd_des_tib %>% dplyr::filter(Match==0)
  
  #
  # Rebuild Top/Bot Calls on Design Sequence::
  #
  fwd_des_tb_tib <- fwd_des_tib %>%
    setTopBot_tib(seqKey="Des_Ref_Seq.x", srdKey = "Strand_TB", topKey = "Top_Sequence",
                  verbose=opt$verbose, tt=pTracker)
  
  fwd_des_csv <- file.path(run$alnDir, paste(opt$runName,"fwd-seqs.csv.gz", sep='-'))
  readr::write_csv(fwd_des_tb_tib, fwd_des_csv)
  
  
  fwd_sel_tib <- fwd_des_tb_tib %>% 
    dplyr::rename(Forward_Sequence=Des_Ref_Seq.x) %>%
    dplyr::mutate(Source_Seq=dplyr::case_when(
      Des_Inf==1 ~ PRB1_N,
      Des_Inf==2 ~ PRB2_N,
      TRUE ~ NA_character_) 
    ) %>%
    dplyr::select(Des_Key,Des_AddU,Des_AddM, 
                  Forward_Sequence,Top_Sequence,Source_Seq, 
                  Strand_SR,Strand_CO,Strand_TB)
  
  man_fwd_ses_col <- c("Probe_ID","Name","U","AlleleA_ProbeSeq","M","AlleleB_ProbeSeq",
                       "Next_Base","Color_Channel","col","Probe_Type",
                       "Strand_FR","Strand_TB","Strand_CO","Infinium_Design","Infinium_Design_Type",
                       "CHR","MAPINFO","Species","Genome_Build","Source_Seq",
                       "Forward_Sequence","Top_Sequence","Rep_Num")
  
  man_fwd_ses_tib <- man_top_sesame_tib %>%
    dplyr::left_join(fwd_sel_tib, 
                     by=c("Probe_ID"="Des_Key", 
                          "U"="Des_AddU", 
                          "M"="Des_AddM", "Strand_CO") ) %>%
    dplyr::mutate(Name=stringr::str_remove(Probe_ID, "_.*$")) %>%
    dplyr::select(dplyr::all_of(man_fwd_ses_col))
  
  #
  # Sesame Controls::
  #
  man_ses_ctl_tib <- 
    readr::read_csv(file.path(par$datDir, "manifest/core/COVIC-C12.manifest.sesame-base.cpg-sorted.csv.gz")) %>%
    dplyr::select(dplyr::any_of(names(man_fwd_ses_tib))) %>% 
    dplyr::filter(stringr::str_starts(Probe_ID, "ct")) %>%
    dplyr::rename(Infinium_Design_Type=Infinium_Design) %>%
    dplyr::mutate(Infinium_Design=as.integer(2))
  
  ses_full_tib <- dplyr::bind_rows(
    man_fwd_ses_tib,man_ses_ctl_tib)
  readr::write_csv(ses_full_tib, run$man_ses_csv)
  
  
  #
  # Genome Studio::
  #
  gs_header_vec <- genome_studio_header(name="Chicago-A5-GRCh37",count=base::nrow(man_fwd_ses_tib))
  gs_controls_tib <- readr::read_csv("/Users/bretbarnes/Documents/data/manifests/raw/manifests/methylation/controls/MethylationEPIC_v-1-0_B4.only-controls-table.csv.gz")
  gs_body_tib <- man_fwd_ses_tib %>% dplyr::rename(
    IlmnID=Probe_ID,
    AddressA_ID=U,
    AddressB_ID=M
  )
  
  readr::write_lines(gs_header_vec, run$man_gsm_csv, append = FALSE)
  readr::write_csv(gs_body_tib, run$man_gsm_csv, append = TRUE, col_names = TRUE)
  readr::write_lines("[Controls],,,,,,,", run$man_gsm_csv, append = TRUE)
  readr::write_csv(gs_controls_tib, run$man_gsm_csv, append = TRUE, col_names = FALSE)
  
  
  if (FALSE) {
    seq_imp_tib %>% dplyr::select(Address, Can_Top) %>% dplyr::filter(!is.na(Can_Top)) %>%
      dplyr::inner_join(exp_tib, by=c("Address"="U"))
    
    # Some quick validation::
    seq_imp_tib %>% dplyr::select(Address, Can_Top) %>% dplyr::filter(!is.na(Can_Top)) %>%
      dplyr::inner_join(exp_tib, by=c("Address"="U")) %>% dplyr::select(Des_Key, Can_Top, Des_Ref_Seq.y) %>% dplyr::mutate(Des_Ref_Seq1=revCmp( shearBrac(Des_Ref_Seq.y))) %>% dplyr::select(Can_Top,Des_Ref_Seq1)
    
  }
}






























# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     5.0 Rebuild Manifest:: ord_tib
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #






# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                             6.0 Add Controls:: 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #






# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        7.0 Convert To Genome Studo:: 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #






# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                 X.X Build improbe Genomic Ranges:: by 50U/49U
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_cgn_mat <- FALSE
if (run_cgn_mat) {
  
  opt$impDir <- "/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput"
  
  imp_prb_dir <- file.path(opt$impDir, "prbDB/merged")
  imp_u49_tsv <- file.path(imp_prb_dir, paste("GRCh36-GRCh37-GRCh38-GRCm10.u49-cgn.seq-sorted.tsv.gz", sep="-") )
  imp_u50_tsv <- file.path(imp_prb_dir, paste("GRCh36-GRCh37-GRCh38-GRCm10.u50-cgn.seq-sorted.tsv.gz", sep="-") )
  
  # imp_u49_tib <- readr::read_tsv(imp_u49_tsv)
  # imp_u50_tib <- readr::read_tsv(imp_u50_tsv)
  
  if (FALSE) {
    # Build/Load imp_bed_grs::
    #  NOTE: imp_bed_tib uses base-0 coordinates
    #
    
    opt$impDir <- "/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput"
    imp_bed_dir <- file.path(opt$impDir, "bed")
    
    gen_vec <- c("GRCh38","GRCh37","GRCh36","GRCm10")
    for (gen in gen_vec) {
      cat(glue::glue("[{par$prgmTag}]: Genome Build={gen}.{RET}"))
      
      imp_bed_tsv <- file.path(imp_bed_dir, paste(gen, "21022021.cgn.sorted.bed.gz", sep="-"))
      imp_grs_rds <- file.path(imp_bed_dir, paste(gen, "21022021.cgn.sorted.grs.rds", sep="-"))
      imp_pos_tsv <- file.path(imp_bed_dir, paste(gen, "21022021.pos.sorted.tsv.gz", sep="-"))
      # imp_pos_rds <- file.path(imp_bed_dir, paste(gen, "21022021.pos.sorted.rds", sep="-"))
      
      cat(glue::glue("[{par$prgmTag}]:{TAB} imp_bed_tsv={imp_bed_tsv}.{RET}"))
      cat(glue::glue("[{par$prgmTag}]:{TAB} imp_pos_tsv={imp_pos_tsv}.{RET}"))
      cat(glue::glue("[{par$prgmTag}]:{TAB} imp_grs_rds={imp_grs_rds}.{RET}"))
      
      if (!file.exists(imp_grs_rds)) {
        
        imp_bed_cols <-
          cols(
            imp_chr = col_character(),
            imp_beg = col_integer(),
            imp_end = col_integer(),
            imp_cgn = col_character(),
            imp_scr = col_character(),
            imp_srd = col_character()
          )
        imp_bed_tib <- 
          readr::read_tsv(imp_bed_tsv, col_names=names(imp_bed_cols$cols), col_types=imp_bed_cols) %>% 
          dplyr::group_by(imp_cgn) %>% 
          dplyr::mutate(imp_rep=dplyr::row_number()) %>% 
          dplyr::ungroup()
        
        imp_bed_grs <-
          GRanges(
            seqnames = Rle(imp_bed_tib$imp_chr),
            strand   = Rle(imp_bed_tib$imp_srd),
            
            cgn=imp_bed_tib$imp_cgn,
            rep=imp_bed_tib$imp_rep,
            
            IRanges(start = imp_bed_tib$imp_beg+1,
                    end   = imp_bed_tib$imp_end+1,
                    names=paste(imp_bed_tib$imp_cgn,imp_bed_tib$imp_rep, sep="_")
            )
          )
        safe_write(imp_bed_grs,"rds",imp_grs_rds, funcTag=par$prgmDir, verbose=opt$verbose)
        
      } else {
        
        if (opt$verbose>=1)
          cat(glue::glue("[{par$prgmTag}]: Loading imp_grs_rds={imp_grs_rds}...{RET}"))
        imp_bed_grs <- readr::read_rds(imp_grs_rds)
        
      }
      
      if (!file.exists(imp_pos_tsv)) {
        
        cgn_bed_tib <- imp_bed_grs %>% # head() %>% 
          as.data.frame() %>% 
          tibble::as_tibble() %>%
          dplyr::rename(Imp_Chr=seqnames, Imp_Pos=start, Imp_Top_Srd=strand, Imp_Cgn=cgn, Imp_Rep=rep) %>%
          dplyr::select(Imp_Chr,Imp_Pos,Imp_Top_Srd,Imp_Cgn,Imp_Rep) %>%
          dplyr::mutate(Imp_Cgn=stringr::str_remove(Imp_Cgn,"^cg") %>% stringr::str_remove("^0+"),
                        Imp_Top_Srd=dplyr::case_when(
                          Imp_Top_Srd=="+" ~ "F",
                          Imp_Top_Srd=="-" ~ "R",
                          TRUE ~ NA_character_)
          ) %>%
          utils::type.convert() %>% 
          dplyr::mutate(across(where(is.factor), as.character) )
        
        print_tib(cgn_bed_tib,par$prgmTag,v=opt$verbose, n=gen)
        safe_write(cgn_bed_tib,"tsv",imp_pos_tsv, funcTag=par$prgmDir, verbose=opt$verbose)
        
      }
    }
    
  }
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                         8.6 Code to be Removed:: 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
