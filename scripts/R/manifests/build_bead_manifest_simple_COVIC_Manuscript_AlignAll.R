
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Genomic Ranges::
suppressWarnings(suppressPackageStartupMessages( base::require("GenomicRanges",quietly=TRUE) ))

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
par$prgmDir <- 'manifests'
par$prgmTag <- 'build_bead_manifest_simple_COVIC_Manuscript_AlignAll'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

par$idat_prefix   <- NULL
par$local_runType <- NULL

# Executables::
opt$Rscript <- NULL

# Run Parameters::
opt$runName   <- NULL

# Directories::
opt$outDir <- NULL
opt$impDir <- NULL
opt$annDir <- NULL
opt$genDir <- NULL
opt$manDir <- NULL

# Required Inputs::
opt$ords <- NULL
opt$mats <- NULL
opt$aqps <- NULL
opt$pqcs <- NULL
opt$ctls <- NULL
opt$idat <- NULL

opt$cpg_s48_tsv <- NULL
opt$cpg_top_tsv <- NULL

opt$cpg_pos_tsv <- NULL
opt$cph_pos_tsv <- NULL
opt$snp_pos_tsv <- NULL

opt$pre_man_csv <- NULL

# opt$cpg_des_csv <- NULL
opt$cph_des_csv <- NULL
opt$snp_des_csv <- NULL

# Platform/Method Options::
opt$genomeBuild <- NULL
opt$platform    <- NULL
opt$version     <- NULL

# Run Options::
opt$fresh   <- FALSE
opt$fixIds  <- FALSE
par$retData <- FALSE

opt$matFormat <- 'new'
opt$ordFormat <- 'old'

opt$matSkip <- 40
opt$ordSkip <- 8
opt$aqpSkip <- 7
opt$pqcSkip <- 7

# Parallel/Cluster Options::
opt$single   <- FALSE
opt$parallel <- TRUE
opt$cluster  <- FALSE

opt$opt_csv  <- NULL
opt$par_csv  <- NULL
opt$time_csv <- NULL
opt$time_org_txt <- NULL

# BSMAP Parameters::
opt$bsmap_exe <- "/Users/bretbarnes/Documents/programs/BSMAPz/bsmapz"
opt$bsmap_opt <- NULL
# opt$bsmap_opt <- "'-s 10 -v 5 -n 1 -r 2 -V 2'"
# opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""

# verbose Options::
opt$verbose <- 3

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
  
  opt$outDir <- file.path(par$topDir, 'scratch')
  locIdatDir <- file.path(par$topDir, 'data/idats')
  opt$impDir <- file.path(par$topDir, 'data/improbe')
  opt$annDir <- file.path(par$topDir, 'data/annotation')
  
  # opt$cpg_s48_tsv <- file.path(opt$impDir, 'designOutput_21092020/seq48U/gz/seq48U-GRCh36-38-10-21092020.unq.noHeader.seq-sorted.tsv.gz')
  opt$cpg_s48_tsv <- file.path(opt$impDir, 'designOutput_21092020/seq48U/un/seq48U-GRCh36-38-10-21092020.unq.noHeader.seq-sorted.tsv')
  
  # The reason we don't use the general form is that specificity for the species may be lost.
  #   I think this is poor form, but using now to meet mouse GS requriments...
  #
  # opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop/GRCh36-GRCh38-GRCm10-21092020.cgnTop.sorted.tsv.gz')
  # opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop/GRCh36-GRCh38-GRCm10-21092020.cgnTop.sorted.tsv')
  
  opt$write_full <- FALSE
  opt$write_base <- FALSE
  
  opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
  opt$bsmap_exe <- "/Users/bretbarnes/Documents/programs/BSMAPz/bsmapz"
  
  opt$manDir    <- "/Users/bretbarnes/Documents/data/manifests"
  opt$genDir    <- "/Users/bretbarnes/Documents/data/iGenomes/Homo_sapiens/NCBI"
  
  #
  # Pre-defined local options runTypes::
  #
  par$local_runType <- 'qcMVP'
  par$local_runType <- 'NZT'
  par$local_runType <- 'COVID'
  par$local_runType <- 'GENK'
  par$local_runType <- 'GRCm38'
  par$local_runType <- 'COVIC'
  
  if (par$local_runType=='COVID') {
    # opt$fresh  <- TRUE
    
    opt$matFormat <- 'old'
    opt$ordFormat <- 'gta'
    opt$matSkip <- 0
    opt$ordSkip <- 15
    
    opt$genomeBuild <- 'COVID'
    opt$platform    <- 'COVID'
    opt$version     <- 'C2'
    
    #
    # TBD:: Add genomic coordinates for virus::
    #
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genomeBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genomeBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/COVID-19_HLA/AQP/COVID-Direct-Detection')
    opt$ords <- paste(
      file.path(par$aqpDir, '371328_CoV_1K_HTS_FinalDesign.design.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20474076_probes.match.gz'),
      sep=',')
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'BS0032777-AQP.txt.gz'),
      sep=',')
    
    opt$pqcs <- NULL
    opt$pqcs <- paste(
      file.path(par$aqpDir, '329922X371395_A_ProductQC.txt.gz'),
      sep=',')
    
    par$idatsTopDir <- file.path(locIdatDir,'idats_COVID-Direct-Set1')
    opt$idat <- paste(
      file.path(par$idatsTopDir, '204756130014'),
      sep=',')
    
  } else if (par$local_runType=='COVIC') {
    opt$matFormat <- 'old'
    opt$ordFormat <- 'old'
    opt$matSkip <- 0
    opt$ordSkip <- 8
    
    opt$genomeBuild <- 'GRCh36'
    opt$genomeBuild <- 'GRCh38'
    opt$genomeBuild <- 'GRCh37'
    
    opt$platform    <- 'COVIC'
    opt$version     <- 'C0'
    opt$version     <- 'B1'
    opt$version     <- 'C2'
    opt$version     <- 'C11'
    opt$version     <- 'C12'
    
    par$idat_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R01C01"
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genomeBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genomeBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    opt$pre_man_csv <- '/Users/bretbarnes/Documents/data/manifests/MethylationEPIC_v-1-0_B2.csv.gz'
    
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
    # opt$pqcs <- paste(
    #   file.path(par$aqpDir, 'BS0032581-AQP.txt.gz'),
    #   sep=',')
    
    par$idatsTopDir <- file.path(locIdatDir,'idats_COVIC-Set1-15052020')
    opt$idat <- paste(
      file.path(par$idatsTopDir, '204500250013'),
      sep=',')
    
  } else if (par$local_runType=='GENK') {
    opt$fresh <- FALSE
    
    opt$matFormat <- 'old'
    opt$ordFormat <- 'old'
    opt$matSkip <- 0
    opt$ordSkip <- 0
    
    opt$genomeBuild <- 'GRCh37'
    opt$genomeBuild <- 'GRCh38'
    opt$platform    <- 'GENK'
    opt$version     <- 'A1'
    opt$version     <- 'A2'
    opt$version     <- 'A3'
    opt$version     <- 'B2'
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genomeBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genomeBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/Genknowme/LS_Epiprofile')
    opt$ords <- paste(
      file.path(par$aqpDir, 'AQP1_NAremoved_GenKnowme_CpG_SNP_order.07082020.csv'),
      file.path(par$aqpDir, 'AQP2_AP_Genknowme_AQP2_replicate_design_file2.csv'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20468029_AQP1_probes.match'),
      file.path(par$aqpDir, '20468029_AQP2_probes.match'),
      sep=',')
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'BS0032678_AQP1-AQP.txt'),
      file.path(par$aqpDir, 'BS0032779_AQP2-AQP.txt'),
      sep=',')
    
    opt$pqcs <- paste(
      file.path(par$aqpDir, '20042793X371678_A_ProductQC_AP.txt'),
      sep=',')
    
    opt$idat <- NULL
  } else if (par$local_runType=='GRCm38') {
    opt$genomeBuild <- 'GRCm38'
    opt$platform    <- 'LEGX'
    
    opt$version     <- 'C0'
    opt$version     <- 'C8'
    opt$version     <- 'C25'
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genomeBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genomeBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    opt$cph_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genomeBuild,'.improbeDesignInput.chn-sorted.tsv.gz') )
    opt$snp_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genomeBuild,'.improbeDesignInput.snp-sorted.tsv.gz') )
    
    opt$snp_des_csv <- file.path(opt$impDir, 'cph-snp-designs/LEGX_SpikeIn_Reorder-SNP-Only.designs.csv.gz')
    opt$cph_des_csv <- file.path(opt$impDir, 'cph-snp-designs/LEGX_SpikeIn_Reorder-CpH-Only.designs.csv.gz')
    
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
    
    par$idatsTopDir <- file.path(locIdatDir,'idats_ILMN_mm10_betaTest_17082020')
    opt$idat <- paste(
      file.path(par$idatsTopDir, '204637490025'),
      sep=',')
    
  } else if (par$local_runType=='NZT') {
    opt$matFormat <- 'old'
    opt$ordFormat <- 'old'
    opt$matSkip <- 0
    opt$ordSkip <- 8
    opt$pqcSkip <- 7
    
    opt$genomeBuild <- 'GRCh36'
    opt$genomeBuild <- 'GRCh38'
    opt$genomeBuild <- 'GRCh37'
    opt$platform    <- 'NZT'
    opt$version     <- 'N0'
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genomeBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genomeBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
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
    
    opt$pqcs <- NULL
    opt$idat <- NULL
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$parallel <- TRUE
  # opt$runName <- paste(opt$genomeBuild,opt$platform,opt$version, sep='-')
  opt$runName <- paste(opt$platform,opt$version, sep='-')
  
  # opt$fresh <- TRUE
  opt$fresh <- FALSE
  
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
    
    make_option(c("--bsmap_exe"), type="character", default=opt$bsmap_exe, 
                help="BSMAP Executable path [default= %default]", metavar="character"),
    make_option(c("--bsmap_opt"), type="character", default=opt$bsmap_opt, 
                help="BSMAP Options [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--pre_man_csv"), type="character", default=opt$pre_man_csv, 
                help="Previously defined manifest [default= %default]", metavar="character"),
    
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
    
    # Pre-defined files (controls)
    make_option(c("--ords"), type="character", default=opt$ords, 
                help="Order files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mats"), type="character", default=opt$mats, 
                help="Match files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqps"), type="character", default=opt$aqps, 
                help="AQP files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--pqcs"), type="character", default=opt$pqcs, 
                help="PQC files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--ctls"), type="character", default=opt$ctls, 
                help="Pre-Defined Infinium Methylation Controls (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--idat"), type="character", default=opt$idat, 
                help="idat directories (comma seperated) [default= %default]", metavar="character"),
    
    make_option(c("--cpg_s48_tsv"), type="character", default=opt$cpg_s48_tsv, 
                help="Seq48U Match file(s) (tab seperated) [default= %default]", metavar="character"),
    make_option(c("--cpg_top_tsv"), type="character", default=opt$cpg_top_tsv, 
                help="Top Sequence CGN Match file(s) (tab seperated) [default= %default]", metavar="character"),
    
    make_option(c("--cpg_pos_tsv"), type="character", default=opt$cpg_pos_tsv, 
                help="CpG Position TSV File [default= %default]", metavar="character"),
    make_option(c("--cph_pos_tsv"), type="character", default=opt$cph_pos_tsv, 
                help="CpH Position TSV File [default= %default]", metavar="character"),
    make_option(c("--snp_pos_tsv"), type="character", default=opt$snp_pos_tsv, 
                help="SNP Position TSV File [default= %default]", metavar="character"),
    
    # make_option(c("--cpg_des_csv"), type="character", default=opt$cpg_des_csv, 
    #             help="CpG Position TSV File [default= %default]", metavar="character"),
    make_option(c("--cph_des_csv"), type="character", default=opt$cph_des_csv, 
                help="CpH Position TSV File [default= %default]", metavar="character"),
    make_option(c("--snp_des_csv"), type="character", default=opt$snp_des_csv, 
                help="SNP Position TSV File [default= %default]", metavar="character"),
    
    # Platform/Method Options::
    make_option(c("--genomeBuild"), type="character", default=opt$genomeBuild, 
                help="Genome Build (e.g. GRCh36, GRCh37, GRCh38, GRCm38) [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform (e.g. HM450, EPIC, LEGX, NZT, COVIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest Version (e.g. B0,B1,B2,B3,B4,C0) [default= %default]", metavar="character"),
    
    # File formats
    make_option(c("--ordFormat"), type="character", default=opt$ordFormat, 
                help="Order File Format [default= %default]", metavar="character"),
    make_option(c("--matFormat"), type="character", default=opt$matFormat, 
                help="Match File Format [default= %default]", metavar="character"),
    
    # Skip Lines
    make_option(c("--matSkip"), type="integer", default=opt$matSkip, 
                help="Match Skip Lines Count [default= %default]", metavar="integer"),
    make_option(c("--ordSkip"), type="integer", default=opt$ordSkip, 
                help="Order Skip Lines Count [default= %default]", metavar="integer"),
    make_option(c("--aqpSkip"), type="integer", default=opt$aqpSkip, 
                help="AQP Skip Lines Count [default= %default]", metavar="integer"),
    make_option(c("--pqcSkip"), type="integer", default=opt$pqcSkip, 
                help="PQC Skip Lines Count [default= %default]", metavar="integer"),
    
    # Run Options::
    make_option(c("--fresh"), action="store_true", default=opt$fresh, 
                help="Boolean variable to run a fresh build [default= %default]", metavar="boolean"),
    make_option(c("--fixIds"), action="store_true", default=opt$fixIds, 
                help="Boolean variable to fix original order ids [default= %default]", metavar="boolean"),
    make_option(c("--retData"), action="store_true", default=opt$retData, 
                help="Developement ONLY Boolean variable to return data for testing [default= %default]", metavar="boolean"),
    
    # Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
    make_option(c("--opt_csv"), type="character", default=opt$opt_csv, 
                help="Unused variable opt_csv [default= %default]", metavar="character"),
    make_option(c("--par_csv"), type="character", default=opt$par_csv, 
                help="Unused variable par_csv [default= %default]", metavar="character"),
    make_option(c("--time_csv"), type="character", default=opt$time_csv, 
                help="Unused variable time_csv [default= %default]", metavar="character"),
    make_option(c("--time_org_txt"), type="character", default=opt$time_org_txt, 
                help="Unused variable time_org_txt [default= %default]", metavar="character"),
    
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
opt_reqs <- c('outDir','impDir','ords','mats','cpg_s48_tsv','cpg_top_tsv',
              'genomeBuild','platform','version','bsmap_exe', # 'bsmap_opt',
              'Rscript','verbose')

par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )

opt <- program_init(name=par$prgmTag,
                    opts=opt, opt_reqs=opt_reqs, 
                    pars=par, par_reqs=par_reqs,
                    libs=TRUE,rcpp=FALSE,
                    verbose=opt$verbose,vt=3,tc=0,tt=NULL)

par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

image_key <- "bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse_Centos"
image_ver <- "v.1.0"
image_ssh <- "run_improbe.sh"
image_str <- glue::glue("{image_key}.{image_ver}")

# Input Definitions::
# if (is.null(opt$ctls)) {
#   opt$ctls <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
# }

ords_vec <- NULL
mats_vec <- NULL
aqps_vec <- NULL
pqcs_vec <- NULL
ctls_vec <- NULL
idat_vec <- NULL
if (!is.null(opt$ords)) ords_vec <- opt$ords %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$mats)) mats_vec <- opt$mats %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$aqps)) aqps_vec <- opt$aqps %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$pqcs)) pqcs_vec <- opt$pqcs %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$ctls)) ctls_vec <- opt$ctls %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$idat)) idat_vec <- opt$idat %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

stopifnot(length(ords_vec)>0)
stopifnot(length(mats_vec)>0)
stopifnot(length(mats_vec)==length(ords_vec))

cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

# Output Definitions::

opt$preDir <- file.path(opt$outDir, 'manifest')
if (!dir.exists(opt$preDir)) dir.create(opt$preDir, recursive=TRUE)

# opt$prdDir <- file.path(opt$outDir, 'product')
# if (!dir.exists(opt$prdDir)) dir.create(opt$prdDir, recursive=TRUE)
# 
# opt$desDir <- file.path(opt$outDir, 'design')
# if (!dir.exists(opt$desDir)) dir.create(opt$desDir, recursive=TRUE)
# 
opt$intDir <- file.path(opt$outDir, 'intersection')
if (!dir.exists(opt$intDir)) dir.create(opt$intDir, recursive=TRUE)

# opt$sumDir <- file.path(opt$outDir, 'summary')
# if (!dir.exists(opt$sumDir)) dir.create(opt$sumDir, recursive=TRUE)

opt$bspDir <- file.path(opt$outDir, 'bspmap')
if (!dir.exists(opt$bspDir)) dir.create(opt$bspDir, recursive=TRUE)

# Define Alignment Genomes::
#
man_fas_pre <- file.path(opt$outDir, "fas", paste(opt$runName,"aln-prbs", sep="_"))
man_prb_fas <- paste0(man_fas_pre,".fa.gz")
man_gen_fas <- file.path(opt$genDir, opt$genomeBuild, "Sequence/WholeGenomeFasta", 
                         paste0(opt$genomeBuild,".genome.fa.gz"))

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Building Output Directories.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Data Collection::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      1.1 Load Pre-built Manifests::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


opt$add_core <- TRUE
if (opt$add_core) {
  man_raw_cor_tib <- 
    loadCoreManifest_COVIC(
      datDir=par$datDir, manDir=opt$manDir,
      verbose=opt$verbose,tc=0,tt=pTracker) %>% 
    dplyr::select(Probe_ID,U,M,Probe_Type,Probe_Type,
                  AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, 
                  chrom_GRCh37,chromStart_GRCh37,chromEnd_GRCh37,FR_GRCh37,
                  nextBase_GRCh37,gene_GRCh37,gene_HGNC_GRCh37,
                  MASK_mapping_GRCh37,MASK_general_GRCh37,Manifest) %>% 
    dplyr::rename(chrom=chrom_GRCh37,chromStart=chromStart_GRCh37,
                  chromEnd=chromEnd_GRCh37,strand_FR=FR_GRCh37,
                  nextBase=nextBase_GRCh37,
                  gene=gene_GRCh37,gene_HGNC=gene_HGNC_GRCh37,
                  MASK_mapping=MASK_mapping_GRCh37,
                  MASK_general=MASK_general_GRCh37) %>%
    dplyr::mutate(strand_FR=dplyr::case_when(
      strand_FR=="F" ~ "+",
      strand_FR=="R" ~ "-",
      TRUE ~ NA_character_),
      chromStart=as.integer(chromStart),
      chromEnd=as.integer(chromEnd),
      nextBase=dplyr::case_when(
        is.na(M) ~ "Both",
        TRUE ~ nextBase)
    )
}

# man_raw_cor_tib %>% head() %>% 
#   manifestToAddress_COVIC(verbose = 3)


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      1.2 New AQP Designs/Results::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_raw_cov_tib <- 
  loadAqpWorkflow_COVIC(
    ords=opt$ords, mats=opt$mats, aqps=opt$aqps,man="C0",
    verbose=opt$verbose,tc=0,tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         1.3 Extract New Overlaps::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# 1. How many man_raw_cor_tib == man_raw_cov_tib by="AlleleA_Probe_Sequence"
#

# man_raw_dup_tib <- man_raw_cov_tib %>%
#   dplyr::filter(AlleleA_Probe_Sequence %in% man_raw_cor_tib$AlleleA_Probe_Sequence) %>%
#   dplyr::filter(AlleleB_Probe_Sequence %in% man_raw_cor_tib$AlleleB_Probe_Sequence) %>%
#   dplyr::select(dplyr::any_of(man_raw_cor_tib %>% names()))
# 
# man_raw_dup_tib <- man_raw_dup_tib %>%
#   dplyr::left_join(
#     dplyr::select(man_raw_cor_tib, AlleleA_Probe_Sequence,!dplyr::any_of(man_raw_dup_tib %>% names())),
#     by="AlleleA_Probe_Sequence") %>%
#   dplyr::select(dplyr::any_of(man_raw_cor_tib %>% names()))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   1.4.0 Build Pre-built Annotation GRS::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# EPIC: B2/B4/C0-Dups
#
# man_sel_all_tib <- 
#   dplyr::bind_rows(man_raw_cor_tib,man_raw_dup_tib) %>%
#   dplyr::rename(Seq_ID=Probe_ID)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 1.5 Extract Novel New Sequences for Alignment::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# man_raw_new_tib <- 
#   dplyr::anti_join(man_raw_cov_tib,man_raw_cor_tib, 
#                  by=c("AlleleA_Probe_Sequence"))

#
# Full Data Set::
#
man_raw_all_tib <- 
  dplyr::bind_rows(
    man_raw_cor_tib %>% dplyr::select(dplyr::any_of(names(man_raw_cov_tib))),
    man_raw_cov_tib
  ) %>% 
  dplyr::select(-Seq_48U) %>%
  dplyr::distinct(U,M, .keep_all=TRUE)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        1.6 Build Address Tibbles::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Core (cor=core, i.e. EPIC B2/B4) Sesame in Address Format::
add_raw_ses_tib <- 
  manifestToAddress_COVIC(
    man_raw_cor_tib, verbose=opt$verbose,tc=0,tt=pTracker) %>%
  dplyr::rename(ses_add=prb_add, ses_cgn=prb_cgn, ses_des=prb_des,
                ses_typ=Probe_Type, prb_ord_seq=ord_seq,
                ses_chr=chrom, ses_pos=chromStart, ses_end=chromEnd, 
                ses_fr=strand_FR, ses_nxb=nextBase, ses_gene=gene, 
                ses_HGNC=gene_HGNC, ses_src=prb_src) %>%
  dplyr::select(ses_add, ses_cgn, ses_des, ses_typ,
                ses_chr, ses_pos, ses_fr, ses_nxb, prb_ord_seq, 
                ses_gene, ses_HGNC, ses_src, MASK_mapping, MASK_general) %>%
  dplyr::mutate(ses_chr=dplyr::case_when(
    ses_chr=="chrM" ~ "chrMT",
    TRUE ~ ses_chr)
  )

#
# All Data (cor=core && cov=covic) in Address Format::
#
add_raw_all_tib <- 
  manifestToAddress_COVIC(
    man_raw_all_tib, verbose=opt$verbose,tc=0,tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      1.7 Write Addresss Fasta Files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_fas_tib <- tibToFas(tib=add_raw_all_tib, key="aln_key", 
                        seq="aln_seq", prefix=man_fas_pre, 
                        verbose=opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       STEP 2:: Alignment:: BSMAAP
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_add_bsp_tsv <- file.path(
  opt$bspDir,"align", paste0(paste(opt$runName,"aln-prbs", sep="_"),"-",
                             opt$genomeBuild,".genome.bsmap.formatted.tsv.gz"))
man_add_bsp_rds <- file.path(
  opt$bspDir,"align", paste0(paste(opt$runName,"aln-prbs", sep="_"),"-",
                             opt$genomeBuild,".genome.bsmap.formatted.rds"))

if (!file.exists(man_add_bsp_tsv)) {
  man_add_bsp_tsv <- bsmapProbeAlign(
    exe=opt$bsmap_exe, fas=man_prb_fas, gen=man_gen_fas, 
    dir=opt$bspDir, opt=opt$bsmap_opt, run=TRUE,
    verbose=opt$verbose,tt=pTracker)
}

if (!file.exists(man_add_bsp_tsv) |
    !file.exists(man_add_bsp_rds) |
    file.mtime(man_add_bsp_rds) < file.mtime(man_add_bsp_tsv) ) {
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading man_add_bsp (TSV)={man_add_bsp_tsv}...{RET}"))
  
  man_add_bsp_tib <- loadBspFormatted_COVIC(
    bsp=man_add_bsp_tsv, src=man_fas_tib, sort=TRUE,
    verbose=opt$verbose,tt=pTracker)
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Building man_add_bsp (RDS)={man_add_bsp_rds}...{RET}"))
  
  man_add_bsp_grs <-bspToGenomicRegion_COVIC(
    bsp=man_add_bsp_tib, rds=man_add_bsp_rds,
    verbose=opt$verbose,tt=pTracker)
  
} else {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading man_add_bsp (RDS)={man_add_bsp_rds}...{RET}"))
  man_add_bsp_grs <- readr::read_rds(man_add_bsp_rds)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#             STEP 3:: Load cg# Database for Genome Build
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Standard Static Inputs::
cgn_pos_db_rds  <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-end-seq.cgn-sorted.rds"
cgn_pos_db_csv  <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-seq.cgn-sorted.csv"

# Run Time Outputs::
add_imp_ses_rds <- file.path(opt$intDir, paste(opt$runName,opt$genomeBuild,"sesame-preproduct.address.GenomicRegions.rds", sep="_"))
add_imp_all_rds <- file.path(opt$intDir, paste(opt$runName,opt$genomeBuild,"improbe-join.pass-address.GenomicRegions.rds", sep="_"))
add_imp_all_csv <- file.path(opt$intDir, paste(opt$runName,opt$genomeBuild,"improbe-join.csv.gz", sep="_"))

opt$writeStep3 <- FALSE # Just a safety during testing to not delete large files...
opt$writeStep4 <- TRUE

if (!file.exists(add_imp_all_rds)) {
  if (!file.exists(add_imp_all_csv)) {
    
    # cgn_pos_db2_csv <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-end-seq.cgn-sorted.csv.gz"
    # cgn_pos_db_csv  <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-seq.cgn-sorted.csv"
    # cgn_pos_db_rds  <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-prb.GenomicRegions.rds"
    
    if (file.exists(cgn_pos_db_rds)) {
      if (opt$verbose>=1)
        cat(glue::glue("[{par$prgmTag}]: Loading GRS RDS={cgn_pos_db_rds}...{RET}"))
      
      cgn_pos_grs <- readr::read_rds(cgn_pos_db_rds)
      
    } else if (file.exists(cgn_pos_db2_csv)) {
      if (opt$verbose>=1)
        cat(glue::glue("[{par$prgmTag}]: Loading Pre-processed CGN CSV={cgn_pos_db2_csv}...{RET}"))
      
      cgn_pos_db_tib <- readr::read_csv(cgn_pos_db2_csv)
      
      # Next: Create GRS and save RDS for future use...
      #
      cgn_pos_grs <- GRanges(
        seqnames = Rle(cgn_pos_db_tib$imp_chr),
        strand=Rle(cgn_pos_db_tib$imp_frs),
        imp_seq1U=cgn_pos_db_tib$imp_seq,
        IRanges(start=cgn_pos_db_tib$imp_pos,
                end=cgn_pos_db_tib$imp_end,
                names=cgn_pos_db_tib$Probe_ID_IMP) )
      
      if (opt$verbose>=1)
        cat(glue::glue("[{par$prgmTag}]: Writing GRS RDS={cgn_pos_db_rds}...{RET}"))
      
      if (opt$writeStep3) readr::write_rds(cgn_pos_grs, cgn_pos_db_rds, compress="gz")    
      
    } else if (file.exists(cgn_pos_db_csv)) {
      
      if (opt$verbose>=1)
        cat(glue::glue("[{par$prgmTag}]: Loading Original CGN CSV={cgn_pos_db_csv}...{RET}"))
      
      cgn_pos_db_col <- cols(
        imp_cgn = col_character(),
        imp_chr = col_character(),
        imp_pos = col_integer(),
        
        imp_frs = col_character(),
        imp_srd = col_character(),
        imp_seq = col_character()
      )
      cgn_pos_db_tib <- readr::read_csv(cgn_pos_db_csv, 
                                        col_names=names(cgn_pos_db_col$cols), 
                                        col_types=cgn_pos_db_col)
      
      # NEW METHOD:: To Be Removed::
      if (FALSE) {
        cgn_pos_grs_tib <- cgn_pos_db_tib %>% 
          dplyr::mutate(imp_frs=stringr::str_replace(imp_frs, "\\+","F") %>% stringr::str_replace("-","R")) %>%
          dplyr::rename(cgn=imp_cgn,chr=imp_chr,pos=imp_pos,FR=imp_frs,prb1U=imp_seq) %>%
          tidyr::separate(imp_srd, into=c("srd", "nxb"), sep=c(2)) %>% 
          tidyr::pivot_wider(names_from="srd", values_from=c("prb1U", "FR","nxb")) %>%
          dplyr::mutate(prb_id=paste(cgn,chr,pos, sep="_"))
        
        cgn_pos_grs <- GRanges(
          seqnames = Rle(cgn_pos_grs_tib$chr),
          # strand=Rle(cgn_pos_db_tib$imp_frs),
          
          prb1U_TC=cgn_pos_grs_tib$prb1U_TC,
          prb1U_TO=cgn_pos_grs_tib$prb1U_TO,
          prb1U_BC=cgn_pos_grs_tib$prb1U_BC,
          prb1U_BO=cgn_pos_grs_tib$prb1U_BO,
          
          FR_TC=cgn_pos_grs_tib$FR_TC,
          FR_TO=cgn_pos_grs_tib$FR_TO,
          FR_BC=cgn_pos_grs_tib$FR_BC,
          FR_BO=cgn_pos_grs_tib$FR_BO,
          
          nxb_TC=cgn_pos_grs_tib$nxb_TC,
          nxb_TO=cgn_pos_grs_tib$nxb_TO,
          nxb_BC=cgn_pos_grs_tib$nxb_BC,
          nxb_BO=cgn_pos_grs_tib$nxb_BO,
          
          IRanges(start=cgn_pos_grs_tib$pos,
                  end=cgn_pos_grs_tib$pos+1,
                  names=cgn_pos_grs_tib$prb_id)
        )
        
        if (opt$verbose>=1)
          cat(glue::glue("[{par$prgmTag}]: Loading GRS RDS={cgn_pos_db_rds}...{RET}"))
        if (opt$writeStep3) readr::write_rds(cgn_pos_grs, cgn_pos_db_rds, compress="gz")    
      }
      
      # Old Method:: To Be Replaced...
      if (FALSE) {
        cgn_pos_db_tib <- cgn_pos_db_tib %>% 
          dplyr::mutate(imp_end=imp_pos+1,
                        Probe_ID_IMP=paste(imp_cgn,imp_srd, imp_chr,imp_pos, sep="_") )
        
        # Pre-processed for faster future loading...
        if (opt$verbose>=1)
          cat(glue::glue("[{par$prgmTag}]: Writing Pre-processed CGN CSV={cgn_pos_db2_csv}...{RET}"))
        
        if (opt$writeStep3) readr::write_csv(cgn_pos_db_tib, cgn_pos_db2_csv)
        
        # Next: Create GRS and save RDS for future use...
        #
        cgn_pos_grs <- GRanges(
          seqnames = Rle(cgn_pos_db_tib$imp_chr),
          strand=Rle(cgn_pos_db_tib$imp_frs),
          imp_seq1U=cgn_pos_db_tib$imp_seq,
          IRanges(start=cgn_pos_db_tib$imp_pos,
                  end=cgn_pos_db_tib$imp_end,
                  names=cgn_pos_db_tib$Probe_ID_IMP) )
        
        if (opt$verbose>=1)
          cat(glue::glue("[{par$prgmTag}]: Loading GRS RDS={cgn_pos_db_rds}...{RET}"))
        if (opt$writeStep3) readr::write_rds(cgn_pos_grs, cgn_pos_db_rds, compress="gz")    
      }
      
      
    } else {
      stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: File does not exist: {cgn_pos_db_csv}.{RET}{RET}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Intersect probes vs. database::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # To Be Removed::
    #
    if (FALSE) {
      int_tib2 <- intersectGranges_COVIC(
        can=head(man_add_bsp_grs, n=100000),
        ref=head(cgn_pos_grs, n=100000),
        verbose=10
      )
      
      int_mat_tib <- int_tib %>% dplyr::mutate(
        imp_prbTC=dplyr::case_when(
          prb_des_CAN=="2" ~ stringr::str_remove(prb1U_TC_REF, "[A-Za-z]$"),
          TRUE ~ prb1U_TC_REF
        ),
        imp_prbTO=dplyr::case_when(
          prb_des_CAN=="2" ~ stringr::str_remove(prb1U_TO_REF, "[A-Za-z]$"),
          TRUE ~ prb1U_TO_REF
        ),
        imp_prbBC=dplyr::case_when(
          prb_des_CAN=="2" ~ stringr::str_remove(prb1U_BC_REF, "[A-Za-z]$"),
          TRUE ~ prb1U_BC_REF
        ),
        imp_prbBO=dplyr::case_when(
          prb_des_CAN=="2" ~ stringr::str_remove(prb1U_BO_REF, "[A-Za-z]$"),
          TRUE ~ prb1U_BO_REF
        ),
        imp_mat_srd=dplyr::case_when(
          prb_des_CAN=="U" & prb_aln_50U_CAN==imp_prbTC ~ "TC",
          prb_des_CAN=="U" & prb_aln_50U_CAN==imp_prbTO ~ "TO",
          prb_des_CAN=="U" & prb_aln_50U_CAN==imp_prbBC ~ "BC",
          prb_des_CAN=="U" & prb_aln_50U_CAN==imp_prbBO ~ "BO",
          
          prb_des_CAN=="M" & prb_aln_50M_CAN==imp_prbTC ~ "TC",
          prb_des_CAN=="M" & prb_aln_50M_CAN==imp_prbTO ~ "TO",
          prb_des_CAN=="M" & prb_aln_50M_CAN==imp_prbBC ~ "BC",
          prb_des_CAN=="M" & prb_aln_50M_CAN==imp_prbBO ~ "BO",
          
          prb_des_CAN=="2" & prb_aln_49U_CAN==imp_prbTC ~ "TC",
          prb_des_CAN=="2" & prb_aln_49U_CAN==imp_prbTO ~ "TO",
          prb_des_CAN=="2" & prb_aln_49U_CAN==imp_prbBC ~ "BC",
          prb_des_CAN=="2" & prb_aln_49U_CAN==imp_prbBO ~ "BO",
          TRUE ~ NA_character_)
      )
    }
    
    #
    # TBD: Will implement mismatch score first...
    # TBD: Load both U/M prb50 for comparisons...
    #
    # imp_aln_scr = [0,1,2, 3] = [ 50U, 50M, 49U, NA ]
    # imp_cgn_scr = [0,1] = [ Match, NA ]
    #
    add_imp_all_tib <- 
      intersectGranges(
        man=man_add_bsp_grs, ref=cgn_pos_grs,
        verbose=opt$verbose, tt=pTracker) %>%
      dplyr::mutate(
        imp_seq2U=stringr::str_remove(imp_seq1U, "[A-Za-z]$"),
        imp_aln_scr=dplyr::case_when(
          prb_aln_50U==imp_seq1U ~ 0,
          prb_aln_50M==imp_seq1U ~ 1,
          prb_aln_49U==imp_seq2U ~ 2,
          TRUE ~ 3) %>% as.integer()
      ) %>%
      # - Filter alignments with sequence mismatch (imp_aln_scr!=3)
      dplyr::filter(imp_aln_scr!=3) %>%
      # - Split imp_srd into TB/CO/Nxb and compare to strand,imp_frs,bsp_srd
      tidyr::separate(Gene, into=c("imp_cgn","imp_srd","imp_chr2","imp_pos2"), sep="_") %>%
      tidyr::separate(imp_srd, into=c("imp_tb", "imp_co", "imp_nxb"), sep=c(1,2)) %>%
      dplyr::mutate(
        chromStrand=chromStrand %>%
          stringr::str_replace("\\+","F") %>% 
          stringr::str_replace("-","R"),
        imp_cgn_scr=dplyr::case_when(
          prb_cgn == imp_cgn ~ 0,
          TRUE ~ 1) %>% as.integer()
      ) %>%
      dplyr::rename(imp_fr=chromStrand)
    
    add_imp_all_sum <- 
      add_imp_all_tib %>%
      dplyr::select(prb_src,prb_des,ends_with("scr")) %>%
      dplyr::group_by_all() %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    
    if (opt$verbose>=1)
      add_imp_all_sum %>% print(n=base::nrow(add_imp_all_sum))
    
    if (opt$writeStep4) readr::write_csv(add_imp_all_tib, add_imp_all_csv)
  } else {
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Loading Prebuilt Manifest/",
                     "Improbe Database Intersection={add_imp_all_csv}...{RET}"))
    add_imp_all_tib <- readr::read_csv(add_imp_all_csv)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Summary BSP Alignments::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  add_imp_all_sum <- add_imp_all_tib %>% 
    dplyr::arrange(imp_aln_scr, bsp_din_scr) %>% 
    dplyr::distinct(prb_add, .keep_all=TRUE) %>%
    dplyr::group_by(prb_src,prb_des,imp_aln_scr,bsp_din_scr,bsp_srd) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  if (opt$verbose>=1)
    add_imp_all_sum %>% print(n=base::nrow(add_imp_all_sum))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Add All Sesame and Previous Manifest Data::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Summary:: Performing Join Later, so no need for this other than a qucik check...
  if (FALSE) {
    add_imp_all_tib %>% dplyr::anti_join(add_raw_ses_tib, by=c("prb_add"="ses_add")) %>% 
      dplyr::group_by(prb_src) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print()
    
    add_imp_all_tib %>% dplyr::anti_join(add_raw_ses_tib, by="prb_ord_seq") %>% 
      dplyr::group_by(prb_src) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print()
    
    add_raw_ses_tib %>% dplyr::anti_join(add_imp_all_tib, by=c("ses_add"="prb_add")) %>% 
      dplyr::group_by(ses_typ,ses_src) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print()
    
    add_raw_ses_tib %>% dplyr::anti_join(add_imp_all_tib, by="prb_ord_seq") %>% 
      dplyr::group_by(ses_typ,ses_src) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print()
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                Build Genome Regions:: Man/imp/Addresses
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  #
  # Sesame/Previous Manifest GRS::
  #
  add_imp_ses_grs <- 
    GRanges(
      seqnames = Rle(add_raw_ses_tib$ses_chr),
      strand=Rle(add_raw_ses_tib$ses_fr),
      
      # improbe cg# database::
      imp_cgn=add_raw_ses_tib$ses_cgn,
      
      imp_fr=add_raw_ses_tib$ses_fr %>% 
        stringr::str_replace("\\+", "F") %>%
        stringr::str_replace("-", "R"),
      imp_tb=NA_character_,
      imp_co=NA_character_, # Could be determined by probe type ch=O, cg/rs=C
      
      imp_chr=add_raw_ses_tib$ses_chr,
      imp_pos=add_raw_ses_tib$ses_pos+1,
      imp_end=add_raw_ses_tib$ses_pos+1,
      imp_nxb=add_raw_ses_tib$ses_nxb,
      imp_src=add_raw_ses_tib$ses_src,
      
      # probe manifest data::
      prb_add=add_raw_ses_tib$ses_add,
      prb_cgn=add_raw_ses_tib$ses_cgn,
      prb_srd=add_raw_ses_tib$ses_fr,  # Same as imp_co above
      prb_des=add_raw_ses_tib$ses_des,
      prb_src=add_raw_ses_tib$ses_src,
      prb_ord=add_raw_ses_tib$prb_ord_seq,
      
      # bsmap alignment data::
      # bsp_tag=add_raw_ses_tib$bsp_tag,
      # bsp_srd=add_raw_ses_tib$bsp_srd,
      # bsp_mis=add_raw_ses_tib$bsp_mis,
      # bsp_gap=add_raw_ses_tib$bsp_gap,
      # 
      # bsp_din_scr=add_raw_ses_tib$bsp_din_scr,
      # bsp_din_ref=add_raw_ses_tib$bsp_din_ref,
      # bsp_din_bsc=add_raw_ses_tib$bsp_din_bsc,
      # bsp_nxb_ref=add_raw_ses_tib$bsp_nxb_ref,
      # bsp_nxb_bsc=add_raw_ses_tib$bsp_nxb_bsc,
      # bsp_ref_seq=add_raw_ses_tib$bsp_ref_seq,
      # 
      # bsp_hit0=add_raw_ses_tib$bsp_hit0,
      # bsp_hit1=add_raw_ses_tib$bsp_hit1,
      # bsp_hit2=add_raw_ses_tib$bsp_hit2,
      # bsp_hit3=add_raw_ses_tib$bsp_hit3,
      # bsp_hit4=add_raw_ses_tib$bsp_hit4,
      # bsp_hit5=add_raw_ses_tib$bsp_hit5,
      
      # Sesame/Previous Data should be handled
      #  in downstream analysis...
      #
      # ses_add=add_raw_ses_tib$ses_add,
      # ses_cgn=add_raw_ses_tib$ses_cgn,
      ses_typ=add_raw_ses_tib$ses_typ,
      # ses_chr=add_raw_ses_tib$ses_chr,
      # ses_pos=add_raw_ses_tib$ses_pos,
      # ses_fr=add_raw_ses_tib$ses_fr,
      # ses_nxb=add_raw_ses_tib$ses_nxb,
      # prb_ord_seq=add_raw_ses_tib$prb_ord_seq,
      # ses_gene=add_raw_ses_tib$ses_gene,
      # ses_HGNC=add_raw_ses_tib$ses_HGNC,
      # ses_src=add_raw_ses_tib$ses_src,
      MASK_mapping=add_raw_ses_tib$MASK_mapping,
      MASK_general=add_raw_ses_tib$MASK_general,
      
      IRanges(start=add_raw_ses_tib$ses_pos,
              end=add_raw_ses_tib$ses_pos+1,
              names=paste(add_raw_ses_tib$ses_cgn,"xxxx",
                          add_raw_ses_tib$ses_add,
                          add_raw_ses_tib$ses_des,
                          add_raw_ses_tib$ses_src, sep="_")
      )
    )
  readr::write_rds(add_imp_ses_grs, add_imp_ses_rds, compress="gz")
  
  #
  # All GRS::
  #
  add_imp_all_grs <- 
    GRanges(
      seqnames = Rle(add_imp_all_tib$seqnames),
      strand=Rle(add_imp_all_tib$strand),
      
      # improbe cg# database::
      imp_cgn=add_imp_all_tib$imp_cgn,
      
      imp_fr=add_imp_all_tib$imp_fr,
      imp_tb=add_imp_all_tib$imp_tb,
      imp_co=add_imp_all_tib$imp_co,
      
      imp_chr=add_imp_all_tib$chrom,
      imp_pos=add_imp_all_tib$chromStart,
      imp_end=add_imp_all_tib$chromEnd,
      imp_nxb=add_imp_all_tib$imp_nxb,
      imp_aln_scr=add_imp_all_tib$imp_aln_scr,
      imp_cgn_scr=add_imp_all_tib$imp_cgn_scr,
      
      # probe manifest data::
      prb_add=add_imp_all_tib$prb_add,
      prb_cgn=add_imp_all_tib$prb_cgn,
      prb_srd=add_imp_all_tib$prb_srd,
      prb_des=add_imp_all_tib$prb_des,
      prb_src=add_imp_all_tib$prb_src,
      prb_ord=add_imp_all_tib$prb_ord_seq,
      
      # bsmap alignment data::
      bsp_tag=add_imp_all_tib$bsp_tag,
      bsp_srd=add_imp_all_tib$bsp_srd,
      bsp_mis=add_imp_all_tib$bsp_mis,
      bsp_gap=add_imp_all_tib$bsp_gap,
      
      bsp_din_scr=add_imp_all_tib$bsp_din_scr,
      bsp_din_ref=add_imp_all_tib$bsp_din_ref,
      bsp_din_bsc=add_imp_all_tib$bsp_din_bsc,
      bsp_nxb_ref=add_imp_all_tib$bsp_nxb_ref,
      bsp_nxb_bsc=add_imp_all_tib$bsp_nxb_bsc,
      bsp_ref_seq=add_imp_all_tib$bsp_ref_seq,
      
      bsp_hit0=add_imp_all_tib$bsp_hit0,
      bsp_hit1=add_imp_all_tib$bsp_hit1,
      bsp_hit2=add_imp_all_tib$bsp_hit2,
      bsp_hit3=add_imp_all_tib$bsp_hit3,
      bsp_hit4=add_imp_all_tib$bsp_hit4,
      bsp_hit5=add_imp_all_tib$bsp_hit5,
      
      # Sesame/Previous Data should be handled
      #  in downstream analysis...
      #
      # ses_add=add_imp_all_tib$ses_add,
      # ses_cgn=add_imp_all_tib$ses_cgn,
      # ses_typ=add_imp_all_tib$ses_typ,
      # ses_chr=add_imp_all_tib$ses_chr,
      # ses_pos=add_imp_all_tib$ses_pos,
      # ses_fr=add_imp_all_tib$ses_fr,
      # ses_nxb=add_imp_all_tib$ses_nxb,
      # prb_ord_seq=add_imp_all_tib$prb_ord_seq,
      # ses_gene=add_imp_all_tib$ses_gene,
      # ses_HGNC=add_imp_all_tib$ses_HGNC,
      # ses_src=add_imp_all_tib$ses_src,
      # MASK_mapping=add_imp_all_tib$MASK_mapping,
      # MASK_general=add_imp_all_tib$MASK_general,
      
      IRanges(start=add_imp_all_tib$start,
              end=add_imp_all_tib$end,
              names=add_imp_all_tib$Seq_ID) )
  
  readr::write_rds(add_imp_all_grs, add_imp_all_rds, compress="gz")
} else {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading Address Sesame/PREMAP RDS={add_imp_ses_rds}...{RET}"))
  add_imp_ses_grs <- readr::read_rds(add_imp_ses_rds)
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading Address BSMAP/IMPROBE RDS={add_imp_all_rds}...{RET}"))
  add_imp_all_grs <- readr::read_rds(add_imp_all_rds)
}

#
# Sesame Tibble Rebuild::
#
add_imp_ses2_tib <- 
  add_imp_ses_grs %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var="Probe_ID") %>% 
  tibble::as_tibble()

add_imp_ses_sum<
  add_imp_ses2_tib %>% 
  dplyr::distinct(prb_add, .keep_all=TRUE) %>%
  dplyr::group_by(prb_src,prb_des,
                  prb_srd) %>% 
  dplyr::summarise(Count=n(), .groups="drop")
if (opt$verbose>=1)
  add_imp_ses2_sum %>% print(n=base::nrow(add_imp_ses2_sum))

#
# All Data Tibble Rebuild::
#
add_imp_all2_tib <- 
  add_imp_all_grs %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var="Probe_ID") %>% 
  tibble::as_tibble()

add_imp_all2_sum <- 
  add_imp_all2_tib %>% 
  dplyr::arrange(imp_aln_scr, bsp_din_scr) %>% 
  dplyr::distinct(prb_add, .keep_all=TRUE) %>%
  dplyr::group_by(prb_src,prb_des,
                  imp_aln_scr,imp_cgn_scr,bsp_din_scr,
                  bsp_srd) %>% 
  dplyr::summarise(Count=n(), .groups="drop")
if (opt$verbose>=1)
  add_imp_all2_sum %>% print(n=base::nrow(add_imp_all2_sum))

# NOT USED...
#
# add_imp_all2_tib %>% 
#   dplyr::inner_join(add_raw_ses_tib, by=c("prb_ord"="prb_ord_seq")) %>%
#   dplyr::mutate(
#     pos_mat=dplyr::case_when(
#       start==ses_pos ~ 0,
#       TRUE ~ 1)
#   ) %>%
#   dplyr::arrange(pos_mat) %>%
#   dplyr::distinct(prb_add)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Annotation Direct:: NCBI
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Double check tango addresses against single IDAT file
#
idat_tib <- loadIdatAddress(prefix="/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R01C01")

ncbi_ann_tsv <- 
  file.path(opt$annDir,opt$genomeBuild, paste(opt$genomeBuild,'ncbi.RefSeqGenes.tsv.gz', sep='.') )

if (!is.null(ncbi_ann_tsv) && file.exists(ncbi_ann_tsv)) {
  # Need to make Unique ID::
  ncbi_fet_tib <- loadNcbiGene_COVIC(file=ncbi_ann_tsv, verbose=10)
  
  # Build gene GRS
  ncbi_ann_grs <- 
    GRanges(
      seqnames = Rle(ncbi_fet_tib$chr),
      strand=Rle(ncbi_fet_tib$srd),
      
      Gene_NCBI=ncbi_fet_tib$gene,
      Gene_NCBI_Simple=stringr::str_remove(ncbi_fet_tib$gene,"-.*$"),
      Tran_NCBI=ncbi_fet_tib$tran,
      Feat_NCBI=ncbi_fet_tib$feat,
      Gene_FR=ncbi_fet_tib$srd,
      
      IRanges(start=ncbi_fet_tib$beg,
              end=ncbi_fet_tib$end,
              names=ncbi_fet_tib$Unique_ID )
    )
  
  if (FALSE) {
    ncbi_ann_tib <- suppressMessages(suppressWarnings( readr::read_tsv(ncbi_ann_tsv) )) %>% 
      dplyr::mutate(
        strand_FR=stringr::str_replace(strand,'\\+',"F") %>% stringr::str_replace("-","R"), 
        Unique_ID=paste(name,chrom,strand_FR,txStart,txEnd, sep="_")
      ) %>%
      dplyr::select(Unique_ID,dplyr::everything()) # %>%
    
    # Build gene GRS
    ncbi_ann_grs <- 
      GRanges(
        seqnames = Rle(ncbi_ann_tib$chrom),
        strand=Rle(ncbi_ann_tib$strand),
        
        Gene_NCBI=ncbi_ann_tib$name2,
        Gene_NCBI_Simple=stringr::str_remove(ncbi_ann_tib$name2,"-.*$"),
        Tran_NCBI=ncbi_ann_tib$name,
        Gene_Bin=ncbi_ann_tib$bin,
        Gene_FR=ncbi_ann_tib$strand_FR,
        
        IRanges(start=ncbi_ann_tib$txStart,
                end=ncbi_ann_tib$txEnd,
                names=ncbi_ann_tib$Unique_ID )
      )
  }
  
  #
  # Run intersection directly:: Sesame + Previous Products
  #
  add_imp_ses_ncbi_tib <- 
    intersectGranges(man=add_imp_ses_grs, 
                     ref=ncbi_ann_grs,
                     verbose=opt$verbose,
                     vt=1,tc=1,tt=pTracker) %>% 
    dplyr::rename(Gene_Chrom=chrom,Gene_Beg=chromStart,
                  Gene_End=chromEnd, Gene_Length=chromLength,
                  Gene_Strand=chromStrand)

  #
  # Run intersection directly:: All
  #
  add_imp_all_ncbi_tib <- 
    intersectGranges(man=add_imp_all_grs, 
                     ref=ncbi_ann_grs,
                     verbose=opt$verbose,
                     vt=1,tc=1,tt=pTracker) %>% 
    dplyr::rename(Gene_Chrom=chrom,Gene_Beg=chromStart,
                  Gene_End=chromEnd, Gene_Length=chromLength,
                  Gene_Strand=chromStrand) %>% 
    dplyr::add_count(prb_add,Gene_NCBI,Feat_NCBI, name="Feat_Count") %>%
    dplyr::arrange(-Feat_Count)
  
  #
  # - Pick best gene hits for add_imp_all_ncbi_tib
  # - left join by ord_seq add_raw_ses_tib & add_imp_all_ncbi_tib$best_hits
  # - remove above from add_imp_all_ncbi_tib
  #
  
  # Remove M address::
  # TBD:: Future analysis keep this values...
  add_cor_ses_tib <- dplyr::filter(add_raw_ses_tib, ses_des != "M")
  
  imp_mat_ses_tib <- 
    add_imp_all_ncbi_tib %>% 
    dplyr::left_join(add_cor_ses_tib, 
                     by=c("prb_ord"="prb_ord_seq"), 
                     suffix=c("_New", "_Ses") ) %>%
    dplyr::mutate(
      pos_delta=dplyr::case_when(
        imp_chr==ses_chr ~ base::abs(imp_pos-ses_pos),
        TRUE ~ NA_real_
      ) %>% as.integer()
    ) %>% 
    dplyr::arrange(pos_delta) %>%
    dplyr::distinct(ses_add, .keep_all=TRUE) %>%
    dplyr::select(dplyr::all_of(names(add_imp_all_ncbi_tib)))
  
  add_ann_sesA_tib <- 
    add_cor_ses_tib %>% 
    dplyr::left_join(imp_mat_ses_tib, 
                     by=c("prb_ord_seq"="prb_ord"),
                     suffix=c("_Seq", "_New") )
  
  add_ann_ses_tib <- bind_rows(
    add_ann_sesA_tib %>% 
      dplyr::filter(!is.na(prb_src) | !is.na(prb_des)) %>% 
      dplyr::filter(ses_add==prb_add & ses_src==prb_src & ses_des==prb_des),
    add_ann_sesA_tib %>% 
      dplyr::filter(is.na(prb_src) | is.na(prb_des))
  )  %>% 
    dplyr::select(-prb_add,-prb_src,-prb_des) %>%
    dplyr::rename(prb_add=ses_add, prb_src=ses_src, prb_des=ses_des)
  
  # Summary::
  add_ann_ses_tib %>% 
    dplyr::group_by(prb_src,prb_des) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  # - remove above from add_imp_all_ncbi_tib
  # - pick best by Feat_Count
  #
  # imp_aln_scr = [0,1,2, 3] = [ 50U, 50M, 49U, NA ]
  # bsp_din_scr = [0,1,2,3, 4] = [ CG, CG, G*, *C, NA ] = [ --, +-, -+, ++, NA ] 
  # imp_cgn_scr = [0,1] = [ Match, NA ]
  #
  add_imp_new_tib <- 
    add_imp_all_ncbi_tib %>% 
    dplyr::anti_join(add_ann_ses_tib, by="prb_add") %>% 
    dplyr::filter(prb_des!="M") %>%
    # dplyr::arrange(imp_aln_scr,bsp_din_scr,-Feat_Count) %>%
    dplyr::arrange(imp_aln_scr,bsp_din_scr,imp_cgn_scr,-Feat_Count) %>%
    dplyr::distinct(prb_add, .keep_all=TRUE)
  
  # Summary::
  add_imp_new_tib %>% 
    dplyr::group_by(prb_src,prb_des) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  #
  # Join all results and check against master::
  #
  
  # Expected U/2
  add_cor_all_tib <- 
    add_raw_all_tib %>% 
    dplyr::filter(prb_des != "M") %>% 
    dplyr::distinct(prb_add, .keep_all=TRUE)
  
  # Summary::
  add_cor_all_tib %>% 
    dplyr::group_by(prb_src,prb_des) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  # Final Output:: Current... NEEDS FIXING...
  add_fin_all_tibA <- dplyr::bind_rows(
    add_ann_ses_tib,add_imp_new_tib)
  
  add_cor_mis_tib <- add_cor_all_tib %>% 
    dplyr::anti_join(add_fin_all_tibA, by="prb_add")
  
  # The missing group to be added to add_fin_all_tibA
  add_fin_mis_tib <- add_imp_all2_tib %>% 
    dplyr::filter(prb_add %in% add_cor_mis_tib$prb_add) %>% 
    dplyr::filter(prb_des != "M") %>% 
    dplyr::distinct(prb_add, .keep_all=TRUE)
  
  add_fin_all_tib <- dplyr::bind_rows(
    add_fin_all_tibA,add_fin_mis_tib) %>% 
    dplyr::filter(prb_add %in% idat_tib$Address)
  
  # Summary::
  add_fin_all_tib %>% 
    dplyr::group_by(prb_src,prb_des) %>% 
    dplyr::summarise(Count=n(), .groups="drop") %>%
    print(n=1000)
  
  #
  # Calculate Errors:: Missing probes both ways::
  #
  add_cor_all_tib %>% 
    dplyr::anti_join(add_fin_all_tib, by="prb_add") %>%
    dplyr::group_by(prb_src,prb_des) %>%
    dplyr::summarise(Count=n(), .groups="drop") %>%
    print(n=1000)
  
  add_fin_all_tib %>% 
    dplyr::anti_join(add_cor_all_tib, by="prb_add") %>%
    dplyr::group_by(prb_src,prb_des,imp_aln_scr,imp_cgn_scr,bsp_din_scr) %>%
    dplyr::summarise(Count=n(), .groups="drop") %>%
    print(n=1000)
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Recombine with Original Manifest::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # man_raw_all_tib
  # man_raw_ful_tib <- dplyr::bind_rows( man_raw_cor_tib, man_raw_cov_tib )
  
  #
  # Build seperately::
  #
  tar_cols <- c("Probe_ID","U","M","AlleleA_Probe_Sequence","AlleleB_Probe_Sequence",
                "Probe_Type","Infinium_Design","DESIGN","COLOR_CHANNEL","col","Next_Base",
                "Probe_Source","Version","Genome_Build",
                "Chromosome","Coordinate","Strand_FR","Strand_TB","Strand_CO",
                "Strand_BSMAP","BSMAP_Tag","Gene_NCBI","Tran_NCBI","Feat_NCBI",
                "Sesame_Gene","Sesame_Gene_HGNC","MASK_mapping","MASK_general")
  
  #
  # EPIC Final::
  #
  man_fin_cor_tib <- man_raw_cor_tib %>%
    dplyr::inner_join(add_fin_all_tib, 
                      by=c("U"="prb_add", "AlleleA_Probe_Sequence"="prb_ord_seq"),
                      suffix=c("","_ADD")) %>%
    dplyr::rename(Chromosome=chrom,
                  Coordinate=chromStart,
                  Strand_TB=imp_tb,
                  Strand_CO=imp_co,
                  Strand_BSMAP=bsp_srd,
                  Sesame_Gene=gene,
                  Sesame_Gene_HGNC=gene_HGNC,
                  BSMAP_Tag=bsp_tag) %>%
    dplyr::mutate(Genome_Build=opt$genomeBuild,
                  Infinium_Design=prb_des,
                  Strand_FR=strand_FR %>% 
                    stringr::str_replace("\\+","F") %>%
                    stringr::str_replace("-","R"),
                  Next_Base=nextBase,
                  Next_Base=dplyr::case_when(
                    Next_Base=="Both" ~ NA_character_,
                    TRUE ~ Next_Base
                  ),
                  DESIGN=dplyr::case_when(
                    !is.na(U) &  is.na(M) ~ "II",
                    !is.na(U) & !is.na(M) ~ "I",
                    TRUE ~ NA_character_
                  ), 
                  COLOR_CHANNEL=dplyr::case_when(
                    prb_des==2 ~ "Both",
                    Next_Base=="C" | Next_Base=="G" ~ "Grn",
                    Next_Base=="A" | Next_Base=="T" ~ "Red",
                    TRUE ~ NA_character_
                  ),
                  col=dplyr::case_when(
                    prb_des==2 ~ NA_character_,
                    Next_Base=="C" | Next_Base=="G" ~ "G",
                    Next_Base=="A" | Next_Base=="T" ~ "R",
                    TRUE ~ NA_character_
                  ),
                  Next_Base=dplyr::case_when(
                    prb_des==1 ~ stringr::str_to_upper(imp_nxb),
                    TRUE ~ NA_character_
                  ),
                  Probe_Source="EPIC",
                  Version=Manifest
    ) %>% dplyr::select(dplyr::any_of(tar_cols))
    
  
  #
  # COVIC Final::
  #
  man_fin_cov_tib <- man_raw_cov_tib %>% 
    dplyr::rename(Probe_ID_Org=Probe_ID) %>%
    dplyr::inner_join(add_fin_all_tib, by=c("U"="prb_add", "AlleleA_Probe_Sequence"="prb_ord")) %>%
    dplyr::select(!dplyr::starts_with("ses_")) %>% 
    dplyr::select(Probe_Type, U,M,
                  AlleleA_Probe_Sequence, AlleleB_Probe_Sequence, 
                  imp_cgn, imp_tb, imp_co, imp_fr, prb_des, 
                  imp_chr, imp_pos, imp_nxb, bsp_tag, bsp_srd, 
                  Gene_NCBI, Tran_NCBI, Feat_NCBI) %>% 
    dplyr::mutate(
      out_srd=paste0(imp_tb,imp_co,prb_des),
      Probe_ID=paste(imp_cgn,out_srd, sep="_"),
      Infinium_Design=prb_des,
      DESIGN=dplyr::case_when(
        !is.na(U) &  is.na(M) ~ "II",
        !is.na(U) & !is.na(M) ~ "I",
        TRUE ~ NA_character_
      ), 
      COLOR_CHANNEL=dplyr::case_when(
        prb_des==2 ~ "Both",
        imp_nxb=="C" | imp_nxb=="G" ~ "Grn",
        imp_nxb=="A" | imp_nxb=="T" ~ "Red",
        TRUE ~ NA_character_
      ),
      col=dplyr::case_when(
        prb_des==2 ~ NA_character_,
        imp_nxb=="C" | imp_nxb=="G" ~ "G",
        imp_nxb=="A" | imp_nxb=="T" ~ "R",
        TRUE ~ NA_character_
      ),
      Next_Base=dplyr::case_when(
        prb_des==1 ~ stringr::str_to_upper(imp_nxb),
        TRUE ~ NA_character_
      ),
      Probe_Source="COVIC",
      Genome_Build=opt$genomeBuild,
      Chromosome=imp_chr,
      Coordinate=imp_pos,
      Strand_FR=imp_fr,
      Strand_TB=imp_tb,
      Strand_CO=imp_co,
      Strand_BSMAP=bsp_srd,
      BSMAP_Tag=bsp_tag,
      Version="C0"
      # Version=opt$version
    ) %>% 
    dplyr::select(dplyr::any_of(tar_cols)) %>% 
    dplyr::group_by(Probe_ID) %>% 
    dplyr::mutate(Probe_ID=paste0(Probe_ID, dplyr::row_number())) %>% 
    dplyr::ungroup()
  
  # man_fin_cov_tib %>% dplyr::distinct(Probe_ID)
  
  #
  # Controls Final::
  #
  epic_ctl_tib <- readr::read_csv(file.path(par$datDir, "manifest/base/EPIC-C0.manifest.sesame-base.cpg-sorted.csv.gz")) %>% 
    dplyr::filter(Probe_Type != "cg") %>% dplyr::filter(Probe_Type != "ch") %>% dplyr::filter(Probe_Type != "rs")
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Output Manifest and Annotation Tables::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  man_fin_all_csv <- file.path(opt$preDir, paste(opt$runName,"manifest.sesame-base.cpg-sorted.csv.gz", sep="."))
  man_fin_all_tib <- dplyr::bind_rows(man_fin_cor_tib,man_fin_cov_tib,epic_ctl_tib) %>%
    dplyr::arrange(Probe_ID)
  readr::write_csv(man_fin_all_tib, man_fin_all_csv)
  
  # man_fin_all_csv <- "/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/dat/manifest/core/COVIC-C12.manifest.sesame-base.cpg-sorted.csv.gz"
  # man_fin_all_tib <- readr::read_csv(man_fin_all_csv)
  
  # man_fin_all_tib %>% dplyr::filter(! U %in% idat_tib$Address)
  # dplyr::bind_rows( 
  #   rdat$cur_list$call_dat %>% dplyr::filter(is.na(raw_betas)) %>% dplyr::inner_join(man_fin_all_tib %>% dplyr::filter(! U %in% idat_tib$Address), by="Probe_ID"),
  #   rdat$cur_list$call_dat %>% dplyr::filter(is.na(raw_betas)) %>% dplyr::inner_join(man_fin_all_tib %>% dplyr::filter(!is.na(M)) %>% dplyr::filter(! U %in% idat_tib$Address), by="Probe_ID")
  # ) %>% dplyr::distinct(Probe_ID)
  
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
