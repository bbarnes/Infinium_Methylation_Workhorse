
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
par$prgmTag <- 'build_decode_manifest'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

par$local_runType <- NULL

# Executables::
opt$Rscript <- NULL

# Run Parameters::
opt$runName   <- NULL

# Directories::
opt$outDir <- NULL
opt$impDir <- NULL
opt$annDir <- NULL

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
  
  opt$outDir <- file.path(par$topDir, 'scratch', par$prgmTag)
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda2-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda3-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/conda_4.6.8/bin/Rscript'
  
  #
  # End of local parameter definitions::
  #
  
  locIdatDir <- file.path(par$topDir, 'data/idats')
  opt$outDir <- file.path(par$topDir, 'scratch')
  opt$impDir <- file.path(par$topDir, 'data/improbe')
  opt$annDir <- file.path(par$topDir, 'data/annotation')

  # opt$cpg_s48_tsv <- file.path(opt$impDir, 'designOutput_21092020/seq48U/gz/seq48U-GRCh36-38-10-21092020.unq.noHeader.seq-sorted.tsv.gz')
  opt$cpg_s48_tsv <- file.path(opt$impDir, 'designOutput_21092020/seq48U/un/seq48U-GRCh36-38-10-21092020.unq.noHeader.seq-sorted.tsv')
  
  # The reason we don't use the general form is that specificity for the species may be lost.
  #   I think this is poor form, but using now to meet mouse GS requriments...
  #
  # opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop/GRCh36-GRCh38-GRCm10-21092020.cgnTop.sorted.tsv.gz')
  # opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop/GRCh36-GRCh38-GRCm10-21092020.cgnTop.sorted.tsv')
  
  opt$fixIds  <- TRUE
  
  opt$write_full <- FALSE
  opt$write_base <- FALSE
  
  #
  # Pre-defined local options runTypes::
  #
  par$local_runType <- 'qcMVP'
  par$local_runType <- 'NZT'
  par$local_runType <- 'COVIC'
  par$local_runType <- 'COVID'
  par$local_runType <- 'GENK'
  par$local_runType <- 'GRCm38'
  
  if (par$local_runType=='COVID') {
    opt$fresh  <- TRUE
    opt$fixIds <- FALSE
    
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
    opt$genomeBuild <- 'GRCh37'
    opt$genomeBuild <- 'GRCh38'
    opt$platform    <- 'COVIC'
    opt$version     <- 'C0'
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genomeBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genomeBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
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
    opt$fresh <- TRUE
    
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
    opt$version     <- 'C11'
    opt$version     <- 'C12'
    opt$version     <- 'C13'
    opt$version     <- 'C14'
    opt$version     <- 'C17'
    opt$version     <- 'C18'
    
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
    # opt$fresh <- TRUE
    opt$matFormat <- 'old'
    opt$ordFormat <- 'old'
    opt$matSkip <- 0
    opt$ordSkip <- 8
    opt$pqcSkip <- 7
    
    opt$genomeBuild <- 'GRCh36'
    opt$genomeBuild <- 'GRCh37'
    opt$genomeBuild <- 'GRCh38'
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
  opt$runName <- paste(opt$genomeBuild,opt$platform,opt$version, sep='-')
  
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
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    
    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    make_option(c("--impDir"), type="character", default=opt$impDir, 
                help="improbe data directory [default= %default]", metavar="character"),
    make_option(c("--annDir"), type="character", default=opt$iannDir, 
                help="Annotation data directory [default= %default]", metavar="character"),
    
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
              'genomeBuild','platform','version','Rscript','verbose')

par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )

opt <- program_init(name=opt$runName,
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

opt$probe_type <- 'cg'
opt$design_key <- 'Seq_ID'
opt$design_seq <- 'Forward_Sequence'
opt$design_seq <- 'Top_Sequence'
opt$design_prb <- 'Probe_Type'
opt$design_srs <- 'TB'
opt$design_cos <- 'CO'

# design_prb_sym <- rlang::sym(opt$design_prb)

# Input Definitions::
if (is.null(opt$ctls)) {
  opt$ctls <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
}

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
opt$manDir <- file.path(opt$outDir, 'manifest')
if (!dir.exists(opt$manDir)) dir.create(opt$manDir, recursive=TRUE)

opt$prdDir <- file.path(opt$outDir, 'product')
if (!dir.exists(opt$prdDir)) dir.create(opt$prdDir, recursive=TRUE)

opt$genDir <- file.path(opt$outDir, 'genomic')
if (!dir.exists(opt$genDir)) dir.create(opt$genDir, recursive=TRUE)

opt$desDir <- file.path(opt$outDir, 'design')
if (!dir.exists(opt$desDir)) dir.create(opt$desDir, recursive=TRUE)

opt$intDir <- file.path(opt$outDir, 'intersection')
if (!dir.exists(opt$intDir)) dir.create(opt$intDir, recursive=TRUE)

opt$sumDir <- file.path(opt$outDir, 'summary')
if (!dir.exists(opt$sumDir)) dir.create(opt$sumDir, recursive=TRUE)

cat(glue::glue("[{par$prgmTag}]: Done. Building Output Directories.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Define Manifest Output Files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
# man_out_key <- paste(opt$genomeBuild,opt$platform,opt$version, sep='-')
man_out_key <- paste(opt$platform,opt$version, sep='-')
gs_swap_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.GenomeStudio.cpg-sorted.csv', sep='.') )
gz_swap_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.GenomeStudio.cpg-sorted.csv.gz', sep='.') )

ann_base_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.annotation-base.pos-sorted.csv.gz', sep='.') )
gen_base_csv <- file.path(opt$genDir, paste(man_out_key,'manifest.genomic-base.pos-sorted.csv.gz', sep='.') )

ses_base_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.sesame-base.cpg-sorted.csv.gz', sep='.') )
pos_base_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.sesame-base.pos-sorted.csv.gz', sep='.') )
ses_mach_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.sesame-mach.cpg-sorted.csv.gz', sep='.') )
ses_epic_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.sesame-epic.cpg-sorted.csv.gz', sep='.') )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Genome Studio Color Codes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Add this to program_init()
#
color_tib <- suppressMessages(suppressWarnings( readr::read_csv(file.path(par$datDir, 'params/GenomeStudioColors.csv')) ))
color_vec <- color_tib %>% dplyr::pull(Color) %>% as.vector()
color_len <- length(color_vec)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         1.0 Process AQP/PQC Data::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$verbose <- 40
opt$verbose <- 3

opt$fixIds <- FALSE
# opt$fresh  <- TRUE
# aqps_vec   <- NULL

man_pqc_dat_tib <- 
  decodeAqpPqcWrapper(
    ord_vec=ords_vec,mat_vec=mats_vec,aqp_vec=aqps_vec,pqc_vec=pqcs_vec,
    platform=opt$platform,version=opt$version,
    ordFormat=opt$ordFormat,ordSkip=opt$ordSkip,
    matFormat=opt$matFormat,matSkip=opt$matSkip,
    aqpSkip=opt$aqpSkip,
    pqcSkip=opt$pqcSkip,
    name=opt$runName,outDir=opt$manDir,origin=opt$time_org_txt,
    fresh=opt$fresh,fixIds=opt$fixIds,full=par$retData,trim=TRUE,
    verbose=opt$verbose,vt=1,tc=0,tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                2.0 Build Probes foreach Type:: RS/CH/CG/etc..
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_prb_list <- NULL
man_raw_list <- man_pqc_dat_tib %>% split(.$Probe_Type)

for (probe_type in rev(names(man_raw_list)) ) {
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Building Probes Types={probe_type}...{RET}"))

  if (probe_type == 'rs') {
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           2.1 Build Probes:: SNP
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

    man_prb_list[[probe_type]] <- NULL
    if (!is.null(opt$snp_des_csv) && file.exists(opt$snp_des_csv) ) {
      man_prb_list[[probe_type]] = 
        adhoc_desToMAN(man=man_pqc_dat_tib, des_csv=opt$snp_des_csv, probe_type=probe_type,
                       name=opt$runName,outDir=opt$desDir,origin=opt$time_org_txt,
                       idsKey=opt$design_key,
                       seqKey=opt$design_seq,prbKey=opt$design_prb,
                       strsSR=opt$design_srs,strsCO=opt$design_cos,
                       parallel=opt$parallel,fresh=opt$fresh,fixIds=FALSE,
                       verbose=opt$verbose,vt=3,tc=1,tt=pTracker)
    }
  } else if (probe_type == 'ch') {
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           2.2 Build Probes:: CpH
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

    man_prb_list[[probe_type]] <- NULL
    if (!is.null(opt$cph_des_csv) && file.exists(opt$cph_des_csv) ) {
      man_prb_list[[probe_type]] =
        adhoc_desToMAN(man=man_pqc_dat_tib, des_csv=opt$cph_des_csv, probe_type=probe_type,
                       name=opt$runName,outDir=opt$desDir,origin=opt$time_org_txt,
                       idsKey=opt$design_key,
                       seqKey=opt$design_seq,prbKey=opt$design_prb,
                       strsSR=opt$design_srs,strsCO=opt$design_cos,
                       parallel=opt$parallel,fresh=opt$fresh,fixIds=FALSE,
                       verbose=opt$verbose,vt=3,tc=1,tt=pTracker)
    }
  } else if (probe_type == 'cg') {
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           2.3 Build Probes:: CpG
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    man_prb_list[[probe_type]] <- NULL
    if (!is.null(opt$cpg_s48_tsv) && file.exists(opt$cpg_s48_tsv) &&
        !is.null(opt$cpg_top_tsv) && file.exists(opt$cpg_top_tsv)) {
      man_prb_list[[probe_type]] =
        clean_manifest_probes(
          tib=man_pqc_dat_tib,
          s48_tsv=opt$cpg_s48_tsv,top_tsv=opt$cpg_top_tsv,
          name=opt$runName,outDir=opt$outDir,origin=opt$time_org_txt,
          
          design_key=opt$design_key,design_seq=opt$design_seq,
          design_prb=opt$design_prb,probe_type=probe_type,
          design_srs=opt$design_srs,design_cos=opt$design_cos,
          parallel=opt$parallel,fresh=opt$fresh,
          
          verbose=opt$verbose,vt=3,tc=0,tt=pTracker)
    }
  } else {
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Probes Types={probe_type} NOT supported!{RET}{RET}"))
  }
}
if (opt$verbose>=1) {
  cat(glue::glue("[{par$prgmTag}]: Done. Building Probes.{RET}{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 3.0 Join All Detected Probes:: SNP/CpH/CpG
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: bindProbeDesignList should use inferFunction...
#
man_pqc_prb_tib <- bindProbeDesignList(
  list=man_prb_list, platform=opt$platform, version=opt$version,
  sumDir=opt$sumDir,del='.',
  verbose=opt$verbose,vt=3,tc=1,tt=pTracker) %>% 
  dplyr::mutate(Assay_Class='Analytical')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#            Extract Missing Probes from Manifest:: i.e. CTL/UNK
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_pqc_unk_tib <- dplyr::anti_join(
  dplyr::distinct(man_pqc_dat_tib, M,U, .keep_all=TRUE),
  dplyr::distinct(man_pqc_prb_tib, M,U, .keep_all=TRUE),
  by=c("M","U")) %>% dplyr::distinct(M,U, .keep_all=TRUE) %>%
  dplyr::rename(Strand_TB=TB, Strand_CO=CO) %>% 
  dplyr::mutate(
    Strand_TB=as.character(Strand_TB),
    Strand_CO=as.character(Strand_CO),
    Strand_TB=dplyr::case_when(
      is.na(Strand_TB) ~ 'N',
      TRUE ~ Strand_TB),
    Strand_CO=dplyr::case_when(
      is.na(Strand_CO) ~ 'N',
      TRUE ~ Strand_CO)
  ) %>%
  dplyr::add_count(Seq_ID,Strand_TB,Strand_CO,Infinium_Design, name='Rep_Max') %>%
  dplyr::group_by(Seq_ID,Strand_TB,Strand_CO,Infinium_Design) %>%
  dplyr::mutate(
    Rep_Cnt=dplyr::row_number(),
    Probe_ID=paste0(Seq_ID,'_',Strand_TB,Strand_CO,Infinium_Design,Rep_Cnt),
    Probe_Class=dplyr::case_when(
      Probe_Type=='cg' ~ 'UK',
      Probe_Type=='ch' ~ 'UK',
      Probe_Type=='mu' ~ 'UK',
      Probe_Type=='rp' ~ 'UK',
      Probe_Type=='rs' ~ 'UK',
      Probe_Type=='BS' ~ 'BS',
      Probe_Type=='NG' ~ 'NG',
      Probe_Type=='NP' ~ 'NP',
      TRUE ~ 'UK'),
    Assay_Class=dplyr::case_when(
      Probe_Type != 'cg' & Probe_Type != 'ch' & Probe_Type != 'mu' & 
        Probe_Type != 'rp' & Probe_Type != 'rs' ~ 'Control',
      TRUE ~ 'Testing'),
    Probe_ID=dplyr::case_when(
      Assay_Class=='Control' ~ paste('ctl',Probe_ID, sep='-'),
      TRUE ~ Probe_ID),
    FN=NA_character_,
    Probe_Source=opt$platform,
    Version=opt$version,
    M=as.integer(M),
    U=as.integer(U),
    Strand_TB=as.character(Strand_TB),
    Strand_CO=as.character(Strand_CO),
    Infinium_Design=as.integer(Infinium_Design),
    Rep_Num=as.integer(Rep_Num),
    Old_Probe_ID=as.character(Old_Probe_ID),
    Di=as.character(Di),
    DS=as.integer(DS),
    HS=as.integer(HS),
    FN=as.character(FN),
    ValidID=FALSE,
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Probe_ID,M,U,DESIGN,COLOR_CHANNEL,col,Next_Base,
                Seq_ID,Probe_Type,Probe_Source,Version,
                Strand_TB,Strand_CO,Infinium_Design,Rep_Num,
                AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                everything()) %>%
  dplyr::arrange(Probe_ID)


# Full Unite::
#  TBD:: Should remove Rep_Num from source...
#
man_pqc_all_col <- names(man_pqc_prb_tib)[names(man_pqc_prb_tib) %in% names(man_pqc_unk_tib)]
man_pqc_prb_col <- c(man_pqc_all_col, 'Top_Sequence')
man_pqc_unk_col <- c(man_pqc_all_col)

# dplyr::select(man_pqc_prb_tib, dplyr::all_of(man_pqc_all_col))
# dplyr::select(man_pqc_unk_tib, dplyr::all_of(man_pqc_all_col))

man_pqc_all_tib <- dplyr::bind_rows(
  dplyr::select(man_pqc_prb_tib, dplyr::all_of(man_pqc_prb_col)),
  dplyr::select(man_pqc_unk_tib, dplyr::all_of(man_pqc_unk_col))
) %>% dplyr::select(-Rep_Num) # %>% dplyr::mutate(Assay_Class='Analytical')

#
# Validation Stats::
#
man_pqc_all_cnt <- man_pqc_all_tib %>% base::nrow()
man_pqc_dat_cnt <- man_pqc_dat_tib %>% base::nrow()
man_dif_all_cnt <- man_pqc_all_cnt - man_pqc_dat_cnt

man_pqc_prb_cnt <- man_pqc_prb_tib %>% base::nrow()
man_pqc_unk_cnt <- man_pqc_unk_tib %>% base::nrow()
man_dif_unk_cnt <- man_pqc_dat_cnt - man_pqc_prb_cnt - man_pqc_unk_cnt

cat(glue::glue("[{par$prgmTag}]: Done; Parsing PQC/All Probes={man_pqc_all_cnt}.{RET}"))
cat(glue::glue("[{par$prgmTag}]: Done; Parsing PQC/Dat Probes={man_pqc_dat_cnt}.{RET}"))
cat(glue::glue("[{par$prgmTag}]: Done; Parsing DiffCnt Probes={man_dif_all_cnt}.{RET}{RET}"))

cat(glue::glue("[{par$prgmTag}]: Done; Parsing PQC/Prb Probes={man_pqc_prb_cnt}.{RET}"))
cat(glue::glue("[{par$prgmTag}]: Done; Parsing Unk/Ctl Probes={man_pqc_unk_cnt}.{RET}"))
cat(glue::glue("[{par$prgmTag}]: Done; Parsing DiffCnt Probes={man_dif_unk_cnt}.{RET}{RET}"))

man_unk_sum_tib <- man_pqc_unk_tib %>% 
  dplyr::group_by(Assay_Class,Probe_Class,Probe_Type,Infinium_Design,AQP) %>% 
  dplyr::summarise(Count=n(), .groups='drop')
print(man_unk_sum_tib,n=base::nrow(man_unk_sum_tib))

man_all_sum_tib <- man_pqc_all_tib %>% 
  dplyr::group_by(Assay_Class,Probe_Class,Probe_Type,Infinium_Design,AQP) %>% 
  dplyr::summarise(Count=n(), .groups='drop')
print(man_all_sum_tib,n=base::nrow(man_all_sum_tib))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#            4.0 Load all Coordinates for Detected Manifest Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: make function()
#

#
# TBD: make pos_tsv into vector input...
#

man_new_list <- NULL
man_pos_list <- NULL
man_new_list <- man_pqc_all_tib %>% split(.$Probe_Class)

for (probe_class in names(man_new_list)) {
  if (opt$verbose>=3)
    cat(glue::glue("[{par$prgmTag}]: Loading Coordinates; probe_class={probe_class}...{RET}"))
  
  pos_key <- paste(opt$runName,probe_class, sep='.')
  pos_col <- 
    cols(Seq_ID=col_character(),Gen_Chr=col_character(),Gen_Pos=col_double() )

  pos_tsv <- NULL
  if (probe_class=='cg') pos_tsv <- opt$cpg_pos_tsv
  if (probe_class=='ch') pos_tsv <- opt$cph_pos_tsv
  if (probe_class=='rs') pos_tsv <- opt$snp_pos_tsv
  if (is.null(pos_tsv)) {
    cat(glue::glue("{RET}[{par$prgmTag}]: Skipping: Unsupported probe_class={probe_class}!{RET}{RET}"))
    next
  }

  man_pos_list[[probe_class]] <- NULL
  man_pos_list[[probe_class]] =
    loadAllGenomicByMAN(
      man=man_new_list[[probe_class]], pos_tsv=pos_tsv, 
      pos_col=pos_col,name=pos_key, outDir=opt$genDir,
      verbose=opt$verbose,vt=3,tc=1,tt=pTracker)
  
  if (opt$verbose>=3)
    man_pos_list[[probe_class]] %>% 
    dplyr::group_by(Genomic_CGN_Count) %>% 
    dplyr::summarise(Genomic_Count_Hist=n(), .groups='drop') %>% print()
  
  if (opt$verbose>=3)
    cat(glue::glue("[{par$prgmTag}]: Done. Loading Coordinates.{RET}{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  4.1 Rejoin all Manifest Mapped Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: add to make function() above::
#
man_pos_tib <- man_pos_list %>% 
  dplyr::bind_rows(.id='Probe_Class') %>% 
  dplyr::rename(chrom=Gen_Chr,chromStart=Gen_Pos,score=Genomic_CGN_Count) %>%
  dplyr::mutate(
    chrom=stringr::str_remove(chrom,'^chr'),
    chrom=paste0('chr',chrom),
    chromEnd=chromStart+1) %>%
  dplyr::group_by(Seq_ID) %>% 
  dplyr::arrange(chrom,chromStart) %>%
  dplyr::select(chrom,chromStart,chromEnd,Seq_ID)

readr::write_csv(man_pos_tib,pos_base_csv)


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                5.0 Build Annotation GRanges:: Genes/Islands::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_ana_tib <- NULL
if (!is.null(opt$annDir) && dir.exists(opt$annDir)) {
  man_ana_tib <- 
    manifestToAnnotation(tib=man_pos_tib, ann=opt$annDir, gen=opt$genomeBuild,
                         verbose=opt$verbose,vt=1,tc=1,tt=pTracker)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              6.0. Collect Remainder:: CTL::Complete HSA
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

std_ctl_tib <- NULL
ctl_csv <- ctls_vec[1]
std_ctl_seq_tsv <- file.path(par$datDir,'manifest/controls/01152015_DarkMatterControls.probe.match.tsv.gz')

hsa_ctl_tib <- NULL
if (!is.null(ctl_csv) && file.exists(ctl_csv) &&
    !is.null(std_ctl_seq_tsv) && file.exists(std_ctl_seq_tsv)) {
  
  hsa_ctl_tib <- format_controls_HSA(file1=std_ctl_seq_tsv, file2=ctl_csv,
                                     verbose=opt$verbose,vt=0,tc=1,tt=pTracker)
  
  hsa_ctl_sum_tib <- hsa_ctl_tib %>% 
    dplyr::group_by(Assay_Class,Probe_Class,Probe_Type,Infinium_Design,AQP) %>% 
    dplyr::summarise(Count=n(), .groups='drop')
  print(hsa_ctl_sum_tib,n=base::nrow(hsa_ctl_sum_tib))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   6.1. Collect Remainder:: CTL::MUS
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Probably don't need this; just use man_pqc_unk_tib
#

# new_ctl_tib <- NULL
# if (par$local_runType=='GRCm38') {
#   new_ctl_tib <- format_controls_New(
#     tib=man_pqc_dat_tib, verbose=opt$verbose,vt=0,tc=1,tt=pTracker)
# }
# 
# new_ctl_tib2 <- NULL
# if (par$local_runType=='GRCm38') {
#   new_ctl_tib2 <- format_controls_Unk(
#     tib=man_pqc_unk_tib, verbose=opt$verbose,vt=0,tc=1,tt=pTracker)
# }

man_prd_all_tib <- dplyr::bind_rows(man_pqc_all_tib,hsa_ctl_tib)
man_prd_sum_tib <- man_prd_all_tib %>% 
  dplyr::group_by(Assay_Class,Probe_Class,Probe_Type,Infinium_Design,AQP) %>% 
  dplyr::summarise(Count=n(), .groups='drop')
print(man_prd_sum_tib,n=base::nrow(man_prd_sum_tib))


#
# QC Checks::
#
# man_prd_all_tib %>% dplyr::group_by(DESIGN,COLOR_CHANNEL,col,Next_Base) %>% dplyr::summarise(Count=n()) %>% as.data.frame()

man_prd_all_tib %>% dplyr::select(Probe_Type, 1:10) %>% split(.$Probe_Type)
man_prd_all_tib %>% dplyr::select(Probe_Type,11:20) %>% split(.$Probe_Type)
man_prd_all_tib %>% dplyr::select(Probe_Type,21:30) %>% split(.$Probe_Type)
man_prd_all_tib %>% dplyr::select(Probe_Type,31:40) %>% split(.$Probe_Type)
man_prd_all_tib %>% dplyr::select(Probe_Type,41:50) %>% split(.$Probe_Type)


#
# Write Sesame Output Full::
#
readr::write_csv(man_prd_all_tib,ses_base_csv)












# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              2.5. Collect Remainder:: CTL::Selected HSA/MUS
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fin_ctl_mat_tib <- NULL
if (par$local_runType=='GRCm38') {
  #
  # Load Manually Selected Human and Mouse Controls for Genome Studio::
  #
  sel_ctl_col <- c('Address','Probe_Type','Control_Color','Probe_ID')
  sel_ctl_csv <- file.path(par$topDir,'data/manifests/tmp/mm10_LEGX_B13.manifest.GenomeStudio.cpg-sorted.clean.controls.csv.csv')
  sel_ctl_tib <- suppressMessages(suppressWarnings( readr::read_csv(sel_ctl_csv, col_names=sel_ctl_col) )) %>%
    dplyr::distinct(Address, .keep_all=TRUE)
  
  ses_sel_ctl_tib <- dplyr::bind_rows(
    new_ctl_tib %>% dplyr::filter(U %in% sel_ctl_tib$Address) %>% 
      dplyr::select(Seq_ID,Infinium_Design,U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,SR_Str,CO_Str,
                    dplyr::everything()),
    
    std_bs1_hum_tib %>% dplyr::filter(U %in% sel_ctl_tib$Address) %>% 
      dplyr::mutate(SR_Str=NA,CO_Str=NA) %>%
      dplyr::select(Seq_ID,Infinium_Design,U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,SR_Str,CO_Str,
                    dplyr::everything()),
    
    std_inf1_hum_tib %>% dplyr::filter(U %in% sel_ctl_tib$Address) %>% 
      dplyr::mutate(SR_Str=NA,CO_Str=NA) %>%
      dplyr::select(Seq_ID,Infinium_Design,U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,SR_Str,CO_Str,
                    dplyr::everything())
  ) %>% 
    dplyr::mutate(Seq_ID=stringr::str_replace_all(Seq_ID,'_','-'),
                  # Seq_ID=paste('ctl',Seq_ID, sep='-'),
                  Seq_ID=paste(Seq_ID, sep='-'),
                  Seq_ID=paste(Seq_ID,Probe_Source, sep='_'))
  
  ses_sel_ctl_cnt <- ses_sel_ctl_tib %>% base::nrow()
  ses_unq_ctl_cnt <- ses_sel_ctl_tib %>% dplyr::distinct(U,M) %>% base::nrow()
  ses_sel_ctl_sum <- ses_sel_ctl_tib %>% dplyr::group_by(Probe_Source,Probe_Type,Infinium_Design) %>% 
    dplyr::summarise(Count=n(), .groups='drop')
  cat(glue::glue("[{par$prgmTag}]: ses_sel_ctl_cnt={ses_sel_ctl_cnt}, ses_unq_ctl_cnt={ses_unq_ctl_cnt}; ses_sel_ctl_sum={RET}"))
  ses_sel_ctl_sum %>% print()
  
  fin_ctl_mat_tib <- ses_sel_ctl_tib %>% 
    dplyr::select(Seq_ID, SR_Str,CO_Str,Infinium_Design, Rep_Num, Probe_Type, U,M,
                  AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, NXB_D, Top_Sequence,
                  dplyr::everything() ) %>% 
    dplyr::arrange(Seq_ID) %>%
    dplyr::mutate(Assay_Class='Control', Seq_ID=paste('ctl',Seq_ID, sep='-'))
  
} else {
  fin_ctl_mat_tib <-
    dplyr::bind_rows(std_inf1_hum_tib,std_bs1_hum_tib) %>%
    dplyr::distinct(M,U, .keep_all=TRUE) %>%
    dplyr::mutate(Seq_ID=Probe_ID,
                  SR_Str=NA,
                  CO_Str=NA,
                  NXB_D=NA,
                  Top_Sequence=NA,
                  Assay_Class='Control',
                  IlmnID=Seq_ID) %>%
    dplyr::select(IlmnID, SR_Str,CO_Str,Infinium_Design, Rep_Num, Probe_Type, U,M,
                  AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, NXB_D, Top_Sequence,
                  dplyr::everything() ) %>% 
    dplyr::arrange(Seq_ID) # %>%
    # dplyr::mutate(Assay_Class='Control', Seq_ID=paste('ctl',Seq_ID, sep='-'))
  
  #
  # Old Method:: Need to re-evaluate::
  #
  if (FALSE) {
    fin_ctl_mat_tib <- dplyr::bind_rows(std_ctl_tib,new_ctl_tib) %>%
      dplyr::distinct(M,U, .keep_all=TRUE) %>%
      dplyr::mutate(Seq_ID=Probe_ID,
                    SR_Str=NA,
                    CO_Str=NA) %>%
      dplyr::select(Seq_ID, SR_Str,CO_Str,Infinium_Design, Rep_Num, Probe_Type, U,M,
                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, NXB_D, Top_Sequence,
                    dplyr::everything() ) %>% 
      dplyr::arrange(Seq_ID) %>%
      dplyr::mutate(Assay_Class='Control', Seq_ID=paste('ctl',Seq_ID, sep='-'))
  }
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               6.1 Intersect improbe Probe Thermo Stats::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  mat_key <- 'Seq_ID'
  man_tsv <- file.path(opt$intDir, paste(opt$runName,'cgn_to_impScr','man',paste0(mat_key,'-sorted'),'tsv', sep='.'))
  int_tsv <- file.path(opt$intDir, paste(opt$runName,'cgn_to_impScr','int',paste0(mat_key,'-sorted'),'tsv.gz', sep='.'))
  
  imp_dat_tsv <- file.path(opt$impDir, 'designOutput_21092020/stats/GRCm10-21092020_improbe-designOutput.cgn-sorted.tsv.gz')
  imp_dat_tsv <- file.path(opt$impDir, 'designOutput_21092020/stats/GRCm10-21092020_improbe-designOutput.cgn-sorted.tsv')
  
  imp_col_tsv <- file.path(opt$impDir, 'designOutput_21092020/stats/GRCm10-21092020_improbe-designOutput.header.tsv')
  imp_col_vec <- suppressMessages(suppressWarnings( readr::read_tsv(imp_col_tsv) )) %>% 
    names() %>% as.vector()
  
  test_max <- 0
  man_scr_tib <- intersect_tsv(
    man=man_pqc_all_tib, 
    man_tsv=man_tsv, imp_tsv=imp_dat_tsv, int_tsv=int_tsv, 
    colA=1,colB=1, mat_vec=c(mat_key), int_col=imp_col_vec, 
    fresh=opt$fresh, max=test_max,
    verbose=opt$verbose+2,tc=1,tt=pTracker)
  
  # man_scr_tib <- readr::read_tsv(int_tsv)
}



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Write Sesame Manifest:: Standard
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                  3.0. Collect Remainder:: CTL/SNP/CpH/CpG::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  fin_ses_tib <- 
    dplyr::bind_rows(
      man_prb_tib %>% dplyr::mutate(Assay_Class='Analytical'),
      fin_ctl_mat_tib ) %>%
    dplyr::distinct(M,U, .keep_all=TRUE) %>%
    dplyr::arrange(IlmnID)
  
  fin_ses_tib$Assay_Class %>% unique()
  
  #
  # Super Quick Fix::
  #
  ctl_only_csv <- '/Users/bretbarnes/Documents/data/CustomContent/Genknowme/manifest/GKME_Sesame_EPIC-Controls-Only.csv'
  ctl_only_tib <- readr::read_csv(ctl_only_csv, col_names=c('Probe_ID','M','U','Infinium_Design'))
  
  out_ses_tib <- dplyr::bind_rows(
    fin_ses_tib %>% dplyr::filter(!stringr::str_starts(IlmnID,'ctl')) %>%
      dplyr::mutate(Infinium_Design=DESIGN),
    ctl_only_tib) %>%
    dplyr::arrange(Probe_ID) %>% 
    dplyr::distinct(M,U, .keep_all=TRUE) %>%
    dplyr::select(Probe_ID:Next_Base,Seq_ID,Probe_Type,Strand_TB,Strand_CO,
                  Infinium_Design,Rep_Num,AlleleA_Probe_Sequence,
                  AlleleB_Probe_Sequence,Top_Sequence)
  
  readr::write_csv(out_ses_tib,ses_base_csv)


  #
  # Not sure the difference between the above and below::
  #
  
  out_ses_tib <- fin_ses_tib %>% 
    dplyr::mutate(Probe_ID=IlmnID) %>%
    dplyr::arrange(Probe_ID) %>% 
    dplyr::distinct(M,U, .keep_all=TRUE) %>%
    dplyr::select(Probe_ID:Next_Base,Seq_ID,Probe_Type,Strand_TB,Strand_CO,
                  Infinium_Design,Rep_Num,AlleleA_Probe_Sequence,
                  AlleleB_Probe_Sequence,Top_Sequence)
  
  readr::write_csv(out_ses_tib,ses_base_csv)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Write Sesame Manifest:: Machine Learning
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Make all mu and rp probes Infinium I:: Future Project
#
if (FALSE) {
  fin_ses_mach_tib <- fin_core_ann_tib %>% 
    dplyr::mutate(Probe_ID=IlmnID) %>%
    dplyr::arrange(Probe_ID) %>% 
    dplyr::distinct(M,U, .keep_all=TRUE) %>% 
    dplyr::select(Probe_ID:Top_Sequence,Assay_Class,MFG_Change_Flagged)
  readr::write_csv(fin_ses_core_tib,ses_mach_csv)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Write Genome Studio Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  fin_core_split <- fin_core_ann_tib %>% split(.$Assay_Class)
  
  #
  # Format Genome Studio Analytical 
  #
  man_gs_col <- c("IlmnID","Name","AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq",
                  "Infinium_Design_Type","Next_Base","Color_Channel","Forward_Sequence","Top_Sequence",
                  "Genome_Build","Genome_Build_UCSC","CHR","MAPINFO","SourceSeq",
                  "Probe_Type","Strand","Strand_TB","Strand_CO","MFG_Change_Flagged",
                  
                  "N_Shelf","N_Shore","CpG_Island","CpG_Island_chrom","CpG_Island_chromStart",
                  "CpG_Island_chromEnd","CpG_Island_length","CpG_Island_cpgNum","CpG_Island_gcNum","CpG_Island_perCpg",
                  "CpG_Island_perGc","CpG_Island_obsExp","S_Shore","S_Shelf")
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Define Genome Studio Headers::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (FALSE) {
    man_gs_csv  <- file.path(par$topDir, 'data/manifests/tmp/mm10_LEGX_B13.manifest.GenomeStudio.cpg-sorted.clean.csv.gz')
    man_head_df <- readr::read_lines(man_gs_csv, n_max = 6) %>% as.data.frame()
    man_head_replace_str <- paste0('_',opt$version,'.csv.gz')
    man_head_df <- man_head_df %>% as_tibble() %>% purrr::set_names(c("Head_Var")) %>% 
      dplyr::mutate(Head_Var=as.character(Head_Var),
                    Head_Var=stringr::str_replace(Head_Var, '_B[0-9]+.csv.gz$', man_head_replace_str) ) %>% 
      dplyr::add_row(Head_Var="[Assay],,,,,,,,,,,,,,,,") %>%
      as.data.frame()
  }
  
  comma_str  <- ',,,,,,,,,,,,,,,,'
  loci_count <- fin_core_split[['Analytical']] %>% base::nrow()
  gs_head_colA <- c('Illumina, Inc.','[Heading]','Descriptor File Name','Assay Format','Date Manufactured','Loci Count','Assay')
  gs_head_colB <- c('','',base::basename(gs_swap_csv),'Infinium 2',as.character( base::date() ),loci_count,'')
  gs_head_colC <- c(comma_str,comma_str,comma_str,comma_str,comma_str,comma_str,comma_str)
  man_head_tib <- tibble::tibble( ValA=gs_head_colA, ValB=gs_head_colB, ValC=gs_head_colC )
  man_head_df  <- man_head_tib %>% as.data.frame()
  ctl_head_df  <- tibble::tibble(Head="[Controls]") %>% tibble::add_column(BL1='',BL2='',BL3='',BL4='') %>% as.data.frame()
  
  gb_ucsc <- NA_character_
  if (opt$genomeBuild=='GRCm38') gb_ucsc <- 'mm10'
  if (opt$genomeBuild=='GRCh38') gb_ucsc <- 'hg38'
  if (opt$genomeBuild=='GRCh37') gb_ucsc <- 'hg19'
  if (opt$genomeBuild=='GRCh36') gb_ucsc <- 'hg18'
  
  man_fin_tib <- fin_core_split[['Analytical']] %>% 
    dplyr::arrange(IlmnID) %>%
    dplyr::distinct(M,U, .keep_all=TRUE) %>% 
    dplyr::mutate(
      Probe_ID=IlmnID,
      SourceSeq=PRB1_D,
      Forward_Sequence=dplyr::case_when(
        FR=='F' & Strand_TB=="T" ~ Top_Sequence,
        FR=='F' & Strand_TB=="B" ~ revCmp(Top_Sequence),
        FR=='R' & Strand_TB=="B" ~ Top_Sequence,
        FR=='R' & Strand_TB=="T" ~ revCmp(Top_Sequence),
        TRUE ~ paste0(paste0(rep('N',60), collapse=''),'[NN]',paste0(rep('N',60), collapse=''))
      ),
      Genome_Build=opt$genomeBuild,
      Genome_Build_UCSC=gb_ucsc
    ) %>%
    # dplyr::rename(IlmnID=Probe_ID, Name=Seq_ID, 
    dplyr::rename(Name=Seq_ID, 
                  AddressA_ID=U,AlleleA_ProbeSeq=AlleleA_Probe_Sequence, 
                  AddressB_ID=M,AlleleB_ProbeSeq=AlleleB_Probe_Sequence,
                  Infinium_Design_Type=DESIGN,
                  Color_Channel=COLOR_CHANNEL,
                  CHR=Gen_Chr,
                  MAPINFO=Gen_Pos,
                  Strand=FR) %>%
    dplyr::select(all_of(man_gs_col))
  
  #
  # Format Genome Studio Controls
  #
  ctl_fin_tib <- fin_core_split[['Control']] %>% 
    dplyr::arrange(Probe_ID) %>%
    dplyr::distinct(M,U, .keep_all=TRUE) %>% 
    dplyr::select(M,U,Control_Group,Probe_ID:Top_Sequence,Control_Group)
  
  #
  # TBD:: For Infinium I designs add appropriate U/C designation to the Probe_ID
  #
  
  #
  # Find this 24637490 human bs1 control is it U?
  #
  ctl_gs_tib <- dplyr::bind_rows(
    ctl_fin_tib %>%
      dplyr::select(U,Control_Group,Probe_ID) %>%
      dplyr::filter(!is.na(U)) %>% 
      dplyr::rename(Address=U) %>%
      dplyr::mutate(Probe_ID=dplyr::case_when(
        Control_Group=='BISULFITE CONVERSION I' ~ stringr::str_replace(Probe_ID, '_([MUSHA]+)$', 'U_\\$1') %>% stringr::str_remove_all('\\\\'),
        TRUE ~ Probe_ID )
      ),
    ctl_fin_tib %>%
      dplyr::select(M,Control_Group,Probe_ID) %>%
      dplyr::filter(!is.na(M)) %>% 
      dplyr::rename(Address=M) %>%
      dplyr::mutate(Probe_ID=dplyr::case_when(
        Control_Group=='BISULFITE CONVERSION I' ~ stringr::str_replace(Probe_ID, '_([MUSHA]+)$', 'M_\\$1') %>% stringr::str_remove_all('\\\\'),
        TRUE ~ Probe_ID )
      )
  ) %>%
    dplyr::arrange(Control_Group) %>% 
    dplyr::inner_join(sel_ctl_tib %>% dplyr::select(Address,Control_Color), by="Address") %>%
    dplyr::select(Address,Control_Group,Control_Color,Probe_ID)
  
  # Sanity Check::
  #  ctl_gs_tib %>% dplyr::filter(Control_Group=='BISULFITE CONVERSION I')
  
  #
  # Write Genome Studio CSV::
  #
  
  write_delim(x=man_head_df, path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=FALSE)
  write_delim(x=man_fin_tib, path=gs_swap_csv, delim=',', col_names=TRUE,  quote_escape=FALSE, na='', append=TRUE)
  write_delim(x=ctl_head_df, path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)
  write_delim(x=ctl_gs_tib,  path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)
  
  cmd <- glue::glue("cat {gs_swap_csv} | perl -pe 's/\"//gi; s/\n/,,,,,\n/;' | gzip -c - > {gz_swap_csv}")
  system(cmd)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
