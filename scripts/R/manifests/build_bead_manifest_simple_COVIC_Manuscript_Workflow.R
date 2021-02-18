
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
par$prgmTag <- 'build_bead_manifest_simple_COVIC_Manuscript'
# par$prgmTag <- 'build_bead_manifest_simple_COVIC_v3'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

par$manDir        <- NULL
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
    opt$version     <- 'C3'
    opt$version     <- 'C4'
    opt$version     <- 'C5'
    opt$version     <- 'C6'
    
    par$manDir <- "/Users/bretbarnes/Documents/data/manifests"
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
    opt$version     <- 'C11'
    opt$version     <- 'C12'
    opt$version     <- 'C13'
    opt$version     <- 'C14'
    opt$version     <- 'C17'
    opt$version     <- 'C18'
    opt$version     <- 'C19'
    opt$version     <- 'C20'
    opt$version     <- 'C21'
    opt$version     <- 'C22'
    opt$version     <- 'C23'
    opt$version     <- 'C24'
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Workflow::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Previous Place to Start::
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                 STEP 1::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Write all probes to fasta::
fas_raw_all_csv <- "/Users/bretbarnes/Documents/data/COVIC/fas/covic-prbU-all.csv.gz"
fas_raw_all_fas <- "/Users/bretbarnes/Documents/data/COVIC/fas/covic-prbU-all.fa.gz"
raw_all_bsp_rds <- "/Users/bretbarnes/Documents/data/COVIC/fas/covic-prbU-all.rds"
  
write_fas <- FALSE
if (write_fas) {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Building FASTA Raw...{RET}"))
  
  # Format EPIC Probes::
  man_pre_cor_tib <- man_raw_cor_tib %>% 
    dplyr::mutate(
      NB=dplyr::case_when(
        is.na(M) ~ "N",
        TRUE ~ nextBase_GRCh37),
      IF=dplyr::case_when(
        is.na(M) ~ 2,
        TRUE ~ 2
      ),
      IF=as.integer(IF),
      Probe_Tag=paste0(FR_GRCh37,NB,"C",IF,"1"),
      Probe_ID=paste(Probe_ID,Probe_Tag, sep='_')
    )
  
  fas_raw_dat_tib <- NULL
  fas_raw_dat_tib <-
    dplyr::bind_rows(man_pre_cor_tib,man_raw_dat_tib) %>%
    dplyr::mutate(FAS_SEQ=AlleleA_Probe_Sequence %>%
                    stringr::str_replace_all("R","A") %>%
                    stringr::str_replace_all("Y","T") %>%
                    stringr::str_to_upper(),
                  Design_Type=dplyr::case_when(
                    !is.na(U) & !is.na(M) ~ 1,
                    !is.na(U) &  is.na(M) ~ 2,
                    TRUE ~ 0),
                  Design_Type=as.integer(Design_Type),
                  FAS_KEY=paste(Probe_ID,Design_Type,Manifest, sep="_")
    ) %>%
    dplyr::select(FAS_KEY,FAS_SEQ,Design_Type,Manifest)
  
  fas_raw_all_vec <- fas_raw_dat_tib %>%
    dplyr::mutate(
      line=paste0('>',FAS_KEY,'\n',FAS_SEQ)
    ) %>%
    dplyr::filter(!is.na(line)) %>% 
    # dplyr::distinct(line) %>%
    dplyr::pull(line)
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Writing DATA Raw={fas_raw_all_csv}...{RET}"))
  readr::write_csv(fas_raw_dat_tib, fas_raw_all_csv)
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Writing FASTA Raw={fas_raw_all_fas}...{RET}"))
  readr::write_lines(fas_raw_all_vec, fas_raw_all_fas)
} else if (file.exists(fas_raw_all_csv)) {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading raw_all_tib={fas_raw_all_csv}...{RET}"))
  
  raw_all_tib <- readr::read_csv(fas_raw_all_csv) %>% # head() %>%
    dplyr::select(FAS_KEY,FAS_SEQ) %>%
    tidyr::separate(FAS_KEY, into=c("prb_cgn","prb_srd","prb_inf","prb_src"), sep="_") %>%
    dplyr::mutate(prb_inf=as.integer(prb_inf)) %>%
    dplyr::rename(prb_seq=FAS_SEQ)
} else {
  stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: File does not exist: {fas_raw_all_csv}.{RET}{RET}"))
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                 STEP 2::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$writeStep2 <- FALSE # Just a safety during testing to not delete large files...

if (file.exists(raw_all_bsp_rds)) {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading Probe GRS RDS={raw_all_bsp_rds}...{RET}"))
  
  raw_all_bsp_grs <- readr::read_rds(raw_all_bsp_rds)
  
} else {

  raw_bsp_tsv <- "/Users/bretbarnes/Documents/data/COVIC/bsmap/COVIC/GRCh37/covic-prbU-all/GRCh37.covic-prbU-all.v5.formatted.bsp.gz"
  raw_bsp_col <- cols(
    # bsp_key = col_character()
    prb_cgn = col_character(),
    prb_srd = col_character(),
    prb_inf = col_integer(),
    prb_src = col_character(),
    
    bsp_seq = col_character(),
    bsp_tag = col_character(),
    
    bsp_chr = col_character(),
    # bsp_cgn = col_character(),
    
    bsp_beg = col_integer(),
    bsp_srd = col_character(),
    bsp_mis = col_integer(),
    bsp_ref = col_character(),
    bsp_gap = col_integer(),
    
    # bsp_str = col_character()
    bsp_hit0 = col_integer(),
    bsp_hit1 = col_integer(),
    bsp_hit2 = col_integer(),
    bsp_hit3 = col_integer(),
    bsp_hit4 = col_integer(),
    bsp_hit5 = col_integer()
  )
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading raw_bsp_tib={raw_bsp_tsv}...{RET}"))
  
  #
  # Load BSP Formatted Data::
  #
  raw_bsp_tib <- readr::read_tsv(
    raw_bsp_tsv, 
    col_names=names(raw_bsp_col$cols),
    col_types=raw_bsp_col) # %>% dplyr::mutate(rvc_seq=revCmp(bsp_seq))
  
  
  # Join Source and Alignment Data::
  #  raw_all_tib %>% dplyr::anti_join(raw_bsp_tib, by=c("prb_cgn","prb_srd","prb_inf","prb_src"))
  raw_all_bsp_tib <- raw_all_tib %>% # head(n=1000) %>%
    dplyr::inner_join(raw_bsp_tib, 
                      by=c("prb_cgn","prb_srd","prb_inf","prb_src")) %>% 
    dplyr::mutate(
      # Generate RevComp for Probe Sequence
      prc_seq=revCmp(prb_seq),
      
      # Confirm Alignment Orientation from Probe Sequence/BSP Seq::
      prb_mat=dplyr::case_when(
        prb_seq==bsp_seq ~ "f",
        prc_seq==bsp_seq ~ "r",
        TRUE ~ NA_character_),
      
      # Extract Expected target CG di-nucleotride from Reference Sequence::
      CG_F1=stringr::str_sub(bsp_ref,-4,-3),
      CG_R1=stringr::str_sub(bsp_ref, 3, 4),
      CG_F2=stringr::str_sub(bsp_ref,-3,-2),
      CG_R2=stringr::str_sub(bsp_ref, 2, 3),
      CG_DIN=dplyr::case_when(
        prb_inf==1 & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ CG_F1,
        prb_inf==1 & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ CG_R1,
        prb_inf==2 & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ CG_F2,
        prb_inf==2 & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ CG_R2,
        TRUE ~ NA_character_
      ) %>% stringr::str_to_upper(),
      
      # Update Correct Genomic CG# Location based on alignment orientation::
      CG_Pos=dplyr::case_when(
        prb_inf==1 & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ bsp_beg +48,
        prb_inf==1 & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ bsp_beg + 0,
        prb_inf==2 & (bsp_srd=="--" | bsp_srd=="++") & prb_mat=="f" ~ bsp_beg +49,
        prb_inf==2 & (bsp_srd=="+-" | bsp_srd=="-+") & prb_mat=="r" ~ bsp_beg - 1,
        TRUE ~ NA_real_
      ) %>% as.integer()
    ) %>% 
    dplyr::group_by(prb_cgn,prb_srd,prb_inf,prb_src) %>% 
    dplyr::mutate(Unique_ID=paste(prb_cgn,prb_srd,prb_inf,prb_src, dplyr::row_number(), sep="_")) %>%
    dplyr::ungroup()
  
  # Quick Validation Summary::
  raw_all_bsp_tib %>% 
    dplyr::group_by(bsp_srd,prb_mat) %>% 
    dplyr::summarise(Count=n(), .groups="drop") %>% print(n=100)
  
  raw_all_bsp_tib %>%
    dplyr::select(prb_cgn,prb_srd,prb_inf,prb_src,prb_seq,bsp_ref,bsp_srd,prb_mat,CG_DIN) %>%
    dplyr::group_by(CG_DIN,prb_inf,bsp_srd,prb_mat) %>%
    dplyr::summarise(Count=n(), .groups="drop") %>% 
    dplyr::arrange(-Count) %>% print(n=100)
  
  #
  # Build BPS GRS::
  #
  raw_all_bsp_grs <- GRanges(
    seqnames = Rle(paste0("chr",raw_all_bsp_tib$bsp_chr)),
    strand=Rle(stringr::str_sub( raw_all_bsp_tib$bsp_srd, 1,1 ) ),
    
    prb_cgn=raw_all_bsp_tib$prb_cgn,
    prb_srd=raw_all_bsp_tib$prb_srd,
    prb_inf=raw_all_bsp_tib$prb_inf,
    prb_src=raw_all_bsp_tib$prb_src,
    
    bsp_seq1U = raw_all_bsp_tib$bsp_seq,
    bsp_ref   = raw_all_bsp_tib$bsp_ref,
    bsp_diNuc = raw_all_bsp_tib$CG_DIN,
    bsp_tag   = raw_all_bsp_tib$bsp_tag,
    bsp_srd   = raw_all_bsp_tib$bsp_srd,
    
    bsp_mis   = raw_all_bsp_tib$bsp_mis,
    bsp_gap   = raw_all_bsp_tib$bsp_gap,
    
    bsp_hit0=raw_all_bsp_tib$bsp_hit0,
    bsp_hit1=raw_all_bsp_tib$bsp_hit1,
    bsp_hit2=raw_all_bsp_tib$bsp_hit2,
    bsp_hit3=raw_all_bsp_tib$bsp_hit3,
    bsp_hit4=raw_all_bsp_tib$bsp_hit4,
    bsp_hit5=raw_all_bsp_tib$bsp_hit5,
    
    IRanges(start=raw_all_bsp_tib$CG_Pos,
            end=raw_all_bsp_tib$CG_Pos+1,
            names=raw_all_bsp_tib$Unique_ID,
            # names=paste(raw_all_bsp_tib$prb_cgn,
            #             raw_all_bsp_tib$prb_srd,
            #             raw_all_bsp_tib$prb_inf,
            #             raw_all_bsp_tib$prb_src,
            #             sep="_")
    )
  )

  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Writing Probe GRS RDS={raw_all_bsp_rds}...{RET}"))
  
  if (opt$writeStep2) readr::write_rds(raw_all_bsp_grs, raw_all_bsp_rds, compress="gz")
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                 STEP 3::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Latest Starting Place::
#
opt$writeStep3 <- FALSE # Just a safety during testing to not delete large files...

cgn_pos_db_rds  <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-end-seq.cgn-sorted.rds"
cgn_pos_db2_csv <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-end-seq.cgn-sorted.csv.gz"
cgn_pos_db_csv  <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-seq.cgn-sorted.csv"

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
    cat(glue::glue("[{par$prgmTag}]: Loading Original CGN CSV={cgn_pos_db2_csv}...{RET}"))
  
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
  
} else {
  stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: File does not exist: {cgn_pos_db_csv}.{RET}{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Intersect probes vs. database::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_imp_int_tib <- intersectGranges(man=raw_all_bsp_grs, ref=cgn_pos_grs,
                                    verbose=opt$verbose, tt=pTracker)


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Analyze probes vs. database Intersection::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Infinium I::
#
# Visual Example::
man_imp_int_tib %>% dplyr::filter(prb_inf==1) %>% head(n=8) %>%
  dplyr::filter(bsp_seq1U==imp_seq1U) %>% as.data.frame()


# Summary Example::
man_imp_int_tib %>% dplyr::filter(prb_inf==1) %>% 
  dplyr::filter(bsp_seq1U==imp_seq1U) %>% 
  dplyr::group_by(prb_inf,prb_src,bsp_diNuc,bsp_tag) %>%
  dplyr::summarise(Count=n(), .groups="drop")

all_sum1_tib <- man_imp_int_tib %>% dplyr::filter(prb_inf==1) %>% 
  dplyr::distinct(bsp_seq1U) %>% 
  dplyr::group_by(prb_inf,prb_src,bsp_diNuc,bsp_tag) %>%
  dplyr::summarise(Count=n(), .groups="drop")

mat_sum1_tib <- man_imp_int_tib %>% dplyr::filter(prb_inf==1) %>% 
  dplyr::filter(bsp_seq1U!=imp_seq1U) %>% 
  dplyr::distinct(bsp_seq1U)  %>% 
  dplyr::group_by(prb_inf,prb_src,bsp_diNuc,bsp_tag) %>%
  dplyr::summarise(Count=n(), .groups="drop")

mis_sum1_tib <- man_imp_int_tib %>% dplyr::filter(prb_inf==1) %>% 
  dplyr::filter(bsp_seq1U==imp_seq1U) %>% 
  dplyr::distinct(bsp_seq1U)  %>% 
  dplyr::group_by(prb_inf,prb_src,bsp_diNuc,bsp_tag) %>%
  dplyr::summarise(Count=n(), .groups="drop")




# Determine Unidentified Loci (i.e. missing from database)
man_imp_int_tib %>% dplyr::filter(prb_inf==1) %>% 
  dplyr::filter(bsp_seq1U==imp_seq1U) %>% 
  dplyr::filter( bsp_seq1U %in% unique(raw_all_bsp_grs$bsp_seq1U) ) %>%
  dplyr::distinct(bsp_seq1U)
  

# This is zero::
man_imp_int_tib %>% dplyr::filter(prb_inf==1) %>% 
  dplyr::filter(!bsp_seq1U %in% unique(raw_all_bsp_grs$bsp_seq1U) ) %>%
  dplyr::distinct(bsp_seq1U)

#
# Infinium II::
#


# Determine Unidentified Loci (i.e. missing from database)
man_imp_int_tib %>% dplyr::filter(prb_inf==2) %>% 
  dplyr::filter( bsp_seq1U %in% unique(raw_all_bsp_grs$bsp_seq1U) ) %>%
  dplyr::distinct(bsp_seq1U)

# This is zero::
man_imp_int_tib %>% dplyr::filter(prb_inf==2) %>% 
  dplyr::filter(!bsp_seq1U %in% unique(raw_all_bsp_grs$bsp_seq1U) ) %>%
  dplyr::distinct(bsp_seq1U)


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
