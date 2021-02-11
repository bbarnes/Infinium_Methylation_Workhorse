
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
par$prgmTag <- 'build_bead_manifest_simple_COVIC_v2'
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
    opt$genomeBuild <- 'GRCh37'
    opt$genomeBuild <- 'GRCh38'
    
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
  opt$runName <- paste(opt$genomeBuild,opt$platform,opt$version, sep='-')
  
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

opt$probe_type <- 'cg'
opt$design_key <- 'Seq_ID'
opt$design_seq <- 'Forward_Sequence'
opt$design_seq <- 'Top_Sequence'
opt$design_prb <- 'Probe_Type'
opt$design_srs <- 'TB'
opt$design_cos <- 'CO'

image_key <- "bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse_Centos"
image_ver <- "v.1.0"
image_ssh <- "run_improbe.sh"
image_str <- glue::glue("{image_key}.{image_ver}")

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
#                     1.0 Clean Tango Build:: COVIC
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (par$local_runType=="COVIC") {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #               1.1 Load Previous Manifests and AQP Results::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

  man_raw_dat_tib <- dplyr::bind_rows(
    # Previous EPIC Designs::
    loadCoreManifest_COVIC(
      datDir=par$datDir, manDir=par$manDir,
      verbose=opt$verbose,tc=0,tt=pTracker),
    
    # New COVIC Designs::
    loadAqpWorkflow_COVIC(
      ords=opt$ords, mats=opt$mats, aqps=opt$aqps,man="C0",
      verbose=opt$verbose,tc=0,tt=pTracker)
  ) %>% 
    dplyr::distinct(U,M, .keep_all=TRUE) %>%
    dplyr::arrange(Probe_ID) %>%
    dplyr::mutate(AQP=1, AQP=as.integer(AQP))
  
  man_add_dat_tib <- 
    manifestToAddress_COVIC(man_raw_dat_tib, verbose=opt$verbose,tc=0,tt=pTracker)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                 1.2 Load All Full IDAT Chip Validation::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  if (!is.null(par$idat_prefix)) {
    idat_tib <- 
      loadIdatAddress(prefix=par$idat_prefix, verbose=opt$verbose,tc=0,tt=pTracker)
    
    man_add_dat_tib <- man_add_dat_tib %>% 
      dplyr::mutate(isIDAT=dplyr::case_when(
        Address %in% idat_tib$Address ~ TRUE,
        TRUE ~ FALSE)
      )
    
    man_add_dat_tib %>% 
      dplyr::group_by(Probe_Type,Probe_Design,Manifest,isIDAT) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print()
  }

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       1.2 Summaries Raw Manifest::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  man_raw_sum_tib <- man_raw_dat_tib %>% 
    dplyr::group_by(Probe_Type,Manifest) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  man_raw_sum_tib %>% print(n=base::nrow(man_raw_sum_tib))
  
  # Only Select a subset of fields:: Skip this for now...
  # man_epic_tib %>% dplyr::select(dplyr::all_of( man_new_raw_tib %>% names() ))

} else {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         1.0 Process AQP/PQC Data::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  base_sel_cols <- c("Probe_ID","M","U","Probe_Type","Next_Base",
                     "AlleleA_Probe_Sequence","AlleleB_Probe_Sequence","Seq_48U",
                     "AQP")
  
  man_pqc_raw_tib <- 
    decodeAqpPqcWrapper(
      ord_vec=ords_vec,mat_vec=mats_vec,aqp_vec=aqps_vec,pqc_vec=pqcs_vec,
      platform=opt$platform,version=opt$version,
      ordFormat=opt$ordFormat,ordSkip=opt$ordSkip,
      matFormat=opt$matFormat,matSkip=opt$matSkip,
      aqpSkip=opt$aqpSkip,
      pqcSkip=opt$pqcSkip,
      name=opt$runName,outDir=opt$manDir,origin=opt$time_org_txt,
      fresh=opt$fresh,fixIds=opt$fixIds,full=par$retData,trim=TRUE,
      verbose=opt$verbose,vt=1,tc=0,tt=pTracker) %>%
    dplyr::select(dplyr::all_of(base_sel_cols))
  
  man_raw_dat_tib <- 
    dplyr::bind_rows(man_pqc_raw_tib,man_pre_dat_tib) %>%
    dplyr::select(dplyr::all_of(base_sel_cols))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                2.0 Build Probes foreach Type:: RS/CH/CG/etc..
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_prb_list <- NULL
man_raw_list <- man_raw_dat_tib %>% split(.$Probe_Type)

for (probe_type in rev(names(man_raw_list)) ) {
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Building Probes Types={probe_type}...{RET}"))
  
  if (probe_type == 'rs') {
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           2.1 Build Probes:: SNP
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    opt$non_s48_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/seq48U/un/SNP-13122020.seq48U.sorted.tsv'
    opt$non_top_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/cgnTop/SNP-13122020.cgnTop.sorted.tsv'
    
    man_prb_list[[probe_type]] <- NULL
    if (!is.null(opt$non_s48_tsv) && file.exists(opt$non_s48_tsv) &&
        !is.null(opt$non_top_tsv) && file.exists(opt$non_top_tsv)) {
      man_prb_list[[probe_type]] =
        clean_manifest_probes(
          tib=man_raw_dat_tib,
          s48_tsv=opt$non_s48_tsv,
          top_tsv=opt$non_top_tsv,
          name=opt$runName,
          outDir=opt$outDir,
          origin=opt$time_org_txt,
          
          design_key=opt$design_key,design_seq=opt$design_seq,
          design_prb=opt$design_prb,probe_type=probe_type,
          design_srs=opt$design_srs,design_cos=opt$design_cos,
          parallel=opt$parallel,fresh=opt$fresh,
          
          verbose=opt$verbose,vt=3,tc=0,tt=pTracker)
    }
    
  } else if (probe_type == 'ch') {
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           2.2 Build Probes:: CpH
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    opt$non_s48_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/seq48U/un/CpH-13122020.seq48U.sorted.tsv'
    opt$non_top_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/cgnTop/CpH-13122020.cgnTop.sorted.tsv'
    
    man_prb_list[[probe_type]] <- NULL
    if (!is.null(opt$non_s48_tsv) && file.exists(opt$non_s48_tsv) &&
        !is.null(opt$non_top_tsv) && file.exists(opt$non_top_tsv)) {
      man_prb_list[[probe_type]] =
        clean_manifest_probes(
          tib=man_raw_dat_tib,
          s48_tsv=opt$non_s48_tsv,
          top_tsv=opt$non_top_tsv,
          name=opt$runName,
          outDir=opt$outDir,
          origin=opt$time_org_txt,
          
          design_key=opt$design_key,design_seq=opt$design_seq,
          design_prb=opt$design_prb,probe_type=probe_type,
          design_srs=opt$design_srs,design_cos=opt$design_cos,
          parallel=opt$parallel,fresh=opt$fresh,
          
          verbose=opt$verbose,vt=3,tc=0,tt=pTracker)
    }

  } else if (probe_type == 'cg') {
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           2.3 Build Probes:: CpG
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    probe_type <- 'cg'
    man_prb_list[[probe_type]] <- NULL
    if (!is.null(opt$cpg_s48_tsv) && file.exists(opt$cpg_s48_tsv) &&
        !is.null(opt$cpg_top_tsv) && file.exists(opt$cpg_top_tsv)) {
      man_prb_list[[probe_type]] =
        clean_manifest_probes(
          tib=man_raw_dat_tib,
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
    next
  }
  
  # Extra QC Checking::
  #
  top_1_cnt <- 0
  top_1_cnt <- man_prb_list[[probe_type]] %>%
    dplyr::filter(!is.na(M),!is.na(U)) %>% base::nrow()
  if (top_1_cnt>0) {
    mat_u_tib <- man_prb_list[[probe_type]] %>% 
      dplyr::filter(AlleleA_Probe_Sequence==PRB1_U_MAT)
    mat_u_cnt <- mat_u_tib %>% base::nrow()
    nan_u_cnt <- mat_u_tib %>% dplyr::filter(is.na(M)) %>% base::nrow()
    if (opt$verbose>=1)
      cat(glue::glue("U: top_1_cnt={top_1_cnt}, mat_u_cnt={mat_u_cnt}, nan_u_cnt={nan_u_cnt}{RET}"))
    
    mat_m_tib <- man_prb_list[[probe_type]] %>% 
      dplyr::filter(AlleleB_Probe_Sequence==PRB1_M_MAT)
    mat_m_cnt <- mat_m_tib %>% base::nrow()
    nan_m_cnt <- mat_m_tib %>% dplyr::filter(is.na(U)) %>% base::nrow()
    if (opt$verbose>=1)
      cat(glue::glue("M: top_1_cnt={top_1_cnt}, mat_m_cnt={mat_m_cnt}, nan_m_cnt={nan_m_cnt}{RET}"))
  }
  
  top_2_cnt <- 0
  top_2_cnt <- man_prb_list[[probe_type]] %>%
    dplyr::filter(is.na(M),!is.na(U)) %>% base::nrow()
  if (top_2_cnt>0) {
    mat_2_tib <- man_prb_list[[probe_type]] %>% 
      dplyr::filter(AlleleA_Probe_Sequence==PRB2_D_MAT)
    mat_2_cnt <- mat_2_tib %>% base::nrow()
    nan_2_cnt <- mat_2_tib %>% dplyr::filter(is.na(M)) %>% base::nrow()
    if (opt$verbose>=1)
      cat(glue::glue("2: top_2_cnt={top_2_cnt}, mat_2_cnt={mat_2_cnt}, nan_2_cnt={nan_2_cnt}{RET}"))
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
# TBD:: Unclear if we really need the function below::
#  - should be cleaned up...
#
# TBD:: The next step should be wrapped up into a previous general function
#   for probe names...
#
# TBD:: The binding should be piped into a general function to add
#   inferred columns
#

man_pqc_prb_tib <- bindProbeDesignList(
  list=man_prb_list, # platform=opt$platform, version=opt$version,
  sumDir=opt$sumDir,del='.',
  verbose=opt$verbose,vt=3,tc=1,tt=pTracker) %>% 
  dplyr::mutate(Assay_Class='Analytical')

core_prb_type_vec <- man_pqc_prb_tib %>% 
  dplyr::distinct(Probe_Type) %>% dplyr::pull(Probe_Type)

# Clear up memory::
# rm(man_prb_list)

# man_pqc_prb_tib %>% dplyr::group_by(Platform, Version) %>% dplyr::summarise(Count=n(), .groups="drop")
# man_pqc_prb_tib %>% dplyr::group_by(Platform, Manifest) %>% dplyr::summarise(Count=n(), .groups="drop")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#            Extract Missing Probes from Manifest:: i.e. CTL/UNK
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_pqc_unk_tib <- NULL
man_pqc_all_tib <- man_pqc_prb_tib %>% 
  tidyr::separate(Probe_ID, into=c("Source_ID","Source_Tag"), 
                  sep="_", remove=FALSE) 

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       0.5 Manifest Summary: S4/C0/B4/B2
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

epic_b2_tag_tib <- dplyr::bind_rows( epic_b4_tib, epic_b2_tib ) %>% 
  dplyr::mutate(Platform="EPIC",Probe_Source="EPIC") %>%
  dplyr::distinct(M,U, .keep_all=TRUE) %>% 
  dplyr::mutate(
    DESIGN=Infinium_Design,
    Infinium_Design=dplyr::case_when(Infinium_Design=="I" ~ "1", 
                                     Infinium_Design=="II" ~ "2", 
                                     TRUE ~ NA_character_), 
    Infinium_Design=as.integer(Infinium_Design) ) %>% 
  dplyr::left_join(man_pqc_all_tib %>% dplyr::select(M,U,PRB1_U_MAT), by=c("M","U"))


epic_b2_sum_tib <- epic_b2_tag_tib %>% 
  dplyr::group_by(Man_Source,Probe_Type,Infinium_Design,
                  Probe_Source,Platform,Version) %>% 
  dplyr::summarise(Count=n(), .groups="drop")

# Quick QC::
epic_b2_ses_mis_cnt <- epic_S4_tib %>% 
  dplyr::anti_join(epic_b2_tag_tib, by=c("M","U")) %>% base::nrow()
cat(glue::glue("[{par$prgmTag}]: epic_b2_ses_mis_cnt={epic_b2_ses_mis_cnt}.{RET}"))


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        0.6 IDATs Stats Summary: S4/C0/B4/B2
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Compare IDATS to Original MAT/AQP and ORD/MAT/AQP Passing::
mis_add_idats_raw_tib <- covic_mat_aqp_raw_tab %>% dplyr::anti_join(idat_tib, by="Address")
mis_add_idats_pas_tib <- covic_mat_aqp_pas_tab %>% dplyr::anti_join(idat_tib, by="Address")
covic_ann_aqp_tab %>% dplyr::anti_join(idat_tib, by="Address")

# Compare new AQP data::
man_pqc_all_tib %>% dplyr::filter(!is.na(U)) %>%
  dplyr::anti_join(idat_tib, by=c("U"="Address"))

# Summary::
man_pqc_all_tib %>% dplyr::filter(!is.na(U)) %>%
  dplyr::anti_join(idat_tib, by=c("U"="Address")) %>% 
  dplyr::group_by(Probe_Type,DESIGN) %>% 
  dplyr::summarise(Count=n(), .groups="drop")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       0.7 Build Label Table: S4/C0/B4/B2
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
# Build labeld M/U Versions::
#
#   RECALCULATE FOR ACCURACY!!!
#
man_ver_all_tib <- dplyr::bind_rows(
  epic_S4_tib %>% dplyr::select(Probe_ID,Probe_Type,U,M) %>% dplyr::mutate(M=as.integer(M), U=as.integer(U), Version="S4"),
  epic_b4_tib %>% dplyr::select(Probe_ID,Probe_Type,U,M) %>% dplyr::mutate(M=as.integer(M), U=as.integer(U), Version="B4"),
  epic_b2_tib %>% dplyr::select(Probe_ID,Probe_Type,U,M) %>% dplyr::mutate(M=as.integer(M), U=as.integer(U), Version="B2"),
  man_pqc_all_tib %>% dplyr::select(Probe_ID,Probe_Type,U,M) %>% dplyr::mutate(M=as.integer(M), U=as.integer(U), Version="C0"),
) %>% dplyr::distinct(M,U, .keep_all=TRUE)

# Compare Labeled Sources vs. Idats:: Tables
man_src_mat_idats_tib <- dplyr::bind_rows(
  man_ver_all_tib %>% 
    dplyr::filter(!is.na(M)) %>% dplyr::inner_join(idat_tib, by=c("M"="Address")),
  man_ver_all_tib %>% 
    dplyr::filter(!is.na(U)) %>% dplyr::inner_join(idat_tib, by=c("U"="Address"))
) %>% dplyr::distinct(M,U, .keep_all=TRUE)

man_src_mis_idats_tib <- dplyr::bind_rows(
  man_ver_all_tib %>% 
    dplyr::filter(!is.na(M)) %>% dplyr::anti_join(idat_tib, by=c("M"="Address")),
  man_ver_all_tib %>% 
    dplyr::filter(!is.na(U)) %>% dplyr::anti_join(idat_tib, by=c("U"="Address"))
) %>% dplyr::distinct(M,U, .keep_all=TRUE)

# Compare Labeled Sources vs. Idats:: Summary
man_src_mat_idats_sum <- man_src_mat_idats_tib %>%
  dplyr::group_by(Version) %>% dplyr::summarise(Count=n(), .groups="drop")
man_src_mis_idats_sum <- man_src_mis_idats_tib %>%
  dplyr::group_by(Version) %>% dplyr::summarise(Count=n(), .groups="drop")


#
# Validation Stats::
#
man_pqc_all_cnt <- man_pqc_all_tib %>% base::nrow()
man_pqc_dat_cnt <- man_raw_dat_tib %>% base::nrow()
man_dif_all_cnt <- man_pqc_all_cnt - man_pqc_dat_cnt

man_pqc_prb_cnt <- man_pqc_prb_tib %>% base::nrow()
man_pqc_unk_cnt <- man_pqc_unk_tib %>% base::nrow()
man_dif_unk_cnt <- man_pqc_dat_cnt - man_pqc_prb_cnt - man_pqc_unk_cnt

cat(glue::glue("[{par$prgmTag}]: Done; Parsing PQC/All Probes={man_pqc_all_cnt}.{RET}"))
cat(glue::glue("[{par$prgmTag}]: Done; Parsing PQC/Dat Probes={man_pqc_dat_cnt}.{RET}"))
cat(glue::glue("[{par$prgmTag}]: Done; Parsing DiffCnt Probes={man_dif_all_cnt}.{RET}{RET}"))

cat(glue::glue("[{par$prgmTag}]: Done; Parsing PQC/Prb Probes={man_pqc_prb_cnt}.{RET}"))

if (! is.null(man_pqc_unk_tib)) {
  cat(glue::glue("[{par$prgmTag}]: Done; Parsing Unk/Ctl Probes={man_pqc_unk_cnt}.{RET}"))
  cat(glue::glue("[{par$prgmTag}]: Done; Parsing DiffCnt Probes={man_dif_unk_cnt}.{RET}{RET}"))
  
  man_unk_sum_tib <- man_pqc_unk_tib %>% 
    dplyr::group_by(Assay_Class,Probe_Class,Probe_Type,Infinium_Design,AQP) %>% 
    dplyr::summarise(Count=n(), .groups='drop')
  print(man_unk_sum_tib,n=base::nrow(man_unk_sum_tib))
}

man_all_sum_tib <- man_pqc_all_tib %>% 
  dplyr::group_by(Assay_Class,Probe_Class,Probe_Type,Infinium_Design,AQP) %>% 
  dplyr::summarise(Count=n(), .groups='drop')
print(man_all_sum_tib,n=base::nrow(man_all_sum_tib))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#            4.0 Load all Coordinates for Detected Manifest Probes::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# NOT NEEDED ANY MORE:: Better method below...
#
if (FALSE) {
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
}

#
# Better method::
#

cgn_pos1_dir <- file.path(opt$outDir, 'intersection/cg')
cgn_pos1_pre <- paste(opt$genomeBuild,par$local_runType,opt$version, sep='-')
cgn_pos1_tsv <- file.path(cgn_pos1_dir, paste(cgn_pos1_pre,'seq48U_to_cgn.int.Seq_48U-sorted.tsv.gz', sep='.'))
out_pos1_tsv <- file.path(cgn_pos1_dir, paste(cgn_pos1_pre,'seq48U_to_cgn.int.cg-sorted.tsv', sep='.'))
out_pos3_tsv <- file.path(cgn_pos1_dir, paste(cgn_pos1_pre,'seq48U_to_cgn.int.improbe-designOutput.cgn-map-seq.cgn-sorted.tsv.gz', sep='.'))

imp_pos3_dir <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020"
cgn_pos3_fns <- paste(opt$genomeBuild,"21092020_improbe-designOutput.cgn-map-seq.cgn-sorted.tsv.gz", sep='-' )
imp_pos3_tsv <- file.path(imp_pos3_dir, cgn_pos3_fns)

# cgn_pos1_tsv <- file.path(opt$outDir, 'intersection/cg/GRCh38-COVIC-C4.seq48U_to_cgn.int.Seq_48U-sorted.tsv.gz')
# out_pos1_tsv <- file.path(opt$outDir, 'intersection/cg/GRCh38-COVIC-C4.seq48U_to_cgn.int.cg-sorted.tsv')
# out_pos3_tsv <- file.path(opt$outDir, 'intersection/cg/GRCh38-COVIC-C4.seq48U_to_cgn.int.GRCh38.improbe-designOutput.cgn-map-seq.cgn-sorted.tsv.gz')
# imp_pos3_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.cgn-sorted.tsv.gz'

if (!file.exists(out_pos1_tsv)) {
  # cmd_1 <- "gzip -dc /Users/bretbarnes/Documents/scratch/build_bead_manifest_simple/GRCh38-COVIC-C4/intersection/cg/GRCh38-COVIC-C4.seq48U_to_cgn.int.Seq_48U-sorted.tsv.gz | sort -k 2,2 > /Users/bretbarnes/Documents/scratch/build_bead_manifest_simple/GRCh38-COVIC-C4/intersection/cg/GRCh38-COVIC-C4.seq48U_to_cgn.int.cg-sorted.tsv"
  int_pos1_cmd <- glue::glue("gzip -dc {cgn_pos1_tsv} | sort -k 2,2 > {out_pos1_tsv}")
  if (opt$verbose>=3)
    cat(glue::glue("[{par$prgmTag}]: Running; cmd={int_pos1_cmd}...{RET}"))
  system(int_pos1_cmd)
}

if (opt$verbose>=3)
  cat(glue::glue("[{par$prgmTag}]: Loading Intersection; out_pos1_tib={out_pos1_tsv}...{RET}"))
out_pos1_col <- cols(
  Seq_48U = col_character(),
  CGN_Imp = col_character(),
  TB = col_character(),
  CO = col_character(),
  M = col_integer(),
  U = col_integer() )
out_pos1_tib <- readr::read_tsv(out_pos1_tsv, col_names=names(out_pos1_col$cols), col_types=out_pos1_col)

if (!file.exists(out_pos3_tsv)) {
  # Build better database::: Combine the two commands below::
  # gzip -dc /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.tsv.gz | cut -f 1,4,5,21-23,45,48 > /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.tsv
  # cat /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.tsv | perl -pe 's/TOP/T/; s/BOT/B/;' | gzip -c - > /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.tsv.gz
  
  int_pos3_cmd <- glue::glue("gzip -dc {imp_pos3_tsv} | join -t $'\t' -11 -22 - {out_pos1_tsv} | gzip -c - > {out_pos3_tsv}")
  if (opt$verbose>=3)
    cat(glue::glue("[{par$prgmTag}]: Running; cmd={int_pos3_cmd}...{RET}"))
  system(int_pos3_cmd)
}

if (opt$verbose>=3)
  cat(glue::glue("[{par$prgmTag}]: Loading Coordinates; out_pos1_tib={out_pos3_tsv}...{RET}"))
imp_pos3_col <- cols(
  Seq_ID = col_character(),
  Chromosome = col_character(),
  Coordinate = col_integer(),
  # FR_DB1 = col_character(),
  # TB_DB1 = col_character(),
  # CO_DB1 = col_character(),
  # NB_DB1 = col_character(),
  FR_IMP = col_character(),
  TB_IMP = col_character(),
  CO_IMP = col_character(),
  NB_IMP = col_character(),
  Seq50U = col_character(),
  
  Seq48U = col_character(),
  TB_DB2 = col_character(),
  CO_DB2 = col_character(),
  M = col_integer(),
  U = col_integer()
)
imp_pos3_tib <- readr::read_tsv(out_pos3_tsv, col_names=names(imp_pos3_col$cols), col_types=imp_pos3_col)





# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             4.1 Match All Probes To CG# Database Coordinates::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
# Data Structures:: 
#   All DES = man_pqc_all_tib
#   All POS = man_pqc_grs_tib
man_pqc_grs_tib <- man_pqc_all_tib %>% 
  dplyr::filter(!is.na(PRB1_U_MAT)) %>%
  dplyr::inner_join(imp_pos3_tib, 
                    by=c("PRB1_U_MAT"="Seq50U"), 
                    suffix=c("_PQC", "_IMP")) %>% 
  dplyr::select(Probe_ID,Chromosome,Coordinate, dplyr::everything()) %>% 
  dplyr::distinct(Probe_ID,Chromosome,Coordinate, .keep_all=TRUE) %>%
  dplyr::rename(chrom=Chromosome,
                chromStart=Coordinate) %>%
  dplyr::mutate(chrom=stringr::str_remove(chrom,'^chr'),
                chrom=paste0('chr',chrom),
                chrom=dplyr::case_when(
                  chrom=="chrMT" ~ "chrM",
                  TRUE ~ chrom
                ),
                chromEnd=chromStart+1) %>%
  dplyr::arrange(chrom,chromStart) %>%
  dplyr::select(chrom,chromStart,chromEnd,Probe_ID, dplyr::everything())

# chrM QC::
# man_pqc_grs_tib %>% dplyr::filter(chrom=="chrM")  # Should be non-zero
# man_pqc_grs_tib %>% dplyr::filter(chrom=="chrMT") # Should be zero

# Collapse down redundant fields::
#
# man_pqc_grs_tib %>% dplyr::filter(Seq_48U != Seq48U) %>% dplyr::select(Seq_ID_IMP,PRB1_U_MAT,Seq_48U,Seq48U)

min_pqc_grs_tib <- 
  man_pqc_all_tib %>% dplyr::select(-Probe_Source,-Version) %>%
  dplyr::left_join(man_ver_all_tib %>% dplyr::rename(M_Lab=M), 
                   by=c("U","Probe_Type"), suffix=c("_PQC","_LAB") ) %>%
  dplyr::filter(!is.na(PRB1_U_MAT)) %>%
  
  dplyr::inner_join(imp_pos3_tib %>% dplyr::distinct(Seq_ID,Chromosome,Coordinate,FR_IMP,TB_IMP,CO_IMP,NB_IMP,Seq50U),
                    by=c("PRB1_U_MAT"="Seq50U"), 
                    suffix=c("_PQC", "_IMP")) %>% 
  # dplyr::select(Probe_ID,Chromosome,Coordinate, dplyr::everything()) %>% 
  # dplyr::distinct(Probe_ID,Chromosome,Coordinate, .keep_all=TRUE) %>%
  dplyr::rename(chrom=Chromosome,
                chromStart=Coordinate) %>%
  dplyr::mutate(chrom=stringr::str_remove(chrom,'^chr'),
                chrom=paste0('chr',chrom),
                chrom=dplyr::case_when(
                  chrom=="chrMT" ~ "chrM",
                  TRUE ~ chrom
                ),
                chromEnd=chromStart+1) %>%
  dplyr::arrange(chrom,chromStart) %>%
  dplyr::select(chrom,chromStart,chromEnd, 
                Seq_ID_IMP,FR_IMP,TB_IMP,CO_IMP,NB_IMP,
                PRB1_U_MAT,
                dplyr::everything())

min_pqc_grs_tib %>% 
  dplyr::group_by(Platform,Version,Probe_Type) %>% 
  dplyr::summarise(Count=n(), .groups="drop")

unq_pqc_cnt_tib <- min_pqc_grs_tib %>% 
  dplyr::add_count(U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                   chrom,chromStart,chromEnd,Seq_ID_IMP,TB_IMP,CO_IMP, name="Unq_Cnt")

unq_pqc_cnt_tib %>% dplyr::filter(Unq_Cnt != 1) %>% 
  dplyr::distinct(U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                  chrom,chromStart,chromEnd,Seq_ID_IMP,TB_IMP,CO_IMP, 
                  Probe_Type,Platform, Version, Unq_Cnt) %>%
  dplyr::group_by(Platform,Version,Probe_Type) %>% 
  dplyr::summarise(Count=n(), .groups="drop")

#
# Below we really only need to check for Tango U matching!!!
#
if (FALSE) {
  # Check for completely missing tango M/U pairs::
  man_pqc_sum_tib <- man_pqc_all_tib %>% 
    dplyr::anti_join(man_pqc_grs_tib, by=c("M"="M_PQC","U"="U_PQC")) %>% 
    dplyr::anti_join(man_pqc_grs_tib, by=c("M"="M_IMP","U"="U_IMP")) %>%
    dplyr::group_by(Probe_Type,DESIGN,Probe_Source,Version) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  
  man_imp_sum_tib <- man_pqc_all_tib %>% 
    dplyr::anti_join(man_pqc_grs_tib, by=c("M"="M_IMP","U"="U_IMP")) %>% 
    # dplyr::anti_join(man_pqc_grs_tib, by=c("M"="M_PQC","U"="U_PQC")) %>%
    dplyr::group_by(Probe_Type,DESIGN,Probe_Source,Version) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  
  # Flip Comparison::
  mis_pqc_grs_tib <- man_pqc_grs_tib %>%
    dplyr::anti_join(man_pqc_all_tib, by=c("M_IMP"="M","U_IMP"="U")) %>%
    # dplyr::anti_join(man_pqc_all_tib, by=c("M_PQC"="M","U_PQC"="U")) %>%
    dplyr::filter(Seq_48U==Seq48U | AlleleA_Probe_Sequence==PRB2_D_MAT) 
  
  man_pqc_all_tib %>% dplyr::inner_join(
    mis_pqc_grs_tib %>% dplyr::select(AlleleA_Probe_Sequence), by="AlleleA_Probe_Sequence") %>%
    dplyr::inner_join(imp_pos3_tib, 
                      by=c("Seq_48U"="Seq48U"),
                      suffix=c("_PQC", "_IMP")) %>%
    dplyr::select(Probe_ID,Chromosome,Coordinate, dplyr::everything()) %>% 
    dplyr::distinct(Probe_ID,Chromosome,Coordinate, .keep_all=TRUE) %>%
    dplyr::rename(chrom=Chromosome,
                  chromStart=Coordinate) %>%
    dplyr::mutate(chrom=stringr::str_remove(chrom,'^chr'),
                  chrom=paste0('chr',chrom),
                  chromEnd=chromStart+1) %>%
    dplyr::arrange(chrom,chromStart) %>%
    dplyr::select(chrom,chromStart,chromEnd,Probe_ID, dplyr::everything()) %>%
    dplyr::select(# chrom,chromStart,chromEnd,
      Probe_ID,PRB1_U_MAT,Seq50U, U_PQC,U_IMP, M_PQC,M_IMP)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 7.1 Extract Predefined Sesame Mappings: S4
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Data Structures:: 
#   All DES = man_pqc_all_tib
#   All POS = man_pqc_grs_tib
#
# TBD:: Need to include man_ver_all_tib labels!!!

epic_S4_grs_tib <- epic_S4_tib %>% 
  dplyr::select(Probe_ID,
                chrom,chromStart,chromEnd,strand,DESIGN,
                U,M,ProbeSeq_A,ProbeSeq_B,
                DESIGN,gene,MASK_mapping,MASK_general) %>% 
  dplyr::arrange(Probe_ID)


# Basic re-ordr of probes based on man_ver_all_tib
man_S4_PQC_grs_intPrb_tib <- epic_S4_grs_tib %>% 
  dplyr::inner_join(man_pqc_grs_tib, 
                    by=c("ProbeSeq_A"="AlleleA_Probe_Sequence", "DESIGN"),
                    suffix=c("_SES", "_PQC") ) %>%
  dplyr::mutate(begDif=chromStart_SES-chromStart_PQC,
                endDif=chromEnd_SES-chromEnd_PQC,
                cgFlag=dplyr::case_when(
                  Probe_ID_SES==Seq_ID_PQC ~ 0,
                  TRUE ~ 1
                ),
                chrFlag=dplyr::case_when(
                  chrom_SES==chrom_PQC ~ 0,
                  TRUE ~ 1
                ),
                absDif=base::abs(begDif)) %>%
  dplyr::arrange(ProbeSeq_A, chrFlag, absDif, cgFlag)

# Visual Inspection::
#
man_S4_PQC_grs_intPrb_tib %>%
  dplyr::select(Probe_ID_SES,Seq_ID_IMP,Seq_ID_PQC,Probe_ID_PQC,
                chrom_SES, chrom_PQC,
                chromStart_SES, chromStart_PQC, begDif,
                # chromEnd_PQC, chromEnd_SES,  endDif,
                ProbeSeq_A,Seq48U,Seq_48U)

# hg38: 861,533 unique ProbeSeq_A
# hg19: 862,799 unique ProbeSeq_A
cur_mat_prb_tib <- man_S4_PQC_grs_intPrb_tib %>% 
  dplyr::filter(begDif == 0 & endDif == 0 & chrom_SES==chrom_PQC) %>% 
  dplyr::distinct(ProbeSeq_A, .keep_all=TRUE)

cur_mat_prb_tib %>%
  dplyr::group_by(chrFlag,absDif,cgFlag,MASK_mapping,MASK_general) %>% 
  dplyr::summarise(Count=n(), .groups="drop")


# Missing Probes: This has to be unique from cur_mat_prb_tib
cur_mis_prb_tib <- man_S4_PQC_grs_intPrb_tib %>% 
  dplyr::filter(! ProbeSeq_A %in% cur_mat_prb_tib$ProbeSeq_A)
  # Old Filtering Method Below:: replaced with line above!!!
  # dplyr::filter(begDif != 0 | endDif != 0 | chrom_SES!=chrom_PQC) %>% 
  # dplyr::arrange(begDif) %>%
  # dplyr::distinct(ProbeSeq_A, .keep_all=TRUE)

cur_mis_prb_tib %>%
  dplyr::group_by(chrFlag,absDif,cgFlag,MASK_mapping,MASK_general) %>% 
  dplyr::summarise(Count=n(), .groups="drop") %>%
  print(n=base::nrow(cur_mis_prb_tib))

# chrM Probes, total: 8
cur_mis_prb_tib %>% 
  dplyr::filter(absDif < 4 & chrFlag==0) %>% 
  dplyr::select(Probe_ID_SES, Probe_ID_PQC, DESIGN, 
                chrFlag,absDif,cgFlag,MASK_mapping,MASK_general, 
                chromStart_SES,chromStart_PQC)

# Example of cgn's not matching (hg38):: 1,701
# Example of cgn's not matching (hg19)::     8
cur_mat_prb_tib %>% 
  dplyr::filter(cgFlag == 1) %>%   
  dplyr::select(Probe_ID_SES,Seq_ID_IMP,Seq_ID_PQC,Probe_ID_PQC,
                chrom_SES, chrom_PQC,
                chromStart_SES, chromStart_PQC, begDif,
                # chromEnd_PQC, chromEnd_SES,  endDif,
                ProbeSeq_A,Seq48U,Seq_48U) %>%
  dplyr::distinct(ProbeSeq_A, .keep_all=TRUE)

if (FALSE) {
  # Check cgn mismatches:: ProbeSeq_A=AAAAAAAAAAATTCTACTACTAAAAACCTCCTCCCCCTCCAAAAATATAC
  epic_b2_tib %>% dplyr::filter(Probe_ID=="cg11769960") %>% dplyr::select(Forward_Sequence)
  # AAGAGGGACGGGAGGAGGGAGATTCTGCTGCTAAGAGCCTCCTCCCCCTCCAAAGATGTG[CG]CTGCGCTTGGGACAAACTTCCAGGTCAAGAAAGTACGGGGGGACAATGGCAACAGCCGCT
  man_pqc_all_tib %>% dplyr::filter(Seq_ID=="cg22471585") %>% dplyr::select(Top_Sequence)
  # AAGAGGGACGGGAGGAGGGAGATTCTGCTGCTAAGAGCCTCCTCCCCCTCCAAAGATGTG[CG]CTGCGCTTGGGACAAACTTCCAGCTCAAGAAAGTACGGGGGGACAATGGCAACAGCCGCT
  
  man_pqc_all_tib %>% dplyr::filter(AlleleA_Probe_Sequence=="AAAAAAAAAAATTCTACTACTAAAAACCTCCTCCCCCTCCAAAAATATAC")
  man_pqc_grs_tib %>% dplyr::filter(AlleleA_Probe_Sequence=="AAAAAAAAAAATTCTACTACTAAAAACCTCCTCCCCCTCCAAAAATATAC")
  
  # How to get any TB/CO
  out_pos1_tib %>% dplyr::filter(CGN_Imp=="cg11769960" | CGN_Imp=="cg22471585")
  
  # Seq_48U                                          CGN_Imp    TB    CO        M        U
  # <chr>                                            <chr>      <chr> <chr> <int>    <int>
  # 1 AAAAAAAAATTCTACTACTAAAAACCTCCTCCCCCTCCAAAAATATAC cg11769960 B     C        NA 59807194
  # 2 AAAAAAAAATTCTACTACTAAAAACCTCCTCCCCCTCCAAAAATATAC cg22471585 B     C        NA 59807194
  
  # Recover TB/CO for cgn mismatches (hg38):: 1,697
  # Recover TB/CO for cgn mismatches (hg19)::     0
  cur_mat_prb_tib %>% 
    dplyr::filter(cgFlag == 1) %>%   
    dplyr::select(Probe_ID_SES,Seq_ID_IMP,Seq_ID_PQC,Probe_ID_PQC,
                  chrom_SES, chrom_PQC,
                  chromStart_SES, chromStart_PQC, begDif,
                  # chromEnd_PQC, chromEnd_SES,  endDif,
                  ProbeSeq_A,Seq48U,Seq_48U) %>% 
    dplyr::inner_join(out_pos1_tib, by="Seq_48U") %>% 
    dplyr::distinct(ProbeSeq_A, .keep_all=TRUE) %>%
    dplyr::filter(Probe_ID_SES==CGN_Imp) %>%
    dplyr::distinct(ProbeSeq_A, .keep_all=TRUE)
  
  # Re-run command with two missing example test case::
  # gzip -dc /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.cgn-sorted.tsv.gz | join -t $'  ' -11 -22 - tmp/two-int-seqs.tsv > tmp/two-int-seqs.results.tsv
  #
  # Re-run single missing example test case::
  # gzip -dc /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.cgn-sorted.tsv.gz | join -t $'  ' -11 -22 - tmp/one-int-seqs.tsv > tmp/one-int-seqs.results.tsv
  #
}






#
# OLD STUFF BELOW::
#




man_S4_PQC_grs_intPrb_tib %>% 
  dplyr::filter(ProbeSeq_A %in% cur_mis_prb_tib$ProbeSeq_A) %>% 
  dplyr::arrange(ProbeSeq_A) %>%
  dplyr::select(Probe_ID_SES,Seq_ID_IMP,Seq_ID_PQC,Probe_ID_PQC,
                chrom_SES, chrom_PQC,
                chromStart_SES, chromStart_PQC, begDif,
                # chromEnd_PQC, chromEnd_SES,  endDif,
                ProbeSeq_A)
                # gene,DESIGN)


epic_S4_grs_tib %>% 
  dplyr::filter(ProbeSeq_A %in% cur_mis_prb_tib$ProbeSeq_A) %>% 
  dplyr::group_by(MASK_mapping,MASK_general) %>% 
  dplyr::summarise(Count=n(), .groups="drop")

man_S4_PQC_grs_intPrb_tib %>% 
  dplyr::filter(ProbeSeq_A %in% cur_mis_prb_tib$ProbeSeq_A) %>% 
  dplyr::distinct(ProbeSeq_A, .keep_all=TRUE) %>%
  dplyr::group_by(MASK_mapping,MASK_general) %>% 
  dplyr::summarise(Count=n(), .groups="drop")

epic_S4_grs_tib %>% 
  dplyr::filter(ProbeSeq_A %in% cur_mat_prb_tib$ProbeSeq_A) %>% 
  dplyr::group_by(MASK_mapping,MASK_general) %>% 
  dplyr::summarise(Count=n(), .groups="drop")


cur_mat_pos_tib <- man_S4_PQC_grs_intPrb_tib %>% 
  # dplyr::distinct(ProbeSeq_A, .keep_all=TRUE) %>%
  dplyr::anti_join(cur_mis_prb_tib, by="ProbeSeq_A") %>% 
  # dplyr::anti_join(cur_mat_prb_tib, by="ProbeSeq_A") %>% 
  dplyr::select(Probe_ID_SES,Seq_ID_IMP,Seq_ID_PQC,Probe_ID_PQC,
                # chrom_SES, chrom_PQC, 
                # chromStart_SES, chromStart_PQC, begDif,
                # chromEnd_PQC, chromEnd_SES,  endDif,
                begDif,gene,
                
                # strand,
                # TB_DB1,TB_DB2,Strand_TB,
                # CO_DB1,CO_DB2,Strand_CO,CO,
                
                U,U_IMP,U_PQC, M,M_IMP,M_PQC,
                ProbeSeq_A,ProbeSeq_B,AlleleB_Probe_Sequence,
                
                DESIGN # DESIGN_SES, DESIGN_PQC,
  )

cur_mat_pos_tib %>% 
  # dplyr::filter(Probe_ID_SES == Seq_ID_IMP) %>%
  dplyr::distinct(Probe_ID_SES,ProbeSeq_A, .keep_all=TRUE)

# BLAT Source Sequece 
#  9 cg00699997   cg00699997 cg00699997 cg00699997_BC11 chr19     chr1               27109958      121743261  -94633303    121743262     27109959  -94633303 I     
# dplyr::filter(Probe_ID_SES=="cg00699997")


man_S4_PQC_grs_intPrb_tib %>% 
  dplyr::select(Probe_ID_SES,Seq_ID_IMP,Seq_ID_PQC,Probe_ID_PQC,
                chrom_SES, chrom_PQC, 
                chromStart_SES, chromStart_PQC, begDif,
                chromEnd_PQC, chromEnd_SES,  endDif,
                
                # strand,
                # TB_DB1,TB_DB2,Strand_TB,
                # CO_DB1,CO_DB2,Strand_CO,CO,
                
                # U,U_IMP,U_PQC, M,M_IMP,M_PQC,
                # ProbeSeq_A,ProbeSeq_B,AlleleB_Probe_Sequence,
                
                DESIGN # DESIGN_SES, DESIGN_PQC,
  )  %>% 
  dplyr::filter(begDif == 0 & endDif == 0) %>% 
  # dplyr::filter(begDif != 0 | endDif != 0) %>% # dplyr::filter(Probe_ID_SES != Seq_ID_PQC & Probe_ID_SES != Seq_ID_IMP)
  dplyr::distinct(ProbeSeq_A)


#  dplyr::distinct(strand,DESIGN, .keep_all=TRUE)

man_S4_PQC_grs_intA_tib <- epic_S4_grs_tib %>% 
  dplyr::inner_join(man_pqc_all_tib, by="U", suffix=c("_SES", "_IMP")) %>% 
  dplyr::inner_join(idat_tib, by=c("U"="Address"))

man_S4_PQC_grs_intS_tib <- epic_S4_grs_tib %>% 
  dplyr::inner_join(man_pqc_all_tib, 
                    by=c("ProbeSeq_A"="AlleleA_Probe_Sequence"),
                    suffix=c("_SES", "_IMP")) %>%
  dplyr::rename(AlleleA_Probe_Sequence=ProbeSeq_A)

man_S4_PQC_grs_intS_tib %>% 
  dplyr::inner_join(idat_tib, by=c("U_IMP"="Address"))
man_S4_PQC_grs_intS_tib %>% 
  dplyr::inner_join(idat_tib, by=c("U_SES"="Address"))

# THIS IS ALL FINE BELOW, Just use Sesame hg38 Probe_IDs
man_S4_PQC_grs_intS_tib %>% 
  dplyr::anti_join(idat_tib, by=c("U_SES"="Address"))

man_S4_PQC_grs_intS_tib %>% 
  dplyr::anti_join(idat_tib, by=c("U_IMP"="Address")) %>%
  dplyr::left_join(man_ver_all_tib, by=c("U_SES"="U"), suffix=c("_MIS","_LAB")) %>%
  dplyr::select(Probe_ID_SES,Probe_ID_IMP, U_SES,U_IMP, Probe_ID,Version_LAB,Probe_Type_LAB)


#
# NOW:: We should use this set to build hg38 coordindates:: man_S4_PQC_grs_intS_tib
#


# Some QC Checks::
if (FALSE) {
  man_pqc_all_tib %>% dplyr::anti_join(man_S4_PQC_grs_tib, by="U")
  man_pqc_all_tib %>% dplyr::anti_join(man_ver_all_tib, by="U")
  
  man_ver_all_tib %>% dplyr::anti_join(man_pqc_all_tib, by="U") %>% 
    dplyr::inner_join(idat_tib, by=c("U"="Address"))
}

# Determine how many other probes can be removed from remainder (non-S4) by
#   comparing ProbeSeqU Sequences::
man_pqc_all_tib %>% dplyr::anti_join(man_S4_PQC_grs_tib, by="U")

man_pqc_all_tib %>% 
  dplyr::anti_join(epic_S4_grs_tib, by=c("AlleleA_Probe_Sequence"="ProbeSeq_A"))


#
# LEFT OFF HERE!!!!
#

# 1.) Extract Sesame hg38 M/U: S4 and format into Genomic Region (grs)

#
# TBD:: Need to pull original and formatted data from source!!!
#  - Add Seq48U to each manifest
#  - Match Seq48U for each manifest against imp_pos3_tib to extract Seq50U
#     this needs to be done based on Infinium I/II
#  - Now remove non-S4 that have the same Seq50U probes
#  - Update all S4 loci with TB/CO and check coordinates
#


man_pqc_new_tib <- man_pqc_all_tib %>% 
  dplyr::anti_join(dplyr::filter(man_ver_all_tib, Version=="S4"), by=c("M","U")) %>%
  dplyr::distinct(U,M, .keep_all=TRUE)

man_pqc_new2_tib <- dplyr::filter(man_ver_all_tib, Version!="S4") %>% 
  dplyr::distinct(U,M, .keep_all=TRUE)

man_pqc_new_tib %>% dplyr::anti_join(man_pqc_new2_tib, by=c("M","U"))
man_pqc_new2_tib %>% dplyr::anti_join(man_pqc_new_tib, by=c("M","U")) %>% 
  dplyr::anti_join(idat_tib, by=c("U"="Address"))


#
# 2.) COVIC Only Analytical Probes::
#
man_cov_dat_tib <- man_pqc_all_tib %>% 
  dplyr::anti_join(epic_b2_tag_tib, by=c("U","M")) %>%
  # dplyr::anti_join(covic_ses_add_tib, by=c("U","M")) %>%
  dplyr::mutate(Platform="EPIC", Probe_Source="COVIC", Version="C4") %>%
  dplyr::rename(AlleleA_ProbeSeq=AlleleA_Probe_Sequence,
                AlleleB_ProbeSeq=AlleleB_Probe_Sequence)

man_cov_dat_tib %>% dplyr::select(dplyr::all_of(fin_epi_cov_ctl_col))

#
# 3.) Add Controls from covic_c0_tib
#
all_ctl_dat_csv <- file.path(par$datDir, 'manifest/core/LEGX-C24.manifest.sesame-base.cpg-sorted.csv.gz')
all_ctl_dat_tib <- readr::read_csv(all_ctl_dat_csv, guess_max = 50000)

man_ctl_dat_tib <- covic_c0_tib %>% 
  dplyr::filter(Probe_Type != "cg") %>%
  dplyr::filter(Probe_Type != "ch") %>%
  dplyr::filter(Probe_Type != "rs")

ctl_cov_dat_tib1 <- all_ctl_dat_tib %>% 
  dplyr::select(-DESIGN,-COLOR_CHANNEL,-col,-Probe_Type,-Probe_Source,-Next_Base,-Version) %>% 
  dplyr::inner_join(man_ctl_dat_tib %>% dplyr::select(-Probe_ID), by=c("M","U"))

ctl_cov_dat_tib <-dplyr::bind_rows(
  ctl_cov_dat_tib1,
  man_ctl_dat_tib %>% dplyr::anti_join(ctl_cov_dat_tib1, by=c("M","U") )
) %>%
  dplyr::mutate(M=as.integer(M), U=as.integer(U),
                Probe_Source="EPIC", Platform="EPIC", Version="B4")

#
# These three need to be joined::
#  - epic_b2_tag_tib,man_cov_dat_tib,man_ctl_dat_tib
#  - ensure we have all original EPIC IDs
#  - ensure all IDs are unique
#  - ensure all tango pairs are unique
#  - cross reference with B2
#  - cut final manifest to get coordinates
#  - write full outputs (TB,CO, Probe_Seq...)
#

fin_epi_cov_ctl_col <- 
  c("Probe_ID","M","U","DESIGN","COLOR_CHANNEL","col","Probe_Type","Next_Base",
    "Probe_Source","Platform","Version",
    "AlleleA_ProbeSeq","AlleleB_ProbeSeq","PRB1_U_MAT","Seq_48U")

fin_epi_cov_ctl_tib <- dplyr::bind_rows(
  epic_b2_tag_tib,man_cov_dat_tib) %>% 
  dplyr::bind_rows(ctl_cov_dat_tib) %>% 
  dplyr::select(dplyr::all_of(fin_epi_cov_ctl_col)) %>%
  dplyr::arrange(Probe_ID)

# Summary::
fin_epi_cov_ctl_tib %>% dplyr::filter(is.na(PRB1_U_MAT)) %>% 
  dplyr::group_by(Probe_Type,Probe_Source,Platform,Version) %>% 
  dplyr::summarise(Count=n(), .groups="drop") %>% 
  dplyr::arrange(-Count)

fin_epi_cov_ctl_tib %>% 
  dplyr::group_by(Probe_Type,Probe_Source,Platform,Version) %>% 
  dplyr::summarise(Count=n(), .groups="drop") %>% 
  dplyr::arrange(-Count)

if (FALSE) {
  # CG# Discrepencies in EPIC::
  #
  # TBD:: Match the discrepency probes::
  #
  man_epi_mis_tib <- man_epi_ses_tib %>% 
    dplyr::filter(Probe_Type_EPI == Probe_Type_PQC) %>% 
    dplyr::filter(Source_ID != Probe_ID_SES) %>% 
    dplyr::select(Probe_ID,Source_ID,Source_Tag,Probe_Type_PQC,Probe_ID_SES,
                  dplyr::everything())
  
  man_pqc_sel_col <- c("Source_ID","Source_Tag", "Probe_ID_V1", "Probe_ID_V2",
                       "Probe_Type_PQC", "Probe_ID_SES",
                       "Seq_ID_PQC","Seq_ID_IMP")
  
  man_cgn_mis_tib <- man_pqc_all_tib %>% 
    dplyr::inner_join(imp_pos3_tib, by=c("PRB1_U_MAT"="Seq50U"), suffix=c("_PQC", "_IMP")) %>% 
    dplyr::select(Probe_ID, Seq_ID_PQC, Seq_ID_IMP) %>% 
    dplyr::inner_join(man_epi_mis_tib, by=c("Seq_ID_PQC"="Seq_ID"), suffix=c("_V1", "_V2")) %>%
    dplyr::select(dplyr::all_of(man_pqc_sel_col))
  
  # 
  # This verifies all the CpG's in EPIC have other names, but all of them are found!!!
  #
  man_epi_mis_tib %>% dplyr::filter(Source_ID %in% man_cgn_mis_tib$Seq_ID_IMP)
  
}

#
# New Version::
#
fin_ses_man_tib <- fin_epi_cov_ctl_tib

#
# Sesame Manifest QC::
#
# This should be zero::
check0_val <- covic_c0_tib %>% 
  dplyr::anti_join(fin_ses_man_tib, by=c("M","U")) %>% 
  dplyr::distinct(M,U, .keep_all=TRUE) %>%
  dplyr::filter(stringr::str_detect(Probe_ID, "_")) %>% 
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ct')) %>%
  base::nrow()
cat(glue::glue("COVIC Validation Checks: 0:{check0_val}{RET}{RET}"))

epi_ses_man_tib <- covic_c0_tib %>% 
  dplyr::anti_join(fin_ses_man_tib, by=c("M","U")) %>%
  dplyr::distinct(M,U, .keep_all=TRUE)

fin_all_man_tib <- dplyr::bind_rows(epi_ses_man_tib,fin_ses_man_tib) %>%
  dplyr::arrange(Probe_ID)

# These should all be the same number::
check1_val <- fin_all_man_tib %>% base::nrow()
check2_val <- fin_all_man_tib %>% dplyr::distinct(M,U) %>% base::nrow()
check3_val <- fin_all_man_tib %>% dplyr::distinct(Probe_ID) %>% base::nrow()
cat(glue::glue("COVIC Validation Checks: 1:{check1_val}, 2:{check2_val}, 3:{check3_val}{RET}{RET}"))

fin_all_man_tib %>% 
  dplyr::group_by(DESIGN,Probe_Type,Probe_Source,Version) %>% 
  dplyr::summarise(Sum_Count=n(), .groups="drop") %>% 
  dplyr::arrange(-Sum_Count) %>%
  print(n=1000)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                  7.1 COVIC Sesame Coordinates Construction::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# man_pqc_all_tib %>% 
#   dplyr::inner_join(imp_pos3_tib, by=c("PRB1_U_MAT"="Seq50U"), suffix=c("_PQC", "_IMP"))
# 
# man_pqc_all_tib %>% 
#   dplyr::inner_join(imp_pos3_tib, by=c("PRB1_U_MAT"="Seq50U"), suffix=c("_PQC", "_IMP")) %>% 
#   dplyr::select(Probe_ID,U_PQC,M_PQC,
#                 AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Chromosome,Coordinate) %>% 
#   dplyr::distinct()

# fin_man_grs_tib <- man_pqc_all_tib %>% 
fin_man_grs_tib <- fin_ses_man_tib %>% 
  dplyr::filter(!is.na(PRB1_U_MAT)) %>%
  dplyr::inner_join(imp_pos3_tib, by=c("PRB1_U_MAT"="Seq50U"), suffix=c("_PQC", "_IMP")) %>% 
  dplyr::select(Probe_ID,Chromosome,Coordinate) %>% 
  dplyr::distinct() %>%
  dplyr::rename(chrom=Chromosome,
                chromStart=Coordinate) %>%
  dplyr::mutate(chrom=stringr::str_remove(chrom,'^chr'),
                chrom=paste0('chr',chrom),
                chromEnd=chromStart+1) %>%
  dplyr::arrange(chrom,chromStart) %>%
  dplyr::select(chrom,chromStart,chromEnd,Probe_ID) %>%
  dplyr::rename(Seq_ID=Probe_ID)
# fin_man_bed_tib <- fin_man_grs_tib %>% dplyr::mutate(Probe_ID=Seq_ID, Seq_ID=stringr::str_remove(Seq_ID, '_.*$'))

fin_man_ana_tib <- NULL
if (!is.null(opt$annDir) && dir.exists(opt$annDir)) {
  fin_man_ana_tib <- 
    manifestToAnnotation(tib=fin_man_grs_tib, 
                         ann=opt$annDir, gen=opt$genomeBuild,
                         verbose=opt$verbose,vt=1,tc=1,tt=pTracker)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 7.2 Write Outputs:: Manifest/BED/Annotation
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Add Genomic Counts to BED file as Score::
#
fin_man_cnt_tib <- fin_man_grs_tib %>% 
  dplyr::add_count(Seq_ID, name="Genomic_Map_Count") %>%
  dplyr::distinct(Seq_ID, Genomic_Map_Count)

fin_man_out_tib <- fin_all_man_tib %>% 
  dplyr::select(-PRB1_U_MAT,-Seq_48U) %>% 
  dplyr::left_join(fin_man_cnt_tib, by=c("Probe_ID"="Seq_ID"))

# Write BED & Annotation Files::
#
fin_all_out_dir <- file.path(par$datDir, 'manifest/covic')
fin_all_man_fns <- paste0( paste('EPIC',opt$version, sep='-'), '.manifest.sesame-base.cpg-sorted.csv.gz')
if (!dir.exists(fin_all_out_dir)) dir.create(fin_all_out_dir, recursive=TRUE)
fin_all_man_csv <- file.path(fin_all_out_dir, fin_all_man_fns)

fin_man_ana_dir <- file.path(par$datDir, 'manifest/annotation/covic')
fin_man_bed_fns <- paste0( paste('EPIC',opt$version, sep='-'), '.bed.gz')
fin_man_ana_fns <- paste0( paste('EPIC',opt$version, sep='-'), '.annotation.csv.gz')

if (!dir.exists(fin_man_ana_dir)) dir.create(fin_man_ana_dir, recursive=TRUE)
fin_man_bed_tsv <- file.path(fin_man_ana_dir,fin_man_bed_fns)
fin_man_ana_csv <- file.path(fin_man_ana_dir,fin_man_ana_fns)

opt$write_data <- TRUE
opt$write_data <- FALSE

# if (opt$write_data) readr::write_csv(fin_all_man_tib, fin_all_man_csv)
# if (opt$write_data) readr::write_csv(fin_ses_man_tib, fin_all_man_csv)
if (opt$write_data) readr::write_csv(fin_man_out_tib, fin_all_man_csv)
if (opt$write_data) readr::write_tsv(fin_man_grs_tib, fin_man_bed_tsv)
if (opt$write_data) readr::write_csv(fin_man_ana_tib, fin_man_ana_csv)






#
#
# Make unique mapping fin_epi_cov_ctl_tib based on Sesame manifest overlap::
#
#   epic_sub_ses_tib
#   fin_epi_cov_ctl_tib



















#
# Something like below::
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  4.0 Rejoin all Manifest Mapped Probes::
#
#  Now for both::
#    Old = base_man_pos_tib
#    New = full_man_pos_tib
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

base_man_grs_tib <- base_man_pos_tib %>% 
  dplyr::rename(Probe_Class=Probe_Type,
                chrom=Gen_Chr,
                chromStart=Gen_Pos) %>%
  dplyr::mutate(chrom=stringr::str_remove(chrom,'^chr'),
                chrom=paste0('chr',chrom),
                chromEnd=chromStart+1) %>%
  dplyr::arrange(chrom,chromStart) %>%
  dplyr::select(chrom,chromStart,chromEnd,Seq_ID)

full_man_grs_tib <- full_man_pos_tib %>%
  dplyr::rename(chrom=Chr_Imp,
                chromStart=Pos_Imp) %>%
  dplyr::mutate(Probe_Class=stringr::str_sub(Seq_ID, 1,2),
                chrom=stringr::str_remove(chrom,'^chr'),
                chrom=paste0('chr',chrom),
                chromEnd=chromStart+1) %>%
  dplyr::arrange(chrom,chromStart) %>%
  dplyr::select(chrom,chromStart,chromEnd,Seq_ID)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                  7.3 COVIC Sesame Annotation Construction::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                5.0.1 Build Annotation GRanges:: Genes/Islands::
#
#  Now for both::
#    Old = base_man_grs_tib
#    New = full_man_grs_tib
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

base_man_ana_tib <- NULL
if (!is.null(opt$annDir) && dir.exists(opt$annDir)) {
  base_man_ana_tib <- 
    manifestToAnnotation(tib=base_man_grs_tib, 
                         ann=opt$annDir, gen=opt$genomeBuild,
                         verbose=opt$verbose,vt=1,tc=1,tt=pTracker)
}

full_man_ana_tib <- NULL
if (!is.null(opt$annDir) && dir.exists(opt$annDir)) {
  full_man_ana_tib <- 
    manifestToAnnotation(tib=full_man_grs_tib, 
                         ann=opt$annDir, gen=opt$genomeBuild,
                         verbose=opt$verbose,vt=1,tc=1,tt=pTracker)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    7.4 COVIC Sesame Validation Step::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# Check Against Previous Results::

# Load to BOX

# Write Email

# Done!!!



# Match Prb1U from Raw Designs to ensure matching for positions::
#   man_pqc_all_tib %>% dplyr::select(Probe_ID, M,U, PRB1_U_MAT) 

#
# 1. Match by PRB1_U_MAT
# 2. Match by both
#

man_pqc_all_tib %>% dplyr::select(Probe_ID,Source_ID,Source_Tag, PRB1_U_MAT, M,U) %>% 
  dplyr::inner_join(imp_pos3_tib, by=c("PRB1_U_MAT"="Seq50U"), suffix=c("_PQC", "_IMP"))


man_pqc_all_tib %>% 
  dplyr::select(Probe_ID,Source_ID,Source_Tag, PRB1_U_MAT, M,U) %>% 
  dplyr::inner_join(imp_pos3_tib, by=c("PRB1_U_MAT"="Seq50U"), suffix=c("_PQC", "_IMP")) %>% 
  dplyr::distinct(PRB1_U_MAT,Chromosome,Coordinate, .keep_all=TRUE) %>%
  
  dplyr::group_by(Source_ID,Seq_ID) %>% dplyr::summarise(Count=n(), .groups='drop') %>% 
  dplyr::arrange(-Count)
  
  
  # dplyr::distinct(Probe_ID,Chromosome,Coordinate, .keep_all=TRUE) %>% 
  dplyr::group_by(Probe_ID) %>% 
  dplyr::summarise(Genomic_CGN_Count=n(), .groups="drop") %>% 
  dplyr::group_by(Genomic_CGN_Count) %>% 
  dplyr::summarise(Genomic_Count_Hist=n(), .groups='drop') %>%
  print(n=1000)









man_pqc_all_tib %>% dplyr::select(Probe_ID, PRB1_U_MAT, M,U) %>% 
  dplyr::inner_join(imp_pos3_tib, by=c("M","U"), suffix=c("_PQC", "_IMP"))


aln_MUs_tib <- man_pqc_all_tib %>% dplyr::select(Probe_ID, PRB1_U_MAT, M,U) %>% 
  dplyr::inner_join(imp_pos3_tib, by=c("PRB1_U_MAT"="Seq50U","M","U"), suffix=c("_PQC", "_IMP"))

aln_MUs_tib %>% dplyr::filter(PRB1_U_MAT==Seq50U)
aln_MUs_tib %>% dplyr::filter(PRB1_U_MAT!=Seq50U)

aln_all_tib <- man_pqc_all_tib %>% dplyr::select(Probe_ID, PRB1_U_MAT, M,U) %>% 
  dplyr::inner_join(imp_pos3_tib, 
                    by=c("PRB1_U_MAT"="Seq50U"), suffix=c("_PQC", "_IMP"))
  # dplyr::filter(Probe_ID != 'cg00001229_TC11') %>% 
  # dplyr::select(Probe_ID,Seq_ID, FR_DB1, TB_DB1, TB_DB2, CO_DB1, CO_DB2)

man_pqc_all_tib %>% dplyr::select(Probe_ID, PRB1_U_MAT, M,U) %>% 
  dplyr::inner_join(imp_pos3_tib, by=c("PRB1_U_MAT"="Seq50U", "M","U"))

#
# NExt Steps::
#   - M/U differences and why
#   - Summary stats on 
#









#
# OLD School Stuff::
#

# OLD Diagnostic Comparison::
#
if (FALSE) {
  
  # Original COVIC Manifest::
  covic_man_csv <- '/Users/bretbarnes/Documents/data/COVIC/transfer/EPIC-C0.manifest.sesame-base.cpg-sorted.csv.gz'
  covic_c0_tib <- readr::read_csv(covic_man_csv) 
  covic_man_tib <- readr::read_csv(covic_man_csv) %>% 
    dplyr::mutate(U=as.integer(U),
                  M=as.integer(M),
                  Design_TypeA=dplyr::case_when(!is.na(M) ~ 'U', TRUE ~ '2' ), 
                  Design_TypeB=dplyr::case_when(!is.na(M) ~ 'M', TRUE ~ '2' ),
                  Probe_Type=stringr::str_sub(Probe_ID, 1,2) ) %>% 
    dplyr::select(Probe_ID, M,U, Design_TypeA, Design_TypeB, Next_Base, 
                  dplyr::everything())
  
  # Seq_ID=stringr::str_remove(Probe_ID, '_.*$'),
  # Probe_Type=stringr::str_sub(Seq_ID, 1,2)
  
  #
  # Previous COVIC Table:: covic_old_man_tib
  #
  covic_old_man_tib <- covic_man_tib %>%
    dplyr::filter(Probe_Type=='cg') %>%
    dplyr::filter( stringr::str_detect(Probe_ID, '_')) %>% 
    dplyr::distinct(M,U, .keep_all=TRUE)
  
  covic_old_man_tab <- dplyr::bind_rows(
    covic_old_man_tib %>% dplyr::select(-M,-Design_TypeB) %>% 
      dplyr::rename(Address=U,Design_Type=Design_TypeA),
    covic_old_man_tib %>% dplyr::select(-U,-Design_TypeA) %>% 
      dplyr::rename(Address=M,Design_Type=Design_TypeB)
  ) %>% 
    dplyr::filter(!is.na(Address) & !is.na(Design_Type)) %>%
    dplyr::arrange(Probe_ID)
  
  #
  # New COVIC Table:: man_pqc_all_tib
  #
  covic_new_man_tib <- man_pqc_all_tib %>%
    dplyr::mutate(U=as.integer(U),
                  M=as.integer(M),
                  Design_TypeA=dplyr::case_when(!is.na(M) ~ 'U', TRUE ~ '2' ), 
                  Design_TypeB=dplyr::case_when(!is.na(M) ~ 'M', TRUE ~ '2' ),
                  Probe_Type=stringr::str_sub(Probe_ID, 1,2) ) %>% 
    dplyr::distinct(M,U, .keep_all=TRUE) %>%
    dplyr::select(Probe_ID, M,U, Design_TypeA, Design_TypeB, Next_Base, 
                  dplyr::everything())
  
  covic_new_man_tab <- dplyr::bind_rows(
    covic_new_man_tib %>% dplyr::select(-M,-Design_TypeB) %>% 
      dplyr::rename(Address=U,Design_Type=Design_TypeA),
    covic_new_man_tib %>% dplyr::select(-U,-Design_TypeA) %>% 
      dplyr::rename(Address=M,Design_Type=Design_TypeB)
  ) %>% 
    dplyr::filter(!is.na(Address) & !is.na(Design_Type)) %>%
    dplyr::arrange(Probe_ID)
  
}

imp_pos3_col <- cols(
  Seq_ID = col_character(),
  Chromosome = col_character(),
  Coordinate = col_integer(),
  Methyl_Allele_FR_Strand = col_character(),
  Methyl_Allele_TB_Strand = col_character(),
  Methyl_Allele_CO_Strand = col_character(),
  Methyl_Next_Base = col_character(),
  UnMethyl_Probe_Sequence = col_character()
)
imp_pos3_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.tsv.gz'
imp_pos3_tib <- readr::read_tsv(imp_pos3_tsv, col_types=imp_pos3_col)

seq48_full_imp_join_tib <- out_pos1_tib %>% 
  dplyr::inner_join(imp_pos3_tib, by=c("CGN_Imp"="Seq_ID", "TB"="Methyl_Allele_TB_Strand", "CO"="Methyl_Allele_CO_Strand"))


# cmd_2 <- "gzip -dc /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/genomic/GRCh38.improbeDesignInput.cgn-sorted.tsv.gz | join -t $'\t' -11 -22 - /Users/bretbarnes/Documents/scratch/build_bead_manifest_simple/GRCh38-COVIC-C4/intersection/cg/GRCh38-COVIC-C4.seq48U_to_cgn.int.cg-sorted.tsv | gzip -c - > /Users/bretbarnes/Documents/scratch/build_bead_manifest_simple/GRCh38-COVIC-C4/intersection/cg/GRCh38-COVIC-C4.seq48U_to_cgn.int.GRCh38.improbeDesignInput.tsv.gz"
imp_pos2_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/genomic/GRCh38.improbeDesignInput.cgn-sorted.tsv.gz'
out_pos2_tsv <- file.path(opt$outDir, 'intersection/cg/GRCh38-COVIC-C4.seq48U_to_cgn.int.GRCh38.improbeDesignInput.tsv.gz')
int_pos2_cmd <- glue::glue("gzip -dc {imp_pos2_tsv} | join -t $'\t' -11 -22 - {out_pos1_tsv} | gzip -c - > {out_pos2_tsv}")
system(int_pos2_cmd)

full_pos_col <- cols(CGN_Imp = col_character(),
                     Chr_Imp = col_character(),
                     Pos_Imp = col_integer(),
                     Seq_48U = col_character(),
                     TB = col_character(),
                     CO = col_character(),
                     M = col_integer(),
                     U = col_integer() )
# full_pos_tsv <- "/Users/bretbarnes/Documents/scratch/build_bead_manifest_simple/GRCh38-COVIC-C4/intersection/cg/GRCh38-COVIC-C4.seq48U_to_cgn.int.GRCh38.improbeDesignInput.tsv.gz"
full_pos_tsv <- out_pos2_tsv
full_pos_tib <- readr::read_tsv(full_pos_tsv, col_names=names(full_pos_col$cols), col_types=full_pos_col)

full_pos_tib %>% dplyr::group_by(M,U) %>% 
  dplyr::summarise(Genomic_CGN_Count=n(), .groups="drop") %>% 
  dplyr::group_by(Genomic_CGN_Count) %>% 
  dplyr::summarise(Genomic_Count_Hist=n(), .groups='drop') %>% print(n=1000)

# First Summary::
full_man_pos_tib <- man_pqc_all_tib %>% 
  dplyr::inner_join(full_pos_tib, by=c("M","U")) %>%
  dplyr::select(Seq_ID,CGN_Imp,Chr_Imp,Pos_Imp) %>%
  dplyr::arrange(Seq_ID) %>%
  dplyr::distinct() %>%
  dplyr::add_count(Seq_ID, name='Genomic_CGN_Count')
full_man_pos_tib %>% 
  dplyr::group_by(Genomic_CGN_Count) %>% 
  dplyr::summarise(Genomic_Count_Hist=n(), .groups='drop') %>% print(n=1000)

# Validation of extra mappings::
full_mat_pos_tib <- full_man_pos_tib %>%
  dplyr::filter(Seq_ID==CGN_Imp) %>%
  dplyr::distinct() %>%
  dplyr::add_count(Seq_ID, name='Genomic_CGN_Count')
full_mat_pos_tib %>% 
  dplyr::group_by(Genomic_CGN_Count) %>% 
  dplyr::summarise(Genomic_Count_Hist=n(), .groups='drop') %>% print(n=1000)

# Original perfect CGN matching::
base_man_pos_tib <- dplyr::bind_rows( man_pos_list, .id = "Probe_Type" )
base_man_pos_tib %>%
  dplyr::group_by(Genomic_CGN_Count) %>% 
  dplyr::summarise(Genomic_Count_Hist=n(), .groups='drop') %>% print(n=1000)


base_vs_full_cnt_tib <- dplyr::inner_join(
  base_man_pos_tib %>% dplyr::select(Seq_ID,Genomic_CGN_Count) %>% dplyr::distinct(),
  full_man_pos_tib %>% dplyr::select(Seq_ID,Genomic_CGN_Count) %>% dplyr::distinct(),
  by="Seq_ID", suffix=c("_base", "_core")
)

ggplot2::ggplot(data=base_vs_full_cnt_tib, aes(x=Genomic_CGN_Count_base, y=Genomic_CGN_Count_core)) +
  ggplot2::geom_point()

#
# TBD:: NEXT:: investigate the biggest differences...
#
base_vs_full_dif_tib <- base_vs_full_cnt_tib %>% 
  dplyr::mutate(CGN_Count_Delta=Genomic_CGN_Count_core-Genomic_CGN_Count_base) %>% 
  dplyr::arrange(-CGN_Count_Delta)

full_man_pos_tib %>% dplyr::filter(Seq_ID=='cg00001229')
full_mat_pos_tib %>% dplyr::filter(Seq_ID=='cg00001229')

# Example of multiple mappings
# gzip -dc data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.tsv.gz | head -n 100000 |  grep -e "^cg00001229" -e "^cg00003967" -e "^cg00004386" -e "^cg00012185" -e "^cg00017083" | cut -f 1,16,21-23 | sort -k 4,5

# TBD NEXT::
#   - Load 13k most important
#   - Get Genes for base and full

# #################################### # ####################################
#
#
#                  TANGO ADDRESS SANITY COUNT CHECK::
#
#
# #################################### # ####################################

top_cgn_csv <- '/Users/bretbarnes/Documents/data/COVIC/transfer/COVID-Controls_EWASsignificantresults_119FIXED.csv'
top_cgn_tib <- readr::read_csv(top_cgn_csv) %>% 
  dplyr::mutate(Seq_ID=stringr::str_remove(probeID, '_.*$')) %>% 
  dplyr::select(Seq_ID, probeID, dplyr::everything())

covic_top_cgn_tib <- top_cgn_tib %>% 
  dplyr::filter(stringr::str_detect(probeID, '_'))

covic_top_cgn_tib %>% dplyr::select(Seq_ID) %>% dplyr::left_join(base_vs_full_dif_tib, by="Seq_ID")

# New Names:: Old Names
# 
# cg25888371_TC21 == cg25888638
# cg25888638
#

# INCONSISTENCIES::
#

man_pqc_all_tib %>% dplyr::anti_join(covic_old_man_tib, by=c("M","U"))
covic_old_man_tib %>% dplyr::anti_join(man_pqc_all_tib, by=c("M","U"))

covic_old_man_tab %>% dplyr::anti_join(covic_new_man_tab, by="Address")
covic_new_man_tab %>% dplyr::anti_join(covic_old_man_tab, by="Address")

covic_mat_aqp_pas_tab %>% dplyr::anti_join(covic_new_man_tab, by="Address")
covic_mat_aqp_pas_tab %>% dplyr::anti_join(covic_old_man_tab, by="Address")

covic_old_man_tab %>% dplyr::anti_join(covic_mat_aqp_pas_tab, by="Address")
covic_new_man_tab %>% dplyr::anti_join(covic_mat_aqp_pas_tab, by="Address")





# #################################### # ####################################
#
#
#          TANGO ADDRESS SANITY COUNT CHECK 2:: Get new/old mapppings::
#
#
# #################################### # ####################################

top_cgn_tib %>% dplyr::select(Seq_ID,probeID) %>%
  dplyr::filter( stringr::str_detect(probeID, '_')) %>%
  dplyr::inner_join(covic_only_man_tib, by="Seq_ID")
  
# Need to break down by tangos and then get old/new name mappings...
# use:: man_pqc_all_tib

covic_only_add_tib <- dplyr::bind_rows(
  covic_only_man_tib %>% dplyr::filter(!is.na(U) & !is.na(M)) %>% 
    dplyr::select(-U) %>%
    dplyr::rename(Address=M) %>%
    # dplyr::select(Probe_ID,Address) %>% 
    dplyr::mutate(Infinium_Type='M'),
  
  covic_only_man_tib %>% dplyr::filter(!is.na(U) & !is.na(M)) %>% 
    dplyr::select(-M) %>%
    dplyr::rename(Address=U) %>%
    # dplyr::select(Probe_ID,Address) %>% 
    dplyr::mutate(Infinium_Type='U'),
  
  covic_only_man_tib %>% dplyr::filter(!is.na(U) &  is.na(M)) %>%  
    dplyr::select(-M) %>%
    dplyr::rename(Address=U) %>%
    # dplyr::select(Probe_ID,Address) %>% 
    dplyr::mutate(Infinium_Type='2')
  
) %>% 
  dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2),
                Seq_ID=stringr::str_remove(Probe_ID, '_.*$')) %>%
  dplyr::select(Seq_ID,Probe_ID,dplyr::everything()) %>%
  dplyr::arrange(Seq_ID, Probe_ID) %>%
  dplyr::filter(!is.na(Address))

# Join with Top List COVIC Only::
# covic_cur_map1_tib <- covic_top_cgn_tib %>% dplyr::select(Seq_ID,probeID) %>%
#   dplyr::inner_join( covic_only_add_tib %>% dplyr::select(Seq_ID,Probe_ID,Address,Next_Base,Infinium_Type),
#                      by="Seq_ID")

covic_top_map_tib <- covic_top_cgn_tib %>% dplyr::select(Seq_ID,probeID) %>%
  dplyr::inner_join( covic_only_add_tib %>% dplyr::select(Seq_ID,Probe_ID,Address,Next_Base,Infinium_Type),
                     by=c("Seq_ID","probeID"="Probe_ID") ) %>%
  dplyr::rename(Seq_ID_Old=Seq_ID, Probe_ID_Old=probeID) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Seq_ID_Old,Probe_ID_Old)


covic_cur_map_tib <- covic_only_add_tib %>% 
  dplyr::select(Seq_ID,Probe_ID,Address,Next_Base,Infinium_Type) %>%
  dplyr::rename(Seq_ID_Old=Seq_ID, Probe_ID_Old=Probe_ID) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Seq_ID_Old,Probe_ID_Old)

old_new_mapM_tib <- covic_cur_map_tib %>% dplyr::inner_join(
  man_pqc_all_tib %>% dplyr::select(Seq_ID, Probe_ID, M) %>% dplyr::filter(!is.na(M)) %>% dplyr::rename(Address=M),
  by="Address")

old_new_mapU_tib <- covic_cur_map_tib %>% dplyr::inner_join(
  man_pqc_all_tib %>% dplyr::select(Seq_ID, Probe_ID, U) %>% dplyr::filter(!is.na(U)) %>% dplyr::rename(Address=U),
  by="Address")

# Now Add Genomic Counts Delta and then Genomic Coordinates 

old_new_mapU_cnt_tib <- old_new_mapU_tib %>% 
  dplyr::inner_join(base_vs_full_dif_tib, by="Seq_ID") %>% 
  dplyr::arrange(-CGN_Count_Delta)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  4.0 Rejoin all Manifest Mapped Probes::
#
#  Now for both::
#    Old = base_man_pos_tib
#    New = full_man_pos_tib
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

base_man_grs_tib <- base_man_pos_tib %>% 
  dplyr::rename(Probe_Class=Probe_Type,
                chrom=Gen_Chr,
                chromStart=Gen_Pos) %>%
  dplyr::mutate(chrom=stringr::str_remove(chrom,'^chr'),
                chrom=paste0('chr',chrom),
                chromEnd=chromStart+1) %>%
  dplyr::arrange(chrom,chromStart) %>%
  dplyr::select(chrom,chromStart,chromEnd,Seq_ID)

full_man_grs_tib <- full_man_pos_tib %>%
  dplyr::rename(chrom=Chr_Imp,
                chromStart=Pos_Imp) %>%
  dplyr::mutate(Probe_Class=stringr::str_sub(Seq_ID, 1,2),
                chrom=stringr::str_remove(chrom,'^chr'),
                chrom=paste0('chr',chrom),
                chromEnd=chromStart+1) %>%
  dplyr::arrange(chrom,chromStart) %>%
  dplyr::select(chrom,chromStart,chromEnd,Seq_ID)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                5.0.1 Build Annotation GRanges:: Genes/Islands::
#
#  Now for both::
#    Old = base_man_grs_tib
#    New = full_man_grs_tib
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

base_man_ana_tib <- NULL
if (!is.null(opt$annDir) && dir.exists(opt$annDir)) {
  base_man_ana_tib <- 
    manifestToAnnotation(tib=base_man_grs_tib, 
                         ann=opt$annDir, gen=opt$genomeBuild,
                         verbose=opt$verbose,vt=1,tc=1,tt=pTracker)
}

full_man_ana_tib <- NULL
if (!is.null(opt$annDir) && dir.exists(opt$annDir)) {
  full_man_ana_tib <- 
    manifestToAnnotation(tib=full_man_grs_tib, 
                         ann=opt$annDir, gen=opt$genomeBuild,
                         verbose=opt$verbose,vt=1,tc=1,tt=pTracker)
}

#
# Check Genes base vs full for top COVIC probes::
#
base_old_new_mapU_cnt_tib <- old_new_mapU_cnt_tib %>%
  dplyr::inner_join(base_man_ana_tib, by="Seq_ID")

base_gene_sum_tib <- base_old_new_mapU_cnt_tib %>% 
  # dplyr::filter(Source=='UCSC') %>% 
  dplyr::filter(Source=='NCBI') %>% 
  dplyr::distinct(Seq_ID,Gene,Feature) %>% 
  dplyr::group_by(Seq_ID,Gene,Feature) %>% 
  dplyr::summarise(Feature_Count=n(), .groups='drop') %>% 
  dplyr::add_count(Seq_ID, name="Loci_Count")

base_gene_sum_tib %>%
  dplyr::arrange(Loci_Count, Seq_ID) %>% 
  print(n=base::nrow(base_gene_sum_tib))

# Full::
full_old_new_mapU_cnt_tib <- old_new_mapU_cnt_tib %>%
  dplyr::inner_join(full_man_ana_tib, by="Seq_ID")

full_gene_sum_tib <- full_old_new_mapU_cnt_tib %>% 
  # dplyr::filter(Source=='UCSC') %>% 
  dplyr::filter(Source=='NCBI') %>% 
  dplyr::distinct(Seq_ID,Gene,Feature) %>% 
  dplyr::group_by(Seq_ID,Gene,Feature) %>% 
  dplyr::summarise(Feature_Count=n(), .groups='drop') %>% 
  dplyr::add_count(Seq_ID, name="Loci_Count")

full_gene_sum_tib %>%
  dplyr::arrange(Loci_Count, Seq_ID) %>% 
  print(n=base::nrow(full_gene_sum_tib))

full_gene_sum_tib %>% dplyr::filter(Seq_ID=='cg00001229')

# Join Summaries::
dplyr::inner_join(
  base_gene_sum_tib %>% dplyr::distinct(Seq_ID,Loci_Count),
  full_gene_sum_tib %>% dplyr::distinct(Seq_ID,Loci_Count),
  by="Seq_ID", suffix=c("_base","_full")
) %>%
  dplyr::mutate(Delta=Loci_Count_full-Loci_Count_base) %>%
  dplyr::arrange(-Delta) %>%
  print(n=100)

full_gene_sum_tib %>% 
  dplyr::distinct(Seq_ID,Gene) %>% 
  dplyr::add_count(Seq_ID, name="Gene_Count") %>% 
  dplyr::distinct(Seq_ID,Gene_Count) %>% 
  dplyr::arrange(-Gene_Count) %>% 
  dplyr::filter(Gene_Count>10) %>% print(n=1000)

#
# TBD Next:: Compare with previous results from Iain...
#








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
  # dplyr::group_by(Seq_ID) %>% 
  dplyr::arrange(chrom,chromStart) %>%
  dplyr::select(chrom,chromStart,chromEnd,Seq_ID) # %>% dplyr::ungroup()

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

ctl_csv <- '/Users/bretbarnes/Documents/tools/dat/manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz'
std_ctl_seq_tsv <- '/Users/bretbarnes/Documents/tools/dat/manifest/controls/01152015_DarkMatterControls.probe.match.tsv.gz'

hsa_ctl_tib <- NULL
if (!is.null(ctl_csv) && file.exists(ctl_csv) &&
    !is.null(std_ctl_seq_tsv) && file.exists(std_ctl_seq_tsv)) {
  
  hsa_ctl_tib <- format_controls_HSA(
    file1=std_ctl_seq_tsv, file2=ctl_csv,
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
#     tib=man_raw_dat_tib, verbose=opt$verbose,vt=0,tc=1,tt=pTracker)
# }
# 
# new_ctl_tib2 <- NULL
# if (par$local_runType=='GRCm38') {
#   new_ctl_tib2 <- format_controls_Unk(
#     tib=man_pqc_unk_tib, verbose=opt$verbose,vt=0,tc=1,tt=pTracker)
# }

man_prd_all_tib <- dplyr::bind_rows(man_pqc_all_tib,hsa_ctl_tib)
man_prd_sum_tib <- man_prd_all_tib %>% 
  dplyr::group_by(Probe_Source,Assay_Class,Probe_Class,Probe_Type,Infinium_Design,AQP) %>% 
  dplyr::summarise(Count=n(), .groups='drop')
print(man_prd_sum_tib,n=base::nrow(man_prd_sum_tib))

#
# QC Checks::
#
man_prd_all_tib %>% dplyr::group_by(DESIGN,COLOR_CHANNEL,col,Next_Base) %>% 
  dplyr::summarise(Count=n(), .groups='drop') %>% as.data.frame()

man_prd_all_tib %>% 
  dplyr::add_count(Probe_ID, name="PID_CNT") %>% 
  dplyr::filter(PID_CNT != 1) %>% 
  dplyr::select(Probe_ID,U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,PID_CNT) %>% 
  dplyr::arrange(Probe_ID)

#
# Write Sesame Output Full::
#
readr::write_csv(man_prd_all_tib,ses_base_csv)


# man_pqc_all_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "ukr") )
# man_prd_all_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "ukr") )
pre_man_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "ukr") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#          Load Previous dat/manifest and compare based on M/U::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

compare_prev_bool <- FALSE
if (compare_prev_bool) {
  
  unq_jon_cols <- c("M","U","AlleleA_Probe_Sequence","AlleleB_Probe_Sequence",
                    "Infinium_Design","AQP")
  
  # Rename derived manifest:: i.e. man_prd_all_tib
  new_man_tib <- man_prd_all_tib
  
  # Load previous manifest::
  pre_man_csv <- file.path(par$datDir,"manifest/core/LEGX-C20.manifest.sesame-base.cpg-sorted.csv.gz")
  pre_man_tib <- suppressMessages(suppressWarnings( readr::read_csv(pre_man_csv) ))
  
  # Known Failed CpH Designs from previous builds::
  pre_mis_tib <- pre_man_tib %>% dplyr::filter(Probe_Class=='ch',Probe_Class=='ch')
  
  #
  # Simple Comparison::
  #
  new_man_tib %>% dplyr::select(dplyr::all_of(base_sel_cols))
  pre_man_tib %>% dplyr::select(dplyr::all_of(base_sel_cols))
  
  min_jon_cols <- c("M","U","AlleleA_Probe_Sequence","AlleleB_Probe_Sequence",
                    "Seq_48U","AQP")
  
  min_jon_tib <- dplyr::inner_join(
    new_man_tib %>% dplyr::select(dplyr::all_of(base_sel_cols)),
    pre_man_tib %>% dplyr::select(dplyr::all_of(base_sel_cols)),
    by=dplyr::all_of(min_jon_cols), suffix=c("_new","_pre")
  )
  
  # min_jon_tib %>% dplyr::filter(Next_Base_new != Next_Base_pre) %>% dplyr::select(starts_with("Probe_ID")) %>% as.data.frame()
  # min_jon_tib %>% dplyr::filter(Probe_ID_new != Probe_ID_pre) %>% dplyr::group_by(Probe_Type_new,Probe_Type_pre) %>% dplyr::summarise(Count=n())
  
  min_sum_tib <- min_jon_tib %>% 
    dplyr::filter(Probe_ID_new != Probe_ID_pre) %>% 
    dplyr::select(Probe_ID_new,Probe_ID_pre,Probe_Type_new,Probe_Type_pre) %>%
    split(.$Probe_Type_new)

  #
  # Concerned about the 100 - 2,465 Probe_Type_new not being cg
  #
  min_sum_tib$ch %>% dplyr::group_by(Probe_Type_new,Probe_Type_pre) %>%
    dplyr::summarise(Count=n())
  
  
  
  #
  #
  # TBD: LEFT OFF HERE::
  #
  #
  
  
  
  
  #
  # Detailed Comparison::
  #
  
  jon_man_tib <- dplyr::inner_join(
    new_man_tib,pre_man_tib, 
    by=dplyr::all_of(unq_jon_cols), suffix=c("_new","_pre")
  )
  
  pre_man_cnt <- pre_man_tib %>% base::nrow()
  new_man_cnt <- new_man_tib %>% base::nrow()
  jon_man_cnt <- jon_man_tib %>% base::nrow()
  
  if (opt$verbose>=1) {
    cat(glue::glue("[{par$prgmTag}]: These should be equal...{RET}"))
    cat(glue::glue("[{par$prgmTag}]: pre_man_cnt={pre_man_cnt}.{RET}"))
    cat(glue::glue("[{par$prgmTag}]: new_man_cnt={new_man_cnt}.{RET}"))
    cat(glue::glue("[{par$prgmTag}]: jon_man_cnt={jon_man_cnt}.{RET}"))
    cat(glue::glue("[{par$prgmTag}]:{RET}"))
  }

  #
  # These should be zero::
  #
  new_mis_cnt <- new_man_tib %>% 
    dplyr::anti_join(pre_man_tib, by=dplyr::all_of(unq_jon_cols)) %>%
    base::nrow()
  pre_mis_cnt <- pre_man_tib %>% 
    dplyr::anti_join(new_man_tib, by=dplyr::all_of(unq_jon_cols)) %>%
    base::nrow()
  
  if (opt$verbose>=1) {
    cat(glue::glue("[{par$prgmTag}]: These should be zero...{RET}"))
    cat(glue::glue("[{par$prgmTag}]: pre_mis_cnt={pre_mis_cnt}.{RET}"))
    cat(glue::glue("[{par$prgmTag}]: new_mis_cnt={new_mis_cnt}.{RET}"))
    cat(glue::glue("[{par$prgmTag}]:{RET}"))
  }

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Load Validation Categories:: CpG/CpH
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  val_dat_col <- c("Chrom","Beg","End","Align_Base","Aln_ID","Align_Score",
                   "Align_Chrom","Align_Beg","Align_Score2","Match_Desc",
                   "Probe_Seq","NM","AS","YD","CN")
  
  #
  # Load:: CpG Count=100
  #
  val_cpg_tsv <- '/Users/bretbarnes/Documents/data/CustomContent/EWAS/from_Wanding_CpH_Lists/true_cpgs.txt'
  val_cpg_tib <- 
    suppressMessages(suppressWarnings( readr::read_tsv(val_cpg_tsv, col_names=val_dat_col) )) %>% 
    dplyr::mutate(Probe_RCS=revCmp(Probe_Seq))
  
  #
  # Load:: CpH Count=2310
  #
  val_cph_tsv <- '/Users/bretbarnes/Documents/data/CustomContent/EWAS/from_Wanding_CpH_Lists/true_cphs.txt'
  val_cph_tib <- 
    suppressMessages(suppressWarnings( readr::read_tsv(val_cph_tsv, col_names=val_dat_col) )) %>% 
    dplyr::mutate(Probe_RCS=revCmp(Probe_Seq))
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #             Summary Validation Comparison:: New/Pre vs CpG/CpH
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  
  #
  # Compare Validated Groups to Previous::
  #
  pre_mat_tib <- dplyr::bind_rows(
    pre_mis_tib %>% dplyr::inner_join(val_cpg_tib, by=c("AlleleA_Probe_Sequence"="Probe_Seq")),
    pre_mis_tib %>% dplyr::inner_join(val_cph_tib, by=c("AlleleA_Probe_Sequence"="Probe_Seq")),
    
    pre_mis_tib %>% dplyr::inner_join(val_cpg_tib, by=c("AlleleB_Probe_Sequence"="Probe_Seq")),
    pre_mis_tib %>% dplyr::inner_join(val_cph_tib, by=c("AlleleB_Probe_Sequence"="Probe_Seq")),
    
    pre_mis_tib %>% dplyr::inner_join(val_cpg_tib, by=c("AlleleA_Probe_Sequence"="Probe_RCS")),
    pre_mis_tib %>% dplyr::inner_join(val_cph_tib, by=c("AlleleA_Probe_Sequence"="Probe_RCS")),
    
    pre_mis_tib %>% dplyr::inner_join(val_cpg_tib, by=c("AlleleB_Probe_Sequence"="Probe_RCS")),
    pre_mis_tib %>% dplyr::inner_join(val_cph_tib, by=c("AlleleB_Probe_Sequence"="Probe_RCS")),
    NULL
  ) %>% 
    dplyr::distinct(M,U, .keep_all=TRUE) %>%
    dplyr::select(Probe_ID,M,U,Aln_ID) %>% 
    dplyr::filter(Probe_ID != Aln_ID)
  
  #
  # Compare Validated Groups to New::
  #
  new_mat_tib <- dplyr::bind_rows(
    new_man_tib %>% dplyr::inner_join(val_cpg_tib, by=c("AlleleA_Probe_Sequence"="Probe_Seq")),
    new_man_tib %>% dplyr::inner_join(val_cph_tib, by=c("AlleleA_Probe_Sequence"="Probe_Seq")),
    
    new_man_tib %>% dplyr::inner_join(val_cpg_tib, by=c("AlleleB_Probe_Sequence"="Probe_Seq")),
    new_man_tib %>% dplyr::inner_join(val_cph_tib, by=c("AlleleB_Probe_Sequence"="Probe_Seq")),
    
    new_man_tib %>% dplyr::inner_join(val_cpg_tib, by=c("AlleleA_Probe_Sequence"="Probe_RCS")),
    new_man_tib %>% dplyr::inner_join(val_cph_tib, by=c("AlleleA_Probe_Sequence"="Probe_RCS")),
    
    new_man_tib %>% dplyr::inner_join(val_cpg_tib, by=c("AlleleB_Probe_Sequence"="Probe_RCS")),
    new_man_tib %>% dplyr::inner_join(val_cph_tib, by=c("AlleleB_Probe_Sequence"="Probe_RCS")),
    NULL
  ) %>% 
    # dplyr::add_count(M,U, name="MU_CNT") %>% dplyr::filter(MU_CNT!=1) %>%
    dplyr::distinct(M,U, .keep_all=TRUE) %>%
    dplyr::select(Probe_ID,M,U,Aln_ID) %>% 
    dplyr::filter(Probe_ID != Aln_ID)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #             Compare Validation Categories to New/Pre:: CpG
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  new_man_tib %>% 
    dplyr::inner_join(val_cpg_tib, by=c("AlleleA_Probe_Sequence"="Probe_RCS")) %>% 
    dplyr::group_by(Probe_Type,Probe_Class,Infinium_Design,AQP) %>% 
    dplyr::summarise(Count=n(), .groups="drop")

  pre_man_tib %>% 
    dplyr::inner_join(val_cpg_tib, by=c("AlleleA_Probe_Sequence"="Probe_RCS")) %>% 
    dplyr::group_by(Probe_Type,Probe_Class,Infinium_Design,AQP) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  jon_man_tib %>% 
    dplyr::inner_join(val_cpg_tib, by=c("AlleleA_Probe_Sequence"="Probe_RCS")) %>% 
    dplyr::group_by(Probe_Type_new,Probe_Class_new,
                    Probe_Type_pre,Probe_Class_pre,
                    Infinium_Design,AQP) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #             Compare Validation Categories to New/Pre:: CpH
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  new_man_tib %>% 
    dplyr::inner_join(val_cph_tib, by=c("AlleleA_Probe_Sequence"="Probe_RCS")) %>% 
    dplyr::group_by(Probe_Type,Probe_Class,Infinium_Design,AQP) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  pre_man_tib %>% 
    dplyr::inner_join(val_cph_tib, by=c("AlleleA_Probe_Sequence"="Probe_RCS")) %>% 
    dplyr::group_by(Probe_Type,Probe_Class,Infinium_Design,AQP) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  jon_man_tib %>% 
    dplyr::inner_join(val_cph_tib, by=c("AlleleA_Probe_Sequence"="Probe_RCS")) %>% 
    dplyr::group_by(Probe_Type_new,Probe_Class_new,
                    Probe_Type_pre,Probe_Class_pre,
                    Infinium_Design,AQP) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  
  #
  # OLD Comparison Code::
  #
  is_relevant <- FALSE
  if (is_relevant) {
    #
    # New/Old Summary Comparison for CpH::
    #
    new_man_tib %>% dplyr::filter(Probe_Type=='ch') %>% 
      dplyr::group_by(Probe_Class,Strand_TB,Strand_CO,Infinium_Design,AQP) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    
    pre_man_tib %>% dplyr::filter(Probe_Type=='ch') %>% 
      dplyr::group_by(Probe_Class,Strand_TB,Strand_CO,Infinium_Design,AQP) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
  }
  
  # TBD:: Compare Changes for each column...
  #
  
}




#
# Write minimal Sesame::
#


#
#
#
#  END HERE...
#
#
#




if (FALSE) {

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

}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
