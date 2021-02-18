
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
  
  man_raw_cor_tib <- 
    loadCoreManifest_COVIC(
      datDir=par$datDir, manDir=par$manDir,
      verbose=opt$verbose,tc=0,tt=pTracker)

  man_raw_dat_tib <- dplyr::bind_rows(
    # Previous EPIC Designs::
    # loadCoreManifest_COVIC(
    #   datDir=par$datDir, manDir=par$manDir,
    #   verbose=opt$verbose,tc=0,tt=pTracker),
    
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
# man_raw_dat_tib <- man_raw_dat_tib %>% dplyr::rename(Probe_ID_SRC=Probe_ID)

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

man_pqc_prb_tib <- bindProbeDesignList(
  list=man_prb_list, # platform=opt$platform, version=opt$version,
  sumDir=opt$sumDir,del='.',
  verbose=opt$verbose,vt=3,tc=1,tt=pTracker) %>% 
  dplyr::mutate(Assay_Class='Analytical')

core_prb_type_vec <- man_pqc_prb_tib %>% 
  dplyr::distinct(Probe_Type) %>% dplyr::pull(Probe_Type)

# Qucik Summary::
man_pqc_prb_tib %>% 
  dplyr::group_by(Manifest,Probe_Type, Infinium_Design) %>% 
  dplyr::summarise(Count=n(), .groups="drop") %>% print()

# Clear up memory::
rm(man_prb_list)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#            Extract Missing Probes from Manifest:: i.e. CTL/UNK
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_pqc_unk_tib <- NULL
man_pqc_all_tib <- man_pqc_prb_tib %>% 
  tidyr::separate(Probe_ID, into=c("Source_ID","Source_Tag"), 
                  sep="_", remove=FALSE) %>% 
  dplyr::filter(!is.na(PRB1_U_MAT))

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
  dplyr::group_by(Manifest,Assay_Class,Probe_Class,Probe_Type,Infinium_Design,AQP) %>% 
  dplyr::summarise(Count=n(), .groups='drop')
print(man_all_sum_tib,n=base::nrow(man_all_sum_tib))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#            4.0 Load all Coordinates for Detected Manifest Probes::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cgn_pos_pre <- paste(par$local_runType,opt$version, sep='-')
cgn_int_dir <- file.path(opt$intDir,'cg')
cgn_out_dir <- file.path(opt$genDir,'cg')
cgn_pos_tsv <- file.path(cgn_int_dir, paste(cgn_pos_pre,'seq48U_to_cgn.int.Seq_48U-sorted.tsv.gz', sep='.'))
out_pos_tsv <- file.path(cgn_out_dir, paste(cgn_pos_pre,'seq48U_to_cgn.int.cg-sorted.tsv', sep='.'))

if (!dir.exists(cgn_int_dir)) dir.create(cgn_int_dir, recursive=TRUE)
if (!dir.exists(cgn_out_dir)) dir.create(cgn_out_dir, recursive=TRUE)

if (!file.exists(out_pos_tsv)) {
  # cmd_1 <- "gzip -dc /Users/bretbarnes/Documents/scratch/build_bead_manifest_simple/GRCh38-COVIC-C4/intersection/cg/GRCh38-COVIC-C4.seq48U_to_cgn.int.Seq_48U-sorted.tsv.gz | sort -k 2,2 > /Users/bretbarnes/Documents/scratch/build_bead_manifest_simple/GRCh38-COVIC-C4/intersection/cg/GRCh38-COVIC-C4.seq48U_to_cgn.int.cg-sorted.tsv"
  int_pos_cmd <- glue::glue("gzip -dc {cgn_pos_tsv} | sort -k 2,2 > {out_pos_tsv}")
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Running; cmd={int_pos_cmd}...{RET}"))
  system(int_pos_cmd)
}

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Loading Intersection; out_pos_tsv={out_pos_tsv}...{RET}"))

out_pos_col <- cols(
  Seq_48U = col_character(),
  CGN_Imp = col_character(),
  TB = col_character(),
  CO = col_character(),
  U = col_integer(),
  M = col_integer()
)
out_pos_tib <- readr::read_tsv(out_pos_tsv, col_names=names(out_pos_col$cols), col_types=out_pos_col)

man_pqc_grs_tibs <- NULL
gbs <- c("GRCh37", "GRCh38")
for (gb in gbs) {
  
  #
  # TBD:: Need to remove man_pqc_all_tib from the function below and replace with a subsetted argument!!!
  #  man_pqc_all_tib %>% dplyr::select(Source_ID,Source_Tag, U,M, AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, Probe_Type,Strand_TB,Strand_CO,Next_Base )
  #
  #  Probe_ID,
  #
  # opt$impDir <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020"
  
  # We can already reduce down to distinct U/M::
  # man_pqc_all_tib %>% dplyr::select(Probe_ID,Source_ID,Source_Tag,U,M,Next_Base,col,Strand_TB,Strand_CO,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence) %>% dplyr::distinct(U,M)
  man_pqc_grs_tibs[[gb]] <- 
    intersectCgnMap_COVIC(tsv=out_pos_tsv, man=man_pqc_all_tib, src=man_pqc_all_tib, build=gb, addBuild=FALSE,
                          datDir="/Users/bretbarnes/Documents/data/improbe/designOutput_21092020", 
                          inpDir=cgn_out_dir, runType=par$local_runType, ver=opt$version,
                          verbose=opt$verbose,tc=1,tt=pTracker)
  # dplyr::mutate(Genome_Build=gb)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           END OF CORE RUN::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


if (FALSE) {

  all_prb_bsp_col <- cols(
    # bsp_key = col_character()
    prb_seq = col_character(),
    prb_inf = col_integer(),
    
    bsp_seq = col_character(),
    bsp_tag = col_character(),
    
    # bsp_chr = col_character()
    bsp_cgn = col_character(),
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
  
  all_prb_bsp_tsv <- "/Users/bretbarnes/Documents/data/COVIC/bsmap/probe/ALL/covic-probes.ALL.s12.v5.g0.p16.n1.r2.Ru.bsp.gz"
  all_prb_bsp_tib <- readr::read_tsv(all_prb_bsp_tsv, 
                                     col_names=names(all_prb_bsp_col$cols), 
                                     col_types=all_prb_bsp_col)

  #
  # Run Mini Test::
  #
  opt$fasDir <- file.path("/Users/bretbarnes/Documents/tmp", "bsmap")
  if (!dir.exists(opt$fasDir)) dir.create(opt$fasDir, recursive=TRUE)
  unlink(list.files(opt$fasDir, full.names = TRUE))

  bsmap_exe <- "/Users/bretbarnes/Documents/programs/BSMAPz/bsmapz"
  test_fas  <- file.path(opt$fasDir, paste("test-covic","prb50U.fa.gz", sep="_"))
  test_csv  <- file.path(opt$fasDir, paste("test-covic","prb50U.csv.gz", sep="_"))
  
  test_ssh  <- file.path(opt$fasDir, paste("test-covic","prb50U.104bp.sh", sep="_"))
  gn37_ssh  <- file.path(opt$fasDir, paste("test-covic","prb50U.gen37.sh", sep="_"))
  gn38_ssh  <- file.path(opt$fasDir, paste("test-covic","prb50U.gen38.sh", sep="_"))
  
  test_bsp  <- file.path(opt$fasDir, paste("test-covic","prb50U.104bp.bsp", sep="_"))
  gn37_bsp  <- file.path(opt$fasDir, paste("test-covic","prb50U.gen37.bsp", sep="_"))
  gn38_bsp  <- file.path(opt$fasDir, paste("test-covic","prb50U.gen38.bsp", sep="_"))
  
  test_bsg  <- file.path(opt$fasDir, paste("test-covic","prb50U.104bp.formatted.bsp.gz", sep="_"))
  gn37_bsg  <- file.path(opt$fasDir, paste("test-covic","prb50U.gen37.formatted.bsp.gz", sep="_"))
  gn38_bsg  <- file.path(opt$fasDir, paste("test-covic","prb50U.gen38.formatted.bsp.gz", sep="_"))
  
  test_ref  <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/cgnTop/GRCh36-GRCh38-GRCm10-21092020.cgnTop.sorted.104bp.fa.gz"
  gn37_ref  <- "/Users/bretbarnes/Documents/data/iGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/WholeGenomeFasta/GRCh37.genome.fa.gz"
  gn38_ref  <- "/Users/bretbarnes/Documents/data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/GRCh38.genome.fa.gz"
  
  #
  # Current Code Below::
  #
  if (FALSE) {
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
    
    #
    # Previous Place to Start::
    #
    
    # Write all probes to fasta::
    fas_raw_all_csv <- "/Users/bretbarnes/Documents/data/COVIC/fas/covic-prbU-all.csv.gz"
    fas_raw_all_fas <- "/Users/bretbarnes/Documents/data/COVIC/fas/covic-prbU-all.fa.gz"
    
    write_fas <- FALSE
    if (write_fas) {
      fas_raw_all_vec <- fas_raw_dat_tib %>%
        dplyr::mutate(
          line=paste0('>',FAS_KEY,'\n',FAS_SEQ)
        ) %>%
        dplyr::filter(!is.na(line)) %>% 
        # dplyr::distinct(line) %>%
        dplyr::pull(line)
      
      readr::write_csv(fas_raw_dat_tib, fas_raw_all_csv)
      readr::write_lines(fas_raw_all_vec, fas_raw_all_fas)
    } else {
      raw_all_tib <- readr::read_csv(fas_raw_all_csv) %>% # head() %>%
        dplyr::select(FAS_KEY,FAS_SEQ) %>%
        tidyr::separate(FAS_KEY, into=c("prb_cgn","prb_srd","prb_inf","prb_src"), sep="_") %>%
        dplyr::mutate(prb_inf=as.integer(prb_inf)) %>%
        dplyr::rename(prb_seq=FAS_SEQ)
    }
    
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
      )
    
    # Quick Validation Summary::
    raw_all_bsp_tib %>% 
      dplyr::group_by(bsp_srd,prb_mat) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print(n=100)
    
    raw_all_bsp_tib %>%
      dplyr::select(prb_cgn,prb_srd,prb_inf,prb_src,prb_seq,bsp_ref,bsp_srd,prb_mat,CG_DIN) %>%
      dplyr::group_by(CG_DIN,prb_inf,bsp_srd,prb_mat) %>%
      dplyr::summarise(Count=n(), .groups="drop") %>% 
      dplyr::arrange(-Count) %>% print(n=100)
    
    # Load True Coordinates from CG# Database::
    #
    # Slow method::
    #  cgn_pos_db_tsv <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-seq.tsv"
    # cgn_pos_db_col <- cols(
    #   Seq_ID     = col_character(),
    #   Chromosome = col_character(),
    #   Coordinate = col_integer(),
    #   
    #   Methyl_Allele_FR_Strand = col_character(),
    #   Methyl_Allele_TB_Strand = col_character(),
    #   Methyl_Allele_CO_Strand = col_character(),
    #   Methyl_Next_Base        = col_character()
    # )
    # cgn_pos_db_tib <- readr::read_tsv(cgn_pos_db_tsv, col_types=cgn_pos_db_col)
    #
    # Old Method Join::
    #  cgn_pos_db_tib %>% dplyr::inner_join(raw_all_bsp_tib, by=c("Chromosome"="bsp_cgn", "Coordinate"="CG_Pos"))
    #

    # This was too slow previously...    
    # if (FALSE) {
    #   cgn_pos_db_tsv <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-seq.cgn-sorted.tsv.gz"
    #   cgn_pos_db_col <- cols(
    #     imp_cgn = col_character(),
    #     imp_chr = col_character(),
    #     imp_pos = col_integer(),
    #     
    #     imp_fr = col_character(),
    #     imp_tb = col_character(),
    #     imp_co = col_character(),
    #     imp_nb = col_character(),
    #     imp_1U = col_character()
    #   )
    #   cgn_pos_db_tib <- readr::read_tsv(cgn_pos_db_tsv, 
    #                                     col_names=names(cgn_pos_db_col$cols), 
    #                                     col_types=cgn_pos_db_col)
    #   
    #   cgn_pos_grs_tib <- cgn_pos_db_tib %>% # head(n=10000) %>%
    #     dplyr::mutate(
    #       strand=dplyr::case_when(
    #         imp_fr=="F" ~ "+",
    #         imp_fr=="R" ~ "-",
    #         TRUE ~ NA_character_
    #       ),
    #       imp_tag=paste0(imp_fr,imp_tb,imp_co,imp_nb),
    #       Probe_ID_IMP=paste(imp_cgn, imp_tag,imp_chr,imp_pos, sep="_"),
    #       imp_chr=paste0("chr",imp_chr)
    #     )
    #   
    #   #
    #   # Next: Create GRS and save RDS for future use...
    #   #
    #   imp_cgn_grs <- GRanges(
    #     seqnames = Rle(cgn_pos_grs_tib$imp_chr),
    #     strand=Rle(cgn_pos_grs_tib$strand),
    #     prb1U=cgn_pos_grs_tib$imp_1U,
    #     fr=cgn_pos_grs_tib$imp_fr,
    #     tb=cgn_pos_grs_tib$imp_tb,
    #     co=cgn_pos_grs_tib$imp_co,
    #     IRanges(start=cgn_pos_grs_tib$imp_pos,
    #             end=cgn_pos_grs_tib$imp_pos+1,
    #             names=cgn_pos_grs_tib$Probe_ID_IMP) )
    # }
    
    #
    # Latest Starting Place::
    #
    cgn_pos_db_csv <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-seq.cgn-sorted.csv"
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
    
    # Preprocessed for faster future loading...
    cgn_pos_db2_csv <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-end-seq.cgn-sorted.csv.gz"
    # readr::write_csv(cgn_pos_db_tib, cgn_pos_db2_csv)
    
    #
    # Next: Create GRS and save RDS for future use...
    #
    # TBD:: Update this::
    cgn_pos_db_rds <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37-21092020_improbe-designOutput.cgn-map-end-seq.cgn-sorted.rds"
    imp_cgn_grs <- GRanges(
      seqnames = Rle(cgn_pos_db_tib$imp_chr),
      strand=Rle(cgn_pos_db_tib$imp_frs),
      imp_seq1U=cgn_pos_db_tib$imp_seq,
      IRanges(start=cgn_pos_db_tib$imp_pos,
              end=cgn_pos_db_tib$imp_end,
              names=cgn_pos_db_tib$Probe_ID_IMP) )
    readr::write_rds(x = imp_cgn_grs, file = cgn_pos_db_rds, compress ="gz")    
    
    
    
    
    
  }

  fas_test_tib <- NULL
  fas_test_tib <- fas_raw_dat_tib %>% 
    dplyr::inner_join( all_prb_bsp_tib, by=c("FAS_SEQ"="prb_seq" ) ) %>% 
    dplyr::distinct(FAS_SEQ, .keep_all=TRUE) %>%
    dplyr::mutate(FAS_KEY=paste(FAS_KEY,bsp_tag, sep="_")) %>%
    dplyr::group_by(Manifest,Design_Type,bsp_tag) %>% 
    dplyr::slice_head(n=3) %>%
    dplyr::ungroup()
  
  fas_test_sum <- NULL
  fas_test_sum <- fas_test_tib %>% 
    dplyr::group_by(Manifest,Design_Type,bsp_tag) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  fas_test_vec <- NULL
  fas_test_vec <- fas_test_tib %>%
    dplyr::mutate(
      line=paste0('>',FAS_KEY,'\n',FAS_SEQ)
    ) %>%
    dplyr::filter(!is.na(line)) %>% 
    dplyr::pull(line)
  
  readr::write_csv(fas_test_tib, test_csv)
  readr::write_lines(fas_test_vec, test_fas)

  test_cmd <- glue::glue("{bsmap_exe} -a {test_fas} -d {test_ref} -o {test_bsp} -s 10 -v 5 -n 1 -r 2 -V 2 -R -u {RET}",
                         "cat {test_bsp} | cut -f 1,2,4-11 | perl -pe 's/_/\t/gi; s/:/\t/gi;' | gzip -c - > {test_bsg}")
  gn37_cmd <- glue::glue("{bsmap_exe} -a {test_fas} -d {gn37_ref} -o {gn37_bsp} -s 10 -v 5 -n 1 -r 2 -V 2 -R -u {RET}",
                         "cat {gn37_bsp} | cut -f 1,2,4-11 | perl -pe 's/_/\t/gi; s/:/\t/gi;' | gzip -c - > {gn37_bsg}")
  gn38_cmd <- glue::glue("{bsmap_exe} -a {test_fas} -d {gn38_ref} -o {gn38_bsp} -s 10 -v 5 -n 1 -r 2 -V 2 -R -u {RET}",
                         "cat {gn38_bsp} | cut -f 1,2,4-11 | perl -pe 's/_/\t/gi; s/:/\t/gi;' | gzip -c - > {gn38_bsg}")

  readr::write_lines(test_cmd, test_ssh)
  readr::write_lines(gn37_cmd, gn37_ssh)
  readr::write_lines(gn38_cmd, gn38_ssh)
  
  system(glue::glue("chmod 777 {test_ssh}"))
  system(glue::glue("chmod 777 {gn37_ssh}"))
  system(glue::glue("chmod 777 {gn38_ssh}"))
  
  #
  # Load Mini Results::
  #
  test_bsg_col <- cols(
    # bsp_key = col_character()
    prb_cgn = col_character(),
    prb_srd = col_character(),
    prb_inf = col_integer(),
    prb_src = col_character(),
    prb_tag = col_character(),
    
    bsp_seq = col_character(),
    bsp_tag = col_character(),
    
    # bsp_chr = col_character()
    bsp_cgn = col_character(),
    
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
  
  #
  # Loda BSP Formatted Data::
  #
  test_bsg_tib <- readr::read_tsv(
    test_bsg, 
    col_names=names(test_bsg_col$cols), 
    col_types=test_bsg_col) %>%
    dplyr::mutate(rvc_seq=revCmp(bsp_seq))
  # test_bsg_tib %>% dplyr::group_by(bsp_tag,prb_tag,prb_src,prb_inf) %>% dplyr::summarise(Count=n(), .groups="drop")
  test_bsg_tib %>% dplyr::group_by(bsp_tag,prb_tag,prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")

  gn37_bsg_tib <- readr::read_tsv(
    gn37_bsg, 
    col_names=names(test_bsg_col$cols), 
    col_types=test_bsg_col) %>%
    dplyr::mutate(rvc_seq=revCmp(bsp_seq))
  # gn37_bsg_tib %>% dplyr::group_by(bsp_tag,prb_tag,prb_src,prb_inf) %>% dplyr::summarise(Count=n(), .groups="drop")
  gn37_bsg_tib %>% dplyr::group_by(bsp_tag,prb_tag,prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
  
  gn38_bsg_tib <- readr::read_tsv(
    gn38_bsg, 
    col_names=names(test_bsg_col$cols), 
    col_types=test_bsg_col) %>%
    dplyr::mutate(rvc_seq=revCmp(bsp_seq))
  # gn38_bsg_tib %>% dplyr::group_by(bsp_tag,prb_tag,prb_src,prb_inf) %>% dplyr::summarise(Count=n(), .groups="drop")
  gn38_bsg_tib %>% dplyr::group_by(bsp_tag,prb_tag,prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
  
  
  #
  # Summarise BSP Alignments::
  #
  test_bsg_tib %>% dplyr::group_by(bsp_tag,prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
  gn37_bsg_tib %>% dplyr::group_by(bsp_tag,prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
  gn38_bsg_tib %>% dplyr::group_by(bsp_tag,prb_src) %>% dplyr::summarise(Count=n(), .groups="drop")
  
  #
  # Investigate CgnTop (test) Alignments:: NM
  #
  test_nm_tib <- test_bsg_tib %>% dplyr::filter(bsp_tag=="NM") %>% 
    dplyr::select(prb_cgn:bsp_tag) %>% dplyr::distinct()

  # gn37_bsg_tib %>% dplyr::filter(bsp_seq %in% test_nm_tib$bsp_seq | bsp_seq %in% test_nm_tib$rvc_seq)
  
  dplyr::bind_rows(
    v1_tib <- dplyr::inner_join(test_nm_tib, gn37_bsg_tib, by=c("bsp_seq","rvc_seq"), suffix=c("_CGN", "_H37") ),
    test_nm_tib %>% dplyr::anti_join(v1_tib, by="bsp_seq") %>%
      dplyr::inner_join(gn37_bsg_tib, by=c("bsp_seq"="rvc_seq"), suffix=c("_CGN", "_H37")) %>%
      dplyr::select(-bsp_tag_H37)
  ) %>% dplyr::arrange(bsp_seq) %>% as.data.frame()
  
  # RBC: CCTTTCTCAACTACAAACAACAACTATAAAAAATAACAAAACAACAAACA
  # PRB: CCTTTCTCAACTACAAACAACAACTATAAAAAATAACAAAACAACAAACA
  # REV: TGTTTGTTGTTTTGTTATTTTTTATAGTTGTTGTTTGTAGTTGAGAAAGG
  # 
  # >cg08065231 FCC21
  #  TOP(refN): AACCCTTTCTCAGCTGCAGGCGACAGCTATGGGAGATGGCGAGGCGACGGA CG CCCACTTCCGCTTCCGGCCAGCAGCTGAAAGGGGTGGGCTGCCCAGCTCGC
  #    TGGGACAGCAACCCTTTCTCAGCTGCAGGCGACAGCTATGGGAGATGGCGAGGCGACGGA[CG]CCCACTTCCGCTTCCGGCCAGCAGCTGAAAGGGGTGGGCTGCCCAGCTCGCACATTTGGC
  
  #  BOT(refN): TTGGGAAAGAGTCGACGTCCGCTGTCGATACCCTCTACCGCTCCGCTGCCT GC GGGTGAAGGCGAAGGCCGGTCGTCGACTTTCCCCACCCGACGGGTCGAGCG
  #    ACCCTGTCGTTGGGAAAGAGTCGACGTCCGCTGTCGATACCCTCTACCGCTCCGCTGCCT[GC]GGGTGAAGGCGAAGGCCGGTCGTCGACTTTCCCCACCCGACGGGTCGAGCGTGTAAACCG
  #  PRB(REV)                                                        TGTTTGTTGTTTTGTTATTTTTTATAGTTGTTGTTTGTAGTTGAGAAAGG
  #      ATTTTGTTGTTGGGAAAGAGTTGATGTTTGTTGTTGATATTTTTTATTGTTTTGTTGTTTGTGGGTGAAGGTGAAGGTTGGTTGTTGATTTTTTTTATTTGATGGGTTGAGTGTGTAAATTG
  #               TTGGGAAAGAGTTGATGTTTGTTGTTGATATTTTTTATTGTTTTGTTGTTTGTGGGTGAAGGTGAAGGTTGGTTGTTGATTTTTTTTATTTGATGGGTTGAGTG
  
  dplyr::bind_rows(
    dplyr::inner_join(test_nm_tib, gn37_bsg_tib, by=c("bsp_seq","rvc_seq"), suffix=c("_CGN", "_H37") ),
    dplyr::inner_join(test_nm_tib, gn37_bsg_tib, by=c("bsp_seq"="rvc_seq", "rvc_seq"="bsp_seq"), suffix=c("_CGN", "_H37"))
  ) %>% dplyr::arrange(prb_cgn_CGN) %>% as.data.frame()
    
  
  test_nm_tib %>% dplyr::inner_join(gn37_bsg_tib, by="bsp_seq", suffix=c("_CGN", "_H37"))
  
  test_nm_tib %>% dplyr::inner_join(gn37_bsg_tib, by="bsp_seq", suffix=c("_CGN", "_H37")) %>% dplyr::filter(prb_cgn_CGN=="cg08065231")
  test_nm_tib %>% dplyr::inner_join(gn37_bsg_tib, by="bsp_seq", suffix=c("_CGN", "_H37")) %>% dplyr::filter(bsp_seq=="CCTTTCTCAACTACAAACAACAACTATAAAAAATAACAAAACAACAAACA")

  
  test_nm_tib %>% dplyr::anti_join(gn37_bsg_tib, by="bsp_seq")
  
  # How do you have a NM in CGN vs. MA/UM in Genomic???
  #  NOTE: This join is done via cgn not probe sequence!!!
  odd_nm_tib <- test_nm_tib %>% dplyr::inner_join(gn37_bsg_tib, by="prb_cgn", suffix=c("_CGN", "_H37"))
  odd_nm_tib %>% 
    dplyr::select(bsp_tag_CGN,bsp_tag_H37) %>% 
    dplyr::group_by_all() %>% dplyr::summarise(Count=n(), .groups="drop")
  
  odd_nm_gn37_tib <- gn37_bsg_tib %>% dplyr::filter(prb_cgn %in% test_nm_tib$prb_cgn)
  
  bscU("AACCCTTTCTCAGCTGCAGGCGACAGCTATGGGAGATGGCGAGGCGACGGACGCCCACTTCCGCTTCCGGCCAGCAGCTGAAAGGGGTGGGCTGCCCAGCTCGC")

  # test_nm_tib::
  # cg08065231 FCC21         1 B2      NM      CCTTTCTCAACTACAAACAACAACTATAAAAAATAACAAAACAACAAACA NM

  # cg08065231 FCC21
  #  TOP(refN): AACCCTTTCTCAGCTGCAGGCGACAGCTATGGGAGATGGCGAGGCGACGGACGCCCACTTCCGCTTCCGGCCAGCAGCTGAAAGGGGTGGGCTGCCCAGCTCGC
  #  BOT(refN): TTGGGAAAGAGTCGACGTCCGCTGTCGATACCCTCTACCGCTCCGCTGCCTGCGGGTGAAGGCGAAGGCCGGTCGTCGACTTTCCCCACCCGACGGGTCGAGCG
  #
  #  5' TGTTTGTTGTTTTGTTATTTTTTATAGTTGTTGTTTGTAGTTGAGAAAGG 3'  revCmpl CCTTTCTCAACTACAAACAACAACTATAAAAAATAACAAAACAACAAACA
  #  3' ACAAACAACAAAACAATAAAAAATATCAACAACAAACATCAACTCTTTCC 5'  revCmpl GGAAAGAGTTGATGTTTGTTGTTGATATTTTTTATTGTTTTGTTGTTTGT
  #
  #                                                             
  #  TOP(bscU): AATTTTTTTTTAGTTGTAGGTGATAGTTATGGGAGATGGTGAGGTGATGGATGTTTATTTTTGTTTTTGGTTAGTAGTTGAAAGGGGTGGGTTGTTTAGTTTGT
  #             TTAAAAAAAAATCAACATCCACTATCAATACCCTCTACCACTCCACTACCTACAAATAAAAACAAAAACCAATCATCAACTTTCCCCACCCAACAAATCAAACA
  #
  #             CACTCAACCCATCAAATAAAAAAAATCAACAACCAACCTTCACCTTCACCCACAAACAACAAAACAATAAAAAATATCAACAACAAACATCAACTCTTTCCCAA
  #  BOT(bscU): GTGAGTTGGGTAGTTTATTTTTTTTAGTTGTTGGTTGGAAGTGGAAGTGGGTGTTTGTTGTTTTGTTATTTTTTATAGTTGTTGTTTGTAGTTGAGAAAGGGTT
  #                                                                TGTTTGTTGTTTTGTTATTTTTTATAGTTGTTGTTTGTAGTTGAGAAAGG
  #  BOT(refN): TTGGGAAAGAGTCGACGTCCGCTGTCGATACCCTCTACCGCTCCGCTGCCTGCGGGTGAAGGCGAAGGCCGGTCGTCGACTTTCCCCACCCGACGGGTCGAGCG
  #
  

  
}


#
# This is it::
#  - Use hg37 cords (chrom_GRCh37,chromStart,chromEnd_GRCh37,FR_GRCh37) for non-COVIC
#    - 
#  - use new positions (man_pqc_grs_tibs$GRCh37) for COVIC
#

# EPIC:: B2/B4
man_raw_org_tib <- man_raw_dat_tib %>% 
  dplyr::select(!dplyr::ends_with("_GRCh38"), -AQP) %>%
  dplyr::filter(!is.na(chrom_GRCh37) & !is.na(chromStart_GRCh37)) %>%
  dplyr::distinct(U,M, .keep_all=TRUE)

man_raw_org_tib %>% dplyr::group_by(Manifest) %>% dplyr::summarise(Count=n(), .groups="drop")
man_pqc_all_tib %>% dplyr::group_by(Manifest) %>% dplyr::summarise(Count=n(), .groups="drop")

# man_pqc_all_tib %>%
#   dplyr::filter(!is.na(chrom_GRCh37) & !is.na(chromStart_GRCh37)) %>%
#   dplyr::mutate(POS_DIF=chromEnd_GRCh37-chromStart_GRCh37) %>%
#   dplyr::group_by(POS_DIF) %>% dplyr::summarise(Count=n(), .groups="drop")
# 
# man_pqc_all_tib %>%
#   dplyr::filter(!is.na(chrom_GRCh37) & !is.na(chromStart_GRCh37)) %>%
#   dplyr::group_by(Manifest,DESIGN) %>% 
#   dplyr::summarise(Count=n(), .groups="drop")

#
# COVIC
#
man_pqc_cov_tib <- 
  man_pqc_grs_tibs$GRCh37 %>% 
  dplyr::filter(is.na(chrom_GRCh37)) %>%
  dplyr::add_count(chrom,chromStart, name="IMP_MAP_Count") %>%
  dplyr::add_count(CGN_IMP, name="IMP_CGN_Count")

man_cov_grs_tib <- man_pqc_cov_tib %>% 
  dplyr::rename(Seq_ID_IMP=CGN_IMP) %>%
  dplyr::mutate(
    strand=dplyr::case_when(
      FR_IMP=="F" ~ "+", FR_IMP=="R" ~ "-", TRUE ~ NA_character_), 
    TB_CO_IMP=paste0(TB_IMP,CO_IMP), 
    Probe_ID_IMP=paste(Seq_ID_IMP,TB_CO_IMP, sep="_")) %>% 
  dplyr::select(chrom,chromStart,chromEnd,strand,Probe_ID_IMP,
                U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Probe_Type,
                Probe_ID,Source_ID,Source_Tag,Next_Base,
                Seq_ID_Uniq,Seq_ID,Strand_TB,Strand_CO,NXB_D,
                DesSeqN,DesBscD) %>%
  dplyr::rename(Probe_ID_SRC=Probe_ID, 
                Probe_ID_DES=Seq_ID_Uniq,
                Seq_ID_DES=Seq_ID,
                Seq_ID=Probe_ID_IMP)

#
# Annotation:: COVIC
#
opt$annDir

tar_ann_vec <- c("NCBI_Genes")
fin_man_ana_tib <- NULL
if (!is.null(opt$annDir) && dir.exists(opt$annDir)) {
  man_cov_ana_tib <- 
    manifestToAnnotation(tib=man_cov_grs_tib,tar=tar_ann_vec,
                         ann=opt$annDir, gen=opt$genomeBuild,
                         verbose=opt$verbose,vt=1,tc=1,tt=pTracker)
  
  
  fin_man_ana_tib <- 
    manifestToAnnotation(tib=fin_man_grs_tib,tar=tar_ann_vec,
                         ann=opt$annDir, gen=opt$genomeBuild,
                         verbose=opt$verbose,vt=1,tc=1,tt=pTracker)
  
}


  
#  How is cgnFLag Calculated:: "CGN_IMP==Seq_ID ~ 0"
#  Where does Seq_ID come from?
man_pqc_cov_tib %>%
  dplyr::distinct(U,M, .keep_all=TRUE) %>%
  dplyr::filter(IMP_CGN_Count != 1) %>% 
  dplyr::group_by(Manifest,cgnFlag,chrFlag,difFlag,srcFlag) %>%
  dplyr::summarise(Count=n(), .groups="drop")


# Should Use these ones:::
man_pqc_cov_tib %>% 
  dplyr::filter(IMP_CGN_Count != 1) %>% 
  dplyr::filter(CGN_IMP==Source_ID) %>% 
  dplyr::distinct(CGN_IMP)




man_pqc_cov_tib %>% 
  dplyr::filter(IMP_CGN_Count != 1) %>%
  dplyr::arrange(-IMP_CGN_Count, CGN_IMP) %>%
  dplyr::select(CGN_IMP,PRB1_U_MAT,IMP_CGN_Count, chrom,chromStart)

man_pqc_cov_tib %>% 
  # dplyr::group_by(Manifest,DESIGN,IMP_CGN_Count) %>% 
  dplyr::group_by(Manifest,DESIGN,IMP_MAP_Count) %>% 
  dplyr::summarise(Count=n(), .groups="drop") %>%
  print(n=1000)

# Find Probes without obvious mappings::
man_pqc_cov_tib %>% 
  dplyr::filter(IMP_CGN_Count != 1) %>%
  dplyr::filter(IMP_MAP_Count != 1) %>%
  dplyr::select(1:10)
  
  



man_pqc_grs_tibs$GRCh37 %>% 
  dplyr::filter(is.na(chrom_GRCh37)) %>% 
  dplyr::distinct(U,M, .keep_all=TRUE) %>%
  dplyr::group_by(Manifest,cgnFlag,chrFlag,difFlag,srcFlag) %>%
  dplyr::summarise(Count=n(), .groups="drop")


# dplyr::select(chrom,chromStart,chromEnd,Probe_ID) %>%
  

man_pqc_grs_tibs$GRCh38 %>% 
  dplyr::filter(is.na(chrom_GRCh38)) %>% 
  dplyr::distinct(U,M, .keep_all=TRUE) %>%
  dplyr::group_by(Manifest,cgnFlag,chrFlag,difFlag,srcFlag) %>%
  dplyr::summarise(Count=n(), .groups="drop")







if (FALSE) {

  

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
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
