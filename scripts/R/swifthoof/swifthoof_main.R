
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Without Warnings for Testing::
#
# base::require("sesame")
# base::require("dbplyr")
# 
# base::require("optparse",quietly=TRUE)
# 
# base::require("tidyverse")
# base::require("plyr")
# base::require("stringr")
# base::require("readr")
# base::require("glue")
# 
# base::require("matrixStats")
# base::require("scales")
#
# base::require("doParallel")

# suppressPackageStartupMessages(base::require() )
# Load sesame:: This causes issues with "ExperimentHub Caching causes a warning"
suppressWarnings(suppressPackageStartupMessages( base::require("sesame") ))
suppressWarnings(suppressPackageStartupMessages( base::require("dbplyr") ))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("plyr")) )
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("readr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

# tidyverse_update(recursive = FALSE, repos = getOption("repos"))
# install.packages("lubridate")

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
par$prgmDir <- 'swifthoof'
par$prgmTag <- paste(par$prgmDir,'main', sep='_')
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

par$retData     <- FALSE

# Executables::
opt$Rscript <- NULL

# Directories::
opt$outDir     <- NULL
opt$idatsDir   <- NULL

# Optional Files::
opt$subManifest  <- FALSE
opt$manifestPath <- 'auto'
opt$auto_sam_csv <- NULL

# Platform/Method Options::
opt$platform  <- NULL
opt$manifest  <- NULL

# Run Options::
opt$fresh        <- FALSE
opt$buildSubDir  <- FALSE
opt$auto_detect  <- FALSE

opt$workflow    <- NULL
opt$manDirName  <- 'core'

# Output Options::
opt$load_idat    <- FALSE
opt$save_idat    <- TRUE

opt$load_sset   <- FALSE
opt$save_sset   <- FALSE

#
# TBD: Add new variable names::
#
opt$write_beta  <- FALSE
opt$write_bsum  <- FALSE

opt$write_pval  <- FALSE
opt$write_psum  <- FALSE

opt$write_sigs  <- FALSE
opt$write_ssum  <- FALSE

opt$write_call  <- FALSE
opt$write_csum  <- FALSE

opt$write_snps  <- TRUE
opt$write_auto  <- FALSE

opt$mask_general <- FALSE

# Threshold Options::
opt$pval <- "pOOBAH,PnegEcdf"
opt$minPval <- "0.1,0.02"
opt$minPerc <- "90,98"

opt$minDeltaBeta <- 0.2

opt$percision_sigs <- 1
opt$percision_beta <- 4
opt$percision_pval <- 6

# Parallel/Cluster Options::
opt$single   <- FALSE
opt$parallel <- FALSE
opt$cluster  <- FALSE

# Plotting Options::
opt$plotSset  <- FALSE
opt$plotCalls <- FALSE
opt$plotAuto  <- FALSE

opt$make_pred <- TRUE

opt$plotFormat <- 'pdf'
opt$plotFormat <- 'png'

opt$dpi <- 72
opt$dpi <- 120

opt$plotMax <- 10000
opt$plotSub <- 5000

opt$time_org_txt <- NULL

# verbose Options::
opt$verbose <- 3

# Set Default Options::
#
def <- opt

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
  
  opt$outDir <- file.path(par$topDir, 'scratch')
  locIdatDir <- file.path(par$topDir, 'data/idats')
  
  opt$buildSubDir  <- FALSE
  opt$auto_detect   <- FALSE

  opt$runName    <- NULL
  opt$platform   <- 'EPIC'
  opt$manifest   <- 'B4'
  opt$platform   <- NULL
  opt$manifest   <- NULL
  
  par$expChipNum <- NULL
  par$expSampNum <- NULL
  
  # Writing Options::
  opt$write_beta  <- FALSE
  opt$write_bsum  <- TRUE
  
  opt$write_pval  <- FALSE
  opt$write_psum  <- TRUE
  
  opt$write_sigs  <- TRUE
  opt$write_ssum  <- TRUE
  
  opt$write_call  <- TRUE
  opt$write_csum  <- TRUE
  
  opt$workflow <- "nd,ind"
  opt$workflow <- 'ind'
  
  opt$auto_detect <- TRUE
  
  par$retData  <- TRUE
  
  opt$single   <- FALSE
  opt$single   <- TRUE
  
  opt$parallel <- TRUE
  opt$parallel <- FALSE
  
  opt$cluster  <- TRUE
  opt$cluster  <- FALSE
  
  opt$manDirName  <- 'base'
  opt$manDirName  <- 'core'
  
  opt$verbose  <- 6
  
  par$local_runType <- 'CORE'
  par$local_runType <- 'EXCBR'
  par$local_runType <- 'GRCm38'
  par$local_runType <- 'COVID'
  par$local_runType <- 'COVIC'
  par$local_runType <- 'GRCm38'
  par$local_runType <- 'DELTA'
  par$local_runType <- 'DKFZ'
  par$local_runType <- 'qcMVP'
  
  opt$fresh <- TRUE
  
  if (par$local_runType=='COVID') {
    opt$runName  <- 'COVID-Direct-Set1'
    par$expChipNum <- '204756130014'
    
    opt$auto_detect <- FALSE
    # opt$single   <- FALSE
    par$retData  <- TRUE
    
    opt$workflow <- "nd,ind"

  } else if (par$local_runType=='COVIC') {
    opt$runName  <- 'COVIC-Set1-15052020'
    par$expChipNum <- '204500250013'
    opt$auto_detect <- TRUE
    opt$workflow <- "ind"
    opt$workflow <- "nd,ind"
    
  } else if (par$local_runType=='GRCm38') {
    opt$runName <- 'MURMETVEP_mm10_betaTest_06082020'
    opt$runName <- 'VanAndel_mm10_betaTest_31082020'
    opt$runName <- 'ILMN_mm10_betaTest_17082020'
    opt$auto_detect <- FALSE
  } else if (par$local_runType=='qcMVP') {
    opt$runName  <- 'CNTL-Samples_VendA_10092020'
    par$expChipNum <- '203962710025'
    par$expSampNum <- '203962710025_R08C01'
    
    par$expChipNum <- "203962710079"
    par$expSampNum <- "203962710079_R08C01"
    
    opt$auto_detect <- TRUE
    opt$dpi <- 72
  } else if (par$local_runType=='DKFZ') {
    opt$runName  <- 'DKFZ'

    opt$auto_detect <- TRUE
    opt$dpi <- 72
  } else if (par$local_runType=='CORE') {
    opt$runName  <- 'BETA-8x1-EPIC-Core'
    par$expChipNum <- '202761400007'

    opt$runName  <- 'ADRN-blood-nonAtopic_EPIC'
    par$expChipNum <- '201125090068'

    opt$runName  <- 'GSE122126_EPIC'
    par$expChipNum <- '202410280180'

    opt$runName  <- 'EPIC-BETA-8x1-CoreCancer'
    par$expChipNum <- '201502830033'
    opt$auto_detect <- TRUE
  } else if (par$local_runType=='EXCBR') {
    opt$runName  <- 'Excalibur-Old-1609202'
    par$expChipNum <- '204076530053'
    par$expChipNum <- '204076530110'
    
    opt$runName  <- 'Excalibur-New-1609202'
    par$expChipNum <- '202915460071'
    opt$auto_detect <- TRUE
  } else if (par$local_runType=='DELTA') {
    opt$runName    <- 'DELTA-8x1-EPIC-Core'
    par$expChipNum <- '203319730022'
    par$expSampNum <- '203319730022_R07C01'
    
    opt$auto_detect <- TRUE
    
    opt$plotSset  <- TRUE
    opt$plotCalls <- TRUE
    opt$plotAuto  <- TRUE
    
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unrecognized local_runType={par$local_runType}.{RET}{RET}"))
  }
  
  opt$save_idat <- TRUE
  opt$load_idat <- TRUE
  
  opt$save_sset <- TRUE
  opt$load_sset <- TRUE

  opt$auto_detect <- FALSE
  opt$auto_detect <- TRUE
  
  opt$idatsDir <- file.path(locIdatDir, paste('idats',opt$runName, sep='_') )
  if (!is.null(par$expChipNum)) opt$idatsDir <- file.path(locIdatDir, paste('idats',opt$runName, sep='_'),  par$expChipNum)
  opt$auto_sam_csv <- file.path(par$datDir, 'ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo_Mean-COVIC-280-NP-ind_negs-0.02.csv.gz')
  
  # opt$outDir <- file.path(par$topDir, 'scratch', par$local_runType, par$prgmTag, opt$runName)
  opt$outDir <- file.path(par$topDir, 'scratch',par$runMode)
  
} else {
  par$runMode    <- 'CommandLine'
  par$exePath <- base::substring(args.dat[grep("--file=", args.dat)], 8)
  
  cat(glue::glue("[{par$prgmTag}]: Local Run par$runMode={par$runMode}.{RET}"))
  
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
    
    # Directories::
    make_option(c("--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    make_option(c("--idatsDir"), type="character", default=opt$idatsDir, 
                help="idats directory [default= %default]", metavar="character"),

    # Optional Files::
    make_option(c("--manifestPath"), type="character", default=opt$manifestPath,
                help="Path to manfifest (CSV) otherwise use dat [default= %default]", metavar="character"),
    make_option(c("--subManifest"), action="store_true", default=opt$subManifest,
                help="Boolean variable to use subset manifest instead of subset. [default= %default]", metavar="boolean"),
    make_option(c("--auto_sam_csv"), type="character", default=opt$auto_sam_csv,
                help="Path to auto detect beta values (CSV) [default= %default]", metavar="character"),
    
    # Platform/Method Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Forced platform [EPIC, 450k, 27k, NZT] otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--manifest"), type="character", default=opt$manifest, 
                help="Forced manifest [B1, B2, B4] otherwise auto-detect [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--fresh"), action="store_true", default=opt$fresh,
                help="Boolean variable to build fresh version of database files [default= %default]", metavar="boolean"),
    make_option(c("--buildSubDir"), action="store_true", default=opt$buildSubDir,
                help="Boolean variable to build subdirectories based on Chip/BeadPool (for R&D purposes) [default= %default]", metavar="boolean"),
    make_option(c("--auto_detect"), action="store_true", default=opt$auto_detect,
                help="Boolean variable to auto detect reference samples. Must provide reference samples. [default= %default]", metavar="boolean"),

    make_option(c("--workflow"), type="character", default=opt$workflow,
                help="Workflows to run i.e. order of operations comma seperated [ raw,ind,ndi,din, etc. ] [default= %default]", metavar="character"),
    make_option(c("--manDirName"), type="character", default=opt$manDirName,
                help="Manifest directory name [default= %default]", metavar="character"),
    
    # Output Options::
    make_option(c("--load_idat"), action="store_true", default=opt$load_idat,
                help="Boolean variable to load existing IDAT from RDS file [default= %default]", metavar="boolean"),
    make_option(c("--save_idat"), action="store_true", default=opt$save_idat,
                help="Boolean variable to write IDAT RDS file [default= %default]", metavar="boolean"),
    
    make_option(c("--load_sset"), action="store_true", default=opt$load_sset,
                help="Boolean variable to load existing Signal Set from RDS file [default= %default]", metavar="boolean"),
    make_option(c("--save_sset"), action="store_true", default=opt$save_sset,
                help="Boolean variable to write Signal Set RDS file [default= %default]", metavar="boolean"),

    #
    # Old Versions to be deleted::
    #
    make_option(c("--write_snps"), action="store_true", default=opt$write_snps,
                help="Boolean variable to write direct and inferred SNP calls file [default= %default]", metavar="boolean"),
    make_option(c("--write_auto"), action="store_true", default=opt$write_auto,
                help="Boolean variable to write Auto-Detection Matricies (Pval/Beta) file [default= %default]", metavar="boolean"),
    make_option(c("--mask_general"), action="store_true", default=opt$mask_general,
                help="Boolean variable to report Sesame masked detection p-value stats [default= %default]", metavar="boolean"),

    #
    # Current Versions::
    #
    make_option(c("--write_beta"), action="store_true", default=opt$write_beta,
                help="Boolean variable to write Beta Set file (CSV) [default= %default]", metavar="boolean"),
    make_option(c("--write_bsum"), action="store_true", default=opt$write_bsum,
                help="Boolean variable to write Beta Set Summary file (CSV) [default= %default]", metavar="boolean"),
    
    make_option(c("--write_pval"), action="store_true", default=opt$write_pval,
                help="Boolean variable to write Pval Set file (CSV) [default= %default]", metavar="boolean"),
    make_option(c("--write_psum"), action="store_true", default=opt$write_psum,
                help="Boolean variable to write Pval Set Summary file (CSV) [default= %default]", metavar="boolean"),
    
    make_option(c("--write_sigs"), action="store_true", default=opt$write_sigs,
                help="Boolean variable to write Signal Set file (CSV) [default= %default]", metavar="boolean"),
    make_option(c("--write_ssum"), action="store_true", default=opt$write_ssum,
                help="Boolean variable to write Signal Set Summary file (CSV) [default= %default]", metavar="boolean"),
    
    make_option(c("--write_call"), action="store_true", default=opt$write_call,
                help="Boolean variable to write Calls (Pval/Beta) file (CSV) [default= %default]", metavar="boolean"),
    make_option(c("--write_csum"), action="store_true", default=opt$write_csum,
                help="Boolean variable to write Calls Summary file (CSV) [default= %default]", metavar="boolean"),
    make_option(c("--make_pred"), action="store_true", default=opt$make_pred,
                help="Boolean variable to make prediction calls (mostly for RD stuff) [default= %default]", metavar="boolean"),
    
    #
    # Threshold Options::
    #
    make_option(c("--pval"), type="character", default=opt$pval, 
                help="Detection p-value methods (commas seperated list)  [default= %default]", metavar="character"),
    make_option(c("--minPval"), type="character", default=opt$minPval, 
                help="Minimum passing detection p-value for each pval method (commas seperated list)  [default= %default]", metavar="character"),
    make_option(c("--minPerc"), type="character", default=opt$minPerc, 
                help="Minimum percentage of loci passing detection p-value for each pval method (commas seperated list)  [default= %default]", metavar="character"),
    
    # Old methods::
    #
    # make_option(c("--minNegPval"), type="double", default=opt$minNegPval, 
    #             help="Minimum passing detection p-value using Negative Controls [default= %default]", metavar="double"),
    # make_option(c("--minOobPval"), type="double", default=opt$minOobPval,
    #             help="Minimum passing detection p-value using Negative Out-Of-Band [default= %default]", metavar="double"),
    # 
    # make_option(c("--minNegPerc"), type="double", default=opt$minNegPerc, 
    #             help="Minimum percentage of loci passing detection p-value using Negative Controls to flag Requeue of sample. [default= %default]", metavar="double"),
    # make_option(c("--minOobPerc"), type="double", default=opt$minOobPerc, 
    #             help="Minimum percentage of loci passing detection p-value using Out-Of-Band to flag Requeue of sample. [default= %default]", metavar="double"),
    
    make_option(c("--minDeltaBeta"), type="double", default=opt$minDeltaBeta,
                help="Minimum passing delta-beta. Used in AutoSampleSheet cacluclations [default= %default]", metavar="double"),
    
    make_option(c("--percision_sigs"), type="integer", default=opt$percision_sigs,
                help="Rounding percision for signal values in calls output files [default= %default]", metavar="double"),
    make_option(c("--percision_beta"), type="integer", default=opt$percision_beta,
                help="Rounding percision for beta values in calls output files [default= %default]", metavar="double"),
    make_option(c("--percision_pval"), type="integer", default=opt$percision_pval,
                help="Rounding percision for detection p-values in calls output files [default= %default]", metavar="double"),
    
    # Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
    # Plotting Options::
    make_option(c("--plotSset"), action="store_true", default=opt$plotSset,
                help="Boolean variable to plot intensity distributions for sset [default= %default]", metavar="boolean"),
    make_option(c("--plotCalls"), action="store_true", default=opt$plotCalls,
                help="Boolean variable to plot detection p-values and beta distributions [default= %default]", metavar="boolean"),
    make_option(c("--plotAuto"), action="store_true", default=opt$plotAuto,
                help="Boolean variable to plot Auto-Detection Matricies (Pval/Beta) file [default= %default]", metavar="boolean"),
    
    make_option(c("--plotFormat"), type="character", default=opt$plotFormat, 
                help="Plotting output format [default= %default]", metavar="character"),
    make_option(c("--dpi"), type="double", default=opt$dpi, 
                help="DPI for plot images Plotting [default= %default]", metavar="double"),
    
    make_option(c("--plotMax"), type="double", default=opt$plotMax, 
                help="Max Sample Display Count for Plotting [default= %default]", metavar="double"),
    make_option(c("--plotSub"), type="double", default=opt$plotSub, 
                help="Sub Sample Display Count for Plotting [default= %default]", metavar="double"),
    
    # make_option(c("--opt_csv"), type="character", default=opt$opt_csv, 
    #             help="Unused variable opt_csv [default= %default]", metavar="character"),
    # make_option(c("--par_csv"), type="character", default=opt$par_csv, 
    #             help="Unused variable par_csv [default= %default]", metavar="character"),
    # make_option(c("--time_csv"), type="character", default=opt$time_csv, 
    #             help="Unused variable time_csv [default= %default]", metavar="character"),
    
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
opt_reqs <- c('outDir','Rscript','verbose')

par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) 
  stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))

for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', 
                         full.names=TRUE, recursive=TRUE)) base::source(sfile)
if (opt$verbose>=0)
  cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form ",
                 "General Source={par$gen_src_dir}!{RET}{RET}") )

# TestCase::
if (is.null(opt$idatsDir) || !dir.exists(opt$idatsDir)) {
  opt$runName  <- "TestCase"
  opt$idatsDir <- file.path(par$datDir,"idats_TestCase")
  
  if (opt$verbose>0) {
    cat(glue::glue("[{par$prgmTag}]: idatsDir does not exist or ",
                   "is missing will use TestData={opt$idatsDir}!!!{RET}"))
    cat(glue::glue("[{par$prgmTag}]: Overriding runName={opt$runName}.{RET}{RET}"))
  }
  
  # stop(glue::glue("{RET}[{par$prgmTag}]: idatsDir={opt$idatsDir} does not exist!!!{RET}{RET}"))
}

opt <- program_init(name=par$prgmTag,
                    opts=opt, opt_reqs=opt_reqs, 
                    pars=par, par_reqs=par_reqs,
                    libs=TRUE,rcpp=FALSE,
                    verbose=opt$verbose,vt=3,tc=0,tt=NULL)

opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
# cat(glue::glue("[{par$prgmTag}]: Output Directory (TOP)={opt$outDir}...{RET}"))

pval_vec     <- splitStrToVec(opt$pval)
min_pval_vec <- splitStrToVec(opt$minPval)
min_perc_vec <- splitStrToVec(opt$minPerc)
workflow_vec <- splitStrToVec(opt$workflow)

# Remove "r/raw" and force "raw" to be first::
workflow_vec <- c("raw",workflow_vec[!workflow_vec %in% c("r","raw")])

cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             Select Chips from idats and/or Target Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
chipPrefixes <- NULL
chipPrefixes <- sesame::searchIDATprefixes(opt$idatsDir)
sampleCounts <- chipPrefixes %>% names() %>% length()

if (is.null(chipPrefixes) || length(chipPrefixes)==0)
  stop(glue::glue("{RET}[{par$prgmTag}]: chipPrefixes is null or length=0!!!{RET}{RET}"))

cat(glue::glue("[{par$prgmTag}]: Found sample counts={sampleCounts}!{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                  Main::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
if (opt$cluster) {
  cat(glue::glue("[{par$prgmTag}]: Launching Chips in Cluster Mode! isSingle={opt$single}"),"\n", sep='')
  
  par$lan_exe <- ''
  par$isLinux <- FALSE
  if (!is.null(par$lixDir1) && length(par$lixDir1)>0 && dir.exists(par$lixDir1)) {
    par$isLinux <- TRUE
    par$lan_exe <- 'qsub -cwd -pe threaded 16 -l excl=true -N'
    if (dir.exists(par$macDir)) stop(glue::glue("[{par$prgmTag}]: Linux/Mac directories exist???{RET}{RET}"))
  }
  chip_list <- prefixesToChipTib(chipPrefixes) %>% split(.$barcode)
  chip_cnts <- chip_list %>% names() %>% length()
  
  cat(glue::glue("[{par$prgmTag}]:{TAB}Cluster Mode; Chip counts={chip_cnts}!{RET}"))
  
  par$shellDir <- file.path(opt$outDir, 'shells')
  if (!dir.exists(par$shellDir)) dir.create(par$shellDir, recursive=TRUE)
  cat(glue::glue("[{par$prgmTag}]:{TAB}Cluster Mode; shellDir={par$shellDir}.{RET}"))
  
  for (chipName in names(chip_list)) { # break }
    runShell <- file.path(par$shellDir, paste0('run_',par$prgmTag,'_',chipName,'.sh'))
    lanShell <- file.path(par$shellDir, paste0('lan_',par$prgmTag,'_',chipName,'.sh'))
    cat(glue::glue("[{par$prgmTag}]:{TAB}{TAB}Cluster Mode; runShell={runShell}.{RET}"))
    
    # Remove Cluster Option
    cmd_bool <- opt %>% bind_rows() %>% gather("Options", "Value") %>% 
      dplyr::filter(Value=='TRUE') %>% 
      dplyr::filter(Options!='cluster') %>% 
      dplyr::select(Options) %>% 
      dplyr::mutate(Options=paste0('--',Options)) %>% 
      dplyr::pull() %>% paste(collapse=" ")
    
    # Add ChipName to idat Directory
    cmd_strs <- opt %>% bind_rows() %>% gather("Options", "Value") %>%
      dplyr::filter(Value!='TRUE' & Value!='FALSE' & Options!='prgmPath') %>%
      dplyr::mutate(Value=case_when(Options=='idatsDir' ~ paste0(Value,'/',chipName), TRUE ~ Value),
                    Value=case_when(Options=='outDir' ~ paste0(Value,'/',chipName), TRUE ~ Value) ) %>%
      dplyr::mutate(Options=paste0('--',Options,'=',Value)) %>%
      dplyr::select(Options) %>% 
      dplyr::pull() %>% paste(collapse=" ")
    
    cmd_full <- paste(opt$Rscript, par$exePath, cmd_bool, cmd_strs,"\n", sep=' ')
    cat(glue::glue("[{par$prgmTag}]:{TAB}{TAB}Writing cmd_full={cmd_full}.{RET}"))
    # readr::write_file(x=cmd_full, file=runShell)
    readr::write_lines(x=cmd_full, file=runShell)
    Sys.chmod(runShell, mode="0777")
    
    # Add cluster execute if avialbel (i.e. linux)
    cmd_lan <- paste0(runShell,RET, sep='')
    if (par$isLinux)
      cmd_lan <- paste(par$lan_exe, paste('imWH',chipName, sep='_'), runShell,RET, sep=' ')
    readr::write_file(cmd_lan, file=lanShell)
    Sys.chmod(lanShell, mode="0777")
    
    # Launch Script
    cat(glue::glue("[{par$prgmTag}]:{TAB}{TAB}Cluster Mode; Launching chip={chipName}, shell={lanShell}.{RET}"))
    system(lanShell)
    
    if (opt$single) break
  }
  cat(glue::glue("[{par$prgmTag}]:{TAB}Cluster Mode; Done.{RET}"))
  
} else {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Load Manifest(s)::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  cat(glue::glue("[{par$prgmTag}]: Launching Samples in Linear Mode! isSingle={opt$single}"),"\n", sep='')

  pTracker <- NULL
  # pTracker <- timeTracker$new(verbose=opt$verbose)
  
  mans <- NULL
  opt$manDir <- file.path(par$datDir, 'manifest',opt$manDirName)
  mans <- getManifestList(path=opt$manifestPath, platform=opt$platform, manifest=opt$manifest, 
                          dir=opt$manDir, verbose=opt$verbose, tt=pTracker)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Load Auto Detection::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  auto_sam_tib <- NULL
  if (opt$auto_detect) {
    if (is.null(opt$auto_sam_csv))
      opt$auto_sam_csv <- file.path(par$datDir, 'ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo_Mean-COVIC-280-NP-ind_negs-0.02.csv.gz')

    cat(glue::glue("[{par$prgmTag}]:{TAB} Loading auto_sam_csv={opt$auto_sam_csv}...{RET}"))
    auto_sam_tib <- suppressMessages(suppressWarnings(readr::read_csv(opt$auto_sam_csv) ))
    cat(glue::glue("[{par$prgmTag}]:{TAB} Done.{RET}"))
  }
  
  # TBD::
  #   - mask pass_perc
  #   - R2/dB returned by Infinium Design
  #
  mask_cpg_vec <- NULL
  mask_cpg_csv <- file.path(par$datDir, 'manifest/mask/sesame-general-mask.cpg.csv.gz')
  if (opt$mask_general) {
    mask_cpg_vec <- suppressMessages(suppressWarnings( 
      readr::read_csv(mask_cpg_csv) )) %>% dplyr::pull(Probe_ID)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Chip::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  chipTimes <- NULL
  sample_cnt <- length(chipPrefixes)
  prefixe_names <- names(chipPrefixes)

  try_str <- ''
  if (opt$parallel) {
    par$funcTag <- 'sesamizeSingleSample-Parallel'
    par$retData <- FALSE
    
    if (opt$verbose>=1) 
      cat(glue::glue("[{par$prgmTag}]: parallelFunc={par$funcTag}: samples={sample_cnt}; ",
                     "num_cores={num_cores}, num_workers={num_workers}, Starting...{RET}"))
    
    par$retData <- FALSE
    # chipTimes <- foreach (prefix=prefixe_names, .inorder=T, .final = function(x) setNames(x, prefixe_names)) %dopar% {
    chipTimes <- foreach (prefix=prefixe_names, .combine = rbind) %dopar% {

      rdat <- NULL
      rdat <- sesamizeSingleSample(prefix=chipPrefixes[[prefix]],
                                   man=mans, ref=auto_sam_tib, opts=opt, defs=def,
                                   mask=mask_cpg_vec,
                                   
                                   pvals=pval_vec,
                                   min_pvals=min_pval_vec,
                                   min_percs=min_perc_vec,
                                   workflows=workflow_vec,
                                   
                                   retData=par$retData,
                                   verbose=opt$verbose, vt=3,tc=1)
      rdat
    }
    cat(glue::glue("[{par$prgmTag}] parallelFunc={par$funcTag}: Done.{RET}{RET}"))
  } else {
    par$funcTag <- 'sesamizeSingleSample-Linear'
    
    cat(glue::glue("[{par$prgmTag}]: linearFunc={par$funcTag}: samples={sample_cnt}.{RET}"))
    
    rdat <- NULL
    for (prefix in prefixe_names) {
      
      # Testing code only::
      #
      if (FALSE && par$runMode=='RStudio' && !is.null(par$expSampNum)) {
        prefix <- par$expSampNum
        opt$single <- TRUE
        opt$fresh  <- TRUE
        opt$fresh  <- FALSE
        
        # workflow_vec <- c('raw')
        # opt$make_pred <- FALSE
      }
      cat(glue::glue("[{par$prgmTag}]: linearFunc={par$funcTag}: Starting; prefix={prefix}...{RET}"))

      rdat <- NULL
      rdat <- sesamizeSingleSample(prefix=chipPrefixes[[prefix]],
                                   man=mans, ref=auto_sam_tib, opts=opt, defs=def,
                                   mask=mask_cpg_vec,
                                   
                                   pvals=pval_vec,
                                   min_pvals=min_pval_vec,
                                   min_percs=min_perc_vec,
                                   workflows=workflow_vec,
                                   
                                   retData=par$retData,
                                   verbose=opt$verbose, vt=3,tc=1)
      
      # rdat1 <- rdat
      # rdat7 <- rdat
      #
      # rdat$ssheet_tib %>% dplyr::select(dplyr::contains("_Requeue_"),
      #                                   dplyr::contains("_pass_perc_") ) %>% as.data.frame()
      # requeueFlag(tib=rdat$cur_list$sums_dat, name = "BLA", csv = "/Users/bretbarnes/Documents/tmp/tmp.requeue.txt", verbose=40)
      #
      # auto_sam_tib %>% dplyr::select(Probe_ID, rdat$ssheet_tib$AutoSample_dB_Key_1)
      # dplyr::select(auto_sam_tib, Probe_ID, rdat$ssheet_tib$AutoSample_dB_Key_1)
      #
      # Generate plot(dB<=0.2) + plot(dB>0.2)
      # Calculate the angle between ablines above
      #
      
      if (FALSE) {
        
        sesame::formatVCF(sset = rdat$open_sset_dat, vcf = file.path(opt$outDir,'tmp.vcf'))
        
        dtmp_tib <- getSsheetDataTab(tib = rdat$ssheet_tib, minOobPval = 0.1, 
                                     minOobPerc = 90, minNegPval = 0.02, 
                                     minNegPerc = 98, minDb = 0.2, verbose = 10)
        print(dtmp_tib, n=base::nrow(dtmp_tib))
        
        # Sentrix_Name
        # Failed_QC
        # Min_Pass_Perc
        rdat$ssheet_tib$Sentrix_Name
        rdat$ssheet_tib$cg_Failed_QC_basic_0
        rdat$ssheet_tib$cg_pass_perc_basic_0
        
        open_req_tib <- requeueFlagOpenSesame(tib=rdat$ssheet_tib, name=rdat$ssheet_tib$Sentrix_Name,
                                              verbose = 10)
        
        #
        #
        # TBD:: Two things::
        #   - Capture r2_val/db_val described below
        #   - Report numerator and denominator in pval perc passing calculations
        #   - Add run_name to Auto Sample Sheet
        #
        
        #
        # rdat$ssheet_tib %>% dplyr::select(dplyr::contains("basic"))
        #
        
        dat_dir <- '/Users/bretbarnes/Documents/scratch/mlk/swifthoof_main/CNTL-Samples_VendA_10092020'
        dat_lst <- list.files(dat_dir, pattern='_AutoSampleSheet.csv.gz$', full.names=TRUE)
        ssh_tib <- lapply(dat_lst, readr::read_csv) %>% dplyr::bind_rows()
        
        ssh_tib %>% 
          dplyr::select(Sentrix_Name,AutoSample_dB_1_Key_1,dplyr::contains('pass_perc')) %>% 
          dplyr::select(Sentrix_Name,AutoSample_dB_1_Key_1,dplyr::starts_with('cg')) %>%
          dplyr::filter(!stringr::str_starts(AutoSample_dB_1_Key_1,'T')) %>%
          dplyr::arrange(AutoSample_dB_1_Key_1,cg_pass_perc_basic_0) %>%
          dplyr::select(!dplyr::contains('_PnegEcdf_')) %>%
          as.data.frame()
        
        #
        # Capture both of these values to demonstrate the value of normalization via Sesame::
        #
        r2_val <- rdat$open_beta_tib %>% 
          tibble::enframe(name="Probe_ID", value="betas") %>%
          dplyr::inner_join(rdat$org_list$beta$beta_dat, 
                            by="Probe_ID", suffix=c("_ref","_can")) %>%
          tibble::column_to_rownames(var="Probe_ID") %>% 
          as.matrix() %>% cor() %>% as_tibble() %>% 
          head(n=1) %>% dplyr::pull(2)
        
        dB_val <- rdat$open_beta_tib %>% 
          tibble::enframe(name="Probe_ID", value="betas") %>%
          dplyr::inner_join(rdat$org_list$beta$beta_dat, 
                            by="Probe_ID", suffix=c("_ref","_can")) %>%
          dplyr::filter(stringr::str_starts(Probe_ID,'cg')) %>% 
          dplyr::mutate(dB=base::abs(betas_ref-betas_can)) %>% 
          dplyr::summarise(pass_perc=cntPer_lte(dB,opt$minDeltaBeta)) %>%
          head(n=1) %>% dplyr::pull(1)
        

        
        
        
        dB_tib <- 
          rdat$open_beta_tib %>% 
          tibble::enframe(name="Probe_ID", value="betas") %>%
          dplyr::inner_join(rdat$org_list$beta$beta_dat, 
                            by="Probe_ID", suffix=c("_ref","_can")) %>%
          dplyr::filter(stringr::str_starts(Probe_ID,'cg')) %>% 
          dplyr::mutate(dB=betas_ref-betas_can)
        
          # dplyr::mutate(dB=base::abs(betas_ref-betas_can))
        
        dB_tib %>% dplyr::summarise(pass_perc=cntPer_lte(dB,opt$minDeltaBeta), pcount=count(dB<0.2),tot_cnt=n(),per2=pcount/tot_cnt)
        
        
        
        # Generate open sset
        #  - Get open_call_tib
        #
        # Get manifest = null stats + requeue
        # Get manifest = EPIC stats
        #

        if (stringr::str_starts(rdat$ssheet_tib$detect_manifest, 'EPIC') ) {
          open_sset_dat <- sesame::openSesame(x=rdat$prefix, what = 'sigset')
          
          open_sum2_tib <- ssetToPassPercSsheet(sset=open_sset_dat, 
                                                man=rdat$sman, min=0.1, per=90,
                                                verbose=10)
          
          open_sum1_tib <- ssetToPassPercSsheet(sset=open_sset_dat,
                                                min=0.1, per=90,
                                                verbose=10)
          
          print(open_sum1_tib)
          print(open_sum2_tib)
          
        }
        
        open_dat %>% dplyr::mutate(Workflow_idx=0) %>% requeueFlag(idx=0)
        
        # Open Sesame Comparison::
        open_sset <- sesame::openSesame(x=rdat$prefix, what = 'sigset')
        open_pval_tib <- open_sset@pval %>% tibble::enframe(name="Probe_ID", value="raw_pvals_pOOBAH")
        
        
        # Add the full filter without Inf I/II spliting...
        open_pval_tib %>% 
          dplyr::filter(stringr::str_starts(Probe_ID,'cg')) %>% 
          dplyr::summarise(pass_perc=cntPer_lte(raw_pvals_pOOBAH, min=0.1), 
                           total_cnt=n(), 
                           pass_cnt=count(raw_pvals_pOOBAH<0.1), 
                           pass_perc2=round(pass_cnt/total_cnt, 3)) %>%
          purrr::set_names(paste('cg',names(.),'basic_0', sep='_'))
        
        # TBD:: Manifest loading function needs to retain original IDs
        #
        open_pval_tib %>% dplyr::left_join(rdat$sman, by="Probe_ID") %>% 
          dplyr::group_by(Probe_Type,Probe_Design) %>% 
          dplyr::summarise(pass_perc=cntPer_lte(raw_pvals_pOOBAH, min=0.1), 
                           Total_Count=n(), Pass_Count=count(raw_pvals_pOOBAH<0.1), 
                           Pass_Perc2=round(Pass_Count/Total_Count, 3))
        
        tmp_outDir <- '/Users/bretbarnes/Documents/tmp'
        open_sum_dat <- ssetToSummary(sset = open_sset, man = rdat$sman, idx=2, 
                                      workflow='open', name='open', outDir=tmp_outDir,
                                      min_percs = c(98), min_pvals = c(0.1), pvals = c('pOOBAH'),
                                      makek_pred=FALSE, fresh=TRUE,
                                      verbose=10)
        
        
        
        val_rdat <- rdat
        
        # v.4.6 values::
        # val_rdat$ssheet_tib$cg_1_pvals_pOOBAH_pass_perc_1 # 87.903
        # val_rdat$ssheet_tib$cg_2_pvals_pOOBAH_pass_perc_1 # 83.219
        
        val_rdat$ssheet_tib$cg_1_pvals_pOOBAH_pass_perc_1
        val_rdat$ssheet_tib$cg_2_pvals_pOOBAH_pass_perc_1
        
        cur_sentrix_name <- rdat$ssheet_tib$Sentrix_Name
        
        old_dir <- '/Users/bretbarnes/Documents/tools/bk/docker-repo/Infinium_Methylation_Workhorse.v.4.6/scratch/swifthoof/swifthoof_main'
        old_dir <- "/Users/bretbarnes/Documents/scratch/RStudio/swifthoof_main/CNTL-Samples_VendA_10092020"
        call_old_csv <- file.path(old_dir, paste(cur_sentrix_name, 'EPIC_B4_raw.call.dat.csv.gz', sep="_"))
        
        call_old_tib <- readr::read_csv(call_old_csv) %>%
          dplyr::left_join(rdat$sman, by="Probe_ID") %>%
          dplyr::arrange(Probe_ID)
        
        # Calculate Pass Percent With/Without NA's
        call_old_tib %>% 
          # dplyr::filter(!is.na(betas)) %>% 
          dplyr::group_by(Probe_Type,Probe_Design) %>%
          dplyr::summarise(Total_Count=n(), 
                           Pass_Perc_Poob=cntPer_lte(pvals_pOOBAH,0.1),
                           Pass_Perc_Negs=cntPer_lte(pvals_PnegEcdf,0.05),
                           .groups="drop")
        
        rdat$new_sset %>% ssetToTib(source="pvals", name="pOOBAH") %>% 
          dplyr::left_join(rdat$sman, by="Probe_ID") %>% 
          dplyr::group_by(Probe_Type,Probe_Design) %>% 
          dplyr::summarise(Total_Count=n(), Pass_Cnt=count(pvals_pOOBAH<=0.1), Pass_Per=Pass_Cnt/Total_Count, Pass_Perc=cntPer_lte(pvals_pOOBAH,0.1), .groups="drop")
        
        
        rdat$cur_list$call_dat %>% dplyr::select(1,2) %>% ssetTibToSummary()
        
        
        
        # Basic Calculation::
        cg_inf1_cnt <- call_old_tib %>% dplyr::filter(Probe_Type=='cg' & Probe_Design==1) %>% base::nrow()
        cg_inf2_cnt <- call_old_tib %>% dplyr::filter(Probe_Type=='cg' & Probe_Design==2) %>% base::nrow()
        
        cg_pas1_cnt <- call_old_tib %>% dplyr::filter(Probe_Type=='cg' & Probe_Design==1 & pvals_pOOBAH <= 0.1) %>% base::nrow()
        cg_pas2_cnt <- call_old_tib %>% dplyr::filter(Probe_Type=='cg' & Probe_Design==2 & pvals_pOOBAH <= 0.1) %>% base::nrow()
        
        cg_pas1_cnt / cg_inf1_cnt
        cg_pas2_cnt / cg_inf2_cnt
        
        #
        # Old analysis
        #
        rdat1 <- rdat
        dB1_tib <- dplyr::inner_join( 
          purrr::set_names(rdat1$org_list$call_dat, c("Probe_ID", "can_poob", "can_negs", "can_beta")),
          dplyr::select(auto_sam_tib, Probe_ID, rdat1$ssheet_tib$AutoSample_dB_Key_1) %>% 
            purrr::set_names(c("Probe_ID","ref_beta")), 
          by="Probe_ID") %>% 
          dplyr::mutate(
            del_beta=can_beta-ref_beta,
            abs_beta=base::abs(can_beta-ref_beta)) %>%
          dplyr::inner_join(rdat$sman, by="Probe_ID") # %>% dplyr::mutate(DESIGN=as.factor(DESIGN))
        
        ind1_tib <- dplyr::inner_join(
          rdat1$cur_list$call_dat %>% dplyr::select(Probe_ID, dplyr::starts_with('ind_')) %>%
          purrr::set_names(c("Probe_ID", "can_poob", "can_negs", "can_beta")),
          
          dplyr::select(auto_sam_tib, Probe_ID, rdat1$ssheet_tib$AutoSample_dB_Key_1) %>% 
            purrr::set_names(c("Probe_ID","ref_beta")), 
          by="Probe_ID") %>% 
          dplyr::mutate(
            del_beta=can_beta-ref_beta,
            abs_beta=base::abs(can_beta-ref_beta)) %>%
          dplyr::inner_join(rdat$sman, by="Probe_ID")
        
        dB1_tib %>%
          ggplot2::ggplot(aes(x=can_poob, y=del_beta, color = DESIGN)) + 
          # ggplot2::geom_point() + 
          # ggplot2::geom_abline() +
          ggplot2::geom_density2d()
        
        ind1_tib %>%
          dplyr::filter(Probe_Type=='cg') %>%
          ggplot2::ggplot(aes(x=can_poob, y=del_beta, color = DESIGN)) + 
          ggplot2::geom_point() +
          ggplot2::geom_abline() +
          # ggplot2::geom_density2d() +
          ggplot2::facet_grid(rows = "DESIGN")
        
        
        
        #
        # Plot histogram of delta-beta split by detection p-value
        #
        dB1_tib %>% 
          # dplyr::filter(can_poob<0.05) %>%
          ggplot2::ggplot(group = DESIGN) + 
          ggplot2::geom_density(aes(x=del_beta))
        
        ggplot2::ggplot(data=dB1_tib) + 
          ggplot2::geom_density(aes(x=del_beta))
          

        
        
        dB7_tib <- dplyr::inner_join( 
          purrr::set_names(rdat7$org_list$call_dat, c("Probe_ID", "can_poob", "can_negs", "can_beta")),
          dplyr::select(auto_sam_tib, Probe_ID, rdat7$ssheet_tib$AutoSample_dB_Key_1) %>% 
            purrr::set_names(c("Probe_ID","ref_beta")), 
          by="Probe_ID") %>% 
          dplyr::mutate(del_beta=base::abs(can_beta-ref_beta))
        
        ggplot2::ggplot(data=dB7_tib, aes(x=can_poob, y=del_beta)) + 
          # ggplot2::geom_point() + 
          ggplot2::geom_density2d() +
          ggplot2::geom_abline()
      }

      cat(glue::glue("[{par$prgmTag}]: linearFunc={par$funcTag}: try_str={try_str}. Done.{RET}{RET}"))
      if (opt$single) break
    }
    cat(glue::glue("[{par$prgmTag}] parallelFunc={par$funcTag}: Done.{RET}{RET}"))
  }
  
  opt_csv <- file.path(opt$outDir, paste(par$prgmTag,'program-options.csv', sep='.') )
  par_csv <- file.path(opt$outDir, paste(par$prgmTag,'program-parameters.csv', sep='.') )
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value") %>% dplyr::arrange(Option)
  par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value") %>% dplyr::arrange(Params)
  
  readr::write_csv(opt_tib, opt_csv)
  readr::write_csv(par_tib, par_csv)
  
  # tim_csv <- file.path(opt$outDir, paste(par$prgmTag,'time-tracker.csv.gz', sep='.') )
  # tim_tib <- pTracker$time %>% dplyr::mutate_if(is.numeric, list(round), 4)
  # readr::write_csv(tim_tib, tim_csv)
  
  # rdat$cur_list$call_dat %>% dplyr::filter(stringr::str_starts(Probe_ID, 'ch'))
  
  # rdat$cur_list$call_dat %>% dplyr::filter(stringr::str_starts(Probe_ID, 'cg')) %>% dplyr::filter(raw_pvals_pOOBAH<=0.1) %>% base::nrow()
  # rdat$cur_list$call_dat %>% dplyr::filter(stringr::str_starts(Probe_ID, 'cg')) %>% base::nrow()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
