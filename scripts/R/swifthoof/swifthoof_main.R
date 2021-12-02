
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Load Core Packages::
suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("sesame",quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("minfi",quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse",quietly=TRUE) ))

# Load Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel",quietly=TRUE) ))

# Load Performance Packages
suppressWarnings(suppressPackageStartupMessages( base::require("profmem",quietly=TRUE) ))

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
par$prgmDir <- 'swifthoof'
par$prgmTag <- paste(par$prgmDir,'main', sep='_')
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

par$retData <- FALSE
par$manDir  <- NULL

# Executables::
opt$Rscript <- NULL

# Directories::
opt$outDir <- NULL
opt$datDir <- NULL

# Platform/Method Options::
opt$manifest <- NULL
opt$platform <- NULL
opt$version  <- NULL
opt$percent  <- NULL

# Optional Files::
opt$subManifest  <- FALSE
opt$auto_sam_csv <- NULL

# Run Options::
opt$fresh        <- FALSE
opt$buildSubDir  <- FALSE
opt$auto_detect  <- FALSE

opt$workflow   <- NULL
opt$manDirPath <- NULL
opt$manDirName <- 'core'
opt$forcedPlat <- NULL
opt$man_suffix <- ".manifest.sesame-base.cpg-sorted.csv.gz"

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

opt$write_call  <- TRUE
opt$write_csum  <- FALSE

opt$write_snps  <- TRUE
opt$write_auto  <- FALSE

opt$mask_general <- FALSE

# Threshold Options::
opt$pval    <- "pOOBAH,PnegEcdf"
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
opt$trackTime    <- FALSE

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
  
  opt$runName  <- NULL
  opt$platform <- NULL
  opt$version  <- NULL
  
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
  
  opt$auto_detect <- TRUE
  
  par$retData  <- TRUE
  
  opt$single   <- FALSE
  opt$single   <- TRUE
  
  opt$parallel <- TRUE
  opt$parallel <- FALSE
  
  opt$cluster  <- TRUE
  opt$cluster  <- FALSE
  
  opt$save_idat <- TRUE
  opt$load_idat <- TRUE
  
  opt$save_sset <- TRUE
  opt$load_sset <- TRUE
  
  opt$manDirPath <- NULL
  opt$manDirName <- 'base'
  opt$manDirName <- 'core'
  
  opt$verbose  <- 3
  
  par$local_runType <- 'VA-MVP-Akesogen_Batch3'
  par$local_runType <- 'VA-MVP-IBX_Batch3'
  
  par$local_runType <- 'CORE'
  par$local_runType <- 'EXCBR'
  par$local_runType <- 'GRCm38'
  par$local_runType <- 'COVID'
  par$local_runType <- 'GRCm38'
  par$local_runType <- 'DELTA-8x1-EPIC-Core'
  par$local_runType <- 'DKFZ'
  par$local_runType <- 'qcMVP'
  par$local_runType <- 'COVIC'
  par$local_runType <- "EPIC-8x1-EM-Sample-Prep"
  par$local_runType <- 'qcMVP2'
  par$local_runType <- 'NA12878'
  par$local_runType <- 'Chicago-Ober-Custom'
  par$local_runType <- 'NZT'
  par$local_runType <- 'COVIC-NZT_23092020'
  
  opt$fresh <- TRUE
  
  opt$auto_sam_csv <- 
    file.path(par$datDir, 'ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo_Mean-COVIC-280-NP-ind_negs-0.02.csv.gz')
  
  opt$auto_detect <- TRUE
  opt$workflow    <- "ind"
  opt$manDirName  <- 'core'
  
  # opt$write_snps  <- TRUE
  opt$write_snps  <- FALSE
  
  opt$single   <- TRUE
  opt$single   <- FALSE
  
  opt$parallel <- FALSE
  opt$parallel <- TRUE
  
  opt$fresh <- TRUE
  opt$trackTime <- TRUE
  
  opt$runName  <- par$local_runType
  
  if (FALSE) {
  } else if (par$local_runType=='NZT' || par$local_runType=="COVIC-NZT_23092020") {
    
    opt$single   <- TRUE
    opt$single   <- FALSE
    opt$parallel <- FALSE
    opt$parallel <- TRUE
    opt$fresh    <- TRUE
    
    # For sub manifest testing::
    opt$platform   <- "NZT"
    opt$version    <- "C2"
    opt$version    <- "B1"
    
    opt$manDirName  <- 'base'
    # opt$manDirPath <- file.path(par$topDir, "data/manifests/methylation/Sesame/NZT")
    
    # opt$auto_sam_csv <- "/Users/bretbarnes/Documents/data/CustomContent/Chicago-Ober-Custom/AutoDetect/v1/AutoSampleDetection_Chicago-Ober-Custom-v1.csv.gz"
    
  } else if (par$local_runType=='Chicago-Ober-Custom') {
    
    opt$single   <- TRUE
    opt$parallel <- FALSE
    opt$fresh    <- TRUE
    
    opt$workflow    <- "d,ind"
    
    # For sub manifest testing::
    opt$platform   <- "Chicago"
    opt$version    <- "S38"
    opt$version    <- "S39"
    opt$manDirPath <- file.path(par$topDir, "data/manifests/methylation/Chicago-Ober-Custom")
    opt$auto_sam_csv <- "/Users/bretbarnes/Documents/data/CustomContent/Chicago-Ober-Custom/AutoDetect/v1/AutoSampleDetection_Chicago-Ober-Custom-v1.csv.gz"
    
  } else if (par$local_runType=='EPIC-8x1-EM-Sample-Prep') {
    
  } else if (par$local_runType=='NA12878') {
    
    opt$single   <- TRUE
    opt$parallel <- FALSE
    opt$fresh    <- TRUE

    # For sub manifest testing::
    opt$platform   <- "Rand1"
    opt$version    <- 20
    opt$percent    <- 10
    # opt$version    <- "S20"
    
    # opt$manDirPath <- file.path(par$topDir, "data/manifests/methylation/McMaster10Kselection", opt$platform)
    # opt$manDirPath <- file.path(par$topDir, "scratch/RStudio/manifest_subset_noob/EPIC-noob-BP4/man", opt$platform)
    opt$manDirPath <- file.path(par$topDir, "scratch/RStudio/manifest_subset_noob/EPIC-noob-BP4/man")
    
    opt$forcedPlat <- "EPIC"
    
  } else if (par$local_runType=='qcMVP2') {
    
    opt$runName  <- 'AKE-Zymogen'
    opt$runName  <- 'AKE-EPIDX'
    
    opt$runName  <- 'IBX-Zymogen'
    opt$runName  <- 'IBX-EPIDX'
    
  } else if (par$local_runType=='COVID') {
    par$expChipNum <- '204756130014'
    
  } else if (par$local_runType=='COVIC') {
    opt$runName  <- 'COVIC-Set1-15052020'
    par$expChipNum <- '204500250013'
    
    opt$manDirName  <- 'core'
    opt$manDirName  <- 'covic'
    
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
    
    par$expChipNum <- "203962710025"
    par$expSampNum <- "203962710025_R08C01"
    
    # Fred's Failure Case::
    par$expChipNum <- "203962710025"
    par$expSampNum <- "203962710025_R01C01"
    
    opt$dpi <- 72
  } else if (par$local_runType=='DKFZ') {
    
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
    
  } else if (par$local_runType=='EXCBR') {
    opt$runName  <- 'Excalibur-Old-1609202'
    par$expChipNum <- '204076530053'
    par$expChipNum <- '204076530110'
    
    opt$runName  <- 'Excalibur-New-1609202'
    par$expChipNum <- '202915460071'
    
  } else if (par$local_runType=='DELTA-8x1-EPIC-Core') {
    par$expChipNum <- '203319730022'
    par$expSampNum <- '203319730022_R07C01'
    
    opt$auto_detect <- TRUE
    
    opt$plotSset  <- TRUE
    opt$plotCalls <- TRUE
    opt$plotAuto  <- TRUE
    
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unrecognized local_runType={par$local_runType}.{RET}{RET}"))
  }
  
  opt$datDir <- file.path(locIdatDir, paste('idats',opt$runName, sep='_') )
  if (!is.null(par$expChipNum)) 
    opt$datDir <- file.path(locIdatDir, paste('idats',opt$runName, sep='_'),  par$expChipNum)
  
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
    make_option(c("-d","--datDir"), type="character", default=opt$datDir, 
                help="List of idats directory(s), commas seperated [default= %default]", metavar="character"),
    
    # Platform/Method Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--manifest"), type="character", default=opt$manifest,
                help="Path to manfifest (CSV) otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Forced platform [EPIC, 450k, 27k, NZT] otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Forced version [B1, B2, B4, etc.] otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--percent"), type="character", default=opt$percent, 
                help="Forced percent [B1, B2, B4, etc.] otherwise auto-detect [default= %default]", metavar="character"),
    
    # Optional Files::
    make_option(c("--subManifest"), action="store_true", default=opt$subManifest,
                help="Boolean variable to use subset manifest instead of subset. [default= %default]", metavar="boolean"),
    make_option(c("--auto_sam_csv"), type="character", default=opt$auto_sam_csv,
                help="Path to auto detect beta values (CSV) [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--fresh"), action="store_true", default=opt$fresh,
                help="Boolean variable to build fresh version of database files [default= %default]", metavar="boolean"),
    make_option(c("--buildSubDir"), action="store_true", default=opt$buildSubDir,
                help="Boolean variable to build subdirectories based on Chip/BeadPool (for R&D purposes) [default= %default]", metavar="boolean"),
    make_option(c("--auto_detect"), action="store_true", default=opt$auto_detect,
                help="Boolean variable to auto detect reference samples. Must provide reference samples. [default= %default]", metavar="boolean"),
    
    make_option(c("--workflow"), type="character", default=opt$workflow,
                help="Workflows to run i.e. order of operations comma seperated [ raw,ind,ndi,din, etc. ] [default= %default]", metavar="character"),
    make_option(c("--manDirPath"), type="character", default=opt$manDirPath,
                help="Manifest directory path otherwise use default dat/manifest directory [default= %default]", metavar="character"),
    make_option(c("--manDirName"), type="character", default=opt$manDirName,
                help="Manifest directory name [default= %default]", metavar="character"),
    make_option(c("--forcedPlat"), type="character", default=opt$forcedPlat,
                help="Forced Manifest name to use in Sesame calculations [default= %default]", metavar="character"),
    make_option(c("--man_suffix"), type="character", default=opt$man_suffix,
                help="Manifest suffix search name. Default is set to predefined Sesame suffix. [default= %default]", metavar="character"),
    
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
    make_option(c("--trackTime"), action="store_true", default=opt$trackTime,
                help="Boolean variable tack run times [default= %default]", metavar="boolean"),
    
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
if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form ",
                 "General Source={par$gen_src_dir}!{RET}{RET}") )

# TBD:: Make sure we load this::
#  library(Rcpp)
#  source("/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/scripts/R/Rcpp/cpgLociVariation.cpp")
#

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

if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             Select Chips from idats and/or Target Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
chipPrefixes <- NULL
chipPrefixes <- sesame::searchIDATprefixes(opt$datDir)
sampleCounts <- chipPrefixes %>% names() %>% length()

if (is.null(chipPrefixes) || length(chipPrefixes)==0)
  stop(glue::glue("{RET}[{par$prgmTag}]: chipPrefixes is null or length=0!!!{RET}{RET}"))

if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Found sample counts={sampleCounts}!{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                  Main::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
if (opt$cluster) {
  if (opt$verbose>0)
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
  
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]:{TAB}Cluster Mode; Chip counts={chip_cnts}!{RET}"))
  
  par$shellDir <- file.path(opt$outDir, 'shells')
  if (!dir.exists(par$shellDir)) dir.create(par$shellDir, recursive=TRUE)
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]:{TAB}Cluster Mode; shellDir={par$shellDir}.{RET}"))
  
  for (chipName in names(chip_list)) { # break }
    runShell <- file.path(par$shellDir, paste0('run_',par$prgmTag,'_',chipName,'.sh'))
    lanShell <- file.path(par$shellDir, paste0('lan_',par$prgmTag,'_',chipName,'.sh'))
    
    if (opt$verbose>0)
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
    
    if (opt$verbose>0)
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
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]:{TAB}{TAB}Cluster Mode; Launching chip={chipName}, shell={lanShell}.{RET}"))
    system(lanShell)
    
    if (opt$single) break
  }
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]:{TAB}Cluster Mode; Done.{RET}"))
  
} else {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Load Manifest(s)::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]: Launching Samples in Linear Mode! isSingle={opt$single}"),"\n", sep='')
  
  pTracker <- NULL
  # pTracker <- timeTracker$new(verbose=opt$verbose)
  
  par$manDir <- NULL
  if (!is.null(opt$manDirPath) && length(opt$manDirPath)>0 && dir.exists(opt$manDirPath)) {
    par$manDir <- opt$manDirPath
  } else if (!is.null(opt$manDirName)) {
    par$manDir <- file.path(par$datDir, 'manifest',opt$manDirName)
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Both manDirPath and manDirName are NULL!!!{RET}{RET}"))
  }
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]: Set manifest search directory={par$manDir}.{RET}"))
  
  if (!is.null(opt$version) && !is.null(opt$percent)) {
    # opt$version <- paste(opt$version, opt$percent, sep='-')
    opt$version <- paste(opt$version, opt$percent, sep='.')
  }
  
  tar_man_tib <- NULL
  tar_man_tib <- get_manifest_list(
    file=opt$manifest, dir=par$manDir,
    platform=opt$platform, version=opt$version,
    suffix=opt$man_suffix,
    verbose=opt$verbose, tt=pTracker)
  
  tar_man_dat <- NULL
  tar_man_dat <- load_manifest_list(
    tib = tar_man_tib, field="path",
    verbose=opt$verbose, tt=pTracker)
  
  # Scratch Space on Manifests::
  if (FALSE) {
    tar_man_tibC <- get_manifest_list(
      file=opt$manifest, dir=file.path(par$datDir, 'manifest',opt$manDirName),
      suffix=opt$man_suffix,
      verbose=opt$verbose, tt=pTracker)
    
    tar_man_datC <- NULL
    tar_man_datC <- load_manifest_list(
      tib = tar_man_tibC, field="path",
      verbose=opt$verbose, tt=pTracker)
    
    # Compare Summaries::  1,081 don't have reference Next_Base/Color_Channel...
    tar_man_dat$`Chicago-S38` %>% 
      dplyr::filter(Probe_Type == "cg") %>%
      dplyr::group_by(Probe_Type,Probe_Design,DESIGN,col,COLOR_CHANNEL) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print(n=1000)
    
    tar_man_datC$`EPIC-B4` %>% 
      dplyr::group_by(Probe_Type,col,COLOR_CHANNEL) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% print(n=1000)

    # Build pseudo Chicago with EPIC controls
    core_man_names_vec <- tar_man_datC$`EPIC-B4` %>% 
      dplyr::filter(Probe_Type=="NEGATIVE") %>% names()
    
    epic_ctl_man_tib <- tar_man_datC$`EPIC-B4` %>% 
      dplyr::filter(Probe_Type != "cg") %>% 
      dplyr::filter(Probe_Type != "ch") %>% 
      dplyr::filter(Probe_Type != "rs") %>%
      dplyr::select(dplyr::all_of(core_man_names_vec)) %>%
      clean_tibble()
    
    chic_cpg_man_tib <- tar_man_dat$`Chicago-S38` %>% 
      dplyr::filter(Probe_Type == "cg") %>%
      dplyr::select(dplyr::all_of(core_man_names_vec)) %>%
      clean_tibble()
    
    chic_out_man_csv <- "/Users/bretbarnes/Documents/data/manifests/methylation/Chicago-Ober-Custom/Chicago-S38.manifest.sesame-base.cpg-sorted.csv.gz"
    chic_out_man_tib <- dplyr::bind_rows(chic_cpg_man_tib,epic_ctl_man_tib)
    
    chic_out_man_tib %>% 
      dplyr::group_by(Probe_Source,Probe_Type,Probe_Design,DESIGN,COLOR_CHANNEL) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>%
      print(n=1000)
    
    readr::write_csv(chic_out_man_tib, chic_out_man_csv)
    
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Load Auto Detection::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  auto_sam_tib <- NULL
  if (opt$auto_detect) {
    if (is.null(opt$auto_sam_csv))
      opt$auto_sam_csv <- file.path(par$datDir, 'ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo_Mean-COVIC-280-NP-ind_negs-0.02.csv.gz')
    
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]:{TAB} Loading auto_sam_csv={opt$auto_sam_csv}...{RET}"))
    auto_sam_tib <- suppressMessages(suppressWarnings(readr::read_csv(opt$auto_sam_csv) ))
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]:{TAB} Done.{RET}{RET}"))
  } else {
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]:{TAB} Will not Auto-Detect Sample.{RET}{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Load Pre-Defined Sesame Mask Probes::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # NOTE:: Pretty sure this is just for EPIC
  #  TBD:: Should load all pre-defined masked CpGs by manifest and pass that
  #   in as a list rather than vector and only mask if the manifest matches...
  mask_cpg_vec <- NULL
  mask_cpg_csv <- file.path(par$datDir, 'manifest/mask/sesame-general-mask.cpg.csv.gz')
  if (opt$mask_general) {
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]:{TAB} Loading mask_cpg_csv={mask_cpg_csv}...{RET}"))
    mask_cpg_vec <- suppressMessages(suppressWarnings( 
      readr::read_csv(mask_cpg_csv) )) %>% dplyr::pull(Probe_ID)
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]:{TAB} Done.{RET}{RET}"))
  } else {
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]:{TAB} Will not Mask Sample CpGs.{RET}{RET}"))
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
    
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]: parallelFunc={par$funcTag}: samples={sample_cnt}; ",
                     "num_cores={num_cores}, num_workers={num_workers}, Starting...{RET}"))
    
    # chipTimes <- foreach (prefix=prefixe_names, .inorder=T, .final = function(x) setNames(x, prefixe_names)) %dopar% {
    chipTimes <- foreach (prefix=prefixe_names, .combine = rbind) %dopar% {
      
      rdat <- NULL
      rdat <- sesamizeSingleSample(prefix=chipPrefixes[[prefix]],
                                   man=tar_man_dat, ref=auto_sam_tib, 
                                   opts=opt, defs=def,
                                   mask=mask_cpg_vec, platform=opt$forcedPlat,
                                   
                                   pvals=pval_vec,
                                   min_pvals=min_pval_vec,
                                   min_percs=min_perc_vec,
                                   workflows=workflow_vec,
                                   
                                   retData=par$retData,
                                   trackTime=opt$trackTime,
                                   verbose=opt$verbose, vt=3,tc=1)
      rdat
    }
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}] parallelFunc={par$funcTag}: Done.{RET}{RET}"))
  } else {
    par$funcTag <- 'sesamizeSingleSample-Linear'
    
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]: linearFunc={par$funcTag}: samples={sample_cnt}.{RET}"))
    
    rdat <- NULL
    for (prefix in prefixe_names) {
      if (opt$verbose>0)
        cat(glue::glue("[{par$prgmTag}]: linearFunc={par$funcTag}: Starting; prefix={prefix}...{RET}"))
      
      # par$retData <- TRUE
      # opt$verbose <- 40
      # ram_performance <- profmem({
      rdat <- NULL
      rdat <- sesamizeSingleSample(prefix=chipPrefixes[[prefix]],
                                   man=tar_man_dat, ref=auto_sam_tib, 
                                   opts=opt, defs=def,
                                   mask=mask_cpg_vec, platform=opt$forcedPlat,
                                   
                                   pvals=pval_vec,
                                   min_pvals=min_pval_vec,
                                   min_percs=min_perc_vec,
                                   workflows=workflow_vec,
                                   
                                   retData=par$retData,
                                   trackTime=opt$trackTime,
                                   verbose=opt$verbose, vt=3,tc=1)
      # })
      
      if (opt$verbose>0)
        cat(glue::glue("[{par$prgmTag}]: linearFunc={par$funcTag}: try_str={try_str}. Done.{RET}{RET}"))
      
      if (opt$single) break
    }
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}] parallelFunc={par$funcTag}: Done.{RET}{RET}"))
  }
  
  opt_csv <- file.path(opt$outDir, paste(par$prgmTag,'program-options.csv', sep='.') )
  par_csv <- file.path(opt$outDir, paste(par$prgmTag,'program-parameters.csv', sep='.') )
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value") %>% dplyr::arrange(Option)
  par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value") %>% dplyr::arrange(Params)
  
  readr::write_csv(opt_tib, opt_csv)
  readr::write_csv(par_tib, par_csv)
}

# Example of perfromance memory usage::
# p <- profmem({
#   x <- integer(1000)
#   Y <- matrix(rnorm(n = 10000), nrow = 100)
# })

if (FALSE) {
  # For sset conversion::
  #  sset_tib <- ssetToTib(sset=rdat$new_sset, source = "sigs", verbose = 10)
  #
  # Current ERROR::
  # 
  # Error in { : task 1 failed - "'x' must have 1 or more non-missing values"
  #   In addition: Warning message:
  #     Unknown or uninitialised column: `percent`. 
    
  ram_performance$bytes %>% as.vector() %>% sum(na.rm=TRUE)
  ram_performance$bytes %>% as.vector() %>% max(na.rm=TRUE)
  ram_performance$bytes %>% as.vector() %>% mean(na.rm=TRUE)
  ram_performance$bytes %>% as.vector() %>% median(na.rm=TRUE)
  which(ram_performance$bytes == ram_performance$bytes %>% as.vector() %>% max(na.rm=TRUE))
  ram_performance$trace[which(ram_performance$bytes == ram_performance$bytes %>% as.vector() %>% max(na.rm=TRUE))]
  
  #
  # Metric for Reference Sample Comparison::
  #
  rdat$cur_list$call_dat %>% 
    dplyr::summarise(raw_avg=mean(abs(raw_dB_ref), na.rm=TRUE), 
                     ind_avg=mean(abs(ind_dB_ref), na.rm=TRUE), 
                     raw_med=median(abs(raw_dB_ref), na.rm=TRUE), 
                     ind_med=median(abs(ind_dB_ref), na.rm=TRUE), 
                     raw_sd=sd(abs(raw_dB_ref), na.rm=TRUE), 
                     ind_sd=sd(abs(ind_dB_ref), na.rm=TRUE)
    )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
