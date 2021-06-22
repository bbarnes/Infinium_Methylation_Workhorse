
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              COVIC Merging and Spliting of Swifthoof Data::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

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
par$prgmDir <- 'analysis'
par$prgmTag <- 'merge_builds'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# Predefined human sample sheet name::
par$humanSampleSheetName <- 'humanSampleSheet.csv'
opt$joinType <- "full"
par$noob_sub <- FALSE

# File Based Parameters::
opt$inputsCsv <- NULL

# Directory Parameters::
opt$outDir  <- NULL
opt$datDir  <- NULL

# Run Parameters::
opt$runName   <- NULL
opt$sampleCsv <- NULL
opt$manifest  <- NULL
opt$findSampleSheet <- FALSE

# Class Parameters::
# Really simple test to make sure we can seperate the sexes...
opt$classVar <- 'Karyotype_0_Call'
opt$classVar <- 'Sample_Name'
opt$classVar <- 'Sample_Class'
opt$workflow <- NULL

opt$select <- FALSE
opt$clean_gta <- FALSE
opt$forceUnq  <- TRUE

# Sample Sheet Parameters::
opt$addSampleName    <- FALSE
opt$addPathsCall     <- TRUE
opt$addPathsSset     <- FALSE

opt$flagDetectPval   <- FALSE
opt$flagSampleDetect <- FALSE
opt$flagRefMatch     <- FALSE

opt$pvalDetectMinKey <- NULL
opt$pvalDetectMinVal <- NULL

# Chip Platform and Version Parameters::
opt$platform <- NULL
opt$version  <- NULL

# Output Format Parameters::
opt$percisionBeta <- 4
opt$percisionPval <- 6

# Parallel/Cluster Options::
opt$execute  <- TRUE
opt$single   <- FALSE
opt$parallel <- FALSE
opt$cluster  <- FALSE

# Make clean output
opt$clean <- FALSE

opt$clean_source <- TRUE

# Executable Paramaters::
opt$Rscript <- NULL

# Verbosity Options::
opt$verbose <- 3

cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

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
  
  opt$addPval     <- TRUE
  opt$buildDml    <- FALSE
  opt$buildDbl    <- FALSE
  opt$buildModels <- FALSE
  
  opt$single   <- FALSE
  opt$single   <- TRUE
  opt$cluster  <- FALSE
  opt$parallel <- FALSE
  
  opt$samplePvalName <- "Poob_Pass_0_Perc"
  opt$samplePvalPerc <- 90
  
  # Loci Level Filtering Parameters::
  opt$lociBetaKey <- "i_beta,ind_beta"
  opt$lociBetaKey <- "ind_beta"
  
  opt$lociPvalKey <- "i_poob,i_negs"
  opt$lociPvalKey <- "i_poob"
  
  opt$lociPvalMin <- "0.02,0.1,0.5,0.9,1.0"
  opt$lociPvalMin <- "0.05,0.1"
  opt$lociPvalMin <- "0.1"
  opt$lociPvalMin <- "0.05"
  
  opt$seeds <- "13,17,42,43,57,61,69"
  opt$seeds <- "13,42"
  opt$seeds <- NULL
  
  opt$featuresCsv <- NULL
  opt$featuresDml <- NULL
  
  opt$classVar <- 'Sample_Class'
  opt$platform <- 'EPIC'
  opt$version  <- 'B4'
  
  opt$clean <- FALSE
  opt$clean <- TRUE
  
  # opt$trainClass <- paste('HELA','JURKAT','MCF7','RAJI', sep=',')
  # opt$trainClass <- paste('Xa','XaXaY','XaXi','XaXiY','XaY', sep=',')
  # opt$trainClass <- paste('XaXi','XaY', sep=',')
  # opt$trainClass <- paste('nSARSCov2', 'pSARSCov2', sep=',')
  
  #
  # Pre-defined local options runTypes::
  #
  par$local_runType <- 'CORE'
  par$local_runType <- 'EXCBR'
  par$local_runType <- 'NZT'
  par$local_runType <- 'COVIC'
  par$local_runType <- 'GENK'
  par$local_runType <- 'COVID'
  par$local_runType <- 'GRCm38'
  par$local_runType <- 'qcMVP'
  par$local_runType <- 'EPIC-8x1-EM-Sample-Prep'
  par$local_runType <- 'Chicago-Ober-Custom'
  par$local_runType <- 'NA12878'
  par$local_runType <- 'qcMVP2'
  
  if (FALSE) {
    
  } else if (par$local_runType=='Chicago-Ober-Custom') {
    opt$runName  <- par$local_runType
    
    opt$workflow    <- "ind"
    
    opt$single   <- FALSE
    # opt$parallel <- TRUE
    
    opt$platform <- NULL
    opt$version  <- NULL
    opt$forceUnq <- FALSE
    
    par$vstr <- ""
    par$vstr <- "v1"
    par$vstr <- "v2"
    
    opt$classVar <- "AutoSample_R2_Key_2"
    
    opt$datDir <- paste(
      # file.path(par$topDir, 'scratch/swifthoof',opt$runName,"Chicago/S38/swifthoof_main"),
      file.path(par$topDir, 'scratch/swifthoof',opt$runName,"Chicago/S38",par$vstr,"swifthoof_main"),
      sep=',')
    
    opt$runName <- paste(opt$runName,par$vstr, sep="-")
    
  } else if (par$local_runType=='EPIC-8x1-EM-Sample-Prep') {
    opt$runName  <- par$local_runType
    
    opt$workflow    <- "ind"
    
    opt$single   <- FALSE
    # opt$parallel <- TRUE
    
    opt$datDir <- paste(
      # file.path(par$topDir, 'scratch/RStudio/swifthoof_main',opt$runName),
      # file.path(par$topDir, 'scratch/RStudio/swifthoof_main',"EPIC-8x1-EM-Sample-Prep.v0"),
      file.path(par$topDir, 'scratch/swifthoof_main',opt$runName),
      file.path(par$topDir, 'scratch/swifthoof_main',"EPIC-8x1-EM-Sample-Prep.v0"),
      sep=',')
    
  } else if (par$local_runType=='NA12878') {
    opt$runName  <- par$local_runType
    
    opt$manifest <- file.path(par$datDir, "manifest/core/EPIC-B4.manifest.sesame-base.cpg-sorted.csv.gz")
    opt$workflow <- "ind"
    
    opt$verbose <- 30
    opt$platform <- NULL
    opt$version  <- NULL
    opt$forceUnq <- FALSE
    
    par$platform <- "Rand3"
    par$platform <- "Rand2"
    par$platform <- "Rand1"
    par$platform <- NULL
    
    # opt$classVar <- "detect_manifest"
    opt$classVar <- "detect_platform"
    opt$classVar <- "detect_version"
    
    opt$single   <- FALSE
    # opt$parallel <- TRUE
    opt$addPathsCall <- TRUE
    opt$addPathsSset <- TRUE
    
    par$mixData <- TRUE
    par$allData <- TRUE
    
    if (!par$mixData) {
      opt$datDir <- paste(
        # file.path(par$topDir, 'scratch/RStudio/swifthoof_main',opt$runName),
        
        # Granularity 0:
        # file.path(par$topDir, 'scratch/noob-sub', par$platform),
        
        # Granularity 1:
        file.path(par$topDir, 'scratch/noob-sub/granular', par$platform),
        sep=',')
      
    } else if (!par$allData) {
      # For individual data::
      opt$datDir <- paste(
        # Granularity 0:
        # file.path(par$topDir, 'scratch/noob-sub', par$platform),
        
        # Granularity 1:
        file.path(par$topDir, 'scratch/noob-sub/granular', par$platform),
        sep=',')
      if (!is.null(par$platform)) opt$runName <- paste(opt$runName,"granular1",par$platform, sep='-')
    } else {
      # For all data::
      opt$datDir <- paste(
        # Granularity 0: Will's sub-sample
        # file.path(par$topDir, 'scratch/noob-sub/granular0/Rand1'),
        # file.path(par$topDir, 'scratch/noob-sub/granular0/Rand2'),
        # file.path(par$topDir, 'scratch/noob-sub/granular0/Rand3'),
        
        # Granularity 1: Will's sub-sample
        # file.path(par$topDir, 'scratch/noob-sub/granular/Rand1'),
        # file.path(par$topDir, 'scratch/noob-sub/granular/Rand2'),
        # file.path(par$topDir, 'scratch/noob-sub/granular/Rand3'),
        # file.path(par$topDir, 'scratch/noob-sub/granular/Rand4'),
        # file.path(par$topDir, 'scratch/noob-sub/granular/Rand5'),
        
        # Rand 1-3 all data::
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand1"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand2"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand3"),
        
        # Rand 1-3 all data:: 0.0, 10.5, 5.10, 20.20, 60.60, 100.100
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand1/0.0"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand1/10.5"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand1/5.10"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand1/20.20"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand1/60.60"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand1/100.100"),
        # 
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand2/0.0"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand2/10.5"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand2/5.10"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand2/20.20"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand2/60.60"),
        # file.path(par$topDir, "scratch/noob-sub/granular2/Rand2/100.100"),
        
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand4/0.0"),
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand4/10.5"),
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand4/5.10"),
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand4/20.20"),
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand4/60.60"),
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand4/100.100"),
        
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand5/0.0"),
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand5/10.5"),
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand5/5.10"),
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand5/20.20"),
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand5/60.60"),
        file.path(par$topDir, "scratch/noob-sub/granular2/Rand5/100.100"),
        
        sep=',')
      
      # opt$runName <- paste(opt$runName,"granular0-1", sep='-')
      # opt$runName <- paste(opt$runName,"granular1.5", sep='-')
      # opt$runName <- paste(opt$runName,"granular0.3", sep='-')
      # opt$runName <- paste(opt$runName,"granular0-1.5", sep='-')
      # opt$runName <- paste(opt$runName,"granular-3.1", sep='-')
      
      # Rand 1-3 all data::
      # opt$runName <- paste(opt$runName,"granular-3-rand1-2-3", sep='-')
      
      # Rand 1-3 all data:: only 
      opt$runName <- paste(opt$runName,"granular-rand1-2.6-test", sep='-')
      
      opt$joinType <- "inner"
    }
    par$noob_sub <- TRUE
    par$noob_sub <- FALSE
    
  } else if (par$local_runType=='COVID') {
    
    opt$flagDetectPval   <- FALSE
    opt$flagSampleDetect <- FALSE
    opt$flagRefMatch     <- FALSE
    
    opt$pvalDetectMinKey <- NULL
    opt$pvalDetectMinVal <- NULL
    
    opt$clean_gta <- TRUE
    
    opt$sampleCsv <- 
      file.path(par$topDir,'data/CustomContent/COVID-19_HLA/data/directDetection/COVID_Direct/mappings/Chip_Result_Mapping.Sample_Class.csv.gz')
    
    opt$trainClass <- paste('nSARSCov2', 'pSARSCov2', sep=',')
    # opt$trainClass <- paste('NEGATIVE', 'POSTIVE', sep=',')
    
    opt$runName  <- 'COVID-Direct-Set1'
    opt$classVar <- 'Sample_Class'
    opt$platform <- 'COVID'
    opt$version  <- 'C1'
    
    opt$datDir <- paste(
      file.path(par$topDir, 'scratch/COVID_All/swifthoof_main',opt$runName),
      sep=',')
    
  } else if (par$local_runType=='COVIC') {
    opt$runName  <- 'COVIC-Set7-06082020'
    opt$runName  <- 'COVIC-Set5-10062020'
    opt$runName  <- 'COVIC-Set1-15052020'
    
    opt$version  <- 'C0'
    
    opt$datDir <- paste(
      file.path(par$topDir, 'scratch/swifthoof_main', opt$runNameA),
      file.path(par$topDir, 'scratch/swifthoof_main', opt$runNameB),
      file.path(par$topDir, 'scratch/swifthoof_main', opt$runNameC),
      sep=',')
    
    opt$trainClass <- paste('nSARSCov2', 'pSARSCov2', sep=',')
    
    # Loci (Feature) Selection Parameters::
    #
    # opt$featuresCsv <- paste( file.path(par$datDir, 'sampleSheets/dmls/Ivana-145.csv.gz'),
    #                           file.path(par$datDir, 'sampleSheets/dmls/Genknowme-2043.csv.gz'),
    #                           file.path(par$datDir, 'sampleSheets/dmls/COVIC-hit.csv.gz'),
    #                           sep=',')
    # opt$featuresDml <- "100"
    
  } else if (par$local_runType=='GRCm38') {
    par$runNameA  <- 'ILMN_mm10_betaTest_17082020'
    par$runNameB  <- 'VanAndel_mm10_betaTest_31082020'
    par$runNameC  <- 'MURMETVEP_mm10_betaTest_06082020'
    
    opt$runName   <- 'mm10_controls'
    
    opt$datDir  <- paste(
      file.path(par$topDir, 'scratch/GRCh38/swifthoof_main', par$runNameA),
      file.path(par$topDir, 'scratch/GRCh38/swifthoof_main', par$runNameB),
      file.path(par$topDir, 'scratch/GRCh38/swifthoof_main', par$runNameC),
      sep=',')
    
    opt$classVar <- 'Sample_Name'
    opt$platform <- 'LEGX'
    opt$version  <- ''
    
    par$titration <- FALSE
    par$titration <- TRUE
    if (par$titration) {
      opt$runName  <- 'mm10-ILS-VAI.Titration'
      opt$trainClass <- paste('T00DZ','T50DZ','T99DZ', sep=',')
    } else {
      opt$runName  <- 'mm10-ILS-VAI.Replicate'
      opt$trainClass <- paste('RepAC','RepS3','RepM1','RepSA', sep=',')
    }
    
    opt$datDir  <- paste(
      file.path(par$topDir,'scratch/merge_builds/LEGX/S1/Sample_Name',opt$runName),
      sep=',')
    
  } else if (par$local_runType=='qcMVP2') {
    # "/Users/bretbarnes/Documents/data/VA_MVP/docker-v.1.11/"
    
    opt$sampleCsv <- file.path(par$topDir, "data/sampleSheets/VA-Controls/Batch3/VA-Controls-Batch3.sampleSheet.csv")
    
    opt$manifest <- file.path(par$datDir, "manifest/core/EPIC-B4.manifest.sesame-base.cpg-sorted.csv.gz")
    opt$workflow <- "ind"
    
    opt$classVar <- 'AutoSample_dB_Key'
    opt$classVar <- 'AutoSample_dB_Key_1'
    
    opt$platform <- 'EPIC'
    opt$version  <- 'B4'
    
    # opt$datDir <- 
    #   do.call(paste, c(as.list(list.files(par$datSrc, full.names = TRUE)), sep=","))
    
    opt$datDir <- paste(
      file.path(par$topDir, "scratch/swifthoof/VA-MVP-Akesogen_Batch3/swifthoof_main"),
      file.path(par$topDir, "scratch/swifthoof/VA-MVP-IBX_Batch3/swifthoof_main"),
      sep=","
    )
    
    opt$runName  <- paste(par$local_runType, sep='_')
    
    opt$addPathsCall <- TRUE
    opt$addPathsSset <- FALSE
    
  } else if (par$local_runType=='qcMVP') {
    # "/Users/bretbarnes/Documents/data/VA_MVP/docker-v.1.11/"
    
    opt$classVar <- 'AutoSample_dB_Key'
    opt$classVar <- 'AutoSample_dB_Key_1'
    
    opt$platform <- 'EPIC'
    opt$version  <- 'B4'
    
    # par$runName1  <- 'CNTL-Samples_VendA_10092020'
    # par$runName2  <- 'CNTL-Samples_VendB_10092020'
    # 
    # par$runName3  <- 'BETA-8x1-EPIC-Core'
    # par$runName4  <- 'DELTA-8x1-EPIC-Core'
    # par$runName5  <- 'DELTA-24x1-EPIC'
    # par$runName6  <- 'BETA-8x1-EPIC-Bad'
    # 
    # par$runName7  <- 'COVIC-Set1-15052020'
    # par$runName8  <- 'COVIC-Set2-31052020'
    # par$runName9  <- 'COVIC-Set3-05062020'
    # par$runName10 <- 'COVIC-Set4-09062020'
    # par$runName11 <- 'COVIC-Set5-10062020'
    # par$runName12 <- 'COVIC-Set7-06082020'
    # par$runName13 <- 'COVIC-Set8-26182020'
    # 
    # par$runName14 <- 'Excalibur-Old-1609202'
    # par$runName15 <- 'Excalibur-New-1609202'
    
    #
    # pass_vec <- c(par$runName1,par$runName2,
    #               par$runName3,par$runName4,par$runName5,
    #               par$runName7,par$runName11,par$runName12,
    #               par$runName15)
    # 
    # opt$datDir  <- paste(
    #   # Known passed chips::
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName1),
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName2),
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName3),
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName4),
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName5),
    #   
    #   # Known failed chips::
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName6),
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName7),
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName8),
    #   
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName9),
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName10),
    # 
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName11),
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName12),
    #   # file.path(par$topDir,'scratch/swifthoof_main',par$runName13),
    #   
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName14),
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runName15),
    #   
    #   sep=',')
    
    par$datSrc <- file.path(par$topDir, "data/VA_MVP/docker-v.1.11")
    opt$datDir <- 
      do.call(paste, c(as.list(list.files(par$datSrc, full.names = TRUE)), sep=","))
    
    opt$runName  <- paste(par$local_runType, sep='_')
    
    opt$addPathsCall <- TRUE
    opt$addPathsSset <- TRUE
    
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$outDir <- file.path(par$topDir, 'scratch',par$runMode)
  
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
    
    # Directory Parameters::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    make_option(c("-d","--datDir"), type="character", default=opt$datDir, 
                help="List of Build Directory(s), commas seperated [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--sampleCsv"), type="character", default=opt$sampleCsv, 
                help="Human provide sample sheet labeling [default= %default]", metavar="character"),
    make_option(c("--manifest"), type="character", default=opt$manifest, 
                help="Human provide manifest [default= %default]", metavar="character"),
    make_option(c("--select"), action="store_true", default=opt$select, 
                help="Boolean variable to only select samples from provided sample sheet [default= %default]", metavar="boolean"),
    
    make_option(c("--findSampleSheet"), action="store_true", default=opt$findSampleSheet,
                help="Boolean variable to search or human provided sample sheet in build directories [default= %default]", metavar="boolean"),
    
    # Chip Platform and Version Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform name (HM50, EPIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest version (B2, B4, C0) [default= %default]", metavar="character"),
    
    # Class Parameters::
    make_option(c("--classVar"), type="character", default=opt$classVar, 
                help="Classification Variable Name [default= %default]", metavar="character"),
    make_option(c("--workflow"), type="character", default=opt$workflow, 
                help="Target Workflow Variable Name [default= %default]", metavar="character"),
    
    # Sample Sheet Parameters::
    make_option(c("--addSampleName"), action="store_true", default=opt$addSampleName, 
                help="Sample Sheet processing to add Auto-SampleNames (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--addPathsCall"), action="store_true", default=opt$addPathsCall, 
                help="Sample Sheet processing to add Calls Local Full Paths (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--addPathsSset"), action="store_true", default=opt$addPathsSset, 
                help="Sample Sheet processing to add Signals Local Full Paths (mostly testing stuff) [default= %default]", metavar="boolean"),
    
    make_option(c("--flagDetectPval"), action="store_true", default=opt$flagDetectPval, 
                help="Sample Sheet processing to add flag for failed detected samples by failed loci percent (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--flagSampleDetect"), action="store_true", default=opt$flagSampleDetect, 
                help="Sample Sheet processing to add flag for failed Auto-Detected Samples (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--flagRefMatch"), action="store_true", default=opt$flagRefMatch, 
                help="Sample Sheet processing to add flag for failed Auto-Detected Methods Agreements (mostly testing stuff) [default= %default]", metavar="boolean"),
    
    make_option(c("--pvalDetectMinKey"), type="character", default=opt$pvalDetectMinKey, 
                help="Sample Sheet processing to Min Detection Pval Key [default= %default]", metavar="character"),
    make_option(c("--pvalDetectMinVal"), type="double", default=opt$pvalDetectMinVal,
                help="Sample Sheet processing to Min Detection Pval Value [default= %default]", metavar="double"),
    
    # Output Format Parameters::
    make_option(c("--percisionBeta"), type="integer", default=opt$percisionBeta,
                help="Rounding percision for beta values in calls output files [default= %default]", metavar="integer"),
    make_option(c("--percisionPval"), type="integer", default=opt$percisionPval,
                help="Rounding percision for detection p-values in calls output files [default= %default]", metavar="integer"),
    
    make_option(c("--clean_gta"), action="store_true", default=opt$clean_gta, 
                help="Boolean variable to remove Genotyping tails from Sentrix Names (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--forceUnq"), action="store_false", default=opt$forceUnq, 
                help="Boolean variable to force unique Sentrix Names [default= %default]", metavar="boolean"),
    make_option(c("--joinType"), type="character", default=opt$joinType, 
                help="Data merging join type (full, inner) [default= %default]", metavar="character"),
    
    # Parallel/Cluster Parameters::
    make_option(c("--execute"), action="store_true", default=opt$execute, 
                help="Boolean variable to shell scripts (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
    # Make clean output
    make_option(c("--clean"), action="store_true", default=opt$clean, 
                help="Boolean variable to clean output directory (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--clean_source"), action="store_true", default=opt$clean_source, 
                help="Boolean variable to clean output directory names for plotting [default= %default]", metavar="boolean"),
    
    # Executable Paramaters::
    make_option(c("--Rscript"), type="character", default=opt$Rscript, 
                help="Rscript path [default= %default]", metavar="character"),
    
    # Verbosity Options::
    make_option(c("-v", "--verbose"), type="integer", default=opt$verbose, 
                help="0-5 (5 is very verbosity) [default= %default]", metavar="integer")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Program Initialization::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_reqs <- c('runMode','prgmTag','scrDir','datDir','exePath')
opt_reqs <- c('outDir','Rscript','verbose')
opt_reqs <- c('outDir','Rscript','verbose','clean',
              'datDir','runName','classVar','workflow',
              'addSampleName','addPathsCall','addPathsSset',
              'flagDetectPval','flagSampleDetect','flagRefMatch',
              # 'platform','version',
              'percisionBeta','percisionPval',
              'execute','single','parallel','cluster')

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

opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

blds_dir_vec  <- opt$datDir %>% 
  str_split(pattern=',', simplify=TRUE) %>% as.vector()

if (is.null(opt$classVar)) opt$classVar <- 'Source_Sample_Name'
class_var <- rlang::sym(opt$classVar)
class_idx <- rlang::sym("Class_Idx")

if (opt$verbose>0) {
  cat(glue::glue("[{par$prgmTag}]: blds_dir_vec={RET}") )
  print(blds_dir_vec)
}

if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)

if (!is.null(opt$platform)) opt$outDir <- file.path(opt$outDir, opt$platform)
if (!is.null(opt$version))  opt$outDir <- file.path(opt$outDir, opt$version)
if (!is.null(opt$classVar)) opt$outDir <- file.path(opt$outDir, opt$classVar)
if (!is.null(opt$workflow)) opt$outDir <- file.path(opt$outDir, opt$workflow)

if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Bulding outDir={opt$outDir}!{RET}") )
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Built; outDir={opt$outDir}!{RET}{RET}") )


if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Loading source manifest={opt$manifest}...{RET}") )
src_man_tib <- suppressMessages(suppressWarnings( readr::read_csv(opt$manifest)))

if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Load Auto Detect Sample Sheets::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

hum_ss_tib  <- NULL
auto_ss_tib <- NULL

for (curDir in blds_dir_vec) {
  
  if (!dir.exists(curDir)) {
    cat(glue::glue("[{par$prgmTag}]:{TAB} Failed to find dir={curDir}, skipping...{RET}") )
    next
  }
  
  if (opt$findSampleSheet) {
    
    cur_hm_csv <- file.path(curDir, par$humanSampleSheetName)
    cat(glue::glue("[{par$prgmTag}]:{TAB}Loading Human Classification; cur_hm_csv={cur_hm_csv}!{RET}") )
    
    cur_hm_tib <- suppressMessages(suppressWarnings( readr::read_csv(cur_hm_csv) ))
    cur_hm_len <- cur_hm_tib %>% base::nrow()
    cat(glue::glue("[{par$prgmTag}]:{TAB}Done. Loading Human Classification; cur_hm_len={cur_hm_len}!{RET}{RET}") )
    # print(cur_hm_tib)
    
    hum_ss_tib <- dplyr::bind_rows(hum_ss_tib,cur_hm_tib)
  }
  
  suffix='AutoSampleSheet.csv.gz'
  
  cur_ss_tib <- NULL
  cur_ss_tib <- 
    loadAutoSampleSheets(dir=curDir, 
                         platform=opt$platform, manifest=opt$version, 
                         workflow=opt$workflow, suffix=suffix,
                         
                         addSampleName=opt$addSampleName, 
                         addPathsCall=opt$addPathsCall, 
                         addPathsSset=opt$addPathsSset,
                         
                         flagDetectPval=opt$flagDetectPval, 
                         flagSampleDetect=opt$flagSampleDetect, 
                         flagRefMatch=opt$flagRefMatch,
                         
                         pvalDetectMinKey=opt$pvalDetectMinKey, 
                         pvalDetectMinVal=opt$pvalDetectMinVal,
                         clean_gta=opt$clean_gta,
                         verbose=opt$verbose,vt=3,tc=1,tt=pTracker)
  
  if (is.null(cur_ss_tib)) {
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]:{TAB} Failed to find any suffix={suffix} ",
                     "in dir={curDir}. Skipping...{RET}") )
    next
  }
  
  cur_ss_tib <- cur_ss_tib %>%
    dplyr::mutate(build_source=base::basename(curDir)) %>% clean_tibble()
  
  # Join builds::
  #
  auto_ss_tib <- auto_ss_tib %>% dplyr::bind_rows(cur_ss_tib) %>% clean_tibble()
  
  cur_ss_len  <- cur_ss_tib %>% base::nrow()
  auto_ss_len <- auto_ss_tib %>% base::nrow()
  
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]:{TAB} Raw Auto Sample Sheet ",
                   "Current={cur_ss_len}, Total={auto_ss_len}.{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Ensure Sentrix_Name is Unique::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$forceUnq) {
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]:{TAB} Only keeping unique Setrix_Name rows!{RET}"))
  
  auto_ss_tib <- auto_ss_tib %>% dplyr::distinct(Sentrix_Name, .keep_all=TRUE)
} else {
  tot_sent_cnt <- auto_ss_tib %>% dplyr::pull(Sentrix_Name) %>% length()
  unq_sent_cnt <- auto_ss_tib %>% dplyr::pull(Sentrix_Name) %>% unique() %>% length()
  
  if (tot_sent_cnt != unq_sent_cnt) {
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]:{TAB} Making Setrix_Name unique!{RET}"))
    
    auto_ss_tib <- auto_ss_tib %>%
      # dplyr::arrange(!!class_var) %>%
      dplyr::mutate(Sentrix_Name=paste(Sentrix_Name,detect_manifest, sep='-'))
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Expand Class Variable if its a Double::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

class_type <- auto_ss_tib %>% 
  dplyr::select(!!class_var) %>% 
  purrr::map_chr(pillar::type_sum) %>% 
  paste(collapse = "_")

if (class_type=="dbl") {
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]:{TAB} Splitting detect_manifest...{RET}"))
  
  auto_ss_tib <- auto_ss_tib %>% 
    tidyr::separate(detect_manifest, 
                    into=c("Sample_Rand","Sample_Per1","Sample_Per2"), 
                    remove=FALSE, convert=TRUE, sep="[.-]") %>% 
    dplyr::mutate(Sample_Perc=paste0("I",Sample_Per1,"_II",Sample_Per2),
                  Sample_Per1_Str=paste0("I",Sample_Per1), Sample_Per2_Str=paste0("II",Sample_Per2)) %>%
    dplyr::select(Sentrix_Name,Sample_Rand,Sample_Perc,Sample_Per1,Sample_Per2,
                  Sample_Per1_Str,Sample_Per2_Str,
                  dplyr::everything())
  
  auto_ss_sum <- auto_ss_tib %>% 
    dplyr::group_by(Sample_Per1,Sample_Per2) %>% 
    # dplyr::group_by(Sample_Rand,Sample_Per1,Sample_Per2) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  auto_ss_sum %>% print(n=base::nrow(auto_ss_sum))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Sample Sheet Summary Check::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Geneknowme Update::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opt$vers_sum1_csv <- file.path(opt$outDir, paste(par$outName, 'detect_Inf-Percent-summary1.csv.gz', sep='_') )
  vers_sum1_tib <- auto_ss_tib %>% 
    # dplyr::mutate(detect_version_int=Sample_Per1_Str) %>%
    dplyr::group_by(Sample_Per1,Sample_Per2) %>%
    dplyr::arrange(Sample_Per1,Sample_Per2) %>%
    dplyr::summarise(Min_PV=min(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE), 
                     Max_PV=max(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     Avg_PV=mean(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE), 
                     Med_PV=median(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     Std_PV=sd(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     TotalN=n(),
                     Min_R2=min(AutoSample_R2_Val_1, na.rm=TRUE), 
                     Max_R2=max(AutoSample_R2_Val_1, na.rm=TRUE),
                     Avg_R2=mean(AutoSample_R2_Val_1, na.rm=TRUE), 
                     Med_R2=median(AutoSample_R2_Val_1, na.rm=TRUE),
                     Std_R2=sd(AutoSample_R2_Val_1, na.rm=TRUE),
                     
                     Min_dB=min(AutoSample_dB_Val_1, na.rm=TRUE), 
                     Max_dB=max(AutoSample_dB_Val_1, na.rm=TRUE),
                     Avg_dB=mean(AutoSample_dB_Val_1, na.rm=TRUE), 
                     Med_dB=median(AutoSample_dB_Val_1, na.rm=TRUE),
                     Std_dB=sd(AutoSample_dB_Val_1, na.rm=TRUE),
                     
                     Min_Age=min(AgeSkinBlood_1, na.rm=TRUE), 
                     Max_Age=max(AgeSkinBlood_1, na.rm=TRUE),
                     Avg_Age=mean(AgeSkinBlood_1, na.rm=TRUE), 
                     Med_Age=median(AgeSkinBlood_1, na.rm=TRUE),
                     Std_Age=sd(AgeSkinBlood_1, na.rm=TRUE),
                     .groups="drop"
    ) %>% dplyr::mutate_if(is.double, base::round, 4)
  vers_sum1_tib %>% print(n=base::nrow(vers_sum1_tib))
  readr::write_csv(vers_sum1_tib, opt$vers_sum1_csv)
  
  opt$vers_sum2_csv <- file.path(opt$outDir, paste(par$outName, 'detect_Inf-Percent-summary2.csv.gz', sep='_') )
  vers_sum2_tib <- auto_ss_tib %>% 
    # dplyr::mutate(detect_version_int=Sample_Per1_Str) %>%
    dplyr::group_by(Sample_Per1,Sample_Per2) %>%
    dplyr::arrange(Sample_Per1,Sample_Per2) %>%
    dplyr::summarise(Min_PV=min(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE), 
                     Max_PV=max(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     Avg_PV=mean(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE), 
                     Med_PV=median(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     Std_PV=sd(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     TotalN=n(),
                     Min_R2=min(AutoSample_R2_Val_2, na.rm=TRUE), 
                     Max_R2=max(AutoSample_R2_Val_2, na.rm=TRUE),
                     Avg_R2=mean(AutoSample_R2_Val_2, na.rm=TRUE), 
                     Med_R2=median(AutoSample_R2_Val_2, na.rm=TRUE),
                     Std_R2=sd(AutoSample_R2_Val_2, na.rm=TRUE),
                     
                     Min_dB=min(AutoSample_dB_Val_2, na.rm=TRUE), 
                     Max_dB=max(AutoSample_dB_Val_2, na.rm=TRUE),
                     Avg_dB=mean(AutoSample_dB_Val_2, na.rm=TRUE), 
                     Med_dB=median(AutoSample_dB_Val_2, na.rm=TRUE),
                     Std_dB=sd(AutoSample_dB_Val_2, na.rm=TRUE),
                     
                     Min_Age=min(AgeSkinBlood_2, na.rm=TRUE), 
                     Max_Age=max(AgeSkinBlood_2, na.rm=TRUE),
                     Avg_Age=mean(AgeSkinBlood_2, na.rm=TRUE), 
                     Med_Age=median(AgeSkinBlood_2, na.rm=TRUE),
                     Std_Age=sd(AgeSkinBlood_2, na.rm=TRUE),
                     .groups="drop"
    ) %>% dplyr::mutate_if(is.double, base::round, 4)
  vers_sum2_tib %>% print(n=base::nrow(vers_sum2_tib))
  readr::write_csv(vers_sum2_tib, opt$vers_sum2_csv)
  
  
  age_dif_tib <- labs_ss_tib %>% dplyr::mutate(AgeSkinBlod_Dif=AgeSkinBlood_1-AgeSkinBlood_2) 
  gg <- ggplot2::ggplot(data=age_dif_tib, aes(Sample_Per1, group=AgeSkinBlod_Dif)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per2)) + ggtitle("Infinium II Bins")
  
  
  
  auto_ss_len <- auto_ss_tib %>% base::nrow()
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]: Done. Raw Auto Sample Sheet; ",
                   "Total={auto_ss_len}.{RET}{RET}"))
  
  # Metrics::
  #  AgeSkinBlood_2
  #
  plot_tib <- auto_ss_tib %>% dplyr::arrange(Sample_Per1,Sample_Per2)
  
  # Violin::
  # Sample_Per1_Str
  # Sample_Per2_Str
  
  # ggplot2::ggplot(data=plot_tib, aes(cg_calls_pass_perc_1, AgeSkinBlood_2)) +
  #   ggplot2::geom_violin() + 
  #   ggplot2::facet_grid(cols = vars(Sample_Per2_Str), rows=vars(Sample_Per1_Str))
  
  ggplot2::ggplot(data=plot_tib, aes(Sample_Per1, AgeSkinBlood_2)) +
    ggplot2::geom_violin() + 
    ggplot2::facet_grid(cols = vars(Sample_Per2))
  
  ggplot2::ggplot(data=plot_tib, aes(Sample_Per2, AgeSkinBlood_2)) +
    ggplot2::geom_violin() + 
    ggplot2::facet_grid(cols = vars(Sample_Per1))
  
  # Heatmap::
  ggplot(data = plot_tib, aes(x=Sample_Per1_Str, y=Sample_Per2_Str, fill=AgeSkinBlood_2)) + 
    geom_tile()
  ggplot(data = plot_tib, aes(x=Sample_Per1_Str, y=Sample_Per2_Str, fill=AutoSample_dB_Val_2)) + 
    geom_tile()
  ggplot(data = plot_tib, aes(x=Sample_Per1_Str, y=Sample_Per2_Str, fill=cg_calls_pass_perc_1)) + 
    geom_tile()
  
  
  ggplot(data = plot_tib, aes(x=Sample_Per1, y=Sample_Per2, fill=AgeSkinBlood_2)) + 
    geom_tile()
  ggplot(data = plot_tib, aes(x=Sample_Per1, y=Sample_Per2, fill=AutoSample_dB_Val_2)) + 
    geom_tile()
  ggplot(data = plot_tib, aes(x=Sample_Per1, y=Sample_Per2, fill=cg_calls_pass_perc_1)) + 
    geom_tile()
  
  # Boxplot::
  box_pdf1 <- file.path(opt$outDir, paste(opt$runName, 'box1.pdf', sep='_') )
  gg <- ggplot2::ggplot(data=plot_tib, aes(Sample_Per1, AgeSkinBlood_2)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per2)) + ggtitle("Infinium II Bins")
  ggplot2::ggsave(box_pdf1,gg)
  
  box_pdf2 <- file.path(opt$outDir, paste(opt$runName, 'box2.pdf', sep='_') )
  gg <- ggplot2::ggplot(data=plot_tib, aes(Sample_Per2, AgeSkinBlood_2)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per1))  + ggtitle("Infinium I Bins")
  ggplot2::ggsave(box_pdf2,gg)
  
  # Metrics::
  #  AutoSample_dB_Val_2
  box_pdf1 <- file.path(opt$outDir, paste(opt$runName, 'AutoSample_dB_Val_2.box1.pdf', sep='_') )
  gg <- ggplot2::ggplot(data=plot_tib, aes(Sample_Per1, AutoSample_dB_Val_2)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per2)) + ggtitle("Infinium II Bins")
  ggplot2::ggsave(box_pdf1,gg)
  
  box_pdf2 <- file.path(opt$outDir, paste(opt$runName, 'AutoSample_dB_Val_2.box2.pdf', sep='_') )
  gg <- ggplot2::ggplot(data=plot_tib, aes(Sample_Per2, AutoSample_dB_Val_2)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per1))  + ggtitle("Infinium I Bins")
  ggplot2::ggsave(box_pdf2,gg)
  
  # Metrics::
  #  cg_calls_pass_perc_1
  box_pdf1 <- file.path(opt$outDir, paste(opt$runName, 'cg_calls_pass_perc_1.box1.pdf', sep='_') )
  gg <- ggplot2::ggplot(data=plot_tib, aes(Sample_Per1, cg_calls_pass_perc_1)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per2)) + ggtitle("Infinium II Bins")
  ggplot2::ggsave(box_pdf1,gg)
  
  box_pdf2 <- file.path(opt$outDir, paste(opt$runName, 'cg_calls_pass_perc_1.box2.pdf', sep='_') )
  gg <- ggplot2::ggplot(data=plot_tib, aes(Sample_Per2, cg_calls_pass_perc_1)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per1))  + ggtitle("Infinium I Bins")
  ggplot2::ggsave(box_pdf2,gg)
  
  
  # Metrics::
  #  AutoSample_dB_Val_2
  box_pdf1 <- file.path(opt$outDir, paste(opt$runName, 'AutoSample_R2_Val_2.box1.pdf', sep='_') )
  gg <- ggplot2::ggplot(data=plot_tib, aes(Sample_Per1, AutoSample_R2_Val_2)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per2)) + ggtitle("Infinium II Bins")
  ggplot2::ggsave(box_pdf1,gg)
  
  box_pdf2 <- file.path(opt$outDir, paste(opt$runName, 'AutoSample_R2_Val_2.box2.pdf', sep='_') )
  gg <- ggplot2::ggplot(data=plot_tib, aes(Sample_Per2, AutoSample_R2_Val_2)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per1))  + ggtitle("Infinium I Bins")
  ggplot2::ggsave(box_pdf2,gg)
  
  
  
  
  
  
  
  # Metrics::
  #  AutoSample_dB_Val_2
  ggplot2::ggplot(data=plot_tib, aes(Sample_Per1, AutoSample_dB_Val_2)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per2))
  
  ggplot2::ggplot(data=plot_tib, aes(Sample_Per2, AutoSample_dB_Val_2)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per1))
  
  # Metrics::
  #  AutoSample_R2_Val_2
  ggplot2::ggplot(data=plot_tib, aes(Sample_Per1, AutoSample_R2_Val_2)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per2))
  
  ggplot2::ggplot(data=plot_tib, aes(Sample_Per2, AutoSample_R2_Val_2)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per1))
  
  # Metrics::
  #  cg_calls_pass_perc_1
  ggplot2::ggplot(data=plot_tib, aes(Sample_Per1, cg_calls_pass_perc_1)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per2))
  
  ggplot2::ggplot(data=plot_tib, aes(Sample_Per2, cg_calls_pass_perc_1)) +
    ggplot2::geom_boxplot() + 
    ggplot2::facet_grid(cols = vars(Sample_Per1))
  
  
  #
  # Density Plots:: AgeSkinBlood_2
  #
  ggplot2::ggplot(data=plot_tib, aes(AgeSkinBlood_2, group=Sample_Per1)) + 
    ggplot2::geom_density(alpha=0.3, aes(fill=Sample_Per1)) +
    ggplot2::facet_grid(cols = vars(Sample_Per2))
  
  ggplot2::ggplot(data=plot_tib, aes(AgeSkinBlood_2, group=Sample_Per2)) + 
    ggplot2::geom_density(alpha=0.3, aes(fill=Sample_Per2)) +
    ggplot2::facet_grid(cols = vars(Sample_Per1))
  
  
  
  
  #  ggplot2::facet_grid(rows = vars(Sample_Per1), cols = vars(Sample_Per2))
  
  ggplot2::ggplot(data=plot_tib, aes(cg_calls_pass_perc_1, AgeSkinBlood_2, color=Sample_Per1)) +
    ggplot2::geom_point() +
    ggplot2::geom_density_2d() +
    ggplot2::facet_grid(cols=vars(Sample_Per2))
  
  ggplot2::ggplot(data=plot_tib, aes(cg_calls_pass_perc_1, AgeSkinBlood_2, color=Sample_Per2)) +
    ggplot2::geom_point() +
    ggplot2::geom_density_2d() +
    ggplot2::facet_grid(cols=vars(Sample_Per1))
  
  
  #  ggplot2::facet_grid(cols=vars(Sample_Per1), rows=vars(Sample_Per2))
  
  
  
  auto_ss_sum <- auto_ss_tib %>% 
    dplyr::group_by(auto_ss_tib$detect_manifest,AutoSample_R2_Key_2) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  auto_ss_sum %>% print(n=base::nrow(auto_ss_sum))
  
  ggplot2::ggplot(data=auto_ss_tib , aes(x=Sample_Perc, y=Sample_Rand, color=AutoSample_dB_Key_2)) +
    ggplot2::geom_boxplot()
  
  ggplot2::ggplot(data=auto_ss_tib %>% dplyr::arrange(Sample_Per1,Sample_Per2), 
                  aes(x=Sample_Per1, y=Sample_Per2, color=AutoSample_dB_Val_2)) +
    # ggplot2::geom_point() +
    ggplot2::geom_density_2d()
  
  
  
  ggplot2::ggplot(data=auto_ss_tib , aes(x=Sample_Rand, y=AutoSample_dB_Val_2, group=Sample_Perc)) +
    ggplot2::geom_boxplot()
  
  ggplot2::ggplot(data=auto_ss_tib %>% dplyr::arrange(Sample_Per1,Sample_Per2), aes(Sample_Per1, y=AutoSample_dB_Val_2, group=Sample_Rand)) +
    ggplot2::geom_boxplot()
  
  ggplot2::ggplot(data=auto_ss_tib %>% dplyr::arrange(Sample_Per1,Sample_Per2), aes(Sample_Per1, y=cg_calls_pass_perc_1, color=Sample_Rand)) +
    ggplot2::geom_boxplot() + ggplot2::facet_grid(rows = vars(Sample_Per2))
  
  ggplot2::ggplot(data=auto_ss_tib %>% dplyr::arrange(Sample_Per1,Sample_Per2), aes(AutoSample_dB_Val_2, y=cg_calls_pass_perc_1, color=Sample_Rand)) +
    ggplot2::geom_boxplot() + ggplot2::facet_grid(rows = vars(Sample_Per1), cols = vars(Sample_Per2))
  
  ggplot2::ggplot(data=auto_ss_tib %>% dplyr::arrange(Sample_Per1,Sample_Per2), aes(AutoSample_dB_Val_2, y=cg_calls_pass_perc_1)) +
    ggplot2::geom_boxplot() + ggplot2::facet_grid(rows = vars(Sample_Per1), cols = vars(Sample_Per2))
  
  
  
  ggplot2::ggplot(data=auto_ss_tib %>% dplyr::arrange(Sample_Per1,Sample_Per2)) +
    ggplot2::geom_density(aes(Sample_Per1), fill="red") +
    ggplot2::geom_density(aes(Sample_Per2), fill="blue")
  
  
  #  ggplot2::facet_grid(rows = vars(Sample_Per1), cols = vars(Sample_Per2))
  
}

#
# Quick fix update calls path for noob_sum::
#
par$chig_sub <- TRUE
par$chig_sub <- FALSE
if (par$noob_sub || par$chig_sub) {
  
  auto_ss_tib <- auto_ss_tib %>% 
    dplyr::mutate(Calls_Path=stringr::str_replace(SampleSheet_Path, "_AutoSampleSheet.csv.gz", "_ind.call.dat.csv.gz"))
  
  # Qucik Validation of Calls Path Update for noob_sum::
  auto_ss_tib %>% 
    dplyr::distinct(Calls_Path) %>% head(n=1) %>% 
    dplyr::pull(Calls_Path) %>% file.exists()
  auto_ss_tib %>% 
    dplyr::distinct(Calls_Path) %>% 
    base::nrow() %>% print()
  
  #
  # Quick FIx for noob_sub
  #
  hum_ss_tib <- auto_ss_tib %>% 
    # dplyr::group_by(Sentrix_Name,detect_platform) %>%
    dplyr::group_by(Sentrix_Name,!!class_var) %>%
    dplyr::mutate(Sample_Class=dplyr::row_number()) %>%
    dplyr::select(Sentrix_Name,Sample_Class,detect_manifest,detect_platform,detect_version) %>%
    dplyr::ungroup() %>%
    # dplyr::arrange(detect_version)
    dplyr::arrange(Sample_Class)
  
  # opt$classVar <- "Sample_Class"
  class_var <- rlang::sym("Sample_Class")
  
}

labs_ss_tib <- NULL
if (!is.null(hum_ss_tib)) {
  cat(glue::glue("[{par$prgmTag}]: Adding default class_var='{opt$classVar}'{RET}") )
  
  if (par$noob_sub) {
    labs_ss_tib <- auto_ss_tib %>% 
      dplyr::left_join(hum_ss_tib, by=c("Sentrix_Name","detect_manifest","detect_platform","detect_version") ) %>% 
      dplyr::arrange(!!class_var) %>%
      dplyr::mutate(Sentrix_Uniq=paste(Sentrix_Name,detect_manifest, sep='-'))
  } else {
    labs_ss_tib <- auto_ss_tib %>% 
      dplyr::left_join(hum_ss_tib, by="Sentrix_Name") %>% 
      dplyr::arrange(!!class_var)
  }
  
} else if (!is.null(opt$sampleCsv) && file.exists(opt$sampleCsv)) {
  
  cat(glue::glue("[{par$prgmTag}]: Loading predfined human classification; sampleCsv='{opt$sampleCsv}'{RET}") )
  
  hum_ss_tib <- suppressMessages(suppressWarnings( readr::read_csv(opt$sampleCsv) )) %>% 
    dplyr::mutate(Sentrix_Name=paste(Sentrix_ID,Sentrix_Position, sep='_')) %>% 
    dplyr::select(Sentrix_Name,Sample_Class, dplyr::everything())
  hum_ss_len <- hum_ss_tib %>% base::nrow()
  cat(glue::glue("[{par$prgmTag}]: Done. Loading Human Classification; hum_ss_len={hum_ss_len}!{RET}{RET}") )
  # print(hum_ss_tib)
  
  # Left Join now that we will force Sample_Class to nSARSCov2 (COVID-) below
  # labs_ss_tib <- auto_ss_tib %>% dplyr::inner_join(hum_ss_tib, by="Sentrix_Name") %>% dplyr::arrange(!!class_var)
  #
  if (opt$select) {
    labs_ss_tib <- auto_ss_tib %>% dplyr::inner_join(hum_ss_tib, by="Sentrix_Name") %>% dplyr::arrange(!!class_var)
  } else {
    labs_ss_tib <- auto_ss_tib %>% dplyr::left_join(hum_ss_tib, by="Sentrix_Name") %>% dplyr::arrange(!!class_var)
  }
} else {
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]: Using Auto Classification; classVar='{opt$classVar}'{RET}") )
  
  labs_ss_tib <- auto_ss_tib
}

labs_ss_len <- labs_ss_tib %>% base::nrow()
if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Done. Joining Human Classification ",
                 "and Auto Sample Sheets; Total={labs_ss_len}.{RET}{RET}"))
print_tib(labs_ss_tib,par$prgmTag, opt$verbose,vt=10,tc=1, n="labs_ss_tib")

# Write Auto Sample Sheet::
if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Writing labs_csv={opt$labs_csv}.{RET}{RET}"))
readr::write_csv(labs_ss_tib, opt$labs_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Import Datasets (Calls)::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par$outName <- opt$runName
if (!is.null(opt$platform)) par$outName <- paste(par$outName, opt$platform, sep='_')
if (!is.null(opt$version))  par$outName <- paste(par$outName, opt$version, sep='_')
if (!is.null(opt$workflow)) par$outName <- paste(par$outName, opt$workflow, sep='_')

opt$beg_file <- file.path(opt$outDir, paste(par$outName, 'load-beg.txt', sep='_') )
opt$end_file <- file.path(opt$outDir, paste(par$outName, 'load-end.txt', sep='_') )
opt$labs_csv <- file.path(opt$outDir, paste(par$outName, 'AutoSampleSheet.csv.gz', sep='_') )
opt$call_csv <- file.path(opt$outDir, paste(par$outName, 'MergedDataFiles.tib.csv.gz', sep='_') )
opt$sset_csv <- file.path(opt$outDir, paste(par$outName, 'MergedSsetFiles.tib.csv.gz', sep='_') )

if (opt$clean) unlink( list.files(opt$outDir, full.names=TRUE) )

pass_time_check <- TRUE
if (file.exists(opt$beg_file) && file.exists(opt$end_file) && file.exists(opt$call_csv)) {
  beg_date <- file.mtime(opt$beg_file)
  end_date <- file.mtime(opt$end_file)
  rawFiles <- list.files(opt$outDir, pattern='raw-data.csv.gz$', full.names=TRUE)
  
  if (length(rawFiles)==0)
    stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Failed to find ANY 'raw-data.csv.gz' files!!!{RET}{RET}"))
  rawFiles <- c(rawFiles,opt$call_csv)
  
  for (file in rawFiles) {
    if (file.mtime(file) - file.mtime(opt$beg_file) > 0) pass_time_check >- FALSE
    if (file.mtime(file) - file.mtime(opt$end_file) > 0) pass_time_check <- FALSE
  }
} else {
  pass_time_check <- FALSE
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Write Outputs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

chipName <- "Sentrix_Name"
if (!pass_time_check) {
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]: Building from scratch...{RET}"))
  unlink( list.files(opt$outDir, full.names=TRUE) )
  
  cmd <- paste('touch',opt$beg_file, sep=' ')
  system(cmd)
  
  # Write Merged Call Files Tibble::
  call_file_tib <- NULL
  if (opt$addPathsCall) {
    call_file_tib <- mergeCallsFromSS(
      ss=labs_ss_tib, max=par$maxTest,
      outName=par$outName, outDir=opt$outDir, 
      chipName=chipName, pathName="Calls_Path", 
      joinNameA="Probe_ID", joinNameB=NULL,
      joinType = opt$joinType,
      verbose=opt$verbose, vt=1, tc=1, tt=pTracker)
    
    safe_write(call_file_tib,"csv",opt$call_csv,verbose=opt$verbose,tt=pTracker)
  }
  

  # Write Merged Sset Files Tibble:: 
  grns_tib <- NULL
  reds_tib <- NULL
  sset_file_tib <- NULL
  if (opt$addPathsSset) {
    # New method for suffix = 'ind.sigs.dat.csv.gz'
    sset_file_tib <- mergeCallsFromSS(
      ss=labs_ss_tib, max=par$maxTest, pre=NULL,
      outName=par$outName, outDir=opt$outDir,
      chipName=chipName, pathName="Ssets_Path", 
      joinNameA="Probe_ID", joinNameB="Probe_Design",
      joinType = opt$joinType,
      verbose=opt$verbose, vt=1, tc=1, tt=pTracker)
    
    safe_write(sset_file_tib,"csv",opt$sset_csv,verbose=opt$verbose,tt=pTracker)
  }
  
  # Only need this for suffix = 'idat.sigs.csv.gz'
  if (FALSE) {
    beta_tib <- NULL
    address_tib <- NULL
    if (!is.null(src_man_tib)) {
      beta_tib <- loadFromFileTib(
        tib=call_file_tib, type="betas", key="Method",
        verbose=opt$verbose, tt=pTracker)
      
      address_tib <- src_man_tib %>%
        dplyr::filter(Probe_ID %in% beta_tib$Probe_ID) %>% 
        dplyr::select(Probe_Design,M,U) %>% tidyr::gather(Probe_Type, Address, -Probe_Design) %>% 
        dplyr::filter(!is.na(Address)) %>% 
        dplyr::arrange(Address) %>%
        clean_tibble()
    }
    
    # Old method for suffix = 'idat.sigs.csv.gz'
    sset_file_tib <- mergeCallsFromSS(
      ss=labs_ss_tib, max=par$maxTest, pre=address_tib,
      outName=par$outName, outDir=opt$outDir,
      chipName=chipName, pathName="Ssets_Path",
      joinNameA="Address", joinNameB=NULL,
      joinType = opt$joinType,
      verbose=opt$verbose, vt=1, tc=1, tt=pTracker)
    
    grns_tib <- loadFromFileTib(
      tib=sset_file_tib, type="Raw_Grn_sig", key="Method",
      verbose=opt$verbose, tt=pTracker)
    reds_tib <- loadFromFileTib(
      tib=sset_file_tib, type="Raw_Red_sig", key="Method",
      verbose=opt$verbose, tt=pTracker)
  }
  
  cmd <- paste('touch',opt$end_file, sep=' ')
  system(cmd)
  
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]: Done. Building from scratch.{RET}{RET}"))
} else {
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]: Build Already Up to Date.{RET}{RET}"))
}

if (FALSE) {
  labs_tib <- labs_ss_tib %>% 
    dplyr::select(Sentrix_Name, dplyr::starts_with("Sample")) %>% 
    dplyr::select(-dplyr::ends_with("Path"))
  
  beta_tib <- NULL
  sigU_tib <- NULL
  sigM_tib <- NULL
  
  beta_tib <- loadFromFileTib(
    tib=call_file_tib, type="betas", key="Method",
    verbose=opt$verbose, tt=pTracker)
  dBrf_tib <- loadFromFileTib(
    tib=call_file_tib, type="dB_ref", key="Method",
    verbose=opt$verbose, tt=pTracker)
  
  sigU_tib <- loadFromFileTib(
    tib=sset_file_tib, type="sig_U", key="Method",
    verbose=opt$verbose, tt=pTracker)
  sigM_tib <- loadFromFileTib(
    tib=sset_file_tib, type="sig_M", key="Method",
    verbose=opt$verbose, tt=pTracker)
  
  
  beta_tab <- beta_tib %>% 
    # dplyr::select(-Probe_ID) %>% 
    tidyr::gather(Sentrix_Name, beta, -Probe_ID) %>%
    dplyr::inner_join(labs_tib, by="Sentrix_Name")

  dBrf_tab <- dBrf_tib %>% 
    # dplyr::select(-Probe_ID) %>% 
    tidyr::gather(Sentrix_Name, dB, -Probe_ID) %>%
    dplyr::inner_join(labs_tib, by="Sentrix_Name")
  
  #
  # CORRECT PLOTTING BELOW FOR BOXPLOT::
  #
  ggplot2::ggplot(data=dBrf_tab, aes(Sample_Perc,log(dB))) +
    ggplot2::geom_boxplot()

  ggplot2::ggplot(data=dBrf_tab, aes(factor(Sample_Per1),log(dB))) +
    ggplot2::geom_boxplot()
  
  ggplot2::ggplot(data=dBrf_tab, aes(factor(Sample_Per2),log(dB))) +
    ggplot2::geom_boxplot()
  
  #
  # These are similar to what worked before...
  #
  ggplot2::ggplot(data=beta_tab) +
    ggplot2::geom_boxplot(aes(factor(Sample_Per1), beta)) +
    ggplot2::facet_grid(cols = vars(Sample_Per2))
    
  ggplot2::ggplot(data=beta_tab) +
    ggplot2::geom_boxplot(aes(factor(Sample_Per2), beta)) +
    ggplot2::facet_grid(cols = vars(Sample_Per1))
  
  
  ggplot2::ggplot(data=beta_tab, aes(x=beta, color=factor(Sample_Perc))) +
    ggplot2::geom_density()
  ggplot2::ggplot(data=beta_tab, aes(x=beta, color=factor(Sample_Per1_Str))) +
    ggplot2::geom_density()
  ggplot2::ggplot(data=beta_tab, aes(x=beta, color=factor(Sample_Per2_Str))) +
    ggplot2::geom_density()
  ggplot2::ggplot(data=beta_tab, aes(x=beta, color=factor(Sample_Rand))) +
    ggplot2::geom_density()
  
  # ggplot2::ggplot(data=beta_plot_tab, aes(x=beta, color=Sample_Perc, group=Sample_Rand)) +
  #   ggplot2::geom_density() +
  #   ggplot2::facet_grid(rows = vars(Sample_Per1_Str), cols = vars(Sample_Per2_Str))
  
  sigs_tab <- dplyr::inner_join(
    sigU_tib %>% tidyr::gather(Sentrix_Name, sigU, -Probe_ID, -Probe_Design),
    sigM_tib %>% tidyr::gather(Sentrix_Name, sigM, -Probe_ID, -Probe_Design),
    by=c("Sentrix_Name","Probe_ID","Probe_Design")
  ) %>%
    dplyr::inner_join(labs_tib, by="Sentrix_Name")

  ggplot2::ggplot(data=sigs_tab) +
    ggplot2::geom_density(aes(x=sigU, color=factor(Sample_Perc))) +
    ggplot2::geom_density(aes(x=sigM, color=factor(Sample_Perc)))
  
  ggplot2::ggplot(data=sigs_tab, aes(x=sigU, y=sigM, color=factor(Sample_Per1_Str))) +
    # ggplot2::geom_point() +
    ggplot2::geom_density_2d()

  
  
  # Only need this for suffix = 'idat.sigs.csv.gz'
  if (FALSE) {
    grns_plot_tab <- grns_tib %>% 
      tidyr::gather(Sentrix_Name, Grn, -Address) %>%
      dplyr::inner_join(address_tib, by="Address") %>%
      dplyr::inner_join(labs_tib, by="Sentrix_Name")
    
    ggplot2::ggplot(data=grns_plot_tab, aes(x=Grn, color=Probe_Type, group=Probe_Design)) +
      ggplot2::geom_density()
    
    ggplot2::ggplot(data=grns_plot_tab, aes(x=Grn, color=Sample_Perc, group=Sentrix_Name)) +
      ggplot2::geom_density()
    
    ggplot2::ggplot(data=grns_plot_tab, aes(x=Grn, color=Sentrix_Name, group=Sample_Perc)) +
      ggplot2::geom_density()
  }
}




if (par$noob_sub) {
  
  # NOTE: The commented code below only looks at Infinium I probe performance::
  #  TBD:: Probably should still write this out, but for both I/II designs...
  #
  # labs_ss_tib %>% 
  #   dplyr::mutate(detect_version_int=stringr::str_remove(detect_version, "^s") %>% as.integer()) %>%
  #   dplyr::group_by(detect_version_int) %>%
  #   dplyr::arrange(detect_version_int) %>%
  #   dplyr::summarise(Min_PV=min(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE), 
  #                    Max_PV=max(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
  #                    Avg_PV=mean(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE), 
  #                    Med_PV=median(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
  #                    Std_PV=sd(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
  #                    TotalN=n(),
  #                    Min_R2=min(AutoSample_R2_1_Val_1, na.rm=TRUE), 
  #                    Max_R2=max(AutoSample_R2_1_Val_1, na.rm=TRUE),
  #                    Avg_R2=mean(AutoSample_R2_1_Val_1, na.rm=TRUE), 
  #                    Med_R2=median(AutoSample_R2_1_Val_1, na.rm=TRUE),
  #                    Std_R2=sd(AutoSample_R2_1_Val_1, na.rm=TRUE),
  #                    
  #                    Min_dB=min(AutoSample_dB_1_Val_1, na.rm=TRUE), 
  #                    Max_dB=max(AutoSample_dB_1_Val_1, na.rm=TRUE),
  #                    Avg_dB=mean(AutoSample_dB_1_Val_1, na.rm=TRUE), 
  #                    Med_dB=median(AutoSample_dB_1_Val_1, na.rm=TRUE),
  #                    Std_dB=sd(AutoSample_dB_1_Val_1, na.rm=TRUE),
  #                    
  #                    Min_Age=min(AgeSkinBlood_1, na.rm=TRUE), 
  #                    Max_Age=max(AgeSkinBlood_1, na.rm=TRUE),
  #                    Avg_Age=mean(AgeSkinBlood_1, na.rm=TRUE), 
  #                    Med_Age=median(AgeSkinBlood_1, na.rm=TRUE),
  #                    Std_Age=sd(AgeSkinBlood_1, na.rm=TRUE),
  #                    .groups="drop"
  #   )
  # 
  # labs_ss_tib %>% 
  #   dplyr::mutate(detect_version_int=stringr::str_remove(detect_version, "^s") %>% as.integer()) %>%
  #   dplyr::group_by(detect_version_int) %>%
  #   dplyr::arrange(detect_version_int) %>%
  #   dplyr::summarise(Min_PV=min(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE), 
  #                    Max_PV=max(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
  #                    Avg_PV=mean(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE), 
  #                    Med_PV=median(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
  #                    Std_PV=sd(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
  #                    TotalN=n(),
  #                    Min_R2=min(AutoSample_R2_1_Val_2, na.rm=TRUE), 
  #                    Max_R2=max(AutoSample_R2_1_Val_2, na.rm=TRUE),
  #                    Avg_R2=mean(AutoSample_R2_1_Val_2, na.rm=TRUE), 
  #                    Med_R2=median(AutoSample_R2_1_Val_2, na.rm=TRUE),
  #                    Std_R2=sd(AutoSample_R2_1_Val_2, na.rm=TRUE),
  #                    
  #                    Min_dB=min(AutoSample_dB_1_Val_2, na.rm=TRUE), 
  #                    Max_dB=max(AutoSample_dB_1_Val_2, na.rm=TRUE),
  #                    Avg_dB=mean(AutoSample_dB_1_Val_2, na.rm=TRUE), 
  #                    Med_dB=median(AutoSample_dB_1_Val_2, na.rm=TRUE),
  #                    Std_dB=sd(AutoSample_dB_1_Val_2, na.rm=TRUE),
  #                    
  #                    Min_Age=min(AgeSkinBlood_2, na.rm=TRUE), 
  #                    Max_Age=max(AgeSkinBlood_2, na.rm=TRUE),
  #                    Avg_Age=mean(AgeSkinBlood_2, na.rm=TRUE), 
  #                    Med_Age=median(AgeSkinBlood_2, na.rm=TRUE),
  #                    Std_Age=sd(AgeSkinBlood_2, na.rm=TRUE),
  #                    .groups="drop"
  #   )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Write Outputs:: detect_version
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opt$vers_sum1_csv <- file.path(opt$outDir, paste(par$outName, 'detect_version-summary1.csv.gz', sep='_') )
  vers_sum1_tib <- labs_ss_tib %>% 
    dplyr::mutate(detect_version_int=detect_version %>%
                    stringr::str_remove("^s") %>%
                    stringr::str_remove("^S") %>% 
                    as.double()) %>%
    dplyr::group_by(detect_version_int) %>%
    dplyr::arrange(detect_version_int) %>%
    dplyr::summarise(Min_PV=min(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE), 
                     Max_PV=max(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     Avg_PV=mean(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE), 
                     Med_PV=median(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     Std_PV=sd(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     TotalN=n(),
                     Min_R2=min(AutoSample_R2_Val_1, na.rm=TRUE), 
                     Max_R2=max(AutoSample_R2_Val_1, na.rm=TRUE),
                     Avg_R2=mean(AutoSample_R2_Val_1, na.rm=TRUE), 
                     Med_R2=median(AutoSample_R2_Val_1, na.rm=TRUE),
                     Std_R2=sd(AutoSample_R2_Val_1, na.rm=TRUE),
                     
                     Min_dB=min(AutoSample_dB_Val_1, na.rm=TRUE), 
                     Max_dB=max(AutoSample_dB_Val_1, na.rm=TRUE),
                     Avg_dB=mean(AutoSample_dB_Val_1, na.rm=TRUE), 
                     Med_dB=median(AutoSample_dB_Val_1, na.rm=TRUE),
                     Std_dB=sd(AutoSample_dB_Val_1, na.rm=TRUE),
                     
                     Min_Age=min(AgeSkinBlood_1, na.rm=TRUE), 
                     Max_Age=max(AgeSkinBlood_1, na.rm=TRUE),
                     Avg_Age=mean(AgeSkinBlood_1, na.rm=TRUE), 
                     Med_Age=median(AgeSkinBlood_1, na.rm=TRUE),
                     Std_Age=sd(AgeSkinBlood_1, na.rm=TRUE),
                     .groups="drop"
    ) %>% dplyr::mutate_if(is.double, base::round, 4)
  readr::write_csv(vers_sum1_tib, opt$vers_sum1_csv)
  vers_sum1_tib %>% print(n=base::nrow(vers_sum1_tib))
  
  opt$vers_sum2_csv <- file.path(opt$outDir, paste(par$outName, 'detect_version-summary2.csv.gz', sep='_') )
  vers_sum2_tib <- labs_ss_tib %>% 
    dplyr::mutate(detect_version_int=detect_version %>%
                    stringr::str_remove("^s") %>%
                    stringr::str_remove("^S") %>% 
                    as.double()) %>%
    dplyr::group_by(detect_version_int) %>%
    dplyr::arrange(detect_version_int) %>%
    dplyr::summarise(Min_PV=min(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE), 
                     Max_PV=max(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     Avg_PV=mean(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE), 
                     Med_PV=median(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     Std_PV=sd(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     TotalN=n(),
                     Min_R2=min(AutoSample_R2_Val_2, na.rm=TRUE), 
                     Max_R2=max(AutoSample_R2_Val_2, na.rm=TRUE),
                     Avg_R2=mean(AutoSample_R2_Val_2, na.rm=TRUE), 
                     Med_R2=median(AutoSample_R2_Val_2, na.rm=TRUE),
                     Std_R2=sd(AutoSample_R2_Val_2, na.rm=TRUE),
                     
                     Min_dB=min(AutoSample_dB_Val_2, na.rm=TRUE), 
                     Max_dB=max(AutoSample_dB_Val_2, na.rm=TRUE),
                     Avg_dB=mean(AutoSample_dB_Val_2, na.rm=TRUE), 
                     Med_dB=median(AutoSample_dB_Val_2, na.rm=TRUE),
                     Std_dB=sd(AutoSample_dB_Val_2, na.rm=TRUE),
                     
                     Min_Age=min(AgeSkinBlood_2, na.rm=TRUE), 
                     Max_Age=max(AgeSkinBlood_2, na.rm=TRUE),
                     Avg_Age=mean(AgeSkinBlood_2, na.rm=TRUE), 
                     Med_Age=median(AgeSkinBlood_2, na.rm=TRUE),
                     Std_Age=sd(AgeSkinBlood_2, na.rm=TRUE),
                     .groups="drop"
    ) %>% dplyr::mutate_if(is.double, base::round, 4)
  readr::write_csv(vers_sum2_tib, opt$vers_sum2_csv)
  vers_sum2_tib %>% print(n=base::nrow(vers_sum2_tib))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Write Outputs:: detect_platform
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opt$plat_sum1_csv <- file.path(opt$outDir, paste(par$outName, 'detect_platform-summary1.csv.gz', sep='_') )
  plat_sum1_tib <- labs_ss_tib %>% 
    dplyr::mutate(detect_platform_int=stringr::str_remove(detect_platform, "^Rand") %>% as.integer()) %>%
    dplyr::group_by(detect_platform_int) %>%
    dplyr::arrange(detect_platform_int) %>%
    dplyr::summarise(Min_PV=min(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE), 
                     Max_PV=max(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     Avg_PV=mean(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE), 
                     Med_PV=median(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     Std_PV=sd(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     TotalN=n(),
                     Min_R2=min(AutoSample_R2_Val_1, na.rm=TRUE), 
                     Max_R2=max(AutoSample_R2_Val_1, na.rm=TRUE),
                     Avg_R2=mean(AutoSample_R2_Val_1, na.rm=TRUE), 
                     Med_R2=median(AutoSample_R2_Val_1, na.rm=TRUE),
                     Std_R2=sd(AutoSample_R2_Val_1, na.rm=TRUE),
                     
                     Min_dB=min(AutoSample_dB_Val_1, na.rm=TRUE), 
                     Max_dB=max(AutoSample_dB_Val_1, na.rm=TRUE),
                     Avg_dB=mean(AutoSample_dB_Val_1, na.rm=TRUE), 
                     Med_dB=median(AutoSample_dB_Val_1, na.rm=TRUE),
                     Std_dB=sd(AutoSample_dB_Val_1, na.rm=TRUE),
                     
                     Min_Age=min(AgeSkinBlood_1, na.rm=TRUE), 
                     Max_Age=max(AgeSkinBlood_1, na.rm=TRUE),
                     Avg_Age=mean(AgeSkinBlood_1, na.rm=TRUE), 
                     Med_Age=median(AgeSkinBlood_1, na.rm=TRUE),
                     Std_Age=sd(AgeSkinBlood_1, na.rm=TRUE),
                     .groups="drop"
    ) %>% dplyr::mutate_if(is.double, base::round, 4)
  readr::write_csv(plat_sum1_tib, opt$plat_sum1_csv)
  plat_sum1_tib %>% print(n=base::nrow(plat_sum1_tib))
  
  opt$plat_sum2_csv <- file.path(opt$outDir, paste(par$outName, 'detect_platform-summary2.csv.gz', sep='_') )
  plat_sum2_tib <- labs_ss_tib %>% 
    dplyr::mutate(detect_platform_int=stringr::str_remove(detect_platform, "^Rand") %>% as.integer()) %>%
    dplyr::group_by(detect_platform_int) %>%
    dplyr::arrange(detect_platform_int) %>%
    dplyr::summarise(Min_PV=min(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE), 
                     Max_PV=max(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     Avg_PV=mean(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE), 
                     Med_PV=median(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     Std_PV=sd(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     TotalN=n(),
                     Min_R2=min(AutoSample_R2_Val_2, na.rm=TRUE), 
                     Max_R2=max(AutoSample_R2_Val_2, na.rm=TRUE),
                     Avg_R2=mean(AutoSample_R2_Val_2, na.rm=TRUE), 
                     Med_R2=median(AutoSample_R2_Val_2, na.rm=TRUE),
                     Std_R2=sd(AutoSample_R2_Val_2, na.rm=TRUE),
                     
                     Min_dB=min(AutoSample_dB_Val_2, na.rm=TRUE), 
                     Max_dB=max(AutoSample_dB_Val_2, na.rm=TRUE),
                     Avg_dB=mean(AutoSample_dB_Val_2, na.rm=TRUE), 
                     Med_dB=median(AutoSample_dB_Val_2, na.rm=TRUE),
                     Std_dB=sd(AutoSample_dB_Val_2, na.rm=TRUE),
                     
                     Min_Age=min(AgeSkinBlood_2, na.rm=TRUE), 
                     Max_Age=max(AgeSkinBlood_2, na.rm=TRUE),
                     Avg_Age=mean(AgeSkinBlood_2, na.rm=TRUE), 
                     Med_Age=median(AgeSkinBlood_2, na.rm=TRUE),
                     Std_Age=sd(AgeSkinBlood_2, na.rm=TRUE),
                     .groups="drop"
    ) %>% dplyr::mutate_if(is.double, base::round, 4)
  readr::write_csv(plat_sum2_tib, opt$plat_sum2_csv)
  plat_sum2_tib %>% print(n=base::nrow(plat_sum2_tib))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #           Write Outputs:: detect_version & detect_platform
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opt$both_sum1_csv <- file.path(opt$outDir, paste(par$outName, 'detect_both-summary1.csv.gz', sep='_') )
  both_sum1_tib <- labs_ss_tib %>% 
    dplyr::mutate(detect_version_int=detect_version %>%
                    stringr::str_remove("^s") %>%
                    stringr::str_remove("^S") %>% 
                    as.double(),
                  detect_platform_int=stringr::str_remove(detect_platform, "^Rand") %>% as.integer()) %>%
    dplyr::arrange(detect_version_int,detect_platform_int) %>%
    dplyr::group_by(detect_version_int,detect_platform_int) %>%
    dplyr::summarise(Min_PV=min(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE), 
                     Max_PV=max(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     Avg_PV=mean(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE), 
                     Med_PV=median(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     Std_PV=sd(cg_pvals_pOOBAH_pass_perc_1, na.rm=TRUE),
                     TotalN=n(),
                     Min_R2=min(AutoSample_R2_Val_1, na.rm=TRUE), 
                     Max_R2=max(AutoSample_R2_Val_1, na.rm=TRUE),
                     Avg_R2=mean(AutoSample_R2_Val_1, na.rm=TRUE), 
                     Med_R2=median(AutoSample_R2_Val_1, na.rm=TRUE),
                     Std_R2=sd(AutoSample_R2_Val_1, na.rm=TRUE),
                     
                     Min_dB=min(AutoSample_dB_Val_1, na.rm=TRUE), 
                     Max_dB=max(AutoSample_dB_Val_1, na.rm=TRUE),
                     Avg_dB=mean(AutoSample_dB_Val_1, na.rm=TRUE), 
                     Med_dB=median(AutoSample_dB_Val_1, na.rm=TRUE),
                     Std_dB=sd(AutoSample_dB_Val_1, na.rm=TRUE),
                     
                     Min_Age=min(AgeSkinBlood_1, na.rm=TRUE), 
                     Max_Age=max(AgeSkinBlood_1, na.rm=TRUE),
                     Avg_Age=mean(AgeSkinBlood_1, na.rm=TRUE), 
                     Med_Age=median(AgeSkinBlood_1, na.rm=TRUE),
                     Std_Age=sd(AgeSkinBlood_1, na.rm=TRUE),
                     .groups="drop"
    ) %>% dplyr::mutate_if(is.double, base::round, 4)
  readr::write_csv(both_sum1_tib, opt$both_sum1_csv)
  both_sum1_tib %>% print(n=base::nrow(both_sum1_tib))
  
  opt$both_sum2_csv <- file.path(opt$outDir, paste(par$outName, 'detect_both-summary2.csv.gz', sep='_') )
  both_sum2_tib <- labs_ss_tib %>% 
    dplyr::mutate(detect_version_int=detect_version %>%
                    stringr::str_remove("^s") %>%
                    stringr::str_remove("^S") %>% 
                    as.double(),
                  detect_platform_int=stringr::str_remove(detect_platform, "^Rand") %>% as.integer()) %>%
    dplyr::arrange(detect_version_int,detect_platform_int) %>%
    dplyr::group_by(detect_version_int,detect_platform_int) %>%
    dplyr::summarise(Min_PV=min(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE), 
                     Max_PV=max(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     Avg_PV=mean(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE), 
                     Med_PV=median(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     Std_PV=sd(cg_pvals_pOOBAH_pass_perc_2, na.rm=TRUE),
                     TotalN=n(),
                     Min_R2=min(AutoSample_R2_Val_2, na.rm=TRUE), 
                     Max_R2=max(AutoSample_R2_Val_2, na.rm=TRUE),
                     Avg_R2=mean(AutoSample_R2_Val_2, na.rm=TRUE), 
                     Med_R2=median(AutoSample_R2_Val_2, na.rm=TRUE),
                     Std_R2=sd(AutoSample_R2_Val_2, na.rm=TRUE),
                     
                     Min_dB=min(AutoSample_dB_Val_2, na.rm=TRUE), 
                     Max_dB=max(AutoSample_dB_Val_2, na.rm=TRUE),
                     Avg_dB=mean(AutoSample_dB_Val_2, na.rm=TRUE), 
                     Med_dB=median(AutoSample_dB_Val_2, na.rm=TRUE),
                     Std_dB=sd(AutoSample_dB_Val_2, na.rm=TRUE),
                     
                     Min_Age=min(AgeSkinBlood_2, na.rm=TRUE), 
                     Max_Age=max(AgeSkinBlood_2, na.rm=TRUE),
                     Avg_Age=mean(AgeSkinBlood_2, na.rm=TRUE), 
                     Med_Age=median(AgeSkinBlood_2, na.rm=TRUE),
                     Std_Age=sd(AgeSkinBlood_2, na.rm=TRUE),
                     .groups="drop"
    ) %>% dplyr::mutate_if(is.double, base::round, 4)
  readr::write_csv(both_sum2_tib, opt$both_sum2_csv)
  both_sum2_tib %>% print(n=base::nrow(both_sum2_tib))
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$opt_csv  <- file.path(opt$outDir, paste(par$prgmTag,'program-options.csv', sep='.') )
opt$par_csv  <- file.path(opt$outDir, paste(par$prgmTag,'program-parameters.csv', sep='.') )
opt$time_csv <- file.path(opt$outDir, paste(par$prgmTag,'time-tracker.csv.gz', sep='.') )

opt_tib  <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
par_tib  <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
time_tib <- pTracker$time %>% dplyr::mutate_if(is.numeric, list(round), 4)

readr::write_csv(opt_tib, opt$opt_csv)
readr::write_csv(par_tib, opt$par_csv)
readr::write_csv(time_tib, opt$time_csv)

sysTime <- Sys.time()
cat(glue::glue("{RET}[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
