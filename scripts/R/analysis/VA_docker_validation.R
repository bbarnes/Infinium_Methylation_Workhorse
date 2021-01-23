
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
par$prgmTag <- 'VA_docker_validation'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# Predefined human sample sheet name::
par$humanSampleSheetName <- 'humanSampleSheet.csv'

# File Based Parameters::
opt$inputsCsv <- NULL

# Directory Parameters::
opt$outDir    <- NULL
opt$buildDir  <- NULL

# Run Parameters::
opt$runName   <- NULL
opt$sampleCsv <- NULL
opt$findSampleSheet <- FALSE

# Class Parameters::
# Really simple test to make sure we can seperate the sexes...
opt$classVar <- 'Karyotype_0_Call'
opt$classVar <- 'Sample_Name'
opt$classVar <- 'Sample_Class'

opt$select <- FALSE
opt$clean_gta <- FALSE

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
  
  if (par$local_runType=='COVID') {
    
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
    
    opt$buildDir <- paste(
      file.path(par$topDir, 'scratch/COVID_All/swifthoof_main',opt$runName),
      sep=',')
    
  } else if (par$local_runType=='COVIC') {
    opt$runName  <- 'COVIC-Set7-06082020'
    opt$runName  <- 'COVIC-Set5-10062020'
    opt$runName  <- 'COVIC-Set1-15052020'
    
    opt$version  <- 'C0'
    
    opt$buildDir <- paste(
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
    
    opt$buildDir  <- paste(
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
    
    opt$mergeDir  <- paste(
      file.path(par$topDir,'scratch/merge_builds/LEGX/S1/Sample_Name',opt$runName),
      sep=',')
    
  } else if (par$local_runType=='qcMVP') {
    opt$classVar <- 'AutoSample_dB_Key'
    opt$classVar <- 'AutoSample_dB_Key_1'
    
    opt$platform <- 'EPIC'
    opt$version  <- 'B4'
    
    par$runName1  <- 'CNTL-Samples_VendA_10092020'
    par$runName2  <- 'CNTL-Samples_VendB_10092020'
    
    par$runName3  <- 'BETA-8x1-EPIC-Core'
    par$runName4  <- 'DELTA-8x1-EPIC-Core'
    par$runName5  <- 'DELTA-24x1-EPIC'
    par$runName6  <- 'BETA-8x1-EPIC-Bad'
    
    par$runName7  <- 'COVIC-Set1-15052020'
    par$runName8  <- 'COVIC-Set2-31052020'
    par$runName9  <- 'COVIC-Set3-05062020'
    par$runName10 <- 'COVIC-Set4-09062020'
    par$runName11 <- 'COVIC-Set5-10062020'
    par$runName12 <- 'COVIC-Set7-06082020'
    par$runName13 <- 'COVIC-Set8-26182020'
    
    par$runName14 <- 'Excalibur-Old-1609202'
    par$runName15 <- 'Excalibur-New-1609202'
    
    # opt$buildDir  <- paste(
    #   file.path(par$topDir,'scratch/swifthoof_main',par$runNameA),
    #   sep=',')
    # opt$runName  <- par$runNameA
    
    # TBD:: Need to turn this into a string and then parse...
    #
    # par$pass_vec <- c(par$runNameA,par$runNameB,par$runNameC,par$runNameD,
    #  par$runNameG)
    # par$fail_vec <- c(par$runNameE,par$runNameF,par$runNameH)
    pass_vec <- c(par$runName1,par$runName2,
                  par$runName3,par$runName4,par$runName5,
                  par$runName7,par$runName11,par$runName12,
                  par$runName15)
    
    opt$buildDir  <- paste(
      # Known passed chips::
      file.path(par$topDir,'scratch/swifthoof_main',par$runName1),
      file.path(par$topDir,'scratch/swifthoof_main',par$runName2),
      file.path(par$topDir,'scratch/swifthoof_main',par$runName3),
      file.path(par$topDir,'scratch/swifthoof_main',par$runName4),
      file.path(par$topDir,'scratch/swifthoof_main',par$runName5),
      
      # Known failed chips::
      file.path(par$topDir,'scratch/swifthoof_main',par$runName6),
      file.path(par$topDir,'scratch/swifthoof_main',par$runName7),
      file.path(par$topDir,'scratch/swifthoof_main',par$runName8),
      
      file.path(par$topDir,'scratch/swifthoof_main',par$runName9),
      file.path(par$topDir,'scratch/swifthoof_main',par$runName10),
      
      file.path(par$topDir,'scratch/swifthoof_main',par$runName11),
      file.path(par$topDir,'scratch/swifthoof_main',par$runName12),
      # file.path(par$topDir,'scratch/swifthoof_main',par$runName13),
      
      file.path(par$topDir,'scratch/swifthoof_main',par$runName14),
      file.path(par$topDir,'scratch/swifthoof_main',par$runName15),
      
      sep=',')
    
    opt$runName  <- paste(par$local_runType, sep='_')
    
    opt$addPathsCall <- TRUE
    opt$addPathsSset <- TRUE
    
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  # opt$outDir <- file.path(par$topDir, 'scratch', par$prgmTag)
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
    make_option(c("-b","--buildDirs"), type="character", default=opt$buildDir, 
                help="List of Build Directory [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--sampleCsv"), type="character", default=opt$sampleCsv, 
                help="Human provide sample sheet labeling [default= %default]", metavar="character"),
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
              'buildDir','runName','classVar',
              'addSampleName','addPathsCall','addPathsSset',
              'flagDetectPval','flagSampleDetect','flagRefMatch',
              'platform','version','percisionBeta','percisionPval',
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






#
# Plots Needed::
#
#   - Open vs. Raw
#   - Log(Pass)
#   - Log(Fail)
#





# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# VA data from v.4.3
prev_ss_csv <- '/Users/bretbarnes/Documents/data/transfer/VA-CNTL-AB-mergedAutoSampleSheet_v.4.3.csv.gz'
prev_ss_tib <- readr::read_csv(prev_ss_csv)

# VA data from cluster::
va18_ss_csv <- '/Users/bretbarnes/Documents/data/transfer/docker.v.1.8.AutoSampleSheet.csv.gz'

docker_dir <- '/Users/bretbarnes/Documents/data/VA_MVP/transfer'

va43_ss_csv <- file.path(docker_dir, 'data.v.4.3', 'docker.v.4.3.AutoSampleSheet.csv.gz')
va19_ss_csv <- file.path(docker_dir, 'data.v.1.9', 'docker.v.1.9.AutoSampleSheet.csv.gz')

va43_ss_tib <- readr::read_csv(va43_ss_csv)
va19_ss_tib <- readr::read_csv(va19_ss_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Validate Auto Sample Sheet Calls/PercPass::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

auto_ss_tib %>% 
  dplyr::select(Sentrix_Name, AutoSample_dB_1_Key_1, 
                cg_calls_pass_perc_1,cg_pvals_pOOBAH_pass_perc_1, 
                AutoSample_dB_Val_1, AutoSample_dB_Val_2 ) %>% 
  dplyr::filter(Sentrix_Name %in% prev_ss_tib$Sentrix_Name) %>%
  dplyr::mutate(Call_Diff=cg_calls_pass_perc_1-cg_pvals_pOOBAH_pass_perc_1) %>%
  dplyr::arrange(-AutoSample_dB_Val_1)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Validate Auto Sample Sheet::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# 203962710025_R01C01 = 
#   cg_calls_pass_perc cg_calls_pass_count cg_calls_total_count cg_calls_min_cutoff cg_calls_metric cg_calls_name
# <dbl>               <int>                <int> <chr>               <chr>           <chr>        
#   1               92.6              799276               862927 0.1                 pvals_pOOBAH    v.1.8
#

# Metrics::
#  - cg_pass_perc_basic_0  # Open Sesame
#  - cg_calls_pass_perc_1  # From Calls File Raw
#  - cg_pvals_pOOBAH_pass_perc_1 # From Calcs Raw
#
auto_df_tib <- auto_ss_tib %>% 
  dplyr::mutate(Call_Calc_Dif=cg_calls_pass_perc_1-cg_pvals_pOOBAH_pass_perc_1)

auto_df_tib <- auto_ss_tib %>% 
  dplyr::mutate(Call_Calc_Dif=cg_calls_pass_perc_1-cg_pass_perc_basic_0)

auto_df_tib <- auto_ss_tib %>% 
  dplyr::mutate(Call_Calc_Dif=cg_calls_pass_perc_1-cg_pvals_pOOBAH_pass_perc_2)

auto_df_tib <- auto_ss_tib %>% 
  dplyr::mutate(Call_Calc_Dif=cg_pvals_pOOBAH_pass_perc_1-cg_pvals_pOOBAH_pass_perc_2)

auto_df_tib$Call_Calc_Dif %>% unique


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Join VA Control & Auto Sample Sheet::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cntl_join_tib <- dplyr::left_join(
  dplyr::select(prev_ss_tib,Sentrix_Name,Requeue_Flag_Oob,Requeue_Flag_Neg,
                cg_1_pvals_pOOBAH_pass_perc,cg_2_pvals_pOOBAH_pass_perc),
  auto_ss_tib, 
  by="Sentrix_Name", suffix=c("_v43", "v_18")
)

dif_csv <- file.path('/Users/bretbarnes/Documents/tmp/dif-43-18.2.tib.csv')
dif_tib <- cntl_join_tib %>% 
  dplyr::mutate(dif1=cg_1_pvals_pOOBAH_pass_perc-cg_1_pass_perc_basic_0, 
                dif2=cg_2_pvals_pOOBAH_pass_perc-cg_2_pass_perc_basic_0) %>% 
  dplyr::arrange(-dif1) %>% 
  dplyr::select(Sentrix_Name,Requeue_Flag_Oob,AutoSample_dB_1_Key_1, dif1, dif2)
dif_tib %>% print(n=100)
readr::write_csv(dif_tib,dif_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Plot Join VA/Auto::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ggplot2::ggplot(data=cntl_join_tib) +
  ggplot2::geom_point(aes(x=cg_1_pvals_pOOBAH_pass_perc, y=cg_1_pass_perc_basic_0, color=Requeue_Flag_Oob)) +
  ggplot2::geom_point(aes(x=cg_2_pvals_pOOBAH_pass_perc, y=cg_2_pass_perc_basic_0, color=Requeue_Flag_Oob))

cntl_plot_tib <- dplyr::bind_rows(
  cntl_join_tib %>% dplyr::select(Sentrix_Name,cg_1_pvals_pOOBAH_pass_perc,cg_1_pass_perc_basic_0,Requeue_Flag_Oob,AutoSample_dB_1_Key_1) %>%
    purrr::set_names(c("Sentrix_Name","Pass_43","Pass_18","Requeue","Sample")) %>% dplyr::mutate(Design='I'),
  cntl_join_tib %>% dplyr::select(Sentrix_Name,cg_2_pvals_pOOBAH_pass_perc,cg_2_pass_perc_basic_0,Requeue_Flag_Oob,AutoSample_dB_1_Key_1) %>%
    purrr::set_names(c("Sentrix_Name","Pass_43","Pass_18","Requeue","Sample")) %>% dplyr::mutate(Design='II'),
)

cntl_plot_tib %>%
  dplyr::mutate(Requeue=as.character(Requeue)) %>%
  dplyr::group_by(Requeue,Design) %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x=Pass_43, y=Pass_18, color=Requeue)) +
  ggplot2::facet_grid(rows=vars(Design),
                      cols=vars(Sample) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Plot All Auto Sample Sheet:: Set Up
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pass_vec <- c("CNTL-Samples_VendA_10092020",
              "CNTL-Samples_VendB_10092020",
              # "BETA-8x1-EPIC-Core",
              "DELTA-8x1-EPIC-Core",
              # "BETA-8x1-EPIC-Bad",
              "COVIC-Set1-15052020",
              # "COVIC-Set2-31052020",
              # "COVIC-Set3-05062020",
              # "COVIC-Set4-09062020",
              "COVIC-Set5-10062020",
              "COVIC-Set7-06082020",
              "COVIC-Set8-26182020",
              # "Excalibur-Old-1609202",
              "Excalibur-New-1609202")

core_tib <- NULL
core_tib <- auto_ss_tib %>%
  dplyr::mutate(
    group_class=dplyr::case_when(build_source %in% pass_vec ~ 'P', TRUE ~ 'F'),
    build_source = 
      stringr::str_remove(build_source, '^CNTL-Samples_') %>% 
      stringr::str_remove('_[0-9]+$') %>% 
      stringr::str_remove('-[0-9]+$') %>%
      stringr::str_remove('-8x1-EPIC-Core$') %>%
      stringr::str_replace('-8x1-EPIC-Bad$', '_Bad') %>%
      stringr::str_replace('Excalibur-New', 'E') %>%
      stringr::str_replace('Excalibur-Old', 'E') %>%
      stringr::str_replace('COVID-Set', 'C') %>%
      stringr::str_replace('COVIC-Set', 'C') %>%
      stringr::str_replace('BETA', 'B') %>%
      stringr::str_replace('DELTA', 'D') %>%
      stringr::str_replace('Vend', 'V') %>%
      # stringr::str_remove('-Set') %>%
      stringr::str_replace_all('-','_')
  ) %>%
  dplyr::filter(! stringr::str_starts(AutoSample_dB_Key_2, 'T')) %>%
  dplyr::group_by(build_source) %>%
  dplyr::mutate(build_source_idx=dplyr::cur_group_id(),
                # build_source=paste(paste0(group_class,build_source_idx),build_source, sep='_'),
                build_source=paste(paste0(build_source_idx),build_source, sep='_'),
                Sample=AutoSample_dB_Key_2) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group_class,build_source_idx,Sample)

# core_tib %>% dplyr::group_by(cg_Failed_QC_basic_0,group_class,build_source_idx,Sample) %>% dplyr::summarise(Call_Cnt=n()) %>% print(n=1000)
# core_tib %>% dplyr::group_by(cg_Failed_QC_basic_0,group_class,build_source_idx) %>% dplyr::summarise(Call_Cnt=n()) %>% print(n=1000)


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Plot All Auto Sample Sheet:: New
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

plot_core_tib <- core_tib %>%
  # dplyr::filter(AutoSample_dB_Key_2 == AutoSample_R2_Key_2) %>%
  # dplyr::select(group_class,build_source,Sample, # build_source_idx,
  dplyr::select(group_class,build_source,Sample, # build_source_idx,
                AutoSample_R2_1_Val_2,AutoSample_R2_2_Val_2,
                AutoSample_dB_1_Val_2,AutoSample_dB_2_Val_2,
                cg_1_pass_perc_basic_0,
                cg_2_pass_perc_basic_0,
                cg_1_pvals_pOOBAH_pass_perc_1,
                cg_2_pvals_pOOBAH_pass_perc_1
  ) %>%
  purrr::set_names(stringr::str_remove_all(names(.), 'AutoSample_')) %>%
  purrr::set_names(stringr::str_replace_all(names(.), '_pvals_', '_')) %>%
  purrr::set_names(stringr::str_replace_all(names(.), '_pass_perc_basic_0', '_Open_pOOBAH')) %>%
  purrr::set_names(stringr::str_remove_all(names(.), '_pass_perc_1')) %>%
  purrr::set_names(stringr::str_remove_all(names(.), '_Val_2'))

plot_core_tib %>% dplyr::group_by(group_class,build_source) %>% 
  dplyr::summarise(Count=n(), .groups="drop")


#
# Direct comparison of Open vs. Raw Passing Percent::
#
dir_gg <- ggplot2::ggplot(data=plot_core_tib) + 
  ggplot2::geom_point(aes(x=cg_1_Open_pOOBAH, y=cg_1_pOOBAH, color=build_source)) + 
  ggplot2::geom_point(aes(x=cg_2_Open_pOOBAH, y=cg_2_pOOBAH, color=build_source)) +
  ggplot2::facet_grid(rows=vars(group_class),
                      cols=vars(Sample) )

dir_png <- file.path(opt$outDir, 'openSesameEPIC-vs-Raw-pOOBAH-passPercentage.png')
ggplot2::ggsave(filename = dir_png, plot = dir_gg)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Plot All Auto Sample Sheet:: Original
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# New Values to test::
#  core_tib$cg_1_pass_perc_basic_0
#  core_tib$cg_2_pass_perc_basic_0
#           cg_1_pass_perc_basic_0

#
# Build Plotting Data:: Original RAW(1) Method
#
plot_stack1_tib <- dplyr::bind_rows(
  dplyr::select(plot_core_tib, group_class,build_source,Sample, R2_1,cg_1_pOOBAH) %>% 
    purrr::set_names(c("Class","Experiment","Sample","Metric","Pval_Perc")) %>%
    dplyr::mutate(Metric_Type="R2",Design_Type="I") %>% 
    dplyr::select(Class,Experiment,Sample,Metric_Type,Design_Type,Metric,Pval_Perc),
  
  dplyr::select(plot_core_tib, group_class,build_source,Sample, R2_2,cg_2_pOOBAH) %>% 
    purrr::set_names(c("Class","Experiment","Sample","Metric","Pval_Perc")) %>%
    dplyr::mutate(Metric_Type="R2",Design_Type="II") %>% 
    dplyr::select(Class,Experiment,Sample,Metric_Type,Design_Type,Metric,Pval_Perc),
  
  dplyr::select(plot_core_tib, group_class,build_source,Sample, dB_1,cg_1_pOOBAH) %>% 
    purrr::set_names(c("Class","Experiment","Sample","Metric","Pval_Perc")) %>%
    dplyr::mutate(Metric_Type="dB",Design_Type="I", Metric=Metric/100) %>% 
    dplyr::select(Class,Experiment,Sample,Metric_Type,Design_Type,Metric,Pval_Perc),
  
  dplyr::select(plot_core_tib, group_class,build_source,Sample, dB_2,cg_2_pOOBAH) %>% 
    purrr::set_names(c("Class","Experiment","Sample","Metric","Pval_Perc")) %>%
    dplyr::mutate(Metric_Type="dB",Design_Type="II", Metric=Metric/100) %>% 
    dplyr::select(Class,Experiment,Sample,Metric_Type,Design_Type,Metric,Pval_Perc)
) %>%
  dplyr::mutate(Pval_Frac=Pval_Perc/100)


plot_stack1_tib %>% 
  dplyr::group_by(Class,Experiment,Metric_Type,Design_Type) %>% 
  dplyr::summarise(Pval_Min=min(Pval_Frac, na.rm=TRUE),
                   Pval_Perc=cntPer_gte(Pval_Frac,0.85), 
                   Metric_Perc=cntPer_gte(Metric,0.98), .groups="drop") %>% print(n=100)

pass_log_gg <- plot_stack1_tib %>% 
  # dplyr::filter(Metric>=0.6) %>%
  dplyr::filter(Class=="P") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x=log(Pval_Frac+1),y=log(Metric+1), color=Sample) ) +
  ggplot2::facet_grid(rows=vars(Class,Experiment),
                      cols=vars(Design_Type,Metric_Type) ) +
  geom_hline(yintercept=log(1.98) ) +
  geom_vline(xintercept=log(1.85) )
pass_log_png <- file.path(opt$outDir, 'pass.canonicalBeta.vs.passPercentage.log.png')
ggplot2::ggsave(filename = pass_log_png, plot = pass_log_gg)

fail_log_gg <- plot_stack1_tib %>% 
  # dplyr::filter(Metric>=0.6) %>%
  dplyr::filter(Class=="F") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x=log(Pval_Frac+1),y=log(Metric+1), color=Sample) ) +
  ggplot2::facet_grid(rows=vars(Class,Experiment),
                      cols=vars(Design_Type,Metric_Type) ) +
  geom_hline(yintercept=log(1.98) ) +
  geom_vline(xintercept=log(1.85) )
fail_log_png <- file.path(opt$outDir, 'fail.canonicalBeta.vs.passPercentage.log.png')
ggplot2::ggsave(filename = fail_log_png, plot = fail_log_gg)


both_log_gg <- plot_stack1_tib %>% 
  # dplyr::filter(Metric>=0.6) %>%
  # dplyr::filter(Class=="P") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x=log(Pval_Frac+1),y=log(Metric+1), color=Sample) ) +
  ggplot2::facet_grid(rows=vars(Class,Experiment),
                      cols=vars(Design_Type,Metric_Type) ) +
  geom_hline(yintercept=log(1.98) ) +
  geom_vline(xintercept=log(1.85) )
both_log_png <- file.path(opt$outDir, 'both.canonicalBeta.vs.passPercentage.log.png')
ggplot2::ggsave(filename = both_log_png, plot = both_log_gg)









#
# Linear Plotting::
#
plot_stack1_tib %>% 
  # dplyr::filter(Metric>=0.6) %>%
  dplyr::filter(Class=="F") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x=Pval_Frac,y=Metric, color=Sample) ) +
  ggplot2::facet_grid(rows=vars(Class,Experiment),
                      cols=vars(Design_Type,Metric_Type) ) +
  geom_hline(yintercept=0.98 ) +
  geom_vline(xintercept=0.85 )

plot_stack1_tib %>% 
  # dplyr::filter(Metric>=0.6) %>%
  dplyr::filter(Class=="P") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x=Pval_Frac,y=Metric, color=Sample) ) +
  ggplot2::facet_grid(rows=vars(Class,Experiment),
                      cols=vars(Design_Type,Metric_Type) ) +
  geom_hline(yintercept=0.98 ) +
  geom_vline(xintercept=0.85 )

#
# Log Versions::
#
plot_stack1_tib %>% 
  # dplyr::filter(Metric>=0.6) %>%
  dplyr::filter(Class=="F") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x=log(Pval_Frac+1),y=log(Metric+1), color=Sample) ) +
  ggplot2::facet_grid(rows=vars(Class,Experiment),
                      cols=vars(Design_Type,Metric_Type) ) +
  geom_hline(yintercept=log(1.98) ) +
  geom_vline(xintercept=log(1.85) )

plot_stack1_tib %>% 
  # dplyr::filter(Metric>=0.6) %>%
  dplyr::filter(Class=="P") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x=log(Pval_Frac+1),y=log(Metric+1), color=Sample) ) +
  ggplot2::facet_grid(rows=vars(Class,Experiment),
                      cols=vars(Design_Type,Metric_Type) ) +
  geom_hline(yintercept=log(1.98) ) +
  geom_vline(xintercept=log(1.85) )




# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Open Method::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Build Plotting Data:: Open01) Method
#
plot_stack0_tib <- dplyr::bind_rows(
  dplyr::select(plot_core_tib, group_class,build_source,Sample, R2_1,cg_1_Open_pOOBAH) %>% 
    purrr::set_names(c("Class","Experiment","Sample","Metric","Pval_Perc")) %>%
    dplyr::mutate(Metric_Type="R2",Design_Type="I") %>% 
    dplyr::select(Class,Experiment,Sample,Metric_Type,Design_Type,Metric,Pval_Perc),
  
  dplyr::select(plot_core_tib, group_class,build_source,Sample, R2_2,cg_2_Open_pOOBAH) %>% 
    purrr::set_names(c("Class","Experiment","Sample","Metric","Pval_Perc")) %>%
    dplyr::mutate(Metric_Type="R2",Design_Type="II") %>% 
    dplyr::select(Class,Experiment,Sample,Metric_Type,Design_Type,Metric,Pval_Perc),
  
  dplyr::select(plot_core_tib, group_class,build_source,Sample, dB_1,cg_1_Open_pOOBAH) %>% 
    purrr::set_names(c("Class","Experiment","Sample","Metric","Pval_Perc")) %>%
    dplyr::mutate(Metric_Type="dB",Design_Type="I", Metric=Metric/100) %>% 
    dplyr::select(Class,Experiment,Sample,Metric_Type,Design_Type,Metric,Pval_Perc),
  
  dplyr::select(plot_core_tib, group_class,build_source,Sample, dB_2,cg_2_Open_pOOBAH) %>% 
    purrr::set_names(c("Class","Experiment","Sample","Metric","Pval_Perc")) %>%
    dplyr::mutate(Metric_Type="dB",Design_Type="II", Metric=Metric/100) %>% 
    dplyr::select(Class,Experiment,Sample,Metric_Type,Design_Type,Metric,Pval_Perc)
) %>%
  dplyr::mutate(Pval_Frac=Pval_Perc/100)


plot_stack0_tib %>% 
  dplyr::group_by(Class,Experiment,Metric_Type,Design_Type) %>% 
  dplyr::summarise(Pval_Min=min(Pval_Frac, na.rm=TRUE),
                   Pval_Perc=cntPer_gte(Pval_Frac,0.85), 
                   Metric_Perc=cntPer_gte(Metric,0.98), .groups="drop") %>% print(n=100)

plot_stack0_tib %>% 
  dplyr::filter(Metric>=0.6) %>%
  # dplyr::filter(Class=="P") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x=log(Pval_Frac+1),y=log(Metric+1), color=Sample) ) +
  ggplot2::facet_grid(rows=vars(Class,Experiment),
                      cols=vars(Design_Type,Metric_Type) ) +
  geom_hline(yintercept=log(1.98) ) +
  geom_vline(xintercept=log(1.85) )










# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Open Sesame Method::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

plot_core1_dB_tib <- core_tib %>%
  dplyr::select(group_class,build_source,Sample,
                cg_basic_dB_2,cg_pass_perc_basic_0) %>%
  dplyr::rename(Metric=cg_basic_dB_2) %>%
  dplyr::mutate(Metric_Type="dB",
                Metric=Metric/100)

plot_core1_r2_tib <- core_tib %>%
  dplyr::select(group_class,build_source,Sample,
                cg_basic_r2_2,cg_pass_perc_basic_0) %>%
  dplyr::rename(Metric=cg_basic_r2_2) %>%
  dplyr::mutate(Metric_Type="r2")

plot_stack1_tib <- dplyr::bind_rows(
  plot_core1_dB_tib,
  plot_core1_r2_tib
) %>% 
  purrr::set_names("Class","Experiment","Sample","Metric","Pval_Frac","Metric_Type") %>%
  dplyr::filter(!is.na(Metric)) %>%
  dplyr::filter(!is.na(Pval_Frac)) %>%
  dplyr::mutate(Pval_Frac=Pval_Frac/100)

#
# Scatter Plot::
#
plot_stack1_tib %>% 
  # dplyr::filter(Metric>=0.9) %>%
  # dplyr::filter(Class=="F") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x=Pval_Frac,y=Metric, color=Sample) ) +
  ggplot2::facet_grid(rows=vars(Class,Experiment),
                      cols=vars(Metric_Type) ) +
  #                    cols=vars(Design_Type,Metric_Type) ) +
  geom_hline(yintercept=0.98 ) +
  geom_vline(xintercept=0.85 )


plot_stack1_tib %>% 
  # dplyr::filter(Metric>=0.6) %>%
  # dplyr::filter(Class=="F") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x=log(Pval_Frac+1),y=log(Metric+1), color=Sample) ) +
  ggplot2::facet_grid(rows=vars(Class,Experiment),
                      cols=vars(Metric_Type) ) +
  #                    cols=vars(Design_Type,Metric_Type) ) +
  geom_hline(yintercept=1.98 ) +
  geom_vline(xintercept=1.85 )

#
# Old plots look better
#

# End of file
