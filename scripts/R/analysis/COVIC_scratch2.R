
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
par$prgmTag <- 'COVIC_scratch2'
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
  
  opt$single   <- FALSE
  opt$single   <- TRUE
  opt$cluster  <- FALSE
  opt$parallel <- FALSE
  
  #
  # Pre-defined local options runTypes::
  #
  par$local_runType <- 'CORE'
  par$local_runType <- 'EXCBR'
  par$local_runType <- 'NZT'
  par$local_runType <- 'GENK'
  par$local_runType <- 'COVID'
  par$local_runType <- 'GRCm38'
  par$local_runType <- 'qcMVP'
  par$local_runType <- 'COVIC'
  
  if (par$local_runType=='COVIC') {
    opt$runName  <- 'COVIC-Set7-06082020'
    opt$runName  <- 'COVIC-Set5-10062020'
    opt$runName  <- 'COVIC-Set1-15052020'
    
    opt$runName  <- 'COVIC_manifest_validation'
    opt$platform <- 'EPIC'
    opt$version  <- 'C0'
    
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
opt_reqs <- c('outDir','Rscript','verbose','clean','runName',
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

# old_man_csv <- '/Users/bretbarnes/Documents/data/COVIC/transfer/EPIC-C0.manifest.sesame-base.cpg-sorted_NOTINOTHER.csv.gz'
# new_man_csv <- '/Users/bretbarnes/Documents/data/COVIC/transfer/GRCh38_COVIC_C0.manifest.sesame-base.pos-sorted_NOTINOTHER.csv.gz'

covic_man_dir <- '/Users/bretbarnes/Documents/data/COVIC/transfer'

cepic_raw_csv <- file.path(covic_man_dir, 'EPIC-C0.manifest.sesame-base.cpg-sorted.csv.gz')
covic_raw_csv <- file.path(covic_man_dir, 'GRCh38_COVIC_C0.manifest.sesame-base.pos-sorted.csv.gz')
nztic_raw_csv <- file.path(covic_man_dir, 'GRCh38_NZT_N0.manifest.sesame-base.pos-sorted.csv.gz')

cepic_raw_tib <- suppressMessages(suppressWarnings( readr::read_csv(cepic_raw_csv) ))
covic_raw_tib <- suppressMessages(suppressWarnings( readr::read_csv(covic_raw_csv) ))
nztic_raw_tib <- suppressMessages(suppressWarnings( readr::read_csv(nztic_raw_csv) ))

cepic_man_tib <- cepic_raw_tib %>% dplyr::select(Probe_ID, M, U, Next_Base, Probe_Type, Probe_Source)
covic_man_tib <- covic_raw_tib %>% dplyr::select(Probe_ID, M, U, Next_Base, Probe_Type, Probe_Source)
nztic_man_tib <- nztic_raw_tib %>% dplyr::select(Probe_ID, M, U, Next_Base, Probe_Type, Probe_Source)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Check Overlap::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# 1. Build a list of unique tangos for each build

cepic_add_vec <- cepic_man_tib %>% dplyr::filter(Probe_Source=='EPIC-C0')

# 1. covic_man_tib vs. covid_man_tib

cepic_man_tib %>% dplyr::inner_join(covic_man_tib, by="U", suffix=c("_C", "_D"))
cepic_man_tib %>% dplyr::inner_join(covic_man_tib, by=c("M","U"), suffix=c("_C", "_D"))

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
