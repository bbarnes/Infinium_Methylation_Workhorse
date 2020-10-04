
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
par$prgmTag <- 'merge_builds_mm10_controls'
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

# Sample Sheet Parameters::
opt$addSampleName    <- FALSE
opt$addPathsCall     <- TRUE
opt$addPathsSigs     <- FALSE

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
  
  par$locRunStr <- 'mm10'
  
  if (par$locRunStr=='mm10') {
    par$runNameA  <- 'ILMN_mm10_betaTest_17082020'
    par$runNameB  <- 'VanAndel_mm10_betaTest_31082020'
    par$runNameC  <- 'MURMETVEP_mm10_betaTest_06082020'
    
    # par$runNameA  <- 'ILMN_mm10_betaTest_17082020.v1'
    # par$runNameB  <- 'VanAndel_mm10_betaTest_31082020'
    # par$runNameC  <- 'MURMETVEP_mm10_betaTest_06082020.v1'
    
    opt$runName   <- 'mm10_controls'
    opt$platform  <- 'LEGX'
    opt$version   <- 'B0'
    
    opt$buildDir  <- paste(
      file.path(par$topDir, 'scratch/swifthoof_main', par$runNameA),
      file.path(par$topDir, 'scratch/swifthoof_main', par$runNameB),
      file.path(par$topDir, 'scratch/swifthoof_main', par$runNameC),
      # file.path(par$topDir, 'scratch/swifthoof_main', par$runNameA, '12x1/BP2'),
      # file.path(par$topDir, 'scratch/swifthoof_main', par$runNameB, '12x1/BP2'),
      # file.path(par$topDir, 'scratch/swifthoof_main', par$runNameC, '12x1/BP2'),
      sep=',')
    
    # opt$sampleCsv <- file.path(par$samDir, 'sampleSheets/BETA-DELTA-EPIC-Core/BETA-DELTA-8x1-EPIC-Core.Sample_Names.SampleSheet.csv.gz')
    
  }
  
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
    make_option(c("--addPathsSigs"), action="store_true", default=opt$addPathsSigs, 
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
#                            Validate Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (is.null(par$runMode) || is.null(par$prgmDir) || is.null(par$prgmTag) || 
    is.null(par$scrDir) || is.null(par$datDir)) {
  
  par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
  par_tib %>% base::print(n=base::nrow(par_tib) )

  if (is.null(par$runMode)) cat(glue::glue("[Usage]: runMode is NULL!!!{RET}"))
  if (is.null(par$prgmDir)) cat(glue::glue("[Usage]: prgmDir is NULL!!!{RET}"))
  if (is.null(par$prgmTag)) cat(glue::glue("[Usage]: prgmTag is NULL!!!{RET}"))
  if (is.null(par$scrDir))  cat(glue::glue("[Usage]: scrDir is NULL!!!{RET}"))
  if (is.null(par$datDir))  cat(glue::glue("[Usage]: darDir is NULL!!!{RET}"))
  base::stop("Null Parameters!\n\n")
}

if (is.null(opt$outDir) || is.null(opt$buildDir) || 
     is.null(opt$runName) || 
     # is.null(opt$sampleCsv) || 
     is.null(opt$classVar) ||
     is.null(opt$addSampleName) || is.null(opt$addPathsCall) || is.null(opt$addPathsSigs) ||
     is.null(opt$flagDetectPval) || is.null(opt$flagSampleDetect) || is.null(opt$flagRefMatch) ||
     # is.null(opt$pvalDetectMinKey) || is.nulll(opt$pvalDetectMinVal) ||
     is.null(opt$platform) || is.null(opt$version) ||
     is.null(opt$percisionBeta) || is.null(opt$percisionPval) ||
     is.null(opt$execute) || is.null(opt$single) || is.null(opt$parallel) || is.null(opt$cluster) ||
     
     is.null(opt$clean) || is.null(opt$Rscript) || is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )

  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  if (is.null(opt$outDir))    cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$buildDir))  cat(glue::glue("[Usage]: buildDirs is NULL!!!{RET}"))
  if (is.null(opt$runName))   cat(glue::glue("[Usage]: runName is NULL!!!{RET}"))
  if (is.null(opt$sampleCsv)) cat(glue::glue("[Usage]: sampleCsv is NULL (Not Required)!!!{RET}"))
  if (is.null(opt$classVar))  cat(glue::glue("[Usage]: class_var is NULL!!!{RET}"))

  if (is.null(opt$addSampleName)) cat(glue::glue("[Usage]: addSampleName is NULL!!!{RET}"))
  if (is.null(opt$addPathsCall))  cat(glue::glue("[Usage]: addPathsCall is NULL!!!{RET}"))
  if (is.null(opt$addPathsSigs))  cat(glue::glue("[Usage]: addPathsSigs is NULL!!!{RET}"))
  
  if (is.null(opt$flagDetectPval))   cat(glue::glue("[Usage]: flagDetectPval is NULL!!!{RET}"))
  if (is.null(opt$flagSampleDetect)) cat(glue::glue("[Usage]: flagSampleDetect is NULL!!!{RET}"))
  if (is.null(opt$flagRefMatch))     cat(glue::glue("[Usage]: flagRefMatch is NULL!!!{RET}"))
  
  # if (is.null(opt$pvalDetectMinKey)) cat(glue::glue("[Usage]: pvalDetectMinKey is NULL!!!{RET}"))
  # if (is.null(opt$pvalDetectMinVal)) cat(glue::glue("[Usage]: pvalDetectMinVal is NULL!!!{RET}"))
  
  if (is.null(opt$platform)) cat(glue::glue("[Usage]: platform is NULL!!!{RET}"))
  if (is.null(opt$version))  cat(glue::glue("[Usage]: version is NULL!!!{RET}"))
  
  if (is.null(opt$percisionBeta)) cat(glue::glue("[Usage]: percisionBeta is NULL!!!{RET}"))
  if (is.null(opt$percisionPval)) cat(glue::glue("[Usage]: percisionPval is NULL!!!{RET}"))
  
  if (is.null(opt$execute))  cat(glue::glue("[Usage]: execute is NULL!!!{RET}"))
  if (is.null(opt$single))   cat(glue::glue("[Usage]: single is NULL!!!{RET}"))
  if (is.null(opt$parallel)) cat(glue::glue("[Usage]: parallel is NULL!!!{RET}"))
  if (is.null(opt$cluster))  cat(glue::glue("[Usage]: cluster is NULL!!!{RET}"))
  
  if (is.null(opt$clean))   cat(glue::glue("[Usage]: clean is NULL!!!{RET}"))
  if (is.null(opt$Rscript)) cat(glue::glue("[Usage]: Rscript is NULL!!!{RET}"))
  if (is.null(opt$verbose)) cat(glue::glue("[Usage]: verbosity is NULL!!!{RET}"))
  base::stop(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
}
par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
if (opt$verbose>=1) par_tib %>% base::print(n=base::nrow(par_tib) )
if (opt$verbose>=1) opt_tib %>% base::print(n=base::nrow(opt_tib) )
cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )

par$prgm_src_dir <- file.path(par$scrDir,par$prgmDir, 'functions')
if (!dir.exists(par$prgm_src_dir)) stop(glue::glue("[{par$prgmTag}]: Program Source={par$prgm_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$prgm_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Program Source={par$prgm_src_dir}!{RET}{RET}") )

# Load All other function methods::
par$man_src_dir <- file.path(par$scrDir, 'manifests/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$man_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$man_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$man_src_dir}!{RET}{RET}") )

par$swt_src_dir <- file.path(par$scrDir, 'swifthoof/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$swt_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$swt_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$swt_src_dir}!{RET}{RET}") )

par$prb_src_dir <- file.path(par$scrDir, 'probe_design/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$prb_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$prb_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$prb_src_dir}!{RET}{RET}") )

par$anl_src_dir <- file.path(par$scrDir, 'analysis/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$anl_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$anl_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$anl_src_dir}!{RET}{RET}") )

cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

blds_dir_vec  <- opt$buildDir %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

if (is.null(opt$classVar)) opt$classVar <- 'Source_Sample_Name'
class_var <- rlang::sym(opt$classVar)
class_idx <- rlang::sym("Class_Idx")

opt <- setLaunchExe(opts=opt, pars=par, verbose=opt$verbose, vt=5,tc=0)

opt$outDir <- file.path(opt$outDir, opt$platform, opt$version, opt$classVar, opt$runName)
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
cat(glue::glue("[{par$prgmTag}]: Built; OutDir={opt$outDir}!{RET}") )

cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing!{RET}{RET}") )
print(blds_dir_vec)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Load Genome Studio Beta Sample Sheets::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

gs_rep_csv <- file.path(par$datDir, 'sampleSheets/mm10/Replicates.ILS-VAI.analytical_SampleSheet.csv.gz')
gs_tit_csv <- file.path(par$datDir, 'sampleSheets/mm10/Titration.ILS-VAI.analytical_SampleSheet.csv.gz')

gs_rep_tib <- suppressMessages(suppressWarnings( readr::read_csv(gs_rep_csv) ))
gs_tit_tib <- suppressMessages(suppressWarnings( readr::read_csv(gs_tit_csv) ))

gs_ss_dir  <- file.path(par$datDir, 'sampleSheets/mm10/GenomeStudio')
gs_ss_csvs <- list.files(gs_ss_dir, pattern='.csv', full.names=TRUE)
gs_ss_tib  <- lapply(gs_ss_csvs, readr::read_csv, skip=7) %>% dplyr::bind_rows() %>% 
  tidyr::separate(Sample_Name, into=c('Laird_ID','Sentrix_Chip','Sentrix_Pos','Species',
                                      'Sample_Str1', 'Sample_Str2','Sample_Str3'), sep='_', remove=FALSE)

default_ng <- 251
ranks_vec <- c('R2', 'R3', 'R4')
input_vec <- c('50ng', '100ng', '250ng', '500ng', '1000ng')
gs_ss_col <- c('Laird_ID','Group','Sentrix_Chip','Sentrix_Pos','Species','Sample_Str1','Sample_Str2',
               'Sample_Well','Sample_Plate','Sample_Group','Pool_ID','Sentrix_ID','Sentrix_Position')

basic_ss_tib <- gs_ss_tib %>% 
  dplyr::filter(! Sentrix_Chip %in% ranks_vec) %>% 
  dplyr::filter(! Sentrix_Chip %in% input_vec) %>% 
  dplyr::mutate(Group=NA) %>% dplyr::select(dplyr::all_of(gs_ss_col) ) %>%
  dplyr::mutate(Input_ng=as.integer(default_ng))

ranks_ss_tib <- gs_ss_tib %>% 
  dplyr::filter(  Sentrix_Chip %in% ranks_vec) %>% 
  dplyr::rename(Group=Sentrix_Chip, Sentrix_Chip=Sentrix_Pos, Sentrix_Pos=Species, 
                Species=Sample_Str1, Sample_Str1=Sample_Str2, Sample_Str2=Sample_Str3) %>%
  dplyr::select(dplyr::all_of(gs_ss_col) ) %>%
  dplyr::mutate(Input_ng=as.integer(default_ng))

input_ss_tib <- gs_ss_tib %>% 
  dplyr::filter(  Sentrix_Chip %in% input_vec) %>% 
  dplyr::rename(Group=Sentrix_Chip, Sentrix_Chip=Sentrix_Pos, Sentrix_Pos=Species, 
                Species=Sample_Str1, Sample_Str1=Sample_Str2, Sample_Str2=Sample_Str3) %>%
  dplyr::select(dplyr::all_of(gs_ss_col) ) %>%
  dplyr::mutate(Input_ng=stringr::str_remove(Group,'ng') %>% as.integer())

input_s1_tib <- input_ss_tib %>% dplyr::filter(! Sentrix_Chip %in% ranks_vec)
input_s2_tib <- input_ss_tib %>% dplyr::filter(  Sentrix_Chip %in% ranks_vec) %>% 
  dplyr::select(-Sentrix_Chip) %>%
  dplyr::rename(Sentrix_Chip=Sentrix_Pos, Sentrix_Pos=Species, 
                Species=Sample_Str1, Sample_Str1=Sample_Str2)

#
# Sample Naming::
#
titration_org_vec <- c("0%MethylMouse","10%Human90%Mouse","10%MethylMouse","100%Human","100%MethylMouse","25%Human75%Mouse","25%MethylMouse","5%MethylMouse",
                       "50%Human50%Mouse","50%MethylMouse","75%Human25%Mouse","75%MethylMouse","90%Human10%Mouse","Mouse" )

titration_new_vec <- c("T00DZ","Z10H90M","T10DZ","Z99H00M","T99DZ","Z25H75M","T25DZ","T05DZ",
                       "Z50H50M","T50DZ","Z75H25M","T75DZ","Z90H10M","Z00H99M" )

titration_hm_vec <-  c(0,10, 0,99, 0,25, 0,0,50, 0,75, 0,90, 0 )
titration_mm_vec <-  c(0,90,10, 0,99,75,25,5,50,50,25,75,10,99 )

titration_map_tib <- tibble::tibble(Org_Sample_Name=titration_org_vec, Sample_Name=titration_new_vec,
                                    HS_Perc=titration_hm_vec, MM_Perc=titration_mm_vec)

#
# Joining of data
#
hum_ss_tib  <- dplyr::bind_rows(basic_ss_tib,ranks_ss_tib,input_s1_tib,input_s2_tib) %>% 
  dplyr::rename(Service_Location=Sample_Str1) %>%
  dplyr::select(-Sample_Str2,-Group) %>% 
  dplyr::mutate(Sentrix_Chip=as.double(Sentrix_Chip)) %>%
  dplyr::filter(Sentrix_Chip==Sentrix_ID) %>% 
  dplyr::left_join(titration_map_tib, by=c("Species"="Org_Sample_Name")) %>% 
  dplyr::select(-Species,-Sentrix_Chip,-Sentrix_Pos) %>%
  tidyr::unite(Sentrix_Name, Sentrix_ID,Sentrix_Position, sep='_') %>% 
  dplyr::select(Sentrix_Name,Service_Location,Laird_ID,Sample_Name,Input_ng,HS_Perc,MM_Perc,everything())

if (FALSE) {
  # 61    20382   4 Replicates     Z00H99M       0      99     9
  # 62    20383   4 Replicates     Z00H99M       0      99     3
  # 63    20384   4 Replicates     Z00H99M       0      99    12
  # 64    20385   4 Replicates     Z00H99M       0      99    12
  rep_vec <- c(20382,20383,20384,20385)
  rep_vec <- c(20382,20384,20385,21026)
  dplyr::bind_rows(basic_ss_tib,ranks_ss_tib,input_s1_tib,input_s2_tib) %>% 
    dplyr::rename(Service_Location=Sample_Str1) %>%
    dplyr::select(-Sample_Str2,-Group) %>% 
    dplyr::mutate(Sentrix_Chip=as.double(Sentrix_Chip)) %>%
    dplyr::filter(Sentrix_Chip==Sentrix_ID) %>% 
    dplyr::left_join(titration_map_tib, by=c("Species"="Org_Sample_Name")) %>% 
    dplyr::select(-Species,-Sentrix_ID,-Sentrix_Pos) %>%
    dplyr::filter(Laird_ID %in% rep_vec) %>% dplyr::group_by(Service_Location,Laird_ID,Sample_Group,Input_ng) %>% dplyr::summarise(Count=n())

  dplyr::bind_rows(basic_ss_tib,ranks_ss_tib,input_s1_tib,input_s2_tib) %>% 
    dplyr::rename(Service_Location=Sample_Str1) %>%
    dplyr::select(-Sample_Str2,-Group) %>% 
    dplyr::mutate(Sentrix_Chip=as.double(Sentrix_Chip)) %>%
    dplyr::filter(Sentrix_Chip==Sentrix_ID) %>% 
    dplyr::select(-Sentrix_Chip) %>%
    dplyr::mutate(Sentrix_Name=paste(Sentrix_ID,Sentrix_Pos, sep='_')) %>%
    dplyr::left_join(titration_map_tib, by=c("Species"="Org_Sample_Name")) %>% 
    dplyr::select(-Species,-Sentrix_ID,-Sentrix_Pos) %>%
    # dplyr::group_by(Service_Location,Sample_Group,Sample_Name,HS_Perc,MM_Perc) %>% 
    dplyr::inner_join(gs_rep_tib, by="Sentrix_Name")
  
  dplyr::bind_rows(basic_ss_tib,ranks_ss_tib,input_s1_tib,input_s2_tib) %>% 
    dplyr::rename(Service_Location=Sample_Str1) %>%
    dplyr::select(-Sample_Str2,-Group) %>% 
    dplyr::mutate(Sentrix_Chip=as.double(Sentrix_Chip)) %>%
    dplyr::filter(Sentrix_Chip==Sentrix_ID) %>% 
    dplyr::select(-Sentrix_Chip) %>%
    dplyr::mutate(Sentrix_Name=paste(Sentrix_ID,Sentrix_Pos, sep='_')) %>%
    dplyr::left_join(titration_map_tib, by=c("Species"="Org_Sample_Name")) %>% 
    dplyr::select(-Species,-Sentrix_ID,-Sentrix_Pos) %>%
    # dplyr::group_by(Service_Location,Sample_Group,Sample_Name,HS_Perc,MM_Perc) %>% 
    dplyr::inner_join(gs_rep_tib, by="Sentrix_Name", suffix=c("_v1","_v2")) %>% dplyr::select(Laird_ID_v1,Sample_Name_v1,Sample_Name_v2,Sample_Group) %>% dplyr::distinct()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Load Auto Detect Sample Sheets::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

auto_ss_tib <- NULL

opt$addPathsSigs <- FALSE
opt$addPathsSigs <- TRUE

opt$verbose <- 30

for (curDir in blds_dir_vec) {
  if (opt$findSampleSheet) {
    
    cur_hm_csv <- file.path(curDir, par$humanSampleSheetName)
    cat(glue::glue("[{par$prgmTag}]:{TAB}Loading Human Classification; cur_hm_csv={cur_hm_csv}!{RET}") )
    
    cur_hm_tib <- suppressMessages(suppressWarnings( readr::read_csv(cur_hm_csv) ))
    cur_hm_len <- cur_hm_tib %>% base::nrow()
    cat(glue::glue("[{par$prgmTag}]:{TAB}Done. Loading Human Classification; cur_hm_len={cur_hm_len}!{RET}{RET}") )
    # print(cur_hm_tib)
    
    hum_ss_tib <- dplyr::bind_rows(hum_ss_tib,cur_hm_tib)
  }
  
  cur_ss_tib <- loadAutoSampleSheets(dir=curDir, platform=opt$platform, manifest=opt$version,
                                     addSampleName=opt$addSampleName,  addPathsCall=opt$addPathsCall, addPathsSigs=opt$addPathsSigs,
                                     flagDetectPval=opt$flagDetectPval, flagSampleDetect=opt$flagSampleDetect, flagRefMatch=opt$flagRefMatch,
                                     pvalDetectMinKey=opt$pvalDetectMinKey, pvalDetectMinVal=opt$pvalDetectMinVal,
                                     
                                     verbose=opt$verbose,tc=2,tt=NULL) %>% 
    dplyr::distinct(Sentrix_Name, .keep_all=TRUE)
  
  auto_ss_tib <- auto_ss_tib %>% dplyr::bind_rows(cur_ss_tib) %>% dplyr::distinct(Sentrix_Name, .keep_all=TRUE)
  cur_ss_len  <- cur_ss_tib %>% base::nrow()
  auto_ss_len <- auto_ss_tib %>% base::nrow()
  
  cat(glue::glue("[{par$prgmTag}]:{TAB} Raw Auto Sample Sheet Current={cur_ss_len}, Total={auto_ss_len}.{RET}"))
  # print(cur_ss_tib)
}
auto_ss_len <- auto_ss_tib %>% base::nrow()
cat(glue::glue("[{par$prgmTag}]: Done. Raw Auto Sample Sheet; Total={auto_ss_len}.{RET}{RET}"))
# print(auto_ss_tib)

# sel_hum_ss_tib <- dplyr::filter(hum_ss_tib, stringr::str_detect(Sample_Group, 'Mix')) %>% 
#   dplyr::inner_join(auto_ss_tib, by="Sentrix_Name")

sel_hum_ss_tib <- hum_ss_tib %>% dplyr::filter( stringr::str_detect(Sample_Group, 'Mix') ) %>% 
  dplyr::filter( (HS_Perc==0 & MM_Perc==99) | (HS_Perc==99 & MM_Perc==0) | (HS_Perc==50 & MM_Perc==50) ) %>%
  dplyr::inner_join(auto_ss_tib, by="Sentrix_Name")

man_360_csv <- file.path(opt$outDir, 'manifest-360.csv.gz')
readr::write_csv(hum_ss_tib, man_360_csv)

max_load <- sel_hum_ss_tib %>% base::nrow()
sigs_list <- sel_hum_ss_tib$Sigs_Path %>% head(n=max_load)
names(sigs_list) <- sel_hum_ss_tib$Sentrix_Name %>% head(n=max_load)

call_list <- sel_hum_ss_tib$Calls_Path %>% head(n=max_load)
names(call_list) <- sel_hum_ss_tib$Sentrix_Name %>% head(n=max_load)

sigs_stack_tib <- lapply(sigs_list, readr::read_csv) %>% dplyr::bind_rows(.id="Sentrix_Name") %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2) )
call_stack_tib <- lapply(call_list, readr::read_csv) %>% dplyr::bind_rows(.id="Sentrix_Name")
detp_stack_tib <- call_stack_tib %>% dplyr::select(Sentrix_Name,Probe_ID,i_negs,i_poob)

opt$pvalMin <- 0.05

#
# Look over control comparisons...
#
hm_ctl_vec <- c('ctl_BS_Conversion_I_','ctl_BS_Conversion_II_','ctl_NP_','ctl_Negative')
mm_ctl_vec <- c('BS','BS','NO','ne')
dp_ctl_vec <- c(0.05,0.05,0.05,1.0)
n1_ctl_vec <- c('II','IG','IG','IG')
n2_ctl_vec <- c('II','IR','IR','IR')

# n1_ctl_vec <- c('IG','II','II','II')
# n2_ctl_vec <- c('IR','II','II','II')

tar_ctl_tib <- tibble::tibble(hm_ctl_key=hm_ctl_vec,
                              mm_ctl_key=mm_ctl_vec,
                              dp_ctl_val=dp_ctl_vec,
                              n1_ctl_val=n1_ctl_vec,
                              n2_ctl_val=n2_ctl_vec)

tar_ctl_len <- tar_ctl_tib %>% base::nrow()

sel_dat_ss_tib  <- dplyr::select(sel_hum_ss_tib, Sentrix_Name, Service_Location, Laird_ID, Sample_Name, HS_Perc, MM_Perc)
prb_int_det_tib <- dplyr::left_join(sigs_stack_tib, detp_stack_tib, by=c("Sentrix_Name","Probe_ID") ) %>%
  dplyr::mutate(
    Probe_Design=dplyr::case_when(
      stringr::str_starts(Probe_ID, 'ctl_BS_Conversion_I_C') ~ 'IG',
      stringr::str_starts(Probe_ID, 'ctl_BS_Conversion_I_M') ~ 'IG',
      stringr::str_starts(Probe_ID, 'ctl_BS_Conversion_I_U') ~ 'IR',
      TRUE ~ Probe_Design )
  ) %>%
  dplyr::filter(!is.na(i_poob)) %>%
  tidyr::pivot_longer(cols=c("M","U"), names_to="Channel", values_to="Intensity") %>% dplyr::filter(!is.na(Intensity))

for (tIdx in c(1:tar_ctl_len)) {
  hm_ctl_key <- tar_ctl_tib[['hm_ctl_key']][tIdx]
  mm_ctl_key <- tar_ctl_tib[['mm_ctl_key']][tIdx]
  dp_ctl_val <- tar_ctl_tib[['dp_ctl_val']][tIdx]
  n1_ctl_val <- tar_ctl_tib[['n1_ctl_val']][tIdx]
  n2_ctl_val <- tar_ctl_tib[['n2_ctl_val']][tIdx]
  
  out_name <- paste(hm_ctl_key,mm_ctl_key,'min-pval',dp_ctl_val,n1_ctl_val,n2_ctl_val, sep='_')

  cat(glue::glue("{RET}[{par$prgmTag}]: Starting out_name={out_name}...{RET}"))
  
  cur_prb_tib <- prb_int_det_tib %>% 
    dplyr::filter(stringr::str_starts(Probe_ID, hm_ctl_key) | Probe_Type==mm_ctl_key) %>%
    dplyr::filter(Probe_Design!=n1_ctl_val) %>% 
    dplyr::filter(Probe_Design!=n2_ctl_val) %>% 
    dplyr::filter(Probe_Type=='ct' | i_poob<=dp_ctl_val)

  # QC Check::
  cur_prb_tib %>% dplyr::group_by(Probe_Design,Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()
  cur_ann_tib <- dplyr::inner_join(sel_dat_ss_tib,cur_prb_tib, by="Sentrix_Name")
    
  # QC Check::
  cur_ann_tib %>% dplyr::group_by(Service_Location,Sample_Name,Probe_Design,Probe_Type) %>% dplyr::summarise(Count=n()) %>% as.data.frame()
  
  cur_ann_pdf <- file.path(opt$outDir, paste(out_name, '.plot-Probe-Intensity.pdf') )
  cur_ann_gg  <- ggplot2::ggplot(data=cur_ann_tib, aes(x=Probe_ID, y=Intensity, color=i_poob)) + 
    ggplot2::geom_point(na.rm = TRUE) +
    # ggplot2::facet_grid(rows = vars(Service_Location,Probe_Design), cols = vars(Sample_Name,Probe_Type,Channel), scales="free_x" ) +
    ggplot2::facet_grid(rows = vars(Service_Location,Probe_Design), cols = vars(Sample_Name,Probe_Type), scales="free_x" ) +
    # ggplot2::facet_grid(rows = vars(Service_Location,Probe_Design), cols = vars(Sample_Name,Probe_Type), scales="free" ) +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
  
  ggplot2::ggsave(cur_ann_pdf,cur_ann_gg)
  # cur_ann_gg
  
  if (mm_ctl_key=='ne') {
    max_hum_int <- 1000
    max_hum_int <- 500
  } else {
    max_hum_int <- cur_ann_tib %>% dplyr::filter(Sample_Name=="Z99H00M") %>% dplyr::filter(Probe_Type=='ct') %>% 
      dplyr::summarise(max=max(Intensity, na.rm=TRUE)) %>% dplyr::pull(max) %>% as.integer()
  }
  
  cur_pval_pdf <- file.path(opt$outDir, paste(out_name, '.plot-DetP-vs-Intensity.pdf') )
  cur_pval_gg  <- ggplot2::ggplot(data=cur_ann_tib, aes(x=i_poob, y=Intensity)) + 
    ggplot2::geom_point(na.rm = TRUE) +
    # ggplot2::geom_density_2d(alpha=0.6, na.rm = TRUE) +
    # ggplot2::facet_grid(rows = vars(Service_Location,Probe_Design), cols = vars(Sample_Name,Probe_Type,Channel), scales="free_x" ) +
    ggplot2::facet_grid(rows = vars(Service_Location,Probe_Design), cols = vars(Sample_Name,Probe_Type), scales="free_x" ) +
    # ggplot2::facet_grid(rows = vars(Service_Location,Probe_Design), cols = vars(Sample_Name,Probe_Type), scales="free" ) +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +
    ylim(0,max_hum_int)

  ggplot2::ggsave(cur_pval_pdf, cur_pval_gg)
  # cur_pval_gg

  cat(glue::glue("{RET}[{par$prgmTag}]: Done. Plotted={cur_pval_pdf}.{RET}{RET}"))
}



















nonp_tib <- dplyr::select(sel_hum_ss_tib, Sentrix_Name, Service_Location, Laird_ID, Sample_Name, HS_Perc, MM_Perc) %>%
  dplyr::inner_join(
    dplyr::mutate(sigs_stack_tib, Probe_Type=stringr::str_sub(Probe_ID, 1,2)) %>% 
      dplyr::filter( stringr::str_starts(Probe_ID, 'ctl_NP_') | stringr::str_starts(Probe_ID, 'NO') ) %>%
      dplyr::left_join(detp_stack_tib, by=c("Sentrix_Name","Probe_ID") ), by="Sentrix_Name") %>% 
  dplyr::filter(!is.na(i_poob)) %>%
  dplyr::filter( Probe_Type=='ct' | i_poob<opt$pvalMin) %>% 
  tidyr::pivot_longer(cols=c("M","U"), names_to="Channel", values_to="Intensity") %>% dplyr::filter(!is.na(Intensity))

top8_tib <- nonp_tib %>% dplyr::group_by(Service_Location,Sample_Name,Probe_Type) %>% dplyr::arrange( Intensity) %>% dplyr::top_n(n=8)
min8_tib <- nonp_tib %>% dplyr::group_by(Service_Location,Sample_Name,Probe_Type) %>% dplyr::arrange(-Intensity) %>% dplyr::top_n(n=8)
all8_tib <- dplyr::bind_rows(top8_tib,min8_tib) %>% dplyr::ungroup() %>% dplyr::distinct()

all8_gg  <- ggplot2::ggplot(data=all8_tib, aes(x=Probe_ID, y=Intensity, color=i_poob)) + 
  ggplot2::geom_point(na.rm = TRUE) +
  ggplot2::facet_grid(rows = vars(Service_Location), cols = vars(Sample_Name,Probe_Type,Channel), scales="free_x" ) +
  theme(axis.text.x = element_text(angle = 45))
all8_gg




max_hum_int <- nonp_tib %>% dplyr::filter(Sample_Name=="Z99H00M") %>% dplyr::filter(Probe_Type=='ct') %>% 
  dplyr::summarise(max=max(Intensity, na.rm=TRUE)) %>% dplyr::pull(max) %>% as.integer()

ctno_tib <- nonp_tib %>% dplyr::filter(Probe_Type=='ct') %>% dplyr::arrange(-Intensity)
ctgg <- ggplot2::ggplot(data=ctno_tib, aes(x=Intensity)) + 
  ggplot2::geom_density(na.rm=TRUE) + 
  ggplot2::facet_grid(rows = vars(Service_Location), cols = vars(Sample_Name,Probe_Type,Channel), scales="free" )

d2gg <- ggplot2::ggplot(data=ctno_tib, aes(x=i_poob, y=Intensity)) + 
  ggplot2::geom_density_2d(na.rm=TRUE) + 
  ggplot2::facet_grid(rows = vars(Service_Location), cols = vars(Sample_Name,Probe_Type,Channel) ) # , scales="free" )

ctgg <- ggplot2::ggplot(data=ctno_tib, aes(x=Probe_ID, y=Intensity, color=i_poob)) + 
  ggplot2::geom_point(na.rm = TRUE) +
  ggplot2::facet_grid(rows = vars(Service_Location), cols = vars(Sample_Name,Probe_Type,Channel), scales="free_x" ) +
  theme(axis.text.x = element_text(angle = 90))


max_hum_int <- 15312
max_hum_int <- 7397
max_hum_int <- 2000

gg <- ggplot2::ggplot(data=nonp_tib, aes(x=Probe_ID, y=Intensity, color=i_poob)) + 
  ggplot2::geom_point(na.rm = TRUE) +
  ggplot2::facet_grid(rows = vars(Service_Location), cols = vars(Sample_Name,Probe_Type,Channel), scales="free_x" ) +
  ylim(0,max_hum_int)

gg



sigs_stack_tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2)) %>% 
  dplyr::filter(Probe_Type=='ct') %>% dplyr::distinct(Probe_ID, .keep_all=TRUE) %>% 
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ctl_Biotin')) %>%
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ctl_Negative')) %>%
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ctl_BS_Conversion')) %>%
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ctl_DNP')) %>%
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ctl_Extension')) %>%
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ctl_GT_Mismatch')) %>%
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ctl_Hyb')) %>%
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ctl_Norm')) %>%
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ctl_Specificity')) %>%
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ctl_Target_Removal')) %>%
  dplyr::filter(! stringr::str_starts(Probe_ID, 'ctl_Restore')) %>%
  dplyr::left_join(detp_stack_tib, by=c("Sentrix_Name","Probe_ID") ) %>%
  as.data.frame()


  

sigs_stack_tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2)) %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n())

sigs_stack_tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2)) %>% dplyr::filter(Probe_Type=='ct')



sigs_stack_tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2)) %>% dplyr::filter(Probe_Type=='NO') %>%
  dplyr::left_join(detp_stack_tib, by=c("Sentrix_Name","Probe_ID") ) %>%
  tidyr::pivot_wider(names_from = "Sentrix_Name", values_from = c("U","M","ind_negs","ind_poob"))


sigs_stack_tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2)) %>% dplyr::filter(Probe_Type=='NO') %>%
  dplyr::left_join(detp_stack_tib, by=c("Sentrix_Name","Probe_ID") ) %>%
  tidyr::pivot_wider(names_from = "Sentrix_Name", values_from = c("U","M","ind_negs","ind_poob")) %>% as.data.frame()



sigs_stack_tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2)) %>% dplyr::filter(Probe_Type=='NO') %>%
  dplyr::left_join(detp_stack_tib, by=c("Sentrix_Name","Probe_ID") ) %>%
  tidyr::pivot_wider(names_from = Sentrix_Name, values_from = c("U","M","ind_negs","ind_poob"), names_glue="{Sentrix_Name}_{.value}", names_sort=TRUE) %>% as.data.frame()



tidyr::pivot_wider(data = sigs_stack_tib, names_from="Probe_Design", values_from=c("U","M"))

sigs_stack_tib %>% tidyr::pivot_wider(
  names_from = c("Sentrix_Name","Probe_ID"), values_from="Probe_Design")

sigs_stack_tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2)) %>% dplyr::filter(Probe_Type=='NO') %>%
  tidyr::pivot_wider(names_from = "Sentrix_Name", values_from = c("U","M"))


sigs_stack_tib %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2)) %>% dplyr::filter(Probe_Type=='NO') %>% tidyr::spread(Sentrix_Name, M,U, -c("Probe_ID","Probe_Design","Probe_Type") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Load Humman Annotation Sample Sheets::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Human Classifcation Sample Sheets are assumed to formatted!!!
#
# TBD::
#  - Assume any unidentified sample is COVID-
#  - Add Nasal Swabs to master sample sheet
#
labs_ss_tib <- NULL
if (! is.null(hum_ss_tib)) {

  labs_ss_tib <- auto_ss_tib %>% dplyr::left_join(hum_ss_tib, by="Sentrix_Name") %>% dplyr::arrange(!!class_var)
  
} else if (!is.null(opt$sampleCsv) && file.exists(opt$sampleCsv)) {
  
  cat(glue::glue("[{par$prgmTag}]: Loading predfined human classification; sampleCsv='{opt$sampleCsv}'{RET}") )
  
  hum_ss_tib <- suppressMessages(suppressWarnings( readr::read_csv(opt$sampleCsv) ))
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
  # stop(glue::glue("[{par$prgmTag}]: Failed to find humman annotation sample sheet={opt$sampleCsv}!!!{RET}{RET}"))
  cat(glue::glue("[{par$prgmTag}]: Using Auto Classification; classVar='{opt$classVar}'{RET}") )
  
  hum_ss_tib <- auto_ss_tib %>% dplyr::select(Sentrix_Name, !!class_var) %>% dplyr::arrange(!!class_var)
  
  # Left Join now that we will force Sample_Class to nSARSCov2 (COVID-) below
  # labs_ss_tib <- auto_ss_tib %>% dplyr::inner_join(hum_ss_tib, by="Sentrix_Name") %>% dplyr::arrange(!!class_var)
  #
  auto_ss_tib <- auto_ss_tib %>% dplyr::select(- !!class_var)
  labs_ss_tib <- auto_ss_tib %>% dplyr::left_join(hum_ss_tib, by="Sentrix_Name") %>% dplyr::arrange(!!class_var)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Join Human and Auto Sample Sheets::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Not sure this is needed yet, but its to force every sample to have
#  a class. Should default to the last classs in opt$trainClass
#
# null_class <- 'nSARSCov2'
# null_class <- rlang::sym(null_class)
# labs_ss_tib <- labs_ss_tib %>% dplyr::mutate(!!class_var := dplyr::case_when(is.na(!!class_var) ~ 'nSARSCov2', TRUE ~ !!class_var) )

labs_ss_len <- labs_ss_tib %>% base::nrow()
cat(glue::glue("[{par$prgmTag}]: Done. Joining Human Classification and Auto Sample Sheets; Total={labs_ss_len}.{RET}{RET}"))
print(labs_ss_tib)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Import Datasets (Calls)::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

outName <- paste(opt$runName, opt$platform, opt$version, sep='_')

opt$beg_file <- file.path(opt$outDir, paste(outName, 'load-beg.txt', sep='_') )
opt$end_file <- file.path(opt$outDir, paste(outName, 'load-end.txt', sep='_') )
opt$labs_csv <- file.path(opt$outDir, paste(outName, 'AutoSampleSheet.csv.gz', sep='_') )
opt$call_csv <- file.path(opt$outDir, paste(outName, 'MergedDataFiles.tib.csv.gz', sep='_') )

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

if (!pass_time_check) {
  cat(glue::glue("[{par$prgmTag}]: Building from scratch...{RET}"))
  unlink( list.files(opt$outDir, full.names=TRUE) )
  
  cmd <- paste('touch',opt$beg_file, sep=' ')
  system(cmd)
  
  call_tib <- mergeCallsFromSS(ss=labs_ss_tib, max=0, outName=outName, outDir=opt$outDir, 
                               verbose=opt$verbose, vt=1, tc=1, tt=pTracker)
  
  # Write Merged Data Files Tibble::
  readr::write_csv(call_tib, opt$call_csv)
  
  cmd <- paste('touch',opt$end_file, sep=' ')
  system(cmd)
  
  cat(glue::glue("[{par$prgmTag}]: Done. Building from scratch.{RET}{RET}"))
} else {
  cat(glue::glue("[{par$prgmTag}]: Build Already Up to Date.{RET}{RET}"))
}

# Write Auto Sample Sheet::
readr::write_csv(labs_ss_tib, opt$labs_csv)

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
