
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              COVIC Merging and Spliting of Swifthoof Data::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages(require("optparse",quietly=TRUE)))

suppressWarnings(suppressPackageStartupMessages(require("plyr")) )
suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )
suppressWarnings(suppressPackageStartupMessages(require("grid")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

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

# Program Parameters::
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'analysis'
par$prgmTag <- 'plot_byExperiment'

# Illumina based directories::
par$macDir  <- '/Users/bbarnes/Documents/Projects/methylation/tools'
par$lixDir  <- '/illumina/scratch/darkmatter/Projects/COVIC'

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
  
  if (dir.exists(par$macDir)) par$topDir <- '/Users/bbarnes/Documents/Projects/methylation/scratch'
  if (dir.exists(par$lixDir)) par$topDir <- '/illumina/scratch/darkmatter/data/scratch'
  if (!dir.exists(par$topDir)) dir.create(par$topDir, recursive=TRUE)
  
  # Default Parameters for local Mac::
  par$runMode    <- args.dat[1]
  par$srcDir     <- file.path(par$macDir, par$codeDir)
  par$scrDir     <- file.path(par$srcDir, 'scripts')
  par$exePath    <- file.path(par$scrDir, 'R', par$prgmDir, paste0(par$prgmTag,'.R'))

  par$prgmTag <- base::sub('\\.R$', '', base::basename(par$exePath))
  par$locPath <- base::dirname(par$exePath)
  par$scrDir  <- base::dirname(base::normalizePath(par$locPath) )
  par$srcDir  <- base::dirname(base::normalizePath(par$scrDir) )
  par$datDir  <- file.path(base::dirname(base::normalizePath(par$srcDir)), 'dat')

  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'

  opt$platform <- 'EPIC'
  opt$version  <- 'B4'
  
  
  grpA <- 'CNTL-Samples_VendA_10092020'
  grpB <- 'CNTL-Samples_VendB_10092020'
  
  opt$runName  <- 'VA_Controls'
  opt$buildDir <- paste(
    file.path('/Users/bbarnes/Documents/Projects/methylation/VA_MVP/scratch/swifthoof_main', grpA),
    file.path('/Users/bbarnes/Documents/Projects/methylation/VA_MVP/scratch/swifthoof_main', grpB),
    sep=',')
  
  opt$outDir <- file.path(par$topDir, par$prgmTag)
  
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
# par$gen_src_dir <- file.path(par$scrDir, 'functions')
# if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
# for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
# cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )
# 
# par$prgm_src_dir <- file.path(par$scrDir,par$prgmDir, 'functions')
# if (!dir.exists(par$prgm_src_dir)) stop(glue::glue("[{par$prgmTag}]: Program Source={par$prgm_src_dir} does not exist!{RET}"))
# for (sfile in list.files(path=par$prgm_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
# cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Program Source={par$prgm_src_dir}!{RET}{RET}") )


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
#                       Load Auto Detect Sample Sheets::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

hum_ss_tib  <- NULL
auto_ss_tib <- NULL

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
    dplyr::distinct(Sentrix_Name, .keep_all=TRUE) %>%
    dplyr::mutate(Experiment_Key=base::basename(curDir))
  
  auto_ss_tib <- auto_ss_tib %>% dplyr::bind_rows(cur_ss_tib) %>% dplyr::distinct(Sentrix_Name, .keep_all=TRUE)
  cur_ss_len  <- cur_ss_tib %>% base::nrow()
  auto_ss_len <- auto_ss_tib %>% base::nrow()
  
  cat(glue::glue("[{par$prgmTag}]:{TAB} Raw Auto Sample Sheet Current={cur_ss_len}, Total={auto_ss_len}.{RET}"))
  # print(cur_ss_tib)
}
auto_ss_len <- auto_ss_tib %>% base::nrow()
cat(glue::glue("[{par$prgmTag}]: Done. Raw Auto Sample Sheet; Total={auto_ss_len}.{RET}{RET}"))
# print(auto_ss_tib)

auto_ss_tib %>% dplyr::mutate(
  Experiment_Key=dplyr::case_when(
    Experiment_Key=='CNTL-Samples_VendA_10092020' ~ 'A',
    Experiment_Key=='CNTL-Samples_VendB_10092020' ~ 'B',
    TRUE ~ NA_character_
  )) %>%
  dplyr::group_by(Experiment_Key,AutoSample_dB_Key) %>% 
  dplyr::summarise(GCT_Avg0=mean(GCT_0_Score, na.rm=TRUE),GCT_Avg1=mean(GCT_1_Score, na.rm=TRUE), GCT_Avg2=mean(GCT_2_Score, na.rm=TRUE),
                   Beta2_Avg=mean(Beta_2_Mean,na.rm=TRUE), 
                   BS1_G_Avg=mean('BISULFITE CONVERSION I_G_avg', na.rm=TRUE), BS1_R_Avg=mean('BISULFITE CONVERSION I_R_avg', na.rm=TRUE),
                   BS2_G_Avg=mean('BISULFITE CONVERSION I_G_avg', na.rm=TRUE), BS2_R_Avg=mean('BISULFITE CONVERSION I_R_avg', na.rm=TRUE) )
                   


# Summary Stats::
auto_ss_tib %>% dplyr::mutate(
  Experiment_Key=dplyr::case_when(
    Experiment_Key=='CNTL-Samples_VendA_10092020' ~ 'A',
    Experiment_Key=='CNTL-Samples_VendB_10092020' ~ 'B',
    TRUE ~ NA_character_
  )) %>%
  dplyr::group_by(Experiment_Key,AutoSample_dB_Key) %>% 
  dplyr::summarise(GCT_Avg0=mean(GCT_0_Score, na.rm=TRUE),GCT_Avg1=mean(GCT_1_Score, na.rm=TRUE), GCT_Avg2=mean(GCT_2_Score, na.rm=TRUE),
                   Beta2_Avg=mean(Beta_2_Mean,na.rm=TRUE), 
                   BS1_G_Avg=mean(`BISULFITE CONVERSION I_G_avg`,  na.rm=TRUE), BS1_R_Avg=mean(`BISULFITE CONVERSION I_R_avg`, na.rm=TRUE),
                   BS2_G_Avg=mean(`BISULFITE CONVERSION II_G_avg`, na.rm=TRUE), BS2_R_Avg=mean(`BISULFITE CONVERSION II_R_avg`, na.rm=TRUE) ) %>% dplyr::arrange(AutoSample_dB_Key)














man_csv <- '/Users/bbarnes/Documents/Projects/methylation/tools/Infinium_Methylation_Workhorse/dat/manifest/base/EPIC-B4.manifest.sesame-base.cpg-sorted.csv.gz'
man_tib <- readr::read_csv(man_csv)


exp_ss_list <- auto_ss_tib %>% split(.$Experiment_Key)

sample_str <- 'HELA'

ss_tibA <- exp_ss_list[[1]] %>% dplyr::filter(AutoSample_dB_Key==sample_str)
ss_tibB <- exp_ss_list[[2]] %>% dplyr::filter(AutoSample_dB_Key==sample_str)

call_tibA <- loadCallFiles(files=ss_tibA$Calls_Path, selKey='Probe_ID', 
                           datKey='ind_beta', minKey='i_poob', minVal=0.1, retMin=FALSE, addSentrix=TRUE,
                           verbose = opt$verbose)
call_tibB <- loadCallFiles(files=ss_tibB$Calls_Path, selKey='Probe_ID', 
                           datKey='ind_beta', minKey='i_poob', minVal=0.1, retMin=FALSE, addSentrix=TRUE,
                           verbose = opt$verbose)

# Change ind_beta -> Beta and i_poob -> Pval???

tibA <- call_tibA %>% purrr::set_names(paste('A',names(.), sep='_') )
names(tibA)[1] <- 'Probe_ID'

tibB <- call_tibB %>% purrr::set_names(paste('B',names(.), sep='_') )
names(tibB)[1] <- 'Probe_ID'

join_prbs <- man_tib %>% dplyr::select(Probe_ID, Probe_Type, DESIGN) %>% dplyr::rename(Design_Type=DESIGN) %>%
  dplyr::right_join(
    join_tib <- dplyr::inner_join(tibA, tibB, by="Probe_ID"), by="Probe_ID")

nameA <- 'A'
nameB <- 'B'

gg <- plotPairs(join_prbs, sample=sample_str, nameA=nameA, nameB=nameB, outDir=opt$outDir,
                probeType='cg', field='ind_beta', field_str='ind_beta', detp='i_poob', # maxCnt=opt$plotMax, 
                minPval=opt$minPval,
                verbose=opt$verbose+30)
# spread=spread, outType=outType, dpi=dpi, format=format,



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Import Datasets (Calls)::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

outName <- paste(opt$runName, opt$platform, opt$version, sep='_')

opt$beg_file <- file.path(opt$outDir, paste(outName, 'load-beg.txt', sep='_') )
opt$end_file <- file.path(opt$outDir, paste(outName, 'load-end.txt', sep='_') )
opt$labs_csv <- file.path(opt$outDir, paste(outName, 'AutoSampleSheet.csv.gz', sep='_') )
opt$call_csv <- file.path(opt$outDir, paste(outName, 'MergedDataFiles.tib.csv.gz', sep='_') )

if (opt$clean) unlink( list.files(opt$outDir, full.names=TRUE) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        New Plotting Effort HERE::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sam_ss_types <- auto_ss_tib %>% split(.$AutoSample_dB_Key)

opt$auto_beta_field <- 'ind_beta'
opt$auto_pval_field <- 'i_poob'
opt$minPval <- 0.1
opt$joinFull <- TRUE
sample <- 'HELA'
# labs_ss_tib <- sam_ss_types[[sample]]
# labs_ss_tib <- sam_ss_types


opt$pvalPassMin <- 96

opt$poobMinPval <- 0.2
opt$negsMinPval <- 0.02

opt$plotMax <- 3
opt$plotSub <- 10000
opt$plotSub <- 5000

opt$plotSpread  <- 'top'
opt$plotSpread  <- 'bot'
opt$plotSpread  <- 'mid'

opt$plotType   <- 'points'
opt$plotType   <- 'density'
opt$plotType   <- 'auto'

opt$plotFormat <- 'pdf'
opt$plotFormat <- 'png'

opt$dpi <- 72
opt$dpi <- 120

man_tib %>% dplyr::select(Probe_ID, Probe_Type, DESIGN)

call_tibs <- loadCallFiles(files = sam_ss_types[[sample]]$Calls_Path, selKey = 'Probe_ID', 
                          datKey = 'ind_beta', minKey = 'i_poob', minVal = 0.1, retMin=TRUE, addSentrix=TRUE,
                          verbose = opt$verbose)



# %>% split(.$Group)

ret_dat <- plotBetaMatrix_bySample(tib=sam_ss_types, sample=sample, field='Beta', minPval=opt$minPval, 
                                   pvalKey=opt$auto_pval_field, betaKey=opt$auto_beta_field, sort='Poob_Pass_0_Perc',
                                   joinFull=opt$joinFull,
                                   manifest=man_tib, outDir=opt$outDir, 
                                   spread=opt$plotSpread, outType=opt$plotType, dpi=opt$dpi, format=opt$plotFormat,
                                   max=2, maxCnt=opt$plotSub,
                                   verbose=opt$verbose+1,vt=0)



call_tib <- mergeCallsFromSS(ss=labs_ss_tib, max=0, outName=outName, outDir=opt$outDir, 
                             verbose=opt$verbose, vt=1, tc=1, tt=pTracker)

beta_csv <- call_tib %>% dplyr::filter(Method=='ind_poob') %>% head(n=1) %>% dplyr::pull(Full_Path)
beta_tib <- readr::read_csv(beta_csv)

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
