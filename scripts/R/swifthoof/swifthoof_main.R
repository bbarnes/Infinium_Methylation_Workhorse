
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))
suppressWarnings(suppressPackageStartupMessages( base::require("grid") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

# Load sesame:: This causes issues with "ExperimentHub Caching causes a warning"
#  suppressWarnings(suppressPackageStartupMessages( base::require("sesame") ))

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
par <- NULL
opt <- NULL

par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'swifthoof'
par$prgmTag <- paste(par$prgmDir,'main', sep='_')
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# Illumina based directories::
par$macDir <- '/Users/bbarnes/Documents/Projects/methylation/tools'
par$lixDir <- '/illumina/scratch/darkmatter'

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
opt$fresh       <- FALSE
opt$buildSubDir <- FALSE
opt$autoDetect  <- FALSE
opt$skipSwap    <- FALSE
opt$workflows   <- NULL

# Output Options::
opt$loadIdat    <- FALSE
opt$saveIdat    <- FALSE

opt$loadSsets   <- FALSE
opt$saveSsets   <- FALSE
opt$saveRawSset <- FALSE
opt$lightFootPrint <- FALSE

opt$addSentrixID <- FALSE
opt$writeSset    <- FALSE
opt$writeSsum    <- FALSE
opt$writeCalls   <- FALSE
opt$writeSsheet  <- FALSE
opt$writeAuto    <- FALSE

opt$addRawCalls <- FALSE

# Reporting Options::
opt$sigs_sum_field <- NULL

# Threshold Options::
opt$minNegPval   <- 0.02
opt$minOobPval   <- 0.1
opt$minNegPerc   <- 98
opt$minOobPerc   <- 90
opt$minDeltaBeta <- 0.2

opt$percisionSigs <- 1
opt$percisionBeta <- 4
opt$percisionPval <- 6

# Parallel/Cluster Options::
opt$single   <- FALSE
opt$parallel <- FALSE
opt$cluster  <- FALSE

# Plotting Options::
opt$plotSset  <- FALSE
opt$plotCalls <- FALSE
opt$plotAuto  <- FALSE

opt$plotFormat <- 'pdf'
opt$plotFormat <- 'png'

opt$dpi <- 72
opt$dpi <- 120

opt$plotMax <- 10000
opt$plotSub <- 5000

# verbose Options::
opt$verbose <- 3

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
  
  opt$single   <- TRUE
  opt$cluster  <- FALSE
  opt$parallel <- FALSE
  opt$single   <- FALSE
  opt$parallel <- TRUE
  opt$parallel <- FALSE
  
  opt$workflows <- 'ind'
  
  opt$buildSubDir  <- FALSE
  opt$autoDetect   <- TRUE
  opt$writeCalls   <- TRUE
  opt$writeSsheet  <- TRUE
  
  par$retData <- TRUE
  
  # isCORE  <- TRUE
  # isCOVIC <- TRUE
  # isCOVIC <- FALSE
  # if (isCOVIC) {
  #   opt$platform   <- 'EPIC'
  #   opt$manifest   <- 'C0'
  #   
  #   opt$platform   <- NULL
  #   opt$manifest   <- NULL
  #   
  #   # Set-1
  #   opt$expRunStr  <- 'idats_COVIC-Set1-15052020'
  #   opt$expChipNum <- '204500250013'
  #   
  # } else if (isCORE) {
  #   opt$platform   <- 'EPIC'
  #   opt$manifest   <- 'B4'
  #   
  #   opt$expRunStr  <- 'idats_BETA-8x1-EPIC-Core'
  #   opt$expChipNum <- '202761400007'
  # 
  #   opt$expRunStr  <- 'idats_ADRN-blood-nonAtopic_EPIC'
  #   opt$expChipNum <- '201125090068'
  #   
  #   opt$expRunStr  <- 'idats_GSE122126_EPIC'
  #   opt$expChipNum <- '202410280180'
  # 
  #   opt$expRunStr  <- 'idats_EPIC-BETA-8x1-CoreCancer'
  #   opt$expChipNum <- '201502830033'
  #   
  # } else {
  #   stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Unsupported pre-defined method! Exiting...{RET}{RET}"))
  # }
  # opt$idatsDir <- file.path('/Users/bbarnes/Documents/Projects/methylation/data/idats', opt$expRunStr, opt$expChipNum)

  locIdatDir <- '/Users/bbarnes/Documents/Projects/methylation/data/idats'
  
  opt$expRunStr  <- 'ReferenceBETA'
  opt$expRunStr  <- 'idats_COVIC-Set1-15052020'
  opt$expChipNum <- '204500250013'
  
  opt$idatsDir <- file.path(locIdatDir, opt$expRunStr)
  opt$idatsDir <- file.path(locIdatDir, opt$expRunStr, opt$expChipNum)
  opt$auto_sam_csv <- file.path(par$datDir, 'ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo_Mean-COVIC-280-NP-ind_negs-0.02.csv.gz')
  
  # mm10
  par$isMouse <- FALSE
  if (par$isMouse) {
    opt$idatsDir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/idats/ILMN_mm10_betaTest_17082020'
    opt$auto_sam_csv <- NULL
    opt$expRunStr  <- 'mm10'
    opt$autoDetect <- FALSE
    opt$cluster  <- FALSE
    opt$single   <- TRUE
    opt$parallel <- FALSE
  }
  par$isMVP <- TRUE
  if (par$isMVP) {
    opt$expRunStr  <- 'CNTL-Samples_VendA_10092020'
    opt$autoDetect <- TRUE
    opt$plotAuto   <- TRUE
    opt$cluster  <- FALSE
    opt$single   <- TRUE
    opt$parallel <- FALSE
    
    opt$dpi <- 72
    opt$idatsDir <- file.path('/Users/bbarnes/Documents/Projects/methylation/VA_MVP/idats',opt$expRunStr)
  }

  opt$verbose <- 3
  opt$verbose <- 6
  opt$verbose <- 60
  
  opt$outDir <- file.path(par$topDir, par$prgmDir, opt$expRunStr)
  
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
    
    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    make_option(c("-i", "--idatsDir"), type="character", default=opt$idatsDir, 
                help="idats directory [default= %default]", metavar="character"),

    # Optional Files::
    make_option(c("--manifestPath"), type="character", default=opt$manifestPath,
                help="Path to manfifest (CSV) otherwise use dat [default= %default]", metavar="character"),
    # make_option(c("--addressPath"), type="character", default=opt$addressPath,
    #             help="Path to address (RDS) otherwise use dat [default= %default]", metavar="character"),
    make_option(c("--subManifest"), action="store_true", default=opt$subManifest,
                help="Boolean variable to use subset manifest instead of subset. [default= %default]", metavar="boolean"),
    make_option(c("--auto_sam_csv"), type="character", default=opt$auto_sam_csv,
                help="Path to auto detect beta values (CSV) [default= %default]", metavar="character"),
    
    # Platform/Method Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Forced platform [EPIC, 450k, 27k, NZT] otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--manifest"), type="character", default=opt$manifest, 
                help="Forced manifest [B1, B2, B4] otherwise auto-detect [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--fresh"), action="store_true", default=opt$fresh,
                help="Boolean variable to build fresh version of database files [default= %default]", metavar="boolean"),
    make_option(c("--buildSubDir"), action="store_true", default=opt$buildSubDir,
                help="Boolean variable to build subdirectories based on Chip/BeadPool (for R&D purposes) [default= %default]", metavar="boolean"),
    make_option(c("--autoDetect"), action="store_true", default=opt$autoDetect,
                help="Boolean variable to auto detect reference samples. Must provide reference samples. [default= %default]", metavar="boolean"),
    make_option(c("--skipSwap"), action="store_true", default=opt$skipSwap,
                help="Boolean variable to skpping appending swap percentages to sample sheet. [default= %default]", metavar="boolean"),
    make_option(c("--workflows"), type="character", default=opt$workflows,
                help="Order of operations comma seperated [ raw,ind,ndi,din ] [default= %default]", metavar="character"),
    # make_option(c("--sampleSheet"), type="character", default=opt$sampleSheet, 
    #             help="Target Sample Sheet containing samples/chips to ONLY analyze [default= %default]", metavar="character"),
    
    # Output Options::
    make_option(c("--loadIdat"), action="store_true", default=opt$loadIdat,
                help="Boolean variable to load existing IDAT from RDS file [default= %default]", metavar="boolean"),
    make_option(c("--saveIdat"), action="store_true", default=opt$saveIdat,
                help="Boolean variable to write IDAT RDS file [default= %default]", metavar="boolean"),
    
    make_option(c("--loadSsets"), action="store_true", default=opt$loadSsets,
                help="Boolean variable to load existing Signal Set from RDS file [default= %default]", metavar="boolean"),
    make_option(c("--saveSsets"), action="store_true", default=opt$saveSsets,
                help="Boolean variable to write Signal Set RDS file [default= %default]", metavar="boolean"),
    make_option(c("--saveRawSset"), action="store_true", default=opt$saveRawSset,
                help="Boolean variable to write Raw Signal Set RDS file [default= %default]", metavar="boolean"),
    make_option(c("--lightFootPrint"), action="store_true", default=opt$lightFootPrint,
                help="Boolean variable to NOT save any RDS files [default= %default]", metavar="boolean"),
    
    make_option(c("--addSentrixID"), action="store_true", default=opt$addSentrixID,
                help="Boolean variable to add Sentrix Name to calls output columns [default= %default]", metavar="boolean"),
    make_option(c("--writeSset"), action="store_true", default=opt$writeSset,
                help="Boolean variable to write Signal Set file [default= %default]", metavar="boolean"),
    make_option(c("--writeSsum"), action="store_true", default=opt$writeSsum,
                help="Boolean variable to write Signal Set Summary file [default= %default]", metavar="boolean"),
    make_option(c("--writeCalls"), action="store_true", default=opt$writeCalls,
                help="Boolean variable to write Calls (Pval/Beta) file [default= %default]", metavar="boolean"),
    make_option(c("--writeSsheet"), action="store_true", default=opt$writeSsheet,
                help="Boolean variable to Sample Sheet file [default= %default]", metavar="boolean"),
    make_option(c("--writeAuto"), action="store_true", default=opt$writeAuto,
                help="Boolean variable to write Auto-Detection Matricies (Pval/Beta) file [default= %default]", metavar="boolean"),
    make_option(c("--addRawCalls"), action="store_true", default=opt$addRawCalls,
                help="Boolean variable to output raw calls [default= %default]", metavar="boolean"),
    
    # Reporting Options::
    make_option(c("--sigs_sum_field"), type="character", default=opt$sigs_sum_field, 
                help="Signal summary field in AutoSampleSheet [default= %default]", metavar="character"),
    
    # Threshold Options::
    make_option(c("--minNegPval"), type="double", default=opt$minNegPval, 
                help="Minimum passing detection p-value using Negative Controls [default= %default]", metavar="double"),
    make_option(c("--minOobPval"), type="double", default=opt$minOobPval,
                help="Minimum passing detection p-value using Negative Out-Of-Band [default= %default]", metavar="double"),
    
    make_option(c("--minNegPerc"), type="double", default=opt$minNegPerc, 
                help="Minimum percentage of loci passing detection p-value using Negative Controls to flag Requeue of sample. [default= %default]", metavar="double"),
    make_option(c("--minOobPerc"), type="double", default=opt$minOobPerc, 
                help="Minimum percentage of loci passing detection p-value using Out-Of-Band to flag Requeue of sample. [default= %default]", metavar="double"),
    
    make_option(c("--minDeltaBeta"), type="double", default=opt$minDeltaBeta,
                help="Minimum passing delta-beta. Used in AutoSampleSheet cacluclations [default= %default]", metavar="double"),
    
    make_option(c("--percisionSigs"), type="integer", default=opt$percisionSigs,
                help="Rounding percision for signal values in calls output files [default= %default]", metavar="double"),
    make_option(c("--percisionBeta"), type="integer", default=opt$percisionBeta,
                help="Rounding percision for beta values in calls output files [default= %default]", metavar="double"),
    make_option(c("--percisionPval"), type="integer", default=opt$percisionPval,
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
    
    # verbose::
    make_option(c("-v", "--verbose"), type="integer", default=opt$verbose, 
                help="0-5 (5 is very verbose) [default= %default]", metavar="integer")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Validate Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (is.null(par$runMode) || is.null(par$prgmTag) || is.null(par$scrDir) || is.null(par$datDir)) {
  
  par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
  par_tib %>% base::print(n=base::nrow(par_tib) )
  
  if (is.null(par$runMode)) cat(glue::glue("[Usage]: runMode is NULL!!!{RET}"))
  if (is.null(par$prgmTag)) cat(glue::glue("[Usage]: prgmTag is NULL!!!{RET}"))
  if (is.null(par$scrDir))  cat(glue::glue("[Usage]: scrDir is NULL!!!{RET}"))
  if (is.null(par$datDir))  cat(glue::glue("[Usage]: darDir is NULL!!!{RET}"))
  base::stop("Null Parameters!\n\n")
}

if (is.null(opt$outDir) || is.null(opt$idatsDir) ||
    is.null(opt$Rscript) ||
    is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )
  
  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )

  if (is.null(opt$outDir))    cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$idatDir))   cat(glue::glue("[Usage]: idatDir is NULL!!!{RET}"))
  if (is.null(opt$Rscript))   cat(glue::glue("[Usage]: Rscript is NULL!!!{RET}"))
  
  if (is.null(opt$verbose))  cat(glue::glue("[Usage]: verbose is NULL!!!{RET}"))
  base::stop(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
}
par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
if (opt$verbose>=1) par_tib %>% base::print(n=base::nrow(par_tib) )
if (opt$verbose>=1) opt_tib %>% base::print(n=base::nrow(opt_tib) )

cat(glue::glue("[{par$prgmTag}]: Done. Validating Options.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )

par$man_src_dir <- file.path(par$scrDir, 'manifests/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$man_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$man_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$man_src_dir}!{RET}{RET}") )

par$prgm_src_dir <- file.path(par$scrDir,par$prgmDir, 'functions')
if (!dir.exists(par$prgm_src_dir)) stop(glue::glue("[{par$prgmTag}]: Program Source={par$prgm_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$prgm_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Program Source={par$prgm_src_dir}!{RET}{RET}") )

cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (!dir.exists(opt$idatsDir)) stop(glue::glue("{RET}[{par$prgmTag}]: idatsDir={opt$idatsDir} does not exist!!!{RET}{RET}"))
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
cat(glue::glue("[{par$prgmTag}]: Output Directory (TOP)={opt$outDir}...{RET}"))

workflows_vec <- NULL
if (!is.null(opt$workflows)) workflows_vec <- opt$workflows %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

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
  if (dir.exists(par$lixDir)) {
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
    readr::write_file(cmd_full, path=runShell)
    Sys.chmod(runShell, mode="0777")
    
    # Add cluster execute if avialbel (i.e. linux)
    cmd_lan <- paste0(runShell,RET, sep='')
    if (par$isLinux)
      cmd_lan <- paste(par$lan_exe, paste('imWH',chipName, sep='_'), runShell,RET, sep=' ')
    readr::write_file(cmd_lan, path=lanShell)
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
  
  pTracker <- timeTracker$new(verbose=opt$verbose)
  
  # opt$outDir <- file.path(opt$outDir, par$prgmTag)
  # if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
  cat(glue::glue("[{par$prgmTag}]:{TAB} Output Directory (SIG)={opt$outDir}...{RET}"))
  
  mans <- NULL
  opt$manDir <- file.path(par$datDir, 'manifest/base')
  mans <- getManifestList(path=opt$manifestPath, platform=opt$platform, manifest=opt$manifest, 
                          dir=opt$manDir, verbose=opt$verbose, tt=pTracker)
  
  auto_opt_tib <- NULL
  auto_ref_tib <- NULL
  auto_can_tib <- NULL
  
  if (opt$autoDetect) {
    cat(glue::glue("[{par$prgmTag}]:{TAB} Loading auto_sam_csv={opt$auto_sam_csv}...{RET}"))
    auto_sam_tib <- suppressMessages(suppressWarnings(readr::read_csv(opt$auto_sam_csv) ))
    cat(glue::glue("[{par$prgmTag}]:{TAB} Done.{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Chip::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  chipTimes <- NULL
  
  # opt_csv  <- file.path(opt$outDir, 'program-options.csv')
  # par_csv  <- file.path(opt$outDir, 'program-parameters.csv')
  # readr::write_csv(opt_tib, opt_csv)
  # readr::write_csv(par_tib, par_csv)
  
  if (opt$parallel) {
    funcTag <- 'sesamizeSingleSample-Parallel'
    par$retData <- FALSE
    
    cat(glue::glue("[{par$prgmTag}]: parallelFunc={funcTag}: num_cores={num_cores}, num_workers={num_workers}, Starting...{RET}"))

    # chipTimes <- foreach (prefix=names(chipPrefixes), .inorder=T, .final = function(x) setNames(x, names(chipPrefixes))) %dopar% {
    chipTimes <- foreach (prefix=names(chipPrefixes), .combine = rbind) %dopar% {
      rdat <- NULL
      try_str <- ''
      rdat = tryCatch({
        try_str <- 'Pass'
        sesamizeSingleSample(prefix=chipPrefixes[[prefix]], man=mans, ref=auto_sam_tib, opt=opt, 
                             retData=par$retData, workflows=workflows_vec, tc=3)
      }, warning = function(w) {
        try_str <- paste('warning',funcTag, sep='-')
        rdat <- NA
      }, error = function(e) {
        try_str <- paste('error',funcTag, sep='-')
        rdat <- NA
      }, finally = {
        try_str <- paste('cleanup',funcTag, sep='-')
        rdat <- NA
      })
      cat(glue::glue("[{par$prgmTag}]: parallelFunc={funcTag}: try_str={try_str}. Done.{RET}{RET}"))
      
      rdat
    }
    
  } else {
    funcTag <- 'sesamizeSingleSample-Linear'

    # opt$skipSwap
    cat(glue::glue("[{par$prgmTag}]: linearFunc={funcTag}: Starting...{RET}"))
    
    for (prefix in names(chipPrefixes)) {
      rdat <- NULL
      try_str <- ''
      rdat = tryCatch({
        try_str <- 'Pass'
        sesamizeSingleSample(prefix=chipPrefixes[[prefix]], man=mans, ref=auto_sam_tib, opt=opt, 
                             retData=par$retData, workflows=workflows_vec, tc=3)
      }, warning = function(w) {
        try_str <- paste('warning',funcTag, sep='-')
        rdat <- NA
      }, error = function(e) {
        try_str <- paste('error',funcTag, sep='-')
        rdat <- NA
      }, finally = {
        try_str <- paste('cleanup',funcTag, sep='-')
        rdat <- NA
      })
      cat(glue::glue("[{par$prgmTag}]: linearFunc={funcTag}: try_str={try_str}. Done.{RET}{RET}"))
      
      if (opt$single) break
    }

    cat(glue::glue("[{par$prgmTag}] parallelFunc={funcTag}: Done.{RET}{RET}"))
  }
  
  opt$opt_csv  <- file.path(opt$outDir, paste(par$prgmTag,'program-options.csv', sep='.') )
  opt$par_csv  <- file.path(opt$outDir, paste(par$prgmTag,'program-parameters.csv', sep='.') )
  opt$time_csv <- file.path(opt$outDir, paste(par$prgmTag,'time-tracker.csv.gz', sep='.') )
  
  opt_tib  <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value") %>% dplyr::arrange(Option)
  par_tib  <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value") %>% dplyr::arrange(Params)
  time_tib <- pTracker$time %>% dplyr::mutate_if(is.numeric, list(round), 4)
  
  readr::write_csv(opt_tib, opt$opt_csv)
  readr::write_csv(par_tib, opt$par_csv)
  readr::write_csv(time_tib, opt$time_csv)

  # pTracker_tib <- pTracker$time %>% dplyr::mutate_if(is.numeric, list(round), 4)
  # time_csv <- file.path(opt$outDir, 'time-tracker.csv.gz')
  # readr::write_csv(pTracker_tib, time_csv)
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
