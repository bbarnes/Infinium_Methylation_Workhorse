
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Load sesame:: This causes issues with "ExperimentHub Caching causes a warning"
suppressWarnings(suppressPackageStartupMessages( base::require("sesame") ))
suppressWarnings(suppressPackageStartupMessages( base::require("dbplyr") ))

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

opt$addSentrixID <- FALSE
opt$writeSset    <- FALSE
opt$writeSsum    <- FALSE
opt$writeCalls   <- FALSE
opt$writeSsheet  <- FALSE
opt$writeAuto    <- FALSE

opt$addRawCalls <- FALSE

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
  
  opt$workflows <- 'ind'
  opt$workflows <- "r,i,ind"

  opt$buildSubDir  <- FALSE
  opt$autoDetect   <- FALSE
  opt$writeCalls   <- TRUE
  opt$writeSsheet  <- TRUE
  
  opt$platform   <- 'EPIC'
  opt$manifest   <- 'B4'
  opt$platform   <- NULL
  opt$manifest   <- NULL
  
  par$expRunStr  <- NULL
  par$expChipNum <- NULL
  
  par$local_runType <- 'CORE'
  par$local_runType <- 'EXCBR'
  par$local_runType <- 'qcMVP'
  par$local_runType <- 'GRCm38'
  par$local_runType <- 'COVID'
  par$local_runType <- 'COVIC'
  
  if (par$local_runType=='COVID') {
    par$expRunStr  <- 'COVID-Direct-Set1'
    par$expChipNum <- '204756130014'
    opt$autoDetect <- FALSE
    
    opt$workflows <- 'i,nd,ndi,ind'
    
  } else if (par$local_runType=='COVIC') {
    par$expRunStr  <- 'COVIC-Set1-15052020'
    par$expChipNum <- '204500250013'
    opt$autoDetect <- TRUE
  } else if (par$local_runType=='GRCm38') {
    par$expRunStr <- 'MURMETVEP_mm10_betaTest_06082020'
    par$expRunStr <- 'VanAndel_mm10_betaTest_31082020'
    par$expRunStr <- 'ILMN_mm10_betaTest_17082020'
    opt$autoDetect <- FALSE
  } else if (par$local_runType=='qcMVP') {
    par$expRunStr  <- 'CNTL-Samples_VendA_10092020'
    par$expRunStr  <- 'CNTL-Samples_VendA_10092020_test'
    opt$autoDetect <- TRUE
    opt$dpi <- 72
  } else if (par$local_runType=='CORE') {
    par$expRunStr  <- 'BETA-8x1-EPIC-Core'
    par$expChipNum <- '202761400007'

    par$expRunStr  <- 'ADRN-blood-nonAtopic_EPIC'
    par$expChipNum <- '201125090068'

    par$expRunStr  <- 'GSE122126_EPIC'
    par$expChipNum <- '202410280180'

    par$expRunStr  <- 'EPIC-BETA-8x1-CoreCancer'
    par$expChipNum <- '201502830033'
    opt$autoDetect <- TRUE
  } else if (par$local_runType=='EXCBR') {
    par$expRunStr  <- 'Excalibur-Old-1609202'
    par$expChipNum <- '204076530053'
    par$expChipNum <- '204076530110'
    
    par$expRunStr  <- 'Excalibur-New-1609202'
    par$expChipNum <- '202915460071'
    opt$autoDetect <- TRUE
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unrecognized local_runType={par$local_runType}.{RET}{RET}"))
  }

  opt$idatsDir <- file.path(locIdatDir, paste('idats',par$expRunStr, sep='_') )
  if (!is.null(par$expChipNum)) opt$idatsDir <- file.path(locIdatDir, paste('idats',par$expRunStr, sep='_'),  par$expChipNum)
  opt$auto_sam_csv <- file.path(par$datDir, 'ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo_Mean-COVIC-280-NP-ind_negs-0.02.csv.gz')
  
  opt$outDir <- file.path(par$topDir, 'scratch', par$prgmTag, par$expRunStr, 'manifest-C1')

  par$retData  <- TRUE
  opt$single   <- TRUE
  opt$parallel <- FALSE
  opt$cluster  <- FALSE
  opt$verbose  <- 3
  
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
opt_reqs <- c('outDir','Rscript','verbose')

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
  opt$manDir <- file.path(par$datDir, 'manifest/base')
  mans <- getManifestList(path=opt$manifestPath, platform=opt$platform, manifest=opt$manifest, 
                          dir=opt$manDir, verbose=opt$verbose, tt=pTracker)
  
  auto_opt_tib <- NULL
  auto_ref_tib <- NULL
  auto_can_tib <- NULL
  
  if (opt$autoDetect) {
    if (is.null(opt$auto_sam_csv))
      opt$auto_sam_csv <- file.path(par$datDir, 'ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo_Mean-COVIC-280-NP-ind_negs-0.02.csv.gz')

    cat(glue::glue("[{par$prgmTag}]:{TAB} Loading auto_sam_csv={opt$auto_sam_csv}...{RET}"))
    auto_sam_tib <- suppressMessages(suppressWarnings(readr::read_csv(opt$auto_sam_csv) ))
    cat(glue::glue("[{par$prgmTag}]:{TAB} Done.{RET}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Chip::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  chipTimes <- NULL
  sample_cnt <- length(chipPrefixes)

  try_str <- ''
  if (opt$parallel) {
    par$funcTag <- 'sesamizeSingleSample-Parallel'
    par$retData <- FALSE
    
    cat(glue::glue("[{par$prgmTag}]: parallelFunc={par$funcTag}: samples={sample_cnt}; ",
                   "num_cores={num_cores}, num_workers={num_workers}, Starting...{RET}"))
    # cat(glue::glue("[{par$prgmTag}]: parallelFunc={par$funcTag}: chipPrefixes={RET}"))
    # print(chipPrefixes)
    
    # chipTimes <- foreach (prefix=names(chipPrefixes), .inorder=T, .final = function(x) setNames(x, names(chipPrefixes))) %dopar% {
    chipTimes <- foreach (prefix=names(chipPrefixes), .combine = rbind) %dopar% {
      rdat <- NULL
      rdat <- sesamizeSingleSample(prefix=chipPrefixes[[prefix]], man=mans, ref=auto_sam_tib, opt=opt, 
                                   retData=par$retData, workflows=workflows_vec, tc=1)
      rdat
    }
    cat(glue::glue("[{par$prgmTag}] parallelFunc={par$funcTag}: Done.{RET}{RET}"))
  } else {
    par$funcTag <- 'sesamizeSingleSample-Linear'
    
    cat(glue::glue("[{par$prgmTag}]: linearFunc={par$funcTag}: samples={sample_cnt}.{RET}"))
    
    rdat <- NULL
    for (prefix in names(chipPrefixes)) {
      cat(glue::glue("[{par$prgmTag}]: linearFunc={par$funcTag}: Starting; prefix={prefix}...{RET}"))
      
      par$retData <- TRUE
      opt$verbose <- 3
      opt$verbose <- 6
      # workflows_vec <- c('r', 'i', 'ind')
      
      rdat <- sesamizeSingleSample(prefix=chipPrefixes[[prefix]], man=mans, ref=auto_sam_tib, opt=opt, 
                                   retData=par$retData, workflows=workflows_vec, tc=1)
      
      idats <- sesamizeSingleSample(prefix=chipPrefixes[[prefix]], man=mans, ref=auto_sam_tib, opt=opt, 
                                    retData=par$retData, workflows=workflows_vec, tc=1)
      
      if (FALSE) {
        
        cur_sset <- rdat$cur_sset
        ses_man_tib <- rdat$sman
        ww <- 3
        cur_workflow <- 'ind'
        
        cur_sset_csv <- NULL
        cur_ssum_csv <- NULL
        
        sum_list <- summarySSET_workflow(
          sset=cur_sset, man=ses_man_tib, idx=ww, workflow=cur_workflow,
          writeSset=opt$writeSset,sset_csv=cur_sset_csv,
          writeSsum=opt$writeSsum,ssum_csv=cur_ssum_csv,
          minNegPval=opt$minNegPval,minOobPval=opt$minOobPval,
          percisionBeta=opt$percisionBeta, 
          percisionPval=opt$percisionPval, 
          percisionSigs=opt$percisionSigs,
          verbose=4)
        
        
        rdat$cur_sset@extra$pvals$pOOBAH %>% head()
        rdat$cur_sset@extra$pvals$pOOBAH <- NULL
        
        new_sset <- rdat$raw_sset
        new_sset@extra$pvals$pOOBAH <- NULL
        new_sset@extra$pvals <- NULL
        new_sset <- pOOBAH2(new_sset, force=TRUE)
        new_sset@extra$pvals$pOOBAH %>% head(n=30)

        two_sset <- rdat$cur_sset
        two_sset@extra$pvals$pOOBAH <- NULL
        two_sset@extra$pvals <- NULL
        two_sset <- pOOBAH2(two_sset, force=TRUE)
        two_sset@extra$pvals$pOOBAH %>% head(n=30)
        
        rdat$raw_sset@extra$pvals$pOOBAH %>% head(n=30)
      }
      
      
      cat(glue::glue("[{par$prgmTag}]: linearFunc={par$funcTag}: try_str={try_str}. Done.{RET}{RET}"))
      if (opt$single) break
    }
    cat(glue::glue("[{par$prgmTag}] parallelFunc={par$funcTag}: Done.{RET}{RET}"))
  }
  
  opt$opt_csv  <- file.path(opt$outDir, paste(par$prgmTag,'program-options.csv', sep='.') )
  opt$par_csv  <- file.path(opt$outDir, paste(par$prgmTag,'program-parameters.csv', sep='.') )
  opt$time_csv <- file.path(opt$outDir, paste(par$prgmTag,'time-tracker.csv.gz', sep='.') )
  
  opt_tib  <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value") %>% dplyr::arrange(Option)
  par_tib  <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value") %>% dplyr::arrange(Params)
  # time_tib <- pTracker$time %>% dplyr::mutate_if(is.numeric, list(round), 4)
  
  readr::write_csv(opt_tib, opt$opt_csv)
  readr::write_csv(par_tib, opt$par_csv)
  # readr::write_csv(time_tib, opt$time_csv)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
