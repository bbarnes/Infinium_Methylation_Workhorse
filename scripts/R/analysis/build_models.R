
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Build COVIC Model From Merged Builds::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Plotting Packages::
suppressWarnings(suppressPackageStartupMessages( base::require("corrplot") ))
suppressWarnings(suppressPackageStartupMessages( base::require("GGally") ))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("plyr")) )
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

# Rcpp::
suppressWarnings(suppressPackageStartupMessages( require("Rcpp") ))

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
par$prgmTag <- 'build_models'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# TBD:: This variable needs to be properly integrated, its a total hack
#  currently...
par$classVar <- NULL
par$joinType <- "left"

# Directory Parameters::
opt$outDir <- NULL
opt$datDir <- NULL

opt$addPval     <- FALSE
opt$buildDml    <- TRUE
opt$buildDbl    <- FALSE
opt$buildModels <- FALSE

# Run Parameters::
opt$runName   <- NULL

# Class Parameters::
# Really simple test to make sure we can seperate the sexes...
opt$classVar <- NULL
opt$workflow <- NULL

opt$trainClass <- NULL
opt$cross_perc_min <- 90

# Sample Level Filtering Parameters::
opt$samplePvalName <- NULL
opt$samplePvalPerc <- NULL

# Loci Level Filtering Parameters::
opt$lociBetaKey <- NULL
opt$lociPvalKey <- NULL
opt$lociPvalMin <- NULL

# Training Parameters::
opt$alphaMin <- 0.0
opt$alphaMax <- 1.0
opt$alphaInc <- 0.2
opt$alphaInc <- 0.5

# Pre-defined Loci (Feature) Selection File Parameters::
opt$featuresCsv <- NULL

# DML-defined Loci (Feature) Selection Size Parameters
opt$featuresDml <- NULL
opt$featuresDbl <- NULL

# Reproducible Seed Parameters::
opt$seeds <- NULL

# Output Format Parameters::
opt$percisionBeta <- 4
opt$percisionPval <- 6

# Plot parameters::
opt$plot_pairs <- FALSE

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
  
  opt$addPval     <- TRUE
  opt$buildDml    <- FALSE
  opt$buildDbl    <- FALSE
  opt$buildModels <- FALSE
  
  opt$buildDml    <- TRUE
  opt$buildDbl    <- TRUE
  
  opt$single   <- TRUE
  opt$single   <- FALSE
  opt$cluster  <- FALSE
  opt$parallel <- FALSE
  
  opt$samplePvalName <- "Poob_Pass_0_Perc"
  opt$samplePvalName <- "cg_pvals_pOOBAH_pass_perc_1"
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
  
  opt$featuresCsv <- NULL
  opt$featuresDml <- NULL
  
  opt$classVar <- 'Sample_Class'
  opt$workflow <- "ind"
  par$platform <- 'EPIC'
  par$version  <- 'B4'
  
  opt$plot_pairs <- FALSE
  
  opt$clean <- FALSE
  opt$clean <- TRUE
  
  # opt$trainClass <- paste('HELA','JURKAT','MCF7','RAJI', sep=',')
  # opt$trainClass <- paste('Xa','XaXaY','XaXi','XaXiY','XaY', sep=',')
  # opt$trainClass <- paste('XaXi','XaY', sep=',')
  # opt$trainClass <- paste('nSARSCov2', 'pSARSCov2', sep=',')
  
  #
  # Pre-defined local options runTypes::
  #
  par$local_runType <- 'qcMVP'
  par$local_runType <- 'CORE'
  par$local_runType <- 'EXCBR'
  par$local_runType <- 'NZT'
  par$local_runType <- 'COVIC'
  par$local_runType <- 'GENK'
  par$local_runType <- 'GRCm38'
  par$local_runType <- 'COVID'
  par$local_runType <- 'EPIC-8x1-EM-Sample-Prep'
  par$local_runType <- 'Chicago-Ober-Custom'
  par$local_runType <- 'NA12878'
  
  if (FALSE) {
    
  } else if (par$local_runType=='Chicago-Ober-Custom') {
    opt$runName  <- par$local_runType
    opt$workflow <- "ind"
    par$platform <- NULL
    par$version  <- NULL
    
    # opt$parallel <- TRUE
    opt$single   <- FALSE
    opt$buildDml <- TRUE
    opt$buildDbl <- TRUE
    
    opt$lociBetaKey <- "betas"
    opt$lociPvalKey <- "pvals_pOOBAH"
    
    opt$classVar <- 'AutoSample_R2_Key_2'
    par$classVar <- 'AutoSample_R2_Key_2'
    
    par$vstr <- ""
    par$vstr <- "v1"
    par$vstr <- "v2"
    opt$runName <- paste(opt$runName,par$vstr, sep="-")
    
    opt$manDirPath <- file.path(par$topDir, "data/manifests/methylation/Chicago-Ober-Custom")

    par$platform <- NULL
    if (!is.null(par$platform)) opt$runName <- paste(opt$runName,par$platform, sep='-')
    
    opt$datDir <- paste(
      file.path(par$topDir,"scratch",par$runMode,"merge_builds",opt$runName,
                # par$platform,par$version,
                # opt$classVar,
                par$classVar,
                opt$workflow),
      sep=',')
    
    opt$trainClass <- paste('ChicagoA1','ChicagoA2', sep=',')
    
    opt$samplePvalPerc <- 50

  } else if (par$local_runType=='EPIC-8x1-EM-Sample-Prep') {
    opt$runName  <- par$local_runType
    opt$classVar <- 'Sample_Class'
    opt$workflow <- "ind"
    par$platform <- 'EPIC'
    par$version  <- 'B4'
    
    # opt$parallel <- TRUE
    opt$single   <- FALSE
    opt$buildDml <- TRUE
    opt$buildDbl <- TRUE

    opt$lociBetaKey <- "betas"
    opt$lociPvalKey <- "pvals_pOOBAH"
    
    # opt$datDir <- paste(
    #   file.path(par$topDir,"scratch","merge_builds_latest/merge_builds",opt$runName,
    #             par$platform,par$version,opt$classVar,opt$workflow),
    #   file.path(par$topDir,"scratch","merge_builds_latest/merge_builds","EPIC-8x1-EM-Sample-Prep.v0",
    #             par$platform,par$version,opt$classVar,opt$workflow),
    #   sep=',')
    
    opt$datDir <- paste(
      file.path(par$topDir, "data/CustomContent/EPIC-8x1-EM-Sample-Prep/docker-v.1.29/merge_builds/EPIC-8x1-EM-Sample-Prep/EPIC/B4"),
      sep=',')
    
    # opt$trainClass <- paste('BS','EM', sep=',')
    opt$trainClass <- paste('BS-FF','BS-PE', 'EM-FF', 'EM-PE', sep=',')
    
  } else if (par$local_runType=='NA12878') {
    opt$runName  <- par$local_runType
    opt$workflow <- "ind"
    par$platform <- NULL
    par$version  <- NULL
    
    # opt$parallel <- TRUE
    opt$single     <- FALSE
    opt$single     <- TRUE
    opt$plot_pairs <- FALSE
    opt$plot_pairs <- TRUE
    
    opt$buildDml <- TRUE
    opt$buildDbl <- TRUE
    opt$buildDml <- FALSE
    opt$buildDbl <- FALSE
    
    opt$lociBetaKey <- "betas"
    opt$lociPvalKey <- "pvals_pOOBAH"
    
    par$joinType <- "inner"
    
    par$classVar <- "detect_platform"
    par$classVar <- "detect_version"
    opt$classVar <- "Sample_Perc"
    
    par$platform <- NULL
    par$platform <- "Rand3"
    par$platform <- "Rand2"
    par$platform <- "Rand1"
    par$platform <- "NA12878-granular-3-rand1"
    par$platform <- "NA12878-granular-3-rand1-2-3"
    
    # if (!is.null(par$platform)) opt$runName <- paste(opt$runName,par$platform, sep='-')

    opt$datDir <- paste(
      file.path(par$topDir, "scratch/RStudio/merge_builds", par$platform,par$classVar,opt$workflow),
      # file.path(par$topDir,"scratch",par$runMode,"merge_builds",opt$runName,
      #           # par$platform,par$version,
      #           # opt$classVar,
      #           par$classVar,
      #           opt$workflow),
      sep=',')
    
    opt$trainClass <- paste('NA12878', sep=',')

  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  # opt$outDir <- file.path(par$topDir, 'scratch', par$prgmTag, par$platform, par$version)
  # opt$outDir <- file.path(par$topDir, 'scratch', par$prgmTag)
  opt$outDir <- file.path(par$topDir, 'scratch', par$runMode)
  
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
                help="List of Merged Swifthoof Build Directory(s), commas seperated [default= %default]", metavar="character"),

    make_option(c("--addPval"), action="store_true", default=opt$addPval, 
                help="Boolean variable to write and return pval matrix [default= %default]", metavar="boolean"),
    
    make_option(c("--buildDbl"), action="store_true", default=opt$buildDbl, 
                help="Boolean variable to build delta beta (needs Rcpp) [default= %default]", metavar="boolean"),
    make_option(c("--buildDml"), action="store_true", default=opt$buildDml, 
                help="Boolean variable to build DML [default= %default]", metavar="boolean"),
    make_option(c("--buildModels"), action="store_true", default=opt$buildModels, 
                help="Boolean variable to build Models [default= %default]", metavar="boolean"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    
    # Chip Platform and Version Parameters::
    # make_option(c("--platform"), type="character", default=opt$platform, 
    #             help="Platform name (HM50, EPIC) [default= %default]", metavar="character"),
    # make_option(c("--version"), type="character", default=opt$version, 
    #             help="Manifest version (B2, B4, C0) [default= %default]", metavar="character"),
    
    # Class Parameters::
    make_option(c("--classVar"), type="character", default=opt$classVar, 
                help="Classification Variable Name [default= %default]", metavar="character"),
    make_option(c("--workflow"), type="character", default=opt$workflow, 
                help="Target Workflow Variable Name [default= %default]", metavar="character"),
    
    make_option(c("--trainClass"), type="character", default=opt$trainClass, 
                help="Training Class Variable Names, comma seperated [default= %default]", metavar="character"),
    make_option(c("--cross_perc_min"), type="double", default=opt$cross_perc_min,
                help="Minimum cross valudation percent to save model [default= %default]", metavar="double"),
    
    # Sample Level Filtering Parameters::
    make_option(c("--samplePvalName"), type="character", default=opt$samplePvalName, 
                help="Pval Method Name to filter Samples for training [default= %default]", metavar="character"),
    make_option(c("--samplePvalPerc"), type="double", default=opt$samplePvalPerc,
                help="Pval Min Percent Passing to filter Sample for training [default= %default]", metavar="double"),
    
    # Loci Level Filtering Parameters::
    make_option(c("--lociBetaKey"), type="character", default=opt$lociBetaKey, 
                help="Loci Beta-Method Name (key) for training. Comma seperated list. [default= %default]", metavar="character"),
    make_option(c("--lociPvalKey"), type="character", default=opt$lociPvalKey, 
                help="Loci Pval-Method Name (key) for filtering for training. Comma seperated list. [default= %default]", metavar="character"),
    make_option(c("--lociPvalMin"), type="character", default=opt$lociPvalMin,
                help="Pval Min Passing for loci filtering for training. Comma seperated list. [default= %default]", metavar="character"),
    
    # Training Parameters::
    make_option(c("--alphaMin"), type="double", default=opt$alphaMin,
                help="Minimum alpha value for Elastic Net training [default= %default]", metavar="double"),
    make_option(c("--alphaMax"), type="double", default=opt$alphaMax,
                help="Maximum alpha value for Elastic Net training [default= %default]", metavar="double"),
    make_option(c("--alphaInc"), type="double", default=opt$alphaInc,
                help="Incremental alpha value for Elastic Net training [default= %default]", metavar="double"),
    
    # Pre-defined Loci (Feature) Selection File Parameters::
    make_option(c("--featuresCsv"), type="character", default=opt$featuresCsv, 
                help="Pre-defined Loci feature selection file for training. Comma seperated list. [default= %default]", metavar="character"),
    
    # DML-defined Loci (Feature) Selection Size Parameters
    make_option(c("--featuresDml"), type="character", default=opt$featuresDml,
                help="DML-defined Loci feature selection counts for training. Comma seperated list. [default= %default]", metavar="character"),
    make_option(c("--featuresDbl"), type="character", default=opt$featuresDbl,
                help="dBL-defined Loci feature selection counts for training. Comma seperated list. [default= %default]", metavar="character"),
    
    # Reproducible Seed Parameters::
    make_option(c("--seeds"), type="character", default=opt$seeds,
                help="List of seeds for machine learning training reproducibility. Comma seperated list. [default= %default]", metavar="character"),
    
    # Output Format Parameters::
    make_option(c("--percisionBeta"), type="integer", default=opt$percisionBeta,
                help="Rounding percision for beta values in calls output files [default= %default]", metavar="integer"),
    make_option(c("--percisionPval"), type="integer", default=opt$percisionPval,
                help="Rounding percision for detection p-values in calls output files [default= %default]", metavar="integer"),
    
    # Plot parameters::
    make_option(c("--plot_pairs"), action="store_true", default=opt$plot_pairs, 
                help="Boolean variable to plot pairs (mostly testing stuff) [default= %default]", metavar="boolean"),
    
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
#                         Program Initialization::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_reqs <- c('runMode','prgmTag','scrDir','datDir','exePath')
opt_reqs <- c('outDir','Rscript','verbose')
opt_reqs <- c('outDir','Rscript','verbose','clean',
              'datDir','runName','classVar','workflow',
              'trainClass', 'cross_perc_min', 
              'samplePvalName', 'samplePvalPerc', 
              'lociBetaKey', 'lociPvalKey', 'lociPvalMin', 
              'alphaMin','alphaMax','alphaInc','seeds', 
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
                    libs=TRUE,rcpp=TRUE,
                    verbose=opt$verbose,vt=3,tc=0,tt=NULL)

opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Preprocessing:: System Params
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Not sure if this is needed anymore...
#
if (FALSE) {
  opt <- setLaunchExe(opts=opt, pars=par, verbose=opt$verbose, vt=5,tc=0)
  if (!opt$isLinux && opt$buildDbl) {
    suppressWarnings(suppressPackageStartupMessages( require("Rcpp") ))
    
    par$sourceCpp <- file.path(par$scrDir, 'R/Rcpp/cpgLociVariation.cpp')
    if (!file.exists(par$sourceCpp)) par$sourceCpp <- file.path(par$scrDir, 'Rcpp/cpgLociVariation.cpp')
    if (!file.exists(par$sourceCpp)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceCpp} does not exist!{RET}"))
    Rcpp::sourceCpp(par$sourceCpp)
    
    cat(glue::glue("[{par$prgmTag}]: Loading Source Files form sourceCpp={par$sourceCpp}!{RET}{RET}") )
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Preprocessing:: General Params
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

mergeDirs_vec   <- opt$datDir %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

lociBetaKey_vec <- opt$lociBetaKey %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
lociPvalKey_vec <- opt$lociPvalKey %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
lociPvalMin_vec <- opt$lociPvalMin %>% str_split(pattern=',', simplify=TRUE) %>% as.double() %>% as.vector()

featuresCsv_vec <- opt$featuresCsv %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
featuresDml_vec <- opt$featuresDml %>% str_split(pattern=',', simplify=TRUE) %>% as.integer() %>% as.vector()

seed_vec <- stringr::str_split(opt$seeds, pattern=',', simplify=TRUE) %>% as.vector() %>% as.integer()

class_var <- rlang::sym(opt$classVar)
class_idx <- rlang::sym("Class_Idx")
exp_var   <- 'Experiment_Key'
exp_sym   <- rlang::sym(exp_var)

if (!is.null(opt$platform)) opt$outDir <- file.path(opt$outDir, opt$platform)
if (!is.null(opt$version))  opt$outDir <- file.path(opt$outDir, opt$version)
if (!is.null(opt$classVar)) opt$outDir <- file.path(opt$outDir, opt$classVar)
if (!is.null(opt$workflow)) opt$outDir <- file.path(opt$outDir, opt$workflow)

# TBD:: par$classVar is a temporary hack, needs to be integrated properly...
#  This should be done now and removed...
#  if (!is.null(par$classVar)) opt$outDir <- file.path(opt$outDir, par$classVar)

if (!dir.exists(opt$outDir))
  dir.create(opt$outDir, recursive=TRUE)

cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Preprocessing:: Load Feature Selection
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cgn_pre_str <- NULL
if (!is.null(featuresCsv_vec) && length(featuresCsv_vec)!=0 ) {
  cgn_pre_str <- featuresCsv_vec %>% base::basename() %>% stringr::str_remove('.csv.gz$') %>% paste(collapse='-')
  cgn_pre_tib <- lapply(featuresCsv_vec, loadCgnFeatures, id="Probe_ID", verbose=opt$verbose) %>% 
    dplyr::bind_rows() %>% dplyr::distinct()
  cgn_pre_len <- cgn_pre_tib %>% base::nrow()
  
  cat(glue::glue("[{par$prgmTag}]: Done. Loading Third Party Features; cgn_pre_len={cgn_pre_len}, cgn_pre_str={cgn_pre_str}!{RET}{RET}") )
} else {
  cat(glue::glue("[{par$prgmTag}]: Done. Did NOT load any Thrid Party features!{RET}{RET}") )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Preprocessing:: File Identification
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

for (betaKey in lociBetaKey_vec) {
  for (pvalKey in lociPvalKey_vec) {
    for (pvalMin in lociPvalMin_vec) {
      cat(glue::glue("[{par$prgmTag}]: Starting; betaKey={betaKey}, pvalKey={pvalKey}, pvalMin={pvalMin}.{RET}{RET}") )
      
      #     }
      #     break
      #   }
      #   break
      # }
      
      opt$clean <- TRUE

      betaStr <- betaKey %>% stringr::str_replace_all('_', '-')
      pvalStr <- paste(pvalKey %>% stringr::str_replace_all('_', '-'), pvalMin, sep='-')
      dirName <- paste(betaStr,pvalStr, sep='_')
      outName <- paste(opt$classVar, opt$runName, dirName, sep='_')
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                     Build Current Output Directory::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      cur_opt_dir <- file.path(opt$outDir, dirName)
      if (!dir.exists(cur_opt_dir)) dir.create(cur_opt_dir, recursive=TRUE)
      
      if (opt$verbose>0)
        cat(glue::glue("[{par$prgmTag}]: Built; cur_opt_dir={cur_opt_dir}!{RET}") )
      
      opt$plotDir <- file.path(cur_opt_dir,'plots')
      if (!dir.exists(opt$plotDir)) dir.create(opt$plotDir, recursive=TRUE)
      
      if (opt$clean) unlink(list.files(cur_opt_dir, full.names=TRUE))
      
      cTracker <- timeTracker$new(verbose=opt$verbose)
      
      # Defined Output files::
      beta_masked_rds <- file.path(cur_opt_dir, paste(outName,'beta_masked_mat.rds', sep='.') )
      pval_masked_rds <- file.path(cur_opt_dir, paste(outName,'pval_masked_mat.rds', sep='.') )
      index_masks_csv <- file.path(cur_opt_dir, paste(outName,'beta_masked_idx.csv.gz', sep='.') )
      class_ss_csv <- file.path(cur_opt_dir, paste(outName,'ClasSampleSheet.sorted.csv.gz', sep='.') )
      
      full_dml_csv <- file.path(cur_opt_dir, paste(outName,'full-dml.csv.gz', sep='.') )
      full_dbl_csv <- file.path(cur_opt_dir, paste(outName,'full-dbl.csv.gz', sep='.') )
      
      dml_opt_csv  <- file.path(cur_opt_dir, paste(outName,'program-options.csv', sep='.') )
      dml_par_csv  <- file.path(cur_opt_dir, paste(outName,'program-parameters.csv', sep='.') )
      dml_time_csv <- file.path(cur_opt_dir, paste(outName,'time-tracker.csv.gz', sep='.') )
      
      # opt$clean <- FALSE
      # opt$clean <- TRUE
      opt$verbose <- 40

      # Classes Variable::
      class_str <- opt$trainClass
      sentrix_name <- "Sentrix_Name"
      
      if (FALSE && stringr::str_starts(opt$runName, "NA12878")) {
        # class_str <- "1,2,3,4,5,6"
        # class_vec <- c(1:44)
        # class_str <- class_vec %>% paste(collapse=",")
        # sentrix_name <- "Sentrix_Uniq"
        # tmp_sam_tib <- readr::read_csv("/Users/bretbarnes/Documents/scratch/RStudio/merge_builds/NA12878-granular0-1.5/detect_platform/ind/NA12878-granular0-1.5_ind_AutoSampleSheet.csv.gz")
        
        # tmp_sam_tib <- readr::read_csv("/Users/bretbarnes/Documents/scratch/RStudio/merge_builds/NA12878-granular-3-rand1-2/detect_version/ind/NA12878-granular-3-rand1-2_ind_betas.dat.csv.gz")
        # tmp_sam_tib <- readr::read_csv("/Users/bretbarnes/Documents/scratch/RStudio/merge_builds/NA12878-granular-3-rand1-2/detect_version/ind/NA12878-granular-3-rand1-2_ind_AutoSampleSheet.csv.gz")
        # class_vec <- tmp_sam_tib %>% dplyr::distinct(detect_version) %>% dplyr::arrange(detect_version) %>% dplyr::pull(detect_version)
        
        tmp_sam_tib <- readr::read_csv("/Users/bretbarnes/Documents/scratch/RStudio/merge_builds/NA12878-granular-3-rand1-2/detect_version/ind/NA12878-granular-3-rand1-2_ind_AutoSampleSheet.csv.gz")
        class_vec <- tmp_sam_tib %>% dplyr::distinct(!!class_var) %>% dplyr::arrange(!!class_var) %>% dplyr::pull(!!class_var)
        class_str <- class_vec %>% paste(collapse=",")
        sentrix_name <- "Sentrix_Name"
      }
      class_str <- NULL
      
      beta_file_tib <- getCallsMatrixFiles(
        betaKey=betaKey,pvalKey=pvalKey,pvalMin=pvalMin, 
        dirs=mergeDirs_vec, cgn=NULL, classes=class_str,
        class_var=class_var, class_idx=class_idx, 
        pval_name=opt$samplePvalName, pval_perc=opt$samplePvalPerc,
        clean=opt$clean,addPval=opt$addPval, 
        sentrix_name=sentrix_name,idKey="Probe_ID", 
        betaName='beta', pvalName='pval', del='.', exp_name=exp_sym,
        beta_rds=beta_masked_rds, pval_rds=pval_masked_rds, 
        ss_csv=class_ss_csv, mask_csv=index_masks_csv,
        sam_suffix="_AutoSampleSheet.csv.gz", 
        dat_suffix="_MergedDataFiles.tib.csv.gz", 
        verbose=opt$verbose, vt=3,tc=1,tt=cTracker)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                    Load Raw and Imputed Sorted Matricies::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # Now load previous results from file::
      if (opt$addPval) pval_select_mat <- 
        loadFromFileTib(tib=beta_file_tib, type="Pval", 
                        verbose=opt$verbose,tc=1,tt=pTracker)
      sampleSheet_tib <- 
        loadFromFileTib(tib=beta_file_tib, type="SampleSheet",
                        verbose=opt$verbose,tc=1,tt=pTracker) %>% 
        dplyr::mutate(Experiment_Key=Experiment_Key %>%
                        stringr::str_remove('-Samples') %>% 
                        stringr::str_remove('-EPIC-Core') %>% 
                        stringr::str_remove('CNTL-'),
                      platformUsed=detect_platform, 
                      platVersUsed=detect_version)
      
      # if ("Sample_Per1" %in% base::names(sampleSheet_tib)) {
      #   sampleSheet_tib <- sampleSheet_tib %>%
      #     dplyr::mutate()
      # } else if ("Sample_Per1" %in% base::names(sampleSheet_tib)) {
      #   
      # } else {
      #   
      # }

      # Fix detect_version (class_var) that have decimal (double) values
      #   - if double add a zero? Split on dot
      if (FALSE) {
        sampleSheet_tib %>% 
          dplyr::mutate(DV=as.character(detect_version)) %>% 
          dplyr::select(Sentrix_Name, DV) %>% 
          as.data.frame()
      }

      index_masks_tib <- loadFromFileTib(tib=beta_file_tib, type="Mask")
      beta_impute_mat <- loadFromFileTib(tib=beta_file_tib, type="Beta")
      beta_masked_mat <- beta_impute_mat
      
      # Rebuild NA beta matrix for DML/dBL calculations::
      pval_na_idx_vec <- index_masks_tib %>% dplyr::arrange(idx) %>% dplyr::distinct(idx) %>% dplyr::pull(idx) %>% as.vector()
      if (!is.null(pval_na_idx_vec) && length(pval_na_idx_vec)!=0) beta_masked_mat[ pval_na_idx_vec ] <- NA
      
      labs_idx_vec <- sampleSheet_tib %>% dplyr::pull(!!class_idx) %>% as.vector()
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                   Plot Pairwise R-squared and DeltaBeta::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      #
      # TBD::
      #  3. Add sample specific opt$samplePvalPerc
      #     - opt$filterSamples <- TRUE
      #     - sampleSheet_tib %>% dplyr::group_by(Experiment_Key) %>% dplyr::summarise(Exp_Pass_avg=mean(!!sam_pval_name_sym, na.rm=TRUE), Exp_Pass_max=max(!!sam_pval_name_sym, na.rm=TRUE))
      #
      #  4. Load signal sets
      #  5. Plot signal sets
      #
      if (opt$plot_pairs) {
        sam_pval_name_sym <- rlang::sym(opt$samplePvalName)
        
        opt$filterSamples <- FALSE
        opt$filterSamples <- TRUE
        # Load Manifest base on Sample Sheet
        #   detect_platform = platformUsed
        
        man_path_csv <- NULL
        man_platform <- sampleSheet_tib %>% dplyr::distinct(platformUsed) %>% head(n=1) %>% dplyr::pull(platformUsed)
        man_version  <- sampleSheet_tib %>% dplyr::distinct(platVersUsed) %>% head(n=1) %>% dplyr::pull(platVersUsed)
        man_pattern  <- paste(paste(man_platform,man_version, sep='-'), 'manifest.sesame-base.cpg-sorted.csv.gz', sep='.')
        man_path_csv <- list.files(file.path(par$datDir, 'manifest/base'), man_pattern, full.names=TRUE) %>% head(n=1)
        if (length(man_path_csv)==0) {
          # opt$manDirPath <- file.path(par$topDir, "data/manifests/methylation/Chicago-Ober-Custom")
          man_path_csv <- list.files(file.path(opt$manDirPath), man_pattern, full.names=TRUE) %>% head(n=1)
          
          if (length(man_path_csv)==0) {
            man_path_csv <- list.files(file.path(par$datDir, 'manifest/base'), pattern = "EPIC", full.names=TRUE) %>% head(n=1)
          }
        }
        ses_man_tib  <- suppressMessages(suppressWarnings( readr::read_csv(man_path_csv) ))
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                   Update Sample Sheet for Plotting::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        # Below are two attempts to set up the Experiment_Key which should be defined by performing each 
        #  merge_builds individually::
        #
        # Update class_var by build_source
        # class_two_str <- 'build_source'
        # class_two_sym <- rlang::sym(class_two_str)
        # sampleSheet_tib <- sampleSheet_tib %>% dplyr::mutate(Experiment_Key=!!class_two_sym) # %>% dplyr::select(Experiment_Key)
        #
        # if (length(grep(class_two_str, names(sampleSheet_tib))) != 0) {
        #   sampleSheet_tib <- sampleSheet_tib %>% 
        #     dplyr::mutate(!!class_var := paste(!!class_two_sym,!!class_var, sep='_')) # %>% dplyr::select(!!class_var) %>% as.data.frame()
        # }
        
        plotSheet_tib <- sampleSheet_tib %>% 
          dplyr::filter(!!sam_pval_name_sym >= opt$samplePvalPerc) %>%
          dplyr::rename(Class_Name=!!class_var) %>% 
          dplyr::group_by(!!exp_sym) %>%
          dplyr::mutate(Class_Int=dplyr::cur_group_id() ) %>% 
          dplyr::arrange(Class_Int) %>% 
          dplyr::mutate(Rank_Chr=rawToChar(as.raw(64+Class_Int[1])) ) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(!!exp_sym,Class_Name) %>%
          dplyr::mutate(Rank_Idx=dplyr::row_number()) %>%
          dplyr::mutate(Plot_Name=paste0(Rank_Chr,'_Rep',Rank_Idx,'_',as.integer(!!sam_pval_name_sym))) %>%
          dplyr::ungroup() %>%
          dplyr::select(Sentrix_Name, Class_Name, !!class_idx, !!exp_sym, Class_Int, !!sam_pval_name_sym, Rank_Idx,Rank_Chr,Plot_Name, everything())
        
        # NOTE:: Previous Versions Re-assigned Sample_Name to equal Class_Name
        # plotSheet_tib <- plotSheet_tib %>% dplyr::mutate(Sample_Name=Class_Name)
        #  plotSheet_tib %>% dplyr::select(Sample_Name, Class_Name)
        # If there is only one sample type (i.e. mix of normal people then set 
        #  Sample_Name to 'Unknown')
        if (grep("Sample_Name", names(plotSheet_tib)) %>% length() == 0) {
          # plotSheet_tib <- plotSheet_tib %>% dplyr::mutate(Sample_Name=Class_Name)
          plotSheet_tib <- plotSheet_tib %>% dplyr::mutate(Sample_Name="Unknown")
        }
        
        # NOTE:: Once we have both Sample_Name and Class_Name we should split
        #  by Sample_Name and then later by Class_Name... Testing....
        #
        # ss_class_list <- plotSheet_tib %>% split(.$Amp_Incubation)
        # ss_class_list <- plotSheet_tib %>% split(.$Class_Name)
        ss_class_list <- plotSheet_tib %>% split(.$Sample_Name)
        ss_class_keys <- names(ss_class_list)
        
        for (cIdx in c(1:length(ss_class_keys))) {
          class_key <- ss_class_keys[[cIdx]]
          ss_plot_tib <- ss_class_list[[class_key]]

          r2_raw_pdf  <- file.path(opt$plotDir, paste('r-squared.all-methods.raw',class_key,'pdf', sep='.'))
          r2_norm_pdf <- file.path(opt$plotDir, paste('r-squared.all-methods.normalized',class_key,'pdf', sep='.'))
          
          chip_vec <- NULL
          if (sentrix_name=="Sentrix_Name") chip_vec <- ss_plot_tib %>% dplyr::pull(Sentrix_Name)
          if (sentrix_name=="Sentrix_Uniq") chip_vec <- ss_plot_tib %>% dplyr::pull(Sentrix_Uniq)
          name_vec <- ss_plot_tib %>% dplyr::pull(Plot_Name)
          
          samp_beta_mat <- beta_masked_mat[, chip_vec]
          colnames(samp_beta_mat) <- name_vec
          
          cor_mat <- cor(samp_beta_mat, use="pairwise.complete.obs", method="pearson")
          cor_min <- min(cor_mat)
          cor_max <- max(cor_mat)
          cor_len <- cor_max - cor_min
          cor_min_rnd <- round(cor_min, 3)
          
          # Raw
          M.raw <- cor_mat
          P.raw <- cor.mtest(M.raw)
          pdf(r2_raw_pdf)
          title.raw <- glue::glue("{RET}{class_key}: R-Squared (min-r2={cor_min_rnd}, sig={pvalMin}) Raw")
          corrplot(M.raw, type="upper", order="hclust", p.mat = P.raw, 
                   sig.level = pvalMin, insig = "blank", title=title.raw,
                   tl.srt=45, tl.cex = 0.5, addrect = 3)
          dev.off()
          
          # Normalized
          M.norm <- (cor_mat - cor_min) / cor_len
          P.norm <- cor.mtest(M.norm)
          pdf(r2_norm_pdf)
          title.norm <- glue::glue("{RET}{class_key}: R-Squared (min-r2={cor_min_rnd},sig={pvalMin}) Normalized")
          corrplot(M.norm, type="upper", order="hclust", p.mat = P.norm, 
                   sig.level = pvalMin, insig = "blank", title=title.norm,
                   tl.srt=45, tl.cex = 0.5, addrect = 3)
          dev.off()
          
          #
          # Run Plot Pairs for each Experiment pairwise-combination::
          #
          # ss_exp_list <- ss_class_list[[class_key]] %>% split(.$Experiment_Key)
          
          ss_exp_list <- ss_class_list[[class_key]] %>% split(.$Class_Name)
          ss_exp_keys <- names(ss_exp_list)
          
          for (idxA in c(1:length(ss_exp_keys))) {
            for (idxB in c(1:length(ss_exp_keys))) {
              if (idxA<idxB) {
                exp_keyA <-ss_exp_keys[idxA]
                exp_keyB <-ss_exp_keys[idxB]
                
                cat(glue::glue("[{par$prgmTag}]:{TAB} Ploting; ",
                               "Class={class_key}: A({idxA})={exp_keyA} vs. ",
                               "B({idxB})={exp_keyB}!{RET}") )
                
                cur_ss_tibA <- ss_exp_list[[exp_keyA]] %>% 
                  dplyr::arrange(desc(!!sam_pval_name_sym))
                cur_ss_tibB <- ss_exp_list[[exp_keyB]] %>% 
                  dplyr::arrange(desc(!!sam_pval_name_sym))
                if (opt$filterSamples) {
                  cur_ss_tibA <- cur_ss_tibA %>% 
                    dplyr::filter(!!sam_pval_name_sym >= opt$samplePvalPerc)
                  cur_ss_tibB <- cur_ss_tibB %>% 
                    dplyr::filter(!!sam_pval_name_sym >= opt$samplePvalPerc)
                }
                
                # Reduce to 3 max for plotting::
                red_ss_tibA <- reduceSortedTib(cur_ss_tibA)
                red_ss_tibB <- reduceSortedTib(cur_ss_tibB) %>% 
                  dplyr::mutate(Plot_Name=stringr::str_replace(Plot_Name,'^A_', 'B_'))
                
                if (sentrix_name=="Sentrix_Uniq") {
                  beta_tibA <- beta_masked_mat[ ,red_ss_tibA$Sentrix_Uniq ] %>% 
                    as.data.frame() %>% purrr::set_names(red_ss_tibA$Plot_Name) %>% 
                    tibble::rownames_to_column(var="Probe_ID") %>% tibble::as_tibble()
                  pval_tibA <- pval_select_mat[ ,red_ss_tibA$Sentrix_Uniq ] %>% 
                    as.data.frame() %>% purrr::set_names(red_ss_tibA$Plot_Name) %>% 
                    tibble::rownames_to_column(var="Probe_ID") %>% tibble::as_tibble()
                  
                  beta_tibB <- beta_masked_mat[ ,red_ss_tibB$Sentrix_Uniq ] %>% 
                    as.data.frame() %>% purrr::set_names(red_ss_tibB$Plot_Name) %>% 
                    tibble::rownames_to_column(var="Probe_ID") %>% tibble::as_tibble()
                  pval_tibB <- pval_select_mat[ ,red_ss_tibB$Sentrix_Uniq ] %>% 
                    as.data.frame() %>% purrr::set_names(red_ss_tibB$Plot_Name) %>% 
                    tibble::rownames_to_column(var="Probe_ID") %>% tibble::as_tibble()                  
                } else {
                  beta_tibA <- beta_masked_mat[ ,red_ss_tibA$Sentrix_Name ] %>% 
                    as.data.frame() %>% purrr::set_names(red_ss_tibA$Plot_Name) %>% 
                    tibble::rownames_to_column(var="Probe_ID") %>% tibble::as_tibble()
                  pval_tibA <- pval_select_mat[ ,red_ss_tibA$Sentrix_Name ] %>% 
                    as.data.frame() %>% purrr::set_names(red_ss_tibA$Plot_Name) %>% 
                    tibble::rownames_to_column(var="Probe_ID") %>% tibble::as_tibble()
                  
                  beta_tibB <- beta_masked_mat[ ,red_ss_tibB$Sentrix_Name ] %>% 
                    as.data.frame() %>% purrr::set_names(red_ss_tibB$Plot_Name) %>% 
                    tibble::rownames_to_column(var="Probe_ID") %>% tibble::as_tibble()
                  pval_tibB <- pval_select_mat[ ,red_ss_tibB$Sentrix_Name ] %>% 
                    as.data.frame() %>% purrr::set_names(red_ss_tibB$Plot_Name) %>% 
                    tibble::rownames_to_column(var="Probe_ID") %>% tibble::as_tibble()
                }
                
                if (par$joinType=="left") {
                  beta_plot_tib <- ses_man_tib %>% 
                    dplyr::select(Probe_ID,Probe_Type,DESIGN) %>% 
                    dplyr::rename(Design_Type=DESIGN) %>%
                    dplyr::left_join(beta_tibA, by="Probe_ID") %>% 
                    dplyr::left_join(beta_tibB, by="Probe_ID")
                  
                  pval_plot_tib <- ses_man_tib %>% 
                    dplyr::select(Probe_ID,Probe_Type,DESIGN) %>% 
                    dplyr::rename(Design_Type=DESIGN) %>%
                    dplyr::left_join(pval_tibA, by="Probe_ID") %>% 
                    dplyr::left_join(pval_tibB, by="Probe_ID")
                } else if (par$joinType=="inner") {
                  beta_plot_tib <- ses_man_tib %>% 
                    dplyr::select(Probe_ID,Probe_Type,DESIGN) %>% 
                    dplyr::rename(Design_Type=DESIGN) %>%
                    dplyr::inner_join(beta_tibA, by="Probe_ID") %>% 
                    dplyr::inner_join(beta_tibB, by="Probe_ID")
                  
                  pval_plot_tib <- ses_man_tib %>% 
                    dplyr::select(Probe_ID,Probe_Type,DESIGN) %>% 
                    dplyr::rename(Design_Type=DESIGN) %>%
                    dplyr::inner_join(pval_tibA, by="Probe_ID") %>% 
                    dplyr::inner_join(pval_tibB, by="Probe_ID")
                  
                } else {
                  cat(glue::glue("[{par$prgmTag}]: ERROR: Unsupported joinType={joinType}!!! Skipping...{RET}{RET}") )
                  next
                }
                
                gg <- plotPairsBeta(
                  beta=beta_plot_tib, pval=pval_plot_tib, 
                  sample=class_key, nameA=exp_keyA, nameB=exp_keyB,
                  outDir=opt$plotDir,
                  probeType='cg', field='Beta', 
                  field_str=betaKey, detp=pvalKey, minPval=pvalMin,
                  format='both', verbose=opt$verbose)
                
                cat(glue::glue("[{par$prgmTag}]:{TAB} Done.{RET}{RET}") )
                
                if (opt$single) break
              }
            }
            if (opt$single) break
          }
          if (opt$single) break
        }
      }
      
      if (!opt$buildDbl && !opt$buildDml) next
      
      # break
      
      # To Do:: MVP
      #
      #  - Plot R-squared in diagnal plot by Sample and by Experiment for Replicates
      #
      #  - Analytical Screen each Experiment with dBL
      #    - Titration
      #    - Replicate
      #
      #  - Summarise AutoSampleSheets:: 4 CTL experiments and DKFZ
      #    - Percent Passing OOB
      #    - Percent Passing NEG
      #    - Average Intensity
      #    - GCT Scores
      #    - Bisulfite I/II Intensities
      #
      # DKFZ = '/Users/bbarnes/Documents/Projects/methylation/VA_MVP/analysis/DKFZ/swifthoof_DKFZ.AutoSampleSheet.csv.gz'
      #
      # SLIDES::
      #
      #   1. Poob vs. Negs
      #      - Previous Poob(0.1, >90%) = >99% for both providers
      #      - DKFZ: Only used Negs, what's their passing rate with Poob?
      #   2. Control Experiments
      #   3. Screening Experiments
      #
      
      # - DKFZ way worse performance vs. Service Provider Controls and Real MVP Samples by both Negs/Poob
      #   - But they have made liget clinical grade actions
      # - Service Provide Control and Internal Control samples have high R-Squared and Delta-Beta Values
      #   Two type of plots: x2 thresholds
      #
      # - Screening Results:: ???
      # 
      # - Reaching out to Horvath, LEGX and UCD regarding docker image
      #
      # - What about beta values?
      #
      opt$isDKFZ <- FALSE
      if (opt$isDKFZ) {
        dkfz_ss_csv <- '/Users/bbarnes/Documents/Projects/methylation/VA_MVP/analysis/DKFZ/swifthoof_DKFZ.AutoSampleSheet.csv.gz'
        dkfz_ss_tib <- readr::read_csv(dkfz_ss_csv)
        
        dkfz_ss_tib %>% dplyr::group_by(platformUsed,platVersUsed,Chip_Format,Bead_Pool) %>% 
          dplyr::summarise(Total_Count=n(), 
                           Pass_Negs_Cnt=count(Negs_Pass_0_Perc>=98, na.rm=TRUE),
                           Pass_Poob_Cnt=count(Poob_Pass_0_Perc>=90, na.rm=TRUE),
                           Pass_Negs_Per=round(100*Pass_Negs_Cnt/Total_Count,2),
                           Pass_Poob_Per=round(100*Pass_Poob_Cnt/Total_Count,2) )
        
        sampleSheet_tib %>% dplyr::group_by(platformUsed,platVersUsed,Chip_Format,Bead_Pool,!!exp_sym) %>% 
          dplyr::summarise(Total_Count=n(), 
                           Pass_Negs_Cnt=count(Negs_Pass_0_Perc>=98, na.rm=TRUE),
                           Pass_Poob_Cnt=count(Poob_Pass_0_Perc>=90, na.rm=TRUE),
                           Pass_Negs_Per=round(100*Pass_Negs_Cnt/Total_Count,2),
                           Pass_Poob_Per=round(100*Pass_Poob_Cnt/Total_Count,2) )
      }
      
      #
      # DKFZ has worse performance, but amazing actionablity
      #
      
      # platformUsed platVersUsed Chip_Format Bead_Pool Total_Count Pass_Negs_Cnt Pass_Poob_Cnt Pass_Negs_Per Pass_Poob_Per
      # <chr>        <chr>        <chr>       <chr>           <int>         <int>         <int>         <dbl>         <dbl>
      # 1 HM450        B2           12x1        BP123            2409          2367          2223          98.3          92.3
      # > 
      
      #
      # CNTL_VendB_10092020 has two non-processed samples (need to investigate)
      #
      
      # platformUsed platVersUsed Chip_Format Bead_Pool Experiment_Key      Total_Count Pass_Negs_Cnt Pass_Poob_Cnt Pass_Negs_Per Pass_Poob_Per
      # <chr>        <chr>        <chr>       <chr>     <chr>                     <int>         <int>         <int>         <dbl>         <dbl>
      # 1 EPIC         B4           8x1         EPIC      BETA-8x1                     16            16            16           100         100  
      # 2 EPIC         B4           8x1         EPIC      CNTL_VendA_10092020          24            24            22           100          91.7
      # 3 EPIC         B4           8x1         EPIC      CNTL_VendB_10092020          22            22            22           100         100  
      # 4 EPIC         B4           8x1         EPIC      DELTA-8x1                    32            32            32           100         100  
      
      
      # Scratch for MVP::
      if (FALSE) {
        sampleSheet_tib %>% dplyr::select(Sentrix_Name,!!class_var,!!class_idx)
        
        mvpA_vec <- c('203962710025', '203962710079', '203962710081')
        mvpB_vec <- c('204229180144', '204229190022', '204229190023')
        
        samp_scr_tib <- lapply(list.files(mergeDirs_vec, pattern = 'SampleSheet.csv.gz', recursive = TRUE, full.names = TRUE), readr::read_csv) %>% 
          dplyr::bind_rows() %>% dplyr::select(Sentrix_Name, Poob_Pass_0_Perc)
        
        mvp_ss_tib <- sampleSheet_tib %>% tidyr::separate(Sentrix_Name, into=c('Sentrix_ID','Sentrix_Pos'), sep='_', remove=FALSE) %>% 
          dplyr::inner_join(samp_scr_tib, by="Sentrix_Name") %>% dplyr::filter(Poob_Pass_0_Perc>opt$samplePvalPerc)
        
        sampleSheetA_tib <- mvp_ss_tib %>% dplyr::filter(Sentrix_ID %in% mvpA_vec) %>% dplyr::mutate(Vender='A')
        sampleSheetB_tib <- mvp_ss_tib %>% dplyr::filter(Sentrix_ID %in% mvpB_vec) %>% dplyr::mutate(Vender='B')
        
        sampleSheetA_sorted_tib <- sampleSheetA_tib %>% dplyr::filter(Class_Idx==0) %>% dplyr::arrange(-Poob_Pass_0_Perc)
        sampleSheetB_sorted_tib <- sampleSheetB_tib %>% dplyr::filter(Class_Idx==0) %>% dplyr::arrange(-Poob_Pass_0_Perc)
        
        pickedA_sentrixAll <- sampleSheetA_sorted_tib %>% dplyr::pull(Sentrix_Name)
        pickedB_sentrixAll <- sampleSheetB_sorted_tib %>% dplyr::pull(Sentrix_Name)
        
        pickedA_sentrix <- pickedA_sentrixAll[c(1,as.integer(length(pickedA_sentrixAll)/2),length(pickedA_sentrixAll))]
        pickedB_sentrix <- pickedB_sentrixAll[c(1,as.integer(length(pickedB_sentrixAll)/2),length(pickedB_sentrixAll))]
        
        sampleSheet_tib <- dplyr::bind_rows(
          dplyr::filter(sampleSheet_tib, Sentrix_Name %in% pickedA_sentrix) %>% 
            dplyr::mutate(AutoSample_dB_Key=paste('A',AutoSample_dB_Key,sep='_'), Class_Idx=Class_Idx+0),
          
          dplyr::filter(sampleSheet_tib, Sentrix_Name %in% pickedB_sentrix) %>% 
            dplyr::mutate(AutoSample_dB_Key=paste('B',AutoSample_dB_Key,sep='_'), Class_Idx=Class_Idx+1)
        ) %>% purrr::set_names(c('Sentrix_Name','Sample_Name','Class_Idx'))
        
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                         Plot R-Squared/DeltaBeta::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # Example for plotting::
      if (FALSE) {
        # Load Auto Sample Sheet from merged builds::
        #
        exp_src_var <- "Exp_Source_Name"
        exp_src_var <- rlang::sym(exp_src_var)
        auto_ss_tib <- NULL
        for (merDir in mergeDirs_vec) {
          exp_source <- base::basename(merDir)
          cur_ss_csv <- list.files(merDir, pattern='_AutoSampleSheet.csv.gz', full.names=TRUE) %>% head(n=1)
          cur_ss_tib <- suppressMessages(suppressWarnings( readr::read_csv(cur_ss_csv) )) %>% 
            dplyr::mutate(!!exp_src_var:=exp_source)
          auto_ss_tib <- auto_ss_tib %>% dplyr::bind_rows(cur_ss_tib)
        }
        auto_ss_tib %>% dplyr::select(Sentrix_Name,!!class_var,!!exp_src_var) %>% 
          dplyr::arrange(!!class_var,!!exp_src_var) %>%
          dplyr::group_by(!!class_var,!!exp_src_var) %>% 
          dplyr::mutate(Class_Idx=dplyr::row_number())
        
        # sampleSheet_tib %>% dplyr::mutate(Class_Int=dplyr::cur_group_id(), !!class_idx:=Class_Int-1 )
        
        # WE really want::
        #   Exp_Idx
        #   Sample_Idx
        #   Class_Idx  = What is going to be tested
        auto_ss_tib %>% dplyr::select(Sentrix_Name,!!class_var,!!exp_src_var) %>% 
          dplyr::arrange(!!class_var,!!exp_src_var) %>%
          dplyr::group_by(!!class_var,!!exp_src_var) %>% 
          dplyr::mutate(Class_Int=dplyr::cur_group_id(), !!class_idx:=Class_Int-1 )
        
        # Load Manifest::
        #
        # par$platform
        
        # Load Detection P-values::
        
        # Only pick two classes::
        picked_sentrix <- sampleSheet_tib %>% 
          dplyr::filter(Class_Idx==0 | Class_Idx==1) %>% 
          dplyr::pull(Sentrix_Name)
        
        useSampleName <- FALSE
        if (useSampleName) {
          picked_sam_vec <- sampleSheet_tib %>% 
            dplyr::filter(Class_Idx==0 | Class_Idx==1) %>% 
            dplyr::group_by(Sample_Name) %>%
            dplyr::mutate(Rep_Num=row_number(), 
                          Rep_Name=paste0(Sample_Name,'_Rep',Rep_Num)) %>% 
            dplyr::pull(Rep_Name) %>% paste('Beta', sep='.')
        } else {
          # NOTE:: Using !!class_var now instead of "Sample_Class"
          #
          picked_sam_vec <- sampleSheet_tib %>% 
            dplyr::filter(Class_Idx==0 | Class_Idx==1) %>% 
            # dplyr::group_by(Sample_Class) %>% 
            dplyr::group_by(!!class_var) %>% 
            dplyr::mutate(Rep_Num=row_number(), 
                          Rep_Name=paste0(!!class_var,'_Rep',Rep_Num)) %>% 
            #              Rep_Name=paste0(Sample_Class,'_Rep',Rep_Num)) %>% 
            dplyr::pull(Rep_Name) %>% paste('Beta', sep='.')
        }
        
        beta_masked_tib <- beta_masked_mat[,picked_sentrix] %>% 
          as.data.frame() %>% 
          purrr::set_names(picked_sam_vec) %>% 
          tibble::rownames_to_column(var="Probe_ID") %>% 
          tibble::as_tibble() %>%
          dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2), 
                        Design_Type=stringr::str_remove(Probe_ID, '^.*_[A-Z][A-Z]') %>% stringr::str_remove('[A-Z][0-9]$'),
                        Design_Type=dplyr::case_when(
                          Design_Type=='1' ~ 'I', Design_Type=='2' ~ 'II', TRUE ~ NA_character_
                        )) %>%
          dplyr::select(Probe_ID,Probe_Type,Design_Type, everything())
        
        # Quick Fix...
        #  TBD:: Use manifest that was loaded above...
        #
        beta_masked_tib$Design_Type <- 'I'
        
        plotDir <- file.path(opt$outDir,'plots')
        if (useSampleName) {
          gg <- plotPairsBeta(beta_masked_tib, sample='T00vs50', 
                              nameA='T00DZ', nameB='T50DZ', outDir=plotDir,
                              probeType='cg', field='Beta', field_str='ind_beta', 
                              detp='i_poob', minPval=pvalMin,
                              format='pdf',
                              verbose=opt$verbose)
        } else {
          beta_masked_tib2 <- beta_masked_tib %>% 
            dplyr::mutate(
              Design_Type=stringr::str_remove(Probe_ID,'^[^_]*_[A-Z][A-Z]') %>% stringr::str_remove('[0-9]$'),
              Design_Type=dplyr::case_when(
                Design_Type=='1' ~ 'I',
                Design_Type=='2' ~ 'II',
                TRUE ~ NA_character_)
            )
          
          gg <- plotPairsBeta(beta_masked_tib2, sample='NEGvsPOS', 
                              nameA='nSARSCov2', nameB='pSARSCov2', outDir=plotDir,
                              probeType='cg', field='Beta', field_str=betaKey, 
                              detp=pvalKey, minPval=pvalMin,
                              format='pdf',
                              verbose=opt$verbose)
        }
        
        # gg <- plotPairsBeta(beta_masked_tib, sample='T00vs50', nameA='T00DZ', nameB='T50DZ', outDir=plotDir,
        #                 probeType='cg', field='Beta', field_str='ind_beta', detp='i_poob', # maxCnt = , minPval=minPval,
        #                 spread=spread, outType=outType, dpi=dpi, format=format,
        #                 verbose=verbose, tc=tc+1)
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                               Build DMLs::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      full_dml_tib <- NULL
      dml_beg_txt  <- paste(full_dml_csv,'begTime.txt', sep='.')
      dml_end_txt  <- paste(full_dml_csv,'endTime.txt', sep='.')
      
      if (opt$buildDml) {
        if (!opt$clean && file.exists(full_dml_csv) &&
            file.exists(dml_beg_txt) && file.exists(dml_end_txt) &&
            file.mtime(dml_beg_txt)  < file.mtime(full_dml_csv) && 
            file.mtime(full_dml_csv) < file.mtime(dml_end_txt) ) {
          
          cat(glue::glue("[{par$prgmTag}]:{TAB} DMLs already exists! Reading; full_dml_csv={full_dml_csv}...{RET}") )
          full_dml_tib <- suppressMessages(suppressWarnings( readr::read_csv(full_dml_csv) ))
          cat(glue::glue("[{par$prgmTag}]:{TAB} Done. Reading Full DMLs{RET}{RET}") )
          
        } else {
          cat(glue::glue("[{par$prgmTag}]:{TAB} Calculating dmls...{RET}") )
          
          if (file.exists(dml_beg_txt))  unlink(dml_beg_txt)
          if (file.exists(full_dml_csv)) unlink(full_dml_csv)
          if (file.exists(dml_end_txt))  unlink(dml_end_txt)
          
          samp_dml_tib <- sampleSheet_tib %>% dplyr::rename(Class_Var = !!class_var)
          full_dml_tib <- dmlsToTib( DML_Local(betas=beta_masked_mat, sample.data=samp_dml_tib, formula = ~Class_Var ), 
                                     verbose=opt$verbose, vt=1, tc=2, tt=cTracker) %>% 
            dplyr::mutate_if( is.double, list(round), opt$percisionPval )
          
          cat(glue::glue("[{par$prgmTag}]:{TAB} Writing Full DML(CSV)={full_dml_csv}...{RET}") )
          system(glue::glue("touch {dml_beg_txt}"))
          Sys.sleep(1)
          readr::write_csv(full_dml_tib,full_dml_csv)
          Sys.sleep(1)
          system(glue::glue("touch {dml_end_txt}"))
          cat(glue::glue("[{par$prgmTag}]:{TAB} Done. Writing Full DMLs{RET}{RET}") )
        }
        rank_dml_tib <- full_dml_tib %>% dplyr::group_by(Probe_ID) %>% 
          dplyr::summarise(Rank_Avg=mean(Rank, na.rm=TRUE)) %>% dplyr::arrange(Rank_Avg)
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                               Build dBLs::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if (opt$buildDbl) {
        # 
        # Suggested Older Install to fix the problem::
        #   packageurl <- "https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.9.900.3.0.tar.gz"
        #   install.packages(packageurl, repos=NULL, type="source")
        #
        # library(Rcpp)
        # library(RcppArmadillo)
        # There is a special Rccp source method...
        # BAD:: source("/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/scripts/R/Rcpp/cpgLociVariation.cpp")
        # GOOD: Rcpp::sourceCpp(par$sourceCpp)
        #       Rcpp::sourceCpp("/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/scripts/R/Rcpp/cpgLociVariation.cpp")
        #
        
        dbl_beg_txt  <- paste(full_dbl_csv,'begTime.txt', sep='.')
        dbl_end_txt  <- paste(full_dbl_csv,'endTime.txt', sep='.')
        
        if (!opt$clean && file.exists(full_dbl_csv) &&
            file.exists(dbl_beg_txt) && file.exists(dbl_end_txt) &&
            file.mtime(dbl_beg_txt) < file.mtime(full_dbl_csv) && 
            file.mtime(dbl_end_txt) > file.mtime(full_dbl_csv) ) {
          
          cat(glue::glue("[{par$prgmTag}]: deltaBetas already existst full_dbl_csv={full_dbl_csv}. Skipping...{RET}") )
          
        } else {
          cat(glue::glue("[{par$prgmTag}]: Calculating deltaBetas...{RET}") )
          cpp.verbose <- 0
          
          #
          # Quick Fix for MVP::
          #
          if (FALSE) {
            mvp_ss_exp_list <- plotSheet_tib %>% dplyr::filter(stringr::str_starts(Class_Name, 'T') ) %>% split(.$Experiment_Key)
            mvp_ss_exp_keys <- names(mvp_ss_exp_list)
            
            for (expKey in mvp_ss_exp_keys) {
              
              sample_cnt_tib <- mvp_ss_exp_list[[expKey]] %>% dplyr::group_by(Class_Name) %>% 
                dplyr::summarise(Count=n()) %>% tibble::as_tibble() %>% dplyr::mutate(Class=as.character(Class_Name) )
              
              sample_ids_vec <- mvp_ss_exp_list[[expKey]]$Sentrix_Name
              sample_key_vec <- sample_cnt_tib %>% dplyr::pull(Class)
              sample_cnt_vec <- sample_cnt_tib %>% dplyr::pull(Count)
              
              beta_tit_mat <- beta_masked_mat[,sample_ids_vec]
              
              full_dbl_tib <- C_crossSampleLociRSquared(beta_tit_mat, sample_cnt_vec, sample_key_vec, cmb=FALSE, verbose=cpp.verbose) %>% 
                as.data.frame() %>% tibble::rownames_to_column(var='Probe_ID') %>% tibble::as_tibble() %>% 
                dplyr::mutate_if(is.double, round, opt$percisionPval)
            }
          }
          
          if (TRUE) {
            # Name by Real Classes::
            #
            sample_cnt_tib <- sampleSheet_tib %>% 
              dplyr::group_by(!!class_var) %>% 
              dplyr::summarise(Count=n()) %>% tibble::as_tibble() %>% dplyr::mutate(Class=as.character(!!class_var) )
            sample_key_vec <- sample_cnt_tib %>% dplyr::pull(Class)
            sample_cnt_vec <- sample_cnt_tib %>% dplyr::pull(Count)
            
            sample_key_vec <- sampleSheet_tib %>% dplyr::pull(!!class_var) %>% unique()
            sample_cnt_vec <- sampleSheet_tib %>% dplyr::group_by(!!class_var) %>% dplyr::summarise(Count=n()) %>% dplyr::pull(Count)
            
            cpp.verbose <- 0
            full_dbl_tib <- C_crossSampleLociRSquared(beta_masked_mat, sample_cnt_vec, sample_key_vec, cmb=TRUE, verbose=cpp.verbose) %>% 
              as.data.frame() %>% tibble::rownames_to_column(var='Probe_ID') %>% tibble::as_tibble() %>% 
              dplyr::mutate_if(is.double, round, opt$percisionPval)
            
            #
            # TBD:: Need to use beta/pval percision on output!!!
            #
            
            system(glue::glue("touch {dbl_beg_txt}"))
            readr::write_csv(full_dbl_tib,full_dbl_csv)
            system(glue::glue("touch {dbl_end_txt}"))
            Sys.sleep(1)
          }
          
          cat(glue::glue("[{par$prgmTag}]: Done. Writing deltaBetas...{RET}{RET}") )
        }
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                       Write Current Params/Options::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      readr::write_csv(opt_tib, dml_opt_csv)
      readr::write_csv(par_tib, dml_par_csv)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                           Build Feature Sets::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if (opt$buildModels) {
        
        # TBD::
        #  - [DONE]: Loop over top feature sizes from DML AND create tag
        #  - [DONE]: Loop over pre-defined features and extract them from DML AND create tag
        #  - [DONE]: Condense top/pre features and sort by Probe_ID
        
        gen_opt_tib <- opt_tib %>% dplyr::filter(
          Option %in% c('lociBetaKey','lociPvalKey','lociPvalMin','samplePvalName', 
                        'samplePvalPerc','platform' ,'version' ,'trainClass','runName')) %>%
          dplyr::mutate(Type="User")
        
        sel_cgn_tib <- tibble::tibble()
        if (!is.null(featuresCsv_vec) && length(featuresCsv_vec)!=0)
          sel_cgn_tib <- full_dml_tib %>% dplyr::filter(Probe_ID %in% cgn_pre_tib$Probe_ID) %>% 
          dplyr::distinct(Probe_ID) %>%
          dplyr::arrange(Probe_ID)
        
        for (dmlSize in featuresDml_vec) {
          
          # TBD:: Add flags for Top/Pre selection source...
          dml_top_str <- paste('dml',dmlSize, sep='-')
          dml_top_tib <- rank_dml_tib %>% head(n=dmlSize) %>% dplyr::bind_rows()
          # dml_top_tib <- full_dml_tib %>% head(n=dmlSize) %>% dplyr::bind_rows()
          
          # Update Current Output Name and Current Output Directory::
          dml_cur_str <- dml_top_str
          if (!is.null(cgn_pre_str) && length(cgn_pre_str)>0)
            dml_cur_str <- paste(dml_top_str,cgn_pre_str, sep='-')
          curName <- paste(outName,dml_cur_str, sep='_')
          
          cur_dml_dir <- file.path(cur_opt_dir, dml_cur_str)
          if (!dir.exists(cur_dml_dir)) dir.create(cur_dml_dir, recursive=TRUE)
          
          # Combine CGN Features::
          cgn_cur_tib <- dplyr::bind_rows(dml_top_tib, sel_cgn_tib) %>% 
            dplyr::distinct(Probe_ID, .keep_all=TRUE) %>% 
            dplyr::arrange(Probe_ID)
          
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          #                Sub-set Beta Matrix by Current CGN Features::
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          
          dml_opt_tib <- gen_opt_tib %>% 
            tibble::add_row( Option="featureSizeDml", Value=as.character(dmlSize), Type="Feature" )
          
          if (!is.null(cgn_pre_str))
            dml_opt_tib <- dml_opt_tib %>% tibble::add_row( Option="featureNamePre", Value=cgn_pre_str, Type="Feature" )
          
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          #                               Build Models::
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          
          beta_impute_matT <- beta_impute_mat %>% t()
          
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          #                       Run Directly or Launch Cluster::
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          
          for (seed_val in seed_vec) {
            cat(glue::glue("[{par$prgmTag}]:{TAB} Starting; seed_val={seed_val}...{RET}") )
            
            seed_str  <- paste('seed',seed_val, sep='-')
            seed_dir  <- file.path(cur_dml_dir, seed_str)
            seed_name <- paste(curName,seed_str, sep='_')
            if (!dir.exists(seed_dir)) dir.create(seed_dir, recursive=TRUE)
            
            if (opt$cluster) {
              
              run_sh  <- file.path(seed_dir, paste(seed_name,'sh', sep='.'))
              add_tib <- tibble::tibble(Option=c("lociBetaKey","lociPvalKey","lociPvalMin","featuresDml"),
                                        Value=c(betaKey,pvalKey,pvalMin,dmlSize))
              
              rm_vec <- c("cluster","lociBetaKey","lociPvalKey","lociPvalMin","featuresDml")
              
              ret_str <- optsToCommand(opts=opt_tib, pre=opt$Rscript, exe=par$exePath, rm=rm_vec, add=add_tib, 
                                       file=run_sh, verbose=opt$verbose,vt=1,tc=1,tt=NULL)
              cat(glue::glue("[{par$prgmTag}]:{TAB}. Wrote shell={run_sh}.{RET}{RET}") )
              
              run_id <- paste0('bm-',seed_val,'-cl')
              cmd <- paste(opt$lanExe,run_id,run_sh, sep=' ')
              if (is.null(opt$lanExe) || stringr::str_length(opt$lanExe)==0) cmd <- run_sh
              
              cat(glue::glue("[{par$prgmTag}]:{TAB}. Launching: cmd={cmd}...{RET}{RET}") )
              sys_ret_val <- base::system(cmd)
              
              if (!sys_ret_val)
                cat(glue::glue("[{par$prgmTag}]: Warning: Bad System Return={sys_ret_val}; cmd='{cmd}'{RET}{RET}"))
              # if (!sys_ret_val)
              #   stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Failed System Command({sys_ret_val}); cmd='{cmd}'{RET}{RET}"))
              
              if (opt$single) break
              
            } else {
              
              fTracker <- timeTracker$new(verbose=opt$verbose)
              cur_time_csv <- file.path(cur_opt_dir, paste(curName,'time-tracker.csv.gz', sep='.') )
              
              files_tib <- NULL
              sed_opt_tib <- dml_opt_tib %>% tibble::add_row( Option="seed", Value=as.character(seed_val), Type="seed" )
              
              file_tib <- glmnetRforestWrapper(
                beta=beta_impute_matT, ss=sampleSheet_tib, cgns=cgn_cur_tib, pars=sed_opt_tib, 
                class_idx=class_idx, seed=seed_val, outName=seed_name, dir=seed_dir, 
                alpha_min=opt$alphaMin, alpha_max=opt$alphaMax, opt$alphaInc, parallel=FALSE,
                crossTrain=TRUE,class_var=class_var,cross_perc_min=opt$cross_perc_min,
                verbose=opt$verbose, vt=1, tc=2, tt=fTracker)
              
              # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
              #                       Write Current Params/Options::
              # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
              
              fTracker_tib <- fTracker$time %>% dplyr::mutate_if(is.numeric, list(round), 4)
              readr::write_csv(fTracker_tib, cur_time_csv)
              
            } # if opt$cluster==TRUE
            
            if (opt$single) break
          } # seed_val
          
          if (opt$single) break
        } # dmlSize
        
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                       Write Current Params/Options::
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      cTracker_tib <- cTracker$time %>% dplyr::mutate_if(is.numeric, list(round), 4)
      readr::write_csv(cTracker_tib, dml_time_csv)
      
      if (opt$single) break
    } # pvalMin
    
    if (opt$single) break
  } # pvalKey
  
  if (opt$single) break
} # betaKey

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker_tib <- pTracker$time %>% dplyr::mutate_if(is.numeric, list(round), 4)

opt_csv  <- file.path(opt$outDir, 'program-options.csv')
par_csv  <- file.path(opt$outDir, 'program-parameters.csv')
time_csv <- file.path(opt$outDir, 'time-tracker.csv.gz')

readr::write_csv(opt_tib, opt_csv)
readr::write_csv(par_tib, par_csv)
readr::write_csv(pTracker_tib, time_csv)

sysTime <- Sys.time()
cat(glue::glue("{RET}[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
