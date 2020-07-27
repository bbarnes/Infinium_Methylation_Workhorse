
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Build COVIC Model From Merged Builds::
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
par$prgmTag <- 'predict_builds'

# Illumina based directories::
par$macDir  <- '/Users/bbarnes/Documents/Projects/methylation/tools'
par$lixDir  <- '/illumina/scratch/darkmatter/Projects/COVIC'

# Directory Parameters::
opt$outDir    <- NULL
opt$mergeDir  <- NULL

# To run many models
opt$modelDir  <- NULL

# To run a single model
opt$model       <- NULL
opt$params      <- NULL
opt$features    <- NULL
opt$sampleSheet <- NULL

# Run Parameters::
opt$runName   <- NULL

# Class Parameters::
opt$classVar <- 'Sample_Class'
opt$trainClass <- NULL
opt$seed_dir   <- NULL

# Loci Level Filtering Parameters::
opt$lociBetaKey <- NULL
opt$lociPvalKey <- NULL
opt$lociPvalMin <- NULL

# Loci (Feature) Selection Parameters::
#  opt$featureCsv <- NULL

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
  
  # Model Dir::
  #  opt$seed_dir <- "seed-13"
  
  # Loci Level Filtering Parameters::
  opt$lociBetaKey <- "i_beta"
  opt$lociPvalKey <- "i_poob"
  opt$lociPvalMin <- 0.9
  opt$lociPvalMin <- 1.0
  opt$lociPvalMin <- "1.0,0.1"
  
  
  opt$classVar <- 'Karyotype_0_Call'
  opt$classVar <- 'Karyotype_1_Call'
  opt$classVar <- 'Sample_Name'
  opt$classVar <- 'Sample_Class'
  
  opt$clean  <- TRUE
  opt$clean  <- FALSE
  opt$single <- TRUE
  
  if (opt$classVar=='Sample_Class') {
    opt$version   <- "C0"
    opt$platform  <- "EPIC"
    
    runNameA  <- "COVIC-Set1-15052020"
    runNameB  <- "COVIC-Set5-10062020"
    
    opt$runName   <-  runNameA
    opt$runName   <-  runNameB
    
    opt$modelDir <- paste(
      file.path(par$topDir, 'build_models',opt$platform,opt$version,opt$classVar,opt$runName,
                'ind-beta_i-poob-1','dml-100-Ivana-145'),
      sep=',')
    
    opt$mergeDir  <- paste( 
      file.path(par$topDir, 'merge_builds', opt$platform, opt$version, opt$classVar, opt$runName),
      sep=',')
    
    opt$trainClass <- paste('nSARSCov2', 'pSARSCov2', sep=',')

  } else if (opt$classVar=='Karyotype_0_Call' || opt$classVar=='Karyotype_1_Call') {
    opt$runNameA  <- "COVIC-Set1-15052020"
    opt$runNameB  <- "COVIC-Set5-10062020"
    
    opt$runName  <- 'COVIC-Set5-10062020'
    
    opt$modelDir <- paste(
      file.path(par$topDir, 'build_models/Sample_Class/COVIC-Set5-10062020/i-beta_i-poob-1/dml-100-Ivana-145',opt$seed_dir),
      sep=',')
    
    opt$version   <- "C0"
    opt$platform  <- "EPIC"
    opt$mergeDir  <- paste( file.path(par$topDir, 'merge_builds', opt$classVar, opt$runNameA, opt$platform, opt$version),
                            file.path(par$topDir, 'merge_builds', opt$classVar, opt$runNameB, opt$platform, opt$version),
                            sep=',')

    # opt$trainClass <- paste('Xa','XaXaY','XaXi','XaXiY','XaY', sep=',')
    opt$trainClass <- paste('XaXi','XaY', sep=',')
  } else if (opt$classVar=='Sample_Name') {
    opt$runName  <- "BETA-DELTA-Decoder"
    opt$modelDir <- paste(
      file.path(par$topDir, 'build_models/Sample_Name',opt$runName,'i-beta_i-poob-0.9/dml-100-Ivana-145',opt$seed_dir ),
      sep=',')
    
    opt$version   <- "B4"
    opt$platform  <- "EPIC"
    
    opt$runNameA  <- 'BETA-8x1-EPIC-Core'
    opt$runNameB  <- 'DELTA-8x1-EPIC-Core'
    opt$runNameC  <- 'BETA-DELTA-Core'
    opt$mergeDir <- paste(
      file.path(par$topDir, 'merge_builds',opt$classVar,opt$runNameC,opt$platform,opt$version ),
      sep=',')

    opt$trainClass <- paste('HELA','JURKAT','MCF7','RAJI', sep=',')
  }
  
  opt$outDir <- file.path(par$topDir, par$prgmTag)
  
} else {
  par$runMode    <- 'CommandLine'
  par$exePath <- base::substring(args.dat[grep("--file=", args.dat)], 8)
  
  par$prgmTag <- base::sub('\\.R$', '', base::basename(par$exePath))
  par$locPath <- base::dirname(par$exePath)
  par$scrDir  <- base::dirname(base::normalizePath(par$locPath) )
  par$srcDir  <- base::dirname(base::normalizePath(par$scrDir) )
  par$datDir  <- file.path(par$srcDir, 'dat')
  
  args.dat <- commandArgs(trailingOnly = TRUE)
  option_list = list(
    
    # Directory Parameters::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    make_option(c("--mergeDir"), type="character", default=opt$mergeDir, 
                help="List of Merged Swifthoof Build Directory(s), commas seperated [default= %default]", metavar="character"),
    
    # To run many models
    make_option(c("--modelDir"), type="character", default=opt$modelDir, 
                help="List of Model Build Directory(s). This will launch many models. [default= %default]", metavar="character"),
    
    # To run a single model
    make_option(c("--model"), type="character", default=opt$model, 
                help="Single Model RDS [default= %default]", metavar="character"),
    make_option(c("--params"), type="character", default=opt$params, 
                help="Single Parameters CSV [default= %default]", metavar="character"),
    make_option(c("--features"), type="character", default=opt$features, 
                help="Single Feature Set for model CSV [default= %default]", metavar="character"),
    make_option(c("--sampleSheet"), type="character", default=opt$sampleSheet, 
                help="Single Sample Sheet for model [default= %default]", metavar="character"),

    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    # make_option(c("--seed_dir"), type="character", default=opt$seed_dir, 
    #             help="Run Name [default= %default]", metavar="character"),
    
    # Class Parameters::
    make_option(c("--classVar"), type="character", default=opt$classVar, 
                help="Classification Variable Name [default= %default]", metavar="character"),
    make_option(c("--trainClass"), type="character", default=opt$trainClass, 
                help="Training Class Variable Names, comma delimited [default= %default]", metavar="character"),
    
    # Sample Level Filtering Parameters::
    # make_option(c("--samplePvalName"), type="character", default=opt$samplePvalName, 
    #             help="Pval Method Name to filter Samples for training [default= %default]", metavar="character"),
    # make_option(c("--samplePvalPerc"), type="double", default=opt$samplePvalPerc,
    #             help="Pval Min Percent Passing to filter Sample for training [default= %default]", metavar="double"),
    
    # Loci Level Filtering Parameters::
    make_option(c("--lociBetaKey"), type="character", default=opt$lociBetaKey,
                help="Loci Beta-Method Name (key) for training, comma delimited [default= %default]", metavar="character"),
    make_option(c("--lociPvalKey"), type="character", default=opt$lociPvalKey, 
                help="Loci Pval-Method Name (key) for filtering for training, comma delimited [default= %default]", metavar="character"),
    make_option(c("--lociPvalMin"), type="character", default=opt$lociPvalMin,
                help="Pval Min Passing for loci filtering for training, comma delimited [default= %default]", metavar="character"),
    
    # Loci (Feature) Selection Parameters::
    # make_option(c("--featureCsv"), type="character", default=opt$featureCsv, 
    #             help="Loci feature selection file for training [default= %default]", metavar="character"),
    
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

if (is.null(par$runMode) || is.null(par$prgmTag) || is.null(par$scrDir) || is.null(par$datDir)) {
  
  par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
  par_tib %>% base::print(n=base::nrow(par_tib) )
  # dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value") %>% as.data.frame() %>% print()
  
  if (is.null(par$runMode)) cat(glue::glue("[Usage]: runMode is NULL!!!{RET}"))
  if (is.null(par$prgmTag)) cat(glue::glue("[Usage]: prgmTag is NULL!!!{RET}"))
  if (is.null(par$scrDir))  cat(glue::glue("[Usage]: scrDir is NULL!!!{RET}"))
  if (is.null(par$datDir))  cat(glue::glue("[Usage]: darDir is NULL!!!{RET}"))
  base::stop("Null Parameters!\n\n")
}

if (is.null(opt$outDir) || is.null(opt$mergeDir) || 
    (is.null(opt$modelDir) ||
     (is.null(opt$model) && is.null(opt$params) && is.null(opt$features) && is.null(opt$sampleSheet) ) ) ||
    is.null(opt$runName) || 
    # is.null(opt$seed_dir) || 
    is.null(opt$classVar) || is.null(opt$trainClass) ||
    # is.null(opt$samplePvalName) || is.null(opt$samplePvalPerc) ||
    is.null(opt$lociBetaKey) || is.null(opt$lociPvalKey) || is.null(opt$lociPvalMin) || 
    # is.null(opt$featureCsv) ||
    is.null(opt$percisionBeta) || is.null(opt$percisionPval) || 
    is.null(opt$execute) || is.null(opt$single) || is.null(opt$parallel) || is.null(opt$cluster) ||
    
    is.null(opt$clean) || is.null(opt$Rscript) || is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )
  # OLD data.frame print method::
  # dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value") %>% as.data.frame() %>% print()
  
  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  if (is.null(opt$outDir))     cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$mergeDir))   cat(glue::glue("[Usage]: mergeDir is NULL!!!{RET}"))
  # if (is.null(opt$modelDir))   cat(glue::glue("[Usage]: modelDir is NULL!!!{RET}"))
  
  (is.null(opt$modelDir) ||
      (is.null(opt$model) && is.null(opt$params) && is.null(opt$features) && is.null(opt$sampleSheet) ) ) {
    cat(glue::glue("[Usage]: modelDir is NULL OR model,params,features and sampleSheet are NULL!!!{RET}"))
  }
  if (is.null(opt$runName))    cat(glue::glue("[Usage]: runName is NULL!!!{RET}"))

  if (is.null(opt$classVar))   cat(glue::glue("[Usage]: classVar is NULL!!!{RET}"))
  if (is.null(opt$trainClass)) cat(glue::glue("[Usage]: trainClass is NULL!!!{RET}"))
  
  # if (is.null(opt$samplePvalName)) cat(glue::glue("[Usage]: samplePvalName is NULL!!!{RET}"))
  # if (is.null(opt$samplePvalPerc)) cat(glue::glue("[Usage]: samplePvalPerc is NULL!!!{RET}"))
  
  if (is.null(opt$lociBetaKey)) cat(glue::glue("[Usage]: lociBetaKey is NULL!!!{RET}"))
  if (is.null(opt$lociPvalKey)) cat(glue::glue("[Usage]: lociPvalKey is NULL!!!{RET}"))
  if (is.null(opt$lociPvalMin)) cat(glue::glue("[Usage]: lociPvalMin is NULL!!!{RET}"))
  # if (is.null(opt$seed_dir)) cat(glue::glue("[Usage]: seed_dir is NULL!!!{RET}"))
  
  # if (is.null(opt$featureCsv)) cat(glue::glue("[Usage]: featureCsv is NULL!!!{RET}"))
  
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Preprocessing:: General
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

mergeDirs_vec   <- opt$mergeDir %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
trainClass_vec  <- opt$trainClass %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

lociBetaKey_vec <- opt$lociBetaKey %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
lociPvalKey_vec <- opt$lociPvalKey %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
lociPvalMin_vec <- opt$lociPvalMin %>% str_split(pattern=',', simplify=TRUE) %>% as.double() %>% as.vector()

# if (is.null(opt$classVar)) opt$classVar <- 'Source_Sample_Name'
class_mat <- rlang::sym(paste(opt$classVar,'CompType', sep='_'))
class_org <- rlang::sym(paste(opt$classVar,'Origin', sep='_'))
class_var <- rlang::sym(opt$classVar)
class_idx <- rlang::sym("Class_Idx")
# opt$samplePvalName <- opt$samplePvalName %>% rlang::sym()

opt <- setLaunchExe(opts=opt, pars=par, verbose=opt$verbose, vt=5,tc=0)

class_tib <- tibble::tibble(Sample_Class=trainClass_vec) %>% 
  dplyr::mutate(Sample_Class=as.factor(Sample_Class), Sample_Idx=as.integer(Sample_Class)-1)
class_nil <- as.character( class_tib$Sample_Class[length(class_tib$Sample_Class)] )

cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#     Preprocessing:: File Identification:: Model/Feature/SS Loading
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Workflow Logic::
#  if modelDir is present then search for files and launch each one
#  else launch current model with each test set
#
if (!is.null(opt$modelDir)) {
  
  search_fn_key <- '.model-files.csv.gz$'
  full_ss_csvs <- list.files(opt$modelDir, pattern=search_fn_key, full.names=TRUE, recursive=TRUE)
  full_ss_cnts <- length(full_ss_csvs)
  cat(glue::glue("[{par$prgmTag}]: full_ss_cnts={full_ss_cnts}.{RET}") )
  
  stopifnot(full_ss_cnts>0)

  opt$outDir <- file.path(opt$outDir, opt$classVar, opt$runName)
  if (!dir.exists(opt$outDir)) dir.create(opt$outDir, opt$classVar, recursive=TRUE)
  
  # opt$shellDir <- file.path(opt$outDir, 'shells')
    
  for (mIdx in c(1:full_ss_cnts)) {
    
    # Gather all files::
    fns_csv <- full_ss_csvs[mIdx]
    stopifnot(file.exists(fns_csv))
    dir_path <- base::dirname(fns_csv)
    
    fns_tib <- suppressMessages(suppressWarnings( readr::read_csv(fns_csv) ))
    sam_csv <- file.path(dir_path, fns_tib %>% dplyr::filter(Type=="SampleSheet") %>% head(n=1) %>% dplyr::pull(File_Name))
    par_csv <- file.path(dir_path, fns_tib %>% dplyr::filter(Type=="Params") %>% head(n=1) %>% dplyr::pull(File_Name))
    fet_csv <- file.path(dir_path, fns_tib %>% dplyr::filter(Type=="Features") %>% head(n=1) %>% dplyr::pull(File_Name))
    mod_rds <- file.path(dir_path, fns_tib %>% dplyr::filter(Type=="Model") %>% head(n=1) %>% dplyr::pull(File_Name))
    
    # Ensure all data exists...
    #
    #  TBD: Skip if files don't exist...
    if (is.null(mod_rds) || length(mod_rds)==0 || !file.exists(mod_rds)) {
      cat(glue::glue("[{par$prgmTag}]:{TAB} Unable to find mod_rds file; dir_path={dir_path}; fns_tib={fns_csv}. Skipping...{RET}") )
      next
    }
    
    stopifnot(file.exists(sam_csv))
    stopifnot(file.exists(par_csv))
    stopifnot(file.exists(par_csv))

    par_tib <- suppressMessages(suppressWarnings( readr::read_csv(par_csv) )) %>% 
      dplyr::mutate(Value_Str=stringr::str_replace_all(Value,',','-'))
    
    cat(glue::glue("[{par$prgmTag}]:{TAB} Found all paths for dir_path={dir_path}.{RET}{RET}") )

    user_dir  <- par_tib %>% dplyr::filter(Type=='User') %>% dplyr::pull(Value_Str) %>% as.vector() %>% stringr::str_c(collapse="_")
    fets_dir  <- par_tib %>% dplyr::filter(Type=='Feature') %>% dplyr::pull(Value_Str) %>% as.vector() %>% stringr::str_c(collapse="_")
    seed_dir  <- stringr::str_c(
      "seed",par_tib %>% dplyr::filter(Type=='seed') %>% dplyr::pull(Value_Str) %>% as.vector() %>% stringr::str_c(collapse="_"),
      sep="-")
    orig_dir  <- stringr::str_c(
      "oSize",par_tib %>% dplyr::filter(Type=='Original') %>% dplyr::pull(Value_Str) %>% as.vector() %>% stringr::str_c(collapse="_"),
      sep="-")
    train_dir <- par_tib %>% dplyr::filter(Type=='train') %>% dplyr::pull(Value_Str) %>% as.vector() %>% stringr::str_c(collapse="_")
    
    for (mergeDir in mergeDirs_vec) {
      if (!dir.exists(mergeDir)) {
        cat(glue::glue("{RET}[{par$prgmTag}]: Warning: mergeDir={mergeDir} does not exist! Skipping...{RET}{RET}"))
        next
      }
      merge_name <- base::basename(mergeDir)
      
      cur_dir <- file.path(opt$outDir,user_dir,fets_dir,seed_dir,orig_dir,train_dir,merge_name)
      if (!dir.exists(cur_dir)) dir.create(cur_dir, recursive=TRUE)
      
      shell_dir <- file.path(cur_dir, "test.model.sh")
      
      add_tib <- tibble::tibble(Option=c("outDir","sampleSheet","params","features","model"),
                                Value=c(cur_dir,sam_csv,par_csv,fet_csv,mod_rds))
      
      rm_vec <- c("modelDir")
      
      cmd <- optsToCommand(opts=opt_tib, exe=par$exePath, rm=rm_vec, add=add_tib, 
                           verbose=opt$verbose,vt=1,tc=1,tt=NULL)
      
      

      break
    }
    break
  }

} else {
  
}

















#
# OLD CODE::
#
if (FALSE) {
  full_ss_csvs <- list.files(opt$modelDir, pattern='.model-SampleSheet.csv.gz$', full.names=TRUE, recursive=TRUE)
  full_ss_cnts <- length(full_ss_csvs)
  cat(glue::glue("[{par$prgmTag}]: full_ss_cnts={full_ss_cnts}.{RET}") )
  
  stopifnot(full_ss_cnts>0)
  
  seed_dirs <- full_ss_csvs %>% base::dirname() %>% base::dirname() %>% unique()
  seed_cnts <- seed_dirs %>% length()
  cat(glue::glue("[{par$prgmTag}]: seed_dirs={seed_cnts}.{RET}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #          Preprocess each Seed Build Directory Individually::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  for (seed_dir in seed_dirs) {
    cur_ss_csvs <- list.files(opt$modelDir, pattern='.model-SampleSheet.csv.gz$', full.names=TRUE, recursive=TRUE)
    train_ss_cnts <- length(cur_ss_csvs)
    cat(glue::glue("[{par$prgmTag}]: train_ss_cnts={train_ss_cnts}.{RET}") )
    stopifnot(train_ss_cnts>0)
    
    break
  }
  
  opt$outDir <- file.path(opt$outDir, opt$classVar, opt$runName, opt$seed_dir)
  if (!dir.exists(opt$outDir)) dir.create(opt$outDir, opt$classVar, recursive=TRUE)
  
  mods <- base::list()
  cgns <- base::list()
  sams <- base::list()
  pars <- base::list()
  for (mIdx in c(1:train_ss_cnts)) {
    train_ss_csv <- cur_ss_csvs[mIdx]
    train_md_rds <- stringr::str_replace(train_ss_csv, '.model-SampleSheet.csv.gz$', '.model.rds')
    train_pr_csv <- stringr::str_replace(train_ss_csv, '.model-SampleSheet.csv.gz$', '.model-params.csv.gz')
    train_cg_csv <- stringr::str_replace(train_ss_csv, '.model-SampleSheet.csv.gz$', '.model-features.csv.gz')
    # train_md_key <- train_md_rds %>% base::basename() %>% stringr::str_remove('.model.rds$') %>% 
    #   stringr::str_replace('^.*_([^_]+)$', "\\$1") %>% stringr::str_remove('\\\\+')
    
    train_md_key <- train_md_rds %>% stringr::str_remove('.rds$')
    if (length(grep(opt$seed_dir, train_md_key))!=1 ) next
    
    stopifnot(file.exists(train_ss_csv))
    stopifnot(file.exists(train_cg_csv))
    stopifnot(file.exists(train_md_rds))
    stopifnot(!is.null(train_md_key))
    stopifnot(length(train_md_key)>0)
    
    mods[[train_md_key]] = readr::read_rds(train_md_rds)
    cgns[[train_md_key]] = suppressMessages(suppressWarnings( readr::read_csv(train_cg_csv) ))
    sams[[train_md_key]] = suppressMessages(suppressWarnings( readr::read_csv(train_ss_csv) ))
    pars[[train_md_key]] = suppressMessages(suppressWarnings( readr::read_csv(train_pr_csv) ))
    
    break
  }
  cgns_max_tib <- cgns %>% dplyr::bind_rows() %>% dplyr::distinct(Probe_ID) %>% dplyr::arrange(Probe_ID)
  
  mods_cnt <- mods %>% names() %>% length()
  sams_cnt <- sams %>% names() %>% length()
  cgns_cnt <- cgns %>% names() %>% length()
  pars_cnt <- pars %>% names() %>% length()
  cgns_max <- cgns_max_tib %>% base::nrow()
  
  cat(glue::glue("[{par$prgmTag}]: Done. models={mods_cnt}, SampleSheets={sams_cnt}, Features={cgns_cnt}, ",
                 "Params={pars_cnt}, cgns_max={cgns_max}; outDir={opt$outDir}.{RET}{RET}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                Load/Filter all Testing Data into Matricies
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  for (betaKey in lociBetaKey_vec) {
    for (pvalKey in lociPvalKey_vec) {
      all_str <- paste(stringr::str_replace_all(betaKey,'_', '-'),
                       stringr::str_replace_all(pvalKey,'_', '-'), sep='_')
      
      all_line_tib <- NULL
      all_line_csv <- file.path(opt$outDir, paste(opt$runName,paste(all_str, 'SamplePrediction.sheet.csv.gz', sep='.'), sep='_' ) )
      for (pvalMin in lociPvalMin_vec) {
        betaStr <- betaKey %>% stringr::str_replace_all('_', '-')
        pvalStr <- paste(pvalKey %>% stringr::str_replace_all('_', '-'), pvalMin, sep='-')
        dirName <- paste(betaStr,pvalStr, sep='_')
        outName <- paste(opt$runName, dirName, sep='_')
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                     Build Current Output Directory::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        cur_opt_dir <- file.path(opt$outDir, dirName)
        if (!dir.exists(cur_opt_dir)) dir.create(cur_opt_dir, recursive=TRUE)
        cat(glue::glue("[{par$prgmTag}]: Built; cur_opt_dir={cur_opt_dir}!{RET}") )
        
        if (opt$clean) unlink(list.files(cur_opt_dir, full.names=TRUE))
        
        cTracker <- timeTracker$new(verbose=opt$verbose)
        
        labs_ss_tib <- NULL
        beta_mat <- NULL
        for (curDir in mergeDirs_vec) {
          
          # Find Sample Sheet
          #
          curr_ss_csv <- findFileByPattern(dir=curDir, patter='_AutoSampleSheet.csv.gz$', max=1, recursive=FALSE, verbose=opt$verbose)
          curr_fn_csv <- curr_ss_csv %>% stringr::str_replace('_AutoSampleSheet.csv.gz$','_MergedDataFiles.tib.csv.gz')
          base_dir <- base::dirname(curr_ss_csv)
          stopifnot(file.exists(curr_ss_csv), file.exists(curr_fn_csv))
          
          cat(glue::glue("[{par$prgmTag}]:{TAB} {outName} Found sample_csv={curr_ss_csv}.{RET}") )
          cat(glue::glue("[{par$prgmTag}]:{TAB} {outName} Found call_table={curr_fn_csv}.{RET}") )
          
          # Load and filter::
          curr_ss_tib <- suppressMessages(suppressWarnings( readr::read_csv(curr_ss_csv) )) # %>% 
          # dplyr::filter(!!opt$samplePvalName > !!opt$samplePvalPerc) %>% 
          # dplyr::filter(!!class_var %in% trainClass_vec) %>%
          # dplyr::mutate(Source_Sample_ID=as.integer(Source_Sample_ID))
          
          calls_path_tib <- suppressMessages(suppressWarnings( readr::read_csv(curr_fn_csv) ))
          
          betas_path_tib <- calls_path_tib %>% dplyr::filter(Method %in% c(betaKey))
          pvals_path_tib <- calls_path_tib %>% dplyr::filter(Method %in% c(pvalKey))
          
          stopifnot(base::nrow(betas_path_tib)==1)
          stopifnot(base::nrow(pvals_path_tib)==1)
          
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          #            Load Beta/Pval And Merge into Previous Matrix::
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          
          beta_csv <- file.path(base_dir, base::basename(betas_path_tib$Full_Path[1]) )
          pval_csv <- file.path(base_dir, base::basename(pvals_path_tib$Full_Path[1]) )
          
          beta_mat <- loadCallsMatrix(betaCSV=beta_csv, pvalCSV=pval_csv, minPval=pvalMin, mat=beta_mat, 
                                      cgn=cgns_max_tib, ss=curr_ss_tib,
                                      verbose=opt$verbose, tc=1, tt=pTracker)
          
          labs_ss_tib <- labs_ss_tib %>% dplyr::bind_rows(curr_ss_tib)
          
          # beta_mat %>% dim() %>% print()
          cat(glue::glue("[{par$prgmTag}]:{TAB}Done. {outName}.{RET}{RET}") )
          
          # break
        }
        if (length(grep("Phase_Num", names(labs_ss_tib)))==0) {
          labs_ss_tib <- labs_ss_tib %>% dplyr::mutate(Phase_Num=1)
        }
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                 Seed-based 3-fold Sample Sheet Partioning::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        # TBD:: 
        #  - Loop over set seeds
        #  - Add partion index to each sample
        #  - Loop over each partition for training (i.e. exclude the thrid to be tested)
        #  - Train on all partitions (i.e. exclude partion 0)
        #  - Possibly test on missing third as validation for later full testing
        #
        
        sort_s1_tib <- labs_ss_tib %>% dplyr::filter(  (!!class_var %in% class_tib$Sample_Class) ) %>%
          dplyr::mutate(!!class_org := !!class_var)
        sort_s2_tib <- labs_ss_tib %>% dplyr::filter(! (!!class_var %in% class_tib$Sample_Class) ) %>%
          dplyr::mutate(!!class_org := !!class_var, 
                        !!class_var := dplyr::case_when(
                          !!class_org=="Xa" ~ "XaXi",
                          TRUE ~ !!class_nil) )
        
        # Try sorting Labeled/Merged Sample Sheet by Sample_Class and then re-order matrix by new order::
        sort_ss_tib <- dplyr::bind_rows(sort_s1_tib,sort_s2_tib) %>% 
          dplyr::mutate(
            # !!class_org := !!class_var,
            !!class_mat := dplyr::case_when(
              !!class_org == !!class_var ~ 'Known',
              !!class_org != !!class_var ~ 'Novel',
              TRUE ~ NA_character_
            )
          ) %>%
          dplyr::arrange(!!class_var) %>% 
          dplyr::mutate(!!class_var := as.factor(!!class_var),
                        !!class_idx := as.integer(!!class_var)-1) %>% 
          dplyr::mutate(!!class_idx := as.integer(!!class_idx) ) %>%
          dplyr::select(Sentrix_Name, !!class_var, !!class_idx, !!class_org, !!class_mat) # %>% as.data.frame()
        
        # QC Sanity Check:: Make sure the new ordering works::
        if (FALSE) {
          cbind(
            beta_mat %>% colnames(),
            sort_ss_tib %>% dplyr::pull(Sentrix_Name),
            beta_mat[ , dplyr::pull(sort_ss_tib, Sentrix_Name)] %>% colnames()
          ) %>% tibble::as_tibble() %>% dplyr::filter(V2==V3)
        }
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                   Build Raw and Imputed Sorted Matricies::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        beta_raw_matT <- beta_mat[ , dplyr::pull(sort_ss_tib, Sentrix_Name) ]
        beta_imp_mat  <- beta_raw_matT %>% as.data.frame() %>% 
          purrr::set_names(paste('X',names(.), sep='_') ) %>% 
          makeX_glmnet_imputeNA(na.impute = TRUE) %>% 
          as.data.frame() %>% purrr::set_names(stringr::str_remove_all(names(.), '^X_')) %>%
          as.matrix() %>% t()
        
        raw_nan_cnt <- which(is.na(beta_raw_matT)) %>% length()
        imp_nan_cnt <- which(is.na(beta_imp_mat)) %>% length()
        if (imp_nan_cnt!=0) stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Failed Imputation: imp_nan_cnt={imp_nan_cnt}, raw_nan_cnt={raw_nan_cnt}{RET}{RET}"))
        cat(glue::glue("[{par$prgmTag}]:{TAB} [{outName}] Imputation Results: imp_nan_cnt={imp_nan_cnt}, raw_nan_cnt={raw_nan_cnt}{RET}") )
        
        # Build Label Data Structures::
        labs_key_vec <- sort_ss_tib %>% dplyr::pull(!!class_var) %>% as.vector()
        labs_idx_vec <- sort_ss_tib %>% dplyr::pull(!!class_idx) %>% as.integer() %>% as.vector()
        beta_labs_df <- sort_ss_tib %>% column_to_rownames(var="Sentrix_Name") %>% as.data.frame()
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                          Predict Data with Models::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        modelNames <- mods %>% names()
        
        # break
        for (curModName in modelNames) {
          # curModName <- "rf-1se-0.1"
          
          curModKey <- curModName %>% base::basename() %>% stringr::str_remove('.model$') %>% 
            stringr::str_replace('^.*_([^_]+)$', "\\$1") %>% stringr::str_remove('\\\\+')
          
          # if (stringr::str_starts(curModName, 'glmnet-') ) {
          if (stringr::str_detect(curModName, 'glmnet-') ) {
            raw_pred_tib <- predGlmnet(mod=mods[[curModName]], data=beta_imp_mat[, cgns[[curModName]]$Probe_ID ], 
                                       labs=labs_idx_vec, name=curModKey, lambda="lambda.1se", verbose=opt$verbose)
            
            # predict(object=mods[[curModName]], newx=beta_imp_mat[, cgns[[curModName]]$Probe_ID ], s="lambda.1se", type="class")
            # predict(object=mods[[curModName]], newx=beta_imp_mat[, cgns[["glmnet-0"]]$Probe_ID ], s=mods[[curModName]]$lambda.1se, type="class")
            # predict(object=mods[[curModName]], newx=beta_imp_mat[, cgns[[curModName]]$Probe_ID ] %>% t(), s=mods[[curModName]]$lambda, type="class")
            
            # } else if (stringr::str_starts(curModName, 'rf-') ) {
          } else if (stringr::str_detect(curModName, 'rforest') ) {
            raw_pred_tib <- predRandomForest(mod=mods[[curModName]], data=beta_imp_mat[, cgns[[curModName]]$Probe_ID ], 
                                             labs=labs_idx_vec, name=curModKey, verbose=opt$verbose)
            
          } else {
            stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Unsupported Model prefix type={curModKey}!!!{RET}{RET}"))
          }
          
          # Accurcy by Sample
          samp_pred_tib  <- raw_pred_tib %>% 
            dplyr::mutate(True_Class=labs_idx_vec) %>%
            
            # dplyr::rename(Sentrix_Name=Full_Sample) %>%
            # dplyr::mutate(Pred_Class=as.integer(Pred_Class), # Pred_Class_Raw=as.integer(Pred_Class_Raw), 
            #               True_Class=labs_idx_vec, Sentrix_Name=as.character(Sentrix_Name)) %>%
            dplyr::left_join(sams[[curModName]], by="Sentrix_Name") %>%
            dplyr::mutate(isTrainSample=dplyr::case_when(
              is.na(!!class_var) ~ 'tests',
              TRUE ~ 'train'),
              Call=dplyr::case_when(
                Pred_Class==True_Class ~ 'TP',
                Pred_Class!=True_Class ~ 'FP',
                TRUE ~ NA_character_
              )
            ) %>% 
            dplyr::select(-!!class_var, -Class_Idx) %>%
            dplyr::left_join(sort_ss_tib, by='Sentrix_Name') %>% 
            dplyr::left_join(dplyr::select(labs_ss_tib, -!!class_var),  by="Sentrix_Name")
          
          # Generate Accuracy Summary::
          #   Phase_Num
          group_vec <- c("class_mat", "isTrainSample", "True_Class", "Call")
          samp_sum_tmp_tib <- samp_pred_tib %>% 
            dplyr::distinct(Sentrix_Name, .keep_all=TRUE) %>% 
            dplyr::group_by(!!class_mat, Phase_Num, isTrainSample, True_Class, Call) %>% 
            dplyr::summarise(Pred_Mode_Count=n()) %>% 
            dplyr::ungroup()
          
          cat(glue::glue("[{par$prgmTag}]:{TAB} [{outName}] samp_sum_tmp_tib=...{RET}") )
          print(samp_sum_tmp_tib)
          
          cat(glue::glue("[{par$prgmTag}]:{TAB} [{outName}] samp_sum_tmp_tib=...{RET}") )
          samp_cnt_tib <- samp_sum_tmp_tib %>% 
            dplyr::group_by(!!class_mat, Phase_Num, isTrainSample, True_Class) %>%
            dplyr::summarise(Total_Count=n(),
                             Pred_Accuracy=round(100*Pred_Mode_Count/Total_Count, 3) )
          print(samp_cnt_tib)
          
          samp_sum_tib <- samp_sum_tmp_tib %>%
            dplyr::group_by(!!class_mat, Phase_Num, isTrainSample, True_Class) %>%
            dplyr::add_tally(wt=Pred_Mode_Count, name="Total_Count") %>% 
            dplyr::mutate(Pred_Accuracy=round(100*Pred_Mode_Count/Total_Count, 3) )
          
          cat(glue::glue("[{par$prgmTag}]:{TAB} [{outName}] samp_sum_tib=...{RET}") )
          print(samp_sum_tib)
          
          # Single line summary::
          # TBD::
          #  - Add additional parameters to output {seed, method, basically-everything-from-params-file... }
          #
          pre_line_tib <- pars[[curModName]] %>% tidyr::unite(Option, Type,Option, sep='_') %>% 
            dplyr::mutate(Option=stringr::str_replace_all(Option, ',',';'),
                          Value=stringr::str_replace_all(Value, ',',';')) %>% 
            tidyr::spread(Option, Value)
          sum_line_tib <- samp_sum_tib %>% dplyr::ungroup() %>% 
            tidyr::gather(Metric, value, -c(!!class_mat, Phase_Num, isTrainSample, True_Class, Call) ) %>% 
            tidyr::unite(key, !!class_mat, Phase_Num, isTrainSample,True_Class,Call,Metric, sep='_') %>% 
            tidyr::spread(key, value)
          
          out_line_tib <- dplyr::bind_cols(pre_line_tib,sum_line_tib)
          
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          #                        Output Prediction Results::
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          
          # Need to add to cur_mod_dir::
          if (!is.null(out_line_tib$Feature_featureSizeDml) &&
              !is.null(out_line_tib$Feature_featureNamePre)) {
            dml_cur_str <- paste('dml',out_line_tib$Feature_featureSizeDml, out_line_tib$Feature_featureNamePre, sep='-')
          } else if (!is.null(out_line_tib$Feature_featureSizeDml)) {
            dml_cur_str <- paste('dml',out_line_tib$Feature_featureSizeDml, sep='-')
          } else if (!is.null(out_line_tib$Feature_featureNamePre)) {
            dml_cur_str <- paste(out_line_tib$Feature_featureNamePre, sep='-')
          } else {
            stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Both Feature_featureSizeDml AND Feature_featureNamePre are NULL!!!{RET}{RET}"))
          }
          seed_str <- paste('seed',out_line_tib$seed_sedd, sep='-')
          model_name <- paste(outName,dml_cur_str,seed_str,curModKey, sep='_')
          
          cur_mod_dir <- file.path(cur_opt_dir, dml_cur_str, seed_str, curModKey)
          if (!dir.exists(cur_mod_dir)) dir.create(cur_mod_dir, recursive=TRUE)
          
          samp_pred_csv <- file.path(cur_mod_dir, paste(model_name,'SamplePrediction.csv.gz', sep='.') )
          samp_sum_csv  <- file.path(cur_mod_dir, paste(model_name,'SamplePrediction.tib.csv.gz', sep='.') )
          out_line_csv  <- file.path(cur_mod_dir, paste(model_name,'SamplePrediction.sheet.csv.gz', sep='.') )
          
          out_line_tib <- out_line_tib %>% dplyr::mutate(Path=out_line_csv)
          
          readr::write_csv(samp_pred_tib, samp_pred_csv)
          readr::write_csv(samp_sum_tib,  samp_sum_csv)
          readr::write_csv(out_line_tib,  out_line_csv)
          
          # Bind all data::
          all_line_tib <- dplyr::bind_rows(all_line_tib, out_line_tib)
          
          if (opt$single) break
        } # curModName
        
        if (opt$single) break
      } # pvalMin
      
      readr::write_csv(all_line_tib,  all_line_csv)
      
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
}

sysTime <- Sys.time()
cat(glue::glue("{RET}[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
