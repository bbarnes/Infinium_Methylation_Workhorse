
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
  opt$single <- FALSE
  opt$single <- TRUE
  opt$execute <-  FALSE
  
  if (FALSE) {
    # cat scratch/predict_builds/Sample_Class/COVIC-Set1-15052020/COVIC-Set1-15052020_nSARSCov2-pSARSCov2_Poob_Pass_0_Perc_96_ind_beta_i_poob_1/1000_Ivana-145/seed-13/oSize-30/rforest_0.5_lambda-1se/test.model.sh | perl -pe 's/ /\n/gi' | perl -pe 's/^--//; s/=(.*)$/="$1"/; print "opt\$"'

    opt$classVar="Sample_Class"
    opt$lociBetaKey="i_beta"
    opt$lociPvalKey="i_poob"
    opt$lociPvalMin="1.0,0.1"
    opt$mergeDir="/Users/bbarnes/Documents/Projects/methylation/scratch/merge_builds/EPIC/C0/Sample_Class/COVIC-Set5-10062020"
    opt$outDir="/Users/bbarnes/Documents/Projects/methylation/scratch/predict_builds"
    opt$percisionBeta="4"
    opt$percisionPval="6"
    opt$Rscript="Rscript"
    opt$runName="COVIC-Set1-15052020"
    opt$trainClass="nSARSCov2,pSARSCov2"
    opt$verbose="3"
    opt$outDir="/Users/bbarnes/Documents/Projects/methylation/scratch/predict_builds/Sample_Class/COVIC-Set1-15052020/COVIC-Set1-15052020_nSARSCov2-pSARSCov2_Poob_Pass_0_Perc_96_ind_beta_i_poob_1/1000_Ivana-145/seed-13/oSize-30/rforest_0.5_lambda-1se"
    opt$sampleSheet="/Users/bbarnes/Documents/Projects/methylation/scratch/build_models/EPIC/C0/Sample_Class/COVIC-Set1-15052020/ind-beta_i-poob-1/dml-1000-Ivana-145/seed-13/alpha-0.5/rforest-0.5-lambda-1se.model-SampleSheet.csv.gz"
    opt$params="/Users/bbarnes/Documents/Projects/methylation/scratch/build_models/EPIC/C0/Sample_Class/COVIC-Set1-15052020/ind-beta_i-poob-1/dml-1000-Ivana-145/seed-13/alpha-0.5/rforest-0.5-lambda-1se.model-params.csv.gz"
    opt$features="/Users/bbarnes/Documents/Projects/methylation/scratch/build_models/EPIC/C0/Sample_Class/COVIC-Set1-15052020/ind-beta_i-poob-1/dml-1000-Ivana-145/seed-13/alpha-0.5/rforest-0.5-lambda-1se.model-features.csv.gz"
    opt$model="/Users/bbarnes/Documents/Projects/methylation/scratch/build_models/EPIC/C0/Sample_Class/COVIC-Set1-15052020/ind-beta_i-poob-1/dml-1000-Ivana-145/seed-13/alpha-0.5/rforest-0.5-lambda-1se.model.rds"
    
  } else if (opt$classVar=='Sample_Class') {
    version   <- "C0"
    platform  <- "EPIC"
    
    runNameA  <- "COVIC-Set1-15052020"
    runNameB  <- "COVIC-Set7-06082020"
    runNameB  <- "COVIC-Set5-10062020"
    
    opt$runName   <-  runNameA
    
    opt$modelDir <- paste(
      file.path(par$topDir, 'build_models',platform,version,opt$classVar,opt$runName,'ind-beta_i-poob-1','dml-1000-Ivana-145'),
      sep=',')
    
    opt$mergeDir  <- paste(
      file.path(par$topDir, 'merge_builds', platform, version, opt$classVar, runNameB),
      sep=',')
    
    opt$trainClass <- paste('nSARSCov2', 'pSARSCov2', sep=',')

  } else if (opt$classVar=='Karyotype_0_Call' || opt$classVar=='Karyotype_1_Call') {
    opt$runNameA  <- "COVIC-Set1-15052020"
    opt$runNameB  <- "COVIC-Set5-10062020"
    
    opt$runName  <- 'COVIC-Set5-10062020'
    
    opt$modelDir <- paste(
      file.path(par$topDir, 'build_models/Sample_Class/COVIC-Set5-10062020/i-beta_i-poob-1/dml-100-Ivana-145',opt$seed_dir),
      sep=',')
    
    version   <- "C0"
    platform  <- "EPIC"
    opt$mergeDir  <- paste( file.path(par$topDir, 'merge_builds', opt$classVar, opt$runNameA, platform, version),
                            file.path(par$topDir, 'merge_builds', opt$classVar, opt$runNameB, platform, version),
                            sep=',')

    # opt$trainClass <- paste('Xa','XaXaY','XaXi','XaXiY','XaY', sep=',')
    opt$trainClass <- paste('XaXi','XaY', sep=',')
  } else if (opt$classVar=='Sample_Name') {
    opt$runName  <- "BETA-DELTA-Decoder"
    opt$modelDir <- paste(
      file.path(par$topDir, 'build_models/Sample_Name',opt$runName,'i-beta_i-poob-0.9/dml-100-Ivana-145',opt$seed_dir ),
      sep=',')
    
    version   <- "B4"
    platform  <- "EPIC"
    
    opt$runNameA  <- 'BETA-8x1-EPIC-Core'
    opt$runNameB  <- 'DELTA-8x1-EPIC-Core'
    opt$runNameC  <- 'BETA-DELTA-Core'
    opt$mergeDir <- paste(
      file.path(par$topDir, 'merge_builds',opt$classVar,opt$runNameC,platform,version ),
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
    
    # Loci Level Filtering Parameters::
    make_option(c("--lociBetaKey"), type="character", default=opt$lociBetaKey,
                help="Loci Beta-Method Name (key) for training, comma delimited [default= %default]", metavar="character"),
    make_option(c("--lociPvalKey"), type="character", default=opt$lociPvalKey, 
                help="Loci Pval-Method Name (key) for filtering for training, comma delimited [default= %default]", metavar="character"),
    make_option(c("--lociPvalMin"), type="character", default=opt$lociPvalMin,
                help="Pval Min Passing for loci filtering for training, comma delimited [default= %default]", metavar="character"),
    
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
    is.null(opt$runName) || 
    is.null(opt$classVar) || is.null(opt$trainClass) ||
    is.null(opt$lociBetaKey) || is.null(opt$lociPvalKey) || is.null(opt$lociPvalMin) || 
    is.null(opt$percisionBeta) || is.null(opt$percisionPval) || 
    is.null(opt$execute) || is.null(opt$single) || is.null(opt$parallel) || is.null(opt$cluster) ||
    
    is.null(opt$clean) || is.null(opt$Rscript) || is.null(opt$verbose) ) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )

  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  if (is.null(opt$outDir))     cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$mergeDir))   cat(glue::glue("[Usage]: mergeDir is NULL!!!{RET}"))

  if (is.null(opt$modelDir) ||
      (is.null(opt$model) && is.null(opt$params) && is.null(opt$features) && is.null(opt$sampleSheet) ) ) {
    cat(glue::glue("[Usage]: modelDir is NULL OR model,params,features and sampleSheet are NULL!!!{RET}"))
    
    if (is.null(opt$modelDir))    cat(glue::glue("[Usage]: modelDir is NULL!!!{RET}"))
    if (is.null(opt$model))       cat(glue::glue("[Usage]: model is NULL!!!{RET}"))
    if (is.null(opt$params))      cat(glue::glue("[Usage]: params is NULL!!!{RET}"))
    if (is.null(opt$features))    cat(glue::glue("[Usage]: features is NULL!!!{RET}"))
    if (is.null(opt$sampleSheet)) cat(glue::glue("[Usage]: sampleSheet is NULL!!!{RET}"))
  }
  if (is.null(opt$runName))    cat(glue::glue("[Usage]: runName is NULL!!!{RET}"))

  if (is.null(opt$classVar))   cat(glue::glue("[Usage]: classVar is NULL!!!{RET}"))
  if (is.null(opt$trainClass)) cat(glue::glue("[Usage]: trainClass is NULL!!!{RET}"))
  
  if (is.null(opt$lociBetaKey)) cat(glue::glue("[Usage]: lociBetaKey is NULL!!!{RET}"))
  if (is.null(opt$lociPvalKey)) cat(glue::glue("[Usage]: lociPvalKey is NULL!!!{RET}"))
  if (is.null(opt$lociPvalMin)) cat(glue::glue("[Usage]: lociPvalMin is NULL!!!{RET}"))

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

class_mat <- rlang::sym(paste(opt$classVar,'CompType', sep='_'))
class_org <- rlang::sym(paste(opt$classVar,'Origin', sep='_'))
class_var <- rlang::sym(opt$classVar)
class_idx <- rlang::sym("Class_Idx")

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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                  Main::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (!is.null(opt$modelDir)) {
  cat(glue::glue("[{par$prgmTag}]: Launching in cluster mode...{RET}"))
  
  opt$outDir <- file.path(opt$outDir, opt$classVar, opt$runName)
  if (!dir.exists(opt$outDir)) dir.create(opt$outDir, opt$classVar, recursive=TRUE)
  cat(glue::glue("[{par$prgmTag}]: outDir={opt$outDir}.{RET}") )
  
  search_fn_key <- '.model-files.csv.gz$'
  full_ss_csvs <- list.files(opt$modelDir, pattern=search_fn_key, full.names=TRUE, recursive=TRUE)
  full_ss_cnts <- length(full_ss_csvs)
  cat(glue::glue("[{par$prgmTag}]: full_ss_cnts={full_ss_cnts}.{RET}") )
  
  stopifnot(full_ss_cnts>0)
  
  prd_cnt <- 0
  for (mIdx in c(1:full_ss_cnts)) {
    
    # Gather all files::
    fns_csv <- full_ss_csvs[mIdx]
    stopifnot(file.exists(fns_csv))
    dir_path <- base::dirname(fns_csv)
    
    sam_csv <- NULL
    par_csv <- NULL
    fet_csv <- NULL
    mod_rds <- NULL
    
    fns_tib <- suppressMessages(suppressWarnings( readr::read_csv(fns_csv) ))
    sam_csv <- file.path(dir_path, fns_tib %>% dplyr::filter(Type=="SampleSheet") %>% head(n=1) %>% dplyr::pull(File_Name))
    par_csv <- file.path(dir_path, fns_tib %>% dplyr::filter(Type=="Params") %>% head(n=1) %>% dplyr::pull(File_Name))
    fet_csv <- file.path(dir_path, fns_tib %>% dplyr::filter(Type=="Features") %>% head(n=1) %>% dplyr::pull(File_Name))
    mod_rds <- file.path(dir_path, fns_tib %>% dplyr::filter(Type=="Model") %>% head(n=1) %>% dplyr::pull(File_Name))
    
    # Ensure all data exists...
    #
    if (is.null(fns_tib) || is.null(sam_csv) || is.null(par_csv) || is.null(fet_csv) || is.null(mod_rds) ) {
      cat(glue::glue("[{par$prgmTag}]:{TAB} Unable to find mod_rds file; dir_path={dir_path}; fns_tib={fns_csv}. Skipping...{RET}") )
      next
    }
    if (is.null(mod_rds) || length(mod_rds)==0 || !file.exists(mod_rds)) {
      cat(glue::glue("[{par$prgmTag}]:{TAB} Unable to find mod_rds file; dir_path={dir_path}; fns_tib={fns_csv}. Skipping...{RET}") )
      next
    }
    stopifnot(file.exists(sam_csv))
    stopifnot(file.exists(par_csv))
    stopifnot(file.exists(par_csv))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Add Cross Validation Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    #
    # TBD:: 
    #  - Incoporate Cross Validation Files
    #  - Split this script into::
    #    - predict_models.R
    #    - load_matrix
    #    - predict_model()
    #
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                         Build Shell Scripts::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
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
      if (!dir.exists(mergeDir))
        stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: mergeDir={mergeDir} does not exist! Skipping...{RET}{RET}"))
      # merge_name <- base::basename(mergeDir)
    }
    # cur_dir <- file.path(opt$outDir,user_dir,fets_dir,seed_dir,orig_dir,train_dir,merge_name)
    
    cur_dir <- file.path(opt$outDir,user_dir,fets_dir,seed_dir,orig_dir,train_dir)
    if (!dir.exists(cur_dir)) dir.create(cur_dir, recursive=TRUE)
    
    sh_path <- file.path(cur_dir, "test.model.sh")
    rm_vec  <- c("modelDir")
    add_tib <- tibble::tibble(Option=c("outDir","sampleSheet","params","features","model"),
                              Value=c(cur_dir,sam_csv,par_csv,fet_csv,mod_rds))
    
    run_sh <- optsToCommand(opts=opt_tib, pre=opt$Rscript, exe=par$exePath, rm=rm_vec, add=add_tib, file=sh_path,
                            verbose=opt$verbose,vt=1,tc=1,tt=NULL)
    
    run_id <- paste0('prd-',prd_cnt,'-cl')
    cmd <- paste(opt$lanExe,run_id,run_sh, sep=' ')
    if (is.null(opt$lanExe) || stringr::str_length(opt$lanExe)==0) cmd <- run_sh
    
    cat(glue::glue("[{par$prgmTag}]:{TAB}. Launching[{prd_cnt}]: cmd={cmd}...{RET}{RET}") )
    if (opt$execute) sys_ret_val <- base::system(cmd)
    
    if (!sys_ret_val)
      cat(glue::glue("[{par$prgmTag}]: Warning: Bad System Return[{prd_cnt}]={sys_ret_val}; cmd='{cmd}'{RET}{RET}"))
    
    prd_cnt <- prd_cnt + 1
    
    # if (opt$single) break
  }
  cat(glue::glue("[{par$prgmTag}]: Done. Launching in cluster mode.{RET}{RET}"))
  
} else {
  cat(glue::glue("[{par$prgmTag}]: Launching in single-job mode...{RET}"))
  
  if (!dir.exists(opt$outDir)) dir.create(opt$outDir, opt$classVar, recursive=TRUE)
  cat(glue::glue("[{par$prgmTag}]: outDir={opt$outDir}.{RET}") )
  
  if (FALSE) {
    # Temp Fix...
    tmp_dir <- '/Users/bbarnes/Documents/Projects/methylation/scratch/tmp'
    tmp_csv <- file.path(tmp_dir, 'tmp_opt_tib.csv')
    readr::write_csv(opt_tib, tmp_csv)
    
    tmp_dir <- '/Users/bbarnes/Documents/Projects/methylation/scratch/tmp'
    tmp_csv <- file.path(tmp_dir, 'tmp_opt_tib.csv')
    
    # Swap for local testing...
    opt_tib1 <- opt_tib
    opt_tib  <- readr::read_csv(tmp_csv)
    opt_tib %>% print(n=22)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #               Preprocessing:: Loading Models/Target CGNs
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  sam_csv <- opt_tib %>% dplyr::filter(Option=='sampleSheet') %>% head(n=1) %>% dplyr::pull(Value)
  par_csv <- opt_tib %>% dplyr::filter(Option=='params') %>% head(n=1) %>% dplyr::pull(Value)
  fet_csv <- opt_tib %>% dplyr::filter(Option=='features') %>% head(n=1) %>% dplyr::pull(Value)
  mod_rds <- opt_tib %>% dplyr::filter(Option=='model') %>% head(n=1) %>% dplyr::pull(Value)
  
  stopifnot(file.exists(sam_csv))
  stopifnot(file.exists(par_csv))
  stopifnot(file.exists(fet_csv))
  stopifnot(file.exists(mod_rds))
  
  sam_tib <- suppressMessages(suppressWarnings( readr::read_csv(sam_csv) ))
  par_tib <- suppressMessages(suppressWarnings( readr::read_csv(par_csv) ))
  fet_tib <- suppressMessages(suppressWarnings( readr::read_csv(fet_csv) ))
  cur_mod <- readr::read_rds(mod_rds)
  
  modName <- par_tib %>% dplyr::filter(Option=='model') %>% head(n=1) %>% dplyr::pull(Value)
  modLamb <- par_tib %>% dplyr::filter(Option=='lambda') %>% head(n=1) %>% dplyr::pull(Value)
  
  modText <- modName
  if (!is.null(modLamb) && length(modLamb)!=0) modText <- paste(modText,modLamb, sep='-')
  
  if (modName!='rforest' && modName!='glmnet')
    stop(glue::glue("[{par$prgmTag}]: ERROR; unsupported model type={modName}; Terminating...{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Preprocessing:: Filtering Test Data
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # opt$single <- FALSE
  opt$single <- TRUE
  
  all_sam_csv <- file.path(opt$outDir, 'combined_performance_samples.csv.gz')
  all_sum_csv <- file.path(opt$outDir, 'combined_performance_summary.csv.gz')
  
  all_sam_tib <- NULL
  all_sum_tib <- NULL
  
  for (betaKey in lociBetaKey_vec) {
    for (pvalKey in lociPvalKey_vec) {
      for (pvalMin in lociPvalMin_vec) {
        
        betaStr <- betaKey %>% stringr::str_replace_all('_', '-')
        pvalStr <- paste(pvalKey %>% stringr::str_replace_all('_', '-'), pvalMin, sep='-')
        dirName <- paste(betaStr,pvalStr,modText, sep='_')
        outName <- paste(opt$classVar, opt$runName, modText, dirName, sep='_')
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                     Build Current Output Directory::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        cur_opt_dir <- file.path(opt$outDir, dirName, modText)
        if (!dir.exists(cur_opt_dir)) dir.create(cur_opt_dir, recursive=TRUE)
        cat(glue::glue("[{par$prgmTag}]: Built; cur_opt_dir={cur_opt_dir}!{RET}") )
        
        if (opt$clean) unlink(list.files(cur_opt_dir, full.names=TRUE))
        
        cTracker <- timeTracker$new(verbose=opt$verbose)
        
        # Defined Output files::
        cur_sam_csv     <- file.path(cur_opt_dir, paste(outName,'method_performance_samples.csv.gz', sep='.') )
        cur_sum_csv     <- file.path(cur_opt_dir, paste(outName,'method_performance_summary.csv.gz', sep='.') )
        class_ss_csv    <- file.path(cur_opt_dir, paste(outName,'ClasSampleSheet.sorted.csv.gz', sep='.') )
        beta_masked_rds <- file.path(cur_opt_dir, paste(outName,'beta_masked_mat.rds', sep='.') )
        index_masks_csv <- file.path(cur_opt_dir, paste(outName,'beta_masked_idx.csv.gz', sep='.') )

        opt$clean <- FALSE
        opt$clean <- TRUE
        beta_file_tib <- getCallsMatrixFiles(
          betaKey=betaKey,pvalKey=pvalKey,pvalMin=pvalMin, dirs=mergeDirs_vec, cgn=fet_tib, classes=opt$trainClass,
          class_var=class_var, class_idx=class_idx, pval_name=NULL, pval_perc=NULL,
          clean=opt$clean, beta_rds=beta_masked_rds, ss_csv=class_ss_csv, mask_csv=index_masks_csv,
          sam_suffix="_AutoSampleSheet.csv.gz$", dat_suffix="_MergedDataFiles.tib.csv.gz", sentrix_name="Sentrix_Name",
          verbose=opt$verbose, vt=3,tc=1,tt=cTracker)
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                   Build Raw and Imputed Sorted Matricies::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        # Now load previous results from file::
        sampleSheet_tib <- loadFromFileTib(tib=beta_file_tib, type="SampleSheet")
        index_masks_tib <- loadFromFileTib(tib=beta_file_tib, type="Mask")
        beta_impute_mat <- loadFromFileTib(tib=beta_file_tib, type="Beta")
        
        # Rebuild NA beta matrix for DML/dBL calculations::
        pval_na_idx_vec <- index_masks_tib %>% dplyr::arrange(idx) %>% dplyr::distinct(idx) %>% dplyr::pull(idx) %>% as.vector()
        beta_masked_mat <- beta_impute_mat
        if (!is.null(pval_na_idx_vec) && length(pval_na_idx_vec)!=0)
          beta_masked_mat[ pval_na_idx_vec ] <- NA
        
        labs_idx_vec <- sampleSheet_tib %>% dplyr::pull(!!class_idx) %>% as.vector()
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                            Make Predictions::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        if (modName=='glmnet') {
          
          cur_pred = predGlmnet(mod=cur_mod, data=t(beta_masked_mat), labs=labs_idx_vec, 
                                name=modName, lambda="lambda.1se", type=type.measure,
                                verbose=opt$verbose,vt=1,tt=pTracker) %>% dplyr::mutate(Group=modText)
          
        } else if (modName=='rforest') {
          
          cur_pred = predRandomForest(mod=cur_mod, data=t(beta_impute_mat), labs=labs_idx_vec, 
                                      name=modName, # lambda="lambda.1se", type=type.measure,
                                      verbose=opt$verbose,vt=1,tt=pTracker) %>% dplyr::mutate(Group=modText)

        } else {
          stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported modName={modName}!!!{RET}{RET}"))
        }
        cur_sam_tib <- predToCalls(pred=cur_pred, labs=labs_idx_vec, pred_lab="Pred_Class",
                                   verbose=opt$verbose,vt=1,tt=pTracker)
        
        cur_sum_tib <- callToSumTib(call=cur_sam_tib, name_lab=modText, true_lab="True_Class",call_lab="Call",
                                    verbose=opt$verbose,vt=1,tt=pTracker)

        # Write current results to local directory::
        readr::write_csv(cur_sam_tib, cur_sam_csv)
        readr::write_csv(cur_sum_tib, cur_sum_csv)
        
        # Add additional variables::
        cur_sam_tib <- cur_sam_tib %>% dplyr::mutate(TestBeta=betaKey, TestPval=pvalKey, TestPvalMin=pvalMin)
        cur_sum_tib <- cur_sum_tib %>% dplyr::mutate(TestBeta=betaKey, TestPval=pvalKey, TestPvalMin=pvalMin)
        
        # Add results to previous summaries::
        all_sam_tib <- all_sam_tib %>% dplyr::bind_rows(cur_sam_tib)
        all_sum_tib <- all_sum_tib %>% dplyr::bind_rows(cur_sum_tib)
        
        cat(glue::glue("[{par$prgmTag}]: Done. triplet=({betaKey},{pvalKey},{pvalMin}).{RET}{RET}"))
        
        if (opt$single) break
      }
      if (opt$single) break
    }
    if (opt$single) break
  }
  
  readr::write_csv(all_sam_tib, all_sam_csv)
  readr::write_csv(all_sum_tib, all_sum_csv)
  
  cat(glue::glue("[{par$prgmTag}]: Done. Launching in single-job mode.{RET}{RET}"))
}


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
