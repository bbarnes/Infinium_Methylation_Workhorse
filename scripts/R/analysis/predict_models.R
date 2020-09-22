
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
# suppressWarnings(suppressPackageStartupMessages(require("grid")) )

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
par$date    <- Sys.Date() %>% as.character()
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'analysis'
par$prgmTag <- 'predict_models'

# Illumina based directories::
par$macDir  <- '/Users/bbarnes/Documents/Projects/methylation/tools'
par$lixDir  <- '/illumina/scratch/darkmatter/Projects/COVIC'

# Directory Parameters::
opt$outDir    <- NULL
opt$mergeDir  <- NULL

# Run Parameters::
opt$runName   <- NULL

opt$build_dir_str <- 'build_models'

# Class Parameters::
# Really simple test to make sure we can seperate the sexes...
opt$classVar   <- NULL
opt$trainClass <- NULL

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
  
  opt$clean  <- FALSE
  opt$single <- TRUE
  
  opt$clean  <- TRUE
  opt$single <- FALSE
  opt$execute <-  FALSE
  
  if (opt$classVar=='Sample_Class') {
    version   <- "C0"
    platform  <- "EPIC"
    
    runNameA  <- "COVIC-Set1-15052020"
    runNameB  <- "COVIC-Set5-10062020"
    runNameC  <- "COVIC-Set7-06082020"
    
    runNameD  <- 'COVIC-Set-17'
    runNameE  <- 'COVIC-Set-15'
    runNameF  <- 'COVIC-Set-57'
    
    sesGroup  <- 'ind-beta_i-poob-1'
    dmlGroup  <- 'dml-1000-Ivana-145'
    seedGroup <- 'seed-13'
    seedGroup <- 'seed-42'

    opt$runName   <-  runNameA
    
    opt$modelDir <- paste(
      # file.path(par$topDir, 'build_models',platform,version,opt$classVar,opt$runName,'ind-beta_i-poob-1'),
      # file.path(par$topDir, 'build_models',platform,version,opt$classVar,opt$runName,sesGroup,dmlGroup,seedGroup),
      file.path(par$topDir, 'build_models',platform,version,opt$classVar,opt$runName,sesGroup),
      sep=',')
    # opt$modelDir <- '/Users/bbarnes/Documents/Projects/methylation/scratch/build_models/EPIC/C0/Sample_Class/COVIC-Set1-15052020/ind-beta_i-poob-1/dml-1000-Ivana-145/seed-42/alpha-1/rforest-1-lambda-min.model.rds'
    
    alphaGroup <- 'alpha-1'
    modelGroup <- 'rforest-1-lambda-min'
    modelFile <- paste(modelGroup,'model.rds', sep='.')
    # opt$modelDir <- file.path(par$topDir, 'build_models',platform,version,opt$classVar,opt$runName,sesGroup,dmlGroup,seedGroup,alphaGroup,modelFile)
    
    opt$mergeDir  <- paste(
      file.path(par$topDir, 'merge_builds', platform, version, opt$classVar, runNameB),
      sep=',')
    
    opt$trainClass <- paste('nSARSCov2', 'pSARSCov2', sep=',')
  }
  
  opt$outDir <- file.path(par$topDir)
  
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
    make_option(c("-m","--mergeDir"), type="character", default=opt$mergeDir, 
                help="List of Merged Swifthoof Build Directory(s), commas seperated [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("-r","--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--build_dir_str"), type="character", default=opt$build_dir_str, 
                help="Assumed build_model sub-directory name. Should not be changed. [default= %default]", metavar="character"),
    
    # To run many models
    make_option(c("--modelDir"), type="character", default=opt$modelDir, 
                help="List of Model Build Directory(s). This will launch many models. [default= %default]", metavar="character"),
    
    # Class Parameters::
    make_option(c("--classVar"), type="character", default=opt$classVar, 
                help="Classification Variable Name [default= %default]", metavar="character"),
    make_option(c("--trainClass"), type="character", default=opt$trainClass, 
                help="Training Class Variable Names, comma seperated [default= %default]", metavar="character"),

    # Loci Level Filtering Parameters::
    make_option(c("--lociBetaKey"), type="character", default=opt$lociBetaKey, 
                help="Loci Beta-Method Name (key) for training. Comma seperated list. [default= %default]", metavar="character"),
    make_option(c("--lociPvalKey"), type="character", default=opt$lociPvalKey, 
                help="Loci Pval-Method Name (key) for filtering for training. Comma seperated list. [default= %default]", metavar="character"),
    make_option(c("--lociPvalMin"), type="character", default=opt$lociPvalMin,
                help="Pval Min Passing for loci filtering for training. Comma seperated list. [default= %default]", metavar="character"),
    
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

if (is.null(par$runMode) || is.null(par$prgmDir) || is.null(par$prgmTag) || is.null(par$scrDir) || is.null(par$datDir)) {
  
  par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
  par_tib %>% base::print(n=base::nrow(par_tib) )
  
  if (is.null(par$runMode)) cat(glue::glue("[Usage]: runMode is NULL!!!{RET}"))
  if (is.null(par$prgmDir)) cat(glue::glue("[Usage]: prgmDir is NULL!!!{RET}"))
  if (is.null(par$prgmTag)) cat(glue::glue("[Usage]: prgmTag is NULL!!!{RET}"))
  if (is.null(par$scrDir))  cat(glue::glue("[Usage]: scrDir is NULL!!!{RET}"))
  if (is.null(par$datDir))  cat(glue::glue("[Usage]: darDir is NULL!!!{RET}"))
  base::stop("Null Parameters!\n\n")
}

if (is.null(opt$outDir) || is.null(opt$mergeDir) || is.null(opt$modelDir) ||
    is.null(opt$runName) || 
    is.null(opt$classVar) || is.null(opt$trainClass) ||
    is.null(opt$lociBetaKey) || is.null(opt$lociPvalKey) || is.null(opt$lociPvalMin) || 
    is.null(opt$percisionBeta) || is.null(opt$percisionPval) || 
    is.null(opt$execute) || is.null(opt$single) || is.null(opt$parallel) || is.null(opt$cluster) ||
    
    is.null(opt$clean) || is.null(opt$Rscript) || is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )

  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  if (is.null(opt$outDir))     cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$mergeDir))   cat(glue::glue("[Usage]: mergeDir is NULL!!!{RET}"))
  if (is.null(opt$modelDir))   cat(glue::glue("[Usage]: modelDir is NULL!!!{RET}"))
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

par$man_src_dir <- file.path(par$scrDir, 'manifests/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$man_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$man_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$man_src_dir}!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Preprocessing:: System Params
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt <- setLaunchExe(opts=opt, pars=par, verbose=opt$verbose, vt=5,tc=0)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Preprocessing:: General Params
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

opt$outDir <- file.path(opt$outDir, par$prgmTag, opt$runName, par$date)
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Preprocessing:: Find all models files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

file_tibs <- NULL
if (file_test("-d", opt$modelDir)) {
  # Write run script without out --cluster
  cat(glue::glue("[{par$prgmTag}]: Launching in cluster mode...{RET}"))
  
  rds_pattern=".model.rds$"
  full_list <- list.files(opt$modelDir, pattern=rds_pattern, full.names=TRUE, recursive=TRUE)
  full_cnts <- length(full_list)
  
  if (full_cnts==0) stop(glue::glue("[{par$prgmTag}]: ERROR: Failed to find any models: full_cnts={full_cnts}; modelDir={opt$modelDir}.{RET}") )
  cat(glue::glue("[{par$prgmTag}]: full_cnts={full_cnts}.{RET}") )
  
  for (mIdx in c(1:full_cnts)) {
    mod_rds  <- full_list[mIdx]
    file_tib <- checkModelFiles(mod_rds=mod_rds, verbose=opt$verbose, vt=4, tc=1, tt=NULL)
    
    if (is.null(file_tib)) {
      cat(glue::glue("[{par$prgmTag}]: ERROR: Failed checkModelFiles() for mod_rds={mod_rds}. Skipping...{RET}") )
      next
    }
    file_tibs <- file_tibs %>% dplyr::bind_rows(file_tib)
    
    if (opt$single) break
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Build Union of all Features::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  all_fet_tib <- lapply( file_tibs$Model_Features, function(x) { suppressMessages(suppressWarnings( readr::read_csv(x) )) %>% dplyr::select(1) }  ) %>% 
    dplyr::bind_rows() %>% dplyr::group_by_all() %>% dplyr::summarise(Count=n()) %>% dplyr::arrange(Probe_ID) # %>% dplyr::arrange(-Count)
  # all_fet_tib %>% dplyr::arrange(-Count) %>% as.data.frame()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Build Minimum Test Output Directory::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  uniq_dir <- file_tibs$Dir %>% unique()
  uniq_cnt <- uniq_dir %>% length()
  unia_max_cnt <- uniq_cnt
  if (opt$verbose>=4) {
    cat(glue::glue("[{par$prgmTag}]: uniq_cnt={uniq_cnt}.{RET}"))
    print(uniq_dir)
  }
  
  uniq_rnd_cnt <- 1
  while(uniq_cnt > 1) {
    uniq_dir <- uniq_dir %>% dirname() %>% unique()
    uniq_cnt <- uniq_dir %>% length()
    if (opt$verbose>=4) {
      cat(glue::glue("[{par$prgmTag}]: uniq_cnt={uniq_cnt}.{RET}"))
      print(uniq_dir)
    }
    
    if (uniq_rnd_cnt > unia_max_cnt) {
      cat(glue::glue("[{par$prgmTag}]: Warning: Reached reduction limit: uniq_cnt={uniq_cnt} > unia_max_cnt={unia_max_cnt}.{RET}"))
      break
    }
    uniq_rnd_cnt <- uniq_rnd_cnt + 1
  }
  
  uniq_out_dir <- NULL
  if (!is.null(opt$build_dir_str) && stringr::str_detect(uniq_dir, paste0('/',opt$build_dir_str,'/') )) {
    uniq_out_dir <- file.path(opt$outDir, stringr::str_remove(uniq_dir, paste0('.*/',opt$build_dir_str,'/') ) )
  } else {
    uniq_out_dir <- file.path(opt$outDir, opt$classVar,opt$runName,dirNameTrain,fetSizeTrain,seedStrTrain,alphaStrTrain,modStr)
  }
  if (!dir.exists(uniq_out_dir)) dir.create(uniq_out_dir, recursive=TRUE)
  cat(glue::glue("[{par$prgmTag}]: uniq_out_dir(cnt={uniq_cnt})={uniq_out_dir}.{RET}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Build Minimum Test Matrix::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  outName <- opt$runName
  uniq_fet_csv    <- file.path(uniq_out_dir, paste(outName,'Feature-histogram.csv.gz', sep='.') )
  class_ss_csv    <- file.path(uniq_out_dir, paste(outName,'ClasSampleSheet.sorted.csv.gz', sep='.') )
  beta_masked_rds <- file.path(uniq_out_dir, paste(outName,'beta_masked_mat.rds', sep='.') )
  index_masks_csv <- file.path(uniq_out_dir, paste(outName,'beta_masked_idx.csv.gz', sep='.') )
  
  #
  # TBD:: Cannot run this step without writing Detection Pvalue file as welll as Beta Matrix...
  #
  if (FALSE) {
    beta_file_tib <- getCallsMatrixFiles(
      betaKey=betaKey,pvalKey=pvalKey,pvalMin=NULL, dirs=mergeDirs_vec, cgn=all_fet_tib, classes=trainClass_vec,
      class_var=class_var, class_idx=class_idx, pval_name=NULL, pval_perc=NULL,
      clean=opt$clean, beta_rds=beta_masked_rds, ss_csv=class_ss_csv, mask_csv=index_masks_csv,
      sam_suffix="_AutoSampleSheet.csv.gz$", dat_suffix="_MergedDataFiles.tib.csv.gz", sentrix_name="Sentrix_Name",
      verbose=opt$verbose, vt=3,tc=1,tt=pTracker)
    
    # Pass File Names file 
  }
  readr::write_csv(all_fet_tib, uniq_fet_csv)
  
} else if (file_test("-f", opt$modelDir) ) {
  # Run module directly::
  cat(glue::glue("[{par$prgmTag}]: Launching in direct single-job...{RET}"))
  
  mod_rds  <- opt$modelDir
  file_tib <- checkModelFiles(mod_rds=mod_rds, verbose=opt$verbose, tc=1, tt=NULL)
  
  if (is.null(file_tib))
    stop(glue::glue("[{par$prgmTag}]: ERROR: Failed checkModelFiles() for mod_rds={mod_rds}. Skipping...{RET}") )
    
  file_tibs <- file_tibs %>% dplyr::bind_rows(file_tib)
  
} else {
  stop(glue::glue("[{par$prgmTag}]: ERROR: modelDir is neither a file or directory: {opt$modelDir}! Exiting...{RET}{RET}") )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Main:: Process Each Model & Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rIdx <- 1
mods_cnts <- base::nrow(file_tibs) %>% as.integer()
for (rIdx in c(1:mods_cnts)) {
  file_tib <- file_tibs[rIdx,]
  cat(glue::glue("[{par$prgmTag}]: Starting model[{rIdx}]={file_tib}...{RET}"))
  
  sam_tib <- suppressMessages(suppressWarnings( readr::read_csv(file_tib$Model_SampleSheet[1]) ))
  par_tib <- suppressMessages(suppressWarnings( readr::read_csv(file_tib$Model_Params[1]) ))
  fet_tib <- suppressMessages(suppressWarnings( readr::read_csv(file_tib$Model_Features[1]) ))
  cur_mod <- readr::read_rds(file_tib$Model_RDS[1])
  
  #
  # Pull Parameters::
  #
  betaKeyTrain <- par_tib %>% dplyr::filter(Option=='lociBetaKey') %>% head(n=1) %>% dplyr::pull(Value)
  pvalKeyTrain <- par_tib %>% dplyr::filter(Option=='lociPvalKey') %>% head(n=1) %>% dplyr::pull(Value)
  pvalMinTrain <- par_tib %>% dplyr::filter(Option=='lociPvalMin') %>% head(n=1) %>% dplyr::pull(Value)
  
  # dml-1000-Ivana-145
  dmlSizeTrain <- par_tib %>% dplyr::filter(Option=='featureSizeDml') %>% head(n=1) %>% dplyr::pull(Value)
  preSizeTrain <- par_tib %>% dplyr::filter(Option=='featureNamePre') %>% head(n=1) %>% dplyr::pull(Value)
  
  fetSizeTrain <- NULL
  if (!is.null(dmlSizeTrain) && length(dmlSizeTrain)!=0) fetSizeTrain <- paste('dml',dmlSizeTrain, sep='-')
  if (!is.null(preSizeTrain) && length(preSizeTrain)!=0) {
    if (is.null(fetSizeTrain)) {
      fetSizeTrain <- preSizeTrain
    } else {
      fetSizeTrain <- paste(fetSizeTrain,preSizeTrain, sep='-')
    }
  }
  
  # Beta/P-value Train Variables and Strings::
  betaStrTrain <- betaKeyTrain %>% stringr::str_replace_all('_', '-')
  pvalStrTrain <- paste(pvalKeyTrain %>% stringr::str_replace_all('_', '-'), pvalMinTrain, sep='-')
  dirNameTrain <- paste(betaStrTrain,pvalStrTrain, sep='_')
  
  # Seed Train Variables and Strings::
  seedValTrain <- par_tib %>% dplyr::filter(Option=='seed') %>% head(n=1) %>% dplyr::pull(Value)
  seedStrTrain <- paste('seed',seedValTrain, sep='-')
  
  # Alpha Train Variables and Strings::
  alphaValTrain <- par_tib %>% dplyr::filter(Option=='alpha') %>% head(n=1) %>% dplyr::pull(Value)
  alphaStrTrain <- paste('alpha',alphaValTrain, sep='-')
  
  # Model Train Variables and Strings::
  modKey <- NULL
  modLab <- NULL
  modStr <- NULL
  modKey <- par_tib %>% dplyr::filter(Option=='model') %>% head(n=1) %>% dplyr::pull(Value)
  modLab <- par_tib %>% dplyr::filter(Option=='lambda') %>% head(n=1) %>% dplyr::pull(Value)
  modStr <- modKey
  if (!is.null(modLab) && length(modLab)!=0) modStr <- paste(modStr,modLab, sep='-')
  
  # if (modKey=='glmnet') next
  if (modKey!='rforest' && modKey!='glmnet')
    stop(glue::glue("[{par$prgmTag}]: ERROR; unsupported model type={modKey}; Terminating...{RET}"))

  #
  # Define additional annotation tib
  #
  add_ano_tib <- tibble::tibble(betaKeyTrain=betaKeyTrain, pvalKeyTrain=pvalKeyTrain, pvalMinTrain=pvalMinTrain,
                                dmlSizeTrain=dmlSizeTrain, preSizeTrain=preSizeTrain, seedValTrain=seedValTrain, 
                                alphaValTrain=alphaValTrain)
  
  #
  # Redefine Output directory based on search directory::
  #
  cur_out_dir <- NULL
  if (!is.null(opt$build_dir_str) && stringr::str_detect(file_tib$Dir, paste0('/',opt$build_dir_str,'/') )) {
    cur_out_dir <- file.path(opt$outDir, stringr::str_remove(file_tib$Dir, paste0('.*/',opt$build_dir_str,'/') ), modStr )
    # file.path(opt$outDir, stringr::str_remove(uniq_dir, paste0('.*/',opt$build_dir_str,'/') ) )
  } else {
    cur_out_dir <- file.path(opt$outDir, opt$classVar,opt$runName,dirNameTrain,fetSizeTrain,seedStrTrain,alphaStrTrain,modStr)
  }
  if (!dir.exists(cur_out_dir)) dir.create(cur_out_dir, recursive=TRUE)
  cat(glue::glue("[{par$prgmTag}]: cur_out_dir={cur_out_dir}.{RET}") )

  full_sum_tib <- NULL
  for (betaKeyTests in lociBetaKey_vec) {
    for (pvalKeyTests in lociPvalKey_vec) {
      for (pvalMinTests in lociPvalMin_vec) {
        
        try_str <- ''
        rdat = tryCatch({
          try_str <- 'Pass'
          # opt$verbose <- 30
          # rdat <- predictModelWrapper(cur_mod=cur_mod,
          predictModelWrapper(cur_mod=cur_mod,
                              betaKey=betaKeyTests, pvalKey=pvalKeyTests, pvalMin=pvalMinTests, ann=add_ano_tib,
                              dir=cur_out_dir, runName=opt$runName, modKey=modKey, modLab=modLab,
                              tests=mergeDirs_vec, cgn=fet_tib, classes=trainClass_vec,
                              classVar=opt$classVar, classIdx=class_idx, pvalName=NULL, pvalPerc=NULL,
                              clean=opt$clean,
                              sam_suffix="_AutoSampleSheet.csv.gz$", dat_suffix="_MergedDataFiles.tib.csv.gz", sentrix_name="Sentrix_Name",
                              verbose=opt$verbose, vt=3,tc=1,tt=pTracker)
        }, warning = function(w) {
          try_str <- paste('warning',par$prgmTag, sep='-')
          rdat <- NA
        }, error = function(e) {
          try_str <- paste('error',par$prgmTag, sep='-')
          rdat <- NA
        }, finally = {
          try_str <- paste('cleanup',par$prgmTag, sep='-')
          rdat <- NA
        })
        
        if (!is.null(rdat)) {
          full_sum_tib <- dplyr::bind_rows(full_sum_tib, rdat)
          cat(glue::glue("[{par$prgmTag}]: parallelFunc={par$prgmTag}: Succesfully tested model: try_str={try_str}.{RET}"))
        } else {
          cat(glue::glue("[{par$prgmTag}]: parallelFunc={par$prgmTag}: ERROR: Model Failed; try_str={try_str}.{RET}"))
        }
        
        if (opt$single) break
      }
      if (opt$single) break
    }
    if (opt$single) break
  }
  cat(glue::glue("[{par$prgmTag}]: Done. model[{rIdx}]={file_tib}.{RET}{RET}"))
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
