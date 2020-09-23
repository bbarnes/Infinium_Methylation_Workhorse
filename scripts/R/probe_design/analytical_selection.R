
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Build COVIC Model From Merged Builds::
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

# Program Parameters::
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'probe_design'
par$prgmTag <- 'analytical_selection'

# Illumina based directories::
par$macDir  <- '/Users/bbarnes/Documents/Projects/methylation/tools'
par$lixDir  <- '/illumina/scratch/darkmatter/Projects/COVIC'

# Directory Parameters::
opt$outDir <- NULL
opt$titDbl <- NULL
opt$repDbl <- NULL

# Run Parameters::
opt$runName   <- NULL

# Class Parameters::
opt$classVar <- NULL

# Manifest File Parameter::
opt$manifest <- NULL

# Sample Level Filtering Parameters::
opt$samplePvalName <- NULL
opt$samplePvalPerc <- NULL

# Loci Level Filtering Parameters::
opt$lociBetaKey <- NULL
opt$lociPvalKey <- NULL
opt$lociPvalMin <- NULL

# Screening Parameters::
opt$huCSS <- 0.2
opt$muCSS <- 0.5
opt$mhCSS <- 0.1

opt$trainClass <- NULL
opt$cross_perc_min <- 80

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
  
  # Pre-defined local options runTypes::
  par$local_runType <- 'mm10'
  
  if (par$local_runType=='mm10') {
    par$topDir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/scratch'

    opt$single   <- TRUE
    opt$single   <- FALSE
    opt$cluster  <- FALSE
    opt$parallel <- FALSE
    
    # Class Parameters::
    opt$classVar <- 'Sample_Name'
    
    par$platform <- 'LEGX'
    par$version  <- 'S1'
    
    opt$lociBetaKey <- "i-beta,ind-beta"
    opt$lociPvalKey <- "i-poob"
    # opt$lociPvalMin <- "1"
    # opt$lociPvalMin <- "0.02,0.1,1"
    opt$lociPvalMin <- "0.02,0.1,0.5,0.9,1"

    # Run Parameters::
    opt$runName   <- 'mm10-ILS-VAI'
    
    par$titName  <- 'mm10-ILS-VAI.Titration'
    # par$trainClass <- paste('T00DZ','T50DZ','T99DZ', sep=',')
    par$repName  <- 'mm10-ILS-VAI.Replicate'
    # par$trainClass <- paste('RepAC','RepS3','RepM1','RepSA', sep=',')
    
    # Directory Parameters::
    par$locDatDir <- file.path(par$topDir, 'build_models',par$platform,par$version,opt$classVar)
    
    # For a single unique pair::
    #
    # par$locParKey <- paste(opt$lociBetaKey, paste(opt$lociPvalKey,opt$lociPvalMin, sep='-'), sep='_')
    # par$locTitDbl <- paste0(opt$classVar,'_',opt$runName,'.Titration_',par$locParKey,'.full-dbl.csv.gz')
    # par$locRepDbl <- paste0(opt$classVar,'_',opt$runName,'.Replicate_',par$locParKey,'.full-dbl.csv.gz')
    # 
    # opt$titDbl <- file.path(par$locDatDir, par$titName, par$locParKey, par$locTitDbl)
    # opt$repDbl <- file.path(par$locDatDir, par$repName, par$locParKey, par$locRepDbl)
    
    opt$titDbl <- file.path(par$locDatDir, par$titName)
    opt$repDbl <- file.path(par$locDatDir, par$repName)

    # Screening Parameters::
    opt$huCSS <- 0.2
    opt$muCSS <- 0.5
    opt$mhCSS <- 0.1
    
    # Manifest File Parameter::
    opt$manifest <- file.path(par$datDir, 'manifest/base/LEGX-S1.manifest.sesame-base.cpg-sorted.csv.gz')

  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$outDir <- file.path(par$topDir, par$prgmTag, par$platform, par$version)
  
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
    make_option(c("-m","--titDbl"), type="character", default=opt$titDbl, 
                help="List of titration DBLs, commas seperated [default= %default]", metavar="character"),
    make_option(c("-m","--repDbl"), type="character", default=opt$repDbl, 
                help="List of replicate DBLs, commas seperated [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    
    # Manifest File Parameter::
    make_option(c("--manifest"), type="character", default=opt$manifest, 
                help="Manifest (Sesame) File [default= %default]", metavar="character"),
    
    # Class Parameters::
    make_option(c("--classVar"), type="character", default=opt$classVar, 
                help="Classification Variable Name [default= %default]", metavar="character"),

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
    
    # Screening Parameters::
    make_option(c("--huCSS"), type="double", default=opt$huCSS,
                help="Minimum H-U CSS [default= %default]", metavar="double"),
    make_option(c("--muCSS"), type="double", default=opt$muCSS,
                help="Minimum M-U CSS [default= %default]", metavar="double"),
    make_option(c("--muCSS"), type="double", default=opt$muCSS,
                help="Minimum M-U CSS [default= %default]", metavar="double"),
    
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

if (is.null(opt$outDir) || is.null(opt$titDbl) || is.null(opt$repDbl) || 
    is.null(opt$runName) || 
    is.null(opt$classVar) || 
    is.null(opt$huCSS) || is.null(opt$muCSS) || is.null(opt$mhCSS) ||
    is.null(opt$manifest) || 
    is.null(opt$percisionBeta) || is.null(opt$percisionPval) || 
    is.null(opt$execute) || is.null(opt$single) || is.null(opt$parallel) || is.null(opt$cluster) ||
    
    is.null(opt$clean) || is.null(opt$Rscript) || is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )

  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  if (is.null(opt$outDir))   cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$titDbl))   cat(glue::glue("[Usage]: titDbl is NULL!!!{RET}"))
  if (is.null(opt$repDbl))   cat(glue::glue("[Usage]: repDbl is NULL!!!{RET}"))
  if (is.null(opt$runName))  cat(glue::glue("[Usage]: runName is NULL!!!{RET}"))
  if (is.null(opt$classVar)) cat(glue::glue("[Usage]: classVar is NULL!!!{RET}"))
  
  if (is.null(opt$huCSS)) cat(glue::glue("[Usage]: huCSS is NULL!!!{RET}"))
  if (is.null(opt$muCSS)) cat(glue::glue("[Usage]: muCSS is NULL!!!{RET}"))
  if (is.null(opt$mhCSS)) cat(glue::glue("[Usage]: mhCSS is NULL!!!{RET}"))
  
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
#                     Preprocessing:: System Params
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$buildDbl <- FALSE
opt <- setLaunchExe(opts=opt, pars=par, verbose=opt$verbose, vt=5,tc=0)
if (!opt$isLinux && opt$buildDbl) {
  suppressWarnings(suppressPackageStartupMessages( require("Rcpp") ))
  
  par$sourceCpp <- file.path(par$scrDir, 'R/Rcpp/cpgLociVariation.cpp')
  if (!file.exists(par$sourceCpp)) par$sourceCpp <- file.path(par$scrDir, 'Rcpp/cpgLociVariation.cpp')
  if (!file.exists(par$sourceCpp)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceCpp} does not exist!{RET}"))
  Rcpp::sourceCpp(par$sourceCpp)
  
  cat(glue::glue("[{par$prgmTag}]: Loading Source Files form sourceCpp={par$sourceCpp}!{RET}{RET}") )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Preprocessing:: General Params
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

# mergeDirs_vec   <- opt$mergeDir %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
lociBetaKey_vec <- opt$lociBetaKey %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
lociPvalKey_vec <- opt$lociPvalKey %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
lociPvalMin_vec <- opt$lociPvalMin %>% str_split(pattern=',', simplify=TRUE) %>% as.double() %>% as.vector()

class_var <- rlang::sym(opt$classVar)
class_idx <- rlang::sym("Class_Idx")

opt$outDir <- file.path(opt$outDir, opt$classVar, opt$runName)
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)

cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Load Sesame Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ses_man_tib <- readr::read_csv(opt$manifest)
cgn_man_tib <- ses_man_tib %>% dplyr::filter(Probe_Type=='cg')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Main Scratch Space::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Options to be added::
#
opt$dB_min <- 0.2
all_sum_tib <- NULL

for (betaKey in lociBetaKey_vec) {
  for (pvalKey in lociPvalKey_vec) {
    for (pvalMin in lociPvalMin_vec) {
      
      cat(glue::glue("[{par$prgmTag}]: Checking: betaKey={betaKey}, pvalKey={pvalKey}, pvalMin={pvalMin}.{RET}"))
      
      
      locParKey <- paste(betaKey, paste(pvalKey,pvalMin, sep='-'), sep='_')
      
      tit_csv <- NULL
      rep_csv <- NULL
      
      tit_csv <- list.files(file.path(opt$titDbl,locParKey), pattern='.full-dbl.csv.gz$', full.names=TRUE) %>% head(n=1)
      rep_csv <- list.files(file.path(opt$repDbl,locParKey), pattern='.full-dbl.csv.gz$', full.names=TRUE) %>% head(n=1)
      
      if (is.null(tit_csv) || length(tit_csv)==0) next
      if (is.null(rep_csv) || length(rep_csv)==0) next
      
      cat(glue::glue("[{par$prgmTag}]: Running: betaKey={betaKey}, pvalKey={pvalKey}, pvalMin={pvalMin}.{RET}"))
      
      tit_tib <- suppressMessages(suppressWarnings( readr::read_csv(tit_csv) ))
      rep_tib <- suppressMessages(suppressWarnings( readr::read_csv(rep_csv) ))
      
      ses_tit_tib <- cgn_man_tib %>% dplyr::inner_join(tit_tib, by="Probe_ID")
      ses_rep_tib <- cgn_man_tib %>% dplyr::inner_join(rep_tib, by="Probe_ID")
      
      # Quick Summary::
      ses_tit_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(PT_Count=n()) %>% print()
      ses_rep_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(PT_Count=n()) %>% print()
      
      #
      # Titration Calling::
      #
      tit_call_tib <- ses_tit_tib %>% dplyr::mutate(
        HU_Call_mu=dplyr::case_when( T00DZ_T50DZ_CSS_mu > opt$huCSS ~ 'A', T00DZ_T50DZ_CSS_mu <= opt$huCSS ~ 'F', TRUE ~ 'N' ),
        MU_Call_mu=dplyr::case_when( T00DZ_T99DZ_CSS_mu > opt$muCSS ~ 'A', T00DZ_T99DZ_CSS_mu <= opt$muCSS ~ 'F', TRUE ~ 'N' ),
        MH_Call_mu=dplyr::case_when( T50DZ_T99DZ_CSS_mu > opt$mhCSS ~ 'A', T50DZ_T99DZ_CSS_mu <= opt$mhCSS ~ 'F', TRUE ~ 'N' )
      )
      
      #
      # Replicate Calling::
      #
      # Get Row stats:: dB_mu
      end_key <- '_beta_sd'
      end_key <- '_dB_q2'
      end_key <- '_dB_mu'
      
      rep_db_mu_tib <- ses_rep_tib %>% dplyr::select(Probe_ID,ends_with(end_key))
      rep_db_mu_mat <- rep_db_mu_tib %>% column_to_rownames(var='Probe_ID') %>% as.matrix()
      
      rep_call_tib <- rep_db_mu_tib %>% 
        dplyr::mutate(
          dB_min_mu=matrixStats::rowMins(rep_db_mu_mat,na.rm=TRUE),
          dB_max_mu=matrixStats::rowMaxs(rep_db_mu_mat,na.rm=TRUE) ) %>% 
        dplyr::mutate( dB_max_Call=dplyr::case_when( dB_max_mu >  opt$dB_min ~ 'F', dB_max_mu <= opt$dB_min ~ 'A', TRUE ~ 'N' ) )

      #
      # Join Manifest with Calls::
      #
      full_call_tib <- cgn_man_tib %>% 
        dplyr::left_join(tit_call_tib %>% dplyr::select(Probe_ID,HU_Call_mu:MH_Call_mu), by="Probe_ID") %>%
        dplyr::left_join(rep_call_tib %>% dplyr::select(Probe_ID,dB_max_Call), by="Probe_ID")
      
      full_tot_cnt <- full_call_tib %>% base::nrow()
      full_sum_tib <- full_call_tib %>% dplyr::select(HU_Call_mu:dB_max_Call) %>% 
        dplyr::mutate( across(everything(), ~replace_na(.x, 'N')) ) %>% dplyr::group_by_all() %>% 
        dplyr::summarise(Call_Cnt=n(),
                         Call_Per=round(100*Call_Cnt/full_tot_cnt,2)) %>% 
        dplyr::arrange(-Call_Cnt) %>% 
        dplyr::mutate(betaKey=betaKey,pvalKey=pvalKey,pvalMin=pvalMin) %>% 
        dplyr::ungroup() %>% dplyr::mutate(Rank=dplyr::row_number())
      
      all_sum_tib <- all_sum_tib %>% dplyr::bind_rows(full_sum_tib)
      
      cat(glue::glue("[{par$prgmTag}]: Done.{RET}{RET}"))
      
      if (opt$single) break
    }
    if (opt$single) break
  }
  if (opt$single) break
}











# Checking for missing data::
#
if (FALSE) {
  cgn_man_tib %>% base::nrow() %>% print()
  ses_tit_tib %>% base::nrow() %>% print()
  ses_rep_tib %>% base::nrow() %>% print()
  
  # Only when the p-value min is set to 1.0 do we have 100% overlap::
  ses_tit_tib %>% dplyr::anti_join(ses_rep_tib, by="Probe_ID") %>% dplyr::select(1:5)
  ses_rep_tib %>% dplyr::anti_join(ses_tit_tib, by="Probe_ID") %>% dplyr::select(1:5)
  
  ses_tit_tib %>% dplyr::filter(! Probe_ID %in% ses_rep_tib$Probe_ID) %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n())
  ses_rep_tib %>% dplyr::filter(! Probe_ID %in% ses_tit_tib$Probe_ID) %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n())
}
# tit_tib <- readr::read_csv(opt$titDbl) %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2) ) %>% split(.$Probe_Type)
# rep_tib <- readr::read_csv(opt$repDbl) %>% dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2) ) %>% split(.$Probe_Type)


#
# Titration Screening::
#

# tit_db_mu_tib <- tit_tib[['cg']] %>% dplyr::select(Probe_ID,ends_with('_CSS_mu'))
# tit_db_mu_mat <- tit_db_mu_tib %>% column_to_rownames(var='Probe_ID') %>% as.matrix()

tot_prb_cnt <- tit_tib[['cg']] %>% base::nrow()

tit_tib[['cg']] %>% dplyr::mutate(
  HU_Call_mu=dplyr::case_when( T00DZ_T50DZ_CSS_mu > opt$huCSS ~ 1, T00DZ_T50DZ_CSS_mu <= opt$huCSS ~ 0, TRUE ~ NA_real_ ),
  MU_Call_mu=dplyr::case_when( T00DZ_T99DZ_CSS_mu > opt$muCSS ~ 1, T00DZ_T99DZ_CSS_mu <= opt$muCSS ~ 0, TRUE ~ NA_real_ ),
  MH_Call_mu=dplyr::case_when( T50DZ_T99DZ_CSS_mu > opt$mhCSS ~ 1, T50DZ_T99DZ_CSS_mu <= opt$mhCSS ~ 0, TRUE ~ NA_real_ )
) %>% 
  dplyr::group_by(Probe_Type,HU_Call_mu,MU_Call_mu,MH_Call_mu) %>% 
  dplyr::summarise(Call_Counts=n(), Call_Perc=round(100*Call_Counts/tot_prb_cnt,1)) %>% dplyr::arrange(-Call_Counts)


#
# Replicate Screening::
#

# Get Row stats:: dB_mu
end_key <- '_dB_q2'
end_key <- '_dB_mu'

rep_db_mu_tib <- rep_tib[['cg']] %>% dplyr::select(Probe_ID,ends_with(end_key))
rep_db_mu_mat <- rep_db_mu_tib %>% column_to_rownames(var='Probe_ID') %>% as.matrix()

rep_db_mu_sum_tib <- rep_db_mu_tib %>% dplyr::mutate(
  dB_min_mu=matrixStats::rowMins(rep_db_mu_mat,na.rm=TRUE),
  dB_max_mu=matrixStats::rowMaxs(rep_db_mu_mat,na.rm=TRUE)
)

rep_db_mu_sum_tib %>% dplyr::summarise(Total_Cnt=n(),
                                       PassMin_Cnt=count(dB_min_mu<=opt$dB_min, na.rm=TRUE),
                                       PassMax_Cnt=count(dB_max_mu<=opt$dB_min, na.rm=TRUE),
                                       PassMin_Per=round(100*PassMin_Cnt/Total_Cnt,3),
                                       PassMax_Per=round(100*PassMax_Cnt/Total_Cnt,3)
)

# Beta SD::
end_key <- '_beta_sd'
rep_beta_sd_tib <- rep_tib[['cg']] %>% dplyr::select(Probe_ID,ends_with(end_key))
rep_beta_sd_mat <- rep_beta_sd_tib %>% column_to_rownames(var='Probe_ID') %>% as.matrix()

rep_beta_sd_sum_tib <- rep_beta_sd_tib %>% dplyr::mutate(
  betaSD_min_mu=matrixStats::rowMins(rep_beta_sd_mat,na.rm=TRUE),
  betaSD_max_mu=matrixStats::rowMaxs(rep_beta_sd_mat,na.rm=TRUE)
)

rep_beta_sd_sum_tib %>% dplyr::summarise(Total_Cnt=n(),
                                       PassMin_Cnt=count(betaSD_min_mu<=opt$dB_min, na.rm=TRUE),
                                       PassMax_Cnt=count(betaSD_max_mu<=opt$dB_min, na.rm=TRUE),
                                       PassMin_Per=round(100*PassMin_Cnt/Total_Cnt,3),
                                       PassMax_Per=round(100*PassMax_Cnt/Total_Cnt,3)
)


#
# Previous CODE::
#

if (FALSE) {
  #
  # Calculate Passing/Failing Probes for Titration::
  #
  tot_prb_cnt <- tit_tib[['cg']] %>% base::nrow()
  
  tit_tib[['cg']] %>% dplyr::mutate(
    HU_Call_mu=dplyr::case_when( T00DZ_T50DZ_CSS_mu > opt$huCSS ~ 1, T00DZ_T50DZ_CSS_mu <= opt$huCSS ~ 0, TRUE ~ NA_real_ ),
    MU_Call_mu=dplyr::case_when( T00DZ_T99DZ_CSS_mu > opt$muCSS ~ 1, T00DZ_T99DZ_CSS_mu <= opt$muCSS ~ 0, TRUE ~ NA_real_ ),
    MH_Call_mu=dplyr::case_when( T50DZ_T99DZ_CSS_mu > opt$mhCSS ~ 1, T50DZ_T99DZ_CSS_mu <= opt$mhCSS ~ 0, TRUE ~ NA_real_ )
  ) %>% 
    dplyr::group_by(Probe_Type,HU_Call_mu,MU_Call_mu,MH_Call_mu) %>% 
    dplyr::summarise(Call_Counts=n(), Call_Perc=round(100*Call_Counts/tot_prb_cnt,1)) %>% dplyr::arrange(-Call_Counts)
  
  
  #
  # Quick plotting::
  #
  
  # Beta Value Distributions::
  #
  ggplot2::ggplot(data=tit_tib[['cg']]) + 
    geom_density(aes(x=T00DZ_beta_mu, fill="T00DZ_beta_mu"), alpha=0.7, na.rm=TRUE ) +
    geom_density(aes(x=T50DZ_beta_mu, fill="T50DZ_beta_mu"), alpha=0.7, na.rm=TRUE ) +
    geom_density(aes(x=T99DZ_beta_mu, fill="T99DZ_beta_mu"), alpha=0.7, na.rm=TRUE )
  
  ggplot2::ggplot(data=tit_tib[['cg']]) + 
    geom_density(aes(x=T00DZ_beta_q2, fill="T00DZ_beta_q2"), alpha=0.7, na.rm=TRUE ) +
    geom_density(aes(x=T50DZ_beta_q2, fill="T50DZ_beta_q2"), alpha=0.7, na.rm=TRUE ) +
    geom_density(aes(x=T99DZ_beta_q2, fill="T99DZ_beta_q2"), alpha=0.7, na.rm=TRUE )
  
  # Cluster Seperation Distribution
  #
  ggplot2::ggplot(data=tit_tib[['cg']]) + 
    geom_density(aes(x=T00DZ_T50DZ_CSS_mu, fill="T00DZ_T50DZ_CSS_mu"), alpha=0.7, na.rm=TRUE ) +
    geom_density(aes(x=T00DZ_T99DZ_CSS_mu, fill="T00DZ_T99DZ_CSS_mu"), alpha=0.7, na.rm=TRUE ) +
    geom_density(aes(x=T50DZ_T99DZ_CSS_mu, fill="T50DZ_T99DZ_CSS_mu"), alpha=0.7, na.rm=TRUE )
  
  ggplot2::ggplot(data=tit_tib[['cg']]) + 
    geom_density(aes(x=T00DZ_T50DZ_CSS_q2, fill="T00DZ_T50DZ_CSS_q2"), alpha=0.7, na.rm=TRUE ) +
    geom_density(aes(x=T00DZ_T99DZ_CSS_q2, fill="T00DZ_T99DZ_CSS_q2"), alpha=0.7, na.rm=TRUE ) +
    geom_density(aes(x=T50DZ_T99DZ_CSS_q2, fill="T50DZ_T99DZ_CSS_q2"), alpha=0.7, na.rm=TRUE )
  
}










if (FALSE) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                  Preprocessing:: Load Feature Selection
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Should Load Manifest Here::
  #
  # cgn_pre_str <- NULL
  # if (!is.null(featuresCsv_vec) && length(featuresCsv_vec)!=0 ) {
  #   cgn_pre_str <- featuresCsv_vec %>% base::basename() %>% stringr::str_remove('.csv.gz$') %>% paste(collapse='-')
  #   cgn_pre_tib <- lapply(featuresCsv_vec, loadCgnFeatures, id="Probe_ID", verbose=opt$verbose) %>% 
  #     dplyr::bind_rows() %>% dplyr::distinct()
  #   cgn_pre_len <- cgn_pre_tib %>% base::nrow()
  #   
  #   cat(glue::glue("[{par$prgmTag}]: Done. Loading Third Party Features; cgn_pre_len={cgn_pre_len}, cgn_pre_str={cgn_pre_str}!{RET}{RET}") )
  # } else {
  #   cat(glue::glue("[{par$prgmTag}]: Done. Did NOT load any Thrid Party features!{RET}{RET}") )
  # }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Preprocessing:: File Identification
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  for (betaKey in lociBetaKey_vec) {
    for (pvalKey in lociPvalKey_vec) {
      for (pvalMin in lociPvalMin_vec) {
        
        #     break }
        #   break }
        # break }
        
        
        betaStr <- betaKey %>% stringr::str_replace_all('_', '-')
        pvalStr <- paste(pvalKey %>% stringr::str_replace_all('_', '-'), pvalMin, sep='-')
        dirName <- paste(betaStr,pvalStr, sep='_')
        outName <- paste(opt$classVar, opt$runName, dirName, sep='_')
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                     Build Current Output Directory::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        cur_opt_dir <- file.path(opt$outDir, dirName)
        if (!dir.exists(cur_opt_dir)) dir.create(cur_opt_dir, recursive=TRUE)
        cat(glue::glue("[{par$prgmTag}]: Built; cur_opt_dir={cur_opt_dir}!{RET}") )
        
        if (opt$clean) unlink(list.files(cur_opt_dir, full.names=TRUE))
        
        cTracker <- timeTracker$new(verbose=opt$verbose)
        
        # Defined Output files::
        beta_masked_rds <- file.path(cur_opt_dir, paste(outName,'beta_masked_mat.rds', sep='.') )
        index_masks_csv <- file.path(cur_opt_dir, paste(outName,'beta_masked_idx.csv.gz', sep='.') )
        class_ss_csv <- file.path(cur_opt_dir, paste(outName,'ClasSampleSheet.sorted.csv.gz', sep='.') )
        
        full_dml_csv <- file.path(cur_opt_dir, paste(outName,'full-dml.csv.gz', sep='.') )
        full_dbl_csv <- file.path(cur_opt_dir, paste(outName,'full-dbl.csv.gz', sep='.') )
        
        dml_opt_csv  <- file.path(cur_opt_dir, paste(outName,'program-options.csv', sep='.') )
        dml_par_csv  <- file.path(cur_opt_dir, paste(outName,'program-parameters.csv', sep='.') )
        dml_time_csv <- file.path(cur_opt_dir, paste(outName,'time-tracker.csv.gz', sep='.') )
        
        beta_file_tib <- getCallsMatrixFiles(
          betaKey=betaKey,pvalKey=pvalKey,pvalMin=pvalMin, dirs=mergeDirs_vec, cgn=NULL, classes=opt$trainClass,
          class_var=class_var, class_idx=class_idx, pval_name=opt$samplePvalName, pval_perc=opt$samplePvalPerc,
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
          if (!opt$isLinux) {
            
            dbl_beg_txt  <- paste(full_dbl_csv,'begTime.txt', sep='.')
            dbl_end_txt  <- paste(full_dbl_csv,'endTime.txt', sep='.')
            
            if (!opt$clean && file.exists(full_dbl_csv) &&
                file.exists(dbl_beg_txt) && file.exists(dbl_end_txt) &&
                file.mtime(dbl_beg_txt) < file.mtime(full_dbl_csv) && 
                file.mtime(dbl_end_txt) > file.mtime(full_dbl_csv) ) {
              
              cat(glue::glue("[{par$prgmTag}]: deltaBetas already existst full_dbl_csv={full_dbl_csv}. Skipping...{RET}") )
              
            } else {
              cat(glue::glue("[{par$prgmTag}]: Calculating deltaBetas...{RET}") )
              
              # Name by Real Classes::
              #
              sample_cnt_tib <- sampleSheet_tib %>% 
                dplyr::group_by(!!class_var) %>% 
                dplyr::summarise(Count=n()) %>% tibble::as_tibble() %>% dplyr::mutate(Class=as.character(!!class_var) )
              sample_key_vec <- sample_cnt_tib %>% dplyr::pull(Class)
              sample_cnt_vec <- sample_cnt_tib %>% dplyr::pull(Count)
              
              dbl_sortMu_key <- paste(paste(sample_key_vec, collapse='_'), 'CSS_mu', sep='_') %>% rlang::sym()
              dbl_sortQ2_key <- paste(paste(sample_key_vec, collapse='_'), 'CSS_q2', sep='_') %>% rlang::sym()
              dbl_sdNeg_key  <- paste(sample_key_vec[1], 'beta_sd', sep='_') %>% rlang::sym()
              dbl_sdPos_key  <- paste(sample_key_vec[2], 'beta_sd', sep='_') %>% rlang::sym()
              
              # Need to pick one:: [ full_dblMu_tib or full_dblQ2_tib ]
              #
              # full_dblMu_tib %>% names() %>% length() =  65 nan
              # full_dblMu_tib %>% names() %>% length() =  80 cmb
              # full_dblMu_tib %>% names() %>% length() = 110 all
              #
              # full_dblMu_tib %>% dplyr::select(Probe_ID, dplyr::ends_with("_CSS_mu"))
              
              cpp.verbose <- 0
              
              # sample_key_vec <- sampleSheet_tib %>% dplyr::pull(Sample_Class) %>% unique()
              # sample_cnt_vec <- sampleSheet_tib %>% dplyr::group_by(Sample_Class) %>% dplyr::summarise(Count=n()) %>% dplyr::pull(Count)
              sample_key_vec <- sampleSheet_tib %>% dplyr::pull(!!class_var) %>% unique()
              sample_cnt_vec <- sampleSheet_tib %>% dplyr::group_by(!!class_var) %>% dplyr::summarise(Count=n()) %>% dplyr::pull(Count)
              
              full_dblMu_tib <- C_crossSampleLociRSquared(beta_masked_mat, sample_cnt_vec, sample_key_vec, cmb=TRUE, verbose=cpp.verbose) %>% 
                as.data.frame() %>% tibble::rownames_to_column(var='Probe_ID') %>% tibble::as_tibble()
              
              # Not used for anything yet::
              #  css_names <- full_dblMu_tib %>% dplyr::select(dplyr::ends_with("_CSS_mu")) %>% names()
              
              #
              # Calculate Passing/Failing Probes for Titration::
              #
              hu_min <- 0.2
              mu_min <- 0.5
              mh_min <- 0.1
              tot_prb_cnt <- full_dblMu_tib %>% base::nrow()
              
              full_dblMu_tib %>% dplyr::mutate(
                Probe_Type=stringr::str_sub(Probe_ID, 1,2),
                HU_Call_mu=dplyr::case_when( T00DZ_T50DZ_CSS_mu > hu_min ~ 1, T00DZ_T50DZ_CSS_mu <= hu_min ~ 0, TRUE ~ NA_real_ ),
                MU_Call_mu=dplyr::case_when( T00DZ_T99DZ_CSS_mu > mu_min ~ 1, T00DZ_T99DZ_CSS_mu <= mu_min ~ 0, TRUE ~ NA_real_ ),
                MH_Call_mu=dplyr::case_when( T50DZ_T99DZ_CSS_mu > mh_min ~ 1, T50DZ_T99DZ_CSS_mu <= mh_min ~ 0, TRUE ~ NA_real_ )
              ) %>% 
                dplyr::group_by(Probe_Type,HU_Call_mu,MU_Call_mu,MH_Call_mu) %>% 
                dplyr::summarise(Call_Counts=n(), Call_Perc=round(100*Call_Counts/PType_Count,1)) %>% dplyr::arrange(-Call_Counts)
              
              
              #
              # Quick plotting::
              #
              
              # Beta Value Distributions::
              #
              ggplot2::ggplot(data=full_dblMu_tib) + 
                geom_density(aes(x=T00DZ_beta_mu, fill="T00DZ_beta_mu"), alpha=0.7, na.rm=TRUE ) +
                geom_density(aes(x=T50DZ_beta_mu, fill="T50DZ_beta_mu"), alpha=0.7, na.rm=TRUE ) +
                geom_density(aes(x=T99DZ_beta_mu, fill="T99DZ_beta_mu"), alpha=0.7, na.rm=TRUE )
              
              ggplot2::ggplot(data=full_dblMu_tib) + 
                geom_density(aes(x=T00DZ_beta_q2, fill="T00DZ_beta_q2"), alpha=0.7, na.rm=TRUE ) +
                geom_density(aes(x=T50DZ_beta_q2, fill="T50DZ_beta_q2"), alpha=0.7, na.rm=TRUE ) +
                geom_density(aes(x=T99DZ_beta_q2, fill="T99DZ_beta_q2"), alpha=0.7, na.rm=TRUE )
              
              # Cluster Seperation Distribution
              #
              ggplot2::ggplot(data=full_dblMu_tib) + 
                geom_density(aes(x=T00DZ_T50DZ_CSS_mu, fill="T00DZ_T50DZ_CSS_mu"), alpha=0.7, na.rm=TRUE ) +
                geom_density(aes(x=T00DZ_T99DZ_CSS_mu, fill="T00DZ_T99DZ_CSS_mu"), alpha=0.7, na.rm=TRUE ) +
                geom_density(aes(x=T50DZ_T99DZ_CSS_mu, fill="T50DZ_T99DZ_CSS_mu"), alpha=0.7, na.rm=TRUE )
              
              ggplot2::ggplot(data=full_dblMu_tib) + 
                geom_density(aes(x=T00DZ_T50DZ_CSS_q2, fill="T00DZ_T50DZ_CSS_q2"), alpha=0.7, na.rm=TRUE ) +
                geom_density(aes(x=T00DZ_T99DZ_CSS_q2, fill="T00DZ_T99DZ_CSS_q2"), alpha=0.7, na.rm=TRUE ) +
                geom_density(aes(x=T50DZ_T99DZ_CSS_q2, fill="T50DZ_T99DZ_CSS_q2"), alpha=0.7, na.rm=TRUE )
              
              
              # TBD::
              #   Now select Probe_ID,css_names[N] %>% bind_rows()
              #
              #
              
              if (FALSE) {
                full_dblMu_tib <- C_crossSampleLociRSquared(beta_masked_mat, sample_cnt_vec, sample_key_vec, cmb=TRUE, verbose=cpp.verbose) %>% 
                  as.data.frame() %>% tibble::rownames_to_column(var='Probe_ID') %>% 
                  tibble::as_tibble() %>% 
                  dplyr::arrange(-!!dbl_sortMu_key) %>% 
                  dplyr::mutate(Rank=row_number()) %>%
                  dplyr::mutate_if(is.double, list(round), opt$percisionPval)
                readr::write_csv(full_dblMu_tib, full_dbl_csv)
              }
              
              if (FALSE) {
                full_dblQ2_tib <- C_crossSampleLociRSquared(beta_masked_mat, sample_cnt_vec, sample_key_vec, verbose=cpp.verbose) %>% 
                  as.data.frame() %>% tibble::rownames_to_column(var='Probe_ID') %>% 
                  tibble::as_tibble() %>% 
                  dplyr::arrange(-!!dbl_sortQ2_key) %>% 
                  dplyr::mutate(Rank=row_number()) %>%
                  dplyr::mutate_if(is.double, list(round), opt$percisionPval)
                readr::write_csv(full_dblQ2_tib, full_dbl_csv)
              }
              
              system(glue::glue("touch {dbl_beg_txt}"))
              readr::write_csv(full_dblMu_tib,full_dbl_csv)
              system(glue::glue("touch {dbl_end_txt}"))
              Sys.sleep(1)
              
              cat(glue::glue("[{par$prgmTag}]: Done. Writing deltaBetas...{RET}{RET}") )
            }
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
