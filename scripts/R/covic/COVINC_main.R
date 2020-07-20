
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              COVINC:: Workflow for building ML COVIC Pipeline::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Data Preprocessing::
#  - Launch latest swifthoof 1.7
#
# Mrege Processed Data::
#  - merge_builds.R
#
# Model Building::
#  - rank_features.R
#  - build_models.R
#
# Predict Builds::
#  - predict_builds.R
#

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

par$prgmTag <- 'COVINC_main'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# Illumina based directories::
par$macDir <- '/Users/bbarnes/Documents/CustomerFacing'
par$lixDir <- '/illumina/scratch/darkmatter/data'

par$retData     <- FALSE

# Executables::
opt$Rscript <- NULL

# Directories::
opt$buildDir   <- NULL
opt$idatsDir   <- NULL

# The Main driver script should check if each step has 
#   already completed.
#
#  1. Swifhoof

#
# Data Preprocessing::
#  - Launch latest swifthoof 1.7
#
# Mrege Processed Data::
#  - merge_builds.R
#
# Model Building::
#  - rank_features.R
#  - build_models.R
#
# Predict Builds::
#  - predict_builds.R
#

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
  
  if (dir.exists(par$macDir)) par$topDir <- par$macDir
  if (dir.exists(par$lixDir)) par$topDir <- par$lixDir
  
  # Default Parameters for local Mac::
  par$runMode    <- args.dat[1]
  par$srcDir     <- file.path(par$topDir, 'workhorse')
  par$scrDir     <- file.path(par$srcDir, 'scripts')
  par$exePath    <- file.path(par$scrDir, 'R', paste0(par$prgmTag,'.R'))
  
  par$prgmTag <- base::sub('\\.R$', '', base::basename(par$exePath))
  par$locPath <- base::dirname(par$exePath)
  par$scrDir  <- base::dirname(base::normalizePath(par$locPath) )
  par$srcDir  <- base::dirname(base::normalizePath(par$scrDir) )
  par$datDir  <- file.path(par$srcDir, 'dat')
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  
  
}


# End of file
