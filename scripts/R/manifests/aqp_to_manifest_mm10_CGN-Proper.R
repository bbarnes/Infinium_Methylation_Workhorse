
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
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
par$prgmDir <- 'manifests'
par$prgmTag <- 'aqp_to_manifest_mm10'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

par$retData <- FALSE

# Executables::
opt$Rscript <- NULL

# Directories::
opt$outDir     <- NULL

# Required Inputs::
opt$ords <- NULL
opt$mats <- NULL
opt$aqps <- NULL
opt$pqcs <- NULL

opt$idats <- NULL

# Pre-defined files (controls)
opt$ctlCSV <- NULL

# Platform/Method Options::
opt$genomeBuild <- NULL
opt$platform    <- NULL
opt$version     <- NULL

# Run Options::
opt$fresh <- FALSE

opt$percisionSigs <- 1
opt$percisionBeta <- 4
opt$percisionPval <- 6

# Parallel/Cluster Options::
opt$single   <- FALSE
opt$parallel <- TRUE
opt$cluster  <- FALSE

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
  
  # Platform/Method Options::
  opt$genomeBuild <- 'mm10'
  opt$platform    <- 'LEGX'
  opt$version     <- 'B3'
  opt$version     <- 'B4'
  opt$version     <- 'S1'
  opt$version     <- 'S2'
  opt$version     <- 'S5'
  opt$version     <- 'S6'
  
  opt$frmt_original <- TRUE
  opt$frmt_original <- FALSE
  
  opt$write_full <- FALSE
  opt$write_base <- FALSE
  
  opt$make_addresss <- FALSE
  
  opt$addControls <- TRUE
  opt$addManifest <- FALSE
  
  #
  # Pre-defined local options runTypes::
  #
  par$local_runType <- 'nzt'
  par$local_runType <- 'covic'
  par$local_runType <- 'mm10'
  
  if (par$local_runType=='mm10') {
    opt$aqpDir <- file.path(par$topDir, 'data/CustomContent/LifeEpigentics/AQP')
    
    opt$ords <- paste(
      file.path(opt$aqpDir, 'orders/Mus_musculus.order_BP1.csv.gz'),
      file.path(opt$aqpDir, 'orders/Mus_musculus.order_BP2.csv.gz'),
      file.path(opt$aqpDir, 'orders/mm10_LEGX_nonCpG_probes.Jan16-2020.order.csv.gz'),
      file.path(opt$aqpDir, 'orders/LEGX_SpikeIn_Reorder-All-06052020.order.withHeader.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(opt$aqpDir, 'BP1/20420178_AQP1_LifeEpigen_BP1.txt.gz'),
      file.path(opt$aqpDir, 'BP2/20420260_AQP1_LifeEpigen_BP2.txt.gz'),
      file.path(opt$aqpDir, 'BP3/20420260_AQP2_LifeEpigen_BP2.txt.gz'),
      file.path(opt$aqpDir, 'BP4/20455357_AQP1_LifeEpigen_BP4.txt.gz'),
      sep=',')
    
    opt$aqps <- paste(
      file.path(opt$aqpDir, 'AQP_Copy/BS0032527-AQP.txt.gz'),
      file.path(opt$aqpDir, 'AQP_Copy/BS0032533-AQP.txt.gz'),
      file.path(opt$aqpDir, 'AQP_Copy/BS0032545-AQP.txt.gz'),
      file.path(opt$aqpDir, 'AQP_Copy/BS0032636-AQP.txt.gz'),
      sep=',')
    
    opt$pqcs <- paste(
      file.path(opt$aqpDir, 'PQC/20042400_A_ProductQC.txt.gz'),
      sep=',')
    
    par$idatsTopDir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/idats/ILMN_mm10_betaTest_17082020'
    opt$idats <- paste(
      file.path(par$idatsTopDir, '204637490002'),
      sep=',')
  } else if (par$local_runType=='covic') {
    
  } else if (par$local_runType=='nzt') {
    
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$parallel <- TRUE
  
  opt$runName <- paste(opt$genomeBuild,opt$platform,opt$version, sep='-')
  opt$outDir <- file.path(par$topDir, 'scratch')
  
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
    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    
    # Pre-defined files (controls)
    make_option(c("--ords"), type="character", default=opt$ords, 
                help="Order files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mats"), type="character", default=opt$mats, 
                help="Match files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqps"), type="character", default=opt$aqps, 
                help="AQP files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--pqcs"), type="character", default=opt$pqcs, 
                help="PQC files (comma seperated) [default= %default]", metavar="character"),
    
    make_option(c("-i", "--idats"), type="character", default=opt$idats, 
                help="idats directory [default= %default]", metavar="character"),
    
    # Required Inputs::
    make_option(c("--ctlCSV"), type="character", default=opt$ctlCSV, 
                help="Standard Pre-Defined Infinium Methylation Controls CSV (no-header) [default= %default]", metavar="character"),
    
    # Platform/Method Options::
    make_option(c("--genomeBuild"), type="character", default=opt$genomeBuild, 
                help="Genome Build (e.g. hg18, hg36, hg19, hg37, hg38, mm10) [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform (e.g. HM450, EPIC, LEGX, NZT, COVIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest Version (e.g. B0,B1,B2,B3,B4,C0) [default= %default]", metavar="character"),
    
    # Executables::
    make_option(c("--Rscript"), type="character", default=opt$Rscript, 
                help="Rscript path [default= %default]", metavar="character"),
    
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

if (is.null(opt$outDir) || 
    is.null(opt$ords) || is.null(opt$mats) || 
    # (is.null(opt$aqps) && is.null(opt$pqcs)) ||
    # is.null(ctlCSV) || # Can be looked up via dat directory...
    is.null(opt$idats) ||
    is.null(opt$genomeBuild) || is.null(opt$platform) || is.null(opt$version) ||
    is.null(opt$Rscript) ||
    is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )
  
  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  
  if (is.null(opt$outDir))    cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  
  if (is.null(opt$ords))   cat(glue::glue("[Usage]: ords is NULL!!!{RET}"))
  if (is.null(opt$mats))   cat(glue::glue("[Usage]: mats is NULL!!!{RET}"))
  if (is.null(opt$aqps))   cat(glue::glue("[Usage]: aqps is NULL [Not required if pqcs defined]!!!{RET}"))
  if (is.null(opt$pqcs))   cat(glue::glue("[Usage]: pqcs is NULL [Not required if aqps defined]!!!{RET}"))
  if (is.null(opt$idats))  cat(glue::glue("[Usage]: idats is NULL [Not Required]!!!{RET}"))
  
  if (is.null(opt$genomeBuild)) cat(glue::glue("[Usage]: genomeBuild is NULL!!!{RET}"))
  if (is.null(opt$platform))    cat(glue::glue("[Usage]: platform is NULL!!!{RET}"))
  if (is.null(opt$version))     cat(glue::glue("[Usage]: version is NULL!!!{RET}"))
  
  if (is.null(opt$Rscript)) cat(glue::glue("[Usage]: Rscript is NULL!!!{RET}"))
  if (is.null(opt$verbose)) cat(glue::glue("[Usage]: verbose is NULL!!!{RET}"))
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

# Input Definitions::
if (is.null(opt$ctlCSV))
  opt$ctlCsv <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')

ord_vec <- NULL
mat_vec <- NULL
aqp_vec <- NULL
pqc_vec <- NULL
if (!is.null(opt$ords)) ord_vec <- opt$ords %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$mats)) mat_vec <- opt$mats %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$aqps)) aqp_vec <- opt$aqps %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$pqcs)) pqc_vec <- opt$pqcs %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

stopifnot(length(ord_vec)>0)
stopifnot(length(mat_vec)>0)
stopifnot(length(mat_vec)==length(ord_vec))

cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

# Output Definitions::
opt$outDir <- file.path(opt$outDir, par$prgmDir, opt$runName)
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)

opt$intDir <- file.path(opt$outDir, 'intersection')
if (!dir.exists(opt$intDir)) dir.create(opt$intDir, recursive=TRUE)

opt$manDir <- file.path(opt$outDir, 'manifest')
if (!dir.exists(opt$manDir)) dir.create(opt$manDir, recursive=TRUE)

opt$desDir <- file.path(opt$outDir, 'design')
if (!dir.exists(opt$desDir)) dir.create(opt$desDir, recursive=TRUE)

cat(glue::glue("[{par$prgmTag}]: Done. Building Output Directories.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Build Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$matFormat <- 'new'

pqc_man_tib <- NULL
aqp_man_tib <- NULL

fix_pqc_man_tib <- NULL
fix_aqp_man_tib <- NULL

# 296724 x 18 - without filtering
# 295796 x 18 - with filtering
if (length(pqc_vec)!=0) {
  pqc_man_tib <- decodeToManifestWrapper(
    ords=ord_vec, mats=mat_vec, pqcs=pqc_vec, aqps=aqp_vec, 
    platform=opt$platform, version=opt$version,
    matFormat=opt$matFormat,
    pqcFormat='pqc',
    full=par$retData, trim=TRUE,
    verbose=opt$verbose,vt=2,tc=0,tt=pTracker)
  fix_pqc_man_tib <- fixOrderProbeIDs(pqc_man_tib, verbose=opt$verbose,vt=3,tc=0,tt=pTracker)
}

# -329221 x 18 - without filtering
# -296174 x 18 - with filtering
if (length(aqp_vec)!=0) {
  aqp_man_tib <- decodeToManifestWrapper(
    ords=ord_vec, mats=mat_vec, pqcs=pqc_vec, aqps=aqp_vec, 
    platform=opt$platform, version=opt$version, 
    matFormat=opt$matFormat,
    pqcFormat='aqp',
    full=par$retData, trim=TRUE,
    verbose=opt$verbose,vt=2,tc=0,tt=pTracker)
  fix_aqp_man_tib <- fixOrderProbeIDs(aqp_man_tib, verbose=opt$verbose,vt=3,tc=0,tt=pTracker)
}

# QC Sanity Checks for AQP/PQC if present::
#
if (!is.null(pqc_man_tib) && !is.null(aqp_man_tib)) {
  aqp_unq_man_tib <- aqp_man_tib %>% dplyr::anti_join(pqc_man_tib, by=c("M","U"))
  pqc_unq_man_tib <- pqc_man_tib %>% dplyr::anti_join(aqp_man_tib, by=c("M","U"))
  
  aqp_list <- rev(aqp_vec)
  names(aqp_list) <- base::basename(aqp_vec) %>% stringr::str_remove('.gz$') %>% stringr::str_remove('.txt$')
  pqc_name <- base::basename(pqc_vec[1]) %>% stringr::str_remove('.gz$') %>% stringr::str_remove('.txt$')
  
  aqp_unq_tib <- lapply(aqp_list, loadPQC, format='aqp', trim=TRUE) %>% dplyr::bind_rows(.id="AQP_Name") %>%
    dplyr::select(Address,Decode_Status,AQP_Name) %>% dplyr::distinct(Address, .keep_all=TRUE) %>% dplyr::arrange(Address)
  pqc_unq_tib <- loadPQC(pqc_vec[1], format='pqc', trim=TRUE) %>% dplyr::mutate(AQP_Name=pqc_name) %>% 
    dplyr::distinct(Address, .keep_all=TRUE) %>% dplyr::arrange(Address)
  
  # Conclusion:: All missing calls from PQC to AQP are 0 at the most recent AQP stage!!!
  #  - This is what we want to see
  pqc_unq_sum_tib <- dplyr::bind_rows(
    pqc_unq_man_tib %>% dplyr::rename(Decode_Status_PQC=Decode_Status_A) %>%
      dplyr::inner_join(aqp_unq_tib, by=c("U"="Address") ) %>% 
      dplyr::select(Decode_Status_PQC,Decode_Status,AQP) %>% dplyr::group_by_all() %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>%
      dplyr::mutate(Allele="U"),
    pqc_unq_man_tib %>% dplyr::rename(Decode_Status_PQC=Decode_Status_A) %>%
      dplyr::inner_join(aqp_unq_tib, by=c("M"="Address") ) %>% 
      dplyr::select(Decode_Status_PQC,Decode_Status,AQP) %>% dplyr::group_by_all() %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>%
      dplyr::mutate(Allele="M") )
  
  pqc_unq_mis_cnt <- pqc_unq_sum_tib %>% dplyr::filter(Decode_Status_PQC != 0 | Decode_Status != 0) %>% base::nrow()
  if (pqc_unq_mis_cnt!=0) {
    print(pqc_unq_sum_tib)
    stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Failed QC Sanity Check; pqc_unq_mis_cnt={pqc_unq_mis_cnt}"))
  }
  
  # Conclusion:: Need to investigate more; take away is that the AQP results will have failures
  #  NOTE: its surprising how other tangos come up often... Implying Biziprobe has a ranked order of probes to use
  aqp_unq_sum_tib <- aqp_unq_man_tib %>% 
    dplyr::rename(Decode_Status_AQP=Decode_Status_A) %>%
    dplyr::inner_join(pqc_unq_tib, by=c("U"="Address") ) %>% 
    dplyr::select(Decode_Status_AQP,Decode_Status,AQP) %>% dplyr::group_by_all() %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
  aqp_unq_mis_cnt <- aqp_unq_sum_tib %>% dplyr::filter(Decode_Status_AQP != 0 & Decode_Status != 0 & Decode_Status != 1 ) %>% base::nrow()
  if (aqp_unq_mis_cnt!=0) {
    print(pqc_unq_sum_tib)
    stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Failed QC Sanity Check; aqp_unq_mis_cnt={aqp_unq_mis_cnt}"))
  }
  cat(glue::glue("[{par$prgmTag}]: Passed all AQP/PQC Sanity Validation; ",
                 "pqc_unq_mis_cnt={pqc_unq_mis_cnt}, aqp_unq_mis_cnt={aqp_unq_mis_cnt}.{RET}"))
  
  fix_man_sum_tib <- dplyr::full_join(
    dplyr::group_by(fix_pqc_man_tib,Probe_Type) %>% dplyr::summarise(Count=n(), .groups="drop"),
    dplyr::group_by(fix_aqp_man_tib,Probe_Type) %>% dplyr::summarise(Count=n(), .groups="drop"),
    by="Probe_Type", suffix=c("_PQC","_AQP")) %>%
    dplyr::mutate(Count_Dif=Count_AQP-Count_PQC)
  fix_man_mis_cnt <- fix_man_sum_tib %>% dplyr::filter(Count_Dif<=0) %>% base::nrow()
  
  if (fix_man_mis_cnt>0) {
    fix_man_sum_tib %>% print()
    stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Failed QC Sanity Check; PQC>AQP; fix_man_sum_tib={fix_man_sum_tib}"))
  }
  cat(glue::glue("[{par$prgmTag}]: Passed all QC Sanity Check; PQC>AQP Validation; ",
                 "fix_man_sum_tib={fix_man_sum_tib}.{RET}{RET}"))
}

# Use the PQC manifest if not use AQP::
#
raw_man_tib <- NULL
if (length(pqc_vec)!=0) {
  raw_man_tib <- fix_pqc_man_tib
} else if (length(aqp_vec)!=0) {
  raw_man_tib <- fix_aqp_man_tib
} else {
  stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Niether PQC or AQP tibble exists!!!}{RET}{RET}"))
}
raw_man_csv <- file.path(opt$manDir, paste(opt$runName,'manifest.raw.tsv.gz', sep='.') )
raw_man_tib <- raw_man_tib %>% dplyr::arrange(Mat_PrbA)
raw_mat_cnt <- raw_man_tib %>% base::nrow()

cat(glue::glue("[{par$prgmTag}]: Writing Raw Manifest(cnt={raw_mat_cnt})={raw_man_csv}...{RET}") )
readr::write_csv(raw_man_tib,raw_man_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               Format Infinium Methylation Standard Controls::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Come back to this later...
#
if (FALSE) {
  ses_ctl_tib <- NULL
  if (!is.null(opt$ctlCsv) && file.exists(opt$ctlCsv)) {
    org_ctl_tib <- suppressMessages(suppressWarnings(
      readr::read_csv(opt$ctlCsv, col_names=c("Address","Probe_Type","COLOR_CHANNEL","Probe_ID"))
    ))
    
    ses_ctl_tib <- org_ctl_tib %>% 
      dplyr::mutate(Probe_ID=stringr::str_replace_all(Probe_ID, ' ','_')) %>% 
      dplyr::mutate(Probe_ID=stringr::str_replace_all(Probe_ID, '-','_')) %>% 
      dplyr::mutate(Probe_ID=stringr::str_replace_all(Probe_ID, '\\(', '')) %>% 
      dplyr::mutate(Probe_ID=stringr::str_replace_all(Probe_ID, '\\)', '')) %>%
      dplyr::rename(U=Address) %>%
      dplyr::mutate(M=NA,DESIGN='II',col=NA,Probe_Source='IM_Controls',Next_Base=NA) %>%
      dplyr::mutate(M=as.double(M), Probe_ID=paste('ctl',Probe_ID, sep='_')) %>%
      dplyr::select(Probe_ID,M,U,DESIGN,COLOR_CHANNEL,col,Probe_Type,Probe_Source,Next_Base) %>%
      dplyr::arrange(Probe_Type,Probe_ID)
  }
  ses_unq_ctl_tib <- dplyr::distinct(ses_ctl_tib, Probe_ID, .keep_all=TRUE)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Get improbe intersection::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# 1. raw_man -> seq48U
# 2. seq48U -> cgn
# 3. cgn -> top
# 4. top -> prb
# 
# 5. raw_man U prb -> prb_man
# 6. prb_man U imp_des
#
# Remainder (CHN, SNP)
#
# Remainder (BS,NO,neg,ct)
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        1. Raw Manifest -> seq48U::
#
#  raw_man_tib[Mat_PrbA,M,U] -> raw_s48_tsv
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Write manifest unmethylated 48-mer sequences to file
raw_s48_tsv <- file.path(opt$intDir, paste(opt$runName,'raw-s48.tsv', sep='.') )
raw_s48_tib <- raw_man_tib %>% dplyr::distinct(Mat_PrbA,M,U) %>% dplyr::arrange(Mat_PrbA)
readr::write_tsv(raw_s48_tib, raw_s48_tsv, col_names=FALSE)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        2. Raw Manifest seq48U U::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Unclear why these don't gzip correctly. 
#  - May need to be sorted
#  - Use locMac for now
#
imp_s48_tsv <- file.path(par$topDir, 'data/improbe/designOutput_21092020/prb48U/prb48U-GRCh36-38-10-21092020.noHeader.unique.prb48U-sorted.tsv.gz')
imp_s48_tsv <- file.path(par$topDir, 'data/improbe/designOutput_21092020/prb48U/prb48U-GRCh36-38-10-21092020.noHeader.unique.prb48U-sorted.locMac.tsv.gz')

# Intersect manifest unmethylated 48-mer sequences with improbe design database seqU48-mers
#
run_join_cmd <- TRUE
run_join_cmd <- FALSE
imp_int_raw_s48_tsv <- file.path(opt$intDir, paste(opt$runName,'imp-s48.intersect.raw-s48.s48-sorted.tsv.gz', sep='.') )
if (!file.exists(imp_int_raw_s48_tsv) || run_join_cmd) {
  join_cmd <- glue::glue("gzip -dc {imp_s48_tsv} | join -t $'\\", "t' -14 -21 - {raw_s48_tsv} | gzip -c -> {imp_int_raw_s48_tsv}")
  cat(glue::glue("[{par$prgmTag}]: Running: cmd='{join_cmd}'...{RET}{RET}") )
  system(join_cmd)
}

# Load intersection results::
#
imp_int_raw_s48_col <- cols(Mat_PrbA = col_character(),
                            Mat_CGN = col_double(),
                            Mat_TB = col_character(),
                            Mat_CO = col_character(),
                            M = col_double(),
                            U = col_double() )

imp_tar_s48_cgn_tib <- 
  readr::read_tsv(imp_int_raw_s48_tsv, col_names=names(imp_int_raw_s48_col$cols), col_types=imp_int_raw_s48_col) %>% 
  dplyr::inner_join(raw_man_tib, by=c("Mat_PrbA","M","U")) %>%
  dplyr::mutate(Seq_CG=paste0( 'cg',stringr::str_pad(Mat_CGN, width=8, side="left", pad="0") ) ) %>%
  dplyr::select(Seq_CG,Mat_PrbA,Mat_CGN,Mat_TB,Mat_CO,M,U,everything() ) %>%
  dplyr::arrange(Seq_CG)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             3. CGN to Top::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Write improbe cgn's to match against top sequence database::
#
imp_tar_s48_cgn_tsv <- file.path(opt$intDir, paste(opt$runName,'imp-s48.intersect.raw-s48.cgn-sorted.tsv', sep='.') )
readr::write_tsv(dplyr::distinct(imp_tar_s48_cgn_tib,Seq_CG),imp_tar_s48_cgn_tsv, col_names=FALSE)

imp_top_tsv <- file.path(par$topDir, 'data/improbe/designOutput_21092020/GRCh36-GRCh38-GRCm10-21092020.cgnTop.sorted.tsv.gz')

run_join_cmd <- TRUE
run_join_cmd <- FALSE
top_int_raw_cgn_tsv <- file.path(opt$intDir, paste(opt$runName,'imp-top.intersect.raw-cgn.cgn-sorted.tsv.gz', sep='.') )
if (!file.exists(top_int_raw_cgn_tsv) || run_join_cmd) {
  join_cmd <- glue::glue("gzip -dc {imp_top_tsv} | join -t $'\\", "t' -11 -21 - {imp_tar_s48_cgn_tsv} | gzip -c -> {top_int_raw_cgn_tsv}")
  cat(glue::glue("[{par$prgmTag}]: Running: cmd='{join_cmd}'...{RET}{RET}") )
  system(join_cmd)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           4. Top to Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$probe_type <- 'cg'
opt$design_key <- 'Seq_ID'
opt$design_seq <- 'Forward_Sequence'
opt$design_seq <- 'Top_Sequence'
opt$design_prb <- 'Probe_Type'
opt$design_srd <- 'TB'

# Load intersection results::
#
top_int_raw_cgn_col <- cols(Seq_ID = col_character(),
                            Top_Sequence = col_character() )

top_int_raw_cgn_tib <- 
  readr::read_tsv(top_int_raw_cgn_tsv, col_names=names(top_int_raw_cgn_col$cols), 
                  col_types=top_int_raw_cgn_col) %>% dplyr::mutate(Probe_Type=opt$probe_type)


# full_prb_des_csv <- 
#   file.path(opt$desDir, paste(opt$runName,opt$design_seq,opt$probe_type,'all-probes.csv.gz', sep='.') )
full_prb_des_rds <-
  file.path(opt$desDir, paste(opt$runName,opt$design_seq,opt$probe_type,'all-probes.rds', sep='.') )

full_prb_des_tib <- NULL
opt$build_probes <- TRUE
opt$build_probes <- FALSE
if (!file.exists(full_prb_des_rds) || opt$build_probes) {
  cat(glue::glue("[{par$prgmTag}]: Building Full Probe Design...{RET}"))
  
  if (FALSE) {
    test_cnt <- 100000
    test10_tib <- tib2prbs(tib=head(top_int_raw_cgn_tib,n=test_cnt), 
                           idsKey=opt$design_key,
                           seqKey=opt$design_seq,
                           prbKey=opt$design_prb,
                           srdStr=opt$design_srd, 
                           parallel=opt$parallel,
                           verbose=opt$verbose,tc=1,tt=pTracker)
    
    test10_tib %>% dplyr::group_by(SR_Str,CO_Str) %>% dplyr::summarise(SRD_Cnt=n())
    test10_csv <- file.path(opt$outDir, paste(opt$runName,opt$design_seq,opt$probe_type,'100000-probes.csv.gz', sep='.') )
    readr::write_csv(test10_tib,test10_csv)
    test10_rds <- file.path(opt$outDir, paste(opt$runName,opt$design_seq,opt$probe_type,'100000-probes.rds', sep='.') )
    readr::write_rds(test10_tib,test10_rds)
    
    test11_tib <- readr::read_csv(test10_csv)
    test11_tib %>% dplyr::group_by(SR_Str,CO_Str) %>% dplyr::summarise(SRD_Cnt=n())
    
    test1r_tib <- readr::read_rds(test10_rds)
    test1r_tib %>% dplyr::group_by(SR_Str,CO_Str) %>% dplyr::summarise(SRD_Cnt=n())
  }
  
  full_prb_des_tib <- tib2prbs(tib=top_int_raw_cgn_tib, 
                               idsKey=opt$design_key,
                               seqKey=opt$design_seq,
                               prbKey=opt$design_prb,
                               srdStr=opt$design_srd, 
                               parallel=opt$parallel,
                               verbose=opt$verbose,tc=1,tt=pTracker)
  
  full_prb_des_tib %>% dplyr::group_by(SR_Str,CO_Str) %>% dplyr::summarise(SRD_Cnt=n())
  
  cat(glue::glue("[{par$prgmTag}]: Writing Full Probe Design File: RDS={full_prb_des_rds}...{RET}"))
  readr::write_rds(full_prb_des_tib,full_prb_des_rds)
  cat(glue::glue("[{par$prgmTag}]: Done. Writing Full Probe Design File.{RET}{RET}"))
  
} else {
  cat(glue::glue("[{par$prgmTag}]: Loading Full Probe Design={full_prb_des_csv}...{RET}"))
  full_prb_des_tib <- readr::read_csv(full_prb_des_csv)
  cat(glue::glue("[{par$prgmTag}]: Done. Loading Full Probe Design.{RET}{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   5. Match Exact Probe Sequences::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

full_prb_des_MAT_tib <- full_prb_des_tib %>% 
  dplyr::mutate(PRB1_U_MAT=stringr::str_to_upper(PRB1_U),
                PRB1_M_MAT=stringr::str_to_upper(PRB1_M),
                PRB2_D_MAT=stringr::str_to_upper(PRB2_D),
                TB_Str=SR_Str)

# Match by Infinium Type::
#
mat_inf1_tib <- dplyr::inner_join(raw_man_tib, full_prb_des_MAT_tib, 
                                  by=c("AlleleA_Probe_Sequence"="PRB1_U_MAT",
                                       "AlleleB_Probe_Sequence"="PRB1_M_MAT"),
                                  suffix=c("_RAW","_DES"))
mat_inf1_tib %>% dplyr::group_by(Probe_Type_RAW,Infinium_Design) %>% 
  dplyr::summarise(Count=n(), .groups='drop') %>% print()

mat_inf2_tib <- dplyr::inner_join(raw_man_tib,full_prb_des_MAT_tib, 
                                  by=c("AlleleA_Probe_Sequence"="PRB2_D_MAT"), 
                                  suffix=c("_RAW","_DES"))
mat_inf2_tib %>% dplyr::group_by(Probe_Type_RAW,Infinium_Design) %>% 
  dplyr::summarise(Count=n(), .groups='drop') %>% print()

# Get Complete Match
#
mat_full_tib <- dplyr::bind_rows(mat_inf1_tib,mat_inf2_tib) %>%
  dplyr::rename(Probe_Type=Probe_Type_RAW, Seq_ID=Seq_ID_DES) %>% 
  add_count(U,M, name='Tango_Count') %>%
  add_count(U,M,Seq_ID, name='Tango_CGN_Count') %>%
  dplyr::select(Seq_ID,Probe_Type,SR_Str,CO_Str,Infinium_Design,
                Tango_CGN_Count,Tango_Count, dplyr::everything())

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   6. Determine Canonical Sequence::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   7. Intersect Full Improbe Designs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #




mat_full_tib %>% dplyr::select(Seq_ID,Seq_ID_RAW,M,U,Probe_Type,Tango_CGN_Count,Tango_Count) %>% dplyr::filter(Seq_ID==Seq_ID_RAW) %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n())
mat_full_tib %>% dplyr::select(Seq_ID,Seq_ID_RAW,M,U,Probe_Type,Tango_CGN_Count,Tango_Count) %>% dplyr::filter(Seq_ID!=Seq_ID_RAW) %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n())


mat_full_tib %>% # dplyr::select(Seq_ID,Seq_ID_RAW,M,U,Probe_Type,Tango_CGN_Count,Tango_Count) %>% 
  dplyr::filter(Seq_ID!=Seq_ID_RAW) %>% dplyr::filter(Probe_Type=='cg') %>% as.data.frame()










mat_inf1_cnt <- mat_inf1_tib %>% base::nrow()
mat_inf2_cnt <- mat_inf2_tib %>% base::nrow()
raw_man_cnt  <- raw_man_tib  %>% base::nrow()

raw_man_tib %>%
  dplyr::group_by(Probe_Type,Infinium_Design) %>% 
  dplyr::summarise(Count=n(), .groups='drop') %>% print()














# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Previous Workflow In Pieces::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  # TBD:: Below are local input files, they should be made options
  # TBD:: Not sure the difference between the first imp_unq_s48.tsv and the locMac one...
  full_des_tsv    <- file.path(par$topDir, 'data/improbe/designOutput_21092020/GRCm10-21092020_improbe-designOutput.tsv.gz')
  imp_unq_s48_tsv <- file.path(par$topDir, 'data/improbe/designOutput_21092020/prb48U/prb48U-GRCh36-38-10-21092020.noHeader.unique.prb48U-sorted.tsv.gz')
  imp_unq_s48_tsv <- file.path(par$topDir, 'data/improbe/designOutput_21092020/prb48U/prb48U-GRCh36-38-10-21092020.noHeader.unique.prb48U-sorted.locMac.tsv.gz')
  
  # Run Time intermediate files
  raw_mat_s48_tsv <- file.path(opt$outDir, paste(opt$runName,'manifest.sesame-base.s48-sorted.tsv', sep='.') )
  raw_int_s48_tsv <- file.path(opt$outDir, paste(opt$runName,'manifest.sesame-base.s48-sorted.full-join.tsv.gz', sep='.') )
  
  # Write manifest unmethylated 48-mer sequences to file
  #
  readr::write_tsv(raw_man_tib, raw_mat_s48_tsv, col_names=FALSE)
  
  # Intersect manifest unmethylated 48-mer sequences with improbe design database seqU48-mers
  #
  run_join_cmd <- TRUE
  run_join_cmd <- FALSE
  if (!file.exists(raw_int_s48_tsv) || run_join_cmd) {
    join_cmd <- glue::glue("gzip -dc {imp_unq_s48_tsv} | join -t $'\\", "t' -14 -27 - {raw_mat_s48_tsv} | gzip -c -> {raw_int_s48_tsv}")
    cat(glue::glue("[{par$prgmTag}]: Running: cmd='{join_cmd}'...{RET}{RET}") )
    system(join_cmd)
  }
  
  # Load intersection results::
  #
  raw_s48_int_col <- c("Mat_PrbA","Mat_CGN","Mat_TB","Mat_CO", names(dplyr::select(raw_man_tib,-Mat_PrbA)))
  raw_s48_int_tib <- readr::read_tsv(raw_int_s48_tsv, col_names=raw_s48_int_col, guess_max=100000) %>% 
    dplyr::mutate(Seq_CGN=stringr::str_remove(Seq_ID, '^[a-zA-Z]+') %>% stringr::str_remove('^0+') %>% as.double()) %>%
    add_count(Mat_CGN,Seq_CGN, name="Paired_Count")
  
  
  
  
  # Write raw s48 matching CGN's to intersect with full improbe database::
  #
  cgn_imp_int_tsv <- file.path(opt$outDir, paste(opt$runName,'improbe-cgn.s48-sorted.full-join.tsv.gz', sep='.') )
  cgn_s48_tar_tsv <- file.path(opt$outDir, paste(opt$runName,'improbe-cgn.s48-sorted.tsv', sep='.') )
  cgn_s48_tar_tib <- raw_s48_int_tib %>% dplyr::mutate(Seq_ID_MAT=paste0('cg',stringr::str_pad(Mat_CGN,8,'left',pad='0')) ) %>% 
    dplyr::distinct(Seq_ID_MAT) %>% dplyr::arrange(Seq_ID_MAT)
  readr::write_tsv(cgn_s48_tar_tib, cgn_s48_tar_tsv)
  
  # Intersect manifest unmethylated 48-mer matching CGN's with improbe design database::
  #
  run_sub_cmd <- TRUE
  run_sub_cmd <- FALSE
  if (!file.exists(cgn_imp_int_tsv) || run_sub_cmd) {
    sub_exe <- '/Users/bretbarnes/Documents/scripts/subset/getSubset.simple.pl'
    sub_cmd <- glue::glue("{sub_exe} -header -t {cgn_s48_tar_tsv} -d {full_des_tsv} | gzip -c -> {cgn_imp_int_tsv}")
    cat(glue::glue("[{par$prgmTag}]: Running: cmd='{sub_cmd}'...{RET}{RET}") )
    system(sub_cmd)
  }
  
  # Full Common::
  full_mat_tib <- loadIMP(file=cgn_imp_int_tsv,verbose=opt$verbose,vt=4,tc=0,tt=pTracker) %>%
    dplyr::mutate(Mat_CGN=stringr::str_remove(Seq_ID,'^cg') %>% stringr::str_remove('^0+') %>% as.double(),
                  Mat_TB=stringr::str_sub(Methyl_Allele_TB_Strand, 1,1),
                  Mat_CO=Methyl_Allele_CO_Strand) %>% 
    dplyr::select(Mat_CGN,Mat_TB,Mat_CO, everything())
  
  # Join Data:: Unclear what this is for...
  #
  # full_int_tib <- full_mat_tib %>% dplyr::inner_join(raw_s48_int_tib, by=c("Mat_CGN","Mat_TB","Mat_CO"), suffix=c("_IMP","_DES"))
  # full_int_tib %>% dplyr::filter(!is.na(Seq_ID_IMP)) %>% dplyr::filter(!is.na(Seq_ID_DES))
  # full_int_tib %>% dplyr::filter(!is.na(Seq_ID_IMP)) %>% dplyr::filter(!is.na(Seq_ID_DES)) %>% dplyr::filter(Seq_ID_IMP == Seq_ID_DES) %>% dplyr::select(Seq_ID_IMP,Seq_ID_DES)
  # full_int_tib %>% dplyr::filter(!is.na(Seq_ID_IMP)) %>% dplyr::filter(!is.na(Seq_ID_DES)) %>% dplyr::filter(Seq_ID_IMP != Seq_ID_DES) %>% dplyr::select(Seq_ID_IMP,Seq_ID_DES)
  
  # Join Data::
  full_prb_inp_tib <- raw_s48_int_tib %>% dplyr::inner_join(full_mat_tib, by=c("Mat_CGN","Mat_TB","Mat_CO"), suffix=c("_DES","_IMP")) %>% 
    dplyr::select(Seq_ID_IMP,Top_Sequence) %>% dplyr::distinct() %>% 
    dplyr::rename(Seq_ID=Seq_ID_IMP) %>% dplyr::arrange(Seq_ID) %>% dplyr::mutate(Probe_Type='cg')
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Rebuild All Probes Designs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  opt$load_probes  <- FALSE
  opt$load_probes  <- TRUE
  opt$probe_type <- 'cg'
  opt$design_key <- 'Seq_ID'
  opt$design_seq <- 'Forward_Sequence'
  opt$design_seq <- 'Top_Sequence'
  opt$design_prb <- 'Probe_Type'
  opt$design_srd <- 'TB'
  
  full_prb_des_csv <- file.path(opt$outDir, paste(opt$runName,opt$design_seq,'all-strands-current-probes.csv.gz', sep='_') )
  
  if (opt$load_probes && file.exists(full_prb_des_csv)) {
    full_prb_des_tib <- readr::read_csv(full_prb_des_csv)
  } else {
    cat(glue::glue("[{par$prgmTag}]: Building Full Probe Design...{RET}"))
    
    opt$parallel <- TRUE
    # test_cnt <- 10
    # test10_tib <- tib2prbs(tib=head(full_prb_inp_tib,n=test_cnt), 
    full_prb_des_tib <- tib2prbs(tib=full_prb_inp_tib, 
                                 idsKey=opt$design_key,
                                 seqKey=opt$design_seq,
                                 prbKey=opt$design_prb,
                                 srdStr=opt$design_srd, 
                                 parallel=opt$parallel,
                                 verbose=opt$verbose,tc=1,tt=pTracker)
    
    full_prb_des_tib %>% dplyr::group_by(SR_Str,CO_Str) %>% dplyr::summarise(SRD_Cnt=n())
    
    cat(glue::glue("[{par$prgmTag}]: Writing Full Probe Design File: {full_prb_des_csv}...{RET}"))
    readr::write_csv(full_prb_des_tib,full_prb_des_csv)
  }
  raw_man_tib %>% dplyr::group_by(Probe_Type,Infinium_Design) %>% dplyr::summarise(Count=n()) %>% print()
  raw_man_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                Build Probe Matching Sequences for Comparison::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  full_prb_des_MAT_tib <- full_prb_des_tib %>% 
    dplyr::mutate(PRB1_U_MAT=stringr::str_to_upper(PRB1_U),
                  PRB1_M_MAT=stringr::str_to_upper(PRB1_M),
                  PRB2_D_MAT=stringr::str_to_upper(PRB2_D),
                  TB_Str=SR_Str)
  
  mat_inf1_tib <- dplyr::inner_join(raw_man_tib, full_prb_des_MAT_tib, 
                                    by=c("AlleleA_Probe_Sequence"="PRB1_U_MAT",
                                         "AlleleB_Probe_Sequence"="PRB1_M_MAT"),
                                    suffix=c("_RAW","_DES"))
  mat_inf1_tib %>% dplyr::group_by(Probe_Type_RAW,Infinium_Design) %>% 
    dplyr::summarise(Count=n(), .groups='drop') %>% print()
  
  mat_inf2_tib <- dplyr::inner_join(raw_man_tib,full_prb_des_MAT_tib, 
                                    by=c("AlleleA_Probe_Sequence"="PRB2_D_MAT"), 
                                    suffix=c("_RAW","_DES"))
  mat_inf2_tib %>% dplyr::group_by(Probe_Type_RAW,Infinium_Design) %>% 
    dplyr::summarise(Count=n(), .groups='drop') %>% print()
  
  raw_man_tib %>%
    dplyr::group_by(Probe_Type,Infinium_Design) %>% 
    dplyr::summarise(Count=n(), .groups='drop') %>% print()
  
  mat_inf1_cnt <- mat_inf1_tib %>% base::nrow()
  mat_inf2_cnt <- mat_inf2_tib %>% base::nrow()
  raw_man_cnt  <- raw_man_tib  %>% base::nrow()
  
  mat_full_tib <- dplyr::bind_rows(mat_inf1_tib,mat_inf2_tib) %>%
    dplyr::rename(Probe_Type=Probe_Type_RAW,
                  Seq_ID=Seq_ID_DES,
    ) %>% 
    add_count(U,M, name='Tango_Count') %>%
    add_count(U,M,Seq_ID)
  
  mat_full_tib %>% dplyr::group_by(Probe_Type,Infinium_Design,AQP) %>% 
    dplyr::summarise(Count=n(), .groups='drop') %>% print()
  
  mat_full_tib %>% dplyr::group_by(Probe_Type,Infinium_Design,AQP,Tango_Count) %>% 
    dplyr::summarise(Count=n(), .groups='drop') %>% print()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                 Determine Canonical Alignments for Matches::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
}

if (FALSE) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Determine Matching and Mismatch Probes::
  #
  # 1. Matching by Tangos
  # 2. Matching by perfect sequence match (not seq48U) *** MATCH BY EVERYTHING ***
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # This is already the full match::
  mat_full_tib %>% dplyr::group_by(Probe_Type,Infinium_Design,AQP) %>% 
    dplyr::summarise(PT_Miss_Count=n(), .groups='drop')
  
  # Full Matching
  mat_man_tib <- dplyr::inner_join(raw_man_tib, mat_full_tib, 
                                   by=c("M","U", "Infinium_Design",) )
  
  
  # Basic Matching
  mat_man_tib <- raw_man_tib %>% dplyr::inner_join(mat_full_tib, by=c("M","U") )
  mat_sum_tib <- mat_man_tib %>% dplyr::group_by(Probe_Type,Infinium_Design,AQP) %>% 
    dplyr::summarise(PT_Miss_Count=n(), .groups='drop')
  
  # Basic Mis Matching
  mis_man_tib <- raw_man_tib %>% dplyr::anti_join(mat_full_tib, by=c("M","U") )
  mis_sum_tib <- mis_man_tib %>% dplyr::group_by(Probe_Type,Infinium_Design,AQP) %>% 
    dplyr::summarise(PT_Miss_Count=n(), .groups='drop')
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Format Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


#
# TBD:: Add back scratch code::
#   - Control format
#   - Output Sesame/GS
#





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
