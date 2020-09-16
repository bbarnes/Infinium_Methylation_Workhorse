
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

# Program Parameters::
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'manifests'
par$prgmTag <- 'aqp_to_manifest_mm10'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# Illumina based directories::
par$macDir <- '/Users/bbarnes/Documents/Projects/methylation/tools'
par$lixDir <- '/illumina/scratch/darkmatter'

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
opt$parallel <- FALSE
opt$cluster  <- FALSE

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
  
  # Platform/Method Options::
  opt$genomeBuild <- 'mm10'
  opt$platform    <- 'LEGX'
  opt$version     <- 'B3'
  opt$version     <- 'B4'
  opt$version     <- 'S1'
  
  opt$frmt_original <- TRUE
  opt$frmt_original <- FALSE
  
  opt$write_full <- FALSE
  opt$write_base <- FALSE
  
  opt$make_addresss <- FALSE
  
  opt$addControls <- TRUE
  opt$addManifest <- FALSE
  
  opt$aqpDir <- file.path('/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics', 'AQP')
  
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
  
  # ord_bp1_csv <- file.path(opt$aqpDir, 'orders/Mus_musculus.order_BP1.csv.gz')
  # ord_bp2_csv <- file.path(opt$aqpDir, 'orders/Mus_musculus.order_BP2.csv.gz')
  # ord_bp3_csv <- file.path(opt$aqpDir, 'orders/mm10_LEGX_nonCpG_probes.Jan16-2020.order.csv.gz')
  # ord_bp4_csv <- file.path(opt$aqpDir, 'orders/LEGX_SpikeIn_Reorder-All-06052020.order.withHeader.csv.gz')
  #
  # mat1_tsv <- file.path(opt$aqpDir, 'BP1/20420178_AQP1_LifeEpigen_BP1.txt.gz')
  # mat2_tsv <- file.path(opt$aqpDir, 'BP2/20420260_AQP1_LifeEpigen_BP2.txt.gz')
  # mat3_tsv <- file.path(opt$aqpDir, 'BP3/20420260_AQP2_LifeEpigen_BP2.txt.gz')
  # mat4_tsv <- file.path(opt$aqpDir, 'BP4/20455357_AQP1_LifeEpigen_BP4.txt.gz')
  #
  # aqp1_tsv <- file.path(opt$aqpDir, 'AQP_Copy/BS0032527-AQP.txt.gz')
  # aqp2_tsv <- file.path(opt$aqpDir, 'AQP_Copy/BS0032533-AQP.txt.gz')
  # aqp3_tsv <- file.path(opt$aqpDir, 'AQP_Copy/BS0032545-AQP.txt.gz')
  # aqp4_tsv <- file.path(opt$aqpDir, 'AQP_Copy/BS0032636-AQP.txt.gz')
  #
  # pqc_tsv <- file.path(opt$aqpDir, 'PQC/20042400_A_ProductQC.txt.gz')
  
  opt$runName <- paste(opt$genomeBuild,opt$platform,opt$version, sep='-')
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

opt$alnDir <- NULL

# Output Definitions::
# opt$outDir <- file.path(opt$outDir, 'manifest/base')
opt$outDir <- file.path(opt$outDir, par$prgmDir, opt$runName)
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)

opt$outName <- paste(opt$genomeBuild,opt$platform,opt$version, sep='_')
ses_man_csv <- file.path(opt$outDir, paste0(opt$outName,'.manifest.sesame-base.cpg-sorted.csv.gz') )
gss_man_csv <- file.path(opt$outDir, paste0(opt$outName,'.manifest.GenomeStudio.cpg-sorted.csv') )
out_add_csv <- file.path(opt$outDir, paste0(opt$outName,'.manifest.address-base.cpg-sorted.csv.gz') )

cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Build Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$frmt_original <- FALSE
isFull <- FALSE

ord_cnt <- length(ord_vec)

full_man_tib <- NULL
if (!is.null(opt$pqcs)) {
  cat(glue::glue("[{par$prgmTag}]: Running with PQC...{RET}") )
  
  for (idx in c(1:ord_cnt)) {
    full_man_tib <- full_man_tib %>% dplyr::bind_rows(
      decodeToManifest(ord=ord_vec[idx], mat=mat_vec[idx], pqc=pqc_vec[1],
                       platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
                       original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=idx, AQP=idx)
    )
  }
  full_man_tib <- full_man_tib %>% dplyr::arrange(-AQP) %>% dplyr::group_by(U) %>% dplyr::slice_head(n=1)
  
} else {
  cat(glue::glue("[{par$prgmTag}]: Running Non-PQC...{RET}") )
  
  full_man_tib <- NULL
  for (idx in c(1:ord_cnt)) {
    full_man_tib <- full_man_tib %>% dplyr::bind_rows(
      decodeToManifest(ord=ord_vec[idx], mat=mat_vec[idx], aqp=aqp_vec[idx],
                       platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
                       original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=idx, AQP=idx)
    )
  }
}

#
# INVESTIGATION PART:: this can be deleted
#
if (FALSE) {
  full_man_tibA <- NULL
  for (idx in c(1:ord_cnt)) {
    full_man_tibA <- full_man_tibA %>% dplyr::bind_rows(
      decodeToManifest(ord=ord_vec[idx], mat=mat_vec[idx], aqp=aqp_vec[idx],
                       platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
                       original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=idx, AQP=idx)
    )
  }
  full_man_tib1 <- full_man_tib %>% dplyr::arrange(-AQP) %>% dplyr::group_by(U) %>% dplyr::slice_head(n=1)
  
  # Check known case::
  full_man_tib  %>% dplyr::filter(U==99706982) %>% as.data.frame()
  full_man_tibA %>% dplyr::filter(U==99706982) %>% as.data.frame()
  full_man_tib1 %>% dplyr::filter(U==99706982) %>% as.data.frame()
  
  # Check summary::
  full_man_tib  %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()
  full_man_tibA %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()
  full_man_tib1 %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()
  
  # Check counts::
  full_man_tib %>% base::nrow()
  full_man_tib %>% dplyr::distinct(U) %>% base::nrow()
  full_man_tib %>% dplyr::distinct(U,AlleleA_Probe_Sequence) %>% base::nrow()
  
  full_man_tibA %>% base::nrow()
  full_man_tibA %>% dplyr::distinct(U) %>% base::nrow()
  full_man_tibA %>% dplyr::distinct(U,AlleleA_Probe_Sequence) %>% base::nrow()
  
  full_man_tib1 %>% base::nrow()
  full_man_tib1 %>% dplyr::distinct(U) %>% base::nrow()
  full_man_tib1 %>% dplyr::distinct(U,AlleleA_Probe_Sequence) %>% base::nrow()
}

#
# Conclusion Use full_man_tib1 to ensure most recent PQC results!!!
#
# Should be zero below::
#  full_man_tib %>% dplyr::add_count(U, name="U_Tango_Count") %>% dplyr::filter(U_Tango_Count!=1) %>% dplyr::arrange(U) 
full_man_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()

mm10_man_tib <- fixOrderProbeIDs(full_man_tib, verbose=opt$verbose)
mm10_man_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               Format Infinium Methylation Standard Controls::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     mm10 CGN Genome Count:: Global Data
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

mm10_cgn_cnt_col <- c('Genomic_CGN_Count', 'Seq_ID')
mm10_cgn_cnt_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-counts.tsv.gz'
mm10_cgn_cnt_tib <- readr::read_tsv(mm10_cgn_cnt_tsv, col_names=mm10_cgn_cnt_col)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Get improbe intersection::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Probe Info with Seq48U:
#  Unix: /Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top-seqU-to-sequ48.sorted-sequ48.tsv.gz 
imp_s48_tsv  <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top-seqU-to-sequ48.sorted-sequ48.tsv.gz'
int_s48_tsv  <- file.path(opt$outDir, 'mm10_LEGX_cp.manifest.sesame-base.s48-sorted.full-join.tsv.gz')
mm10_s48_tsv <- file.path(opt$outDir, 'mm10_LEGX_cp.manifest.sesame-base.s48-sorted.tsv')
mm10_s48_tib <- mm10_man_tib %>% dplyr::arrange(Mat_PrbA)
readr::write_tsv(mm10_s48_tib, mm10_s48_tsv, col_names=FALSE)

# Probe Info with Seq48U:
#  gzip -dc /Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top-seqU-to-sequ48.sorted-sequ48.tsv.gz
#
# FULL JOIN COMMAND:: Single Line
#  gzip -dc /Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top-seqU-to-sequ48.sorted-sequ48.tsv.gz | join -t $'\t' -14 -27 - /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-cp/mm10_LEGX_cp.manifest.sesame-base.s48-sorted.tsv | gzip -c -> /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-cp/mm10_LEGX_cp.manifest.sesame-base.s48-sorted.full-join.tsv.gz
#
# FULL JOIN COMMAND:: Split Line
# gzip -dc /Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top-seqU-to-sequ48.sorted-sequ48.tsv.gz | \
#  join -t $'\t' -14 -27 - /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-cp/mm10_LEGX_cp.manifest.sesame-base.s48-sorted.tsv | \
#  gzip -c -> /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-cp/mm10_LEGX_cp.manifest.sesame-base.s48-sorted.full-join.tsv.gz

run_join_cmd <- TRUE
run_join_cmd <- FALSE
if (run_join_cmd) {
  join_cmd <- glue::glue("gzip -dc {imp_s48_tsv} | join -t $'\\", "t' -14 -27 - {mm10_s48_tsv} | gzip -c -> {int_s48_tsv}")
  cat(glue::glue("[{par$prgmTag}]: Running: cmd='{join_cmd}'...{RET}{RET}") )
  system(join_cmd)
}

# For completely missing ("off") source files::
#
#  unix: /illumina/scratch/darkmatter/Projects/LifeEpigenetics/data/designInputs/*
#
#  rs = /Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/Redesign/data/SNP/selected_SNP_probes.bed
#

mm10_s48_int_col <- c("Mat_PrbA",'Mat_CGN', 'Mat_TB', 'Mat_CO',
                      "Seq_ID","ID_FR","ID_TB","ID_CO","ID_PD","Infinium_Design","Mat_Prb",
                      "Probe_ID","M","U","DESIGN","COLOR_CHANNEL","col",
                      "Probe_Type","Probe_Source","Next_Base",
                      "AlleleA_Probe_Sequence","AlleleB_Probe_Sequence","Normalization_Bin",
                      "Address_A_Seq","QC_A_Action","Address_B_Seq","QC_B_Action","BP","AQP","Rep_Max","Rep_Num",
                      "AlleleA_Probe_Length","AlleleB_Probe_Length","Old_Probe_ID","Di","DS","HS","FN")
mm10_s48_int_tib <- readr::read_tsv(int_s48_tsv, col_names=mm10_s48_int_col) %>% 
  add_count(Mat_CGN,Seq_ID, name="Paired_Count") %>% 
  dplyr::left_join(mm10_cgn_cnt_tib, by=c("Mat_CGN"="Seq_ID"))


# Split by matching and mismatch Seq_ID
#
mm10_s48_mat_tib <- mm10_s48_int_tib %>% dplyr::filter(Mat_CGN==Seq_ID) %>% dplyr::arrange(Seq_ID)
mm10_s48_mis_tib <- mm10_s48_int_tib %>% dplyr::filter(Mat_CGN!=Seq_ID) %>% dplyr::arrange(Seq_ID)
#  AND make a list of missing probes from original data::
mm10_s48_off_tib <- mm10_man_tib %>% dplyr::anti_join(mm10_s48_int_tib, by="Mat_PrbA")

# QC Counts Matching::
mm10_s48_int_tib %>% base::nrow()
mm10_s48_mat_tib %>% base::nrow()
mm10_s48_mis_tib %>% base::nrow()
mm10_s48_off_tib %>% base::nrow()

# QC: Check type distributions::
#  mat = cg
#  mis = rp,mu
mm10_s48_mat_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Grp_Count=n())
mm10_s48_mis_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Grp_Count=n())
mm10_s48_off_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Grp_Count=n())

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Manifest Improbe Matching:: CGN ONLY
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Write improbe lookup file::
mm10_s48_int_mat_cgn_tsv <- file.path(opt$outDir, 'mm10_s48_int_mat_cgn_tsv.gz')
mm10_s48_int_mat_cgn_tib <- mm10_s48_mat_tib %>% dplyr::distinct(Mat_CGN) %>% dplyr::arrange(Mat_CGN) %>%
  dplyr::rename(Seq_ID=Mat_CGN)
readr::write_tsv(mm10_s48_int_mat_cgn_tib, mm10_s48_int_mat_cgn_tsv)

# /Users/bbarnes/Documents/Projects/scripts/subset/getSubset.simple.pl -header -CO C -t /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-cp/mm10_s48_int_mat_cgn_tsv.gz -d /Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.tsv.gz | head
#  > /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-cp/mm10_s48_int_mat.improbe.tsv

mm10_s48_int_mat_imp_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.LEGX-mat.tsv.gz'
mm10_s48_int_mat_imp_tib <- read_tsv(mm10_s48_int_mat_imp_tsv) %>% 
  dplyr::mutate(Imp_FR=Methyl_Allele_FR_Strand,
                Imp_TB=stringr::str_sub(Methyl_Allele_TB_Strand, 1,1),
                Imp_CO=Methyl_Allele_CO_Strand)

#
# Join improbe data to designs::
#
mm10_man_mat_imp_tib <- mm10_s48_mat_tib %>% 
  dplyr::left_join(mm10_s48_int_mat_imp_tib, by=c("Seq_ID", "Probe_Type","Mat_TB"="Imp_TB", "Mat_CO"="Imp_CO")) %>% 
  dplyr::mutate(Inter_Type='T') %>% tidyr::unite(Probe_ID_Suffix, Mat_TB,Mat_CO,Infinium_Design,Inter_Type,Rep_Num, sep='', remove=FALSE) %>% 
  tidyr::unite(Probe_ID, Seq_ID,Probe_ID_Suffix, sep='_')

mm10_man_mat_ses_tib <- mm10_man_mat_imp_tib %>%
  dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, AlleleA_Probe_Sequence,AlleleB_Probe_Sequence)

# QC Counts Matching::
mm10_man_mat_imp_tib %>% base::nrow()
mm10_s48_mat_tib %>% base::nrow()
mm10_man_mat_imp_tib %>% dplyr::distinct(Probe_ID) %>% base::nrow()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Manifest Improbe Matching:: Off ONLY
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Split ch/mu/rp analysis
mm10_s48_off_types <- mm10_s48_off_tib %>% split(.$Probe_Type)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Manifest Improbe Matching:: rs ONLY
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Design Source::
#  rs = /Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/Redesign/data/SNP/selected_SNP_probes.bed

#
# TBD:: Add back previous SNPs!!!!!       dplyr::filter(TB ===== 'N' & !is.na(TB) )         Missing most SNPs now...
#
mm10_man_newRS_ses_tib <- mm10_s48_off_types[['rs']] %>% dplyr::filter(TB != 'N' & !is.na(TB) ) %>%
  dplyr::mutate(Inter_Type='T') %>% 
  tidyr::unite(Probe_ID_Suffix, TB,CO,Infinium_Design,Inter_Type,Rep_Num, sep='', remove=FALSE) %>% 
  tidyr::unite(Probe_ID, Seq_ID,Probe_ID_Suffix, sep='_') %>%
  dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, AlleleA_Probe_Sequence,AlleleB_Probe_Sequence)

mm10_man_newRS_ses_tib %>% base::nrow()
mm10_man_newRS_ses_tib %>% dplyr::distinct(Probe_ID) %>% base::nrow()


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Sesame Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# current = B2
outName <- paste(opt$platform,opt$version, sep='-')
outFile <- paste(outName,'manifest.sesame-base.cpg-sorted.csv.gz', sep='.')
mm10_ses_man_cgn_out_csv <- file.path(opt$outDir, outFile)
mm10_ses_man_cgn_git_csv <- file.path(par$datDir,'manifest/base',outFile)
mm10_ses_man_cgn_tib <- dplyr::bind_rows(mm10_man_mat_ses_tib,
                                         mm10_man_newRS_ses_tib,
                                         ses_unq_ctl_tib) %>% dplyr::arrange(Probe_ID)
readr::write_csv(mm10_ses_man_cgn_tib,mm10_ses_man_cgn_out_csv)
readr::write_csv(mm10_ses_man_cgn_tib,mm10_ses_man_cgn_git_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Manifest Improbe Matching:: Non ONLY
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Split ch/mu/rp analysis
mm10_s48_mis_types <- mm10_s48_mis_tib %>% split(.$Probe_Type)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Manifest Improbe Matching:: mu ONLY
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Conclusion:: mu probes can be converted to cg and added back to mat_tib with extra cg's added as extra field...
#
clean_mu_all_tib <- mm10_s48_mis_types[['mu']] %>% dplyr::mutate(Seq_ID=stringr::str_replace(Seq_ID, '^mu', 'cg'))
clean_mu_mat_tib <- clean_mu_all_tib %>% dplyr::filter(Mat_CGN==Seq_ID) %>% dplyr::add_count(Seq_ID,Mat_TB,Mat_CO,Infinium_Design, name="ID_Count")
clean_mu_mis_tib <- clean_mu_all_tib %>% dplyr::filter(! Seq_ID %in% clean_mu_mat_tib$Seq_ID)

clean_mu_mat_tib %>% dplyr::filter(ID_Count!=1) %>% as.data.frame()

# Write improbe lookup file::
mm10_s48_int_mis_cgn_tsv <- file.path(opt$outDir, 'mm10_s48_int_mis_cgn_tsv.gz')
mm10_s48_int_mis_cgn_tib <- clean_mu_all_tib %>% dplyr::distinct(Mat_CGN) %>% dplyr::arrange(Mat_CGN) %>%
  dplyr::rename(Seq_ID=Mat_CGN)
# readr::write_tsv(mm10_s48_int_mis_cgn_tib, mm10_s48_int_mis_cgn_tsv)








# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               Quick Look at Auto Sample Sheets from CG Only::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# mm10_hum_ss_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/Laird-IDs_to_SampleNames_basic.csv'
# mm10_hum_ss_tib <- readr::read_csv(mm10_hum_ss_csv) %>% purrr::set_names(c('Laird_ID', 'Sentrix_Pos', 'Cell_Type', 'Sample_Class'))

# Possible mismatch: 20364 = 20382
RepAC_ID <- 20364
RepAC_ID <- 20382

sam_map_tib <- tibble::tibble(
  Laird_ID   =c( RepAC_ID,  20384,  20385,  21026,   20010,   20015,   20012),
  Sample_Name=c('RepAC','RepS3','RepM1','RepSA', 'T00DZ', 'T50DZ', 'T99DZ')
)

tit_map_tib <- tibble::tibble(
  Laird_ID   =c( 20010,   20015,   20012),
  Sample_Name=c('T00DZ', 'T50DZ', 'T99DZ')
)

rep_map_tib <- tibble::tibble(
  Laird_ID   =c( RepAC_ID,  20384,  20385,  21026),
  Sample_Name=c('RepAC','RepS3','RepM1','RepSA')
)

# Load ILS Sample Sheet::
#
ils_map_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/ILS-Mouse_Methylation_samplesheet.csv'
ils_map_tib <- readr::read_csv(ils_map_csv) %>%
  dplyr::mutate(Laird_ID=as.integer(stringr::str_remove(Sample_ID, '_.*$')),
                Sentrix_ID=SentrixBarcode_A,
                Sentrix_Pos=SentrixPosition_A,
                Sentrix_Name=paste(Sentrix_ID,Sentrix_Pos, sep='_')
  ) %>% 
  dplyr::select(Laird_ID,Sentrix_Name,Sentrix_ID,Sentrix_Pos) %>% dplyr::mutate(Lab='ILS')

# Load VAI Sample Sheet::
#
vai_map_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/VAI_SampleSheetPlate1-2.no-header.csv'
vai_map_tib <- readr::read_csv(vai_map_csv) %>% 
  dplyr::rename(Laird_Str=Sample_Name, Sentrix_Pos=Sentrix_Position) %>% 
  dplyr::mutate(Laird_ID=stringr::str_remove(Laird_Str, '_.*$') %>% as.integer(),
                Sentrix_Name=paste(Sentrix_ID,Sentrix_Pos, sep='_')) %>% 
  dplyr::select(Laird_ID,Sentrix_Name,Sentrix_ID,Sentrix_Pos) %>% dplyr::mutate(Lab='VAI')

#
# Build beta full Sample Sheet::
#

# Full Analytical Sample Sheet::
beta_ss_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/analytical/ILS-VAI.analytical_SampleSheet.csv.gz'
beta_ss_tib <- dplyr::bind_rows(
  dplyr::inner_join(ils_map_tib, sam_map_tib, by="Laird_ID"),
  dplyr::inner_join(vai_map_tib, sam_map_tib, by="Laird_ID") ) %>%
  dplyr::select(Sentrix_Name,Sample_Name, everything())
readr::write_csv(beta_ss_tib,beta_ss_csv)

# Titration Analytical Sample Sheet::
beta_tit_ss_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/analytical/Titration.ILS-VAI.analytical_SampleSheet.csv.gz'
beta_tit_ss_tib <- beta_ss_tib %>% dplyr::filter(Laird_ID %in% tit_map_tib$Laird_ID)
readr::write_csv(beta_tit_ss_tib,beta_tit_ss_csv)

# Replicate Analytical Sample Sheet::
beta_rep_ss_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/analytical/Replicates.ILS-VAI.analytical_SampleSheet.csv.gz'
beta_rep_ss_tib <- beta_ss_tib %>% dplyr::filter(Laird_ID %in% rep_map_tib$Laird_ID)
readr::write_csv(beta_rep_ss_tib,beta_rep_ss_csv)

# cp /Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/sampleSheets/betaTest/analytical/* /Users/bbarnes/Documents/Projects/methylation/tools/Infinium_Methylation_Workhorse/dat/sampleSheets/mm10/


#
# Done with Sample Sheets!!!
#

#
#
# Need to write Titration and Replicates Sample Sheets seperately...
#
#

#
# Actual Idats and plots::
#
ils_dat_dir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/scratch/swifthoof_main/ILMN_mm10_betaTest_17082020'
ils_ss_tib  <- loadAutoSampleSheets(dir=ils_dat_dir, verbose=opt$verbose) %>% dplyr::mutate(Lab='ILS')

van_dat_dir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/scratch/swifthoof_main/VanAndel_mm10_betaTest_31082020'
van_ss_tib  <- loadAutoSampleSheets(dir=van_dat_dir, verbose=opt$verbose) %>% dplyr::mutate(Lab='Van')

fox_dat_dir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/scratch/swifthoof_main/MURMETVEP_mm10_betaTest_06082020'
fox_ss_tib  <- loadAutoSampleSheets(dir=fox_dat_dir, verbose=opt$verbose) %>% dplyr::mutate(Lab='FOX')

mm10_ss_tib <- dplyr::bind_rows(ils_ss_tib, van_ss_tib, fox_ss_tib) %>% dplyr::group_by(Lab)
ggplot2::ggplot(data=mm10_ss_tib, aes(x=Beta_2_Mean, y=Poob_Pass_0_Perc, color=Lab)) + ggplot2::geom_point() # + ylim(80,100)
ggplot2::ggplot(data=mm10_ss_tib, aes(x=Beta_2_Mean, y=Poob_Pass_0_Perc, color=Lab)) + ggplot2::geom_point() + ylim(75,100)


#
# Evidence For Missing Replicate Group::
#
ils_ss_tib %>% dplyr::inner_join(beta_ss_tib, by="Sentrix_Name") %>% dplyr::group_by(Sample_Name) %>% dplyr::summarise(SGroup_Count=n())

#
# Current Merged Data::
#
lifeDir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics'

mm10_tit_ss_merge_csv <- file.path(lifeDir, 'scratch/merge_builds/LEGX/S1/Sample_Name/mm10-ILS-VAI.Titration/mm10-ILS-VAI.Titration_LEGX_S1_AutoSampleSheet.csv.gz')
mm10_tit_ss_merge_tib <- readr::read_csv(mm10_tit_ss_merge_csv)

mm10_rep_ss_merge_csv <- file.path(lifeDir, 'scratch/merge_builds/LEGX/S1/Sample_Name/mm10-ILS-VAI.Replicates/mm10-ILS-VAI.Replicates_LEGX_S1_AutoSampleSheet.csv.gz')
mm10_rep_ss_merge_tib <- readr::read_csv(mm10_rep_ss_merge_csv)





















# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Source and Design Search::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  # Objectives::
  #
  #  1. Gather all mouse idats tangos from a few recent idats
  #
  #  2. Gather all public product improbe GGNs and designs 
  #     - Add Horvath Manifests
  #     - Add NZT Manifest
  #     - Add Excalibur Manifest
  #     - Add Custom Manifests (Genknowme, TruDX, Univ Chicago)
  #  3. Gather all mouse methylation improbe CGNs and designs
  #
  #  4. Comapre Previous Products with Mouse for CGN -> Top Integrity
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Identify Address Not Found in IDATs::
  #
  #  1. Gather all mouse idats tangos from a few recent idats
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  par$compare_idats <- FALSE
  if (par$compare_idats) {
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                         Load Tangos From IDATS::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    idat_vec <- NULL
    if (!is.null(opt$idats)) idat_vec <- opt$idats %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
    
    max_idat_load <- 3
    prefix_list <- sesame::searchIDATprefixes(idat_vec[1], recursive=TRUE, use.basename=TRUE)
    idat_add_tib <- NULL
    for (pIdx in c(1:length(prefix_list))) {
      idat_cur_dat <- prefixToIdat(prefix=prefix_list[[pIdx]], verbose=opt$verbose)
      idat_add_tib <- idat_add_tib %>% dplyr::bind_rows(idat_cur_dat$sig %>% dplyr::select(Address))
      
      if (pIdx >= max_idat_load) break
    }
    adds_cnt_tib <- idat_add_tib %>% dplyr::group_by(Address) %>% dplyr::summarise(Count=n())
    idat_add_vec <- adds_cnt_tib %>% dplyr::pull(Address)
    
    # Sanity Check
    adds_mis_cnt <- adds_cnt_tib %>% dplyr::filter(Count != max_idat_load) %>% base::nrow()
    stopifnot(adds_mis_cnt==0)
    
    add_mis_man_tib <- dplyr::bind_rows(
      mm10_man_tib %>% dplyr::filter(! U %in% idat_add_vec),
      mm10_man_tib %>% dplyr::filter(!is.na(M)) %>% dplyr::filter(! M %in% idat_add_vec)
    )
    add_mis_ctl_tib <- dplyr::bind_rows(
      ses_ctl_tib %>% dplyr::filter(! U %in% idat_add_vec),
      ses_ctl_tib %>% dplyr::filter(!is.na(M)) %>% dplyr::filter(! M %in% idat_add_vec)
    )
    
    # List of Probe IDs with missing tangos in idats from Wanding
    add_mis_cgn_wand_tsv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests/tango-issues/probe_with_tango_issue.tsv.gz'
    add_mis_cgn_wand_tib <- readr::read_tsv(add_mis_cgn_wand_tsv)
    
    #
    # CONCLUSION: I don't see any missing addresses: Only checked ILMN so far...
    #  TBD:: Check IDATs from other BETA sites
    #
    
    cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Identify Address Not Found in IDATs!{RET}{RET}") )
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #           Preprocessing:: Improbe Previous Product Designs
  #
  #  2. Gather all public product improbe GGNs and designs 
  #     - Add Horvath Manifests
  #     - Add NZT Manifest
  #     - Add Excalibur Manifest
  #     - Add Custom Manifests (Genknowme, TruDX, Univ Chicago)
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  imp_prd_tsv <- '/Users/bbarnes/Documents/Projects/methylation/EWAS/improbe/EPIC-reorder.improbe-design.tsv.gz'
  imp_prd_tib <- suppressMessages(suppressWarnings( readr::read_tsv(imp_prd_tsv) )) %>% 
    dplyr::mutate(Min_Final_Score=pmin(Methyl_Final_Score,UnMethyl_Final_Score),
                  Methyl_Allele_TB_Strand=stringr::str_sub(Methyl_Allele_TB_Strand, 1,1))
  
  #
  # Investigate Integrety of CGN -> Top_Sequence in Previous Products::
  #
  
  # Build Summary:: i.e. 27k vs. 450k/EPIC
  imp_prd_tib %>% dplyr::group_by(Genome_Build) %>% dplyr::summarise(Count=n()) %>% print()
  
  # Make sure CGN -> Top is unique::
  prd_uniq_cgn_top_tib <- imp_prd_tib %>% dplyr::select(Seq_ID, Top_Sequence) %>% dplyr::distinct()
  
  # Check Top_Sequences all have a Count == 1
  prd_top_cnt_tib <- prd_uniq_cgn_top_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Count=n())
  
  # There should be 7 and they should be from the 27k:: i.e. Genome_Build == 36
  prd_top_mis_tib <- prd_top_cnt_tib %>% dplyr::filter(Count != 1)
  prd_bad_cgn_tib <- imp_prd_tib %>% dplyr::filter(Top_Sequence %in% prd_top_mis_tib$Top_Sequence) %>% 
    dplyr::select(Seq_ID, Genome_Build) %>% dplyr::distinct()
  prb_bad_cgn_cnt <- prd_bad_cgn_tib %>% base::nrow()
  
  # Check CGNs all have a Count == 1
  prd_cgn_cnt_tib <- prd_uniq_cgn_top_tib %>% dplyr::group_by(Seq_ID) %>% dplyr::summarise(Count=n())
  
  # This should be zero::
  prd_cgn_mis_tib <- prd_cgn_cnt_tib %>% dplyr::filter(Count != 1)
  
  cat(glue::glue("[{par$prgmTag}]: Done. Investigating Previous Product CGN -> Top_Sequence Designs; prb_bad_cgn_cnt={prb_bad_cgn_cnt}!{RET}{RET}") )
  
  #
  # Investigate uniqueness of probes::
  #
  imp_unq_tib <- imp_prd_tib %>% dplyr::distinct(Seq_ID, Methyl_Allele_CO_Strand, Methyl_Allele_TB_Strand, UnMethyl_Probe_Sequence)
  imp_mis_tib <- imp_unq_tib %>% dplyr::group_by(UnMethyl_Probe_Sequence) %>% dplyr::summarise(Seq_Count=n())
  imp_ann_tib <- imp_unq_tib %>% dplyr::left_join(imp_mis_tib, by="UnMethyl_Probe_Sequence")
  
  # imp_ann_tib %>% dplyr::filter(Seq_Count!=1, Methyl_Allele_CO_Strand=='C')
  imp_bad_prb_tib <- imp_ann_tib %>% dplyr::filter(Seq_Count!=1, Methyl_Allele_CO_Strand=='C') %>% dplyr::group_by(Seq_ID) %>% 
    dplyr::summarise(CGN_Count=n()) %>% dplyr::filter(CGN_Count!=1)
  prd_bad_prb_tib <- imp_prd_tib %>% dplyr::filter(Seq_ID %in% imp_bad_prb_tib$Seq_ID) %>% dplyr::distinct(Seq_ID, Genome_Build)
  prd_bad_prb_cnt <- prd_bad_prb_tib %>% base::nrow()
  
  cat(glue::glue("[{par$prgmTag}]: Done. Investigating Previous Product Non-Unique Prb_Seq_U; prd_bad_prb_cnt={prd_bad_prb_cnt}!{RET}{RET}") )

  #
  # Join all Failed Probes due to CGN->Top or UnMethyl_Probe_Sequeunce Integrity Issues::
  #
  prd_fin_bad_tib <- dplyr::bind_rows( prd_bad_cgn_tib, prd_bad_prb_tib ) %>% dplyr::distinct()
  prd_fin_bad_cnt <- prd_fin_bad_tib %>% base::nrow()

  cat(glue::glue("[{par$prgmTag}]: Done. Investigating Previous Product Integrity Failures; prd_fin_bad_cnt={prd_fin_bad_cnt}!{RET}{RET}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                Preprocessing:: Improbe Mouse Designs
  #
  #  3. Gather all mouse methylation improbe CGNs and designs
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # mm10_imp_AL_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top-seqU.tsv.gz'
  # mm10_imp_AL_tib <- readr::read_tsv(mm10_imp_AL_tsv)
  
  # Ensure we have the correct matching database between our designs, our design database (v1) and the full database (AL)
  mm10_imp_v1_tsv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/data/dropbox/merged_with_raw_ordered_withHeader.tsv.gz'
  mm10_imp_v1_tib <- suppressMessages(suppressWarnings( readr::read_tsv(mm10_imp_v1_tsv) )) %>% 
    dplyr::mutate(
      Seq_ID = stringr::str_remove_all(Seq_ID, ' '),
      Mat_Prb1=stringr::str_sub(UnMethyl_Probe_Sequence, 2,50), 
      Mat_Prb2=stringr::str_sub(UnMethyl_Probe_Sequence, 2,49) )
  
  # Non Traditional Probes::
  #  mm10_imp_v1_tib %>% dplyr::filter(Customer_Probe_IU!=UnMethyl_Probe_Sequence) 
  
  #
  # Make sure both V1 have CGN->TOP Integrity::
  #
  mm10_uniq_cgn_top_tib <- mm10_imp_v1_tib %>% dplyr::distinct(Seq_ID,Top_Sequence)
  mm10_uniq_top_cnt_tib <- mm10_uniq_cgn_top_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Top_Count=n())
  mm10_uniq_top_mis_cnt <- mm10_uniq_top_cnt_tib %>% dplyr::filter(Top_Count != 1) %>% base::nrow()
  
  cat(glue::glue("[{par$prgmTag}]: Done. Investigating Mouse Product Integrity Failures; ",
                 "mm10_uniq_top_mis_cnt={mm10_uniq_top_mis_cnt}!{RET}{RET}") )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                Preprocessing:: Improbe Mouse Designs
  #
  #  4. Comapre Previous Products with Mouse for CGN -> Top Integrity
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  prd_mm10_uniq_cgn_top_tib <- dplyr::bind_rows( prd_uniq_cgn_top_tib, mm10_uniq_cgn_top_tib ) %>% dplyr::distinct()
  prd_mm10_uniq_top_cnt_tib <- prd_mm10_uniq_cgn_top_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Top_Count=n())
  
  # Validation:: Only loss of integrity comes from 27k::
  prd_mm10_uniq_top_min_tib <- prd_mm10_uniq_top_cnt_tib %>% dplyr::filter(Top_Count != 1) %>% 
    dplyr::inner_join(imp_prd_tib %>% dplyr::select(Seq_ID,Top_Sequence,Genome_Build), by="Top_Sequence") %>% dplyr::distinct(Seq_ID, Genome_Build)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                Preprocessing:: Update mouse CGN names::
  #
  #  5. Use improbe mouse designs to generate new proper CGNs
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  #
  # LAST POINT BEFORE RESTART: TBD:
  #  - Validate designs for all probes on mm10_man_tib
  #
  par$tbd_step <- FALSE
  if (par$tbd_step) {
    mm10_man_v1_matPrbCgn_tib <- dplyr::bind_rows(
      mm10_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb1", "Seq_ID")),
      mm10_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb2", "Seq_ID")) )

    mm10_man_v1_matPrb_tib <- dplyr::bind_rows(
      mm10_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb1")),
      mm10_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb2")) )
    
    # Can't Anit-join until we remove the non Seq_ID matches::
    # mm10_man_tib %>% dplyr::anti_join(mm10_man_v1_matPrb_tib, by='')
    mm10_man_v1_matPrb_tib 
    
    
    # Probes with non-match Seq_ID's
    mm10_man_v1_matPrb_tib %>% dplyr::filter(Seq_ID.x != Seq_ID.y) %>% dplyr::select(Seq_ID.x, Seq_ID.y, Top_Sequence, Mat_Prb)
    
    mm10_man_v1_matPrb_tib %>% dplyr::filter(Seq_ID.x != Seq_ID.y) %>% dplyr::select(Seq_ID.x, Seq_ID.y) %>% as.data.frame()
    
  }
  
  #
  # NEXT STEP:: 
  #
  #  - Check missing probes from basic sort above
  #
  #  - Perl PrbSeqU_48 and sort by it from full design file.
  #  - Intersect by Mat_Prb
  #  - Checking missing probes
  #  - Identify RP/MU probes and add their values
  #
  mm10_des_AL_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top-seqU-to-sequ48.tsv.gz'
  mm10_des_AL_tib <- readr::read_tsv(mm10_des_AL_tsv)
  

  #
  # Full Unique:: 
  #
  par$full_cross_verification <- FALSE
  if (par$full_cross_verification) {
    #
    # Takes too much memory locally::
    #  mm10_imp_AL_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top-seqU.tsv.gz'
    #  mm10_imp_AL_tib <- readr::read_tsv(mm10_imp_AL_tsv)
    #
    # Next step: cross validate all CGN/TOP against all manifest designs...
    #  Unique Count = 20239851
    #
    mm10_imp_AL_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top.uniq.tsv'
    mm10_imp_AL_rds <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.cgn-top.uniq.rds'
    # mm10_imp_AL_tib <- readr::read_tsv(mm10_imp_AL_tsv)
    # readr::write_rds(mm10_imp_AL_tib, mm10_imp_AL_rds,  compress="gz")
    mm10_imp_AL_tib <- readr::read_rds(mm10_imp_AL_rds)
    
    prd_mm10_full_cgn_top_tib <- dplyr::bind_rows( prd_uniq_cgn_top_tib, mm10_uniq_cgn_top_tib, mm10_imp_AL_tib ) %>% dplyr::distinct()
    prd_mm10_full_top_cnt_tib <- prd_mm10_full_cgn_top_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Top_Count=n())
    prd_mm10_full_top_cnt_csv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.prd_mm10_full_top_cnt.csv.gz'
    readr::write_csv(prd_mm10_full_top_cnt_tib, prd_mm10_full_top_cnt_csv)
    
    # Same ones as 27k::
    prd_mm10_full_top_cnt_tib %>% dplyr::filter(Top_Count != 1)

    # Conclusion: Cross Validation Validated!
  }
  
}




















if (FALSE) {
  
  #
  # Cannot trust below becuase coordinate system is off for columns:: Also, the above seems to work...
  #
  if (FALSE) {
    
    #
    # TBD:: Repeat above with AL improbe mouse designs for integrity check...
    #  4. Comapre Previous Products with Mouse for CGN -> Top Integrity
    #
    
    #
    # Complete mouse improbe design::
    #
    mm10_imp_org_csv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/complicated-header.csv'
    mm10_imp_org_col <- suppressMessages(suppressWarnings( readr::read_csv(mm10_imp_org_csv) )) %>% names() %>% as.vector()
    # mm10_imp_org_col %>% length() = 92
    mm10_imp_ext_col <- c('Add_Forward_Sequence','Add_Genome_Build','Add_Chromosome','Add_Position','Add_Top_Sequence',
                          'Add_Probe_Seq_M','Add_FR','Add_TB','Add_CO','Add_Score_M','Add_CpgCnt','Add_CpgPos','Add_Probe_Seq_U',
                          'Add_Score_U')
    mm10_imp_AL_col <- c(mm10_imp_org_col,mm10_imp_ext_col)
    # mm10_imp_AL_col %>% length() = 106
    
    # mm10_imp_AL_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests/joined.designed.dat.head.csv'
    mm10_imp_AL_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests/joined.designed.dat.csv.gz'
    mm10_imp_AL_tib <- suppressMessages(suppressWarnings( readr::read_csv(mm10_imp_AL_csv, col_names=mm10_imp_AL_col) )) %>% 
      dplyr::mutate(
        Mat_Prb1=stringr::str_sub(Add_Probe_Seq_U, 2,50), 
        Mat_Prb2=stringr::str_sub(Add_Probe_Seq_U, 2,49) )
  }
  
  # Join manifest probes and improbe based on Allele_A probe
  if (FALSE) {
    # v1 improbe::
    mm10_man_v1_mat_tib <- dplyr::bind_rows(
      mm10_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb1")),
      mm10_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb2")) )
    
    # AL improbe::
    mm10_man_AL_mat_tib <- dplyr::bind_rows(
      mm10_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb1")),
      mm10_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb2")) )
    
    
    #
    # Check improbe design uniqueness by CGN->TOP Integrity::
    #
    mm10_imp_v1_tib %>% base::nrow() # 285604
    mm10_imp_AL_tib %>% base::nrow() # 1108632
    
    mm10_unq_v1_cgTop_tib <- mm10_imp_v1_tib %>% dplyr::distinct(Seq_ID, Top_Sequence)
    mm10_unq_AL_cgTop_tib <- mm10_imp_AL_tib %>% dplyr::distinct(Probe_ID, Top_Sequence) %>% dplyr::rename(Seq_ID=Probe_ID)
    
    # Failed in translation: Probe_ID=='cg00554227'
    
    #
    # By the numbers::
    #
    mm10_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb1")) %>% base::nrow()
    mm10_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb1")) %>% base::nrow()
    
    mm10_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb2")) %>% base::nrow()
    mm10_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb2")) %>% base::nrow()
    
    
    #
    # Probe IDs that don't match::
    #
    mm10_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb1")) %>% dplyr::filter(Seq_ID.x != Seq_ID.y) %>% as.data.frame()
    mm10_man_tib %>% dplyr::inner_join(mm10_imp_v1_tib, by=c("Mat_Prb"="Mat_Prb2")) %>% dplyr::filter(Seq_ID.x != Seq_ID.y) %>% as.data.frame()
    
    
  }
  
  #
  # Left off checking the other mouse designs here::
  #

  # improbe_AL_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/mm10/mm10.improbeDesignOutput.core-v1.tsv.gz'
  # imp_des_AL_tib <- readr::read_tsv(improbe_AL_tsv)
  # improbe_AL_cols <- c('Probe_ID', 'Chrom1', 'Pos1', 'Pos2', 'Chip_Str1', 
  #                      'Customer_Score', 'Customer_Underlying_CpG_Count', 'Customer_Pos', 'Customer_Design_Char', 'Customer_Dist1', 'Customer_Dist2',
  #                      'Customer_Long_ID', 'Customer_ID_CO', 'Customer_Prb_A', 'Customer_Prb_B', 'Customer_Design_Type'
  # )
  
  mm10_imp_org_csv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/complicated-header.csv'
  mm10_imp_org_col <- suppressMessages(suppressWarnings( readr::read_csv(mm10_imp_org_csv) )) %>% names() %>% as.vector()
  # mm10_imp_org_col %>% length() = 92
  mm10_imp_ext_col <- c('Add_Forward_Sequence','Add_Genome_Build','Add_Chromosome','Add_Position','Add_Top_Sequence',
                    'Add_Probe_Seq_M','Add_FR','Add_TB','Add_CO','Add_Score_M','Add_CpgCnt','Add_CpgPos','Add_Probe_Seq_U',
                    'Add_Score_U')
  
  mm10_imp_AL_col <- c(mm10_imp_org_col,mm10_imp_ext_col)
  # mm10_imp_AL_col %>% length() = 106
  
  # mm10_imp_AL_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests/joined.designed.dat.head.csv'
  mm10_imp_AL_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests/joined.designed.dat.csv.gz'
  mm10_imp_AL_tib <- suppressMessages(suppressWarnings( readr::read_csv(mm10_imp_AL_csv, col_names=mm10_imp_AL_col) )) %>% 
    dplyr::mutate(
      Mat_Prb1=stringr::str_sub(UnMethyl_Probe_Sequence, 2,50), 
      Mat_Prb2=stringr::str_sub(UnMethyl_Probe_Sequence, 2,49) )
  
  # Join manifest probes and improbe based on Allele_A probe
  # mm10_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb1"))
  # mm10_man_tib %>% dplyr::inner_join(mm10_imp_AL_tib, by=c("Mat_Prb"="Mat_Prb2"))
  # mm10_man_tib %>% dplyr::full_join(mm10_imp_AL_tib, by=c("Seq_ID"="Probe_ID") ) %>% 
  #   dplyr::select(Seq_ID, PD, Address_A_Seq, Mat_Prb1, Mat_Prb2)

  
    
  #
  # Make sure both V1/AL have CGN->TOP Integrity::
  #
  unq_imp_v1_tib <- imp_des_v1_tib %>% dplyr::distinct(Seq_ID,Top_Sequence)
  top_cnt_v1_tib <- unq_imp_v1_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Top_Count=n())

  unq_imp_AL_tib <- improbe_AL_tib %>% dplyr::distinct(Probe_ID,Top_Sequence)
  top_cnt_AL_tib <- unq_imp_AL_tib %>% dplyr::group_by(Top_Sequence) %>% dplyr::summarise(Top_Count=n())
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                Preprocessing:: Improbe Mouse Designs
  #
  #  4. Comapre Previous Products with Mouse for CGN -> Top Integrity
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #              Dirty and Quick Selection of Optimal Content::
  #
  #                           EXAMPLE CODE BELOW::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (FALSE) {
    man_cn_des_all_tib <- dplyr::left_join(man_sel_cn_tib,
                                           imp_des_tib %>% dplyr::filter(Methyl_Allele_CO_Strand=='C') %>% dplyr::select(-Probe_Type),
                                           by=c("IlmnID"="Seq_ID"), suffix=c("_Man", "_Des")
    ) %>% dplyr::mutate(
      Mat_Prb1_Man=stringr::str_sub(AlleleA_ProbeSeq, 2,50) %>% stringr::str_replace_all('R', 'A'), 
      Mat_Prb2_Man=stringr::str_sub(AlleleA_ProbeSeq, 3,50) %>% stringr::str_replace_all('R', 'A'), 
      Mat_Prb1_Des=stringr::str_sub(UnMethyl_Probe_Sequence, 2,50), 
      Mat_Prb2_Des=stringr::str_sub(UnMethyl_Probe_Sequence, 2,49) ) %>% dplyr::filter(! IlmnID %in% man_27k_tib$IlmnID)
    
    mat_des1_tib <- man_cn_des_all_tib %>% dplyr::filter(Mat_Prb1_Man==Mat_Prb1_Des & Forward_Sequence_Man == Forward_Sequence_Des)
    mat_des2_tib <- man_cn_des_all_tib %>% dplyr::filter(Mat_Prb2_Man==Mat_Prb2_Des & Forward_Sequence_Man == Forward_Sequence_Des)
    
    man_cn_des_mat_tib <- dplyr::bind_rows( mat_des1_tib,mat_des2_tib ) %>% 
      dplyr::arrange(IlmnID) %>% dplyr::rename(Forward_Sequence=Forward_Sequence_Des) %>% 
      dplyr::select(-Forward_Sequence_Man)
    
    man_cn_des_mat_tib %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise(Count=n()) %>% print()
  }
  
  
}











# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Simplified Output Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

out_man_tib <- mm10_man_tib %>% 
  dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, 
                AlleleA_Probe_Sequence,AlleleB_Probe_Sequence) %>%
  dplyr::mutate(Version=opt$tar_version)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Extract Probe CGN Positions for cg & mu::
#                              NOT USED YET!
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  cgn_key_tib <- out_man_tib %>% dplyr::filter(Probe_Type=='mu' | Probe_Type=='cg') %>% 
    dplyr::mutate(CGN=stringr::str_remove(Probe_ID,'_.*$')) %>% 
    dplyr::mutate(CGN=stringr::str_replace(CGN,'^mu','cg')) %>% 
    dplyr::select(CGN,Probe_ID,Probe_Type) %>% dplyr::arrange(CGN)
  
  cgn_key_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Old Control Probe Loading::
#                              NOT USED YET!
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  opt$manifestPath <- file.path(opt$manDir, 'HM450-B2.manifest.sesame-base.cpg-sorted.csv.gz')
  if (is.null(opt$manifestPath))
    opt$manifestPath <- file.path(opt$manDir, paste0(opt$platform,'-',opt$manifest,'.manifest.sesame-base.cpg-sorted.csv.gz') )
  
  full_org_man_tib <- loadManifestSource(opt$manifestPath, verbose=opt$verbose,tt=NULL)
  org_man_tib <- loadManifestSource(opt$manifestPath, verbose=opt$verbose,tt=NULL) %>%
    dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base) %>% 
    dplyr::mutate(Version=opt$manifest)
  org_ctl_tib <- org_man_tib %>% dplyr::filter(Probe_Type!='cg', Probe_Type!='ch', Probe_Type!='rs') %>% dplyr::arrange(Probe_ID)
  
  if (opt$addControls) out_man_tib <- dplyr::bind_rows(out_man_tib, org_ctl_tib) %>% dplyr::arrange(Probe_ID, DESIGN)
  if (opt$addManifest) out_man_tib <- dplyr::bind_rows(out_man_tib, org_man_tib) %>% dplyr::arrange(Probe_ID, DESIGN)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Sesame Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  ses_neg_tib <- ses_ctl_tib %>% dplyr::filter(Probe_Type=='NEGATIVE')
  # ses_man_tib <- dplyr::bind_rows(out_man_tib) %>% 
  # ses_man_tib <- dplyr::bind_rows(out_man_tib,ses_ctl_tib) %>% 
  ses_man_tib <- dplyr::bind_rows(out_man_tib,ses_neg_tib) %>% 
    # dplyr::select(-AlleleA_Probe_Sequence,-AlleleB_Probe_Sequence) %>%
    dplyr::arrange(Probe_ID)
  readr::write_csv(ses_man_tib, ses_man_csv)
  
  # Local Mac Copy Command For Testing::
  #  cp /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-cp/mm10_LEGX_cp.manifest.sesame-base.cpg-sorted.csv.gz tools/Infinium_Methylation_Workhorse/dat/manifest/base/LEGX-B0..manifest.sesame-base.cpg-sorted.csv.gz
  #
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Genome Studio Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loci_cnt <- out_man_tib %>% base::nrow()
gs_head_tib <- tibble::tibble(Key=character(), Value=character())
gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="Illumina, Inc.", Value="")
gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="[Heading]", Value="")
gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="Descriptor File Name", Value=paste0('Infinium_Methylation_',opt$outName,'.csv.gz'))
gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="Assay Format", Value="Infinium 2")
gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="Date Manufactured", Value=as.character(Sys.Date()))
gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="Loci Count", Value=as.character(loci_cnt))
gs_head_tib <- gs_head_tib %>% tibble::add_row(Key="[Assay]", Value="")

#  IlmnID,Name,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,
#   Infinium_Design_Type,Next_Base,Color_Channel,Forward_Sequence,
#   Genome_Build,CHR,MAPINFO,SourceSeq,Strand
#
#   Missing: Source_Seq
#
src_seq <- as.character(paste(rep('N',50), collapse='') )
gs_body_tib <- out_man_tib %>% # head(n=1000) %>%
  dplyr::rename(IlmnID=Probe_ID,
                AddressA_ID=U,AlleleA_ProbeSeq=AlleleA_Probe_Sequence,
                AddressB_ID=M,AlleleB_ProbeSeq=AlleleB_Probe_Sequence,
                Infinium_Design_Type=DESIGN,Color_Channel=col) %>%
  dplyr::mutate(Name=IlmnID,Forward_Sequence='N',Genome_Build=opt$genomeBuild,
                CHR='chr1',MAPINFO=1,Strand='F',Source_Seq=!!src_seq) %>%
  dplyr::select(IlmnID,Name,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,
                Infinium_Design_Type,Next_Base,Color_Channel,Forward_Sequence,
                Genome_Build,CHR,MAPINFO,Source_Seq,Strand,Source_Seq,Probe_Type,Probe_Source) %>%
  dplyr::arrange(IlmnID)
ctl_head_df <- tibble::tibble(Head="[Controls]")

# Silly addition of commas to the end of the controls section...
#
org_ctl_tib <- org_ctl_tib %>% dplyr::mutate(Probe_ID=paste0(Probe_ID,",,,,,,"))

# Silly incosistent padding of Tango Addresses..
gs_body_tib <- gs_body_tib %>% dplyr::mutate(
  AddressA_ID=stringr::str_pad(string=AddressA_ID, pad='0', side="left", width=10), 
  AddressB_ID=stringr::str_pad(string=AddressB_ID, pad='0', side="left", width=10),
  Color_Channel=dplyr::case_when(
    Color_Channel=='R' ~ 'Red',
    Color_Channel=='G' ~ 'Grn',
    TRUE ~ ''
  )
)

# Write Genome Studio CSV::
write_delim(x=gs_head_tib, path=gss_man_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=FALSE)
write_delim(x=gs_body_tib, path=gss_man_csv, delim=',', col_names=TRUE,  quote_escape=FALSE, na='', append=TRUE)
write_delim(x=ctl_head_df, path=gss_man_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)
write_delim(x=org_ctl_tib, path=gss_man_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)

# Follow up silly command to remove quotes::
#  gzip -dc /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-B4/mm10_LEGX_B4.manifest.GenomeStudio.cpg-sorted.csv.gz | perl -pe 's/\"//gi' | gzip -c - > /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-B4/mm10_LEGX_B4.manifest.GenomeStudio.cpg-sorted.clean.csv.gz

# gss_man_zip <- paste(gss_man_csv,'gz', sep='.')
# cmd <- glue::glue("cat {gss_man_csv} | perl -pe 's/\"//gi' | gzip -c - > {gss_man_zip}")
# system(cmd)

# Follow up silly command to remove quotes::
#  cat /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-B5/mm10_LEGX_B5.manifest.GenomeStudio.cpg-sorted.csv | perl -pe 's/\"//gi' | gzip -c - > /Users/bbarnes/Documents/Projects/methylation/scratch/manifests/mm10-LEGX-B5/mm10_LEGX_B5.manifest.GenomeStudio.cpg-sorted.clean.csv.gz

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Write mm10 Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$make_addresss) {
  out_add_tib <- NULL
  out_add_tib <- dplyr::bind_rows(
    out_man_tib %>% dplyr::filter(!is.na(M)) %>% dplyr::select(Probe_ID, M, col, DESIGN, Probe_Type) %>% dplyr::rename(Address=M),
    out_man_tib %>% dplyr::filter(!is.na(U)) %>% dplyr::select(Probe_ID, U, col, DESIGN, Probe_Type) %>% dplyr::rename(Address=U)
  ) %>% dplyr::rename(Man_Col=col, Design_Type=DESIGN)
  
  readr::write_csv(out_add_tib, out_add_csv)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Load Annotation/Alignment Results::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  ann_cnames <- c('Probe_ID', 'Probe_ID_A', 'Probe_Seq_A', 'Probe_ID_B', 'Probe_Seq_B', 
                  'Normalization_Bin', 'Probe_Chrom', 'Probe_Beg', 'Probe_End', 'Strand_FR', 'Probe_Cnt',
                  'Reg_Type', 'Transcript', 'Gene_Name', 'Reg_Beg', 'Reg_End', 'Reg_Idx', 
                  'Gene_Chrom', 'Gene_Beg', 'Gene_End', 'Gene_Score')
  ann1_tsv <- file.path(opt$alnDir, 'COVID_EPIC_Round1.COVID-Genes.intersection.tsv')
  ann1_tib <- readr::read_tsv(file=ann1_tsv, col_names=ann_cnames)
  
  gene_ann_tib <- ann1_tib %>% dplyr::rename(Aln_Probe_ID=Probe_ID) %>% 
    dplyr::mutate(Probe_ID=stringr::str_remove(Probe_ID_A, '_A$')) %>% 
    tidyr::separate(Probe_ID, into=c("CG","FR","TB","CO","DT","Hits"), sep='_', remove=FALSE) %>%
    dplyr::select(Probe_ID, CG, FR, TB, CO, DT, Hits, everything())
  
  hits_ann_tib <- gene_ann_tib %>% dplyr::select(Probe_ID:DT, Probe_Chrom, Probe_Beg, Probe_End, Strand_FR) %>% dplyr::distinct() %>% 
    dplyr::add_count(Probe_ID, name='Aln_Count') %>% dplyr::arrange(-Aln_Count) %>% dplyr::distinct(Probe_ID, Aln_Count)
  
  gene_man_tib <- full_man_tib %>% dplyr::inner_join(gene_ann_tib, by="Probe_ID", suffix=c("_Man", "_Ann"))
  misA_man_tib <- full_man_tib %>% dplyr::anti_join(gene_ann_tib, by="Probe_ID", suffix=c("_Man", "_Ann"))
  misB_man_tib <- gene_ann_tib %>% dplyr::anti_join(full_man_tib, by="Probe_ID", suffix=c("_Man", "_Ann"))
  
  # full_man_tib %>% dplyr::filter(stringr::str_starts( 'cg00065388' ) )
  
  
  #
  # Quick look up for Fasial for chr18 HLA region coverage "ESTIMATE"
  #
  #  hg19: chr5 140164938 140864474
  #  mm10: chr18 36930065 37815275
  #
  #
  # mm_tib <- readr::read_tsv('/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/to-order/designW.Final_Dec-10-2019.cpg-only/Mus_musculus.annotation.genomic.sorted-cgn.chr-pos.tsv', col_names = c("Probe_ID",'Pos'))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Genome Studio Controls Swap::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

  color_tib <- suppressMessages(suppressWarnings( readr::read_csv(file.path(par$datDir, 'params/GenomeStudioColors.csv')) ))
  color_vec <- color_tib %>% dplyr::pull(Color) %>% as.vector()
  color_len <- length(color_vec)
  
  man_gs_csv  <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests/mm10-LEGX-B5/mm10_LEGX_B5.manifest.GenomeStudio.cpg-sorted.clean.csv.gz'
  man_gs_list <- loadManifestGenomeStudio(file = man_gs_csv, addSource = TRUE, normalize = FALSE, verbose = 20)
  man_head_df <- readr::read_lines(man_gs_csv, n_max = 6) %>% as.data.frame()
  ctl_head_df <- tibble::tibble(Head="[Controls]") %>% tibble::add_column(BL1='',BL2='',BL3='',BL4='') %>% as.data.frame()
  
  man_gs_list$ctl %>% dplyr::group_by(Control_Group) %>% dplyr::summarise(Count=n())
  man_gs_list$man %>% dplyr::group_by(Probe_Type, Infinium_Design_Type) %>% dplyr::summarise(Type_Count=n())
  
  man_ct_tib <- man_gs_list$man %>% dplyr::filter(Probe_Type=='BS' | Probe_Type=='NO')
  man_gs_tib <- man_gs_list$man %>% dplyr::filter(Probe_Type!='BS') %>% dplyr::filter(Probe_Type!='NO')
  ctl_gs_tib <- man_gs_list$ctl
  
  ctl_cols <- ctl_gs_tib %>% dplyr::select(1:4) %>% names()
  
  # Naming investigation::
  ctl_gs_tib %>% dplyr::filter(Control_Group=="BISULFITE CONVERSION I") %>% dplyr::arrange(Control_Name) %>% dplyr::select(-Control_Color, -Man_Source) %>% as.data.frame()  
  ctl_gs_tib %>% dplyr::filter(Control_Group=="BISULFITE CONVERSION II") %>% dplyr::arrange(Control_Name) %>% dplyr::select(-Control_Color, -Man_Source) %>% as.data.frame()
  ctl_gs_tib %>% dplyr::filter(Control_Group=="NON-POLYMORPHIC") %>% dplyr::arrange(Control_Name) %>% dplyr::select(-Control_Color, -Man_Source) %>% as.data.frame()
  
  # Make new controls::
  new_ct_list <- man_ct_tib %>% dplyr::arrange(Probe_Type) %>%
    dplyr::mutate(
      Control_Group=dplyr::case_when(
        Probe_Type=='BS' & Infinium_Design_Type=='I'  ~ "BISULFITE CONVERSION I",
        Probe_Type=='BS' & Infinium_Design_Type=='II' ~ "BISULFITE CONVERSION II",
        Probe_Type=='NO' ~ 'NON-POLYMORPHIC',
        TRUE ~ NA_character_
      ),
      DiNuc=dplyr::case_when(
        Probe_Type=='BS' ~ stringr::str_replace(IlmnID, '^.*-([ACTG][ACTG])-.*$', '\\$1') %>% stringr::str_remove_all('\\\\'),
        Probe_Type=='NO' ~ stringr::str_replace(IlmnID, '^.*_([ACTG][ACTG])_.*$', '\\$1') %>% stringr::str_remove_all('\\\\'),
        TRUE ~ NA_character_
      ),
      N1=stringr::str_sub(DiNuc, 1,1),
      N2=stringr::str_sub(DiNuc, 2,2),
      
      Control_ID=paste(IlmnID, stringr::str_sub(Color_Channel, 1,1), sep='_')
    ) %>% split(.$Infinium_Design_Type)
    # dplyr::select(IlmnID,DiNuc,N1,N2,Control_Group,Control_ID,Color_Channel,Probe_Type,Infinium_Design_Type)
  
  # Notes:: 
  #  BISULFITE CONVERSION I; All M's get 'SkyBlue' the C's get a unique color
  #     N2.B == A(G) => C
  #     N2.B == T(R) -> U
  #     N1.A == C    -> M
  #
  bs1_both_tib <- new_ct_list$I %>% dplyr::mutate(
    Row_Idx=dplyr::row_number() + 100,
    Row_Str=Row_Idx,
    Col_Idx=Row_Idx %% (color_len-100),
    
    Control_Name_A=paste0('BS Conversion I M',Row_Str),
    Control_Name_B=dplyr::case_when(
      Color_Channel=='Grn' ~ paste0('BS Conversion I-C',Row_Str),
      Color_Channel=='Red' ~ paste0('BS Conversion I-U',Row_Str),
      TRUE ~ NA_character_
    ),
    Control_ColorA='Red',
    Control_ColorB=color_vec[Col_Idx+1]
  )
  
  # %>% dplyr::select(Control_Name_A,Control_Name_B,
  #                     IlmnID,DiNuc,N1,N2,Control_Group,Row,Control_ID,
  #                     Color_Channel,Probe_Type,Infinium_Design_Type) %>% dplyr::select(1,2) %>% as.data.frame()

  fin_bs1_tib <- dplyr::bind_rows(
    bs1_both_tib %>% dplyr::select(AddressA_ID,Control_Group,Control_ColorA,Control_Name_A) %>% purrr::set_names(ctl_cols),
    bs1_both_tib %>% dplyr::select(AddressB_ID,Control_Group,Control_ColorB,Control_Name_B) %>% purrr::set_names(ctl_cols)
  ) %>% dplyr::arrange(Control_Name)

  fin_bs2_tib <- new_ct_list$II %>% 
    dplyr::filter(Control_Group=='BISULFITE CONVERSION II') %>%
    dplyr::mutate(Row_Idx=row_number()+100,
                  Row_Str=Row_Idx,
                  Control_Color='Blue',
                  Control_Name=paste('BS Conversion II',Row_Str, sep='-')
    ) %>% 
    dplyr::select(AddressA_ID,Control_Group,Control_Color,Control_Name) %>% 
    purrr::set_names(ctl_cols) %>%
    dplyr::arrange(Control_Name)
  
  fin_non_tib <- new_ct_list$II %>% 
    dplyr::filter(Control_Group=='NON-POLYMORPHIC') %>%
    dplyr::mutate(Row_Idx=row_number()+100,
                  Row_Str=Row_Idx,
                  Control_Color='Blue',
                  Control_Name=paste0('NP (',N1,') ',Row_Str)
    ) %>% 
    dplyr::select(AddressA_ID,Control_Group,Control_Color,Control_Name) %>% 
    purrr::set_names(ctl_cols) %>%
    dplyr::arrange(Control_Name)
  
  fin_hum_tib <- ctl_gs_tib %>% dplyr::select(1:4) %>% 
    purrr::set_names(ctl_cols) %>% dplyr::arrange(Control_Name) %>%
    dplyr::mutate(Address=as.double(Address))
  
  fin_ctl_tib <- dplyr::bind_rows(fin_bs1_tib,fin_bs2_tib,fin_non_tib,fin_hum_tib)
  
  #
  # Investigat Duplicates::
  #
  add_cnt_tib <- dplyr::bind_rows(
    dplyr::select(man_gs_tib,AddressA_ID) %>% dplyr::rename(Address=AddressA_ID), 
    dplyr::select(man_gs_tib,AddressB_ID) %>% dplyr::rename(Address=AddressB_ID),
    dplyr::select(fin_ctl_tib,Address)
  ) %>% dplyr::filter(!is.na(Address)) %>%
    dplyr::group_by(Address) %>% dplyr::summarise(Add_Dup_Count=n()) %>% dplyr::arrange(-Add_Dup_Count)
    
  add_tib_tib <- add_cnt_tib %>% dplyr::filter(Add_Dup_Count>1)
  
  man_fin_tib <- man_gs_tib %>% dplyr::filter(!AddressA_ID %in% add_tib_tib$Address) %>% dplyr::filter(!AddressB_ID %in% add_tib_tib$Address)
  ctl_fin_tib <- fin_ctl_tib %>% dplyr::filter(! Address %in% add_tib_tib$Address) %>% 
    tibble::add_column(BL1='',BL2='',BL3='',BL4='',BL5='',BL6='')
    # tibble::add_column(BL1='',BL2='',BL3='',BL4='',BL5='',BL6='',BL7='',BL8='')
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Write Control Swapped Manifest::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

  gs_man_ver  <- 8
  gs_man_ver  <- 9
  gs_swap_dir <- file.path('/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/manifests', paste0('mm10-LEGX-B',gs_man_ver))
  gs_swap_csv <- file.path(gs_swap_dir, paste0('mm10_LEGX_B',gs_man_ver,'.manifest.GenomeStudio.cpg-sorted.clean.csv') )
  gz_swap_csv <- file.path(gs_swap_dir, paste0('mm10_LEGX_B',gs_man_ver,'.manifest.GenomeStudio.cpg-sorted.clean.csv.gz') )
  
  if (!dir.exists(gs_swap_dir)) dir.create(gs_swap_dir, recursive=TRUE)

  man_head_replace_str <- paste0('_B',gs_man_ver,'.csv.gz')
  man_head_df <- man_head_df %>% as_tibble() %>% purrr::set_names(c("Head_Var")) %>% 
    dplyr::mutate(Head_Var=as.character(Head_Var),
                  Head_Var=stringr::str_replace(Head_Var, '_B[0-9]+.csv.gz$', man_head_replace_str) ) %>%
    as.data.frame()
  
  # Write Genome Studio CSV::
  write_delim(x=man_head_df, path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=FALSE)
  write_delim(x=man_fin_tib, path=gs_swap_csv, delim=',', col_names=TRUE,  quote_escape=FALSE, na='', append=TRUE)
  write_delim(x=ctl_head_df, path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)
  write_delim(x=ctl_fin_tib, path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)
  
  cmd <- glue::glue("cat {gs_swap_csv} | perl -pe 's/\"//gi' | gzip -c - > {gz_swap_csv}")
  # system(cmd)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))


# End of file
