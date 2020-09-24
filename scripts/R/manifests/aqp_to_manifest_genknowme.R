
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
par <- NULL
opt <- NULL

par$runMode <- ''
par$macDir1 <- NULL
par$macDir2 <- NULL
par$lixDir1 <- NULL

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
opt$parallel <- FALSE
opt$cluster  <- FALSE

# verbose Options::
opt$verbose <- 3

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Parse Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
args.dat <- commandArgs(trailingOnly = FALSE)
if (args.dat[1]=='RStudio') {
  # Illumina based directories::
  par$macDir1 <- '/Users/bbarnes/Documents/Projects/methylation/tools'
  par$macDir2 <- '/Users/bretbarnes/Documents/tools'
  par$lixDir1 <- '/illumina/scratch/darkmatter'
  
  par$runMode    <- args.dat[1]
  cat(glue::glue("[{par$prgmTag}]: Local Run args.dat[1]={args.dat[1]}.{RET}"))
  cat(glue::glue("[{par$prgmTag}]: Local Run     runMode={par$runMode}.{RET}"))
  
  if (dir.exists(par$macDir1)) par$topDir <- '/Users/bbarnes/Documents/Projects/methylation'
  if (dir.exists(par$macDir2)) par$topDir <- '/Users/bretbarnes/Documents'
  if (!dir.exists(par$topDir)) dir.create(par$topDir, recursive=TRUE)
  
  if (dir.exists(par$macDir1)) par$macDir <- par$macDir1
  if (dir.exists(par$macDir2)) par$macDir <- par$macDir2
  
  # Default Parameters for local Mac::
  par$srcDir     <- file.path(par$macDir, par$codeDir)
  par$scrDir     <- file.path(par$srcDir, 'scripts')
  par$exePath    <- file.path(par$scrDir, 'R', par$prgmDir, paste0(par$prgmTag,'.R'))
  
  par$prgmTag <- base::sub('\\.R$', '', base::basename(par$exePath))
  par$locPath <- base::dirname(par$exePath)
  par$scrDir  <- base::dirname(base::normalizePath(par$locPath) )
  par$srcDir  <- base::dirname(base::normalizePath(par$scrDir) )
  par$datDir  <- file.path(base::dirname(base::normalizePath(par$srcDir)), 'dat')
  
  opt$outDir <- file.path(par$topDir, 'scratch')
  locIdatDir <- file.path(par$topDir, 'data/idats')
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  
  # Platform/Method Options::
  opt$genomeBuild <- 'hg38'
  opt$platform    <- 'GENK'
  opt$version     <- 'B3'
  opt$version     <- 'B4'
  opt$version     <- 'A0'
  
  opt$frmt_original <- TRUE
  opt$frmt_original <- FALSE
  
  opt$write_full <- FALSE
  opt$write_base <- FALSE
  
  opt$make_addresss <- FALSE
  
  opt$addControls <- TRUE
  opt$addManifest <- FALSE
  
  # opt$aqpDir <- file.path('/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics', 'AQP')
  opt$aqpDir <- '/Users/bbarnes/Documents/Projects/methylation/CustomContent/Genknowme/LS_Epiprofile'
  
  opt$ords <- paste(
    file.path(opt$aqpDir, 'AQP1_NAremoved_GenKnowme_CpG_SNP_order.07082020.csv'),
    file.path(opt$aqpDir, 'AQP2_AP_Genknowme_AQP2_replicate_design_file2.csv'),
    sep=',')
  
  opt$mats <- paste(
    file.path(opt$aqpDir, '20468029_AQP1_probes.match'),
    file.path(opt$aqpDir, '20468029_AQP2_probes.match'),
    sep=',')
  
  opt$aqps <- paste(
    file.path(opt$aqpDir, 'BS0032678_AQP1-AQP.txt'),
    file.path(opt$aqpDir, 'BS0032779_AQP2-AQP.txt'),
    sep=',')
  
  opt$pqcs <- paste(
    file.path(opt$aqpDir, '20042793X371678_A_ProductQC_AP.txt'),
    sep=',')
  
  par$idatsTopDir <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/idats/ILMN_mm10_betaTest_17082020'
  opt$idats <- paste(
    file.path(par$idatsTopDir, '204637490002'),
    sep=',')
  opt$idats <- NULL

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
    # is.null(opt$idats) ||
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
opt$frmt_original <- TRUE
isFull <- FALSE

ord_cnt <- length(ord_vec)

full_man_tib <- NULL
if (!is.null(opt$pqcs)) {
  cat(glue::glue("[{par$prgmTag}]: Running with PQC...{RET}") )
  
  for (idx in c(1:ord_cnt)) {
    full_man_tib <- full_man_tib %>% dplyr::bind_rows(
      decodeToManifest(ord=ord_vec[idx], mat=mat_vec[idx], pqc=pqc_vec[1], ordSkip=0, matSkip=0, 
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
# INVESTIGATION PART:: this can be deleted now
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

ses_man_tib <- fixOrderProbeIDs(full_man_tib, verbose=opt$verbose)
ses_man_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()

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
#                           Map Back to Designs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# This is just Seq_ID
# ses_man_tib %>% dplyr::mutate(CGN=stringr::str_remove(Probe_ID, '_.*$'))

EPIC_tab_csv <- '/Users/bbarnes/Documents/Projects/manifests/methylation/MethylationEPIC_v-1-0_B4.core.cpg-only.table.csv.gz'
EPIC_tab_tib <- suppressMessages(suppressWarnings( readr::read_csv(EPIC_tab_csv) ))

ses_man_tib %>% dplyr::left_join(EPIC_tab_tib, by=c("Seq_ID"="Name"), suffix=c("_Ses", "_EPIC"))

ses_man_tib %>% dplyr::left_join(EPIC_tab_tib, by=c("Seq_ID"="Name"), suffix=c("_Ses", "_EPIC")) %>% dplyr::filter(is.na(IlmnID))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Simplified Output Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

out_man_tib <- ses_man_tib %>% 
  dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, 
                AlleleA_Probe_Sequence,AlleleB_Probe_Sequence) %>%
  dplyr::mutate(Version=opt$tar_version)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Sesame Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

out_man_tib <- dplyr::bind_rows(out_man_tib,ses_ctl_tib) %>% dplyr::arrange(Probe_ID)
readr::write_csv(out_man_tib, ses_man_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))


# End of file
