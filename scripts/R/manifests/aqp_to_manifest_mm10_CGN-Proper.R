
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

# Executables::
opt$Rscript <- NULL

# Directories::
opt$outDir <- NULL
opt$impDir <- NULL

# Required Inputs::
opt$ords <- NULL
opt$mats <- NULL
opt$aqps <- NULL
opt$pqcs <- NULL
opt$ctls <- NULL
opt$idat <- NULL

# Platform/Method Options::
opt$genomeBuild <- NULL
opt$platform    <- NULL
opt$version     <- NULL

# Run Options::
opt$fresh <- FALSE
par$retData <- FALSE
opt$matFormat <- 'new'

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
  opt$version     <- 'S6'
  opt$version     <- 'C0'
  opt$version     <- 'C1'
  opt$version     <- 'C2'
  opt$version     <- 'C4'
  opt$version     <- 'C5'
  opt$version     <- 'C6'
  
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
    opt$idat <- paste(
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
  opt$impDir <- file.path(par$topDir, 'data/improbe')
  
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
    # Executables::
    make_option(c("--Rscript"), type="character", default=opt$Rscript, 
                help="Rscript path [default= %default]", metavar="character"),
    
    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    make_option(c("--impDir"), type="character", default=opt$impDir, 
                help="improbe data directory [default= %default]", metavar="character"),
    
    # Pre-defined files (controls)
    make_option(c("--ords"), type="character", default=opt$ords, 
                help="Order files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mats"), type="character", default=opt$mats, 
                help="Match files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqps"), type="character", default=opt$aqps, 
                help="AQP files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--pqcs"), type="character", default=opt$pqcs, 
                help="PQC files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--ctls"), type="character", default=opt$ctls, 
                help="Pre-Defined Infinium Methylation Controls (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--idat"), type="character", default=opt$idat, 
                help="idat directories (comma seperated) [default= %default]", metavar="character"),
    
    # Platform/Method Options::
    make_option(c("--genomeBuild"), type="character", default=opt$genomeBuild, 
                help="Genome Build (e.g. hg18, hg36, hg19, hg37, hg38, mm10) [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform (e.g. HM450, EPIC, LEGX, NZT, COVIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest Version (e.g. B0,B1,B2,B3,B4,C0) [default= %default]", metavar="character"),
    
    # Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
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

if (is.null(opt$outDir) || is.null(opt$impDir) ||
    is.null(opt$ords) || is.null(opt$mats) || 
    is.null(opt$genomeBuild) || is.null(opt$platform) || is.null(opt$version) ||
    is.null(opt$Rscript) ||
    is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )
  
  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  
  if (is.null(opt$outDir))    cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$impDir))    cat(glue::glue("[Usage]: impDir is NULL!!!{RET}"))
  
  if (is.null(opt$ords)) cat(glue::glue("[Usage]: ords is NULL!!!{RET}"))
  if (is.null(opt$mats)) cat(glue::glue("[Usage]: mats is NULL!!!{RET}"))
  if (is.null(opt$aqps)) cat(glue::glue("[Usage]: aqps is NULL [Not required if pqcs defined]!!!{RET}"))
  if (is.null(opt$pqcs)) cat(glue::glue("[Usage]: pqcs is NULL [Not required if aqps defined]!!!{RET}"))
  if (is.null(opt$ctls)) cat(glue::glue("[Usage]: ctls is NULL [Not required]!!!{RET}"))
  if (is.null(opt$idat)) cat(glue::glue("[Usage]: idat is NULL [Not Required]!!!{RET}"))
  
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

opt$probe_type <- 'cg'
opt$design_key <- 'Seq_ID'
opt$design_seq <- 'Forward_Sequence'
opt$design_seq <- 'Top_Sequence'
opt$design_prb <- 'Probe_Type'
opt$design_srs <- 'TB'
opt$design_cos <- 'CO'

design_prb_sym <- rlang::sym(opt$design_prb)


# Input Definitions::
if (is.null(opt$ctls)) {
  opt$ctls <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
}


ords_vec <- NULL
mats_vec <- NULL
aqps_vec <- NULL
pqcs_vec <- NULL
ctls_vec <- NULL
idat_vec <- NULL
if (!is.null(opt$ords)) ords_vec <- opt$ords %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$mats)) mats_vec <- opt$mats %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$aqps)) aqps_vec <- opt$aqps %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$pqcs)) pqcs_vec <- opt$pqcs %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$ctls)) ctls_vec <- opt$ctls %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$idat)) idat_vec <- opt$idat %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

stopifnot(length(ords_vec)>0)
stopifnot(length(mats_vec)>0)
stopifnot(length(mats_vec)==length(ords_vec))

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

opt$prdDir <- file.path(opt$outDir, 'product')
if (!dir.exists(opt$prdDir)) dir.create(opt$prdDir, recursive=TRUE)

cat(glue::glue("[{par$prgmTag}]: Done. Building Output Directories.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Genome Studio Color Codes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

color_tib <- suppressMessages(suppressWarnings( readr::read_csv(file.path(par$datDir, 'params/GenomeStudioColors.csv')) ))
color_vec <- color_tib %>% dplyr::pull(Color) %>% as.vector()
color_len <- length(color_vec)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Process AQP/PQC Data::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$fresh <- TRUE
opt$verbose <- 4

man_raw_tib <- decodeAqpPqcWrapper(
  ord_vec=ords_vec,mat_vec=mats_vec,aqp_vec=aqps_vec,pqc_vec=pqcs_vec, 
  platform=opt$platform,version=opt$version,matFormat=opt$matFormat,
  name=opt$runName,outDir=opt$manDir,fresh=opt$fresh,full=par$retData,trim=TRUE,
  verbose=opt$verbose,vt=3,tc=0,tt=pTracker)

if (FALSE) {
  man_raw_tib %>% dplyr::filter(U==59792834) %>% as.data.frame()
  man_raw_tib %>% dplyr::filter(M==59792834) %>% as.data.frame()
  
  man_raw_tib %>% base::nrow()
  man_raw_tib %>% dplyr::group_by(U,M) %>% base::nrow()
  man_raw_tib %>% dplyr::group_by(U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence) %>% base::nrow()
  
  add_cnt_tib <- dplyr::bind_rows(
    dplyr::select(man_raw_tib, U) %>% dplyr::rename(Address=U) %>% dplyr::filter(!is.na(Address)),
    dplyr::select(man_raw_tib, M) %>% dplyr::rename(Address=M) %>% dplyr::filter(!is.na(Address))
  ) %>% dplyr::add_count(Address, name="Address_Count")
  
  add_fail_cnt <- add_cnt_tib %>% dplyr::filter(Address_Count!=1)
  
  man_raw_tib %>% dplyr::filter(U %in% add_fail_cnt$Address) %>% dplyr::group_by(Rep_Max) %>% dplyr::summarise(Rep_Count=n(), .groups='drop')
  man_raw_tib %>% dplyr::filter(M %in% add_fail_cnt$Address) %>% dplyr::group_by(Rep_Max) %>% dplyr::summarise(Rep_Count=n(), .groups='drop')
  
  man_raw_tib %>% dplyr::filter(U %in% add_fail_cnt$Address)
  man_raw_tib %>% dplyr::filter(M %in% add_fail_cnt$Address) %>% dplyr::filter(Rep_Max==1)
}

opt$fresh <- FALSE

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
#
#                        1. RawMan -> seq48U::
#                        2. seq48U -> CGN/TB/CO::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imp_s48_tsv <- file.path(opt$impDir, 'designOutput_21092020/seq48U/seq48U-GRCh36-38-10-21092020.unq.noHeader.seq-sorted.csv.gz')
imp_s48_tsv <- file.path(opt$impDir, 'designOutput_21092020/seq48U/un/seq48U-GRCh36-38-10-21092020.unq.noHeader.seq-sorted.tsv')

imp_s48_tib <- seq48U_to_cgn(tib=man_raw_tib, imp_tsv=imp_s48_tsv,
                             name=opt$runName,outDir=opt$intDir,
                             mat_key='Mat_PrbA', fresh=opt$fresh,
                             colA=4,colB=1,
                             verbose=opt$verbose,tc=0,tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             3. CGN to Top::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# This is for all species::
if (FALSE) {
  imp_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/GRCh36-GRCh38-GRCm10-21092020.cgnTop.sorted.tsv.gz')
  imp_top_tib <- cgn_to_topSeq(tib=imp_s48_tib, imp_tsv=imp_top_tsv,
                               name=opt$runName,outDir=opt$intDir,
                               mat_key='CGN_Imp', fresh=opt$fresh,
                               colA=1,colB=1,
                               verbose=opt$verbose,tc=0,tt=pTracker) %>%
    dplyr::mutate(!!design_prb_sym := opt$probe_type)
}

# Best to use species specific::
imp_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop/GRCm10-21092020.cgnTop.sorted.tsv.gz')
imp_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop/GRCm10-21092020.cgnTop.sorted.tsv')
imp_top_tib <- cgn_to_topSeq(tib=imp_s48_tib, imp_tsv=imp_top_tsv,
                             name=opt$runName,outDir=opt$intDir,
                             mat_key='CGN_Imp', fresh=opt$fresh,
                             colA=1,colB=1,
                             verbose=opt$verbose,tc=0,tt=pTracker) %>%
  dplyr::mutate(!!design_prb_sym := opt$probe_type)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           4. Top to Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_prb_rds <-
  file.path(opt$desDir, paste(opt$runName,opt$design_seq,opt$probe_type,'mat-probes.rds', sep='.') )
man_full_rds <-
  file.path(opt$desDir, paste(opt$runName,opt$design_seq,opt$probe_type,'all-probes.rds', sep='.') )

full_prb_des_tib <- NULL
if (!file.exists(man_prb_rds) || opt$fresh) {
  cat(glue::glue("[{par$prgmTag}]: Building Full Probe Design...{RET}"))
  
  test_max <- 0
  full_prb_des_tib <- 
    desSeq_to_prbs(tib=imp_top_tib, 
                   idsKey=opt$design_key,
                   seqKey=opt$design_seq,prbKey=opt$design_prb,
                   strsSR=opt$design_srs,strsCO=opt$design_cos,
                   parallel=opt$parallel, max=test_max,
                   verbose=opt$verbose,tc=1,tt=pTracker)
  
  full_prb_des_tib %>% dplyr::group_by(SR_Str,CO_Str) %>% 
    dplyr::summarise(SRD_Cnt=n(), .groups='drop') %>% print()
  
  # Remove:: only supports intermediates...
  cat(glue::glue("[{par$prgmTag}]: Writing Full Probe Match Design File: RDS={man_full_rds}...{RET}"))
  readr::write_rds(full_prb_des_tib,man_full_rds)
  cat(glue::glue("[{par$prgmTag}]: Done. Writing Full Probe Design File.{RET}{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        5. Match to Manifest Probes::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  man_prb_tib <- 
    man_join_prbs(man_raw_tib,full_prb_des_tib,
                  man_mat1U_key="AlleleA_Probe_Sequence", prb_mat1U_key="PRB1_U_MAT",
                  man_mat1M_key="AlleleB_Probe_Sequence", prb_mat1M_key="PRB1_M_MAT", 
                  man_mat2D_key="AlleleA_Probe_Sequence", prb_mat2D_key="PRB2_D_MAT",
                  verbose=opt$verbose+10,tc=1,tt=pTracker)
  
  cat(glue::glue("[{par$prgmTag}]: Writing Join Probe Match Design File: RDS={man_prb_rds}...{RET}"))
  readr::write_rds(man_prb_tib,man_prb_rds)
  cat(glue::glue("[{par$prgmTag}]: Done. Writing Join Probe Design File.{RET}{RET}"))
} else {
  cat(glue::glue("[{par$prgmTag}]: Loading Full Match Probe Design; RDS={man_prb_rds}...{RET}"))
  man_prb_tib <- readr::read_rds(man_prb_rds)
  cat(glue::glue("[{par$prgmTag}]: Done. Loading Join Probe Design.{RET}{RET}"))
}

if (FALSE) {
  man_prb_tib %>% dplyr::filter(U==59792834) %>% as.data.frame()
  man_prb_tib %>% dplyr::filter(M==59792834) %>% as.data.frame()
  
  man_prb_tib %>% base::nrow()
  man_prb_tib %>% dplyr::group_by(U,M) %>% base::nrow()
  man_prb_tib %>% dplyr::group_by(U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence) %>% base::nrow()
  
  add_cnt_tib <- dplyr::bind_rows(
    dplyr::select(man_prb_tib, U) %>% dplyr::rename(Address=U) %>% dplyr::filter(!is.na(Address)),
    dplyr::select(man_prb_tib, M) %>% dplyr::rename(Address=M) %>% dplyr::filter(!is.na(Address))
  ) %>% dplyr::add_count(Address, name="Address_Count")
  
  add_fail_cnt <- add_cnt_tib %>% dplyr::filter(Address_Count!=1)
  
  man_prb_tib %>% dplyr::filter(U %in% add_fail_cnt$Address) %>% dplyr::group_by(Rep_Max) %>% dplyr::summarise(Rep_Count=n(), .groups='drop')
  man_prb_tib %>% dplyr::filter(M %in% add_fail_cnt$Address) %>% dplyr::group_by(Rep_Max) %>% dplyr::summarise(Rep_Count=n(), .groups='drop')
  
  man_prb_tib %>% dplyr::filter(U %in% add_fail_cnt$Address)
  man_prb_tib %>% dplyr::filter(M %in% add_fail_cnt$Address) %>% dplyr::filter(Rep_Max==1)
}



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               6.1 Intersect improbe Genomic Coordinates::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# - This needs to be moved up in order. Basically after manifest generation...
# - Also, needs to be saved with Genomic_CGN_Count added...
#

#
# Turn into function()
#
imp_pos_col <- cols(Seq_ID = col_character(),
                    Gen_Chr = col_character(),
                    Gen_Pos = col_double() )
imp_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic/GRCm10.improbeDesignInput.cgn-sorted.tsv.gz')

if (FALSE) {
  opt$cntDir  <- file.path(opt$outDir,'genomic')
  if (!dir.exists(opt$cntDir)) dir.create(opt$cntDir, recursive=TRUE)
  gen_cnt_csv <-
    file.path(opt$cntDir, paste(opt$runName,'genomic-counts.csv.gz', sep='.') )
  man_pos_tib <- 
    readr::read_tsv(imp_pos_tsv, col_names=names(imp_pos_col$cols), col_types=imp_pos_col) %>%
    dplyr::add_count(Seq_ID, name='Genomic_CGN_Count')
  readr::write_csv(man_pos_tib,gen_cnt_csv)
}

man_pos_tib <- 
  readr::read_tsv(imp_pos_tsv, col_names=names(imp_pos_col$cols), col_types=imp_pos_col) %>%
  dplyr::filter(Seq_ID %in% man_prb_tib$Seq_ID) %>%
  dplyr::add_count(Seq_ID, name='Genomic_CGN_Count')

#
# NOTES::
#  man_pos_tib = Genomic Position Table
#  
man_gen_tib <- man_prb_tib %>% dplyr::inner_join(
  dplyr::distinct(man_pos_tib, Seq_ID, .keep_all=TRUE), by="Seq_ID")

#
# Validation Summary of expected Genomic_CGN_Counts for matching designs!!!
#
man_gen_tib %>% 
  dplyr::group_by(Probe_Type,Genomic_CGN_Count,Infinium_Design) %>% 
  dplyr::summarise(Count=n(), .groups='drop') %>% as.data.frame()

#
# Need to understand Anti-join as well
#  - Should be zero
#
man_nog_tib <- man_prb_tib %>% dplyr::anti_join(
  dplyr::distinct(man_pos_tib, Seq_ID, .keep_all=TRUE), by="Seq_ID")

man_nog_tib %>%
  dplyr::group_by(Probe_Type,Infinium_Design,AQP) %>% 
  dplyr::summarise(Count=n(), .groups='drop')



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               6.2 Intersect improbe Probe Thermo Stats::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  mat_key <- 'Seq_ID'
  man_tsv <- file.path(opt$intDir, paste(opt$runName,'cgn_to_impScr','man',paste0(mat_key,'-sorted'),'tsv', sep='.'))
  int_tsv <- file.path(opt$intDir, paste(opt$runName,'cgn_to_impScr','int',paste0(mat_key,'-sorted'),'tsv.gz', sep='.'))
  
  imp_dat_tsv <- file.path(opt$impDir, 'designOutput_21092020/stats/GRCm10-21092020_improbe-designOutput.cgn-sorted.tsv.gz')
  imp_dat_tsv <- file.path(opt$impDir, 'designOutput_21092020/stats/GRCm10-21092020_improbe-designOutput.cgn-sorted.tsv')
  
  imp_col_tsv <- file.path(opt$impDir, 'designOutput_21092020/stats/GRCm10-21092020_improbe-designOutput.header.tsv')
  imp_col_vec <- suppressMessages(suppressWarnings( readr::read_tsv(imp_col_tsv) )) %>% 
    names() %>% as.vector()
  
  test_max <- 0
  man_scr_tib <- intersect_tsv(
    man=man_prb_tib, 
    man_tsv=man_tsv, imp_tsv=imp_dat_tsv, int_tsv=int_tsv, 
    colA=1,colB=1, mat_vec=c(mat_key), int_col=imp_col_vec, 
    fresh=opt$fresh, max=test_max,
    verbose=opt$verbose+2,tc=1,tt=pTracker)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   6.OLD Intersect Full Improbe Designs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Make function = man_join_improbe()
#

if (FALSE) {
  
  full_des_tsv <- file.path(opt$impDir, 'designOutput_21092020/GRCm10-21092020_improbe-designOutput.tsv.gz')
  imp_full_des_tsv <- file.path(opt$impDir, 'designOutput_21092020/GRCm10-21092020_improbe-designOutput.tsv.gz')
  man_uniq_ids_tsv <- file.path(opt$intDir, paste(opt$runName,'imp-cgn-srd-uniq-mat.cgn-sorted.tsv.gz', sep='.') )
  man_uniq_ids_tib <- dplyr::distinct(man_prb_tib,Seq_ID,SR_Str,CO_Str) %>% 
    dplyr::mutate(Probe_ID=paste(Seq_ID,paste0(SR_Str,CO_Str), sep='_')) %>%
    dplyr::select(Probe_ID) %>% dplyr::arrange(Probe_ID)
  readr::write_tsv(man_uniq_ids_tib,man_uniq_ids_tsv)
  
  # Intersect manifest unmethylated 48-mer matching CGN's with improbe design database::
  #
  run_sub_cmd <- TRUE
  run_sub_cmd <- FALSE
  imp_int_mat_cgn_tsv <- file.path(opt$intDir, paste(opt$runName,'imp-all.intersect.mat-cgn-srd.cgn-srd-sorted.tsv.gz', sep='.') )
  if (!file.exists(imp_int_mat_cgn_tsv) || run_sub_cmd) {
    # sub_exe <- '/Users/bretbarnes/Documents/scripts/subset/getSubset.improbe_TB-CO.pl'
    sub_exe <- file.path(par$srcDir, 'scripts/perl/improbe/getSubset.improbe_TB-CO.pl')
    sub_cmd <- glue::glue("{sub_exe} -header -t {man_uniq_ids_tsv} -d {full_des_tsv} | gzip -c -> {imp_int_mat_cgn_tsv}")
    cat(glue::glue("[{par$prgmTag}]: Running: cmd='{sub_cmd}'...{RET}{RET}") )
    system(sub_cmd)
  }
  
  # Load Matched improbe designs::
  #
  imp_int_mat_cgn_tib <- loadIMP(imp_int_mat_cgn_tsv, verbose=opt$verbose,tc=1,tt=pTracker) %>% 
    dplyr::mutate(SR_Str=stringr::str_sub(Methyl_Allele_TB_Strand,1,1),
                  CO_Str=Methyl_Allele_CO_Strand) %>%
    dplyr::rename(Probe_Type_IMP=Probe_Type)
  
  imp_int_mat_cgn_tib %>% dplyr::group_by(SR_Str,CO_Str) %>% 
    dplyr::summarise(Count=n(), .groups='drop') %>% print()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            7. Format Output::
  #
  #  - Determine Canonical Sequence
  #  - Split Manifest & Alignment
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Determine Matched Probes::
  mat_full_imp_tib <- mat_full_tib %>% 
    dplyr::inner_join(imp_int_mat_cgn_tib, by=c("Seq_ID","SR_Str","CO_Str")) %>% 
    dplyr::arrange(-Tango_CGN_Count) %>% 
    dplyr::distinct(U,M,SR_Str,CO_Str,Infinium_Design,Probe_Type, .keep_all=TRUE)
  mat_full_imp_tib %>% dplyr::group_by(Probe_Type,AQP) %>% 
    dplyr::summarise(Count=n(), .groups='drop')
  
  # Determine Missed Probes::
  mis_full_imp_tib <- mat_full_tib %>% dplyr::anti_join(imp_int_mat_cgn_tib, by=c("Seq_ID","SR_Str","CO_Str"))
  mis_full_imp_tib %>% dplyr::group_by(Probe_Type,AQP) %>% 
    dplyr::summarise(Count=n(), .groups='drop')
  
}



#
#
# TBD:: Need to join these two:: fin_full_tib & mat_full_imp_tib
#   - Chromosome Positions...
#

# fin_full_tib
# mat_full_imp_tib
#
# fin_full_tib %>% dplyr::inner_join(mat_full_imp_tib, by=c("M","U"), suffix=c("_Raw","_Des"))






# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          2.0. Remainder:: CpH/SNP/CTL::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Instead Check ALL individually::
non_imp_man_tib <- man_raw_tib

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Rebuilding Remainder with ALL:: SNP/CpH::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           2.1. Remainder:: SNP::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

snp_prb_des_csv <- file.path(opt$impDir, 'cph-snp-designs/LEGX_SpikeIn_Reorder-SNP-Only.designs.csv.gz')
snp_prb_des_tib <- readr::read_csv(snp_prb_des_csv) %>%
  dplyr::mutate(PRB1_U_MAT=stringr::str_to_upper(PRB1_U),
                PRB1_M_MAT=stringr::str_to_upper(PRB1_M),
                PRB2_D_MAT=stringr::str_to_upper(PRB2_D) )

# Match by ID::
#
snp_top_des_tib <- dplyr::inner_join(
  non_imp_man_tib, snp_prb_des_tib, by=c("Seq_ID"),
  suffix=c("_RAW","_DES")) %>% 
  dplyr::select(Seq_ID, IUPAC_Forward_Sequence,Probe_Type_RAW) %>%
  dplyr::rename(Top_Sequence=IUPAC_Forward_Sequence,Probe_Type=Probe_Type_RAW) %>% 
  dplyr::distinct()

test_max <- 0
snp_full_prb_des_tib <- desSeq_to_prbs(
  tib=snp_top_des_tib, 
  idsKey=opt$design_key,
  seqKey=opt$design_seq,prbKey=opt$design_prb,
  strsSR=opt$design_srs,strsCO=opt$design_cos,
  parallel=opt$parallel, max=test_max,
  verbose=opt$verbose,tc=1,tt=pTracker)

# snp_full_prb_des_tib <- tib2prbs(tib=snp_top_des_tib, 
#                                  idsKey=opt$design_key,
#                                  seqKey=opt$design_seq,
#                                  prbKey=opt$design_prb,
#                                  srdStr=opt$design_srd, 
#                                  parallel=opt$parallel,
#                                  verbose=opt$verbose,tc=1,tt=pTracker) %>%
#   dplyr::mutate(PRB1_U_MAT=stringr::str_to_upper(PRB1_U),
#                 PRB1_M_MAT=stringr::str_to_upper(PRB1_M),
#                 PRB2_D_MAT=stringr::str_to_upper(PRB2_D) )

#
# Match by Infinium I(1):: This produces nothing... No Infinum I CH 
#
snp_mat_inf1_tib <- 
  dplyr::inner_join(non_imp_man_tib, snp_full_prb_des_tib,
                    by=c("AlleleA_Probe_Sequence"="PRB1_U_MAT",
                         "AlleleB_Probe_Sequence"="PRB1_M_MAT"),
                    suffix=c("_RAW","_DES"))

snp_mat_inf1_tib %>% dplyr::group_by(Probe_Type_RAW,Infinium_Design) %>%
  dplyr::summarise(Count=n(), .groups='drop') %>% print()

#
# Match by Infinium II(2)::
#
snp_mat_inf2_tib <- 
  dplyr::inner_join(non_imp_man_tib,snp_full_prb_des_tib, 
                    by=c("AlleleA_Probe_Sequence"="PRB2_D_MAT"), 
                    suffix=c("_RAW","_DES"))

snp_mat_inf2_tib %>% dplyr::group_by(Probe_Type_RAW,Infinium_Design) %>% 
  dplyr::summarise(Count=n(), .groups='drop') %>% print()


snp_mat_tib <- 
  dplyr::bind_rows(snp_mat_inf1_tib,snp_mat_inf2_tib) %>%
  dplyr::rename(Seq_ID=Seq_ID_RAW,
                Probe_Type=Probe_Type_RAW) %>%
  dplyr::mutate(Seq_ID=stringr::str_replace_all(Seq_ID, ':','-') %>%
                  stringr::str_replace('^rs-', 'chr') %>% 
                  stringr::str_replace_all(' ',''))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           2.2. Remainder:: CpH::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# NOTE:: May have to use rds files instead of csv for some weird reason...
#
cph_prb_des_csv <- file.path(opt$impDir, 'cph-snp-designs/LEGX_SpikeIn_Reorder-CpH-Only.designs.csv.gz')
cph_prb_des_tib <- readr::read_csv(cph_prb_des_csv) %>%
  dplyr::mutate(PRB1_U_MAT=stringr::str_to_upper(PRB1_U),
                PRB1_M_MAT=stringr::str_to_upper(PRB1_M),
                PRB2_D_MAT=stringr::str_to_upper(PRB2_D) )

# Match by ID::
#
cph_top_des_tib <- dplyr::inner_join(
  non_imp_man_tib, cph_prb_des_tib, by=c("Seq_ID"),
  suffix=c("_RAW","_DES")) %>% 
  dplyr::select(Seq_ID, IUPAC_Forward_Sequence,Probe_Type_RAW) %>%
  dplyr::rename(Top_Sequence=IUPAC_Forward_Sequence,Probe_Type=Probe_Type_RAW) %>% 
  dplyr::distinct()

cph_full_prb_des_tib <-desSeq_to_prbs(
  tib=cph_top_des_tib, 
  idsKey=opt$design_key,
  seqKey=opt$design_seq,prbKey=opt$design_prb,
  strsSR=opt$design_srs,strsCO=opt$design_cos,
  parallel=opt$parallel, max=test_max,
  verbose=opt$verbose,tc=1,tt=pTracker)

# cph_full_prb_des_tib <- tib2prbs(tib=cph_top_des_tib, 
#                                  idsKey=opt$design_key,
#                                  seqKey=opt$design_seq,
#                                  prbKey=opt$design_prb,
#                                  srdStr=opt$design_srd, 
#                                  parallel=opt$parallel,
#                                  verbose=opt$verbose,tc=1,tt=pTracker) %>%
#   dplyr::mutate(PRB1_U_MAT=stringr::str_to_upper(PRB1_U),
#                 PRB1_M_MAT=stringr::str_to_upper(PRB1_M),
#                 PRB2_D_MAT=stringr::str_to_upper(PRB2_D) )

#
# Match by Infinium I(1):: This produces nothing... No Infinum I CH 
#
cph_mat_inf1_tib <- NULL
# if (FALSE) {
#   cph_mat_inf1_tib <- dplyr::inner_join(non_imp_man_tib, cph_full_prb_des_tib,
#                                         by=c("AlleleA_Probe_Sequence"="PRB1_U_MAT",
#                                              "AlleleB_Probe_Sequence"="PRB1_M_MAT"),
#                                         suffix=c("_RAW","_DES"))
#   cph_mat_inf1_tib %>% dplyr::group_by(Probe_Type_RAW,Infinium_Design) %>%
#     dplyr::summarise(Count=n(), .groups='drop') %>% print()
# }

#
# Match by Infinium II(2)::
#
cph_mat_inf2_tib <- dplyr::inner_join(
  non_imp_man_tib,cph_full_prb_des_tib, 
  by=c("AlleleA_Probe_Sequence"="PRB2_D_MAT"), 
  suffix=c("_RAW","_DES"))

cph_mat_inf2_tib %>% 
  dplyr::group_by(Probe_Type_RAW,Infinium_Design) %>% 
  dplyr::summarise(Count=n(), .groups='drop') %>% print()

cph_mat_tib <- 
  dplyr::bind_rows(cph_mat_inf1_tib,cph_mat_inf2_tib) %>%
  dplyr::rename(Seq_ID=Seq_ID_RAW,
                Probe_Type=Probe_Type_RAW) %>% 
  dplyr::mutate(Seq_ID=paste0('ch',stringr::str_remove(Seq_ID, '^ch') %>% 
                                stringr::str_pad(width=9, side='left', pad='0')) )

# cph_mat_tib %>% dplyr::select(Seq_ID, SR_Str,CO_Str,Infinium_Design, Rep_Num, Probe_Type, U,M,
#                               AlleleB_Probe_Sequence,AlleleB_Probe_Sequence, NXB_D, Top_Sequence,
#                               dplyr::everything()) %>% dplyr::arrange(Seq_ID)
#
# CpH Probes can be unique by Top Sequence
#
# if (FALSE) {
#   join_snp_cph_tib %>% dplyr::filter(Probe_Type_RAW=='ch') %>% 
#     dplyr::distinct(SR_Str, CO_Str,Top_Sequence, AlleleA_Probe_Sequence, AlleleB_Probe_Sequence) %>%
#     dplyr::pull(Top_Sequence)
#   
#   join_snp_cph_tib %>% dplyr::filter(Probe_Type_RAW=='rs') %>% 
#     dplyr::distinct(SR_Str, CO_Str,Top_Sequence, AlleleA_Probe_Sequence, AlleleB_Probe_Sequence) %>% 
#     dplyr::pull(Top_Sequence)
# }












# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   2.3. Collect Remainder:: CTL::MUS
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

raw_ctl_tib <- man_raw_tib %>% 
  dplyr::filter(Probe_Type != 'cg') %>% 
  dplyr::filter(Probe_Type != 'ch') %>% 
  dplyr::filter(Probe_Type != 'rs') %>% 
  dplyr::filter(Probe_Type != 'mu') %>% 
  dplyr::filter(Probe_Type != 'rp') %>% 
  dplyr::filter(Probe_Type == 'BS' | Probe_Type == 'NO' | Probe_Type == 'ne') %>% 
  dplyr::distinct(M,U, .keep_all=TRUE,Top_Sequence=NA) %>% 
  dplyr::mutate(SR_Str=FR,CO_Str=CO,Rep_Num=NA,Probe_Source='MUS',NXB_D=NA) %>% 
  dplyr::arrange(Probe_Type,Infinium_Design) %>%
  dplyr::mutate(
    Control_Group=dplyr::case_when(
      Probe_Type=='BS' & Infinium_Design==1 ~ "BISULFITE CONVERSION I",
      Probe_Type=='BS' & Infinium_Design==2 ~ "BISULFITE CONVERSION II",
      Probe_Type=='NO' ~ 'NON-POLYMORPHIC',
      Probe_Type=='ne' ~ 'NEGATIVE',
      TRUE ~ NA_character_
    ),
    Control_Group_Str=stringr::str_replace_all(Control_Group,' ','_') %>% 
      stringr::str_replace_all('-','_'),
    DiNuc=dplyr::case_when(
      Probe_Type=='BS' ~ stringr::str_replace(Seq_ID, '^.*-([ACTG][ACTG])-.*$', '\\$1') %>% 
        stringr::str_remove_all('\\\\'),
      Probe_Type=='NO' ~ stringr::str_replace(Seq_ID, '^.*_([ACTG][ACTG])_.*$', '\\$1') %>% 
        stringr::str_remove_all('\\\\'),
      TRUE ~ NA_character_
    ),
    N1=stringr::str_sub(DiNuc, 1,1),
    N2=stringr::str_sub(DiNuc, 2,2),
    Last_BaseA=stringr::str_sub(AlleleA_Probe_Sequence,-1),
    Last_BaseB=stringr::str_sub(AlleleB_Probe_Sequence,-1)
  ) %>% 
  dplyr::group_by(Probe_Type,Infinium_Design) %>%
  dplyr::mutate(
    Row_Idx=dplyr::row_number() + 100,
    Row_Str=Row_Idx,
    Col_Idx=Row_Idx %% (color_len-100)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Control_Name=dplyr::case_when(
      # Probe_Type=='BS' & Infinium_Design==1 & Last_BaseA=='G' ~ paste0(Control_Group_Str,'_U',Row_Str),
      # Probe_Type=='BS' & Infinium_Design==1 & Last_BaseA=='A' ~ paste0(Control_Group_Str,'_C',Row_Str),
      Probe_Type=='BS' & Infinium_Design==1 ~ paste0(Control_Group_Str,'_',Row_Str),
      Probe_Type=='BS' & Infinium_Design==2 ~ paste0(Control_Group_Str,'_',Row_Str),
      Probe_Type=='NO' & Infinium_Design==2 ~ paste0(Control_Group_Str,'_',Row_Str),
      Probe_Type=='ne' & Infinium_Design==2 ~ paste0(Control_Group_Str,'_',Row_Str),
      TRUE ~ NA_character_ ),
    Seq_ID_Org=Seq_ID,
    Seq_ID=Control_Name
  ) %>%
  dplyr::select(Seq_ID,SR_Str,CO_Str,Infinium_Design,Rep_Num,Probe_Type,U,M, 
                dplyr::everything())

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              2.4. Collect Remainder:: CTL::Complete HSA
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

std_ctl_tib <- NULL
ctl_csv <- ctls_vec[1]
std_ctl_seq_tsv <- file.path(par$datDir,'manifest/controls/01152015_DarkMatterControls.probe.match.tsv.gz')

if (!is.null(ctl_csv) && file.exists(ctl_csv) &&
    !is.null(std_ctl_seq_tsv) && file.exists(std_ctl_seq_tsv)) {
  
  # TBD:: This should be transferred to 
  std_ctl_seq_tib <- dplyr::inner_join(
    suppressMessages(suppressWarnings(readr::read_tsv(std_ctl_seq_tsv) )) %>% 
      dplyr::mutate(Address=stringr::str_remove(address_name, '^1') %>% as.integer()),
    suppressMessages(suppressWarnings(
      readr::read_csv(ctl_csv, col_names=c("Address","Probe_Type","COLOR_CHANNEL","Probe_ID")) )),
    by="Address") %>%
    dplyr::select(Address,Probe_Type,COLOR_CHANNEL,Probe_ID,probe_id,sequence, everything()) %>%
    dplyr::rename(Design_ID=probe_id) %>% dplyr::distinct(Address, .keep_all=TRUE) %>%
    dplyr::select(-type_b,-bo_seq,-address_name)
  
  std_ctl_tib <- std_ctl_seq_tib %>% 
    dplyr::mutate(PIDX=Probe_ID %>%
                    stringr::str_replace('^.*[^0-9]([0-9]+)$', '\\$1') %>% 
                    stringr::str_remove_all('\\\\')) %>%
    dplyr::mutate(Probe_ID=Probe_ID %>%
                    stringr::str_replace_all(' ', '_') %>%
                    stringr::str_replace_all('-', '_') %>% 
                    stringr::str_replace_all('\\(', '') %>% 
                    stringr::str_replace_all('\\)', '')) %>%
    dplyr::rename(U=Address) %>%
    dplyr::mutate(M=NA,DESIGN='II',col=NA,Probe_Source='HSA',Next_Base=NA) %>%
    dplyr::mutate(M=as.double(M), Probe_ID=paste('ctl',Probe_ID, sep='_')) %>%
    dplyr::select(Probe_ID,M,U,DESIGN,COLOR_CHANNEL,col,Probe_Type,Probe_Source,Next_Base,
                  dplyr::everything()) %>%
    dplyr::arrange(Probe_Type,Probe_ID) %>%
    dplyr::distinct(M,U, .keep_all=TRUE) %>% 
    dplyr::mutate(
      Design_Base_ID=stringr::str_remove(Design_ID, '_[AB]$'),
      Design_Base_AB=stringr::str_replace(Design_ID, '^.*_([AB])$','\\$1') %>%
        stringr::str_remove_all('\\\\'),
      Control_Group=Probe_Type,
      Control_Group_Str=stringr::str_replace_all(Control_Group,' ','_') %>% 
        stringr::str_replace_all('-','_'),
      Probe_Type=dplyr::case_when(
        Control_Group=="BISULFITE CONVERSION I"  ~ 'BS',
        Control_Group=="BISULFITE CONVERSION II" ~ 'BS',
        Control_Group=="EXTENSION"       ~ 'EX',
        Control_Group=="HYBRIDIZATION"   ~ 'HB',
        Control_Group=="NEGATIVE"        ~ 'ne',
        Control_Group=="NON-POLYMORPHIC" ~ 'NO',
        
        Control_Group=="NORM_A"          ~ 'NA',
        Control_Group=="NORM_C"          ~ 'NC',
        Control_Group=="NORM_G"          ~ 'NG',
        Control_Group=="NORM_T"          ~ 'NT',
        
        Control_Group=="RESTORATION"     ~ 'RE',
        Control_Group=="SPECIFICITY I"   ~ 'SP',
        Control_Group=="SPECIFICITY II"  ~ 'SP',
        
        Control_Group=="TARGET REMOVAL"  ~ 'TR',
        TRUE ~ NA_character_
      )
    )
  
  std_bsU_tib <- std_ctl_tib %>% dplyr::filter(stringr::str_detect(Probe_ID,'ctl_BS_Conversion_I_U')) %>%
    dplyr::rename(Address=U,AlleleA_Probe_Sequence=sequence) %>% 
    dplyr::select(-Probe_ID,-Design_ID,-M,-COLOR_CHANNEL,-col,-DESIGN,-Next_Base)
  
  std_bsC_tib <- std_ctl_tib %>% dplyr::filter(stringr::str_detect(Probe_ID,'ctl_BS_Conversion_I_C')) %>%
    dplyr::rename(Address=U,AlleleB_Probe_Sequence=sequence) %>% 
    dplyr::select(-Probe_ID,-Design_ID,-M,-COLOR_CHANNEL,-col,-DESIGN,-Next_Base)
  
  # Unique Addresses to Remove from Standard Controls
  std_bs1_add_vec <- dplyr::bind_rows(std_bsU_tib,std_bsC_tib) %>% 
    dplyr::distinct(Address) %>% dplyr::pull(Address)
  
  # Bisulfite Conversion I Probes for Standard Human Controls
  std_bs1_hum_tib <- 
    dplyr::inner_join(std_bsU_tib,std_bsC_tib, 
                      by=c("Design_Base_ID","Probe_Type","Probe_Source",
                           "Control_Group","Control_Group_Str","PIDX"), 
                      suffix=c("_U", "_C") ) %>%
    dplyr::rename(U=Address_U,M=Address_C) %>% 
    dplyr::mutate(Seq_ID=paste0(Control_Group_Str,'_',PIDX),
                  Infinium_Design=1,Rep_Num=1) %>%
    dplyr::select(Seq_ID,Infinium_Design,Rep_Num,Probe_Type,U,M,
                  AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                  dplyr::everything())
  
  #
  # Build remainder of Standard EPIC Controls::
  #
  std_inf1_hum_tib <- std_ctl_tib %>% 
    dplyr::filter(!U %in% std_bs1_add_vec) %>% 
    dplyr::rename(AlleleA_Probe_Sequence=sequence) %>%
    dplyr::mutate(
      AlleleB_Probe_Sequence=NA,
      DiNuc=NA_character_,
      N1=NA_character_,
      N2=NA_character_
    ) %>%
    dplyr::mutate(
      Infinium_Design=1,
      Seq_ID=paste0(Control_Group_Str,'_',PIDX) %>%
        stringr::str_replace_all(' ', '_') %>%
        stringr::str_replace_all('-', '_') %>% 
        stringr::str_replace_all('\\(', '') %>% 
        stringr::str_replace_all('\\)', ''))
}

# These are not compatiable in format yet::
all_ctl_tib <- dplyr::bind_rows(std_ctl_tib,raw_ctl_tib) %>%
  dplyr::distinct(M,U, .keep_all=TRUE)

# Standard EPIC Probes need to be joined for Infinium I probes

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              2.5. Collect Remainder:: CTL::Selected HSA/MUS
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Load Manually Selected Human and Mouse Controls for Genome Studio::
#
sel_ctl_col <- c('Address','Probe_Type','Control_Color','Probe_ID')
sel_ctl_csv <- file.path(par$topDir,'data/manifests/tmp/mm10_LEGX_B13.manifest.GenomeStudio.cpg-sorted.clean.controls.csv.csv')
sel_ctl_tib <- suppressMessages(suppressWarnings( readr::read_csv(sel_ctl_csv, col_names=sel_ctl_col) )) %>%
  dplyr::distinct(Address, .keep_all=TRUE)

ses_sel_ctl_tib <- dplyr::bind_rows(
  raw_ctl_tib %>% dplyr::filter(U %in% sel_ctl_tib$Address) %>% 
    dplyr::select(Seq_ID,Infinium_Design,U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,SR_Str,CO_Str,
                  dplyr::everything()),
  
  std_bs1_hum_tib %>% dplyr::filter(U %in% sel_ctl_tib$Address) %>% 
    dplyr::mutate(SR_Str=NA,CO_Str=NA) %>%
    dplyr::select(Seq_ID,Infinium_Design,U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,SR_Str,CO_Str,
                  dplyr::everything()),
  
  std_inf1_hum_tib %>% dplyr::filter(U %in% sel_ctl_tib$Address) %>% 
    dplyr::mutate(SR_Str=NA,CO_Str=NA) %>%
    dplyr::select(Seq_ID,Infinium_Design,U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,SR_Str,CO_Str,
                  dplyr::everything())
) %>% 
  dplyr::mutate(Seq_ID=stringr::str_replace_all(Seq_ID,'_','-'),
                # Seq_ID=paste('ctl',Seq_ID, sep='-'),
                Seq_ID=paste(Seq_ID, sep='-'),
                Seq_ID=paste(Seq_ID,Probe_Source, sep='_'))

ses_sel_ctl_cnt <- ses_sel_ctl_tib %>% base::nrow()
ses_unq_ctl_cnt <- ses_sel_ctl_tib %>% dplyr::distinct(U,M) %>% base::nrow()
ses_sel_ctl_sum <- ses_sel_ctl_tib %>% dplyr::group_by(Probe_Source,Probe_Type,Infinium_Design) %>% 
  dplyr::summarise(Count=n(), .groups='drop')
cat(glue::glue("[{par$prgmTag}]: ses_sel_ctl_cnt={ses_sel_ctl_cnt}, ses_unq_ctl_cnt={ses_unq_ctl_cnt}; ses_sel_ctl_sum={RET}"))
ses_sel_ctl_sum %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  3.0. Collect Remainder:: CTL/SNP/CpH/CpG::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fin_ctl_mat_tib <- 
  ses_sel_ctl_tib %>% dplyr::select(Seq_ID, SR_Str,CO_Str,Infinium_Design, Rep_Num, Probe_Type, U,M,
                                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, NXB_D, Top_Sequence,
                                    dplyr::everything()
                                    ) %>% dplyr::arrange(Seq_ID) %>%
  dplyr::mutate(Assay_Class='Control', Seq_ID=paste('ctl',Seq_ID, sep='-'))

fin_snp_mat_tib <-
  snp_mat_tib %>% dplyr::select(Seq_ID, SR_Str,CO_Str,Infinium_Design, Rep_Num, Probe_Type, U,M,
                                AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, NXB_D, Top_Sequence,
                                dplyr::everything()
                                ) %>% dplyr::arrange(Seq_ID) %>%
  dplyr::mutate(Assay_Class='Analytical', Probe_Source='MUS')

fin_cph_mat_tib <-
  cph_mat_tib %>% dplyr::select(Seq_ID, SR_Str,CO_Str,Infinium_Design, Rep_Num, Probe_Type, U,M,
                                AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, NXB_D, Top_Sequence,
                                dplyr::everything()
                                ) %>% dplyr::arrange(Seq_ID) %>%
  dplyr::mutate(Assay_Class='Analytical', Probe_Source='MUS')

fin_full_tib <- man_gen_tib
fin_cpg_mat_tib <- 
  fin_full_tib %>% dplyr::select(Seq_ID, SR_Str,CO_Str,Infinium_Design, Rep_Num, Probe_Type, U,M,
                                 AlleleA_Probe_Sequence, AlleleB_Probe_Sequence, NXB_D, Top_Sequence,
                                 dplyr::everything()
                                 ) %>% dplyr::arrange(Seq_ID) %>%
  dplyr::mutate(Assay_Class='Analytical', Probe_Source='MUS')


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                  3.1. Clean Up Remainder:: CTL/SNP/CpH/CpG::
#
#  - Make unique Tango Addresses (M,U)
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Uniaue Tango Address Pairs::
#
fin_core_all_tib <- NULL

opt$skip_controls <- FALSE
if (opt$skip_controls) {
  prb_core_all_tib <-
    dplyr::bind_rows(fin_snp_mat_tib,fin_cph_mat_tib,fin_cpg_mat_tib)
  prb_core_unq_tib <- prb_core_all_tib %>% 
    dplyr::distinct(M,U, .keep_all=TRUE) %>%
    dplyr::add_count(M,U, name='Tango_Pair_Count')
  
  prb_core_unq_tib %>% dplyr::group_by(Probe_Type,Infinium_Design,AQP) %>%
    dplyr::summarise(Count=n(), .groups='drop') %>% as.data.frame()
  prb_core_unq_tib %>% dplyr::filter(Tango_Pair_Count!=1) %>% dplyr::arrange(U)
  man_raw_tib %>% 
    dplyr::anti_join(prb_core_unq_tib, by=c("M","U")) %>% 
    dplyr::anti_join(fin_ctl_mat_tib,  by=c("M","U")) %>%
    dplyr::group_by(Probe_Type,AQP) %>% dplyr::summarise(Mis_Count=n())
} else {
  
  # Doesn't Work::
  #
  # fin_core_all_tib <- 
  #   dplyr::bind_rows(fin_ctl_mat_tib,fin_snp_mat_tib,fin_cph_mat_tib,fin_cpg_mat_tib)
  
  # Woks::
  #  dplyr::bind_rows(fin_cpg_mat_tib %>% dplyr::select(1:34,37:82), dplyr::bind_rows(fin_ctl_mat_tib,fin_snp_mat_tib,fin_cph_mat_tib) )
  
  fin_core_all_tib <-
    dplyr::bind_rows(fin_ctl_mat_tib,fin_snp_mat_tib,fin_cph_mat_tib, dplyr::select(fin_cpg_mat_tib,1:34,37:82))
  
  fin_core_unq_tib <- fin_core_all_tib %>% dplyr::distinct(U,M, .keep_all=TRUE) %>%
    dplyr::add_count(M,U, name='Tango_Pair_Count')
  fin_core_unq_tib %>% dplyr::group_by(Probe_Type,Infinium_Design,AQP) %>% 
    dplyr::summarise(Count=n(), .groups='drop') %>% as.data.frame()
  fin_core_mis_cnt <- fin_core_unq_tib %>% dplyr::filter(Tango_Pair_Count!=1) %>% base::nrow()
}























#
# Clean Up Fields to match Sesame::
#
# fin_core_ses_tib <- prb_core_unq_tib %>% 

fin_core_ses_tib <- fin_core_unq_tib %>% 
  dplyr::rename(Probe_ID2=Probe_ID) %>% 
  dplyr::mutate(
    Probe_ID=dplyr::case_when(
      Assay_Class=='Analytical' ~ paste0(Seq_ID,'_',SR_Str,CO_Str,Infinium_Design,Rep_Num),
      Assay_Class=='Control'    ~ Seq_ID,
      TRUE ~ NA_character_
    ),
    DESIGN=dplyr::case_when(
      Infinium_Design==1 ~ 'I',
      Infinium_Design==2 ~ 'II',
      TRUE ~ NA_character_),
    COLOR_CHANNEL=dplyr::case_when(
      Infinium_Design==2 ~ 'Both',
      stringr::str_to_upper(NXB_D)=='A' ~ 'Red',
      stringr::str_to_upper(NXB_D)=='T' ~ 'Red',
      stringr::str_to_upper(NXB_D)=='C' ~ 'Grn',
      stringr::str_to_upper(NXB_D)=='G' ~ 'Grn',
      TRUE ~ NA_character_),
    col=dplyr::case_when(
      Infinium_Design==2 ~ NA_character_,
      stringr::str_to_upper(NXB_D)=='A' ~ 'R',
      stringr::str_to_upper(NXB_D)=='T' ~ 'R',
      stringr::str_to_upper(NXB_D)=='C' ~ 'G',
      stringr::str_to_upper(NXB_D)=='G' ~ 'G',
      TRUE ~ NA_character_),
    Next_Base=dplyr::case_when(
      Infinium_Design==2 ~ NA_character_,
      Infinium_Design==1 ~ stringr::str_to_upper(NXB_D),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::rename(Strand_TB=SR_Str,Strand_CO=CO_Str) %>%
  dplyr::select(Probe_ID,M,U,DESIGN,COLOR_CHANNEL,col,Probe_Type, Probe_Source,Next_Base,
                Seq_ID,Strand_TB,Strand_CO,Infinium_Design,Rep_Num,
                AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Top_Sequence, everything())

#
# General Stats::
#
fin_core_ses_tib %>% dplyr::group_by(Probe_Type,DESIGN,AQP) %>% 
  dplyr::summarise(Count=n(), .groups='drop') %>% 
  print(n=base::nrow(fin_core_ses_tib))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  4.0. Add Filtering Field to all Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

css_csv <- file.path(par$topDir, 'tmp/LifeEpigentics/scratch/build_models/LEGX/S1/Sample_Name/mm10-ILS-VAI.Titration/ind-beta_i-poob-1/Sample_Name_mm10-ILS-VAI.Titration_ind-beta_i-poob-1.full-dbl.csv.gz')
css_tib <- readr::read_csv(css_csv) %>%
  dplyr::filter(!stringr::str_starts(Probe_ID,'ctl_')) %>%
  dplyr::mutate(Seq_ID=stringr::str_remove(Probe_ID,'_.*$'),
                Strand_TB=stringr::str_remove(Probe_ID,'^.*_') %>% stringr::str_replace('([A-Z]).*$','\\$1') %>% stringr::str_remove_all('\\\\'),
                Strand_CO=stringr::str_remove(Probe_ID,'^.*_') %>% stringr::str_replace('[A-Z]([A-Z]).*$','\\$1') %>% stringr::str_remove_all('\\\\'),
                Infinium_Design=stringr::str_remove(Probe_ID,'^.*_') %>% stringr::str_replace('[A-Z][A-Z]([12]).*$','\\$1') %>% 
                  stringr::str_remove_all('\\\\') %>% as.integer()
  ) %>% 
  dplyr::select(Seq_ID, Strand_TB,Strand_CO,Infinium_Design, everything()) %>%
  dplyr::distinct(Seq_ID, Strand_TB,Strand_CO,Infinium_Design, .keep_all=TRUE)

fin_core_css_tib <- 
  dplyr::left_join(fin_core_ses_tib,css_tib, 
                   by=c("Seq_ID","Strand_TB","Strand_CO","Infinium_Design"),
                   suffix=c("_MAN","_CSS")) %>%
  dplyr::rename(Probe_ID=Probe_ID_MAN)

#
# TBD:: Fail all non-cpgs::
#
fin_core_tag_tib <- fin_core_css_tib %>% 
  dplyr::mutate(
    CSS_Fail=dplyr::case_when(
      T00DZ_T99DZ_CSS_mu>=0.5 & T00DZ_T50DZ_CSS_mu>=0.2 & T50DZ_T99DZ_CSS_mu>=0.1 ~ 'FALSE',
      T00DZ_T99DZ_CSS_mu<0.5 | T00DZ_T50DZ_CSS_mu<0.2 | T50DZ_T99DZ_CSS_mu<0.1 ~ 'TRUE',
      TRUE ~ NA_character_
    ),
    MFG_Change_Flagged=dplyr::case_when(
      Probe_Type != 'cg' ~ 'TRUE',
      CSS_Fail == 'FALSE' ~ 'FALSE',
      CSS_Fail == 'TRUE' ~ 'TRUE',
      is.na(CSS_Fail) ~ 'TRUE',
      TRUE ~ NA_character_
    )
  )
fin_core_tag_tib %>% dplyr::group_by(CSS_Fail) %>% 
  dplyr::summarise(MFG_Count=n(), .groups='drop')
fin_core_tag_tib %>% dplyr::group_by(Probe_Type,MFG_Change_Flagged) %>% 
  dplyr::summarise(MFG_Count=n(), .groups='drop')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Write Genome Studio Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# gs_man_ver  <- 0
# gs_man_let  <- 'C'
# gs_swap_dir <- file.path(opt$outDir,  paste0('mm10-LEGX-',gs_man_let,gs_man_ver))
# gs_swap_csv <- file.path(gs_swap_dir, paste0('mm10_LEGX_',gs_man_let,gs_man_ver,'.manifest.GenomeStudio.cpg-sorted.clean.csv') )
# gz_swap_csv <- file.path(gs_swap_dir, paste0('mm10_LEGX_',gs_man_let,gs_man_ver,'.manifest.GenomeStudio.cpg-sorted.clean.csv.gz') )
# ses_man_csv <- file.path(gs_swap_dir, paste0('mm10_LEGX_',gs_man_let,gs_man_ver,'.manifest.Sesame.cpg-sorted.clean.csv.gz') )

opt$version
gs_swap_dir <- file.path(opt$prdDir,  paste0('mm10-LEGX-',opt$version))
gs_swap_csv <- file.path(gs_swap_dir, paste0('mm10_LEGX_',opt$version,'.manifest.GenomeStudio.cpg-sorted.clean.csv') )
gz_swap_csv <- file.path(gs_swap_dir, paste0('mm10_LEGX_',opt$version,'.manifest.GenomeStudio.cpg-sorted.clean.csv.gz') )
ses_man_csv <- file.path(gs_swap_dir, paste0('mm10_LEGX_',opt$version,'.manifest.Sesame.cpg-sorted.clean.csv.gz') )

if (!dir.exists(gs_swap_dir)) dir.create(gs_swap_dir, recursive=TRUE)

man_gs_csv  <- file.path(par$topDir, 'data/manifests/tmp/mm10_LEGX_B13.manifest.GenomeStudio.cpg-sorted.clean.csv.gz')
man_head_df <- readr::read_lines(man_gs_csv, n_max = 6) %>% as.data.frame()

# man_gs_list <- loadManifestGenomeStudio(file = man_gs_csv, addSource = TRUE, normalize = FALSE, verbose = 20)
ctl_head_df <- tibble::tibble(Head="[Controls]") %>% tibble::add_column(BL1='',BL2='',BL3='',BL4='') %>% as.data.frame()

# man_head_replace_str <- paste0('_B',gs_man_ver,'.csv.gz')
man_head_replace_str <- paste0('_',opt$version,'.csv.gz')
man_head_df <- man_head_df %>% as_tibble() %>% purrr::set_names(c("Head_Var")) %>% 
  dplyr::mutate(Head_Var=as.character(Head_Var),
                Head_Var=stringr::str_replace(Head_Var, '_B[0-9]+.csv.gz$', man_head_replace_str) ) %>%
  as.data.frame()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Write Sesame Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fin_ses_core_csv <- file.path(opt$prdDir)
fin_ses_core_tib <- fin_core_tag_tib %>% dplyr::arrange(Probe_ID) %>% 
  dplyr::distinct(M,U, .keep_all=TRUE) %>% 
  dplyr::select(Probe_ID:Top_Sequence,Assay_Class,MFG_Change_Flagged)
readr::write_csv(fin_ses_core_tib,ses_man_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Write Genome Studio Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fin_core_split <- fin_core_tag_tib %>% split(.$Assay_Class)

#
# Format Genome Studio Analytical 
#
man_gs_col <- c("IlmnID","Name","AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq",
                "Infinium_Design_Type","Next_Base","Color_Channel","Forward_Sequence","Top_Sequence",
                "Genome_Build","Genome_Build_NCBI","CHR","MAPINFO","SourceSeq",
                "Strand","Strand_TB","Strand_CO")

man_fin_tib <- fin_core_split[['Analytical']] %>% 
  dplyr::arrange(Probe_ID) %>%
  dplyr::distinct(M,U, .keep_all=TRUE) %>% 
  dplyr::mutate(
    SourceSeq=PRB1_D,
    Forward_Sequence=dplyr::case_when(
      FR=='F' & Strand_TB=="T" ~ Top_Sequence,
      FR=='F' & Strand_TB=="B" ~ revCmp(Top_Sequence),
      FR=='R' & Strand_TB=="B" ~ Top_Sequence,
      FR=='R' & Strand_TB=="T" ~ revCmp(Top_Sequence),
      TRUE ~ paste0(paste0(rep('N',60), collapse=''),'[NN]',paste0(rep('N',60), collapse=''))
    ),
    Genome_Build="mm10",
    Genome_Build_NCBI="GRCm10"
  ) %>%
  dplyr::rename(IlmnID=Probe_ID, Name=Seq_ID, 
                AddressA_ID=U,AlleleA_ProbeSeq=AlleleA_Probe_Sequence, 
                AddressB_ID=M,AlleleB_ProbeSeq=AlleleB_Probe_Sequence,
                Infinium_Design_Type=DESIGN,
                Color_Channel=COLOR_CHANNEL,
                CHR=Gen_Chr,
                MAPINFO=Gen_Pos,
                Strand=FR) %>%
  dplyr::select(all_of(man_gs_col))

#
# Format Genome Studio Controls
#
ctl_fin_tib <- fin_core_split[['Control']] %>% 
  dplyr::arrange(Probe_ID) %>%
  dplyr::distinct(M,U, .keep_all=TRUE) %>% 
  dplyr::select(M,U,Control_Group,Probe_ID:Top_Sequence,Control_Group)

#
# TBD:: For Infinium I designs add appropriate U/C designation to the Probe_ID
#

#
# Find this 24637490 human bs1 control is it U?
#
ctl_gs_tib <- dplyr::bind_rows(
  ctl_fin_tib %>%
    dplyr::select(U,Control_Group,Probe_ID) %>%
    dplyr::filter(!is.na(U)) %>% 
    dplyr::rename(Address=U) %>%
    dplyr::mutate(Probe_ID=dplyr::case_when(
      Control_Group=='BISULFITE CONVERSION I' ~ stringr::str_replace(Probe_ID, '_([MUSHA]+)$', 'U_\\$1') %>% stringr::str_remove_all('\\\\'),
      TRUE ~ Probe_ID )
    ),
  ctl_fin_tib %>%
    dplyr::select(M,Control_Group,Probe_ID) %>%
    dplyr::filter(!is.na(M)) %>% 
    dplyr::rename(Address=M) %>%
    dplyr::mutate(Probe_ID=dplyr::case_when(
      Control_Group=='BISULFITE CONVERSION I' ~ stringr::str_replace(Probe_ID, '_([MUSHA]+)$', 'M_\\$1') %>% stringr::str_remove_all('\\\\'),
      TRUE ~ Probe_ID )
    )
) %>%
  dplyr::arrange(Control_Group) %>% 
  dplyr::inner_join(sel_ctl_tib %>% dplyr::select(Address,Control_Color), by="Address") %>%
  dplyr::select(Address,Control_Group,Control_Color,Probe_ID)

# Sanity Check::
#  ctl_gs_tib %>% dplyr::filter(Control_Group=='BISULFITE CONVERSION I')

#
# Write Genome Studio CSV::
#

write_delim(x=man_head_df, path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=FALSE)
write_delim(x=man_fin_tib, path=gs_swap_csv, delim=',', col_names=TRUE,  quote_escape=FALSE, na='', append=TRUE)
write_delim(x=ctl_head_df, path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)
write_delim(x=ctl_gs_tib,  path=gs_swap_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)

cmd <- glue::glue("cat {gs_swap_csv} | perl -pe 's/\"//gi; s/\n/,,,,,\n/;' | gzip -c - > {gz_swap_csv}")
system(cmd)

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
