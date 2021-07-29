
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages( base::require("GenomicRanges",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("plyr")) )
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringi") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("data.table") ))
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
BNG <- "|"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Define Default Params and Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
opt <- NULL
par <- NULL
run <- NULL

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
par$prgmDir <- 'workhorse'
par$prgmTag <- 'workhorse_main_dev_latest'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))


# Executables::
opt$Rscript <- NULL

# BSMAP Parameters::
opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
opt$bsmap_opt <- "\"-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R\""
opt$bsmap_exe <- "/Users/bretbarnes/Documents/tools/programs/BSMAPz/bsmapz"

# Run Parameters::
opt$runName    <- NULL
opt$Species    <- NULL

# Null Place Holders::
opt$cpg_top_tsv <- NULL
opt$cpg_pos_tsv <- NULL
opt$cph_pos_tsv <- NULL
opt$snp_pos_tsv <- NULL
opt$ord_des_csv <- NULL
opt$canonical_csv <- NULL

# Directories::
opt$outDir  <- NULL
opt$impDir  <- NULL
opt$annDir  <- NULL
opt$gen_dir  <- NULL
opt$manDir  <- NULL

# Manufacturing Files:: Required
opt$ords <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib
opt$mats <- NULL
opt$aqps <- NULL

opt$mans <- NULL  # Manifest of probes to be added
opt$noob <- NULL  # Manifest of probes to be noob masked
opt$vcfs <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib
opt$beds <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib
opt$snps <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib

# Manufacturing Info:: Required
opt$bpns <- NULL
opt$aqpn <- NULL

# Pre-defined files (Controls & IDAT Validation):: Optional
opt$gs_ctl_csv  <- NULL
opt$ses_ctl_csv <- NULL

opt$idat <- NULL

# Platform/Method Options::
opt$genome_build <- NULL
opt$platform <- NULL
opt$version  <- NULL

# Process Parallel/Cluster Parameters::
opt$single   <- FALSE
opt$parallel <- TRUE
opt$cluster  <- FALSE

# Run-time Files
opt$opt_csv  <- NULL
opt$par_csv  <- NULL
opt$time_csv <- NULL
opt$time_org_txt <- NULL

# Run-time Options::
opt$fresh  <- FALSE
opt$reload <- FALSE

# verbose Options::
opt$verbose <- 3

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Order/Match/AQP/PQC Expected Columns::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: This is copied in the code, so it should be removed from here
#   or moved to a config file... Maybe a single r-data structure (RDS)
#   for manifest parameter defaults is the way to go...
#
par_cols <- list()
par_cols$ord <- 
  cols(
    Assay_Design_Id        = col_character(),
    AlleleA_Probe_Id       = col_character(),
    AlleleA_Probe_Sequence = col_character(),
    AlleleB_Probe_Id       = col_character(),
    AlleleB_Probe_Sequence = col_character(),
    Normalization_Bin      = col_character()
  )

par_cols$mat <- 
  cols(
    Plate    = col_character(),
    Row      = col_character(),
    Col      = col_integer(),
    Address  = col_integer(),
    Mod5     = col_character(),
    Sequence = col_character(),
    Mod3     = col_character(),
    Comments = col_character()
  )

par_cols$ma2 <- 
  cols(
    address_names = col_integer(),
    probe_id      = col_character(),
    sequence      = col_character(),
    type_b        = col_character(),
    address_name  = col_integer(),
    bo_seq        = col_character()
  )

par_cols$aqp <- 
  cols(
    Address           = col_integer(),
    Decode_Status     = col_integer(),
    Decode_Error_Code = col_integer(),
    Decode_Score      = col_integer(),
    Func_Status       = col_integer(),
    Func_Error_Code   = col_integer(),
    QC_Action         = col_integer()
  )

par_cols$pqc <- 
  cols(
    Address      = col_integer(),
    Status       = col_integer(),
    Eval_Code    = col_integer(),
    Average_Rep  = col_integer(),
    Expected_Rep = col_integer()
  )

par$ord_col <- par_cols$ord$cols %>% names()

par$mat_col <- par_cols$mat$cols %>% names()
par$ma2_col <- par_cols$ma2$cols %>% names()

par$aqp_col <- par_cols$aqp$cols %>% names()
par$pqc_col <- par_cols$pqc$cols %>% names()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Local Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_source_files = function(dir, verbose, funcTag="load_source_files") {
  gen_src_dir <- file.path(dir, 'functions')
  if (!dir.exists(gen_src_dir))
    stop(glue::glue("[{funcTag}]: General Source={gen_src_dir} ",
                    "does not exist!!!{RET}{RET}"))
  
  for (sfile in list.files(path=gen_src_dir, pattern='.R$', 
                           full.names=TRUE, recursive=TRUE)) base::source(sfile)
  if (verbose>0)
    cat(glue::glue("[{funcTag}]: Done. Loading Source Files form ",
                   "General Source={gen_src_dir}!{RET}{RET}") )
  
  gen_src_dir
}

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
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda2-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda3-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/conda_4.6.8/bin/Rscript'
  
  #
  # End of local parameter definitions::
  #
  
  opt$outDir  <- file.path(par$topDir, 'scratch')
  opt$impDir  <- file.path(par$topDir, 'data/improbe')
  opt$annDir  <- file.path(par$topDir, 'data/annotation')
  opt$manDir  <- file.path(par$topDir, 'data/manifests')
  opt$gen_dir  <- file.path(par$topDir, 'data/iGenomes/Homo_sapiens/NCBI')
  # opt$idatDir <- file.path(par$topDir, 'data/idats')
  
  opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
  opt$bsmap_opt <- "\"-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R\""
  opt$bsmap_exe <- "/Users/bretbarnes/Documents/tools/programs/BSMAPz/bsmapz"
  
  # Pre-defined local options runTypes::
  #
  par$local_runType <- NULL
  par$local_runType <- 'COVIC'
  par$local_runType <- 'GSA'
  par$local_runType <- 'HM450'
  par$local_runType <- 'TruDx'
  par$local_runType <- 'EWAS'
  par$local_runType <- 'GRCm10'
  par$local_runType <- 'NZT'
  par$local_runType <- 'Chicago'
  par$local_runType <- 'McMaster10Kselection'
  
  opt$parallel <- TRUE
  
  opt$verbose <- 10
  opt$verbose <- 3
  opt$verbose <- 5
  
  opt$fresh  <- TRUE
  opt$fresh  <- FALSE
  opt$reload <- TRUE
  
  opt$run_improbe    <- TRUE
  opt$build_manifest <- TRUE
  
  if (FALSE) {
    
  } else if (par$local_runType=='McMaster10Kselection') {
    opt$genome_build <- 'GRCh37'
    opt$platform <- 'MCM'
    opt$version  <- 'v2'
    opt$Species  <- "Human"
    
    par$aqpDir <- file.path(par$topDir, "data/CustomContent/McMaster/McMaster10Kselection/AQP.v2")
    opt$ords <- paste(
      file.path(par$aqpDir, 'McMaster_CpG_DesignFile_v4.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20532820_probes1.match.gz'),
      file.path(par$aqpDir, '20532820_probes2.match.gz'),
      sep=',')
    
    opt$aqps <- paste(
      # file.path(par$aqpDir, 'BS0033057-AQP1.txt.gz'),
      # file.path(par$aqpDir, 'BS0033090-AQP2.txt.gz'),
      file.path(par$aqpDir, '20051339_A_ProductQC.txt.gz'),
      sep=',')
    
    if (FALSE) {
      opt$aqps <- paste(
        file.path(par$aqpDir, 'BS0033057-AQP1.txt.gz'),
        file.path(par$aqpDir, 'BS0033090-AQP2.txt.gz'),
        # file.path(par$aqpDir, '20051339_A_ProductQC.txt.gz'),
        sep=',')
    }
    
    opt$noob <- paste(
      file.path(par$topDir, "data/CustomContent/transfer/updated_manifest.csv.gz"),
      sep=',')
    
    opt$bpns <- paste(1, sep=",")
    opt$aqpn <- paste(1, sep=",")
    
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$runName <- paste(par$local_runType,opt$platform,opt$version,opt$genome_build, sep='-')
  
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
    
    make_option(c("--bsmap_opt"), type="character", default=opt$bsmap_opt, 
                help="BSMAP Options [default= %default]", metavar="character"),
    make_option(c("--bsmap_exe"), type="character", default=opt$bsmap_exe, 
                help="BSMAP Executable path [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--Species"), type="character", default=opt$Species, 
                help="Target Species Name [default= %default]", metavar="character"),
    make_option(c("--pre_man_csv"), type="character", default=opt$pre_man_csv, 
                help="Previously defined manifest [default= %default]", metavar="character"),
    
    make_option(c("--cpg_pos_tsv"), type="character", default=opt$cpg_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--cpg_top_tsv"), type="character", default=opt$cpg_top_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--cph_pos_tsv"), type="character", default=opt$cph_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--snp_pos_tsv"), type="character", default=opt$snp_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--ord_des_csv"), type="character", default=opt$ord_des_csv, 
                help="Original design file used for canonical position selection. [default= %default]", metavar="character"),
    make_option(c("--canonical_csv"), type="character", default=opt$canonical_csv, 
                help="Canonical CGN names from previous products. [default= %default]", metavar="character"),
    
    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    
    make_option(c("--impDir"), type="character", default=opt$impDir, 
                help="improbe data directory [default= %default]", metavar="character"),
    make_option(c("--annDir"), type="character", default=opt$iannDir, 
                help="Annotation data directory [default= %default]", metavar="character"),
    make_option(c("--gen_dir"), type="character", default=opt$gen_dir, 
                help="Genomic data directory [default= %default]", metavar="character"),
    make_option(c("--manDir"), type="character", default=opt$manDir, 
                help="Pre-built Manifest data directory [default= %default]", metavar="character"),
    
    # Manufacturing Files:: Required
    make_option(c("--ords"), type="character", default=opt$ords, 
                help="Order file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mats"), type="character", default=opt$mats, 
                help="Match (format 1 or 2) file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqps"), type="character", default=opt$aqps, 
                help="AQP file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mans"), type="character", default=opt$mans,
                help="Manifest CSV file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--noob"), type="character", default=opt$noob, 
                help="Manifest of probes to be added and noob-masked CSV file(s) (comma seperated) [default= %default]", metavar="character"),
    
    # Not fully supported yet::
    make_option(c("--vcfs"), type="character", default=opt$vcfs, 
                help="Target Design VCF file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--beds"), type="character", default=opt$beds, 
                help="Target Design Coordinate BED file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--snps"), type="character", default=opt$snps, 
                help="Target Design Coordinate SNP file(s) (comma seperated) [default= %default]", metavar="character"),
    
    # Manufacturing Info:: Required
    make_option(c("--bpns"), type="character", default=opt$bpns, 
                help="Bead Pool Numbers (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqpn"), type="character", default=opt$aqpn, 
                help="AQP Numbers (comma seperated) [default= %default]", metavar="character"),
    
    # Pre-defined files (Controls & IDAT Validation):: Optional
    make_option(c("--gs_ctl_csv"), type="character", default=opt$gs_ctl_csv, 
                help="Pre-Defined Infinium Methylation Controls (Genome Studio) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--ses_ctl_csv"), type="character", default=opt$ses_ctl_csv, 
                help="Pre-Defined Infinium Methylation Controls (Sesame) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--idat"), type="character", default=opt$idat, 
                help="idat directories (comma seperated) [default= %default]", metavar="character"),
    
    # Platform/Method Options::
    make_option(c("--genome_build"), type="character", default=opt$genome_build, 
                help="Genome Build (e.g. GRCh36, GRCh37, GRCh38, GRCm38) [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform (e.g. HM450, EPIC, LEGX, NZT, COVIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest Version (e.g. B0,B1,B2,B3,B4,C0) [default= %default]", metavar="character"),
    
    # Process Parallel/Cluster Parameters::
    make_option(c("--build_manifest"), action="store_true", default=opt$build_manifest, 
                help="Boolean variable to build basic manifest [default= %default]", metavar="boolean"),
    make_option(c("--run_improbe"), action="store_true", default=opt$run_improbe, 
                help="Boolean variable to run improbe [default= %default]", metavar="boolean"),
    
    # Process Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
    # Run-time Files::
    make_option(c("--opt_csv"), type="character", default=opt$opt_csv, 
                help="Unused variable opt_csv [default= %default]", metavar="character"),
    make_option(c("--par_csv"), type="character", default=opt$par_csv, 
                help="Unused variable par_csv [default= %default]", metavar="character"),
    make_option(c("--time_csv"), type="character", default=opt$time_csv, 
                help="Unused variable time_csv [default= %default]", metavar="character"),
    make_option(c("--time_org_txt"), type="character", default=opt$time_org_txt, 
                help="Unused variable time_org_txt [default= %default]", metavar="character"),
    
    # Run=time Options::
    make_option(c("--fresh"), action="store_true", default=opt$fresh, 
                help="Boolean variable to run a fresh build [default= %default]", metavar="boolean"),
    make_option(c("--reload"), action="store_true", default=opt$reload, 
                help="Boolean variable reload intermediate files (for testing). [default= %default]", metavar="boolean"),
    
    # verbose::
    make_option(c("-v", "--verbose"), type="integer", default=opt$verbose, 
                help="0-5 (5 is very verbose) [default= %default]", metavar="integer")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Program Initialization::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_reqs <- c('runMode','prgmTag','scrDir','datDir','exePath')
opt_reqs <- c('outDir','impDir','Species',
              'genome_build','platform','version','bsmap_exe',
              'Rscript','verbose')

par$gen_src_dir <- load_source_files(dir=par$scrDir, verbose=opt$verbose)

opt <- program_init(name=par$prgmTag,
                    opts=opt, opt_reqs=opt_reqs, 
                    pars=par, par_reqs=par_reqs,
                    libs=TRUE,rcpp=FALSE,
                    verbose=opt$verbose,vt=3,tc=0,tt=NULL)

par_tib <- par %>%
  unlist(recursive = TRUE) %>%
  dplyr::bind_rows() %>% 
  tidyr::gather("Params", "Value")

opt_tib <- opt %>%
  dplyr::bind_rows() %>% 
  tidyr::gather("Option", "Value")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Parse List Options
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run <- NULL
pTracker <- timeTracker$new()

run$image_key <- "bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse_Centos"
run$image_ver <- "v.1.25"
run$doc_shell <- "run_improbe.sh"
run$doc_image <- glue::glue("{run$image_key}.{run$image_ver}")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Pre-defined Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Checking pre-defined files.{RET}"))

# Define Pre-built improbe directories and files::
#   - Using split files now instead of two single large files...
run$cgn_seq_dir <- 
  file.path(opt$impDir, "scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49-split")

stopifnot(dir.exists(run$cgn_seq_dir))

run$cgn_bed_dir <- file.path(opt$impDir, "scratch/cgnDB/dbSNP_Core4/design-input/min")
run$cgn_bed_tsv <- file.path(run$cgn_bed_dir, paste(opt$genome_build,"cgn.min.txt.gz", sep="."))
run$canonical_csv <- file.path(par$datDir, "manifest/cgnDB/canonical.cgn-top-grp.csv.gz")

stopifnot(dir.exists(run$cgn_bed_dir))
stopifnot(file.exists(run$cgn_bed_tsv))
stopifnot(file.exists(run$canonical_csv))

# Change these to options::
if (is.null(opt$gs_ctl_csv))
  opt$gs_ctl_csv  <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
if (is.null(opt$ses_ctl_csv))
  opt$ses_ctl_csv <- file.path(par$topDir, "data/manifests/methylation/Sesame/EPIC-B4-BP4.manifest.sesame-base.controls-only.csv.gz")

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Checking pre-defined files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Intermediate Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Defining Intermediate Run Time Files...{RET}"))

#
# MOST declaration of intermediate directories is no longer needed
#   They build themselves....
#

# Manifest Directory::
run$man_dir <- file.path(opt$outDir, 'man')
run$aqp_man_csv  <- file.path(run$man_dir, paste(opt$runName,"aqp-pass.manifest-sesame.csv.gz", sep="."))

# Annotation Directory::
#
run$ann_dir <- file.path(opt$outDir, 'ann')
run$ann_int_csv  <- file.path(run$ann_dir,paste(opt$runName,'cpg-pass.annotation.csv.gz', sep='.'))

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Defining Intermediate Run Time Files.{RET}{RET}"))

#
#
# Genomes Manifest Generation Parameters::
#  - TBD:: These should be moved to some config file...
#
#

# For accumilating failed probes and why::
error_ledgar <- NULL

run$gen_nrec <- 0
run$gen_key <- "Genome_Key"

# Field (key) Parameters:: general
run$del     <- "_"
run$ids_key <- "Prb_Key"
run$unq_key <- "Prb_Key_Unq"

run$add_key <- "Address"
run$din_key <- "Ord_Din"
run$des_key <- "Ord_Des"
run$map_key <- "Ord_Map"
run$prb_key <- "Ord_Prb"

run$bsp_srd <- "Bsp_FR"
run$bsp_cos <- "Bsp_CO"
run$pos_key <- "Bsp_Pos"
run$chr_key <- "Bsp_Chr"

run$Cgn_Int <- "Cgn_Int"
run$Can_Cgn <- "Can_Cgn"
run$Ord_Cgn <- "Ord_Cgn"
run$Bsp_Cgn <- "Bsp_Cgn"
run$Imp_Cgn <- "Imp_Cgn"

run$out_col <- c(run$ids_key, run$add_key, run$des_key,
                 run$din_key, run$map_key, run$prb_key)
run$unq_col <- c(run$din_key, run$map_key, run$Cgn_Int)

# Default run parameters by workflow::
run$bsp_full   <- FALSE
run$bsp_sort   <- TRUE
run$bsp_light  <- TRUE
run$bsp_merge  <- FALSE
run$cgn_merge  <- FALSE
run$cgn_join   <- "inner"
run$seq_suffix <- "probe-subseq"
run$seq_idxA   <- 1
run$seq_idxB   <- 1
run$seq_pattern_U <- "-probe_U49_cgn-table.tsv.gz"
run$seq_pattern_M <- "-probe_M49_cgn-table.tsv.gz"

# Field Parameters:: s-improbe
run$ext_seq="Ext_Forward_Seq"
run$iup_seq="Iupac_Forward_Sequence"
run$imp_seq="Forward_Sequence"

# Field Parameters:: r-improbe
run$srsplit <- TRUE
run$srd_key <- "Bsp_FR"
run$srd_str <- "FR"

run$cosplit <- TRUE
run$cos_key <- "Bsp_CO"
run$cos_str <- "CO"

run$re_load <- TRUE

opt$build_manifest <- FALSE
opt$run_improbe    <- FALSE
run$imp_level <- 3
run$add_inf <- TRUE

run$join_new_vec <- c()
run$join_old_vec <- c()

run$c_improbe <- TRUE


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#            0.0 Load any pre-defined Standard Manifest to be added::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (!is.null(opt$mans)) {
  #
  # TBD:: User should be able to rebuild an existing or old manifest,
  #  or add manifests together...
  #
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#          0.1 Load any pre-defined Noob-Masked Manifest to be added::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# noob_ctl_tib <- NULL
# if (!is.null(opt$noob)) {
#   noob_ctl_tib <- noob_mask(noob_csv = opt$noob, 
#                             ctl_csv = opt$ses_ctl_csv, 
#                             verbose=opt$verbose, tt=pTracker)
# }

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      0.2 Load any other pre-defined data
#                        dbCGN, improbe, imGenomes data::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imGenome_tib <- load_imGenomes_table(dir = opt$gen_dir,
                                     genome_build = opt$genome_build, 
                                     ret_list = FALSE,
                                     verbose = opt$verbose, tt=pTracker)


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   1.0 AQP Address Manifest Workflow: 
#                           Order/Match/AQP/PQC
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker$addFile(opt$time_org_txt)

ord_tib <- NULL
ord_tib <- 
  aqp_mapping_workflow(ord = opt$ords,
                       mat = opt$mats,
                       aqp = opt$aqps,
                       
                       prb_key = run$prb_key,
                       add_key = run$add_key,
                       des_key = run$des_key,
                       din_key = run$din_key,
                       ids_key = run$ids_key,
                       
                       out_dir = opt$outDir,
                       out_col = run$out_col,
                       run_tag = opt$runName,
                       re_load = run$re_load,
                       pre_tag = pTracker$file_vec,
                       
                       verbose=opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    2.0 Align All Probe Sequence:: BSMAP
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsp_tib <- bsp_mapping_workflow(ref_fas = NULL,
                                ref_tib = imGenome_tib,
                                can_fas = NULL,
                                can_tib = ord_tib,
                                
                                cgn_src = run$cgn_bed_tsv,
                                
                                ids_key = run$ids_key,
                                unq_key = run$unq_key,
                                prb_key = run$prb_key,
                                des_key = run$des_key,
                                din_key = run$din_key,
                                
                                join_key  = run$ids_key,
                                join_type = "inner",
                                
                                sort    = run$bsp_sort,
                                full    = run$bsp_full,
                                merge   = run$bsp_merge,
                                
                                light   = run$bsp_light,
                                reload  = opt$reload,
                                retData = FALSE,
                                
                                bsp_exe = opt$bsmap_exe,
                                bsp_opt = opt$bsmap_opt,
                                
                                out_dir = opt$outDir,
                                out_col = run$out_col,
                                run_tag = opt$runName,
                                re_load = run$re_load,
                                pre_tag = pTracker$file_vec,
                                
                                verbose=opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#         3.0 Intersect Sequences Address and improbe:: U49/M49
#                         CGN Mapping Workflow()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# NOTE: McMaster10Kselection = 822s
seq_tib <- 
  seq_mapping_workflow(ord_tib = ord_tib,
                       
                       seq_dir   = run$cgn_seq_dir,
                       pattern_u = run$seq_pattern_U, 
                       pattern_m = run$seq_pattern_M,
                       
                       prb_key = run$prb_key,
                       add_key = run$add_key,
                       des_key = run$des_key,
                       din_key = run$din_key,
                       ids_key = run$ids_key,
                       
                       prefix = opt$runName,
                       suffix = run$seq_suffix, 
                       
                       idxA = run$seq_idxA,
                       idxB = run$seq_idxB,
                       
                       reload   = opt$reload,
                       parallel = opt$parallel,
                       
                       del = run$del,
                       
                       out_dir = opt$outDir,
                       out_col = run$out_col,
                       run_tag = opt$runName,
                       re_load = run$re_load,
                       pre_tag = pTracker$file_vec,
                       
                       verbose=opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       4.0 Analyze and Assign Cgn:: 
#                      CGN-Map/BSMAP/dbGCGN look-up
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cgn_tib <- 
  cgn_mapping_workflow(ord_tib = ord_tib,
                       bsp_tib = bsp_tib,
                       seq_tib = seq_tib,
                       
                       ids_key = run$ids_key,
                       des_key = run$des_key,
                       din_key = run$din_key,
                       map_key = run$map_key,
                       
                       Cgn_Int = run$Cgn_Int,
                       Can_Cgn = run$Can_Cgn,
                       Ord_Cgn = run$Ord_Cgn,
                       Bsp_Cgn = run$Bsp_Cgn,
                       Imp_Cgn = run$Imp_Cgn,
                       
                       can_csv = run$canonical_csv,
                       
                       join    = run$cgn_join,
                       merge   = run$cgn_merge,
                       retData = FALSE,
                       
                       out_dir = opt$outDir,
                       out_col = run$out_col,
                       unq_col = run$unq_col,
                       run_tag = opt$runName,
                       re_load = run$re_load,
                       pre_tag = pTracker$file_vec,
                       
                       verbose=opt$verbose, tt=pTracker)



print(ord_tib)
print(bsp_tib)
print(seq_tib)
print(cgn_tib)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
# Ending the first of this script here and developing the second half in a 
#   completely speparte script. This is due to some naming convention 
#   updates. Its just easier to keep them separated!!!
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   5.0 Probe Design Validation via imGenome:: 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
#
# The following sections have been moved to seperate script: 
#
#  Probe Design/Validation/Sequence Extraction::
#    c_improbe_workflow()
#    r_improbe_workflow()
#    s_improbe_workflow()
#
#  Manifest Final Generation()
#    build_manifest()
#    build_sesame_manifest()
#    build_genome_studio_manifest()
#    build_minfi_manifest()
#
#  Auxilary Annotation::
#
#    add_bed_annotaitons()
#
#  Sesame Cross Validation::
#  Docker Image Support::
#
# New Script: "Infinium_Methylation_Workhorse/scripts/R/workhorse/workhorse_main_dev_latest.part2.R"
#
#

#
# Build individual parts of the template sequence:: 
#
#                                            iupac
#     up01.up02...up11.up12.up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn12.dn11.dn02.dn01
#                                       Nxb [  C   G  ] Nxb
#
# Probe Design Formulas::
#
# Inf1C                               Nxb60* up61------------------dn58
# Inf2C                                     ext61* dn61-------------------dn12
#
# Inf1O                up12------------------up61 Nxb61*
# Inf2O           up11------------------up60 ext61*
#
#


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
