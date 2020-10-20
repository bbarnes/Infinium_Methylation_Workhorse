
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

# Docker Image Parameters::
par$docker_image <- 'bbarnesimdocker/im_workhorse:improbe.v1.1'
par$docker_shell <- 'run_improbe.sh'

# Program Parameters::
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'probe_design'
par$prgmTag <- 'improbe_main_non-CpG'
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

opt$snp_src_tsv <- NULL

# opt$cpg_s48_tsv <- NULL
# opt$cpg_top_tsv <- NULL
# 
# opt$cpg_pos_tsv <- NULL
# opt$cph_pos_tsv <- NULL
# opt$snp_pos_tsv <- NULL

# Platform/Method Options::
# opt$genomeBuild <- NULL
# opt$platform    <- NULL
# opt$version     <- NULL

# Run Options::
opt$fresh <- FALSE
par$retData <- FALSE

opt$matFormat <- 'new'
opt$ordFormat <- 'old'

opt$matSkip <- 40
opt$ordSkip <- 8
opt$pqcSkip <- 7

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
  
  locIdatDir <- file.path(par$topDir, 'data/idats')
  opt$outDir <- file.path(par$topDir, 'scratch')
  opt$impDir <- file.path(par$topDir, 'data/improbe')
  
  opt$cpg_s48_tsv <- file.path(opt$impDir, 'designOutput_21092020/seq48U/gz/seq48U-GRCh36-38-10-21092020.unq.noHeader.seq-sorted.tsv.gz')
  opt$cpg_s48_tsv <- file.path(opt$impDir, 'designOutput_21092020/seq48U/un/seq48U-GRCh36-38-10-21092020.unq.noHeader.seq-sorted.tsv')
  
  opt$write_full <- FALSE
  opt$write_base <- FALSE
  
  #
  # Pre-defined local options runTypes::
  #
  par$local_runType <- 'GENK'
  par$local_runType <- 'NZT'
  par$local_runType <- 'COVIC'
  par$local_runType <- 'GRCm38'

  opt$snp_mus_src_tsv <- '/Users/bretbarnes/Documents/tmp/LifeEpigentics/data/20190228_input_files/snp_input.txt.gz'
  opt$snp_hsa_src_tsv <- '/Users/bretbarnes/Documents/data/manifests/raw/manifests/methylation/rs-repair/rs.swapDatabase.txt'
  
  opt$fresh <- TRUE
  opt$parallel <- TRUE
  opt$runName  <- par$prgmTag
  
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
    
    make_option(c("--cpg_s48_tsv"), type="character", default=opt$cpg_s48_tsv, 
                help="Seq48U Match file(s) (tab seperated) [default= %default]", metavar="character"),
    make_option(c("--cpg_top_tsv"), type="character", default=opt$cpg_top_tsv, 
                help="Top Sequence CGN Match file(s) (tab seperated) [default= %default]", metavar="character"),
    
    make_option(c("--cpg_pos_tsv"), type="character", default=opt$cpg_pos_tsv, 
                help="CpG Position TSV File [default= %default]", metavar="character"),
    make_option(c("--cph_pos_tsv"), type="character", default=opt$cph_pos_tsv, 
                help="CpH Position TSV File [default= %default]", metavar="character"),
    make_option(c("--snp_pos_tsv"), type="character", default=opt$snp_pos_tsv, 
                help="SNP Position TSV File [default= %default]", metavar="character"),
    
    # Platform/Method Options::
    # make_option(c("--genomeBuild"), type="character", default=opt$genomeBuild, 
    #             help="Genome Build (e.g. GRCh36, GRCh37, GRCh38, GRCm38) [default= %default]", metavar="character"),
    # make_option(c("--platform"), type="character", default=opt$platform, 
    #             help="Platform (e.g. HM450, EPIC, LEGX, NZT, COVIC) [default= %default]", metavar="character"),
    # make_option(c("--version"), type="character", default=opt$version, 
    #             help="Manifest Version (e.g. B0,B1,B2,B3,B4,C0) [default= %default]", metavar="character"),
    
    # File formats
    make_option(c("--ordFormat"), type="character", default=opt$ordFormat, 
                help="Order File Format [default= %default]", metavar="character"),
    make_option(c("--matFormat"), type="character", default=opt$matFormat, 
                help="Match File Format [default= %default]", metavar="character"),
    
    # Skip Lines
    make_option(c("--matSkip"), type="integer", default=opt$matSkip, 
                help="Match Skip Lines Count [default= %default]", metavar="integer"),
    make_option(c("--ordSkip"), type="integer", default=opt$ordSkip, 
                help="Order Skip Lines Count [default= %default]", metavar="integer"),
    make_option(c("--pqcSkip"), type="integer", default=opt$pqcSkip, 
                help="AQP/PQC Skip Lines Count [default= %default]", metavar="integer"),
    
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
    is.null(opt$snp_src_tsv) ||
    # is.null(opt$ords) || is.null(opt$mats) || 
    # is.null(opt$cpg_s48_tsv) || is.null(opt$cpg_top_tsv) ||
    # is.null(opt$genomeBuild) || is.null(opt$platform) || is.null(opt$version) ||
    is.null(opt$Rscript) ||
    is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )
  
  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  
  if (is.null(opt$outDir))    cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$impDir))    cat(glue::glue("[Usage]: impDir is NULL!!!{RET}"))
  
  if (is.null(opt$snp_src_tsv)) cat(glue::glue("[Usage]: snp_src_tsv is NULL!!!{RET}"))
  
  if (is.null(opt$ords)) cat(glue::glue("[Usage]: ords is NULL!!!{RET}"))
  if (is.null(opt$mats)) cat(glue::glue("[Usage]: mats is NULL!!!{RET}"))
  if (is.null(opt$aqps)) cat(glue::glue("[Usage]: aqps is NULL [Not required if pqcs defined]!!!{RET}"))
  if (is.null(opt$pqcs)) cat(glue::glue("[Usage]: pqcs is NULL [Not required if aqps defined]!!!{RET}"))
  if (is.null(opt$ctls)) cat(glue::glue("[Usage]: ctls is NULL [Not required]!!!{RET}"))
  if (is.null(opt$idat)) cat(glue::glue("[Usage]: idat is NULL [Not Required]!!!{RET}"))
  
  if (is.null(opt$cpg_s48_tsv)) cat(glue::glue("[Usage]: cpg_s48_tsv is NULL!!!{RET}"))
  if (is.null(opt$cpg_top_tsv)) cat(glue::glue("[Usage]: cpg_top_tsv is NULL!!!{RET}"))
  
  # if (is.null(opt$genomeBuild)) cat(glue::glue("[Usage]: genomeBuild is NULL!!!{RET}"))
  # if (is.null(opt$platform))    cat(glue::glue("[Usage]: platform is NULL!!!{RET}"))
  # if (is.null(opt$version))     cat(glue::glue("[Usage]: version is NULL!!!{RET}"))
  
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

par$anl_src_dir <- file.path(par$scrDir, 'annotation/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$anl_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$anl_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$anl_src_dir}!{RET}{RET}") )

cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

# opt$probe_type <- 'cg'
# opt$design_key <- 'Seq_ID'
# opt$design_seq <- 'Forward_Sequence'
# opt$design_seq <- 'Top_Sequence'
# opt$design_prb <- 'Probe_Type'
# opt$design_srs <- 'TB'
# opt$design_cos <- 'CO'
# 
# design_prb_sym <- rlang::sym(opt$design_prb)

# Input Definitions::
if (is.null(opt$ctls)) {
  opt$ctls <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
}

# ords_vec <- NULL
# mats_vec <- NULL
# aqps_vec <- NULL
# pqcs_vec <- NULL
# ctls_vec <- NULL
# idat_vec <- NULL
# if (!is.null(opt$ords)) ords_vec <- opt$ords %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
# if (!is.null(opt$mats)) mats_vec <- opt$mats %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
# if (!is.null(opt$aqps)) aqps_vec <- opt$aqps %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
# if (!is.null(opt$pqcs)) pqcs_vec <- opt$pqcs %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
# if (!is.null(opt$ctls)) ctls_vec <- opt$ctls %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
# if (!is.null(opt$idat)) idat_vec <- opt$idat %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
# 
# stopifnot(length(ords_vec)>0)
# stopifnot(length(mats_vec)>0)
# stopifnot(length(mats_vec)==length(ords_vec))
# 
# cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

# Output Definitions::
opt$outDir <- file.path(opt$outDir, par$prgmDir, opt$runName)
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)

# opt$manDir <- file.path(opt$outDir, 'manifest')
# if (!dir.exists(opt$manDir)) dir.create(opt$manDir, recursive=TRUE)
# 
# opt$prdDir <- file.path(opt$outDir, 'product')
# if (!dir.exists(opt$prdDir)) dir.create(opt$prdDir, recursive=TRUE)
# 
# 
# opt$genDir <- file.path(opt$outDir, 'genomic')
# if (!dir.exists(opt$genDir)) dir.create(opt$genDir, recursive=TRUE)

cat(glue::glue("[{par$prgmTag}]: Done. Building Output Directories.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Define Manifest Output Files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# man_out_key <- paste(opt$genomeBuild,opt$platform,opt$version, sep='_')
# gs_swap_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.GenomeStudio.cpg-sorted.csv', sep='.') )
# gz_swap_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.GenomeStudio.cpg-sorted.csv.gz', sep='.') )
# 
# ses_base_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.sesame-base.cpg-sorted.csv.gz', sep='.') )
# pos_base_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.sesame-base.pos-sorted.csv.gz', sep='.') )
# ses_mach_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.sesame-mach.cpg-sorted.csv.gz', sep='.') )
# ses_epic_csv <- file.path(opt$prdDir, paste(man_out_key,'manifest.sesame-epic.cpg-sorted.csv.gz', sep='.') )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Genome Studio Color Codes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

color_tib <- suppressMessages(suppressWarnings( readr::read_csv(file.path(par$datDir, 'params/GenomeStudioColors.csv')) ))
color_vec <- color_tib %>% dplyr::pull(Color) %>% as.vector()
color_len <- length(color_vec)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Historical Investigation:: SNP
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD::
#     1. Write improbe input tsv
#     2. Run Docker improbe.v1.1 on improbe input tsv
#     3. Build s48U
#     4. Build pos tsv
#     5. Reuse top tsv
#     6. Test clean_manifest_probes(SNP)
#

#
# MUS SNPs Top='/Users/bretbarnes/Documents/tmp/LifeEpigentics/data/20190228_input_files/snp_input.txt.gz'
#   opt$snp_mus_src_tsv='/Users/bretbarnes/Documents/tmp/LifeEpigentics/data/20190228_input_files/snp_input.txt.gz'
#   opt$snp_hsa_src_tsv='/Users/bretbarnes/Documents/data/manifests/raw/manifests/methylation/rs-repair/rs.swapDatabase.txt'
#
# MUS CpHs Top='/Users/bretbarnes/Documents/tmp/LifeEpigentics/data/20190228_input_files/cph_input.txt.gz'
#
# HSA SNPs Prb='data/manifests/raw/manifests/methylation/HumanMethylation450_15017482_v.1.2.rs-only.csv.gz'
# HSA SNPs Top='data/manifests/raw/manifests/methylation/rs-repair/rs.swapDatabase.txt'
#
# HSA CpHs Prb='data/manifests/raw/manifests/methylation/HumanMethylation450_15017482_v.1.2.ch.csv'
# HSA CpHs Top='data/manifests/raw/manifests/methylation/ch-repair/HumanMethylation450_15017482_v.1.2.ch.csv'
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Historical Investigation:: SNP
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$curDir <- file.path(opt$outDir,'snp')
if (!dir.exists(opt$curDir)) dir.create(opt$curDir, recursive=TRUE)
if (opt$fresh) list.files(opt$curDir, full.names=TRUE) %>% unlink()

opt$prdDir <- file.path(opt$curDir,'product')
if (!dir.exists(opt$prdDir)) dir.create(opt$prdDir, recursive=TRUE)
if (opt$fresh) list.files(opt$prdDir, full.names=TRUE) %>% unlink()

genomeBuilds <- NULL
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Historical Investigation:: SNP-MUS
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

snp_mus_list <- loadGRCm38_SNP(
  file=opt$snp_mus_src_tsv,outDir=opt$curDir,
  verbose=opt$verbose,vt=1,tc=0,tt=pTracker)

genomeBuild='GRCm38'
genomeBuilds <- paste(genomeBuilds,genomeBuild, sep='-')
snp_mus_imp_tib <- improbe_docker(
  dir=opt$curDir,
  file=base::basename(snp_mus_list$imp), name=genomeBuild,
  image=par$docker_image, shell=par$docker_shell,
  verbose=opt$verbose,vt=1,tc=0,tt=pTracker)

snp_mus_prb_tib <- desSeq_to_prbs(
  snp_mus_list$src, 
  idsKey="Seq_ID",seqKey="IUPAC_Sequence",prbKey="Probe_Type",
  strsSR="FR",strsCO="CO",
  verbose=opt$verbose,vt=1,tc=0,tt=pTracker)

snp_mus_all_tib <- 
  dplyr::inner_join(snp_mus_prb_tib,snp_mus_imp_tib,
                    by=c("Seq_ID","SR_Str"="Methyl_Allele_FR_Strand","CO_Str"="Methyl_Allele_CO_Strand")) %>%
  dplyr::mutate(Seq48U=stringr::str_sub(UnMethyl_Probe_Sequence, 2,49)) %>% 
  dplyr::mutate(Chromosome=as.character(Chromosome),
                FR_Str=SR_Str,
                SR_Str=stringr::str_sub(TB_Strand,1,1))

# Per Species Outputs::
#
snp_mus_pos_tsv <- file.path(opt$prdDir,paste(genomeBuild,'improbe-pos.snp.tsv.gz',sep='.'))
snp_mus_pos_tib <- snp_mus_all_tib %>% 
  dplyr::select(Seq_ID,Chromosome,Coordinate) %>%
  dplyr::distinct() %>% dplyr::arrange(Seq_ID)
readr::write_tsv(snp_mus_pos_tib,snp_mus_pos_tsv)

snp_mus_top_tsv <- file.path(opt$prdDir,paste(genomeBuild,'improbe-top.snp.tsv.gz',sep='.'))
snp_mus_top_tib <- snp_mus_all_tib %>% 
  dplyr::select(Seq_ID,Top_Sequence) %>% 
  dplyr::distinct() %>% dplyr::arrange(Seq_ID)
readr::write_tsv(snp_mus_top_tib,snp_mus_top_tsv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Historical Investigation:: SNP-HSA
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

snp_hsa_list <- loadGRCh37_SNP(
  file=opt$snp_hsa_src_tsv,outDir=opt$curDir,
  verbose=opt$verbose,vt=1,tc=0,tt=pTracker)

genomeBuild='GRCh37'
genomeBuilds <- paste(genomeBuilds,genomeBuild, sep='-')
snp_hsa_imp_tib <- improbe_docker(
  dir=opt$curDir,
  file=base::basename(snp_hsa_list$imp), name=genomeBuild,
  image=par$docker_image, shell=par$docker_shell,
  verbose=opt$verbose,vt=1,tc=0,tt=pTracker)

snp_hsa_prb_tib <- desSeq_to_prbs(
  snp_hsa_list$src, 
  idsKey="Seq_ID",seqKey="IUPAC_Sequence",prbKey="Probe_Type",
  strsSR="FR",strsCO="CO",
  verbose=opt$verbose,vt=1,tc=0,tt=pTracker)

snp_hsa_all_tib <- 
  dplyr::inner_join(snp_hsa_prb_tib,snp_hsa_imp_tib,
                    by=c("Seq_ID","SR_Str"="Methyl_Allele_FR_Strand","CO_Str"="Methyl_Allele_CO_Strand")) %>%
  dplyr::mutate(Seq48U=stringr::str_sub(UnMethyl_Probe_Sequence, 2,49)) %>%
  dplyr::mutate(Chromosome=as.character(Chromosome),
                FR_Str=SR_Str,
                SR_Str=stringr::str_sub(TB_Strand,1,1))

# Per Species Outputs::
#
snp_hsa_pos_tsv <- file.path(opt$prdDir,paste(genomeBuild,'improbe-pos.snp.tsv.gz',sep='.'))
snp_hsa_pos_tib <- snp_hsa_all_tib %>% 
  dplyr::select(Seq_ID,Chromosome,Coordinate) %>%
  dplyr::distinct() %>% dplyr::arrange(Seq_ID)
readr::write_tsv(snp_hsa_pos_tib,snp_hsa_pos_tsv)

snp_hsa_top_tsv <- file.path(opt$prdDir,paste(genomeBuild,'improbe-top.snp.tsv.gz',sep='.'))
snp_hsa_top_tib <- snp_hsa_all_tib %>% 
  dplyr::select(Seq_ID,Top_Sequence) %>% 
  dplyr::distinct() %>% dplyr::arrange(Seq_ID)
readr::write_tsv(snp_hsa_top_tib,snp_hsa_top_tsv)


#
# We need to match based on Seq_ID, FR, CO then, build Match Seq (Seq48U)...
#

#
# Produce Required Files::
#
#   - opt$cpg_s48_tsv -> opt$snp_s48_tsv
#
# By Species::
#   - opt$cpg_top_tsv -> opt$snp_top_tsv
#   - opt$cpg_pos_tsv -> opt$snp_pos_tsv
#
#   - opt$cpg_pos_tsv -> opt$snp_ann_bed
#                      +
#   - opt$cpg_pos_tsv -> opt$snp_scr_tsv
#

genomeBuilds <- genomeBuilds %>% stringr::str_remove('^-')

snp_all_tib <- dplyr::bind_rows(snp_hsa_all_tib,snp_mus_all_tib)

snp_all_s48_tsv <- file.path(opt$prdDir,paste(genomeBuilds,'improbe-s48.snp.tsv.gz',sep='.'))
snp_all_s48_tib <- snp_all_tib %>% 
  dplyr::select(Seq_ID,SR_Str,CO_Str,Seq48U) %>% 
  dplyr::rename(Strand_SR=SR_Str,Strand_CO=CO_Str) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Seq48U)
readr::write_tsv(snp_all_s48_tib,snp_all_s48_tsv)

#
# MUS SNPs Top='/Users/bretbarnes/Documents/tmp/LifeEpigentics/data/20190228_input_files/snp_input.txt.gz'
#
# MUS CpHs Top='/Users/bretbarnes/Documents/tmp/LifeEpigentics/data/20190228_input_files/cph_input.txt.gz'
#
# HSA SNPs Prb='data/manifests/raw/manifests/methylation/HumanMethylation450_15017482_v.1.2.rs-only.csv.gz'
# HSA SNPs Top='data/manifests/raw/manifests/methylation/rs-repair/rs.swapDatabase.txt'
#
# HSA CpHs Prb='data/manifests/raw/manifests/methylation/HumanMethylation450_15017482_v.1.2.ch.csv'
# HSA CpHs Top='data/manifests/raw/manifests/methylation/ch-repair/HumanMethylation450_15017482_v.1.2.ch.csv'
#


#
# tag=cgn_[TB][CO]
#
# TBD:: Split into
#   - s48U,tag,Probe_Type,topSeq
#   - chr,beg,end,tag,[FR],genome_build
#   - tag,[SeqU1,SeqM1],[scores]
#
#   - chr,beg,end,tag,[FR],genome_build,fwdSeq
#




if (FALSE) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #             Verify all probes in manifest with proper names::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  #
  # SNP:: may need to test sidx=2|1 & plen=50 for SNPs...
  #
  opt$fresh <- TRUE
  
  sidx <- 0
  sidx <- 2
  pmax <- 50
  man_raw_tib0 <- decodeAqpPqcWrapper(
    ord_vec=ords_vec,mat_vec=mats_vec,aqp_vec=aqps_vec,pqc_vec=pqcs_vec, 
    platform=opt$platform,version=opt$version,
    ordFormat=opt$ordFormat,ordSkip=opt$ordSkip,
    matFormat=opt$matFormat,matSkip=opt$matSkip,
    pqcSkip=opt$pqcSkip,
    sidx=sidx, plen=pmax-2,
    name=opt$runName,outDir=paste(opt$manDir,sidx,sep='_'),
    fresh=opt$fresh,full=par$retData,trim=TRUE,
    verbose=opt$verbose,vt=1,tc=0,tt=pTracker) %>%
    dplyr::arrange(Seq_ID)
  
  raw_rep_tib <- manifestCheckSummary(man_raw_tib0, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)
  
  opt$probe_type <- 'rs'
  opt$snp_s48_tsv <- '/Users/bretbarnes/Documents/scratch/probe_design/improbe_main_non-CpG/snp/product/GRCm38-GRCh37.improbe-s48.snp.tsv.gz'
  opt$snp_top_tsv <- '/Users/bretbarnes/Documents/scratch/probe_design/improbe_main_non-CpG/snp/product/GRCm38.improbe-top.snp.tsv.gz'
  
  man_snp_prb_tib0 <- clean_manifest_probes(
    tib=man_raw_tib0,s48_tsv=opt$snp_s48_tsv,top_tsv=opt$snp_top_tsv,
    name=opt$runName,outDir=paste(opt$outDir,sidx,sep='_'),
    
    design_key=opt$design_key,
    design_seq=opt$design_seq,
    design_prb=opt$design_prb,
    probe_type=opt$probe_type,
    design_srs=opt$design_srs,
    design_cos=opt$design_cos,
    parallel=opt$parallel,fresh=opt$fresh,
    
    verbose=opt$verbose+10,vt=1,tc=0,tt=pTracker)
  
  #
  # Standard::
  #
  
  opt$probe_type <- 'rs'
  opt$snp_s48_tsv <- '/Users/bretbarnes/Documents/scratch/probe_design/improbe_main_non-CpG/snp/product/GRCm38-GRCh37.improbe-s48.snp.tsv.gz'
  opt$snp_top_tsv <- '/Users/bretbarnes/Documents/scratch/probe_design/improbe_main_non-CpG/snp/product/GRCm38.improbe-top.snp.tsv.gz'
  
  man_snp_prb_tib <- clean_manifest_probes(
    tib=man_raw_tib,
    s48_tsv=opt$snp_s48_tsv,top_tsv=opt$snp_top_tsv,
    name=opt$runName,outDir=opt$outDir,
    
    design_key=opt$design_key,
    design_seq=opt$design_seq,
    design_prb=opt$design_prb,
    probe_type=opt$probe_type,
    design_srs=opt$design_srs,
    design_cos=opt$design_cos,
    parallel=opt$parallel,fresh=opt$fresh,
    
    verbose=opt$verbose+10,vt=1,tc=0,tt=pTracker)
  
  
  
  #
  # TBD:: Format all CpH and SNP files to allow clean_manifest_probes()
  #
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Historical Investigation:: SNP
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  #
  # TBD::
  #     1. Write improbe input tsv
  #     2. Run Docker improbe.v1.1 on improbe input tsv
  #     3. Build s48U
  #     4. Build pos tsv
  #     5. Reuse top tsv
  #     6. Test clean_manifest_probes(SNP)
  #
  
  snp_src_col <- c('Seq_ID','Probe_Type_Str','Design_Org_Seq',
                   'Genome_Build','Chromosome','Coordinate','AlleleA','AlleleB',
                   'IntA','IntB','Ann_Str')
  snp_src_tsv <- '/Users/bretbarnes/Documents/tmp/LifeEpigentics/data/20190228_input_files/snp_input.txt.gz'
  
  snp_inp_dir <- file.path(opt$impDir,'designOutput_21092020/snp')
  snp_inp_fon <- paste(opt$genomeBuild,'improbeDesignInput.tsv.gz', sep='.')
  snp_inp_tsv <- file.path(snp_inp_dir, snp_inp_fon)
  
  snp_inp_tib <- 
    suppressMessages(suppressWarnings( readr::read_tsv(snp_src_tsv, col_names=snp_src_col) )) %>% 
    dplyr::mutate(
      Iupac_Allele=mapDIs(paste0(AlleleA,AlleleB)),
      Pre_Seq=stringr::str_remove(Design_Org_Seq, '\\[.*$'),
      Pos_Seq=stringr::str_remove(Design_Org_Seq, '^.*\\]'),
      Pos_Nuc=stringr::str_sub(Pos_Seq,1,1),
      Pos_Seq=stringr::str_sub(Pos_Seq,2),
      Sequence=stringr::str_to_upper(paste0(Pre_Seq,'[CG]',Pos_Seq,'N')),
      IUPAC_Sequence=paste0(Pre_Seq,'[',Iupac_Allele,Pos_Nuc,']',Pos_Seq,'N'),
      Chromosome=stringr::str_remove(Chromosome,'^chr'),
      Genome_Build=dplyr::case_when(
        Genome_Build=='hg18' ~ "GRCh36",
        Genome_Build=='hg19' ~ "GRCh37",
        Genome_Build=='hg38' ~ "GRCh38",
        Genome_Build=='mm10' ~ "GRCm38",
        TRUE ~ Genome_Build
      ),
      CpG_Island="FALSE"
    ) %>% 
    dplyr::select(Seq_ID,Sequence,Genome_Build,Chromosome,Coordinate,
                  CpG_Island,IUPAC_Sequence) %>%
    dplyr::arrange(Seq_ID)
  
  readr::write_tsv(snp_inp_tib,snp_inp_tsv)
  
  
  #
  # OLD STUFF::
  #
  # snp_top_tib <- readr::read_csv(opt$snp_des_csv) %>% 
  #   dplyr::select(Seq_ID, IUPAC_Forward_Sequence, Probe_Type) %>% 
  #   dplyr::mutate(Seq_ID=stringr::str_replace(Seq_ID,':','-')) %>%
  #   dplyr::distinct()
  # 
  # snp_pos_col <- 
  #   cols(Seq_ID=col_character(),Gen_Chr=col_character(),Gen_Pos=col_double() )
  # snp_pos_tib <- readr::read_tsv(opt$snp_pos_tsv, 
  #                                col_names=names(snp_pos_col$cols), 
  #                                col_types=snp_pos_col)
  # 
  # snp_pos_tib %>% dplyr::anti_join(snp_top_tib, by="Seq_ID")
  # snp_top_tib %>% dplyr::anti_join(snp_pos_tib, by="Seq_ID")
  # snp_top_tib %>% dplyr::anti_join(snp_org_tib, by="Seq_ID")
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Historical Investigation:: CpH
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  readr::read_csv(opt$cph_des_csv) %>% 
    dplyr::select(Seq_ID, IUPAC_Forward_Sequence, Probe_Type) %>% 
    # dplyr::group_by_all() %>% dplyr::summarise(Count=n())
    dplyr::distinct()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #             Rebuild full manifest with updated clean probes::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  man_prb_tib <- dplyr::bind_rows(  
    man_cpg_prb_tib %>%
      dplyr::distinct(M,U, .keep_all=TRUE) %>%
      dplyr::add_count(Seq_ID,SR_Str,CO_Str,Infinium_Design, name='Rep_Max') %>%
      dplyr::group_by(Seq_ID,SR_Str,CO_Str,Infinium_Design) %>%
      dplyr::mutate(
        Rep_Cnt=dplyr::row_number(),
        IlmnID=paste0(Seq_ID,'_',SR_Str,CO_Str,Infinium_Design,Rep_Cnt),
        ValidID=TRUE
      ) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(IlmnID, dplyr::everything()),
    # man_unk_tib,
    NULL
  ) %>% dplyr::arrange(IlmnID) %>%
    dplyr::select(IlmnID,ValidID,dplyr::everything())
}



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
