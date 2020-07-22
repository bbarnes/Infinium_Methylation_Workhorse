
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

par$prgmDir <- 'manifests'
par$prgmTag <- paste(par$prgmDir,'aqp_to_manifest_mm10', sep='_')
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# Illumina based directories::
par$macDir <- '/Users/bbarnes/Documents/Projects/methylation/tools'
par$lixDir <- '/illumina/scratch/darkmatter'

par$retData <- FALSE

# Executables::
opt$Rscript <- NULL

# Directories::
opt$outDir     <- NULL
opt$datDir     <- NULL

# Required Inputs::
opt$ords <- NULL
opt$mats <- NULL
opt$aqps <- NULL
opt$pqcs <- NULL

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
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  
  # Default Parameters for local Mac::
  par$runMode    <- args.dat[1]
  par$srcDir     <- file.path(par$macDir, 'Infinium_Methylation_Workhorse')
  par$scrDir     <- file.path(par$srcDir, 'scripts')
  par$exePath    <- file.path(par$scrDir, 'R', par$prgmDir, paste0(par$prgmTag,'.R'))
  
  par$prgmTag <- base::sub('\\.R$', '', base::basename(par$exePath))
  par$locPath <- base::dirname(par$exePath)
  par$scrDir  <- base::dirname(base::normalizePath(par$locPath) )
  par$srcDir  <- base::dirname(base::normalizePath(par$scrDir) )
  par$datDir  <- file.path(base::dirname(base::normalizePath(par$srcDir)), 'dat')
  
  # Platform/Method Options::
  opt$genomeBuild <- 'mm10'
  opt$platform    <- 'LEGX'
  opt$version     <- 'B3'
  
  opt$frmt_original <- TRUE
  
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
    (is.null(opt$aqps) && is.null(opt$pqcs)) ||
    # is.null(ctlCSV) ||
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
  if (is.null(opt$ctlCSV)) cat(glue::glue("[Usage]: ctlCSV is NULL [Not Required]!!!{RET}"))
  
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
cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )

par$man_src_dir <- file.path(par$scrDir, 'manifests/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$man_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$man_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$man_src_dir}!{RET}{RET}") )

par$prgm_src_dir <- file.path(par$scrDir,par$prgmDir, 'functions')
if (!dir.exists(par$prgm_src_dir)) stop(glue::glue("[{par$prgmTag}]: Program Source={par$prgm_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$prgm_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Program Source={par$prgm_src_dir}!{RET}{RET}") )

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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Build Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$frmt_original <- FALSE
isFull <- FALSE

ord_cnt <- length(ord_vec)

full_man_tib <- NULL
if (!is.null(opt$pqcs)) {
  for (idx in c(1:ord_cnt)) {
    full_man_tib <- full_man_tib %>% dplyr::bind_rows(
      decodeToManifest(ord=ord_vec[idx], mat=mat_vec[idx], pqc=pqc_vec[1],
                       platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
                       original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=idx, AQP=idx)
    )
  }
  # pqcs_man_tib <- dplyr::bind_rows(
  #   decodeToManifest(ord=ord_bp1_csv, mat=mat1_tsv, pqc=pqc_tsv,
  #                    platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
  #                    original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=1, AQP=1),
  #   decodeToManifest(ord=ord_bp2_csv, mat=mat2_tsv, pqc=pqc_tsv, 
  #                    platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
  #                    original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=2, AQP=1),
  #   decodeToManifest(ord=ord_bp3_csv, mat=mat3_tsv, pqc=pqc_tsv, 
  #                    platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
  #                    original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=2, AQP=2),
  #   decodeToManifest(ord=ord_bp4_csv, mat=mat4_tsv, pqc=pqc_tsv, 
  #                    platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
  #                    original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=3, AQP=1)
  # )
} else {
  for (idx in c(1:ord_cnt)) {
    full_man_tib <- full_man_tib %>% dplyr::bind_rows(
      decodeToManifest(ord=ord_vec[idx], mat=mat_vec[idx], pqc=aqp_vec[idx],
                       platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
                       original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=idx, AQP=idx)
    )
  }
  # aqps_man_tib <- dplyr::bind_rows(
  #   decodeToManifest(ord=ord_bp1_csv, mat=mat1_tsv, aqp=aqp1_tsv,
  #                    platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
  #                    original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=1, AQP=1),
  #   decodeToManifest(ord=ord_bp2_csv, mat=mat2_tsv, aqp=aqp2_tsv, 
  #                    platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
  #                    original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=2, AQP=1),
  #   decodeToManifest(ord=ord_bp3_csv, mat=mat3_tsv, aqp=aqp3_tsv, 
  #                    platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
  #                    original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=2, AQP=2),
  #   decodeToManifest(ord=ord_bp4_csv, mat=mat4_tsv, aqp=aqp4_tsv, 
  #                    platform=opt$platform, version=opt$version, full=isFull, cleanAdds=TRUE,
  #                    original=opt$frmt_original, verbose=opt$verbose) %>% dplyr::mutate(BP=3, AQP=1)
  # )
}
full_man_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()

# Generate Unique Probe Stats::
full_man_tib %>% dplyr::select(Probe_ID) %>% dplyr::group_by(Probe_ID) %>% dplyr::mutate(Rep=row_number()) %>% dplyr::ungroup() %>% 
  dplyr::group_by(Rep) %>% dplyr::summarise(Count=n()) %>% print()

# Create Unique Names for duplicates::
uniq_man_tib <- full_man_tib %>% dplyr::group_by(Probe_ID) %>% 
  dplyr::mutate(Rep=paste0('r',row_number()) ) %>% 
  tidyr::unite(Probe_ID, Probe_ID,Rep, sep='_')

# Generate Unique Probe Stats::
# TBD:: Only add r# if rep > 1
uniq_man_tib %>% dplyr::select(Probe_ID) %>% dplyr::group_by(Probe_ID) %>% dplyr::mutate(Rep=row_number()) %>% dplyr::ungroup() %>% 
  dplyr::group_by(Rep) %>% dplyr::summarise(Count=n()) %>% print()

out_man_tib <- uniq_man_tib %>% dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, 
                                              AlleleA_Probe_Sequence,AlleleB_Probe_Sequence) %>%
  dplyr::mutate(Version=opt$tar_version)

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
  
} else {
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
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Sesame Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ses_man_tib <- dplyr::bind_rows(out_man_tib,ses_ctl_tib) %>% 
  # dplyr::select(-AlleleA_Probe_Sequence,-AlleleB_Probe_Sequence) %>%
  dplyr::arrange(Probe_ID)
readr::write_csv(ses_man_tib, ses_man_csv)

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
gs_body_tib <- out_man_tib %>% # head(n=1000) %>%
  dplyr::rename(IlmnID=Probe_ID,
                AddressA_ID=U,AlleleA_ProbeSeq=AlleleA_Probe_Sequence,
                AddressB_ID=M,AlleleB_ProbeSeq=AlleleB_Probe_Sequence,
                Infinium_Design_Type=DESIGN,Color_Channel=col) %>%
  dplyr::mutate(Name=IlmnID,Forward_Sequence='N',Genome_Build=opt$genomeBuild,CHR='chr1',MAPINFO=1,Source_Seq='N',Strand='F') %>%
  dplyr::select(IlmnID,Name,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,
                Infinium_Design_Type,Next_Base,Color_Channel,Forward_Sequence,Genome_Build,CHR,MAPINFO,Source_Seq,Strand) %>%
  dplyr::arrange(IlmnID)
ctl_head_df <- tibble::tibble(Head="[Controls]")

# Write Genome Studio CSV::
write_delim(x=gs_head_tib, path=gss_man_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=FALSE)
write_delim(x=gs_body_tib, path=gss_man_csv, delim=',', col_names=TRUE,  quote_escape=FALSE, na='', append=TRUE)
write_delim(x=ctl_head_df, path=gss_man_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)
write_delim(x=org_ctl_tib, path=gss_man_csv, delim=',', col_names=FALSE, quote_escape=FALSE, na='', append=TRUE)

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
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))


# End of file
