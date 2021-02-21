
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Genomic Ranges::
suppressWarnings(suppressPackageStartupMessages( base::require("GenomicRanges",quietly=TRUE) ))

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
par$prgmTag <- 'aqp_to_manifest'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))


# Executables::
opt$Rscript <- NULL

# BSMAP Parameters::
opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
opt$bsmap_exe <- "/Users/bretbarnes/Documents/programs/BSMAPz/bsmapz"

# Run Parameters::
opt$runName    <- NULL
opt$idatPrefix <- NULL

# Directories::
opt$outDir <- NULL
opt$impDir <- NULL
opt$annDir <- NULL
opt$genDir <- NULL
opt$manDir <- NULL

# Required Inputs::
opt$ords <- NULL
opt$mats <- NULL
opt$aqps <- NULL
opt$pqcs <- NULL
opt$ctls <- NULL
opt$idat <- NULL

# Platform/Method Options::
opt$genBuild <- NULL
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
opt$fresh   <- FALSE

# verbose Options::
opt$verbose <- 3


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Order/Match/AQP/PQC Expected Columns::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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
par$aqp_col <- par_cols$aqp$cols %>% names()
par$pqc_col <- par_cols$pqc$cols %>% names()

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
  opt$genDir  <- file.path(par$topDir, 'data/iGenomes/Homo_sapiens/NCBI')
  opt$idatDir <- file.path(par$topDir, 'data/idats')
  
  opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
  opt$bsmap_exe <- "/Users/bretbarnes/Documents/programs/BSMAPz/bsmapz"
  
  # Pre-defined local options runTypes::
  #
  par$local_runType <- NULL
  par$local_runType <- 'NZT'
  par$local_runType <- 'COVIC'
  par$local_runType <- 'GRCm10'

  if (par$local_runType=='COVIC') {
    opt$matFormat <- 'old'
    opt$ordFormat <- 'old'
    opt$matSkip <- 0
    opt$ordSkip <- 8
    
    opt$genBuild <- 'GRCh36'
    opt$genBuild <- 'GRCh38'
    opt$genBuild <- 'GRCh37'
    
    opt$platform    <- 'COVIC'
    opt$version     <- 'C0'
    opt$version     <- 'B1'
    opt$version     <- 'C2'
    opt$version     <- 'C11'
    opt$version     <- 'C12'
    
    opt$idatPrefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R01C01"
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    opt$pre_man_csv <- '/Users/bretbarnes/Documents/data/manifests/MethylationEPIC_v-1-0_B2.csv.gz'
    
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/COVID-19_HLA/AQP/COVIC-Host-Immune-Detection')
    opt$ords <- paste(
      file.path(par$aqpDir, 'COVID_EPIC_Round1.03172020.unique.order_AP.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20447043_probes.match.gz'),
      sep=',')
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'BS0032581-AQP.txt.gz'),
      sep=',')
    
    opt$pqcs <- NULL
    # opt$pqcs <- paste(
    #   file.path(par$aqpDir, 'BS0032581-AQP.txt.gz'),
    #   sep=',')
    
    par$idatsTopDir <- file.path(opt$idatDir,'idats_COVIC-Set1-15052020')
    opt$idat <- paste(
      file.path(par$idatsTopDir, '204500250013'),
      sep=',')
    
  } else if (par$local_runType=='GRCm10') {
    opt$platform <- 'LEGX'
    opt$version  <- 'C0'
    opt$version  <- 'C8'
    opt$version  <- 'C25'
    opt$genBuild <- 'GRCm38'
    opt$genBuild <- 'GRCm10'
    
    opt$genDir  <- file.path(par$topDir, 'data/iGenomes/Mus_musculus/NCBI')
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
    opt$cph_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.chn-sorted.tsv.gz') )
    opt$snp_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.snp-sorted.tsv.gz') )
    
    # opt$snp_des_csv <- file.path(opt$impDir, 'cph-snp-designs/LEGX_SpikeIn_Reorder-SNP-Only.designs.csv.gz')
    # opt$cph_des_csv <- file.path(opt$impDir, 'cph-snp-designs/LEGX_SpikeIn_Reorder-CpH-Only.designs.csv.gz')
    
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/LifeEpigentics/AQP')
    opt$ords <- paste(
      file.path(par$aqpDir, 'orders/Mus_musculus.order_BP1.csv.gz'),
      file.path(par$aqpDir, 'orders/Mus_musculus.order_BP2.csv.gz'),
      file.path(par$aqpDir, 'orders/mm10_LEGX_nonCpG_probes.Jan16-2020.order.csv.gz'),
      file.path(par$aqpDir, 'orders/LEGX_SpikeIn_Reorder-All-06052020.order.withHeader.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, 'BP1/20420178_AQP1_LifeEpigen_BP1.txt.gz'),
      file.path(par$aqpDir, 'BP2/20420260_AQP1_LifeEpigen_BP2.txt.gz'),
      file.path(par$aqpDir, 'BP3/20420260_AQP2_LifeEpigen_BP2.txt.gz'),
      file.path(par$aqpDir, 'BP4/20455357_AQP1_LifeEpigen_BP4.txt.gz'),
      sep=',')
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'AQP_Copy/BS0032527-AQP.txt.gz'),
      file.path(par$aqpDir, 'AQP_Copy/BS0032533-AQP.txt.gz'),
      file.path(par$aqpDir, 'AQP_Copy/BS0032545-AQP.txt.gz'),
      file.path(par$aqpDir, 'AQP_Copy/BS0032636-AQP.txt.gz'),
      sep=',')
    
    opt$pqcs <- paste(
      file.path(par$aqpDir, 'PQC/20042400_A_ProductQC.txt.gz'),
      sep=',')
    
    par$idatsTopDir <- file.path(opt$idatDir,'idats_ILMN_mm10_betaTest_17082020')
    opt$idat <- paste(
      file.path(par$idatsTopDir, '204637490025'),
      sep=',')
    
  } else if (par$local_runType=='NZT') {
    opt$genBuild <- 'GRCh36'
    opt$genBuild <- 'GRCh38'
    opt$genBuild <- 'GRCh37'
    opt$platform    <- 'NZT'
    opt$version     <- 'N0'
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/NZT/decode')
    opt$ords <- paste(
      file.path(par$aqpDir, 'selected.order1.csv.gz'),
      file.path(par$aqpDir, 'selected.order2.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20297484_probes.match1.tsv.gz'),
      file.path(par$aqpDir, '20297484_probes.match2.tsv.gz'),
      sep=',')
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'BS0031918-AQP1.txt.gz'),
      file.path(par$aqpDir, 'BS0032272-AQP2.txt.gz'),
      sep=',')
    
    opt$pqcs <- NULL
    opt$idat <- NULL
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$parallel <- TRUE
  opt$runName <- paste(opt$platform,opt$version,opt$genBuild, sep='-')
  
  # opt$fresh <- TRUE
  opt$fresh <- FALSE
  
  opt$verbose <- 3
  
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
    make_option(c("--pre_man_csv"), type="character", default=opt$pre_man_csv, 
                help="Previously defined manifest [default= %default]", metavar="character"),

    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    
    make_option(c("--impDir"), type="character", default=opt$impDir, 
                help="improbe data directory [default= %default]", metavar="character"),
    make_option(c("--annDir"), type="character", default=opt$iannDir, 
                help="Annotation data directory [default= %default]", metavar="character"),
    make_option(c("--genDir"), type="character", default=opt$genDir, 
                help="Genomic data directory [default= %default]", metavar="character"),
    make_option(c("--manDir"), type="character", default=opt$manDir, 
                help="Pre-built Manifest data directory [default= %default]", metavar="character"),
    
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
    make_option(c("--opt$genBuild"), type="character", default=opt$genBuild, 
                help="Genome Build (e.g. GRCh36, GRCh37, GRCh38, GRCm38) [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform (e.g. HM450, EPIC, LEGX, NZT, COVIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest Version (e.g. B0,B1,B2,B3,B4,C0) [default= %default]", metavar="character"),
    

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
opt_reqs <- c('outDir','impDir','ords','mats',
              'genBuild','platform','version','bsmap_exe', # 'bsmap_opt',
              'Rscript','verbose')

par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )

opt <- program_init(name=par$prgmTag,
                    opts=opt, opt_reqs=opt_reqs, 
                    pars=par, par_reqs=par_reqs,
                    libs=TRUE,rcpp=FALSE,
                    verbose=opt$verbose,vt=3,tc=0,tt=NULL)

par_tib <- par %>%
  unlist(recursive = TRUE) %>%
  dplyr::bind_rows() %>% 
  tidyr::gather("Params", "Value")
opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

# image_key <- "bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse_Centos"
# image_ver <- "v.1.0"
# image_ssh <- "run_improbe.sh"
# image_str <- glue::glue("{image_key}.{image_ver}")

# Manifest Control Defaults::
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
#
opt$preDir <- file.path(opt$outDir, 'manifest')
if (!dir.exists(opt$preDir)) dir.create(opt$preDir, recursive=TRUE)

opt$intDir <- file.path(opt$outDir, 'intersection')
if (!dir.exists(opt$intDir)) dir.create(opt$intDir, recursive=TRUE)

opt$bspDir <- file.path(opt$outDir, 'bspmap')
if (!dir.exists(opt$bspDir)) dir.create(opt$bspDir, recursive=TRUE)

# Define Alignment Genomes::
#
man_fas_pre <- file.path(opt$outDir, "fas", paste(opt$runName,"aln-prbs", sep="_"))
man_prb_fas <- paste0(man_fas_pre,".fa.gz")
man_gen_fas <- file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
                         paste0(opt$genBuild,".genome.fa.gz"))

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Building Output Directories.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       1.0 Data Collection:: Ord/Mat/AQP/PQC
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD::
#
#  - Create Function for 1.1 - 1.3
#  - Pass in ord,mat,aqp,pqc strings and pars into vectors in function
#
#  - Add 2.0 Data Merging into 1.0 Data Collection

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     1.1.1 Data Collection:: Order Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_ord_tib <- 
  lapply(ords_vec, load_manifestBuildFile, 
         field=par$ord_col[1], cols=par$ord_col,
         verbose=opt$verbose,tt=pTracker) %>% 
  dplyr::bind_rows(.id="ord_aqp")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    1.1.2 Gather (Stack):: Order Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_ord_tib <- 
  format_ORD(man_ord_tib, verbose=opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    1.2.0 Data Collection:: Match Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
add_mat_tib <- 
  lapply(mats_vec, load_manifestBuildFile, 
         field=par$mat_col[1], cols=par$mat_col,
         verbose=opt$verbose,tt=pTracker) %>% 
  dplyr::bind_rows(.id="mat_aqp") %>%
  format_MAT(trim=TRUE, verbose=opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     1.3.1 Data Collection:: AQP Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_aqp_tib <- 
  lapply(aqps_vec, load_manifestBuildFile, 
         field=par$aqp_col[1], cols=par$aqp_col,
         verbose=opt$verbose,tt=pTracker) %>% 
  dplyr::bind_rows(.id="aqp_idx") %>%
  format_AQP(sort=TRUE,filt=TRUE,
             verbose=opt$verbose,tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     1.3.2 Data Collection:: PQC Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_pqc_tib <- 
  lapply(pqcs_vec, load_manifestBuildFile, 
         field=par$pqc_col[1], cols=par$pqc_col,
         verbose=opt$verbose,tt=pTracker) %>% 
  dplyr::bind_rows(.id="aqp_idx") %>%
  format_AQP(sort=TRUE,filt=TRUE,
             verbose=opt$verbose,tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    2.0 Data Merging:: Ord/Mat/AQP/PQC
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     2.2 Join All Stack:: Ord/Mat/AQP/PQC
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# AQC Join::
#
add_aqp_pas_tib <- 
  dplyr::inner_join(
    add_mat_tib, add_aqp_tib,
    by=c("mat_add"="aqp_add")
  ) %>% 
  dplyr::select(mat_add,mat_aqp,mat_seq) %>% 
  dplyr::distinct()

# QC: This should be zero
#
add_aqp_mis_cnt <- 
  add_aqp_pas_tib %>% 
  dplyr::add_count(mat_add, name="Add_Cnt") %>% 
  dplyr::filter(Add_Cnt != 1 ) %>% 
  dplyr::arrange(mat_add) %>%
  base::nrow()
cat(glue::glue("[{par$prgmTag}]: add_aqp_mis_cnt={add_aqp_mis_cnt}{RET}{RET}"))


# PQC Join::
#
add_pqc_pas_tib <- 
  dplyr::inner_join(
    add_mat_tib, add_pqc_tib,
    by=c("mat_add"="aqp_add")
  ) %>% 
  dplyr::select(mat_add,mat_aqp,mat_seq) %>% 
  dplyr::distinct()

# QC: This should be zero
#
add_pqc_mis_cnt <- 
  add_pqc_pas_tib %>% 
  dplyr::add_count(mat_add, name="Add_Cnt") %>% 
  dplyr::filter(Add_Cnt != 1 ) %>% 
  dplyr::arrange(mat_add) %>%
  base::nrow()
cat(glue::glue("[{par$prgmTag}]: add_pqc_mis_cnt={add_pqc_mis_cnt}{RET}{RET}"))


#
# Done::This should be moved into format_ORD()
#
if (FALSE) {

  add_ord_unq_tib <- 
    dplyr::select(add_ord_tib, ord_id,ord_des,ord_aqp,ord_seq) %>% 
    dplyr::arrange(ord_seq, -ord_aqp, ord_des) %>% 
    dplyr::add_count(ord_seq, name="Prb_Cnt") %>%
    dplyr::distinct(ord_seq, .keep_all=TRUE)
  
  add_ord_unq_tib %>% dplyr::filter(Prb_Cnt!=1)
  add_ord_unq_tib %>% dplyr::distinct(ord_seq)
  
}

# Final PQC Address Data::
add_pqc_all_tib <-
  dplyr::inner_join(
    add_ord_tib,add_pqc_pas_tib,
    by=c("ord_seq"="mat_seq")) %>%
  dplyr::arrange(ord_seq) %>%
  dplyr::select(mat_add, mat_aqp, ord_seq,
                ord_id, ord_des, ord_aqp) %>%
  dplyr::rename(prb_add=mat_add,prb_des=ord_des)

# Summary::
add_pqc_all_tib %>%
  dplyr::group_by(mat_aqp,prb_des) %>%
  dplyr::summarise(Count=n(), .groups="drop") %>%
  print(n=100)

# Double Check::
# add_pqc_all_tib <-
#   dplyr::inner_join(
#     add_pqc_pas_tib,add_ord_tib,
#     by=c("mat_seq"="ord_seq")) %>%
#   dplyr::arrange(mat_seq) %>%
#   dplyr::select(mat_add, mat_aqp, mat_seq, 
#                 ord_id, ord_des, ord_aqp)
# 
# add_pqc_all_tib %>% 
#   dplyr::group_by(mat_aqp,ord_des) %>% 
#   dplyr::summarise(Count=n(), .groups="drop") %>%
#   print(n=100)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    3.0 Align All Probe Sequence:: BSMAP
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# aln_key=ord_id %>%
#   stringr::str_replace("BSC_","bsc-") %>%
#   stringr::str_replace("NON_","non-") %>%
#   stringr::str_replace("neg_","neg-") %>%
#   stringr::str_replace_all("\\.","-") %>%
#   stringr::str_replace("_","-") %>%
#   stringr::str_replace_all("_","")
#
# add_pqc_all_tib2 %>% 
#   dplyr::filter(!stringr::str_starts(ord_id, "cg")) %>%
#   dplyr::filter(!stringr::str_starts(ord_id, "ch")) %>%
#   dplyr::filter(!stringr::str_starts(ord_id, "rs")) %>%
#   dplyr::filter(!stringr::str_starts(ord_id, "mu")) %>%
#   dplyr::filter(!stringr::str_starts(ord_id, "rp")) %>%
#   dplyr::filter(!stringr::str_starts(ord_id, "BSC")) %>%
#   dplyr::filter(!stringr::str_starts(ord_id, "neg")) %>%
#   dplyr::filter(!stringr::str_starts(ord_id, "NON"))
#   dplyr::mutate(ord_cgn=stringr::str_remove(ord_id, "_.*$"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       STEP 2:: Alignment:: BSMAAP
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_add_bsp_tsv <- file.path(
  opt$bspDir,"align", paste0(paste(opt$runName,"aln-prbs", sep="_"),"-",
                             opt$genBuild,".genome.bsmap.formatted.tsv.gz"))
man_add_bsp_rds <- file.path(
  opt$bspDir,"align", paste0(paste(opt$runName,"aln-prbs", sep="_"),"-",
                             opt$genBuild,".genome.bsmap.formatted.rds"))

# Add Fasta File::
add_raw_all_tib <- 
  add_pqc_all_tib %>%
  dplyr::mutate(aln_seq=ord_seq %>%
                  stringr::str_replace_all("R","A") %>%
                  stringr::str_replace_all("Y","T"),
                aln_rev=revCmp(aln_seq),
                aln_key=paste(prb_add,prb_des, sep="_")
  )

# Write Fasta File::
man_fas_tib <- 
  tibToFas(tib=add_raw_all_tib, 
           key="aln_key",seq="aln_seq", prefix=man_fas_pre, 
           verbose=opt$verbose, tt=pTracker)

if (!file.exists(man_add_bsp_tsv)) {
  man_add_bsp_tsv <- bsmapProbeAlign(
    exe=opt$bsmap_exe, fas=man_prb_fas, gen=man_gen_fas, 
    dir=opt$bspDir, opt=NULL, run=TRUE,
    # dir=opt$bspDir, opt=opt$bsmap_opt, run=TRUE,
    verbose=opt$verbose,tt=pTracker)
}

if (!file.exists(man_add_bsp_tsv) |
    !file.exists(man_add_bsp_rds) |
    file.mtime(man_add_bsp_rds) < file.mtime(man_add_bsp_tsv) ) {
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading man_add_bsp (TSV)={man_add_bsp_tsv}...{RET}"))
  
  man_add_bsp_tib <- loadBspFormatted(
    bsp=man_add_bsp_tsv, src=dplyr::mutate(man_fas_tib, prb_cgn=as.character(prb_add)), sort=TRUE,
    # bsp=man_add_bsp_tsv, src=man_fas_tib, sort=TRUE,
    verbose=opt$verbose,tt=pTracker)
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Building man_add_bsp (RDS)={man_add_bsp_rds}...{RET}"))
  
  man_add_bsp_grs <-bspToGenomicRegion(
    bsp=man_add_bsp_tib, rds=man_add_bsp_rds,
    verbose=opt$verbose+10,tt=pTracker)
  
} else {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading man_add_bsp (RDS)={man_add_bsp_rds}...{RET}"))
  man_add_bsp_grs <- readr::read_rds(man_add_bsp_rds)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#               4.0 Annotate All Probe Alignments:: CG# Database
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# New Method::
#
#  1. Build cpg_pos_rds (add Top => FR) [ cgn, chr, pos, Top, F/R ]
#  2. cgn/pos <- Intersect(man_add_bsp_grs,cpg_pos_rds)
#  3.0 - Fix Pre-loading issues...
#
#
#  3.1 - NEXT:: Load cpg_top_tsv
#
#
#
#
#  4. Pull all cpgs
#
# opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genBuild,'-21092020.cgnTop.sorted.tsv.gz') )
# opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )

# opt$cph_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.chn-sorted.tsv.gz') )
# opt$snp_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.snp-sorted.tsv.gz') )


imp_top_col <- col <- cols(
  imp_cgn = col_character(),
  imp_top = col_character()
)

imp_pos_rds <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCm10-21092020_improbe-designOutput.cgn-pos-srd.rds"
imp_top_tsv <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/cgnTop/GRCh36-GRCh38-GRCm10-21092020.cgnTop.sorted.tsv.gz"

imp_pos_grs <- readr::read_rds(imp_pos_rds)
imp_top_tib <- readr::read_tsv(imp_top_tsv, 
                               col_names=names(imp_top_col$cols), 
                               col_types=imp_top_col)

# cur_pos_db_grs <- write_impGenomeGRS(genBuild="GRCm10", verbose=10, tt=pTracker)
# cur_pos_db_grs <- write_impGenomeGRS(genBuild="GRCh37", verbose=10, tt=pTracker)
# cur_pos_db_grs <- write_impGenomeGRS(genBuild="GRCh38", verbose=10, tt=pTracker)

#
# NEXT:: Test intersect( man_add_bsp_grs,cur_pos_db_grs )
#
add_imp_all_tib <- 
  intersectGranges(
    man=man_add_bsp_grs, ref=imp_pos_grs,
    verbose=opt$verbose, tt=pTracker)

add_imp_all_tib

#
# Others::
#
#   1. Extract Forward Sequence from BSMAP Alignment...
#
# seqnames   start     end width strand              ord_id  prb_add  prb_cgn prb_des                                        prb_ord_seq
# 20700992_U_1       chr1 3005998 3005999     2      +  cg36602742_F_T_C_I 20700992 20700992       U TTATAAACTTCTCTACAAAACCCAAAACATCACTAACCCTAAATAATTCA
# 28608858_M_1       chr1 3005998 3005999     2      +  cg36602742_F_T_C_I 28608858 28608858       M TTATAAACTTCTCTACAAAACCCAAAACATCACTAACCCTAAATAATTCG
#
#  beg=3005998-3005938; end=3006059-3005999
#  CMD="samtools faidx /Users/bretbarnes/Documents/data/iGenomes/Mus_musculus/NCBI/GRCm10/Sequence/WholeGenomeFasta/GRCm10.genome.fa 1:3005938-3006059"
#  >1:3005938-3006059
#  GATTGGACTTACAATCCATTCAATAACAACAAAAAGTCATGCTTGGTCCTGAAAACCTAA[CG]AACTACCCAGGGCTAGTGATGTCCTGGGTCTTGTAGAGAAGCCCACAACTTTTACTTTAA
#
#   2. Add CH/RS via extraction method above::
#
#   3. Simplify non-improbe design code to run on full cg# database
#


#
# New Method Takes too much memory::
#
if (FALSE) {

  cgn_pos_db_rds  <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCm10-21092020_improbe-designOutput.cgn-map-prb.rds"
  cgn_pos_db_tsv  <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCm10-21092020_improbe-designOutput.cgn-map-prb.tsv.gz"
  
  if (file.exists(cgn_pos_db_tsv)) cgn_pos_db_tib <- readr::read_tsv(cgn_pos_db_tsv)
  # 6:49 to 7:16
  
  if (!file.exists(cgn_pos_db_rds)) {
    
    cgn_pos_db_grs <- 
      GRanges(
        seqnames = Rle(cgn_pos_db_tib$chrChromosome),
        # strand=Rle(stringr::str_sub( bsp$bsp_srd, 1,1 ) ),
        
        imp_cg=cgn_pos_db_tib$Seq_ID,
        imp_fr=cgn_pos_db_tib$Methyl_Allele_FR_Strand,
        imp_tb=cgn_pos_db_tib$Methyl_Allele_TB_Strand,
        imp_co=cgn_pos_db_tib$Methyl_Allele_CO_Strand,
        imp_nb=cgn_pos_db_tib$Methyl_Next_Base,
        
        imp_seq1U=cgn_pos_db_tib$UnMethyl_Probe_Sequence,
        imp_seq1M=cgn_pos_db_tib$Methyl_Probe_Sequence,
        
        IRanges(start=cgn_pos_db_tib$Coordinate,
                end=cgn_pos_db_tib$Coordinate+1,
                names=paste(cgn_pos_db_tib$Seq_ID,
                            paste0(cgn_pos_db_tib$Methyl_Allele_TB_Strand,
                                   cgn_pos_db_tib$Methyl_Allele_CO_Strand),
                            paste0(cgn_pos_db_tib$Methyl_Allele_FR_Strand,
                                   imp_nb=cgn_pos_db_tib$Methyl_Next_Base),
                            cgn_pos_db_tib$chrChromosome,
                            cgn_pos_db_tib$Coordinate,
                            sep="_")
        )
      )
    
  }
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     5.0 Rebuild Manifest:: ord_tib
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #






# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                             6.0 Add Controls:: 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #






# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        7.0 Convert To Genome Studo:: 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #




# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
