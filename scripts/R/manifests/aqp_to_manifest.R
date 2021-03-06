
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

# Null Place Holders::
opt$cpg_top_tsv <- NULL
opt$cpg_pos_tsv <- NULL
opt$cph_pos_tsv <- NULL
opt$snp_pos_tsv <- NULL

# Directories::
opt$outDir  <- NULL
opt$impDir  <- NULL
opt$annDir  <- NULL
opt$genDir  <- NULL
opt$manDir  <- NULL

# Manufacturing Files:: Required
opt$ords <- NULL
opt$mat1 <- NULL
opt$mat2 <- NULL
opt$aqps <- NULL
opt$pqcs <- NULL

# Manufacturing Info:: Required
opt$bpns <- NULL
opt$aqpn <- NULL
opt$pqcn <- NULL

# Pre-defined files (Controls & IDAT Validation):: Optional
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
  # opt$idatDir <- file.path(par$topDir, 'data/idats')
  
  opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
  opt$bsmap_exe <- "/Users/bretbarnes/Documents/programs/BSMAPz/bsmapz"
  
  # Pre-defined local options runTypes::
  #
  par$local_runType <- NULL
  par$local_runType <- 'NZT'
  par$local_runType <- 'COVIC'
  par$local_runType <- 'Chicago'
  par$local_runType <- 'GRCm10'
  
  if (par$local_runType=='Chicago') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'EPIC'
    opt$version  <- 'A1'
    
    opt$idat   <- NULL
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/UnivChicago/latest')
    
    opt$ords <- paste(
      file.path(par$aqpDir, 'UofChicago-A_A_Array-CpG-order-FINAL.csv'),
      sep=',')
    
    opt$mat1 <- NULL
    
    opt$mat2 <- paste(
      file.path(par$aqpDir, '20504790_probes.match.tsv'),
      sep=',')
    
    opt$aqps <- NULL
    
    opt$pqcs <- paste(
      file.path(par$aqpDir, '329922X374054_A_ProductQC.txt'),
      sep=',')
    
  } else if (par$local_runType=='COVIC') {
    opt$genBuild <- 'GRCh36'
    opt$genBuild <- 'GRCh38'
    opt$genBuild <- 'GRCh37'
    
    opt$platform    <- 'COVIC'
    opt$version     <- 'C0'
    opt$version     <- 'B1'
    opt$version     <- 'C2'
    opt$version     <- 'C11'
    opt$version     <- 'C12'
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    opt$pre_man_csv <- '/Users/bretbarnes/Documents/data/manifests/MethylationEPIC_v-1-0_B2.csv.gz'
    
    opt$idat  <- NULL
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/COVID-19_HLA/AQP/COVIC-Host-Immune-Detection')
    opt$ords <- paste(
      file.path(par$aqpDir, 'COVID_EPIC_Round1.03172020.unique.order_AP.csv.gz'),
      sep=',')
    
    opt$mat1 <- paste(
      file.path(par$aqpDir, '20447043_probes.match.gz'),
      sep=',')
    
    opt$mat2 <- NULL
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'BS0032581-AQP.txt.gz'),
      sep=',')
    
    opt$pqcs <- NULL
    # opt$pqcs <- paste(
    #   file.path(par$aqpDir, 'BS0032581-AQP.txt.gz'),
    #   sep=',')
    
    opt$bpns <- paste(1, sep=",")
    opt$aqpn <- paste(1, sep=",")
    opt$pqcn <- NULL
    
  } else if (par$local_runType=='GRCm10') {
    opt$platform <- 'LEGX'
    opt$version  <- 'C0'
    opt$version  <- 'C8'
    opt$version  <- 'C25'
    opt$version  <- 'C26'
    opt$version  <- 'C27'
    opt$version  <- 'C28'
    opt$version  <- 'C29'
    opt$version  <- 'C30'
    
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
    
    opt$mat1 <- paste(
      file.path(par$aqpDir, 'BP1/20420178_AQP1_LifeEpigen_BP1.txt.gz'),
      file.path(par$aqpDir, 'BP2/20420260_AQP1_LifeEpigen_BP2.txt.gz'),
      file.path(par$aqpDir, 'BP3/20420260_AQP2_LifeEpigen_BP2.txt.gz'),
      file.path(par$aqpDir, 'BP4/20455357_AQP1_LifeEpigen_BP4.txt.gz'),
      sep=',')
    
    opt$mat2 <- NULL
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'AQP_Copy/BS0032527-AQP.txt.gz'),
      file.path(par$aqpDir, 'AQP_Copy/BS0032533-AQP.txt.gz'),
      file.path(par$aqpDir, 'AQP_Copy/BS0032545-AQP.txt.gz'),
      file.path(par$aqpDir, 'AQP_Copy/BS0032636-AQP.txt.gz'),
      sep=',')
    
    opt$pqcs <- paste(
      file.path(par$aqpDir, 'PQC/20042400_A_ProductQC.txt.gz'),
      sep=',')
    
    opt$bpns <- paste(1,2,2,3, sep=",")
    opt$aqpn <- paste(1,1,2,1, sep=",")
    opt$pqcn <- paste(1, sep=",")
    
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
    
    opt$mat1 <- paste(
      file.path(par$aqpDir, '20297484_probes.match1.tsv.gz'),
      file.path(par$aqpDir, '20297484_probes.match2.tsv.gz'),
      sep=',')
    
    opt$mat2 <- NULL
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'BS0031918-AQP1.txt.gz'),
      file.path(par$aqpDir, 'BS0032272-AQP2.txt.gz'),
      sep=',')
    
    opt$bpns <- paste(1,2, sep=",")
    opt$aqpn <- paste(1,2, sep=",")
    opt$pqcn <- NULL
    
    opt$pqcs <- NULL
    opt$idat <- NULL
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$parallel <- TRUE
  opt$runName <- paste(par$local_runType,opt$platform,opt$version,opt$genBuild, sep='-')
  
  # opt$fresh <- TRUE
  opt$fresh <- FALSE
  
  opt$verbose <- 10
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
    
    make_option(c("--cpg_pos_tsv"), type="character", default=opt$cpg_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--cpg_top_tsv"), type="character", default=opt$cpg_top_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--cph_pos_tsv"), type="character", default=opt$cph_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--snp_pos_tsv"), type="character", default=opt$snp_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    
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
    
    # Manufacturing Files:: Required
    make_option(c("--ords"), type="character", default=opt$ords, 
                help="Order file(s) (comma seperated) [default= %default]", metavar="character"),
    
    make_option(c("--mat1"), type="character", default=opt$mat1, 
                help="Match (format 1) file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mat2"), type="character", default=opt$mat2, 
                help="Match (format 2) file(s) (comma seperated) [default= %default]", metavar="character"),
    
    make_option(c("--aqps"), type="character", default=opt$aqps, 
                help="AQP file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--pqcs"), type="character", default=opt$pqcs, 
                help="PQC file(s) (comma seperated) [default= %default]", metavar="character"),
    
    # Manufacturing Info:: Required
    make_option(c("--bpns"), type="character", default=opt$bpns, 
                help="Bead Pool Numbers (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqpn"), type="character", default=opt$aqpn, 
                help="AQP Numbers (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--pqcn"), type="character", default=opt$pqcn, 
                help="PQC Numbers; should be 1 (comma seperated) [default= %default]", metavar="character"),
    
    # Pre-defined files (Controls & IDAT Validation):: Optional
    make_option(c("--ctls"), type="character", default=opt$ctls, 
                help="Pre-Defined Infinium Methylation Controls (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--idat"), type="character", default=opt$idat, 
                help="idat directories (comma seperated) [default= %default]", metavar="character"),
    
    # Platform/Method Options::
    make_option(c("--genBuild"), type="character", default=opt$genBuild, 
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
opt_reqs <- c('outDir','impDir','ords',
              'genBuild','platform','version','bsmap_exe', # 'bsmap_opt',
              'Rscript','verbose')

par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )

# print(opt)
# print(opt$genBuild)

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
#                    Pre-processing:: Parse List Options
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
mat1_vec <- NULL
mat2_vec <- NULL
aqps_vec <- NULL
pqcs_vec <- NULL
if (!is.null(opt$ords)) ords_vec <- opt$ords %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$mat1)) mat1_vec <- opt$mat1 %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$mat2)) mat2_vec <- opt$mat2 %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$aqps)) aqps_vec <- opt$aqps %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$pqcs)) pqcs_vec <- opt$pqcs %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

ords_len <- length(ords_vec)
mat1_len <- length(mat1_vec)
mat2_len <- length(mat2_vec)
stopifnot(ords_len>0)
stopifnot(mat1_len>0 || mat2_len>0)
stopifnot(ords_len==mat1_len || ords_len==mat2_len)

bpns_vec <- NULL
aqpn_vec <- NULL
pqcn_vec <- NULL
if (!is.null(opt$bpns)) bpns_vec <- opt$bpns %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$aqpn)) aqpn_vec <- opt$aqpn %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$pqcn)) pqcn_vec <- opt$pqcn %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

# TBD:: Should throw an error if pqcn_vec > 1...

if (is.null(bpns_vec)) bpns_vec <- c(1:ords_len)
if (is.null(aqpn_vec)) aqpn_vec <- c(1:ords_len)
# if (is.null(pqcn_vec)) pqcn_vec <- c(1:ords_len)

man_info_tib <- 
  tibble::tibble(Dat_BPN=bpns_vec, Dat_AQP=aqpn_vec) %>%
  dplyr::mutate(Dat_IDX=dplyr::row_number()) %>% 
  dplyr::select(Dat_IDX, dplyr::everything()) %>% 
  utils::type.convert()

ctls_vec <- NULL
idat_vec <- NULL
if (!is.null(opt$ctls)) ctls_vec <- opt$ctls %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
if (!is.null(opt$idat)) idat_vec <- opt$idat %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Pre-processing:: Run Time:: Output Directories
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run <- NULL

run$manDir <- file.path(opt$outDir, 'man')
if (!dir.exists(run$manDir)) dir.create(run$manDir, recursive=TRUE)

run$addDir <- file.path(opt$outDir, 'add')
if (!dir.exists(run$addDir)) dir.create(run$addDir, recursive=TRUE)

run$intDir <- file.path(opt$outDir, 'int')
if (!dir.exists(run$intDir)) dir.create(run$intDir, recursive=TRUE)

run$fasDir <- file.path(opt$outDir, 'fas')
if (!dir.exists(run$fasDir)) dir.create(run$fasDir, recursive=TRUE)

run$alnDir <- file.path(opt$outDir, 'aln')
if (!dir.exists(run$alnDir)) dir.create(run$alnDir, recursive=TRUE)

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Building Output Directories.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Pre-processing:: Run Time:: Intermediate Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Define Run Time:: Alignment Genome
run$man_gen_fas <- file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
                             paste0(opt$genBuild,".genome.fa.gz"))

# Defined Run Time:: Intermediate Files
run$man_fun_csv <- file.path(run$manDir, paste(opt$runName,"functional-sesame.manifest.csv.gz", sep="_"))

run$add_prb_fas <- file.path(run$fasDir, paste(opt$runName, "aln-seq.fa.gz",  sep='.') )
run$add_dat_csv <- file.path(run$addDir, paste(opt$runName, "add_dat.csv.gz", sep='.') )

run$add_u49_tsv <- file.path(run$intDir, paste(opt$runName, "map-u49.tsv.gz", sep='.') )
run$add_u50_tsv <- file.path(run$intDir, paste(opt$runName, "map-u50.tsv.gz", sep='.') )

run$add_prb_bsp <- file.path(run$alnDir, paste(opt$runName, "bsp",  sep='.') )

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Defining Run Time Files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       1.0 Data Collection:: Ord/Mat/AQP/PQC
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     1.0.0 Data Collection:: IDAT Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_idat_tib <- NULL
if (!is.null(opt$idat)) {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading idats from dir={opt$idat}.{RET}"))
  
  grn_idat_dat <- 
    loadIdat(prefix=file.path(opt$idat, paste(base::basename(opt$idat),"R01C01", sep="_")), 
             col="Grn", gzip=TRUE, verbose=opt$verbose, tt=pTracker )
  
  add_idat_tib <- 
    grn_idat_dat$Quants %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="Address") %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(Address=as.integer(Address)) %>%
    dplyr::distinct(Address)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     1.1.0 Data Collection:: Order Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  is_test_cases <- FALSE
  if (is_test_cases) {
    ver <- 10
    guess_file_del(ords_vec[1], verbose=ver)
    guess_file_del(mat1_vec[1], verbose=ver)
    # guess_file_del(mat2_vec[1], verbose=ver)
    guess_file_del(aqps_vec[1], verbose=ver)
    guess_file_del(pqcs_vec[1], verbose=ver)
    
    
    ver=1
    ord_guess_tib <- guess_man_file(ords_vec[1], verbose=ver)
    mat_guess_tib <- guess_man_file(mat1_vec[1], verbose=ver)
    aqp_guess_tib <- guess_man_file(aqps_vec[1], verbose=ver)
    pqc_guess_tib <- guess_man_file(pqcs_vec[1], verbose=ver)
    
    
    ver=10
    ord_dat_tib <- load_man_file(ords_vec[1], verbose=ver)
    mat_dat_tib <- load_man_file(mat1_vec[1], verbose=ver)
    aqp_dat_tib <- load_man_file(aqps_vec[1], verbose=ver)
    pqc_dat_tib <- load_man_file(pqcs_vec[1], verbose=ver)
    
    ver=3
    idx=1
    ord_dat_tib <- load_man_file(ords_vec[1], idx=idx, verbose=ver)
    mat_dat_tib <- load_man_file(mat1_vec[1], idx=idx, verbose=ver)
    aqp_dat_tib <- load_man_file(aqps_vec[1], idx=idx, verbose=ver)
    pqc_dat_tib <- load_man_file(pqcs_vec[1], idx=idx, verbose=ver)
  }
  
  #
  # [Done]: TBD:: Add Formating for::
  #   mat: uc(Sequence), Fix(Address)
  #   aqs: Fix(Address)
  #
  # [Half]: TBD:: Add skipping null inputs (i.e. mat1/mat2)
  #   Solution:: Try passing in both mat1/mat2 as mats now...
  #
  # [Done]: TBD:: Pipe in vectors to make it look more seemeless...
  #
  # TBD:: Pass in all pre-defined cols (val, sel, key)
  #
  
  # This method requires removing idx option
  #  - Probably better option for simplicity...
  ords_dat_tib <- ords_vec %>%
    lapply(load_man_file,
           verbose=opt$verbose,tt=pTracker) %>% 
    dplyr::bind_rows(.id="Dat_IDX") %>%
    utils::type.convert() %>%
    dplyr::mutate(across(where(is.factor),  as.character) ) %>%
    dplyr::left_join(man_info_tib, by="Dat_IDX")
  
  # Latest Order Summary::
  ords_sum_tib <- ords_dat_tib %>% 
    dplyr::group_by(Dat_IDX,Dat_BPN,Dat_AQP,Dat_PQC, Ord_Din,Ord_Des,Ord_Col) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  ords_sum_tib %>% print(n=base::nrow(ords_sum_tib))
  
  #
  # TBD::: Build general summaries with file and screen outputs...
  #
  mat1_dat_tib <- mat1_vec %>%
    lapply(load_man_file,
           verbose=opt$verbose,tt=pTracker) %>% 
    dplyr::bind_rows(.id="Dat_IDX") %>%
    utils::type.convert() %>%
    dplyr::mutate(across(where(is.factor),  as.character) ) %>%
    dplyr::left_join(man_info_tib, by="Dat_IDX")

  mat2_dat_tib <- mat2_vec %>%
    lapply(load_man_file,
           verbose=opt$verbose,tt=pTracker) %>% 
    dplyr::bind_rows(.id="Dat_IDX") %>%
    utils::type.convert() %>%
    dplyr::mutate(across(where(is.factor),  as.character) ) %>%
    dplyr::left_join(man_info_tib, by="Dat_IDX")
  
  aqps_dat_tib <- aqps_vec %>%
    lapply(load_man_file,
           verbose=opt$verbose,tt=pTracker) %>% 
    dplyr::bind_rows(.id="Dat_IDX") %>%
    utils::type.convert() %>%
    dplyr::mutate(across(where(is.factor),  as.character) ) %>%
    dplyr::left_join(man_info_tib, by="Dat_IDX")
  
  #
  # NOTE:: There's only one PQC, so no need to join with man_info_tib
  #
  pqcs_dat_tib <- pqcs_vec %>%
    lapply(load_man_file,
           verbose=opt$verbose,tt=pTracker) %>% 
    dplyr::bind_rows(.id="Dat_PQC") %>%
    utils::type.convert() %>%
    dplyr::mutate(across(where(is.factor),  as.character) )

  
  
  
  
  
  
    
  # TBD:: Slower, but possibily better method
  #   - Maybe, maybe look into this later...
  # purrr::imap(ords_vec, function(x,y) { 
  #   load_man_file(x,y, verbose=opt$verbose,tt=pTracker) } ) %>%
  #   dplyr::bind_rows()
  
}


# man_ord_tib <- NULL
# man_ord_tib <-
#   lapply(ords_vec, load_manifestBuildFile,
#          field=par$ord_col[1], cols=par$ord_col,
#          verbose=opt$verbose,tt=pTracker) %>%
#   dplyr::bind_rows(.id="IDX")
# 
# We need the mappings::
#  man_ord_tib %>% dplyr::filter(!is.na(AlleleB_Probe_Id))
#

add_ord_tib <- NULL
add_ord_tib <- 
  lapply(ords_vec, load_manifestBuildFile, 
         field=par$ord_col[1], cols=par$ord_col,
         verbose=opt$verbose,tt=pTracker) %>% 
  dplyr::bind_rows(.id="IDX") %>%
  format_ORD(idx_key="IDX", uniq=TRUE,
             verbose=opt$verbose, tt=pTracker) %>%
  dplyr::left_join(man_info_tib, by="IDX") %>% 
  dplyr::mutate(prb_type=stringr::str_sub(ord_id, 1,2)) %>%
  dplyr::select( dplyr::any_of(
    c( base::names(man_info_tib),
       "ord_id","prb_type","prb_des","prb_col","prb_seq","prb_par" )) ) %>%
  dplyr::distinct()

add_ord_sum <- add_ord_tib %>%
  dplyr::group_by( dplyr::across( 
    c(base::names(man_info_tib),"prb_type","prb_des","prb_col") ) ) %>% 
  dplyr::summarise(Count=n(), .groups="drop")
add_ord_sum %>% print(n=base::nrow(add_ord_sum))

# all_ord_tib <- NULL
# all_ord_tib <- 
#   lapply(ords_vec, load_manifestBuildFile, 
#          field=par$ord_col[1], cols=par$ord_col,
#          verbose=opt$verbose,tt=pTracker) %>% 
#   dplyr::bind_rows(.id="IDX") %>%
#   format_ORD(idx_key="IDX", uniq=FALSE,
#              verbose=opt$verbose, tt=pTracker) %>%
#   dplyr::left_join(man_info_tib, by="IDX")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    1.2.0 Data Collection:: Match Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_mat_tib <- NULL
if (!is.null(mat1_vec)) {
  
  add_mat_tib <- 
    lapply(mat1_vec, load_manifestBuildFile, 
           field=par$mat_col[1], cols=par$mat_col,
           verbose=opt$verbose,tt=pTracker) %>% 
    dplyr::bind_rows(.id="IDX") %>%
    format_MAT(idx_key="IDX",uniq=TRUE,trim=TRUE, 
               verbose=opt$verbose, tt=pTracker) %>%
    dplyr::left_join(man_info_tib, by="IDX") %>% 
    dplyr::select( dplyr::any_of( 
      c( base::names(man_info_tib), "Address","prb_seq","tan_seq" )) ) %>%
    dplyr::distinct()
  
  add_mat_sum <- add_mat_tib %>%
    dplyr::group_by( dplyr::across( base::names(man_info_tib) )) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  add_mat_sum %>% print(n=base::nrow(add_mat_sum))
  
  # all_mat_tib <- NULL
  # all_mat_tib <- 
  #   lapply(mat1_vec, load_manifestBuildFile, 
  #          field=par$mat_col[1], cols=par$mat_col,
  #          verbose=opt$verbose,tt=pTracker) %>% 
  #   dplyr::bind_rows(.id="IDX") %>%
  #   format_MAT(idx_key="IDX",uniq=FALSE,trim=TRUE, 
  #              verbose=opt$verbose, tt=pTracker) %>%
  #   dplyr::left_join(man_info_tib, by="IDX")
  
} else if (!is.null(mat2_vec)) {
  
  add_mat_tib <- 
    lapply(mat2_vec, load_manifestBuildFile, 
           field=par$ma2_col[1], cols=par$ma2_col,
           verbose=opt$verbose,tt=pTracker) %>% 
    dplyr::bind_rows(.id="IDX") %>%
    format_MAT(idx_key="IDX",uniq=TRUE,trim=TRUE, 
               verbose=opt$verbose, tt=pTracker) %>%
    dplyr::left_join(man_info_tib, by="IDX") %>% 
    dplyr::select( dplyr::any_of( 
      c( base::names(man_info_tib), "Address","prb_seq","tan_seq" )) ) %>%
    dplyr::distinct()
  
  add_mat_sum <- add_mat_tib %>%
    dplyr::group_by( dplyr::across( base::names(man_info_tib) )) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  add_mat_sum %>% print(n=base::nrow(add_mat_sum))
  
} else {
  stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Neither mats or mat2 defined!!!{RET}{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             1.2.1 Data Collection:: Remove Missing IDAT Tangos
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (!is.null(add_idat_tib)) {
  if (opt$verbose>=1) {
    cat(glue::glue("[{par$prgmTag}]: Removing tangos missing in IDATs; Beg={RET}"))
    print(add_mat_tib)
  }
  
  add_mat_tib <- 
    dplyr::filter(Address %in% add_idat_tib$Address)
  
  if (opt$verbose>=1) {
    cat(glue::glue("[{par$prgmTag}]: Removing tangos missing in IDATs; End={RET}"))
    print(add_mat_tib)
    cat("\n")
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     1.3.1 Data Collection:: AQP Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_aqp_tib <- NULL
if (!is.null(aqps_vec)) {
  add_aqp_tib <- 
    lapply(aqps_vec, load_manifestBuildFile, 
           field=par$aqp_col[1], cols=par$aqp_col,
           verbose=opt$verbose,tt=pTracker) %>% 
    dplyr::bind_rows(.id="IDX") %>%
    format_AQP(idx_key="IDX",uniq=TRUE,filt=TRUE,
               verbose=opt$verbose,tt=pTracker) %>%
    dplyr::left_join(man_info_tib, by="IDX")
  
  add_aqp_sum <- add_aqp_tib %>% 
    dplyr::group_by( dplyr::across( base::names(man_info_tib) )) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  add_aqp_sum %>% print(n=base::nrow(add_aqp_sum))
  
  # all_aqp_tib <- NULL
  # all_aqp_tib <- 
  #   lapply(aqps_vec, load_manifestBuildFile, 
  #          field=par$aqp_col[1], cols=par$aqp_col,
  #          verbose=opt$verbose,tt=pTracker) %>% 
  #   dplyr::bind_rows(.id="IDX") %>%
  #   format_AQP(idx_key="IDX",uniq=FALSE,filt=FALSE,
  #              verbose=opt$verbose,tt=pTracker) %>%
  #   dplyr::left_join(man_info_tib, by="IDX")
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     1.3.2 Data Collection:: PQC Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_pqc_tib <- NULL
if (!is.null(pqcs_vec)) {
  add_pqc_tib <- 
    lapply(pqcs_vec, load_manifestBuildFile, 
           field=par$pqc_col[1], cols=par$pqc_col,
           verbose=opt$verbose,tt=pTracker) %>% 
    dplyr::bind_rows(.id="PQC") %>%
    format_AQP(idx_key="PQC",uniq=TRUE,filt=TRUE,
               verbose=opt$verbose,tt=pTracker) %>%
    dplyr::distinct(Address)
  
  # add_pqc_sum <- add_pqc_tib %>%
  #   dplyr::group_by(PQC) %>%
  #   dplyr::summarise(Count=n(), .groups="drop")
  # add_pqc_sum %>% print(n=base::nrow(add_pqc_sum))
  # 
  # all_pqc_tib <- NULL
  # all_pqc_tib <- 
  #   lapply(pqcs_vec, load_manifestBuildFile, 
  #          field=par$pqc_col[1], cols=par$pqc_col,
  #          verbose=opt$verbose,tt=pTracker) %>% 
  #   dplyr::bind_rows(.id="PQC") %>%
  #   format_AQP(idx_key="PQC",uniq=FALSE,filt=FALSE,
  #              verbose=opt$verbose,tt=pTracker)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    2.0 Data Merging:: Ord/Mat/AQP/PQC
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_pas_tib <- NULL
if (!is.null(add_pqc_tib)) {
  add_pas_pqc_mat_tib <- 
    add_mat_tib %>% 
    dplyr::filter(Address %in% add_pqc_tib$Address) %>%
    dplyr::inner_join(add_ord_tib, by=c("prb_seq","IDX","BPN","AQP","PQC") )
  
  add_pas_pqc_mat_sum <- 
    add_pas_pqc_mat_tib %>%
    dplyr::group_by( dplyr::across( c(prb_type,prb_des,prb_col,base::names(man_info_tib)) )) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  add_pas_pqc_mat_sum %>% 
    print(n=base::nrow(add_pas_pqc_mat_sum))
  
  # QC: This should be zero
  add_pas_pqc_mat_cnt <- 
    add_pas_pqc_mat_tib %>% 
    dplyr::add_count(Address, name="Add_Cnt") %>% 
    dplyr::filter(Add_Cnt != 1 ) %>% 
    base::nrow()
  cat(glue::glue("[{par$prgmTag}]: add_pas_pqc_mat_cnt={add_pas_pqc_mat_cnt}{RET}{RET}"))
  
  add_pas_tib <- add_pas_pqc_mat_tib
} else if (!is.null(add_aqp_tib)) {
  add_pas_aqp_mat_tib <- 
    add_mat_tib %>% 
    dplyr::filter(Address %in% add_aqp_tib$Address) %>%
    dplyr::inner_join(add_ord_tib, by=c("prb_seq","IDX","BPN","AQP","PQC") )
  
  add_pas_aqp_mat_sum <- 
    add_pas_aqp_mat_tib %>%
    dplyr::group_by( dplyr::across( c(prb_type,prb_des,base::names(man_info_tib)) )) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  add_pas_aqp_mat_sum %>% 
    print(n=base::nrow(add_pas_aqp_mat_sum))
  
  # QC: This should be zero
  add_pas_aqp_mat_cnt <- 
    add_pas_aqp_mat_tib %>% 
    dplyr::add_count(Address, name="Add_Cnt") %>% 
    dplyr::filter(Add_Cnt != 1 ) %>% 
    base::nrow()
  cat(glue::glue("[{par$prgmTag}]: add_pas_aqp_mat_cnt={add_pas_aqp_mat_cnt}{RET}{RET}"))
  
  add_pas_tib <- add_pas_aqp_mat_tib
} else {
  stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: Both AQP and PQC are length zero!!!{RET}{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  2.1 Functional Manifest Generation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: 
#  - Clear up extra validation
#  - Write Functional Manifest with ord_id's
#  - Write Function 
#     addressToManifest(add=add_pas_tib, outDir=run$manDir, name=paste(opt$runName,"functional", sep="_"), validate=TRUE, verbose=opt$verbose, tt=pTracker)
#
man_fun_vec <- c("ord_id","prb_type","prb_col", "prb_par",
                 dplyr::all_of(base::names(man_info_tib) ) )

man_fun_tib <- add_pas_tib %>%
  addressToManifest(des="prb_des", join=man_fun_vec,
                    csv=run$man_fun_csv, 
                    validate=TRUE, 
                    verbose=opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    3.0 Align All Probe Sequence:: BSMAP
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_pas_fas_tib <-
  tib_to_fasta(
    tib=add_pas_tib, 
    prb_key="prb_seq",
    add_key="Address", des_key="prb_des", type_key="prb_type",
    prb_fas=run$add_prb_fas, dat_csv=run$add_dat_csv,
    u49_tsv=run$add_u49_tsv, u50_tsv=run$add_u50_tsv,
    verbose=opt$verbose+10, tt=pTracker)

# total alignments = 109367
run$add_prb_bspz <- 
  run_bsmap(
    exe=opt$bsmap_exe, 
    fas=run$add_prb_fas, 
    gen=run$man_gen_fas,
    bsp=run$add_prb_bsp,
    # dir=run$alnDir, 
    opt=NULL, lan=NULL, run=TRUE,
    verbose=opt$verbose,tt=pTracker)










add_pas_bsp_tib <- NULL
add_pas_bsp_grs <- NULL
# if (!file.exists(man_add_bsp_tsv) ||
#     !file.exists(man_add_bsp_rds) ||
#     file.mtime(man_add_bsp_rds) < file.mtime(man_add_bsp_tsv) ) {
if (TRUE) {
  
  add_pas_bsp_tibI <-
    join_bsmap(
      add=add_pas_fas_tib,
      file=run$add_prb_bspz,
      join_key="aln_key",
      sort=TRUE,
      verbose=opt$verbose+10,tt=pTracker)
  
  add_pas_bsp_tibL <-
    join_bsmap(
      add=add_pas_fas_tib,
      file=run$add_prb_bspz,
      join_key="aln_key", join_type="left",
      sort=TRUE,
      verbose=opt$verbose+10,tt=pTracker)

  add_pas_bsp_tibR <-
    join_bsmap(
      add=add_pas_fas_tib,
      file=run$add_prb_bspz,
      join_key="aln_key", join_type="right",
      sort=TRUE,
      verbose=opt$verbose+10,tt=pTracker)

  add_pas_bsp_tibF <-
    join_bsmap(
      add=add_pas_fas_tib,
      file=run$add_prb_bspz,
      join_key="aln_key", join_type="full",
      sort=TRUE,
      verbose=opt$verbose+10,tt=pTracker)
  
  # Structure for grs::
  # add_pas_bsp_grs
  
  
  # Example of two stage: Load -> join::
  #
  # bsp_tib <- 
  #   load_bsmap(
  #     bsp=run$add_prb_bspz,
  #     sort=TRUE,
  #     verbose=opt$verbose+10,tt=pTracker)
  # 
  # add_pas_bsp_tib <- 
  #   join_bsmap(
  #     add=add_pas_fas_tib, 
  #     bsp=bsp_tib,
  #     join_key="aln_key",
  #     sort=TRUE,
  #     verbose=opt$verbose+10,tt=pTracker)
}











#
# Old Code Below::
#
if (FALSE) {

  if (!file.exists(man_add_bsp_tsv) ||
      !file.exists(man_add_bsp_rds) ||
      file.mtime(man_add_bsp_rds) < file.mtime(man_add_bsp_tsv) ) {
    
    if (FALSE) {
      man_add_bsp_tib <- NULL
      man_add_bsp_grs <- NULL
      
      man_add_bsp_tib <- loadBspFormatted(
        bsp=man_add_bsp_tsv, src=man_fas_tib, sort=TRUE,
        # bsp=man_add_bsp_tsv, src=dplyr::mutate(man_fas_tib, prb_cgn=as.character(prb_add)), sort=TRUE,
        # bsp=man_add_bsp_tsv, src=man_fas_tib, sort=TRUE,
        verbose=opt$verbose,tt=pTracker)
      
      if (opt$verbose>=1)
        cat(glue::glue("[{par$prgmTag}]: Building man_add_bsp (RDS)={man_add_bsp_rds}...{RET}"))
      
      man_add_bsp_grs <- bspToGenomicRegion(
        bsp=man_add_bsp_tib, rds=man_add_bsp_rds,
        verbose=opt$verbose,tt=pTracker)
      
    }
    
  } else {
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Loading man_add_bsp (RDS)={man_add_bsp_rds}...{RET}"))
    man_add_bsp_grs <- readr::read_rds(man_add_bsp_rds)
  }
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       STEP 3:: Annotate:: by 50U/49U
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_cgn_mat <- TRUE
if (run_cgn_mat) {
  
  imp_prb49U_tsv <- "/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/prbDbs/GRCh36-38.mm10.cgn-srd-seq49U.prb-sorted.tsv.gz"
  imp_prb50U_tsv <- "/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/prbDbs/GRCh36-38.mm10.cgn-srd-seq50U.prb-sorted.tsv.gz"
  
  add_prb49U_tsv <- file.path(opt$outDir,"intersect", paste(opt$runName,"add-prb49U.seq-sorted.tsv", sep="_"))
  add_prb50U_tsv <- file.path(opt$outDir,"intersect", paste(opt$runName,"add-prb50U.seq-sorted.tsv", sep="_"))
  
  int_prb49U_tsv <- file.path(opt$outDir,"intersect", paste(opt$runName,"imp-add-prb49U.seq-intersect.tsv.gz", sep="_"))
  int_prb50U_tsv <- file.path(opt$outDir,"intersect", paste(opt$runName,"imp-add-prb50U.seq-intersect.tsv.gz", sep="_"))
  
  #
  #
  # Big Change:: 
  #  - add_pas_tib should replace man_add_bsp_tib below::
  #  - Rational need to use all passing sequences not just the ones that aligned to the genome of interest...
  #
  #  Tasks::
  #   - Add prb49U/prb50U need to be defined well above (probably at fasta file generation:: man_fas_tib)
  #   - This will take a little restructuring...
  #
  #  Convert 
  #
  
  # Previous Conversion of Genomic Regions Alignment::
  use_grs_mat <- FALSE
  if (use_grs_mat) {
    # Address BSP tibble instead of GRS::
    man_add_bsp_tib <- 
      man_add_bsp_grs %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var="unq_add_id") %>% 
      tibble::as_tibble() %>%
      dplyr::arrange(prb_aln_50U) %>%
      dplyr::distinct(prb_aln_50U, .keep_all=TRUE)
  }
  
  #
  # TBD:: These two steps below need to be performed during the fasta file generation
  #
  
  # Build/Write prb49U
  add_prb49U_tib <- 
    man_add_bsp_tib %>%
    dplyr::filter(prb_des=="2") %>%
    dplyr::mutate(
      add_prb49U=stringr::str_sub(prb_aln_50U, 2)
    ) %>%
    dplyr::select(unq_add_id,add_prb49U)
  
  if (!file.exists(imp_prb49U_tsv) ||
      !file.exists(add_prb49U_tsv) ||
      file.mtime(imp_prb49U_tsv) > file.mtime(add_prb49U_tsv) ) {
    
    readr::write_tsv(
      add_prb49U_tib, add_prb49U_tsv, col_names=FALSE)
  }
  
  # Build/Write prb50U
  add_prb50U_tib <- 
    man_add_bsp_tib %>%
    dplyr::filter(prb_des!="2") %>%
    dplyr::select(unq_add_id,prb_aln_50U)
  
  if (!file.exists(imp_prb50U_tsv) ||
      !file.exists(add_prb50U_tsv) ||
      file.mtime(imp_prb50U_tsv) > file.mtime(add_prb50U_tsv) ) {
    
    readr::write_tsv(
      add_prb50U_tib, add_prb50U_tsv, col_names=FALSE)
  }
  
  # TBD:: 
  #  Candidate::
  #  - Don't need to add ord_id,prb_des
  #  - Split files by aln_prb50U,aln_prb49U
  #
  #  Reference::
  #  - Combine cgn_int,srd_int
  #  - Split files by aln_prb50U,aln_prb49U
  #    i.e. split this file: "/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/GRCh36-21022021_improbe-designOutput.cgn-seq50U.unq-cgn-srd.seq50U-sorted.short-fields.tsv.gz"
  #    [done]::
  #
  # opt$fresh <- TRUE
  
  int_seq_col <- cols(
    mat_seq = col_character(),
    mat_cgn = col_integer(),
    mat_sr2 = col_integer(),
    mat_can = col_character()
  )
  
  if (opt$fresh ||
      !file.exists(imp_prb49U_tsv) ||
      !file.exists(add_prb49U_tsv) ||
      !file.exists(int_prb49U_tsv) ||
      file.mtime(imp_prb49U_tsv) > file.mtime(add_prb49U_tsv) ||
      file.mtime(add_prb49U_tsv) > file.mtime(int_prb49U_tsv) ) {
    
    int_49U_cmd = glue::glue("gzip -dc {imp_prb49U_tsv} | join -t $'\t' -13 -22 - {add_prb49U_tsv} | gzip -c - > {int_prb49U_tsv}")
    int_49U_ret <- system(int_49U_cmd)
  }
  
  int_prb49U_tib <- 
    readr::read_tsv(int_prb49U_tsv,
                    col_names=names(int_seq_col$cols),
                    col_types=int_seq_col)
  
  # Split into compareable table::
  int_mat49U_tib <- int_prb49U_tib %>% 
    tidyr::separate(mat_can, into=c("can_add","can_des","can_rep"), sep="_") %>% 
    utils::type.convert() %>% 
    dplyr::mutate(across(where(is.factor),  as.character) )
  
  
  if (opt$fresh ||
      !file.exists(imp_prb50U_tsv) ||
      !file.exists(add_prb50U_tsv) ||
      !file.exists(int_prb50U_tsv) ||
      file.mtime(imp_prb50U_tsv) > file.mtime(add_prb50U_tsv) ||
      file.mtime(add_prb50U_tsv) > file.mtime(int_prb50U_tsv) ) {
    
    int_50U_cmd = glue::glue("gzip -dc {imp_prb50U_tsv} | join -t $'\t' -13 -22 - {add_prb50U_tsv} | gzip -c - > {int_prb50U_tsv}")
    int_50U_ret <- system(int_50U_cmd)
  }
  
  int_prb50U_tib <- 
    readr::read_tsv(int_prb50U_tsv,
                    col_names=names(int_seq_col$cols),
                    col_types=int_seq_col)
  
  # Split into comparable table::
  int_mat50U_tib <- int_prb50U_tib %>% 
    tidyr::separate(mat_can, into=c("can_add","can_des","can_rep"), sep="_") %>% 
    utils::type.convert() %>% 
    dplyr::mutate(across(where(is.factor),  as.character) )
  
  
  # Summary::
  int_mat50U_tib %>% 
    # dplyr::add_count(can_add, name="Comb_Cnt") %>% 
    dplyr::add_count(mat_cgn, name="Comb_Cnt") %>% 
    # dplyr::add_count(mat_cgn,can_add, name="Comb_Cnt") %>% 
    dplyr::group_by(Comb_Cnt) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  
}




# Add Extra Infinum II Column::
#
# cmd="gzip -dc /Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/GRCh36-21022021_improbe-designOutput.cgn-seq50U.unq-cgn-srd.seq50U-sorted.short-fields.tsv.gz | perl -pe 's/\n$//; @d=split("\t",$_); s/^.*$//; print "$d[0]\t$d[1]\t$d[2]\t".substr($d[2],0,length($d[2])-1)."\n"; ' | gzip -c -> /Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/GRCh36-21022021_improbe-designOutput.cgn-seq50U.unq-cgn-srd.seq50U-sorted.both-fields.tsv.gz"
# Script=/Users/bretbarnes/Documents/tools/shells/improbe/sort-prb50U-all.sh
#

#
# Run Basic Probe seq50U intersection::
#  - This needs to be cleaned up for seq49U Infinium II matching...
# cmd="gzip -dc /Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/GRCh36-21022021_improbe-designOutput.cgn-seq50U.unq-cgn-srd.seq50U-sorted.short-fields.tsv.gz | join -t $'\t' -13 -24 - /Users/bretbarnes/Documents/scratch/aqp_to_manifest/LEGX-C25-GRCm10/manifest/LEGX-C25-GRCm10_aln-prbs-GRCm10.UniqueID-AlnSeq.tsv > /Users/bretbarnes/Documents/scratch/aqp_to_manifest/LEGX-C25-GRCm10/intersection/intersect-seq50U.tsv"
#
# cmd="gzip -dc /Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/GRCh36-21022021_improbe-designOutput.cgn-seq50U.unq-cgn-srd.seq50U-sorted.both-fields.tsv.gz | join -t $'\t' -13 -24 - /Users/bretbarnes/Documents/scratch/aqp_to_manifest/LEGX-C25-GRCm10/manifest/LEGX-C25-GRCm10_aln-prbs-GRCm10.UniqueID-AlnSeq.tsv > /Users/bretbarnes/Documents/scratch/aqp_to_manifest/LEGX-C25-GRCm10/intersection/intersect-seq50U.tsv"
# cmd="gzip -dc /Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/GRCh36-21022021_improbe-designOutput.cgn-seq50U.unq-cgn-srd.seq50U-sorted.both-fields.tsv.gz | join -t $'\t' -14 -25 - /Users/bretbarnes/Documents/scratch/aqp_to_manifest/LEGX-C25-GRCm10/manifest/LEGX-C25-GRCm10_aln-prbs-GRCm10.UniqueID-AlnSeq.tsv > /Users/bretbarnes/Documents/scratch/aqp_to_manifest/LEGX-C25-GRCm10/intersection/intersect-seq49U.tsv"

# Run Intersection::

# man_add_prb_tsv2   <- "/Users/bretbarnes/Documents/scratch/aqp_to_manifest/LEGX-C25-GRCm10/manifest/LEGX-C25-GRCm10_aln-prbs-GRCm10.UniqueID-AlnSeq.tsv"

# cgn_ann_seq50U_tsv <- "/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/GRCh36-21022021_improbe-designOutput.cgn-seq50U.unq-cgn-srd.seq50U-sorted.both-fields.tsv.gz"
# cmd_49U=glue::glue("gzip -dc {cgn_ann_seq50U_tsv} | join -t $'\t' -14 -25 - {man_add_prb_tsv} | cut -f 1,2,3,5,6,7,8 | gzip -c - > {man_add_int49U_tsv}")
# cmd_50U=glue::glue("gzip -dc {cgn_ann_seq50U_tsv} | join -t $'\t' -13 -24 - {man_add_prb_tsv} | cut -f 1,2,3,5,6,7,8 | gzip -c - > {man_add_int50U_tsv}")

if (FALSE) {
  
  if (FALSE) {
    #
    #  NOTE: This is just a quick fix for the Horvath Sesame Array::
    #    - This should be moved to a seperate script. 
    #
    man_hov_csv <- "/Users/bretbarnes/Documents/data/manifests/HorvathMammal40-A2.sesame-base.cpg-sorted.csv.gz"
    man_hov_tib <- 
      readr::read_csv( file.path(par$datDir, "manifest/core/HorvathMammal40-A1.sesame-base.cpg-sorted.csv.gz") ) %>% 
      dplyr::mutate(
        DESIGN_COLOR=dplyr::case_when(
          col=="R" ~ "Red", 
          col=="G" ~ "Grn", 
          TRUE ~ NA_character_),
        Probe_Type=stringr::str_sub(Probe_ID, 1,2),
        Next_Base=dplyr::case_when(
          col=="R" ~ "T",
          col=="G" ~ "C",
          TRUE ~ NA_character_
        ),
        DESIGN=dplyr::case_when(
          !is.na(U) & !is.na(M) ~ "I",
          !is.na(U) &  is.na(M) ~ "II",
          TRUE ~ NA_character_
        ),
        Probe_Design=dplyr::case_when(
          !is.na(U) & !is.na(M) ~ 1,
          !is.na(U) &  is.na(M) ~ 2,
          TRUE ~ NA_real_
        ),
        Probe_Source="HorvathMammal40"
      )
    readr::write_csv(man_hov_tib, man_hov_csv)
    
    
    # DESIGN_COLOR,col,Probe_Type,Probe_Source,Next_Base,Probe_Design
  }
  
  # Load First Round:: II
  
  int_seq49U_col <- cols(
    mat_49U = col_character(),
    
    imp_cgn = col_integer(),
    imp_srd = col_integer(),
    
    # mat_50U = col_character(),
    # mat_49U = col_character(),
    
    unq_id    = col_character(),
    ord_id    = col_character(),
    prb_des   = col_character(),
    aln_prb50U = col_character()
  )
  
  int_seq49U_tsv <- "/Users/bretbarnes/Documents/scratch/aqp_to_manifest/LEGX-C26-GRCm10/intersection/intersect-seq49U.tsv.gz"
  int_seq49U_tib <- readr::read_tsv(int_seq49U_tsv,
                                    col_names=names(int_seq49U_col$cols),
                                    col_types=int_seq49U_col)
  
  # Load Second Round:: I
  int_seq50U_tsv <- "/Users/bretbarnes/Documents/scratch/aqp_to_manifest/LEGX-C25-GRCm10/intersection/intersect-seq50U.tsv.gz"
  int_seq50U_tib <- readr::read_tsv(int_seq50U_tsv)
  
  # Build Address Data (add = prb)
  # Align 
  #  - Need to improve alignment summary
  # Annotate
  #  - Load seq50U
  #  - Parse imp48U
  #  - Parse prb48U
  #  - Join by imp48U/prb48
  #  - Add Canonical cg# flag
  # Redesign
  #  - Extract 122mer (ref/snp)
  #  - Redesign
  #  - SNP affect summary
  
  imp_cgn_prb_tsv <- "/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/GRCh36-21022021_improbe-designOutput.cgn-seq50U.unq-cgn-srd-sorted.short-fields.tsv.gz"
  imp_cgn_prb_rds <- "/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/GRCh36-21022021_improbe-designOutput.cgn-seq50U.unq-cgn-srd-sorted.short-fields.rds"
  
  if (file.exists(imp_cgn_prb_rds)) {
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Loading imp_cgn_prb_rds (RDS)={imp_cgn_prb_rds}...{RET}"))
    
    imp_cgn_prb_tib <- readr::read_rds(imp_cgn_prb_rds)
  } else {
    if (file.exists(imp_cgn_prb_tsv)) {
      if (opt$verbose>=1)
        cat(glue::glue("[{par$prgmTag}]: Loading imp_cgn_prb_tsv (TSV)={imp_cgn_prb_tsv}...{RET}"))
      
      imp_time <- system.time({
        imp_cgn_prb_col <- cols(
          imp_cgn = col_integer(),
          imp_srd = col_integer(),
          imp_seq = col_character()
        )
        
        imp_cgn_prb_tib <- 
          readr::read_tsv(imp_cgn_prb_tsv,
                          col_names=names(imp_cgn_prb_col$cols),
                          col_types=imp_cgn_prb_col)
        
        if (opt$verbose>=1)
          cat(glue::glue("[{par$prgmTag}]: Writing imp_cgn_prb_rds (RDS)={imp_cgn_prb_rds}...{RET}"))
        
        # readr::write_rds(imp_cgn_prb_tib, imp_cgn_prb_rds, compress="gz")
      })
    } else {
      cat(glue::glue("ERROR: Failed to find imp_cgn_prb_tsv={imp_cgn_prb_tsv}!!!{RET}{RET}"))
    }
  }
  
  add_bsp_imp_tib <- man_add_bsp_tib2 %>% 
    dplyr::left_join(imp_cgn_prb_tib, by=c("aln_prb"="imp_seq"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #               4.0 Annotate All Probe Alignments:: CG# Database
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  gen_map_tsv <- "/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designInput/GRCm10.cgn-pos-fwd.base.tsv.gz"
  gen_map_tib <- write_impTopGenomeGRS(genBuild = opt$genBuild, verbose = 10)
  
  # Load genome for sub-string of design sequence::
  #
  tmp_fas_dat <- Biostrings::readDNAStringSet(filepath = run$man_gen_fas, format = "fasta", nrec = 2)
  
  #
  # Tasks::
  #  - Extract (sub-string) forward sequence from NCBI Reference Genome
  #  - Convert designOutput.cgn-pos-fwd.tsv.gz to GR Range
  #    Unique_ID=[cg,srd,chr,pos], chr=paste("chr",chr), top_seq, start=pos, width=2
  #  - Validate forward sequence
  #  - Extract (sub-string) forward sequence from dbSNP151 Reference Genome
  #  
  #
  
  man_add_bsp_tib %>% dplyr::filter(seqnames=="chr10") %>% head() %>% tail(n=3) %>% as.data.frame()
  # cg28108104_F_T_C_II
  # 3102674
  bsp_pos=3102674
  
  # How to extract forward sequence::
  subseq(tmp_fas_dat[["10"]], start=bsp_pos - 60, width = 60)
  subseq(tmp_fas_dat[["10"]], start=bsp_pos +  0, width = 2)
  subseq(tmp_fas_dat[["10"]], start=bsp_pos +  2, width = 60)
  
  fwd_seq <- subseq(tmp_fas_dat[["10"]], start=bsp_pos - 60, width = 122) %>% as.character()
  rev_seq <- fwd_seq %>% revCmp()
  
  fwd_brac_seq <- fwd_seq %>% addBrac()
  rev_brac_seq <- rev_seq %>% addBrac()
  
  
  # cgn_pos_db_tsv <- file.path("/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designOutput/GRCh37-21022021_improbe-designOutput.cgn-pos-fwd.tsv.gz")
  
  cur_pos_db_grs <- write_impTopGenomeGRS(genBuild=opt$genBuild, verbose=10, tt=pTracker)
  cur_pos_db_grs %>% dplyr::filter(imp_chr==10 && imp_pos==bsp_pos)
  
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
  
  # imp_pos_rds <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCm10-21092020_improbe-designOutput.cgn-pos-srd.rds"
  imp_pos_dir <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020"
  imp_pos_rds <- file.path(imp_pos_dir, paste(opt$genBuild, "21092020_improbe-designOutput.cgn-pos-srd.rds", sep='-'))
  imp_pos_grs <- readr::read_rds(imp_pos_rds)
  
  imp_top_tsv <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/cgnTop/GRCh36-GRCh38-GRCm10-21092020.cgnTop.sorted.tsv.gz"
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
  # seqnames   start     end width strand              ord_id  prb_add  prb_cgn prb_des                                        prb_prb_seq
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
  
  
}



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
