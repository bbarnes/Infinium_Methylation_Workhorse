
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
par$prgmDir <- 'analysis'
par$prgmTag <- 'COVIC_analyze_signals'
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
  par$local_runType <- 'GRCm10'
  par$local_runType <- 'COVIC'
  
  if (par$local_runType=='COVIC') {
    opt$genBuild <- 'GRCh36'
    opt$genBuild <- 'GRCh38'
    opt$genBuild <- 'GRCh37'
    
    opt$platform    <- 'COVIC'
    opt$version     <- 'C0'
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
    opt$version  <- 'C26'
    
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

opt$rawDir <- file.path(opt$outDir, "raw")
if (!dir.exists(opt$rawDir)) dir.create(opt$rawDir, recursive=TRUE)
ssh_pas_csv <- file.path(opt$rawDir, paste(opt$runName,"Labled_Sample_Sheets.csv.gz", sep="_"))

data_dir <- "/Users/bretbarnes/Documents/data/COVIC/COVIC-Signals"

if (!file.exists(ssh_pas_csv)) {

  ssh_lab_csv <- "/Users/bretbarnes/Documents/data/sampleSheets/raw/sampleSheets/annotation/Human-Classification_COVID_Count-921_AnnotatedOnlySampleSheet.csv.gz"
  ssh_lab_tib <- readr::read_csv(ssh_lab_csv)
  
  ssh_csv_list <- list.files(path=data_dir, pattern="AutoSampleSheet.csv.gz$", full.names=TRUE)
  
  ssh_all_tib <- lapply(ssh_csv_list, readr::read_csv) %>% dplyr::bind_rows()
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Merge Human Classes::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ssh_hum_tib <- 
    ssh_lab_tib %>% 
    dplyr::select(Sentrix_Name:Tissue_Source) %>% 
    dplyr::inner_join(ssh_all_tib, by="Sentrix_Name")
  
  ssh_hum_sum <- 
    ssh_hum_tib %>% 
    dplyr::group_by(Sample_Class,COVID_Status,Tissue_Source) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  ssh_hum_sum %>% print(n=base::nrow(ssh_hum_sum))
  
  # Only Passing Human COVID+/- Samples moving forward::
  ssh_pas_tib <- 
    ssh_hum_tib %>% 
    dplyr::filter(Sample_Class=="nSARSCov2" | Sample_Class=="pSARSCov2") %>%
    dplyr::filter(Tissue_Source=="Whole_Blood") %>%
    dplyr::group_by(Sample_Class) %>% 
    dplyr::mutate(Class_Rank=dplyr::row_number(), Class_Tag=paste(Sample_Class,Class_Rank, sep="_")) %>%
    dplyr::ungroup()
  readr::write_csv(ssh_pas_tib, ssh_pas_csv)
  
} else {
  ssh_pas_tib <- readr::read_csv(ssh_pas_csv)
}

ssh_pas_sum <- 
  ssh_pas_tib %>% 
  dplyr::group_by(Sample_Class,COVID_Status,Tissue_Source) %>% 
  dplyr::summarise(Count=n(), .groups="drop")
ssh_pas_sum %>% print(n=base::nrow(ssh_pas_sum))


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Load Signal Data::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sig_met <- "ind"
raw_dat_rds <- file.path(opt$rawDir, paste(opt$runName,sig_met,"signal-data.rds", sep="_"))

ssh_pat  <- "AutoSampleSheet.csv.gz"
ssh_csvs <- list.files(data_dir, pattern=ssh_pat, full.names=TRUE)

ssh_cnt <- ssh_pas_tib %>% base::nrow()
# ssh_cnt <- 30

sig_tab <- NULL
sig_row <- NULL

if (!file.exists(raw_dat_rds)) {
  
  for (ii in c(1:ssh_cnt)) {
    sentrix_name <- ssh_pas_tib$Sentrix_Name[ii]
    sample_class <- ssh_pas_tib$Class_Tag[ii]
    # sample_class <- ssh_pas_tib$Sample_Class[ii]
    
    sig_csv <- file.path(data_dir, paste(
      paste(sentrix_name,opt$platform,opt$version,sig_met, sep="_"),
      "sigs.dat.csv.gz", sep=".")
    )
    
    cat(glue::glue("[{par$prgmTag}]: Reading sig_csv({ii}/{sample_class})={sig_csv}{RET}") )
    sig_tib <- suppressMessages(suppressWarnings( readr::read_csv(sig_csv) )) %>%
      dplyr::filter(!is.na(sig_U)) %>%
      dplyr::filter(!is.na(sig_M))
    
    # Probably should normalize each group by their type...
    # sig_tib %>% dplyr::group_by(Probe_Design) %>% dplyr::summarise(M_min=min(sig_M),M_max=max(sig_M),U_min=min(sig_U),U_max=max(sig_U))
    
    cur_min_U <- base::min(sig_tib$sig_U)
    cur_min_M <- base::min(sig_tib$sig_M)
    
    cur_max_U <- base::max(sig_tib$sig_U)
    cur_max_M <- base::max(sig_tib$sig_M)
    
    cat(glue::glue("[{par$prgmTag}]: Building transform; class={sample_class}; ",
                   "min={cur_min_U}/{cur_min_M}, ",
                   "max={cur_max_U}/{cur_max_M}).{RET}{RET}") )
    
    cur_tib <- sig_tib %>% # head() %>% 
      dplyr::mutate(
        sig_U=(sig_U-cur_min_U)/cur_max_U,
        sig_M=(sig_M-cur_min_M)/cur_max_M
      ) %>%
      dplyr::mutate(Probe_Design=stringr::str_remove(Probe_Design, "^[IO]")) %>% 
      dplyr::select(-Probe_Type) %>% 
      tidyr::unite(Probe_ID, Probe_ID,Probe_Design, sep="_") %>% 
      dplyr::rename(M=sig_M,U=sig_U) %>% 
      tidyr::pivot_longer(cols=c("M","U")) %>% 
      tidyr::unite(Probe_ID, Probe_ID,name, sep="") %>% 
      dplyr::arrange(Probe_ID) %>% 
      # tibble::column_to_rownames(var="Probe_ID") %>% 
      # purrr::set_names(c(dplyr::all_of(c("Probe_ID",sample_class) ) ))
      purrr::set_names(c(dplyr::all_of(c("Probe_ID",sample_class) ) ))
    
    if (is.null(sig_row)) {
      sig_row <- cur_tib
    } else {
      sig_row <- dplyr::full_join(sig_row, cur_tib, by="Probe_ID")
    }
    
    cat(glue::glue("[{par$prgmTag}]: Done class={sample_class}; ",
                   "min={cur_min_U}/{cur_min_M}, ",
                   "max={cur_max_U}/{cur_max_M}).{RET}{RET}") )
  }
  
  if (FALSE && !file.exists(raw_dat_rds))
    readr::write_rds(sig_row, raw_dat_rds, compress="gz")
  
} else {
  cat(glue::glue("[{par$prgmTag}]: Loading Experiment RDS={raw_dat_rds}...{RET}"))
  sig_row <- readr::read_rds(raw_dat_rds)
}

sig_mat <- sig_row %>% tibble::column_to_rownames(var="Probe_ID") %>% as.matrix() %>% t()


sig_df  <- sig_mat %>% as.data.frame() %>% 
  tibble::rownames_to_column(var="Class") %>% 
  dplyr::mutate(Class=stringr::str_remove(Class, "\\..*$"))

sig_mat[1:3,1:3] %>% print()
sig_row_tib <- sig_mat %>% rownames() %>% 
  tibble::enframe(name = "Row", value = "Class") %>% 
  tidyr::separate(Class, into=c("Class_Name", "Class_Rank"), sep="_")
sig_row_sum <- sig_row_tib %>% 
  dplyr::group_by(Class_Name) %>% 
  dplyr::summarise(Count=n(), .groups="drop")
sig_row_sum %>% print(n=base::nrow(sig_row_sum))

# Set Row Names in Matrix to Class::
rownames(sig_mat) <- sig_row_tib$Class_Name

# Get CLass Vector::
sig_class_vec <- sig_row_tib %>% 
  dplyr::mutate(Class_Fact=as.factor(Class_Name)) %>% 
  dplyr::pull(Class_Fact) %>% as.integer()

test_df <- sig_mat[,1:30] %>% as.data.frame() %>% dplyr::mutate(Class=sig_class_vec)
rf_classifier = randomForest(Class ~ ., data=test_df, ntree=100, mtry=2, importance=TRUE)


# Simple Test::
#  rf_classifier = randomForest(x=sig_mat, y=sig_class_vec)

t_len <- 50
rf_classifier = randomForest(Class ~ ., data=sig_df[,1:t_len], ntree=100, mtry=2, importance=TRUE)



rf_classifier = randomForest(x=sig_mat[ ,1:t_len], y=sig_class_vec)
predict( rf_classifier, sig_mat[,1:t_len]) %>% tibble::enframe() %>% 
  dplyr::arrange(name) %>% print(n=100)


# Get 1/3 data sets
sig_len <- sig_mat %>% rownames() %>% length()
idx1_vec <- which(1:sig_len %% 3 == 1)
idx2_vec <- which(1:sig_len %% 3 == 2)
idx3_vec <- which(1:sig_len %% 3 == 0)

# sig_mat[idx1_vec, ] %>% dim()

sig_mat_12 <- rbind(
  sig_mat[idx1_vec, ],
  sig_mat[idx2_vec, ]
) %>% dim()






# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Random Forest Example::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Set random seed to make results reproducible:
# set.seed(17)
# Calculate the size of each of the data sets:
data_set_size <- floor(nrow(iris)/2)
# Generate a random sample of "data_set_size" indexes
indexes <- sample(1:nrow(iris), size = data_set_size)
# Assign the data to the correct sets
training <- iris[indexes,]
validation1 <- iris[-indexes,]


#import the package
library(randomForest)
# Perform training:
rf_classifier = randomForest(Species ~ ., data=training, ntree=100, mtry=2, importance=TRUE)


# Bret Comparison::
dat_mat <- training[, 1:4]
key_vec <- training[, 5] %>% as.factor() %>% as.integer()

# randomForest(Species ~ ., data=training, ntree=100, mtry=2, importance=TRUE)
randomForest(x=dat_mat, y=key_vec, ntree=100, mtry=2, importance=TRUE)


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Split Dataset::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# sigs_tib$pSARSCov2 %>% head(n=20) %>% tibble::column_to_rownames(var="Probe_ID") %>% as.matrix() %>% t()



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           OLD CODE TO BE DELTED::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# Old method::
if (FALSE) {
  
  if (FALSE) {
    cur_tib <- sig_tib %>% 
      dplyr::mutate(
        sig_U=(sig_U-cur_min_U)/cur_max_U,
        sig_M=(sig_M-cur_min_M)/cur_max_M
      ) %>%
      dplyr::mutate(Probe_Design=stringr::str_remove(Probe_Design, "^[IO]")) %>% 
      dplyr::select(-Probe_Type) %>% 
      tidyr::unite(Probe_ID, Probe_ID,Probe_Design, sep="_") %>% 
      dplyr::rename(M=sig_M,U=sig_U) %>% 
      tidyr::pivot_longer(cols=c("M","U")) %>% 
      tidyr::unite(Probe_ID, Probe_ID,name, sep="") %>% 
      dplyr::arrange(Probe_ID) %>% 
      tibble::column_to_rownames(var="Probe_ID") %>% 
      #
      # Need to break here and join by full_join()
      #
      t() %>% 
      as.data.frame() %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(Class := sample_class) %>% 
      dplyr::select(Class, dplyr::everything())
    
    sig_tab <- sig_tab %>% dplyr::bind_rows(cur_tib)
  }
  
  cur_col <- c("Probe_ID",
               paste("U",sentrix_name, sep="_"),
               paste("M",sentrix_name, sep="_"))
  
  cur_tib <- sig_tib %>% 
    dplyr::mutate(
      Probe_Class=dplyr::case_when(
        Probe_Design=="OR" ~ 0,  # Unmethylated (Out-of-bound)
        Probe_Design=="OG" ~ 1,  # Methylated (Out-of-bound)
        Probe_Design=="2"  ~ 2,  # Type 2
        Probe_Design=="IR" ~ 3,  # Unmethylated (In-bound)
        Probe_Design=="IG" ~ 4,  # Methylated (In-bound)
        TRUE ~ 5
      ) %>% as.integer(),
      Probe_ID=paste(Probe_ID,Probe_Class, sep="_")
    ) %>% 
    dplyr::select(Probe_ID,sig_U,sig_M) %>%
    dplyr::filter(!is.na(sig_U)) %>%
    dplyr::filter(!is.na(sig_M))
  
  cur_min_U <- base::min(cur_tib$sig_U)
  cur_min_M <- base::min(cur_tib$sig_M)
  
  cur_max_U <- base::max(cur_tib$sig_U)
  cur_max_M <- base::max(cur_tib$sig_M)
  
  cur_tib <- cur_tib %>% dplyr::mutate(
    sig_U=(sig_U-cur_min_U)/cur_max_U,
    sig_M=(sig_M-cur_min_M)/cur_max_M
  ) %>% purrr::set_names(cur_col)
  
  if (is.null(sigs_tib[[sample_class]])) {
    sigs_tib[[sample_class]] <- cur_tib
  } else {
    sigs_tib[[sample_class]] <- sigs_tib[[sample_class]] %>% 
      dplyr::full_join(cur_tib, by="Probe_ID")
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
