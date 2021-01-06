
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
par$prgmDir <- 'probe_design'
par$prgmTag <- 'improbe_EWAS_scratch'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

par$local_runType <- NULL

par$cpg_wat_int_csv <- NULL

# Executables::
opt$Rscript <- NULL

# Run Parameters::
opt$runName   <- NULL

# Directories::
opt$outDir <- NULL
opt$impDir <- NULL
opt$annDir <- NULL

# Required Inputs::
opt$ords <- NULL
opt$mats <- NULL
opt$aqps <- NULL
opt$pqcs <- NULL
opt$ctls <- NULL
opt$idat <- NULL

opt$cpg_s48_tsv <- NULL
opt$cpg_top_tsv <- NULL

opt$cpg_pos_tsv <- NULL
opt$cph_pos_tsv <- NULL
opt$snp_pos_tsv <- NULL

opt$pre_man_csv <- NULL

# opt$cpg_des_csv <- NULL
opt$cph_des_csv <- NULL
opt$snp_des_csv <- NULL

# Platform/Method Options::
opt$genomeBuild <- NULL
opt$platform    <- NULL
opt$version     <- NULL

# Run Options::
opt$fresh   <- FALSE
opt$fixIds  <- FALSE
par$retData <- FALSE

opt$matFormat <- 'new'
opt$ordFormat <- 'old'

opt$matSkip <- 40
opt$ordSkip <- 8
opt$aqpSkip <- 7
opt$pqcSkip <- 7

# Parallel/Cluster Options::
opt$single   <- FALSE
opt$parallel <- TRUE
opt$cluster  <- FALSE

opt$opt_csv  <- NULL
opt$par_csv  <- NULL
opt$time_csv <- NULL
opt$time_org_txt <- NULL

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
  opt$impDir <- file.path(par$topDir, 'data/improbe')
  opt$annDir <- file.path(par$topDir, 'data/annotation')
  
  # opt$cpg_s48_tsv <- file.path(opt$impDir, 'designOutput_21092020/seq48U/gz/seq48U-GRCh36-38-10-21092020.unq.noHeader.seq-sorted.tsv.gz')
  opt$cpg_s48_tsv <- file.path(opt$impDir, 'designOutput_21092020/seq48U/un/seq48U-GRCh36-38-10-21092020.unq.noHeader.seq-sorted.tsv')
  
  # The reason we don't use the general form is that specificity for the species may be lost.
  #   I think this is poor form, but using now to meet mouse GS requriments...
  #
  # opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop/GRCh36-GRCh38-GRCm10-21092020.cgnTop.sorted.tsv.gz')
  # opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop/GRCh36-GRCh38-GRCm10-21092020.cgnTop.sorted.tsv')
  
  opt$write_full <- FALSE
  opt$write_base <- FALSE
  
  #
  # Pre-defined local options runTypes::
  #
  par$local_runType <- 'qcMVP'
  par$local_runType <- 'NZT'
  par$local_runType <- 'COVID'
  par$local_runType <- 'GENK'
  par$local_runType <- 'COVIC'
  par$local_runType <- 'GRCm38'
  par$local_runType <- 'EWAS'
  
  if (par$local_runType=='EWAS') {
    opt$genomeBuild <- 'GRCh36'
    opt$genomeBuild <- 'GRCh38'
    opt$genomeBuild <- 'GRCh37'
    opt$platform    <- 'EWAS'
    opt$version     <- 'E0'
    opt$version     <- 'E1'
    
    opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genomeBuild,'-21092020.cgnTop.sorted.tsv') )
    opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genomeBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
    opt$pre_man_csv <- paste(
      '/Users/bretbarnes/Documents/data/manifests/MethylationEPIC_v-1-0_B2.csv.gz',
      '/Users/bretbarnes/Documents/data/manifests/HumanMethylation450_15017482_v.1.2.csv.gz',
      '/Users/bretbarnes/Documents/data/manifests/HumanMethylation27_270596_v.1.2.csv.gz',
      sep=',')
    
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$parallel <- TRUE
  opt$runName <- paste(opt$genomeBuild,opt$platform,opt$version, sep='-')
  
  par$cpg_wat_int_csv <- file.path(opt$outDir, paste(opt$runName,"Waterland_CORSIVS_int.tsv.gz", sep="_") )
  
  opt$fresh <- FALSE
  opt$fresh <- TRUE
  
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
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    
    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    make_option(c("--impDir"), type="character", default=opt$impDir, 
                help="improbe data directory [default= %default]", metavar="character"),
    make_option(c("--annDir"), type="character", default=opt$iannDir, 
                help="Annotation data directory [default= %default]", metavar="character"),
    
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
    
    # make_option(c("--cpg_des_csv"), type="character", default=opt$cpg_des_csv, 
    #             help="CpG Position TSV File [default= %default]", metavar="character"),
    make_option(c("--cph_des_csv"), type="character", default=opt$cph_des_csv, 
                help="CpH Position TSV File [default= %default]", metavar="character"),
    make_option(c("--snp_des_csv"), type="character", default=opt$snp_des_csv, 
                help="SNP Position TSV File [default= %default]", metavar="character"),
    
    # Platform/Method Options::
    make_option(c("--genomeBuild"), type="character", default=opt$genomeBuild, 
                help="Genome Build (e.g. GRCh36, GRCh37, GRCh38, GRCm38) [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform (e.g. HM450, EPIC, LEGX, NZT, COVIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest Version (e.g. B0,B1,B2,B3,B4,C0) [default= %default]", metavar="character"),
    
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
    make_option(c("--aqpSkip"), type="integer", default=opt$aqpSkip, 
                help="AQP Skip Lines Count [default= %default]", metavar="integer"),
    make_option(c("--pqcSkip"), type="integer", default=opt$pqcSkip, 
                help="PQC Skip Lines Count [default= %default]", metavar="integer"),
    
    # Run Options::
    make_option(c("--fresh"), action="store_true", default=opt$fresh, 
                help="Boolean variable to run a fresh build [default= %default]", metavar="boolean"),
    make_option(c("--fixIds"), action="store_true", default=opt$fixIds, 
                help="Boolean variable to fix original order ids [default= %default]", metavar="boolean"),
    make_option(c("--retData"), action="store_true", default=opt$retData, 
                help="Developement ONLY Boolean variable to return data for testing [default= %default]", metavar="boolean"),
    
    # Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
    make_option(c("--opt_csv"), type="character", default=opt$opt_csv, 
                help="Unused variable opt_csv [default= %default]", metavar="character"),
    make_option(c("--par_csv"), type="character", default=opt$par_csv, 
                help="Unused variable par_csv [default= %default]", metavar="character"),
    make_option(c("--time_csv"), type="character", default=opt$time_csv, 
                help="Unused variable time_csv [default= %default]", metavar="character"),
    make_option(c("--time_org_txt"), type="character", default=opt$time_org_txt, 
                help="Unused variable time_org_txt [default= %default]", metavar="character"),
    
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
opt_reqs <- c('outDir','impDir','cpg_s48_tsv','cpg_top_tsv',
              'genomeBuild','platform','version','Rscript','verbose')

par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )

opt <- program_init(name=par$prgmTag,
                    opts=opt, opt_reqs=opt_reqs, 
                    pars=par, par_reqs=par_reqs,
                    libs=TRUE,rcpp=FALSE,
                    verbose=opt$verbose,vt=3,tc=0,tt=NULL)

par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

image_key <- "bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse_Centos"
image_ver <- "v.1.0"
image_ssh <- "run_improbe.sh"
image_str <- glue::glue("{image_key}.{image_ver}")

# opt$probe_type <- 'cg'
# opt$design_key <- 'Seq_ID'
# opt$design_seq <- 'Forward_Sequence'
# opt$design_seq <- 'Top_Sequence'
# opt$design_prb <- 'Probe_Type'
# opt$design_srs <- 'TB'
# opt$design_cos <- 'CO'

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Preprocessing:: Current EWAS+EPIC-V2
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pre_ewas_dir <- "/Users/bretbarnes/Documents/data/CustomContent/EWAS/orders/round1"
pre_epic_dir <- "/Users/bretbarnes/Documents/data/CustomContent/EPIC_v2/11102020"

pre_ewas_pat <- "*order.csv.gz$"
pre_epic_pat <- "bp_[0-9]+.tsv.gz$"

pre_ewas_csvs <- list.files(pre_ewas_dir, pattern=pre_ewas_pat, full.names=TRUE)
pre_epic_tsvs <- list.files(pre_epic_dir, pattern=pre_epic_pat, full.names=TRUE)

# pre_ewas_tib <- lapply(pre_ewas_csvs, readr::read_csv) %>% dplyr::bind_rows() %>%
#   dplyr::select(Assay_Design_Id,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Normalization_Bin) %>%
#   tidyr::separate(Assay_Design_Id, into=c("Seq_ID","Probe_Desc"), sep="_")
# pre_epic_tib <- lapply(pre_epic_tsvs, readr::read_tsv) %>% dplyr::bind_rows() %>%
#   dplyr::select(Assay_Design_Id,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Normalization_Bin) %>%
#   tidyr::separate(Assay_Design_Id, into=c("Seq_ID","Probe_Desc"), sep="_")

# Total Order Probes: 926,093 as of Jan-1-2021
pre_order_csv <- "/Users/bretbarnes/Documents/data/CustomContent/EWASv1_EPICv2_current_reorder.csv.gz"
if (file.exists(pre_order_csv)) {
  pre_order_tib <- readr::read_csv(pre_order_csv)
} else {
  pre_order_tib <- dplyr::bind_rows(
    lapply(pre_ewas_csvs, readr::read_csv) %>% dplyr::bind_rows() %>% 
      dplyr::mutate(Order_Pool="EWASv1"),
    lapply(pre_epic_tsvs, readr::read_tsv) %>% dplyr::bind_rows() %>% 
      dplyr::mutate(Order_Pool="EPICv2")
  ) %>%
    dplyr::select(Assay_Design_Id,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                  Normalization_Bin,Order_Pool) %>%
    tidyr::separate(Assay_Design_Id, into=c("Seq_ID","Probe_Desc"), sep="_") %>%
    dplyr::arrange(Seq_ID)
  
  readr::write_csv(pre_order_tib, pre_order_csv)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             3. CGN to Top::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Load All Method:: for high meomory machines
#
all_fwd_col <- cols(Seq_ID       = col_character(),
                    Sequence     = col_character(),
                    Genome_Build = col_character(),
                    Chromosome   = col_character(),
                    Coordinate   = col_integer(),
                    CpG_Island   = col_logical()
)
all_fwd_tsv <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh37.improbeDesignInput.tsv.gz"
all_fwd_rds <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/rds/GRCh37.improbeDesignInput.rds"
all_grs_rds <- "/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/rds/GRCh37.improbeDesignInput.grs.rds"

# Concerned about this function below::
# all_top_files <- c(all_top_tsv,all_top_rds,all_grs_rds)
# stp_tib <- check_timeStamps(name="load_full_top",outDir=opt$outDir,
#                             origin=NULL,
#                             files=all_top_files,
#                             verbose=opt$verbose,tt=pTracker)


#
#
#  Switch to  par$cpg_wat_int_csv instead of cpg_wat_int_csv
#
#

all_fwd_grs <- NULL
if (file.exists(all_grs_rds) && 
    file.mtime(all_grs_rds) >= file.mtime(all_fwd_rds) ) {
  
  #
  # Huge Temp Fix For Fast Processing::
  #
  if (!file.exists(par$cpg_wat_int_csv)) {
    cat(glue::glue("Loading: {all_grs_rds}...{RET}{RET}"))
    load_fwd_grs_time <- system.time({
      all_fwd_grs <- readr::read_rds(all_grs_rds)
    })
    pTracker$addTime(time=load_fwd_grs_time, key="load_fwd_grs_time")
    
    stopifnot(!is.null(all_fwd_grs))
  }
  
} else {
  
  #
  # Skipping all of this for now since the above should be available::
  #
  if (FALSE) {

    if (file.exists(all_fwd_rds) && 
        file.mtime(all_fwd_rds) >= file.mtime(all_fwd_tsv) ) {
      
    } else {
      
      # Load original forward cpg tsv
      load_fwd_tsv_time <- system.time({
        all_fwd_tib <- readr::read_tsv(all_fwd_tsv,col_types=all_fwd_col) %>%
          dplyr::mutate(Genome_Build=as.factor(Genome_Build),
                        Strand=as.factor("+"), 
                        Chromosome=as.factor(paste0("chr",Chromosome)), 
                        Rep=dplyr::row_number()) %>% 
          tidyr::unite(Seq_ID, Seq_ID,Rep, sep="_")
      })
      
      # Write original forward cpg RDS::
      write_fwd_cpg_time <- system.time({
        readr::write_rds(all_fwd_tib, all_fwd_rds, compress="gz")
      })
      
      # Load/Build/Write Genomic Ranges
      #
      build_fwd_grs_tibe <- system.time({
        all_fwd_grs <- GRanges(
          seqnames=Rle(all_fwd_tib$Chromosome), strand=Rle(all_fwd_tib$Strand), 
          Sequence=all_fwd_tib$Sequence, Genome_Build=all_fwd_tib$Genome_Build,
          IRanges(start=all_fwd_tib$Coordinate, end=all_fwd_tib$Coordinate+1, names=all_fwd_tib$Seq_ID) )
      })
      
      # Write original forward Genomic Ranges RDS::
      write_fwd_grs_time <- system.time({
        if (!file.exists(all_grs_rds)) readr::write_rds(all_fwd_grs, all_grs_rds, compress="gz")
      })
    }
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              New Designs::
#
#  1. Build Pre-BED Designs
#
#  2. Load Waterland_CORSIVS Coordinate New Design Files
#    - S2 (Additional file 2),"Annotated list of all genomic bins within the 39,424 CoRSIVs (Unfiltered)",,,,
#    - S3,"Annotated list of all genomic bins within the 9,926 CoRSIVs (Filtered for number of CpGs and Interindividual Range)",,,,
#  3. Intersect 1 vs. 2
#
#  4. Load Cancer_Driver_Mutation_Calling New Coordinate Design Files
#  5. Intersect 1 vs. 4
#
#  4. Load New Low-Coverage Gene Files
#  5. Load Gene Coordinates
#  6. Intersect 1 vs. 5
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# 1. --- Build Pre-BED Designs ----
#
#  NOT NEEDED ANYMORE SINCE WE CAN LOAD EVERYTHING INTO MEMORY!!!
#
if (FALSE) {
  gen_name_col <- cols(IlmnID = col_character(),
                       chrom = col_character(),
                       start = col_integer()
  )
  
  # cpg_gen_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/genomic/GRCh38.improbeDesignInput.cgn-sorted.tsv.gz'
  cpg_gen_tsv <- '/Users/bretbarnes/Documents/data/improbe/designOutput_21092020/genomic/GRCh37.improbeDesignInput.cgn-sorted.tsv.gz'
  cpg_gen_tib <- readr::read_tsv(cpg_gen_tsv, 
                                 col_names=names(gen_name_col$cols),
                                 col_types=gen_name_col) %>% 
    dplyr::mutate(chrom=paste0("chr",chrom), strand="+", Rep=dplyr::row_number()) %>% 
    tidyr::unite(IlmnID, IlmnID,Rep, sep="_")
  
  cpg_gen_grs <- GRanges(
    seqnames=Rle(cpg_gen_tib$chrom), strand=Rle(cpg_gen_tib$strand),
    IRanges(start=cpg_gen_tib$start, end=cpg_gen_tib$start, names=cpg_gen_tib$IlmnID) )
}

#
#  2. Load New Waterland_CORSIVS Coordinate Design Files
#
#    - S2 (Additional file 2),"Annotated list of all genomic bins within the 39,424 CoRSIVs (Unfiltered)",,,,
# new_wat_csv <- file.path(new_cpg_dir, 'Waterland_CORSIVS/Waterland_CORSIV_Genomic_Coordinates_S2-Table_1-13059_2019_1708_MOESM2_ESM.csv.gz')
#
#    - S3,"Annotated list of all genomic bins within the 9,926 CoRSIVs (Filtered for number of CpGs and Interindividual Range)",,,,
#
new_cpg_dir <- '/Users/bretbarnes/Documents/data/CustomContent/EWAS/CORSIVs_and_others/Not_Already_Designed_Content'
new_wat_csv <- file.path(new_cpg_dir, 'Waterland_CORSIVS/Waterland_CORSIV_Genomic_Coordinates_S3-Table_1.csv.gz')
if (opt$version=='E0') {
  new_wat_tib <- readr::read_csv(new_wat_csv) %>% 
    dplyr::select(Uniq_ID,USCS_Coordinates_CoRSIV) %>%
    dplyr::distinct() %>%
    tidyr::separate(USCS_Coordinates_CoRSIV, into=c("chrom","range"), sep=":") %>%
    tidyr::separate(range, into=c("beg","end"), sep="-") %>%
    dplyr::mutate(beg=as.integer(beg), end=as.integer(end), strand="+")
  
  new_wat_grs <- GRanges(
    seqnames=Rle(new_wat_tib$chrom), strand=Rle(new_wat_tib$strand),
    IRanges(start=new_wat_tib$beg, end=new_wat_tib$end, names=new_wat_tib$Uniq_ID) )

} else if (opt$version=='E1') {
  #
  # What we should use::
  #
  #  tmp_tib %>% dplyr::select(Uniq_ID,USCS_Coordinates_CoRSIV,CpG.Count,UCSC.Coordinates)
  #
  new_wat_tib <- readr::read_csv(new_wat_csv) %>% 
    dplyr::select(Uniq_ID,CpG.Count,UCSC.Coordinates,USCS_Coordinates_CoRSIV) %>%
    dplyr::rename(CpG_Count=CpG.Count) %>%
    tidyr::separate(UCSC.Coordinates, into=c("chrom","range"), sep=":") %>%
    tidyr::separate(range, into=c("beg","end"), sep="-") %>%
    dplyr::mutate(beg=as.integer(beg), end=as.integer(end), strand="+",
                  Rep=dplyr::row_number()) %>%
    tidyr::unite(Uniq_ID, Uniq_ID,Rep, sep="_")
    
  new_wat_grs <- GRanges(
    seqnames=Rle(new_wat_tib$chrom), strand=Rle(new_wat_tib$strand),
    CpG_Count=new_wat_tib$CpG_Count, 
    USCS_Coordinates_CoRSIV=new_wat_tib$USCS_Coordinates_CoRSIV,
    IRanges(start=new_wat_tib$beg, end=new_wat_tib$end, names=new_wat_tib$Uniq_ID) )
} else {
  stop(glue::glue("[{par$prgmTag}]: Unsupported option version={opt$version}!{RET}{RET}"))
}

#
# 3. Intersect all_fwd_grs with new_wat_grs::
#
#  TBD:: Move some of the renaming code into the function: intersectGranges()
#

cpg_wat_int_tib <- NULL
if (!is.null(par$cpg_wat_int_csv) & file.exists(par$cpg_wat_int_csv)) {
  cpg_wat_int_tib <- readr::read_csv(par$cpg_wat_int_csv)
} else {
  cpg_wat_int_tib <- intersectGranges(all_fwd_grs,new_wat_grs, 
                                      verbose=opt$verbose, tt=pTracker) %>% 
    tidyr::separate(Seq_ID, into=c("Seq_ID", "Seq_Idx"), sep="_") %>%
    dplyr::mutate(CpG_Island=FALSE,Probe_Type="cg") %>%
    dplyr::rename(Chromosome=seqnames, Coordinate=start) %>%
    dplyr::select(Seq_ID,Sequence,Genome_Build,Chromosome,Coordinate,CpG_Island,
                  Probe_Type,dplyr::everything()) %>%
    dplyr::rename(Gene_Beg=chromStart,Gene_End=chromEnd) %>% 
    dplyr::select(-c(end,width,strand,chrom,chromLength,chromStrand))
  
  if (!is.null(par$cpg_wat_int_csv))
    readr::write_csv(cpg_wat_int_tib,par$cpg_wat_int_csv)
}

#
# Remove some large data structures::
#
if (!is.null(all_fwd_grs)) rm(all_fwd_grs)

#
# 4. Build all designs::
#
des_all_wat_tib <- cpg_wat_int_tib %>% # head() %>%
  improbe_design_all(ptype="cg", outDir=opt$outDir, 
                     gen="all", image=image_str, shell=image_ssh, 
                     seqKey="Sequence",strsSR="FR",reduce_imp=TRUE,
                     verbose=opt$verbose, tt=pTracker)

# - Remove Loci Already in Current Re-order::
# - Basic Filtering::
# - [after alignment] Add Group Annodation::
#
# Write FASTA!!! and do Alignment for screening!!!
#
des_pas_wat_tib <- des_all_wat_tib %>% 
  dplyr::filter(!Seq_ID %in% pre_order_tib$Seq_ID) %>%
  dplyr::filter(Probe_Score_Min>=0.3) %>% 
  dplyr::mutate(Probe_Score_Min=as.integer(Probe_Score_Min*100),
                Next_Base_Inf1=NXB_D,Next_Base_Inf2=CPN_D,
                DT1=1,DT2=2) %>%
  # dplyr::select(Seq_ID,Strand_TB,Strand_CO,Probe_Score_Min, 
  #               Next_Base_Inf1,Next_Base_Inf2,DT1,DT2,
  #               Underlying_CpG_Count,Strand_SR,Underlying_CpG_Min_Dist,
  #               Seq_48U_1, Seq_48U_2) %>%
  tidyr::unite(Probe_Desc1, Strand_TB,DT1,Strand_CO,Probe_Score_Min,Next_Base_Inf1,
               Underlying_CpG_Count,Strand_SR,Underlying_CpG_Min_Dist, 
               sep="", remove=FALSE) %>%
  tidyr::unite(Probe_Desc2, Strand_TB,DT2,Strand_CO,Probe_Score_Min,Next_Base_Inf2,
               Underlying_CpG_Count,Strand_SR,Underlying_CpG_Min_Dist, 
               sep="", remove=FALSE) %>%
  tidyr::unite(Seq_ID1, Seq_ID,Probe_Desc1, sep="_", remove=FALSE) %>%
  tidyr::unite(Seq_ID2, Seq_ID,Probe_Desc2, sep="_", remove=FALSE)


#
# Write FASTA:: Run BSMAP!!!
#
aln_fas <- file.path(opt$outDir, paste0(opt$runName,"aln.fas.gz"))
aln_tib <- des_pas_wat_tib %>% 
  dplyr::select(Seq_ID1,Seq_48U_1,Seq_ID2,Seq_48U_2) %>%
  dplyr::mutate(line=paste0('>',Seq_ID1,'\n',Seq_48U_1,'\n',
                            '>',Seq_ID2,'\n',Seq_48U_2)) %>% 
  dplyr::pull(line)

readr::write_lines(x=aln_tib, file=aln_fas)

#
# Screen Alignment Results::
#
aln_bsp_col <- cols(Bsp_ID       = col_character(),
                    Bsp_Seq      = col_character(),
                    Bsp_Qual     = col_character(),
                    Bsp_Map      = col_character(),
                    Bsp_Chr      = col_character(),
                    Bsp_Pos      = col_integer(),
                    Bsp_Srd      = col_character(),
                    Bsp_Mis_Cnt  = col_integer(),
                    Bsp_Aln_Seq  = col_character(),
                    Bsp_Gap_Cnt  = col_integer(),
                    Bsp_Mis_Str  = col_character()
)

aln_bsp_tsv <- "/Users/bretbarnes/Documents/data/CustomContent/EWAS/alignment/hg19/bsmap/cpg/cpg-hg19.s12.v5.g0.p16.n1.r2.Ru.bsp.gz"
aln_bsp_tib <- readr::read_tsv(aln_bsp_tsv,
                               col_names=names(aln_bsp_col$cols),
                               col_types=aln_bsp_col) %>%
  dplyr::select(-Bsp_Qual) %>% 
  tidyr::separate(Bsp_ID, into=c("Seq_ID","Probe_Desc"), sep="_", remove=FALSE)

aln_bsp_list <- aln_bsp_tib %>% split(.$Bsp_Map)

#
# Passing Signle Alignments::
#
sig_pas_wat_tib <- des_pas_wat_tib %>% 
  dplyr::filter(Seq_ID1 %in% aln_bsp_list$UM$Bsp_ID & 
                  Seq_ID2 %in% aln_bsp_list$UM$Bsp_ID)

#  dplyr::distinct(Seq_ID)


#
# Annotation Joining::
# - [after alignment] Add Group Annodation::
#
des_ann_wat_tib <- cpg_wat_int_tib %>% 
  dplyr::inner_join(sig_pas_wat_tib, 
                    by=c("Seq_ID","Genome_Build","Chromosome","Coordinate",
                         "Sequence"="Forward_Sequence","Probe_Type")) %>%
  dplyr::add_count(Gene, name="Gene_Loci_Count") %>% 
  dplyr::select(Gene,Gene_Loci_Count,Seq_ID,
                Strand_SR,Strand_TB,Strand_CO,Probe_Score_Min) %>%
  dplyr::arrange(Gene,-Probe_Score_Min)

# This gives you the top 1 for each gene::
sel1_ann_wat_tib <- des_ann_wat_tib %>% dplyr::group_by(Gene) %>% 
  dplyr::top_n(n=1, wt=Probe_Score_Min) %>% dplyr::ungroup() %>% 
  dplyr::distinct(Gene, .keep_all=TRUE) %>% 
  dplyr::add_count(Gene, name="Gene_Count") %>% 
  dplyr::filter(Gene_Count == 1)

# How many genes are left???
#

#
# Questions::
#  1. Missing any previous loci?
#  2. Ion Torrent Sites? Rishi?
#  3. 
#


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           2.3 Build Probes:: CpG
#
# From: build_bead_manifest_simple.R
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  man_prb_list <- NULL
  
  probe_type <- 'cg'
  man_prb_list[[probe_type]] <- NULL
  if (!is.null(opt$cpg_s48_tsv) && file.exists(opt$cpg_s48_tsv) &&
      !is.null(opt$cpg_top_tsv) && file.exists(opt$cpg_top_tsv)) {
    man_prb_list[[probe_type]] =
      clean_manifest_probes(
        tib=man_raw_dat_tib,
        s48_tsv=opt$cpg_s48_tsv,top_tsv=opt$cpg_top_tsv,
        name=opt$runName,outDir=opt$outDir,origin=opt$time_org_txt,
        
        design_key=opt$design_key,design_seq=opt$design_seq,
        design_prb=opt$design_prb,probe_type=probe_type,
        design_srs=opt$design_srs,design_cos=opt$design_cos,
        parallel=opt$parallel,fresh=opt$fresh,
        
        verbose=opt$verbose,vt=3,tc=0,tt=pTracker)
  }
}




#
#  4. Load Cancer_Driver_Mutation_Calling New Coordinate Design Files
#
new_can_csv <- file.path(new_cpg_dir, 'Cancer_Driver_Mutation_Calling/IonTorrent_Cancer_Hotspot_Panel_v2.csv.gz')
new_can_tib <- readr::read_csv(new_can_csv)

#
# 4. Load New Low-Coverage Gene Files
#
new_pro_csv <- file.path(new_cpg_dir, 'Promoter_Methylation_Cancer_Associated_Genes/Genes_Whose_CpG_Islands_Should_Be_More_Completely_Covered_111720.csv.gz')
new_cnv_csv <- file.path(new_cpg_dir, 'Enhanced_CNV_Calling_for_CNS_Tumor_Diagnostics/Genes_Whose_Exons_Should_Be_More_Completely_Covered_111720.csv.gz')


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Pre Designs::
#
#  1. Load Target CpG Designs
#  2. Load Manifests CpG Designs
#     - Load Ordered CpGs from Cluster
#  3. Intersect Target and Manifest CpGs
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# 1. Load Target CpG Designs::
#
pre_cpg_dir  <- '/Users/bretbarnes/Documents/data/CustomContent/EWAS/CORSIVs_and_others/Already_Designed_Content'
pre_cpg_tsv1 <- file.path(pre_cpg_dir, 'Database_Associations/EWAS_Catalog_03-07-2019.txt.gz')
pre_cpg_tsv2 <- file.path(pre_cpg_dir, 'Database_Associations/EWAS_Atlas_associations.tsv.gz')
pre_cpg_csv3 <- file.path(pre_cpg_dir, 'CpGs_for_Age_Labs/cpg_sites_of_interest_Age_Labs.csv.gz')

pre_cpg_tib1 <- readr::read_tsv(pre_cpg_tsv1) %>% 
  dplyr::rename(IlmnID=CpG) %>%
  dplyr::select(IlmnID) %>% 
  dplyr::distinct(IlmnID) %>% 
  dplyr::mutate(Request_Group="Pre_Design", Request_Name="EWAS_Catalog")

pre_cpg_tib2 <- readr::read_tsv(pre_cpg_tsv2) %>%
  dplyr::rename(IlmnID=`Probe id`) %>%
  dplyr::select(IlmnID) %>% 
  dplyr::distinct(IlmnID) %>% 
  dplyr::mutate(Request_Group="Pre_Design", Request_Name="EWAS_Atlas")

pre_cpg_tib3 <- readr::read_csv(pre_cpg_csv3) %>% 
  dplyr::rename(IlmnID=ProbeID) %>%
  dplyr::select(IlmnID) %>% 
  dplyr::distinct(IlmnID) %>% 
  dplyr::mutate(Request_Group="Pre_Design", Request_Name="Age_Labs")

pre_all_tib <- pre_cpg_tib1 %>% 
  dplyr::bind_rows(pre_cpg_tib2) %>% 
  dplyr::bind_rows(pre_cpg_tib3)

pre_unq_tib <- pre_all_tib %>% 
  dplyr::distinct(IlmnID, .keep_all=TRUE) %>%
  dplyr::rename(Assay_Design_Id=IlmnID)

#
# 2.0 Intersect Target and Pre-Ordered CpGs::
#

pre_unq_order_mis_tib <- pre_unq_tib %>% 
  dplyr::anti_join(pre_order_tib, by=c("Assay_Design_Id"="Seq_ID"))

pre_unq_order_mis_sum <- pre_unq_order_mis_tib %>% 
  dplyr::group_by(Request_Group,Request_Name) %>% 
  dplyr::summarise(Count=n())

if (opt$verbose>=1) {
  cat(glue::glue("[{par$prgmTag}]: pre_unq_order_mis_sum={RET}"))
  pre_unq_order_mis_sum %>% print(n=base::nrow(pre_unq_order_mis_sum))
}


#
# 2.1 Load Manifests CpG Designs::
#     - The goal here is identify where the missing probes came from::
#
sel_man_col <- c("Assay_Design_Id", "AlleleA_ProbeSeq","AlleleB_ProbeSeq",
                 "Man_Source", "Next_Base", "Man_Source")
pre_man_vec <- splitStrToVec(opt$pre_man_csv)

pre_man_tib <- NULL
for (ii in c(1:length(pre_man_vec))) {
  cur_man_tib <- pre_man_vec[ii] %>%
    loadManifestGenomeStudio(addSource=TRUE, normalize=TRUE, retType="man", 
                             verbose=opt$verbose, tt=pTracker) %>%
    dplyr::rename(Assay_Design_Id=IlmnID) %>%
    dplyr::select(dplyr::all_of(sel_man_col))
  
  pre_man_tib <- pre_man_tib %>% dplyr::bind_rows(cur_man_tib)
}

pre_ord_man_tib <- pre_man_tib %>% 
  dplyr::mutate(
    Normalization_Bin=dplyr::case_when(
      is.na(Next_Base) ~ 'C', 
      Next_Base == 'A' | Next_Base == 'T' ~ 'A', 
      Next_Base == 'C' | Next_Base == 'G' ~ 'B', 
      TRUE ~ NA_character_)
  ) %>% 
  dplyr::mutate(AlleleA_Probe_Id=paste(Assay_Design_Id,"A", sep="_"),
                AlleleB_Probe_Id=paste(Assay_Design_Id,"B", sep="_")) %>%
  dplyr::select(dplyr::all_of(sel_man_col), dplyr::everything()) %>%
  dplyr::select(Assay_Design_Id,
                AlleleA_Probe_Id,AlleleA_ProbeSeq,
                AlleleB_Probe_Id,AlleleB_ProbeSeq,
                Normalization_Bin, dplyr::everything() )

org_man_unq_tib <- pre_ord_man_tib %>% 
  dplyr::distinct(Assay_Design_Id, .keep_all=TRUE)


#
# 3. Intersect Target and Manifest CpGs::
#
pre_unq_tib %>% dplyr::filter( Assay_Design_Id %in% org_man_unq_tib$IlmnID)
pre_unq_tib %>% dplyr::filter(!Assay_Design_Id %in% org_man_unq_tib$IlmnID)

#
# Match the missing probes::
#
pre_unq_order_mis_tib %>% dplyr::filter(IlmnID %in% pre_man_tib$IlmnID)
pre_man_tib %>% dplyr::filter(IlmnID %in% pre_unq_order_mis_tib$IlmnID)


pre_man_tib %>% 
  dplyr::filter(IlmnID %in% pre_unq_order_mis_tib$IlmnID) %>% 
  dplyr::distinct(IlmnID, .keep_all=TRUE)










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
