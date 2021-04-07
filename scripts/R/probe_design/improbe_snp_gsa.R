
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Should improve after new cgnDB is built...
#
# use Digest::SHA qw(sha1_hex);
# my $var = 123;
# my $sha1_hash = sha1_hex($var);
# print $sha1_hash;
#
# Need to add this command to all tibble functions::
#
#    select(where(~sum(!is.na(.x)) > 0)) %>%
#    utils::type.convert() %>% 
#    dplyr::mutate(across(where(is.factor), as.character) )
#

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
par$prgmTag <- 'improbe_snp_gsa'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))


# Executables::
opt$Rscript <- NULL

# BSMAP Parameters::
opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
opt$bsmap_exe <- "/Users/bretbarnes/Documents/programs/BSMAPz/bsmapz"

# Run Parameters::
opt$runName    <- NULL
opt$Species    <- NULL

# Null Place Holders::
opt$cpg_top_tsv <- NULL
opt$cpg_pos_tsv <- NULL
opt$cph_pos_tsv <- NULL
opt$snp_pos_tsv <- NULL
opt$org_des_tsv <- NULL

# Directories::
opt$outDir  <- NULL
opt$impDir  <- NULL
opt$annDir  <- NULL
opt$genDir  <- NULL
opt$manDir  <- NULL

# Manufacturing Files:: Required
opt$ords <- NULL
opt$mats <- NULL
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

#
# Data Structures to be pre-defined
#
org_des_tib <- NULL
add_org_tib <- NULL
cgn_bed_tib <- NULL

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
  par$local_runType <- 'GSA'
  
  if (par$local_runType=='GSA') {
    opt$genBuild <- 'GRCh38'
    opt$genBuild <- 'GRCh37'
    
    opt$platform <- 'EPIC'
    opt$version  <- 'A2'

    opt$Species <- "Human"
    
    opt$idat   <- NULL
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/UnivChicago/latest')
    
    opt$ords <- paste(
      file.path(par$aqpDir, 'UofChicago-A_A_Array-CpG-order-FINAL.csv'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20504790_probes.match.tsv'),
      sep=',')
    
    opt$aqps <- NULL
    
    opt$pqcs <- paste(
      file.path(par$aqpDir, '329922X374054_A_ProductQC.txt'),
      sep=',')
    
    # opt$org_des_tsv <- file.path(par$topDir, "data/CustomContent/UnivChicago/improbe_input/CpGs_UnivChicago_alldesigns_55860sites.cgn-pos-srd-prbs.tsv.gz")
    opt$org_des_tsv <- NULL
    
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
    make_option(c("--org_des_tsv"), type="character", default=opt$org_des_tsv, 
                help="Original design file used for canonical position selection. [default= %default]", metavar="character"),
    
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
    make_option(c("--mats"), type="character", default=opt$mats, 
                help="Match (format 1 or 2) file(s) (comma seperated) [default= %default]", metavar="character"),
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
opt_reqs <- c('outDir','impDir','ords','Species',
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

stamp_vec <- NULL
stamp_vec <- c(stamp_vec,opt$time_org_txt)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Parse List Options
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

image_key <- "bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse_Centos"
image_ver <- "v.1.0"
image_ver <- "v.1.11"
image_ssh <- "run_improbe.sh"
image_str <- glue::glue("{image_key}.{image_ver}")

# Manifest Control Defaults::
if (is.null(opt$ctls)) {
  opt$ctls <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
}

# ords_vec <- NULL
# mats_vec <- NULL
# aqps_vec <- NULL
# pqcs_vec <- NULL
# if (!is.null(opt$ords)) ords_vec <- opt$ords %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
# if (!is.null(opt$mats)) mats_vec <- opt$mats %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
# if (!is.null(opt$aqps)) aqps_vec <- opt$aqps %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
# if (!is.null(opt$pqcs)) pqcs_vec <- opt$pqcs %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

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
#                Pre-processing:: Run Time:: Ref Alignment Genome
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

gen_ref_dat <- NULL
run$gen_ref_fas <- 
  file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
            paste0(opt$genBuild,".genome.fa.gz"))

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Loading gen_ref_fas={run$gen_ref_fas}...{RET}"))

gen_ref_dat <- 
  Biostrings::readDNAStringSet(filepath = run$gen_ref_fas, format = "fasta") # , nrec = 2)

gen_ref_tab <- gen_ref_dat %>% names() %>% 
  stringr::str_remove(" .*$") %>% 
  stringr::str_remove("^chr") %>%
  tibble::tibble() %>% 
  purrr::set_names("Chrom_Char") %>% 
  dplyr::mutate(Idx=dplyr::row_number(),
                Chrom_Char=paste0("chr",Chrom_Char) )
print(gen_ref_tab)

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Loading gen_ref_fas={run$gen_ref_fas}{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Pre-processing:: Run Time:: SNP IUPAC Genome
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

gen_snp_dat <- NULL
run$gen_snp_fas <- 
  file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
            paste0(opt$genBuild,".dbSNP151-genome.fa.gz"))

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Loading gen_snp_fas={run$gen_snp_fas}...{RET}"))

if (!is.null(run$gen_snp_fas) && file.exists(run$gen_snp_fas)) {
  gen_snp_dat <- 
    Biostrings::readDNAStringSet(filepath = run$gen_snp_fas, format = "fasta") # , nrec = 2)
}

gen_snp_tab <- gen_snp_dat %>% names() %>% 
  stringr::str_remove(" .*$") %>% 
  stringr::str_remove("^chr") %>%
  tibble::tibble() %>% 
  purrr::set_names("Chrom_Char") %>% 
  dplyr::mutate(Idx=dplyr::row_number(),
                Chrom_Char=paste0("chr",Chrom_Char) )
print(gen_snp_tab)

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Loading gen_snp_fas={run$gen_snp_fas}{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Pre-processing:: Run Time:: Intermediate Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Defined Run Time:: Intermediate Files
# run$add_pas_csv <- file.path(run$addDir, paste(opt$runName,"address-pass.csv.gz", sep="."))
# run$man_fun_csv <- file.path(run$manDir, paste(opt$runName,"functional-sesame.manifest.csv.gz", sep="."))
# run$add_fin_csv <- file.path(run$manDir, paste(opt$runName,"final-address-aligned.csv.gz", sep="."))
# 
# run$man_mas_csv <- file.path(run$manDir, paste(opt$runName,"final-master.manifest.csv.gz", sep="."))
# run$man_ses_csv <- file.path(run$manDir, paste(opt$runName,"final-sesame.manifest.csv.gz", sep="."))
# run$man_gsm_csv <- file.path(run$manDir, paste(opt$runName,"final-GenomeStudio.manifest.csv.gz", sep="."))

run$add_prb_fas <- file.path(run$fasDir, paste(opt$runName, "aln-seq.fa.gz",  sep='.') )
run$add_dat_csv <- file.path(run$addDir, paste(opt$runName, "add_dat.csv.gz", sep='.') )
run$add_u49_tsv <- file.path(run$intDir, paste(opt$runName, "map-u49.tsv", sep='.') )
run$add_m49_tsv <- file.path(run$intDir, paste(opt$runName, "map-u50.tsv", sep='.') )

run$add_prb_bsp  <- file.path(run$alnDir, paste(opt$runName, "bsp",  sep='.') )
run$add_prb_bspz <- paste(run$add_prb_bsp, 'tsv.gz', sep='.')

#
# run$int_u49_tsv <- file.path(run$intDir, paste(opt$runName, "int-u49.tsv.gz", sep='.') )
# run$int_m49_tsv <- file.path(run$intDir, paste(opt$runName, "int-m49.tsv.gz", sep='.') )
# run$int_seq_tsv <- file.path(run$intDir, paste(opt$runName, "int-seq-imp.tsv.gz", sep='.') )
# 
# run$add_pas_bsp_csv <- file.path(run$alnDir, paste(opt$runName, "add_pas_bsp.csv.gz",  sep='.') )
# run$add_pas_grs_rds <- file.path(run$alnDir, paste(opt$runName, "add_pas_grs.rds",  sep='.') )

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Defining Run Time Files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: SNP/GSA Raw Manifests
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Loading Manifest Files...{RET}"))

man_450_csv <- file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz")
man_450_tib <- loadManifestGenomeStudio(
  file = man_450_csv, retType = "man", verbose = opt$verbose, tt=pTracker)

man_gsa_csv <- file.path(opt$manDir, "genotyping/GSA-24v2-0_A1.csv.gz")
raw_gsa_tib <- loadManifestGenomeStudio(
  file = man_gsa_csv, retType = "man", verbose = opt$verbose, tt=pTracker)

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Loading Manifest Files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: GSA Raw -> Manifests
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_gsa_tib <- raw_gsa_tib %>%
  dplyr::mutate(
    SNP_Top_Sequence=stringr::str_to_upper(TopGenomicSeq),
    UP_T=SNP_Top_Sequence %>%
      stringr::str_remove("\\[.*$"),
    QD_T=SNP_Top_Sequence %>% 
      stringr::str_remove("^.*\\[") %>%
      stringr::str_remove("\\].*$") %>%
      stringr::str_remove("\\/"),
    QA_T=QD_T %>% stringr::str_sub(1,1),
    QB_T=QD_T %>% stringr::str_sub(2,2),
    QT_T=dplyr::case_when(
      QD_T=="AG" | QD_T=="GA" ~ "R",
      QD_T=="CT" | QD_T=="TC" ~ "Y",
      
      QD_T=="CG" | QD_T=="GC" ~ "S",
      QD_T=="AT" | QD_T=="TA" ~ "W",
      
      QD_T=="GT" | QD_T=="TG" ~ "K",
      QD_T=="AC" | QD_T=="CA" ~ "M",
      
      TRUE ~ NA_character_
    ),
    # IlmnStrand/SourceStrand/RefStrand
    #
    # TOP/TOP/+   => Fwd (do na)
    # TOP/TOP/-   => Rev (do rc)
    #
    # TOP/BOT/+   => Fwd
    # TOP/BOT/-   => Rev
    #
    # BOT/TOP/+   => Rev
    # BOT/TOP/-   => Fwd
    #
    # BOT/BOT/+   => Rev
    # BOT/BOT/-   => Fwd
    #
    QI_T=dplyr::case_when(
      IlmnStrand=="TOP" & SourceStrand=="TOP" & RefStrand=="-" ~ cmpl(QT_T),
      IlmnStrand=="TOP" & SourceStrand=="BOT" & RefStrand=="-" ~ cmpl(QT_T),
      IlmnStrand=="BOT" & SourceStrand=="TOP" & RefStrand=="+" ~ cmpl(QT_T),
      IlmnStrand=="BOT" & SourceStrand=="BOT" & RefStrand=="+" ~ cmpl(QT_T),
      
      TRUE ~ QT_T
    ),
    DN_T=SNP_Top_Sequence %>% 
      stringr::str_remove("^.*\\]"),
    NB_T=DN_T %>%
      stringr::str_sub(1,1),
    DN_T=DN_T %>%
      stringr::str_sub(2)
  ) %>%
  dplyr::select(Name,Probe_Type,
                Chromosome,Coordinate,
                UP_T,
                QD_T,QA_T,QB_T,QI_T,
                NB_T,
                DN_T,
                dplyr::everything())

#
# TBD:: Mark Targets with CG[I]CG
#

man_gsa_tib %>% 
  dplyr::filter(Probe_Type=="rs") %>%
  dplyr::select(Name,Probe_Type,
                Chromosome,Coordinate,
                UP_T,
                QD_T,QA_T,QB_T,QI_T,
                NB_T,
                DN_T)

cur_gsa_tib <- man_gsa_tib %>% 
  dplyr::filter(Probe_Type=="rs")

# Target known SNP designs from 450k/EPIC::
#
# cur_gsa_tib <- man_gsa_tib %>% 
#   dplyr::filter(Probe_Type=="rs") %>% 
#   dplyr::filter(Name %in% man_450_tib$IlmnID)


#
# TBD:: Extract Forward_Sequence from coordinates
#

dna_gsa_fwd_tib <- 
  dna_to_template(tib = cur_gsa_tib, dna = gen_ref_dat, map = gen_ref_tab, iupac = "QI_T",
                  add_flank=TRUE,
                  verbose = opt$verbose, tt=pTracker)

dna_gsa_fwd_tib %>% dplyr::select(Name,Des_Din,Chromosome,Coordinate,Fwd_Seq:Des_Seq) %>% as.data.frame()

# Test Cases::
# dna_gsa_fwd_tib %>% dplyr::select(Name,Chromosome,Coordinate,Fwd_Seq:Des_Seq) %>% dplyr::filter(stringr::str_starts(Name,"rs1001737")) %>% as.data.frame()
# dna_gsa_fwd_tib %>% dplyr::select(Name,Chromosome,Coordinate,Fwd_Seq:Des_Seq) %>% dplyr::filter(stringr::str_starts(Name,"rs10159416")) %>% as.data.frame()

snp_gsa_fwd_tib <- 
  dna_to_template(tib = cur_gsa_tib, dna = gen_snp_dat, map = gen_snp_tab, iupac = "QI_T",
                  add_flank=TRUE,
                  verbose = opt$verbose, tt=pTracker)

snp_gsa_fwd_tib %>% dplyr::select(Name,Des_Din,Chromosome,Coordinate,Fwd_Seq:Des_Seq) %>% as.data.frame()

# Build Probes:: DNA
#
des_dna_gsa_list <- NULL
dna_gsa_fwd_list <- dna_gsa_fwd_tib %>% split(.$Des_Din)
for (ptype in names(dna_gsa_fwd_list)) {

  des_dna_gsa_list[[ptype]] <- 
    desSeq_to_prbs(
      dna_gsa_fwd_list[[ptype]],
      ids_key="Name",
      seq_key="Des_Seq",
      prb_key="Des_Din",
      strsSR="FR",
      strsCO="CO",
      addMatSeq=TRUE, parallel=TRUE, # max = 10,
      verbose=opt$verbose, tt=pTracker)
}

# Save Previous Work::
# des_dna_gsa_rds <- file.path(par$topDir, "data/manifests/genotyping/GSA-24v2-0_A1.probe-des-dna.v1.rds")
# readr::write_rds(des_dna_gsa_list, des_dna_gsa_rds, compress="gz")

# Add Top/Bot Calls::
#  -  Code below works, but not vlaidated...
#
# des_dna_gsa_cg_tb_tib <-
#   des_dna_gsa_list$cg %>%
#   dplyr::filter(PRB1_U_MAT!=PRB1_M_MAT) %>%
#   dplyr::distinct(Name, .keep_all=TRUE) %>%
#   setTopBot_tib(seqKey="DesSeqN", srdKey = "SRD_TB", topKey = "Strand_TB")

#
# Validate Previous 11 known SNPs:: DONE!!!
#
des_dna_gsa_list$rs %>%
  dplyr::inner_join(man_450_tib, by=c("Name"="IlmnID")) %>%
  dplyr::filter(
    PRB2_D_MAT==AlleleA_ProbeSeq |
      (PRB1_M_MAT==AlleleA_ProbeSeq & PRB1_U_MAT==AlleleB_ProbeSeq)) %>% 
  dplyr::select(Name,DesSeqN) %>% 
  dplyr::mutate(DesSeqN=addBrac(DesSeqN))


# Build Probes:: SNP
#
des_snp_gsa_list <- NULL
snp_gsa_fwd_list <- snp_gsa_fwd_tib %>% split(.$Des_Din)
for (ptype in names(snp_gsa_fwd_list)) {
  
  des_snp_gsa_list[[ptype]] <- 
    desSeq_to_prbs(
      snp_gsa_fwd_list[[ptype]],
      ids_key="Name",
      seq_key="Des_Seq",
      prb_key="Des_Din",
      strsSR="FR",
      strsCO="CO",
      addMatSeq=TRUE, parallel=TRUE, # max = 10,
      verbose=opt$verbose, tt=pTracker)
}

# Save Previous Work::
# des_snp_gsa_rds <- file.path(par$topDir, "data/manifests/genotyping/GSA-24v2-0_A1.probe-des-snp.v1.rds")
# readr::write_rds(des_snp_gsa_list, des_snp_gsa_rds, compress="gz")











# dna_gsa_fwd_tib %>% dplyr::select(Fwd_Seq)
# stringr::str_sub(dna_gsa_fwd_tib$Fwd_Seq, 62,62) <- dna_gsa_fwd_tib$QI_T; dna_gsa_fwd_tib %>% dplyr::select(Fwd_Seq)

#
# Convert SNP_Top_Sequence => Forward_Sequence
#
# TOP/TOP/+   => Fwd (do na)
# TOP/TOP/-   => Rev (do rc)
#
# TOP/BOT/+   => Fwd
# TOP/BOT/-   => Rev
#
# BOT/TOP/+   => Rev
# BOT/TOP/-   => Fwd
#
# BOT/BOT/+   => Rev
# BOT/BOT/-   => Fwd
#

dna_gsa_fwd_tib %>% head() %>% dplyr::select(Name,IlmnStrand,SourceStrand,RefStrand,SNP_Top_Sequence,SourceSeq,Fwd_Seq) %>% as.data.frame()

cur_gsa_tib %>% dplyr::group_by(IlmnStrand,SourceStrand,RefStrand,QD_T,QI_T) %>% dplyr::filter(!stringr::str_detect(QD_T, "-")) %>% dplyr::summarise(Count=n(), .groups="drop") %>% as.data.frame()


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             SNP Validation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Known Methylation SNPs:: Manifest
#  - NOTE:: Need to swap U/M Order...
#
swp_450_snp_tib <- man_450_tib %>% 
  dplyr::filter(Probe_Type=="rs") %>%
  clean_tibble()

man_450_snp_tib <- man_450_tib %>% 
  dplyr::filter(Probe_Type=="rs") %>%
  clean_tibble() %>%
  dplyr::mutate(
    TMP_ADD_U=AddressA_ID,
    TMP_PRB_U=AlleleA_ProbeSeq,
    
    # TMP_ADD_M=AddressB_ID,
    # TMP_PRB_M=AlleleB_ProbeSeq,
    
    AddressA_ID=dplyr::case_when(
      Infinium_Design_Type=="I" ~ AddressB_ID,
      TRUE ~ AddressA_ID
    ),
    AlleleA_ProbeSeq=dplyr::case_when(
      Infinium_Design_Type=="I" ~ AlleleB_ProbeSeq,
      TRUE ~ AlleleA_ProbeSeq
    ),
    
    AddressB_ID=dplyr::case_when(
      Infinium_Design_Type=="I" ~ TMP_ADD_U,
      TRUE ~ AddressB_ID
    ),
    AlleleB_ProbeSeq=dplyr::case_when(
      Infinium_Design_Type=="I" ~ TMP_PRB_U,
      TRUE ~ AlleleB_ProbeSeq
    )
  ) %>%
  dplyr::select(-TMP_ADD_U,-TMP_PRB_U)

# Viz Check Swap::
# swp_450_snp_tib %>% dplyr::select(1:6) %>% print()
# man_450_snp_tib %>% dplyr::select(1:6) %>% print()

#
# Known Methylation SNPs:: Address
#
add_450_snp_tib <- man_450_snp_tib %>% 
  man_to_add(verbose=opt$verbose, tt=pTracker)

fas_450_snp_tib <-
  add_to_fas(
    tib=add_450_snp_tib, 
    prb_key="Man_Prb",
    add_key="Man_Add", des_key="Man_Des", type_key="Man_Din",
    prb_fas=run$add_prb_fas, dat_csv=run$add_dat_csv,
    u49_tsv=run$add_u49_tsv, m49_tsv=run$add_m49_tsv,
    verbose=opt$verbose, tt=pTracker)

#
# Align SNPs to get genomic pos -> forward sequence
#
run$add_prb_bspz <- 
  run_bsmap(
    exe=opt$bsmap_exe, 
    fas=run$add_prb_fas, 
    gen=run$gen_ref_fas,
    bsp=run$add_prb_bsp,
    opt=NULL, lan=NULL, run=TRUE,
    verbose=opt$verbose,tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  3.3 Join Address and Alignment Data:: BSMAP
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsp_imp_tib <-
  join_bsmap(
    add=fas_450_snp_tib,
    bed=cgn_bed_tib, org=add_org_tib,
    file=run$add_prb_bspz,
    join_key="Aln_Key", 
    prb_des_key="Man_Des", prb_din_key="Man_Din",
    join_type="inner",
    sort=TRUE,
    verbose=opt$verbose,tt=pTracker)

#
# Join GSA/450
#
man_sel_snp_tib <- man_gsa_tib %>% 
  dplyr::inner_join(man_450_snp_tib, 
                    by=c("Name","Probe_Type"), 
                    suffix=c("_GSA","_450")) %>% 
  dplyr::distinct(Name, .keep_all = TRUE)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    3.7 Extract 122mer & SNP Probes::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Build Target List::
#   - Basically change chr/pos variable names...
#
tar_add_list <-
  bsp_imp_tib %>%
  dplyr::rename(Imp_Chr=Bsp_Chr,Imp_Pos=Bsp_Pos) %>%
  dplyr::arrange(Imp_Chr,Imp_Pos) %>%
  split(.$Imp_Chr)

fwd_des_tib <- NULL
fet_chroms <- names(tar_add_list)
for (chr_str in fet_chroms) {
  cat(glue::glue("[{par$prgmTag}]: chr={chr_str}.{RET}"))
  
  if (is.null(tar_add_list[[chr_str]])) next
  
  chr_idx <- gen_snp_tab %>% 
    dplyr::filter(Chrom_Char==chr_str) %>% 
    head(n=1) %>% pull(Idx) %>% as.integer()
  # print(chr_idx)
  
  cur_add_tib <- tar_add_list[[chr_str]]
  bsp_begs <- cur_add_tib$Imp_Pos - 60
  bsp_ends <- bsp_begs + 122 - 1
  # print(bsp_begs)
  
  ref_seqs <- stringr::str_sub( as.character(gen_ref_dat[[chr_idx]]), bsp_begs, bsp_ends) %>% addBrac()
  snp_seqs <- stringr::str_sub( as.character(gen_snp_dat[[chr_idx]]), bsp_begs, bsp_ends) %>% addBrac()
  # print(snp_seqs)
  
  cur_add_tib <- 
    cur_add_tib %>% 
    dplyr::mutate(
      Gen_Ref_Seq=ref_seqs,
      Gen_Snp_Seq=snp_seqs
    )
  
  #
  # REF::
  #
  des_ref_tib <- 
    desSeq_to_prbs(
      cur_add_tib,
      ids_key="Aln_Key",
      seq_key="Gen_Ref_Seq",
      prb_key="Man_Din",
      strsSR="FR",
      strsCO="CO",
      addMatSeq=TRUE, parallel=TRUE, # max = 10,
      verbose=opt$verbose, tt=pTracker)
  
  ref_join_tib <- NULL
  ref_join_tib <- 
    cur_add_tib %>% 
    dplyr::inner_join(des_ref_tib, 
                      by=c("Aln_Key","Man_Din") ) %>% 
    dplyr::mutate(Match=dplyr::case_when(
      Man_Des=="U" & Man_Prb==PRB1_U_MAT ~ 0,
      Man_Des=="M" & Man_Prb==PRB1_M_MAT ~ 0,
      Man_Des=="2" & Man_Prb==PRB2_D_MAT ~ 0,
      TRUE ~ 1),
      Fas_Src="Ref"
    ) %>% dplyr::filter(Match==0)
  
  #
  # SNP::
  #
  des_snp_tib <- 
    desSeq_to_prbs(
      cur_add_tib,
      ids_key="Aln_Key",
      seq_key="Gen_Snp_Seq",
      prb_key="Man_Din",
      strsSR="FR",
      strsCO="CO",
      addMatSeq=TRUE, parallel=TRUE, # max = 10,
      verbose=opt$verbose, tt=pTracker)
  
  snp_join_tib <- NULL
  snp_join_tib <- 
    cur_add_tib %>% 
    dplyr::inner_join(des_snp_tib, 
                      by=c("Aln_Key","Man_Din") ) %>% 
    dplyr::mutate(Match=dplyr::case_when(
      Man_Des=="U" & Man_Prb==PRB1_U_MAT ~ 0,
      Man_Des=="M" & Man_Prb==PRB1_M_MAT ~ 0,
      Man_Des=="2" & Man_Prb==PRB2_D_MAT ~ 0,
      TRUE ~ 1),
      Fas_Src="Snp"
    ) %>% dplyr::filter(Match==0)
  
  if (FALSE) {
    #
    # Run improbe full design for thermodynamic scores::
    #
    opt$desDir <- file.path(opt$outDir, "des")
    if (!dir.exists(opt$desDir)) dir.create(opt$desDir, recursive = TRUE)
    ref_des_tsv <- file.path(opt$desDir, paste(opt$genBuild,chr_str,"improbe-ref-design-input.tsv", sep='-'))
    
    ref_des_tib <- tibble::tibble(
      Seq_ID=cur_des_tib$Des_Key,
      Sequence=cur_des_tib$Des_Ref_Seq,
      Genome_Build=cur_des_tib$Des_Bld,
      Chromosome=cur_des_tib$Des_Chr,
      Coordinate=cur_des_tib$Des_Pos,
      CpG_Island="FALSE"
    )
    readr::write_tsv(ref_des_tib,ref_des_tsv)
    
    # Not sure why this isn't running...
    ref_imp_des_tib <- 
      improbe_docker(dir = opt$desDir, file = ref_des_tsv, 
                     name = paste(opt$genBuild,chr_str, sep="-"), 
                     image = image_str, shell = image_ssh,
                     verbose = opt$verbose, tt=pTracker)
  }
  
  fwd_des_tib <- 
    bind_rows(fwd_des_tib, 
              ref_join_tib,
              snp_join_tib)
  # print(fwd_des_tib)
  
  cat(glue::glue("[{par$prgmTag}]: Done.{RET}{RET}"))
  
  # break
}
# fwd_des_tib %>% dplyr::filter(Match==0)






# TBD:: Might need to revCmp()...
#
# TBD:: Compare the GSA vs. Ref/Snp 122mers.
#      - Add the IUPAC Code in the []..
#
#
# fwd_des_tib$Man_PID
# man_sel_snp_tib

# fwd_des_tib %>% dplyr::inner_join(man_sel_snp_tib, by=c("Man_PID"="Name"))
tmp_des_tib <- fwd_des_tib %>% 
  dplyr::mutate(
    DesSeqN=stringr::str_to_upper(DesSeqN) %>% addBrac(),
    GEN_Pre_Seq=stringr::str_remove(DesSeqN, "\\[.*$"),
    GEN_Din_Seq=DesSeqN %>% 
      stringr::str_remove("^.*\\[") %>%
      stringr::str_remove("\\].*$"),
    GEN_Sin_Nuc=stringr::str_sub(GEN_Din_Seq,1,1),
    GEN_Pos_Nuc=stringr::str_sub(GEN_Din_Seq,2,1),
    GEN_Pos_Seq=DesSeqN %>% 
      stringr::str_remove("^.*\\]") %>%
      stringr::str_sub(2)
  ) %>%
  dplyr::inner_join(man_sel_snp_tib, by=c("Man_PID"="Name")) %>% 
  dplyr::select(Man_PID,Man_Des,Bsp_Srd,CPN_D,
                Imp_Chr,Imp_Pos,Fas_Src,
                DesSeqN,SNP_Top_Sequence) %>% 
  dplyr::mutate(DesSeqN=addBrac(DesSeqN))
  # dplyr::filter(Fas_Src=="Snp")

tmp_des_tib %>% as.data.frame() %>% print()

# TBD:: Need to compare pre/mid/pos sequence with IUPAC
#















ret_list <- NULL
ref_list <- des_ref_tib %>% 
  dplyr::left_join(cur_add_tib, by=c("Aln_Key","Man_Din")) %>% 
  dplyr::mutate(
    Match=dplyr::case_when(
      Man_Des=="U" & Man_Prb==PRB1_U_MAT ~ 0,
      Man_Des=="M" & Man_Prb==PRB1_M_MAT ~ 0,
      Man_Des=="2" & Man_Prb==PRB2_D_MAT ~ 0,
      TRUE ~ 1),
    Mat_Src="ref"
  ) %>%
  dplyr::arrange(Man_PID,Strand_SR,Strand_CO) %>%
  split(.[["Man_Des"]])

# ref_list[["U"]] %>% dplyr::select(Man_PID,Strand_SR,Strand_CO,Aln_Key,Man_Din,Man_Des,Man_Din, Man_Prb,PRB1_U_MAT, Match)
# ref_list[["M"]] %>% dplyr::select(Man_PID,Strand_SR,Strand_CO,Aln_Key,Man_Din,Man_Des,Man_Din, Man_Prb,PRB1_M_MAT, Match)

# TBD:: Need to Swap in the SNP!!!!
ref_list[["U"]] %>% dplyr::select(Man_PID,Strand_SR,Strand_CO,Aln_Key,Man_Din,Man_Des,Man_Din, Man_Prb,PRB1_U_MAT, Match)
ref_list[["M"]] %>% dplyr::select(Man_PID,Strand_SR,Strand_CO,Aln_Key,Man_Din,Man_Des,Man_Din, Man_Prb,PRB1_M_MAT, Match)
ref_list[["2"]] %>% dplyr::select(Man_PID,Strand_SR,Strand_CO,Aln_Key,Man_Din,Man_Des,Man_Din, Man_Prb,PRB2_D_MAT, Match)

#   Man_PID   Strand_SR Strand_CO Aln_Key       Man_Din Man_Des Man_Prb                                            PRB1_*_MAT                                         Match
# 5 rs3936238 F         C         42802473_U_rs rs      U       AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCC AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCT     1
# 5 rs3936238 F         C         29627504_M_rs rs      M       AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCT AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCT     0
#
#
# 5 rs3936238 F         C         29627504_U_rs rs      U       AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCT AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCT     0
# 5 rs3936238 F         C         42802473_M_rs rs      M       AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCC AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCT     1

snp_list <- NULL
snp_list <- des_snp_tib %>% 
  dplyr::left_join(cur_add_tib, by=c("Aln_Key","Man_Din")) %>% 
  dplyr::mutate(
    Match=dplyr::case_when(
      Man_Des=="U" & Man_Prb==PRB1_U_MAT ~ 0,
      Man_Des=="M" & Man_Prb==PRB1_M_MAT ~ 0,
      Man_Des=="2" & Man_Prb==PRB2_D_MAT ~ 0,
      TRUE ~ 1),
    Mat_Src="snp"
  ) %>%
  dplyr::arrange(Man_PID,Strand_SR,Strand_CO) %>%
  split(.[["Man_Des"]])

snp_list[["U"]] %>% dplyr::select(Man_PID,Strand_SR,Strand_CO,Aln_Key,Man_Din,Man_Des,Man_Din, Man_Prb,PRB1_U_MAT, Match)
snp_list[["M"]] %>% dplyr::select(Man_PID,Strand_SR,Strand_CO,Aln_Key,Man_Din,Man_Des,Man_Din, Man_Prb,PRB1_M_MAT, Match)
snp_list[["2"]] %>% dplyr::select(Man_PID,Strand_SR,Strand_CO,Aln_Key,Man_Din,Man_Des,Man_Din, Man_Prb,PRB2_D_MAT, Match)

#   Man_PID   Strand_SR Strand_CO Aln_Key       Man_Din Man_Des Man_Prb                                            PRB1_M_MAT                                         Match
# 5 rs3936238 F         C         42802473_U_rs rs      U       AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCC AATACCCTACAATATACAAAAACAACAMTTCCTATAYAAAACTACTACCT     1
# 5 rs3936238 F         C         29627504_M_rs rs      M       AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCT AATACCCTACAATATACAAAAACAACAMTTCCTATAYAAAACTACTACCC     1
#
# 5 rs3936238 F         C         29627504_U_rs rs      U       AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCT AATACCCTACAATATACAAAAACAACAMTTCCTATAYAAAACTACTACCT     1
# 5 rs3936238 F         C         42802473_M_rs rs      M       AATACCCTACAATATACAAAAACAACACTTCCTATACAAAACTACTACCC AATACCCTACAATATACAAAAACAACAMTTCCTATAYAAAACTACTACCC     1

mat_add_tib <-
  dplyr::bind_rows(
    dplyr::bind_rows(ref_list),
    dplyr::bind_rows(snp_list)
  ) %>% 
  dplyr::filter(Match==0)

mat_add_tib %>% 
  dplyr::group_by(Man_PID,Aln_Key,Mat_Src) %>% 
  dplyr::summarise(Count=n(), .groups="drop") %>%
  dplyr::arrange(Man_PID)



























# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     4.0 Write Sesame & Genome Studio Manifest::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  fwd_des_csv <- file.path(run$alnDir, paste(opt$runName,"fwd-seqs.csv.gz", sep='-'))
  readr::write_csv(fwd_des_tib, fwd_des_csv)
  
  
  fwd_sel_tib <- fwd_des_tib %>% 
    dplyr::rename(Forward_Sequence=Des_Ref_Seq.x) %>%
    dplyr::mutate(Source_Seq=dplyr::case_when(
      Des_Inf==1 ~ PRB1_N,
      Des_Inf==2 ~ PRB2_N,
      TRUE ~ NA_character_) 
    ) %>%
    dplyr::select(Des_Key,Des_AddU,Des_AddM, Forward_Sequence,Source_Seq, Strand_SR,Strand_CO)
  
  man_fwd_ses_col <- c("Probe_ID","Name","U","AlleleA_ProbeSeq","M","AlleleB_ProbeSeq",
                       "Next_Base","Color_Channel","col","Probe_Type",
                       "Strand_FR","Strand_TB","Strand_CO","Infinium_Design","Infinium_Design_Type",
                       "CHR","MAPINFO","Species","Genome_Build","Source_Seq","Forward_Sequence","Rep_Num")
  man_fwd_ses_tib <- man_top_sesame_tib %>%
    dplyr::left_join(fwd_sel_tib, 
                     by=c("Probe_ID"="Des_Key", 
                          "U"="Des_AddU", 
                          "M"="Des_AddM", "Strand_CO") ) %>%
    dplyr::mutate(Name=stringr::str_remove(Probe_ID, "_.*$")) %>%
    dplyr::select(dplyr::all_of(man_fwd_ses_col))
  
  #
  # Sesame Controls::
  #
  man_ses_ctl_tib <- 
    readr::read_csv(file.path(par$datDir, "manifest/core/COVIC-C12.manifest.sesame-base.cpg-sorted.csv.gz")) %>%
    dplyr::select(dplyr::any_of(names(man_fwd_ses_tib))) %>% 
    dplyr::filter(stringr::str_starts(Probe_ID, "ct")) %>%
    dplyr::rename(Infinium_Design_Type=Infinium_Design) %>%
    dplyr::mutate(Infinium_Design=as.integer(2))
  
  ses_full_tib <- dplyr::bind_rows(
    man_fwd_ses_tib,man_ses_ctl_tib)
  readr::write_csv(ses_full_tib, run$man_ses_csv)
  
  
  #
  # Genome Studio::
  #
  gs_header_vec <- genome_studio_header(name="Chicago-A5-GRCh37",count=base::nrow(man_fwd_ses_tib))
  gs_controls_tib <- readr::read_csv("/Users/bretbarnes/Documents/data/manifests/raw/manifests/methylation/controls/MethylationEPIC_v-1-0_B4.only-controls-table.csv.gz")
  gs_body_tib <- man_fwd_ses_tib %>% dplyr::rename(
    IlmnID=Probe_ID,
    AddressA_ID=U,
    AddressB_ID=M
  )
  
  readr::write_lines(gs_header_vec, run$man_gsm_csv, append = FALSE)
  readr::write_csv(gs_body_tib, run$man_gsm_csv, append = TRUE, col_names = TRUE)
  readr::write_lines("[Controls],,,,,,,", run$man_gsm_csv, append = TRUE)
  readr::write_csv(gs_controls_tib, run$man_gsm_csv, append = TRUE, col_names = FALSE)
  
  
  if (FALSE) {
    seq_imp_tib %>% dplyr::select(Address, Can_Top) %>% dplyr::filter(!is.na(Can_Top)) %>%
      dplyr::inner_join(exp_tib, by=c("Address"="U"))
    
    # Some quick validation::
    seq_imp_tib %>% dplyr::select(Address, Can_Top) %>% dplyr::filter(!is.na(Can_Top)) %>%
      dplyr::inner_join(exp_tib, by=c("Address"="U")) %>% dplyr::select(Des_Key, Can_Top, Des_Ref_Seq.y) %>% dplyr::mutate(Des_Ref_Seq1=revCmp( shearBrac(Des_Ref_Seq.y))) %>% dplyr::select(Can_Top,Des_Ref_Seq1)
    
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file

