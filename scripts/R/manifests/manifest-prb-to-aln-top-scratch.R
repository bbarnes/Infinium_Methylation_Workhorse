
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# This builds the offical Canonical CGN->TOP->GRP prereq for cgnDB creation!!!
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
par$prgmDir <- 'manifests'
par$prgmTag <- 'manifest-prb-to-aln-top'
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
  par$local_runType <- 'All'

  if (par$local_runType=='All') {
    # opt$genBuild <- paste('GRCh38','GRCh37','GRCm10', sep=',')
    opt$genBuild <- 'GRCh37'
    
    opt$platform <- 'EPIC'
    opt$version  <- 'A1'

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
opt_reqs <- c('outDir','impDir','ords',
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
image_ssh <- "run_improbe.sh"
image_str <- glue::glue("{image_key}.{image_ver}")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Pre-processing:: Run Time:: Output Directories
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run <- list()

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


# Manifest Control Defaults::
if (is.null(opt$ctls)) {
  opt$ctls <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Pre-processing:: Run Time:: Intermediate Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Define Run Time:: Ref Alignment Genome
run$gen_ref_fas <- 
  file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
            paste0(opt$genBuild,".genome.fa.gz"))

# Define Run Time:: SNP IUPAC Genome
run$gen_snp_fas <- 
  file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
            paste0(opt$genBuild,".dbSNP151-genome.fa.gz"))

# Defined Run Time:: Intermediate Files
run$add_pas_csv <- file.path(run$addDir, paste(opt$runName,"address-pass.csv.gz", sep="."))
run$man_fun_csv <- file.path(run$manDir, paste(opt$runName,"functional-sesame.manifest.csv.gz", sep="."))
run$man_ses_csv <- file.path(run$manDir, paste(opt$runName,"final-sesame.manifest.csv.gz", sep="."))
run$add_fin_csv <- file.path(run$manDir, paste(opt$runName,"final-address-aligned.csv.gz", sep="."))

run$add_prb_fas <- file.path(run$fasDir, paste(opt$runName, "aln-seq.fa.gz",  sep='.') )
run$add_dat_csv <- file.path(run$addDir, paste(opt$runName, "add_dat.csv.gz", sep='.') )

run$add_u49_tsv <- file.path(run$intDir, paste(opt$runName, "map-u49.tsv", sep='.') )
run$add_m49_tsv <- file.path(run$intDir, paste(opt$runName, "map-u50.tsv", sep='.') )

run$int_u49_tsv <- file.path(run$intDir, paste(opt$runName, "int-u49.tsv.gz", sep='.') )
run$int_m49_tsv <- file.path(run$intDir, paste(opt$runName, "int-m49.tsv.gz", sep='.') )
run$int_seq_tsv <- file.path(run$intDir, paste(opt$runName, "int-seq-imp.tsv.gz", sep='.') )

run$add_prb_bsp  <- file.path(run$alnDir, paste(opt$runName, "bsp",  sep='.') )
run$add_prb_bspz <- paste(run$add_prb_bsp, 'tsv.gz', sep='.')

run$add_pas_bsp_csv <- file.path(run$alnDir, paste(opt$runName, "add_pas_bsp.csv.gz",  sep='.') )
run$add_pas_grs_rds <- file.path(run$alnDir, paste(opt$runName, "add_pas_grs.rds",  sep='.') )

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Defining Run Time Files.{RET}{RET}"))

#
#
# Load All Manifest GSM/Ses and Full 
#
#
gsm_files <- list.files( file.path(opt$manDir, "GenomeStudio"), pattern=".csv.gz", full.names=TRUE)
gsm_ranks <- tibble::tibble(
  Man_Source=c(
    "HumanMethylation450_15017482_v.1.2",
    "MethylationEPIC_v-1-0_B2",
    "HumanMethylation27_270596_v.1.2",
    "Mpanel1.LEGX.CombinedManifestEPIC.manifest",
    "HorvathMammal40-A2"
  ),
  Man_Rank=c(1,2,3,4,5) %>% as.integer()
)

# We want unique(probe, din, des)
#  Aln_Seq=strip-degen-bases...
#   Make this look like a manifest
#

gsm_prb_tib <- lapply(gsm_files, loadManifestGenomeStudio,
       addSource=TRUE,normalize=TRUE,retType="man", 
       verbose=opt$verbose+10,tt=pTracker) %>% 
  dplyr::bind_rows()

gsm_sel_tib <- gsm_ranks %>% 
  dplyr::right_join(gsm_prb_tib, by="Man_Source") %>% 
  dplyr::arrange(Man_Rank) %>% 
  dplyr::distinct(AlleleA_ProbeSeq,AlleleB_ProbeSeq,Probe_Type, .keep_all=TRUE)

gsm_sel_tib %>%
  dplyr::filter(Infinium_Design=="I" )
  dplyr::group_by(Probe_Type) %>% 
  dplyr::summarise(Count=n(), .groups="drop")





#
# Load All Human Control Sequences
#  data/Infinium_Methylation_Workhorse/dat/manifest/controls/Infinium_Methylation_Controls_1983_full.csv.gz
#

hum_ctl_csv <- file.path(par$topDir, "data/Infinium_Methylation_Workhorse/dat/manifest/controls/Infinium_Methylation_Controls_1983_full.csv.gz")
hum_ctl_tib <- suppressMessages(suppressWarnings( readr::read_csv(hum_ctl_csv) )) %>%
  dplyr::filter(!is.na(Probe_Seq)) %>%
  dplyr::mutate(Address=as.integer(U),Ord_Din=Probe_Type,Ord_Des="2",
                Ord_Prb=stringr::str_to_upper(Probe_Seq),
                Aln_Prb=Ord_Prb,
                Aln_Rev=revCmp(Aln_Prb),
                Aln_Key=paste(Address,Ord_Din,Ord_Des, sep="_"),
                fas_line=paste0(">",Aln_Key,"\n",Ord_Prb)
  )

readr::write_lines(hum_ctl_tib$fas_line, run$add_prb_fas, append = FALSE)

# Split genBuild
# Add full genome path
# Run alignment
# Extract 122-mers
# 

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    3.2 Align All Probe Sequence:: BSMAP
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

stamp_vec <- c(stamp_vec, run$add_prb_bspz)
if (opt$fresh || !valid_time_stamp(stamp_vec)) {
  
  run$add_prb_bspz <- 
    run_bsmap(
      exe=opt$bsmap_exe, 
      fas=run$add_prb_fas, 
      gen=run$gen_ref_fas,
      bsp=run$add_prb_bsp,
      opt=NULL, lan=NULL, run=TRUE,
      verbose=opt$verbose,tt=pTracker)
  
} else {
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Using existing file add_prb_bspz={run$add_prb_bspz}...{RET}"))
}

# load_bsmap(run$add_prb_bspz)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  3.3 Join Address and Alignment Data:: BSMAP
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsp_imp_tib <- NULL
stamp_vec <- c(stamp_vec, run$add_pas_bsp_csv)
if (opt$fresh || !valid_time_stamp(stamp_vec)) {
  
  bsp_imp_tib <-
    join_bsmap(
      add=hum_ctl_tib, # bed=cgn_bed_tib,
      file=run$add_prb_bspz,
      join_key="Aln_Key", prb_des="Ord_Des",
      join_type="inner",
      sort=TRUE,
      verbose=opt$verbose,tt=pTracker)
  
  safe_write(bsp_imp_tib,"csv",run$add_pas_bsp_csv, funcTag=par$prgmTag, 
             verbose=opt$verbose)
  
} else {
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Loading add_pas_bsp_csv={run$add_pas_bsp_csv}...{RET}"))
  
  bsp_imp_tib <- suppressMessages(suppressWarnings(
    readr::read_csv(run$add_pas_bsp_csv, guess_max=100000) )) %>%
    utils::type.convert() %>% 
    dplyr::mutate(across(where(is.factor),  as.character) )
  
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    3.7 Extract 122mer & SNP Probes::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Define Run Time:: Ref Alignment Genome
gen_ref_dat <- NULL
run$gen_ref_fas <- 
  file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
            paste0(opt$genBuild,".genome.fa.gz"))
gen_ref_dat <- 
  Biostrings::readDNAStringSet(filepath = run$gen_ref_fas, format = "fasta") # , nrec = 2)


# Define Run Time:: SNP IUPAC Genome
gen_snp_dat <- NULL
run$gen_snp_fas <- 
  file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
            paste0(opt$genBuild,".dbSNP151-genome.fa.gz"))
if (!is.null(run$gen_snp_fas) && file.exists(run$gen_snp_fas)) {
  gen_snp_dat <- 
    Biostrings::readDNAStringSet(filepath = run$gen_snp_fas, format = "fasta") # , nrec = 2)
}

gen_ref_tab <- gen_ref_dat %>% names() %>% 
  stringr::str_remove(" .*$") %>% 
  stringr::str_remove("^chr") %>%
  tibble::tibble() %>% 
  purrr::set_names("Chrom_Char") %>% 
  dplyr::mutate(Idx=dplyr::row_number(),
                Chrom_Char=paste0("chr",Chrom_Char) )
print(gen_ref_tab)

gen_snp_tab <- gen_snp_dat %>% names() %>% 
  stringr::str_remove(" .*$") %>% 
  stringr::str_remove("^chr") %>%
  tibble::tibble() %>% 
  purrr::set_names("Chrom_Char") %>% 
  dplyr::mutate(Idx=dplyr::row_number(),
                Chrom_Char=paste0("chr",Chrom_Char) )
print(gen_snp_tab)

tar_add_list <- 
  bsp_imp_tib %>% 
  dplyr::distinct(Address,Ord_Des,Ord_Din, Bsp_Chr,Bsp_Pos) %>% 
  dplyr::arrange(Bsp_Chr,Bsp_Pos) %>%
  split(.$Bsp_Chr)

fwd_des_tib <- NULL
fet_chroms <- names(tar_add_list)
for (chr_str in fet_chroms) {
  cat(glue::glue("[{par$prgmTag}]: chr={chr_str}.{RET}"))
  
  if (is.null(tar_add_list[[chr_str]])) next
  
  chr_idx <- gen_snp_tab %>% 
    dplyr::filter(Chrom_Char==chr_str) %>% 
    head(n=1) %>% pull(Idx) %>% as.integer()
  print(chr_idx)
  
  bsp_begs <- tar_add_list[[chr_str]]$Bsp_Pos - 60
  bsp_ends <- bsp_begs + 122
  # print(bsp_begs)
  
  ref_seqs <- stringr::str_sub( as.character(gen_ref_dat[[chr_idx]]), bsp_begs, bsp_ends) %>% addBrac()
  snp_seqs <- stringr::str_sub( as.character(gen_snp_dat[[chr_idx]]), bsp_begs, bsp_ends) %>% addBrac()
  # print(snp_seqs)
  
  #
  # TBD:: Should determine Des_Din based on Des_Din + [XX]
  #
  
  cur_des_tib <- 
    tibble::tibble(
      Des_Key=paste(tar_add_list[[chr_str]]$Address,
                    tar_add_list[[chr_str]]$Ord_Des,
                    tar_add_list[[chr_str]]$Ord_Din, sep="_"),
      Des_Bld=opt$genBuild,
      Des_Chr=chr_str,
      Des_Pos=tar_add_list[[chr_str]]$Imp_Pos,
      Des_Din=tar_add_list[[chr_str]]$Ord_Din,
      Des_FR="F",
      Des_CO="C",
      
      Des_Ref_Seq=ref_seqs,
      Des_Snp_Seq=snp_seqs
    ) %>% 
    dplyr::mutate(
      Des_Chr=as.character(Des_Chr),
      Des_Pos=as.integer(Des_Pos)
    )
  
  # TBD:: Should Split based on Des_Din (cg, ch, rs, etc...)
  #
  des_ref_tib <- 
    desSeq_to_prbs(cur_des_tib, 
                   idsKey = "Des_Key", seqKey = "Des_Ref_Seq", prbKey = "Des_Din", # prbKey = "cg",
                   strsSR = "Des_FR",strsCO = "Des_CO", 
                   addMatSeq = TRUE, parallel = TRUE, max = 10,
                   verbose = opt$verbose, tt = pTracker)
  
  des_snp_tib <- 
    desSeq_to_prbs(cur_des_tib, 
                   idsKey = "Des_Key", seqKey = "Des_Snp_Seq", prbKey = "Des_Din", # prbKey = "cg",
                   strsSR = "Des_FR",strsCO = "Des_CO", 
                   addMatSeq = TRUE, parallel = TRUE, max = 10,
                   verbose = opt$verbose, tt = pTracker)
  
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
  
  
  
  fwd_des_tib <- bind_rows(fwd_des_tib, cur_des_tib)
  # print(fwd_des_tib)
  
  cat(glue::glue("[{par$prgmTag}]: Done.{RET}{RET}"))
  
  break
}
















if (FALSE) {
  ords_vec <- NULL
  mats_vec <- NULL
  aqps_vec <- NULL
  pqcs_vec <- NULL
  if (!is.null(opt$ords)) ords_vec <- opt$ords %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  if (!is.null(opt$mats)) mats_vec <- opt$mats %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  if (!is.null(opt$aqps)) aqps_vec <- opt$aqps %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  if (!is.null(opt$pqcs)) pqcs_vec <- opt$pqcs %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  
  ords_len <- length(ords_vec)
  mats_len <- length(mats_vec)
  stopifnot(ords_len>0)
  stopifnot(mats_len>0)
  stopifnot(ords_len==mats_len)
  
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
  #                Pre-processing:: Run Time:: Intermediate Files
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Define Run Time:: Alignment Genome
  run$man_gen_fas <- file.path(opt$genDir, opt$genBuild, "Sequence/WholeGenomeFasta",
                               paste0(opt$genBuild,".genome.fa.gz"))
  
  # Defined Run Time:: Intermediate Files
  run$add_pas_csv <- file.path(run$addDir, paste(opt$runName,"address-pass.csv.gz", sep="_"))
  run$man_fun_csv <- file.path(run$manDir, paste(opt$runName,"functional-sesame.manifest.csv.gz", sep="_"))
  
  run$add_prb_fas <- file.path(run$fasDir, paste(opt$runName, "aln-seq.fa.gz",  sep='.') )
  run$add_dat_csv <- file.path(run$addDir, paste(opt$runName, "add_dat.csv.gz", sep='.') )
  
  run$add_u49_tsv <- file.path(run$intDir, paste(opt$runName, "map-u49.tsv", sep='.') )
  run$add_u50_tsv <- file.path(run$intDir, paste(opt$runName, "map-u50.tsv", sep='.') )
  
  run$int_u49_tsv <- file.path(run$intDir, paste(opt$runName, "int-u49.tsv.gz", sep='.') )
  run$int_u50_tsv <- file.path(run$intDir, paste(opt$runName, "int-u50.tsv.gz", sep='.') )
  run$int_seq_tsv <- file.path(run$intDir, paste(opt$runName, "int-seq-imp.tsv.gz", sep='.') )
  
  run$add_prb_bsp  <- file.path(run$alnDir, paste(opt$runName, "bsp",  sep='.') )
  run$add_prb_bspz <- paste(run$add_prb_bsp, 'tsv.gz', sep='.')
  
  run$add_pas_bsp_csv <- file.path(run$alnDir, paste(opt$runName, "add_pas_bsp.csv.gz",  sep='.') )
  run$add_pas_grs_rds <- file.path(run$alnDir, paste(opt$runName, "add_pas_grs.rds",  sep='.') )
  
  # run$add_join_bsp_csv <- file.path(run$alnDir, paste(opt$runName, "add-join-bsp.csv.gz",  sep='.') )
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Done. Defining Run Time Files.{RET}{RET}"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #                       1.0 Validate Pre-built:: CGN -> TOP
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  pre_cgn_top_csv  <- "/Users/bretbarnes/Documents/data/improbe/dbCGN/pre-assigned/pre-assigned.cgnTop.hash.csv.gz"
  mm10_ses_man_csv <- file.path(opt$manDir, "LEGX-C24A.manifest.sesame-base.cpg-sorted.csv.gz")
  
  pre_cgn_top_col <-
    cols(
      CGN  = col_character(),
      TOP  = col_character(),
      SRC  = col_character()
    )
  
  pre_cgn_top_tib  <- 
    readr::read_csv(pre_cgn_top_csv,
                    # col_names=names(pre_cgn_top_col$cols), 
                    col_types=pre_cgn_top_col)  %>% 
    dplyr::filter(stringr::str_starts(CGN, "[0-9]")) %>%
    utils::type.convert() %>% 
    dplyr::mutate(across(where(is.factor),  as.character) )
  
  pre_cgn_top_sum <- pre_cgn_top_tib %>%
    dplyr::group_by(SRC) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  pre_cgn_top_sum %>% 
    dplyr::arrange(-Count) %>%
    print(n=base::nrow(pre_cgn_top_sum))
  
  mm10_ses_man_tib <- readr::read_csv(mm10_ses_man_csv)
  
  mm10_unq_man_tib <- mm10_ses_man_tib %>% 
    dplyr::select(Probe_ID,Probe_Type,Probe_Class,Top_Sequence) %>% 
    dplyr::filter(!is.na(Top_Sequence)) %>% 
    dplyr::filter(Probe_Type=='cg') %>%
    tidyr::separate(Probe_ID, into=c("CGN","SRD_Str"), sep="_", remove=TRUE) %>%
    dplyr::mutate(
      Top_Sequence=shearBrac(Top_Sequence),
      CGN=stringr::str_remove(CGN,"^cg") %>% stringr::str_remove("^0+") %>% as.integer()
    ) %>%
    dplyr::distinct(CGN,Probe_Type,Probe_Class,Top_Sequence) %>%
    dplyr::add_count(CGN, name="CGN_Count") %>%
    dplyr::add_count(CGN,Top_Sequence, name="CGN_TOP_Count")
  
  mm10_all_man_tib <- mm10_ses_man_tib %>% 
    dplyr::select(Probe_ID,Probe_Type,Probe_Class,Top_Sequence) %>% 
    dplyr::filter(!is.na(Top_Sequence)) %>% 
    dplyr::filter(Probe_Type=='cg') %>%
    tidyr::separate(Probe_ID, into=c("CGN","SRD_Str"), sep="_", remove=TRUE) %>%
    dplyr::mutate(
      Top_Sequence=shearBrac(Top_Sequence),
      CGN=stringr::str_remove(CGN,"^cg") %>% stringr::str_remove("^0+") %>% as.integer()
    ) %>%
    # dplyr::distinct(CGN,Probe_Type,Probe_Class,Top_Sequence) %>%
    dplyr::add_count(CGN, name="CGN_Count") %>%
    dplyr::add_count(CGN,Top_Sequence, name="CGN_TOP_Count")
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Compare:: Pre vs. mm10
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  mm10_join_tib <- dplyr::inner_join(pre_cgn_top_tib, mm10_unq_man_tib, by="CGN")
  mm10_anti_tib <- dplyr::anti_join(pre_cgn_top_tib, mm10_unq_man_tib, by="CGN")
  
  # This is good::
  mm10_join_tib %>% dplyr::filter(TOP != Top_Sequence)
  
  # This is good::
  dplyr::anti_join(mm10_unq_man_tib, pre_cgn_top_tib, by="CGN")
  
  mm10_join_sum <- mm10_join_tib %>%
    dplyr::group_by(SRC,Probe_Type,Probe_Class) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  mm10_join_sum %>% print(n=base::nrow(mm10_join_sum))
  
  mm10_anti_sum <- mm10_anti_tib %>%
    dplyr::group_by(SRC) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  mm10_anti_sum %>% print(n=base::names(mm10_anti_sum))
  
  #
  # Question:: What are these missing values???
  #
  mm10_anti_tib %>% dplyr::filter(SRC=="mm10-LEGX-S2") %>%
    dplyr::inner_join(mm10_unq_man_tib, 
                      by=c("TOP"="Top_Sequence"), 
                      suffix=c("_PRE", "_MAN"))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Compare:: Pre vs. mm10
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Modifications::
  #  - Remove "univ-chicago"
  #  - Rename mm10_unq_man_tib/-mm10_unq_man_tib SRC
  
  fin_cgn_top_tib <- pre_cgn_top_tib %>% 
    dplyr::filter(SRC!="univ-chicago") %>%
    dplyr::mutate(
      SRC=dplyr::case_when(
        SRC=="mm10-LEGX-S2" & CGN %in% mm10_unq_man_tib$CGN ~ "mm10-LEGX-C24A",
        TRUE ~ SRC
      )
    ) %>% 
    dplyr::filter(!stringr::str_detect(TOP,"N")) %>%
    dplyr::arrange(CGN)
  
  #
  # Investigate TOP Seqs with N's
  #
  # CGN TOP                                                                                                                        SRC                               
  # <int> <chr>                                                                                                                      <chr>                             
  # 1   718858 NATTTTAAAATACCCAGCTCCACCCCTTCCTGTTAGGCTTTCGCGTGTCGCAGCTGTGCACGCTGATTGGTCCTCTGCTGGCCAATCACCACTGCACTTCATGACGGCTGTAGTTTTCAAAA HumanMethylation450_15017482_v.1.2
  # 2 10890917 ACAAAAAATATGCTACCAGGGAATTTTTTGTTTTCTAACTAAAGATTCAATCTGGCCAAGCGCAGTGGCCCACACCTGTAATCCCAGCATTTTGGGAGGTGGAGGTGGGAGGATCNNNNNNN HumanMethylation450_15017482_v.1.2
  # 3 12857372 TAAACCGAGAGCTTCTCCAGGTCCCGAAGTTAAGAGTAAATCTGTAGACATTTATCTTATCGCTGCGAAGGGAAACACATCCCAGATGAATTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNN MethylationEPIC_v-1-0_B2          
  # 4 13279585 TTATTTATTAGCCACAAGGGAACTCTTTTTTTTCAAGTTCTTAATCAGAGCACTGGTCATCGTTCCCTGGAGGTGAATCCTGATTATTCATAAGACAAACCTGAATTCNNNNNNNNNNNNNN HumanMethylation27_270596_v.1.2   
  # 5 15312899 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGAATTCATCTGGGATGTGTTTCCCTTCGCAGCGATAAGATAAATGTCTACAGATTTACTCTTAACTTCGGGACCTGGAGAAGCTCTCG MethylationEPIC_v-1-0_B2          
  #
  # Investigate BAD CGN's::
  bad_cgn_vec <- c(3458191, 11037148, 11718090, 15092802, 15408454, 16390856, 23509027)
  fin_cgn_top_tib %>% dplyr::filter(CGN %in% bad_cgn_vec)
  
  fin_cgn_top_sum <- fin_cgn_top_tib %>%
    dplyr::group_by(SRC) %>%
    dplyr::summarise(Count=n(), .groups="drop")
  fin_cgn_top_sum %>% 
    dplyr::arrange(-Count) %>%
    print(n=base::nrow(fin_cgn_top_sum))
  
  fin_cgn_top_csv <- file.path("/Users/bretbarnes/Documents/data/improbe/dbCGN/pre-assigned/canonical-assignment.cgn-top-grp.csv.gz")
  readr::write_csv(fin_cgn_top_tib, fin_cgn_top_csv)
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
