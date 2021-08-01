
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Global Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM  <- ","
TAB  <- "\t"
RET  <- "\n"
RET2 <- "\n\n"
BNG  <- "|"
BRK  <- paste0("# ",
               paste(rep("-----",6),collapse=" "),"|",
               paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Define Default Params and Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt <- NULL  # List of user options
par <- NULL  # List of static program parameters (NOT accessible to user)
run <- NULL  # List of non-static run-time program parameters

par$date    <- Sys.Date() %>% as.character()
par$runMode <- ''
par$maxTest <- NULL

# Default local Mac/sd-isilon directories for ease of use::
par$macDir1 <- '/Users/bbarnes/Documents/Projects/methylation'
par$macDir2 <- '/Users/bretbarnes/Documents'
par$lixDir1 <- '/illumina/scratch/darkmatter'
par$lixDir  <- '/illumina/scratch/darkmatter'

# Program Name Parameters::
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'workhorse'
par$prgmTag <- 'workhorse_main_dev_latest'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Run Time Version Options:: 
#                       Platform, Genome Build, etc
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$run_name     <- NULL
opt$platform     <- NULL
opt$version      <- NULL
opt$genome_build <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Run Time User Input Directories:: 
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$out_dir <- NULL

opt$ord_dir <- NULL
opt$mat_dir <- NULL
opt$aqp_dir <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Run Time User Input Files:: 
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$ord_csv <- NULL
opt$mat_tsv <- NULL
opt$aqp_tsv <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Run Time User Input Executable(s):: 
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$Rscript   <- NULL
opt$bsmap_opt <- NULL
opt$bsmap_exe <- NULL
opt$align_chroms <- FALSE

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Pre-defined Static Data Directories:: 
#            improbe, Annotation, Genomic, Manifest, Validation Idats
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$imp_dir  <- NULL
opt$ann_dir  <- NULL
opt$gen_dir  <- NULL
opt$man_dir  <- NULL
opt$idat_dir <- NULL

opt$cgn_seq_dir <- NULL
opt$cgn_bed_dir <- NULL

opt$canonical_cgn_dir <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Pre-defined Static External File Options:: 
#                   Manifest, Controls, Design Coordinates
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$sesame_manfiest_dat   <- NULL
opt$sesame_manifest_csv   <- NULL
opt$genome_manifest_csv   <- NULL

opt$sesame_controls_csv   <- NULL
opt$genome_controls_csv   <- NULL
opt$noob_controls_csv     <- NULL

opt$source_coordinate_csv <- NULL
opt$canonical_cgn_csv     <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Run Time File Options:: Time Stamps
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
opt$time_org_txt <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Run Time Mode Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$single    <- FALSE
opt$parallel  <- FALSE
opt$cluster   <- FALSE

opt$trackTime <- NULL
opt$fresh     <- FALSE
opt$reload    <- FALSE

opt$verbose   <- 3

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Order/Match/AQP/PQC Expected Columns::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# NEW SOLUTION FOR FILE HEADERS::
#  - Store all file headers in a single rData (rds) file. They're 
#    minature, but required for files without headers and they 
#    massively speed up reading of very large data files since
#    variable types can be expected!
#
# TBD:: This is copied in the code, so it should be removed from here
#   or moved to a config file... Maybe a single r-data structure (RDS)
#   for manifest parameter defaults is the way to go...
#

# par_cols <- list()
# par_cols$ord <- 
#   cols(
#     Assay_Design_Id        = col_character(),
#     AlleleA_Probe_Id       = col_character(),
#     AlleleA_Probe_Sequence = col_character(),
#     AlleleB_Probe_Id       = col_character(),
#     AlleleB_Probe_Sequence = col_character(),
#     Normalization_Bin      = col_character()
#   )
# 
# par_cols$mat <- 
#   cols(
#     Plate    = col_character(),
#     Row      = col_character(),
#     Col      = col_integer(),
#     Address  = col_integer(),
#     Mod5     = col_character(),
#     Sequence = col_character(),
#     Mod3     = col_character(),
#     Comments = col_character()
#   )
# 
# par_cols$ma2 <- 
#   cols(
#     address_names = col_integer(),
#     probe_id      = col_character(),
#     sequence      = col_character(),
#     type_b        = col_character(),
#     address_name  = col_integer(),
#     bo_seq        = col_character()
#   )
# 
# par_cols$aqp <- 
#   cols(
#     Address           = col_integer(),
#     Decode_Status     = col_integer(),
#     Decode_Error_Code = col_integer(),
#     Decode_Score      = col_integer(),
#     Func_Status       = col_integer(),
#     Func_Error_Code   = col_integer(),
#     QC_Action         = col_integer()
#   )
# 
# par_cols$pqc <- 
#   cols(
#     Address      = col_integer(),
#     Status       = col_integer(),
#     Eval_Code    = col_integer(),
#     Average_Rep  = col_integer(),
#     Expected_Rep = col_integer()
#   )
# 
# par$ord_col <- par_cols$ord$cols %>% names()
# 
# par$mat_col <- par_cols$mat$cols %>% names()
# par$ma2_col <- par_cols$ma2$cols %>% names()
# 
# par$aqp_col <- par_cols$aqp$cols %>% names()
# par$pqc_col <- par_cols$pqc$cols %>% names()

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
  
  opt$out_dir  <- file.path(par$topDir, 'scratch')
  opt$imp_dir  <- file.path(par$topDir, 'data/improbe')
  opt$ann_dir  <- file.path(par$topDir, 'data/annotation')
  opt$man_dir  <- file.path(par$topDir, 'data/manifests')
  opt$gen_dir  <- file.path(par$topDir, 'data/iGenomes/Homo_sapiens/NCBI')
  opt$idat_dir <- file.path(par$topDir, 'data/idats')
  
  # opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
  # opt$bsmap_opt <- "\"-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R\""
  
  opt$bsmap_opt <- "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
  opt$bsmap_exe <- "/Users/bretbarnes/Documents/tools/programs/BSMAPz/bsmapz"
  opt$cgn_seq_dir <- 
    file.path(opt$imp_dir, "scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49-split")
  opt$cgn_bed_dir <- 
    file.path(opt$imp_dir, "scratch/cgnDB/dbSNP_Core4/design-input/min")
  
  opt$canonical_cgn_dir <- file.path(par$datDir, "manifest/cgnDB")
  opt$canonical_cgn_csv <- "canonical.cgn-top-grp.csv.gz"
  
  # opt$genome_controls_csv <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
  # opt$sesame_manifest_csv <- file.path(par$topDir,"data/manifests/methylation/Sesame/EPIC-B4-BP4.manifest.sesame-base.controls-only.csv.gz")
  
  # Pre-defined local options runTypes::
  #
  par$local_runType <- NULL
  par$local_runType <- 'COVIC'
  par$local_runType <- 'GSA'
  par$local_runType <- 'HM450'
  par$local_runType <- 'TruDx'
  par$local_runType <- 'GRCm10'
  par$local_runType <- 'NZT'
  par$local_runType <- 'EWAS'
  par$local_runType <- 'EPIC_v2'
  par$local_runType <- 'Chicago'
  par$local_runType <- 'McMaster10Kselection'
  
  opt$parallel <- TRUE
  
  opt$verbose <- 10
  opt$verbose <- 3
  opt$verbose <- 5
  
  opt$fresh  <- TRUE
  opt$fresh  <- FALSE
  opt$reload <- TRUE
  
  if (FALSE) {
    
  } else if (par$local_runType=='EPIC_v2') {
    opt$genome_build <- 'GRCh37'
    opt$platform <- 'EPIC'
    opt$version  <- 'v4'
    
    opt$sesame_manifest_dat <- "EPIC.hg19.manifest,HM450.hg19.manifest"
    genome_manifest_dir <- file.path(par$topDir, "data/manifests/methylation/GenomeStudio")
    opt$genome_manifest_csv <- paste(
      file.path(genome_manifest_dir, "MethylationEPIC_v-1-0_B4-Beadpool_ID.csv.gz"),
      file.path(genome_manifest_dir, "HumanMethylation450_15017482_v.1.2.csv.gz"),
      sep = ","
    )
    
  } else if (par$local_runType=='McMaster10Kselection') {
    opt$genome_build <- 'GRCh37'
    opt$platform <- 'MCM'
    opt$version  <- 'v3'
    
    opt$ord_dir <- file.path(par$topDir, "data/CustomContent/McMaster/McMaster10Kselection/AQP.v2")
    opt$mat_dir <- file.path(par$topDir, "data/CustomContent/McMaster/McMaster10Kselection/AQP.v2")
    opt$aqp_dir <- file.path(par$topDir, "data/CustomContent/McMaster/McMaster10Kselection/AQP.v2")
    
    opt$ord_csv <- paste('McMaster_CpG_DesignFile_v4.csv.gz',sep=',')
    opt$mat_tsv <- paste('20532820_probes1.match.gz',
                         '20532820_probes2.match.gz', sep=',')
    opt$aqp_tsv <- paste('20051339_A_ProductQC.txt.gz', sep=',')
    
    if (FALSE) {
      opt$aqp_tsv <- paste('BS0033057-AQP1.txt.gz',
                           'BS0033090-AQP2.txt.gz',
                           #'20051339_A_ProductQC.txt.gz',
                           sep=',')
    }
    
    opt$noob <- paste(
      file.path(par$topDir, "data/CustomContent/transfer/updated_manifest.csv.gz"),
      sep=',')
    
  } else if (par$local_runType=='Chicago') {
    opt$genome_build <- 'GRCh38'
    opt$genome_build <- 'GRCh37'
    opt$platform <- 'EPIC'
    opt$version  <- 'B3'
    
    opt$idat_dir <- file.path(opt$idat_dir, "idats_Chicago-Ober-Custom")
    par$ord_dir  <- file.path(par$topDir, 'data/CustomContent/UnivChicago/latest')
    par$mat_dir  <- file.path(par$topDir, 'data/CustomContent/UnivChicago/latest')
    par$aqp_dir  <- file.path(par$topDir, 'data/CustomContent/UnivChicago/latest')
    
    opt$ord_csv <- paste('UofChicago-A_A_Array-CpG-order-FINAL.csv', sep=',')
    opt$mat_tsv <- paste('20504790_probes.match.tsv', sep=',')
    opt$aqp_tsv <- paste('329922X374054_A_ProductQC.txt', sep=',')
    
    # Use this later in the process for picking coordinates::
    par$ord_pos_csv <- file.path(par$topDir, "data/CustomContent/UnivChicago/improbe_input/CpGs_UnivChicago_alldesigns_55860sites.cgn-pos-srd-prbs.tsv.gz")
    
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$run_name <- paste(par$local_runType,opt$platform,opt$version,opt$genome_build, sep='-')
  
} else {
  par$runMode    <- 'CommandLine'
  par$exePath <- base::substring(args.dat[grep("--file=", args.dat)], 8)
  
  par$prgmTag <- base::sub('\\.R$', '', base::basename(par$exePath))
  par$locPath <- base::dirname(par$exePath)
  par$scrDir  <- base::dirname(base::normalizePath(par$locPath) )
  par$srcDir  <- base::dirname(base::normalizePath(par$scrDir) )
  par$datDir  <- file.path(base::dirname(base::normalizePath(par$srcDir)), 'dat')
  
  args.dat <- base::commandArgs(trailingOnly = TRUE)
  option_list = list(
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Run Time Version Options:: 
    #                       Platform, Genome Build, etc
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Run Parameters::
    optparse::make_option(
      c("--run_name"), type="character", default=opt$run_name, 
      help=paste0("Run Name [default= %default]"), 
      metavar="character"),
    
    # Platform/Method Options::
    optparse::make_option(
      c("--platform"), type="character", default=opt$platform, 
      help=paste0("Platform (e.g. HM450, EPIC, LEGX, NZT, ",
                  "COVIC) [default= %default]"), 
      metavar="character"),
    optparse::make_option(
      c("--version"), type="character", default=opt$version, 
      help=paste0("Manifest Version (e.g. B0,B1,B2,B3,B4,C0) ",
                  "[default= %default]"), 
      metavar="character"),
    optparse::make_option(
      c("--genome_build"), type="character", default=opt$genome_build, 
      help=paste0("Genome Build (e.g. GRch36, GRCh37, GRCh38, GRCm38) ",
                  "[default= %default]"), 
      metavar="character"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Run Time User Input Directories:: 
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    optparse::make_option(
      c("--out_dir"), type="character", default=opt$out_dir, 
      help=paste0("Output directory [default= %default]"), 
      metavar="character"),
    
    # Manufacturing Files:: Required
    optparse::make_option(
      c("--ord_dir"), type="character", default=opt$ord_dir, 
      help=paste0("Order directories. Either a single directory or one-to-one ",
                  "pairing with order files (comma seperated list) ",
                  "[default= %default]"), 
      metavar="character"),
    optparse::make_option(
      c("--mat_dir"), type="character", default=opt$mat_dir, 
      help=paste0("Biziprobe Match directories. Either a single directory or ",
                  "one-to-one pairing with match files (comma seperated list) ",
                  "[default= %default]"), 
      metavar="character"),
    optparse::make_option(
      c("--aqp_dir"), type="character", default=opt$aqp_dir, 
      help=paste0("AQP/PQC directories. Either a single directory or ",
                  "one-to-one pairing with AQP/PQC files ",
                  "(comma seperated list) ",
                  "[default= %default]"), 
      metavar="character"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Run Time User Input Files:: 
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Manufacturing Files:: Required
    optparse::make_option(
      c("--ord_csv"), type="character", default=opt$ord_csv,
      help=paste0("Order file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),
    optparse::make_option(
      c("--mat_tsv"), type="character", default=opt$mat_tsv, 
      help=paste0("Biziprobe Match file(s) (comma seperated list) ",
                  "[default= %default]"), 
      metavar="character"),
    optparse::make_option(
      c("--aqp_tsv"), type="character", default=opt$aqp_tsv, 
      help=paste0("AQP/PQC file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Run Time User Input Executable(s):: 
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    optparse::make_option(
      c("--Rscript"), type="character", default=opt$Rscript, 
      help=paste0("Rscript path [default= %default]"), 
      metavar="character"),
    
    optparse::make_option(
      c("--bsmap_opt"), type="character", default=opt$bsmap_opt, 
      help=paste0("BSMAP Options [default= %default]"),
      metavar="character"),
    optparse::make_option(
      c("--bsmap_exe"), type="character", default=opt$bsmap_exe, 
      help=paste0("BSMAP Executable path [default= %default]"),
      metavar="character"),
    optparse::make_option(
      c("--align_chroms"), action="store_true", default=opt$align_chroms, 
      help=paste0("Boolean flag to align against individual chromosomes. ",
                  "This provides more alignments than the best hit. ",
                  "[default= %default]"),
      metavar="boolean"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Pre-defined Static Data Directories:: 
    #            improbe, Annotation, Genomic, Manifest, Validation Idats
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    optparse::make_option(
      c("--imp_dir"), type="character", default=opt$imp_dir, 
      help=paste0("improbe data directory [default= %default]"),
      metavar="character"),
    optparse::make_option(
      c("--ann_dir"), type="character", default=opt$ann_dir, 
      help=paste0("Annotation data directory [default= %default]"),
      metavar="character"),
    optparse::make_option(
      c("--gen_dir"), type="character", default=opt$gen_dir, 
      help=paste0("Genomic data directory [default= %default]"),
      metavar="character"),
    optparse::make_option(
      c("--man_dir"), type="character", default=opt$man_dir, 
      help=paste0("Pre-built Manifest data directory [default= %default]"),
      metavar="character"),
    
    # Validation existing idats directory to confirm Addresses against::
    optparse::make_option(
      c("--idat_dir"), type="character", default=opt$idat_dir, 
      help=paste0("Validation existing idats directory ",
                  "to confirm Addresses against. ",
                  "CSV file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),
    
    # Pre-defined directory with files containing coordinates to cg numbers::
    optparse::make_option(
      c("--cgn_seq_dir"), type="character", default=opt$cgn_seq_dir, 
      help=paste0("Pre-defined directory with files containing coordinates ",
                  "to cg numbers. [default= %default]"),
      metavar="character"),
    
    # Pre-defined directory with files containing coordinates to cg numbers::
    optparse::make_option(
      c("--cgn_bed_dir"), type="character", default=opt$cgn_bed_dir, 
      help=paste0("Pre-defined directory with files containing coordinates ",
                  "to cg numbers. [default= %default]"), 
      metavar="character"),

    # Pre-defined directory with file containing canonical cg number assignments::
    optparse::make_option(
      c("--canonical_cgn_dir"), type="character", default=opt$canonical_cgn_dir, 
      help=paste0("Pre-defined canonical cg-numbers file used for cg number ",
                  "resolution assignment. Directory path not file name(s)! ",
                  "CSV file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Pre-defined Static External File Options:: 
    #                   Manifest, Controls, Design Coordinates
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Pre-defined manifest(s) to be re-built and/or added to new manifest
    #  from Sesame Repo::
    optparse::make_option(
      c("--sesame_manfiest_dat"), type="character", 
      default=opt$sesame_manfiest_dat,
      help=paste0("Sesame Manifest(s) to be re-built and/or added to ",
                  "new manifest from Sesame Repo. ",
                  "Example = 'HM450.hg19.manifest,EPIC.hg19.manifest",
                  "CSV file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),
    
    # Pre-defined manifest(s) to be re-built and/or added to new manifest::
    optparse::make_option(
      c("--sesame_manifest_csv"), type="character", 
      default=opt$sesame_manifest_csv,
      help=paste0("Sesame Manifest(s) to be re-built and/or added ",
                  "to new manifest. Probe Seq required! ",
                  "CSV file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),
    optparse::make_option(
      c("--genome_manifest_csv"), type="character", 
      default=opt$genome_manifest_csv,
      help=paste0("Genome Studio Manifest(s) to be re-built and/or ",
                  "added to new manifest. Probe Seq required! ",
                  "CSV file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),
    
    # Pre-defined manifest control(s) to be added to new manifest::
    optparse::make_option(
      c("--sesame_controls_csv"), type="character", 
      default=opt$sesame_controls_csv, 
      help=paste0("Sesame Pre-defined manifest control(s)  ",
                  "to be added to new manifest. ",
                  "CSV file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),
    optparse::make_option(
      c("--genome_controls_csv"), type="character", 
      default=opt$genome_controls_csv, 
      help=paste0("Genome Studio Pre-defined manifest control(s)  ",
                  "to be added to new manifest. ",
                  "CSV file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),
    
    # Pre-defined noob-masked control(s) to be added to new manifest::
    optparse::make_option(
      c("--noob_controls_csv"), type="character", 
      default=opt$noob_controls_csv, 
      help=paste0("Noob-Masked Pre-defined control(s) ",
                  "to be added to new manifest. ",
                  "CSV file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),
    
    # Original source design file used for canonical position selection::
    optparse::make_option(
      c("--source_coordinate_csv"), type="character", 
      default=opt$source_coordinate_csv, 
      help=paste0("Original source design file used for canonical ",
                  "position selection. ",
                  "CSV file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),
    
    # Pre-defined file names with canonical cg number assignments::
    optparse::make_option(
      c("--canonical_cgn_csv"), type="character", 
      default=opt$canonical_cgn_csv, 
      help=paste0("Pre-defined canonical cg-numbers file used for cg number ",
                  "resolution assignment. File(s) name, not path(s)! ",
                  "CSV file(s) (comma seperated list) [default= %default]"), 
      metavar="character"),
    
    # optparse::make_option(c("--cpg_pos_tsv"), type="character", default=opt$cpg_pos_tsv, 
    #             help="Null value for passing arguments [default= %default]", metavar="character"),
    # optparse::make_option(c("--cpg_top_tsv"), type="character", default=opt$cpg_top_tsv, 
    #             help="Null value for passing arguments [default= %default]", metavar="character"),
    # optparse::make_option(c("--cph_pos_tsv"), type="character", default=opt$cph_pos_tsv, 
    #             help="Null value for passing arguments [default= %default]", metavar="character"),
    # optparse::make_option(c("--snp_pos_tsv"), type="character", default=opt$snp_pos_tsv, 
    #             help="Null value for passing arguments [default= %default]", metavar="character"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Run Time File Options:: Time Stamps
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    optparse::make_option(
      c("--time_org_txt"), type="character", default=opt$time_org_txt, 
      help=paste0("Unused variable time_org_txt [default= %default]"),
      metavar="character"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Run Time Mode Options::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Process Parallel/Cluster Parameters::
    optparse::make_option(
      c("--single"), action="store_true", default=opt$single, 
      help=paste0("Boolean variable to run a single sample on a single-core ",
                  "[default= %default]"),
      metavar="boolean"),
    optparse::make_option(
      c("--parallel"), action="store_true", default=opt$parallel, 
      help=paste0("Boolean variable to run parallel on multi-core ",
                  "[default= %default]"), 
      metavar="boolean"),
    optparse::make_option(
      c("--cluster"), action="store_true", default=opt$cluster,
      help=paste0("Boolean variable to run jobs on cluster by chip ",
                  "[default= %default]"),
      metavar="boolean"),
    
    # Run=time Options::
    optparse::make_option(
      c("--trackTime"), action="store_true", default=opt$trackTime,
      help=paste0("Boolean variable tack run times [default= %default]"),
      metavar="boolean"),
    optparse::make_option(
      c("--fresh"), action="store_true", default=opt$fresh, 
      help=paste0("Boolean variable to run a fresh build [default= %default]"),
      metavar="boolean"),
    optparse::make_option(
      c("--reload"), action="store_true", default=opt$reload, 
      help=paste0("Boolean variable reload intermediate files (for testing). ",
                  "[default= %default]"),
      metavar="boolean"),
    
    # Verbosity level::
    optparse::make_option(
      c("-v", "--verbose"), type="integer", default=opt$verbose, 
      help=paste0("Verbosity level: 0-5 (5 is very verbose) ",
                  "[default= %default]"), 
      metavar="integer")
  )
  opt_parser = optparse::OptionParser(option_list=option_list)
  opt = optparse::parse_args(opt_parser)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Program Initialization::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_reqs <- c('runMode','prgmTag','scrDir','datDir','exePath')
opt_reqs <- c('out_dir','imp_dir',
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

pTracker <- timeTracker$new()
pTracker$addFile(opt$time_org_txt)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Pre-processing:: 
#                      Pre-defined & Intermediate Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Checking pre-defined files...{RET}"))

# For error ledger for accumulating failed probes and why::
error_ledger <- NULL

run <- NULL
run <- get_run_defaults(fresh = opt$fresh,
                        genome_build = opt$genome_build,
                        cgn_seq_dir  = opt$cgn_seq_dir,
                        cgn_bed_dir  = opt$cgn_bed_dir,
                        canonical_cgn_dir = opt$canonical_cgn_dir,
                        canonical_cgn_csv = opt$canonical_cgn_csv,
                        verbose = opt$verbose )


if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Checking pre-defined files.{RET2}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        0.0 Validate AQP Inputs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

valid_files <- valid_aqp_inputs(ord_dir = opt$ord_dir, ord_csv = opt$ord_csv,
                                mat_dir = opt$mat_dir, mat_tsv = opt$mat_tsv,
                                aqp_dir = opt$aqp_dir, aqp_tsv = opt$aqp_tsv,
                                verbose = opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#          0.1 Load any pre-defined Noob-Masked Manifest to be added::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# noob_ctl_tib <- NULL
# if (!is.null(opt$noob)) {
#   noob_ctl_tib <- noob_mask(noob_csv = opt$noob, 
#                             ctl_csv = opt$sesame_manifest_csv, 
#                             verbose=opt$verbose, tt=pTracker)
# }

#
# TBD:: User should be able to rebuild an existing or old manifest,
#  or add manifests together...
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      0.2 Load any other pre-defined data
#                        dbCGN, improbe, imGenomes data::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imGenome_tib <- load_imGenomes_table(dir = opt$gen_dir,
                                     genome_build = opt$genome_build, 
                                     ret_list = FALSE, 
                                     load_chroms = opt$align_chroms,
                                     verbose = opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   1.0 AQP Address Manifest Workflow: 
#                           Order/Match/AQP/PQC
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ord_tib <- NULL
ord_tib <- 
  aqp_mapping_workflow(ord = valid_files$ord_str,
                       mat = valid_files$mat_str,
                       aqp = valid_files$aqp_str,
                       
                       prb_key = run$prb_key,
                       add_key = run$add_key,
                       des_key = run$des_key,
                       din_key = run$din_key,
                       ids_key = run$ids_key,
                       
                       out_dir = opt$out_dir,
                       out_col = run$out_col,
                       run_tag = opt$run_name,
                       re_load = run$re_load,
                       pre_tag = pTracker$file_vec,
                       
                       verbose=opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    2.0 Align All Probe Sequence:: BSMAP
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Plot normalized intensity x individual chr hits x binding energy
#

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
                                
                                out_dir = opt$out_dir,
                                out_col = run$out_col,
                                run_tag = opt$run_name,
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
                       
                       seq_dir   = opt$cgn_seq_dir,
                       pattern_u = run$seq_pattern_U, 
                       pattern_m = run$seq_pattern_M,
                       
                       prb_key = run$prb_key,
                       add_key = run$add_key,
                       des_key = run$des_key,
                       din_key = run$din_key,
                       ids_key = run$ids_key,
                       
                       prefix = opt$run_name,
                       suffix = run$seq_suffix, 
                       
                       idxA = run$seq_idxA,
                       idxB = run$seq_idxB,
                       
                       reload   = opt$reload,
                       parallel = opt$parallel,
                       
                       out_dir = opt$out_dir,
                       out_col = run$out_col,
                       run_tag = opt$run_name,
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
                       
                       can_csv = run$canonical_cgn_csv,
                       
                       join    = run$cgn_join,
                       merge   = run$cgn_merge,
                       retData = FALSE,
                       
                       out_dir = opt$out_dir,
                       out_col = run$out_col,
                       unq_col = run$unq_col,
                       run_tag = opt$run_name,
                       re_load = run$re_load,
                       pre_tag = pTracker$file_vec,
                       
                       verbose=opt$verbose, tt=pTracker)

if (FALSE) {
  print(ord_tib)
  print(bsp_tib)
  print(seq_tib)
  print(cgn_tib)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
# Ending the first of this script here and developing the second half in a 
#   completely separate script. This is due to some naming convention 
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
