
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages(require("optparse",quietly=TRUE)))

suppressWarnings(suppressPackageStartupMessages(require("plyr")) )
suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )
# suppressWarnings(suppressPackageStartupMessages(require("grid")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

# Manifest RDS Required Packages
# suppressWarnings(suppressPackageStartupMessages(require("GenomicRanges")) )
# suppressWarnings(suppressPackageStartupMessages(require("GenomeInfoDb")) )

# Fasta file reading Packages::
suppressWarnings(suppressPackageStartupMessages(require("Biostrings")) )

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
# Program Parameters::
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'analysis'
par$prgmTag <- 'sample_sheet_scratch'

# Illumina based directories::
par$macDir  <- '/Users/bbarnes/Documents/Projects/methylation/tools'
par$lixDir  <- '/illumina/scratch/darkmatter/Projects/COVIC'

par$improbe_exe <- '/illumina/scratch/darkmatter/bin/improbe'
par$tan_file <- '/illumina/scratch/darkmatter/dat/Tango_A_or_B_11mer_s1.dat'
par$mer_file <- '/illumina/scratch/darkmatter/dat/human-36.1-methyl-and-unmethyl-13mer-s3-for-infinium-methylation.dat'
par$bsp_exe  <- '/illumina/scratch/methylation/software/bsmap-2.90/bsmap'
par$bow_exe  <- '/illumina/thirdparty/bowtie2/bowtie2-2.2.2/bowtie2'

# Directory Parameters::
opt$outDir    <- NULL

# Run Parameters::
opt$runName <- NULL

# Parallel/Cluster Options::
opt$execute  <- TRUE
opt$single   <- FALSE
opt$parallel <- FALSE
opt$cluster  <- FALSE

# Make clean output
opt$clean <- FALSE

# Executable Paramaters::
opt$Rscript <- NULL

# Verbosity Options::
opt$verbose <- 3

cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Parse Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
args.dat <- commandArgs(trailingOnly = FALSE)
if (args.dat[1]=='RStudio') {
  
  if (dir.exists(par$macDir)) par$topDir <- '/Users/bbarnes/Documents/Projects/methylation/scratch'
  if (dir.exists(par$lixDir)) par$topDir <- '/illumina/scratch/darkmatter/data/scratch'
  if (!dir.exists(par$topDir)) dir.create(par$topDir, recursive=TRUE)
  
  # Default Parameters for local Mac::
  par$runMode    <- args.dat[1]
  par$srcDir     <- file.path(par$macDir, par$codeDir)
  par$scrDir     <- file.path(par$srcDir, 'scripts')
  par$exePath    <- file.path(par$scrDir, 'R', par$prgmDir, paste0(par$prgmTag,'.R'))
  
  par$prgmTag <- base::sub('\\.R$', '', base::basename(par$exePath))
  par$locPath <- base::dirname(par$exePath)
  par$scrDir  <- base::dirname(base::normalizePath(par$locPath) )
  par$srcDir  <- base::dirname(base::normalizePath(par$scrDir) )
  par$datDir  <- file.path(base::dirname(base::normalizePath(par$srcDir)), 'dat')
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  
  opt$platform  <- 'EPIC'
  opt$version   <- 'SARS-CoV-2'
  opt$build     <- 'MN908947'
  
  opt$fasta   <- '/Users/bbarnes/Documents/Projects/iGenomes/COVID-19/nCoV_Wuhan_Sequence_MN908947.3.fa.gz'
  opt$runName <- base::basename(opt$fasta) %>% stringr::str_remove('\\.gz$') %>% stringr::str_remove('\\.fa')
  opt$runName <- 'COVIC'
  
  opt$aln_dir <- '/Users/bbarnes/Documents/Projects/methylation/scratch/small.index/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/align'
  opt$aln_dir <- '/Users/bbarnes/Documents/Projects/methylation/scratch/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/n5/aln/align'
  opt$aln_dir <- '/Users/bbarnes/Documents/Projects/methylation/scratch/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/n20000/aln/align/bowtie'
  
  opt$max <- 500
  opt$parallel <- TRUE
  
  des_dir  <- '/Users/bbarnes/Documents/Projects/COVID-19_HLA/data/directDetection/snps'
  opt$des1_csv <- paste(file.path(des_dir, '370992_SARS-CoV-2_probes_F2BT_BEST.score.csv.gz'),
                        file.path(des_dir, '371240_SARS-CoV-2_probes_F2BT_OTHER.score.csv.gz'),
                        sep=',')
  opt$des2_csv <- paste(file.path(des_dir,'370986_SARS-CoV-2_probes_BEST.score.csv.gz'),
                        file.path(des_dir,'371241_SARS-CoV-2_probes_OTHER.score.csv.gz'),
                        sep=',')
  
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
    
    # Directory Parameters::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),

    make_option(c("--max"), type="integer", default=opt$max, 
                help="Max files to process [default= %default]", metavar="integer"),
    
    # Chip Platform and Version Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform name (HM50, EPIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest version (B2, B4, C0) [default= %default]", metavar="character"),
    make_option(c("--build"), type="character", default=opt$build, 
                help="Manifest build (hg19, hg38) [default= %default]", metavar="character"),
    
    # Parallel/Cluster Parameters::
    make_option(c("--execute"), action="store_true", default=opt$execute, 
                help="Boolean variable to shell scripts (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
    # Make clean output
    make_option(c("--clean"), action="store_true", default=opt$clean, 
                help="Boolean variable to clean output directory (mostly testing stuff) [default= %default]", metavar="boolean"),
    
    # Executable Paramaters::
    make_option(c("--Rscript"), type="character", default=opt$Rscript, 
                help="Rscript path [default= %default]", metavar="character"),
    
    # Verbosity Options::
    make_option(c("-v", "--verbose"), type="integer", default=opt$verbose, 
                help="0-5 (5 is very verbosity) [default= %default]", metavar="integer")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Validate Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (is.null(par$runMode) || is.null(par$prgmDir) || is.null(par$prgmTag) || 
    is.null(par$scrDir) || is.null(par$datDir)) {
  
  par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
  par_tib %>% base::print(n=base::nrow(par_tib) )
  
  if (is.null(par$runMode)) cat(glue::glue("[Usage]: runMode is NULL!!!{RET}"))
  if (is.null(par$prgmDir)) cat(glue::glue("[Usage]: prgmDir is NULL!!!{RET}"))
  if (is.null(par$prgmTag)) cat(glue::glue("[Usage]: prgmTag is NULL!!!{RET}"))
  if (is.null(par$scrDir))  cat(glue::glue("[Usage]: scrDir is NULL!!!{RET}"))
  if (is.null(par$datDir))  cat(glue::glue("[Usage]: darDir is NULL!!!{RET}"))
  base::stop("Null Parameters!\n\n")
}

if (is.null(opt$outDir) || 
    is.null(opt$runName) || is.null(opt$aln_dir) ||
    is.null(opt$platform) || is.null(opt$version) || is.null(opt$build) ||
    
    is.null(opt$minPrbScore) || is.null(opt$minCpgRank) || is.null(opt$minScrRank) ||
    is.null(opt$strandCO) || is.null(opt$pickBest) ||
    
    is.null(opt$execute) || is.null(opt$single) || is.null(opt$parallel) || is.null(opt$cluster) ||
    
    is.null(opt$clean) || is.null(opt$Rscript) || is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )
  
  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  if (is.null(opt$outDir))    cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$buildDir))  cat(glue::glue("[Usage]: buildDirs is NULL!!!{RET}"))
  if (is.null(opt$runName))   cat(glue::glue("[Usage]: runName is NULL!!!{RET}"))
  # if (is.null(opt$fasta))     cat(glue::glue("[Usage]: fasta is NULL!!!{RET}"))
  if (is.null(opt$aln_dir))   cat(glue::glue("[Usage]: aln_dir is NULL!!!{RET}"))
  
  if (is.null(opt$des1_csv)) cat(glue::glue("[Usage]: des1_csv is NULL (not required)!!!{RET}"))
  if (is.null(opt$des2_csv)) cat(glue::glue("[Usage]: des2_csv is NULL (not required)!!!{RET}"))
  
  if (is.null(opt$platform)) cat(glue::glue("[Usage]: platform is NULL!!!{RET}"))
  if (is.null(opt$version))  cat(glue::glue("[Usage]: version is NULL!!!{RET}"))
  if (is.null(opt$build))    cat(glue::glue("[Usage]: build is NULL!!!{RET}"))
  
  if (is.null(opt$minPrbScore))  cat(glue::glue("[Usage]: minPrbScore is NULL!!!{RET}"))
  if (is.null(opt$minCpgRank))   cat(glue::glue("[Usage]: minCpgRank is NULL!!!{RET}"))
  if (is.null(opt$minScrRank))   cat(glue::glue("[Usage]: minScrRank is NULL!!!{RET}"))
  if (is.null(opt$strandCO))     cat(glue::glue("[Usage]: strandCO is NULL!!!{RET}"))
  if (is.null(opt$pickBest))     cat(glue::glue("[Usage]: pickBest is NULL!!!{RET}"))
  
  if (is.null(opt$execute))  cat(glue::glue("[Usage]: execute is NULL!!!{RET}"))
  if (is.null(opt$single))   cat(glue::glue("[Usage]: single is NULL!!!{RET}"))
  if (is.null(opt$parallel)) cat(glue::glue("[Usage]: parallel is NULL!!!{RET}"))
  if (is.null(opt$cluster))  cat(glue::glue("[Usage]: cluster is NULL!!!{RET}"))
  
  if (is.null(opt$clean))   cat(glue::glue("[Usage]: clean is NULL!!!{RET}"))
  if (is.null(opt$Rscript)) cat(glue::glue("[Usage]: Rscript is NULL!!!{RET}"))
  if (is.null(opt$verbose)) cat(glue::glue("[Usage]: verbosity is NULL!!!{RET}"))
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

par$prgm_src_dir <- file.path(par$scrDir,par$prgmDir, 'functions')
if (!dir.exists(par$prgm_src_dir)) stop(glue::glue("[{par$prgmTag}]: Program Source={par$prgm_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$prgm_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Program Source={par$prgm_src_dir}!{RET}{RET}") )

par$man_src_dir <- file.path(par$scrDir, 'manifests/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$man_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$man_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$man_src_dir}!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# pTracker <- timeTracker$new(verbose=opt$verbose)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Defined Outputs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# opt$outDir <- file.path(opt$outDir, par$prgmTag, opt$platform, opt$version, opt$build, opt$runName)
# if (!is.null(opt$max)) opt$outDir <- file.path(opt$outDir, paste0('n',opt$max) )
# if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
# cat(glue::glue("[{par$prgmTag}]: Built; OutDir={opt$outDir}!{RET}") )
# 
# if (opt$clean) list.files(opt$outDir, full.names=TRUE) %>% unlink()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Load Manifests::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sam_dir <- '/Users/bbarnes/Documents/Projects/methylation/sampleSheets'

# old_ss_csv <- '/Users/bbarnes/Documents/CustomerFacing/sampleSheets/annotation/Human-Classification_COVID_Count-656_AnnotatedOnlySampleSheet.csv'
old_ss_csv <- file.path(sam_dir, 'annotation/Human-Classification_COVID_Count-656_AnnotatedOnlySampleSheet.csv')
new_ss_csv <- '/Users/bbarnes/Documents/Projects/COVID-19_HLA/data/phenotype_data/Finalpd_basic.clean.csv'

old_ss_tib <- suppressMessages(suppressWarnings( readr::read_csv(old_ss_csv) ))
new_ss_tib <- suppressMessages(suppressWarnings( readr::read_csv(new_ss_csv) ))

add_ss_tib <- new_ss_tib %>% dplyr::mutate(
  Sample_Class=dplyr::case_when(
    covid_status=='Positive' ~ 'pSARSCov2', covid_status=='Negative' ~ 'nSARSCov2', TRUE ~ NA_character_),
  COVID_Status=dplyr::case_when(
    covid_status=='Positive' ~ '+', covid_status=='Negative' ~ '-', TRUE ~ NA_character_),
  Tissue_Source='Whole_Blood',
  Source_Sample_ID=as.character(sample_id)
  ) %>% 
  dplyr::rename(Sample_Plate=sample_plate) %>% 
  dplyr::select(Sentrix_Name,Sample_Class,COVID_Status,Source_Sample_ID,Tissue_Source,Sample_Plate)

mis_cnt <- dplyr::inner_join(add_ss_tib,old_ss_tib, by="Sentrix_Name") %>% dplyr::filter(Sample_Class.x != Sample_Class.y) %>% base::nrow()

out_ss_csv <- file.path(sam_dir, 'annotation/Human-Classification_COVID_Count-921_AnnotatedOnlySampleSheet.csv.gz')
out_ss_tib <- dplyr::bind_rows(
  dplyr::anti_join(add_ss_tib,old_ss_tib, by="Sentrix_Name"),
  old_ss_tib) %>% dplyr::arrange(Sentrix_Name)

readr::write_csv(out_ss_tib, out_ss_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sysTime <- Sys.time()
cat(glue::glue("{RET}[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
