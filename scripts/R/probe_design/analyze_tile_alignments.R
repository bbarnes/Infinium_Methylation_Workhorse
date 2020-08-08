
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
suppressWarnings(suppressPackageStartupMessages(require("grid")) )

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
par$prgmDir <- 'probe_design'
par$prgmTag <- 'analyze_tile_alignments'

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
opt$fasta   <- NULL
opt$aln_dir <- NULL
opt$genome  <- NULL

opt$des1_csv <- NULL
opt$des2_csv <- NULL

# Probe Filter Parameters::
opt$minPrbScore <- 0.3
opt$minCpgRank  <- "3,2,1,0"
opt$minScrRank  <- "0.6,0.5,0.4,0.3"
opt$strandCO    <- 'C,O'
opt$pickBest    <- FALSE

opt$max <- NULL

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
    make_option(c("--fasta"), type="character", default=opt$fasta, 
                help="Whole Genome Fasta [default= %default]", metavar="character"),
    make_option(c("--genome"), type="character", default=opt$genome, 
                help="Genome Fasta to align probes against [default= %default]", metavar="character"),
    make_option(c("--aln_dir"), type="character", default=opt$aln_dir, 
                help="Alignment Directory [default= %default]", metavar="character"),

    make_option(c("--des1_csv"), type="character", default=opt$des1_csv, 
                help="SNP Designs Ininium I (comma seperate list) [default= %default]", metavar="character"),
    make_option(c("--des2_csv"), type="character", default=opt$des2_csv, 
                help="SNP Designs Ininium II (comma seperate list) [default= %default]", metavar="character"),
    
    make_option(c("--max"), type="integer", default=opt$max, 
                help="Max files to process [default= %default]", metavar="integer"),
    
    # Chip Platform and Version Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform name (HM50, EPIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest version (B2, B4, C0) [default= %default]", metavar="character"),
    make_option(c("--build"), type="character", default=opt$build, 
                help="Manifest build (hg19, hg38) [default= %default]", metavar="character"),
    
    # Probe Filter Parameters::
    make_option(c("--minPrbScore"), type="double", default=opt$minPrbScore,
                help="Minimum Probe Score [default= %default]", metavar="double"),
    make_option(c("--minCpgRank"), type="character", default=opt$minCpgRank, 
                help="Infinium I/II Probe Underlying CpG Count Rank String (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--minScrRank"), type="character", default=opt$minScrRank, 
                help="Infinium I/II Probe Score Rank String (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--strandCO"), type="character", default=opt$strandCO, 
                help="Target Strand Converted/Opposite to design (C,O) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--pickBest"), action="store_true", default=opt$pickBest, 
                help="Boolean variable to only return the best scoring probe for each strand [default= %default]", metavar="boolean"),
    
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

pTracker <- timeTracker$new(verbose=opt$verbose)

des1_csv_vec <- NULL
des2_csv_vec <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Defined Outputs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$outDir <- file.path(opt$outDir, par$prgmTag, opt$platform, opt$version, opt$build, opt$runName)
if (!is.null(opt$max)) opt$outDir <- file.path(opt$outDir, paste0('n',opt$max) )
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
cat(glue::glue("[{par$prgmTag}]: Built; OutDir={opt$outDir}!{RET}") )

if (opt$clean) list.files(opt$outDir, full.names=TRUE) %>% unlink()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Load Probe Designs:: cDNA
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (!is.null(opt$des1_csv) && !is.null(opt$des2_csv) ) {
  
  des1_csv_vec <- opt$des1_csv %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  des2_csv_vec <- opt$des2_csv %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

  snp_des1_tib <- suppressMessages(suppressWarnings( lapply(des1_csv_vec, readr::read_csv, skip=15) )) %>%
    dplyr::bind_rows()
  snp_des2_tib <- suppressMessages(suppressWarnings( lapply(des2_csv_vec, readr::read_csv, skip=15) )) %>%
    dplyr::bind_rows()
  
  # Look at Allele Cases::
  #  Design Studio Always forces AT -> CG SNPs
  # snp_des1_tib %>% dplyr::mutate(AlleleA=stringr::str_sub(AlleleA_Probe_Sequence, stringr::str_length(AlleleA_Probe_Sequence)),
  #                                AlleleB=stringr::str_sub(AlleleB_Probe_Sequence, stringr::str_length(AlleleB_Probe_Sequence)) ) %>% 
  #   dplyr::select(Sequence, AlleleA_Probe_Sequence,AlleleA, AlleleB_Probe_Sequence,AlleleB) %>% dplyr::group_by(AlleleA,AlleleB) %>% dplyr::summarise(Count=n())
  #
  # Missing Alignment:: Always CG/AT...
  #   snp_des1_tib %>% dplyr::filter(AlleleA_Probe_Sequence=='CTCTGCAAAACAGCTGAGGTGATAGAGGTTTGTGGTGGTTGGTAAAGAAT') %>% as.data.frame()
  
  # Need to preprocess and format before joining...
  #  - [Done]: Strands need to match!
  snp_des_tib <- dplyr::inner_join(
    snp_des1_tib %>% dplyr::select(Ilmn_Id,Locus_Name,Sequence,Coordinate,Final_Score,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence) %>%
      tidyr::separate(Ilmn_Id, into=c("IDA","IDB"), sep='-', remove=FALSE) %>%
      tidyr::separate(IDB, into=c("IDN","TB","FR","Assay_Index"), sep='_', remove=TRUE) %>%
      dplyr::rename(Prb_Seq_IA=AlleleA_Probe_Sequence,Prb_Seq_IB=AlleleB_Probe_Sequence) %>%
      dplyr::select(-IDA, -IDN),
    snp_des2_tib %>% dplyr::select(Ilmn_Id,Locus_Name,Sequence,Coordinate,Final_Score,AlleleA_Probe_Sequence) %>%
      tidyr::separate(Ilmn_Id, into=c("IDA","IDB"), sep='-', remove=FALSE) %>%
      tidyr::separate(IDB, into=c("IDN","TB","FR","Assay_Index"), sep='_', remove=TRUE) %>%
      dplyr::rename(Prb_Seq_II=AlleleA_Probe_Sequence) %>%
      dplyr::select(-IDA, -IDN),
    by=c("Sequence","Coordinate","TB","FR"), suffix=c("_I","_II")
  ) %>% dplyr::arrange(Coordinate,TB,FR)

  # Add Score Probes by # of G's
  #  
  # Its lack of G::
  # gzip -dc /Users/bbarnes/Documents/Projects/manifests/methylation/HumanMethylation450_15017482_v.1.2.rs-only.csv.gz | grep G
  #
  snp_des_tib <- snp_des_tib %>% 
    dplyr::mutate(G_Count_IA=stringr::str_length(Prb_Seq_IA) - stringr::str_length(stringr::str_remove_all(Prb_Seq_IA, 'G')) ) %>% 
    dplyr::mutate(G_Count_IB=stringr::str_length(Prb_Seq_IB) - stringr::str_length(stringr::str_remove_all(Prb_Seq_IB, 'G')) ) %>% 
    dplyr::mutate(G_Count_II=stringr::str_length(Prb_Seq_II) - stringr::str_length(stringr::str_remove_all(Prb_Seq_II, 'G')) )
  
  # Score Summary::
  #  snp_des_tib %>% dplyr::arrange(G_Count_II) %>% dplyr::select(Ilmn_Id_II,Prb_Seq_II,G_Count_IA,G_Count_IB,G_Count_II) %>% print()
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Load Alignments::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

snp_pattern <- 'tile_main_EPIC_SARS-CoV-2_MN908947_COVIC.snp'
cgn_pattern <- 'tile_main_EPIC_SARS-CoV-2_MN908947_COVIC.cgn'

snp_aln_files <- list.files(opt$aln_dir, pattern=snp_pattern, full.names=TRUE)
cgn_aln_files <- list.files(opt$aln_dir, pattern=cgn_pattern, full.names=TRUE)

snp_aln_cnts <- snp_aln_files %>% length()
cgn_aln_cnts <- cgn_aln_files %>% length()
cat(glue::glue("[{par$prgmTag}]: Found snp={snp_aln_cnts}, cgn={cgn_aln_cnts}{RET}"))

if (!is.null(opt$max)) {
  snp_aln_files <- snp_aln_files %>% head(opt$max)
  cgn_aln_files <- cgn_aln_files %>% head(opt$max)
  
  snp_aln_cnts <- snp_aln_files %>% length()
  cgn_aln_cnts <- cgn_aln_files %>% length()
  cat(glue::glue("[{par$prgmTag}]: Found snp={snp_aln_cnts}, cgn={cgn_aln_cnts}{RET}"))
}

#
# Found forced SNPs::
#
art_snp_md_vec <- c("MD:Z:49T0", "MD:Z:0T49", "MD:Z:0A49", "MD:Z:49A0",
                    "MD:Z:49C0", "MD:Z:0C49", "MD:Z:0G49", "MD:Z:49G0")

all_aln_tib <- NULL
if (opt$parallel) {
  cat(glue::glue("[{par$prgmTag}]: Loading Alignments in parallel mode; num_cores={num_cores}, num_workers={num_workers}{RET}"))
  
  all_aln_tib = foreach (snp_aln_file=snp_aln_files, .combine = rbind) %dopar% {
    loadProbeAlignBowtieInfI(sam=snp_aln_file, reduced=TRUE, filtered=TRUE, flipSeq=TRUE,
                             verbose=opt$verbose,tt=pTracker)
  }

} else {
  cat(glue::glue("[{par$prgmTag}]: Loading Alignments in linear mode...{RET}"))
  
  for (snp_aln_file in snp_aln_files) {
    snp_sam_tib <- loadProbeAlignBowtieInfI(sam=snp_aln_file, reduced=TRUE, filtered=TRUE, flipSeq=TRUE,
                                            verbose=opt$verbose,tt=pTracker)
    all_aln_tib <- all_aln_tib %>% dplyr::bind_rows(snp_sam_tib)
  }
}
all_aln_len <- all_aln_tib %>% base::nrow()
cat(glue::glue("[{par$prgmTag}]: Done. Loading Alignments; all_aln_len={all_aln_len}.{RET}"))

# Add MD Class for Infinium I Probe Pairs::
#  TBD:: Add SNP Classes (A,C,T,G) from art_snp_md_vec + Original SNP; i.e. A/C, A/G, etc. 
#
art_snp_md_vec <- c("MD:Z:49T0", "MD:Z:0T49", "MD:Z:0A49", "MD:Z:49A0",
                    "MD:Z:49C0", "MD:Z:0C49", "MD:Z:0G49", "MD:Z:49G0")

ann_aln_tib <- all_aln_tib %>% dplyr::mutate(
  MD_Class_A=dplyr::case_when(
    MD_IA=="MD:Z:50" ~ 'P',
    FLAG==16 & MD_IA=='MD:Z:0A49' ~ 'A',
    FLAG==0  & MD_IA=='MD:Z:49A0' ~ 'A',
    FLAG==16 & MD_IA=='MD:Z:0C49' ~ 'C',
    FLAG==0  & MD_IA=='MD:Z:49C0' ~ 'C',
    FLAG==16 & MD_IA=='MD:Z:0G49' ~ 'G',
    FLAG==0  & MD_IA=='MD:Z:49G0' ~ 'G',
    FLAG==16 & MD_IA=='MD:Z:0T49' ~ 'T',
    FLAG==0  & MD_IA=='MD:Z:49T0' ~ 'T',
    TRUE ~ 'U'),
  MD_Class_B=dplyr::case_when(
    MD_IB=="MD:Z:50" ~ 'P',
    FLAG==16 & MD_IB=='MD:Z:0A49' ~ 'A',
    FLAG==0  & MD_IB=='MD:Z:49A0' ~ 'A',
    FLAG==16 & MD_IB=='MD:Z:0C49' ~ 'C',
    FLAG==0  & MD_IB=='MD:Z:49C0' ~ 'C',
    FLAG==16 & MD_IB=='MD:Z:0G49' ~ 'G',
    FLAG==0  & MD_IB=='MD:Z:49G0' ~ 'G',
    FLAG==16 & MD_IB=='MD:Z:0T49' ~ 'T',
    FLAG==0  & MD_IB=='MD:Z:49T0' ~ 'T',
    TRUE ~ 'U')
)

# Old Method:: Not-Specific::
# ann_aln_tib <- all_aln_tib %>% dplyr::mutate(
#   MD_Class_A=dplyr::case_when(
#     MD_IA=="MD:Z:50" ~ 'PER',
#     MD_IA %in% art_snp_md_vec ~ 'SNP',
#     TRUE ~ 'UND'),
#   MD_Class_B=dplyr::case_when(
#     MD_IB=="MD:Z:50" ~ 'PER',
#     MD_IB %in% art_snp_md_vec ~ 'SNP',
#     TRUE ~ 'UND')
# )


# Summary by class
#  ann_aln_tib %>% dplyr::select(starts_with('MD_Class')) %>% dplyr::group_by_all() %>% dplyr::summarise(CNT=n()) %>% dplyr::arrange(-CNT) %>% as.data.frame()  
#  ann_aln_tib %>% dplyr::filter(MD_Class_A=='SNP' & MD_Class_B=='SNP') %>% dplyr::select(FLAG, MD_IA, MD_IB, SEQ_IA, SEQ_IB)
#  ann_aln_tib %>% dplyr::filter(MD_Class_A=='SNP' & MD_Class_B=='PER') %>% dplyr::select(FLAG, MD_IA, MD_IB, SEQ_IA, SEQ_IB)

# Gather Stats::
#
ann_sum_tib <- dplyr::full_join(
  ann_aln_tib %>% dplyr::group_by(QNAME, MD_Class_A) %>% dplyr::summarise(MD_A_Count=n()) %>% 
    tidyr::spread(MD_Class_A, MD_A_Count, fill=0) %>% dplyr::ungroup() %>% 
    dplyr::mutate(TOT=PER+SNP+UND, PER_perc=round(100*PER/TOT,1), SNP_perc=round(100*SNP/TOT,1), UND_perc=round(100*UND/TOT,1)),
  
  ann_aln_tib %>% dplyr::group_by(QNAME, MD_Class_B) %>% dplyr::summarise(MD_B_Count=n()) %>% 
    tidyr::spread(MD_Class_B, MD_B_Count, fill=0) %>% dplyr::ungroup() %>% 
    dplyr::mutate(TOT=PER+SNP+UND, PER_perc=round(100*PER/TOT,1), SNP_perc=round(100*SNP/TOT,1), UND_perc=round(100*UND/TOT,1)),
  by="QNAME", suffix=c("_A", "_B")
)

# Filter Design Pairs Groups::
#  - Remove Underlying SNPs
#  - Select Low G Probes
#  - Select SNP Probes (SNP/PER or SNP/SNP)
#
ann_sum_tib %>% dplyr::filter(UND_A<=0 & UND_B <=0) %>% dplyr::select(-QNAME,-PER_A,-PER_B) %>% 
  dplyr::group_by_all() %>% dplyr::summarise(Count=n()) %>% dplyr::arrange(-Count) %>% as.data.frame()
  
ann_sum_tib %>% dplyr::filter(UND_A<=2 & UND_B <=2) %>% dplyr::arrange(-SNP_A)

# Add sequences back::
#  - seq_aln_tib <- all_aln_tib %>% dplyr::distinct(QNAME,SEQ_IA,SEQ_IB)
#  - seq_aln_tib+ann_sum_tib
# Join Design Data::
#  - des_aln_mat_tib <- dplyr::inner_join(seq_aln_tib, snp_des_tib, by=c("SEQ_IA"="Prb_Seq_IA", "SEQ_IB"="Prb_Seq_IB") )
#





# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Build Match Descriptor Summary::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

mds_sum_tib <- all_aln_tib %>% dplyr::group_by(QNAME,MD_IA,MD_IB) %>% dplyr::summarise(Count=n()) %>% 
  # dplyr::filter(Count!=11) %>% 
  dplyr::arrange(QNAME)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Output Merged Alignments/Summary::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sum_mds_csv <- file.path(opt$outDir, 'snp-matchDesc-merge.csv.gz')
sum_aln_csv <- file.path(opt$outDir, 'snp-alignment-merge.csv.gz')

cat(glue::glue("[{par$prgmTag}]: Writing {sum_mds_csv}..."))
readr::write_csv(mds_sum_tib, sum_mds_csv)

cat(glue::glue("[{par$prgmTag}]: Writing {sum_aln_csv}..."))
readr::write_csv(all_aln_tib, sum_aln_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Make Selectionn::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  # Max Count = 3200
  
  # Cluster Data::
  #
  datDir <- '/Users/bbarnes/Documents/Projects/methylation/scratch/analyze_tile_alignments/EPIC/SARS-CoV-2/MN908947/COVIC/n20'
  datDir <- '/Users/bbarnes/Documents/Projects/methylation/scratch/analyze_tile_alignments/EPIC/SARS-CoV-2/MN908947/COVIC/n300'
  snp_aln_csv <- file.path(datDir, 'snp-alignment-merge.csv.gz')
  snp_sum_csv <- file.path(datDir, 'snp-matchDesc-merge.csv.gz')
  
  all_aln_tib <- readr::read_csv(snp_aln_csv)
  mds_sum_tib <- readr::read_csv(snp_sum_csv)
  
  # Probe-Design/Probe-Alignment Matching::
  #
  seq_aln_tib <- all_aln_tib %>% dplyr::distinct(QNAME,SEQ_IA,SEQ_IB)
  des_aln_mat_tib <- dplyr::inner_join(seq_aln_tib, snp_des_tib, by=c("SEQ_IA"="Prb_Seq_IA", "SEQ_IB"="Prb_Seq_IB") )
  des_aln_mis_tib <- dplyr::anti_join(seq_aln_tib, snp_des_tib, by=c("SEQ_IA"="Prb_Seq_IA", "SEQ_IB"="Prb_Seq_IB") )
  
  # Alignment Analysis::
  #
  per_snp_md_vec <- c('MD:Z:50')
  art_snp_md_vec <- c("MD:Z:49T0", "MD:Z:0T49", "MD:Z:0A49", "MD:Z:49A0",
                      "MD:Z:49C0", "MD:Z:0C49", "MD:Z:0G49", "MD:Z:49G0")
  und_snp_md_vec <- c(per_snp_md_vec, art_snp_md_vec)
  
  # Get Full Genome Counts
  #
  # all_cnt_tib <- mds_sum_tib %>% dplyr::ungroup() %>% dplyr::arrange(QNAME) %>% dplyr::group_by(QNAME) %>% dplyr::summarise(Total_Count=n())
  all_cnt_tib <- all_aln_tib %>% dplyr::ungroup() %>% dplyr::select(QNAME) %>%
    dplyr::arrange(QNAME) %>% dplyr::group_by(QNAME) %>% dplyr::summarise(Total_Count=n())
  
  # Alignment Distribution
  #  all_cnt_tib %>% dplyr::group_by(Total_Count) %>% dplyr::summarise(cnt=n())
  
  # Get Sub Category Counts::
  #
  # TBD:: Add a new Class Definition
  #
  #  mds_sum_tib %>% dplyr::mutate(Count_IA_PER=dplyr::case_when( MD_IA %in% per_snp_md_vec ~ Count,TRUE ~ 0.0 ) )
  #
  # Reference Infinium IA:
  #
  per_IA_cnt_tib <- mds_sum_tib %>% dplyr::ungroup() %>% dplyr::filter(  MD_IA %in% per_snp_md_vec) %>% 
    dplyr::select(QNAME,Count) %>% dplyr::distinct() %>% dplyr::rename(Count_IA_PER=Count)
  
  art_IA_cnt_tib <- mds_sum_tib %>% dplyr::ungroup() %>% dplyr::filter(  MD_IA %in% art_snp_md_vec) %>% 
    dplyr::select(QNAME,Count) %>% dplyr::distinct() %>% dplyr::rename(Count_IA_SNP=Count)
  
  und_IA_cnt_tib <- mds_sum_tib %>% dplyr::ungroup() %>% dplyr::filter(! MD_IA %in% und_snp_md_vec) %>% 
    dplyr::select(QNAME,Count) %>% dplyr::distinct() %>% dplyr::rename(Count_IA_UND=Count)
  
  # Alternate Infinium IB:
  #
  per_IB_cnt_tib <- mds_sum_tib %>% dplyr::ungroup() %>% dplyr::filter(  MD_IB %in% per_snp_md_vec) %>% 
    dplyr::select(QNAME,Count) %>% dplyr::distinct() %>% dplyr::rename(Count_IB_PER=Count)
  
  art_IB_cnt_tib <- mds_sum_tib %>% dplyr::ungroup() %>% dplyr::filter(  MD_IB %in% art_snp_md_vec) %>% 
    dplyr::select(QNAME,Count) %>% dplyr::distinct() %>% dplyr::rename(Count_IB_SNP=Count)
  
  und_IB_cnt_tib <- mds_sum_tib %>% dplyr::ungroup() %>% dplyr::filter(! MD_IB %in% und_snp_md_vec) %>% 
    dplyr::select(QNAME,Count) %>% dplyr::distinct() %>% dplyr::rename(Count_IB_UND=Count)

  # Join the class counts and makes sure they add up::
  #
  join_cnt_tib <- all_cnt_tib %>% 
    dplyr::left_join(per_IA_cnt_tib, by="QNAME") %>%
    dplyr::left_join(art_IA_cnt_tib, by="QNAME") %>%
    dplyr::left_join(und_IA_cnt_tib, by="QNAME") %>%
    dplyr::left_join(per_IB_cnt_tib, by="QNAME") %>%
    dplyr::left_join(art_IB_cnt_tib, by="QNAME") %>% 
    dplyr::left_join(und_IB_cnt_tib, by="QNAME") %>% 
    base::replace(is.na(.), 0) %>%
    dplyr::mutate(Sum_IA=Count_IA_PER+Count_IA_SNP+Count_IA_UND, Sum_IB=Count_IB_PER+Count_IB_SNP+Count_IB_UND)
    
  
  # Generate summarys::
  #
  join_cnt_tib %>% dplyr::arrange(-Count_IA_SNP)
  
  join_cnt_tib %>% dplyr::select(-QNAME, -Total_Count) %>% dplyr::group_by_all() %>% 
    dplyr::summarise(Count=n()) %>% dplyr::arrange(-Count) %>% as.data.frame()
  
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sysTime <- Sys.time()
cat(glue::glue("{RET}[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
