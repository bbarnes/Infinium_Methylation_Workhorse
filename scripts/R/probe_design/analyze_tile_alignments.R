
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
  
  aln_dir <- '/Users/bbarnes/Documents/Projects/methylation/scratch/small.index/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/align'

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
#                               Defined Outputs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$outDir <- file.path(opt$outDir, par$prgmTag, opt$platform, opt$version, opt$build, opt$runName)
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
cat(glue::glue("[{par$prgmTag}]: Built; OutDir={opt$outDir}!{RET}") )

if (opt$clean) list.files(opt$outDir, full.names=TRUE) %>% unlink()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Load Probe Designs:: cDNA
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  snp_des1_csv <- '~/Downloads/370992_SARS-CoV-2_probes_F2BT_BEST.score.csv.gz'
  snp_des2_csv <- '~/Downloads/370986_SARS-CoV-2_probes_BEST.score.csv.gz'
  
  snp_des1_tib <- readr::read_csv(snp_des1_csv, skip=15)
  snp_des2_tib <- readr::read_csv(snp_des2_csv, skip=15)
  
  snp_des_tib <- dplyr::inner_join(
    snp_des1_tib %>% dplyr::select(Locus_Name,Coordinate,Sequence_Orientation,Final_Score,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence) %>%
      dplyr::rename(Prb_Seq_IA=AlleleA_Probe_Sequence,Prb_Seq_IB=AlleleB_Probe_Sequence),
    snp_des2_tib %>% dplyr::select(Locus_Name,Coordinate,Sequence_Orientation,Final_Score,AlleleA_Probe_Sequence) %>%
      dplyr::rename(Prb_Seq_II=AlleleA_Probe_Sequence),
    by=c("Coordinate","Sequence_Orientation"), suffix=c("_I","_II")
  ) %>% 
    dplyr::rename(Strand_FR=Sequence_Orientation,Score_I=Final_Score_I,Score_II=Final_Score_II) %>%
    dplyr::mutate(Strand_FR=case_when(Strand_FR=="FORWARD" ~ 'F', Strand_FR=="REVERSE" ~ 'R', TRUE ~ NA_character_)) %>%
    dplyr::select(Locus_Name_I,Locus_Name_II,Coordinate,Strand_FR,Prb_Seq_IA,Prb_Seq_IB,Prb_Seq_II,Score_I,Score_II, everything())
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

sam_col_vec <- c('QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL',
                 'AS', 'XN', 'XM', 'XO', 'XG', 'NM', 'MD', 'YT')

# snp_covid_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file = opt$snp_covid_sam, col_names=sam_col_vec, comment='@') ))
# snp_ncbi_tib  <- suppressMessages(suppressWarnings( readr::read_tsv(file = opt$snp_ncbi_sam, col_names=sam_col_vec, comment='@') ))
#
# snp_covid_tib %>% dplyr::group_by(CIGAR) %>% dplyr::summarise(Count=n())
# snp_ncbi_tib %>% dplyr::group_by(CIGAR) %>% dplyr::summarise(Count=n())

# Its lack of G::
# gzip -dc /Users/bbarnes/Documents/Projects/manifests/methylation/HumanMethylation450_15017482_v.1.2.rs-only.csv.gz | grep G
#
# snp_ncbi_tib %>% dplyr::mutate(G_Count=50 - stringr::str_length(stringr::str_remove_all(SEQ, 'G')) ) %>% dplyr::arrange(G_Count) %>% dplyr::select(QNAME,SEQ,G_Count)

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

sum_aln_tib <- NULL
for (snp_aln_file in snp_aln_files) {
  snp_raw_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file=snp_aln_file, col_names=sam_col_vec, comment='@') )) %>% 
    dplyr::select(QNAME:POS,SEQ,MD) %>%
    dplyr::filter(FLAG==0 | FLAG==16) %>%
    dplyr::mutate(SEQ=case_when(
      FLAG==16 ~ revCmp(SEQ),
      TRUE ~ SEQ
    ))
  
  # Split by Infinium I Probe Design::
  snp_sam_1A_tib <- snp_raw_tib %>% dplyr::filter(stringr::str_ends(QNAME,'_IA')) %>% dplyr::mutate(QNAME=stringr::str_remove(QNAME, '_IA$'))
  snp_sam_1B_tib <- snp_raw_tib %>% dplyr::filter(stringr::str_ends(QNAME,'_IB')) %>% dplyr::mutate(QNAME=stringr::str_remove(QNAME, '_IB$'))
  
  # Skipping Probe Design Join Step for Now::
  snp_sam_tib <- dplyr::inner_join(snp_sam_1A_tib,snp_sam_1B_tib, by=c("QNAME","FLAG","RNAME","POS"), suffix=c("_IA", "_IB")) # %>%
  #  dplyr::inner_join(snp_des_tib, by=c("SEQ_IA"="Prb_Seq_IA", "SEQ_IB"="Prb_Seq_IB") )
  
  MD_IA_cnt <- snp_sam_tib %>% dplyr::filter(MD_IA %in% art_snp_md_vec) %>% base::nrow()
  MD_IB_cnt <- snp_sam_tib %>% dplyr::filter(MD_IB %in% art_snp_md_vec) %>% base::nrow()
  cat(glue::glue("[{par$prgmTag}]: IA/IB = {MD_IA_cnt}, {MD_IB_cnt}{RET}"))
  
  # View SNP probes::
  # snp_sam_tib %>% dplyr::filter(MD_IA %in% art_snp_md_vec) %>% dplyr::select(QNAME,FLAG,RNAME,POS, MD_IA,MD_IB) %>% print()
  
  # Bind all alignments::
  sum_aln_tib <- sum_aln_tib %>% dplyr::bind_rows(snp_sam_tib)
}

mds_sum_tib <- sum_aln_tib %>% dplyr::group_by(QNAME,MD_IA,MD_IB) %>% dplyr::summarise(Count=n()) %>% dplyr::filter(Count!=11) %>% dplyr::arrange(QNAME)

sum_mds_csv <- file.path(opt$outDir, 'snp-matchDesc-merge.csv.gz')
sum_aln_csv <- file.path(opt$outDir, 'snp-alignment-merge.csv.gz')

cat(glue::glue("[{par$prgmTag}]: Writing {sum_mds_csv}..."))
readr::write_csv(mds_sum_tib, sum_mds_csv)

cat(glue::glue("[{par$prgmTag}]: Writing {sum_aln_csv}..."))
readr::write_csv(sum_aln_tib, sum_aln_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sysTime <- Sys.time()
cat(glue::glue("{RET}[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
