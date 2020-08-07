
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
par$prgmTag <- 'tile_main'

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

opt$genome  <- NULL

# Probe Filter Parameters::
opt$minPrbScore <- 0.3
opt$minCpgRank  <- "3,2,1,0"
opt$minScrRank  <- "0.6,0.5,0.4,0.3"
opt$strandCO    <- 'C,O'
opt$pickBest    <- FALSE

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
    is.null(opt$runName) || is.null(opt$fasta) ||
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
  if (is.null(opt$fasta))     cat(glue::glue("[Usage]: fasta is NULL!!!{RET}"))

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

strandCO_vec <- opt$strandCO %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
cpgRank_vec  <- opt$minCpgRank %>% str_split(pattern=',', simplify=TRUE) %>% as.integer() %>% as.vector()
scrRank_vec  <- opt$minScrRank %>% str_split(pattern=',', simplify=TRUE) %>% as.double() %>% as.vector()
genAlign_vec <- opt$genome %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

print(genAlign_vec)

opt <- setLaunchExe(opts=opt, pars=par, verbose=opt$verbose, vt=5,tc=0)

fas_files_vec  <- opt$fasta %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

fas_dat <- Biostrings::readDNAStringSet(fas_files_vec[1])
fas_tib <- tibble::tibble(Chrom=names(fas_dat), Sequence=paste(fas_dat),
                          Seq_Length=stringr::str_length(Sequence),
                          platform=opt$platform, version=opt$version, build=opt$build)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Defined Outputs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$outDir <- file.path(opt$outDir, par$prgmTag, opt$platform, opt$version, opt$build, opt$runName)
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
cat(glue::glue("[{par$prgmTag}]: Built; OutDir={opt$outDir}!{RET}") )

if (opt$clean) list.files(opt$outDir, full.names=TRUE) %>% unlink()

out_des_str <- paste(par$prgmTag, opt$platform, opt$version, opt$build, opt$runName, sep='_')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                  Main::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$des_seq_len <- 122

opt$pre_beg_pos <- 1
opt$pre_end_pos <- opt$pre_beg_pos + 59

opt$mid_beg_pos <- opt$pre_end_pos + 1
opt$mid_end_pos <- opt$mid_beg_pos + 0

opt$pos_beg_pos <- opt$mid_end_pos + 1
opt$pos_end_pos <- opt$des_seq_len

beg_vec <- seq(from=1, to=fas_tib$Seq_Length, by=1)
end_vec <- beg_vec+opt$des_seq_len
snp_vec <- beg_vec+(opt$mid_beg_pos-1)

all_des_tib <- tibble::tibble(
  Chromosome=opt$build,
  Coordinate=snp_vec,
  Start=snp_vec,
  End=snp_vec+1,
  Strand='F',
  Forward_Sequence_Raw=stringr::str_sub(fas_tib$Sequence[1], start=beg_vec, end=end_vec),
  Forward_Sequence_Len=stringr::str_length(Forward_Sequence_Raw),
  
  Pre_Seq=stringr::str_sub(Forward_Sequence_Raw, opt$pre_beg_pos,opt$pre_end_pos),
  Ref_Nuc=stringr::str_sub(Forward_Sequence_Raw, opt$mid_beg_pos,opt$mid_end_pos),
  Nxt_Nuc=stringr::str_sub(Forward_Sequence_Raw, opt$mid_end_pos+1,opt$mid_end_pos+1),
  Alt_Nuc=dplyr::case_when( Ref_Nuc=='A' ~ 'G', Ref_Nuc=='T' ~ 'C', Ref_Nuc=='C' ~ 'T', Ref_Nuc=='G' ~ 'A', TRUE ~ NA_character_),
  Di_Nuc=paste0(Ref_Nuc,Nxt_Nuc),
  
  Pos_NUC_Seq=stringr::str_sub(Forward_Sequence_Raw, opt$pos_beg_pos,opt$pos_end_pos),
  Pos_CGN_Seq=stringr::str_sub(Forward_Sequence_Raw, opt$pos_beg_pos+1,opt$pos_end_pos),
  
  Forward_Sequence_SNP_Des=paste0(Pre_Seq,'[',Ref_Nuc,'/',Alt_Nuc,']',Pos_NUC_Seq),
  Forward_Sequence_CGN_Org=paste0(Pre_Seq,'[',Ref_Nuc,Nxt_Nuc,']',Pos_CGN_Seq),
  Forward_Sequence_CGN_Des=paste0(Pre_Seq,'[','C','G',']',Pos_CGN_Seq),
  platform=opt$platform, version=opt$version, Genome_Build=opt$build,genome=opt$runName,
  Seq_ID=paste0(Genome_Build,'.',Coordinate,'_',Di_Nuc),
  
  # Build SNP Forward Probes::
  Prb_SNP_F_IA = paste0(stringr::str_sub(Pre_Seq, stringr::str_length(Pre_Seq) - 48), Ref_Nuc ),
  Prb_SNP_F_IB = paste0(stringr::str_sub(Pre_Seq, stringr::str_length(Pre_Seq) - 48), Alt_Nuc ),
  Prb_SNP_F_II = paste0(stringr::str_sub(Pre_Seq, stringr::str_length(Pre_Seq) - 49), '' ),
  
  # Build SNP Reverse Probes::
  Prb_SNP_R_IA = revCmp( paste0(Ref_Nuc,stringr::str_sub(Pos_NUC_Seq, 1,49) ) ),
  Prb_SNP_R_IB = revCmp( paste0(Alt_Nuc,stringr::str_sub(Pos_NUC_Seq, 1,49) ) ),
  Prb_SNP_R_II = revCmp( paste0('',stringr::str_sub(Pos_NUC_Seq, 1,50) ) )
  
) %>% dplyr::filter(Forward_Sequence_Len==opt$des_seq_len+1) %>% 
  dplyr::select(Seq_ID,Genome_Build,Chromosome,Coordinate, everything())


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Build cDNA Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Read cDNA Designs::
prb_snp_fas <- NULL
prb_dna_best_csv <- NULL
prb_dna_other_csv <- NULL
# prb_dna_best_csv  <- '/Users/bbarnes/Documents/Projects/COVID-19_HLA/data/directDetection/370845_SARS-CoV-2_probes_BEST.score.csv'
# prb_dna_other_csv <- '/Users/bbarnes/Documents/Projects/COVID-19_HLA/data/directDetection/370846_SARS-CoV-2_probes_OTHER.score.csv'

if (opt$isLinux) {
  prb_dna_best_csv  <- 'designs/370845_SARS-CoV-2_probes_BEST.score.csv'
  prb_dna_other_csv <- 'designs/370846_SARS-CoV-2_probes_OTHER.score.csv'
}

if (file.exists(prb_dna_best_csv) && file.exists(prb_dna_other_csv)) {
  prb_dna_best_tib  <- suppressMessages(suppressWarnings( readr::read_csv(file=prb_dna_best_csv, skip=15) ))
  prb_dna_other_tib <- suppressMessages(suppressWarnings( readr::read_csv(file=prb_dna_other_csv, skip=15) ))
  
  # Join Score Files::
  prb_dna_tib <- dplyr::bind_rows(prb_dna_best_tib,prb_dna_other_tib) %>% 
    dplyr::select(Ilmn_Id,Locus_Name,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Assay_Type,Final_Score,Sequence) %>% dplyr::distinct()

  # Merge Score Data with Full Design Data::
  #
  prb_dna_tib %>% dplyr::inner_join(all_des_tib, by=c("Sequence"="Forward_Sequence_SNP_Des")) %>%
    dplyr::filter(AlleleA_Probe_Sequence==Prb_SNP_F_IA & AlleleB_Probe_Sequence==Prb_SNP_F_IB) 
  
  prb_dna_tib %>% dplyr::inner_join(all_des_tib, by=c("Sequence"="Forward_Sequence_SNP_Des")) %>%
    dplyr::filter(AlleleA_Probe_Sequence==Prb_SNP_R_IA & AlleleB_Probe_Sequence==Prb_SNP_R_IB) 
  
  snp_des_tib <- prb_dna_tib %>% dplyr::inner_join(all_des_tib, by=c("Sequence"="Forward_Sequence_SNP_Des")) %>%
    dplyr::filter(
      (AlleleA_Probe_Sequence==Prb_SNP_F_IA & AlleleB_Probe_Sequence==Prb_SNP_F_IB) |
        (AlleleA_Probe_Sequence==Prb_SNP_R_IA & AlleleB_Probe_Sequence==Prb_SNP_R_IB) ) %>%
    dplyr::mutate(Ilmn_Id=stringr::str_replace(Ilmn_Id,'_','.'))
  
  # Now we have scores and Sequences::
  #  - [DONE] Validate uniqueness
  #  - Generate new probe names to match CGN
  #  - Write Fasta
  
  #  - Align (bowtie) to hg19,hg38,COVID,COVID-NCBI
  #  - Use match descriptor from COVID-NCBI alignments to:
  #    - Remove probes with underlying SNPs
  #    - Identify probes with target SNPs
  #
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Write Fasta File::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  prb_snp_tib <- dplyr::bind_rows(
    snp_des_tib %>% dplyr::select(Ilmn_Id,AlleleA_Probe_Sequence) %>% dplyr::mutate(Fas_Id=paste(Ilmn_Id,'IA', sep='_'), line=paste0('>',Fas_Id,'\n',AlleleA_Probe_Sequence) ),
    snp_des_tib %>% dplyr::select(Ilmn_Id,AlleleB_Probe_Sequence) %>% dplyr::mutate(Fas_Id=paste(Ilmn_Id,'IB', sep='_'), line=paste0('>',Fas_Id,'\n',AlleleB_Probe_Sequence) )
  ) %>% dplyr::select(Ilmn_Id, line) %>% dplyr::arrange(Ilmn_Id)
  
  prb_snp_str <- prb_snp_tib %>% dplyr::pull(line)
  
  # This needs to be split, but we can ignore it for now because alignment 
  #  can be captured with Infinium I probes, which are the only ones we have
  #  scores for at the moment anyways...
  # snp_des_tib %>% dplyr::select(Ilmn_Id,Prb_SNP_F_II,Prb_SNP_R_II) %>% 
  #   dplyr::mutate(Ilmn_Id=paste(Ilmn_Id,'II', sep='_'), 
  #                 lineF=paste0('>',Ilmn_Id,'\n',Prb_SNP_F_II),
  #                 lineR=paste0('>',Ilmn_Id,'\n',Prb_SNP_R_II) ) %>% 
  #   dplyr::pull(lineF,lineR)

  prb_snp_fas  <- file.path(opt$outDir, paste0(out_des_str,'.snp.fa.gz'))
  readr::write_lines(x=prb_snp_str, path=prb_snp_fas)
  bow_name <- prb_snp_fas %>% stringr::str_remove('.fa.gz$')
  
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Build CGN Tibbles::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cgn_des_tib <- all_des_tib %>% dplyr::select(Seq_ID,Forward_Sequence_CGN_Des,Genome_Build,Chromosome,Coordinate) %>%
  dplyr::rename(Sequence=Forward_Sequence_CGN_Des) %>%
  dplyr::mutate(CpG_Island="FALSE")

cgn_org_tib <- all_des_tib %>% dplyr::select(Seq_ID,Forward_Sequence_CGN_Org,Genome_Build,Chromosome,Coordinate) %>%
  dplyr::rename(Sequence=Forward_Sequence_CGN_Org) %>%
  dplyr::mutate(CpG_Island="FALSE")

if (opt$verbose>4) {
  all_des_tib$Forward_Sequence_Raw %>% head(n=3) %>% print()
  all_des_tib$Forward_Sequence_SNP_Des %>% head(n=3) %>% print()
  all_des_tib$Forward_Sequence_CGN_Org %>% head(n=3) %>% print()
  all_des_tib$Forward_Sequence_CGN_Des %>% head(n=3) %>% print()
  
  all_des_tib %>% head(n=3) %>% as.data.frame()
  
  all_des_tib %>% dplyr::select(Seq_ID, Ref_Nuc, Nxt_Nuc, Di_Nuc, Forward_Sequence_Raw) %>% tail()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Write Outputs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

all_des_csv <- file.path(opt$outDir, paste(out_des_str, 'tiled-all-dat.csv.gz', sep='_') )
cgn_des_tsv <- file.path(opt$outDir, paste(out_des_str, 'tiled-cgn-des.tsv.gz', sep='_') )
cgn_org_tsv <- file.path(opt$outDir, paste(out_des_str, 'tiled-cgn-org.tsv.gz', sep='_') )

readr::write_csv(all_des_tib,all_des_csv)
readr::write_tsv(cgn_des_tib,cgn_des_tsv)
readr::write_tsv(cgn_org_tib,cgn_org_tsv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Run improbe design::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imp_out_tsv <- file.path(opt$outDir, paste(out_des_str, 'tiled-cgn-improbe-designs.tsv.gz', sep='_') )
imp_out_log <- file.path(opt$outDir, paste(out_des_str, 'tiled-cgn-improbe-designs.log', sep='_') )

if (opt$isLinux) {
  if (! file.exists(par$tan_file))
    stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: On linux and can't find tangos={par$tan_file}!{RET}{RET}"))
  if (! file.exists(par$mer_file))
    stop(glue::glue("{RET}[{par$prgmTag}]: ERROR: On linux and can't find 13-mer={par$mer_file}!{RET}{RET}"))

  shell_dir <- file.path(opt$outDir, 'shells')
  if (!dir.exists(shell_dir)) dir.create(shell_dir, recursive=TRUE)
  shell_file <- file.path(shell_dir, 'run_improbe.sh')
  
  cmd <- paste(
    'gzip -dc',cgn_des_tsv,'|',
    par$improbe_exe,
    '-oASPE -tBOTH -cBoth',
    '-n', par$mer_file,
    '-a', par$tan_file,
    '-V - 2>',imp_out_log,
    '| gzip -c - >',imp_out_tsv,
    sep=' '
  )
  
  readr::write_lines(x=cmd, path=shell_file, append=FALSE)
  Sys.chmod(paths=shell_file, mode="0777")
  base::system(shell_file)
}

if (!is.null(imp_out_tsv) & file.exists(imp_out_tsv)) {
  # Load designs and extract original di-nucleotide
  des_ord_tib <- suppressMessages(suppressWarnings( readr::read_tsv(imp_out_tsv) )) %>% 
    dplyr::mutate(Di_Nuc=stringr::str_to_upper( stringr::str_remove(Seq_ID,'^.*_') ),
                  Probe_Type=dplyr::case_when(Di_Nuc=='CG' ~ 'cg', TRUE ~ 'ch' ) )

  des_scr_tib <- des_ord_tib %>% dplyr::select(Seq_ID,Forward_Sequence,Genome_Build,Chromosome,Coordinate,
                                UnMethyl_Probe_Sequence,Methyl_Probe_Sequence, Probe_Type,
                                UnMethyl_Final_Score,Methyl_Final_Score,Methyl_Underlying_CpG_Count,
                                Methyl_Allele_FR_Strand,Methyl_Allele_TB_Strand,Methyl_Allele_CO_Strand) %>% 
    dplyr::rename(Prb_Seq_IU_IMP=UnMethyl_Probe_Sequence,Prb_Seq_IM_IMP=Methyl_Probe_Sequence,
                  Prb_Scr_IU=UnMethyl_Final_Score,Prb_Scr_IM=Methyl_Final_Score,CpgCnt=Methyl_Underlying_CpG_Count,
                  FR_Str=Methyl_Allele_FR_Strand,TB_Str=Methyl_Allele_TB_Strand,CO_Str=Methyl_Allele_CO_Strand) %>%
    dplyr::mutate(TB_Str=stringr::str_sub(TB_Str,1,1), Probe_ID=paste(Seq_ID, FR_Str,TB_Str,CO_Str, sep='_'), Prb_Scr_Min=pmin(Prb_Scr_IU,Prb_Scr_IM)) %>%
    dplyr::select(Seq_ID,Probe_ID,everything()) %>%
    dplyr::filter(Prb_Scr_Min>=opt$minPrbScore) %>% 
    dplyr::filter(CO_Str %in% strandCO_vec) %>% 
    dplyr::filter(Prb_Scr_Min>=opt$minPrbScore) %>%
    dplyr::mutate(Design_Type=dplyr::case_when(
      CpgCnt==cpgRank_vec[1] & Prb_Scr_Min>=scrRank_vec[1] ~ 'II',
      CpgCnt==cpgRank_vec[2] & Prb_Scr_Min>=scrRank_vec[2] ~ 'II',
      CpgCnt==cpgRank_vec[3] & Prb_Scr_Min>=scrRank_vec[3] ~ 'II',
      CpgCnt==cpgRank_vec[1] & Prb_Scr_Min>=scrRank_vec[1] ~ 'II',
      TRUE ~ 'I'
    )) %>% dplyr::arrange(-Prb_Scr_Min) %>% dplyr::distinct(Seq_ID, .keep_all=TRUE)
  
  des_scr_tib %>% dplyr::group_by(Design_Type) %>% dplyr::summarise(Count=n()) %>% print()
  
  if (!opt$pickBest) {
    des_scr_tib <- des_scr_tib %>% dplyr::filter(Design_Type=='II') %>% 
      dplyr::mutate(Design_Type='I') %>% dplyr::bind_rows(des_scr_tib) %>% dplyr::arrange(Seq_ID)
  }
  des_scr_tib %>% dplyr::group_by(Design_Type) %>% dplyr::summarise(Count=n()) %>% print()
  
  des_fwd_tib <- des_scr_tib %>% dplyr::select(Seq_ID, Forward_Sequence, Probe_Type) %>%
    dplyr::distinct(Seq_ID, .keep_all=TRUE)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Design All Probes::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  new_prb_tib <- tib2prbs(tib=des_fwd_tib, idsKey="Seq_ID", prbKey="Probe_Type", seqKey="Forward_Sequence", verbose=opt$verbose)
  new_prb_tib %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                  Match Probes by Strand to Top Picks::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  sel_prb_tib <- dplyr::inner_join(new_prb_tib,des_scr_tib, by=c("Seq_ID","FR_Str","CO_Str","Probe_Type","Forward_Sequence") ) 
  sel_prb_tib %>% dplyr::group_by(Design_Type) %>% dplyr::summarise(Count=n()) %>% print()

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Format Order File::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  new_ord1_tib <- prbs2order(sel_prb_tib, verbose=opt$verbose) %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(Valid_Design_Bool) %>% 
    dplyr::select(Assay_Design_Id:Normalization_Bin) %>% dplyr::filter(Normalization_Bin!='C')
  
  new_ord2_tib <- prbs2order(sel_prb_tib, verbose=opt$verbose) %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(Valid_Design_Bool) %>% 
    dplyr::select(Assay_Design_Id:Normalization_Bin) %>% dplyr::filter(Normalization_Bin=='C')
  
  new_ord_tib <- dplyr::bind_rows(new_ord1_tib,new_ord2_tib) %>%
    dplyr::arrange(Assay_Design_Id) %>%
    dplyr::mutate(
      AlleleB_Probe_Id=dplyr::case_when(
        is.na(AlleleB_Probe_Id) ~ '',
        TRUE ~ AlleleB_Probe_Id
      ),
      AlleleB_Probe_Sequence=dplyr::case_when(
        is.na(AlleleB_Probe_Sequence) ~ '',
        TRUE ~ AlleleB_Probe_Sequence
      )
    ) %>% dplyr::distinct(AlleleA_Probe_Id,AlleleA_Probe_Sequence, .keep_all=TRUE)
  
  order2stats(new_ord_tib) %>% print()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         Write Final Order File::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  new_ord_csv <- file.path(opt$outDir, paste(opt$runName,'order.csv.gz', sep='.'))
  readr::write_csv(new_ord_tib, new_ord_csv)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Write Fasta File::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

  prb_cgn_fas  <- ordToFas(tib=new_ord_tib, dir=opt$outDir, name=out_des_str, verbose=opt$verbose)
  bsp_name <- prb_cgn_fas %>% stringr::str_remove('.fa.gz$')

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Build Alignments::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  gen_cnt <- length(genAlign_vec)
  for (gen_idx in c(1:gen_cnt)) {
    gen_path <- genAlign_vec[gen_idx]
    cat(glue::glue("[{par$prgmTag}]: Ready to launch alignments; genomes count={gen_cnt}; gen_path={gen_path}...{RET}"))
    
    if (dir.exists(gen_path)) {
      cat(glue::glue("[{par$prgmTag}]:{TAB} Searching directory for genomes: gen_path={gen_path}...{RET}"))
      # gen_paths <- list.files(gen_path, pattern='.fa[.gz]$', full.names=TRUE)
      gen_paths <- list.files(gen_path, pattern='.fa.gz$', full.names=TRUE)
      
      gen_cnts <- 0
      for (gpath in gen_paths) {
        if (file.exists(gpath)) {
          cat(glue::glue("[{par$prgmTag}]:{TAB}{TAB}Ready to launch alignments; genomes count={gen_cnts}/{gen_cnt}; gen_path={gen_path}...{RET}"))
          
          bow_tsv <- bowtieProbeAlign(exe=par$bow_exe, fas=prb_snp_fas, gen=gpath, dir=opt$outDir,
                                      verbose=opt$verbose,vt=5,tc=1,tt=pTracker)
          
          bsp_tsv <- bsmapProbeAlign(exe=par$bsp_exe, fas=prb_cgn_fas, gen=gpath, dir=opt$outDir,
                                     verbose=opt$verbose,vt=5,tc=1,tt=pTracker)
        }
        
        gen_cnts <- gen_cnts + 1
        # if (gen_cnts > 10) break
      }
    } else if (file.exists(gen_path)) {
      cat(glue::glue("[{par$prgmTag}]:{TAB} Using Genome Path={gen_path}...{RET}"))
      bow_tsv <- bowtieProbeAlign(exe=par$bow_exe, fas=prb_snp_fas, gen=gen_path, dir=opt$outDir,
                                  verbose=opt$verbose,vt=5,tc=1,tt=pTracker)
      
      bsp_tsv <- bsmapProbeAlign(exe=par$bsp_exe, fas=prb_cgn_fas, gen=gen_path, dir=opt$outDir,
                                 verbose=opt$verbose,vt=5,tc=1,tt=pTracker)

    }
    

  }
  cat(glue::glue("[{par$prgmTag}]: Done launching alignments.{RET}{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$opt_csv  <- file.path(opt$outDir, paste(par$prgmTag,'program-options.csv', sep='.') )
opt$par_csv  <- file.path(opt$outDir, paste(par$prgmTag,'program-parameters.csv', sep='.') )
opt$time_csv <- file.path(opt$outDir, paste(par$prgmTag,'time-tracker.csv.gz', sep='.') )

opt_tib  <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
par_tib  <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
# time_tib <- pTracker$time %>% dplyr::mutate_if(is.numeric, list(round), 4)

readr::write_csv(opt_tib, opt$opt_csv)
readr::write_csv(par_tib, opt$par_csv)
# readr::write_csv(time_tib, opt$time_csv)

sysTime <- Sys.time()
cat(glue::glue("{RET}[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
