
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

run_tib <- NULL

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
    dplyr::mutate(G_Count_II=stringr::str_length(Prb_Seq_II) - stringr::str_length(stringr::str_remove_all(Prb_Seq_II, 'G')) ) %>%
    dplyr::mutate(G_Count_Max=pmax(G_Count_IA,G_Count_IB,G_Count_II))

  org_des_cnt <- snp_des_tib %>% base::nrow()
  cat(glue::glue("[{par$prgmTag}]: Design Studio: org_des_cnt={org_des_cnt}{RET}"))
  # Filter on Probe (Pair) Score::
  #
  snp_des_tib <- snp_des_tib %>% dplyr::filter(
    (Final_Score_I >= 0.3 & Final_Score_II >= 0.3 & G_Count_Max==1) |
    (Final_Score_I >= 0.4 & Final_Score_II >= 0.4 & G_Count_Max>1)
  )
    
  scr_des_cnt <- snp_des_tib %>% base::nrow()
  cat(glue::glue("[{par$prgmTag}]: Design Studio: scr_des_cnt={scr_des_cnt}{RET}"))
  
  # Score Summary::
  #  snp_des_tib %>% dplyr::arrange(G_Count_II) %>% dplyr::select(Ilmn_Id_II,Prb_Seq_II,G_Count_IA,G_Count_IB,G_Count_II) %>% print()
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Format and Summarize Alignments::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# reduced_col <- c('Ref_Genome','Ref_Pos','Des_Type','TB','FR','Des_Order','Probe_Type','ILMN_Id','FLAG','Aln_Type')
# reduced_tsv <- '/Users/bbarnes/Documents/Projects/methylation/scratch/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/n20000/aln/reduced-bowtie.50000000.tsv.gz'

reduced_col <- c('ILMN_Id', 'Ref_Genome','Ref_Pos','Des_Type','TB','FR','Des_Order','Probe_Type','FLAG','Aln_Type','Uniq_Count')
# reduced_tsv <- '/Users/bbarnes/Documents/Projects/methylation/scratch/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/n20000/aln/align-summay.bowtie.40000.tsv.gz'
reduced_tsv <- '/Users/bbarnes/Documents/Projects/methylation/scratch/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/n20000/cluster/align-summay.bowtie.50000000.tsv.gz'
# reduced_tsv <- '/Users/bbarnes/Documents/Projects/methylation/scratch/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/n20000/cluster/align-summay.bowtie.100000000.tsv.gz'
# reduced_tsv <- '/Users/bbarnes/Documents/Projects/methylation/scratch/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/n20000/cluster/align-summay.bowtie.1000000000.tsv.gz'
# reduced_tsv <- '/Users/bbarnes/Documents/Projects/methylation/scratch/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC/n20000/cluster/align-summay.bowtie.10000000000.tsv.gz'

reduced_tib <- suppressMessages(suppressWarnings( readr::read_tsv(reduced_tsv, col_names=reduced_col) ))

# Only look at Infinium I A/B for now...
# Also skip; 
#   IA_0_NM:i:0,
#   IA_0_NM:i:1
#  red_sum_tib <- reduced_tib %>% dplyr::mutate(Aln_Type=dplyr::case_when(stringr::str_starts(Aln_Type, 'NM:i:') ~ 'X', TRUE ~ Aln_Type )) %>%
red_sum_tib <- reduced_tib %>% 
  dplyr::filter(!(!is.na(Aln_Type) & stringr::str_starts(Aln_Type, 'NM:')) ) %>%
  dplyr::filter(Probe_Type!='II') %>%
  dplyr::mutate(Aln_Type=dplyr::case_when(stringr::str_starts(Aln_Type, 'NM:i:') ~ 'X', TRUE ~ Aln_Type )) %>%
  dplyr::mutate(Prb_Key=stringr::str_remove(ILMN_Id, '_I[IAB]$')) %>%
  dplyr::mutate(Aln_Type=dplyr::case_when(
    FLAG!=0 & FLAG!=16 ~ 'N',
    is.na(Aln_Type) ~ 'N',
    Aln_Type != 'A' & Aln_Type != 'C' & Aln_Type != 'G' & Aln_Type != 'T' & Aln_Type != 'P' & Aln_Type != 'U' & Aln_Type != 'X' ~ 'Z',
    TRUE ~ Aln_Type)
  ) %>% # dplyr::group_by(Aln_Type) %>% dplyr::summarise(Count=n()) %>%
  tidyr::unite(Aln_Key, Probe_Type,Aln_Type, sep='_') %>% 
  dplyr::select(Prb_Key, Aln_Key, Uniq_Count) %>% # dplyr::group_by(Aln_Key) %>% dplyr::summarise(Count=n())
  tidyr::spread(Aln_Key, Uniq_Count) %>%
  base::replace(is.na(.), 0)
run_tib$red_sum <- red_sum_tib %>% base::nrow()

red_tot_vec <- red_sum_tib %>% tibble::column_to_rownames('Prb_Key') %>% as.matrix() %>% matrixStats::rowSums2()  

#
# Add Stat for summed mutation needed for filtering::
#  IA_0_S = IA_0_A+IA_0_C+IA_0_G+IA_0_T
#  TA_Max = max(IA_0_A, IA_0_C, IA_0_G, IA_0_T)
#  TA_Max_Per = TA_Max/Aln_Tot_Cnt
#  TA_Alt_Per = 1 - TA_Max_Per
#

# TBD:: Need to builld a fake tibble to join with at the end that includes all columns...


red_ann_tib <- red_sum_tib %>%
  dplyr::mutate(
    Tot_Cnt=red_tot_vec,
    
    # Non Calls::
    NON_SumA=IA_N,
    NON_SumB=IB_N,
    
    NON_Max=pmax(NON_SumA,NON_SumB),
    NON_Min=pmin(NON_SumA,NON_SumB),
    NON_Sum=NON_SumA+NON_SumB,

    NON_PercA=round(100*NON_SumA/(Tot_Cnt/2), 1),
    NON_PercB=round(100*NON_SumB/(Tot_Cnt/2), 1),
    NON_PercM=round(100*NON_Max/(Tot_Cnt/2), 1),
    NON_PercS=round((100*NON_Max+NON_Min)/(Tot_Cnt/2), 1),
    
    NON_Max_med=median(NON_Max), NON_Max_avg=mean(NON_Max), NON_Max_sd=sd(NON_Max),
    NON_Min_med=median(NON_Min), NON_Min_avg=mean(NON_Min), NON_Min_sd=sd(NON_Min),
    NON_Sum_med=median(NON_Sum), NON_Sum_avg=mean(NON_Sum), NON_Sum_sd=sd(NON_Sum),
    
    # Und Calls::
    #
    UND_SumA=IA_U,
    UND_SumB=IB_U,
    
    UND_Max=pmax(UND_SumA,UND_SumB),
    UND_Min=pmin(UND_SumA,UND_SumB),
    UND_Sum=UND_SumA+UND_SumB,
    
    UND_PercA=round(100*UND_SumA/(Tot_Cnt/2), 1),
    UND_PercB=round(100*UND_SumB/(Tot_Cnt/2), 1),
    UND_PercM=round(100*UND_Max/(Tot_Cnt/2), 1),
    UND_PercS=round((100*UND_Max+UND_Min)/(Tot_Cnt/2), 1),
    
    UND_Max_med=median(UND_Max), UND_Max_avg=mean(UND_Max), UND_Max_sd=sd(UND_Max),
    UND_Min_med=median(UND_Min), UND_Min_avg=mean(UND_Min), UND_Min_sd=sd(UND_Min),
    UND_Sum_med=median(UND_Sum), UND_Sum_avg=mean(UND_Sum), UND_Sum_sd=sd(UND_Sum),
    
    # ALt Calls::
    #
    SNP_SumA=IA_A+IA_C+IA_G+IA_T,
    SNP_SumB=IB_A+IB_C+IB_G+IB_T,
    
    SNP_Max=pmax(SNP_SumA,SNP_SumB),
    SNP_Min=pmin(SNP_SumA,SNP_SumB),
    SNP_Sum=SNP_SumA+SNP_SumB,
    
    SNP_PercA=round(100*SNP_SumA/(Tot_Cnt/2), 1),
    SNP_PercB=round(100*SNP_SumB/(Tot_Cnt/2), 1),
    SNP_PercM=round(100*SNP_Max/(Tot_Cnt/2), 1),
    SNP_PercS=round((100*SNP_Max+SNP_Min)/(Tot_Cnt/2), 1),
    
    SNP_Max_med=median(SNP_Max), SNP_Max_avg=mean(SNP_Max), SNP_Max_sd=sd(SNP_Max),
    SNP_Min_med=median(SNP_Min), SNP_Min_avg=mean(SNP_Min), SNP_Min_sd=sd(SNP_Min),
    SNP_Sum_med=median(SNP_Sum), SNP_Sum_avg=mean(SNP_Sum), SNP_Sum_sd=sd(SNP_Sum),
    
    # Ref Calls::
    #
    REF_SumA=IA_P,
    REF_SumB=IB_P,
    
    REF_Max=pmax(REF_SumA,REF_SumB),
    REF_Min=pmin(REF_SumA,REF_SumB),
    REF_Sum=REF_SumA+REF_SumB,
    
    REF_PercA=round(100*REF_SumA/(Tot_Cnt/2), 1),
    REF_PercB=round(100*REF_SumB/(Tot_Cnt/2), 1),
    REF_PercM=round(100*REF_Max/(Tot_Cnt/2), 1),
    REF_PercS=round((100*REF_Max+REF_Min)/(Tot_Cnt/2), 1),

    REF_Max_med=median(REF_Max), REF_Max_avg=mean(REF_Max), REF_Max_sd=sd(REF_Max),
    REF_Min_med=median(REF_Min), REF_Min_avg=mean(REF_Min), REF_Min_sd=sd(REF_Min),
    REF_Sum_med=median(REF_Sum), REF_Sum_avg=mean(REF_Sum), REF_Sum_sd=sd(REF_Sum),

    Tot_Cnt=red_tot_vec
  )

#
# Real filtering methods::
#

#  NON_Max > NON_Max_avg+NON_Max_sd
opt$NON_val <- 0.1
opt$NON_val <- 1.0
rem_non_tib <- red_ann_tib %>% dplyr::filter(NON_Max >= NON_Max_avg+(opt$NON_val*NON_Max_sd) ) %>% 
  dplyr::select(Prb_Key,Tot_Cnt,NON_Max,NON_Max_med,NON_Max_avg,NON_Max_sd )
# red_ann_tib %>% dplyr::filter(NON_Min <= pmax(0,NON_Min_avg-(opt$NON_val*NON_Min_sd) ) ) %>% dplyr::select(Prb_Key,NON_Min,NON_Min_med,NON_Min_avg,NON_Min_sd )
run_tib$rem_non <- rem_non_tib %>% base::nrow()

#  UND_Max[A|B] > UND_Max_avg+UND_Max_sd
opt$UND_val <- 0.1
opt$UND_val <- 1.0
rem_und_tib <- red_ann_tib %>% dplyr::filter(UND_Max >= UND_Max_avg+(opt$UND_val*UND_Max_sd) ) %>% 
  dplyr::select(Prb_Key,Tot_Cnt,UND_Max,UND_Max_med,UND_Max_avg,UND_Max_sd )
# red_ann_tib %>% dplyr::filter(UND_Min <= pmax(0,UND_Min_avg-(opt$UND_val*UND_Min_sd) ) ) %>% dplyr::select(Prb_Key,UND_Min,UND_Min_med,UND_Min_avg,UND_Min_sd )
run_tib$rem_und <- rem_und_tib %>% base::nrow()

rem_all_tib <- dplyr::bind_rows(
  dplyr::distinct(rem_non_tib,Prb_Key),
  dplyr::distinct(rem_und_tib,Prb_Key) ) %>% 
  dplyr::group_by(Prb_Key) %>% dplyr::summarise(Count=n())
run_tib$rem_all <- rem_all_tib %>% base::nrow()


#
# Selection::
#

#  Ref_Max > Ref_Max_avg+Aln_Ref_sd
#  Ref_Min < Ref_Min_avg-Aln_Ref_sd

opt$ref_val <- 0.1
opt$ref_val <- 0.3
opt$ref_val <- 0.4

sel_ref_tib <- red_ann_tib %>% dplyr::filter(REF_Max >= REF_Max_avg+(opt$ref_val*REF_Max_sd) ) %>% 
  dplyr::filter(REF_SumA==0 | REF_SumB==0) %>%
  dplyr::select(Prb_Key,Tot_Cnt,REF_SumA,REF_SumB,REF_Max,REF_Max_med,REF_Max_avg,REF_Max_sd ) %>% dplyr::arrange(-REF_Max)
run_tib$sel_ref <- sel_ref_tib %>% base::nrow()

#
# TBD:: Add REF_Ratio::
#


#
#  Aln_SNP_Perc[A|B] > Aln_SNP_avg+Aln_SNP_sd
#

# You really want the SNP_Sum vs. REF_Sum::
#
opt$SNP_val <- 0.4
sel_SNP_Max_tib <- red_ann_tib %>% dplyr::filter(SNP_Max >= SNP_Max_avg+(opt$SNP_val*SNP_Max_sd) ) %>% 
  dplyr::select(Prb_Key,Tot_Cnt,SNP_SumA,SNP_SumB,SNP_Max,SNP_Max_med,SNP_Max_avg,SNP_Max_sd ) %>% dplyr::arrange(-SNP_Max)
run_tib$sel_SNP_Max <- sel_SNP_Max_tib %>% base::nrow()

opt$SNP_val <- 0.4
sel_SNP_Sum_tib <- red_ann_tib %>% dplyr::filter(SNP_Sum >= SNP_Sum_avg+(opt$SNP_val*SNP_Sum_sd) ) %>% 
  dplyr::select(Prb_Key,Tot_Cnt,SNP_SumA,SNP_SumB,
                SNP_Sum,SNP_Sum_med,SNP_Sum_avg,SNP_Sum_sd,
                SNP_Max,SNP_Max_med,SNP_Max_avg,SNP_Max_sd,REF_Sum ) %>% 
  dplyr::mutate(SNP_Ratio=SNP_Sum/REF_Sum) %>% 
  dplyr::arrange(-SNP_Ratio) %>% 
  dplyr::select(SNP_Ratio, everything())

sel_SNP_Sum_tib <- sel_SNP_Sum_tib %>% dplyr::filter(SNP_Ratio>1)
run_tib$sel_SNP_Sum <- sel_SNP_Sum_tib %>% base::nrow()

run_tib %>% dplyr::bind_rows()

# align-summay.bowtie.50000000.tsv.gz
#
# red_sum rem_non rem_und rem_all sel_ref sel_SNP_Max sel_SNP_Sum
# <int>   <int>   <int>   <int>   <int>       <int>       <int>
# 59564    4368     657    4999    9250        9328         237
#

#
# Run all data::
#
# Screen by N/U Probes
#
# Check G-1/2 vs. probe scores
#
# Select top_g-probes
# Select sel_SNP_sum_tib
# Select sel_REF_sum_tib
#



#
# align-summay.bowtie.50000000.tsv.gz
# red_sum rem_non rem_und rem_all sel_ref sel_SNP_1 sel_SNP_tib_2
# <int>   <int>   <int>   <int>   <int>     <int>         <int>
# 59564    4368     657    4999     9250     9328           237
#
# align-summay.bowtie.100000000.tsv.gz
# red_sum rem_non rem_und rem_all sel_ref sel_SNP_1 sel_SNP_tib_2
# <int>   <int>   <int>   <int>   <int>     <int>         <int>
# 59564    3751     814    4488     8708     8832           332
#
# align-summay.bowtie.1000000000.tsv.gz
# red_sum rem_non rem_und rem_all sel_ref sel_SNP_1 sel_SNP_tib_2
# <int>   <int>   <int>   <int>   <int>     <int>         <int>
# 59564    3320    1373    4612    581       618            96
# 
# align-summay.bowtie.10000000000.tsv.gz
# red_sum rem_non rem_und rem_all sel_ref sel_SNP_1 sel_SNP_tib_2
# <int>   <int>   <int>   <int>   <int>     <int>         <int>
# 59564    3200    1230    4382     0         0            19
#

# red_ann_tib %>% dplyr::mutate(SNP_Ratio=SNP_Sum/REF_Sum) %>% dplyr::select(SNP_Ratio, SNP_Sum,REF_Sum)

# red_ann_tib %>% dplyr::filter(SNP_Min <= pmax(0,SNP_Min_avg-(opt$SNP_val*SNP_Min_sd) ) ) %>% dplyr::select(Prb_Key,Tot_Cnt,SNP_Min,SNP_Min_med,SNP_Min_avg,SNP_Min_sd )












#
# Below is general filtering methods::
#

# Filtering::
#  Aln_NON_Perc[A|B] > Aln_NON_avg-Aln_NON_sd
#  Aln_UND_Perc[A|B] > Aln_UND_avg-Aln_UND_sd
#
opt$NON_val <- 0.1
opt$NON_val <- 1
red_ann_tib %>% dplyr::filter(NON_Max >= NON_Max_avg+(opt$NON_val*NON_Max_sd) ) %>% dplyr::select(NON_Max,NON_Max_med,NON_Max_avg,NON_Max_sd )
red_ann_tib %>% dplyr::filter(NON_Min <= pmax(0,NON_Min_avg-(opt$NON_val*NON_Min_sd) ) ) %>% dplyr::select(NON_Min,NON_Min_med,NON_Min_avg,NON_Min_sd )

opt$UND_val <- 0.1
red_ann_tib %>% dplyr::filter(UND_Max >= UND_Max_avg+(opt$UND_val*UND_Max_sd) ) %>% dplyr::select(UND_Max,UND_Max_med,UND_Max_avg,UND_Max_sd )
red_ann_tib %>% dplyr::filter(UND_Min <= pmax(0,UND_Min_avg-(opt$UND_val*UND_Min_sd) ) ) %>% dplyr::select(UND_Min,UND_Min_med,UND_Min_avg,UND_Min_sd )

# Selection::
#  Aln_REF_Perc[A|B] > Aln_REF_avg+Aln_REF_sd
#  Aln_SNP_Perc[A|B] > Aln_SNP_avg+Aln_SNP_sd
#
opt$REF_val <- 0.1
red_ann_tib %>% dplyr::filter(REF_Max >= REF_Max_avg+(opt$ref_val*REF_Max_sd) ) %>% dplyr::select(REF_Max,REF_Max_med,REF_Max_avg,REF_Max_sd )
red_ann_tib %>% dplyr::filter(REF_Min <= pmax(0,REF_Min_avg-(opt$ref_val*REF_Min_sd) ) ) %>% dplyr::select(REF_Min,REF_Min_med,REF_Min_avg,REF_Min_sd )

opt$SNP_val <- 0.1
red_ann_tib %>% dplyr::filter(SNP_Max >= SNP_Max_avg+(opt$snp_val*SNP_Max_sd) ) %>% dplyr::select(SNP_Max,SNP_Max_med,SNP_Max_avg,SNP_Max_sd )
red_ann_tib %>% dplyr::filter(SNP_Min <= pmax(0,SNP_Min_avg-(opt$SNP_val*SNP_Min_sd) ) ) %>% dplyr::select(SNP_Min,SNP_Min_med,SNP_Min_avg,SNP_Min_sd )












# red_ann_tib <- red_sum_tib %>%
#   dplyr::mutate(
#     Aln_Tot_Cnt=red_tot_vec,
#     
#     # Ref Calls::
#     #
#     Aln_Ref_SumA=IA_P,
#     Aln_Ref_SumB=IB_P,
#     
#     Aln_Ref_Max=pmax(Aln_Ref_SumA,Aln_Ref_SumB),
#     Aln_Ref_Min=pmin(Aln_Ref_SumA,Aln_Ref_SumB),
#     
#     Aln_Ref_PercA=round(100*Aln_Ref_SumA/(Aln_Tot_Cnt/2), 1),
#     Aln_Ref_PercB=round(100*Aln_Ref_SumB/(Aln_Tot_Cnt/2), 1),
#     Aln_Ref_PercM=round(100*Aln_Ref_Max/(Aln_Tot_Cnt/2), 1),
#     
#     Aln_Ref_Max_med=median(Aln_Ref_Max), Aln_Ref_Max_avg=mean(Aln_Ref_Max), Aln_Ref_Max_sd=sd(Aln_Ref_Max),
#     Aln_Ref_Min_med=median(Aln_Ref_Min), Aln_Ref_Min_avg=mean(Aln_Ref_Min), Aln_Ref_Min_sd=sd(Aln_Ref_Min),
#     
#     # ALt Calls::
#     #
#     Aln_SNP_SumA=IA_A+IA_C+IA_G+IA_T,
#     Aln_SNP_SumB=IB_A+IB_C+IB_G+IB_T,
#     
#     Aln_SNP_Max=pmax(Aln_SNP_SumA,Aln_SNP_SumB),
#     Aln_SNP_Min=pmin(Aln_SNP_SumA,Aln_SNP_SumB),
#     
#     Aln_SNP_PercA=round(100*Aln_SNP_SumA/(Aln_Tot_Cnt/2), 1),
#     Aln_SNP_PercB=round(100*Aln_SNP_SumB/(Aln_Tot_Cnt/2), 1),
#     Aln_SNP_PercM=round(100*Aln_SNP_Max/(Aln_Tot_Cnt/2), 1),
#     
#     # Und Calls::
#     #
#     Aln_UND_SumA=IA_U,
#     Aln_UND_SumB=IB_U,
#     
#     Aln_UND_Max=pmax(Aln_UND_SumA,Aln_UND_SumB),
#     Aln_UND_Min=pmin(Aln_UND_SumA,Aln_UND_SumB),
#     
#     Aln_UND_PercA=round(100*Aln_UND_SumA/(Aln_Tot_Cnt/2), 1),
#     Aln_UND_PercB=round(100*Aln_UND_SumB/(Aln_Tot_Cnt/2), 1),
#     Aln_UND_PercM=round(100*Aln_UND_Max/(Aln_Tot_Cnt/2), 1),
#     
#     # Non Calls::
#     Aln_NON_SumA=IA_N,
#     Aln_NON_SumB=IB_N,
#     
#     Aln_NON_Max=pmax(Aln_NON_SumA,Aln_NON_SumB),
#     Aln_NON_Min=pmin(Aln_NON_SumA,Aln_NON_SumB),
#     
#     Aln_NON_PercA=round(100*Aln_NON_SumA/(Aln_Tot_Cnt/2), 1),
#     Aln_NON_PercB=round(100*Aln_NON_SumB/(Aln_Tot_Cnt/2), 1),
#     Aln_NON_PercM=round(100*Aln_NON_Max/(Aln_Tot_Cnt/2), 1),
#     
#     
#     # Over all stats for comparison later...
#     #
#     Aln_Ref_med=median(Aln_Ref_PercM), Aln_Ref_avg=mean(Aln_Ref_PercM), Aln_Ref_sd=sd(Aln_Ref_PercM),
#     Aln_Per_med=median(Aln_Ref_Max),  Aln_Per_avg=mean(Aln_Ref_Max),  Aln_Per_sd=sd(Aln_Ref_Max),
#     
#     IA_N_med=median(IA_N), IA_N_avg=mean(IA_N), IA_N_sd=sd(IA_N),
#     IB_N_med=median(IB_N), IB_N_avg=mean(IB_N), IB_N_sd=sd(IB_N),
#     IA_U_med=median(IA_U), IA_U_avg=mean(IA_U), IA_U_sd=sd(IA_U),
#     IB_U_med=median(IB_U), IB_U_avg=mean(IB_U), IB_U_sd=sd(IB_U)
#   )




# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Join/Filter Probes by Design Score::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Join Probe Alignment Stats with Probe Design Information::
#  NOT BY II:
sel_des_tib <- red_ann_tib %>% dplyr::inner_join(snp_des_tib, by=c("Prb_Key"="Ilmn_Id_I") ) 
sel_des_tib %>% base::nrow() %>% print()
run_tib$sel_des_prbFilt <- sel_des_tib %>% base::nrow()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              Filter Probes with Low Global Genomes Alignment::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Filtering::
#  - Remove Probes with > average NA % + sd [Poor Global Virus Alignment]
sel_des_tib <- sel_des_tib %>% dplyr::filter(IA_N<IA_N_med+IA_N_sd & IB_N<IB_N_med+IB_N_sd)
sel_des_tib %>% base::nrow() %>% print()
run_tib$sel_des_lowAln <- sel_des_tib %>% base::nrow()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Filter Probes with Underlying SNPs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Filtering::
#  - Remove Probes with > average U  % - sd [Too Many Underlying Virus SNPs]
sel_des_tib <- sel_des_tib %>% dplyr::filter(IA_U <IA_U_med+IA_U_sd & IB_U <IB_U_med+IB_U_sd)
sel_des_tib %>% base::nrow() %>% print()
run_tib$sel_des_uSNPs <- sel_des_tib %>% base::nrow()

#
# Now calculate Selection Stats::
#

#
# Original Selection Parameters::
#
# opt$g_max <- 3
# opt$a_val <- 2.4
# opt$p_val <- 4

# New Parameters after earlier filtering::

opt$g_max <- 3
opt$a_per <- 0.78
opt$a_val <- 2.4
opt$a_val <- 0.25
opt$p_val <- 0.25

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Select Bisulfite Conversion:: low G content
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: Look at score of probes with 1 G::
#
#    sel_des_tib %>% dplyr::group_by(G_Count_Max) %>% dplyr::summarise(Count=n())
sel_des_gmax_tib <- sel_des_tib %>% dplyr::filter(G_Count_Max<=opt$g_max) %>% dplyr::arrange(G_Count_Max)
sel_des_gmax_tib %>% base::nrow() %>% print()
sel_des_gmax_tib %>% dplyr::select(Prb_Key,G_Count_Max) %>% head(n=3) %>% print()
sel_des_gmax_tib %>% dplyr::select(Prb_Key,G_Count_Max) %>% tail(n=3) %>% print()
sel_des_gmax_tib %>% dplyr::group_by(G_Count_Max) %>% dplyr::summarise(Count=n()) %>% print()
run_tib$sel_des_gmax <- sel_des_gmax_tib %>% base::nrow()

#
#
# The Two Below Need to be Switched in Logic; calling opposite things...
#
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Select High Viral SNP Detection
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Sum of IA(ACTG)

sel_des_sper_tib <- sel_des_tib %>% dplyr::filter(Aln_Ref_PercM<opt$a_per) %>% dplyr::arrange(Aln_Ref_PercM)
sel_des_sper_tib %>% base::nrow() %>% print()
sel_des_sper_tib %>% dplyr::select(Prb_Key,Aln_Ref_PercM,Aln_Ref_med,Aln_Ref_avg,Aln_Ref_sd) %>% head(n=3) %>% print()
sel_des_sper_tib %>% dplyr::select(Prb_Key,Aln_Ref_PercM,Aln_Ref_med,Aln_Ref_avg,Aln_Ref_sd) %>% tail(n=3) %>% print()
run_tib$sel_des_sper <- sel_des_sper_tib %>% base::nrow()

sel_des_smax_tib <- sel_des_tib %>% dplyr::filter(Aln_Ref_PercM <Aln_Ref_avg-(opt$a_val*Aln_Ref_sd) ) %>% dplyr::arrange(Aln_Ref_PercM)
sel_des_smax_tib %>% base::nrow() %>% print()
sel_des_smax_tib %>% dplyr::select(Prb_Key,Aln_Ref_PercM,Aln_Ref_med,Aln_Ref_avg,Aln_Ref_sd) %>% head(n=3) %>% print()
sel_des_smax_tib %>% dplyr::select(Prb_Key,Aln_Ref_PercM,Aln_Ref_med,Aln_Ref_avg,Aln_Ref_sd) %>% tail(n=3) %>% print()
run_tib$sel_des_smax <- sel_des_smax_tib %>% base::nrow()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Select Highly Conserved Viral Probes
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sel_des_pmax_tib <- sel_des_tib %>% dplyr::filter(IA_P==0 | IB_P==0) %>%
  dplyr::filter(Aln_Ref_Max <Aln_Per_avg-(opt$p_val*Aln_Per_sd) ) %>% dplyr::arrange(Aln_Ref_Max)
sel_des_pmax_tib %>% base::nrow() %>% print()
sel_des_pmax_tib %>% dplyr::select(Prb_Key,Aln_Ref_Max,Aln_Per_med,Aln_Per_avg,Aln_Per_sd) %>% head(n=3) %>% print()
sel_des_pmax_tib %>% dplyr::select(Prb_Key,Aln_Ref_Max,Aln_Per_med,Aln_Per_avg,Aln_Per_sd) %>% tail(n=3) %>% print()
run_tib$sel_des_pmax <- sel_des_pmax_tib %>% base::nrow()




# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Join Selected Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Join alll selections::
sel_bind_tib <- dplyr::bind_rows(sel_des_gmax_tib,sel_des_smax_tib,sel_des_pmax_tib)
sel_bind_tib %>% dplyr::distinct(Prb_Key) %>% base::nrow() %>% print()
run_tib$sel_bind <- sel_bind_tib %>% base::nrow()

# 2336 = align-summay.bowtie.50000000.tsv.gz
# 2172 = align-summay.bowtie.100000000.tsv.gz
# 1586 = align-summay.bowtie.1000000000.tsv.gz

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               Plot Probe Density by Genome Coordinates::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ggplot2::ggplot(data=sel_des_gmax_tib, aes(x=Coordinate)) + ggplot2::geom_density()
ggplot2::ggplot(data=sel_des_smax_tib, aes(x=Coordinate)) + ggplot2::geom_density()
ggplot2::ggplot(data=sel_des_pmax_tib, aes(x=Coordinate)) + ggplot2::geom_density()



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                 OLD CODE::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                          Search for Alignments::
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
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Load Alignments::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
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
}

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
