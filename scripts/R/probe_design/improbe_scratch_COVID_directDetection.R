

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
par$prgmDir <- 'probe_design'
par$prgmTag <- 'improbe_main'

# Illumina based directories::
par$macDir  <- '/Users/bbarnes/Documents/CustomerFacing'
par$lixDir  <- '/illumina/scratch/darkmatter/Projects/COVIC'

# Directory Parameters::
opt$outDir <- NULL
opt$topDir <- '/Users/bbarnes/Documents/Projects/methylation/'
opt$datDir <- file.path(opt$topDir, 'scratch/tile_main/EPIC/SARS-CoV-2/MN908947/COVIC')

opt$snp_covid_sam <- file.path(opt$datDir, 'tile_main_EPIC_SARS-CoV-2_MN908947_COVIC.snp.nCoV_Wuhan_Sequence_MN908947.3.tsv.gz')
opt$snp_ncbi_sam  <- file.path(opt$datDir, 'tile_main_EPIC_SARS-CoV-2_MN908947_COVIC.snp.ncbi.tsv.gz')

opt$snp_sam <- file.path('/Users/bbarnes/Documents/Projects/methylation/scratch/tile_main_EPIC_SARS-CoV-2_MN908947_COVIC.snp-LC528232.1_11042bowtie.sam.gz')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Load Alignments::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sam_col_vec <- c('QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL',
                 'AS', 'XN', 'XM', 'XO', 'XG', 'NM', 'MD', 'YT')

# snp_covid_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file = opt$snp_covid_sam, col_names=sam_col_vec, comment='@') ))
# snp_ncbi_tib  <- suppressMessages(suppressWarnings( readr::read_tsv(file = opt$snp_ncbi_sam, col_names=sam_col_vec, comment='@') ))

snp_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file = opt$snp_sam, col_names=sam_col_vec, comment='@') ))

snp_covid_tib %>% dplyr::group_by(CIGAR) %>% dplyr::summarise(Count=n())
snp_ncbi_tib %>% dplyr::group_by(CIGAR) %>% dplyr::summarise(Count=n())

# Its lack of G::
# gzip -dc /Users/bbarnes/Documents/Projects/manifests/methylation/HumanMethylation450_15017482_v.1.2.rs-only.csv.gz | grep G
#
# snp_ncbi_tib %>% dplyr::mutate(G_Count=50 - stringr::str_length(stringr::str_remove_all(SEQ, 'G')) ) %>% dplyr::arrange(G_Count) %>% dplyr::select(QNAME,SEQ,G_Count)




# End of file
