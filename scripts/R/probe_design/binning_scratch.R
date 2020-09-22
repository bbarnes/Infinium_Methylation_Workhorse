
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
suppressWarnings(suppressPackageStartupMessages(require("GenomicRanges")) )
suppressWarnings(suppressPackageStartupMessages(require("GenomeInfoDb")) )


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
#                              Main Scratch::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ordDir <- '/Users/bbarnes/Documents/Projects/manifests'

ord_cols <- c('Seq_ID', 'Sequence', 'Genome_Build', 'Chromosome', 'Coordinate', 'CpG_Island')
ordA_tsv <- file.path(ordDir, 'methylation/designInput.27k.tsv.gz')
ordB_tsv <- file.path(ordDir, 'methylation/designInput.450k-EPIC-B2.tsv.gz')
ordC_tsv <- file.path(ordDir, 'methylation/designInput.EPIC-B2.tsv.gz')

ordA_tib <- suppressMessages(suppressWarnings( readr::read_tsv(ordA_tsv, col_names=ord_cols) ))
ordB_tib <- suppressMessages(suppressWarnings( readr::read_tsv(ordB_tsv, col_names=ord_cols) ))
ordC_tib <- suppressMessages(suppressWarnings( readr::read_tsv(ordC_tsv, col_names=ord_cols) ))

ord_bind_tib <- dplyr::bind_rows(
  ordC_tib,
  ordB_tib,
  ordA_tib
)
ord_uniq_tib <- ord_bind_tib %>% dplyr::distinct(Seq_ID, .keep_all=TRUE)
ord_uniq_tsv <- file.path('/Users/bbarnes/Documents/Projects/methylation/EWAS/improbe', 'EPIC-reorder.improbe-input.tsv.gz')

readr::write_tsv(ord_uniq_tib, ord_uniq_tsv)

elly_tib <- readr::read_csv('/Users/bbarnes/Documents/Projects/methylation/ElysiumHealth/targets/EH-cg.txt.gz')
legx_tib <- readr::read_csv('/Users/bbarnes/Documents/Projects/methylation/EWAS/data/LEGX/epicplus-important-probes.csv') %>% 
  tidyr::separate(probe_name, into=c('Probe_ID', 'IF', 'FR', 'CO', 'RP', 'GN'), sep='_')

dplyr::right_join(
  elly_tib %>% dplyr::distinct(Probe_ID),
  legx_tib %>% dplyr::distinct(Probe_ID),
  by="Probe_ID"
)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Load Manifests::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

manDir   <- '/Users/bbarnes/Documents/Projects/manifests/methylation/Sesame/hg38'
manA_tsv <- file.path(manDir, 'HM27.hg38.manifest.gencode.v22.tsv.gz')
manB_tsv <- file.path(manDir, 'HM450.hg38.manifest.gencode.v22.tsv.gz')
manC_tsv <- file.path(manDir, 'EPIC.hg38.manifest.gencode.v22.tsv.gz')

manA_tib <- suppressMessages(suppressWarnings( readr::read_tsv(manA_tsv) ))
manB_tib <- suppressMessages(suppressWarnings( readr::read_tsv(manB_tsv) ))
manC_tib <- suppressMessages(suppressWarnings( readr::read_tsv(manC_tsv) ))


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              Squash All Current Design Seqs to Uniq:: hg19
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

des_seq_tsv  <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/improbeOutput.37.cgn-topSeq.tsv.gz'
unq_seq_tsv  <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/improbeOutput.37.cgn-topSeq.uniq.tsv.gz'
des_seq_tibc <- readr::read_tsv(des_seq_tsv)
des_unq_tib <- des_seq_tib %>% dplyr::distinct()
des_unq_tib <- readr::write_tsv(des_unq_tib,unq_seq_tsv)


# End of file
