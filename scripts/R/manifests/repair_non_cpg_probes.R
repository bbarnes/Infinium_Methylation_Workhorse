
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# suppressPackageStartupMessages(base::require() )
# Load sesame:: This causes issues with "ExperimentHub Caching causes a warning"
suppressWarnings(suppressPackageStartupMessages( base::require("sesame") ))
suppressWarnings(suppressPackageStartupMessages( base::require("dbplyr") ))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("plyr")) )
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("readr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

# tidyverse_update(recursive = FALSE, repos = getOption("repos"))
# install.packages("lubridate")

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
par$prgmTag <- 'repair_non_cpg_probes'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

par$retData     <- FALSE

# Executables::
opt$Rscript <- NULL

# Directories::
opt$outDir     <- NULL

# Platform/Method Options::
opt$platform  <- NULL
opt$manifest  <- NULL

# Run Options::
opt$fresh        <- FALSE

opt$time_org_txt <- NULL

# verbose Options::
opt$verbose <- 3

# Set Default Options::
#
def <- opt

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
  
  opt$outDir <- file.path(par$topDir, 'scratch', par$prgmTag)
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda2-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda3-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/conda_4.6.8/bin/Rscript'
  
  #
  # End of local parameter definitions::
  #
  
  opt$fresh  <- TRUE
  opt$outDir <- file.path(par$topDir, 'scratch')
  
} else {
  par$runMode    <- 'CommandLine'
  par$exePath <- base::substring(args.dat[grep("--file=", args.dat)], 8)
  
  cat(glue::glue("[{par$prgmTag}]: Local Run par$runMode={par$runMode}.{RET}"))
  
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
    
    # Directories::
    make_option(c("--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    
    # Optional Files::
    make_option(c("--manifestPath"), type="character", default=opt$manifestPath,
                help="Path to manfifest (CSV) otherwise use dat [default= %default]", metavar="character"),
    
    # Platform/Method Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Forced platform [EPIC, 450k, 27k, NZT] otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--manifest"), type="character", default=opt$manifest, 
                help="Forced manifest [B1, B2, B4] otherwise auto-detect [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--fresh"), action="store_true", default=opt$fresh,
                help="Boolean variable to build fresh version of database files [default= %default]", metavar="boolean"),
    
    make_option(c("--time_org_txt"), type="character", default=opt$time_org_txt, 
                help="Unused variable time_org_txt [default= %default]", metavar="character"),
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
opt_reqs <- c('outDir','Rscript','verbose')

par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) 
  stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))

for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', 
                         full.names=TRUE, recursive=TRUE)) base::source(sfile)
if (opt$verbose>=0)
  cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form ",
                 "General Source={par$gen_src_dir}!{RET}{RET}") )

opt <- program_init(name=par$prgmTag,
                    opts=opt, opt_reqs=opt_reqs, 
                    pars=par, par_reqs=par_reqs,
                    libs=TRUE,rcpp=FALSE,
                    verbose=opt$verbose,vt=3,tc=0,tt=NULL)

opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker  <- timeTracker$new(verbose=opt$verbose)

image_key <- "bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse_Centos"
image_ver <- "v.1.0"
image_ssh <- "run_improbe.sh"
image_str <- glue::glue("{image_key}.{image_ver}")

#
# Sesame Manifest Loading::
#
ses_19_tib <- sesameData::sesameDataGet("HM450.hg19.manifest") %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var="Seq_ID") %>%
  tibble::as_tibble() %>%
  dplyr::mutate(Probe_Type=stringr::str_sub(Seq_ID, 1,2),
                Genome_Build="37",
                seqnames=stringr::str_remove(seqnames,'^chr'),
                Strand_FR=dplyr::case_when(
                  strand=='+' ~ 'F', strand=='-' ~ 'R', TRUE ~ NA_character_)
  ) %>%
  dplyr::rename(AlleleA_ProbeSeq=ProbeSeq_A,
                AlleleB_ProbeSeq=ProbeSeq_B,
                Chromosome=seqnames,
                Coordinate=start) %>%
  dplyr::select(Seq_ID,
                AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                Genome_Build,Chromosome,Coordinate,Strand_FR,
                Probe_Type) %>% 
  dplyr::mutate(AlleleB_ProbeSeq=dplyr::case_when(
    stringr::str_length(AlleleB_ProbeSeq)==0 ~ NA_character_, 
    TRUE ~ AlleleB_ProbeSeq)) %>%
  dplyr::mutate(Genome_Build="GRCh37")

ses_hsa_snp_tib <- ses_19_tib %>% dplyr::filter(Probe_Type=='rs')
ses_hsa_cph_tib <- ses_19_tib %>% dplyr::filter(Probe_Type=='ch')

ses_mus_man_csv <- 
  file.path(par$datDir, 'manifest/core/LEGX-C20.manifest.sesame-base.cpg-sorted.csv.gz')
ses_mus_man_tib <- 
  suppressMessages(suppressWarnings( readr::read_csv(ses_mus_man_csv) )) %>%
  dplyr::mutate(Genome_Build="GRCm38")

ses_mus_snp_tib <- ses_mus_man_tib %>% dplyr::filter(Probe_Type=='rs')
ses_mus_cph_tib <- ses_mus_man_tib %>% dplyr::filter(Probe_Type=='ch')

#
#
# TBD:: Add minimalization fields and inferred fields from Sesame...
#
#

sidx=2 
plen=50

idx1 <- sidx
len1 <- plen - 1
idx2 <- sidx + 1
len2 <- plen

prb_seq_tib <- dplyr::bind_rows(
  ses_hsa_snp_tib %>%
    dplyr::select(Seq_ID,AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                  Genome_Build,Probe_Type) %>%
    dplyr::distinct(AlleleA_ProbeSeq,AlleleB_ProbeSeq, .keep_all=TRUE),
  
  ses_hsa_cph_tib %>%
    dplyr::select(Seq_ID,AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                  Genome_Build,Probe_Type) %>%
    dplyr::distinct(AlleleA_ProbeSeq,AlleleB_ProbeSeq, .keep_all=TRUE),
  
  ses_mus_snp_tib %>% 
    dplyr::rename(AlleleA_ProbeSeq=AlleleA_Probe_Sequence,
                  AlleleB_ProbeSeq=AlleleB_Probe_Sequence) %>%
    dplyr::select(Seq_ID,AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                  Genome_Build,Probe_Type) %>%
    dplyr::distinct(AlleleA_ProbeSeq,AlleleB_ProbeSeq, .keep_all=TRUE),
  
  ses_mus_cph_tib %>% 
    dplyr::rename(AlleleA_ProbeSeq=AlleleA_Probe_Sequence,
                  AlleleB_ProbeSeq=AlleleB_Probe_Sequence) %>%
    dplyr::select(Seq_ID,AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                  Genome_Build,Probe_Type) %>%
    dplyr::distinct(AlleleA_ProbeSeq,AlleleB_ProbeSeq, .keep_all=TRUE)
) %>% 
  dplyr::mutate(
    AlleleA_ProbeSeq=dplyr::case_when(
      stringr::str_length(AlleleA_ProbeSeq)==0 ~ NA_character_,
      TRUE ~ AlleleA_ProbeSeq
    ),
    AlleleB_ProbeSeq=dplyr::case_when(
      stringr::str_length(AlleleB_ProbeSeq)==0 ~ NA_character_,
      TRUE ~ AlleleB_ProbeSeq
    ),
    Infinium_Design=dplyr::case_when(
      !is.na(AlleleA_ProbeSeq) & !is.na(AlleleB_ProbeSeq) ~ 1,
      !is.na(AlleleA_ProbeSeq) &  is.na(AlleleB_ProbeSeq) ~ 2,
      TRUE ~ NA_real_),
    Seq_48U=dplyr::case_when(
      Infinium_Design==1 ~ stringr::str_sub(AlleleA_ProbeSeq, idx1,len1),
      Infinium_Design==2 ~ stringr::str_sub(AlleleA_ProbeSeq, idx2,len2),
      TRUE ~ NA_character_
    ) %>% 
      stringr::str_to_upper() %>% 
      stringr::str_replace_all('R', 'A') %>% stringr::str_replace_all('Y', 'T'),
    Mat_48U=stringr::str_sub(AlleleA_ProbeSeq,2)
  )

# Summary Count::
prb_seq_tib %>% 
  dplyr::add_count(Seq_ID, name="Seq_ID_Count") %>% 
  # dplyr::filter(Seq_ID_Count!=1) %>% 
  dplyr::group_by(Genome_Build,Probe_Type,Infinium_Design,Seq_ID_Count) %>% 
  dplyr::summarise(Count=n(), .groups="drop") %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Load Data:: MUS-CpHs
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

probe_type <- "ch"
genome_ver <- "GRCm38"
build_dir  <- file.path(opt$outDir, probe_type, genome_ver)
if (!dir.exists(build_dir)) dir.create(build_dir, recursive=TRUE)

# Mouse Data::
mus_cph_des_tsv <- '/Users/bretbarnes/Documents/data/CustomContent/LifeEpigentics/data/20190228_input_files/cph_input.txt.gz'
mus_cph_pre_csv <- '/Users/bretbarnes/Documents/data/improbe/cph-snp-designs/LEGX_SpikeIn_Reorder-CpH-Only.designs.csv.gz'
mus_cph_bed_tsv <- '/Users/bretbarnes/Documents/data/CustomContent/LifeEpigentics/Redesign/data/CpH/selected_CpH_probes.bed'

mus_cph_bed_tib <- readr::read_tsv(
  mus_cph_bed_tsv, 
  col_names=c("Chromosome","Coordinate","Pos_End","Direction","Full_ID")) %>% 
  dplyr::mutate(Seq_ID=stringr::str_remove(Full_ID, '_.*$')) %>% 
  dplyr::select(Seq_ID,Chromosome,Coordinate)

mus_cph_pre_tib <- 
  suppressMessages(suppressWarnings( readr::read_csv(mus_cph_pre_csv) )) %>%
  dplyr::mutate(Seq_ID=stringr::str_replace_all(Seq_ID, regex("\\W+"), ""),
                Source=1) %>% 
  dplyr::select(Seq_ID,IUPAC_Forward_Sequence,Source) %>% 
  dplyr::rename(IUPAC_Sequence=IUPAC_Forward_Sequence) %>%
  dplyr::inner_join(mus_cph_bed_tib, by="Seq_ID") %>%
  dplyr::distinct()

mus_cph_des_tib <- 
  suppressMessages(suppressWarnings( readr::read_tsv(mus_cph_des_tsv) )) %>%
  dplyr::mutate(Seq_ID=stringr::str_replace_all(Seq_ID, regex("\\W+"), ""),
                Source=2) %>% 
  dplyr::select(Seq_ID,Sequence,Source,Chromosome,Coordinate) %>% 
  dplyr::rename(IUPAC_Sequence=Sequence) %>%
  dplyr::distinct()

mus_cph_tib <- dplyr::bind_rows(mus_cph_des_tib, mus_cph_pre_tib) %>% 
  dplyr::arrange(Seq_ID,Source) %>%
  dplyr::add_count(Seq_ID, name="ID_Rep_Count") %>% 
  dplyr::distinct(Seq_ID, .keep_all=TRUE) %>%
  addDesignSeqCG(seq="IUPAC_Sequence",add="Sequence", din="DiNuc",
                 verbose=opt$verbose, vt=1,tc=1,tt=pTracker) %>%
  dplyr::mutate(Genome_Build=genome_ver,
                CpG_Island="FALSE",
                Chromosome=stringr::str_remove(Chromosome,'^chr'),
                Probe_Type=probe_type) %>%
  dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island,
                DiNuc,Probe_Type,IUPAC_Sequence)

# mus_cph_des_dat <- improbe_design_all(
#   tib=mus_cph_tib, ptype=probe_type, outDir=build_dir, 
#   gen=genome_ver, image=image_str, shell=image_ssh, parse_din=TRUE,
#   verbose=opt$verbose, vt=1,tc=1,tt=pTracker)

#
# Re-joing::
#
# prb_seq_tib %>% dplyr::inner_join(
#   dplyr::mutate(mus_cph_des_dat$iup, Mat_48U=stringr::str_sub(PRB2_D_MAT,2)),
#   by=c("Seq_ID","Mat_48U","Probe_Type")) %>% 
#   dplyr::group_by(Probe_Type,Strand_SR,Strand_CO,Infinium_Design) %>% 
#   dplyr::summarise(Count=n(), .groups="drop")


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Load Data:: HSA-CpHs
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

probe_type <- "ch"
genome_ver <- "GRCh37"
build_dir  <- file.path(opt$outDir, probe_type, genome_ver)
if (!dir.exists(build_dir)) dir.create(build_dir, recursive=TRUE)

hsa_cph_csv <- "/Users/bretbarnes/Documents/data/manifests/raw/manifests/methylation/ch-repair/HumanMethylation450_15017482_v.1.2.ch-noAnnotation.csv.gz"
hsa_cph_tib <- suppressMessages(suppressWarnings( readr::read_csv(hsa_cph_csv) )) %>%
  addDesignSeqCG(seq="Forward_Sequence", add="Sequence", din="DiNuc", verbose=4) %>% 
  dplyr::select(IlmnID,DiNuc,Sequence) %>%
  dplyr::rename(Seq_ID=IlmnID) %>%
  dplyr::inner_join(
    dplyr::select(ses_hsa_cph_tib, Seq_ID,Chromosome,Coordinate),
    by="Seq_ID"
  ) %>%
  dplyr::mutate(Genome_Build=genome_ver,
                CpG_Island="FALSE",
                Probe_Type=probe_type) %>%
  replaceDesignSeqCG(seq="Sequence", add="IUPAC_Sequence", nuc="DiNuc",
                     verbose=opt$verbose, vt=1,tc=1,tt=pTracker) %>%
  dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island,
                DiNuc,Probe_Type,IUPAC_Sequence)

# hsa_cph_des_dat <- improbe_design_all(
#   tib=hsa_cph_tib, ptype=probe_type, outDir=build_dir, 
#   gen=genome_ver, image=image_str, shell=image_ssh, parse_din=TRUE,
#   verbose=opt$verbose, vt=1,tc=1,tt=pTracker)

#
# Re-joing::
#
# prb_seq_tib %>% dplyr::inner_join(
#   dplyr::mutate(hsa_cph_des_dat$iup, Mat_48U=stringr::str_sub(PRB2_D_MAT,2)),
#   by=c("Seq_ID","Mat_48U","Probe_Type")) %>% 
#   dplyr::group_by(Probe_Type,Strand_SR,Strand_CO,Infinium_Design) %>% 
#   dplyr::summarise(Count=n(), .groups="drop")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Load Data:: MUS-SNPs
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

probe_type <- "rs"
genome_ver <- "GRCm38"
build_dir  <- file.path(opt$outDir, probe_type, genome_ver)
if (!dir.exists(build_dir)) dir.create(build_dir, recursive=TRUE)

# Mouse Data::
mus_snp_bed_tsv <- '/Users/bretbarnes/Documents/data/CustomContent/LifeEpigentics/Redesign/data/SNP/selected_SNP_probes.bed'
mus_snp_pre_csv <- '/Users/bretbarnes/Documents/data/improbe/cph-snp-designs/LEGX_SpikeIn_Reorder-SNP-Only.designs.csv.gz'

mus_snp_bed_tib <- suppressMessages(suppressWarnings( readr::read_tsv(
  mus_snp_bed_tsv,col_names=c("Chromosome","Coordinate","Pos_End","Direction","Full_ID") ) ) ) %>% 
  dplyr::mutate(Seq_ID=stringr::str_remove(Full_ID, '_.*$')) %>% 
  dplyr::mutate(Genome_Build=genome_ver, 
                Probe_Type=probe_type,
                CpG_Island="FALSE") %>% 
  dplyr::select(Seq_ID,Genome_Build,Chromosome,Coordinate,CpG_Island,
                dplyr::everything())

mus_snp_pre_tib <- 
  suppressMessages(suppressWarnings( readr::read_csv(mus_snp_pre_csv) )) %>%
  dplyr::select(Seq_ID,IUPAC_Forward_Sequence) %>%
  dplyr::rename(IUPAC_Sequence=IUPAC_Forward_Sequence) %>%
  addDesignSeqCG(seq="IUPAC_Sequence",add="Sequence", din="DiNuc",
                 verbose=opt$verbose, vt=1,tc=1,tt=pTracker)

mus_snp_tib <- mus_snp_pre_tib %>% 
  dplyr::inner_join(mus_snp_bed_tib, by="Seq_ID") %>% 
  dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island,
                DiNuc,Probe_Type,IUPAC_Sequence)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Load Data:: HSA-SNPs
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

probe_type <- "rs"
genome_ver <- "GRCh37"
build_dir  <- file.path(opt$outDir, probe_type, genome_ver)
if (!dir.exists(build_dir)) dir.create(build_dir, recursive=TRUE)

hsa_snp_tsv <- '/Users/bretbarnes/Documents/data/manifests/raw/manifests/methylation/rs-repair/rs_swap_improbe_input.tsv.gz'
hsa_snp_tib <- suppressMessages(suppressWarnings( readr::read_tsv(hsa_snp_tsv) )) %>%
  dplyr::select(Seq_ID,Sequence) %>%
  tidyr::separate(Seq_ID, into=c("Seq_ID","DiNuc"), sep='_') %>%
  dplyr::inner_join(
    dplyr::select(ses_hsa_snp_tib, Seq_ID,Chromosome,Coordinate),
    by="Seq_ID"
  ) %>%
  dplyr::mutate(Genome_Build=genome_ver,
                CpG_Island="FALSE",
                Probe_Type=probe_type) %>%
  replaceDesignSeqCG(seq="Sequence", add="IUPAC_Sequence", nuc="DiNuc",
                     verbose=opt$verbose, vt=1,tc=1,tt=pTracker) %>%
  dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island,
                DiNuc,Probe_Type,IUPAC_Sequence)

if (FALSE) {
  # Match real probes to IUPAC designs:: SNP
  # mat1_snp_tib <- ses_hsa_snp_tib %>% 
  mat1_snp_tib <- prb_seq_tib %>% 
    dplyr::inner_join(hsa_snp_des_dat$iup, by=c("Seq_ID","Probe_Type","AlleleA_ProbeSeq"="PRB1_M_MAT")) %>%
    dplyr::distinct(Seq_ID, .keep_all=TRUE) %>%
    dplyr::mutate(Infinium_Design=1)
  
  # mat2_snp_tib <- ses_hsa_snp_tib %>% 
  mat2_snp_tib <- prb_seq_tib %>% 
    dplyr::inner_join(hsa_snp_des_dat$iup, by=c("Seq_ID","Probe_Type","AlleleA_ProbeSeq"="PRB2_D_MAT")) %>%
    dplyr::distinct(Seq_ID, .keep_all=TRUE) %>%
    dplyr::mutate(Infinium_Design=2)
  
  mat_snp_tib <- dplyr::bind_rows(mat1_snp_tib,mat2_snp_tib)
  
  mat_snp_tib %>%
    dplyr::group_by(Probe_Type,Strand_SR,Strand_CO,Infinium_Design) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Join all Unique All Probes Designs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

all_des_tib <- dplyr::bind_rows(
  dplyr::distinct(hsa_cph_tib),
  dplyr::distinct(hsa_snp_tib),
  dplyr::distinct(mus_cph_tib),
  dplyr::distinct(mus_snp_tib) )
  
all_des_tib %>% base::nrow() %>% print()
all_des_tib %>% dplyr::distinct(Seq_ID) %>% base::nrow() %>% print()
all_des_tib %>% dplyr::distinct(Sequence) %>% base::nrow() %>% print()
all_des_tib %>% dplyr::distinct(Genome_Build) %>% base::nrow() %>% print()
all_des_tib %>% dplyr::distinct(Genome_Build,Chromosome,Coordinate) %>% base::nrow() %>% print()
all_des_tib %>% dplyr::distinct(IUPAC_Sequence) %>% base::nrow() %>% print()

all_des_tib %>% dplyr::add_count(Seq_ID, name="SeqID_Cnt") %>% dplyr::filter(SeqID_Cnt != 1)
all_des_tib %>% dplyr::add_count(Genome_Build,Chromosome,Coordinate, name="SeqID_Cnt") %>% dplyr::filter(SeqID_Cnt != 1)

all_des_tib %>% dplyr::group_by(Probe_Type,Genome_Build) %>% dplyr::summarise(Count=n(), .grops="drop")

#
# Split Designs by type and genome build::
#   - build Seq_48U
#   - remove/classify junk designs
#

ptype_fwd_dat <- all_des_tib %>% split(.$Probe_Type)

ptype_des_dat <- NULL
for (ptype in names(ptype_fwd_dat)) {
  cur_out_dir <- file.path(opt$outDir, ptype)
  cur_des_tib <- ptype_fwd_dat[[ptype]]
  cur_des_dat <- improbe_design_all(
    tib=cur_des_tib, ptype=ptype,
    outDir=cur_out_dir, gen="all",
    image=image_str, shell=image_ssh, 
    parse_din=FALSE,
    # parse_din=TRUE,
    verbose=opt$verbose, vt=1,tc=1,tt=pTracker)
  
  ptype_des_dat[[ptype]] = cur_des_dat
}

prb_seq_tib %>% dplyr::inner_join(
  ptype_des_dat$rs$imp, 
  by=c("Mat_48U"="Seq_48U_1") )

prb_seq_tib %>% dplyr::inner_join(
  ptype_des_dat$rs$imp, 
  by=c("Mat_48U"="Seq_48U_2") )

prb_seq_tib %>% dplyr::inner_join(
  ptype_des_dat$rs$imp, 
  by=c("Seq_48U"="Seq_48U_1") )

prb_seq_tib %>% dplyr::inner_join(
  ptype_des_dat$rs$imp, 
  by=c("Seq_48U"="Seq_48U_2") )




prb_seq_tib %>% dplyr::inner_join(
  ptype_des_dat$ch$imp, 
  by=c("Mat_48U"="Seq_48U_1") )

prb_seq_tib %>% dplyr::inner_join(
  ptype_des_dat$ch$imp, 
  by=c("Mat_48U"="Seq_48U_2") )

prb_seq_tib %>% dplyr::inner_join(
  ptype_des_dat$ch$imp, 
  by=c("Seq_48U"="Seq_48U_1") )

prb_seq_tib %>% dplyr::inner_join(
  ptype_des_dat$ch$imp, 
  by=c("Seq_48U"="Seq_48U_2") )



# mus_cph_des_dat <- improbe_design_all(
#   tib=all_des_tib, ptype=probe_type, outDir=build_dir,
#   gen=genome_ver, image=image_str, shell=image_ssh, parse_din=TRUE,
#   verbose=opt$verbose, vt=1,tc=1,tt=pTracker)


# all_des_tib %>% 
#   dplyr::inner_join(prb_seq_tib, by=c("Seq_ID","Genome_Build","Probe_Type") ) %>%
#   dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island,
#                 DiNuc,Probe_Type,Infinium_Design,
#                 IUPAC_Sequence,AlleleA_ProbeSeq,AlleleB_ProbeSeq,
#                 Seq_48U,Mat_48U,dplyr::everything())

# prb_seq_tib %>% dplyr::inner_join(
#   dplyr::mutate(hsa_cph_des_dat$iup, Mat_48U=stringr::str_sub(PRB2_D_MAT,2)),
#   by=c("Seq_ID","Mat_48U","Probe_Type")) %>% 
#   dplyr::group_by(Probe_Type,Strand_SR,Strand_CO,Infinium_Design) %>% 
#   dplyr::summarise(Count=n(), .groups="drop")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             improbe design:: All
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  hsa_cph_tib %>% dplyr::inner_join(prb_seq_tib, by=c("Seq_ID","Genome_Build","Probe_Type") ) %>%
    dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island,
                  DiNuc,Probe_Type,Infinium_Design,
                  IUPAC_Sequence,AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                  Seq_48U,Mat_48U,dplyr::everything())
  
  hsa_snp_tib %>% dplyr::inner_join(prb_seq_tib, by=c("Seq_ID","Genome_Build","Probe_Type") ) %>%
    dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island,
                  DiNuc,Probe_Type,Infinium_Design,
                  IUPAC_Sequence,AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                  Seq_48U,Mat_48U,dplyr::everything())
  
  mus_cph_tib %>% dplyr::inner_join(prb_seq_tib, by=c("Seq_ID","Genome_Build","Probe_Type") ) %>%
    dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island,
                  DiNuc,Probe_Type,Infinium_Design,
                  IUPAC_Sequence,AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                  Seq_48U,Mat_48U,dplyr::everything())
  
  mus_snp_tib %>% dplyr::inner_join(prb_seq_tib, by=c("Seq_ID","Genome_Build","Probe_Type") ) %>%
    dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island,
                  DiNuc,Probe_Type,Infinium_Design,
                  IUPAC_Sequence,AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                  Seq_48U,Mat_48U,dplyr::everything())
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
