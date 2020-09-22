
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Manifest Colllection Scratch::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))
suppressWarnings(suppressPackageStartupMessages( base::require("grid") ))

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
par <- NULL
opt <- NULL

# Program Parameters::
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'manifests'
par$prgmTag <- 'manifest_scratch'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# Illumina based directories::
par$macDir <- '/Users/bbarnes/Documents/Projects/methylation/tools'
par$lixDir <- '/illumina/scratch/darkmatter'

par$retData <- FALSE

# Executables::
opt$Rscript <- NULL

# Directories::
opt$outDir     <- NULL

# Required Inputs::
opt$ords <- NULL
opt$mats <- NULL
opt$aqps <- NULL
opt$pqcs <- NULL

opt$idats <- NULL

# Pre-defined files (controls)
opt$ctlCSV <- NULL

# Platform/Method Options::
opt$genomeBuild <- NULL
opt$platform    <- NULL
opt$version     <- NULL

# Run Options::
opt$fresh <- FALSE

opt$percisionSigs <- 1
opt$percisionBeta <- 4
opt$percisionPval <- 6

# Parallel/Cluster Options::
opt$single   <- FALSE
opt$parallel <- FALSE
opt$cluster  <- FALSE

# verbose Options::
opt$verbose <- 3

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
    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    
    # Pre-defined files (controls)
    make_option(c("--ords"), type="character", default=opt$ords, 
                help="Order files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mats"), type="character", default=opt$mats, 
                help="Match files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqps"), type="character", default=opt$aqps, 
                help="AQP files (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--pqcs"), type="character", default=opt$pqcs, 
                help="PQC files (comma seperated) [default= %default]", metavar="character"),
    
    make_option(c("-i", "--idats"), type="character", default=opt$idats, 
                help="idats directory [default= %default]", metavar="character"),
    
    # Required Inputs::
    make_option(c("--ctlCSV"), type="character", default=opt$ctlCSV, 
                help="Standard Pre-Defined Infinium Methylation Controls CSV (no-header) [default= %default]", metavar="character"),
    
    # Platform/Method Options::
    make_option(c("--genomeBuild"), type="character", default=opt$genomeBuild, 
                help="Genome Build (e.g. hg18, hg36, hg19, hg37, hg38, mm10) [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform (e.g. HM450, EPIC, LEGX, NZT, COVIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest Version (e.g. B0,B1,B2,B3,B4,C0) [default= %default]", metavar="character"),
    
    # Executables::
    make_option(c("--Rscript"), type="character", default=opt$Rscript, 
                help="Rscript path [default= %default]", metavar="character"),
    
    # verbose::
    make_option(c("-v", "--verbose"), type="integer", default=opt$verbose, 
                help="0-5 (5 is very verbose) [default= %default]", metavar="integer")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Validate Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (is.null(par$runMode) || is.null(par$prgmTag) || is.null(par$scrDir) || is.null(par$datDir)) {
  
  par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
  par_tib %>% base::print(n=base::nrow(par_tib) )
  
  if (is.null(par$runMode)) cat(glue::glue("[Usage]: runMode is NULL!!!{RET}"))
  if (is.null(par$prgmTag)) cat(glue::glue("[Usage]: prgmTag is NULL!!!{RET}"))
  if (is.null(par$scrDir))  cat(glue::glue("[Usage]: scrDir is NULL!!!{RET}"))
  if (is.null(par$datDir))  cat(glue::glue("[Usage]: darDir is NULL!!!{RET}"))
  base::stop("Null Parameters!\n\n")
}

if (is.null(opt$outDir) || 
    is.null(opt$Rscript) ||
    is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )
  
  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  
  if (is.null(opt$outDir))    cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  
  if (is.null(opt$Rscript)) cat(glue::glue("[Usage]: Rscript is NULL!!!{RET}"))
  if (is.null(opt$verbose)) cat(glue::glue("[Usage]: verbose is NULL!!!{RET}"))
  base::stop(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
}
par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
if (opt$verbose>=1) par_tib %>% base::print(n=base::nrow(par_tib) )
if (opt$verbose>=1) opt_tib %>% base::print(n=base::nrow(opt_tib) )

cat(glue::glue("[{par$prgmTag}]: Done. Validating Options.{RET}{RET}"))

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

# Load All other function methods::
par$man_src_dir <- file.path(par$scrDir, 'manifests/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$man_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$man_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$man_src_dir}!{RET}{RET}") )

par$swt_src_dir <- file.path(par$scrDir, 'swifthoof/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$swt_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$swt_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$swt_src_dir}!{RET}{RET}") )

par$prb_src_dir <- file.path(par$scrDir, 'probe_design/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$prb_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$prb_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$prb_src_dir}!{RET}{RET}") )

par$anl_src_dir <- file.path(par$scrDir, 'analysis/functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: Manifest Source={par$anl_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$anl_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form Manifest Source={par$anl_src_dir}!{RET}{RET}") )

cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Manifests::
#  - [Done] HM27
#  - [Done] HM450
#  - [Done] EPIC (B0/B2/B4)
#
#  - Mpanel1
#
#  - NZT (not needed yet)
#  - [Done] Excalibur (15k, 5k)
#
#  - [Done] Horvath
#  - mm10
#
#  - Geneknowme
#  - University Chicago
#

# Special Checks for 27k::
#
hm27_man_csv <- '/Users/bbarnes/Documents/Projects/manifests/current/HumanMethylation27_270596_v.1.2.csv.gz'
hm27_man_tib <- loadManifestGenomeStudio(hm27_man_csv, addSource=TRUE, retType='man', normalize=TRUE, verbose=opt$verbose) %>% 
  dplyr::mutate(Seq_ID=IlmnID, Forward_Sequence=Top_Sequence) %>% dplyr::select(Seq_ID, Probe_Type, Man_Source,Forward_Sequence) 

hm27_cgn_tib <- hm27_man_tib %>% 
  dplyr::filter(Probe_Type=='cg') %>%
  dplyr::filter(!is.na(Forward_Sequence)) %>%
  setTopBot_tib(seqKey='Forward_Sequence', srdKey='Forward_StrandTB', topKey='Top_Sequence', verbose=opt$verbose)

hm27_uniq_cgnTop_tib <- hm27_cgn_tib %>% 
  dplyr::select(Seq_ID, Man_Source, Top_Sequence) %>% dplyr::distinct() %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count')

hm27_uniq_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
hm27_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1)

# Special Checks for 450k:: CG
#  - cgn
#  - chn
hm45_man_csv <- '/Users/bbarnes/Documents/Projects/manifests/current/HumanMethylation450_15017482_v.1.2.csv.gz'
hm45_man_tib <- loadManifestGenomeStudio(hm45_man_csv, addSource=TRUE, retType='man', normalize=TRUE, verbose=opt$verbose) %>%
  dplyr::mutate(Seq_ID=IlmnID) %>% dplyr::select(Seq_ID, Probe_Type, Man_Source,Forward_Sequence)

hm45_cgn_tib <- hm45_man_tib %>% 
  dplyr::filter(Probe_Type=='cg') %>%
  dplyr::filter(!is.na(Forward_Sequence)) %>%
  setTopBot_tib(seqKey='Forward_Sequence', srdKey='Forward_StrandTB', topKey='Top_Sequence', verbose=opt$verbose)

hm45_uniq_cgnTop_tib <- hm45_cgn_tib %>% 
  dplyr::select(Seq_ID, Man_Source, Top_Sequence) %>% dplyr::distinct() %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count')

hm45_uniq_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
hm45_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1)

# Special Checks for 450k:: CH
#  - cgn
#  - chn
if (FALSE) {
  hm45_chn_tib <- hm45_man_tib %>% 
    dplyr::filter(Probe_Type=='ch') %>% head() %>%
    dplyr::filter(!is.na(Forward_Sequence)) %>%
    dplyr::mutate(Forward_Sequence_CH=Forward_Sequence) %>%
    dplyr::mutate(Forward_Sequence=stringr::str_replace(Forward_Sequence, '\\[[A-Z][A-Z]\\]', "[CG]") ) %>% 
    setTopBot_tib(seqKey='Forward_Sequence', srdKey='Forward_StrandTB', topKey='Top_Sequence', verbose=opt$verbose)
  
  hm45_uniq_chnTop_tib <- hm45_chn_tib %>% 
    dplyr::select(Seq_ID, Man_Source, Top_Sequence) %>% dplyr::distinct() %>%
    dplyr::add_count(Seq_ID, name='CGN_Count') %>%
    dplyr::add_count(Top_Sequence, name='Top_Count')
  
  hm45_uniq_chnTop_tib %>% dplyr::filter(CGN_Count!=1)
  hm45_uniq_chnTop_tib %>% dplyr::filter(Top_Count!=1)
}

# Special Checks for EPIC:: CG
#  - cgn
#  - chn
hm85_man_csv <- '/Users/bbarnes/Documents/Projects/manifests/current/MethylationEPIC_v-1-0_B2.csv.gz'
hm85_man_tib <- loadManifestGenomeStudio(hm85_man_csv, addSource=TRUE, retType='man', normalize=TRUE, verbose=opt$verbose) %>%
  dplyr::mutate(Seq_ID=IlmnID) %>% dplyr::select(Seq_ID, Probe_Type, Man_Source, Forward_Sequence)

hm85_cgn_tib <- hm85_man_tib %>% 
  dplyr::filter(Probe_Type=='cg') %>%
  dplyr::filter(!is.na(Forward_Sequence)) %>%
  setTopBot_tib(seqKey='Forward_Sequence', srdKey='Forward_StrandTB', topKey='Top_Sequence', verbose=opt$verbose)

hm85_uniq_cgnTop_tib <- hm85_cgn_tib %>% 
  dplyr::select(Seq_ID, Man_Source, Top_Sequence) %>% dplyr::distinct() %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count')

hm85_uniq_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
hm85_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1)

# Special Checks for EPIC:: CH
#  - cgn
#  - chn
if (FALSE) {
  hm85_chn_tib <- hm85_man_tib %>% 
    dplyr::filter(Probe_Type=='ch') %>% head() %>%
    dplyr::filter(!is.na(Forward_Sequence)) %>%
    dplyr::mutate(Forward_Sequence_CH=Forward_Sequence) %>%
    dplyr::mutate(Forward_Sequence=stringr::str_replace(Forward_Sequence, '\\[[A-Z][A-Z]\\]', "[CG]") ) %>% 
    setTopBot_tib(seqKey='Forward_Sequence', srdKey='Forward_StrandTB', topKey='Top_Sequence', verbose=opt$verbose)
  
  hm85_uniq_chnTop_tib <- hm85_chn_tib %>% 
    dplyr::select(Seq_ID, Man_Source, Top_Sequence) %>% dplyr::distinct() %>%
    dplyr::add_count(Seq_ID, name='CGN_Count') %>%
    dplyr::add_count(Top_Sequence, name='Top_Count')
  
  hm85_uniq_chnTop_tib %>% dplyr::filter(CGN_Count!=1)
  hm85_uniq_chnTop_tib %>% dplyr::filter(Top_Count!=1)
}

#
# Reduce human core arrays::
#

hm27_uniq_cgnTop_tib
hm45_uniq_cgnTop_tib
hm85_uniq_cgnTop_tib

# Mpanel1::
#
mpan_man_csv <- '/Users/bbarnes/Documents/Projects/manifests/current/Mpanel1.LEGX.CombinedManifestEPIC.manifest.csv.gz'
mpan_man_tib <- loadManifestGenomeStudio(mpan_man_csv, addSource=TRUE, retType='man', normalize=TRUE, verbose=opt$verbose) %>%
  dplyr::mutate(Seq_ID=IlmnID) %>% dplyr::select(Seq_ID, Probe_Type, Man_Source, Forward_Sequence)

mpan_cgn_tib <- mpan_man_tib %>% 
  dplyr::filter(Probe_Type=='cg') %>%
  dplyr::filter(!is.na(Forward_Sequence)) %>%
  setTopBot_tib(seqKey='Forward_Sequence', srdKey='Forward_StrandTB', topKey='Top_Sequence', verbose=opt$verbose)

mpan_uniq_cgnTop_tib <- mpan_cgn_tib %>% 
  dplyr::select(Seq_ID, Man_Source, Top_Sequence) %>% dplyr::distinct() %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count')

mpan_uniq_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
mpan_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1)


# Excalibur (15k)
#
excb_ann_csv <- '/Users/bbarnes/Documents/Projects/methylation/Excalibur/orders/order.annotation.csv.gz'
excb_ann_tib <- suppressMessages(suppressWarnings( readr::read_csv(excb_ann_csv) )) %>% 
  dplyr::filter(!is.na(Forward_Sequence)) %>%
  setTopBot_tib(seqKey='Forward_Sequence', srdKey='Forward_StrandTB', topKey='Top_Sequence', verbose=opt$verbose) %>% 
  dplyr::mutate(Seq_ID=stringr::str_remove(Assay_Design_Id, '_.*$'))

excb_uniq_cgnTop_tib <- excb_ann_tib %>% 
  dplyr::select(Seq_ID, Top_Sequence) %>% dplyr::distinct() %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count') %>% 
  dplyr::mutate(Man_Source="Excalibur15k") %>% dplyr::select(Seq_ID,Man_Source,Top_Sequence, everything())

excb_uniq_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
excb_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1)

# Horvath Mamalian::
#
horv_man_csv <- '/Users/bbarnes/Documents/Projects/manifests/current/HorvathMammMeth.csv.gz'
horv_man_tib <- loadManifestGenomeStudio(file=horv_man_csv, addSource=TRUE, retType='man', verbose=opt$verbose) %>%
  dplyr::filter(!is.na(Forward_Sequence)) %>%
  setTopBot_tib(seqKey='Forward_Sequence', srdKey='Forward_StrandTB', topKey='Top_Sequence', verbose=opt$verbose) %>% 
  dplyr::mutate(Seq_ID=stringr::str_remove(IlmnID, '_.*$'))

horv_uniq_cgnTop_tib <- horv_man_tib %>% 
  dplyr::select(Seq_ID, Top_Sequence) %>% dplyr::distinct() %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count') %>% 
  dplyr::mutate(Man_Source="Horvath40k") %>% dplyr::select(Seq_ID,Man_Source,Top_Sequence, everything())

horv_uniq_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
horv_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1)

# mm10 forced:: cg+mu
#
mm10_uniq_cgnTop_csv <- '/Users/bbarnes/Documents/Projects/manifests/current/mm10-LEGX-S2.forced-defined-cgnToTopSeq.csv.gz'
mm10_uniq_cgnTop_tib <- suppressMessages(suppressWarnings( readr::read_csv(mm10_uniq_cgnTop_csv) )) %>% 
  dplyr::mutate(Man_Source="mm10-LEGX-S2") %>% dplyr::select(Seq_ID,Man_Source,Top_Sequence, everything())

mm10_uniq_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
mm10_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1)

# University of Chicago::
#
#  - Univ Chicago
#    - Make sure methylation/CustomContent/UnivChicago/improbe_input/CpGs_UnivChicago_alldesigns_55860sites.tsv.gz covers
#       all probes in methylation/CustomContent/UnivChicago/jira-data/CpGs_UnivChicago_alldesigns_55860sites.order.csv
#       and add them to the contigency check...
#
chig_des_tsv <- '/Users/bbarnes/Documents/Projects/methylation/CustomContent/UnivChicago/improbe_input/CpGs_UnivChicago_alldesigns_55860sites.tsv.gz'
chig_ord_csv <- '/Users/bbarnes/Documents/Projects/methylation/CustomContent/UnivChicago/jira-data/CpGs_UnivChicago_alldesigns_55860sites.order.csv'

chig_des_tib <- suppressMessages(suppressWarnings( readr::read_tsv(chig_des_tsv) ))
chig_ord_tib <- suppressMessages(suppressWarnings( readr::read_csv(chig_ord_csv) )) %>% dplyr::mutate(Seq_ID=stringr::str_remove(Assay_Design_Id, '_.*$'))

chig_uniq_cgnTop_tib <- chig_ord_tib %>% dplyr::inner_join(chig_des_tib, by="Seq_ID") %>% 
  dplyr::select(Seq_ID, Top_Sequence) %>% dplyr::distinct() %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count') %>% 
  dplyr::mutate(Man_Source="univ-chicago") %>% dplyr::select(Seq_ID,Man_Source,Top_Sequence, everything())
chig_uniq_cgnTop_tib %>% base::nrow()
chig_ord_tib %>% dplyr::distinct(Seq_ID) %>% base::nrow()

# EPIC B0:: 
#   /Users/bbarnes/Documents/Projects/manifests/current/MethylationEPIC_v-1-0_B0.forced-defined-cgnToTopSeq.tsv.gz
#
epi0_man_tsv <- '/Users/bbarnes/Documents/Projects/manifests/current/MethylationEPIC_v-1-0_B0.forced-defined-cgnToTopSeq.tsv.gz'
epi0_man_tib <- suppressMessages(suppressWarnings( readr::read_tsv(epi0_man_tsv) )) %>%
  dplyr::mutate(Man_Source='MethylationEPIC_v-1-0_B0', Probe_Type=stringr::str_sub(Seq_ID,1,2)) %>%
  dplyr::rename(Forward_Sequence=Top_Sequence) %>%
  dplyr::select(Seq_ID, Probe_Type, Man_Source, Forward_Sequence)

epi0_cgn_tib <- epi0_man_tib %>% 
  dplyr::filter(Probe_Type=='cg') %>%
  dplyr::filter(!is.na(Forward_Sequence)) %>%
  setTopBot_tib(seqKey='Forward_Sequence', srdKey='Forward_StrandTB', topKey='Top_Sequence', verbose=opt$verbose)

epi0_uniq_cgnTop_tib <- epi0_cgn_tib %>% 
  dplyr::select(Seq_ID, Man_Source, Top_Sequence) %>% dplyr::distinct() %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count')

epi0_uniq_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
epi0_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1)



#
# Combine all Current arrays::
#
hm45_uniq_cgnTop_tib
hm85_uniq_cgnTop_tib

excb_uniq_cgnTop_tib
horv_uniq_cgnTop_tib
mm10_uniq_cgnTop_tib
chig_uniq_cgnTop_tib

hm27_uniq_cgnTop_tib
epi0_uniq_cgnTop_tib

#
# Validate All::
#  - add EPIC_B0
#  - write corrections: (mm10, chig)
#  - write black-list cgns from 27k
#
#  - add COVIC ???
#  - add NZT ???
#
all_full_cgnTop_tib <- dplyr::bind_rows(hm45_uniq_cgnTop_tib,hm85_uniq_cgnTop_tib,
                                        excb_uniq_cgnTop_tib,horv_uniq_cgnTop_tib,
                                        mm10_uniq_cgnTop_tib,chig_uniq_cgnTop_tib,
                                        hm27_uniq_cgnTop_tib,epi0_uniq_cgnTop_tib) %>%
  dplyr::select(Seq_ID,Man_Source,Top_Sequence)

all_tops_cgnTop_tib <- all_full_cgnTop_tib %>% 
  dplyr::distinct(Seq_ID, .keep_all=TRUE) %>%
  dplyr::distinct(Top_Sequence, .keep_all=TRUE) %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count') %>% 
  dplyr::add_count(Seq_ID,Top_Sequence, name="CGN_Top_Count")

all_tops_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
all_tops_cgnTop_tib %>% dplyr::filter(Top_Count!=1)
all_tops_cgnTop_tib %>% dplyr::filter(CGN_Top_Count!=1)

all_tops_cgnTop_tib %>% dplyr::group_by(Man_Source) %>% dplyr::summarise(Src_Cnt=n()) %>% dplyr::arrange(-Src_Cnt)

fin_tops_cgnTop_tib <- all_tops_cgnTop_tib %>% 
  dplyr::select(Seq_ID, Man_Source, Top_Sequence) %>% 
  dplyr::mutate(CGN=stringr::str_remove(Seq_ID,'cg') %>% stringr::str_remove('^0+'),
                Top_Sequence=shearBrac(Top_Sequence)) %>%
  dplyr::rename(TOP=Top_Sequence,SRC=Man_Source) %>%
  dplyr::select(CGN,TOP,SRC)

fin_tops_cgnTop_csv <- file.path(opt$outDir, par$prgmTag, 'pre-assigned.cgnTop.hash.csv.gz')
readr::write_csv(fin_tops_cgnTop_tib, fin_tops_cgnTop_csv)


# Pre epi0::
# Man_Source                         Src_Cnt
# 1 HumanMethylation450_15017482_v.1.2  482421
# 2 MethylationEPIC_v-1-0_B2            413743
# 3 mm10-LEGX-S2                        279051
# 4 univ-chicago                         37043
# 5 Horvath40k                           29769
# 6 HumanMethylation27_270596_v.1.2       1135

# Pos epi0::
# Man_Source                         Src_Cnt
# 1 HumanMethylation450_15017482_v.1.2  482421
# 2 MethylationEPIC_v-1-0_B2            413743
# 3 mm10-LEGX-S2                        279051
# 4 univ-chicago                         37043
# 5 MethylationEPIC_v-1-0_B0             36810
# 6 Horvath40k                           29769
# 7 HumanMethylation27_270596_v.1.2       1135





# Old Testing Method without ordering::
#
all_tops_cgnTop_tib <- all_full_cgnTop_tib %>% 
  dplyr::distinct(Seq_ID,Top_Sequence) %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count')

all_tops_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
all_tops_cgnTop_tib %>% dplyr::filter(Top_Count!=1)

# Failures::
all_fail_uniq_cgnTop_tib <- all_tops_cgnTop_tib %>% dplyr::filter(CGN_Count!=1 | Top_Count!=1)

hm45_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% all_fail_uniq_cgnTop_tib$Top_Sequence)
hm85_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% all_fail_uniq_cgnTop_tib$Top_Sequence)
epi0_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% all_fail_uniq_cgnTop_tib$Top_Sequence)

excb_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% all_fail_uniq_cgnTop_tib$Top_Sequence)
horv_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% all_fail_uniq_cgnTop_tib$Top_Sequence)
mm10_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% all_fail_uniq_cgnTop_tib$Top_Sequence)
chig_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% all_fail_uniq_cgnTop_tib$Top_Sequence)
hm27_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% all_fail_uniq_cgnTop_tib$Top_Sequence)


# Remove 27k failures::
#  - fix 27k failures
#    - write black-list cgns from 27k
#
# Conclusion:: Pick 7 and black list 7
#
fx27_full_cgnTop_tib <- dplyr::bind_rows(hm27_uniq_cgnTop_tib,hm45_uniq_cgnTop_tib,hm85_uniq_cgnTop_tib) %>%
  dplyr::select(Seq_ID,Man_Source,Top_Sequence)

fx27_uniq_cgnTop_tib <- fx27_full_cgnTop_tib %>% 
  dplyr::distinct(Seq_ID,Top_Sequence) %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count')

fx27_uniq_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
fx27_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1)

fx27_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1) %>% dplyr::inner_join(fx27_full_cgnTop_tib, by="Top_Sequence")

fx27_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1) %>% dplyr::distinct(Top_Sequence) %>%
  dplyr::left_join(fx27_full_cgnTop_tib, by="Top_Sequence") %>% 
  dplyr::distinct(Top_Sequence,Man_Source)
  # dplyr::group_by(Man_Source) %>% dplyr::summarise(Src_Count=n())

no27_full_cgnTop_tib %>% dplyr::inner_join(
  dplyr::filter(fx27_uniq_cgnTop_tib,Top_Count!=1) %>% dplyr::distinct(Seq_ID),
  by="Seq_ID"
)

#
# Current Steps::
#
#  - fix 27k failures
#    - write black-list cgns from 27k
#  - add COVIC
#  - add EPIC_B0
#  - Univ Chicago
#    - Make sure methylation/CustomContent/UnivChicago/improbe_input/CpGs_UnivChicago_alldesigns_55860sites.tsv.gz covers
#       all probes in methylation/CustomContent/UnivChicago/jira-data/CpGs_UnivChicago_alldesigns_55860sites.order.csv
#       and add them to the contigency check...
#  - validate Geneknowme
#  - order current cgnTop tibs by importance and distinct on cgn order
#    - Horvath before mm10 & Univ-Chicago
#    - Last 27k
#

#
# Validate all non-27k CGN->Top assignments::
#
no27_full_cgnTop_tib <- dplyr::bind_rows(hm45_uniq_cgnTop_tib,hm85_uniq_cgnTop_tib,
                                         excb_uniq_cgnTop_tib,horv_uniq_cgnTop_tib,mm10_uniq_cgnTop_tib,
                                         chig_uniq_cgnTop_tib) %>%
  dplyr::select(Seq_ID,Man_Source,Top_Sequence)

no27_uniq_cgnTop_tib <- no27_full_cgnTop_tib %>% 
  dplyr::distinct(Seq_ID,Top_Sequence) %>%
  dplyr::add_count(Seq_ID, name='CGN_Count') %>%
  dplyr::add_count(Top_Sequence, name='Top_Count')

no27_uniq_cgnTop_tib %>% dplyr::filter(CGN_Count!=1)
no27_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1)

# Failures::
no27_fail_uniq_cgnTop_tib <- no27_uniq_cgnTop_tib %>% dplyr::filter(Top_Count!=1)

hm45_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% no27_fail_uniq_cgnTop_tib$Top_Sequence)
hm85_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% no27_fail_uniq_cgnTop_tib$Top_Sequence)
excb_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% no27_fail_uniq_cgnTop_tib$Top_Sequence)
horv_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% no27_fail_uniq_cgnTop_tib$Top_Sequence)
mm10_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% no27_fail_uniq_cgnTop_tib$Top_Sequence)
chig_uniq_cgnTop_tib %>% dplyr::filter(Top_Sequence %in% no27_fail_uniq_cgnTop_tib$Top_Sequence)

#
# Resolve Horvath and mm10::
#
# Conclusion:: Need to use the Horvath Assignments!!!
#
if (FALSE) {
  opt$outDir <- file.path(opt$outDir, par$prgmTag)
  if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
  
  err_csv <- file.path(opt$outDir, 'error-cgnToTopSeq.csv.gz')
  err_tib <- dplyr::inner_join(
    dplyr::filter(horv_uniq_cgnTop_tib,Top_Sequence %in% no27_fail_uniq_cgnTop_tib$Top_Sequence) %>% dplyr::select(Seq_ID,Man_Source,Top_Sequence),
    dplyr::filter(mm10_uniq_cgnTop_tib,Top_Sequence %in% no27_fail_uniq_cgnTop_tib$Top_Sequence) %>% dplyr::select(Seq_ID,Man_Source,Top_Sequence),
    by="Top_Sequence", suffix=c("_horv","_mm10"))
  readr::write_csv(err_tib, err_csv)
}



  
# Geneknowme::
#
genk_aqpDir <- '/Users/bbarnes/Documents/Projects/methylation/CustomContent/Genknowme/LS_Epiprofile'
genk_ords  <- paste(
  file.path(genk_aqpDir, 'AQP1_NAremoved_GenKnowme_CpG_SNP_order.07082020.csv'),
  file.path(genk_aqpDir, 'AQP2_AP_Genknowme_AQP2_replicate_design_file2.csv'),
  sep=',')
genk_ord_vec <- genk_ords %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

genk_ord_tib <- lapply(genk_ord_vec, readr::read_csv) %>% 
  dplyr::bind_rows() %>% dplyr::mutate(Seq_ID=stringr::str_remove(Assay_Design_Id, '_.*$')) %>% 
  dplyr::filter(stringr::str_starts(Seq_ID, 'cg') )

# Just ensure all cgs are in core_uniq_cgnTop_tib
genk_ord_tib %>% dplyr::filter(Seq_ID %in% all_full_cgnTop_tib$Seq_ID)

#
# TBD:: Next ensure by seq48U matching::
#

# OLD::
#  genk_ord_tib %>% dplyr::filter(Seq_ID %in% core_uniq_cgnTop_tib$Seq_ID)
#
# genk_ord_tib %>% dplyr::filter(Seq_ID %in% core_uniq_cgnTop_tib$Seq_ID)
# genk_ord_tib %>% dplyr::filter(!Seq_ID %in% core_uniq_cgnTop_tib$Seq_ID) %>% 
#   dplyr::filter(Seq_ID %in% horv_uniq_cgnTop_tib$Seq_ID)



# NZT::
#
# nzt_full_man_csv <- file.path(par$datDir, 'manifest/full/NZT-B1.manifest.sesame-full.cpg-sorted.csv.gz')
# nzt_full_add_csv <- file.path(par$datDir, 'manifest/full/NZT-B1.manifest.address-full.cpg-sorted.csv.gz')
# 
# nzt_full_man_tib <- suppressMessages(suppressWarnings( readr::read_csv(nzt_full_man_csv) ))
# nzt_full_add_tib <- suppressMessages(suppressWarnings( readr::read_csv(nzt_full_add_csv) ))


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))


# End of file
