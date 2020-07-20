
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
opt$outDir    <- NULL
opt$buildDir  <- NULL

# Run Parameters::
opt$runName   <- NULL
opt$sampleCsv <- NULL

# Class Parameters::
# Really simple test to make sure we can seperate the sexes...
opt$classVar <- 'Karyotype_0_Call'
opt$classVar <- 'Sample_Name'
opt$classVar <- 'Sample_Class'

# Sample Sheet Parameters::
opt$addSampleName    <- TRUE
opt$addPathsCall     <- TRUE
opt$addPathsSigs     <- FALSE

opt$flagDetectPval   <- FALSE
opt$flagSampleDetect <- FALSE
opt$flagRefMatch     <- FALSE

opt$pvalDetectMinKey <- NULL
opt$pvalDetectMinVal <- NULL

# Chip Platform and Version Parameters::
opt$platform <- NULL
opt$version  <- NULL

# Output Format Parameters::
opt$percisionBeta <- 4
opt$percisionPval <- 6

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
  
  if (dir.exists(par$macDir)) par$topDir <- par$macDir
  if (dir.exists(par$lixDir)) par$topDir <- par$lixDir
  
  # Default Parameters for local Mac::
  par$runMode    <- args.dat[1]
  par$srcDir     <- file.path(par$topDir, 'workhorse')
  par$scrDir     <- file.path(par$srcDir, 'scripts')
  par$exePath    <- file.path(par$scrDir, 'R', par$prgmDir, paste0(par$prgmTag,'.R'))
  
  par$prgmTag <- base::sub('\\.R$', '', base::basename(par$exePath))
  par$locPath <- base::dirname(par$exePath)
  par$scrDir  <- base::dirname(base::normalizePath(par$locPath) )
  par$srcDir  <- base::dirname(base::normalizePath(par$scrDir) )
  par$datDir  <- file.path(par$srcDir, 'dat')
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  
  opt$platform  <- 'EPIC'
  opt$version   <- 'C0'
  
  opt$classVar <- 'Karyotype_0_Call'
  opt$classVar <- 'Karyotype_1_Call'
  opt$classVar <- 'Sample_Name'
  opt$classVar <- 'Sample_Class'
  
  opt$clean  <- TRUE
  opt$clean  <- FALSE
  
  if (opt$classVar=='Sample_Class') {
    opt$runName1  <- 'COVIC-Set1-15052020'
    opt$runName5  <- 'COVIC-Set5-10062020'
    opt$runName   <- opt$runName1
    
    opt$buildDir  <- paste(
      file.path(par$topDir, 'builds/swifthoof_main', opt$runName1),
      # file.path(par$topDir, 'builds/swifthoof_main', opt$runName5),
      sep=',')
    
    opt$sampleCsv <- file.path(par$topDir, 'sampleSheets/annotation/Human-Classification_COVID_Count-656_AnnotatedMultiSampleSheet.csv')
    opt$trainClass <- paste('nSARSCov2', 'pSARSCov2', sep=',')
    
  } else if (opt$classVar=='Karyotype_0_Call' || opt$classVar=='Karyotype_1_Call') {
    opt$runNameA  <- "COVIC-Set1-15052020"
    opt$runNameB  <- "COVIC-Set5-10062020"
    
    opt$runName  <- 'COVIC-Set1-5'
    opt$runName  <- opt$runNameB
    opt$runName  <- opt$runNameA
    
    opt$mergeDir <- paste(
      file.path(par$topDir, 'merge_builds',opt$classVar,opt$runName,opt$platform,opt$version ),
      sep=',')
    
    # opt$trainClass <- paste('Xa','XaXaY','XaXi','XaXiY','XaY', sep=',')
    opt$trainClass <- paste('XaXi','XaY', sep=',')
    
  } else if (opt$classVar=='Sample_Name') {
    opt$runNameA  <- 'DELTA-8x1-EPIC-Core'
    opt$runNameB  <- 'BETA-8x1-EPIC-Core'
    opt$runName   <- 'BETA-DELTA-8x1-EPIC-Core'
    
    opt$version   <- 'B4'
    
    opt$buildDir  <- paste(
      file.path(par$topDir, 'builds/swifthoof_main', opt$runNameA),
      file.path(par$topDir, 'builds/swifthoof_main', opt$runNameB),
      sep=',')
    
    # readr::write_csv( auto_ss_tib %>% dplyr::select(Sentrix_Name, Auto_Sample_Name) %>% dplyr::rename(Sample_Name=Auto_Sample_Name),
    #                   "/Users/bbarnes/Documents/CustomerFacing/sampleSheets/BETA/BETA-8x1-EPIC-Core.Sample_Names.SampleSheet.csv.gz")
    # opt$sampleCsv <- file.path(par$topDir, 'sampleSheets/BETA', paste(opt$runName, 'Sample_Names.SampleSheet.csv.gz', sep='.'))
    
    # readr::write_csv( auto_ss_tib %>% dplyr::select(Sentrix_Name, Auto_Sample_Name) %>% dplyr::rename(Sample_Name=Auto_Sample_Name),
    #                   "/Users/bbarnes/Documents/CustomerFacing/sampleSheets/DELTA/DELTA-8x1-EPIC-Core.Sample_Names.SampleSheet.csv.gz")
    # opt$sampleCsv <- file.path(par$topDir, 'sampleSheets/DELTA', paste(opt$runName, 'Sample_Names.SampleSheet.csv.gz', sep='.'))
    opt$sampleCsv <- file.path(par$topDir, 'sampleSheets/BETA-DELTA-EPIC-Core/BETA-DELTA-8x1-EPIC-Core.Sample_Names.SampleSheet.csv.gz')
    
    opt$trainClass <- paste('HELA','JURKAT','MCF7','RAJI', sep=',')
  }
  opt$outDir <- file.path(par$topDir, par$prgmTag)
  
} else {
  par$runMode    <- 'CommandLine'
  par$exePath <- base::substring(args.dat[grep("--file=", args.dat)], 8)
  
  par$prgmTag <- base::sub('\\.R$', '', base::basename(par$exePath))
  par$locPath <- base::dirname(par$exePath)
  par$scrDir  <- base::dirname(base::normalizePath(par$locPath) )
  par$srcDir  <- base::dirname(base::normalizePath(par$scrDir) )
  par$datDir  <- file.path(par$srcDir, 'dat')
  
  args.dat <- commandArgs(trailingOnly = TRUE)
  option_list = list(
    
    # Directory Parameters::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    make_option(c("--buildDirs"), type="character", default=opt$buildDir, 
                help="List of Build Directory [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--sampleCsv"), type="character", default=opt$sampleCsv, 
                help="Human provide sample sheet labeling [default= %default]", metavar="character"),
    
    # Chip Platform and Version Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform name (HM50, EPIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest version (B2, B4, C0) [default= %default]", metavar="character"),
    
    # Class Parameters::
    make_option(c("--classVar"), type="character", default=opt$classVar, 
                help="Classification Variable Name [default= %default]", metavar="character"),
    
    # Sample Sheet Parameters::
    make_option(c("--addSampleName"), action="store_true", default=opt$addSampleName, 
                help="Sample Sheet processing to add Auto-SampleNames (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--addPathsCall"), action="store_true", default=opt$addPathsCall, 
                help="Sample Sheet processing to add Calls Local Full Paths (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--addPathsSigs"), action="store_true", default=opt$addPathsSigs, 
                help="Sample Sheet processing to add Signals Local Full Paths (mostly testing stuff) [default= %default]", metavar="boolean"),
    
    make_option(c("--flagDetectPval"), action="store_true", default=opt$flagDetectPval, 
                help="Sample Sheet processing to add flag for failed detected samples by failed loci percent (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--flagSampleDetect"), action="store_true", default=opt$flagSampleDetect, 
                help="Sample Sheet processing to add flag for failed Auto-Detected Samples (mostly testing stuff) [default= %default]", metavar="boolean"),
    make_option(c("--flagRefMatch"), action="store_true", default=opt$flagRefMatch, 
                help="Sample Sheet processing to add flag for failed Auto-Detected Methods Agreements (mostly testing stuff) [default= %default]", metavar="boolean"),
    
    make_option(c("--pvalDetectMinKey"), type="character", default=opt$pvalDetectMinKey, 
                help="Sample Sheet processing to Min Detection Pval Key [default= %default]", metavar="character"),
    make_option(c("--pvalDetectMinVal"), type="double", default=opt$pvalDetectMinVal,
                help="Sample Sheet processing to Min Detection Pval Value [default= %default]", metavar="double"),
    
    # Output Format Parameters::
    make_option(c("--percisionBeta"), type="integer", default=opt$percisionBeta,
                help="Rounding percision for beta values in calls output files [default= %default]", metavar="integer"),
    make_option(c("--percisionPval"), type="integer", default=opt$percisionPval,
                help="Rounding percision for detection p-values in calls output files [default= %default]", metavar="integer"),
    
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

if (is.null(par$runMode) || is.null(par$prgmTag) || is.null(par$scrDir) || is.null(par$datDir)) {
  
  par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
  par_tib %>% base::print(n=base::nrow(par_tib) )
  
  if (is.null(par$runMode)) cat(glue::glue("[Usage]: runMode is NULL!!!{RET}"))
  if (is.null(par$prgmTag)) cat(glue::glue("[Usage]: prgmTag is NULL!!!{RET}"))
  if (is.null(par$scrDir))  cat(glue::glue("[Usage]: scrDir is NULL!!!{RET}"))
  if (is.null(par$datDir))  cat(glue::glue("[Usage]: darDir is NULL!!!{RET}"))
  base::stop("Null Parameters!\n\n")
}

if (is.null(opt$outDir) || is.null(opt$buildDir) || 
    is.null(opt$runName) || is.null(opt$sampleCsv) || 
    is.null(opt$classVar) ||
    is.null(opt$addSampleName) || is.null(opt$addPathsCall) || is.null(opt$addPathsSigs) ||
    is.null(opt$flagDetectPval) || is.null(opt$flagSampleDetect) || is.null(opt$flagRefMatch) ||
    # is.null(opt$pvalDetectMinKey) || is.nulll(opt$pvalDetectMinVal) ||
    is.null(opt$platform) || is.null(opt$version) ||
    is.null(opt$percisionBeta) || is.null(opt$percisionPval) ||
    is.null(opt$execute) || is.null(opt$single) || is.null(opt$parallel) || is.null(opt$cluster) ||
    
    is.null(opt$clean) || is.null(opt$Rscript) || is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )
  
  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  if (is.null(opt$outDir))    cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$buildDir))  cat(glue::glue("[Usage]: buildDirs is NULL!!!{RET}"))
  if (is.null(opt$runName))   cat(glue::glue("[Usage]: runName is NULL!!!{RET}"))
  if (is.null(opt$sampleCsv)) cat(glue::glue("[Usage]: sampleCsv is NULL!!!{RET}"))
  if (is.null(opt$classVar))  cat(glue::glue("[Usage]: class_var is NULL!!!{RET}"))
  
  if (is.null(opt$addSampleName)) cat(glue::glue("[Usage]: addSampleName is NULL!!!{RET}"))
  if (is.null(opt$addPathsCall))  cat(glue::glue("[Usage]: addPathsCall is NULL!!!{RET}"))
  if (is.null(opt$addPathsSigs))  cat(glue::glue("[Usage]: addPathsSigs is NULL!!!{RET}"))
  
  if (is.null(opt$flagDetectPval))   cat(glue::glue("[Usage]: flagDetectPval is NULL!!!{RET}"))
  if (is.null(opt$flagSampleDetect)) cat(glue::glue("[Usage]: flagSampleDetect is NULL!!!{RET}"))
  if (is.null(opt$flagRefMatch))     cat(glue::glue("[Usage]: flagRefMatch is NULL!!!{RET}"))
  
  # if (is.null(opt$pvalDetectMinKey)) cat(glue::glue("[Usage]: pvalDetectMinKey is NULL!!!{RET}"))
  # if (is.null(opt$pvalDetectMinVal)) cat(glue::glue("[Usage]: pvalDetectMinVal is NULL!!!{RET}"))
  
  if (is.null(opt$platform)) cat(glue::glue("[Usage]: platform is NULL!!!{RET}"))
  if (is.null(opt$version))  cat(glue::glue("[Usage]: version is NULL!!!{RET}"))
  
  if (is.null(opt$percisionBeta)) cat(glue::glue("[Usage]: percisionBeta is NULL!!!{RET}"))
  if (is.null(opt$percisionPval)) cat(glue::glue("[Usage]: percisionPval is NULL!!!{RET}"))
  
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
#                      Example SNP Data for GenKnowMe::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Workflow::
#
#  1. Format improbe inpute
#  2. Run improbe
#  3. Extract design scores from improbe-output
#  4. De-novo design probe sequecnes and match with design scores
#  5. Write order/annotation file
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         General Order CpGs Only::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$topDir <- '/Users/bbarnes/Documents/Projects/methylation/CustomContent'
opt$topDir <- '/Users/bbarnes/Documents/CustomerFacing/CustomContent'
ord_out_dir <- file.path(opt$topDir, 'orders')
if (!dir.exists(ord_out_dir)) dir.create(ord_out_dir, recursive=TRUE)

opt$runName <- 'CpGs_UnivChicago_alldesigns_55860sites'
fin_sum_csv <- file.path(opt$topDir, paste(opt$runName, 'order.csv', sep='.') )
fin_ord_csv <- file.path(opt$topDir, paste(opt$runName, 'order-summary.csv', sep='.') )

des_ord_tsv <- file.path(opt$topDir, 'data', paste(opt$runName,'ForwardSequence.tsv.gz', sep='.') )
des_ord_tib <- suppressMessages(suppressWarnings( readr::read_tsv(des_ord_tsv)) ) %>% 
  dplyr::mutate(Probe_Type='cg')

new_prb_tib <- tib2prbs(tib=des_ord_tib, idsKey="Seq_ID", prbKey="Probe_Type", seqKey="Forward_Sequence", verbose=opt$verbose)
new_ord_tib <- prbs2order(new_prb_tib, verbose=opt$verbose) %>% 
  dplyr::bind_rows() %>% 
  dplyr::filter(Valid_Design_Bool) %>% dplyr::select(Assay_Design_Id:Normalization_Bin)














# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         GenKnowMe Order CpGs Only::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fin_out_dir <- '/Users/bbarnes/Documents/Projects/methylation/CustomContent/Genknowme/orders/final.07072020'
if (!dir.exists(ord_out_dir)) dir.create(ord_out_dir, recursive=TRUE)

fin_sum_csv <- file.path(fin_out_dir, 'GenKnowme_CpG_SNP_summary.07072020.csv')
fin_ord_csv <- file.path(fin_out_dir, 'GenKnowme_CpG_SNP_order.07072020.csv')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         GenKnowMe Order CpGs Only::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cpg_ord_csv <- '/Users/bbarnes/Documents/Projects/methylation/CustomContent/Genknowme/orders/v2/GenKnowme_New_CpG-ONLY_v2_DesignFile.csv'
cpg_ord_tib <- suppressMessages(suppressWarnings( readr::read_csv(cpg_ord_csv)) ) %>% 
  dplyr::select(Assay_Design_Id:Normalization_Bin) %>%
  dplyr::filter(Assay_Design_Id!='cg02613386')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      GenKnowMe Order New SNPs Only::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

snp_new_csv <- '/Users/bbarnes/Documents/Projects/methylation/CustomContent/Genknowme/367503_Genknowme_rsIDs_v1-1.score_noHeader.csv'

#  New 9 SNPs
new_snp_vec <- c("rs11583993", "rs12073947", "rs1938454", "rs1985313", "rs2009774", "rs2741308", "rs3757351", "rs76293297", "rs8070212")
snp_new_tib <- suppressMessages(suppressWarnings( readr::read_csv(snp_new_csv) )) %>% 
  dplyr::mutate(Name=stringr::str_remove(Locus_Name, "_.*$")) %>%
  dplyr::filter(Name %in% new_snp_vec) %>%
  dplyr::rename(Seq_ID=Name, Genome_Build=Genome_Build_Version, Forward_Sequence=Sequence) %>%
  dplyr::mutate(
    SNP=paste0("[", Forward_Sequence %>% stringr::str_remove("^.*\\[") %>% stringr::str_remove("\\].*$"), "]" ),
    Forward_Sequence=stringr::str_to_upper(Forward_Sequence),
    SNP_IUPAC=dplyr::case_when(
      SNP=="[A/G]" | SNP=="[G/A]" ~ MAPDi("AG"),
      SNP=="[A/C]" | SNP=="[C/A]" ~ MAPDi("AC"),
      SNP=="[C/T]" | SNP=="[T/C]" ~ MAPDi("CT"),
      SNP=="[G/T]" | SNP=="[T/G]" ~ MAPDi("GT"),
      TRUE ~ NA_character_),
    Forward_Sequence=dplyr::case_when(
      Sequence_Orientation=='FORWARD' ~ Forward_Sequence,
      Sequence_Orientation=='REVERSE' ~ revCmp(Forward_Sequence),
      TRUE ~ NA_character_
    ),
    preSeq=Forward_Sequence %>% stringr::str_remove("\\[.*$"),
    posSeq=Forward_Sequence %>% stringr::str_remove("^.*\\]") %>% stringr::str_sub(2),
    nxtNuc=Forward_Sequence %>% stringr::str_replace("^[^]]+[]]([A-Z]).*$", "\\$1") %>% stringr::str_remove("\\\\"),
    Forward_Sequence=paste0(preSeq,"[",SNP_IUPAC,nxtNuc,"]",posSeq,"N"),
    diNUC=paste0(SNP_IUPAC,nxtNuc),
    Probe_Type='rs'
  ) %>%
  dplyr::select(Seq_ID, Forward_Sequence, Genome_Build, Chromosome, Coordinate, Probe_Type, diNUC, everything() )
  
snp_man_cnt <- snp_new_tib %>% base::nrow()
# snp_new_tib %>% print(n=snp_man_cnt)
# snp_new_tib %>% head(n=2) %>% as.data.frame()

new_prb_tib <- tib2prbs(tib=snp_new_tib, idsKey="Seq_ID", prbKey="Probe_Type", seqKey="Forward_Sequence", verbose=opt$verbose)
new_ord_tib <- prbs2order(new_prb_tib, verbose=opt$verbose) %>% 
  dplyr::bind_rows() %>% 
  dplyr::filter(Valid_Design_Bool) %>% dplyr::select(Assay_Design_Id:Normalization_Bin)

new_all_tib <- prbs2order(new_prb_tib, verbose=opt$verbose) %>% dplyr::bind_rows()

# Summary::
new_sum_tib <- new_ord_tib %>% dplyr::mutate(Name=stringr::str_remove(Assay_Design_Id, "_.*$")) %>% dplyr::group_by(Name) %>% dplyr::summarise(Count=n())

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                HM450k Original SNPs (64) Formatting Only::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

snp_val_csv <- '/Users/bbarnes/Documents/Projects/manifests/methylation/HumanMethylation450_15017482_v.1.2.rs-only.csv.gz'
snp_val_tib <- suppressMessages(suppressWarnings( readr::read_csv(snp_val_csv) ))

org_ord_tib <- snp_val_tib %>% # dplyr::select(Name, AlleleA_ProbeSeq, AlleleB_ProbeSeq, everything()) %>% 
  dplyr::rename(Assay_Design_Id=Name,AlleleA_Probe_Sequence=AlleleA_ProbeSeq,AlleleB_Probe_Sequence=AlleleB_ProbeSeq) %>%
  dplyr::mutate(AlleleA_Probe_Id=paste(Assay_Design_Id,"A", sep="_"), 
                AlleleB_Probe_Id=paste(Assay_Design_Id,"B", sep="_"),
                Normalization_Bin=dplyr::case_when(
                  is.na(Next_Base) ~ 'C',
                  Next_Base=='A' | Next_Base=='T' ~ 'A',
                  Next_Base=='C' | Next_Base=='G' ~ 'B',
                  TRUE ~ NA_character_
                ),
                AlleleB_Probe_Id=dplyr::case_when(
                  is.na(AlleleB_Probe_Sequence) ~ '',
                  TRUE ~ AlleleB_Probe_Id
                ),
                AlleleB_Probe_Sequence=dplyr::case_when(
                  is.na(AlleleB_Probe_Sequence) ~ '',
                  TRUE ~ AlleleB_Probe_Sequence
                )
  ) %>% dplyr::select(Assay_Design_Id, AlleleA_Probe_Id, AlleleA_Probe_Sequence, AlleleB_Probe_Id, AlleleB_Probe_Sequence, Normalization_Bin)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Join and Sort Final Order::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fin_ord_tib <- dplyr::bind_rows(
  cpg_ord_tib,
  new_ord_tib,
  org_ord_tib
) %>% dplyr::distinct()


# cpg_ord_tib %>% dplyr::filter(Assay_Design_Id=='cg02613386')
# new_ord_tib %>% dplyr::filter(Assay_Design_Id=='cg02613386')
# org_ord_tib %>% dplyr::filter(Assay_Design_Id=='cg02613386')

fin_ord_tib %>% base::nrow() %>% print()
fin_ord_tib %>% dplyr::distinct(Assay_Design_Id) %>% base::nrow() %>% print()
fin_ord_tib %>% dplyr::distinct(AlleleA_Probe_Id) %>% base::nrow() %>% print()
fin_ord_tib %>% dplyr::distinct(AlleleA_Probe_Sequence) %>% base::nrow() %>% print()
fin_ord_tib %>% dplyr::distinct(AlleleB_Probe_Id) %>% base::nrow() %>% print()
fin_ord_tib %>% dplyr::distinct(AlleleB_Probe_Sequence) %>% base::nrow() %>% print()
fin_ord_tib %>% dplyr::distinct(Normalization_Bin) %>% base::nrow() %>% print()

# fin_ord_tib %>% dplyr::group_by(Assay_Design_Id) %>% dplyr::summarise(Count=n()) %>% dplyr::filter(Count!=1)

fin_sum_tib <- order2stats(fin_ord_tib)

readr::write_csv(fin_sum_tib, fin_sum_csv)
readr::write_csv(fin_ord_tib, fin_ord_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                   DONE::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #










validation_check <- FALSE
if (validation_check) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #         HM450k Original SNPs (64) AND GenKnowMe Order New SNPs Only::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  snp_man_csv <- '/Users/bbarnes/Documents/Projects/methylation/improbe/data/manifests/GenKnowMe.snp.selected.csv'
  snp_man_tib <- suppressMessages(suppressWarnings( readr::read_csv(snp_man_csv) )) %>% 
    dplyr::select(Name, SNP, SourceSeq, GenomeBuild, Chr, MapInfo, RefStrand, SourceStrand) %>%
    dplyr::rename(Seq_ID=Name, Genome_Build=GenomeBuild, Chromosome=Chr, Coordinate=MapInfo) %>%
    dplyr::mutate(
      SourceSeq=stringr::str_to_upper(SourceSeq),
      SNP_IUPAC=dplyr::case_when(
        SNP=="[A/G]" | SNP=="[G/A]" ~ MAPDi("AG"),
        SNP=="[A/C]" | SNP=="[C/A]" ~ MAPDi("AC"),
        SNP=="[C/T]" | SNP=="[T/C]" ~ MAPDi("CT"),
        SNP=="[G/T]" | SNP=="[T/G]" ~ MAPDi("GT"),
        TRUE ~ NA_character_),
      # FR=dplyr::case_when(RefStrand=='+' ~ 'F', RefStrand=='-' ~ 'R', TRUE ~ NA_character_),
      TB_Source=dplyr::case_when(RefStrand=='+' & SourceStrand=='TOP' ~ 'TOP', 
                                 RefStrand=='+' & SourceStrand=='BOT' ~ 'BOT',
                                 
                                 RefStrand=='-' & SourceStrand=='TOP' ~ 'BOT',
                                 RefStrand=='-' & SourceStrand=='BOT' ~ 'TOP',
                                 TRUE ~ NA_character_),
      Forward_Sequence=dplyr::case_when(
        RefStrand=='+' ~ SourceSeq,
        RefStrand=='-' ~ revCmp(SourceSeq),
        TRUE ~ NA_character_
      ),
      preSeq=Forward_Sequence %>% stringr::str_remove("\\[.*$"),
      posSeq=Forward_Sequence %>% stringr::str_remove("^.*\\]") %>% stringr::str_sub(2),
      nxtNuc=Forward_Sequence %>% stringr::str_replace("^[^]]+[]]([A-Z]).*$", "\\$1") %>% stringr::str_remove("\\\\"),
      Forward_Sequence=paste0(preSeq,"[",SNP_IUPAC,nxtNuc,"]",posSeq,"N"),
      diNUC=paste0(SNP_IUPAC,nxtNuc),
      Probe_Type='rs'
    ) %>% dplyr::filter(!is.na(SNP)) %>%
    dplyr::select(Seq_ID, Forward_Sequence, Genome_Build, Chromosome, Coordinate, Probe_Type, diNUC, everything()) %>% 
    dplyr::distinct(Seq_ID, .keep_all=TRUE)
  
  snp_man_cnt <- snp_man_tib %>% base::nrow()
  # snp_man_tib %>% print(n=snp_man_cnt)
  # snp_man_tib %>% head(n=2) %>% as.data.frame()
  
  prb_tib <- tib2prbs(tib=snp_man_tib, idsKey="Seq_ID", prbKey="Probe_Type", seqKey="Forward_Sequence", verbose=opt$verbose)
  ord_tib <- prbs2order(prb_tib, verbose=opt$verbose) %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(Valid_Design_Bool) # %>% dplyr::select(Assay_Design_Id:Normalization_Bin)
  
  all_tib <- prbs2order(prb_tib, verbose=opt$verbose) %>% dplyr::bind_rows()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                 Validate Designs With Previous Designs::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ord_mat_tib <- ord_tib %>% 
    dplyr::mutate(AlleleA_Probe_Sequence=stringr::str_to_upper(AlleleA_Probe_Sequence),
                  AlleleB_Probe_Sequence=stringr::str_to_upper(AlleleB_Probe_Sequence),
                  Name=stringr::str_remove(Assay_Design_Id, "_.*$"))
  
  all_mat_tib <- all_tib %>% 
    dplyr::mutate(AlleleA_Probe_Sequence=stringr::str_to_upper(AlleleA_Probe_Sequence),
                  AlleleB_Probe_Sequence=stringr::str_to_upper(AlleleB_Probe_Sequence),
                  Name=stringr::str_remove(Assay_Design_Id, "_.*$"))
  
  #
  # Almost full match:: Infinium II
  #
  snp_val_cnt2 <- snp_val_tib %>% dplyr::filter(Infinium_Design_Type=="II") %>% dplyr::distinct(Name, .keep_all=TRUE) %>% base::nrow()
  
  ord_mat2_tib <- snp_val_tib %>% dplyr::filter(Infinium_Design_Type=="II") %>% 
    dplyr::inner_join(ord_mat_tib, by="Name") %>% 
    dplyr::filter(AlleleA_ProbeSeq==AlleleA_Probe_Sequence) %>%
    dplyr::select(Assay_Design_Id,AlleleA_ProbeSeq,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence)
  
  all_mat2_tib <- snp_val_tib %>% dplyr::filter(Infinium_Design_Type=="II") %>% 
    dplyr::inner_join(all_mat_tib, by="Name") %>% 
    dplyr::filter(AlleleA_ProbeSeq==AlleleA_Probe_Sequence) %>%
    dplyr::select(Assay_Design_Id,AlleleA_ProbeSeq,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, everything())
  
  # These four may be junk::
  #   rs2125573_F_C_II
  #   rs264581_F_C_II
  #   rs472920_F_C_II
  #   rs3818562_F_C_II
  non_mat2_tib <- all_mat2_tib %>% dplyr::anti_join(ord_mat2_tib, by="Assay_Design_Id") %>% as.data.frame()
  
  #
  # Need to make sure for Infinium I designs A/T -> AlleleA and C/G -> AlleleB
  #
  snp_val_cnt1 <- snp_val_tib %>% dplyr::filter(Infinium_Design_Type=="I") %>% dplyr::distinct(Name, .keep_all=TRUE) %>% base::nrow()
  
  ord_mat1_tib <- snp_val_tib %>% dplyr::filter(Infinium_Design_Type=="I") %>% 
    dplyr::inner_join( dplyr::filter(ord_mat_tib, !is.na(AlleleB_Probe_Sequence)), by="Name") %>%
    dplyr::filter(AlleleA_ProbeSeq==AlleleB_Probe_Sequence & AlleleB_ProbeSeq==AlleleA_Probe_Sequence)
  
  # These four may be junk::
  #   rs715359_R_C_I
  #   rs3936238_F_C_I
  #   rs3936238_R_C_I
  #   rs1520670_R_C_I
  ord_mat1_tib <- snp_val_tib %>% dplyr::filter(Infinium_Design_Type=="I") %>% 
    dplyr::inner_join( dplyr::filter(ord_mat_tib, !is.na(AlleleB_Probe_Sequence)), by="Name") %>%
    dplyr::filter(!stringr::str_detect(Assay_Design_Id, '_O_') ) %>%
    dplyr::filter(AlleleA_ProbeSeq!=AlleleB_Probe_Sequence) %>%
    dplyr::select(Assay_Design_Id,AlleleA_ProbeSeq,AlleleB_Probe_Sequence)
  
  all_mat1_tib <- snp_val_tib %>% dplyr::filter(Infinium_Design_Type=="I") %>% 
    dplyr::inner_join( dplyr::filter(all_mat_tib, !is.na(AlleleB_Probe_Sequence)), by="Name") %>%
    dplyr::filter(!stringr::str_detect(Assay_Design_Id, '_O_') ) %>%
    dplyr::filter(AlleleA_ProbeSeq!=AlleleB_Probe_Sequence) %>%
    dplyr::select(Assay_Design_Id,AlleleA_ProbeSeq,AlleleB_Probe_Sequence)
  
  non_mat1_tib <- all_mat1_tib %>% dplyr::anti_join(ord_mat1_tib, by="Assay_Design_Id") %>% as.data.frame()
  
  
  # Check anit-join...
  
  
  
  
  
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Out of Date Workflow Below::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

  if (FALSE) {
    tar_snp_csv <- '/Users/bbarnes/Documents/Projects/methylation/CustomContent/Genknowme/367503_Genknowme_rsIDs_v1-1.score_noHeader.csv'
    tar_snp_tib <- readr::read_csv(tar_snp_csv)
    
    #
    # Add SNP IUPAC Code and format Design Sequence::
    #
    man_snp_tib <- tar_snp_tib %>% dplyr::mutate(
      Probe_ID=Locus_Name,
      PRB_TB=stringr::str_remove_all(stringr::str_replace(Ilmn_Id, '^.*_([TB])_.*$', '\\$1'), '\\\\'),
      # Probe_ID=stringr::str_remove(Ilmn_Id, '-[0-9]+'),
      # Probe_ID=stringr::str_remove(Probe_ID, '_[0-9]+$'),
      Sequence=stringr::str_to_upper(Sequence),
      Design_Prefix_Seq=stringr::str_remove_all(stringr::str_replace(Sequence, '^(.*)\\[.*$', '\\$1'), '\\\\'),
      Allele_A_Nuc=stringr::str_remove_all(stringr::str_replace(Sequence, '^.*\\[([A-Z]).*$', '\\$1'), '\\\\'),
      Allele_B_Nuc=stringr::str_remove_all(stringr::str_replace(Sequence, '^.*\\[[A-Z]/([A-Z]).*$', '\\$1'), '\\\\'),
      Design_Suffix_Seq=stringr::str_remove_all(stringr::str_replace(Sequence, '^.*\\](.*)$', '\\$1'), '\\\\'),
      IUPAC_Nuc=case_when(
        Allele_A_Nuc=='A' & Allele_B_Nuc=='A' ~ 'A',
        Allele_A_Nuc=='C' & Allele_B_Nuc=='C' ~ 'C',
        Allele_A_Nuc=='G' & Allele_B_Nuc=='G' ~ 'G',
        Allele_A_Nuc=='T' & Allele_B_Nuc=='T' ~ 'T',
        
        Allele_A_Nuc=='A' & Allele_B_Nuc=='G' ~ 'R',
        Allele_A_Nuc=='T' & Allele_B_Nuc=='C' ~ 'Y',
        Allele_A_Nuc=='C' & Allele_B_Nuc=='G' ~ 'S',
        Allele_A_Nuc=='A' & Allele_B_Nuc=='T' ~ 'W',
        Allele_A_Nuc=='T' & Allele_B_Nuc=='G' ~ 'K',
        Allele_A_Nuc=='A' & Allele_B_Nuc=='C' ~ 'M',
        
        Allele_A_Nuc=='G' & Allele_B_Nuc=='A' ~ 'R',
        Allele_A_Nuc=='C' & Allele_B_Nuc=='T' ~ 'Y',
        Allele_A_Nuc=='G' & Allele_B_Nuc=='C' ~ 'S',
        Allele_A_Nuc=='T' & Allele_B_Nuc=='A' ~ 'W',
        Allele_A_Nuc=='G' & Allele_B_Nuc=='T' ~ 'K',
        Allele_A_Nuc=='C' & Allele_B_Nuc=='A' ~ 'M',
        
        TRUE ~ NA_character_
      ),
      Next_Base=stringr::str_sub(Design_Suffix_Seq, 1,1),
      Design_Suffix_Seq=stringr::str_sub(Design_Suffix_Seq, 2),
      DesSeqN=paste0(Design_Prefix_Seq,'[',IUPAC_Nuc,Next_Base,']',Design_Suffix_Seq,'N'),
      PRB_DES='rs'
    ) %>% dplyr::distinct(Probe_ID, .keep_all=TRUE) %>%
      dplyr::select(Probe_ID, PRB_DES, PRB_TB, DesSeqN, IUPAC_Nuc, Design_Prefix_Seq, Allele_A_Nuc, Allele_B_Nuc, Next_Base, Design_Suffix_Seq, Sequence)  
    
    # man_snp_tib %>% dplyr::select(Probe_ID, DesSeqN, IUPAC_Nuc, Design_Prefix_Seq, Allele_A_Nuc,Allele_B_Nuc, Next_Base, Design_Suffix_Seq, Sequence) %>% as.data.frame()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      Calcluate Probes on All Strands::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    des_seq_F_C <- man_snp_tib %>% dplyr::mutate(DesSeqN=shearBrac(DesSeqN) ) %>% dplyr::mutate(FR=TRUE, CO=TRUE)
    des_seq_R_C <- des_seq_F_C %>% dplyr::mutate(FR=!FR,CO=CO, DesSeqN=revCmp(DesSeqN) )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Build All Design Strands::
    bsc_tibs <- NULL
    
    # BSC-Forward-Converted::
    bsc_tibs$FC <- des_seq_F_C %>% dplyr::mutate(
      DesBscU = bscUs(DesSeqN),
      DesBscM = bscMs(DesSeqN),
      DesBscD = bscDs(DesSeqN) )
    
    # BSC-Foward-Opposite::
    bsc_tibs$FO <- bsc_tibs$FC %>% dplyr::mutate(
      FR=FR,CO=!CO, 
      # DesSeqN=revCmp(DesSeqN),
      DesBscU=revCmp(DesBscU),
      DesBscM=revCmp(DesBscM),
      DesBscD=revCmp(DesBscD) )
    
    # BSC-Reverse-Converted::
    bsc_tibs$RC <- des_seq_R_C %>% dplyr::mutate(
      DesBscU = bscUs(DesSeqN),
      DesBscM = bscMs(DesSeqN),
      DesBscD = bscDs(DesSeqN) )
    
    # BSC-Reverse-Opposite::
    bsc_tibs$RO <- bsc_tibs$RC %>% dplyr::mutate(
      FR=FR,CO=!CO, 
      # DesSeqN=revCmp(DesSeqN),
      DesBscU=revCmp(DesBscU),
      DesBscM=revCmp(DesBscM),
      DesBscD=revCmp(DesBscD) )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Compare Manifest with Calculated Probes::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    prb_tibs <- NULL
    prb_tibs <- foreach (srd=names(bsc_tibs), .combine=rbind) %dopar% {
      lapply(split(bsc_tibs[[srd]], bsc_tibs[[srd]]$PRB_DES), desAllPrbs) %>% dplyr::bind_rows()
    } %>% dplyr::mutate(PRB_FR=case_when(FR ~ 'F', TRUE ~ 'R'), 
                        PRB_CO=case_when(CO ~ 'C', TRUE ~ 'O'), 
                        Probe_ID=paste(Probe_ID, PRB_FR, PRB_TB, PRB_CO, sep='_'))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Designs To Orders:: I
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    inf1_ord_tib <- prb_tibs %>% dplyr::filter(PRB1_M!=PRB1_U) %>% dplyr::bind_rows() %>%
      dplyr::mutate(Probe_ID=paste(Probe_ID,'I', sep='_')) %>%
      dplyr::rename(Assay_Design_Id=Probe_ID, AlleleA_Probe_Sequence=PRB1_U, AlleleB_Probe_Sequence=PRB1_M) %>% 
      dplyr::mutate(AlleleA_Probe_Id=paste(Assay_Design_Id,'A',sep='_'),
                    AlleleB_Probe_Id=case_when(!is.na(AlleleB_Probe_Sequence) ~ paste(Assay_Design_Id,'B',sep='_'),
                                               TRUE ~ NA_character_),
                    # AlleleB_Probe_Sequence=AlleleB_Probe_Sequence,
                    Normalization_Bin=case_when(is.na(AlleleB_Probe_Sequence) ~ 'C',
                                                NXB_M=='A' | NXB_M=='T' ~ 'A',
                                                NXB_M=='a' | NXB_M=='t' ~ 'A',
                                                NXB_M=='C' | NXB_M=='G' ~ 'B',
                                                NXB_M=='c' | NXB_M=='g' ~ 'B',
                                                TRUE ~ NA_character_)
      ) %>% 
      dplyr::select(Assay_Design_Id,AlleleA_Probe_Id,AlleleA_Probe_Sequence,AlleleB_Probe_Id,AlleleB_Probe_Sequence,Normalization_Bin) %>% 
      dplyr::mutate(Last_Base_A=stringr::str_sub(AlleleA_Probe_Sequence,50), Last_Base_B=stringr::str_sub(AlleleB_Probe_Sequence,50)) %>% 
      dplyr::filter(stringr::str_to_upper(Last_Base_A) != stringr::str_to_upper(Last_Base_B)) 
    
    # Fix Last Base Inconsistencies::
    #
    #  inf1_ord_tib %>% dplyr::mutate(Last_Base_A=stringr::str_sub(AlleleA_Probe_Sequence,50), Last_Base_B=stringr::str_sub(AlleleB_Probe_Sequence,50)) %>% dplyr::select(Assay_Design_Id,Last_Base_A, Last_Base_B) %>% group_by(Last_Base_B) %>% dplyr::summarise(Count=n())
    
    inf1_ord_tib %>% dplyr::select(1) %>% dplyr::distinct() %>% base::nrow()
    inf1_ord_tib %>% dplyr::select(2) %>% dplyr::distinct() %>% base::nrow()
    inf1_ord_tib %>% dplyr::select(3) %>% dplyr::distinct() %>% base::nrow()
    
    inf1_ord_tib %>% dplyr::select(4) %>% dplyr::distinct() %>% base::nrow()
    inf1_ord_tib %>% dplyr::select(5) %>% dplyr::distinct() %>% base::nrow()
    inf1_ord_tib %>% dplyr::select(6) %>% dplyr::distinct() %>% base::nrow()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Designs To Orders:: II
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    inf2_ord_tib <- prb_tibs %>% dplyr::filter(PRB1_M!=PRB1_U) %>% dplyr::bind_rows() %>%
      dplyr::mutate(Probe_ID=paste(Probe_ID,'II', sep='_')) %>%
      dplyr::rename(Assay_Design_Id=Probe_ID, AlleleA_Probe_Sequence=PRB2_D) %>% 
      dplyr::mutate(AlleleA_Probe_Id=paste(Assay_Design_Id,'A',sep='_'),
                    AlleleB_Probe_Id=NA,
                    AlleleB_Probe_Sequence=NA,
                    Normalization_Bin='C'
      ) %>% dplyr::select(Assay_Design_Id,AlleleA_Probe_Id,AlleleA_Probe_Sequence,AlleleB_Probe_Id,AlleleB_Probe_Sequence,Normalization_Bin)
    
    inf2_ord_tib %>% dplyr::select(1) %>% dplyr::distinct() %>% base::nrow()
    inf2_ord_tib %>% dplyr::select(2) %>% dplyr::distinct() %>% base::nrow()
    inf2_ord_tib %>% dplyr::select(3) %>% dplyr::distinct() %>% base::nrow()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Write Order File::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ord_tib <- dplyr::bind_rows(inf1_ord_tib, inf2_ord_tib) %>% dplyr::arrange(Assay_Design_Id)
    
    ord_tib %>% dplyr::select(1) %>% dplyr::distinct() %>% base::nrow()
    ord_tib %>% dplyr::select(2) %>% dplyr::distinct() %>% base::nrow()
    ord_tib %>% dplyr::select(3) %>% dplyr::distinct() %>% base::nrow()
    
    ord_tib %>% dplyr::select(4) %>% dplyr::distinct() %>% base::nrow()
    ord_tib %>% dplyr::select(5) %>% dplyr::distinct() %>% base::nrow()
    ord_tib %>% dplyr::select(Normalization_Bin) %>% dplyr::group_by(Normalization_Bin) %>% dplyr::summarise(Bin_Count=n())
    
    readr::write_csv(ord_tib, ord_snp_csv)
  }
}

# End of file
