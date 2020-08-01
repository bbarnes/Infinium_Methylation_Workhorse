
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
opt$improbeTSV <- NULL

# Chip Platform and Version Parameters::
opt$platform <- NULL
opt$version  <- NULL

# Probe Filter Parameters::
opt$minPrbScore <- 0.3
opt$minCpgRank  <- "3,2,1,0"
opt$minScrRank  <- "0.6,0.5,0.4,0.3"
opt$strandCO    <- 'C,O'
opt$pickBest    <- FALSE

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
  par$datDir  <- file.path(base::dirname(base::normalizePath(par$srcDir)), 'dat')
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  
  opt$platform  <- 'EPIC'
  opt$version   <- 'C0'

  opt$clean  <- TRUE
  opt$clean  <- FALSE
  
  opt$runName <- 'CpGs_UnivChicago_alldesigns_55860sites'
  opt$improbeTSV  <- '/Users/bbarnes/Documents/CustomerFacing/CustomContent/data/CpGs_UnivChicago_alldesigns_55860sites.tsv.gz'
  
  # fin_sum_csv <- file.path(opt$outDir, paste(opt$runName, 'order.csv', sep='.') )
  # fin_ord_csv <- file.path(opt$outDir, paste(opt$runName, 'order-summary.csv', sep='.') )
  # des_ord_tsv <- file.path(opt$outDir, 'data', paste(opt$runName,'ForwardSequence.tsv.gz', sep='.') )
  
  opt$outDir <- file.path(par$topDir, par$prgmTag, 'CustomContent')

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
    # make_option(c("--buildDirs"), type="character", default=opt$buildDir, 
    #             help="List of Build Directory [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--improbeTSV"), type="character", default=opt$improbeTSV, 
                help="Improbe Design Output (TSV) [default= %default]", metavar="character"),
    
    # Chip Platform and Version Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform name (HM50, EPIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest version (B2, B4, C0) [default= %default]", metavar="character"),
    
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

if (is.null(opt$outDir) || 
    is.null(opt$runName) || is.null(opt$improbeTSV) ||
    is.null(opt$platform) || is.null(opt$version) ||
    
    is.null(opt$minPrbScore) || is.null(opt$minCpgRank) || is.null(opt$minScrRank) ||
    is.null(opt$strandCO) || is.null(opt$pickBest) ||
    
    is.null(opt$percisionBeta) || is.null(opt$percisionPval) ||
    is.null(opt$execute) || is.null(opt$single) || is.null(opt$parallel) || is.null(opt$cluster) ||
    
    is.null(opt$clean) || is.null(opt$Rscript) || is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )
  
  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  if (is.null(opt$outDir))     cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$buildDir))   cat(glue::glue("[Usage]: buildDirs is NULL!!!{RET}"))
  if (is.null(opt$runName))    cat(glue::glue("[Usage]: runName is NULL!!!{RET}"))
  if (is.null(opt$improbeTSV)) cat(glue::glue("[Usage]: improbeTSV is NULL!!!{RET}"))
  
  if (is.null(opt$platform)) cat(glue::glue("[Usage]: platform is NULL!!!{RET}"))
  if (is.null(opt$version))  cat(glue::glue("[Usage]: version is NULL!!!{RET}"))

  if (is.null(opt$minPrbScore))  cat(glue::glue("[Usage]: minPrbScore is NULL!!!{RET}"))
  if (is.null(opt$minCpgRank))   cat(glue::glue("[Usage]: minCpgRank is NULL!!!{RET}"))
  if (is.null(opt$minScrRank))   cat(glue::glue("[Usage]: minScrRank is NULL!!!{RET}"))
  if (is.null(opt$strandCO))     cat(glue::glue("[Usage]: strandCO is NULL!!!{RET}"))
  if (is.null(opt$pickBest))     cat(glue::glue("[Usage]: pickBest is NULL!!!{RET}"))
  
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

strandCO_vec <- opt$strandCO %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
cpgRank_vec  <- opt$minCpgRank %>% str_split(pattern=',', simplify=TRUE) %>% as.integer() %>% as.vector()
scrRank_vec  <- opt$minScrRank %>% str_split(pattern=',', simplify=TRUE) %>% as.double() %>% as.vector()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         General Order CpGs Only::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ord_out_dir <- file.path(opt$outDir, 'orders')
if (!dir.exists(ord_out_dir)) dir.create(ord_out_dir, recursive=TRUE)

des_ord_tib <- suppressMessages(suppressWarnings( readr::read_tsv(opt$improbeTSV)) ) %>% 
  dplyr::mutate(Probe_Type='cg')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Filter and Assign Optimal Design Type::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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

des_fwd_tib <- des_scr_tib %>% dplyr::select(Seq_ID, Forward_Sequence, Probe_Type) %>%
  dplyr::distinct(Seq_ID, .keep_all=TRUE)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Design All Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

new_prb_tib <- tib2prbs(tib=des_fwd_tib, idsKey="Seq_ID", prbKey="Probe_Type", seqKey="Forward_Sequence", verbose=opt$verbose)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Match Probes by Strand to Top Picks::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sel_prb_tib <- dplyr::inner_join(new_prb_tib,des_scr_tib, by=c("Seq_ID","FR_Str","CO_Str") ) 
sel_prb_tib %>% dplyr::group_by(Design_Type) %>% dplyr::summarise(Count=n()) %>% print()

sel_prb1_tib <- sel_prb_tib %>% dplyr::filter(Design_Type=='I')
sel_prb2_tib <- sel_prb_tib %>% dplyr::filter(Design_Type=='II')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Format Order File::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

new_ord1_tib <- prbs2order(sel_prb1_tib, verbose=opt$verbose) %>% 
  dplyr::bind_rows() %>% 
  dplyr::filter(Valid_Design_Bool) %>% 
  dplyr::select(Assay_Design_Id:Normalization_Bin) %>% dplyr::filter(Normalization_Bin!='C')

new_ord2_tib <- prbs2order(sel_prb2_tib, verbose=opt$verbose) %>% 
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
  )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Write Final Order File::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

new_ord_csv <- file.path(opt$outDir, paste(opt$runName,'order.csv.gz', sep='.'))
readr::write_csv(new_ord_tib, new_ord_csv)

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
