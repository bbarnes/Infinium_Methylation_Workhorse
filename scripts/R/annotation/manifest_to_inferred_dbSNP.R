

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Load sesame:: This causes issues with "ExperimentHub Caching causes a warning"
suppressWarnings(suppressPackageStartupMessages( base::require("sesame") ))
suppressWarnings(suppressPackageStartupMessages( base::require("dbplyr") ))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("plyr")) )
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))

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
par$prgmDir <- 'annotation'
par$prgmTag <- 'manifest_to_inferred_dbSNP'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

par$retData     <- FALSE

# Executables::
opt$Rscript <- NULL

# Directories::
opt$outDir     <- NULL
opt$idatsDir   <- NULL

# Optional Files::
opt$subManifest  <- FALSE
opt$manifestPath <- 'auto'
opt$auto_sam_csv <- NULL

# Platform/Method Options::
opt$platform  <- NULL
opt$manifest  <- NULL

# Run Options::
opt$fresh       <- FALSE
opt$buildSubDir <- FALSE
opt$autoDetect  <- FALSE
opt$skipSwap    <- FALSE
opt$workflows   <- NULL

# Output Options::
opt$loadIdat    <- FALSE
opt$saveIdat    <- FALSE

opt$loadSsets   <- FALSE
opt$saveSsets   <- FALSE
opt$saveRawSset <- FALSE

opt$addSentrixID <- FALSE
opt$writeSset    <- FALSE
opt$writeSsum    <- FALSE
opt$writeCalls   <- FALSE
opt$writeSsheet  <- FALSE
opt$writeAuto    <- FALSE

opt$addRawCalls <- FALSE

# Threshold Options::
opt$minNegPval   <- 0.02
opt$minOobPval   <- 0.1
opt$minNegPerc   <- 98
opt$minOobPerc   <- 90
opt$minDeltaBeta <- 0.2

opt$percisionSigs <- 1
opt$percisionBeta <- 4
opt$percisionPval <- 6

# Parallel/Cluster Options::
opt$single   <- FALSE
opt$parallel <- FALSE
opt$cluster  <- FALSE

# Plotting Options::
opt$plotSset  <- FALSE
opt$plotCalls <- FALSE
opt$plotAuto  <- FALSE

opt$plotFormat <- 'pdf'
opt$plotFormat <- 'png'

opt$dpi <- 72
opt$dpi <- 120

opt$plotMax <- 10000
opt$plotSub <- 5000

opt$opt_csv  <- NULL
opt$par_csv  <- NULL
opt$time_csv <- NULL
opt$time_org_txt <- NULL

# verbose Options::
opt$verbose <- 3

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
  
  opt$outDir <- file.path(par$topDir, 'scratch')
  locIdatDir <- file.path(par$topDir, 'data/idats')
  
  opt$workflows <- 'ind'
  
  opt$buildSubDir  <- FALSE
  opt$autoDetect   <- FALSE
  opt$writeCalls   <- TRUE
  opt$writeSsheet  <- TRUE
  
  opt$platform   <- 'EPIC'
  opt$manifest   <- 'B4'
  opt$platform   <- NULL
  opt$manifest   <- NULL
  
  par$expRunStr  <- NULL
  par$expChipNum <- NULL
  
  par$local_runType <- 'CORE'
  par$local_runType <- 'EXCBR'
  par$local_runType <- 'qcMVP'
  par$local_runType <- 'GRCm38'
  par$local_runType <- 'COVID'
  par$local_runType <- 'COVIC'
  
  if (par$local_runType=='COVID') {
    par$expRunStr  <- 'COVID-Direct-Set1'
    par$expChipNum <- '204756130014'
    opt$autoDetect <- FALSE
    
    opt$workflows <- 'i,nd,ndi,ind'
    
  } else if (par$local_runType=='COVIC') {
    par$expRunStr  <- 'COVIC-Set1-15052020'
    par$expChipNum <- '204500250013'
    opt$autoDetect <- TRUE
  } else if (par$local_runType=='GRCm38') {
    par$expRunStr <- 'MURMETVEP_mm10_betaTest_06082020'
    par$expRunStr <- 'VanAndel_mm10_betaTest_31082020'
    par$expRunStr <- 'ILMN_mm10_betaTest_17082020'
    opt$autoDetect <- FALSE
  } else if (par$local_runType=='qcMVP') {
    par$expRunStr  <- 'CNTL-Samples_VendA_10092020'
    par$expRunStr  <- 'CNTL-Samples_VendA_10092020_test'
    opt$autoDetect <- TRUE
    opt$dpi <- 72
  } else if (par$local_runType=='CORE') {
    par$expRunStr  <- 'BETA-8x1-EPIC-Core'
    par$expChipNum <- '202761400007'
    
    par$expRunStr  <- 'ADRN-blood-nonAtopic_EPIC'
    par$expChipNum <- '201125090068'
    
    par$expRunStr  <- 'GSE122126_EPIC'
    par$expChipNum <- '202410280180'
    
    par$expRunStr  <- 'EPIC-BETA-8x1-CoreCancer'
    par$expChipNum <- '201502830033'
    opt$autoDetect <- TRUE
  } else if (par$local_runType=='EXCBR') {
    par$expRunStr  <- 'Excalibur-Old-1609202'
    par$expChipNum <- '204076530053'
    par$expChipNum <- '204076530110'
    
    par$expRunStr  <- 'Excalibur-New-1609202'
    par$expChipNum <- '202915460071'
    opt$autoDetect <- TRUE
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unrecognized local_runType={par$local_runType}.{RET}{RET}"))
  }
  
  opt$idatsDir <- file.path(locIdatDir, paste('idats',par$expRunStr, sep='_') )
  if (!is.null(par$expChipNum)) opt$idatsDir <- file.path(locIdatDir, paste('idats',par$expRunStr, sep='_'),  par$expChipNum)
  opt$auto_sam_csv <- file.path(par$datDir, 'ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo_Mean-COVIC-280-NP-ind_negs-0.02.csv.gz')
  
  opt$outDir <- file.path(par$topDir, 'scratch', par$prgmTag, par$expRunStr, 'manifest-C1')
  
  par$retData  <- TRUE
  opt$single   <- TRUE
  opt$parallel <- FALSE
  opt$cluster  <- FALSE
  opt$verbose  <- 3
  
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
    make_option(c("--idatsDir"), type="character", default=opt$idatsDir, 
                help="idats directory [default= %default]", metavar="character"),
    
    # Optional Files::
    make_option(c("--manifestPath"), type="character", default=opt$manifestPath,
                help="Path to manfifest (CSV) otherwise use dat [default= %default]", metavar="character"),
    # make_option(c("--addressPath"), type="character", default=opt$addressPath,
    #             help="Path to address (RDS) otherwise use dat [default= %default]", metavar="character"),
    make_option(c("--subManifest"), action="store_true", default=opt$subManifest,
                help="Boolean variable to use subset manifest instead of subset. [default= %default]", metavar="boolean"),
    make_option(c("--auto_sam_csv"), type="character", default=opt$auto_sam_csv,
                help="Path to auto detect beta values (CSV) [default= %default]", metavar="character"),
    
    # Platform/Method Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Forced platform [EPIC, 450k, 27k, NZT] otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--manifest"), type="character", default=opt$manifest, 
                help="Forced manifest [B1, B2, B4] otherwise auto-detect [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--fresh"), action="store_true", default=opt$fresh,
                help="Boolean variable to build fresh version of database files [default= %default]", metavar="boolean"),
    make_option(c("--buildSubDir"), action="store_true", default=opt$buildSubDir,
                help="Boolean variable to build subdirectories based on Chip/BeadPool (for R&D purposes) [default= %default]", metavar="boolean"),
    make_option(c("--autoDetect"), action="store_true", default=opt$autoDetect,
                help="Boolean variable to auto detect reference samples. Must provide reference samples. [default= %default]", metavar="boolean"),
    make_option(c("--skipSwap"), action="store_true", default=opt$skipSwap,
                help="Boolean variable to skpping appending swap percentages to sample sheet. [default= %default]", metavar="boolean"),
    make_option(c("--workflows"), type="character", default=opt$workflows,
                help="Order of operations comma seperated [ raw,ind,ndi,din ] [default= %default]", metavar="character"),
    # make_option(c("--sampleSheet"), type="character", default=opt$sampleSheet, 
    #             help="Target Sample Sheet containing samples/chips to ONLY analyze [default= %default]", metavar="character"),
    
    # Output Options::
    make_option(c("--loadIdat"), action="store_true", default=opt$loadIdat,
                help="Boolean variable to load existing IDAT from RDS file [default= %default]", metavar="boolean"),
    make_option(c("--saveIdat"), action="store_true", default=opt$saveIdat,
                help="Boolean variable to write IDAT RDS file [default= %default]", metavar="boolean"),
    
    make_option(c("--loadSsets"), action="store_true", default=opt$loadSsets,
                help="Boolean variable to load existing Signal Set from RDS file [default= %default]", metavar="boolean"),
    make_option(c("--saveSsets"), action="store_true", default=opt$saveSsets,
                help="Boolean variable to write Signal Set RDS file [default= %default]", metavar="boolean"),
    make_option(c("--saveRawSset"), action="store_true", default=opt$saveRawSset,
                help="Boolean variable to write Raw Signal Set RDS file [default= %default]", metavar="boolean"),
    
    make_option(c("--addSentrixID"), action="store_true", default=opt$addSentrixID,
                help="Boolean variable to add Sentrix Name to calls output columns [default= %default]", metavar="boolean"),
    make_option(c("--writeSset"), action="store_true", default=opt$writeSset,
                help="Boolean variable to write Signal Set file [default= %default]", metavar="boolean"),
    make_option(c("--writeSsum"), action="store_true", default=opt$writeSsum,
                help="Boolean variable to write Signal Set Summary file [default= %default]", metavar="boolean"),
    make_option(c("--writeCalls"), action="store_true", default=opt$writeCalls,
                help="Boolean variable to write Calls (Pval/Beta) file [default= %default]", metavar="boolean"),
    make_option(c("--writeSsheet"), action="store_true", default=opt$writeSsheet,
                help="Boolean variable to Sample Sheet file [default= %default]", metavar="boolean"),
    make_option(c("--writeAuto"), action="store_true", default=opt$writeAuto,
                help="Boolean variable to write Auto-Detection Matricies (Pval/Beta) file [default= %default]", metavar="boolean"),
    make_option(c("--addRawCalls"), action="store_true", default=opt$addRawCalls,
                help="Boolean variable to output raw calls [default= %default]", metavar="boolean"),
    
    # Threshold Options::
    make_option(c("--minNegPval"), type="double", default=opt$minNegPval, 
                help="Minimum passing detection p-value using Negative Controls [default= %default]", metavar="double"),
    make_option(c("--minOobPval"), type="double", default=opt$minOobPval,
                help="Minimum passing detection p-value using Negative Out-Of-Band [default= %default]", metavar="double"),
    
    make_option(c("--minNegPerc"), type="double", default=opt$minNegPerc, 
                help="Minimum percentage of loci passing detection p-value using Negative Controls to flag Requeue of sample. [default= %default]", metavar="double"),
    make_option(c("--minOobPerc"), type="double", default=opt$minOobPerc, 
                help="Minimum percentage of loci passing detection p-value using Out-Of-Band to flag Requeue of sample. [default= %default]", metavar="double"),
    
    make_option(c("--minDeltaBeta"), type="double", default=opt$minDeltaBeta,
                help="Minimum passing delta-beta. Used in AutoSampleSheet cacluclations [default= %default]", metavar="double"),
    
    make_option(c("--percisionSigs"), type="integer", default=opt$percisionSigs,
                help="Rounding percision for signal values in calls output files [default= %default]", metavar="double"),
    make_option(c("--percisionBeta"), type="integer", default=opt$percisionBeta,
                help="Rounding percision for beta values in calls output files [default= %default]", metavar="double"),
    make_option(c("--percisionPval"), type="integer", default=opt$percisionPval,
                help="Rounding percision for detection p-values in calls output files [default= %default]", metavar="double"),
    
    # Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
    # Plotting Options::
    make_option(c("--plotSset"), action="store_true", default=opt$plotSset,
                help="Boolean variable to plot intensity distributions for sset [default= %default]", metavar="boolean"),
    make_option(c("--plotCalls"), action="store_true", default=opt$plotCalls,
                help="Boolean variable to plot detection p-values and beta distributions [default= %default]", metavar="boolean"),
    make_option(c("--plotAuto"), action="store_true", default=opt$plotAuto,
                help="Boolean variable to plot Auto-Detection Matricies (Pval/Beta) file [default= %default]", metavar="boolean"),
    
    make_option(c("--plotFormat"), type="character", default=opt$plotFormat, 
                help="Plotting output format [default= %default]", metavar="character"),
    make_option(c("--dpi"), type="double", default=opt$dpi, 
                help="DPI for plot images Plotting [default= %default]", metavar="double"),
    
    make_option(c("--plotMax"), type="double", default=opt$plotMax, 
                help="Max Sample Display Count for Plotting [default= %default]", metavar="double"),
    make_option(c("--plotSub"), type="double", default=opt$plotSub, 
                help="Sub Sample Display Count for Plotting [default= %default]", metavar="double"),
    
    make_option(c("--opt_csv"), type="character", default=opt$opt_csv, 
                help="Unused variable opt_csv [default= %default]", metavar="character"),
    make_option(c("--par_csv"), type="character", default=opt$par_csv, 
                help="Unused variable par_csv [default= %default]", metavar="character"),
    make_option(c("--time_csv"), type="character", default=opt$time_csv, 
                help="Unused variable time_csv [default= %default]", metavar="character"),
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
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )

opt <- program_init(name=par$prgmTag,
                    opts=opt, opt_reqs=opt_reqs, 
                    pars=par, par_reqs=par_reqs,
                    libs=TRUE,rcpp=FALSE,
                    verbose=opt$verbose,vt=3,tc=0,tt=NULL)

par_tib <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")

if (TRUE) {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #            EPIC hg19 Infinium I Inferred SNP Table Construction::
  #
  #               Should move this scratch to somewhere else::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # 450k Genome Studio::
  #
  # hm45_man_csv <- '/Users/bretbarnes/Documents/data/manifests/HumanMethylation450_15017482_v.1.2.csv.gz'
  # hm45_man_lst <- loadManifestGenomeStudio(file=hm45_man_csv, verbose=opt$verbose+4)
  # hm45_man_tib <- hm45_man_lst$man %>% dplyr::select(IlmnID,Genome_Build,CHR,MAPINFO,Strand)
  # hm45_ctl_tib <- hm45_man_lst$ctl %>% dplyr::mutate(Address=as.double(Address))
  
  #
  # SNP Table
  #
  epicI_ses_grs <- sesameDataPullVariantAnno_InfiniumI(platform = "EPIC")
  # epicR_ses_grs <- sesameDataPullVariantAnno_SNP()
  
  epicI_ses_tib <- epicI_ses_grs %>% 
    as.data.frame() %>% tibble::rownames_to_column(var="IlmnID") %>% tibble::as_tibble()
  
  # Five Examples by hand::
  #  gs_ses_csv <- '/Users/bretbarnes/Documents/tools/notes/inferred-snps.GS-vs-Ses.csv'
  #  gs_ses_tib <- readr::read_csv(gs_ses_csv)
  #
  #
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                              Preprocessing::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # epic_man_csv <- file.path(par$datDir, 'manifest/base/EPIC-B4.manifest.sesame-base.cpg-sorted.csv.gz')
  epic_man_csv <- '/Users/bretbarnes/Documents/data/manifests/MethylationEPIC_v-1-0_B4.csv.gz'
  epic_man_lst <- loadManifestGenomeStudio(file=epic_man_csv, verbose=opt$verbose+4)
  epic_ctl_tib <- epic_man_lst$ctl %>% dplyr::mutate(Address=as.double(Address))
  epic_man_tib <- epic_man_lst$man %>% 
    dplyr::select(IlmnID,Genome_Build,CHR,MAPINFO,Strand,Next_Base,Infinium_Design_Type) %>%
    dplyr::filter(Infinium_Design_Type=='I') %>% 
    dplyr::mutate(
      POS=dplyr::case_when(
        Strand=='R' ~ MAPINFO+2,
        Strand=='F' ~ MAPINFO-1,
        TRUE ~ NA_real_
      )
    ) %>% dplyr::rename(Chrom=CHR) %>%
    dplyr::arrange(Chrom,MAPINFO)
  
  #
  # Validation Shown Below::
  #
  join_man_tib <- epic_man_tib %>% 
    dplyr::inner_join(epicI_ses_tib, by="IlmnID", suffix=c("_Man", "_Ses")) %>%
    dplyr::mutate(REF_CON=dplyr::case_when(
      strand=='+' & REF=='G' ~ 'A',
      strand=='-' & REF=='C' ~ 'T',
      TRUE ~ REF)
    ) %>% dplyr::mutate(seqnames=stringr::str_remove(seqnames, 'chr'))
  
  join_sum_tib <- join_man_tib %>% dplyr::filter(Chrom==seqnames) %>% 
    dplyr::mutate(Pos_Dif=MAPINFO-start) %>% 
    dplyr::group_by(Pos_Dif,Strand,strand,Next_Base,REF_CON) %>% 
    dplyr::summarise(Count=n(), .groups='drop')
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                                  Main::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  db151_hs37_col <- cols(CHR = col_character(),
                         POS = col_double(),
                         SNP = col_character(),
                         REF = col_character(),
                         ALT = col_character())
  
  db151_hs37_pat <- 'All_20180423.snp-5.chr.*Stsv.gz'
  db151_hs37_dir <- '/Users/bretbarnes/Documents/data/annotation/GRCh37'
  db151_hs37_fns <- list.files(path = db151_hs37_dir, pattern = db151_hs37_pat, full.names = TRUE)
  
  all_int_tib <- NULL
  for (file in db151_hs37_fns) {
    cat(glue::glue("Loading: file={file}{RET}"))
    db151_hs37_tib <- readr::read_tsv(file, col_names=names(db151_hs37_col$cols), col_types=db151_hs37_col)
    
    cur_int_tib <- db151_hs37_tib %>% 
      dplyr::inner_join(epic_man_tib, by="POS") %>%
      dplyr::filter(CHR==Chrom)
    
    # Validation::
    cur_val_tib <- cur_int_tib %>% 
      dplyr::inner_join(epicI_ses_tib, by="IlmnID", suffix=c("_Man", "_Ses"))
    
    cur_val_cnt <- cur_val_tib %>% base::nrow()
    cur_mat_cnt <- cur_val_tib %>% dplyr::filter(SNP==rs) %>% base::nrow()
    
    all_int_tib <- all_int_tib %>% dplyr::bind_rows(cur_int_tib)
    all_int_cnt <- all_int_tib %>% base::nrow()
    
    cat(glue::glue("Validation: {cur_val_cnt} vs {cur_mat_cnt}; total={all_int_cnt}{RET}"))
  }
  
  out_dat <- all_int_tib
  out_csv <- '/Users/bretbarnes/Documents/data/annotation/GRCh37/current-inferred-dbSNP151.hg19.clean.csv.gz'
  readr::write_csv(out_dat,out_csv)
  
  # all_int_tib %>% dplyr::mutate(REF_Base=dplyr::case_when())
  
}

# End of file
