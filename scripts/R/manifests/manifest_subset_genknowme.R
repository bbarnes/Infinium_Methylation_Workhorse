
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Load Core Packages::
suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("sesame",quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("minfi",quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse",quietly=TRUE) ))

# Load Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel",quietly=TRUE) ))

# Load Performance Packages
suppressWarnings(suppressPackageStartupMessages( base::require("profmem",quietly=TRUE) ))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Global Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ses_mm10_man_tib <- sesameData::sesameDataGet("MM285.mm10.manifest") %>% 
#   as.data.frame() %>% 
#   tibble::as_tibble()
# 
# ses_mm10_man_csv <- "/Users/bretbarnes/Documents/tools/Infinium_Methylation_Workhorse/dat/manifest/sesame_src/MM285.mm10.manifest.csv.gz"
# readr::write_csv(ses_mm10_man_tib, ses_mm10_man_csv)

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM <- ","
TAB <- "\t"
RET <- "\n"
BNG <- "|"

# mm10_man_grs = sesameDataGet('MM285.mm10.manifest')
# mm10_man_tib <- mm10_man_grs %>% as.data.frame() %>% rownames_to_column(var="Probe_ID") %>% tibble::as_tibble()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Define Default Params and Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
opt <- NULL
par <- NULL
run <- NULL

# Illumina based directories::
par$date    <- Sys.Date() %>% as.character()
par$runMode <- ''
par$maxTest <- NULL
par$macDir1 <- '/Users/bbarnes/Documents/Projects/methylation'
par$macDir2 <- '/Users/bretbarnes/Documents'
par$lixDir1 <- '/illumina/scratch/darkmatter'
par$lixDir  <- '/illumina/scratch/darkmatter'
par$retData <- FALSE
par$manDir  <- NULL

# Program Parameters::
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'manifests'
par$prgmTag <- 'manifest_subset_genknowme'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# Executables::
opt$Rscript <- NULL

# Directories::
opt$outDir <- NULL
opt$datDir <- NULL

# Platform/Method Options::
opt$manifest <- NULL
opt$platform <- NULL
opt$version  <- NULL

# Run Options::
opt$fresh       <- FALSE
opt$workflow    <- NULL
opt$manDirPath  <- NULL
opt$manDirName  <- 'core'
opt$subManifest <- NULL
opt$man_suffix  <- ".manifest.sesame-base.cpg-sorted.csv.gz"

# Threshold Options::
opt$rand <- NULL
opt$per1 <- NULL
opt$per2 <- NULL
opt$type <- NULL

# Parallel/Cluster Options::
opt$single   <- FALSE
opt$parallel <- FALSE
opt$cluster  <- FALSE

# Time Tracking::
opt$time_org_txt <- NULL
opt$trackTime    <- FALSE

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
  
  opt$manDirPath <- NULL
  opt$manDirName <- 'base'
  opt$manDirName <- 'core'

  par$local_runType <- 'EPIC-noob-BP4'
  opt$runName  <- par$local_runType
  
  if (FALSE) {
  } else if (par$local_runType=='EPIC-noob-BP4') {
    
    # For sub manifest testing::
    opt$platform   <- "EPIC"
    opt$version    <- "B4"

    # opt$per2 <- "0,5,10,15,20,30,60,90,100"

    opt$rand <- 5
    opt$per1 <- "0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100"
    opt$per2 <- "0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100"
    opt$type <- "cg"
    
    opt$manifest <- 
      file.path(par$topDir, "data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B4-Beadpool_ID.csv.gz")
    
    par$ses_man_dir <- file.path(par$topDir, "data/manifests/methylation/Sesame")
    opt$subManifest <- file.path(par$ses_man_dir, "EPIC-B4-BP4.manifest.sesame-base.cpg-sorted.csv.gz")
    
    opt$verbose <- 40
    
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unrecognized local_runType={par$local_runType}.{RET}{RET}"))
  }
  
  opt$outDir <- file.path(par$topDir, 'scratch',"test",par$runMode)
  opt$outDir <- file.path(par$topDir, 'scratch',par$runMode)
  
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
    make_option(c("-d","--datDir"), type="character", default=opt$datDir, 
                help="List of idats directory(s), commas seperated [default= %default]", metavar="character"),
    
    # Platform/Method Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--manifest"), type="character", default=opt$manifest,
                help="Path to manfifest (CSV) otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Forced platform [EPIC, 450k, 27k, NZT] otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Forced version [B1, B2, B4, etc.] otherwise auto-detect [default= %default]", metavar="character"),
    
    # Optional Files::
    make_option(c("--subManifest"), type="character", default=opt$subManifest,
                help="Path to manifest to subset (CSV) [default= %default]", metavar="character"),

    # Run Parameters::
    make_option(c("--fresh"), action="store_true", default=opt$fresh,
                help="Boolean variable to build fresh version of database files [default= %default]", metavar="boolean"),

    make_option(c("--manDirPath"), type="character", default=opt$manDirPath,
                help="Manifest directory path otherwise use default dat/manifest directory [default= %default]", metavar="character"),
    make_option(c("--manDirName"), type="character", default=opt$manDirName,
                help="Manifest directory name [default= %default]", metavar="character"),
    make_option(c("--man_suffix"), type="character", default=opt$man_suffix,
                help="Manifest suffix search name. Default is set to predefined Sesame suffix. [default= %default]", metavar="character"),
    
    #
    # Threshold Options::
    #
    make_option(c("--rand"), type="integer", default=opt$rand, 
                help="Number of randomization rounds  [default= %default]", metavar="integer"),
    make_option(c("--per1"), type="character", default=opt$per1, 
                help="Percent subsets Infinuum I (commas seperated list)  [default= %default]", metavar="character"),
    make_option(c("--per2"), type="character", default=opt$per2, 
                help="Percent subsets Infinuum II (commas seperated list)  [default= %default]", metavar="character"),
    make_option(c("--type"), type="character", default=opt$type, 
                help="Probe type(s) to subset (commas seperated list)  [default= %default]", metavar="character"),
    
    # Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
    make_option(c("--time_org_txt"), type="character", default=opt$time_org_txt, 
                help="Unused variable time_org_txt [default= %default]", metavar="character"),
    make_option(c("--trackTime"), action="store_true", default=opt$trackTime,
                help="Boolean variable tack run times [default= %default]", metavar="boolean"),
    
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

par_reqs <- c('runMode','prgmTag','scrDir','exePath')
opt_reqs <- c('outDir','Rscript','verbose',
              'rand','per1','per2','type')

par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) 
  stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))

for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', 
                         full.names=TRUE, recursive=TRUE)) base::source(sfile)
if (opt$verbose>0)
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

if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
cat(glue::glue("[{par$prgmTag}]: Output Directory (TOP)={opt$outDir}...{RET}"))

run$manDir <- file.path(opt$outDir, 'man')
if (!dir.exists(run$manDir)) dir.create(run$manDir, recursive=TRUE)

cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Suboptimal Selection for 20k:: Open DMAPs
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  
  maxG_cnt <- 5000
  maxR_cnt <- 5000
  max2_cnt <- 10000
  
  dmap_cols <- cols(
    tango    = col_integer(),
    color    = col_integer()
  )
  
  dmap_tan_tsv <- "/Users/bretbarnes/Documents/data/CustomContent/Genknowme/dmaps/15073382_A_MRGPOOL,_EPIC_15054510_UsedCodes.txt"
  dmap_tan_tib <- suppressMessages(suppressWarnings(
    readr::read_tsv(dmap_tan_tsv, col_names=names(dmap_cols$cols), col_types=dmap_cols) ))

  ses_man_dir <- file.path(par$topDir, "data/manifests/methylation/Sesame")
  sel_aqp_dir <- file.path(par$topDir, "data/CustomContent/McMaster/McMaster10Kselection/AQP.v1")
  
  sel_csv <- file.path(ses_man_dir, "EPIC-B4-BP4.manifest.sesame-base.cpg-sorted.csv.gz")
  pre_csv <- file.path(sel_aqp_dir, "Rand3-S0.060.manifest.sesame-base.cpg-sorted.csv.gz")

  # Remove Non-CpGs
  sel_tib <- suppressMessages(suppressWarnings(readr::read_csv(sel_csv) )) %>%
    dplyr::filter(stringr::str_starts(Probe_ID, "cg"))
  pre_tib <- suppressMessages(suppressWarnings(readr::read_csv(pre_csv) )) %>%
    dplyr::filter(stringr::str_starts(Probe_ID, "cg"))
  
  # Remove pre's in sel
  #
  sel_tib <- sel_tib %>% dplyr::filter(! Probe_ID %in% pre_tib$Probe_ID)
  
  
  # Split both sel and pre by G/R/2
  #
  pre_list <- pre_tib %>% dplyr::mutate(Prb_Des=dplyr::case_when(
    is.na(col) ~ "2", TRUE ~ col)
  ) %>% split(.$Prb_Des)
  sel_list <- sel_tib %>% dplyr::mutate(Prb_Des=dplyr::case_when(
    is.na(col) ~ "2", TRUE ~ col)
  ) %>% split(.$Prb_Des)
  
  
  # Tally up all of pre: preG_cnt, preR_cnt, pre2_cnt
  #
  pre2_cnt <- pre_list[["2"]] %>% base::nrow()
  preG_cnt <- pre_list[["G"]] %>% base::nrow()
  preR_cnt <- pre_list[["R"]] %>% base::nrow()
  
  # Pick the the top X type from sel to satisfy
  #
  #  maxG_cnt - preG_cnt = 0
  #  maxR_cnt - preR_cnt = 0
  #  max2_cnt - pre2_cnt = 0
  #
  mis2_cnt <- max2_cnt - pre2_cnt
  misR_cnt <- maxR_cnt - preR_cnt
  misG_cnt <- maxG_cnt - preG_cnt

  new2_tib <- sel_list[["2"]] %>% head(n=mis2_cnt)
  newG_tib <- sel_list[["G"]] %>% head(n=misG_cnt)
  newR_tib <- sel_list[["R"]] %>% head(n=misR_cnt)
  
  fin_sel_tib  <- dplyr::bind_rows(pre_tib,new2_tib,newG_tib,newR_tib) %>%
    clean_tibble()
  fin_sel_list <- fin_sel_tib %>% dplyr::mutate(Prb_Des=dplyr::case_when(
    is.na(col) ~ "2", TRUE ~ col)
  ) %>% split(.$Prb_Des)
  
  fin_tan_tib <- dplyr::bind_rows(
    fin_sel_tib %>% dplyr::filter(!is.na(M)) %>% dplyr::select(M) %>% dplyr::rename(tango=M),
    fin_sel_tib %>% dplyr::filter(!is.na(U)) %>% dplyr::select(U) %>% dplyr::rename(tango=U)
  ) %>% dplyr::arrange(tango)
  
  dmap_out_tib <- dmap_tan_tib %>% dplyr::filter(tango %in% fin_tan_tib$tango)
  
  dmap_out_dir <- "/Users/bretbarnes/Documents/data/CustomContent/Genknowme/dmaps/open-dmaps20k"
  dmap_man_csv <- file.path(dmap_out_dir, "open-dmaps20k.manifest.csv.gz")
  dmap_out_tsv <- file.path(dmap_out_dir, "open-dmaps20k.unlock.tsv.gz")
  
  readr::write_csv(fin_sel_tib,dmap_man_csv)
  readr::write_csv(dmap_out_tib,dmap_out_tsv)
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                  Main::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (!is.null(opt$subManifest) && file.exists(opt$subManifest)) {
  
  epic_man_tib <- 
    loadManifestGenomeStudio(opt$manifest, addSource = TRUE, normalize = TRUE, 
                             retType = "man", verbose = opt$verbose) %>% 
    clean_tibble()
  
  h450_man_csv <- file.path(par$topDir, "data/manifests/methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz")
  h450_man_tib <- loadManifestGenomeStudio(h450_man_csv, addSource = TRUE, normalize = TRUE, 
                                           retType = "man", verbose = opt$verbose) %>% 
    clean_tibble()
  
  # Customer Manifest:
  gkm_man_csv <- file.path(par$topDir, "data/CustomContent/Genknowme/manifest/GKME_EPIC_Overlap_Sesame.csv")
  gkm_man_tib <- suppressMessages(
    suppressWarnings(readr::read_csv(gkm_man_csv))) %>% clean_tibble()
  
  # 80 To be added from BP4::
  e80_man_csv <- file.path(par$topDir, "data/CustomContent/Genknowme/dat/add-on-list-n80.cpg-sorted.csv.gz")
  e80_man_tib <- suppressMessages(suppressWarnings(
    readr::read_csv(e80_man_csv, col_names=c("Probe_ID")) ))
  
  # 11k Previously selected from BP4::
  k11_man_csv <- file.path(par$topDir, "data/CustomContent/McMaster/McMaster10Kselection/Rand3-S0.060.manifest.sesame-base.cpg-sorted.csv.gz")
  k11_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(k11_man_csv) ))

  # Now we're missing extra BP4 to meet a total of 10k Inf2 and 5k Inf1_Grn/Red  
  bp4_prb_csv <- file.path(par$topDir, "data/CustomContent/Genknowme/manifest/bp4/bp4.rand3.cpg-sorted.csv")
  bp4_prb_tib <- suppressMessages(suppressWarnings(
    readr::read_csv(bp4_prb_csv, col_names=c("Probe_ID")) ))
  bp4_man_tib <- epic_man_tib %>% 
    dplyr::filter(IlmnID %in% bp4_prb_tib$Probe_ID)

  # bp4_man_tib %>% dplyr::filter(IlmnID %in% k11_man_tib$Probe_ID)
  # bp4_man_tib %>% dplyr::filter(IlmnID %in% gkm_man_tib$Probe_ID)
  
  # This should be zero::
  k11_sel_mis_cnt <- k11_man_tib %>% 
    dplyr::filter(!stringr::str_starts(Probe_ID, "ctl")) %>%
    dplyr::anti_join(bp4_man_tib, by=c("Probe_ID"="IlmnID")) %>%
    base::nrow()
  cat(glue::glue("[{par$prgmTag}]: k11_sel_mis_cnt={k11_sel_mis_cnt}; Exp==0.{RET}"))
  
  # This should be the mismatch of BP4 vs. e80 probes:: non-zero
  e80_bp4_mis_cnt <- e80_man_tib %>% 
    dplyr::filter(!stringr::str_starts(Probe_ID, "ctl")) %>%
    dplyr::anti_join(bp4_man_tib, by=c("Probe_ID"="IlmnID")) %>%
    base::nrow()
  cat(glue::glue("[{par$prgmTag}]: e80_bp4_mis_cnt={e80_bp4_mis_cnt}; Exp==0.{RET}"))
  
  # These should be the number of probes not found in BP4:: non zero
  gkm_sel_mis_cnt <- gkm_man_tib %>% 
    dplyr::filter(!stringr::str_starts(Probe_ID, "ctl")) %>%
    dplyr::anti_join(bp4_man_tib, by=c("Probe_ID"="IlmnID")) %>%
    base::nrow()
  cat(glue::glue("[{par$prgmTag}]: gkm_sel_mis_cnt={gkm_sel_mis_cnt}; Exp!=0.{RET}"))

  # This should be the mismatch of EPIC probes:: should be zero
  gkm_epi_mis_cnt <- gkm_man_tib %>% 
    dplyr::filter(!stringr::str_starts(Probe_ID, "ctl")) %>%
    dplyr::anti_join(epic_man_tib, by=c("Probe_ID"="IlmnID")) %>%
    base::nrow()
  cat(glue::glue("[{par$prgmTag}]: gkm_epi_mis_cnt={gkm_epi_mis_cnt}; Exp==0.{RET}"))
  
  # This should be the mismatch of 450k probes:: non-zero
  gkm_450_mis_cnt <- gkm_man_tib %>% 
    dplyr::filter(!stringr::str_starts(Probe_ID, "ctl")) %>%
    dplyr::anti_join(h450_man_tib, by=c("Probe_ID"="IlmnID")) %>%
    base::nrow()
  cat(glue::glue("[{par$prgmTag}]: gkm_450_mis_cnt={gkm_450_mis_cnt}; Exp!=0.{RET}"))

  # Make sure BP4 does not contain tangos found in gkm_man_tib
  gkm_450_misA_cnt <- gkm_man_tib %>% 
    dplyr::filter(!is.na(U)) %>%
    dplyr::inner_join(h450_man_tib, by=c("U"="AddressA_ID")) %>%
    base::nrow()
  cat(glue::glue("[{par$prgmTag}]: gkm_450_misA_cnt={gkm_450_misA_cnt}; Exp==0.{RET}"))
  
  # Make sure BP4 does not contain tangos found in gkm_man_tib
  gkm_450_misB_cnt <- gkm_man_tib %>% 
    dplyr::filter(!is.na(M)) %>%
    dplyr::inner_join(h450_man_tib, by=c("U"="AddressB_ID")) %>%
    base::nrow()
  cat(glue::glue("[{par$prgmTag}]: gkm_450_misB_cnt={gkm_450_misB_cnt}; Exp==0.{RET}"))
  
  #
  # TBD:: Combine all requested tangos
  #   - Add non k11 bp4 probes
  #
}

if (FALSE) {

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Load Manifest(s)::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  cat(glue::glue("[{par$prgmTag}]: Launching Samples in Linear Mode! ",
                 "isSingle={opt$single}"),"\n", sep='')
  
  pTracker <- NULL
  pTracker <- timeTracker$new()
  
  tar_man_dat <- base::list()
  
  par$manDir <- NULL
  if (!is.null(opt$manDirPath) && length(opt$manDirPath)>0 && dir.exists(opt$manDirPath)) {
    par$manDir <- opt$manDirPath
  } else if (!is.null(opt$manDirName)) {
    par$manDir <- file.path(par$datDir, 'manifest',opt$manDirName)
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Both manDirPath and manDirName are NULL!!!{RET}{RET}"))
  }
  if (opt$verbose>0)
    cat(glue::glue("[{par$prgmTag}]: Set manifest search directory={par$manDir}.{RET}"))
  
  tar_man_tib <- NULL
  tar_man_tib <- get_manifest_list(
    file=NULL, dir=par$manDir,
    platform=opt$platform, version=opt$version,
    suffix=opt$man_suffix,
    verbose=opt$verbose, tt=pTracker)
  
  tar_man_dat <- NULL
  tar_man_dat <- load_manifest_list(
    tib = tar_man_tib, field="path",
    verbose=opt$verbose, tt=pTracker)
  
  if (!is.null(opt$manifest)) {
    
    opt$man_key <- "Unknown"
    if (!is.null(opt$platform) && !is.null(opt$version))
      opt$man_key <- paste(opt$platform,opt$version)
    
    # tar_man_dat[[opt$man_key]] <-
    tmp_man_dat <- 
      loadManifestGenomeStudio(
        file=opt$manifest, addSource=TRUE, # normalize=TRUE, 
        verbose=opt$verbose, tt=pTracker)
    
    tmp_man_dat$man %>% 
      dplyr::group_by(Beadpool_ID, Infinium_Design_Type) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    
    tmp_man_dat$man %>% 
      dplyr::group_by(Probe_Type, Infinium_Design_Type) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    
    # tar_man_dat$`EPIC-B4` %>% dplyr::anti_join(tmp_man_dat$man, by=c(Probe_ID="IlmnID")) %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n())
    # tar_man_dat$`EPIC-B4` %>% dplyr::anti_join(tmp_man_dat$man, by=c("U"="AddressA_ID")) %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n())
    # tmp_man_dat$man %>% dplyr::anti_join(tar_man_dat$`EPIC-B4`, by=c("AddressA_ID"="U")) %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n())
    
    ses_bead_tib <- tar_man_dat$`EPIC-B4` %>% 
      dplyr::inner_join(dplyr::select(tmp_man_dat$man, AddressA_ID, Beadpool_ID),
                        by=c("U"="AddressA_ID"))
    ses_bead_sum <- ses_bead_tib %>%
      dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n())
    ses_bead_sum %>% print(n=base::nrow(ses_bead_sum))
    
    ses_ctl_tib <- tar_man_dat$`EPIC-B4` %>% 
      dplyr::filter(stringr::str_starts(Probe_ID, "ct")) %>% 
      dplyr::mutate(Beadpool_ID="CTL")
    
    ses_bead_dat <- ses_bead_tib %>% split(.$Beadpool_ID)
    
    bead_pools <- names(ses_bead_dat)
    for (bp in bead_pools) {
      if (opt$verbose>0)
        cat(glue::glue("[{par$prgmTag}]:{TAB} Current Bead Pool={bp}.{RET}"))
      
      cur_ses_tib <- dplyr::bind_rows(ses_bead_dat[[bp]], ses_ctl_tib) %>%
        clean_tibble()
      
      cur_out_fns <- paste(opt$platform, opt$version, bp, sep='-')
      cur_out_fns <- paste0(cur_out_fns, opt$man_suffix)
      cur_out_csv <- file.path(par$ses_man_dir,cur_out_fns)
      
      if (opt$verbose>0)
        cat(glue::glue("[{par$prgmTag}]:{TAB} Writing manifest={cur_out_csv}.{RET}"))
      
      readr::write_csv(cur_ses_tib, cur_out_csv)
      
      if (opt$verbose>0)
        cat(glue::glue("[{par$prgmTag}]:{TAB} Done. Writing manifest={cur_out_csv}.{RET}{RET}"))
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Load Required Probe_IDs (CpGs)::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # ses_dat_list <- sesameData::sesameDataList()
    # 
    # ses_sex_dat <- sesameData::sesameDataGet("sex.inference")
    # ses_det_dat <- sesameData::sesameDataGet("detection.stats")
    # ses_ref_dat <- sesameData::sesameDataGet("ref.methylation")
    # ses_leu_dat <- sesameData::sesameDataGet("leukocyte.betas")
    
    ses_eth_dat <- sesameData::sesameDataGet("ethnicity.inference")
    ses_eth_tib <- dplyr::bind_rows(
      tibble::tibble(Probe_ID=ses_eth_dat$ccs.probes),
      tibble::tibble(Probe_ID=ses_eth_dat$rs.probes)
    ) %>%
      dplyr::mutate(Inference_Name="Ethnicity")
    
    ses_age_dat <- sesameData::sesameDataGet("age.inference")
    ses_age_tib <- dplyr::bind_rows(
      ses_age_dat$SkinBlood %>% as.data.frame() %>% tibble::as_tibble() %>%
        dplyr::rename(Probe_ID=CpGmarker) %>%
        dplyr::filter(Probe_ID!="intercept") %>%
        dplyr::mutate(Inference_Name="SkinBlood_Clock"),
      ses_age_dat$PhenoAge  %>% as.data.frame() %>% tibble::as_tibble() %>%
        dplyr::rename(Probe_ID=CpGmarker) %>%
        dplyr::filter(Probe_ID!="intercept") %>%
        dplyr::mutate(Inference_Name="PhenoAge_Clock")
    )
    
    snp_inf_tib <- sesameData::sesameDataPullVariantAnno_InfiniumI() %>% 
      as.data.frame() %>% rownames_to_column(var="Probe_ID") %>% 
      tibble::as_tibble()
    snp_dir_tib <- sesameData::sesameDataPullVariantAnno_SNP() %>% 
      as.data.frame() %>% rownames_to_column(var="Probe_ID") %>% 
      tibble::as_tibble()
    
    ses_all_tib <- dplyr::bind_rows(
      ses_eth_tib,
      ses_age_tib,
      snp_inf_tib,
      snp_dir_tib
    )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Subset Manifest(s)::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Get Full Manifest and defined required set and control set
    ses_base_tib <- ses_bead_tib %>% 
      dplyr::filter(Probe_ID %in% ses_all_tib$Probe_ID) %>%
      dplyr::bind_rows(ses_ctl_tib) %>% 
      dplyr::distinct(Probe_ID, .keep_all=TRUE) %>% 
      dplyr::arrange(Probe_ID) %>%
      clean_tibble()
    ses_base_sum <- ses_base_tib %>% 
      dplyr::group_by(Probe_Type,Probe_Design) %>% 
      dplyr::summarise(Count=n(), .groups="drop")
    ses_base_sum %>% print(n=base::nrow(ses_base_sum))
    
    if (!is.null(opt$subManifest) && file.exists(opt$subManifest)) {
      man_base_name <- base::basename(opt$subManifest)
      man_base_vecs <- man_base_name %>% 
        stringr::str_remove(opt$man_suffix) %>% 
        stringr::str_split(pattern='-', simplify=TRUE) %>% 
        as.vector()
      
      cur_platform <- man_base_vecs[1]
      cur_version  <- man_base_vecs[2]
      cur_beadpool <- man_base_vecs[3]
      
      if (opt$verbose>0)
        cat(glue::glue("[{par$prgmTag}]: Current ",
                       "platform={cur_platform}, ",
                       "version={cur_version},",
                       "beadpool={cur_beadpool}",
                       "{RET}"))
      
      # Load SubSet Manifest and remove controls and Base Probes::
      top_man_tib <- 
        suppressMessages(suppressWarnings(readr::read_csv(opt$subManifest))) %>%
        dplyr::filter(!stringr::str_starts(Probe_ID,"ct")) %>%
        dplyr::anti_join(ses_base_tib, by="Probe_ID") %>%
        clean_tibble()
      top_man_sum <- top_man_tib %>% 
        dplyr::group_by(Probe_Type,Probe_Design) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      top_man_sum %>% print(n=base::nrow(top_man_sum))
      
      for (rand in c(1:opt$rand)) {
        if (opt$verbose>0)
          cat(glue::glue("[{par$prgmTag}]:{TAB}Current Randomization={rand}...{RET}"))
        
        # Set Randomization with seed::
        set.seed(rand)
        rand_man_tib <- top_man_tib[sample(base::nrow(top_man_tib)), ]
        rand_man_tib %>% print()
        
        inf1_man_tib <- rand_man_tib %>% dplyr::filter(Probe_Design==1)
        inf2_man_tib <- rand_man_tib %>% dplyr::filter(Probe_Design==2)
        
        for (per1 in per1_vec) {
          frac1_cnt <-
            base::as.integer(base::as.double(as.integer(per1)/100) * base::nrow(inf1_man_tib))
          sub1_man_tib <- inf1_man_tib %>% head(n=frac1_cnt)
          
          for (per2 in per2_vec) {
            # cur_man_fns <- paste0(paste(paste0("Rand",rand),per1,per2, sep='-'), opt$man_suffix)
            perc_str    <- paste(per1,per2, sep='.')
            cur_man_fns <- paste0(paste(paste0("Rand",rand),perc_str, sep='-'), opt$man_suffix)
            cur_man_csv <- file.path(run$manDir, cur_man_fns)
            
            frac2_cnt <- 
              base::as.integer(base::as.double(as.integer(per2)/100) * base::nrow(inf2_man_tib))
            sub2_man_tib <- inf2_man_tib %>% head(n=frac2_cnt)
            
            cur_man_tib <- dplyr::bind_rows(
              sub1_man_tib,sub2_man_tib,ses_base_tib
            ) %>% 
              dplyr::distinct(Probe_ID, .keep_all=TRUE) %>%
              dplyr::arrange(Probe_ID)
            
            cur_man_sum <- cur_man_tib %>% 
              dplyr::group_by(Beadpool_ID,Probe_Type,Probe_Design) %>%
              dplyr::summarise(Count=n(), .groups="drop")
            cur_man_sum %>% print(n=base::nrow(cur_man_sum))
            # print_tib(t = cur_man_sum, f=par$prgmTag, v=opt$verbose, n="cur_man_sum")
            
            if (opt$verbose>0)
              cat(glue::glue("[{par$prgmTag}]:{TAB}{TAB}Writing Manifest CSV={cur_man_csv}...{RET}"))
            readr::write_csv(cur_man_tib, cur_man_csv)
          }
          
        }
        if (opt$verbose>0)
          cat(glue::glue("[{par$prgmTag}]:{TAB}Done. Randomization={rand}{RET}{RET}"))
      }
      
      if (opt$verbose>0)
        cat(glue::glue("[{par$prgmTag}]: Done.{RET}{RET}"))
    }
    
  }

}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
