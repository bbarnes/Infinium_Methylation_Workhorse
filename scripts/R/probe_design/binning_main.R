
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
#                      Define Default Params and Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
opt <- NULL
par <- NULL

# Program Parameters::
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'probe_design'
par$prgmTag <- 'binning_main'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# Illumina based directories::
par$macDir  <- '/Users/bbarnes/Documents/Projects/methylation/tools'
par$lixDir <- '/illumina/scratch/darkmatter'

# Directory Parameters::
opt$outDir    <- NULL

# Run Parameters::
opt$runName   <- NULL
opt$manifest  <- NULL
# opt$sampleCsv <- NULL

# Chip Platform and Version Parameters::
opt$platform <- NULL
opt$version  <- NULL
opt$build    <- NULL

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
  opt$version   <- 'B4'
  opt$build     <- 'hg38'
  
  opt$runName <- Sys.Date() %>% as.character()
  
  par$manDir <- '/Users/bbarnes/Documents/Projects/manifests/methylation'
  opt$manifest <- paste(
    '/Users/bbarnes/Documents/Projects/manifests/methylation/MethylationEPIC_v-1-0_B2.csv.gz',
    '/Users/bbarnes/Documents/Projects/manifests/methylation/HumanMethylation450_15017482_v.1.2.csv.gz',
    '/Users/bbarnes/Documents/Projects/manifests/methylation/HumanMethylation27_270596_v.1.2.csv.gz',
    sep=','
  )
  
  par$annDir <- file.path(par$manDir, 'Sesame')
  opt$ann_gene  <- paste(
    file.path(par$annDir, opt$build, 'EPIC.hg38.manifest.gencode.v22.tsv.gz'),
    file.path(par$annDir, opt$build, 'HM450.hg38.manifest.gencode.v22.tsv.gz'),
    file.path(par$annDir, opt$build, 'HM27.hg38.manifest.gencode.v22.tsv.gz'),
    sep=','
  )
  
  build <- 'hg19'
  opt$ann_tfbs  <- paste(
    file.path(par$annDir, build, 'EPIC.hg19.ENCODE.TFBS.tsv.gz'),
    file.path(par$annDir, build, 'HM450.hg19.ENCODE.TFBS.tsv.gz'),
    file.path(par$annDir, build, 'HM27.hg19.ENCODE.TFBS.tsv.gz'),
    sep=','
  )
  

  opt$outDir <- file.path(par$topDir)
  
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

    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--manifests"), type="character", default=opt$manifests, 
                help="Human provide manifests to bin [default= %default]", metavar="character"),
    
    make_option(c("--sampleCsv"), type="character", default=opt$sampleCsv, 
                help="Human provide sample sheet labeling [default= %default]", metavar="character"),
    
    # Chip Platform and Version Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform name (HM50, EPIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest version (B2, B4, C0) [default= %default]", metavar="character"),
    make_option(c("--build"), type="character", default=opt$build, 
                help="Manifest build (hg19, hg38) [default= %default]", metavar="character"),
    
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
    is.null(opt$runName) || is.null(opt$manifest) ||
    is.null(opt$platform) || is.null(opt$version) || is.null(opt$build) ||
    
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
  if (is.null(opt$manifest))  cat(glue::glue("[Usage]: manifest is NULL!!!{RET}"))
  # if (is.null(opt$sampleCsv)) cat(glue::glue("[Usage]: sampleCsv is NULL!!!{RET}"))

  if (is.null(opt$platform)) cat(glue::glue("[Usage]: platform is NULL!!!{RET}"))
  if (is.null(opt$version))  cat(glue::glue("[Usage]: version is NULL!!!{RET}"))
  if (is.null(opt$build))    cat(glue::glue("[Usage]: build is NULL!!!{RET}"))
  
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
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
pTracker <- timeTracker$new(verbose=opt$verbose)

opt <- setLaunchExe(opts=opt, pars=par, verbose=opt$verbose, vt=5,tc=0)

opt$outDir <- file.path(opt$outDir, par$prgmTag, opt$runName, opt$build)
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
cat(glue::glue("[{par$prgmTag}]: Built; OutDir={opt$outDir}!{RET}") )

cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing!{RET}{RET}") )













# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               Preprocessing:: Genome Studio Manifests
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

des_cols <- c('Seq_ID', 'Sequence', 'Genome_Build', 'Chromosome', 'Coordinate', 'CpG_Island')

man_gs_vec <- opt$manifest %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
man_gs_tib <- lapply(man_gs_vec, loadManifestGenomeStudio, 
                     normalize=TRUE, addSource=TRUE, retType='man', 
                     verbose=opt$verbose, vt=1,tc=1,tt=pTracker) %>% dplyr::bind_rows()

# Split by design types::
man_cg_tib <- man_gs_tib %>% dplyr::filter(Probe_Type=='cg') %>%
  dplyr::distinct(IlmnID, .keep_all=TRUE) %>% 
  dplyr::mutate(Forward_Sequence=dplyr::case_when(
    is.na(Forward_Sequence) ~ Top_Sequence,
    TRUE ~ Forward_Sequence
  ))

man_ch_tib <- man_gs_tib %>% dplyr::filter(Probe_Type=='ch') %>%
  dplyr::distinct(IlmnID, .keep_all=TRUE) %>% 
  dplyr::filter(!is.na(Forward_Sequence))

# All Combined together::
man_cn_tib <- dplyr::bind_rows( man_cg_tib, man_ch_tib ) %>% dplyr::arrange(IlmnID)

# Need to remove 27k probes that do not overlap with 450k/EPIC for now::
#
man_gs_list <- man_gs_tib %>% split(.$Man_Source)
man_27k_tib <- man_gs_list$HumanMethylation27_270596_v.1.2 %>% 
  dplyr::filter(! IlmnID %in% man_gs_list$HumanMethylation450_15017482_v.1.2$IlmnID) %>%
  dplyr::filter(! IlmnID %in% man_gs_list$`MethylationEPIC_v-1-0_B2`$IlmnID)

cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Genome Studio Manifests!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Preprocessing:: Improbe Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: This should actually read the manifests directly and write the improbe input files
#  - FOR NOW: We'll use the pre-computed versions...
#
if (FALSE) {
  opt$ordDir <- '/Users/bbarnes/Documents/Projects/manifests'
  opt$impDir <- file.path(opt$outDir, 'improbe')
  if (!dir.exists(opt$impDir)) dir.create(opt$impDir, recursive=TRUE)
  
  ord_cols <- c('Seq_ID', 'Sequence', 'Genome_Build', 'Chromosome', 'Coordinate', 'CpG_Island')
  imp_inpA_tsv <- file.path(opt$ordDir, 'methylation/designInput.27k.tsv.gz')
  imp_inpB_tsv <- file.path(opt$ordDir, 'methylation/designInput.450k-EPIC-B2.tsv.gz')
  imp_inpC_tsv <- file.path(opt$ordDir, 'methylation/designInput.EPIC-B2.tsv.gz')
  
  imp_inpA_tib <- suppressMessages(suppressWarnings( readr::read_tsv(imp_inpA_tsv, col_names=ord_cols) ))
  imp_inpB_tib <- suppressMessages(suppressWarnings( readr::read_tsv(imp_inpB_tsv, col_names=ord_cols) ))
  imp_inpC_tib <- suppressMessages(suppressWarnings( readr::read_tsv(imp_inpC_tsv, col_names=ord_cols) ))
  
  imp_inp_bind_tib <- dplyr::bind_rows(
    imp_inpC_tib,
    imp_inpB_tib,
    imp_inpA_tib
  )
  imp_inp_uniq_tib <- dplyr::distinct(imp_inp_bind_tib, Seq_ID, .keep_all=TRUE)
  imp_inp_uniq_tsv <- file.path(opt$impDir, 'EPIC-reorder.improbe-input.tsv.gz')
  
  readr::write_tsv(imp_inp_uniq_tib, imp_inp_uniq_tsv)
  
  cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Order Files!{RET}{RET}") )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Preprocessing:: Improbe Designs
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imp_des_tsv <- '/Users/bbarnes/Documents/Projects/methylation/EWAS/improbe/EPIC-reorder.improbe-design.tsv.gz'
imp_des_tib <- suppressMessages(suppressWarnings( readr::read_tsv(imp_des_tsv) )) %>% 
  dplyr::mutate(Min_Final_Score=pmin(Methyl_Final_Score,UnMethyl_Final_Score),
                Methyl_Allele_TB_Strand=stringr::str_sub(Methyl_Allele_TB_Strand, 1,1))

imp_unq_tib <- imp_des_tib %>% dplyr::distinct(Seq_ID, Methyl_Allele_CO_Strand, Methyl_Allele_TB_Strand, UnMethyl_Probe_Sequence)
imp_mis_tib <- imp_unq_tib %>% dplyr::group_by(UnMethyl_Probe_Sequence) %>% dplyr::summarise(Seq_Count=n())
imp_ann_tib <- imp_unq_tib %>% dplyr::left_join(imp_mis_tib, by="UnMethyl_Probe_Sequence")

imp_ann_tib %>% dplyr::filter(Seq_Count!=1, Methyl_Allele_CO_Strand=='C')


cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Improbe Designs!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Preprocessing:: Content Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

trdx_csv <- '/Users/bbarnes/Documents/Projects/methylation/Content/TruDx/missing-13cpgs.csv'
grim_csv <- '/Users/bbarnes/Documents/Projects/methylation/Content/Steve_Horvath/GrimAgeCpGs.csv'
horv_csv <- '/Users/bbarnes/Documents/Projects/methylation/Content/Steve_Horvath/datMiniAnnotation3.csv'
genk_csv <- '/Users/bbarnes/Documents/Projects/methylation/CustomContent/Genknowme/orders/final.07072020/GenKnowme_CpG_SNP_order.07072020.csv'
elly_csv <- '/Users/bbarnes/Documents/Projects/methylation/ElysiumHealth/targets/EH-cg.txt.gz'
legx_csv <- '/Users/bbarnes/Documents/Projects/methylation/EWAS/data/LEGX/epicplus-important-probes.csv'

trdx_tib <- suppressMessages(suppressWarnings( readr::read_csv(trdx_csv) ))
grim_tib <- suppressMessages(suppressWarnings( readr::read_csv(grim_csv) )) %>% dplyr::rename(Probe_ID=Probe)
horv_tib <- suppressMessages(suppressWarnings( readr::read_csv(horv_csv) )) %>% dplyr::rename(Probe_ID=Name)
genk_tib <- suppressMessages(suppressWarnings( readr::read_csv(genk_csv) )) %>% dplyr::rename(Probe_ID=Assay_Design_Id) %>% dplyr::distinct(Probe_ID)
elly_tib <- suppressMessages(suppressWarnings( readr::read_csv(elly_csv) ))
legx_tib <- suppressMessages(suppressWarnings( readr::read_csv(legx_csv) )) %>% 
  dplyr::mutate(Probe_Type=stringr::str_sub(probe_name, 1,2)) %>%
  tidyr::separate(probe_name, into=c('Probe_ID', 'IF', 'FR', 'CO', 'RP', 'GN'), sep='_')

# Getting Extra Fill-in from Phantom::
phat_tib <- loadManifestGenomeStudio(man_gs_vec[2], retType="man") %>% dplyr::filter(!is.na(Phantom)) %>% dplyr::rename(Probe_ID=IlmnID)
# phat_tib %>% dplyr::select(Probe_ID, Phantom) %>% print()

ann_all_tib <- dplyr::bind_rows(
  genk_tib %>% dplyr::distinct(Probe_ID),
  grim_tib %>% dplyr::distinct(Probe_ID),
  horv_tib %>% dplyr::distinct(Probe_ID),
  elly_tib %>% dplyr::distinct(Probe_ID),
  legx_tib %>% dplyr::distinct(Probe_ID),
  phat_tib %>% dplyr::distinct(Probe_ID)
) %>% 
  dplyr::group_by(Probe_ID) %>% dplyr::summarise(Evd_Count=n()) %>% dplyr::arrange(-Evd_Count) %>%
  dplyr::mutate(
    isGenK=dplyr::case_when( Probe_ID %in% genk_tib$Probe_ID ~ TRUE, TRUE ~ FALSE ),
    isGrim=dplyr::case_when( Probe_ID %in% grim_tib$Probe_ID ~ TRUE, TRUE ~ FALSE ),
    isHorv=dplyr::case_when( Probe_ID %in% horv_tib$Probe_ID ~ TRUE, TRUE ~ FALSE ),
    isElly=dplyr::case_when( Probe_ID %in% elly_tib$Probe_ID ~ TRUE, TRUE ~ FALSE ),
    isLEGX=dplyr::case_when( Probe_ID %in% legx_tib$Probe_ID ~ TRUE, TRUE ~ FALSE ),
    isPHAT=dplyr::case_when( Probe_ID %in% phat_tib$Probe_ID ~ TRUE, TRUE ~ FALSE )
    
  )

ann_all_tib %>% dplyr::group_by(Evd_Count) %>% dplyr::summarise(Group_Count=n()) %>% print()
man_sel_cn_tib <- dplyr::inner_join(man_cn_tib, ann_all_tib, by=c("IlmnID"="Probe_ID") )

cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Content Files!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#         Preprocessing:: Sesame Manifests (Annotation): Genes
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$pre_ann_dir <- '/Users/bbarnes/Documents/Projects/manifests/methylation/Sesame/pre-built'

ann_gene_all_csv <- file.path(opt$pre_ann_dir, 'hg38.manifest.gencode.v22.full.csv.gz')
ann_tran_cgn_csv <- file.path(opt$pre_ann_dir, 'hg38.manifest.gencode.v22.tran-cgn-table.csv.gz')
ann_tran_sum_csv <- file.path(opt$pre_ann_dir, 'hg38.manifest.gencode.v22.tran-summary.csv.gz')
ann_gene_sum_csv <- file.path(opt$pre_ann_dir, 'hg38.manifest.gencode.v22.gene-summary.csv.gz')
ann_uniq_sum_csv <- file.path(opt$pre_ann_dir, 'hg38.manifest.gencode.v22.uniq-summary.csv.gz')

if (!file.exists(ann_gene_all_csv) || !file.exists(ann_tran_sum_csv) ||
    !file.exists(ann_gene_sum_csv) || !file.exists(ann_uniq_sum_csv)) {
  
  ann_gene_all_vec <- opt$ann_gene %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
  ann_gene_all_tib <- lapply(ann_gene_all_vec, loadManifestSource, addSource=TRUE, verbose=opt$verbose, vt=1,tc=1,tt=pTracker) %>%
    dplyr::bind_rows()
  readr::write_csv(ann_gene_all_tib, ann_gene_all_csv)
  
  ann_tran_cgn_tib <- ann_gene_all_tib %>% dplyr::select(probeID,transcriptIDs) %>% dplyr::filter(!is.na(transcriptIDs)) %>% 
    dplyr::mutate(trans= stringr::str_split(transcriptIDs, pattern=';', simplify=TRUE) ) %>% dplyr::select(-transcriptIDs) %>% 
    as.data.frame() %>% tidyr::gather(Group, Gene, -probeID) %>% dplyr::select(-Group) %>% tibble::tibble() %>%
    dplyr::filter(!is.na(Gene))
  readr::write_csv(ann_tran_cgn_tib,ann_tran_cgn_csv)
  
  ann_tran_sum_tib <- ann_tran_cgn_tib %>% dplyr::group_by(Gene) %>% dplyr::summarise(Count=n())
  readr::write_csv(ann_tran_sum_tib,ann_tran_sum_csv)
  
  ann_gene_sum_tib <- ann_gene_all_tib %>% dplyr::select(probeID,geneNames) %>% dplyr::filter(!is.na(geneNames)) %>% 
    dplyr::mutate(trans= stringr::str_split(geneNames, pattern=';', simplify=TRUE) ) %>% dplyr::select(-geneNames) %>% 
    as.data.frame() %>% tidyr::gather(Group, Gene, -probeID) %>% dplyr::select(-Group) %>% tibble::tibble() %>%
    dplyr::filter(!is.na(Gene)) %>% dplyr::group_by(Gene) %>% dplyr::summarise(Count=n())
  readr::write_csv(ann_gene_sum_tib,ann_gene_sum_csv)
  
  ann_uniq_sum_tib <- ann_gene_all_tib %>% dplyr::select(probeID,genesUniq) %>% dplyr::filter(!is.na(genesUniq)) %>%
    dplyr::mutate(trans= stringr::str_split(genesUniq, pattern=';', simplify=TRUE) ) %>% dplyr::select(-genesUniq) %>% 
    as.data.frame() %>% tidyr::gather(Group, Gene, -probeID) %>% dplyr::select(-Group) %>% tibble::tibble() %>%
    dplyr::filter(!is.na(Gene)) %>% dplyr::group_by(Gene) %>% dplyr::summarise(Count=n())
  readr::write_csv(ann_uniq_sum_tib,ann_uniq_sum_csv)
  
} else {
  ann_gene_all_tib <- suppressMessages(suppressWarnings( readr::read_csv(ann_gene_all_csv) ))
  if (FALSE) {
    ann_tran_cgn_tib <- suppressMessages(suppressWarnings( readr::read_csv(ann_tran_cgn_csv) ))
    ann_tran_sum_tib <- suppressMessages(suppressWarnings( readr::read_csv(ann_tran_sum_csv) ))
    ann_gene_sum_tib <- suppressMessages(suppressWarnings( readr::read_csv(ann_gene_sum_csv) ))
    ann_uniq_sum_tib <- suppressMessages(suppressWarnings( readr::read_csv(ann_uniq_sum_csv) ))
  }
}
cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Sesame Manifests (Annotation::Genes)!{RET}{RET}") )

#
# Summarize by:: transcriptIDs,geneNames,genesUniq
#






# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#         Preprocessing:: Sesame Manifests (Annotation): TFBS
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  pre_ann_tfbs_csv <- file.path(opt$pre_ann_dir, 'hg19.ENCODE.TFBS.counts.csv.gz')
  if (!file.exists(pre_ann_tfbs_csv)) {
    ann_tfbs_vec <- opt$ann_tfbs %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
    ann_tfbs_tib <- lapply(ann_tfbs_vec, loadManifestSource, addSource=TRUE, verbose=opt$verbose, vt=1,tc=1,tt=pTracker) %>%
      dplyr::bind_rows()
    ann_tfbs_cnt_tib <- ann_tfbs_tib %>% dplyr::group_by(probeID,TF,loc_summit) %>% dplyr::summarise(Chip_Count=n())
    readr::write_csv(ann_tfbs_cnt_tib, pre_ann_tfbs_csv)
  } else {
    ann_tfbs_cnt_tib <- readr::read_csv(pre_ann_tfbs_csv)
  }
  cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Sesame Manifests (Annotation::TFBS)!{RET}{RET}") )
}

#
# TBD: May want to save summary annotation files...
#

#
# TBD: Need to rebuild probes from scratch (mostly to fix Infinium II designs)
#

#
# TBD: Add annotation to binning::
#  - Bin A: Add genes with low counts
#  - Bin B: Add break into bins based on gene counts/TFBS
#


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              Dirty and Quick Selection of Optimal Content::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_cn_des_all_tib <- dplyr::left_join(man_sel_cn_tib,
                                       imp_des_tib %>% dplyr::filter(Methyl_Allele_CO_Strand=='C') %>% dplyr::select(-Probe_Type),
                                       by=c("IlmnID"="Seq_ID"), suffix=c("_Man", "_Des")
) %>% dplyr::mutate(
  Mat_Prb1_Man=stringr::str_sub(AlleleA_ProbeSeq, 2,50) %>% stringr::str_replace_all('R', 'A'), 
  Mat_Prb2_Man=stringr::str_sub(AlleleA_ProbeSeq, 3,50) %>% stringr::str_replace_all('R', 'A'), 
  Mat_Prb1_Des=stringr::str_sub(UnMethyl_Probe_Sequence, 2,50), 
  Mat_Prb2_Des=stringr::str_sub(UnMethyl_Probe_Sequence, 2,49) ) %>% dplyr::filter(! IlmnID %in% man_27k_tib$IlmnID)

mat_des1_tib <- man_cn_des_all_tib %>% dplyr::filter(Mat_Prb1_Man==Mat_Prb1_Des & Forward_Sequence_Man == Forward_Sequence_Des)
mat_des2_tib <- man_cn_des_all_tib %>% dplyr::filter(Mat_Prb2_Man==Mat_Prb2_Des & Forward_Sequence_Man == Forward_Sequence_Des)

man_cn_des_mat_tib <- dplyr::bind_rows( mat_des1_tib,mat_des2_tib ) %>% 
  dplyr::arrange(IlmnID) %>% dplyr::rename(Forward_Sequence=Forward_Sequence_Des) %>% 
  dplyr::select(-Forward_Sequence_Man)

man_cn_des_mat_tib %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise(Count=n()) %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Final Dirty Quick Selection::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$binMax <- 40000

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Partition Preperation:: bit manuel
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

binA_tib <- man_cn_des_mat_tib %>% 
  dplyr::filter(isGenK | isGrim | isHorv) %>% 
  dplyr::distinct(IlmnID, Infinium_Design, .keep_all=TRUE) %>%
  dplyr::bind_rows(
    man_cn_des_mat_tib %>% dplyr::filter(Infinium_Design=='II') %>% arrange(Min_Final_Score) %>% 
      dplyr::filter(Min_Final_Score < 0.30) %>% dplyr::mutate(Infinium_Design='I')
  )
binA_cnt <- base::nrow(binA_tib)
binA_sum <- binA_tib %>% dplyr::group_by(Probe_Type,Infinium_Design) %>% dplyr::summarise(Count=n())
sumA_tib <- manToBeadSummary(binA_tib)

binB_tib <- man_cn_des_mat_tib %>% 
  dplyr::anti_join(binA_tib, by="IlmnID") %>% 
  dplyr::filter(isElly & isLEGX) %>%
  dplyr::filter(Infinium_Design=='I') %>% 
  dplyr::arrange(Min_Final_Score) %>% dplyr::distinct(IlmnID, Infinium_Design, .keep_all=TRUE)
binB_cnt <- base::nrow(binB_tib)
binB_sum <- binB_tib %>% dplyr::group_by(Probe_Type,Infinium_Design) %>% dplyr::summarise(Count=n())
sumB_tib <- manToBeadSummary(binB_tib)

binC_tib <- man_cn_des_mat_tib %>% 
  dplyr::anti_join(binA_tib, by="IlmnID") %>% 
  dplyr::anti_join(binB_tib, by="IlmnID") %>% 
  dplyr::filter(isElly | isLEGX) %>%
  dplyr::arrange(Min_Final_Score) %>% dplyr::distinct(IlmnID, Infinium_Design, .keep_all=TRUE)
binC_cnt <- base::nrow(binC_tib)
binC_sum <- binC_tib %>% dplyr::group_by(Probe_Type,Infinium_Design) %>% dplyr::summarise(Count=n())
sumC_tib <- manToBeadSummary(binC_tib)

binD_tib <- man_cn_des_mat_tib %>% 
  dplyr::anti_join(binA_tib, by="IlmnID") %>% 
  dplyr::anti_join(binB_tib, by="IlmnID") %>% 
  dplyr::anti_join(binC_tib, by="IlmnID") %>% 
  dplyr::arrange(Min_Final_Score) %>% dplyr::distinct(IlmnID, Infinium_Design, .keep_all=TRUE)
binD_cnt <- base::nrow(binD_tib)
binD_sum <- binD_tib %>% dplyr::group_by(Probe_Type,Infinium_Design) %>% dplyr::summarise(Count=n())
sumD_tib <- manToBeadSummary(binD_tib)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Partition Data into Bins::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bin_all_tib <- dplyr::bind_rows(binA_tib,binB_tib,binC_tib,binD_tib) %>%
  dplyr::mutate(Bead_Cost=as.integer( dplyr::case_when(Infinium_Design=='II' ~ 1, Infinium_Design=='I' ~ 2, TRUE ~ NA_real_) ),
                  Bead_Cost_Agg=cumsum(Bead_Cost))
cur_prb_cnt <- bin_all_tib %>% tail(n=1) %>% dplyr::pull(Bead_Cost_Agg) %>% as.integer()
num_bin_cnt <- as.integer( cur_prb_cnt / opt$binMax )
max_prb_cnt <- num_bin_cnt * opt$binMax
# rem_prb_cnt <- cur_prb_cnt - max_prb_cnt

bin_sel_tib <- bin_all_tib %>% dplyr::filter(Bead_Cost_Agg<=max_prb_cnt)
bin_all_cnt <- bin_sel_tib %>% base::nrow()

bin_grp_tib <- bin_sel_tib %>% dplyr::select(IlmnID,Infinium_Design)
bin_grp_cnt <- bin_grp_tib %>% dplyr::distinct() %>% base::nrow()

# Ensure everything was unique... This could be done earlier...
stopifnot(bin_all_cnt==bin_grp_cnt)

sum_all_tib <- bin_sel_tib %>% manToBeadSummary()
sum_bin_cnt <- as.integer(sum_all_tib$Beads/opt$binMax)

chunk_tmp_tib <- dplyr::bind_rows(
  dplyr::mutate(bin_grp_tib, Row=dplyr::row_number()),
  dplyr::mutate(bin_grp_tib, Row=dplyr::row_number()) %>% dplyr::filter(Infinium_Design=='I')
) %>% dplyr::arrange(Row)

par_tag_tib <- split(chunk_tmp_tib, factor(sort(rank( row.names(chunk_tmp_tib) ) %% num_bin_cnt))) %>% 
  dplyr::bind_rows(.id="Partition") %>% 
  dplyr::distinct(IlmnID,Infinium_Design, .keep_all=TRUE)

par_prb_tib <- bin_all_tib %>% dplyr::inner_join(par_tag_tib, by=c("IlmnID","Infinium_Design"))

# Partition Validation Summary Table::
par_val_tib <- par_prb_tib %>% dplyr::group_by(Partition,Infinium_Design) %>% 
  dplyr::summarise(Count=n()) %>% tidyr::spread(Infinium_Design,Count) %>% 
  dplyr::mutate(Total=I+I+II)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Build and Write Partitions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_prb_list <- par_prb_tib %>% split(.$Partition)
for (part in names(par_prb_list)) {
  cat(glue::glue("[par$prgmTag]: Building partition={part}...{RET}"))

  # Allow for dual design targeting::
  tar_prb1_tib <- par_prb_list[[part]] %>% dplyr::filter(Infinium_Design=='I') %>% dplyr::select(IlmnID)
  tar_prb2_tib <- par_prb_list[[part]] %>% dplyr::filter(Infinium_Design=='II') %>% dplyr::select(IlmnID)
  
  bin_tib <- par_prb_list[[part]] %>% dplyr::distinct(Forward_Sequence, .keep_all=TRUE)
  prb_tib <- tib2prbs(tib=bin_tib, idsKey="IlmnID", prbKey="Probe_Type", seqKey="Forward_Sequence", 
                      verbose=opt$verbose+10, tt=pTracker) %>% 
    dplyr::mutate(CMP1_U=stringr::str_to_upper(PRB1_U),
                  CMP1_M=stringr::str_to_upper(PRB1_M),
                  CMP2_D=stringr::str_to_upper(PRB2_D),
                  DesSeqN_Brac=addBrac(DesSeqN))
  
  # Unique, Join and Build Final Key::
  sel_prb_tib <- dplyr::bind_rows(
    bin_tib %>% dplyr::inner_join(prb_tib, by=c("IlmnID","Probe_Type","AlleleA_ProbeSeq"="CMP1_U", "AlleleB_ProbeSeq"="CMP1_M"), 
                                   suffix=c("_Man", "_Prb") ),
    bin_tib %>% dplyr::inner_join(prb_tib, by=c("IlmnID","Probe_Type","AlleleA_ProbeSeq"="CMP2_D"), 
                                   suffix=c("_Man", "_Prb") )
  ) %>% dplyr::mutate(SRD_Str=paste0(Methyl_Allele_TB_Strand,Methyl_Allele_CO_Strand),
                      Seq_ID_Uniq=paste(IlmnID,SRD_Str, sep='_') )
  
  # sel_prb1_tib <- sel_prb_tib %>% dplyr::filter(Infinium_Design=='I')
  # sel_prb2_tib <- sel_prb_tib %>% dplyr::filter(Infinium_Design=='II')
  
  sel_prb1_tib <- sel_prb_tib %>% dplyr::filter(IlmnID %in% tar_prb1_tib$IlmnID)
  sel_prb2_tib <- sel_prb_tib %>% dplyr::filter(IlmnID %in% tar_prb2_tib$IlmnID)
  
  outIdx <- as.integer(part) + 1
  outName <- paste('EPIC-reorder.partition',outIdx, sep='-')
  ann_tib <- formatReorderEPIC(sel_prb1_tib,sel_prb2_tib, dir=opt$outDir, name=outName, verbose=opt$verbose)
  
  ann_sum <- suppressMessages(suppressWarnings( readr::read_csv( file.path(opt$outDir, 'EPIC-reorder.partition-1.summary.csv.gz') ) ))
  ann_cnt <- manToBeadSummary(ann_tib)
  
  
  cat(glue::glue("[par$prgmTag]: Done. Building partition={part}.{RET}{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$opt_csv  <- file.path(opt$outDir, paste(par$prgmTag,'program-options.csv', sep='.') )
opt$par_csv  <- file.path(opt$outDir, paste(par$prgmTag,'program-parameters.csv', sep='.') )
opt$time_csv <- file.path(opt$outDir, paste(par$prgmTag,'time-tracker.csv.gz', sep='.') )

opt_tib  <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
par_tib  <- dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value")
time_tib <- pTracker$time %>% dplyr::mutate_if(is.numeric, list(round), 4)

readr::write_csv(opt_tib, opt$opt_csv)
readr::write_csv(par_tib, opt$par_csv)
readr::write_csv(time_tib, opt$time_csv)

sysTime <- Sys.time()
cat(glue::glue("{RET}[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
