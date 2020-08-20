
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
  opt$annotation  <- paste(
    file.path(par$manDir, opt$build, 'EPIC.hg38.manifest.gencode.v22.tsv.gz'),
    file.path(par$manDir, opt$build, 'HM450.hg38.manifest.gencode.v22.tsv.gz'),
    file.path(par$manDir, opt$build, 'HM27.hg38.manifest.gencode.v22.tsv.gz'),
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



#
# Probe Designs::
#
des_cg_prb_tib <- tib2prbs(tib=head(man_cg_tib, n=40), idsKey="IlmnID", prbKey="Probe_Type", seqKey="Forward_Sequence", 
                           verbose=opt$verbose+10, tt=pTracker)
des_ch_prb_tib <- tib2prbs(tib=head(man_ch_tib, n=40), idsKey="IlmnID", prbKey="Probe_Type", seqKey="Forward_Sequence", 
                           verbose=opt$verbose+10, tt=pTracker)

des_cg_prb_tib <- tib2prbs(tib=man_cg_tib, idsKey="IlmnID", prbKey="Probe_Type", seqKey="Forward_Sequence", 
                           verbose=opt$verbose, tt=pTracker)
des_ch_prb_tib <- tib2prbs(tib=man_ch_tib, idsKey="IlmnID", prbKey="Probe_Type", seqKey="Forward_Sequence", 
                           verbose=opt$verbose, tt=pTracker)



#
# Improbe design input formats::
#
des_cg_tib <- man_cg_tib %>% 
  dplyr::select(IlmnID,Forward_Sequence,Genome_Build,Chromosome,Coordinate,Probe_Type) # %>% dplyr::mutate(CpG_Island="FALSE") %>% purrr::set_names(des_cols)
des_ch_tib <- man_ch_tib %>%
  dplyr::select(IlmnID,Forward_Sequence,Genome_Build,Chromosome,Coordinate,Probe_Type) # %>% dplyr::mutate(CpG_Island="FALSE") %>% purrr::set_names(des_cols)


if (FALSE) {
  #
  # Known Foward Sequence/CG# Mappings...
  #
  
  unq_seq_tsv <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/improbeOutput.37.cgn-topSeq.uniq.tsv.gz'
  des_unq_tib <- readr::read_tsv(unq_seq_tsv)
  
  man_gs_tib %>% dplyr::anti_join(des_unq_tib, by=c("IlmnID"="Seq_ID") ) %>% dplyr::group_by(Probe_Type) %>% dplyr::summarise(Count=n()) %>% print()
  
  man_des_lf_tib <- man_gs_tib %>% dplyr::left_join(des_unq_tib,  by=c("IlmnID"="Seq_ID"), suffix=c("_MAN", "_DES") )
  man_des_in_tib <- man_gs_tib %>% dplyr::inner_join(des_unq_tib, by=c("IlmnID"="Seq_ID"), suffix=c("_MAN", "_DES") )
}

if (FALSE) {
  cgn_gs_cnt <- man_gs_tib %>% dplyr::group_by(IlmnID) %>% dplyr::summarise(Count=n()) %>% dplyr::arrange(-Count)
  cgn_gs_cnt %>% dplyr::group_by(Count) %>% dplyr::summarise(Group_Count=n())
  
  #  man_gs_tib %>% dplyr::group_by(Man_Source) %>% dplyr::summarise(Count=n())
  m27_gs_cnt_tib <- man_gs_tib %>% dplyr::filter(Man_Source=='HumanMethylation27_270596_v.1.2') %>%
    dplyr::left_join(cgn_gs_cnt, by="IlmnID") %>% dplyr::select(IlmnID, Count) %>% dplyr::arrange(Count)
  m27_gs_cnt_tib %>% dplyr::group_by(Count) %>% dplyr::summarise(Group_Count=n())
  #
  # 1156 unique Methyl 27k probes...
  #
  
  #
  # Template to add Top_Sequence to 450k/EPIC
  #
  man_23_gs_tib <- man_gs_tib %>% dplyr::filter(Man_Source!='HumanMethylation27_270596_v.1.2') %>%
    setTopBot_tib(seqKey="Forward_Sequence", srdKey="Fwd_Strand_TB",
                  verbose=opt$verbose+20, tt=pTracker)
  
  man_23_gs_tib %>% dplyr::group_by(Fwd_Strand_TB) %>% dplyr::summarise(Count=n())
  
  man_23_gs_tib2 <- man_23_gs_tib %>% dplyr::filter(Fwd_Strand_TB==-1) %>% 
    setTopBot_tib(seqKey="Forward_Sequence", srdKey="Fwd_Strand_TB2", verbose=opt$verbose+20, tt=pTracker)
  
  man_23_gs_tib2 %>% dplyr::group_by(Fwd_Strand_TB,Fwd_Strand_TB2) %>% dplyr::summarise(Count=n())
  
  man_23_gs_tib3 <- man_gs_tib %>% dplyr::filter(Man_Source!='HumanMethylation27_270596_v.1.2') %>%
    dplyr::filter(IlmnID %in% man_23_gs_tib2$IlmnID) %>%
    setTopBot_tib(seqKey="Forward_Sequence", srdKey="Fwd_Strand_TB3", verbose=opt$verbose+20, tt=pTracker)
  
  des_uniq_tib <- des_seq_tib %>% dplyr::distinct()
}

#
# TBD::
#  - Understand why TB Calling is acting weird
#  - Design all probes based from R code
#  - Match manifest designs to Rcode designs (InfI short by a 5' base)
#

#
# End goal is normalized manifest joined with designs
#

#
# Example:: How to Add Forward TOP/BOT Strand...
#
# new2_tib <- setTopBot_tib(tib=man2_tib, seqKey="Forward_Sequence", srdKey="Fwd_Strand_TB",
#                           verbose=opt$verbose, tt=pTracker)
#

cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Genome Studio Manifests!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Preprocessing:: Improbe Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: This should actually read the manifests directly and write the improbe input files
#  - FOR NOW: We'll use the pre-computed versions...
#
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Preprocessing:: Improbe Designs
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #



cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Improbe Designs!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Preprocessing:: Content Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

genk_tib <- readr::read_csv('/Users/bbarnes/Documents/Projects/methylation/CustomContent/Genknowme/orders/final.07072020/GenKnowme_CpG_SNP_order.07072020.csv')
elly_tib <- readr::read_csv('/Users/bbarnes/Documents/Projects/methylation/ElysiumHealth/targets/EH-cg.txt.gz')
legx_tib <- readr::read_csv('/Users/bbarnes/Documents/Projects/methylation/EWAS/data/LEGX/epicplus-important-probes.csv') %>% 
  tidyr::separate(probe_name, into=c('Probe_ID', 'IF', 'FR', 'CO', 'RP', 'GN'), sep='_')

dplyr::right_join(
  elly_tib %>% dplyr::distinct(Probe_ID),
  legx_tib %>% dplyr::distinct(Probe_ID),
  by="Probe_ID"
)

cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Content Files!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               Preprocessing:: Sesame Manifests (Annotation)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_ses_vec <- opt$annotation %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
man_ses_tib <- lapply(man_ses_vec, loadManifestSource, addSource=TRUE, verbose=opt$verbose, vt=1,tc=1,tt=pTracker) %>%
  dplyr::bind_rows()

cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing::Sesame Manifests (Annotation)!{RET}{RET}") )


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
