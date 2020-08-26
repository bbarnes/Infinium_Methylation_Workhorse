
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               COVIC Scratch::
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
par$date    <- Sys.Date() %>% as.character()
par$codeDir <- 'Infinium_Methylation_Workhorse'
par$prgmDir <- 'analysis'
par$prgmTag <- 'COVIC_scratch'

cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}; Date={par$date}.{RET}{RET}"))

# Illumina based directories::
par$macDir  <- '/Users/bbarnes/Documents/Projects/methylation/tools'
par$lixDir <- '/illumina/scratch/darkmatter'

# Directory Parameters::
opt$outDir <- NULL
opt$modDir <- NULL

# Run Parameters::
opt$runName <- NULL

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
  
  opt$runName    <- 'COVIC-Set7-06082020'
  opt$trainClass <- paste('nSARSCov2', 'pSARSCov2', sep=',')
  
  # opt$modDir <- file.path(par$topDir, 'build_models/EPIC/C0/Sample_Class/COVIC-Set7-06082020/ind-beta_i-poob-1')
  opt$modDir <- file.path(par$topDir, 'build_models/EPIC/C0/Sample_Class', opt$runName)
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
    make_option(c("-m", "--modDir"), type="character", default=opt$modDir, 
                help="Models directory [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("-r","--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),

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

if (is.null(opt$outDir) || is.null(opt$modDir) ||
    is.null(opt$runName) || 
    is.null(opt$execute) || is.null(opt$single) || is.null(opt$parallel) || is.null(opt$cluster) ||
    
    # is.null(opt$Rscript) || 
    is.null(opt$clean) || is.null(opt$verbose)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  
  opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")
  opt_tib %>% base::print(n=base::nrow(opt_tib) )
  
  cat(glue::glue("{RET}[Usage]: Missing arguements!!!{RET}{RET}") )
  if (is.null(opt$outDir))  cat(glue::glue("[Usage]: outDir is NULL!!!{RET}"))
  if (is.null(opt$modDir))  cat(glue::glue("[Usage]: modDir is NULL!!!{RET}"))
  if (is.null(opt$runName)) cat(glue::glue("[Usage]: runName is NULL!!!{RET}"))

  if (is.null(opt$execute))  cat(glue::glue("[Usage]: execute is NULL!!!{RET}"))
  if (is.null(opt$single))   cat(glue::glue("[Usage]: single is NULL!!!{RET}"))
  if (is.null(opt$parallel)) cat(glue::glue("[Usage]: parallel is NULL!!!{RET}"))
  if (is.null(opt$cluster))  cat(glue::glue("[Usage]: cluster is NULL!!!{RET}"))
  
  if (is.null(opt$clean))   cat(glue::glue("[Usage]: clean is NULL!!!{RET}"))
  # if (is.null(opt$Rscript)) cat(glue::glue("[Usage]: Rscript is NULL!!!{RET}"))
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

trainClass_vec <- opt$trainClass %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

opt$outDir <- file.path(opt$outDir, par$prgmTag, opt$runName, par$date)
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
cat(glue::glue("[{par$prgmTag}]: Built; OutDir={opt$outDir}!{RET}") )

cat(glue::glue("[{par$prgmTag}]: Done. Preprocessing!{RET}{RET}") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                  Main::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cross_pat <- 'crossValidation-calls-summary.csv.gz'

cross_dirs <- list.files(opt$modDir, pattern=cross_pat, recursive=TRUE, full.names=TRUE) %>% 
  dirname() %>% unique() %>% dirname() %>% unique() %>% dirname() %>% unique()

print(cross_dirs)

cross_dat_csv <- file.path(opt$outDir, 'cross_validation.dat.csv.gz')
cross_dat_tib <- NULL

for (cross_dir in cross_dirs) {
  
  sam_group <- cross_dir %>% dirname() %>% dirname() %>% basename()
  ses_group <- cross_dir %>% dirname() %>% basename()
  dml_group <- cross_dir %>% basename()
  
  cat(glue::glue("[{par$prgmTag}]: Train={sam_group}, Filter={ses_group}, DML={dml_group}; dir={cross_dir}.{RET}"))
  cross_csvs <- list.files(cross_dir, pattern=cross_pat, recursive=TRUE, full.names=TRUE)

  for (ii in c(1:length(cross_csvs))) {
    cur_csv <- cross_csvs[ii]
    cur_tib <- suppressMessages(suppressWarnings( readr::read_csv(cur_csv) )) %>% 
      tidyr::separate(CrossGroup, into=c('ML_Method', 'Cross_Partition'), sep='-Part') %>%
      dplyr::mutate(Training_Group={sam_group}, Filter_Group={ses_group}, DML_Group=dml_group)
    
    # print(cur_csv)
    # print(cur_tib)
    
    cross_dat_tib <- cross_dat_tib %>% dplyr::bind_rows(cur_tib)
    
    # break
  }
  
  # cross_sum_tib <- cross_dat_tib %>% 
  #   group_by(Training_Group,Filter_Group,DML_Group, ML_Method,True_Class,Call) %>%
  #   dplyr::summarise(Group_Count_Sum=sum(Group_Count), 
  #                    Total_Count_Group=sum(Total_Count),
  #                    Percent_Avg=mean(Group_Perc)) %>% 
  #   dplyr::ungroup() %>% 
  #   dplyr::mutate(True_Class_Str=dplyr::case_when(
  #     !is.na(True_Class) ~ trainClass_vec[True_Class+1],
  #     TRUE ~ NA_character_))
  # print(cross_sum_tib)
  # full_sum_tib <- dplyr::bind_rows(full_sum_tib, cross_sum_tib)
  
  # curr_sum_tib <- cross_sum_tib %>% 
  #   tidyr::unite(True_Class_Str, True_Class_Str,Call, sep='_') %>%
  #   dplyr::select(-True_Class, -Group_Count_Sum, -Total_Count_Group) %>%
  #   tidyr::spread(True_Class_Str,Percent_Avg)
  # print(curr_sum_tib)
  # 
  # full_sum_tib <- full_sum_tib %>% dplyr::bind_rows(curr_sum_tib)
}
cat("\n")
readr::write_csv(cross_dat_tib, cross_dat_csv)
cat(glue::glue("[{par$prgmTag}]: Done Loaeing Data!{RET}{RET}"))

cross_sum_csv <- file.path(opt$outDir, 'cross_validation.summary.csv.gz')
cross_sum_tib <- cross_dat_tib %>% 
  group_by(Training_Group,Filter_Group,DML_Group, ML_Method,True_Class,Call) %>%
  dplyr::summarise(Group_Count_Sum=sum(Group_Count), 
                   Total_Count_Group=sum(Total_Count),
                   Percent_Avg=mean(Group_Perc)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(True_Class_Str=dplyr::case_when(
    !is.na(True_Class) ~ trainClass_vec[True_Class+1],
    TRUE ~ NA_character_))
print(cross_sum_tib)

readr::write_csv(cross_sum_tib, cross_sum_csv)
cat(glue::glue("[{par$prgmTag}]: Done Building Summary!{RET}{RET}"))

if (FALSE) {
  cat("\n\n\n\n")
  
  cross_cur_csv <- '/Users/bbarnes/Documents/Projects/methylation/fromCluster/COVIC_scratch/COVIC-Set5-10062020/2020-08-26/cross_validation.dat.csv.gz'
  cross_cur_csv <- '/Users/bbarnes/Documents/Projects/methylation/fromCluster/COVIC_scratch/COVIC-Set-17/2020-08-26/cross_validation.dat.csv.gz'
  cross_cur_csv <- '/Users/bbarnes/Documents/Projects/methylation/fromCluster/COVIC_scratch/COVIC-Set1-15052020/2020-08-26/cross_validation.dat.csv.gz'
  cross_dat_tib <- readr::read_csv(cross_cur_csv)
  
  cross_sum_tib <- cross_dat_tib %>% 
    group_by(Training_Group,Filter_Group,DML_Group, ML_Method,True_Class,Call) %>%
    dplyr::summarise(Total_Group_Count=n(),
                     Group_Count_Sum=sum(Group_Count), 
                     Total_Count_Group=sum(Total_Count),
                     Percent_SD=sd(Group_Perc),
                     Percent_Avg=mean(Group_Perc),
                     Percent_Med=median(Group_Perc),
                     Percent_Min=min(Group_Perc),
                     Percent_Max=max(Group_Perc)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(True_Class_Str=dplyr::case_when(
      !is.na(True_Class) ~ trainClass_vec[True_Class+1],
      TRUE ~ NA_character_)) %>%
    tidyr::unite(True_Class_Str, True_Class_Str,Call, sep='_') %>%
    dplyr::select(-True_Class)
    
  
  print(cross_sum_tib)
  
  
  all_sum_tib <- cross_sum_tib %>% 
    dplyr::ungroup() %>% group_by(Training_Group,Filter_Group,DML_Group, ML_Method) %>%
    dplyr::select(-Group_Count_Sum, -Total_Count_Group) %>%
    tidyr::spread(True_Class_Str,Percent_Avg) %>% 
    dplyr::arrange(-nSARSCov2_TP,-pSARSCov2_TP)
  print(all_sum_tib)
  
  # Rscript tools/Infinium_Methylation_Workhorse/scripts/R/analysis/COVIC_scratch.R -o scratch -m /Users/bbarnes/Documents/Projects/methylation/scratch/build_models/EPIC/C0/Sample_Class/COVIC-Set7-06082020 -r COVIC-Set7-06082020
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