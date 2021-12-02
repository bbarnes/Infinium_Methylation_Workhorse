
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Global Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM  <- ","
TAB  <- "\t"
RET  <- "\n"
RET2 <- "\n\n"
BNG  <- "|"
BRK  <- paste0("# ",
               paste(rep("-----",6),collapse=" "),"|",
               paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Define Default Params and Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt <- NULL  # List of user options
par <- NULL  # List of static program parameters (NOT accessible to user)
run <- NULL  # List of non-static run-time program parameters

par$date    <- Sys.Date() %>% as.character()
par$run_mode <- ''
par$maxTest <- NULL

# Default local Mac/sd-isilon directories for ease of use::
par$macDir1 <- '/Users/bbarnes/Documents/Projects/methylation'
par$macDir2 <- '/Users/bretbarnes/Documents'
par$lixDir1 <- '/illumina/scratch/darkmatter'
par$lixDir  <- '/illumina/scratch/darkmatter'

# Program Name Parameters::
par$code_dir <- 'Infinium_Methylation_Workhorse'
par$prgm_dir <- 'workhorse'
par$prgm_tag <- 'workhorse_main_EPIC_EWAS_Content'
cat(glue::glue("[{par$prgm_tag}]: Starting; {par$prgm_tag}.{RET2}"))

opt$verbose <- 3

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Local Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_source_files = function(pars, verbose = 0, funcTag = "load_source_files") {
  
  pars$prgm_tag   <- base::basename( stringr::str_remove( pars$exe_path, ".R$"))
  pars$prgm_dir   <- base::dirname( base::normalizePath( pars$exe_path ) )
  pars$source_dir <- base::dirname( pars$prgm_dir )
  pars$script_dir <- base::dirname( pars$source_dir )
  pars$dat_dir    <- 
    file.path(base::dirname(base::normalizePath(pars$script_dir)), 'dat')
  
  if (verbose>0) {
    cat(glue::glue("[{funcTag}]:   exe_path = {pars$exe_path}.{RET}"))
    cat(glue::glue("[{funcTag}]:   prgm_tag = {pars$prgm_tag}.{RET}"))
    cat(glue::glue("[{funcTag}]:   prgm_dir = {pars$prgm_dir}.{RET}"))
    cat(glue::glue("[{funcTag}]: source_dir = {pars$source_dir}.{RET}"))
    cat(glue::glue("[{funcTag}]: script_dir = {pars$script_dir}.{RET}"))
    cat(glue::glue("[{funcTag}]:    dat_dir = {pars$dat_dir}.{RET}"))
    cat(glue::glue("[{funcTag}]:{RET}"))
  }
  
  func_dir <- file.path(pars$source_dir, 'functions')
  if (!dir.exists(func_dir)) {
    stop(glue::glue("{RET}[{funcTag}]: ERROR: General Source {func_dir} ",
                    "does NOT exist!!!{RET2}"))
    return(NULL)
  }
  
  for (sfile in list.files(path=func_dir, pattern='.R$', full.names=TRUE)) 
    base::source(sfile)
  
  if (verbose>0)
    cat(glue::glue("[{funcTag}]: Done. Loading Source Files from ",
                   "General Source = '{func_dir}'!{RET}") )
  
  opts_source <- 
    file.path(pars$prgm_dir, 'functions', paste0(pars$prgm_tag,"_options.R"))
  if (!file.exists(opts_source)) {
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Program Options Source File ",
                    "{opts_source} does NOT exist!{RET2}"))
    return(NULL)
  }
  base::source(opts_source)
  if (verbose>0)
    cat(glue::glue("[{funcTag}]: Done. Loading Options Source File for ",
                   "Program ({pars$prgm_tag}) Source = '{opts_source}'!{RET2}"))
  
  pars
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Parse Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

args.dat <- commandArgs(trailingOnly = FALSE)
if (args.dat[1]=='RStudio') {
  
  if (dir.exists(par$macDir1)) par$top_dir <- par$macDir1
  if (dir.exists(par$macDir2)) par$top_dir <- par$macDir2
  if (dir.exists(par$lixDir1)) par$top_dir <- par$lixDir1
  if (dir.exists(par$lixDir))  par$top_dir <- par$lixDir

  par$run_mode    <- args.dat[1]
  cat(glue::glue("[{par$prgm_tag}]: Local args.dat[1]={args.dat[1]}.{RET}"))
  cat(glue::glue("[{par$prgm_tag}]: Local    run_mode={par$run_mode}.{RET}"))
  cat(glue::glue("[{par$prgm_tag}]: Local      top_dir={par$top_dir}.{RET}"))
  
  # Default Parameters for local Mac::
  par$exe_path <- file.path(par$top_dir, 'tools', par$code_dir, 'scripts/R', 
                            par$prgm_dir, paste0(par$prgm_tag,'.R') )
  par <- load_source_files(pars = par, verbose = opt$verbose)
  opt <- program_default_options(verbose = opt$verbose)
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda2-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda3-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/conda_4.6.8/bin/Rscript'
  
  #
  # End of local parameter definitions::
  #
  
  opt$out_dir  <- file.path(par$top_dir, 'scratch')
  opt$imp_dir  <- file.path(par$top_dir, 'data/improbe')
  opt$ann_dir  <- file.path(par$top_dir, 'data/annotation')
  opt$man_dir  <- file.path(par$top_dir, 'data/manifests')
  opt$gen_dir  <- file.path(par$top_dir, 'data/imGenomes/Homo_sapiens/NCBI')
  opt$idat_dir <- file.path(par$top_dir, 'data/idats')
  
  opt$bsmap_opt <- "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
  opt$bsmap_exe <- "/Users/bretbarnes/Documents/tools/programs/BSMAPz/bsmapz"
  opt$cgn_seq_dir <- 
    file.path(opt$imp_dir, "scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49-split")
  opt$cgn_bed_dir <- 
    file.path(opt$imp_dir, "scratch/cgnDB/dbSNP_Core4/design-input/min")
  
  opt$canonical_cgn_dir <- file.path(par$dat_dir, "manifest/cgnDB")
  opt$canonical_cgn_csv <- "canonical.cgn-top-grp.csv.gz"
  
  # opt$genome_controls_csv <- file.path(par$dat_dir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
  # opt$sesame_manifest_csv <- file.path(par$top_dir,"data/manifests/methylation/Sesame/EPIC-B4-BP4.manifest.sesame-base.controls-only.csv.gz")
  
  # Pre-defined local options runTypes::
  #
  par$local_runType <- NULL
  par$local_runType <- 'COVIC'
  par$local_runType <- 'GSA'
  par$local_runType <- 'HM450'
  par$local_runType <- 'TruDx'
  par$local_runType <- 'GRCm10'
  par$local_runType <- 'NZT'
  par$local_runType <- 'Chicago'
  par$local_runType <- 'McMaster10Kselection'
  par$local_runType <- 'EPIC_v2'
  par$local_runType <- 'EWAS'
  
  opt$parallel <- TRUE
  
  opt$verbose <- 10
  opt$verbose <- 3
  opt$verbose <- 5

  opt$fresh  <- TRUE
  opt$fresh  <- FALSE
  opt$reload <- TRUE
  
  par$align  <- NULL
  
  if (FALSE) {
    
  } else if (par$local_runType=='EPIC_v2' ||
             par$local_runType=='EWAS') {
    
    par$align <- "p"
    par$align <- "F"

    map_col <- c("Bead_Pool", "MN", "Match_Num", "Bead_Pool_Name", "Bucket_Name", "Order_Path", "AQP1_Num", "AQP2_Num")
    top_dir <- "/Users/bretbarnes/Documents/data/CustomContent/EPIC_v2/AQP-Files-EPIC_and_EWAS-Content"
    ord_dir <- file.path(top_dir, "order")
    mat_dir <- file.path(top_dir, "match")
    aqp_dir <- file.path(top_dir, "aqp")
    
    EWAS_bucket <- c("EWAS_01", "EWAS_02", "EPICv2_EWAS_1-7_SI")
    EPIC_bucket <- c("EWAS_01", "EWAS_02")
    
    map_csv <- file.path(top_dir, "EPICv2EWAS_screening_master-Input-for-Elisa-V4.csv")
    map_tib <- readr::read_csv(map_csv) %>% 
      purrr::set_names(map_col) %>% 
      tibble::as_tibble() %>% 
      tidyr::pivot_longer(cols=c(AQP1_Num,AQP2_Num), names_to="AQP_Num", values_to = "AQP_Name") %>% 
      dplyr::filter(!is.na(AQP_Name)) %>% 
      dplyr::filter(stringr::str_starts(AQP_Name, "BS")) %>%
      dplyr::mutate(AQP_Num=AQP_Num %>%
                      stringr::str_remove("^AQP") %>% 
                      stringr::str_remove("_Num$") %>% 
                      as.integer(),
                    Order_Path=Order_Path %>% 
                      stringr::str_remove("^.*\\\\") %>% 
                      paste0(ord_dir,'/',.,".gz"),
                    Match_Path=dplyr::case_when(
                      AQP_Num == 2 ~ paste0(mat_dir,"/AQP2-",Match_Num,"_probes.match.gz"),
                      AQP_Num == 1 ~ paste0(mat_dir,"/",Match_Num,"_probes.match.gz"),
                      TRUE ~ NA_character_
                    ),
                    AQP_Path = paste0(aqp_dir,"/",AQP_Name,".txt.gz"),
                    # AQP_Path=dplyr::case_when(
                    #   AQP2_Num != "-" ~ paste0(aqp_dir,"/",AQP2_Num,".txt.gz"),
                    #   TRUE ~ paste0(aqp_dir,"/",AQP1_Num,".txt.gz")
                    # ),
                    Order_File_Name = base::basename(Order_Path),
                    Match_File_Name = base::basename(Match_Path),
                    AQP_File_Name = base::basename(AQP_Path) ) %>%
      dplyr::select(Bead_Pool:Bucket_Name,AQP_Num,AQP_Name,
                    Order_File_Name,Match_File_Name,AQP_File_Name,
                    Order_Path,Match_Path,AQP_Path) %>%
      dplyr::arrange(AQP_Num, Bead_Pool)

    lapply(map_tib$Order_Path, file.exists) %>% cbind() %>% as.vector() %>% unique()
    lapply(map_tib$Match_Path, file.exists) %>% cbind() %>% as.vector() %>% unique()
    lapply(map_tib$AQP_Path, file.exists) %>% cbind() %>% as.vector() %>% unique()
    
    # epic_map_tib <- map_tib %>% dplyr::filter(!Bucket_Name %in% EPIC_bucket)
    # ewas_map_tib <- map_tib %>% dplyr::filter( Bucket_Name %in% EWAS_bucket)
    epic_map_tib <- map_tib %>% dplyr::filter(stringr::str_detect(Bucket_Name, "EPIC") )
    ewas_map_tib <- map_tib %>% dplyr::filter(stringr::str_detect(Bucket_Name, "EWAS") )
    
    if (par$local_runType=='EPIC_v2') {
      # EPIC v2:: 
      opt$genome_build <- 'GRCh37'
      opt$platform <- 'EPIC_v2'
      opt$version  <- '1'
      opt$version <- paste0(par$align,opt$version)
      
      opt$ord_dir <- ord_dir
      opt$mat_dir <- mat_dir
      opt$aqp_dir <- aqp_dir
      
      opt$ord_csv <- paste(epic_map_tib$Order_File_Name, collapse = ",")
      opt$mat_tsv <- paste(epic_map_tib$Match_File_Name, collapse = ",")
      opt$aqp_tsv <- paste(epic_map_tib$AQP_File_Name, collapse = ",")
      
    } else if (par$local_runType=='EWAS') {
      opt$genome_build <- 'GRCh37'
      opt$platform <- 'EWAS'
      opt$version  <- '1'
      opt$version <- paste0(par$align,opt$version)
      
      opt$ord_dir <- ord_dir
      opt$mat_dir <- mat_dir
      opt$aqp_dir <- aqp_dir
      
      opt$ord_csv <- paste(ewas_map_tib$Order_File_Name, collapse = ",")
      opt$mat_tsv <- paste(ewas_map_tib$Match_File_Name, collapse = ",")
      opt$aqp_tsv <- paste(ewas_map_tib$AQP_File_Name, collapse = ",")
    }

    opt$sesame_manifest_dat <- "EPIC.hg19.manifest,HM450.hg19.manifest"
    genome_manifest_dir <- file.path(par$top_dir, "data/manifests/methylation/GenomeStudio")
    opt$genome_manifest_csv <- paste(
      file.path(genome_manifest_dir, "MethylationEPIC_v-1-0_B4-Beadpool_ID.csv.gz"),
      file.path(genome_manifest_dir, "HumanMethylation450_15017482_v.1.2.csv.gz"),
      sep = ","
    )
    
  } else if (par$local_runType=='McMaster10Kselection') {
    opt$genome_build <- 'GRCh37'
    opt$platform <- 'MCM'
    opt$version  <- 'v3'
    
    opt$ord_dir <- file.path(par$top_dir, "data/CustomContent/McMaster/McMaster10Kselection/AQP.v2/order")
    opt$mat_dir <- file.path(par$top_dir, "data/CustomContent/McMaster/McMaster10Kselection/AQP.v2/match")
    opt$aqp_dir <- file.path(par$top_dir, "data/CustomContent/McMaster/McMaster10Kselection/AQP.v2/aqp")
    
    opt$ord_csv <- paste('McMaster_CpG_DesignFile_v4.csv.gz',sep=',')
    opt$mat_tsv <- paste('20532820_probes1.match.gz',
                         '20532820_probes2.match.gz', sep=',')
    opt$aqp_tsv <- paste('20051339_A_ProductQC.txt.gz', sep=',')
    
    if (FALSE) {
      opt$aqp_tsv <- paste('BS0033057-AQP1.txt.gz',
                           'BS0033090-AQP2.txt.gz',
                           # '20051339_A_ProductQC.txt.gz',
                           sep=',')
    }
    opt$noob <- paste(
      file.path(par$top_dir, "data/CustomContent/transfer/updated_manifest.csv.gz"),
      sep=',')
    
  } else if (par$local_runType=='Chicago') {
    opt$genome_build <- 'GRCh38'
    opt$genome_build <- 'GRCh37'
    opt$platform <- 'EPIC'
    opt$version  <- 'B3'
    
    opt$idat_dir <- file.path(opt$idat_dir, "idats_Chicago-Ober-Custom")
    par$ord_dir  <- file.path(par$top_dir, 'data/CustomContent/UnivChicago/latest')
    par$mat_dir  <- file.path(par$top_dir, 'data/CustomContent/UnivChicago/latest')
    par$aqp_dir  <- file.path(par$top_dir, 'data/CustomContent/UnivChicago/latest')
    
    opt$ord_csv <- paste('UofChicago-A_A_Array-CpG-order-FINAL.csv', sep=',')
    opt$mat_tsv <- paste('20504790_probes.match.tsv', sep=',')
    opt$aqp_tsv <- paste('329922X374054_A_ProductQC.txt', sep=',')
    
    # Use this later in the process for picking coordinates::
    par$ord_pos_csv <- file.path(par$top_dir, "data/CustomContent/UnivChicago/improbe_input/CpGs_UnivChicago_alldesigns_55860sites.cgn-pos-srd-prbs.tsv.gz")
    
  } else {
    stop(glue::glue("{RET}[{par$prgm_tag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$run_name <- paste(par$local_runType,opt$platform,opt$version,opt$genome_build, sep='-')
  
} else {
  
  par$run_mode <- 'CommandLine'
  par$exe_path <- args.dat[grep("--file=", args.dat)] %>%
    stringr::str_remove("^.*=") %>% head( n = 1)
  # par <- load_source_files(pars = par, verbose = opt$verbose)
  par <- load_source_files(pars = par, verbose = opt$verbose)
  
  opt  <- program_options(verbose = opt$verbose)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Program Initialization::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_reqs <- c('run_mode','prgm_tag','source_dir','dat_dir','exe_path')
opt_reqs <- c('out_dir','imp_dir',
              'genome_build','platform','version','bsmap_exe',
              'Rscript','verbose')

opt <- program_init(name=par$prgm_tag,
                    opts=opt, opt_reqs=opt_reqs, 
                    pars=par, par_reqs=par_reqs,
                    libs=TRUE,rcpp=FALSE,
                    verbose=opt$verbose,vt=3,tc=0,tt=NULL)

par_tib <- par %>%
  unlist(recursive = TRUE) %>%
  dplyr::bind_rows() %>% 
  tidyr::gather("Params", "Value")

opt_tib <- opt %>%
  dplyr::bind_rows() %>% 
  tidyr::gather("Option", "Value")

pTracker <- timeTracker$new()
pTracker$addFile(opt$time_org_txt)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Pre-processing:: 
#                      Pre-defined & Intermediate Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgm_tag}]: Checking pre-defined files...{RET}"))

# For error ledger for accumulating failed probes and why::
error_ledger <- NULL

run <- NULL
run <- get_run_defaults(fresh = opt$fresh,
                        genome_build = opt$genome_build,
                        cgn_seq_dir  = opt$cgn_seq_dir,
                        cgn_bed_dir  = opt$cgn_bed_dir,
                        canonical_cgn_dir = opt$canonical_cgn_dir,
                        canonical_cgn_csv = opt$canonical_cgn_csv,
                        verbose = opt$verbose )

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgm_tag}]: Done. Checking pre-defined files.{RET2}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        0.0 Validate AQP Inputs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

valid_files <- NULL
valid_files <- valid_aqp_inputs(ord_dir = opt$ord_dir, ord_csv = opt$ord_csv,
                                mat_dir = opt$mat_dir, mat_tsv = opt$mat_tsv,
                                aqp_dir = opt$aqp_dir, aqp_tsv = opt$aqp_tsv,
                                verbose = opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#          0.1 Load any pre-defined Noob-Masked Manifest to be added::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# noob_ctl_tib <- NULL
# if (!is.null(opt$noob)) {
#   noob_ctl_tib <- noob_mask(noob_csv = opt$noob, 
#                             ctl_csv = opt$sesame_manifest_csv, 
#                             verbose=opt$verbose, tt=pTracker)
# }

#
# TBD:: User should be able to rebuild an existing or old manifest,
#  or add manifests together...
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      0.2 Load any other pre-defined data
#                        dbCGN, improbe, imGenomes data::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

imGenome_tib <- NULL
imGenome_tib <- load_imGenomes_table(dir = opt$gen_dir,
                                     genome_build = opt$genome_build, 
                                     ret_list = FALSE, 
                                     load_chroms = opt$align_chroms,
                                     verbose = opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#            0.0 Load any pre-defined Standard Manifest to be added::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (!is.null(opt$sesame_manifest_dat)) {
  
  sesame_address_list <- get_file_list(files=opt$sesame_manifest_dat, 
                                       alpha_numeric = TRUE, del = COM)
  
  found_all <- TRUE
  for (sesame_key in names(sesame_address_list)) {
    out_csv <- file.path(opt$out_dir, "Sesame_Address", paste0(sesame_key,"address.csv.gz"))
    if (!file.exists(out_csv)) found_all <- FALSE
  }
  
  if (!found_all) {
    
    sesame_address_dat  <- lapply(sesame_address_list, load_sesame_repo_address,
                                  add_decoy = TRUE,
                                  add_masks = TRUE,
                                  verbose=opt$verbose, tt=pTracker)
    
    sesame_comp_tib <- manifest_column_agreement(
      sesame_address_dat[[1]], sesame_address_dat[[2]],
      verbose = opt$verbose, tt = pTracker
    )
    
    for (sesame_key in names(sesame_address_list)) {
      out_csv <- file.path(opt$out_dir, "Sesame_Address", paste0(sesame_key,"address.csv.gz"))
      
      safe_write(sesame_address_dat[[sesame_key]], file = out_csv, done = TRUE,
                 verbose = opt$verbose, tt = pTracker)
    }
  } else {
    sesame_address_dat <- NULL
    for (sesame_key in names(sesame_address_list)) {
      out_csv <- file.path(opt$out_dir, "Sesame_Address", paste0(sesame_key,"address.csv.gz"))
      sesame_address_dat[[sesame_key]] <- 
        safe_read(out_csv, verbose = opt$verbose, tt = pTracker)
    }
  }
}

if (FALSE && !is.null(opt$genome_manifest_csv)) {
  
  genome_manifest_list <- get_file_list(files=opt$genome_manifest_csv,
                                        trim = c(".csv.gz"), 
                                        alpha_numeric = TRUE, del = COM)
  
  found_all <- TRUE
  for (genome_key in names(genome_manifest_list)) {
    out_csv <- file.path(opt$out_dir, "Genome_Studio_Address", paste0(genome_key,"address.csv.gz"))
    if (!file.exists(out_csv)) found_all <- FALSE
  }
  
  if (!found_all) {
    genome_manifest_dat <- lapply(genome_manifest_list, load_genome_studio_address,
                                  load_clean     = TRUE,
                                  load_controls  = TRUE,
                                  write_clean    = TRUE,
                                  overwrite      = TRUE, 
                                  add_annotation = TRUE,
                                  ret_data       = FALSE,
                                  verbose = opt$verbose, tt = pTracker)
    
    genome_studio_comp_tib <- manifest_column_agreement(
      genome_manifest_dat[[1]], genome_manifest_dat[[2]],
      verbose = opt$verbose, tt = pTracker
    )
    
    for (genome_key in names(genome_manifest_list)) {
      out_csv <- file.path(opt$out_dir, "Genome_Studio_Address", paste0(genome_key,"address.csv.gz"))
      
      safe_write(genome_manifest_dat[[genome_key]], file = out_csv, done = TRUE,
                 verbose = opt$verbose, tt = pTracker)
    }
  } else {
    genome_manifest_dat <- NULL
    for (genome_key in names(genome_manifest_list)) {
      out_csv <- file.path(opt$out_dir, "Genome_Studio_Address", paste0(genome_key,"address.csv.gz"))
      genome_manifest_dat[[genome_key]] <- 
        safe_read(out_csv, verbose = opt$verbose, tt = pTracker)
    }
  }
  
}

if (FALSE) {

  #
  # TBD:: Print non-matching columns to find missing matches!
  #
  all_manifest_dat   <- c(sesame_address_dat, genome_manifest_dat)
  all_manifest_len   <- all_manifest_dat %>% length()
  all_manifest_names <- names(all_manifest_dat)
  
  comparison_tib <- NULL
  for (ii in c(1:all_manifest_len)) {
    src_a <- all_manifest_names[ii]
    
    for (jj in c(1:all_manifest_len)) {
      
      if (ii < jj) {
        src_b <- all_manifest_names[jj]
        
        comp_tib <- manifest_column_agreement(
          all_manifest_dat[[src_a]],
          all_manifest_dat[[src_b]],
          src_a, src_b,
          verbose = opt$verbose, tt = pTracker )
        
        comparison_tib <- comparison_tib %>%
          dplyr::bind_rows(comp_tib)
      }
    }
  }
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   1.0 AQP Address Manifest Workflow: 
#                           Order/Match/AQP/PQC
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ord_tib <- NULL
ord_tib <- 
  aqp_mapping_workflow(ord = valid_files$ord_str,
                       mat = valid_files$mat_str,
                       aqp = valid_files$aqp_str,
                       
                       prb_key = run$prb_key,
                       add_key = run$add_key,
                       des_key = run$des_key,
                       din_key = run$din_key,
                       ids_key = run$ids_key,
                       
                       out_dir = opt$out_dir,
                       out_col = run$out_col,
                       run_tag = opt$run_name,
                       re_load = run$re_load,
                       pre_tag = pTracker$file_vec,
                       
                       verbose=opt$verbose, tt=pTracker)

# List of Infinium I Swaped Prb_Des::
ord_flip_sum <- ord_tib %>% 
  dplyr::filter( (Ord_Des=="U" & stringr::str_ends(Ord_Prb, "G") | 
                    (Ord_Des == "M" & stringr::str_ends(Ord_Prb, "A")) ) ) %>% 
  dplyr::group_by(Ord_Idx) %>% dplyr::summarise(COunt=n(), .groups = "drop")

ord_fix_tib <- ord_tib %>% 
  dplyr::mutate(Ord_Des_Orign=Ord_Des,
                Ord_Des=dplyr::case_when(
                  Ord_Des_Orign=="U" & stringr::str_ends(Ord_Prb, "G") ~ "M",
                  Ord_Des_Orign=="M" & stringr::str_ends(Ord_Prb, "A") ~ "U",
                  TRUE ~ Ord_Des
                ))

ord_fix_tib %>% 
  dplyr::filter( (Ord_Des=="U" & stringr::str_ends(Ord_Prb, "G") | 
                    (Ord_Des == "M" & stringr::str_ends(Ord_Prb, "A")) ) ) %>% 
  dplyr::group_by(Ord_Idx) %>% dplyr::summarise(COunt=n(), .groups = "drop")

#
# Steps::
#  - [Done]: Write out intermediates above::
#  - Intersect with Sesame EPIC
#  - Intersect remainder with Sesame 450k
#  - Bsp/Seq/Cgn remainder
#
#  - Import Annotation
#  - Sesame EPIC GRS %in% Annotation
#

sesame_address_all_dat <- dplyr::bind_rows(
  sesame_address_dat[["EPIC_hg19_manifest"]],
  sesame_address_dat[["HM450_hg19_manifest"]] ) %>%
  dplyr::distinct(Probe_ID, Prb_Des, .keep_all = TRUE)

ses_ord_inn <- sesame_address_all_dat %>% 
  dplyr::inner_join(ord_fix_tib %>% dplyr::mutate(Ord_Key_Cgn=stringr::str_remove(Ord_Key,"-.*$")), 
                    by=c("Probe_ID"= "Ord_Key_Cgn",
                         "Prb_Seq" = "Ord_Prb",
                         "Prb_Des" = "Ord_Des",
                         "Prb_Din" = "Ord_Din"),
                    suffix=c("_ses", "_ord")) %>%
  dplyr::rename(Address=Address_ord) %>%
  dplyr::distinct(Address, .keep_all=TRUE) %>%
  dplyr::mutate(Ord_Unq_Key=paste(Ord_Key,Address,sep="_"))
# dplyr::distinct(Ord_Key,Address,Prb_Seq, .keep_all=TRUE)

ses_epic_pos_grs <- tib_to_grs(sesame_address_dat[["EPIC_hg19_manifest"]], ids_key="Probe_ID", verbose = 10)
ses_450k_pos_grs <- tib_to_grs(sesame_address_dat[["HM450_hg19_manifest"]], ids_key="Probe_ID", verbose = 10)
ses_ord_inn_grs  <- tib_to_grs(ses_ord_inn, 
                               ids_key="Ord_Unq_Key", 
                               chr_key = "Chromosome", verbose = 10)

ord_ses_ant_both <- ord_fix_tib %>% 
  dplyr::mutate( Ord_Key_Cgn=stringr::str_remove(Ord_Key,"-.*$") ) %>%
  dplyr::anti_join(sesame_address_all_dat, 
                   by=c("Ord_Des"="Prb_Des",
                        "Ord_Key_Cgn"="Probe_ID",
                        "Ord_Prb"="Prb_Seq",
                        "Ord_Din"="Prb_Din") ) 
# %>% dplyr::select(Ord_Key,Ord_Des,Ord_Din, Ord_Prb,Prb_Seq) %>% dplyr::filter(Ord_Prb != Prb_Seq)

ord_ses_ant_both %>%
  dplyr::group_by(Ord_Des,Ord_Din,Ord_Idx) %>% 
  dplyr::summarise(Count=n(), .groups = "drop") %>% 
  dplyr::arrange(-Count) %>% print(n=100)

if (!is.null(par$align)) {
  if (par$align=="F") {
    ord_tib <- ord_fix_tib
  } else if (par$align=="p") {
    ord_tib <- ord_ses_ant_both
  } else {
    # Do nothing for now... this shouldn't happen...
  }
}

if (FALSE) {
  # Quick Summary on Tri-fecta probes::
  ses_ord_inn %>% 
    dplyr::group_by(MASK_typeINextBaseSwitch) %>% 
    dplyr::summarise(Count=n(), .groups = "drop") %>% print()
  
  #
  # Unique to EPIC_v2:: by CGN Match
  #
  ord_ses_ant_cgn <- ord_fix_tib %>% 
    dplyr::mutate( Ord_Key_Cgn=stringr::str_remove(Ord_Key,"-.*$") ) %>%
    dplyr::anti_join(sesame_address_all_dat, 
                     by=c("Ord_Des"="Prb_Des",
                          "Ord_Key_Cgn"="Probe_ID",
                          # "Ord_Prb"="Prb_Seq",
                          "Ord_Din"="Prb_Din") ) # %>% dplyr::select(Ord_Key,Ord_Des,Ord_Din, Ord_Prb,Prb_Seq) %>% dplyr::filter(Ord_Prb != Prb_Seq)
  ord_ses_ant_cgn %>% 
    dplyr::group_by(Ord_Des,Ord_Din,Ord_Idx) %>% 
    dplyr::summarise(Count=n(), .groups = "drop") %>% 
    dplyr::arrange(-Count) %>% print(n=100)
  
  #
  # Unique to EPIC_v2:: by Prb Match
  #
  ord_ses_ant_prb <- ord_fix_tib %>% 
    dplyr::mutate( Ord_Key_Cgn=stringr::str_remove(Ord_Key,"-.*$") ) %>%
    dplyr::anti_join(sesame_address_all_dat, 
                     by=c("Ord_Des"="Prb_Des",
                          # "Ord_Key_Cgn"="Probe_ID",
                          "Ord_Prb"="Prb_Seq",
                          "Ord_Din"="Prb_Din") ) # %>% dplyr::select(Ord_Key,Ord_Des,Ord_Din, Ord_Prb,Prb_Seq) %>% dplyr::filter(Ord_Prb != Prb_Seq)
  ord_ses_ant_prb %>%
    dplyr::group_by(Ord_Des,Ord_Din,Ord_Idx) %>% 
    dplyr::summarise(Count=n(), .groups = "drop") %>% 
    dplyr::arrange(-Count) %>% print(n=100)
  
  #
  # Unique to EPIC_v2:: by Both Match
  #
  ord_ses_ant_both <- ord_fix_tib %>% 
    dplyr::mutate( Ord_Key_Cgn=stringr::str_remove(Ord_Key,"-.*$") ) %>%
    dplyr::anti_join(sesame_address_all_dat, 
                     by=c("Ord_Des"="Prb_Des",
                          "Ord_Key_Cgn"="Probe_ID",
                          "Ord_Prb"="Prb_Seq",
                          "Ord_Din"="Prb_Din") ) # %>% dplyr::select(Ord_Key,Ord_Des,Ord_Din, Ord_Prb,Prb_Seq) %>% dplyr::filter(Ord_Prb != Prb_Seq)
  ord_ses_ant_both %>%
    dplyr::group_by(Ord_Des,Ord_Din,Ord_Idx) %>% 
    dplyr::summarise(Count=n(), .groups = "drop") %>% 
    dplyr::arrange(-Count) %>% print(n=100)
  
  #
  # Inner Join::CGN
  #
  ord_ses_inn_cgn <- ord_fix_tib %>% 
    dplyr::mutate( Ord_Key_Cgn=stringr::str_remove(Ord_Key,"-.*$") ) %>%
    dplyr::inner_join(sesame_address_all_dat, 
                      by=c("Ord_Des"="Prb_Des",
                           "Ord_Key_Cgn"="Probe_ID",
                           # "Ord_Prb"="Prb_Seq",
                           "Ord_Din"="Prb_Din") ) %>% 
    dplyr::select(Ord_Key,Ord_Des,Ord_Din, Ord_Prb,Prb_Seq) %>% dplyr::filter(Ord_Prb != Prb_Seq)
  
  #
  # Inner Join::Prb
  #
  ord_ses_inn_prb <- ord_fix_tib %>% 
    dplyr::mutate( Ord_Key_Cgn=stringr::str_remove(Ord_Key,"-.*$") ) %>%
    dplyr::inner_join(sesame_address_all_dat, 
                      by=c("Ord_Des"="Prb_Des",
                           # "Ord_Key_Cgn"="Probe_ID",
                           "Ord_Prb"="Prb_Seq",
                           "Ord_Din"="Prb_Din") ) %>% 
    dplyr::select(Ord_Key,Ord_Des,Ord_Din, Ord_Key_Cgn,Probe_ID) %>% dplyr::filter(Ord_Key_Cgn != Probe_ID)
  
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    2.0 Align All Probe Sequence:: BSMAP
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Plot normalized intensity x individual chr hits x binding energy
#

bsp_tib <- NULL
bsp_tib <- bsp_mapping_workflow(ref_fas = NULL,
                                ref_tib = imGenome_tib,
                                can_fas = NULL,
                                can_tib = ord_tib,
                                
                                cgn_src = run$cgn_bed_tsv,
                                
                                ids_key = run$ids_key,
                                unq_key = run$unq_key,
                                prb_key = run$prb_key,
                                des_key = run$des_key,
                                din_key = run$din_key,
                                
                                join_key  = run$ids_key,
                                join_type = "inner",
                                
                                sort    = run$bsp_sort,
                                full    = run$bsp_full,
                                merge   = run$bsp_merge,
                                
                                light   = run$bsp_light,
                                reload  = opt$reload,
                                retData = FALSE,
                                
                                bsp_exe = opt$bsmap_exe,
                                bsp_opt = opt$bsmap_opt,
                                
                                out_dir = opt$out_dir,
                                out_col = run$out_col,
                                run_tag = opt$run_name,
                                re_load = run$re_load,
                                pre_tag = pTracker$file_vec,
                                
                                verbose=opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#         3.0 Intersect Sequences Address and improbe:: U49/M49
#                         CGN Mapping Workflow()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: TOP Priority:: small table [cgn_int, chr_let, position] cgn_int sorted
#   Load into memory and intersect with seq cgns (analogous to bsp position,
#   mapping).
# TBD:: Update cgn_mapping_workflow rank with seq position map AND 
#   bsp cgn/position agreement!
#
seq_tib <- NULL
seq_tib <- 
  seq_mapping_workflow(ord_tib = ord_tib,
                       
                       seq_dir   = opt$cgn_seq_dir,
                       pattern_u = run$seq_pattern_U, 
                       pattern_m = run$seq_pattern_M,
                       
                       prb_key = run$prb_key,
                       add_key = run$add_key,
                       des_key = run$des_key,
                       din_key = run$din_key,
                       ids_key = run$ids_key,
                       
                       prefix = opt$run_name,
                       suffix = run$seq_suffix, 
                       
                       idxA = run$seq_idxA,
                       idxB = run$seq_idxB,
                       
                       reload   = opt$reload,
                       parallel = opt$parallel,
                       
                       out_dir = opt$out_dir,
                       out_col = run$out_col,
                       run_tag = opt$run_name,
                       re_load = run$re_load,
                       pre_tag = pTracker$file_vec,
                       
                       verbose=opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       4.0 Analyze and Assign Cgn:: 
#                      CGN-Map/BSMAP/dbGCGN look-up
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cgn_tib <- NULL
cgn_tib <- 
  cgn_mapping_workflow(ord_tib = ord_tib,
                       bsp_tib = bsp_tib,
                       seq_tib = seq_tib,
                       
                       ids_key = run$ids_key,
                       des_key = run$des_key,
                       din_key = run$din_key,
                       map_key = run$map_key,
                       
                       Cgn_Int = run$Cgn_Int,
                       Can_Cgn = run$Can_Cgn,
                       Ord_Cgn = run$Ord_Cgn,
                       Bsp_Cgn = run$Bsp_Cgn,
                       Imp_Cgn = run$Imp_Cgn,
                       
                       can_csv = run$canonical_cgn_csv,
                       
                       join    = run$cgn_join,
                       merge   = run$cgn_merge,
                       retData = FALSE,
                       
                       out_dir = opt$out_dir,
                       out_col = run$out_col,
                       unq_col = run$unq_col,
                       run_tag = opt$run_name,
                       re_load = run$re_load,
                       pre_tag = pTracker$file_vec,
                       
                       verbose=opt$verbose, tt=pTracker)

if (FALSE) {
  print(ord_tib)
  print(bsp_tib)
  print(seq_tib)
  print(cgn_tib)
}


#
# Initial Code for Resolving Alignments::
#   - Should start with a full alignment of ord_tib...
#
if (FALSE) {

  bsp_cgn_col <- intersect( names(bsp_tib), names(cgn_tib) )
  
  bsp_cgn_inn <- 
    dplyr::inner_join(bsp_tib, cgn_tib, 
                      by=bsp_cgn_col,
                      suffix=c("_bsp","_cgn") 
    ) %>%
    dplyr::mutate(
      Cgn_Mat_Scr=dplyr::case_when(
        Ord_Key_Cgn == Cgn_Str ~ 0,
        TRUE ~ 1 )
    )
  
  bsp_cgn_inn %>% dplyr::group_by(Cgn_Mat_Scr, Ord_Des, Ord_Din, Cgn_Tag) %>%
    dplyr::summarise(Count=n(), .groups = "drop")
  
  #
  # bind_rows(ses_ord_inn, bsp_cgn_inn)
  #   - First map keys::
  #
  
  dplyr::rename(Chromosome=Bsp_Chr, Coordinate=Bsp_Pos)
  
  common_cols <- intersect( names(ses_ord_inn), names(bsp_cgn_inn) )
  
  ses_ord_inn_grs  <- tib_to_grs(ses_ord_inn, 
                                 ids_key="Ord_Unq_Key", 
                                 chr_key = "Chromosome", verbose = 10)
  
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
# Ending the first of this script here and developing the second half in a 
#   completely separate script. This is due to some naming convention 
#   updates. Its just easier to keep them separated!!!
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   5.0 Probe Design Validation via imGenome:: 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$genome_build=="GRCh37" || opt$genome_build=="GRCh38") {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #              6.1 Load Annotations:: EPIC_CORE/UPDATE_CORE/CHROM_HMM
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Could loop over everything and make subdirectories for EPIC,UPDATE,HMM
  #
  # epic_ann_file <- list.dirs(epic_ann_path, full.names = TRUE)[-1] 
  # epic_dir_list <- as.list(epic_ann_file)
  # names(epic_dir_list) <- base::basename(epic_ann_file)
  
  # Candidates:: 
  can_grs <- NULL
  can_grs$ses_epic_mask0 <- ses_epic_pos_grs
  can_grs$ses_450k_mask0 <- ses_450k_pos_grs
  can_grs$ord_epic_inn_mask0 <- ses_ord_inn_grs

  # Set Up Target
  core_ann_path <- file.path(par$top_dir, "scratch/annotation_to_workhorse_bed/methylation-Human-GRCh37-v2")
  # epic_ann_path <- file.path(core_ann_path, "EPIC_CORE/UCSC")
  epic_ann_path <- file.path(core_ann_path, "EPIC_CORE")
  epic_ann_file <- list.files(epic_ann_path,pattern=".bed.gz$",full.names=TRUE, recursive = TRUE)
  epic_fns_list <- as.list(epic_ann_file)
  names(epic_fns_list) <- base::basename(epic_ann_file) %>% stringr::str_remove(".bed.gz")
  
  
  # TBD:: Load GRS all from list of files with lapply...
  epic_grs_list <- lapply(epic_fns_list, ann_to_grs, verbose=opt$verbose,tt=pTracker)
  
  # tmp_tib <- ann_to_grs(epic_fns_list[[1]], verbose = 10)
  
  can_int_list <- NULL
  can_grs_keys <- names(can_grs)
  for (can_grs_key in can_grs_keys) {
    can_key <- "IlmnID"
    cur_int_list <- NULL
    cur_int_list <- c(cur_int_list,
                       lapply(epic_grs_list, intersect_GRS, can=can_grs[[can_grs_key]], 
                              ref_key=NULL,ref_col=NULL,ref_prefix=NULL,ref_red=TRUE,
                              can_key=can_key,can_col=can_key,can_prefix=NULL, 
                              verbose=opt$verbose, tt=pTracker) )
    
    can_int_list[[can_grs_key]] <- cur_int_list
  }
  
  # Now write each annotation to run$ann_dir
  for (name in names(epic_int_list)) {
    out_csv <- file.path(run$ann_dir, paste(opt$runName,name,'annotation.csv.gz', sep='-'))
    
    # TBD:: The writing function should be moved to cgn_mapping_workflow()
    # TBD:: Generate summary coverage stats. This should be done in the function
    #  as well
    
    safe_write(epic_int_list[[name]],"csv",out_csv, funcTag=par$prgmTag,
               verbose=opt$verbose)
  }
  
}

#
#
# The following sections have been moved to seperate script: 
#
#  Probe Design/Validation/Sequence Extraction::
#    c_improbe_workflow()
#    r_improbe_workflow()
#    s_improbe_workflow()
#
#  Manifest Final Generation()
#    build_manifest()
#    build_sesame_manifest()
#    build_genome_studio_manifest()
#    build_minfi_manifest()
#
#  Auxilary Annotation::
#
#    add_bed_annotaitons()
#
#  Sesame Cross Validation::
#  Docker Image Support::
#
# New Script: "Infinium_Methylation_Workhorse/scripts/R/workhorse/workhorse_main_dev_latest.part2.R"
#
#

#
# Build individual parts of the template sequence:: 
#
#                                            iupac
#     up01.up02...up11.up12.up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn12.dn11.dn02.dn01
#                                       Nxb [  C   G  ] Nxb
#
# Probe Design Formulas::
#
# Inf1C                               Nxb60* up61------------------dn58
# Inf2C                                     ext61* dn61-------------------dn12
#
# Inf1O                up12------------------up61 Nxb61*
# Inf2O           up11------------------up60 ext61*
#
#


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
