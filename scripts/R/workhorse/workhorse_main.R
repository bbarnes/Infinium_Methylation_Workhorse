
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
par$prgm_tag <- 'workhorse_main'
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
  opt$gen_dir  <- file.path(par$top_dir, 'data/iGenomes/Homo_sapiens/NCBI')
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
  par$local_runType <- 'EWAS'
  par$local_runType <- 'EPIC_v2'
  
  opt$parallel <- TRUE
  
  opt$verbose <- 10
  opt$verbose <- 3
  opt$verbose <- 5
  opt$verbose <- 200
  
  opt$fresh  <- TRUE
  opt$fresh  <- FALSE
  opt$reload <- TRUE
  
  if (FALSE) {
    
  } else if (par$local_runType=='EPIC_v2' ||
             par$local_runType=='EWAS') {
    
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
      dplyr::filter(!is.na(AQP1_Num)) %>%
      dplyr::mutate(Order_Path=Order_Path %>% 
                      stringr::str_remove("^.*\\\\") %>% 
                      paste0(ord_dir,'/',.,".gz"),
                    Match_Path=dplyr::case_when(
                      AQP2_Num != "-" ~ paste0(mat_dir,"/AQP2-",Match_Num,"_probes.match.gz"),
                      TRUE ~ paste0(mat_dir,"/",Match_Num,"_probes.match.gz")
                    ),
                    AQP_Path=dplyr::case_when(
                      AQP2_Num != "-" ~ paste0(aqp_dir,"/",AQP2_Num,".txt.gz"),
                      TRUE ~ paste0(aqp_dir,"/",AQP1_Num,".txt.gz")
                    ),
                    Order_File_Name = base::basename(Order_Path),
                    Match_File_Name = base::basename(Match_Path),
                    AQP_File_Name = base::basename(AQP_Path) )
    
    lapply(map_tib$Order_Path, file.exists) %>% cbind() %>% as.vector() %>% unique()
    lapply(map_tib$Match_Path, file.exists) %>% cbind() %>% as.vector() %>% unique()
    lapply(map_tib$AQP_Path, file.exists) %>% cbind() %>% as.vector() %>% unique()
    
    epic_map_tib <- map_tib %>% dplyr::filter(!Bucket_Name %in% EPIC_bucket)
    ewas_map_tib <- map_tib %>% dplyr::filter( Bucket_Name %in% EWAS_bucket)
    
    if (par$local_runType=='EPIC_v2') {
      # EPIC v2:: 
      opt$genome_build <- 'GRCh37'
      opt$platform <- 'EPIC_v2'
      opt$version  <- 'v1'
      
      opt$ord_dir <- ord_dir
      opt$mat_dir <- mat_dir
      opt$aqp_dir <- aqp_dir
      
      opt$ord_csv <- paste(epic_map_tib$Order_File_Name, collapse = ",")
      opt$mat_tsv <- paste(epic_map_tib$Match_File_Name, collapse = ",")
      opt$aqp_tsv <- paste(epic_map_tib$AQP_File_Name, collapse = ",")
      
    } else if (par$local_runType=='EWAS') {
      opt$genome_build <- 'GRCh37'
      opt$platform <- 'EWAS'
      opt$version  <- 'v1'
      
      opt$ord_dir <- ord_dir
      opt$mat_dir <- mat_dir
      opt$aqp_dir <- aqp_dir
      
      opt$ord_csv <- paste(ewas_map_tib$Order_File_Name, collapse = ",")
      opt$mat_tsv <- paste(ewas_map_tib$Match_File_Name, collapse = ",")
      opt$aqp_tsv <- paste(ewas_map_tib$AQP_File_Name, collapse = ",")
    }

    # opt$sesame_manifest_dat <- "EPIC.hg19.manifest,HM450.hg19.manifest"
    # genome_manifest_dir <- file.path(par$top_dir, "data/manifests/methylation/GenomeStudio")
    # opt$genome_manifest_csv <- paste(
    #   file.path(genome_manifest_dir, "MethylationEPIC_v-1-0_B4-Beadpool_ID.csv.gz"),
    #   file.path(genome_manifest_dir, "HumanMethylation450_15017482_v.1.2.csv.gz"),
    #   sep = ","
    # )
    
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
