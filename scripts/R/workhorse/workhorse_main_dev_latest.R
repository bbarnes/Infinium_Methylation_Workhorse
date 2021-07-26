
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Big Science::
#
# Producing noob masked probe ids: digest::digest2int("3", 0L) %% 28000000
#  - library(digest)
#  - digest::digest2int("aasdfdsfadssdssdad234sdf", 0L) %% 100000000 %>% stringr::str_length()
#  - need to run a quick scan test scan, but probably not. Just need an 8
#    digit number. 
#  - need to explain this to the team tomorrow.
#
# Sequencing -> Discovery -> Array Dominace (Can't compete)
# COVID: Sample Prep -> Gary; do you know about Direct Detection
#   - Existing array
#   - GC Content difference with infectious cells
#   - Normalization (Noob) applied to both arrays and sequencing
#   - SPIT F-ING FIRE!!!
#
# VACINE PASSPORT REPORT CARD!!!
#
# Should improve after new cgnDB is built...
#
# use Digest::SHA qw(sha1_hex);
# my $var = 123;
# my $sha1_hash = sha1_hex($var);
# print $sha1_hash;
#

rm(list=ls(all=TRUE))

# Genomic Ranges::
suppressWarnings(suppressPackageStartupMessages( base::require("GenomicRanges",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("plyr")) )
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringi") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("data.table") ))
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
BNG <- "|"

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
par$prgmDir <- 'workhorse'
par$prgmTag <- 'workhorse_main_dev_latest'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))


# Executables::
opt$Rscript <- NULL

# BSMAP Parameters::
opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
opt$bsmap_opt <- "\"-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R\""
opt$bsmap_exe <- "/Users/bretbarnes/Documents/tools/programs/BSMAPz/bsmapz"

# Run Parameters::
opt$runName    <- NULL
opt$Species    <- NULL

# Null Place Holders::
opt$cpg_top_tsv <- NULL
opt$cpg_pos_tsv <- NULL
opt$cph_pos_tsv <- NULL
opt$snp_pos_tsv <- NULL
opt$ord_des_csv <- NULL
opt$canonical_csv <- NULL

# Directories::
opt$outDir  <- NULL
opt$impDir  <- NULL
opt$annDir  <- NULL
opt$genDir  <- NULL
opt$manDir  <- NULL

# Manufacturing Files:: Required
opt$ords <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib
opt$mats <- NULL
opt$aqps <- NULL

opt$mans <- NULL  # Manifest of probes to be added
opt$noob <- NULL  # Manifest of probes to be noob masked
opt$vcfs <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib
opt$beds <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib
opt$snps <- NULL  # -> prb_tib -> bsp_tib -> fwd_tib -> des_tib -> man_tib

# Manufacturing Info:: Required
opt$bpns <- NULL
opt$aqpn <- NULL

# Pre-defined files (Controls & IDAT Validation):: Optional
opt$gs_ctl_csv  <- NULL
opt$ses_ctl_csv <- NULL

opt$idat <- NULL

# Platform/Method Options::
opt$genBuild <- NULL
opt$platform <- NULL
opt$version  <- NULL

# Process Parallel/Cluster Parameters::
opt$single   <- FALSE
opt$parallel <- TRUE
opt$cluster  <- FALSE

# Run-time Files
opt$opt_csv  <- NULL
opt$par_csv  <- NULL
opt$time_csv <- NULL
opt$time_org_txt <- NULL

# Run-time Options::
opt$fresh  <- FALSE
opt$reload <- FALSE

# verbose Options::
opt$verbose <- 3

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Order/Match/AQP/PQC Expected Columns::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_cols <- list()
par_cols$ord <- 
  cols(
    Assay_Design_Id        = col_character(),
    AlleleA_Probe_Id       = col_character(),
    AlleleA_Probe_Sequence = col_character(),
    AlleleB_Probe_Id       = col_character(),
    AlleleB_Probe_Sequence = col_character(),
    Normalization_Bin      = col_character()
  )

par_cols$mat <- 
  cols(
    Plate    = col_character(),
    Row      = col_character(),
    Col      = col_integer(),
    Address  = col_integer(),
    Mod5     = col_character(),
    Sequence = col_character(),
    Mod3     = col_character(),
    Comments = col_character()
  )

par_cols$ma2 <- 
  cols(
    address_names = col_integer(),
    probe_id      = col_character(),
    sequence      = col_character(),
    type_b        = col_character(),
    address_name  = col_integer(),
    bo_seq        = col_character()
  )

par_cols$aqp <- 
  cols(
    Address           = col_integer(),
    Decode_Status     = col_integer(),
    Decode_Error_Code = col_integer(),
    Decode_Score      = col_integer(),
    Func_Status       = col_integer(),
    Func_Error_Code   = col_integer(),
    QC_Action         = col_integer()
  )

par_cols$pqc <- 
  cols(
    Address      = col_integer(),
    Status       = col_integer(),
    Eval_Code    = col_integer(),
    Average_Rep  = col_integer(),
    Expected_Rep = col_integer()
  )

par$ord_col <- par_cols$ord$cols %>% names()

par$mat_col <- par_cols$mat$cols %>% names()
par$ma2_col <- par_cols$ma2$cols %>% names()

par$aqp_col <- par_cols$aqp$cols %>% names()
par$pqc_col <- par_cols$pqc$cols %>% names()

opt$build_manifest <- FALSE
opt$run_improbe   <- FALSE
par$load_ann <- FALSE

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Local Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_source_files = function(dir, verbose, funcTag="load_source_files") {
  gen_src_dir <- file.path(dir, 'functions')
  if (!dir.exists(gen_src_dir))
    stop(glue::glue("[{funcTag}]: General Source={gen_src_dir} ",
                    "does not exist!!!{RET}{RET}"))
  
  for (sfile in list.files(path=gen_src_dir, pattern='.R$', 
                           full.names=TRUE, recursive=TRUE)) base::source(sfile)
  if (verbose>0)
    cat(glue::glue("[{funcTag}]: Done. Loading Source Files form ",
                   "General Source={gen_src_dir}!{RET}{RET}") )
  
  gen_src_dir
}

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
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda2-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/Anaconda3-2019.10-Linux-x86_64/bin/Rscript'
  if (dir.exists(par$lixDir1)) opt$Rscript <- '/illumina/scratch/darkmatter/thirdparty/conda_4.6.8/bin/Rscript'
  
  #
  # End of local parameter definitions::
  #
  
  opt$outDir  <- file.path(par$topDir, 'scratch')
  opt$impDir  <- file.path(par$topDir, 'data/improbe')
  opt$annDir  <- file.path(par$topDir, 'data/annotation')
  opt$manDir  <- file.path(par$topDir, 'data/manifests')
  opt$genDir  <- file.path(par$topDir, 'data/iGenomes/Homo_sapiens/NCBI')
  # opt$idatDir <- file.path(par$topDir, 'data/idats')
  
  opt$bsmap_opt <- "\"-s 10 -v 5 -n 1 -r 2 -V 2\""
  opt$bsmap_opt <- "\"-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R\""
  opt$bsmap_exe <- "/Users/bretbarnes/Documents/tools/programs/BSMAPz/bsmapz"
  
  # Pre-defined local options runTypes::
  #
  par$local_runType <- NULL
  par$local_runType <- 'COVIC'
  par$local_runType <- 'GSA'
  par$local_runType <- 'HM450'
  par$local_runType <- 'TruDx'
  par$local_runType <- 'EWAS'
  par$local_runType <- 'GRCm10'
  par$local_runType <- 'NZT'
  par$local_runType <- 'Chicago'
  par$local_runType <- 'NZT'
  par$local_runType <- 'McMaster10Kselection'
  
  opt$parallel <- TRUE
  
  opt$verbose <- 10
  opt$verbose <- 3
  opt$verbose <- 5
  
  opt$fresh  <- TRUE
  opt$fresh  <- FALSE
  opt$reload <- TRUE
  
  opt$run_improbe    <- TRUE
  opt$build_manifest <- TRUE
  par$load_ann       <- TRUE
  
  if (FALSE) {
    
  } else if (par$local_runType=='McMaster10Kselection') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'MCM'
    opt$version  <- 'v9'
    opt$Species  <- "Human"
    
    #
    # TBD: Rebuild GMAIL with all the new settings in: 
    #  data/CustomContent/McMaster/McMaster10Kselection/AQP.v1
    # NOTE: We're going to run all three versions::
    #
    #  v1 : Check results, then modularize sesame manifest generation
    #  v0 : 
    #  v00: 
    #
    # TBD: Chicago: with opt$ord_des_csv set to default canonical
    #
    # Below are some examples of how things don't match up...
    # ord_tib <- readr::read_csv(opt$ords) %>% 
    #   dplyr::mutate(Ord_Key=Assay_Design_Id %>% stringr::str_remove("_.*$"))
    # man_tib <- readr::read_csv(opt$mans) %>% 
    #   dplyr::mutate(Ord_Key=Probe_ID %>% stringr::str_remove("_.*$"))
    # 
    # ord_tib %>% dplyr::filter(Ord_Key %in% man_tib$Ord_Key)
    # man_tib %>% dplyr::filter(Ord_Key %in% ord_tib$Ord_Key)
    # 
    # cur_tib <- aqp_add_tib %>% 
    #   dplyr::mutate(Ord_Key=Ord_Key %>% stringr::str_remove("_.*$"))
    # 
    # ord_tib %>% dplyr::filter(Ord_Key %in% cur_tib$Ord_Key)
    # man_tib %>% dplyr::filter(Ord_Key %in% cur_tib$Ord_Key)
    # 
    # cur_tib %>% dplyr::filter(Ord_Key %in% ord_tib$Ord_Key)
    # cur_tib %>% dplyr::filter(Ord_Key %in% man_tib$Ord_Key)
    #
    
    par$aqpDir <- file.path(par$topDir, "data/CustomContent/McMaster/McMaster10Kselection/AQP.v2")
    opt$ords <- paste(
      file.path(par$aqpDir, 'McMaster_CpG_DesignFile_v4.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20532820_probes1.match.gz'),
      file.path(par$aqpDir, '20532820_probes2.match.gz'),
      sep=',')
    
    opt$aqps <- paste(
      # file.path(par$aqpDir, 'BS0033057-AQP1.txt.gz'),
      # file.path(par$aqpDir, 'BS0033090-AQP2.txt.gz'),
      file.path(par$aqpDir, '20051339_A_ProductQC.txt.gz'),
      sep=',')
    
    if (FALSE) {
      opt$aqps <- paste(
        file.path(par$aqpDir, 'BS0033057-AQP1.txt.gz'),
        file.path(par$aqpDir, 'BS0033090-AQP2.txt.gz'),
        # file.path(par$aqpDir, '20051339_A_ProductQC.txt.gz'),
        sep=',')
    }
    
    opt$noob <- paste(
      file.path(par$topDir, "data/CustomContent/transfer/updated_manifest.csv.gz"),
      sep=',')
    
    # noob-mask demo::
    if (FALSE) {
      
      # opt$noob <- paste(
      #   file.path(par$aqpDir, "Rand3-S0.060.manifest.sesame-base.cpg-sorted.csv.gz"),
      #   sep=',')
      
      noob_csv <- file.path(par$topDir, "data/CustomContent/McMaster/McMaster10Kselection/AQP.v1",
                            "Rand3-S0.060.manifest.sesame-base.cpg-sorted.csv.gz")
      noob_tib <- suppressMessages(suppressWarnings(readr::read_csv(noob_csv) ))
      head(noob_tib)
      
      mask_tib <- noob_mask_manifest(tib = noob_tib, field = "Probe_ID", verbose = opt$verbose)
      head(mask_tib)
      
    }
    
    opt$bpns <- paste(1, sep=",")
    opt$aqpn <- paste(1, sep=",")
    
  } else if (par$local_runType=='EWAS') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'EWAS'
    opt$version  <- 'A1'
    opt$version  <- 'A2'
    opt$Species  <- "Human"
    
    opt$idat   <- NULL
    opt$ords   <- NULL
    # opt$ords   <- paste(
    #   "/Users/bretbarnes/Documents/data/CustomContent/EWAS/orders/round1/EPIC-reorder.partition-*.order.csv.gz",
    #   sep=','
    # )
    new_files <- NULL
    tsv_files <- list.files(file.path(par$topDir, "data/CustomContent/EPIC_v2/11102020"), pattern=".tsv.gz$", full.names=TRUE)
    csv_files <- c(
      list.files(file.path(par$topDir, "data/CustomContent/EWAS/orders/round1"), pattern=".order.csv.gz$", full.names=TRUE),
      list.files(file.path(par$topDir, "data/CustomContent/EWAS/orders/round2"), pattern=".order.csv.gz$", full.names=TRUE)
    )
    
    if (FALSE) {
      new_dir <- file.path(par$topDir, "data/CustomContent/EPIC_v2/11102020/csv")
      if (!dir.exists(new_dir)) dir.create(new_dir, recursive=TRUE)
      
      for (tsv in tsv_files) {
        cur_nam <- base::basename(tsv) %>% 
          stringr::str_remove('.gz') %>% 
          stringr::str_remove('.tsv') %>%
          stringr::str_remove("^bp_")
        
        cur_nam <- paste0("9",cur_nam)
        
        out_nam <- paste("EPIC-reorder.partition",cur_nam, sep='-')
        out_csv <- file.path(new_dir, paste(out_nam,'order.csv.gz', sep='.'))
        cur_tib <- suppressMessages(suppressWarnings( readr::read_tsv(tsv, guess_max = 50000) ))
        
        cat(glue::glue("[Pre-process]: name={cur_nam}; inp_csv={tsv}{RET}"))
        
        readr::write_csv(cur_tib, out_csv)
        
        cat(glue::glue("[Pre-process]: name={cur_nam}; out_csv={out_csv}{RET}{RET}"))
        # print(cur_tib)
      }
      new_files <- list.files(new_dir, pattern=".csv.gz$", full.names=TRUE)
    }
    
    opt$ords <- c(csv_files)
    par$ordn <- opt$ords %>% base::basename() %>% 
      stringr::str_remove("^EPIC-reorder.partition-") %>%
      stringr::str_remove(".order.csv.gz$")
    
    # opt$mans   <- paste(
    #   file.path(opt$manDir, "genotyping/GSA-24v2-0_A1.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz"),
    # file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz"),
    # sep=',')
    
    # opt$vcfs   <- paste(
    # "/Users/bretbarnes/Documents/data/annotation/dbSNP/dbSNP-151/GRCh37/common_all_20180423.snps.csv.gz",
    # "/Users/bretbarnes/Documents/data/annotation/dbSNP/dbSNP-151/GRCh37/All_20180423.vcf.gz",
    #  sep=',')
    
    # opt$snps   <- paste(
    #   "/Users/bretbarnes/Documents/data/CustomContent/TruDx/SNPs/TruDx_target_SNPs.csv",
    #   sep=',')
    
  } else if (par$local_runType=='TruDx') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'GSA'
    opt$version  <- 'A4'
    opt$Species  <- "Human"
    
    # opt$mans   <- paste(
    #   file.path(opt$manDir, "genotyping/GSA-24v2-0_A1.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz"),
    # file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz"),
    # sep=',')
    
    opt$vcfs   <- paste(
      # "/Users/bretbarnes/Documents/data/annotation/dbSNP/dbSNP-151/GRCh37/common_all_20180423.snps.csv.gz",
      "/Users/bretbarnes/Documents/data/annotation/dbSNP/dbSNP-151/GRCh37/All_20180423.vcf.gz",
      sep=',')
    opt$beds   <- NULL
    opt$snps   <- paste(
      "/Users/bretbarnes/Documents/data/CustomContent/TruDx/SNPs/TruDx_target_SNPs.csv",
      sep=',')
    
  } else if (par$local_runType=='HM450') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'GSA'
    opt$version  <- 'A2'
    opt$Species  <- "Human"
    
    opt$mans   <- paste(
      #   file.path(opt$manDir, "genotyping/GSA-24v2-0_A1.csv.gz"),
      #   file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz"),
      file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz"),
      #   file.path(opt$manDir, "methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz"),
      sep=',')
    
  } else if (par$local_runType=='GSA') {
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'GSA'
    opt$version  <- 'A2'
    opt$Species  <- "Human"
    
    # opt$mans   <- paste(
    #   file.path(opt$manDir, "genotyping/GSA-24v2-0_A1.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz"),
    #   file.path(opt$manDir, "methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz"),
    #   sep=',')
    
  } else if (par$local_runType=='Chicago') {
    opt$genBuild <- 'GRCh38'
    opt$genBuild <- 'GRCh37'
    
    opt$platform <- 'EPIC'
    opt$version  <- 'A1'
    opt$version  <- 'A2'
    opt$version  <- 'A3'
    opt$version  <- 'A4'
    opt$version  <- 'A5'
    opt$version  <- 'A6'
    opt$version  <- 'A7'
    opt$version  <- 'A8'
    opt$version  <- 'A9'
    opt$version  <- 'A10'
    opt$version  <- 'A11'
    opt$version  <- 'A12'
    opt$version  <- 'A13'
    opt$version  <- 'A14'
    opt$version  <- 'A15'
    
    opt$Species  <- "Human"
    
    opt$idat   <- NULL
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/UnivChicago/latest')
    
    opt$ords <- paste(
      file.path(par$aqpDir, 'UofChicago-A_A_Array-CpG-order-FINAL.csv'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20504790_probes.match.tsv'),
      sep=',')
    
    opt$aqps <- paste(
      file.path(par$aqpDir, '329922X374054_A_ProductQC.txt'),
      sep=',')
    
    opt$bpns <- paste(1, sep=",")
    opt$aqpn <- paste(1, sep=",")
    
    # Use this later in the process for picking coordinates::
    #
    # opt$ord_des_csv <- file.path(par$topDir, "data/CustomContent/UnivChicago/improbe_input/CpGs_UnivChicago_alldesigns_55860sites.cgn-pos-srd-prbs.tsv.gz")
    
    par$ord_pos_csv <- file.path(par$topDir, "data/CustomContent/UnivChicago/improbe_input/CpGs_UnivChicago_alldesigns_55860sites.cgn-pos-srd-prbs.tsv.gz")
    
  } else if (par$local_runType=='COVIC') {
    opt$genBuild <- 'GRCh36'
    opt$genBuild <- 'GRCh38'
    opt$genBuild <- 'GRCh37'
    
    opt$platform    <- 'COVIC'
    opt$version     <- 'C0'
    opt$version     <- 'B1'
    opt$version     <- 'C2'
    opt$version     <- 'C11'
    opt$version     <- 'C12'
    
    opt$Species <- "Human"
    
    # opt$pre_man_csv <- '/Users/bretbarnes/Documents/data/manifests/MethylationEPIC_v-1-0_B2.csv.gz'
    
    opt$idat  <- NULL
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/COVID-19_HLA/AQP/COVIC-Host-Immune-Detection')
    opt$ords <- paste(
      file.path(par$aqpDir, 'COVID_EPIC_Round1.03172020.unique.order_AP.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, '20447043_probes.match.gz'),
      sep=',')
    
    opt$aqps <- paste(
      file.path(par$aqpDir, 'BS0032581-AQP.txt.gz'),
      sep=',')
    
    opt$bpns <- paste(1, sep=",")
    opt$aqpn <- paste(1, sep=",")
    
  } else if (par$local_runType=='GRCm10') {
    opt$platform <- 'LEGX'
    opt$version  <- 'C0'
    opt$version  <- 'C30'
    opt$version  <- 'C31'
    opt$version  <- 'C32'
    opt$version  <- 'C0'
    opt$version  <- 'P0'
    
    opt$Species <- "Mouse"
    
    opt$genBuild <- 'GRCm38'
    opt$genBuild <- 'GRCm10'
    
    opt$genDir  <- file.path(par$topDir, 'data/iGenomes/Mus_musculus/NCBI')
    
    # This file exists, but not what we want now::
    # opt$ord_des_csv <- file.path(par$topDir, "data/CustomContent/LifeEpigentics/data/dropbox/merged_with_raw_ordered.cgn-pos-srd-prbs.tsv.gz")
    
    
    # opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genBuild,'-21092020.cgnTop.sorted.tsv') )
    # opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
    # opt$cph_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.chn-sorted.tsv.gz') )
    # opt$snp_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.snp-sorted.tsv.gz') )
    # opt$snp_des_csv <- file.path(opt$impDir, 'cph-snp-designs/LEGX_SpikeIn_Reorder-SNP-Only.designs.csv.gz')
    # opt$cph_des_csv <- file.path(opt$impDir, 'cph-snp-designs/LEGX_SpikeIn_Reorder-CpH-Only.designs.csv.gz')
    
    par$aqpDir <- file.path(par$topDir, 'data/CustomContent/LifeEpigentics/AQP')
    opt$ords <- paste(
      file.path(par$aqpDir, 'orders/Mus_musculus.order_BP1.csv.gz'),
      file.path(par$aqpDir, 'orders/Mus_musculus.order_BP2.csv.gz'),
      file.path(par$aqpDir, 'orders/mm10_LEGX_nonCpG_probes.Jan16-2020.order.csv.gz'),
      file.path(par$aqpDir, 'orders/LEGX_SpikeIn_Reorder-All-06052020.order.withHeader.csv.gz'),
      sep=',')
    
    opt$mats <- paste(
      file.path(par$aqpDir, 'BP1/20420178_AQP1_LifeEpigen_BP1.txt.gz'),
      file.path(par$aqpDir, 'BP2/20420260_AQP1_LifeEpigen_BP2.txt.gz'),
      file.path(par$aqpDir, 'BP3/20420260_AQP2_LifeEpigen_BP2.txt.gz'),
      file.path(par$aqpDir, 'BP4/20455357_AQP1_LifeEpigen_BP4.txt.gz'),
      sep=',')
    
    if (opt$version == 'P0') {
      opt$aqps <- paste(
        file.path(par$aqpDir, 'PQC/20042400_A_ProductQC.txt.gz'),
        sep=',')
      opt$aqpn <- paste(1, sep=",")
      
    } else {
      
      opt$aqps <- paste(
        file.path(par$aqpDir, 'AQP_Copy/BS0032527-AQP.txt.gz'),
        file.path(par$aqpDir, 'AQP_Copy/BS0032533-AQP.txt.gz'),
        file.path(par$aqpDir, 'AQP_Copy/BS0032545-AQP.txt.gz'),
        file.path(par$aqpDir, 'AQP_Copy/BS0032636-AQP.txt.gz'),
        sep=',')
      opt$aqpn <- paste(1,1,2,1, sep=",")
    }
    
    opt$bpns <- paste(1,2,2,3, sep=",")
    
  } else if (par$local_runType=='NZT') {
    opt$genBuild <- 'GRCh36'
    opt$genBuild <- 'GRCh38'
    opt$genBuild <- 'GRCh37'
    opt$platform <- 'NZT'
    opt$Species  <- "Human"
    opt$version  <- 'N0'
    opt$version  <- 'C4'
    
    # These don't do anything yet...
    # opt$cpg_top_tsv <- file.path(opt$impDir, 'designOutput_21092020/cgnTop',  paste0(opt$genBuild,'-21092020.cgnTop.sorted.tsv') )
    # opt$cpg_pos_tsv <- file.path(opt$impDir, 'designOutput_21092020/genomic', paste0(opt$genBuild,'.improbeDesignInput.cgn-sorted.tsv.gz') )
    
    par$combined <- FALSE
    par$combined <- TRUE
    
    if (par$combined) {
      
      par$aqpDir <- file.path(par$topDir, 'data/CustomContent/NZT/decode/combined')
      opt$ords <- paste(
        file.path(par$aqpDir, 'order-12-combined.csv.gz'),
        sep=',')
      
      opt$mats <- paste(
        file.path(par$aqpDir, 'match-12-combined.match.tsv.gz'),
        sep=',')
      
      opt$aqps <- paste(
        file.path(par$aqpDir, 'aqp-12-combined.txt.gz'),
        sep=',')
      
      opt$bpns <- paste(1, sep=",")
      opt$aqpn <- paste(1, sep=",")
      
    } else {
      par$aqpDir <- file.path(par$topDir, 'data/CustomContent/NZT/decode')
      opt$ords <- paste(
        file.path(par$aqpDir, 'selected.order1.sixCols.csv.gz'),
        file.path(par$aqpDir, 'selected.order2.sixCols.csv.gz'),
        sep=',')
      
      opt$mats <- paste(
        file.path(par$aqpDir, '20297484_probes.match1.tsv.gz'),
        file.path(par$aqpDir, '20297484_probes.match2.tsv.gz'),
        sep=',')
      
      opt$aqps <- paste(
        file.path(par$aqpDir, 'BS0031918-AQP1.txt.gz'),
        file.path(par$aqpDir, 'BS0032272-AQP2.txt.gz'),
        sep=',')
      
      opt$bpns <- paste(1,2, sep=",")
      opt$aqpn <- paste(1,2, sep=",")
    }
    
    opt$idat <- NULL
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  opt$runName <- paste(par$local_runType,opt$platform,opt$version,opt$genBuild, sep='-')
  
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
    # Executables::
    make_option(c("--Rscript"), type="character", default=opt$Rscript, 
                help="Rscript path [default= %default]", metavar="character"),
    
    make_option(c("--bsmap_opt"), type="character", default=opt$bsmap_opt, 
                help="BSMAP Options [default= %default]", metavar="character"),
    make_option(c("--bsmap_exe"), type="character", default=opt$bsmap_exe, 
                help="BSMAP Executable path [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--Species"), type="character", default=opt$Species, 
                help="Target Species Name [default= %default]", metavar="character"),
    make_option(c("--pre_man_csv"), type="character", default=opt$pre_man_csv, 
                help="Previously defined manifest [default= %default]", metavar="character"),
    
    make_option(c("--cpg_pos_tsv"), type="character", default=opt$cpg_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--cpg_top_tsv"), type="character", default=opt$cpg_top_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--cph_pos_tsv"), type="character", default=opt$cph_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--snp_pos_tsv"), type="character", default=opt$snp_pos_tsv, 
                help="Null value for passing arguments [default= %default]", metavar="character"),
    make_option(c("--ord_des_csv"), type="character", default=opt$ord_des_csv, 
                help="Original design file used for canonical position selection. [default= %default]", metavar="character"),
    make_option(c("--canonical_csv"), type="character", default=opt$canonical_csv, 
                help="Canonical CGN names from previous products. [default= %default]", metavar="character"),
    
    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    
    make_option(c("--impDir"), type="character", default=opt$impDir, 
                help="improbe data directory [default= %default]", metavar="character"),
    make_option(c("--annDir"), type="character", default=opt$iannDir, 
                help="Annotation data directory [default= %default]", metavar="character"),
    make_option(c("--genDir"), type="character", default=opt$genDir, 
                help="Genomic data directory [default= %default]", metavar="character"),
    make_option(c("--manDir"), type="character", default=opt$manDir, 
                help="Pre-built Manifest data directory [default= %default]", metavar="character"),
    
    # Manufacturing Files:: Required
    make_option(c("--ords"), type="character", default=opt$ords, 
                help="Order file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mats"), type="character", default=opt$mats, 
                help="Match (format 1 or 2) file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqps"), type="character", default=opt$aqps, 
                help="AQP file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--mans"), type="character", default=opt$mans,
                help="Manifest CSV file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--noob"), type="character", default=opt$noob, 
                help="Manifest of probes to be added and noob-masked CSV file(s) (comma seperated) [default= %default]", metavar="character"),
    
    # Not fully supported yet::
    make_option(c("--vcfs"), type="character", default=opt$vcfs, 
                help="Target Design VCF file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--beds"), type="character", default=opt$beds, 
                help="Target Design Coordinate BED file(s) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--snps"), type="character", default=opt$snps, 
                help="Target Design Coordinate SNP file(s) (comma seperated) [default= %default]", metavar="character"),
    
    # Manufacturing Info:: Required
    make_option(c("--bpns"), type="character", default=opt$bpns, 
                help="Bead Pool Numbers (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--aqpn"), type="character", default=opt$aqpn, 
                help="AQP Numbers (comma seperated) [default= %default]", metavar="character"),
    
    # Pre-defined files (Controls & IDAT Validation):: Optional
    make_option(c("--gs_ctl_csv"), type="character", default=opt$gs_ctl_csv, 
                help="Pre-Defined Infinium Methylation Controls (Genome Studio) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--ses_ctl_csv"), type="character", default=opt$ses_ctl_csv, 
                help="Pre-Defined Infinium Methylation Controls (Sesame) (comma seperated) [default= %default]", metavar="character"),
    make_option(c("--idat"), type="character", default=opt$idat, 
                help="idat directories (comma seperated) [default= %default]", metavar="character"),
    
    # Platform/Method Options::
    make_option(c("--genBuild"), type="character", default=opt$genBuild, 
                help="Genome Build (e.g. GRCh36, GRCh37, GRCh38, GRCm38) [default= %default]", metavar="character"),
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform (e.g. HM450, EPIC, LEGX, NZT, COVIC) [default= %default]", metavar="character"),
    make_option(c("--version"), type="character", default=opt$version, 
                help="Manifest Version (e.g. B0,B1,B2,B3,B4,C0) [default= %default]", metavar="character"),
    
    # Process Parallel/Cluster Parameters::
    make_option(c("--build_manifest"), action="store_true", default=opt$build_manifest, 
                help="Boolean variable to build basic manifest [default= %default]", metavar="boolean"),
    make_option(c("--run_improbe"), action="store_true", default=opt$run_improbe, 
                help="Boolean variable to run improbe [default= %default]", metavar="boolean"),
    
    # Process Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
    # Run-time Files::
    make_option(c("--opt_csv"), type="character", default=opt$opt_csv, 
                help="Unused variable opt_csv [default= %default]", metavar="character"),
    make_option(c("--par_csv"), type="character", default=opt$par_csv, 
                help="Unused variable par_csv [default= %default]", metavar="character"),
    make_option(c("--time_csv"), type="character", default=opt$time_csv, 
                help="Unused variable time_csv [default= %default]", metavar="character"),
    make_option(c("--time_org_txt"), type="character", default=opt$time_org_txt, 
                help="Unused variable time_org_txt [default= %default]", metavar="character"),
    
    # Run=time Options::
    make_option(c("--fresh"), action="store_true", default=opt$fresh, 
                help="Boolean variable to run a fresh build [default= %default]", metavar="boolean"),
    make_option(c("--reload"), action="store_true", default=opt$reload, 
                help="Boolean variable reload intermediate files (for testing). [default= %default]", metavar="boolean"),
    
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
opt_reqs <- c('outDir','impDir','Species',
              'genBuild','platform','version','bsmap_exe',
              'Rscript','verbose')

par$gen_src_dir <- load_source_files(dir=par$scrDir, verbose=opt$verbose)

opt <- program_init(name=par$prgmTag,
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Parse List Options
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run <- NULL
pTracker <- timeTracker$new()

run$image_key <- "bbarnesimdocker/im_workhorse:Infinium_Methylation_Workhorse_Centos"
run$image_ver <- "v.1.25"
run$doc_shell <- "run_improbe.sh"
run$doc_image <- glue::glue("{run$image_key}.{run$image_ver}")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Pre-defined Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Checking pre-defined files.{RET}"))

# Find Genome FASTA Files::
#
run$gen_dir <- file.path(opt$genDir, opt$genBuild,"Sequence/WholeGenomeFasta")
# bsc_ref_pat <- paste0(opt$genBuild,".genome.[FR]C[MUD].fa.gz$")
# bsc_snp_pat <- paste0(opt$genBuild,".genome.dbSNP-151.iupac.[FR]C[MUD].fa.gz$")

# Genome:: Reference
#   - Ref
#   - SNP (dbSNP-151)
# run$gen_ref_fas <- file.path(run$gen_dir, paste0(opt$genBuild,".genome.fa.gz"))
# stopifnot(file.exists(run$gen_ref_fas))
# 
# run$gen_snp_fas <- file.path(run$gen_dir, paste0(opt$genBuild,".dbSNP151-genome.fa.gz"))
# stopifnot(file.exists(run$gen_snp_fas))

all_gen_pat <- paste0(opt$genBuild,".genome.*.fa.gz$")
all_gen_tib <- list.files(run$gen_dir, pattern=all_gen_pat, full.names=TRUE) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("Path")) %>%
  dplyr::mutate(Base_Name=base::basename(Path) %>% stringr::str_remove(".fa.gz$")) %>%
  dplyr::mutate(Unq_ID=stringr::str_remove(Base_Name, paste(opt$genBuild,"genome.",sep=".")), 
                Unq_ID=stringr::str_replace(Unq_ID,"^F",paste(opt$genBuild,"NCBI.dna.F", sep='.') ), 
                Unq_ID=stringr::str_replace(Unq_ID,"^R",paste(opt$genBuild,"NCBI.dna.R", sep='.') ),
                Unq_ID=stringr::str_replace(
                  Unq_ID, paste(opt$genBuild,"genome$", sep='.'), paste(opt$genBuild,"NCBI.dna.FCN", sep='.') ),
                
                
                Unq_ID=stringr::str_replace(Unq_ID,"dbSNP-151.iupac$", "dbSNP-151.iupac.FCN"), 
                Unq_ID=stringr::str_replace(Unq_ID,"dbSNP-151.iupac",paste0(opt$genBuild,".dbSNP-151.snp"))
  ) %>%
  dplyr::select(-Base_Name) %>% 
  tidyr::separate(Unq_ID, into=c("Genome_Build","Source","Alphabet","Genome_Key"), sep="\\.") %>%
  tidyr::separate(Genome_Key, into=c("Strand_FR","Strand_CO","Strand_BSC"), sep=c(1,2), remove=FALSE) %>%
  dplyr::mutate(Genome_Key=paste(Genome_Key,Alphabet, sep="_"))


# Genome:: Bisulfite Converted
#   - Ref (U/M/D)
#   - SNP (dbSNP: U/M/D)
#
# run$bsc_ref_fas <- get_file_list(run$gen_dir, pattern=bsc_ref_pat, 
#                                  trim=c(opt$genBuild,".genome.",".fa.gz"),
#                                  verbose=opt$verbose)
# run$bsc_snp_fas <- get_file_list(run$gen_dir, pattern=bsc_snp_pat, 
#                                  trim=c(opt$genBuild,".genome.",".fa.gz"),
#                                  verbose=opt$verbose)

# Define Pre-built improbe directories and files::

run$imp_prb_dir <- 
  file.path(opt$impDir, "scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49-split")

# run$imp_U49_tsvs <- get_file_list(run$imp_prb_dir, 
#                                   pattern="_U49_cgn-table.tsv.gz", 
#                                   trim=c("^probe_","_U49_cgn-table.tsv.gz"),
#                                   verbose=opt$verbose)
# 
# run$imp_M49_tsvs <- get_file_list(run$imp_prb_dir, 
#                                   pattern="_M49_cgn-table.tsv.gz", 
#                                   trim=c("^probe_","_M49_cgn-table.tsv.gz"),
#                                   verbose=opt$verbose)

# run$imp_prb_dir <- file.path(opt$impDir, "scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49")
# run$imp_u49_tsv <- file.path(run$imp_prb_dir, paste("probe_U49_cgn-table.csv.gz", sep="-") )
# run$imp_m49_tsv <- file.path(run$imp_prb_dir, paste("probe_M49_cgn-table.csv.gz", sep="-") )

# stopifnot(dir.exists(run$imp_prb_dir))
# stopifnot(file.exists(run$imp_u49_tsv))
# stopifnot(file.exists(run$imp_m49_tsv))


run$cgn_bed_dir <- file.path(opt$impDir, "scratch/cgnDB/dbSNP_Core4/design-input/min")
run$cgn_bed_tsv <- file.path(run$cgn_bed_dir, paste(opt$genBuild,"cgn.min.txt.gz", sep="."))
run$canonical_csv <- file.path(par$datDir, "manifest/cgnDB/canonical-assignment.cgn-top-grp.csv.gz")

stopifnot(dir.exists(run$cgn_bed_dir))
stopifnot(file.exists(run$cgn_bed_tsv))
stopifnot(file.exists(run$canonical_csv))

if (is.null(opt$gs_ctl_csv))
  opt$gs_ctl_csv  <- file.path(par$datDir,'manifest/controls/Infinium_Methylation_Controls_15_1983_manifest.noHeader.csv.gz')
if (is.null(opt$ses_ctl_csv))
  opt$ses_ctl_csv <- file.path(par$topDir, "data/manifests/methylation/Sesame/EPIC-B4-BP4.manifest.sesame-base.controls-only.csv.gz")

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Checking pre-defined files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Intermediate Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Defining Intermediate Run Time Files...{RET}"))

# Genomes Parameters::
run$gen_nrec <- 0
run$gen_key <- "Genome_Key"

# Field Parameters:: general
run$del     <- "_"
run$ids_key <- "Prb_Key"
run$unq_key <- "Prb_Key_Unq"

run$add_key <- "Address"
run$din_key <- "Ord_Din"
run$des_key <- "Ord_Des"
run$prb_key <- "Ord_Prb"
run$pos_key <- "Bsp_Pos"
run$chr_key <- "Bsp_Chr"

# Field Parameters:: s-improbe
run$ext_seq="Ext_Forward_Seq"
run$iup_seq="Iupac_Forward_Sequence"
run$imp_seq="Forward_Sequence"

# Field Parameters:: r-improbe
run$srsplit <- TRUE
run$srd_key <- "Bsp_FR"
run$srd_str <- "FR"

run$cosplit <- TRUE
run$cos_key <- "Bsp_CO"
run$cos_str <- "CO"

run$seq_dir <- file.path(opt$outDir, "seq")


# Order Directory::
run$ord_dir <- file.path(opt$outDir, 'ord')
run$ord_suffix  <- "aqp-pass.probe-order"
run$aqp_ord_csv <- 
  file.path(run$ord_dir, paste(opt$runName,run$ord_suffix,"csv.gz", sep='.') )

# CGN/Sequence Intersection Directory::
run$int_dir <- file.path(opt$outDir, 'int')
run$int_suffix  <- "aqp-pass.probe-subseq"
run$int_seq_tsv <- 
  file.path(run$int_dir, paste(opt$runName,run$int_suffix,"tsv.gz", sep='.') )

# BSMAP Alignment Directory::
run$bsp_dir <- file.path(opt$outDir, 'bsp')
run$bsp_suffix  <- "aqp-pass.probe-bsmap-cgn"
run$bsp_cgn_csv  <- file.path(run$bsp_dir, paste(opt$runName,run$bsp_suffix,"csv.gz",  sep='.') )

# run$aqp_prb_bsp  <- file.path(run$bsp_dir, paste(opt$runName,"aqp-pass.address.bsp",  sep='.') )
# run$aqp_bsp_tsv  <- file.path(run$bsp_dir, paste(opt$runName,"aqp-pass.address-bsp.tsv.gz",  sep='.') )
# run$bsp_cgn_csv  <- file.path(run$bsp_dir, paste(opt$runName,"aqp-pass.bsp-cgn.csv.gz",  sep='.') )

# Improbe Design Directory::
run$imp_dir <- file.path(opt$outDir, 'imp')
run$imp_inp_suffix  <- "aqp-pass.improbe-input"
run$imp_out_suffix  <- "aqp-pass.improbe-output"

# run$imp_seq_csv  <- file.path(run$imp_dir, paste(opt$runName, 'improbe-sequence.csv.gz', sep='.') )
# run$imp_inp_tsv  <- file.path(run$imp_dir, paste(opt$runName, 'improbe-inputs.tsv.gz', sep='.') )
# run$imp_des_tsv  <- file.path(run$imp_dir, paste(opt$runName, 'improbe-designOutput.tsv.gz', sep='.') )
# run$imp_fin_tsv  <- file.path(run$imp_dir, paste(opt$runName, 'improbe-designOutput.clean.tsv.gz', sep='.') )

# Manifest Directory::
run$man_dir <- file.path(opt$outDir, 'man')
run$aqp_man_csv  <- file.path(run$man_dir, paste(opt$runName,"aqp-pass.manifest-sesame.csv.gz", sep="."))

# Annotation Directory::
run$ann_dir <- file.path(opt$outDir, 'ann')
run$ann_int_csv  <- file.path(run$ann_dir,paste(opt$runName,'cpg-pass.annotation.csv.gz', sep='.'))

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Defining Intermediate Run Time Files.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Define Global Data Structures::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

aqp_add_tib <- NULL
aqp_seq_tib <- NULL
aqp_cgn_tib <- NULL
aqp_imp_tib <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#            0.0 Load any pre-defined Standard Manifest to be added::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (!is.null(opt$mans)) {
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#          0.1 Load any pre-defined Noob-Masked Manifest to be added::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# noob_ctl_tib <- NULL
# if (!is.null(opt$noob)) {
#   noob_ctl_tib <- noob_mask(noob_csv = opt$noob, 
#                             ctl_csv = opt$ses_ctl_csv, 
#                             verbose=opt$verbose, tt=pTracker)
# }

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#          0.2 Load any pre-defined dbCGN or improbe data::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$build_manifest) {
  
  stamp_vec <- NULL
  stamp_vec <- c(stamp_vec,opt$time_org_txt)
  
  pTracker$addFile(opt$time_org_txt)
  
  par$retData <- TRUE
  par$retData <- FALSE
  # opt$verbose <- 100
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                1.1 New Manifest Workflow: Design/Match/AQP
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (opt$fresh || !valid_time_stamp(stamp_vec)) {
    
    aqp_ord_tib <- 
      aqp_address_workflow(ord = opt$ords,
                           mat = opt$mats,
                           aqp = opt$aqps,

                           prefix = opt$runName,
                           suffix = run$ord_suffix,

                           prb_key = run$prb_key,
                           add_key = run$add_key,
                           des_key = run$des_key,
                           din_key = run$din_key, 
                           ids_key = run$ids_key,
                           
                           del = run$del,
                           
                           out_dir = opt$outDir,
                           run_tag = opt$runName,
                           re_load = TRUE,
                           pre_tag = pTracker$file_vec,
                           verbose=opt$verbose, tt=pTracker)

    
  } else {
    aqp_ord_tib <- safe_read(
      run$aqp_ord_csv, funcTag="aqp-ord", clean=TRUE, guess_max=100000,
      verbose=opt$verbose,tt=pTracker)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #         3.4 Intersect Sequences Address and improbe:: U49/M49
  #                         CGN Mapping Workflow()
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # NOTE: McMaster10Kselection = 822s
  
  stamp_vec <- c(stamp_vec,run$int_seq_tsv)
  

  if (opt$fresh || !valid_time_stamp(stamp_vec)) {
    
    aqp_seq_tib <- 
      cgn_mapping_workflow(tib = aqp_ord_tib, 
                           dir = run$imp_prb_dir, 
                           
                           pattern_u = "-probe_U49_cgn-table.tsv.gz", 
                           pattern_m = "-probe_M49_cgn-table.tsv.gz", 
                           
                           prb_key = run$prb_key,
                           add_key = run$add_key,
                           des_key = run$des_key,
                           din_key = run$din_key,
                           ids_key = run$ids_key,

                           out    = run$int_dir,
                           prefix = opt$runName,
                           suffix = run$int_suffix, 
                           
                           idxA = 1,
                           idxB = 1,
                           
                           reload   = opt$reload,
                           parallel = opt$parallel,
                           
                           del = run$del,
                           
                           out_dir = opt$outDir,
                           run_tag = opt$runName,
                           re_load = TRUE,
                           pre_tag = pTracker$file_vec,
                           verbose=opt$verbose, tt=pTracker)

  } else {
    aqp_seq_tib <- safe_read(
      run$int_seq_tsv, funcTag="aqp-cgn", clean=TRUE, 
      verbose=opt$verbose,tt=pTracker)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    3.2 Align All Probe Sequence:: BSMAP
  #                  3.3 Join Address and Alignment Data:: BSMAP
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # TBD:: Return bsp data and 
  stamp_vec <- c(stamp_vec, run$bsp_cgn_csv)
  if (opt$fresh || !valid_time_stamp(stamp_vec)) {
    
    aqp_bsp_tib <- bsp_mapping_workflow(ref_fas = NULL,
                                        ref_tib = all_gen_tib,
                                        can_fas = NULL,
                                        can_tib = aqp_ord_tib,
                                        seq_tib = aqp_seq_tib,
                                        
                                        cgn_src = run$cgn_bed_tsv,
                                        canonical = run$canonical_csv,
                                        
                                        ids_key = run$ids_key,
                                        unq_key = run$unq_key,
                                        prb_key = run$prb_key,
                                        des_key = run$des_key,
                                        din_key = run$din_key,
                                        
                                        join_key  = run$ids_key,  # Use to be "Aln_Key"
                                        join_type = "inner",

                                        sort  = TRUE,
                                        full  = FALSE,
                                        merge = FALSE,
                                        light = TRUE,
                                        reload = opt$reload,
                                        
                                        bsp_exe = opt$bsmap_exe,
                                        bsp_opt = opt$bsmap_opt,
                                        
                                        out_dir = opt$outDir,
                                        run_tag = opt$runName,
                                        re_load = TRUE,
                                        pre_tag = pTracker$file_vec,
                                        verbose=opt$verbose, tt=pTracker)
    
    
  } else {
    aqp_bsp_tib <- safe_read(
      run$bsp_cgn_csv, funcTag="aqp-bsp", clean=TRUE, guess_max=100000,
      verbose=opt$verbose,tt=pTracker)
  }

  #
  # Workflow:: AQP Manifest Preperation
  #
  #   aqp_address_workflow()
  #   cgn_mapping_workflow()
  #   bsp_mapping_workflow()
  #
  # Workflow:: Probe Design
  #
  #   prb_designs_workflow()
  #     - s_improbe_workflow()
  #     - c_improbe_workflow()
  #     - r_improbe_workflow()
  #
  # Workflow:: Annotation Mapping
  #
  #   prb_annotation_workflow()
  #
  # Workflow:: Manifest Generation
  #
  #   build_sesame_manifest()
  #   build_genome_studio_manifest()
  #
  
  # TBD:: NEXT::
  #
  #   - Minimize aqp_bsp_tib as input to below::
  #
  # aqp_prb_tib2 <- aqp_prb_tib
  aqp_prb_tib <- prb_designs_workflow(
    tib = aqp_bsp_tib %>% 
      dplyr::select(run$unq_key, run$ids_key, run$add_key, 
                    run$des_key, run$din_key, 
                    run$srd_key, run$cos_key,
                    run$pos_key, run$chr_key), 
    max = 0,
    
    out_dir  = run$seq_dir,
    run_name = opt$runName,
    
    imp_level  = 3,
    imp_prefix = opt$runName,
    imp_inp_suffix = run$imp_inp_suffix,
    imp_out_suffix = run$imp_out_suffix,
    
    # imp_inp_tsv = run$imp_inp_tsv,
    # imp_des_tsv = run$imp_des_tsv,
    # imp_fin_tsv = run$imp_fin_tsv,
    # imp_seq_csv = run$imp_seq_csv,
    
    # Genomes Parameters::
    gen_bld  = opt$genBuild,
    gen_nrec = run$gen_nrec,
    gen_key  = run$gen_key,
    gen_tib  = all_gen_tib,
    
    # Field Parameters:: general
    # ids_key = run$ids_key, 
    ids_key = run$unq_key, 
    din_key = run$din_key, 
    pos_key = run$pos_key,
    chr_key = run$chr_key,
    
    # Field Parameters:: s-improbe
    ext_seq = run$ext_seq,
    iup_seq = run$iup_seq,
    imp_seq = run$imp_seq,
    
    # Field Parameters:: r-improbe
    srsplit = run$srsplit,
    srd_key = run$srd_key,
    srd_str = run$srd_str,
    
    cosplit = run$cosplit,
    cos_key = run$cos_key,
    cos_str = run$cos_str,
    
    # Docker Parameters::
    doc_image = run$doc_image,
    doc_shell = run$doc_shell,
    
    join     = FALSE,
    join_new = c("Aln_Key_Unq","Bsp_Chr","Bsp_Pos","Bsp_FR","Bsp_CO"),
    join_old = c("Seq_ID","Chromosome","Coordinate","Strand_FR","Strand_CO"),
    
    subset   = TRUE,
    sub_cols = NULL,
    
    # reload=opt$reload,
    reload  = TRUE,
    retData = TRUE,
    # retData=FALSE,
    
    parallel=opt$parallel,
    # parallel=FALSE,
    
    r_improbe = TRUE,
    s_improbe = FALSE,
    c_improbe = TRUE,
    
    add_flanks = TRUE,
    add_matseq = TRUE,
    
    verbose=opt$verbose, tt=pTracker
  )
  
}











if (FALSE) {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         4.0 improbe fwd design::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  stamp_vec <- c(stamp_vec,
                 run$imp_inp_tsv,
                 run$imp_des_tsv,
                 run$imp_fin_tsv)
  
  if (opt$fresh || !valid_time_stamp(stamp_vec)) {
    
    # TBD:: Add a safty check for Aln_Key_Unq > 15 characters!!!
    # TBD:: Add converted genomes and substring probe checking!!!
    # TBD:: Add r-improbe designer
    #
    
    #    aqp_imp_list <- imp_designs_workflow(
    aqp_imp_tib_pre <- aqp_imp_tib
    aqp_imp_tib <- imp_designs_workflow(
      tib=aqp_cgn_tib, # max=20,
      
      imp_inp_tsv=run$imp_inp_tsv,
      imp_des_tsv=run$imp_des_tsv,
      imp_fin_tsv=run$imp_fin_tsv,
      imp_seq_csv=run$imp_seq_csv,
      
      # Genomes::
      gen_bld=opt$genBuild,
      
      gen_ref_fas=run$gen_ref_fas,
      bsc_ref_fas=run$bsc_ref_fas,
      gen_snp_fas=run$gen_snp_fas,
      bsc_snp_fas=run$bsc_snp_fas,
      
      # Field Parameters::
      ids_key="Aln_Key_Unq",
      din_key="Ord_Din",
      pos_key="Bsp_Pos",
      chr_key="Bsp_Chr",
      
      srsplit=TRUE,
      srd_key="Bsp_FR",
      cosplit=TRUE,
      cos_key="Bsp_CO",
      
      # Docker Parameters::
      run_name=opt$runName,
      doc_image=image_str,
      doc_shell=image_ssh,
      
      join_new=c("Aln_Key_Unq","Bsp_Chr","Bsp_Pos","Bsp_FR","Bsp_CO"),
      join_old=c("Seq_ID","Chromosome","Coordinate","Strand_FR","Strand_CO"),
      
      subset=TRUE,
      sub_cols=NULL,
      
      # reload=opt$reload,
      reload=FALSE,
      retData=TRUE,
      # retData=FALSE,
      
      parallel=opt$parallel,
      # parallel=FALSE,
      r_improbe=TRUE,
      s_improbe=FALSE,
      
      add_flanks=TRUE,
      add_matseq=TRUE,
      
      verbose=opt$verbose, tt=pTracker)
    
  } else {
    
    aqp_imp_tib <- safe_read(
      run$imp_fin_tsv, funcTag="aqp-imp", clean=TRUE,
      verbose=opt$verbose,tt=pTracker)
    
    if (FALSE) {
      imp_fmt_tib <- aqp_imp_tib %>% purrr::set_names(
        c("Imp_Cgn","Imp_Bld","Imp_Chr","Imp_Pos","Imp_Fwd","Imp_Top",
          "Imp_PrbT","Imp_PrbU","Imp_PrbM",
          "Imp_FR","Imp_TB","Imp_CO","Imp_Nxt",
          "Imp_ScrU","Imp_ScrM",
          "Imp_CpgCnt","Imp_CpgDis","Imp_Scr","Imp_NxbScr",
          "Imp_ScrMin","Imp_Inf") )
    }
  }
  
  # Qucik QC::
  #  aqp_imp_tib %>% dplyr::filter(Inf_Type != 0) %>% dplyr::distinct(Seq_ID)
  #  aqp_imp_tib %>% dplyr::group_by(Ord_Des,Ord_Din) %>% dplyr::summarise(Count=n(), .groups = "drop")
  #
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#                              END OF ROUND 1::
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# Load Genome::
#
chr_key <- "Bsp_Chr"
gen_ref_fas <- run$gen_ref_fas
nrec <- 1


dna_dat <- load_genome(file=gen_ref_fas, nrec=nrec,
                       chr_key=chr_key, ret_map=TRUE,
                       verbose=opt$verbose)

prb_seqs <- min_bsp_tib %>% 
  dplyr::mutate(Prb_Beg=dplyr::case_when(
    Bsp_FR=="F" & Bsp_CO=="O" & Ord_Des=="U" ~ Bsp_Pos + 0.0,
    Bsp_FR=="F" & Bsp_CO=="O" & Ord_Des=="M" ~ Bsp_Pos + 0.0,
    Bsp_FR=="F" & Bsp_CO=="O" & Ord_Des=="2" ~ Bsp_Pos + 1.0,
    TRUE ~ NA_real_ ) ) %>%
  dplyr::filter(Bsp_FR=="F") %>% 
  dplyr::filter(Bsp_CO=="O") %>% 
  parse_genomic_seqs(seq = as.character(dna_dat$seqs[1]), 
                     srd_str = "F",
                     pos_key = "Prb_Beg",
                     ups_len = 0,
                     seq_len = 50,
                     pad_len = 0,
                     chr_str = "chr1") %>% revCmp() 

dplyr::mutate(Prb_Seq=prb_seqs)

#
# - fas_to_seq() reference should be moved into bsp workflow
# - setTopBot_tib() on reference seqs 
# - validate TB calls against improbe
#

#
# - Select Probe Order Design from s_improbe
# - Select U/M to compare against improbe
#

#
#  N => r_improbe()
# !N => s_improbe()
#

sel_cols <- c("Aln_Key_Unq", 
              "Ord_Des","Ord_Din", 
              "Bsp_Chr","Bsp_Pos","Bsp_Cgn","Bsp_FR","Bsp_CO","Ord_Prb")

all_prb_tib <- NULL
min_bsp_tib <- aqp_cgn_tib %>% dplyr::select(dplyr::all_of(sel_cols))
all_gen_list <- all_gen_tib %>% split(f=all_gen_tib$Genome_Key)
#for (ii in c(1:base::nrow(all_gen_tib))) {
for (gen_key in names(all_gen_list)) {
  # gen_key  <- all_gen_tib[[]]$Genome_Key
  gen_fas  <- all_gen_list[[gen_key]]$Path
  gen_FR   <- all_gen_list[[gen_key]]$Strand_FR
  gen_CO   <- all_gen_list[[gen_key]]$Strand_CO
  gen_BSC  <- all_gen_list[[gen_key]]$Strand_BSC
  gen_prep <- all_gen_list[[gen_key]]$Alphabet
  
  if (gen_prep=="N") r_improbe_val <- TRUE
  if (gen_prep=="N") s_improbe_val <- FALSE
  if (gen_prep!="N") r_improbe_val <- FALSE
  if (gen_prep!="N") s_improbe_val <- TRUE
  
  cat(glue::glue("[{par$prgmTag}]: Genome({gen_key})={gen_fas}...{RET}"))
  
  gen_key <- "RCM_dna"
  # all_prb_tib[[gen_key]] <- fas_to_seq(tib=min_bsp_tib,
  RCM_dna_revCmp_tib <- fas_to_seq(tib=min_bsp_tib, srd_str = "R",
                                   gen_bld=paste(opt$genBuild,gen_key, sep="_"),
                                   
                                   gen_ref_fas=gen_fas,
                                   # imp_tsv=imp_inp_tsv,
                                   # seq_csv=imp_seq_csv,
                                   
                                   build=c("Prb1C","Prb2C","Prb1O","Prb2O"),
                                   
                                   ids_key="Aln_Key_Unq",
                                   din_key="Ord_Din",
                                   tar_din="rs",
                                   pos_key="Bsp_Pos",
                                   chr_key="Bsp_Chr",
                                   
                                   # subset=subset,
                                   # sub_cols=sub_cols,
                                   
                                   reload=FALSE,
                                   retData=FALSE,
                                   parallel=opt$parallel,
                                   r_improbe=r_improbe_val,
                                   s_improbe=s_improbe_val,
                                   add_flanks=TRUE,
                                   
                                   verbose=opt$verbose+1,tt=pTracker)
  
  cat(glue::glue("[{par$prgmTag}]: Genome({gen_key})={gen_fas}...{RET}{RET}"))
}




pairwiseAlignment(FCN_r_improbe_tib$Prb_1M %>% revCmp(), FCN_r_improbe_tib$DesBscM %>% stringr::str_to_upper(), type="local-global")

all_prb_tib$RCM_dna %>% dplyr::filter(Ord_Des=="M") %>% dplyr::select(Aln_Key_Unq, Bsp_FR, Bsp_CO, Ord_Prb, Ord_Prb, Prb1C, Temp, Iupac_Forward_Sequence) %>% head() %>% as.data.frame()
RCM_dna_revCmp_tib %>% dplyr::filter(Ord_Des=="M") %>% dplyr::select(Aln_Key_Unq, Bsp_FR, Bsp_CO, Ord_Prb, Ord_Prb, Prb1C, Temp, Iupac_Forward_Sequence) %>% head() %>% as.data.frame()

FCN_r_improbe_tib <- r_improbe(tib=aqp_imp_tib$s_ref, ids_key="Aln_Key_Unq", seq_key="Forward_Sequence", din_key="Ord_Din", srsplit = TRUE, srd_key = "Bsp_FR", cosplit = TRUE, cos_key = "Bsp_CO", parallel = TRUE, add_matseq = TRUE, verbose = opt$verbose, tt=pTracker)

print_prbs(tib = FCN_r_improbe_tib, 
           tar_des = "cg", 
           ids_key = "Aln_Key_Unq", 
           prb_key = "Ord_Prb",
           des_key = "Ord_Des", 
           din_key = "Ord_Din",
           outDir = file.path(opt$outDir),
           plotName = "FCN_dna",
           verbose = opt$verbose+10)

FCN_r_improbe_tib %>% head(n=1) %>%
  dplyr::rename(
    StrandFR=Strand_FR,
    StrandCO=Strand_CO,
  ) %>%
  prbsToStr(pr="cg", verbose = opt$verbose+5)

all_prb_tib$FCM_iupac %>% 
  cmpInfII_MisMatch(fieldA="Ord_Prb", 
                    fieldB="Prb1C", 
                    verbose=opt$verbose) %>% 
  dplyr::group_by(Ord_Des, Ord_Din, Bsp_FR, Bsp_CO, Man_MisMatch) %>%
  dplyr::summarise(Count=n(), .groups="drop")








aqp_top_tib <- setTopBot_tib(aqp_imp_tib$s_ref, 
                             seqKey="Forward_Sequence", 
                             srdKey="Bsp_TB2", 
                             verbose=opt$verbose,tt=pTracker)

aqp_top_tib <- setTopBot_tib(aqp_imp_tib$s_ref, 
                             seqKey="Forward_Sequence", 
                             srdKey="Bsp_TB2", 
                             verbose=opt$verbose,tt=pTracker)

imp_top_tib <- setTopBot_tib(aqp_imp_tib$i_imp, 
                             seqKey="Forward_Sequence_imp", 
                             srdKey="Bsp_TB2",
                             verbose=opt$verbose,tt=pTracker)

all_top_tib <- setTopBot_tib(aqp_imp_tib$i_imp,
                             seqKey="Top_Sequence", 
                             srdKey="Bsp_TB2",
                             verbose=opt$verbose,tt=pTracker)

#
# Design Steps::
#  gen_ref_fas
#    - Don't build probes (--build_prbs=FALSE)
#
#  r-improbe(Forward_Sequence=Imp_Temp_Seq)
#    - Compare Prb1_[U/M] against improbe [U/M]
#    - Compare by Ord_Des (2/U/M) against Ord_Prb
#
#  s-improbe (substring improbe) [BSC U/M/D]
#    - Compare Prb1_[U/M] against improbe [U/M]
#    - Compare by Ord_Des (2/U/M) against Ord_Prb
#
# SNP Check::
#    - Pos:Iupac (Include Next Base: pos=0)
#
# Update Probe_ID
#    - rs/ch database
#    - mu = multiple zero mismatch hits
#    - ma = multiple non-zero mismatch hits
#    - um = un-paired Infinium I probes
#
# Extension/Color Summary::
#    - Extension Summary (Cpg, Nxb)
#    - Color Summary (Red, Grn)
#
# Clean-Up Steps::
#  - Remove/Rename Temp_Seq
#  - Remove *_Len
#  - Only build probes when nescessary
#

tmp_join_vec <- c("Aln_Key_Unq", "Address","Ord_Des", "Ord_Din")

tmp_join_tib <- dplyr::inner_join(
  aqp_imp_list$fwd, aqp_imp_list$ret, 
  by=tmp_join_vec,
  suffix=c("_fwd","_imp") )






# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       1.0 Write improbe input::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: 
#  - Compare designs from r-improbe to substr-BSC-genomes!!!
#

nrec <- 2

fwd_seq_tsv <- file.path(par$topDir, "tmp/imp.fwd-seq.tsv.gz")
snp_seq_tsv <- file.path(par$topDir, "tmp/imp.snp-seq.tsv.gz")
fwd_des_tsv <- file.path(par$topDir, "tmp/test.improbe-designOutput.tsv.gz")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   0.0 Get Order Probes for Comparison::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ord_imp_tib <- aqp_imp_tib %>% 
  dplyr::mutate(Seq_ID=Aln_Key_Unq,
                Len=stringr::str_length(Seq_ID))
# Sanity Check length(Seq_ID) <= 15::
ord_imp_tib %>% dplyr::select(Seq_ID, Len) %>% dplyr::arrange(-Len)

ord_prb_tib <- ord_imp_tib %>% 
  dplyr::select(Seq_ID,Ord_Des,Ord_Din, 
                Bsp_FR,Strand_Ref_FR,Bsp_CO, Ord_Prb) %>%
  dplyr::rename(Strand_FR=Strand_Ref_FR,
                Strand_CO=Bsp_CO)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       1.0 Get Templates from Genome::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fwd_seq_tib <- ord_imp_tib %>% 
  dplyr::arrange(Bsp_Chr,Bsp_Pos) %>%
  fas_to_seq(fas=run$gen_ref_fas, file=fwd_seq_tsv, 
             name="Seq_ID", din="Ord_Din", 
             gen=opt$genBuild,
             chr1="Bsp_Chr", pos="Bsp_Pos",
             nrec=nrec, verbose=opt$verbose)

# snp_seq_org <- snp_seq_tib
# snp_seq_org2 <- snp_seq_tib

snp_seq_tib <- ord_imp_tib %>% 
  dplyr::arrange(Bsp_Chr,Bsp_Pos) %>%
  fas_to_seq(fas=run$gen_snp_fas, file=snp_seq_tsv, 
             name="Seq_ID", din="Ord_Din", 
             gen=opt$genBuild,
             chr1="Bsp_Chr", pos="Bsp_Pos",
             nrec=nrec, verbose=opt$verbose)

# Compare probes against designs::
#
fwd_seq_tib %>% dplyr::filter(Bsp_FR=="F", Bsp_CO=="C") %>% 
  dplyr::select(Aln_Key_Unq, Probe_Seq_T,Prb1_FC)
snp_seq_tib %>% dplyr::filter(Bsp_FR=="F", Bsp_CO=="C") %>% 
  dplyr::select(Aln_Key_Unq, Probe_Seq_T,Prb1_FC)

#
# Get Converted Genomes::
#

cur_gen_dir <- file.path(opt$genDir, opt$genBuild,"Sequence/WholeGenomeFasta")
ref_gen_pat <- paste0(opt$genBuild,".genome.[FR]C[MUD].fa.gz$")
bsc_gen_pat <- paste0(opt$genBuild,".genome.dbSNP-151.iupac.[FR]C[MUD].fa.gz$")

ref_bsc_files <- list.files(cur_gen_dir, pattern=ref_gen_pat, full.names=TRUE)
snp_bsc_files <- list.files(cur_gen_dir, pattern=bsc_gen_pat, full.names=TRUE)

#
# Validate Comparison of difference (i.e. contains SNPs)
#
if (FALSE) {
  fwdSnp_seq_tib <- dplyr::inner_join(fwd_seq_tib,snp_seq_tib, 
                                      by="Seq_ID",
                                      suffix=c("_fwd","_snp"))
  
  fwdSnp_seq_tib %>% 
    dplyr::filter(Fwd_Temp_Seq_fwd != Fwd_Temp_Seq_snp) %>% 
    dplyr::select(Fwd_Temp_Seq_fwd,Fwd_Temp_Seq_snp)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         2.0 improbe/docker::fwd
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ret_val <- 
  run_improbe_docker(
    file=fwd_seq_tsv, 
    name="test", image=image_str, shell=image_ssh,
    verbose=opt$verbose)

fwd_ides_tib <-
  load_improbe_design(
    file=fwd_des_tsv, join=NULL, out=NULL,
    level=3, add_inf=TRUE,
    verbose=opt$verbose)
fwd_iprb_tib <- fwd_ides_tib %>% 
  dplyr::select(Seq_ID,Strand_Ref_FR,Strand_TB,
                Strand_CO,Probe_Seq_U,Probe_Seq_M) %>% 
  dplyr::rename(Strand_FR=Strand_Ref_FR) %>%
  dplyr::left_join(dplyr::select(ord_imp_tib, Seq_ID,Ord_Des,Ord_Din),
                   by="Seq_ID", suffix=c("_fwd","_ord")) %>%
  dplyr::select(Seq_ID,Ord_Des,Ord_Din, dplyr::everything())

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           3.0 r-improbe::fwd
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fwd_rdes_tib <- desSeq_to_prbs(
  tib=fwd_seq_tib, 
  ids_key="Seq_ID", seq_key="Imp_Temp_Seq", prb_key="Ord_Din", 
  strsSR="FR", strsCO="CO", 
  parallel=TRUE, verbose=opt$verbose )
fwd_rprb_tib <- fwd_rdes_tib %>% 
  dplyr::select(Seq_ID, Strand_SR,Strand_CO, 
                PRB1_U_MAT,PRB1_M_MAT,PRB2_D_MAT) %>%
  dplyr::rename(Strand_FR=Strand_SR) %>%
  dplyr::left_join(dplyr::select(ord_imp_tib, Seq_ID,Ord_Des,Ord_Din),
                   by="Seq_ID", suffix=c("_fwd","_ord")) %>%
  dplyr::select(Seq_ID,Ord_Des,Ord_Din, dplyr::everything())


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           4.0 r-improbe::snp
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

snp_rdes_tib <- desSeq_to_prbs(
  tib=snp_seq_tib,
  ids_key="Seq_ID", seq_key="Iup_Temp_Seq", prb_key="Ord_Din", 
  strsSR="FR", strsCO="CO", 
  parallel=TRUE, verbose=opt$verbose )
snp_rprb_tib <- snp_rdes_tib %>% 
  dplyr::select(Seq_ID, Strand_SR,Strand_CO, 
                PRB1_U_MAT,PRB1_M_MAT,PRB2_D_MAT) %>%
  dplyr::rename(Strand_FR=Strand_SR) %>%
  dplyr::left_join(dplyr::select(ord_imp_tib, Seq_ID,Ord_Des,Ord_Din),
                   by="Seq_ID", suffix=c("_fwd","_ord")) %>%
  dplyr::select(Seq_ID,Ord_Des,Ord_Din, dplyr::everything())

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           4.0 Compare Results::
#
#  - fwd_iprb_tib vs. [fwd_rprb_tib, snp_rprb_tib]
#  - ord_prb_tib  vs. [fwd_iprb_tib, fwd_rprb_tib, snp_rprb_tib]
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Functionize Intersect/Comparison Methods::
#

bys_vec  <- c("Seq_ID","Ord_Des","Ord_Din","Strand_FR","Strand_CO")
grp_vec  <- c("Ord_Des","Ord_Din","Strand_FR","Strand_CO")

ref_keys <- c("Probe_Seq_U","Probe_Seq_M")
can_keys <- c("PRB1_U_MAT","PRB1_M_MAT")

ir_fwd_cmp <- compare_probes(fwd_iprb_tib,fwd_rprb_tib, 
                             ref_keys=ref_keys, can_keys=can_keys,
                             by=bys_vec, grp=grp_vec,
                             retData=TRUE,
                             verbose=opt$verbose)

ref_keys <- c("PRB1_U_MAT","PRB1_M_MAT","PRB2_D_MAT")
can_keys <- c("PRB1_U_MAT","PRB1_M_MAT","PRB2_D_MAT")

ii_fwdSnp_cmp <- compare_probes(fwd_rprb_tib,snp_rprb_tib,
                                ref_keys=ref_keys, can_keys=can_keys,
                                by=bys_vec, grp=grp_vec,
                                retData=TRUE, pivot=TRUE,
                                verbose=opt$verbose)

ref_keys <- c("Ord_Prb")
can_keys <- c("Probe_Seq")

snp_rprb_tab <- dplyr::bind_rows(
  dplyr::filter(snp_rprb_tib,Ord_Des=="U") %>% dplyr::select(-PRB1_M_MAT,-PRB2_D_MAT) %>% dplyr::rename(Probe_Seq=PRB1_U_MAT),
  dplyr::filter(snp_rprb_tib,Ord_Des=="M") %>% dplyr::select(-PRB1_U_MAT,-PRB2_D_MAT) %>% dplyr::rename(Probe_Seq=PRB1_M_MAT),
  dplyr::filter(snp_rprb_tib,Ord_Des=="2") %>% dplyr::select(-PRB1_U_MAT,-PRB1_M_MAT) %>% dplyr::rename(Probe_Seq=PRB2_D_MAT)
)

oi_ordSnp_cmp <- compare_probes(ord_prb_tib,snp_rprb_tab,
                                ref_keys=ref_keys, can_keys=can_keys,
                                by=bys_vec, grp=grp_vec,
                                retData=TRUE, pivot=TRUE,
                                verbose=opt$verbose)



snp_rprb_tib_lc <- snp_rprb_tib %>% dplyr::mutate(PRB1_U_MAT=stringr::str_to_lower(PRB1_U_MAT),PRB1_M_MAT=stringr::str_to_lower(PRB1_M_MAT),PRB2_D_MAT=stringr::str_to_lower(PRB2_D_MAT))

mapply(list.string.diff, fwd_rprb_tib$PRB1_U_MAT,snp_rprb_tib$PRB1_U_MAT)

mapply(list.string.diff, fwd_rprb_tib$PRB1_U_MAT,snp_rprb_tib$PRB1_U_MAT)

# Pivot for human viewing::
ii_fwdSnp_cmp$cmp %>% 
  # dplyr::filter(Ord_Des=="M" & Man_MisMatch==3) %>% head() %>% 
  tidyr::pivot_longer(cols=c("Probe_A","Probe_B"), 
                      names_to="Prb_Source", values_to="Probe_Seq") %>% 
  as.data.frame()





# r-improbe:: improbe_design_all()
#

#
# TBD::
#
#  Manifest Generation::
#    - Code Clean Up
#       - [done] bsp_mapping_workflow()
#
#    - Calculate extension/color distribution
#    - Extract BSC Top/Probe-Design
#
#    - Join by position
#    - Add masked-controls from source
#
#  Annotation::
#    - Incorporate
#    - Annotation Summary
#    - Implement ftp for missing files
#
#  Cluster::
#    - Transfer all files to cluster
#



# TBD::
#  - Try to join Singles::
#  - Add Singles
#
opt$multi_unique <- TRUE
opt$multi_unique <- FALSE

if (opt$multi_unique) {
  aqp_imp_tib <- aqp_imp_tib %>% 
    dplyr::distinct(Address,Bsp_Chr,Bsp_Pos, .keep_all=TRUE) %>%
    dplyr::add_count(Address, name="Bsp_Hits") %>%
    dplyr::mutate(
      Cgn_New=paste0(Ord_Din,stringr::str_remove(Cgn_Str,"^[^0-9]+")),
      Cgn_Str=dplyr::case_when(
        Bsp_Hits>1 ~ Cgn_New,
        TRUE ~ Cgn_Str
      )
    )
}

#
# Fix rs/ch probes::
#
aqp_imp_tib <- aqp_imp_tib %>% 
  dplyr::mutate(
    Cgn_Str=dplyr::case_when(
      Ord_Din=="rs" | Ord_Din=="ch" ~ Ord_Key,
      TRUE ~ Cgn_Str
    )
  )

aqp_imp_tib %>%
  dplyr::group_by(Ord_Des,Ord_Din) %>% 
  dplyr::summarise(Count=n(), .groups = "drop") %>%
  print()

imp_des_list <- aqp_imp_tib %>% split(.$Ord_Des)

inf1_tib <- dplyr::inner_join(imp_des_list[["U"]],imp_des_list[["M"]],
                              by=c("Ord_Key","Cgn","Ord_Din",
                                   "Cgn_Str", "Bsp_Chr","Bsp_Pos",
                                   "Top_Sequence","Forward_Sequence",
                                   "Cpg_Cnt",
                                   "Genome_Build"),
                              suffix=c("_U","_M") )
inf1_tib %>%
  dplyr::group_by(Ord_Des_U,Ord_Des_M,Ord_Din) %>% 
  dplyr::summarise(Count=n(), .groups = "drop") %>%
  print()

inf2_tib <- imp_des_list[["2"]]
inf2_tib %>%
  dplyr::group_by(Ord_Des,Ord_Din) %>% 
  dplyr::summarise(Count=n(), .groups = "drop") %>%
  print()

#
# TBD:: aqp_seq_tib needs to split Infinium II probes into U/M by degenerate
#   bases to distinguish alignments...
#

prb1_dat <- aqp_to_sesame1(inf1_tib, isMU=TRUE, retData=TRUE, verbose=opt$verbose)
prb2_dat <- aqp_to_sesame2(inf2_tib, isMU=TRUE, retData=TRUE, verbose=opt$verbose)

prb2_dat$ses_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "rs"))
prb2_dat$ses_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "ch"))

#
#
# New Merging::
#  - ses_man_tib <- bind(sesN_tib's)
#  - aux_pos_tib <- bind(allN_tib's) => extract coordinates
#
#

# aux_pos_tib <- dplyr::bind_rows(all1_tib, all2_tib) %>% 
aux_pos_tib <- dplyr::bind_rows(prb1_dat$all_tib, prb2_dat$all_tib) %>% 
  dplyr::arrange(CHR,MAPINFO) %>% 
  dplyr::select(Probe_ID, CHR,MAPINFO, Strand_FR,Strand_CO, Species,Genome_Build)

aux_pos_name <- paste0(par$local_runType,"-",opt$version,".alignments.csv.gz")
aux_pos_csv  <- file.path(run$man_dir, aux_pos_name)
safe_write(x=aux_pos_tib,file=aux_pos_csv, verbose=opt$verbose)

# dplyr::select(Probe_ID,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
#               Next_Base,Color_Channel,Col,Probe_Type,
#               Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,Rep_Num,
#               CHR,MAPINFO,Species,Genome_Build,
#               Source_Seq,Forward_Sequence,Top_Sequence,
#               Underlying_CpG_Count,Underlying_CpG_Min_Dist,
#               # Alt_Cgn_Count_U,Alt_Cgn_Count_M, 
#               Bsp_Tag_U,Bsp_Tag_M)

#
# Need to reduce Sesame manifest and move alignments to aux files...
#
ses_man_tib <- 
  dplyr::bind_rows(prb1_dat$ses_tib,prb2_dat$ses_tib) %>% 
  dplyr::arrange(CHR,MAPINFO) %>% 
  dplyr::select(Probe_ID,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                Next_Base,Color_Channel,Col,Probe_Type,
                Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,Rep_Num,
                CHR,MAPINFO,Species,Genome_Build,
                Source_Seq,Forward_Sequence,Top_Sequence,
                Underlying_CpG_Count,Underlying_CpG_Min_Dist,
                # Alt_Cgn_Count_U,Alt_Cgn_Count_M, 
                Bsp_Tag_U,Bsp_Tag_M)

ses_man_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "rs"))
ses_man_tib %>% dplyr::filter(stringr::str_starts(Probe_ID, "ch"))

#
# So many singletons::
#
if (FALSE) {
  aqp_imp_tib %>% 
    dplyr::anti_join(ses_man_tib, by=c("Address"="U")) %>% 
    dplyr::anti_join(ses_man_tib %>% dplyr::filter(!is.na(M)), 
                     by=c("Address"="M")) %>% 
    dplyr::group_by(Ord_Des,Ord_Din) %>% 
    dplyr::summarise(Count=n(), .groups = "drop")
  #
  # Make Them ALL Infinium I Probes::
  #
  prb0_dat <- aqp_imp_tib %>% 
    dplyr::anti_join(ses_man_tib, by=c("Address"="U")) %>% 
    dplyr::anti_join(ses_man_tib %>% dplyr::filter(!is.na(M)), 
                     by=c("Address"="M")) %>% 
    aqp_to_sesame(isMU=TRUE, retData=TRUE, verbose=opt$verbose)
}

#
# Add Controls and Write Manifest::
#
ses_man_tib <- ses_man_tib %>% 
  # dplyr::bind_rows(prb0_dat$ses_tib) %>%
  dplyr::rename(col=Col) %>%
  dplyr::bind_rows(dplyr::rename(ses_ctl_tib,Color_Channel=COLOR_CHANNEL)) %>%
  clean_tibble()

ses_man_name <- paste0(par$local_runType,"-",opt$version,".manifest.sesame-base.cpg-sorted.csv.gz")
ses_man_name <- paste0(par$local_runType,"-",opt$version,".manifest.sesame-base.11k-controls.cpg-sorted.csv.gz")
ses_man_csv  <- file.path(run$man_dir, ses_man_name)
safe_write(x=ses_man_tib,file=ses_man_csv, verbose=opt$verbose)

if (FALSE) {
  ses_red_tib <- ses_man_tib %>%
    dplyr::select(Probe_ID,Name,U,AlleleA_ProbeSeq,M,AlleleB_ProbeSeq,
                  Next_Base,Color_Channel,col,Probe_Type,
                  Strand_FR,Strand_TB,Strand_CO,Infinium_Design_Type,Rep_Num,
                  CHR,MAPINFO,Species,Genome_Build,
                  Source_Seq,Forward_Sequence,Top_Sequence)
  
  ses_man_name <- paste0(par$local_runType,"-",opt$version,".manifest.sesame-base.cpg-sorted.csv.gz")
  ses_man_name <- paste0(par$local_runType,"-",opt$version,".manifest.sesame-base.11k-controls.cpg-sorted.csv.gz")
  ses_man_csv  <- file.path(run$man_dir, ses_man_name)
  safe_write(x=ses_red_tib,file=ses_man_csv, verbose=opt$verbose)
}





#
# Validate Against Order File::
#
old_man_csv <- "/Users/bretbarnes/Documents/data/manifests/methylation/Chicago-Ober-Custom.original/Chicago-S38.manifest.sesame-base.cpg-sorted.csv"
old_man_tib <- safe_read(old_man_csv)

ses_man_tib %>% dplyr::filter(!is.na(U) & ! U %in% old_man_tib$U) %>% as.data.frame()
ses_man_tib %>% dplyr::filter(!is.na(M) & ! M %in% old_man_tib$M) %>% as.data.frame()

ses_mis_tib <- dplyr::bind_rows(
  old_man_tib %>% dplyr::filter(!is.na(U) & ! U %in% ses_man_tib$U),
  old_man_tib %>% dplyr::filter(!is.na(M) & ! M %in% ses_man_tib$M)
) %>% dplyr::distinct() %>% 
  dplyr::mutate(
    Infinium_Design_Type=dplyr::case_when(
      Infinium_Design_Type==1 ~ "I",
      Infinium_Design_Type==2 ~ "II",
      TRUE ~ NA_character_)
  ) %>% clean_tibble()


#
# Validation that all the unpaired probes have failing partners!!!
#
if (FALSE) {
  # These seem to be singletons::
  #
  misM_tib <- aqp_add_tib %>% 
    dplyr::filter(Ord_Des=="M") %>% 
    dplyr::filter(! Address %in% ses_man_tib$M)
  misU_tib <- aqp_add_tib %>% 
    dplyr::filter(Ord_Des!="M") %>% 
    dplyr::filter(! Address %in% ses_man_tib$U)
  
  # No Overlap::
  #
  misM_tib %>% dplyr::inner_join(misU_tib, by="Ord_Cgn")
  
  mat_vec  <- splitStrToVec(opt$mats)
  mats_tib <- dplyr::bind_rows(load_aqp_files(mat_vec[1]),load_aqp_files(mat_vec[2]))
  
  mat_mis1_tib <- dplyr::bind_rows(
    mats_tib %>% dplyr::filter(Mat_Prb %in% misM_tib$Ord_Prb),
    mats_tib %>% dplyr::filter(Mat_Prb %in% misU_tib$Ord_Prb)
  )
  pqc_tib %>% dplyr::filter(Address %in% mat_mis1_tib$Address) %>% print(n=10000)
  
  # The Pairs all failed!!!
  mat_mis2_tib <- dplyr::bind_rows(
    mats_tib %>% dplyr::filter(Mat_Prb %in% misM_tib$Ord_Par),
    mats_tib %>% dplyr::filter(Mat_Prb %in% misU_tib$Ord_Par)
  )
  pqc_tib %>% dplyr::filter(Address %in% mat_mis2_tib$Address) %>% print(n=10000)
  
  pass_mis_csv <- file.path(par$topDir, "data/CustomContent/McMaster/manifest/passing-but-missing.10072021.txt")
  pass_mis_tib <- readr::read_csv(pass_mis_csv)
  
  mis1_tib <- aqp_add_tib %>% 
    dplyr::mutate(Ord_Key2=stringr::str_remove(Ord_Key,"[-_:].*$")) %>% 
    dplyr::filter(Ord_Key2 %in% pass_mis_tib$Probe_ID) %>% 
    dplyr::distinct(Ord_Key, .keep_all = TRUE)
  
  ord_tib <- load_aqp_file( opt$ords )
  mis2_tib <- ord_tib %>% 
    dplyr::filter(Ord_Key %in% pass_mis_tib$Probe_ID) %>% 
    dplyr::distinct(Ord_Key, .keep_all = TRUE)
  
  mis2_tib %>% dplyr::filter(Ord_Key %in% mis1_tib$Ord_Key)
  mis1_tib %>% dplyr::filter(Ord_Key %in% mis2_tib$Ord_Key)
  
  aqp_add_tib %>% dplyr::filter(Ord_Prb %in% mis2_tib$Ord_Prb)
  
  mis2_tib %>% dplyr::inner_join(mats_tib, by=c("Ord_Prb"="Mat_Prb"))
  
  load_aqp_files(mat_vec)
  
}


#
# Mapping to original coordinates with Genomic Regions::
#
imp_trim_tib <- NULL
if (FALSE) {
  
  pos_cols <- 
    cols(
      Cgn    = col_integer(),
      Chr    = col_character(),
      Pos    = col_integer(),
      PrbA   = col_character(),
      ExtA   = col_character(),
      Srd_FR = col_character(),
      Srd_TB = col_character(),
      Srd_CO = col_character(),
      Nxb    = col_character(),
      PrbB   = col_character(),
      ExtB   = col_character()
    )
  
  # ord_pos_tib <- safe_read(par$ord_pos_csv, verbose=opt$verbose)
  ord_pos_tib <- readr::read_tsv(par$ord_pos_csv, 
                                 col_names = names(pos_cols$cols), 
                                 col_types=pos_cols) %>%
    dplyr::mutate(Chr=paste0("chr",Chr)) %>% 
    dplyr::distinct(Chr,Pos, .keep_all = TRUE) %>%
    dplyr::mutate(Rank=dplyr::row_number(),
                  Unq_Key=paste(Chr,Pos,Rank,sep="_"))
  
  buf_len <- 50
  ord_grs =
    GenomicRanges::GRanges(
      seqnames=Rle(ord_pos_tib$Chr),
      # strand=Rle(ret_tib$srd),
      
      PrbA=ord_pos_tib$PrbA,
      
      IRanges(start=ord_pos_tib$Pos-buf_len,
              end=ord_pos_tib$Pos+buf_len,
              names=ord_pos_tib$Unq_Key)
    )
  
  imp_des_list <- aqp_imp_tib %>% split(.$Ord_Des)
  
  imp_grs =
    GenomicRanges::GRanges(
      seqnames=Rle(aqp_imp_tib$Bsp_Chr),
      # strand=Rle(ret_tib$srd),
      
      Cgn_Str=aqp_imp_tib$Cgn_Str,
      Aln_Key=aqp_imp_tib$Aln_Key,
      Ord_Key=aqp_imp_tib$Ord_Key,
      Cpg_Pos=aqp_imp_tib$Bsp_Pos,
      
      IRanges(start=aqp_imp_tib$Bsp_Pos,
              width = 2,
              names=aqp_imp_tib$Aln_Key_Unq)
    )
  
  imp_int_tib <- intersect_GRS(ord_grs,imp_grs, 
                               can_prefix="Imp", 
                               can_key="Unq_Key", 
                               verbose=opt$verbose)
  
  imp_trim_tib <- aqp_imp_tib %>% dplyr::filter(Aln_Key %in% imp_int_tib$Imp_Aln_Key)
  
  imp_des_list <- NULL
  if (!is.null(imp_trim_tib)) {
    imp_des_list <- imp_trim_tib %>% split(.$Ord_Des)
  } else {
    imp_des_list <- aqp_imp_tib %>% split(.$Ord_Des)
  }
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#                              END OF ROUND 2::
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Current thoughts are true joining:: Two ways::
#   A: by coordinate and orientation
#      - TBD: Should add coordinate cgn look up!!
#      - TBD: Use the extension distribution some how
#   C: by cgn in case alignment failed
#
# This should all be tested with a small multi-unique set like NZT!
#
# Both bsp_cgn_tib and aqp_bsp_tib need to be joined seperately
#  with aqp_add_tib
#
if (FALSE) {
  #   A: by coordinate and orientation
  aqp_bsp_tibA <- aqp_bsp_tib %>% 
    dplyr::select(Address,Bsp_Chr,Bsp_Pos,Bsp_Din_Bsc,Bsp_Nxb_Bsc,
                  Bsp_Tag,Bsp_Srd,Ord_Key) %>%
    dplyr::arrange(Bsp_Chr,Bsp_Pos)
  
  aqp_bsp_tibB <- aqp_bsp_tib %>% 
    dplyr::select(Address,Bsp_Chr,Bsp_Pos,Bsp_Din_Bsc,Bsp_Nxb_Bsc,
                  Bsp_Tag,Bsp_Srd,Ord_Key) %>%
    dplyr::arrange(Bsp_Chr,Bsp_Pos)
  
  aqp_bsp_inn <- aqp_bsp_tibA %>% 
    dplyr::left_join(aqp_bsp_tibB, by=c("Bsp_Chr","Bsp_Pos"), suffix=c("_A", "_B")) %>% 
    dplyr::distinct() %>%
    dplyr::filter(Address_A!=Address_B)
  
  
  aqp_bsp_pas <- dplyr::bind_rows(
    aqp_bsp_inn %>% dplyr::filter(Ord_Key_A==Ord_Key_B),
    aqp_bsp_inn %>% dplyr::filter(Ord_Key_A!=Ord_Key_B)
  )
  aqp_bsp_bad <- dplyr::anti_join(
    aqp_bsp_inn,
    aqp_bsp_pas, 
    by=c("Address_A","Address_B")
  )
  
  # Reload the original order file::
  ord_tib <- load_aqp_file(opt$ords)
  mat_tib <- load_aqp_file(opt$mats)
  aqp_tib <- load_aqp_file(opt$aqps)
  aqp_tib %>% dplyr::group_by(Decode_Status) %>% dplyr::summarise(Count=n(), .groups="drop")
  
  bsp_tsv <- file.path(par$topDir, "scratch/workhorse_main_dev_latest/NZT-NZT-C0-GRCh37/aln/NZT-NZT-C0-GRCh37.aqp-pass.address.bsp.tsv.gz")
  bsp_tib <- load_bsmap(bsp_tsv, verbose = opt$verbose)
  
  bed_cols <- 
    cols(
      Imp_Chr = col_character(),
      Imp_Pos = col_integer(),
      imp_End = col_integer(),
      Imp_Cgn = col_integer(),
      Imp_Bld = col_character(),
      Imp_Srd = col_character()
    )
  
  imp_bed <- file.path(par$topDir, "data/improbe/scratch/cgnDB/dbSNP_Core4/design-input/GRCh37.cgn.bed.gz")
  imp_tib <- readr::read_tsv(imp_bed, col_names=names(bed_cols), col_types=bed_cols)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        6.0 Annotate Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$genBuild=="GRCh37" || opt$genBuild=="GRCh38") {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #              6.1 Load Annotations:: EPIC_CORE/UPDATE_CORE/CHROM_HMM
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Could loop over everything and make subdirectories for EPIC,UPDATE,HMM
  #
  # epic_ann_file <- list.dirs(epic_ann_path, full.names = TRUE)[-1] 
  # epic_dir_list <- as.list(epic_ann_file)
  # names(epic_dir_list) <- base::basename(epic_ann_file)
  
  core_ann_path <- file.path(par$topDir, "scratch/annotation_to_workhorse_bed/methylation-Human-GRCh37-v2")
  epic_ann_path <- file.path(core_ann_path, "EPIC_CORE/UCSC")
  epic_ann_file <- list.files(epic_ann_path,pattern=".bed.gz$",full.names=TRUE)
  epic_fns_list <- as.list(epic_ann_file)
  names(epic_fns_list) <- base::basename(epic_ann_file) %>% stringr::str_remove(".bed.gz")
  
  
  # TBD:: Load GRS all from list of files with lapply...
  epic_grs_list <- lapply(epic_fns_list, ann_to_grs, verbose=opt$verbose,tt=pTracker)
  
  can_key <- "IlmnID"
  epic_int_list <- NULL
  epic_int_list <- c(epic_int_list,
                     lapply(epic_grs_list, intersect_GRS, can=man_pos_grs, 
                            ref_key=NULL,ref_col=NULL,ref_prefix=NULL,ref_red=TRUE,
                            can_key=can_key,can_col=can_key,can_prefix=NULL, 
                            verbose=opt$verbose, tt=pTracker)
  )
  
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       6.0 Annotation Conformation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
# if (opt$genBuild=="GRCh37" || opt$genBuild=="GRCh38") {
if (FALSE) {
  
  core_anno_dir <- file.path(par$topDir, "data/annotation", opt$genBuild)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       6.2 Load Annotation:: EPIC_CORE
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # epic_ann_trim <- c(".bed.gz$",".sorted$",".intersect",".sorted",".formatted")
  # epic_ann_trim <- c(".sorted.merged.bed.gz",
  #                    ".sorted.intersect.bed.gz",
  #                    ".intersect.bed.gz",
  #                    ".formatted.sorted.bed.gz")
  
  epic_ann_trim <- c(".bed.gz$")
  epic_ann_path <- file.path(core_anno_dir, "EPIC_CORE")
  epic_ann_file <- list.dirs(epic_ann_path, full.names = TRUE)[-1] 
  epic_dir_list <- as.list(epic_ann_file)
  names(epic_dir_list) <- base::basename(epic_ann_file)
  
  epic_int_list <- NULL
  for (source in names(epic_dir_list)) {
    file_list <- get_file_list(
      dir=epic_dir_list[[source]], pattern=epic_ann_trim[1],
      trim=epic_ann_trim,
      verbose=opt$verbose)
    
    epic_out_path <- file.path(core_anno_dir, "EPIC_CORE_CLEAN")
    epic_tib_list <- lapply(file_list, grs=FALSE, load_epic_anno, source=source,
                            out=epic_out_path,
                            verbose=opt$verbose)
    
    next
    if (FALSE) {
      epic_grs_list <- lapply(file_list, grs=TRUE, load_epic_anno, source=source, 
                              verbose=opt$verbose)
      epic_int_list <- c(epic_int_list,
                         lapply(epic_grs_list, intersect_GRS, can=man_pos_grs, 
                                can_key="IlmnID", ref_prefix=NULL,
                                verbose=opt$verbose, tt=pTracker)
      )
    }
  }
  
  # Functionalize: Chrom_HMM and EPIC Annotation lists
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       6.3 Load Annotation:: NCBI/UCSC
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ncib_gene_tsv <- file.path(core_anno_dir, "NCBI", paste(opt$genBuild,"ncbi.RefSeqGenes.tsv.gz", sep='.'))
  ucsc_gene_tsv <- file.path(core_anno_dir, "UCSC", paste(opt$genBuild,"ucsc.knownGene.tsv.gz", sep='.'))
  ucsc_cpgs_tsv <- file.path(core_anno_dir, "UCSC", paste(opt$genBuild,"ucsc.CpG-Islands.tsv.gz", sep='.'))
  
  ncbi_gene_tib <- load_ncbi_gene(file=ncib_gene_tsv, verbose=opt$verbose, tt=pTracker)
  ucsc_gene_tib <- load_ucsc_gene(file=ucsc_gene_tsv, verbose=opt$verbose, tt=pTracker)
  ucsc_cpgs_tib <- load_ucsc_cpgs(file=ucsc_cpgs_tsv, verbose=opt$verbose, tt=pTracker)
  
  ncbi_gene_grs <- load_ncbi_gene(file=ncib_gene_tsv, grs=TRUE, verbose=opt$verbose, tt=pTracker)
  ucsc_gene_grs <- load_ucsc_gene(file=ucsc_gene_tsv, grs=TRUE, verbose=opt$verbose, tt=pTracker)
  ucsc_cpgs_grs <- load_ucsc_cpgs(file=ucsc_cpgs_tsv, grs=TRUE, verbose=opt$verbose, tt=pTracker)
  
  ann_int_list <- NULL
  # NCBI Gene Comparison::
  #
  ref_pre_str <- "NCBI_Gene"
  ann_key_str <- ref_pre_str
  ref_pre_str <- NULL
  # ncbi_gene_int_tib <- 
  ann_int_list[[ann_key_str]] <-
    intersect_GRS(can=man_pos_grs, ref=ncbi_gene_grs, 
                  can_key="IlmnID", ref_prefix=ref_pre_str, 
                  verbose=opt$verbose, tt=pTracker)
  
  if (FALSE) {
    ncbi_gene_crs_tib <- 
      ann_int_list[[ann_key_str]] %>% 
      dplyr::inner_join(ses_ann_tab, by="IlmnID")
    ncbi_gene_crs_sum <- ncbi_gene_crs_tib %>% 
      dplyr::group_by(class,Value_Str) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% 
      dplyr::arrange(-Count)
    ncbi_gene_crs_sum %>% print(n=base::nrow(ncbi_gene_crs_sum))
  }
  
  # UCSC Gene Comparison::
  #
  ref_pre_str <- "UCSC_Gene"
  ann_key_str <- ref_pre_str
  ref_pre_str <- NULL
  # ucsc_gene_int_tib <- 
  ann_int_list[[ann_key_str]] <-
    intersect_GRS(can=man_pos_grs, ref=ucsc_gene_grs, 
                  can_key="IlmnID", ref_prefix=ref_pre_str, 
                  verbose=opt$verbose, tt=pTracker)
  
  if (FALSE) {
    ucsc_gene_crs_tib <- 
      ann_int_list[[ann_key_str]] %>% 
      dplyr::inner_join(ses_ann_tab, by="IlmnID")
    ucsc_gene_crs_sum <- ucsc_gene_crs_tib %>% 
      dplyr::group_by(class,Value_Str) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% 
      dplyr::arrange(-Count)
    ucsc_gene_crs_sum %>% print(n=base::nrow(ucsc_gene_crs_sum))
  }
  
  # UCSC Islands Comparison::
  #
  ref_pre_str <- "UCSC_Islands"
  ann_key_str <- ref_pre_str
  ref_pre_str <- NULL
  # ucsc_cpgs_int_tib <- 
  ann_int_list[[ann_key_str]] <-
    intersect_GRS(can=man_pos_grs, ref=ucsc_cpgs_grs, 
                  can_key="IlmnID", ref_prefix=ref_pre_str, 
                  verbose=opt$verbose, tt=pTracker)
  
  if (FALSE) {
    ucsc_cpgs_crs_tib <- 
      ann_int_list[[ann_key_str]] %>% 
      dplyr::inner_join(ses_ann_tab, by="IlmnID")
    
    ucsc_cpgs_crs_sum <- ucsc_cpgs_crs_tib %>% 
      dplyr::group_by(class,Value_Str) %>% 
      dplyr::summarise(Count=n(), .groups="drop") %>% 
      dplyr::arrange(-Count)
    ucsc_cpgs_crs_sum %>% print(n=base::nrow(ucsc_cpgs_crs_sum))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       6.4 Load Annotation:: Chrom HMM
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  hmm_ann_drn <- paste("GRCh36",opt$genBuild, sep='-')
  hmm_ann_dir <- file.path(par$topDir, "data/annotation/liftOver/chrom_hmm/wgEncodeBroadHmm/ucsc_liftover_main",hmm_ann_drn)
  hmm_ann_fns <- list.files(hmm_ann_dir, pattern=".map.bed.gz$", full.names=TRUE, recursive=FALSE)
  
  hmm_cols <- 
    cols(
      chr    = col_character(),
      beg    = col_integer(),
      end    = col_integer(),
      class  = col_character(),
      val1   = col_integer(),
      val2   = col_character(),
      val3   = col_integer(),
      val4   = col_integer(),
      val5   = col_character()
    )
  
  hmm_sufix <- ".map.bed.gz"
  hmm_names <- hmm_ann_fns %>% base::basename() %>% 
    stringr::str_remove(hmm_sufix) %>% 
    stringr::str_remove("-.*$") %>% 
    stringr::str_remove("^wgEncodeBroadHmm") %>% 
    stringr::str_remove("HMM$")
  
  hmm_list <- NULL
  for (ii in c(1:length(hmm_ann_fns))) {
    name <- hmm_names[ii]
    hmm_list[[name]] <- hmm_ann_fns[ii]
  }
  
  # lapply(hmm_ann_fns, suppressMessages(suppressWarnings(readr::read_tsv)))
  hmm_dat_list <- lapply(hmm_list, 
                         readr::read_tsv, 
                         col_names=names(hmm_cols$cols), 
                         col_types=hmm_cols)
  
  for (samp in base::names(hmm_dat_list)) {
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]: Sample={samp}...{RET}"))
    
    cur_key <- paste(samp,"hmm", sep="_")
    cur_grp <- paste(cur_key,"class", sep="_")
    cur_grp_sym <- rlang::sym(cur_grp)
    cur_name <- paste(opt)
    
    cur_tib <- hmm_dat_list[[samp]] %>%
      dplyr::arrange(chr,beg) %>% 
      dplyr::group_by(class) %>% 
      dplyr::mutate(class=stringr::str_replace_all(class, '_','-'), 
                    Class_Rank=dplyr::row_number(), 
                    Uniq_Id=paste(class,Class_Rank, sep="_")) %>%
      dplyr::ungroup()
    
    cur_grs <- 
      GenomicRanges::GRanges(
        seqnames=Rle(cur_tib$chr),
        # strand=Rle(ses_pos_tib$srd),
        # Probe_Type=ses_pos_tib$Probe_Type,
        name=paste(cur_tib$chr,cur_tib$beg,cur_tib$end, sep='-'),
        name2=NA_character_,
        class=cur_tib$class,
        source="Chrom-HMM",
        tissue=samp,
        rank=cur_tib$Class_Rank,
        
        IRanges(start=cur_tib$beg,
                end=cur_tib$end,
                names=cur_tib$Uniq_Id)
      )
    
    if (opt$verbose>0) {
      cat(glue::glue("[{par$prgmTag}]: Sample={samp}; cur_grs={RET}"))
      print(cur_grs)
    }
    
    # Comparison::
    #
    ref_pre_str <- cur_key
    ann_key_str <- cur_key
    ref_pre_str <- NULL
    
    # cur_int_tib <- 
    ann_int_list[[ann_key_str]] <-
      intersect_GRS(can=man_pos_grs, ref=cur_grs,
                    can_key="IlmnID", ref_prefix=ref_pre_str, 
                    verbose=opt$verbose, tt=pTracker)
    
    if (FALSE) {
      
      cur_crs_tib <- 
        ann_int_list[[ann_key_str]] %>% 
        dplyr::inner_join(ses_ann_tab, by="IlmnID")
      cur_crs_sum <- cur_crs_tib %>% 
        # dplyr::group_by(!!cur_grp_sym,Value_Str) %>% 
        dplyr::group_by(class,Value_Str) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>% 
        dplyr::arrange(-Count)
      cur_crs_sum %>% print(n=base::nrow(cur_crs_sum))
      
    }
    
    if (opt$verbose>0)
      cat(glue::glue("[{par$prgmTag}]: Done. Sample={samp}.{RET}{RET}"))
    
    # break
  }
}











#
#  Validation of Manifest/BSP Coordinates via Sesame
#
par$validateSesame <- FALSE
if (par$validateSesame) {
  
  #
  # Load Sesame equivlent and compare coordinates and names by sequence
  #
  par$sesBuild <- NULL
  if (opt$genBuild=="GRCh37") par$sesBuild <- "hg19"
  if (opt$genBuild=="GRCh38") par$sesBuild <- "hg18"
  ses_man_grs <- sesameData::sesameDataGet(paste("EPIC",par$sesBuild,"manifest", sep='.'))
  ses_man_tib <- ses_man_grs %>% as.data.frame() %>% 
    tibble::rownames_to_column(var="Ses_Cgn") %>%
    tibble::as_tibble()
  
  #
  # Compare aqp_bsp_tib vs. ses_man_tib
  #
  ses_pos_col <- c("Ses_Cgn", "seqnames", "start", "end", "probeBeg", "probeEnd", "strand", 
                   "designType", "probeType")
  ses_pos_colA <- c(ses_pos_col,"ProbeSeq_A")
  ses_pos_colB <- c(ses_pos_col,"ProbeSeq_B")
  
  ses_pos_tab <- dplyr::bind_rows(
    ses_man_tib %>% 
      dplyr::filter(designType=="I") %>%
      dplyr::select(dplyr::all_of(ses_pos_colA)) %>% 
      dplyr::rename(Ord_Prb=ProbeSeq_A) %>%
      dplyr::mutate(Ord_Prb=stringr::str_to_upper(Ord_Prb),
                    Ord_Des="U") %>%
      tibble::as_tibble(),
    ses_man_tib %>% 
      dplyr::filter(designType=="I") %>%
      dplyr::select(dplyr::all_of(ses_pos_colB)) %>% 
      dplyr::rename(Ord_Prb=ProbeSeq_B) %>%
      dplyr::mutate(Ord_Prb=stringr::str_to_upper(Ord_Prb),
                    Ord_Des="M") %>%
      tibble::as_tibble(),
    ses_man_tib %>% 
      dplyr::filter(designType=="II") %>%
      dplyr::select(dplyr::all_of(ses_pos_colA)) %>% 
      dplyr::rename(Ord_Prb=ProbeSeq_A) %>%
      dplyr::mutate(Ord_Prb=stringr::str_to_upper(Ord_Prb),
                    Ord_Des="2") %>%
      tibble::as_tibble()
  ) %>% 
    dplyr::rename(Ord_Din=probeType) %>%
    clean_tibble() %>%
    dplyr::select(Ses_Cgn,Ord_Des,Ord_Din,Ord_Prb, dplyr::everything())
  
  #
  # TBD:: We should use the manifest from above rather than raw aqp_bsp_tib
  #   results...
  #
  aqp_bsp_tib %>% 
    dplyr::inner_join(ses_pos_tab, by=c("Ord_Prb","Ord_Des","Ord_Din")) %>%
    dplyr::mutate(Pos_Dif=Bsp_Pos-start) %>% 
    dplyr::filter(Pos_Dif!=0) %>% 
    dplyr::filter(Bsp_Tag!="UM") %>%
    dplyr::group_by(Ord_Des,Ord_Din,
                    # Aqp_Idx,
                    Bsp_Nxb_Bsc,
                    Bsp_Din_Bsc,
                    # CG_F1,CG_R1,
                    Bsp_Tag) %>% 
    dplyr::summarise(Count=n(), .groups="drop") %>%
    print(n=1000)
  
  #
  # Matching by manifest::
  #
  ses_man_inn_tib <- add_cgn_imp_bsp_man %>% 
    dplyr::mutate(Ord_Prb=stringr::str_to_upper(Ord_Prb_U),
                  Ord_Des=Ord_Des_U,
                  Bsp_Chr=Chromosome_U,
                  Bsp_Pos=Coordinate_U,
                  Bsp_Tag=Bsp_Tag_U,
                  Aqp_Idx=Aqp_Idx_U,
                  Bsp_Din_Scr=Bsp_Din_Scr_U) %>%
    dplyr::select(Ord_Key,Ord_Des,Ord_Din,Ord_Prb,Bsp_Chr,Bsp_Pos,Bsp_Tag,
                  Aqp_Idx,Bsp_Din_Scr) %>%
    dplyr::inner_join(ses_pos_tab, by=c("Ord_Prb","Ord_Des","Ord_Din"))
  
  # Match Cases::
  #
  ses_man_mat_tib <- ses_man_inn_tib %>% 
    dplyr::mutate(Pos_Dif=Bsp_Pos-start) %>% 
    dplyr::filter(Pos_Dif==0) %>% 
    dplyr::select(Bsp_Din_Scr,Ord_Key,Ses_Cgn,Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx,
                  Ord_Prb,Pos_Dif,Bsp_Chr,seqnames)
  ses_man_mat_sum <- ses_man_mat_tib %>% 
    # dplyr::group_by(Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx) %>% 
    dplyr::group_by(Bsp_Tag,Bsp_Din_Scr,Ord_Des,Ord_Din,Aqp_Idx) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  ses_man_mat_sum %>% print(n=base::nrow(ses_man_mat_sum))
  
  # Miss Cases::
  #
  ses_man_mis_tib <- ses_man_inn_tib %>% 
    dplyr::mutate(Pos_Dif=Bsp_Pos-start) %>% 
    dplyr::filter(Pos_Dif!=0) %>% 
    dplyr::select(Bsp_Din_Scr,Ord_Key,Ses_Cgn,Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx,
                  Ord_Prb,Pos_Dif,Bsp_Chr,seqnames)
  ses_man_mis_sum <- ses_man_mis_tib %>% 
    # dplyr::group_by(Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx) %>% 
    dplyr::group_by(Bsp_Tag,Bsp_Din_Scr,Ord_Des,Ord_Din,Aqp_Idx) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  ses_man_mis_sum %>% print(n=base::nrow(ses_man_mis_sum))
  
  
  #
  # Probe Order Sequence match to check position:: BSP vs. Sesame::
  #
  ses_bsp_inn_tib <- aqp_bsp_tib %>% 
    dplyr::select(Ord_Key,Ord_Des,Ord_Din,Ord_Prb,Bsp_Chr,Bsp_Pos,Bsp_Tag,
                  Aqp_Idx,Bsp_Din_Scr) %>%
    dplyr::inner_join(ses_pos_tab, by=c("Ord_Prb","Ord_Des","Ord_Din"))
  
  # Match Cases::
  #
  ses_bsp_mat_tib <- ses_bsp_inn_tib %>% 
    dplyr::mutate(Pos_Dif=Bsp_Pos-start) %>% 
    dplyr::filter(Pos_Dif==0) %>% 
    dplyr::select(Bsp_Din_Scr,Ord_Key,Ses_Cgn,Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx,
                  Ord_Prb,Pos_Dif,Bsp_Chr,seqnames)
  ses_bsp_mat_sum <- ses_bsp_mat_tib %>% 
    # dplyr::group_by(Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx) %>% 
    dplyr::group_by(Bsp_Tag,Bsp_Din_Scr,Ord_Des,Ord_Din,Aqp_Idx) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  ses_bsp_mat_sum %>% print(n=base::nrow(ses_bsp_mat_sum))
  
  # Miss Cases::
  #
  ses_bsp_mis_tib <- ses_bsp_inn_tib %>% 
    dplyr::mutate(Pos_Dif=Bsp_Pos-start) %>% 
    dplyr::filter(Pos_Dif!=0) %>% 
    dplyr::select(Bsp_Din_Scr,Ord_Key,Ses_Cgn,Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx,
                  Ord_Prb,Pos_Dif,Bsp_Chr,seqnames)
  ses_bsp_mis_sum <- ses_bsp_mis_tib %>% 
    # dplyr::group_by(Ord_Des,Ord_Din,Bsp_Tag,Aqp_Idx) %>% 
    dplyr::group_by(Bsp_Tag,Bsp_Din_Scr,Ord_Des,Ord_Din,Aqp_Idx) %>% 
    dplyr::summarise(Count=n(), .groups="drop")
  ses_bsp_mis_sum %>% print(n=base::nrow(ses_bsp_mis_sum))
  
  #
  # BSP:: Missing Probes???
  #
  aqp_bsp_tib %>% 
    dplyr::select(Ord_Key,Ord_Des,Ord_Din,Ord_Prb,Bsp_Chr,Bsp_Pos,Aqp_Idx) %>%
    dplyr::anti_join(ses_pos_tab, by=c("Ord_Prb","Ord_Des","Ord_Din")) %>% 
    dplyr::add_count(Ord_Prb, name="Aln_Cnt") %>%
    dplyr::group_by(Aqp_Idx,Aln_Cnt,Ord_Des,Ord_Din) %>%
    dplyr::summarise(Count=n(), .groups="drop") %>% print(n=10000)
  
  #
  # Compare aqp_seq_tib vs. ses_man_tib
  #
  aqp_bsp_tib %>% dplyr::filter(Ord_Key %in% aqp_seq_tib$Address) %>%
    dplyr::select(Bsp_Seq,Aln_Prb)
  
  aqp_bsp_tib %>% dplyr::filter(!Bsp_Seq %in% aqp_seq_tib$Aln_Prb)
  aqp_bsp_tib %>% dplyr::filter(!Aln_Prb %in% aqp_seq_tib$Aln_Prb)
  
  aqp_bsp_tib %>% dplyr::select(Ord_Des:Aln_Key) %>% dplyr::arrange(Aln_Prb)
  aqp_seq_tib %>% dplyr::select(Address:Aln_Prb)
}


#
# All srd probe extraction::
#  - Seems to be working, will validate later...
#  - Not really needed at this moment, but useful...
#
if (FALSE) {
  
  all_srd_prb_tib <- aqp_bsp_tib %>% head(n=10) %>%
    bed_to_prbs(fas=run$gen_ref_fas, 
                din="Ord_Din", cgn="Aln_Key",
                chr="Bsp_Chr", pos="Bsp_Pos",
                nrec = 1,
                verbose=opt$verbose+10, tt=pTracker)
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
