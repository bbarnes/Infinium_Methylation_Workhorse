
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
run <- NULL

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

#
# TBD:: This is copied in the code, so it should be removed from here
#   or moved to a config file... Maybe a single r-data structure (RDS)
#   for manifest parameter defaults is the way to go...
#
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
  par$local_runType <- 'NZT'
  par$local_runType <- 'Chicago'
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
    opt$version  <- 'v1'
    opt$Species  <- "Human"
    
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
    opt$Species  <- "Human"
    opt$version  <- 'B1'
    
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
    opt$platform <- 'COVIC'
    opt$Species  <- "Human"
    opt$version  <- 'C0'
    
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
    opt$genBuild <- 'GRCm38'
    opt$genBuild <- 'GRCm10'
    opt$Species  <- "Mouse"
    opt$version  <- 'M0'
    
    opt$genDir <- file.path(par$topDir, 'data/iGenomes/Mus_musculus/NCBI')
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

# Genomic:: References (Bisulfite and Non-Converted) are loaded into a table now
#   - Ref
#   - SNP (dbSNP-151)
# TBD:: Functionalize this...
run$gen_dir <- file.path(opt$genDir, opt$genBuild,"Sequence/WholeGenomeFasta")
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

# Define Pre-built improbe directories and files::
#   - Using split files now instead of two single large files...
run$cgn_seq_dir <- 
  file.path(opt$impDir, "scratch/cgnDB/dbSNP_Core4/design-output/prbs-p49-split")

stopifnot(dir.exists(run$cgn_seq_dir))

run$cgn_bed_dir <- file.path(opt$impDir, "scratch/cgnDB/dbSNP_Core4/design-input/min")
run$cgn_bed_tsv <- file.path(run$cgn_bed_dir, paste(opt$genBuild,"cgn.min.txt.gz", sep="."))
run$canonical_csv <- file.path(par$datDir, "manifest/cgnDB/canonical-assignment.cgn-top-grp.csv.gz")

stopifnot(dir.exists(run$cgn_bed_dir))
stopifnot(file.exists(run$cgn_bed_tsv))
stopifnot(file.exists(run$canonical_csv))

# Change these to options::
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

#
# MOST declaration of intermediate directories is no longer needed
#   They build themselves....
#

# Manifest Directory::
run$man_dir <- file.path(opt$outDir, 'man')
run$aqp_man_csv  <- file.path(run$man_dir, paste(opt$runName,"aqp-pass.manifest-sesame.csv.gz", sep="."))

# Annotation Directory::
#
run$ann_dir <- file.path(opt$outDir, 'ann')
run$ann_int_csv  <- file.path(run$ann_dir,paste(opt$runName,'cpg-pass.annotation.csv.gz', sep='.'))

if (opt$verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: Done. Defining Intermediate Run Time Files.{RET}{RET}"))

#
#
# Genomes Manifest Generation Parameters::
#  - TBD:: These should be moved to some config file...
#
#
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

run$int_suffix <- "probe-subseq"

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

opt$build_manifest <- FALSE
opt$run_improbe    <- FALSE
par$load_ann       <- FALSE

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#            0.0 Load any pre-defined Standard Manifest to be added::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (!is.null(opt$mans)) {
  #
  # TBD:: User should be able to rebuild an existing or old manifest,
  #  or add manifests together...
  #
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


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   1.0 AQP Address Manifest Workflow: 
#                           Order/Match/AQP/PQC
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker$addFile(opt$time_org_txt)

ord_tib <- NULL
ord_tib <- 
  aqp_address_workflow(ord = opt$ords,
                       mat = opt$mats,
                       aqp = opt$aqps,
                       
                       prb_key = run$prb_key,
                       add_key = run$add_key,
                       des_key = run$des_key,
                       din_key = run$din_key,
                       ids_key = run$ids_key,
                       
                       out_dir = opt$outDir,
                       run_tag = opt$runName,
                       re_load = TRUE,
                       pre_tag = pTracker$file_vec,
                       
                       verbose=opt$verbose+100, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    2.0 Align All Probe Sequence:: BSMAP
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsp_tib <- bsp_mapping_workflow(ref_fas = NULL,
                                ref_tib = all_gen_tib,
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
                                
                                sort    = TRUE,
                                full    = FALSE,
                                merge   = FALSE,
                                light   = TRUE,
                                reload  = opt$reload,
                                retData = FALSE,
                                
                                bsp_exe = opt$bsmap_exe,
                                bsp_opt = opt$bsmap_opt,
                                
                                out_dir = opt$outDir,
                                run_tag = opt$runName,
                                re_load = TRUE,
                                pre_tag = pTracker$file_vec,
                                
                                verbose=opt$verbose, tt=pTracker)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#         3.0 Intersect Sequences Address and improbe:: U49/M49
#                         CGN Mapping Workflow()
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# NOTE: McMaster10Kselection = 822s
seq_tib <- 
  seq_mapping_workflow(ord_tib = ord_tib,
                       
                       seq_dir   = run$cgn_seq_dir,
                       pattern_u = "-probe_U49_cgn-table.tsv.gz", 
                       pattern_m = "-probe_M49_cgn-table.tsv.gz", 
                       
                       prb_key = run$prb_key,
                       add_key = run$add_key,
                       des_key = run$des_key,
                       din_key = run$din_key,
                       ids_key = run$ids_key,
                       
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       4.0 Analyze and Assign Cgn:: 
#                      CGN-Map/BSMAP/dbGCGN look-up
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cgn_tib <- 
  cgn_mapping_workflow(ord_tib = ord_tib,
                       bsp_tib = bsp_tib,
                       seq_tib = seq_tib,

                       ids_key = ids_key,
                       bsp_csv = bsp_key,
                       can_csv = run$canonical_csv,
                       merge   = run$merge,
                       
                       out_dir = opt$outDir,
                       run_tag = opt$runName,
                       re_load = TRUE,
                       pre_tag = pTracker$file_vec,
                       
                       verbose=opt$verbose, tt=pTracker)

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

run$re_load <- TRUE

aqp_prb_tib <- prb_designs_workflow(
  tib = aqp_bsp_tib$bsp_join %>% 
    dplyr::select(run$unq_key, run$ids_key, run$add_key, 
                  run$des_key, run$din_key, 
                  run$srd_key, run$cos_key,
                  run$pos_key, run$chr_key, 
                  run$prb_key, "Bsp_Str", "Bsp_Tag", "Bsp_Prb_Dir"), 
  max = 0,
  
  # out_dir = run$seq_dir,
  out_dir = opt$outDir,
  run_tag = opt$runName, 
  pre_tag = pTracker$file_vec, 
  re_load = TRUE,
  
  imp_level  = 3,
  
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
  
  subset   = FALSE,
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
  add_probes = TRUE,
  add_matseq = TRUE,
  
  verbose=opt$verbose, tt=pTracker
)

#
#
# MAJOR TBD:: ADD FDC/dna/snp/rcp during sequence extraction/probe-design...
#
#


# Probe Design Scratch::
#
if (FALSE) {
  
  files_to_tib = function(dir, pat) {
    
    list.files(dir, pattern = pat, full.names = TRUE) %>% 
      tibble::as_tibble() %>% purrr::set_names("Path") %>%  
      dplyr::mutate(Base_Name=base::basename(Path) %>% 
                      stringr::str_remove(".[ct]sv.gz$")) %>%
      tidyr::separate(Base_Name, into=c("Run_Tag","improbe_type"), sep="\\.",
                      remove = FALSE) %>% 
      dplyr::mutate(Alphabet=stringr::str_remove(Run_Tag, "^.*_"),
                    Bsc_Tag=stringr::str_remove(Run_Tag, "^.*-") %>%
                      stringr::str_remove("_.*$")) %>% 
      tidyr::separate(Bsc_Tag, into=c("Strand_FR","Strand_CO", "Strand_BSC"),
                      sep=c(1,2), remove=FALSE) %>% 
      dplyr::mutate(improbe_type=improbe_type %>% 
                      stringr::str_remove("_improbe_workflow$"),
                    Key=paste(improbe_type, Alphabet, Bsc_Tag, sep='_')) %>%
      dplyr::select(Key, improbe_type, Alphabet, Bsc_Tag, 
                    Strand_FR, Strand_CO, Strand_BSC, Path)
  }
  
  c_improbe_dir <- file.path(opt$outDir, "c_improbe_workflow")
  c_improbe_pat <- paste0("c_improbe_workflow.tsv.gz$")
  # c_improbe_pat <- paste0("c_improbe_workflow.csv.gz$")
  
  r_improbe_dir <- file.path(opt$outDir, "r_improbe_workflow")
  r_improbe_pat <- paste0(".r_improbe_workflow.csv.gz$")
  # r_improbe_pat <- paste0("._dna.csv.gz$")
  
  s_improbe_dir <- file.path(opt$outDir, "s_improbe_workflow")
  s_improbe_pat <- paste0(".s_improbe_workflow.csv.gz$")
  
  file_tab <- dplyr::bind_rows(
    files_to_tib(c_improbe_dir,c_improbe_pat),
    files_to_tib(r_improbe_dir,r_improbe_pat),
    files_to_tib(r_improbe_dir,r_improbe_pat),
    files_to_tib(s_improbe_dir,s_improbe_pat)
  )
  
  #
  # c_U vs. r_U  (dna: Bsc/Ref)
  # c_M vs. r_M  (dna: Bsc/Ref)
  # c_U vs. s_U  (Bsc: Bsc/Ref)
  # c_M vs. s_M  (Bsc: Bsc/Ref)
  #
  # r_U vs. s_U
  #
  # dna vs. snp   (r)
  #
  #  s vs.   c   (Ref-dna) - Sanity Check: Are you getting the same coordinates???
  #  s vs.   r   (Ref-dna) - Sanity Check: Are you getting the same coordinates???
  #
  file_cnt <- file_tab %>% base::nrow()
  file_vec <- c(1:file_cnt)
  
  file_list <- file_tab %>% split(.$Key)
  
  # lapply(list, function)
  
  # Load files::
  for (i in file_vec) {
    file_dat[i] <- safe_read
  }
  
  for (i in file_vec) {
    for (j in file_vec) {
      compare_seqs(i,j, file_tab, file_dat)
    }
  }
  
  
  ord_des_tib <- aqp_bsp_tib$bsp_join %>% 
    dplyr::select(run$unq_key, run$ids_key, run$add_key, 
                  run$des_key, run$din_key, 
                  run$srd_key, run$cos_key,
                  run$pos_key, run$chr_key, 
                  run$prb_key, "Bsp_Str", "Bsp_Tag", "Bsp_Prb_Dir")
  
  r_improbe_fwd_tib <- dplyr::select(ord_des_tib, Prb_Key_Unq, Ord_Din, Bsp_FR) %>% 
    dplyr::distinct() %>% dplyr::right_join(c_imp_tib, by=c("Prb_Key_Unq"="Seq_ID")) %>%
    dplyr::select(Prb_Key_Unq, Forward_Sequence, Ord_Din, Bsp_FR) %>%
    r_improbe(ids_key = "Prb_Key_Unq", seq_key = "Forward_Sequence", din_key = "Ord_Din", srd_str = "FR", verbose = 40)
  
  
  #
  #
  # Head-to-head test:: old implementation vs. new::
  #
  #
  r_improbe_testset <- dplyr::select(ord_des_tib, Prb_Key_Unq, Ord_Din) %>% dplyr::distinct() %>% 
    dplyr::right_join(safe_read(file_tab %>% dplyr::filter(Alphabet=="dna" & improbe_type=="c") %>% dplyr::pull(Path)), 
                      by=c("Prb_Key_Unq"="Seq_ID")) %>% 
    dplyr::filter(Strand_FR=="F", Strand_CO=="C") %>% 
    dplyr::select(Prb_Key_Unq,Forward_Sequence,Ord_Din,Strand_FR,Strand_CO)
  
  r_improbe_old_tib <- r_improbe_testset %>%
    desSeq_to_prbs(idsKey = "Prb_Key_Unq", 
                   seqKey = "Forward_Sequence", 
                   prbKey = "Ord_Din", 
                   strsSR = "FR", 
                   strsCO = "CO", 
                   parallel = FALSE, 
                   verbose  = 3, tt=pTracker)
  
  r_improbe_new_tib <- r_improbe_testset %>%
    r_improbe(ids_key = "Prb_Key_Unq",
              seq_key = "Forward_Sequence",
              din_key = "Ord_Din",
              srd_str = "FR",
              cos_str = "CO",
              parallel = FALSE, 
              verbose  = 3, 
              tt=pTracker)
  
  r_improbe_new_tib <- r_improbe_testset %>%
    r_improbe(ids_key = "Prb_Key_Unq",
              seq_key = "Forward_Sequence",
              din_key = "Ord_Din",
              srsplit = TRUE,
              srd_key = "Strand_FR",
              srd_str = "FR",
              cosplit = TRUE,
              cos_key = "Strand_CO",
              cos_str = "CO",
              parallel = TRUE,
              verbose = 30, 
              tt=pTracker)
  
  
  c_imp_tib <- safe_read(file_tab %>% dplyr::filter(Alphabet=="dna" & improbe_type=="c") %>% dplyr::pull(Path)) %>%
    dplyr::select(Seq_ID,Strand_Ref_FR,Strand_TB,Strand_CO, Probe_Seq_U,Probe_Seq_M)
  
  s_imp_tib_FCM <- readr::read_csv(file_tab %>% dplyr::filter(Alphabet=="dna" & improbe_type=="s" & Strand_CO=="C", Strand_FR=="F" & Strand_BSC=="M") %>% dplyr::pull(Path)) %>% 
    setTopBot_tib(seqKey="Forward_Sequence", srdKey="Raw_TB", topKey = "Top_Seq") %>%
    dplyr::select(Raw_TB,Bsp_FR, Top_Seq, starts_with("Prb"))
  
  r_imp_tib <- readr::read_csv(file_tab %>% dplyr::filter(Alphabet=="dna" & improbe_type=="r") %>% dplyr::pull(Path)) %>%
    dplyr::select(r_imp_tib, Seq_ID,Forward_Sequence,Strand_SR,Strand_CO, PRB1_U_MAT,PRB1_M_MAT,PRB2_D_MAT) %>%
    setTopBot_tib(seqKey="Forward_Sequence", srdKey="Raw_TB", topKey = "Top_Seq") %>%
    dplyr::distinct()
  
  
  cs_FCM_inn <- dplyr::inner_join(c_imp_tib,s_imp_tib_FCM, by=c("Seq_ID"="Prb_Key_Unq") ) %>% 
    dplyr::filter(Probe_Seq_M == Prb1C)
  cs_FCM_inn %>%
    dplyr::group_by(Strand_Ref_FR,Strand_TB,Strand_CO, Bsp_FR,Raw_TB) %>%
    dplyr::summarise(Count=n(), .groups = "drop")
  
  rs_FCM_inn <- dplyr::inner_join(r_imp_tib,s_imp_tib_FCM, by=c("Seq_ID"="Prb_Key_Unq") ) %>%
    dplyr::filter(PRB1_M_MAT == Prb1C)
  rs_FCM_inn %>%
    dplyr::group_by(Strand_Ref_FR,Strand_TB,Strand_CO, Bsp_FR,Raw_TB) %>%
    dplyr::summarise(Count=n(), .groups = "drop")
  
  
}






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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#                              END OF ROUND 2::
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


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
