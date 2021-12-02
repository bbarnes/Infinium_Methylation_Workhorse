
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
opt$gen_dir  <- NULL
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
opt$genome_build <- NULL
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
  opt$gen_dir  <- file.path(par$topDir, 'data/iGenomes/Homo_sapiens/NCBI')
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
  
  if (FALSE) {
    
  } else if (par$local_runType=='McMaster10Kselection') {
    opt$genome_build <- 'GRCh37'
    opt$platform <- 'MCM'
    opt$version  <- 'v2.5'
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
    opt$genome_build <- 'GRCh37'
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
    opt$genome_build <- 'GRCh37'
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
    opt$genome_build <- 'GRCh37'
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
    opt$genome_build <- 'GRCh37'
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
    opt$genome_build <- 'GRCh38'
    opt$genome_build <- 'GRCh37'
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
    opt$genome_build <- 'GRCh36'
    opt$genome_build <- 'GRCh38'
    opt$genome_build <- 'GRCh37'
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
    opt$genome_build <- 'GRCm38'
    opt$genome_build <- 'GRCm10'
    opt$Species  <- "Mouse"
    opt$version  <- 'M0'
    
    opt$gen_dir <- file.path(par$topDir, 'data/iGenomes/Mus_musculus/NCBI')
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
    opt$genome_build <- 'GRCh36'
    opt$genome_build <- 'GRCh38'
    opt$genome_build <- 'GRCh37'
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
  
  opt$runName <- paste(par$local_runType,opt$platform,opt$version,opt$genome_build, sep='-')
  
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
    make_option(c("--gen_dir"), type="character", default=opt$gen_dir, 
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
    make_option(c("--genome_build"), type="character", default=opt$genome_build, 
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
              'genome_build','platform','version','bsmap_exe',
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

# For accumilating failed probes and why::
error_ledgar <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      0.2 Load any other pre-defined data
#                        dbCGN, improbe, imGenomes data::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   5.0 Probe Design Validation via imGenome:: 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
#
# NOTE:: Should load files from the first half of this program directly from
#   from CSV files!!!
#
#
bsp_dir <- file.path(opt$outDir, "../McMaster10Kselection-MCM-v2-GRCh37")
bsp_csv <- list.files(file.path(bsp_dir, "bsp_mapping_workflow"), pattern="GRCh37.bsp_mapping_workflow.csv.gz$", full.names = TRUE)
bsp_tib <- safe_read(bsp_csv)

imp_tsv <- file.path(
  par$topDir,"scratch/workhorse_main_dev_latest/McMaster10Kselection-MCM-v2.5-GRCh37",
  "c_improbe_workflow/McMaster10Kselection-MCM-v2.5-GRCh37-FCN_dna.c_improbe_workflow.improbe-designOutput.tsv.gz")
imp_tib <- readr::read_tsv(imp_tsv)

large_tsv <- "/Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-input/min/GRCh37.cgn.min.txt.gz"
large_tib <- 

# print(ord_tib)
print(bsp_tib)
# print(seq_tib)
# print(cgn_tib)

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
#         Initialize Probe Data Structures and Design Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: Functionalize this::
fr1_vec <- c("F","R")
fr2_vec <- c("R","F")
tbs_vec <- c("T","B")
cos_vec <- c("C","O")
des_vec <- c("1","2")

probe_offset_tib <- 
  expand.grid(tbs_vec, cos_vec, des_vec) %>% 
  tibble::as_tibble() %>% 
  purrr::set_names("TB", "CO", "Des") %>% 
  dplyr::mutate(Offset=c( 0,-49,-50,1,1,0,-51,2),
                Length=rep(50,8),
                NxbPos=c(-1,  1,-50,0,0,1,  0,1) ) %>%
  dplyr::mutate(FR1 = rep(fr1_vec, 4), 
                FR2 = rep(fr2_vec, 4),
                Srd_Key = paste0(TB,  CO, Des),
                FR1_Key = paste0(FR1, CO, Des),
                FR2_Key = paste0(FR2, CO, Des)) %>%
  dplyr::select(Srd_Key, FR1_Key, FR2_Key, FR1, FR2, 
                dplyr::everything())
print(probe_offset_tib)

template_offset_tib <- 
  tibble::tibble(
    dns_pad =  -2,
    dns_len = -60,
    din_len =   2,
    ups_len =  60,
    ups_pad =   2,
    tmp_len = 
      base::abs(dns_len) +
      din_len + ups_len,
    ext_len = 
      base::abs(dns_pad) +
      base::abs(dns_len) + 
      base::abs(din_len) +
      base::abs(ups_len) +
      base::abs(ups_pad) )

print(template_offset_tib)

template_flank_len <- 60
template_flank_ext <- 2
dinucleotide_len   <- 2

upstream_offset  <- template_flank_len + template_flank_ext
dnstream_offset  <- template_flank_len + template_flank_ext + dinucleotide_len
template_ext_len <- upstream_offset + dnstream_offset

run$ord_ids_key <- "Prb_Key_Unq"
run$ord_des_key <- "Ord_Des"
run$ord_din_key <- "Ord_Din"
run$ord_inf_key <- "Ord_Inf"

run$bsp_frs_key <- "Bsp_FR"
run$bsp_cos_key <- "Bsp_CO"
run$bsp_tbs_key <- "Bsp_TB"
run$bsp_chr_key <- "Bsp_Chr"
run$bsp_pos_key <- "Bsp_Pos"

run$imp_chr_key <- "Chromosome"
run$imp_pos_key <- "Coordinate"
run$imp_level   <- 3
run$call_inf    <- TRUE

run$des_srd_key <- "Srd_Key"
run$tmp_tbs_key <- "Temp_Strand_TB"
run$tmp_fwd_seq <- "Forward_Sequence"
run$tmp_top_seq <- "Top_Sequence"

run$ord_prb_seq <- "Ord_Prb"
run$ord_prb_aln <- "Aln_Prb"

# This is for new probes::
run$prb_seq <- "Probe_Sequence"

ord_ids_sym <- rlang::sym(run$ord_ids_key)
ord_des_sym <- rlang::sym(run$ord_des_key)
ord_din_sym <- rlang::sym(run$ord_din_key)
ord_inf_sym <- rlang::sym(run$ord_inf_key)

bsp_frs_sym <- rlang::sym(run$bsp_frs_key)
bsp_cos_sym <- rlang::sym(run$bsp_cos_key)
bsp_tbs_sym <- rlang::sym(run$bsp_tbs_key)
bsp_chr_sym <- rlang::sym(run$bsp_chr_key)
bsp_pos_sym <- rlang::sym(run$bsp_pos_key)

tmp_tbs_sym <- rlang::sym(run$tmp_tbs_key)
tmp_fwd_sym <- rlang::sym(run$tmp_fwd_seq)
tmp_top_sym <- rlang::sym(run$tmp_top_seq)

probe_cols <- c(run$ord_ids_key, 
                run$ord_des_key, run$ord_din_key,
                run$bsp_frs_key, run$bsp_cos_key,
                run$bsp_chr_key, run$bsp_pos_key,
                run$ord_prb_seq, run$ord_prb_aln )

probes <- bsp_tib %>% 
  dplyr::select( dplyr::all_of( probe_cols) ) %>% 
  dplyr::mutate( !!ord_inf_sym := 
                   dplyr::case_when(
                     !!ord_des_sym == "2" ~ 2,
                     !!ord_des_sym == "U" ~ 1,
                     !!ord_des_sym == "M" ~ 1,
                     TRUE ~ NA_real_
                   ) %>% as.integer(),
                 chr=stringr::str_remove(!!bsp_chr_sym, "^chr"),
                 template_beg = !!bsp_pos_sym + 
                   template_offset_tib$dns_len,
                 template_len = 
                   base::abs(template_offset_tib$dns_len) +
                   base::abs(template_offset_tib$din_len) +
                   base::abs(template_offset_tib$ups_len),
                 template_end = !!bsp_pos_sym + 
                   template_offset_tib$din_len +
                   template_offset_tib$ups_len,
                 ext_temp_beg = !!bsp_pos_sym +
                   template_offset_tib$dns_len +
                   template_offset_tib$dns_pad,
                 ext_temp_len = 
                   base::abs(template_offset_tib$dns_pad) +
                   base::abs(template_offset_tib$dns_len) +
                   base::abs(template_offset_tib$din_len) +
                   base::abs(template_offset_tib$ups_len) +
                   base::abs(template_offset_tib$ups_pad),
                 ext_temp_end = !!bsp_pos_sym + 
                   template_offset_tib$din_len +
                   template_offset_tib$ups_len +
                   template_offset_tib$ups_pad)

stopifnot(probes$template_end - probes$template_beg == probes$template_len)
stopifnot(probes$ext_temp_end - probes$ext_temp_beg == probes$ext_temp_len)

probes_list <- probes %>% split(.$chr)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   5.1 Pass Over all imGenomes Builds:: 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Genomic:: References (Bisulfite and Non-Converted) are loaded into a table now
#   - Ref
#   - SNP (dbSNP-151)
#
imGenome_tib <- load_imGenomes_table(dir = opt$gen_dir,
                                     genome_build = opt$genome_build, 
                                     ret_list = FALSE,
                                     verbose = opt$verbose, tt=pTracker)

imGenome_list <- imGenome_tib %>% split(f=imGenome_tib$Genome_Key)

# Implement once we're on the cluster
# if (opt$parallel) {
# } else {
# }

for (imGenome in imGenome_tib$Genome_Key) {
  
  imGenome_dat <- imGenome_list[[imGenome]][1,]
  imGenome_Fas <- imGenome_dat %>% dplyr::pull("Path")
  imGenome_Src <- imGenome_dat %>% dplyr::pull("Genome_Source")
  imGenome_GnB <- imGenome_dat %>% dplyr::pull("Genome_Version")
  
  imGenome_FR  <- imGenome_dat %>% dplyr::pull("Genome_Strand_FR")
  imGenome_CO  <- imGenome_dat %>% dplyr::pull("Genome_Strand_CO")
  imGenome_BSC <- imGenome_dat %>% dplyr::pull("Genome_Strand_BSC")
  imGenome_Dna <- imGenome_dat %>% dplyr::pull("Genome_Alphabet")
  
  chrom_seqs   <- Biostrings::readDNAStringSet(filepath=imGenome_Fas, 
                                               format="fasta")
  chrom_names  <- intersect( names(probes_list), names(chrom_seqs) )
  
  
  # if (opt$parallel) {
  # } else {
  # }
  
  for (chr in chrom_names) {
    if (opt$verbose>=0) cat(glue::glue("Parsing {chr}...{RET}"))
    
    chr_probes_cnt <- probes_list[[chr]] %>% base::nrow()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Extended Template Sequences::
    #               Meant for Tri-fecta Probe Identification/Discovery::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    opt$check_trifecta_flanks <- FALSE
    if (opt$check_trifecta_flanks) {
      # Note: This can be done on any genome input. 
      
      ext_template_ranges <- 
        IRanges::IRanges( start = probes_list[[chr]]$ext_temp_beg, 
                          width = probes_list[[chr]]$ext_temp_len,
                          names = probes_list[[chr]] %>% dplyr::pull(run$ord_ids_key) %>% 
                            paste("ext_fwd_template", sep="_") )
      
      ext_fwd_template_seqs <- 
        Biostrings::extractAt( chrom_seqs[[chr]], at = ext_template_ranges )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Standard Forward Template Sequences::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    fwd_template_ranges <- 
      IRanges::IRanges( start = probes_list[[chr]]$template_beg, 
                        width = probes_list[[chr]]$template_len,
                        names = probes_list[[chr]] %>% dplyr::pull(run$ord_ids_key) )
    
    fwd_template_seqs <- 
      Biostrings::extractAt( chrom_seqs[[chr]], at = fwd_template_ranges )
    
    fwd_template_tib <- as.character( fwd_template_seqs ) %>% cbind() %>% 
      tibble::as_tibble( rownames = run$ord_ids_key ) %>% 
      purrr::set_names( run$ord_ids_key, run$tmp_fwd_seq )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Only Call Top/Bot on 
    #                         Forward/non-snp/non-bsc 
    #                           Reference Templates
    #
    #       Only c_improbe on Forward/non-snp/non-bsc Reference Templates
    #
    #            Both Top/Bot calling and Thermodynamic Calculations 
    #               assume there there is only four bases. 
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Important Points about Strand_TB at this point::
    #   - Validate c_improbe vs. r_improbe set_topbot_tib() function
    #   - Since this genome always happens first the Strand_TB should be 
    #     passed on in each round.
    #   - It needs to be clearly calculated which strand designation goes
    #     to the template and which goes to the probe that aligns to the 
    #     correct template. 
    
    if (imGenome_Dna=="dna") {
      
      if (imGenome_BSC=="N") {
        
        fwd_top_cols <- c(run$des_srd_key,
                          run$ord_ids_key, run$ord_des_key, run$ord_din_key,
                          run$bsp_frs_key, run$bsp_tbs_key, run$bsp_cos_key,
                          run$bsp_chr_key, run$bsp_pos_key,
                          run$ord_prb_seq, run$ord_prb_aln,
                          run$tmp_tbs_key, run$tmp_fwd_seq, run$tmp_top_seq )
        
        key_cols <- c( run$bsp_tbs_key, run$bsp_cos_key, run$ord_inf_key )
        
        top_template_tib <- fwd_template_tib %>% 
          dplyr::mutate( !!run$tmp_fwd_seq := add_brac( !!tmp_fwd_sym ) ) %>%
          set_topbot_tib( seq_key = run$tmp_fwd_seq, 
                          top_col = run$tmp_tbs_key, 
                          top_key = run$tmp_top_seq ) %>%
          dplyr::right_join( probes_list[[chr]], by=c( run$ord_ids_key ) ) %>%
          dplyr::mutate(
            !!bsp_tbs_sym := dplyr::case_when(
              !!bsp_frs_sym=="F" & !!bsp_cos_sym=="C" ~ !!tmp_tbs_sym,
              !!bsp_frs_sym=="F" & !!bsp_cos_sym=="O" ~ !!tmp_tbs_sym,
              !!bsp_frs_sym=="R" & !!bsp_cos_sym=="C" ~ cmpl_TB(!!tmp_tbs_sym),
              !!bsp_frs_sym=="R" & !!bsp_cos_sym=="O" ~ cmpl_TB(!!tmp_tbs_sym),
              TRUE ~ NA_character_
            ) ) %>% 
          tidyr::unite( !!run$des_srd_key, dplyr::all_of( key_cols ), 
                        sep='', remove=FALSE ) %>%
          dplyr::select( dplyr::all_of( fwd_top_cols ), dplyr::everything() )
        
        # Warning message caused by topBot code::
        #
        # Warning message:
        #   The `x` argument of `as_tibble.matrix()` must have unique 
        #     column names if `.name_repair` is omitted as of tibble 2.0.0.
        # Using compatibility `.name_repair`.
        # This warning is displayed once every 8 hours.
        # Call `lifecycle::last_warnings()` to see where this warning was 
        #   generated. 
        #
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #
        #                Infinium Methylation Probe Design Methods::
        #             c-improbe = c++ improbe (traditional) via docker image
        #                    Includes Thermodynamic Calculations
        #                   Only Designs Infinium I U/M Probes 
        #
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        #
        # TBD:: Probe Strand TB and Template Sequence Strand TB!!!
        #
        
        # 597
        # Add Ord_Prb...
        # top_template_tib %>% dplyr::inner_join(c_imp_tib, by=c("Prb_Key_Unq"="Seq_ID", dplyr::all_of(intersect(names(top_template_tib), names(c_imp_tib) ) ) ) ) %>% dplyr::select()
        
        # dplyr::bind_rows(
        #   dplyr::select(c_imp_tib, Probe_Seq_U,Imp_U49,Scr_U, dplyr::everything()) %>% purrr::set_names("Prb_Seq", "Prb_P49", "Prb_Scr"),
        #   dplyr::select(c_imp_tib, Probe_Seq_U,Imp_U49,Scr_U, dplyr::everything()) %>% purrr::set_names("Prb_Seq", "Prb_P49", "Prb_Scr")
        # )
        # 
        # dplyr::bind_rows(
        #   dplyr::select(c_imp_tib, Probe_Seq_U,Imp_U49,Scr_U) %>% purrr::set_names("Prb_Seq", "Prb_P49", "Prb_Scr"),
        #   dplyr::select(c_imp_tib, Probe_Seq_M,Imp_M49,Scr_M) %>% purrr::set_names("Prb_Seq", "Prb_P49", "Prb_Scr")
        # )
        
        
        # dplyr::inner_join( c_imp_tib,top_template_tib, 
        #                    
        #                    by=c( "Seq_ID"=run$ord_ids_key, 
        #                          "Chromosome"=run$chr_key, 
        #                          "Coordinate"=run$pos_key, 
        #                          "Forward_Sequence", 
        #                          "Top_Sequence", 
        #                          "Ord_Des", 
        #                          "Ord_Din"),
        #                    suffix=c( "_c_imp", "_ord" ) ) %>% 
        #   
        #   dplyr::select(Seq_ID, dplyr::contains("Strand"),
        #                 Probe_Seq_U,Probe_Seq_M,
        #                 Ord_Prb, Aln_Prb, 
        #                 Bsp_FR, Bsp_CO) %>% 
        #   dplyr::select(-dplyr::starts_with("Genome_Strand") ) %>% 
        #   dplyr::filter(Strand_CO==Bsp_CO & Strand_Ref_FR==Bsp_FR)
        
        # This one makes sense!!!!!
        #
            # c_imp_tib %>% 
            #   dplyr::filter(Inf_Type==1) %>% 
            #   dplyr::mutate(Scr_M=as.character(Scr_M), 
            #                 Scr_U=as.character(Scr_U)) %>% 
            #   dplyr::select(Seq_ID, Chromosome, Coordinate,
            #                 Probe_Seq_U,Probe_Seq_M,
            #                 Scr_U,Scr_M,
            #                 Imp_U49, Imp_M49) %>% 
            #   dplyr::rename(Probe_U=Probe_Seq_U, 
            #                 Probe_M=Probe_Seq_M) %>% 
            #   tidyr::pivot_longer(cols = c("Probe_U","Probe_M"), 
            #                       names_to = c("MUD"), 
            #                       values_to="Probes", 
            #                       names_prefix = "Probe_") %>%
            #   tidyr::pivot_longer(cols = c("Scr_U","Scr_M"), 
            #                       names_to = c("MUD1"),
            #                       values_to = "Scores") %>%
            #   tidyr::pivot_longer(cols = c("Imp_U49","Imp_M49"), 
            #                       names_to = c("MUD2"),
            #                       values_to="Probes_Aln49")
            # 
        
        #
        # - Generalized Comparison Method
        # - build all s_improbe & r_improbe
        # - Gneralized Summary Stats (alignment, nxt/ext base, SNPs, color balance)
        # - Implement name prefix optimization cg -> mu, etc.
        # - Support singletons (unpaired reads)
        # - Support swifthoof annotation coverage summary files
        # - Look into Minfi support
        # - VCF Output???
        #
        
        if (run$c_improbe) {
          
          run$re_load <- FALSE
          
          # TBD:: Validate that the Top Sequence Calculations are the same!
          #
          # TBD:: Option to return a table rather than tibble, i.e. one probe
          #   per line!
          #
          
          c_imp_tib <- NULL
          c_imp_tib <- top_template_tib %>% 
            dplyr::select( run$ord_ids_key, run$tmp_fwd_seq, 
                           run$bsp_chr_key, run$bsp_pos_key) %>% 
            dplyr::rename( Chromosome=run$bsp_chr_key, 
                           Coordinate=run$bsp_pos_key ) %>%
            c_improbe_workflow(imGenome = imGenome_dat, # Add genome source row::
                               
                               ids_key = run$ord_ids_key,
                               fwd_seq = run$tmp_fwd_seq,                              
                               pos_key = run$imp_pos_key,  # run$bsp_pos_key,
                               chr_key = run$imp_chr_key,  # run$bsp_chr_key,
                               
                               doc_image = run$doc_image,
                               doc_shell = run$doc_shell,
                               
                                                         # Flag to tell the prg
                                                         #  to make an Infinium 
                                                         #  I/II call based on 
                               call_inf = run$call_inf,  #  scores/cpg-counts
                               
                               outlevel = run$imp_level, # Output verbosity level
                               
                               re_join  = run$rejoin,
                               new_join = run$join_new_vec,
                               old_join = run$join_old_vec,
                               
                               reload = run$reload,
                               
                               out_dir = opt$outDir,
                               run_tag = paste(opt$runName,imGenome,sep='-'),
                               re_load = run$re_load,
                               pre_tag = pTracker$file_vec,
                               end_str = 'tsv.gz',
                               
                               verbose=opt$verbose, tt=pTracker)
          # verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
        }
      }
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #
      #                Infinium Methylation Probe Design Methods::
      #                        r-improbe re-implemented
      #                    Allows All Probe Designs:: cg,ch,rs 
      #
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      if (r_improbe) {
        
        # If dna calculation Top Sequence
        # If snp do NOT calculation Top Sequence
        
        top_col <- "Strand_TB"
        top_key <- "Top_Sequence"
        add_topseq <- FALSE
        if (cur_gen_tib$Gen_Alphabet=="dna")
          add_topseq <- TRUE
        
        r_imp_tib <- NULL
        r_imp_tib <- 
          r_improbe_workflow(tib = s_imp_tib,
                             
                             # Add genome source row::
                             gen_tib = cur_gen_tib,
                             
                             ids_key = ids_key,
                             seq_key = imp_seq,
                             din_key = din_key,
                             
                             top_col = top_col,
                             top_key = top_key,
                             
                             srsplit = srsplit,
                             srd_key = "Strand_FR",
                             srd_str = srd_str,
                             
                             cosplit = cosplit,
                             cos_key = "Strand_CO",
                             cos_str = cos_str,
                             
                             ups_len = ups_len,
                             seq_len = seq_len,
                             
                             # subset   = subset,
                             # sub_cols = sub_cols,
                             
                             reload     = reload,
                             parallel   = parallel,
                             add_matseq = add_matseq,
                             add_topseq = add_topseq,
                             
                             out_dir = out_dir,
                             run_tag = paste(run_tag,cur_gen_key,sep='-'),
                             re_load = TRUE,
                             pre_tag = tt$file_vec,
                             
                             verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
        
      }
      
      
      
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                          Extract Probes Sequences:: 
    #                           Targeted or Design All
    #                               s_improbe()
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    cur_probe_tib <- NULL
    
    cur_offset_tib <- 
      dplyr::right_join( probe_offset_tib, fwd_top_tib, by=run$des_srd_key ) %>% 
      split(f=dplyr::pull(., run$des_srd_key) )
    
    for ( key in names(cur_offset_tib) ) {
      
      # - Can't pull probe sequences from dna, must be bsc!!!
      #
      # - We can run r_improbe & c_improbe
      #
      
      probe_ranges <-
        IRanges::IRanges( 
          start = 
            cur_offset_tib[[key]] %>% dplyr::pull( run$bsp_pos_key ) +
            cur_offset_tib[[key]] %>% dplyr::pull( Offset ),
          width = cur_offset_tib[[key]] %>% dplyr::pull( Length ),
          names = cur_offset_tib[[key]] %>% dplyr::pull( run$ord_ids_key ) )
      
      probe_seqs <- 
        Biostrings::extractAt( chrom_seqs[[chr]], at = probe_ranges ) %>%
        as.character() %>% cbind() %>% 
        tibble::as_tibble( rownames = run$ord_ids_key ) %>% 
        purrr::set_names( run$ord_ids_key, run$prb_seq ) %>% 
        dplyr::right_join(cur_offset_tib[[key]], by=c(run$ord_ids_key) )
      
      key_cols <- c( run$ord_ids_key, run$bsp_chr_key, run$bsp_pos_key,
                     run$bsp_frs_key, run$bsp_tbs_key, run$bsp_cos_key, 
                     run$ord_inf_key )
      
      c_dna_tib %>% dplyr::select(Seq_ID, Chromosome, Coordinate, 
                                  Forward_Sequence, Top_Sequence, 
                                  Probe_Seq_U, Probe_Seq_M, 
                                  Strand_FR, Strand_TB, Strand_CO, 
                                  Ord_Din, Ord_Des ) %>% 
        dplyr::inner_join(probe_seqs , by=c("Seq_ID" = "Prb_Key_Unq", 
                                            "Chromosome" = "Bsp_Chr",
                                            "Coordinate" = "Bsp_Pos",
                                            "Strand_FR"  = "Bsp_FR",
                                            "Strand_CO"  = "Bsp_CO") ) %>%
        dplyr::mutate(
          Match=dplyr::case_when(
            Ord_Des=="U" & Probe_Seq_U==Probe_Sequence ~ 0,
            Ord_Des=="U" & Probe_Seq_U==Probe_Sequence ~ 1,
            Ord_Des=="M" & Probe_Seq_M==Probe_Sequence ~ 0,
            Ord_Des=="M" & Probe_Seq_M==Probe_Sequence ~ 1,
          )
        )
      
      break
    }
    
    break
  }
  
  break
}

template_ranges <- 
  IRanges::IRanges( start = template_beg_vec, 
                    width = template_ext_len,
                    names = probes[[chr]]$Prb_Key )

template_beg_vec <- dplyr::pull( probes[[chr]], pos_key ) - upstream_offset
template_ranges <- 
  IRanges::IRanges( start = template_beg_vec, 
                    width = template_ext_len,
                    names = probes[[chr]]$Prb_Key )

begs <- NULL
ends <- NULL
begs <- tib %>% dplyr::pull(!!pos_sym) - ups_len - pad_len
ends <- begs + seq_len - 1 + (2 * pad_len)


seq_vec <- stringr::str_sub( as.character(seq), begs, ends ) %>%
  stringr::str_to_upper()
ret_cnt <- seq_vec %>% length()





# - DMAP validation???
# - Add idat validation
#
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


gen_cnt <- prb_designs_workflow(
  tib = bsp_tib %>% 
    dplyr::select(run$unq_key, run$ids_key, run$add_key, 
                  run$des_key, run$din_key, 
                  run$srd_key, run$cos_key,
                  run$pos_key, run$chr_key, 
                  run$prb_key, "Bsp_Str", "Bsp_Tag", "Bsp_Prb_Dir"), 
  max = 0,
  
  out_dir = opt$outDir,
  run_tag = opt$runName, 
  pre_tag = pTracker$file_vec, 
  re_load = TRUE,
  
  imp_level  = 3,
  
  # Genomes Parameters::
  gen_bld  = opt$genome_build,
  gen_nrec = run$gen_nrec,
  gen_key  = run$gen_key,
  gen_tib  = imGenome_tib,
  
  # Field Parameters:: general
  ids_key = run$unq_key,
  des_key = run$des_key, 
  din_key = run$din_key,
  
  # bsp_srd = run$bsp_srd,
  # bsp_cos = run$bsp_cos,
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
  
  # subset   = FALSE,
  # sub_cols = NULL,
  
  # reload=opt$reload,
  reload  = TRUE,
  retData = TRUE,
  # retData=FALSE,
  
  parallel=opt$parallel,
  # parallel=FALSE,
  
  r_improbe = TRUE,
  s_improbe = TRUE,
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

tmp_csv <- "/Users/bretbarnes/Documents/scratch/workhorse_main_dev_latest/McMaster10Kselection-MCM-v1-GRCh37/pre/s_improbe_workflow/McMaster10Kselection-MCM-v1-GRCh37-FCM_dna.s_improbe_workflow.csv.gz"
s_tib <- safe_read(tmp_csv, clean = FALSE)

# Probe Design Scratch::
#
if (FALSE) {
  
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
    files_to_tib(s_improbe_dir,s_improbe_pat)
  ) %>% dplyr::distinct()
  
  file_cnt <- file_tab %>% base::nrow()
  file_vec <- c(1:file_cnt)
  
  improbe_file_list <- file_tab %>% split(.$Key)
  
  # Check out c_improbe/r_improbe::
  c_dna_tib <- safe_read(improbe_file_list[["c_dna_FCN"]]$Path, 
                         verbose=opt$verbose)
  
  r_dna_tib <- safe_read(improbe_file_list[["r_dna_FCN"]]$Path, 
                         verbose=opt$verbose)
  
  r_snp_tib <- safe_read(improbe_file_list[["r_snp_FCN"]]$Path, 
                         verbose=opt$verbose)
  
  top_tib <- r_dna_tib %>% head() %>% set_topbot_tib(seq_key="Forward_Sequence", top_col="Strand_TB", top_key="Top_Sequence", verbose = 1000)
  
  top_tib %>% dplyr::select(Prb_Key_Unq_srd, Strand_TB, Top_Sequence)
  
  c_dna_tib %>% dplyr::inner_join(top_tib, by=c("Seq_ID"="Prb_Key_Unq", "Forward_Sequence", "Top_Sequence", "Gen_Genome_Build", "Gen_Source","Gen_Alphabet","Genome_Key", "Gen_Strand_FR", "Gen_Strand_CO", "Gen_Strand_BSC", "Strand_TB", "Ord_Des", "Ord_Din", "Probe_Seq_U"="Prb_1U", "Probe_Seq_M"="Prb_1M")) 
  
  # Load files::
  prb_des_tib <- NULL
  for (i in file_vec) {
    file_key <- file_tab[i,]$Key
    file_imp <- file_tab[i,]$improbe_type
    file_fns <- file_tab[i,]$Path
    
    if ((file_tab[i,]$improbe_type=="s" & 
         file_tab[i,]$Alphabet=="dna" & 
         file_tab[i,]$Strand_BSC != "N")) {
      
      cur_gen_tib <- file_tab[i,] %>% 
        dplyr::select(-Path) %>%
        dplyr::rename(Gen_Strand_FR=Strand_FR,
                      Gen_Strand_CO=Strand_CO,
                      Gen_improbe_type=improbe_type)
      
      cat(glue::glue("[{par$prgmTag}]: Loading {file_key}={file_fns}...{RET}"))
      
      cur_tib <- safe_read(file_fns, verbose=opt$verbose) %>% 
        dplyr::mutate(Tmp_FR = as.character(Tmp_FR))
      
      prb_des_tib <- prb_des_tib %>% dplyr::bind_rows(
        bind_cols(
          cur_gen_tib,
          dplyr::bind_rows(
            cur_tib %>% dplyr::select(-Prb_1O, -Prb_2O) %>% dplyr::mutate(Des_Srd="C"),
            cur_tib %>% dplyr::select(-Prb_1C, -Prb_2C) %>% dplyr::mutate(Des_Srd="O")
          )
        )
      )
      
      cat(glue::glue("[{par$prgmTag}]: Done.{RET2}"))
      
    }
  }
  
  prb_des_tib %>% dplyr::select(Prb_Key_Unq, Key, Gen_Strand_FR,Gen_Strand_CO,Strand_BSC, Strand_FR,Strand_CO, Tmp_FR,Tmp_CO,Tmp_Bsc) %>% dplyr::filter(Gen_Strand_CO == Strand_CO)
  
  # - Select s + dna
  # - bind_rows
  # - pivot longer [Prb_1C,Prb_2C, Prb_1O,Prb_2O]
  # - Select Rows where Tmp_FR==BSP_FR & Tmp_CO==Bsp_CO
  #
  # - Compare U/M against c-improbe
  tar_key_vec <- file_tab %>% 
    dplyr::filter(improbe_type=="s" & Alphabet=="dna" & Strand_BSC != "N") %>% 
    dplyr::pull(Key)
  
  s_dna_tib <- NULL
  for (key in tar_key_vec) {
    s_dna_tib <- s_dna_tib %>% 
      dplyr::bind_rows(
        prb_des_list[[4]] %>% dplyr::select(-Prb_1O, -Prb_2O) %>% dplyr::mutate(Des_Srd="C"),
        prb_des_list[[4]] %>% dplyr::select(-Prb_1C, -Prb_2C) %>% dplyr::mutate(Des_Srd="O"),
      )
  }
  
  
  
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
  
  # Add improbe_type = [c,r,s]
  # Add Tmp_FR for everything!!!
  # s:: rename Strand_FR=BSP_FR, Strand_CO=CO
  # r:: add Ord_Des
  
  # c_improbe:: NOTE: we can check topBot function here with Forward_Sequence
  #
  c_imp_tib <- data_des[["c_dna_FCN"]] %>% 
    dplyr::select(Seq_ID,Ord_Des,Ord_Din, 
                  Probe_Seq_U,Probe_Seq_M, 
                  Strand_TB,Strand_CO,Strand_FR,Strand_Ref_FR,
                  Tmp_FR,Tmp_CO,Tmp_BC,Tmp_AB, improbe_type,
                  Forward_Sequence, Top_Sequence) %>% 
    dplyr::mutate(Tmp_FR="F")  %>% 
    setTopBot_tib(seqKey = "Forward_Sequence",  srdKey = "Auto_Strand_TB", 
                  topKey = "Auto_Top_Sequence", verbose = opt$verbose)
  
  
  r_imp_tib <- data_des[["r_dna_FCN"]] %>% 
    dplyr::select(Prb_Key_Unq,Ord_Des,Ord_Din, 
                  Prb_1U,Prb_1M,Prb_2D, 
                  Strand_FR,Strand_CO, 
                  Tmp_FR,Tmp_CO,Tmp_BC,Tmp_AB, improbe_type,
                  Forward_Sequence) %>% 
    dplyr::mutate(Tmp_FR="F") %>% 
    setTopBot_tib(seqKey = "Forward_Sequence",  srdKey = "Auto_Strand_TB", 
                  topKey = "Auto_Top_Sequence", verbose = opt$verbose)
  
  # 30 Probe Sequence Match Failures::
  dplyr::inner_join(c_imp_tib,r_imp_tib, 
                    by=c("Seq_ID"="Prb_Key_Unq", 
                         "Ord_Des","Ord_Din",
                         "Strand_TB"="Auto_Strand_TB","Strand_CO", 
                         "Tmp_FR","Tmp_CO","Tmp_BC","Tmp_AB", 
                         "Probe_Seq_U"="Prb_1U", 
                         "Probe_Seq_M"="Prb_1M" ) )
  
  # - Compare Top Calculations::
  #
  dplyr::inner_join(c_imp_tib,r_imp_tib, 
                    by=c("Seq_ID"="Prb_Key_Unq", 
                         "Ord_Des","Ord_Din","Strand_FR",
                         "Auto_Strand_TB","Strand_CO", 
                         "Tmp_FR","Tmp_CO","Tmp_BC","Tmp_AB", 
                         "Auto_Top_Sequence", "Forward_Sequence")) %>%
    dplyr::filter(Top_Sequence==Auto_Top_Sequence)
  
  #
  # Recombine s_improbe:: two-steps
  #  NOTE:: This should actually be done in the s_improbe function...
  #
  file_dat <- file_tab %>% 
    dplyr::filter(improbe_type=='s') %>% 
    select(Key:Strand_BSC) %>% 
    dplyr::arrange(Alphabet) %>%
    split(.$Key)
  
  data_des[[4]] %>% dplyr::mutate(
    Tmp_FR=dplyr::case_when(
      purrr::is_logical(Tmp_FR) ~ "F",
      TRUE ~ Tmp_FR
    )
  )
  
  for (i in c(1:17)) {
    data_des[[i]] %>% names() %>% length() %>% print()
  }
  
  data_des[["s_dna_FCM"]] %>% 
    dplyr::select(Prb_Key_Unq, Ord_Des,Ord_Din, 
                  Strand_FR,Strand_CO, 
                  Prb_1C,Prb_2C,Prb_1O,Prb_2O, 
                  Tmp_FR,Tmp_CO,Tmp_BC,Tmp_AB, improbe_type) %>% 
    dplyr::mutate(Tmp_FR="F")
  
  data_des[["s_dna_FCM"]] %>% 
    dplyr::select(Prb_Key_Unq, Ord_Des,Ord_Din, 
                  Strand_FR,Strand_CO, 
                  Prb1C,Prb2C,Prb1O,Prb2O, 
                  Tmp_FR,Tmp_CO,Tmp_BC,Tmp_AB, improbe_type) %>% 
    dplyr::mutate(Tmp_FR="F") %>% 
    dplyr::rename(Prb_1MC=Prb1C, Prb_2MC=Prb2C, Prb_1MO=Prb1O, Prb_2MO=Prb2O)
  
  
  
  file_dat[[4]] %>% 
    dplyr::select(Prb_Key_Unq, Ord_Des,Ord_Din, 
                  Strand_FR,Strand_CO, 
                  Prb1C,Prb2C,Prb1O,Prb2O, 
                  Tmp_FR,Tmp_CO,Tmp_BC,Tmp_AB, improbe_type) %>% 
    dplyr::mutate(Tmp_FR="F")
  
  
  for (i in file_vec) {
    for (j in file_vec) {
      if (i<=j) {
        
        
        frA <- file_tab[i, ]$Strand_FR
        frB <- file_tab[j, ]$Strand_FR
        
        coA <- file_tab[i, ]$Strand_CO
        coB <- file_tab[j, ]$Strand_CO
        
        if (frA==frB && coA==coB)
          compare_seqs(i,j, file_tab, file_dat)
        
        # break
      }
    }
    break
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

if (opt$genome_build=="GRCh37" || opt$genome_build=="GRCh38") {
  
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
# if (opt$genome_build=="GRCh37" || opt$genome_build=="GRCh38") {
if (FALSE) {
  
  core_anno_dir <- file.path(par$topDir, "data/annotation", opt$genome_build)
  
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
  
  ncib_gene_tsv <- file.path(core_anno_dir, "NCBI", paste(opt$genome_build,"ncbi.RefSeqGenes.tsv.gz", sep='.'))
  ucsc_gene_tsv <- file.path(core_anno_dir, "UCSC", paste(opt$genome_build,"ucsc.knownGene.tsv.gz", sep='.'))
  ucsc_cpgs_tsv <- file.path(core_anno_dir, "UCSC", paste(opt$genome_build,"ucsc.CpG-Islands.tsv.gz", sep='.'))
  
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
  hmm_ann_drn <- paste("GRCh36",opt$genome_build, sep='-')
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
  if (opt$genome_build=="GRCh37") par$sesBuild <- "hg19"
  if (opt$genome_build=="GRCh38") par$sesBuild <- "hg18"
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
