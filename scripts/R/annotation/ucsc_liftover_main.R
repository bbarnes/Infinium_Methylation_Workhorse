
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Genomic Ranges::
suppressWarnings(suppressPackageStartupMessages( base::require("GenomicRanges",quietly=TRUE) ))

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
par$prgmTag <- 'ucsc_liftover_main'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))

# Executables::
opt$Rscript  <- NULL
opt$lift_exe <- "/Users/bretbarnes/Documents/programs/ucsc/liftOver"

# Run Parameters:: Required
opt$runName  <- NULL
opt$oldBuild <- NULL
opt$newBuild <- NULL

# Directories:: Required
opt$outDir  <- NULL
opt$annDir  <- NULL

# Input Files:: Required
opt$chain    <- NULL
opt$chainKey <- NULL

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
opt$fresh   <- FALSE

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

  # Pre-defined local options runTypes::
  #
  par$local_runType <- NULL
  par$local_runType <- "EWAS"
  
  if (par$local_runType=='EWAS') {

    # Human::
    opt$oldBuild <- 'GRCh36'
    
    opt$newBuild <- 'GRCh37'
    opt$chainKey <- 'hg18ToHg19'

    opt$newBuild <- 'GRCh38'
    opt$chainKey <- 'hg18ToHg38'
    
    opt$liftSub  <- 'wgEncodeBroadHmm'

    # Mouse::
    opt$oldBuild <- 'GRCm9'
    opt$newBuild <- 'GRCm10'
    opt$chainKey <- 'mm9ToMm10'
    opt$liftSub  <- 'bed'

    opt$chain    <- file.path(par$topDir, "programs/ucsc/data/chain", paste(opt$chainKey, "over.chain.gz", sep='.'))
    opt$runName  <- paste(opt$oldBuild,opt$newBuild, sep="-")

    # opt$oldDir  <- file.path(opt$annDir, opt$oldBuild, "chrom_hmm/wgEncodeBroadHmm")
    # opt$outDir  <- file.path(opt$annDir, "liftOver", "chrom_hmm/wgEncodeBroadHmm")
    opt$oldDir  <- file.path(opt$annDir, opt$oldBuild, "chrom_hmm",opt$liftSub)
    opt$outDir  <- file.path(opt$annDir, "liftOver", "chrom_hmm",opt$liftSub)
    
  } else {
    stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported pre-options local type: local_runType={par$local_runType}!{RET}{RET}"))
  }
  
  # opt$fresh <- TRUE
  opt$fresh <- FALSE
  
  opt$verbose <- 10
  opt$verbose <- 3
  
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
    
    make_option(c("--lift_exe"), type="character", default=opt$lift_exe, 
                help="UCSC Lift Over Executable path [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),

    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    make_option(c("--oldDir"), type="character", default=opt$oldDir, 
                help="Directory to lift from (old directory) [default= %default]", metavar="character"),
    
    # Manufacturing Files:: Required
    make_option(c("--oldBuild"), type="character", default=opt$oldBuild, 
                help="Old Build Name [default= %default]", metavar="character"),
    make_option(c("--newBuild"), type="character", default=opt$newBuild, 
                help="New Build Name [default= %default]", metavar="character"),
    
    make_option(c("--chain"), type="character", default=opt$chain, 
                help="UCSC Lift Over Chain Path [default= %default]", metavar="character"),
    make_option(c("--chainKey"), type="character", default=opt$chainKey, 
                help="UCSC Lift Over Chain Key [default= %default]", metavar="character"),

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
opt_reqs <- c('outDir','oldDir','oldBuild','newBuild',
              'lift_exe','chain','chainKey',
              'Rscript','verbose')

par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir)) stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))
for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)) base::source(sfile)
cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form General Source={par$gen_src_dir}!{RET}{RET}") )

# print(opt)
# print(opt$genBuild)

opt <- program_init(name=par$prgmTag,
                    opts=opt, opt_reqs=opt_reqs, 
                    pars=par, par_reqs=par_reqs,
                    libs=TRUE,rcpp=FALSE,
                    verbose=opt$verbose,vt=3,tc=0,tt=NULL)

par_tib <- par %>%
  unlist(recursive = TRUE) %>%
  dplyr::bind_rows() %>% 
  tidyr::gather("Params", "Value")
opt_tib <- dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value")

stamp_vec <- NULL
stamp_vec <- c(stamp_vec,opt$time_org_txt)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Parse List Options
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

pTracker <- timeTracker$new(verbose=opt$verbose)

# List target bed files::
old_files <- list.files(opt$oldDir, pattern=".bed.gz", full.names=TRUE)

opt$rdsDir <- file.path(opt$annDir,opt$newBuild,"rds")
if (!dir.exists(opt$rdsDir)) dir.create(opt$rdsDir, recursive=TRUE)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                  Main::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rds_dir <- file.path(opt$annDir, opt$newBuild, "rds")
if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive=TRUE)

if (verbose>=1)
  cat(glue::glue("[{par$prgmTag}]: rds_dir={rds_dir}{RET}") )

ncbi_gene_bed <- file.path(opt$annDir,opt$newBuild, paste(opt$newBuild,"ncbi.RefSeqGenes.tsv.gz", sep='.'))
ucsc_gene_bed <- file.path(opt$annDir,opt$newBuild, paste(opt$newBuild,"ucsc.knownGene.tsv.gz", sep='.'))
ucsc_cpgs_bed <- file.path(opt$annDir,opt$newBuild, paste(opt$newBuild,"ucsc.CpG-Islands.tsv.gz", sep='.'))

bind_all_tib <- NULL
opt$bind_all <- TRUE
if (opt$bind_all) {
  ncbi_gene_tib <- NULL
  ucsc_gene_tib <- NULL
  ucsc_cpgs_tib <- NULL
  
  ncbi_gene_tib <- load_ncbi_gene(file=ncbi_gene_bed, verbose=opt$verbose)
  ucsc_gene_tib <- load_ucsc_gene(file=ucsc_gene_bed, verbose=opt$verbose)
  ucsc_cpgs_tib <- load_ucsc_cpgs(file=ucsc_cpgs_bed, verbose=opt$verbose)
  
  bind_all_tib <- bind_all_tib %>%
    dplyr::bind_rows(ncbi_gene_tib) %>%
    dplyr::bind_rows(ucsc_gene_tib) %>%
    dplyr::bind_rows(ucsc_cpgs_tib) %>%
    utils::type.convert() %>% 
    dplyr::mutate(across(where(is.factor),  as.character) )

  if (opt$verbose>=1) {
    cat(glue::glue("[{par$prgmTag}]: bind_all_tib={RET}") )
    print(bind_all_tib)
  }
} else {
  ncbi_gene_grs <- NULL
  ucsc_gene_grs <- NULL
  ucsc_cpgs_grs <- NULL
  
  ncbi_gene_grs <- load_ncbi_gene(file=ncbi_gene_bed, grs=TRUE, verbose=opt$verbose)
  ucsc_gene_grs <- load_ucsc_gene(file=ucsc_gene_bed, grs=TRUE, verbose=opt$verbose)
  ucsc_cpgs_grs <- load_ucsc_cpgs(file=ucsc_cpgs_bed, grs=TRUE, verbose=opt$verbose)
  
  ncbi_gene_rds <- file.path(opt$rdsDir, paste(opt$newBuild,"ncbi.RefSeqGenes.rds", sep='.'))
  ucsc_gene_rds <- file.path(opt$rdsDir, paste(opt$newBuild,"ucsc.knownGene.rds", sep='.'))
  ucsc_cpgs_rds <- file.path(opt$rdsDir, paste(opt$newBuild,"ucsc.CpG-Islands.rds", sep='.'))
  
  readr::write_rds(ncbi_gene_grs, ncbi_gene_rds, compress="gz")
  readr::write_rds(ucsc_gene_grs, ucsc_gene_rds, compress="gz")
  readr::write_rds(ucsc_cpgs_grs, ucsc_cpgs_rds, compress="gz")
}

tar_cnt <- length(old_files)
for (ii in c(1:tar_cnt)) {
  old_file <- old_files[ii]
  old_name <- basename(old_file) %>% 
    stringr::str_remove('.gz$') %>%
    stringr::str_remove('.bed$')
  new_file <- file.path(opt$outDir, paste0(old_name,'-lift-',opt$chainKey,'.map.bed'))
  non_file <- file.path(opt$outDir, paste0(old_name,'-lift-',opt$chainKey,'.non.bed'))
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Processing({ii}) old_file={old_file}...{RET}"))
  
  cmd <- glue::glue("{opt$lift_exe} {old_file} {opt$chain} {new_file} {non_file}")
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: cmd='{cmd}'.{RET}"))
  system(cmd)
  
  cmd <- glue::glue("gzip -f {new_file}")
  system(cmd)
  cmd <- glue::glue("gzip -f {non_file}")
  system(cmd)
  
  rds_file <- file.path(opt$rdsDir, paste0(opt$newBuild,'.',old_name,'-lift-',opt$chainKey,'.map.rds'))
  new_file <- paste(new_file, 'gz', sep='.')
  
  bed_cols <-
    cols(
      chr  = col_character(),
      beg  = col_integer(),
      end  = col_integer(),
      key  = col_character(),
      scr  = col_integer(),
      srd  = col_character(),
      beg2 = col_integer(),
      end2 = col_integer(),
      cols = col_character()
    )
  
  cur_tib <- 
    readr::read_tsv(new_file, 
                    col_names=names(bed_cols$cols), 
                    col_types=bed_cols) %>%
    dplyr::mutate(
      class_str=stringr::str_replace_all(key, '\\/',"_"),
      class=stringr::str_remove(class_str, '_.*$') %>% as.integer(),
      name2=stringr::str_remove(class, '^[0-9]+_')) %>%
    dplyr::group_by(class) %>%
    dplyr::mutate(rank=dplyr::row_number(),
                  source="Chrom_HMM",
                  source2=opt$oldBuild,
                  name=paste(old_name, sep='_'),
                  unq_key=paste(name2,rank, sep='_')) %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::all_of(names(bind_all_tib))) %>%
    dplyr::mutate(class=as.character(class))

  if (opt$bind_all) {

    bind_all_tib <- bind_all_tib %>%
      dplyr::bind_rows(cur_tib)
    
  } else {
    # Create basic Genomic Regions RDS as well::
    bed_grs <- 
      GRanges(
        seqnames = Rle(bed_tib$chr),
        # strand = Rle(stringr::str_sub(bed_tib$bsp_srd, 1,1)),
        
        Name=bed_tib$Name,
        Name2=bed_tib$Name2,
        Class=bed_tib$Class,
        
        IRanges(start = bed_tib$beg,
                end   = bed_tib$end,
                names=paste(bed_tib$Name)
        )
      )
    
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Writing RDS={rds_file}...{RET}"))
    readr::write_rds(bed_grs,rds_file, compress="gz")
  }
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Done.{RET}{RET}"))
  
  # break
}

if (opt$bind_all) {
  
  core_rds <- file.path(rds_dir, paste(opt$newBuild,"core-annotation.rds", sep='.'))
  core_csv <- file.path(rds_dir, paste(opt$newBuild,"core-annotation.csv.gz", sep='.'))
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Writing core_csv={core_csv}...{RET}"))
  readr::write_csv(bind_all_tib, core_csv)
  
  if (opt$verbose>=1)
    cat(glue::glue("[{par$prgmTag}]: Writing core_rds={core_rds}...{RET}"))
  readr::write_rds(bind_all_tib, core_rds, compress="gz")
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
