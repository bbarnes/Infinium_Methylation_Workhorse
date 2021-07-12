
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Big Science::
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
par$prgmDir <- 'annotation'
par$prgmTag <- 'annotation_to_workhorse_bed'
cat(glue::glue("[{par$prgmTag}]: Starting; {par$prgmTag}.{RET}{RET}"))


# Executables::
opt$Rscript <- NULL

# Run Parameters::
opt$runName    <- NULL
opt$Species    <- NULL

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
  opt$parallel <- TRUE
  
  opt$verbose <- 10
  opt$verbose <- 3
  
  opt$fresh <- FALSE
  opt$fresh <- TRUE
  
  opt$genBuild <- 'GRCh37'
  opt$platform <- 'methylation'
  opt$Species  <- "Human"
  opt$version  <- 'v1'
  opt$version  <- 'v2'
  
  opt$runName <- paste(opt$platform,opt$Species,opt$genBuild,opt$version, sep='-')
  
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
    
    # Run Parameters::
    make_option(c("--runName"), type="character", default=opt$runName, 
                help="Run Name [default= %default]", metavar="character"),
    make_option(c("--Species"), type="character", default=opt$Species, 
                help="Target Species Name [default= %default]", metavar="character"),
    make_option(c("--pre_man_csv"), type="character", default=opt$pre_man_csv, 
                help="Previously defined manifest [default= %default]", metavar="character"),

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
opt_reqs <- c('outDir','Species',
              'genBuild','platform','version',
              'Rscript','verbose')

par$gen_src_dir <- file.path(par$scrDir, 'functions')
if (!dir.exists(par$gen_src_dir))
  stop(glue::glue("[{par$prgmTag}]: General Source={par$gen_src_dir} does not exist!{RET}"))

for (sfile in list.files(path=par$gen_src_dir, pattern='.R$', 
                         full.names=TRUE, recursive=TRUE)) base::source(sfile)
if (opt$verbose>=0)
  cat(glue::glue("[{par$prgmTag}]: Done. Loading Source Files form ",
                 "General Source={par$gen_src_dir}!{RET}{RET}") )

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

pTracker <- timeTracker$new()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Pre-processing:: Run Time:: Output Directories
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run <- NULL

# run$manDir <- file.path(opt$outDir, 'man')
# if (!dir.exists(run$manDir)) dir.create(run$manDir, recursive=TRUE)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       6.0 Annotation Conformation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
if (opt$genBuild=="GRCh37" || opt$genBuild=="GRCh38") {
  
  core_anno_dir <- file.path(par$topDir, "data/annotation", opt$genBuild)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       6.2 Load Annotation:: EPIC_CORE
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # epic_ann_trim <- c(".bed.gz$",".sorted$",".intersect",".sorted",".formatted")
  # epic_ann_trim <- c(".sorted.merged.bed.gz",
  #                    ".sorted.intersect.bed.gz",
  #                    ".intersect.bed.gz",
  #                    ".formatted.sorted.bed.gz")

  #
  # TBD::
  #  - Fix Evidence field [0-9]+
  #  - Fix gene names: Needs parsing
  #
  epic_ann_trim <- c(".bed.gz$")
  epic_ann_path <- file.path(core_anno_dir, "EPIC_CORE_CLEAN")
  epic_ann_file <- list.dirs(epic_ann_path, full.names = TRUE)[-1] 
  epic_dir_list <- as.list(epic_ann_file)
  names(epic_dir_list) <- base::basename(epic_ann_file)
  
  opt$verbose <- 3
  epic_tib_list <- NULL
  epic_int_list <- NULL
  for (source in names(epic_dir_list)) {
    file_list <- get_file_list(
      dir=epic_dir_list[[source]], pattern=epic_ann_trim[1],
      trim=epic_ann_trim,
      verbose=opt$verbose)
    
    epic_out_path <- file.path(opt$outDir, "EPIC_CORE")
    epic_tib_list <- c(epic_tib_list,
                       lapply(file_list, grs=FALSE, load_epic_anno, source=source,
                              out=epic_out_path,
                              verbose=opt$verbose)
    )
    
    # This is for actually running the code this builds::
    #
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
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       6.3 Load Annotation:: NCBI/UCSC
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  updat_out_path <- file.path(opt$outDir, "UPDATED_CORE")
  
  ncib_gene_tsv <- file.path(core_anno_dir, "NCBI", paste(opt$genBuild,"ncbi.RefSeqGenes.tsv.gz", sep='.'))
  ucsc_gene_tsv <- file.path(core_anno_dir, "UCSC", paste(opt$genBuild,"ucsc.knownGene.tsv.gz", sep='.'))
  ucsc_cpgs_tsv <- file.path(core_anno_dir, "UCSC", paste(opt$genBuild,"ucsc.CpG-Islands.tsv.gz", sep='.'))
  
  ncbi_gene_tib <- load_ncbi_gene(file=ncib_gene_tsv, out=updat_out_path, verbose=opt$verbose, tt=pTracker)
  ucsc_gene_tib <- load_ucsc_gene(file=ucsc_gene_tsv, out=updat_out_path, verbose=opt$verbose, tt=pTracker)
  ucsc_cpgs_tib <- load_ucsc_cpgs(file=ucsc_cpgs_tsv, out=updat_out_path, verbose=opt$verbose, tt=pTracker)
  
  print(epic_tib_list$ucscRefGene)
  print(ncbi_gene_tib)
  print(ucsc_gene_tib)
  print(ucsc_cpgs_tib)
  print(epic_tib_list$cpgIslands)
  
}

#
# This is for actually comparing overlaps against manifests::
#
if (FALSE) {
  
  if (FALSE) {
    opt$verbose <- 10
    
    ncbi_gene_grs <- load_ncbi_gene(file=ncib_gene_tsv, grs=TRUE, verbose=opt$verbose, tt=pTracker)
    ucsc_gene_grs <- load_ucsc_gene(file=ucsc_gene_tsv, grs=TRUE, verbose=opt$verbose, tt=pTracker)
    ucsc_cpgs_grs <- load_ucsc_cpgs(file=ucsc_cpgs_tsv, grs=TRUE, verbose=opt$verbose, tt=pTracker)
    
    # names(ucsc_gene_grs@elementMetadata) <- c("name","name2","class","source","tissue","rank")
    
    #
    # Lets actually test some of these
    #
    # can_cols <- c("name","name2","class","source","tissue","rank")
    
    can_pre_str <- "UCSC"
    ref_pre_str <- "NCBI"
    ucsc_ncbi_int <- intersect_GRS(ref=ncbi_gene_grs,can=ucsc_gene_grs, 
                                    
                                   ref_key=NULL,ref_col=NULL,ref_prefix=ref_pre_str,
                                   can_key=NULL,can_col=NULL,can_prefix=can_pre_str, 
                                   # ref_key=NULL,ref_col=NULL,ref_prefix=NULL,
                                   # can_key=NULL,can_col=NULL,can_prefix=NULL, 
                                   
                                   # can_key="ucsc_name", all_can=TRUE,
                                   # can_prefix=can_pre_str,ref_prefix=ref_pre_str, 
                                   verbose=opt$verbose, tt=pTracker)
    
    
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
    
    #
    # This is for actually comparing overlaps against manifests::
    #
    if (FALSE) {
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
  # Compare seq_cgn_tib vs. ses_man_tib
  #
  aqp_bsp_tib %>% dplyr::filter(Ord_Key %in% seq_cgn_tib$Address) %>%
    dplyr::select(Bsp_Seq,Aln_Prb)
  
  aqp_bsp_tib %>% dplyr::filter(!Bsp_Seq %in% seq_cgn_tib$Aln_Prb)
  aqp_bsp_tib %>% dplyr::filter(!Aln_Prb %in% seq_cgn_tib$Aln_Prb)
  
  aqp_bsp_tib %>% dplyr::select(Ord_Des:Aln_Key) %>% dplyr::arrange(Aln_Prb)
  seq_cgn_tib %>% dplyr::select(Address:Aln_Prb)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 3.4.3 Bind BSMAP & Seq-Match into Table::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


#
# Compare seq_cgn_tib vs. imp_des_tib for TB/CO strand names
#

#
# Compare aqp_ann_tib vs. ses_man_tib
#   NOTE: Done below doesn't exist yet...
#


#
# improbe matching testing code below::
#
if (FALSE) {
  
  if (FALSE) {
    
    imp_ext_tib <- imp_des_tib %>% 
      # dplyr::select(Seq_ID,Probe_Seq_U,Probe_Seq_M) %>% 
      dplyr::mutate(Aln_U49=stringr::str_sub(Probe_Seq_U, 1,49), 
                    Aln_M49=stringr::str_sub(Probe_Seq_M, 1,49))
    
    aqp_cgn_vec <- aqp_add_tib %>% dplyr::mutate(Ord_Cgn=stringr::str_remove(Ord_Key, "-.*$")) %>% dplyr::pull(Ord_Cgn) %>% unique()
    
    #
    # Matching only "Mat_Prb" Plus
    #
    aqp_add_des_tib6 <- dplyr::bind_rows(
      
      # UC::
      aqp_add_tib %>% 
        dplyr::filter(Ord_Des=="U") %>% 
        dplyr::mutate(Mat_Prb=Aln_Prb) %>%
        # head(n=2) %>% 
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=Probe_Seq_U), 
                          by=c("Mat_Prb")),
      
      # UO::
      aqp_add_tib %>% 
        dplyr::filter(Ord_Des=="U") %>% 
        dplyr::mutate(Mat_Prb=stringr::str_sub(Aln_Prb, 2)) %>%
        # head(n=2) %>% 
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=stringr::str_sub(Aln_U49,1)),
                          by=c("Mat_Prb")),
      
      # MC::
      aqp_add_tib %>% 
        dplyr::filter(Ord_Des=="M") %>% 
        dplyr::mutate(Mat_Prb=Aln_Prb) %>%
        # head(n=2) %>% 
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=Probe_Seq_M), 
                          by=c("Mat_Prb")),
      
      # MO::
      aqp_add_tib %>% 
        dplyr::filter(Ord_Des=="M") %>% 
        dplyr::mutate(Mat_Prb=stringr::str_sub(Aln_Prb, 2)) %>%
        # head(n=2) %>% 
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=stringr::str_sub(Aln_U49,1)),
                          by=c("Mat_Prb")),
      
      # 2C::
      aqp_add_tib %>% 
        dplyr::filter(Ord_Des=="2") %>% 
        dplyr::mutate(Mat_Prb=Aln_P49) %>%
        # head(n=2) %>% 
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=Aln_U49), 
                          by=c("Mat_Prb")),
      
      # 2O::
      aqp_add_tib %>%
        dplyr::filter(Ord_Des=="2") %>%
        dplyr::mutate(Mat_Prb=stringr::str_sub(Aln_Prb, 1,49)) %>%
        # head(n=2) %>%
        dplyr::inner_join(imp_ext_tib %>% dplyr::mutate(Mat_Prb=Aln_U49),
                          by=c("Mat_Prb"))
      
    )
    aqp_add_tib$Ord_Prb %>% unique() %>% length()
    aqp_add_des_tib6$Ord_Prb %>% unique() %>% length()
    aqp_add_tib %>% dplyr::filter( Ord_Prb %in% aqp_add_des_tib6$Ord_Prb) %>% dplyr::distinct(Ord_Prb) %>% base::nrow()
    aqp_add_tib %>% dplyr::filter(!Ord_Prb %in% aqp_add_des_tib6$Ord_Prb) %>% dplyr::distinct(Ord_Prb) %>% base::nrow()
    
    aqp_add_tib %>% dplyr::filter(!Ord_Prb %in% aqp_add_des_tib6$Ord_Prb) %>%
      dplyr::mutate(Ord_Srd=stringr::str_remove(Ord_Key, "^.*-")) %>%
      dplyr::group_by(Ord_Idx,Ord_Srd,Ord_Des) %>%
      dplyr::summarise(Count=n(), .groups="drop") %>% print(n=1000)
    
    aqp_add_tib %>% dplyr::filter(Ord_Prb %in% aqp_add_des_tib6$Ord_Prb) %>%
      dplyr::mutate(Ord_Srd=stringr::str_remove(Ord_Key, "^.*-")) %>%
      dplyr::group_by(Ord_Idx,Ord_Srd,Ord_Des) %>%
      dplyr::summarise(Count=n(), .groups="drop") %>% print(n=1000)
    
    # aqp_add_des_tib6 %>% dplyr::mutate(Ord_Cgn=stringr::str_remove(Ord_Key, "-.*$")) %>% dplyr::pull(Ord_Cgn) %>% unique()
    aqp_add_des_tib6 %>% dplyr::mutate(Ord_Cgn=stringr::str_remove(Ord_Key, "-.*$")) %>%
      dplyr::filter(!Ord_Cgn %in% aqp_cgn_vec)
    
    #
    #
    # CONCLUSION: WE GET ALL THE CGN's!!!
    #
    #
  }
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

if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Binding all annotation...{RET}"))

ann_int_tib <- ann_int_list %>% 
  dplyr::bind_rows() %>%
  dplyr::arrange(IlmnID)

ann_int_sum <- ann_int_tib %>% 
  dplyr::group_by(source,class,tissue) %>% 
  dplyr::summarise(Count=n(), .groups="drop")
ann_int_sum %>% print(n=base::nrow(ann_int_sum))

if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Writing all annotation(CSV)={run$ann_int_csv}...{RET}"))
readr::write_csv(ann_int_tib, run$ann_int_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, verbose=opt$verbose,vt=1,tc=1,tt=pTracker)

# End of file
