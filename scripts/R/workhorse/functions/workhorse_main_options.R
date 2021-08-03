
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("optparse", quietly=TRUE) ))
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

COM  <- ","
TAB  <- "\t"
RET  <- "\n"
RET2 <- "\n\n"
BNG  <- "|"
BRK  <- paste0("# ",
               paste(rep("-----",6),collapse=" "),"|",
               paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='template_func') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}   funcTag={funcTag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Get Program Options Defaults::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

program_default_options = function(verbose=3, vt=6,tc=1,tt=NULL,
                                   funcTag='program_default_options') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  
  etime   <- 0
  ret_cnt <- 0
  
  opts <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Run Time Version Options::
  #                       Platform, Genome Build, etc
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opts$run_name     <- NULL
  opts$platform     <- NULL
  opts$version      <- NULL
  opts$genome_build <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Run Time User Input Directories::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opts$out_dir <- NULL
  
  opts$ord_dir <- NULL
  opts$mat_dir <- NULL
  opts$aqp_dir <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Run Time User Input Files::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opts$ord_csv <- NULL
  opts$mat_tsv <- NULL
  opts$aqp_tsv <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Run Time User Input Executable(s)::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # opts$Rscript   <- NULL
  # opts$bsmap_opt <- "-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
  # opts$bsmap_dir <- NULL
  # opts$bsmap_exe <- NULL
  # opts$align_chroms <- FALSE
  # 
  # # Set Docker Defaults::
  # bsmap_dir <- '/repo/bsmap-2.90'
  # if ( dir.exists(bsmap_dir ) ) opts$bsmap_dir <- bsmap_dir
  # bsmap_exe <- file.path(bsmap_dir, 'bsmap')
  # if ( file.exists(bsmap_exe) )  opts$bsmap_exe <- bsmap_exe
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Pre-defined Static Data Directories::
  #            improbe, Annotation, Genomic, Manifest, Validation Idats
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # opts$imp_dir  <- NULL
  # opts$ann_dir  <- NULL
  # opts$gen_dir  <- NULL
  # opts$man_dir  <- NULL
  # opts$idat_dir <- NULL
  # 
  # opts$tag_seq_dir <- NULL
  # opts$tag_map_dir <- NULL
  # opts$bsp_map_dir <- NULL
  # 
  # opts$canonical_cgn_dir <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                  Pre-defined Static External File Options::
  #                   Manifest, Controls, Design Coordinates
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # opts$sesame_manfiest_dat   <- NULL
  # opts$sesame_manifest_csv   <- NULL
  # opts$genome_manifest_csv   <- NULL
  # 
  # opts$sesame_controls_csv   <- NULL
  # opts$genome_controls_csv   <- NULL
  # opts$noob_controls_csv     <- NULL
  # 
  # opts$source_coordinate_csv <- NULL
  # 
  # opts$tag_map_tsv <- "GRCh37.chr-pos-srd.slim.cgn-sorted.txt.gz"
  # opts$bsp_map_tsv <- "GRCh37.chr-pos-srd.slim.pos-sorted.txt.gz"
  # 
  # opts$canonical_cgn_csv <- "canonical.cgn-top-grp.csv.gz"
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Run Time File Options:: Time Stamps
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  opts$time_org_txt <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                          Run Time Mode Options::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  opts$single    <- FALSE
  opts$parallel  <- FALSE
  opts$cluster   <- FALSE
  
  opts$trackTime <- NULL
  opts$fresh     <- FALSE
  opts$reload    <- FALSE
  
  opts$verbose   <- verbose
  
  ret_cnt <- opts %>% length()
  
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  opts
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Get Program Options List::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

program_options = function(verbose=3, vt=3,tc=1,tt=NULL,
                           funcTag='program_options') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")

  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  # if (verbose>=vt+2) {
  #   cat(glue::glue("{RET}"))
  #   cat(glue::glue("{mssg} Function Parameters::{RET}"))
  #   cat(glue::glue("{mssg}   funcTag={funcTag}.{RET}"))
  #   cat(glue::glue("{RET}"))
  # }
  # 
  etime   <- 0
  ret_cnt <- 0
  # 
  opt <- program_default_options()
  
  option_list = list(
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Run Time Version Options::
    #                       Platform, Genome Build, etc
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Run Parameters::
    optparse::make_option(
      c("--run_name"), type="character", default=opt$run_name,
      help=paste0("Run Name [default= %default]"),
      metavar="character"),
    
    # Platform/Method Options::
    optparse::make_option(
      c("--platform"), type="character", default=opt$platform,
      help=paste0("Platform (e.g. HM450, EPIC, LEGX, NZT, ",
                  "COVIC) [default= %default]"),
      metavar="character"),
    optparse::make_option(
      c("--version"), type="character", default=opt$version,
      help=paste0("Manifest Version (e.g. B0,B1,B2,B3,B4,C0) ",
                  "[default= %default]"),
      metavar="character"),
    optparse::make_option(
      c("--genome_build"), type="character", default=opt$genome_build,
      help=paste0("Genome Build (e.g. GRch36, GRCh37, GRCh38, GRCm38) ",
                  "[default= %default]"),
      metavar="character"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Run Time User Input Directories::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    optparse::make_option(
      c("--out_dir"), type="character", default=opt$out_dir,
      help=paste0("Output directory [default= %default]"),
      metavar="character"),
    
    # Manufacturing Files:: Required
    optparse::make_option(
      c("--ord_dir"), type="character", default=opt$ord_dir,
      help=paste0("Order directories. Either a single directory or one-to-one ",
                  "pairing with order files (comma seperated list) ",
                  "[default= %default]"),
      metavar="character"),
    optparse::make_option(
      c("--mat_dir"), type="character", default=opt$mat_dir,
      help=paste0("Biziprobe Match directories. Either a single directory or ",
                  "one-to-one pairing with match files (comma seperated list) ",
                  "[default= %default]"),
      metavar="character"),
    optparse::make_option(
      c("--aqp_dir"), type="character", default=opt$aqp_dir,
      help=paste0("AQP/PQC directories. Either a single directory or ",
                  "one-to-one pairing with AQP/PQC files ",
                  "(comma seperated list) ",
                  "[default= %default]"),
      metavar="character"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Run Time User Input Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Manufacturing Files:: Required
    optparse::make_option(
      c("--ord_csv"), type="character", default=opt$ord_csv,
      help=paste0("Order file(s) (comma seperated list) [default= %default]"),
      metavar="character"),
    optparse::make_option(
      c("--mat_tsv"), type="character", default=opt$mat_tsv,
      help=paste0("Biziprobe Match file(s) (comma seperated list) ",
                  "[default= %default]"),
      metavar="character"),
    optparse::make_option(
      c("--aqp_tsv"), type="character", default=opt$aqp_tsv,
      help=paste0("AQP/PQC file(s) (comma seperated list) [default= %default]"),
      metavar="character"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Run Time User Input Executable(s)::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # optparse::make_option(
    #   c("--Rscript"), type="character", default=opt$Rscript,
    #   help=paste0("Rscript path [default= %default]"),
    #   metavar="character"),
    # 
    # optparse::make_option(
    #   c("--bsmap_opt"), type="character", default=opt$bsmap_opt,
    #   help=paste0("BSMAP Options [default= %default]"),
    #   metavar="character"),
    # optparse::make_option(
    #   c("--bsmap_dir"), type="character", default=opt$bsmap_dir,
    #   help=paste0("BSMAP Executable directory path [default= %default]"),
    #   metavar="character"),
    # optparse::make_option(
    #   c("--bsmap_exe"), type="character", default=opt$bsmap_exe,
    #   help=paste0("BSMAP Executable file name [default= %default]"),
    #   metavar="character"),
    # optparse::make_option(
    #   c("--align_chroms"), action="store_true", default=opt$align_chroms,
    #   help=paste0("Boolean flag to align against individual chromosomes. ",
    #               "This provides more alignments than the best hit. ",
    #               "[default= %default]"),
    #   metavar="boolean"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Pre-defined Static Data Directories::
    #            improbe, Annotation, Genomic, Manifest, Validation Idats
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # optparse::make_option(
    #   c("--imp_dir"), type="character", default=opt$imp_dir,
    #   help=paste0("improbe data directory [default= %default]"),
    #   metavar="character"),
    # # optparse::make_option(
    # #   c("--ann_dir"), type="character", default=opt$ann_dir,
    # #   help=paste0("Annotation data directory [default= %default]"),
    # #   metavar="character"),
    # optparse::make_option(
    #   c("--gen_dir"), type="character", default=opt$gen_dir,
    #   help=paste0("Genomic data directory [default= %default]"),
    #   metavar="character"),
    # optparse::make_option(
    #   c("--man_dir"), type="character", default=opt$man_dir,
    #   help=paste0("Pre-built Manifest data directory [default= %default]"),
    #   metavar="character"),
    
    # Validation existing idats directory to confirm Addresses against::
    # optparse::make_option(
    #   c("--idat_dir"), type="character", default=opt$idat_dir,
    #   help=paste0("Validation existing idats directory ",
    #               "to confirm Addresses against. ",
    #               "CSV file(s) (comma seperated list) [default= %default]"),
    #   metavar="character"),
    
    # Pre-defined directory with files containing cg numbers to coordinates::
    #
    # optparse::make_option(
    #   c("--tag_seq_dir"), type="character", default=opt$tag_seq_dir,
    #   help=paste0("Pre-defined directory with files containing tag probe ",
    #               "sequences to cg numbers. [default= %default]"),
    #   metavar="character"),
    # 
    # optparse::make_option(
    #   c("--tag_map_dir"), type="character", default=opt$tag_map_dir,
    #   help=paste0("Pre-defined directory with files containing cg numbers ",
    #               "to coordinates. [default= %default]"),
    #   metavar="character"),
    # 
    # # Pre-defined directory with files containing coordinates to cg numbers::
    # optparse::make_option(
    #   c("--tag_map_tsv"), type="character", default=opt$tag_map_tsv,
    #   help=paste0("Pre-defined file names containing cg number mappings ",
    #               "to coordinates. [default= %default]"),
    #   metavar="character"),
    
    
    # Pre-defined directory with files containing coordinates to cg numbers::
    #
    # optparse::make_option(
    #   c("--bsp_map_dir"), type="character", default=opt$bsp_map_dir,
    #   help=paste0("Pre-defined directory with files containing coordinates ",
    #               "to cg numbers. [default= %default]"),
    #   metavar="character"),
    
    # Pre-defined directory with files containing coordinates to cg numbers::
    # optparse::make_option(
    #   c("--bsp_map_tsv"), type="character", default=opt$bsp_map_tsv,
    #   help=paste0("Pre-defined file names containing coordinates mappings ",
    #               "to cg numbers. [default= %default]"),
    #   metavar="character"),

    # Pre-defined directory with file containing canonical cg number assignments::
    #
    # optparse::make_option(
    #   c("--canonical_cgn_dir"), type="character", default=opt$canonical_cgn_dir,
    #   help=paste0("Pre-defined canonical cg-numbers file used for cg number ",
    #               "resolution assignment. Directory path not file name(s)! ",
    #               "CSV file(s) (comma seperated list) [default= %default]"),
    #   metavar="character"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Pre-defined Static External File Options::
    #                   Manifest, Controls, Design Coordinates
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Pre-defined manifest(s) to be re-built and/or added to new manifest
    #  from Sesame Repo::
    # optparse::make_option(
    #   c("--sesame_manfiest_dat"), type="character",
    #   default=opt$sesame_manfiest_dat,
    #   help=paste0("Sesame Manifest(s) to be re-built and/or added to ",
    #               "new manifest from Sesame Repo. ",
    #               "Example = 'HM450.hg19.manifest,EPIC.hg19.manifest",
    #               "CSV file(s) (comma seperated list) [default= %default]"),
    #   metavar="character"),
    # 
    # # Pre-defined manifest(s) to be re-built and/or added to new manifest::
    # optparse::make_option(
    #   c("--sesame_manifest_csv"), type="character",
    #   default=opt$sesame_manifest_csv,
    #   help=paste0("Sesame Manifest(s) to be re-built and/or added ",
    #               "to new manifest. Probe Seq required! ",
    #               "CSV file(s) (comma seperated list) [default= %default]"),
    #   metavar="character"),
    # optparse::make_option(
    #   c("--genome_manifest_csv"), type="character",
    #   default=opt$genome_manifest_csv,
    #   help=paste0("Genome Studio Manifest(s) to be re-built and/or ",
    #               "added to new manifest. Probe Seq required! ",
    #               "CSV file(s) (comma seperated list) [default= %default]"),
    #   metavar="character"),
    
    # Pre-defined manifest control(s) to be added to new manifest::
    # optparse::make_option(
    #   c("--sesame_controls_csv"), type="character",
    #   default=opt$sesame_controls_csv,
    #   help=paste0("Sesame Pre-defined manifest control(s)  ",
    #               "to be added to new manifest. ",
    #               "CSV file(s) (comma seperated list) [default= %default]"),
    #   metavar="character"),
    # optparse::make_option(
    #   c("--genome_controls_csv"), type="character",
    #   default=opt$genome_controls_csv,
    #   help=paste0("Genome Studio Pre-defined manifest control(s)  ",
    #               "to be added to new manifest. ",
    #               "CSV file(s) (comma seperated list) [default= %default]"),
    #   metavar="character"),
    # 
    # # Pre-defined noob-masked control(s) to be added to new manifest::
    # optparse::make_option(
    #   c("--noob_controls_csv"), type="character",
    #   default=opt$noob_controls_csv,
    #   help=paste0("Noob-Masked Pre-defined control(s) ",
    #               "to be added to new manifest. ",
    #               "CSV file(s) (comma seperated list) [default= %default]"),
    #   metavar="character"),
    
    # Original source design file used for canonical position selection::
    # optparse::make_option(
    #   c("--source_coordinate_csv"), type="character",
    #   default=opt$source_coordinate_csv,
    #   help=paste0("Original source design file used for canonical ",
    #               "position selection. ",
    #               "CSV file(s) (comma seperated list) [default= %default]"),
    #   metavar="character"),
    
    # Pre-defined file names with canonical cg number assignments::
    # optparse::make_option(
    #   c("--canonical_cgn_csv"), type="character",
    #   default=opt$canonical_cgn_csv,
    #   help=paste0("Pre-defined canonical cg-numbers file used for cg number ",
    #               "resolution assignment. File(s) name, not path(s)! ",
    #               "CSV file(s) (comma seperated list) [default= %default]"),
    #   metavar="character"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Run Time File Options:: Time Stamps
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    optparse::make_option(
      c("--time_org_txt"), type="character", default=opt$time_org_txt,
      help=paste0("Unused variable time_org_txt [default= %default]"),
      metavar="character"),
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Run Time Mode Options::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Process Parallel/Cluster Parameters::
    optparse::make_option(
      c("--single"), action="store_true", default=opt$single,
      help=paste0("Boolean variable to run a single sample on a single-core ",
                  "[default= %default]"),
      metavar="boolean"),
    optparse::make_option(
      c("--parallel"), action="store_true", default=opt$parallel,
      help=paste0("Boolean variable to run parallel on multi-core ",
                  "[default= %default]"),
      metavar="boolean"),
    optparse::make_option(
      c("--cluster"), action="store_true", default=opt$cluster,
      help=paste0("Boolean variable to run jobs on cluster by chip ",
                  "[default= %default]"),
      metavar="boolean"),
    
    # Run=time Options::
    optparse::make_option(
      c("--trackTime"), action="store_true", default=opt$trackTime,
      help=paste0("Boolean variable tack run times [default= %default]"),
      metavar="boolean"),
    optparse::make_option(
      c("--fresh"), action="store_true", default=opt$fresh,
      help=paste0("Boolean variable to run a fresh build [default= %default]"),
      metavar="boolean"),
    optparse::make_option(
      c("--reload"), action="store_true", default=opt$reload,
      help=paste0("Boolean variable reload intermediate files (for testing). ",
                  "[default= %default]"),
      metavar="boolean"),
    
    # Verbosity level::
    optparse::make_option(
      c("-v", "--verbose"), type="integer", default=opt$verbose,
      help=paste0("Verbosity level: 0-5 (5 is very verbose) ",
                  "[default= %default]"),
      metavar="integer")
  )
  
  opt_parse <- optparse::OptionParser(option_list = option_list)
  opt = optparse::parse_args(opt_parse)
  
  ret_cnt <- opt %>% length()
  
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  opt
}

# End of file
