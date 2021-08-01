
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Workflow Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Simple place to store run defaults::
#
get_run_defaults = function(ver = "1.0",
                            fresh = FALSE,
                            
                            genome_build,
                            cgn_seq_dir,
                            cgn_bed_dir,
                            canonical_cgn_dir,
                            canonical_cgn_csv,
                            
                            verbose=0,vt=3,tc=1,tt=NULL,
                            funcTag='get_run_defaults') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}                 ver={ver}.{RET}"))
    cat(glue::glue("{mssg}               fresh={fresh}.{RET}"))
    cat(glue::glue("{mssg}        genome_build={genome_build}.{RET}"))
    cat(glue::glue("{mssg}         cgn_seq_dir={cgn_seq_dir}.{RET}"))
    cat(glue::glue("{mssg}         cgn_bed_dir={cgn_bed_dir}.{RET}"))
    cat(glue::glue("{mssg}   canonical_cgn_dir={canonical_cgn_dir}.{RET}"))
    cat(glue::glue("{mssg}   canonical_cgn_csv={canonical_cgn_csv}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  etime   <- 0
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  stopifnot( dir.exists( cgn_seq_dir ) )
  stopifnot( dir.exists( cgn_bed_dir ) )
  stopifnot( dir.exists( canonical_cgn_dir ) )

  if (ver == "1.0") {
    # ret_tib <- tibble::tribble(
    #   ~name,     ~val,      ~func,
    #   "ids_key", "Prb_Key", "all",
    #   
    # )

    # Field (key) Parameters:: general
    ret_dat$ids_key <- "Prb_Key"
    ret_dat$unq_key <- "Prb_Key_Unq"
    
    ret_dat$add_key <- "Address"
    ret_dat$din_key <- "Ord_Din"
    ret_dat$des_key <- "Ord_Des"
    ret_dat$map_key <- "Ord_Map"
    ret_dat$prb_key <- "Ord_Prb"
    
    ret_dat$bsp_srd <- "Bsp_FR"
    ret_dat$bsp_cos <- "Bsp_CO"
    ret_dat$pos_key <- "Bsp_Pos"
    ret_dat$chr_key <- "Bsp_Chr"
    
    ret_dat$Cgn_Int <- "Cgn_Int"
    ret_dat$Can_Cgn <- "Can_Cgn"
    ret_dat$Ord_Cgn <- "Ord_Cgn"
    ret_dat$Bsp_Cgn <- "Bsp_Cgn"
    ret_dat$Imp_Cgn <- "Imp_Cgn"
    
    ret_dat$out_col <- c(ret_dat$ids_key, ret_dat$add_key, ret_dat$des_key,
                         ret_dat$din_key, ret_dat$map_key, ret_dat$prb_key)
    ret_dat$unq_col <- c(ret_dat$din_key, ret_dat$map_key, ret_dat$Cgn_Int)
    
    # Default run parameters by workflow::
    ret_dat$bsp_full   <- FALSE
    ret_dat$bsp_sort   <- TRUE
    ret_dat$bsp_light  <- TRUE
    ret_dat$bsp_merge  <- FALSE
    ret_dat$bsp_suffix <- "cgn.min.txt.gz"
    
    ret_dat$cgn_merge  <- FALSE
    ret_dat$cgn_join   <- "inner"

    ret_dat$seq_idxA      <- 1
    ret_dat$seq_idxB      <- 1
    ret_dat$seq_suffix    <- "probe-subseq"
    ret_dat$seq_pattern_U <- "-probe_U49_cgn-table.tsv.gz"
    ret_dat$seq_pattern_M <- "-probe_M49_cgn-table.tsv.gz"
    
    ret_dat$cgn_bed_tsv <- 
      file.path(cgn_bed_dir, paste(genome_build, ret_dat$bsp_suffix, sep="."))
    ret_dat$canonical_cgn_csv <- file.path(canonical_cgn_dir, canonical_cgn_csv)
    
    stopifnot( file.exists( ret_dat$cgn_bed_tsv ) )
    stopifnot( file.exists( ret_dat$canonical_cgn_csv ) )
    
    ret_dat$re_load <- TRUE
    if (fresh) ret_dat$re_load <- FALSE
  } else {
    stop(glue::glue("{RET}{mssg} ERROR: Only supports version 1.0 NOT ",
                    "{ver}!{RET2}"))
    return(ret_tib)
  }

  ret_key <- glue::glue("ret-FIN({funcTag})")
  ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  ret_dat
}

# End of file
