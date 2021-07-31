
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("plyr")) )
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))

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
#                             Noob Probe_ID Masking::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_sesame_repo_address = function(
  name,
  
  add_decoy = FALSE,
  add_masks = FALSE,
  normalize = TRUE,
  old_cols = c("Probe_ID", "seqnames", "start", "end", "strand", "designType", 
                "channel", "nextBase", "nextBaseRef", "probeType", "gene", "
                 gene_HGNC"),
  new_cols = c("Probe_ID", "Chromosome", "Coordinate", "CoordinateG", 
                "Strand_FR", "Infinium_Design_Type", "Color", "Prb_Nxb", 
                "Prb_Nxb_Ref", "Prb_Din", "gene", "gene_HGNC"),
  
  verbose=0,vt=3,tc=1,tt=NULL,
  funcTag='load_sesame_repo_manifest') {
  
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
    
    # Cache Sesame Manifest Data::
    data_key <- name %>% stringr::str_remove("\\..*$")
    sesameData::sesameDataCache(data_key)
    
    man_tib <- sesameData::sesameDataGet( name ) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var="Probe_ID") %>% 
      tibble::as_tibble()
    man_key <- glue::glue("man-tib({funcTag})")
    man_cnt <- print_tib(man_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=man_key)
    
    # Split Probes into A/B
    prb_key <- "A"
    oldA_cols <- 
      man_tib %>% 
      dplyr::select(dplyr::ends_with(paste0("_",prb_key))) %>% 
      dplyr::select(-dplyr::starts_with("wDecoy_")) %>% 
      names()
    newA_cols <- oldA_names %>% 
      stringr::str_remove(paste0("_",prb_key,"$")) %>% 
      stringr::str_to_title() %>%
      paste("Prb",., sep="_") %>%
      stringr::str_replace("Prb_Address", "Address") %>%
      stringr::str_replace("Prb_Probeseq", "Prb_Seq") %>%
      stringr::str_replace("Prb_Chrm", "Prb_Chr")
    
    mask_cols <- NULL
    if (add_masks) man_tib %>% dplyr::select(dplyr::starts_with("MASK_"))
    
    selA_cols <- c(old_cols, oldA_cols, mask_cols)
    
    
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
#                             Manifest I/O::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# epic_gst_dat <- 
#   load_genome_studio_manifest(file = opt$genome_manifest_csv,
#                               
#                               load_clean    = TRUE,
#                               load_controls = TRUE,
#                               write_clean   = TRUE,
#                               overwrite     = TRUE,
#                               ret_data      = TRUE,
#                               
#                               verbose = opt$verbose, tt = pTracker)
#

load_genome_studio_manifest = function(file,
                                       
                                       load_clean    = TRUE,
                                       load_controls = FALSE,
                                       write_clean   = TRUE,
                                       overwrite     = FALSE,
                                       ret_data      = FALSE,
                                       
                                       verbose=0, vt=3,tc=1,tt=NULL,
                                       funcTag='load_genome_studio_manifest') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}          file={file}.{RET}"))
    cat(glue::glue("{mssg}    load_clean={load_clean}.{RET}"))
    cat(glue::glue("{mssg}   write_clean={write_clean}.{RET}"))
    cat(glue::glue("{mssg}     overwrite={overwrite}.{RET}"))
    cat(glue::glue("{mssg}      ret_data={ret_data}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  probes_tib   <- NULL
  controls_tib <- NULL
  
  stime <- base::system.time({
    
    # Set Column Repair Prefix::
    col_prefix <- "Val_"
    if (stringr::str_detect(file,"EPIC")) col_prefix <- "HM"
    
    gs_control_names <- 
      c("Address","Control_Group","GS_Color","Control_Type")
    
    analytical_suffix_csv <- ".analytical.csv.gz"
    analytical_suffix_rds <- ".analytical.rds"
    control_suffix_csv    <- ".controls.csv.gz"
    
    # if (!R.utils::isGzipped(file, ".gz"))
    #   file <- R.utils::gzip(file) %>% as.character()
    if (!stringr::str_ends(file,".gz")) {
      system(glue::glue("gzip {file}"))
      file <- paste0(file,".gz")
    }
    
    if (verbose>=vt+3) {
      cat(glue::glue("{mssg} Setting suffix variables::{RET}"))
      cat(glue::glue("{mssg}       control_suffix_csv=",
                     "{control_suffix_csv}.{RET}"))
      cat(glue::glue("{mssg}    analytical_suffix_csv=",
                     "{analytical_suffix_csv}.{RET}"))
      cat(glue::glue("{mssg}    analytical_suffix_rds=",
                     "{analytical_suffix_rds}.{RET}"))
      cat(glue::glue("{RET}"))
    }
    
    # If the input file is the analytical file load it!
    if (stringr::str_ends(file, analytical_suffix_csv)) {
      clean_prb_csv <- file
      clean_col_rds <- clean_basename %>%
        paste(analytical_suffix_csv)
      clean_prb_done <- paste(clean_prb_csv,'done.txt', sep='.')
      
      if (verbose>=vt+2)
        cat(glue::glue("{mssg} Found clean manifest! Will load clean ",
                       "manifest={clean_prb_csv}...{RET}"))
      
      if (file.exists(clean_col_rds) &&
          file.mtime(clean_prb_csv) <= file.mtime(clean_col_rds)) {
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Using column header={clean_col_rds}...{RET}"))
        
        probes_type_cols <- readr::read_rds(clean_col_rds)
        probes_tib  <- readr::read_csv(clean_prb_csv,
                                       col_names = names(probes_type_cols$cols),
                                       col_types = probes_type_cols)
      } else {
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Will build column header...{RET}"))
        
        probes_tib <- safe_read( clean_prb_csv, clean = TRUE, 
                                 verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
        
        # Get Column Types String::
        col_types_str <- sapply(probes_tib, typeof) %>% 
          tibble::as_tibble(rownames = "Key") %>% 
          dplyr::mutate(value=stringr::str_sub(value, 1,1)) %>% 
          dplyr::pull(value) %>% paste0(collapse = '')
        
        # Build cols() object for future loading::
        probes_type_cols <- spec(readr::read_csv(
          readr::format_csv(probes_tib), col_types = col_types_str) )
        
        if (write_clean) {
          if (verbose>=vt+2) 
            cat(glue::glue("{mssg} Writing clean manifest cols RDS...{RET}"))
          
          safe_write( probes_type_cols, file = clean_col_rds, done = TRUE,
                      verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
        }
      }
      
      if (load_controls) {
        controls_csv <- file %>% 
          stringr::str_replace(analytical_suffix_csv, control_suffix_csv)
        
        if (!file.exists(controls_csv)) {
          cat(glue::glue("{mssg} Warning: No controls file present. Try ",
                         "rebuilding genome studio manifest from scratch.{RET}"))
        } else {
          if (verbose>=vt+2)
            cat(glue::glue("{mssg} Loading controls={controls_csv}...{RET}"))
          
          controls_tib <- 
            suppressMessages(suppressWarnings( readr::read_csv( controls_csv )))
        }
      }
      
    } else {
      
      # Define clean formatted files::
      clean_basename <- file %>% 
        stringr::str_remove(".gz$") %>%
        stringr::str_remove(".csv$")
      
      clean_prb_csv <- clean_basename %>%
        paste0(analytical_suffix_csv)
      clean_col_rds <- clean_basename %>%
        paste0(analytical_suffix_rds)
      clean_prb_done <- paste(clean_prb_csv,'done.txt', sep='.')
      controls_csv   <- clean_basename %>%
        paste0(control_suffix_csv)
      
      if (verbose>=vt+3) {
        cat(glue::glue("{mssg} Setting output variables::{RET}"))
        cat(glue::glue("{mssg}    clean_basename={clean_basename}.{RET}"))
        cat(glue::glue("{mssg}      controls_csv={controls_csv}.{RET}"))
        cat(glue::glue("{mssg}     clean_col_rds={clean_col_rds}.{RET}"))
        cat(glue::glue("{mssg}     clean_prb_csv={clean_prb_csv}.{RET}"))
        cat(glue::glue("{mssg}    clean_prb_done={clean_prb_done}.{RET}"))
        cat(glue::glue("{RET}"))
      }
      
      if (load_clean && file.exists(clean_prb_csv) &&
          file.mtime(file) <= file.mtime(clean_prb_csv) &&
          file.mtime(clean_prb_csv) <= file.mtime(clean_prb_done) ) {
        
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Loading clean manifest={clean_prb_csv}...{RET}"))
        
        if (file.exists(clean_col_rds) &&
            file.mtime(clean_prb_csv) <= file.mtime(clean_col_rds)) {
          if (verbose>=vt+2)
            cat(glue::glue("{mssg} Using column header={clean_col_rds}...{RET}"))
          
          probes_type_cols <- readr::read_rds(clean_col_rds)
          
          probes_tib  <- readr::read_csv(clean_prb_csv,
                                         col_names = names(probes_type_cols$cols),
                                         col_types = probes_type_cols)
        } else {
          if (verbose>=vt+2) 
            cat(glue::glue("{mssg} Will build column header...{RET}"))
          
          probes_tib <- safe_read( clean_prb_csv, clean = TRUE, 
                                   verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
          
          # Get Column Types String::
          col_types_str <- sapply(probes_tib, typeof) %>% 
            tibble::as_tibble(rownames = "Key") %>% 
            dplyr::mutate(value=stringr::str_sub(value, 1,1)) %>% 
            dplyr::pull(value) %>% paste0(collapse = '')
          
          # Build cols() object for future loading::
          probes_type_cols <- spec(readr::read_csv(
            readr::format_csv(probes_tib), col_types = col_types_str) )
          
          if (write_clean) {
            if (verbose>=vt+2) 
              cat(glue::glue("{mssg} Writing clean manifest cols RDS...{RET}"))
            
            safe_write( probes_type_cols, file = clean_col_rds, done = TRUE,
                        verbose=verbose, vt=vt+1,tc=tc+1,tt=tt )
          }
        }
        
        # Load Controls if requested
        if (load_controls) {
          if (!file.exists(controls_csv))
            cat(glue::glue("{mssg} Warning: No controls file present. Try ",
                           "rebuilding genome studio manifest from scratch.{RET}"))
          else
            controls_tib <- 
              suppressMessages(suppressWarnings( readr::read_csv( controls_csv )))
        }
        
      } else {
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Will load raw manifest={file}...{RET}"))
        
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Will build column header...{RET}"))
        
        lines_vec <- readr::read_lines(file)
        if (ret_data) ret_dat$lines <- lines_vec
        
        # Get Analytical and Control Start & End Indexes::
        beg_idx <- head(
          which( stringr::str_starts(head(lines_vec, n=20), "IlmnID" ) ), n=1)
        con_idx <- tail(
          which( stringr::str_starts(lines_vec, "\\[Controls\\]" ) ), n=1)
        end_idx <- lines_vec %>% length()
        
        if (verbose>=vt+2)
          cat(glue::glue("{mssg} Setting Indexes=({beg_idx}, {con_idx}, ",
                         "{end_idx}).{RET}"))
        
        # Get Analytical Column Names::
        probes_name_cols <- 
          lines_vec[beg_idx] %>% 
          stringr::str_split(pattern = ",", simplify = TRUE) %>% 
          as.vector() %>% 
          stringr::str_replace("^([0-9])", paste0(col_prefix, "\\$1") ) %>%
          stringr::str_remove_all("\\\\")
        
        probes_name_cnt <- probes_name_cols %>% length()
        if (verbose>=vt+3) {
          name1_str <- probes_name_cols[1]
          name2_str <- probes_name_cols[probes_name_cnt]
          out_str <- glue::glue("({name1_str} ... {name2_str})")
          cat(glue::glue("{mssg} Extracted Probe Names Vec({probes_name_cnt})",
                         " = {out_str}.{RET2}"))
        }
        
        # Build Analytical Probes Tibble
        probes_tib <- 
          lines_vec[c( (beg_idx+1):(con_idx-1) ) ] %>% # head() %>%
          tibble::as_tibble() %>% 
          tidyr::separate(value, into=c(probes_name_cols), sep=',') %>%
          clean_tibble()
        ret_key <- glue::glue("probes-tib({funcTag})")
        ret_cnt <- print_tib(probes_tib,funcTag, verbose,vt=vt+4,tc=tc+1, 
                             n=ret_key)
        
        # Get Column Types String::
        col_types_str <- sapply(probes_tib, typeof) %>% 
          tibble::as_tibble(rownames = "Key") %>% 
          dplyr::mutate(value=stringr::str_sub(value, 1,1)) %>% 
          dplyr::pull(value) %>% paste0(collapse = '')
        
        # Build cols() object for future loading::
        probes_type_cols <- spec(readr::read_csv(
          readr::format_csv(probes_tib), col_types = col_types_str) )
        if (verbose>=vt+4) {
          cat(glue::glue("{mssg} probes_type_cols={RET}"))
          probes_type_cols %>% print()
        }
        
        # Build Controls Tibble::
        controls_tib <- 
          lines_vec[c( (con_idx+1):(end_idx) ) ] %>%
          stringr::str_remove(",+$") %>%
          stringr::str_remove(",AVG$") %>%
          tibble::as_tibble() %>% 
          tidyr::separate(value, into=c(gs_control_names), sep=',') %>% 
          clean_tibble()
        ret_key <- glue::glue("controls-tib({funcTag})")
        ret_cnt <- print_tib(controls_tib,funcTag, verbose,vt=vt+4,tc=tc+1,
                             n=ret_key)
        
        if (ret_data) ret_dat$probes   <- probes_tib
        if (ret_data) ret_dat$controls <- controls_tib
        
        if (ret_data) ret_dat$probes_key_str <- probes_name_cols
        if (ret_data) ret_dat$probes_col_str <- col_types_str
        if (ret_data) ret_dat$probes_cols    <- probes_type_cols
        
        # Write clean files if requested::
        if (write_clean) {
          if (verbose>=vt+2) 
            cat(glue::glue("{mssg} Writing clean manifest files...{RET}"))
          
          out_tag <- glue::glue("write-controls-CSV({funcTag})")
          safe_write( controls_tib, file = controls_csv, done = TRUE, 
                      funcTag = out_tag, verbose=verbose, vt=vt+3,tc=tc+1,tt=tt)
          
          out_tag <- glue::glue("write-probes-CSV({funcTag})")
          safe_write( probes_tib, file = clean_prb_csv, done = TRUE,
                      funcTag = out_tag, verbose=verbose, vt=vt+3,tc=tc+1,tt=tt)
          
          out_tag <- glue::glue("write-probes-col-RDS({funcTag})")
          safe_write( probes_type_cols, file = clean_col_rds, done = TRUE,
                      funcTag = out_tag, verbose=verbose, vt=vt+3,tc=tc+1,tt=tt)
        }
      }
    }
    
    if (load_controls && is.null(controls_tib)) {
      ret_dat$probes   <- probes_tib
      ret_dat$controls <- controls_tib
      ret_data <- TRUE
    }
    
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(probes_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  if (ret_data) return(ret_dat)
  
  probes_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Noob Probe_ID Masking::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# changing default mod=100000000 to mod=999999999
# changing default prefix="nb" to prefix="cg9999"
# NOTE: 5 9's is optimal because improbe can only handle 15 characters
#  as input 2 character (cg) + 8 digits + 2 character ([TB][CO]) = 12
# However, we can reduce characters to a single leter 1+8 = 9 
#  Using 5 would put us at 14 leaving one charcter to expand 8 digits to 9
#
# None of that really matters since we'll never use these probes in improbe
#  The issue is running into your growing cg# (8 to 9 to 10) space. 
#  We'll make 5 9's
# Wait why don't we just make cgX
noob_mask = function(x, seed=21L, mod=100000000, prefix="cgBK",
                     funcTag='noob_mask') {
  
  # Clean input::
  x <- x %>% stringr:: str_remove_all("[^0-9]+") %>% as.integer()
  
  hash <- digest::digest2int(as.character(x), seed) %% mod
  
  # Waiting to see if this ever fails
  m_len <- stringr::str_length(format(mod, scientific = FALSE))
  h_len <- stringr::str_length(hash)
  
  if (h_len>=m_len) {
    stop(glue::glue("{RET}[{funcTag}]: ERROR: ",
                    "h_len({h_len}) >= m_len({m_len})!!!{RET}",
                    "hash={hash}{RET}",
                    "mod={mod}{RET}{RET}"))
    return(NULL)
  }
  
  # Now remake cg style number::
  hash <- paste0(prefix,
                 stringr::str_pad(hash, width=m_len-1,side="left", pad="0"))
  
  hash
}

noob_mask_manifest = function(tib,
                              key = "Probe_ID", 
                              out = NULL, 
                              prefix = "cg",
                              
                              verbose=0, vt=3,tc=1,tt=NULL,
                              funcTag='noob_mask_manifest') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}      key={key}.{RET}"))
    cat(glue::glue("{mssg}      out={out}.{RET}"))
    cat(glue::glue("{mssg}   prefix={prefix}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (is.null(out)) out <- key
    key_sym <- rlang::sym(key)
    out_sym <- rlang::sym(out)
    
    noob_vec <- tib %>% 
      dplyr::pull(key_sym) %>% 
      lapply( noob_mask) %>% unlist()
    
    ret_tib <- tib %>% dplyr::mutate(!!out_sym := noob_vec)
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2{tabs}{BRK}{RET2}"))
  
  ret_tib
}

# End of file
