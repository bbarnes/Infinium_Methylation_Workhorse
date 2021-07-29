
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           imGenomes Functions::
#            Similar to iGenomes, but for infinium methylation 
#                          Hence the lower-case im....
#
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
#                            Locate and Load all 
#                             Normal, SNP and 
#                   Pre-Bisulfite Converted Genomes (bsc)::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_imGenomes_table = function(dir, genome_build, ret_list = TRUE,
                                verbose=0,vt=3,tc=1,tt=NULL,
                                funcTag='load_imGenomes_table') {
  
  tabs <- paste0(rep(TAB, tc), collapse='')
  mssg <- glue::glue("[{funcTag}]:{tabs}")
  
  if (verbose>=vt) cat(glue::glue("{mssg} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("{mssg} Function Parameters::{RET}"))
    cat(glue::glue("{mssg}            dir={dir}.{RET}"))
    cat(glue::glue("{mssg}       ret_list={ret_list}.{RET}"))
    cat(glue::glue("{mssg}   genome_build={genome_build}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_dat <- NULL
  
  stime <- base::system.time({
    
    gen_dir <- file.path(dir, genome_build,"Sequence/WholeGenomeFasta")
    gen_pattern <- paste0(genome_build,".genome.*.fa.gz$")
    
    if (verbose>=vt) {
      cat(glue::glue("{mssg} Using Genome Build={genome_build}...{RET}"))
      cat(glue::glue("{mssg} Using Genome Build Directory={gen_dir}...{RET}"))
      cat(glue::glue("{mssg} Searching with Pattern={gen_pattern}...{RET}"))
    }
    fas_list <- list.files(gen_dir, pattern=gen_pattern, full.names=TRUE)
    fas_count <- fas_list %>% length()
    
    if (verbose>=vt) {
      cat(glue::glue("{mssg} Found {fas_count} Fasta File(s)={RET}"))
      print(fas_list)
    }

    ret_tib <- fas_list %>% 
      tibble::as_tibble() %>% 
      purrr::set_names(c("Path")) %>%
      dplyr::mutate(Base_Name=base::basename(Path) %>% stringr::str_remove(".fa.gz$")) %>%
      dplyr::mutate(Unq_ID=stringr::str_remove(Base_Name, paste(genome_build,"genome.",sep=".")), 
                    Unq_ID=stringr::str_replace(Unq_ID,"^F",paste(genome_build,"NCBI.dna.F", sep='.') ), 
                    Unq_ID=stringr::str_replace(Unq_ID,"^R",paste(genome_build,"NCBI.dna.R", sep='.') ),
                    Unq_ID=stringr::str_replace(
                      Unq_ID, paste(genome_build,"genome$", sep='.'), paste(genome_build,"NCBI.dna.FCN", sep='.') ),
                    Unq_ID=stringr::str_replace(Unq_ID,"dbSNP-151.iupac$", "dbSNP-151.iupac.FCN"), 
                    Unq_ID=stringr::str_replace(Unq_ID,"dbSNP-151.iupac",paste0(genome_build,".dbSNP-151.snp")) ) %>%
      dplyr::select(-Base_Name)

    ret_key <- glue::glue("mid-formatting({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)

    ret_tib <- ret_tib %>% 
      tidyr::separate(Unq_ID, into=c("Genome_Build","Source","Alphabet","Genome_Key"), sep="\\.") %>%
      tidyr::separate(Genome_Key, into=c("Strand_FR","Strand_CO","Strand_BSC"), sep=c(1,2), remove=FALSE) %>%
      dplyr::mutate(Genome_Key=paste(Genome_Key,Alphabet, sep="_")) %>%
      dplyr::mutate(
        Alphabet_Int=dplyr::case_when(
          Alphabet=="dna" ~ 0,
          Alphabet=="snp" ~ 1,
          TRUE ~ 3,
        ) %>% as.integer(),
        Strand_BSC_Int=dplyr::case_when(
          Strand_BSC=="N" ~ 0,
          Strand_BSC=="U" ~ 1,
          Strand_BSC=="M" ~ 2,
          Strand_BSC=="D" ~ 3,
          TRUE ~ 4,
        ) %>% as.integer(),
        Strand_CO_Int=dplyr::case_when(
          Strand_CO=="C" ~ 0,
          Strand_CO=="O" ~ 1,
          TRUE ~ 2,
        ) %>% as.integer(),
        Strand_FR_Int=dplyr::case_when(
          Strand_FR=="F" ~ 0,
          Strand_FR=="R" ~ 1,
          TRUE ~ 2
        ) %>% as.integer()
      ) %>% dplyr::arrange(Alphabet_Int,
                           Strand_BSC_Int,
                           Strand_CO_Int,
                           Strand_FR_Int)
    
    if (ret_list) ret_dat <- ret_tib %>% split(f=ret_tib$Genome_Key)

    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt=vt+4,tc=tc+1, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "{mssg} Done; Count={ret_cnt}; elapsed={etime}.{RET2}{tabs}{BRK}{RET2}"))
  
  if (ret_list) return(ret_dat)
  
  ret_tib
}







# End of file
