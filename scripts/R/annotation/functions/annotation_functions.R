
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Annotation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Genomic Ranges::
suppressWarnings(suppressPackageStartupMessages( base::require("GenomicRanges",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("plyr")) )
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Method::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadUcscIslandGR = function(file,
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'template_func'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; file={file}...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
    
    #
    # Build GRange Data Structures:: CpG Islands
    #
    ret_grs$island <- GRanges(seqnames = Rle(dat_tib$chrom), 
                              IRanges(start=dat_tib$chromStart, end=dat_tib$chromEnd, names= dat_tib$name) )
    
    ret_grs$nshore <- GRanges(seqnames = Rle(dat_tib$chrom), 
                              IRanges(start=dat_tib$chromStart-2000, end=dat_tib$chromStart, names= dat_tib$name) )
    
    ret_grs$nshelf <- GRanges(seqnames = Rle(dat_tib$chrom), 
                              IRanges(start=dat_tib$chromStart-4000, end=dat_tib$chromStart-2000, names= dat_tib$name) )
    
    ret_grs$sshore <- GRanges(seqnames = Rle(dat_tib$chrom), 
                              IRanges(start=dat_tib$chromEnd, end=dat_tib$chromEnd+2000, names= dat_tib$name) )
    
    ret_grs$sshelf <- GRanges(seqnames = Rle(dat_tib$chrom), 
                              IRanges(start=dat_tib$chromEnd+2000, end=dat_tib$chromEnd+4000, names= dat_tib$name) )
    
    ret_cnt <- ret_grs %>% names %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_grs
}



# End of file
