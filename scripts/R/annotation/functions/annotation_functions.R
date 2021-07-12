
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
#                               dbSNP Method::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

mutate_chrom_seq = function(seq, pos, val,
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutate_chrom_seq'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_seq <- NULL
  stime <- system.time({
    
    pos_len <- length(pos)
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Position Count={pos_len}.{RET}"))
    
    ret_seq <- seq %>%
      stringi::stri_sub_all_replace(from=pos, to=pos, value=val) %>%
      Biostrings::DNAString()
    ret_len <- ret_seq %>% length()
    
    # ret_cnt <- ret_tib %>% base::nrow()
    # ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    ret_cnt <- ret_seq %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_seq
}

load_dbSNP_vcf = function(vcf, file=NULL, 
                          fresh=FALSE, unique=FALSE,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_dbSNP_vcf'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_sum <- NULL
  stime <- system.time({
    
    done_txt  <- paste0(file,'.done.txt')
    stamp_vec <- c(vcf,file,done_txt)
    
    if (fresh || !valid_time_stamp(stamp_vec)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Building clean snp_csv from={vcf}...{RET}"))
      
      vcf_sel <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
      vcf_col <- c("Chromosome","Coordinate","Seq_ID","AlleleA_Str","AlleleB_Str","Qual","Filter","Info")
      
      ret_tib <- fread(vcf, select=vcf_sel) %>%
        tibble::as_tibble() %>%
        purrr::set_names(vcf_col) %>% 
        dplyr::filter(stringr::str_detect(AlleleA_Str, "^[A-Z]$")) %>%
        dplyr::filter(
          stringr::str_detect(AlleleB_Str, "^[A-Z]$") | 
            stringr::str_detect(AlleleB_Str, "^[A-Z],[A-Z]$") |
            stringr::str_detect(AlleleB_Str, "^[A-Z],[A-Z],[A-Z]$")) %>%
        dplyr::mutate(
          AlleleB_Str1=stringr::str_remove(AlleleB_Str, ",.*$"),
          AlleleB_Str2=dplyr::case_when(
            stringr::str_detect(AlleleB_Str, "^[A-Z],[A-Z],[A-Z]$") ~ stringr::str_remove(AlleleB_Str,",[A-Z]$"),
            TRUE ~ AlleleB_Str),
          
          AlleleC_Iup1=paste0(AlleleA_Str, stringr::str_remove_all(AlleleB_Str1, ",")) %>% 
            mapDIs(),
          AlleleC_Iup2=paste0(AlleleA_Str, stringr::str_remove_all(AlleleB_Str2, ",")) %>% 
            mapDIs(),
          
          Chromosome=paste0("chr",stringr::str_remove(Chromosome, "^chr")),
          AlleleB_Len=stringr::str_remove_all(AlleleB_Str,",") %>% 
            stringr::str_length()
        ) %>% 
        dplyr::arrange(Chromosome,Coordinate, -AlleleB_Len)
      
      # TBD:: Add VC filtering like below::
      # ret_tib %>% dplyr::filter(stringr::str_detect(Info,"VC=SNV"))
      
      if (unique) ret_tib <- ret_tib %>%
        dplyr::distinct(Chromosome,Coordinate, .keep_all=TRUE)
      
      if (!is.null(file)) {
        safe_write(x=ret_tib,type="csv",file=file,done=TRUE,funcTag=funcTag, 
                   verbose=verbose,vt=vt,tc=tc,append=FALSE)
      }
    } else {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Loading snp_csv={file}...{RET}"))
      
      ret_tib <- suppressMessages(suppressWarnings( readr::read_csv(file) ))
    }
    
    if (verbose>=vt+4) {
      # Summary by chromosome coverage::
      ret_sum <- ret_tib %>% 
        dplyr::group_by(Chromosome) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      ret_sum %>% print(n=base::nrow(ret_sum))
      
      # Summary by SNP Allele coverage::
      ret_sum <- ret_tib %>% 
        dplyr::mutate(Probe_Type=stringr::str_sub(Seq_ID, 1,2)) %>% 
        dplyr::group_by(AlleleA_Str,AlleleB_Str2,AlleleC_Iup) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      ret_sum %>% print(n=base::nrow(ret_sum))
    }
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Manifest Intersection Method::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

manifestToAnnotation = function(tib, ann, gen, csv=NULL, key="Seq_ID",
                                tar=c("UCSC_Islands","UCSC_Genes","NCBI_Genes"),
                                verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'manifestToAnnotation'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; key={key}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Annotation   = {ann}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr} Genome Build = {gen}{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr} Targets={RET}"))
      tar %>% print()
      cat("\n")
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      4.2 Build Genomic Region Set::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # key_sym <- rlang::sym(key)
    
    # Ensure we have unique key fields::
    tib <- makeFieldUnique(
      tib=tib, field=key, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} tib={RET}"))
      print(tib)
      cat("\n")
    }
    
    man_grs <- GenomicRanges::GRanges(
      seqnames = Rle(tib$chrom), 
      IRanges(start=tib$chromStart, 
              end=tib$chromEnd, 
              names=tib$Seq_ID) )
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} man_grs={RET}"))
      print(man_grs)
      cat("\n")
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                5.1 Intersect Mapped Probes with Annotation::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # NCBI Genes::
    ncbi_mat_tib <- NULL
    if ("NCBI_Genes" %in% tar) {
      ncbi_ann_tsv <- file.path(ann,gen,paste(gen,'ncbi.RefSeqGenes.tsv.gz', sep='.') )
      
      if (!is.null(ncbi_ann_tsv) && file.exists(ncbi_ann_tsv)) {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} ncbi_ann_tsv={ncbi_ann_tsv}...{RET}"))
        
        ncbi_ann_tib <- suppressMessages(suppressWarnings( readr::read_tsv(ncbi_ann_tsv) ))
        colnames(ncbi_ann_tib)[1] <- stringr::str_remove(colnames(ncbi_ann_tib)[1], '^#')
        
        ncbi_ref_grs <- loadNcbiGeneGR(file=ncbi_ann_tsv, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        ncbi_mat_tib <- intersectGranges(man=man_grs,ref=ncbi_ref_grs,
                                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
          dplyr::left_join(dplyr::select(ncbi_ann_tib,name,name2), by=c("Gene"="name")) %>% 
          dplyr::rename(Transcript=Gene, Gene=name2) %>% 
          dplyr::select(Seq_ID, Gene,Transcript,dplyr::everything()) %>%
          dplyr::mutate(Source="NCBI")
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. ncbi.{RET}{RET}"))
      } else {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Failed to find NCBI_Genes={ncbi_ann_tsv}...{RET}"))
      }
    } else {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Will skip NCBI_Genes...{RET}"))
    }
    
    # UCSC Genes::
    gene_mat_tib <- NULL
    if ("UCSC_Genes" %in% tar) {
      gene_ann_tsv <- file.path(ann,gen,paste(gen,'ucsc.knownGene.tsv.gz', sep='.') )
      if (!is.null(gene_ann_tsv) && file.exists(gene_ann_tsv)) {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} gene_ann_tsv={gene_ann_tsv}...{RET}"))
        
        gene_ann_tib <- suppressMessages(suppressWarnings( readr::read_tsv(gene_ann_tsv) )) %>% 
          dplyr::rename(name2=proteinID)
        colnames(gene_ann_tib)[1] <- stringr::str_remove(colnames(gene_ann_tib)[1], '^#')
        
        # proteinID
        gene_ref_grs <- loadUcscGeneGR(file=gene_ann_tsv, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        gene_mat_tib <- intersectGranges(man=man_grs,ref=gene_ref_grs,
                                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
          dplyr::left_join(dplyr::select(gene_ann_tib,name,name2), by=c("Gene"="name")) %>% 
          dplyr::rename(Transcript=Gene, Gene=name2) %>% 
          dplyr::select(Seq_ID, Gene,Transcript,dplyr::everything()) %>%
          dplyr::mutate(Source="UCSC")
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. ucsc.{RET}{RET}"))
      } else {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Failed to find UCSC_Genes={ucsc_ann_tsv}...{RET}"))
      }
    } else {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Will skip UCSC_Genes...{RET}"))
    }
    
    # CpG Islands::
    cpgs_mat_tib <- NULL
    if ("UCSC_Islands" %in% tar) {
      cpgs_ann_tsv <- file.path(ann,gen,paste(gen,'ucsc.CpG-Islands.tsv.gz', sep='.') )
      if (!is.null(cpgs_ann_tsv) && file.exists(cpgs_ann_tsv)) {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} cpgs_ann_tsv={cpgs_ann_tsv}...{RET}"))
        
        cpgs_ann_tib <- suppressMessages(suppressWarnings( readr::read_tsv(cpgs_ann_tsv) )) %>% dplyr::select(-1)
        cpgs_ref_grs <- loadUcscCpgsGR(file=cpgs_ann_tsv, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        cpgs_mat_tib <- intersectGranges(man=man_grs,ref=cpgs_ref_grs,
                                         verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
          dplyr::mutate(Source="UCSC")
        
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. cpgs.{RET}{RET}"))
      } else {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Failed to find UCSC_Islands={cpgs_ann_tsv}...{RET}"))
      }
    } else {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Will skip UCSC_Islands...{RET}"))
    }
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. Mapping Annotation.{RET}{RET}"))
    
    ret_tib <- dplyr::bind_rows(ncbi_mat_tib,gene_mat_tib,cpgs_mat_tib) %>% 
      dplyr::arrange(chrom,chromStart,chromEnd) %>%
      dplyr::mutate(Seq_Rep=stringr::str_remove(Seq_ID, '^.*_'),
                    Seq_ID=stringr::str_remove(Seq_ID, '_[0-9]+$'))
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      print(ret_tib)
    }
    
    if (!is.null(csv)) {
      outDir <- base::dirname(csv)
      if (verbose>=vt) cat(glue::glue("[{par$prgmTag}]: Writing annotation mappings={csv}...{RET}"))
      readr::write_csv(ret_tib, csv)
      if (verbose>=vt) cat(glue::glue("[{par$prgmTag}]: Done. Writing annotation mappings.{RET}{RET}"))
    }
    
    # man_ana_key <- names(ret_tib)[1]
    # man_ana_beg <- names(ret_tib)[2]
    # man_ana_beg <- names(ret_tib)[2]
    # 
    # man_ana_col <- ret_tib %>% dplyr::select(-1) %>% names()
    # man_all_col <- ret_tib %>% names()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Standard EPIC Annotation Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_epic_anno = function(file,grs=FALSE, source="EPIC", tissue="All",
                          out=NULL,
                          sel_col=c("chr","beg","end","srd",
                                    "name","name2","class",
                                    "source","tissue","rank",
                                    "unq_key","evidence"),
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='load_epic_anno') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  val_col <- 
    cols(
      chr  = col_character(),
      beg  = col_integer(),
      end  = col_integer(),
      srd  = col_character(),
      unq_key  = col_character(),
      name2 = col_character()
    )
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data({source})={file}...{RET}"))
    
    file_type <- guess_file_del(file)

    if (file_type==COM) {
      ret_tib <- 
        readr::read_csv(file, col_names=names(val_col$cols), col_types=val_col)
    } else if (file_type==TAB) {
      ret_tib <- 
        readr::read_tsv(file, col_names=names(val_col$cols), col_types=val_col)
    } else if (file_type==" ") {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: file_type=Single-Space!!!{RET}{RET}"))
      return(ret_tib)
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: file_type=NULL!!!{RET}{RET}"))
      return(ret_tib)
    }
    ret_key <- glue::glue("ret-raw({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    #
    # Determine the data size field 'name2'
    #
    dat_len <- ret_tib %>% 
      dplyr::pull(name2) %>% head(n=1) %>%
      stringr::str_split(pattern="[-:]", simplify = TRUE) %>% 
      as.vector() %>% length()
    if (verbose>=vt+4)
      cat(glue::glue("[{funcTag}]:{tabsStr} dat_len={dat_len}.{RET}"))
    
    ret_tib <- ret_tib %>% # head(4) %>%
        # dplyr::select(chr,beg,end,srd,name2) %>%
        dplyr::mutate(name2=stringr::str_replace(name2, ".evd=", ":"),
                      dat_len=dat_len)
    ret_key <- glue::glue("ret-spt({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)

    # name2=transcript/chr-beg-end
    # name =gene
    #
    if (dat_len==1) {
      # class = name2 (DMR)
      # name  = chr-beg-end
      # name2 = name
      # evidence = 0 [default]
      #
      # [DONE]: Files formatted
      
      ret_tib <- ret_tib %>%
        dplyr::mutate(class=name2,
                      name2=paste(chr,beg,end, sep="-"),
                      name=paste(chr,beg,end, sep="-"),
                      evidence=as.integer(0)
        )

    } else if (dat_len==4) {
      # class = n1
      # name  = chr-beg-end
      # name2 = n2-n3-n4
      # source = UCSC
      # evidence = n5
      #   - Fix Evidence field [0-9]+
      #
      # QUESTION: collapse source from UCSC_Island to just UCSC
      #   NOTE: evidence in this case is CpG Count:
      #      Leave that for technical writing and marketing to explain...
      #
      #
      # [Done]: File modifications and directory collapse
      #
      all_col <- c("n1","n2","n3","n4","n5")

      ret_tib <- ret_tib %>%
        tidyr::separate(name2, into=all_col, sep="[:]") %>%
        dplyr::mutate(class=n1,
                      name=paste(chr,beg,end, sep="-"),
                      name2=paste(n2,n3,n4, sep="-"),
                      evidence=as.integer(n5)
        )
      
      
    } else if (dat_len==5) {
      # This should be true:
      #   all_tib %>% filter(Len==5 & (chr!=n3 | beg!=n4 | end!=n5))
      #
      # Source and Class will take some file manipulation
      #   1. TFBS_PeakSeq-based_Peaks -> TFBS_PeakSeq_based_Peaks
      #
      # class = n1 Dnase1,TFBS.../Enhacners
      # name  = n3-n4-n5
      # name2 = name
      # source = ENCODE/PHANTOM
      # evidence = n6
      #   - Fix Evidence field [0-9]+
      #
      # [Done]: Fixed TFBS name to change - to _
      # [Done]: DNase's
      # [Done]: Phantom5
      #
      all_col <- c("n1","n2","n3","n4","n5","n6")
      
      ret_tib <- ret_tib %>%
        tidyr::separate(name2, into=all_col, sep="[:]") %>%
        dplyr::mutate(class=n1,
                      name=paste(n3,n4,n5, sep="-"),
                      name2=name,
                      evidence=as.integer(n6)
        )

    } else if (dat_len==6) {
      # class = n1
      # name  = n4
      # name2 = n3
      # evidence = n6
      all_col <- c("n1","n2","n3","n4","n5","n6")
      
      ret_tib <- ret_tib %>%
        tidyr::separate(name2, into=all_col, sep="[:]") %>%
        dplyr::mutate(class=n1,
                      name=n4,
                      name2=n3,
                      evidence=as.integer(0)
        )
    } else if (dat_len==7) {
      # class = n1
      # name  = n4
      # name2 = n5
      # evidence = n6
      all_col <- c("n1","n2","n3","n4","n5","n6","n7")
      
      ret_tib <- ret_tib %>%
        tidyr::separate(name2, into=all_col, sep="[:]") %>%
        dplyr::mutate(class=n1,
                      name=n4,
                      name2=n3,
                      evidence=as.integer(n7)
        )
    } else {
      cat("Unsupported length=",dat_len,"\n")
    }
    ret_tib <- ret_tib %>% 
      dplyr::mutate(source=source,tissue=tissue) %>%
      dplyr::group_by(class) %>%
      dplyr::mutate(rank=dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(chr,beg,end) %>%
      dplyr::select(dplyr::all_of(sel_col))
    
    if (!is.null(out)) {
      out_file <- file.path(out, source, base::basename(file))
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing; out={out_file}{RET}"))

      safe_write(x=ret_tib,type="tsv",file=out_file,done=TRUE,funcTag=funcTag, 
                 verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }   
    
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    if (grs) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
      ret_tib = 
        GenomicRanges::GRanges(
          seqnames=Rle(ret_tib$chr), 
          strand=Rle(ret_tib$srd),
          
          name=ret_tib$name,
          name2=ret_tib$name2,
          class=ret_tib$class,
          source=ret_tib$source,
          tissue=ret_tib$tissue,
          rank=ret_tib$rank,
          evidence=ret_tib$evidence,
          
          IRanges(start=ret_tib$beg, 
                  end=ret_tib$end, 
                  names=ret_tib$unq_key)
        )
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; ",
                   "Return Count={ret_cnt}; ",
                   "elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ",
                   "----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Standard NCBI Annnotation Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_ncbi_gene = function(file,grs=FALSE,
                          source="NCBI",tissue="All",out=NULL,
                          sel_col=c("chr","beg","end","srd",
                                    "name","name2","class",
                                    "source","tissue","rank",
                                    "unq_key","evidence"),
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='load_ncbi_gene') {
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
    colnames(dat_tib)[1] <- stringr::str_remove(colnames(dat_tib)[1], '^#')
    dat_tib <- dat_tib %>% dplyr::distinct(name, .keep_all=TRUE) # %>% dplyr::mutate(transcript=name, name=name2)
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Raw Data={RET}"))
    if (verbose>=vt+4) print(dat_tib)
    
    dat_tib <- dat_tib %>% 
      dplyr::mutate(
        gene_tss=dplyr::case_when(
          strand=='+' ~ txStart,
          strand=='-' ~ txEnd,
          TRUE ~ NA_real_
        ),
        
        tss_200_beg=dplyr::case_when(
          strand=='+' ~ gene_tss-200,
          strand=='-' ~ gene_tss,
          TRUE ~ NA_real_
        ),
        tss_200_end=dplyr::case_when(
          strand=='+' ~ gene_tss,
          strand=='-' ~ gene_tss+200,
          TRUE ~ NA_real_
        ),
        
        tss_1500_beg=dplyr::case_when(
          strand=='+' ~ gene_tss-1500,
          strand=='-' ~ gene_tss+200,
          TRUE ~ NA_real_
        ),
        tss_1500_end=dplyr::case_when(
          strand=='+' ~ gene_tss-200,
          strand=='-' ~ gene_tss+1500,
          TRUE ~ NA_real_
        ),
      )
    
    ret_tib <- 
      dplyr::bind_rows(
        tibble::tibble(chr=dat_tib$chrom,
                       beg=dat_tib$txStart,
                       end=dat_tib$txEnd,
                       srd=dat_tib$strand,
                       name2=dat_tib$name2,
                       name=dat_tib$name,
                       class="GeneBody",
                       source=source,
                       tissue=tissue
        ),
        
        tibble::tibble(chr=dat_tib$chrom,
                       beg=dat_tib$tss_200_beg,
                       end=dat_tib$tss_200_end,
                       srd=dat_tib$strand,
                       name2=dat_tib$name2,
                       name=dat_tib$name,
                       class="TSS200",
                       source=source,
                       tissue=tissue
        ),
        
        tibble::tibble(chr=dat_tib$chrom,
                       beg=dat_tib$tss_1500_beg,
                       end=dat_tib$tss_1500_end,
                       srd=dat_tib$strand,
                       name2=dat_tib$name2,
                       name=dat_tib$name,
                       class="TSS1500",
                       source=source,
                       tissue=tissue
        )
      ) %>% 
      dplyr::group_by(class) %>%
      dplyr::mutate(rank=dplyr::row_number(),
                    unq_key=paste(class,rank, sep="_"),
                    evidence=as.integer(0)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(chr,beg,end) %>%
      dplyr::select(dplyr::all_of(sel_col))
    
    if (!is.null(out)) {
      out_file <- file.path(out, source, base::basename(file))
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing; out={out_file}{RET}"))
      
      safe_write(x=ret_tib,type="tsv",file=out_file,done=TRUE,funcTag=funcTag, 
                 verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    if (grs) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
      ret_tib = 
        GenomicRanges::GRanges(
          seqnames=Rle(ret_tib$chr), 
          strand=Rle(ret_tib$srd),
          
          # name=ret_tib$tran,
          # name2=ret_tib$gene,
          name=ret_tib$name,
          name2=ret_tib$name2,
          class=ret_tib$class,
          source=ret_tib$source,
          tissue=ret_tib$tissue,
          rank=ret_tib$rank,
          
          IRanges(start=ret_tib$beg, 
                  end=ret_tib$end, 
                  names=ret_tib$unq_key)
        )
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Standard UCSC Annotation Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_ucsc_gene = function(file,grs=FALSE,
                          source="UCSC",tissue="All",out=NULL,
                          sel_col=c("chr","beg","end","srd",
                                    "name","name2","class",
                                    "source","tissue","rank",
                                    "unq_key","evidence"),
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='load_ucsc_gene') {
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
    colnames(dat_tib)[1] <- stringr::str_remove(colnames(dat_tib)[1], '^#')
    dat_tib <- dat_tib %>% dplyr::distinct(name, .keep_all=TRUE)
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Raw Data={RET}"))
    if (verbose>=vt+4) print(dat_tib)
    
    dat_tib <- dat_tib %>% 
      dplyr::mutate(
        gene_tss=dplyr::case_when(
          strand=='+' ~ txStart,
          strand=='-' ~ txEnd,
          TRUE ~ NA_real_
        ),
        
        tss_200_beg=dplyr::case_when(
          strand=='+' ~ gene_tss-200,
          strand=='-' ~ gene_tss,
          TRUE ~ NA_real_),
        tss_200_end=dplyr::case_when(
          strand=='+' ~ gene_tss,
          strand=='-' ~ gene_tss+200,
          TRUE ~ NA_real_),
        
        tss_1500_beg=dplyr::case_when(
          strand=='+' ~ gene_tss-1500,
          strand=='-' ~ gene_tss+200,
          TRUE ~ NA_real_),
        tss_1500_end=dplyr::case_when(
          strand=='+' ~ gene_tss-200,
          strand=='-' ~ gene_tss+1500,
          TRUE ~ NA_real_)
      )
    
    ret_tib <- 
      dplyr::bind_rows(
        tibble::tibble(chr=dat_tib$chrom,
                       beg=dat_tib$txStart,
                       end=dat_tib$txEnd,
                       srd=dat_tib$strand,
                       name2=dat_tib$proteinID,
                       name=dat_tib$name,
                       class="GeneBody",
                       source=source,
                       tissue=tissue
        ),
        
        tibble::tibble(chr=dat_tib$chrom,
                       beg=dat_tib$tss_200_beg,
                       end=dat_tib$tss_200_end,
                       srd=dat_tib$strand,
                       name2=dat_tib$proteinID,
                       name=dat_tib$name,
                       class="TSS200",
                       source=source,
                       tissue=tissue
        ),
        
        tibble::tibble(chr=dat_tib$chrom,
                       beg=dat_tib$tss_1500_beg,
                       end=dat_tib$tss_1500_end,
                       srd=dat_tib$strand,
                       name2=dat_tib$proteinID,
                       name=dat_tib$name,
                       class="TSS1500",
                       source=source,
                       tissue=tissue
        )
      ) %>% 
      dplyr::group_by(class) %>%
      dplyr::mutate(rank=dplyr::row_number(),
                    unq_key=paste(class,rank, sep="_"),
                    evidence=as.integer(0)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(chr,beg,end) %>%
      dplyr::select(dplyr::all_of(sel_col))
    
    if (!is.null(out)) {
      out_file <- file.path(out, source, base::basename(file))
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing; out={out_file}{RET}"))
      
      safe_write(x=ret_tib,type="tsv",file=out_file,done=TRUE,funcTag=funcTag, 
                 verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)

    if (grs) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
      ret_tib = 
        GenomicRanges::GRanges(
          seqnames=Rle(ret_tib$chr), 
          strand=Rle(ret_tib$srd),
          
          name=ret_tib$name,
          name2=ret_tib$name2,
          class=ret_tib$class,
          source=ret_tib$source,
          tissue=ret_tib$tissue,
          rank=ret_tib$rank,
          
          IRanges(start=ret_tib$beg, 
                  end=ret_tib$end, 
                  names=ret_tib$unq_key)
        )
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

load_ucsc_cpgs = function(file, grs=FALSE,
                          source="UCSC",tissue="All",out=NULL,
                          sel_col=c("chr","beg","end","srd",
                                    "name","name2","class",
                                    "source","tissue","rank",
                                    "unq_key","evidence"),
                          verbose=0,vt=3,tc=1,tt=NULL,
                          funcTag='load_ucsc_cpgs') {
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) )) %>% 
      dplyr::mutate(name2=paste0(chrom,'-',chromStart,'-',chromEnd))
    
    ret_tib <- dplyr::bind_rows(
      
      tibble::tibble(chr=dat_tib$chrom,
                     beg=dat_tib$chromStart-4000,
                     end=dat_tib$chromStart-2000,
                     srd="+",
                     name2=dat_tib$name2,
                     class="NShelf",
                     source=source,
                     tissue=tissue
      ),
      
      tibble::tibble(chr=dat_tib$chrom,
                     beg=dat_tib$chromStart-2000,
                     end=dat_tib$chromStart,
                     srd="+",
                     name2=dat_tib$name2,
                     class="NShore",
                     source=source,
                     tissue=tissue
      ),
      
      tibble::tibble(chr=dat_tib$chrom,
                     beg=dat_tib$chromStart,
                     end=dat_tib$chromEnd,
                     srd="+",
                     name2=dat_tib$name2,
                     class="Island",
                     source=source,
                     tissue=tissue
      ),
      
      tibble::tibble(chr=dat_tib$chrom,
                     beg=dat_tib$chromEnd,
                     end=dat_tib$chromEnd+2000,
                     srd="+",
                     name2=dat_tib$name2,
                     class="SShore",
                     source=source,
                     tissue=tissue
      ),
      
      tibble::tibble(chr=dat_tib$chrom,
                     beg=dat_tib$chromEnd+2000,
                     end=dat_tib$chromEnd+4000,
                     srd="+",
                     name2=dat_tib$name2,
                     class="SShelf",
                     source=source,
                     tissue=tissue
      )
    ) %>% 
      dplyr::group_by(class) %>%
      dplyr::mutate(rank=dplyr::row_number(),
                    unq_key=paste(class,rank, sep="_"),
                    name=paste(chr,beg,end, sep="-"),
                    evidence=as.integer(0)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(chr,beg,end) %>%
      dplyr::select(dplyr::all_of(sel_col))
    
    if (!is.null(out)) {
      out_file <- file.path(out, source, base::basename(file))
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing; out={out_file}{RET}"))
      
      safe_write(x=ret_tib,type="tsv",file=out_file,done=TRUE,funcTag=funcTag, 
                 verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    ret_key <- glue::glue("ret-fin({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
    
    if (grs) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
      ret_tib <-
        GenomicRanges::GRanges(
          seqnames=Rle(ret_tib$chr), 
          strand=Rle(ret_tib$srd),
          
          name=ret_tib$name,
          class=ret_tib$class,
          name2=NA_character_,
          source=ret_tib$source,
          tissue=ret_tib$tissue,
          rank=ret_tib$rank,
          
          IRanges(start=ret_tib$beg, 
                  end=ret_tib$end, 
                  names=ret_tib$unq_key)
        )
    }
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Standard UCSC Loading Methods:: OLD CODE
#                           SHOULD BE DELETED!!
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadNcbiGeneGR = function(file,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadNcbiGeneGR'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
    colnames(dat_tib)[1] <- stringr::str_remove(colnames(dat_tib)[1], '^#')
    dat_tib <- dat_tib %>% dplyr::distinct(name, .keep_all=TRUE) # %>% dplyr::mutate(transcript=name, name=name2)
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Raw Data={RET}"))
    if (verbose>=vt+4) print(dat_tib)
    
    dat_tib <- dat_tib %>% 
      dplyr::mutate(
        gene_tss=dplyr::case_when(
          strand=='+' ~ txStart,
          strand=='-' ~ txEnd,
          TRUE ~ NA_real_
        ),
        
        tss_200_beg=dplyr::case_when(
          strand=='+' ~ gene_tss-200,
          strand=='-' ~ gene_tss,
          TRUE ~ NA_real_
        ),
        tss_200_end=dplyr::case_when(
          strand=='+' ~ gene_tss,
          strand=='-' ~ gene_tss+200,
          TRUE ~ NA_real_
        ),
        
        tss_1500_beg=dplyr::case_when(
          strand=='+' ~ gene_tss-1500,
          strand=='-' ~ gene_tss+200,
          TRUE ~ NA_real_
        ),
        tss_1500_end=dplyr::case_when(
          strand=='+' ~ gene_tss-200,
          strand=='-' ~ gene_tss+1500,
          TRUE ~ NA_real_
        ),
      )
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
    ret_grs$tss_1500 <- 
      GenomicRanges::GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand), # seqinfo=dat_tib$name2,
                             IRanges(start=dat_tib$tss_1500_beg, end=dat_tib$tss_1500_end, names=dat_tib$name) )
    ret_grs$tss_200 <- 
      GenomicRanges::GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand), # seqinfo=dat_tib$name2,
                             IRanges(start=dat_tib$tss_200_beg, end=dat_tib$tss_200_end, names=dat_tib$name) )
    ret_grs$tss_body <- 
      GenomicRanges::GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand), # seqinfo=dat_tib$name2,
                             IRanges(start=dat_tib$txStart, end=dat_tib$txEnd, names=dat_tib$name) )
    
    ret_cnt <- ret_grs %>% names %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_grs
}

loadUcscGeneGR = function(file,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadUcscGeneGR'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
    colnames(dat_tib)[1] <- stringr::str_remove(colnames(dat_tib)[1], '^#')
    dat_tib <- dat_tib %>% dplyr::distinct(name, .keep_all=TRUE)
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Raw Data={RET}"))
    if (verbose>=vt+4) print(dat_tib)
    
    dat_tib <- dat_tib %>% 
      dplyr::mutate(
        gene_tss=dplyr::case_when(
          strand=='+' ~ txStart,
          strand=='-' ~ txEnd,
          TRUE ~ NA_real_
        ),
        
        tss_200_beg=dplyr::case_when(
          strand=='+' ~ gene_tss-200,
          strand=='-' ~ gene_tss,
          TRUE ~ NA_real_
        ),
        tss_200_end=dplyr::case_when(
          strand=='+' ~ gene_tss,
          strand=='-' ~ gene_tss+200,
          TRUE ~ NA_real_
        ),
        
        tss_1500_beg=dplyr::case_when(
          strand=='+' ~ gene_tss-1500,
          strand=='-' ~ gene_tss+200,
          TRUE ~ NA_real_
        ),
        tss_1500_end=dplyr::case_when(
          strand=='+' ~ gene_tss-200,
          strand=='-' ~ gene_tss+1500,
          TRUE ~ NA_real_
        ),
      )
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
    ret_grs$tss_1500 <- 
      GenomicRanges::GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand),
                             IRanges(start=dat_tib$tss_1500_beg, end=dat_tib$tss_1500_end, names=dat_tib$name) )
    ret_grs$tss_200 <- 
      GenomicRanges::GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand),
                             IRanges(start=dat_tib$tss_200_beg, end=dat_tib$tss_200_end, names=dat_tib$name) )
    ret_grs$tss_body <- 
      GenomicRanges::GRanges(seqnames=Rle(dat_tib$chrom), strand=Rle(dat_tib$strand),
                             IRanges(start=dat_tib$txStart, end=dat_tib$txEnd, names=dat_tib$name) )
    
    ret_cnt <- ret_grs %>% names %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_grs
}

loadUcscCpgsGR = function(file,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadUcscCpgsGR'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) )) %>% 
      dplyr::mutate(name=paste0(chrom,'-',chromStart,'-',chromEnd))
    
    # Build GRange Data Structures:: CpG Islands
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
    ret_grs$NShelf <- GenomicRanges::GRanges(seqnames=Rle(dat_tib$chrom), 
                                             IRanges(start=dat_tib$chromStart-4000, end=dat_tib$chromStart-2000, names=dat_tib$name) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Built North Shelves.{RET}"))
    
    ret_grs$NShore <- GenomicRanges::GRanges(seqnames=Rle(dat_tib$chrom), 
                                             IRanges(start=dat_tib$chromStart-2000, end=dat_tib$chromStart, names=dat_tib$name) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Built North Shores.{RET}"))
    
    ret_grs$Island <- GenomicRanges::GRanges(seqnames=Rle(dat_tib$chrom), 
                                             IRanges(start=dat_tib$chromStart, end=dat_tib$chromEnd, names=dat_tib$name) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Built Islands.{RET}"))
    
    ret_grs$SShore <- GenomicRanges::GRanges(seqnames=Rle(dat_tib$chrom), 
                                             IRanges(start=dat_tib$chromEnd, end=dat_tib$chromEnd+2000, names=dat_tib$name) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Built South Shores.{RET}"))
    
    ret_grs$SShelf <- GenomicRanges::GRanges(seqnames=Rle(dat_tib$chrom), 
                                             IRanges(start=dat_tib$chromEnd+2000, end=dat_tib$chromEnd+4000, names=dat_tib$name) )
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Built South Shelves.{RET}"))
    
    ret_cnt <- ret_grs %>% names %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_grs
}



# End of file
