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

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Manifest Coordinate Intersection Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

impBuildSummary = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'impBuildSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- tib %>% 
      dplyr::distinct(U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                      chrom,chromStart,chromEnd,CGN_IMP,TB_IMP,CO_IMP, 
                      Probe_Type,Manifest,Unq_Cnt, .keep_all=TRUE) %>%
      dplyr::group_by(Manifest,Probe_Type,chrFlag,difFlag,cgnFlag,srcFlag,
                      Infinium_Design) %>%
      #                srcFlagU,srcFlagM) %>% 
      dplyr::summarise(Count=n(), .groups="drop")

    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt) {
      cat(glue::glue("[Summary]: sum_tib({ret_cnt})={RET}"))
      ret_tib %>% print(n=base::nrow(ret_tib))
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# tsv=out_pos1_tsv
# datDir="/Users/bretbarnes/Documents/data/improbe/designOutput_21092020"
# inpDir=cgn_pos1_dir
# runType=par$local_runType
# ver=opt$version
# build="GRCh37"
intersectCgnMap_COVIC = function(tsv, man, src=NULL,
                                 build, datDir, inpDir, runType, ver,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'intersectCgnMap_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # TBD:: Make dir part of par$posDir
    #
    cgn_pos_fns <- paste(build,"21092020_improbe-designOutput.cgn-map-seq.cgn-sorted.tsv.gz", sep='-' )
    cgn_pos_pre <- paste(build,runType,ver, sep='-')
    imp_pos_tsv <- file.path(datDir, cgn_pos_fns)
    out_pos_tsv <- file.path(inpDir, paste(cgn_pos_pre,'seq48U_to_cgn.int.improbe-designOutput.cgn-map-seq.cgn-sorted.tsv.gz', sep='.'))
    
    if (!file.exists(out_pos_tsv)) {
      # Build better database::: Combine the two commands below::
      # gzip -dc /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.tsv.gz | cut -f 1,4,5,21-23,45,48 > /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.tsv
      # cat /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.tsv | perl -pe 's/TOP/T/; s/BOT/B/;' | gzip -c - > /Users/bretbarnes/Documents/data/improbe/designOutput_21092020/GRCh38-21092020_improbe-designOutput.cgn-map-seq.tsv.gz
      
      int_pos_cmd <- glue::glue("gzip -dc {imp_pos_tsv} | join -t $'\t' -11 -22 - {tsv} | gzip -c - > {out_pos_tsv}")
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]: Running; cmd={int_pos_cmd}...{RET}"))
      system(int_pos_cmd)
    }
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: Loading Coordinates; out_pos_tsv={out_pos_tsv}...{RET}"))
    imp_pos_col <- cols(
      CGN_IMP = col_character(),
      chrom = col_character(),
      chromStart = col_integer(),
      FR_IMP = col_character(),
      TB_IMP = col_character(),
      CO_IMP = col_character(),
      NB_IMP = col_character(),
      Seq50U = col_character(),
      
      Seq48U = col_character(),
      TB_DB2 = col_character(),
      CO_DB2 = col_character(),
      U_DB2 = col_integer(),
      M_DB2 = col_integer()
    )
    imp_pos_tib <- 
      readr::read_tsv(out_pos_tsv, 
                      col_names=names(imp_pos_col$cols), 
                      col_types=imp_pos_col) %>% 
      # dplyr::distinct(CGN_IMP,chrom,chromStart,FR_IMP,TB_IMP,CO_IMP,NB_IMP,Seq50U) %>%
      dplyr::mutate(chrom=stringr::str_remove(chrom,'^chr'),
                    chrom=paste0('chr',chrom),
                    chromEnd=chromStart+1)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #             4.1 Match All Probes To CG# Database Coordinates::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- man %>%
      dplyr::inner_join(imp_pos_tib,
                        by=c("PRB1_U_MAT"="Seq50U"), 
                        suffix=c("_PQC", "_IMP")) %>% 
      dplyr::arrange(chrom,chromStart) %>%
      dplyr::select(chrom,chromStart,chromEnd, 
                    CGN_IMP,FR_IMP,TB_IMP,CO_IMP,NB_IMP,
                    PRB1_U_MAT,
                    dplyr::everything()) %>% 
      dplyr::add_count(U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                       chrom,chromStart,chromEnd,CGN_IMP,TB_IMP,CO_IMP, 
                       name="Unq_Cnt")
    
    if (build=="GRCh37") {
      ret_tib <- ret_tib %>%
        dplyr::mutate(
          begDif=chromStart-chromStart_GRCh37,
          endDif=chromEnd-chromEnd_GRCh37,
          absDif=base::abs(begDif),
          cgnFlag=dplyr::case_when(
            is.na(CGN_IMP) ~ -2,
            is.na(Seq_ID) ~ -1,
            CGN_IMP==Seq_ID ~ 0,
            TRUE ~ 1
          ),
          chrFlag=dplyr::case_when(
            is.na(chrom) ~ -2,
            is.na(chrom_GRCh37) ~ -1,
            chrom==chrom_GRCh37 ~ 0,
            TRUE ~ 1
          ),
          difFlag=dplyr::case_when(
            is.na(absDif) ~ -1,
            absDif==0 ~ 0,
            TRUE ~ 1
          )
        )
    } else if (build=="GRCh38") {
      ret_tib <- ret_tib %>%
        dplyr::mutate(
          begDif=chromStart-chromStart_GRCh38,
          endDif=chromEnd-chromEnd_GRCh38,
          absDif=base::abs(begDif),
          cgnFlag=dplyr::case_when(
            is.na(CGN_IMP) ~ -2,
            is.na(Seq_ID) ~ -1,
            CGN_IMP==Seq_ID ~ 0,
            TRUE ~ 1
          ),
          chrFlag=dplyr::case_when(
            is.na(chrom) ~ -2,
            is.na(chrom_GRCh38) ~ -1,
            chrom==chrom_GRCh38 ~ 0,
            TRUE ~ 1
          ),
          difFlag=dplyr::case_when(
            is.na(absDif) ~ -1,
            absDif==0 ~ 0,
            TRUE ~ 1
          )
        )
    } else {
      stop(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: Unsupported build type={build}!{RET}{RET}"))
      return(NULL)
    }
    
    ret_tib <- ret_tib %>%
      dplyr::mutate(
        begDif=as.integer(begDif),
        endDif=as.integer(endDif),
        absDif=as.integer(absDif),
        cgnFlag=as.integer(cgnFlag),
        chrFlag=as.integer(chrFlag),
        difFlag=as.integer(difFlag)
      )
    
    if (!is.null(src)) {
      ret_tib <- ret_tib %>%
        dplyr::mutate(
          srcFlagU=dplyr::case_when(
            U %in% src$U ~ 0,
            TRUE ~ 1),
          srcFlagM=dplyr::case_when(
            is.na(M) ~ 0,
            M %in% src$M ~ 0,
            TRUE ~ 1),
          srcFlagU=as.integer(srcFlagU),
          srcFlagM=as.integer(srcFlagM),
          srcFlag=dplyr::case_when(
            srcFlagU==0 & srcFlagM==0 ~ 0,
            TRUE ~ 1
          ),
          srcFlag=as.integer(srcFlag)
        ) %>%
        dplyr::arrange(AlleleA_Probe_Sequence, chrFlag, difFlag, cgnFlag, absDif, 
                       srcFlagM, srcFlagU, srcFlagM)
    } else {
      ret_tib <- ret_tib %>%
        dplyr::arrange(AlleleA_Probe_Sequence, chrFlag, difFlag, cgnFlag, absDif)
    }

    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt) {
      sum1_tib <- ret_tib %>% 
        dplyr::group_by(Manifest,Probe_Type) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      cat(glue::glue("[{funcTag}]: sum1_tib({ret_cnt})={RET}"))
      sum1_tib %>% print(n=base::nrow(sum1_tib))

      # For hg19 & hg38 this has nrows = 0
      mis_cnt <- ret_tib %>% dplyr::filter(Unq_Cnt != 1) %>% base::nrow()
      cat(glue::glue("[{funcTag}]: mis_cnt({ret_cnt})={mis_cnt}.{RET}"))

      sum2_tib <- ret_tib %>% 
        dplyr::distinct(U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                        chrom,chromStart,chromEnd,CGN_IMP,TB_IMP,CO_IMP, 
                        Probe_Type,Manifest,Unq_Cnt, .keep_all=TRUE) %>%
        dplyr::group_by(Manifest,Probe_Type,chrFlag,difFlag,cgnFlag,Infinium_Design) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      cat(glue::glue("[{funcTag}]: sum2_tib({ret_cnt})={RET}"))
      sum2_tib %>% print(n=base::nrow(sum2_tib))
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Standard Manifest Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

manifestToAddress_COVIC = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'manifestToAddress_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- dplyr::bind_rows(
      
      tib %>% dplyr::filter(!is.na(U)) %>%
        dplyr::filter(is.na(M)) %>% 
        dplyr::select(U, Probe_ID,Probe_Type,Manifest) %>%
        dplyr::rename(Address=U) %>%
        dplyr::mutate(Probe_Design="2"),
      
      tib %>% dplyr::filter(!is.na(U)) %>%
        dplyr::filter(!is.na(M)) %>% 
        dplyr::select(U, Probe_ID,Probe_Type,Manifest) %>%
        dplyr::rename(Address=U) %>%
        dplyr::mutate(Probe_Design="U"),
      
      tib %>% dplyr::filter(!is.na(U)) %>%
        dplyr::filter(!is.na(M)) %>% 
        dplyr::select(M, Probe_ID,Probe_Type,Manifest) %>%
        dplyr::rename(Address=M) %>%
        dplyr::mutate(Probe_Design="M")
    ) %>%
      dplyr::mutate(Address=as.integer(Address)) %>%
      dplyr::arrange(Address) %>%
      dplyr::select(Address,Probe_ID,Probe_Type,Probe_Design,Manifest)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Standard Manifest Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadCoreManifest_COVIC = function(datDir, manDir,
                                  verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadCoreManifest_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                1.0 Load Sesame hg37/hg38 Manifest:: S4
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ses_hg37_tib <- 
      loadSesameManifest_COVIC(build="GRCh37", tag=TRUE, 
                               verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    ses_hg38_tib <-
      loadSesameManifest_COVIC(build="GRCh38", tag=TRUE,
                               verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)

    ret_tib <- dplyr::inner_join(
      ses_hg37_tib,ses_hg38_tib, 
      by=c("Probe_ID","Probe_Type","U","M",
           "AlleleA_Probe_Sequence","AlleleB_Probe_Sequence")) %>%
      dplyr::mutate(Manifest="B4") %>%
      dplyr::select(Probe_ID,Probe_Type,U,M,
                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Manifest,
                    dplyr::everything())
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Sesame(hg37/hg38) manifest={RET}"))
      print(ret_tib)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       1.1 Load EPIC Manifest:: B2
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    epic_csv <- file.path(manDir, 'MethylationEPIC_v-1-0_B2.csv.gz')
    epic_tib <- 
      loadManifestGenomeStudio(file = epic_csv, addSource = TRUE, 
                               normalize = TRUE, retType = "man", 
                               verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>% 
      dplyr::rename(Probe_ID=IlmnID, 
                    chrom_GRCh37=Chromosome,
                    chromStart_GRCh37=Coordinate,
                    FR_GRCh37=Strand_FR,
                    nextBase_GRCh37=Next_Base,
                    AlleleA_Probe_Sequence=AlleleA_ProbeSeq,
                    AlleleB_Probe_Sequence=AlleleB_ProbeSeq,
                    U=AddressA_ID, M=AddressB_ID) %>%
      dplyr::mutate(U=as.integer(U),
                    M=as.integer(M),
                    chrom_GRCh37=paste0("chr",chrom_GRCh37),
                    chromEnd_GRCh37=chromStart_GRCh37+1,
                    Manifest="B2") %>%
      dplyr::select(Probe_ID,Probe_Type,U,M,
                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                    Manifest,chrom_GRCh37,chromStart_GRCh37,chromEnd_GRCh37,
                    FR_GRCh37,nextBase_GRCh37)
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} EPIC(B2) manifest={RET}"))
      print(epic_tib)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       1.2 Bind Manifests:: C0/B2
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ret_tib <- ret_tib %>% 
      dplyr::bind_rows(dplyr::anti_join(epic_tib, ret_tib, by="U") ) %>% 
      addSeq48U(field="AlleleA_Probe_Sequence") %>%
      dplyr::mutate(
        Seq_48U=dplyr::case_when( is.na(M) ~ Seq_48U_2,TRUE ~ Seq_48U_1 )
      ) %>%
      dplyr::select(-Seq_48U_1,-Seq_48U_2) %>%
      dplyr::select(Probe_ID,Probe_Type,U,M,
                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Seq_48U,
                    dplyr::everything())
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       1.X Load COVIC Manifest:: C0
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Load Original COVIC (EPIC_C0) manifest and compare to idat_join_tib
    # covic_c0_csv <- file.path(datDir, 'manifest/base/EPIC-C0.manifest.sesame-base.cpg-sorted.csv.gz')
    # covic_c0_tib <- readr::read_csv(covic_c0_csv) %>% 
    #   dplyr::mutate(M=as.integer(M), U=as.integer(U))
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Core Manifests({ret_cnt})={RET}{RET}"))
      print(ret_tib)
    }
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Core Summary({ret_cnt})={RET}{RET}"))
      ret_tib %>% 
        dplyr::group_by(Probe_Type,Manifest) %>% 
        dplyr::summarise(Count=n(), .groups="drop") %>% print()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Sesame Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadSesameManifest_COVIC = function(build, tag=FALSE, del="_",
                                    verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadSesameManifest_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    genomeBuild_UCSC <- NA_character_
    if (build=='GRCm38') genomeBuild_UCSC <- 'mm10'
    if (build=='GRCh38') genomeBuild_UCSC <- 'hg38'
    if (build=='GRCh37') genomeBuild_UCSC <- 'hg19'
    if (build=='GRCh36') genomeBuild_UCSC <- 'hg18'
    
    if (is.na(genomeBuild_UCSC)) {
      stop(glue::glue("{RET}[{funcTag}]: genomeBuild_UCSC is null!{RET}{RET}"))
      return(ret_tib)
    }
    
    ses_fns <- paste('EPIC',genomeBuild_UCSC,'manifest', sep='.')
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading Sesame Manifest={ses_fns}...{RET}"))
    }
    
    ret_tib <- sesameData::sesameDataGet(ses_fns) %>%
      as.data.frame() %>% 
      tibble::rownames_to_column(var="Probe_ID") %>%
      tibble::as_tibble() %>%
      dplyr::rename(Probe_Type=probeType, 
                    U=address_A, M=address_B,
                    AlleleA_Probe_Sequence=ProbeSeq_A,
                    AlleleB_Probe_Sequence=ProbeSeq_B,
                    chrom=seqnames, chromStart=start, chromEnd=end, 
                    Design_Type=designType) %>% 
      dplyr::mutate(U=as.integer(U), M=as.integer(M),
                    chrom=as.character(chrom),
                    chromStart=as.integer(chromStart),
                    chromEnd=as.integer(chromEnd),
                    FR=dplyr::case_when(
                      strand=='+' ~ "F",
                      strand=='-' ~ "R",
                      TRUE ~ NA_character_
                    ),
                    AlleleB_Probe_Sequence=dplyr::case_when(
                      stringr::str_length(AlleleB_Probe_Sequence)==0 ~ NA_character_, 
                      TRUE ~ AlleleB_Probe_Sequence)
      ) %>%
      dplyr::distinct(U,M, .keep_all=TRUE) %>% 
      dplyr::select(Probe_ID,Probe_Type, 
                    U,M,AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, 
                    chrom,chromStart,chromEnd,
                    FR,nextBase, 
                    gene,gene_HGNC,MASK_mapping,MASK_general)
    
    if (tag) {
      ret_tib <- ret_tib %>% purrr::set_names(paste(names(.),build, sep=del))
      
      colnames(ret_tib)[1] = "Probe_ID"
      colnames(ret_tib)[2] = "Probe_Type"
      colnames(ret_tib)[3] = "U"
      colnames(ret_tib)[4] = "M"
      colnames(ret_tib)[5] = "AlleleA_Probe_Sequence"
      colnames(ret_tib)[6] = "AlleleB_Probe_Sequence"
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} Core Manifests({ret_cnt})={RET}{RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           COVIC Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadAqpWorkflow_COVIC = function(ords, mats, aqps, man=NULL,
                                 verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadAqpWorkflow_COVIC'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ords_vec <- ords %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
    mats_vec <- mats %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
    aqps_vec <- aqps %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()
    
    # Check ORD Data::
    covic_ord_tib <- suppressMessages(suppressWarnings( 
      readr::read_csv(ords_vec[1], skip = 8, guess_max = 50000) )) %>%
      dplyr::rename(Order_ID=Assay_Design_Id) %>%
      dplyr::mutate(AlleleA_Probe_Sequence=stringr::str_to_upper(AlleleA_Probe_Sequence),
                    AlleleB_Probe_Sequence=stringr::str_to_upper(AlleleB_Probe_Sequence),
                    Design_TypeA=dplyr::case_when(!is.na(AlleleB_Probe_Sequence) ~ 'U', TRUE ~ '2' ), 
                    Design_TypeB=dplyr::case_when(!is.na(AlleleB_Probe_Sequence) ~ 'M', TRUE ~ '2' ) ) %>% 
      dplyr::select(-AlleleA_Probe_Id,-AlleleB_Probe_Id,-Normalization_Bin) %>%
      dplyr::arrange(Order_ID) %>%
      dplyr::distinct(AlleleA_Probe_Sequence,AlleleB_Probe_Sequence, .keep_all=TRUE) %>%
      dplyr::mutate(Order_ID=stringr::str_remove(Order_ID, '_[0-9]+$') %>%
                      stringr::str_replace('_','-') %>% 
                      stringr::str_replace('II', '2') %>% 
                      stringr::str_replace('I', '1') %>%
                      stringr::str_remove_all('_') %>%
                      stringr::str_replace('-','_') )
    
    # Build flatten ORD table::
    #   - Non Unique = 10,537; NEW => 10,462
    covic_ord_tab <- dplyr::bind_rows(
      covic_ord_tib %>% dplyr::select(-AlleleB_Probe_Sequence,-Design_TypeB) %>% 
        dplyr::rename(Probe_Seq=AlleleA_Probe_Sequence,Design_Type=Design_TypeA),
      covic_ord_tib %>% dplyr::select(-AlleleA_Probe_Sequence,-Design_TypeA) %>% 
        dplyr::rename(Probe_Seq=AlleleB_Probe_Sequence,Design_Type=Design_TypeB),
    ) %>% dplyr::filter(!is.na(Probe_Seq) & !is.na(Design_Type)) %>% 
      dplyr::arrange(Order_ID) %>%
      dplyr::distinct(Probe_Seq,Design_Type, .keep_all=TRUE)    
    
    # Resolve Duplicates:: INCOMPLETE
    covic_ord_tab %>% 
      dplyr::add_count(Probe_Seq, name="Seq_Count") %>% 
      dplyr::filter(Seq_Count!=1) %>% 
      dplyr::arrange(-Seq_Count,Probe_Seq)
    
    # Check MAT Data::
    covic_mat_tib <- suppressMessages(suppressWarnings( 
      readr::read_tsv(mats_vec[1]) )) %>%
      dplyr::rename(Match_ID=probe_id,
                    Address=address_name,
                    Probe_Seq=sequence) %>%
      dplyr::mutate(Address=stringr::str_remove(Address, '^1') %>% 
                      stringr::str_remove('^0+') %>% as.integer(),
                    Probe_Seq=stringr::str_to_upper(Probe_Seq)) %>%
      dplyr::select(Match_ID,Probe_Seq,Address) %>%
      dplyr::arrange(Address) %>%
      dplyr::distinct(Address, .keep_all=TRUE) %>%
      dplyr::mutate(Match_ID=stringr::str_replace(Match_ID, '_','-') %>%
                      stringr::str_replace('II', '2') %>% 
                      stringr::str_replace('I', '1') %>%
                      stringr::str_remove_all('_') %>%
                      stringr::str_replace('-','_') )
    
    # Check PQC Data::
    #
    covic_aqp_tib <- loadPQC(file = aqps_vec[1], format = 'aqp', verbose = 4) %>% 
      dplyr::mutate(Address=stringr::str_remove(Address, '^1') %>% 
                      stringr::str_remove('^0+') %>% as.integer(),
                    Decode_Status=as.integer(Decode_Status)) %>%
      dplyr::select(Address,Decode_Status) %>%
      dplyr::arrange(Address) %>%
      dplyr::distinct(Address, .keep_all=TRUE)
    
    # Join MAT & AQP::
    covic_mat_aqp_raw_tab <- covic_mat_tib %>% 
      dplyr::inner_join(covic_aqp_tib, by="Address") %>% 
      dplyr::distinct(Probe_Seq,Address, .keep_all=TRUE)
    
    covic_mat_aqp_pas_tab <- covic_mat_aqp_raw_tab %>%
      dplyr::filter(!is.na(Decode_Status) & Decode_Status==0) %>%
      dplyr::distinct(Probe_Seq,Address, .keep_all=TRUE)
    
    # Join all data:: Decode Pass
    covic_ann_aqp_tab <- covic_ord_tab %>% 
      dplyr::inner_join(covic_mat_aqp_pas_tab, by="Probe_Seq") %>% 
      dplyr::select(-Decode_Status)
    
    # Rebuild Pairs:: covic_ord_tib U covic_ann_aqp_tab
    #
    ret_tib <- dplyr::bind_rows(
      covic_ord_tib %>% 
        dplyr::filter(is.na(AlleleB_Probe_Sequence)) %>%
        dplyr::inner_join(
          covic_ann_aqp_tab %>% 
            dplyr::rename(Order_ID_A=Order_ID, Match_ID_A=Match_ID, U=Address), 
          by=c("AlleleA_Probe_Sequence"="Probe_Seq",
               "Design_TypeA"="Design_Type")
        ) %>% dplyr::filter(!is.na(AlleleA_Probe_Sequence) & !is.na(U)),
      
      covic_ord_tib %>% 
        dplyr::filter(!is.na(AlleleB_Probe_Sequence)) %>%
        dplyr::inner_join(
          covic_ann_aqp_tab %>% 
            dplyr::rename(Order_ID_A=Order_ID, Match_ID_A=Match_ID, U=Address), 
          by=c("AlleleA_Probe_Sequence"="Probe_Seq",
               "Design_TypeA"="Design_Type")
        ) %>%
        dplyr::left_join(
          covic_ann_aqp_tab %>% 
            dplyr::rename(Order_ID_B=Order_ID, Match_ID_B=Match_ID, M=Address), 
          by=c("AlleleB_Probe_Sequence"="Probe_Seq",
               "Design_TypeB"="Design_Type")
        ) %>% dplyr::filter(!is.na(AlleleA_Probe_Sequence) & !is.na(U),
                            !is.na(AlleleB_Probe_Sequence) & !is.na(M))
    ) %>% dplyr::distinct(M,U, .keep_all=TRUE) %>% 
      dplyr::select(Order_ID,M,U,
                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,
                    Design_TypeA,Design_TypeB) %>%
      dplyr::rename(Probe_ID=Order_ID) %>%
      dplyr::mutate(Probe_Type=stringr::str_sub(Probe_ID, 1,2))  %>%
      addSeq48U(field="AlleleA_Probe_Sequence") %>%
      dplyr::mutate(
        Seq_48U=dplyr::case_when( is.na(M) ~ Seq_48U_2,TRUE ~ Seq_48U_1 )
      ) %>%
      dplyr::select(Probe_ID,Probe_Type,U,M,
                    AlleleA_Probe_Sequence,AlleleB_Probe_Sequence,Seq_48U)
    
    if (!is.null(man)) ret_tib <- ret_tib %>% dplyr::mutate(Manifest=man)

    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} AQP({ret_cnt})={RET}{RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        IDAT Samples Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# loadIdatAddress(prefix="/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R01C01")
loadIdatAddress = function(prefix,
                           verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'loadIdatAddress'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # idat_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R01C01"
    grn_dat <- illuminaio::readIDAT(file=paste(prefix,"Grn.idat.gz", sep="_"), what="all")
    ret_tib <- grn_dat$Quants %>% as.data.frame() %>% 
      tibble::rownames_to_column(var="Address") %>% 
      dplyr::mutate(Address=as.integer(Address)) %>%
      tibble::as_tibble() %>% 
      dplyr::filter(!is.na(Address)) %>% 
      dplyr::distinct(Address) %>% 
      dplyr::arrange(Address)
    
    if (FALSE) {
      idat1_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R01C01"
      idat1_tib <- prefixToIdat(prefix = idat1_prefix, verbose=opt$verbose)
      
      # idat2_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R02C01"
      # idat2_tib <- prefixToIdat(prefix = idat2_prefix, verbose=opt$verbose)
      # 
      # idat3_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R03C01"
      # idat3_tib <- prefixToIdat(prefix = idat3_prefix, verbose=opt$verbose)
      # 
      # idat4_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R04C01"
      # idat4_tib <- prefixToIdat(prefix = idat4_prefix, verbose=opt$verbose)
      # 
      # idat5_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R05C01"
      # idat5_tib <- prefixToIdat(prefix = idat5_prefix, verbose=opt$verbose)
      # 
      # idat6_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R06C01"
      # idat6_tib <- prefixToIdat(prefix = idat6_prefix, verbose=opt$verbose)
      # 
      # idat7_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R07C01"
      # idat7_tib <- prefixToIdat(prefix = idat7_prefix, verbose=opt$verbose)
      
      idat8_prefix <- "/Users/bretbarnes/Documents/data/idats/idats_COVIC-Set1-15052020/204500250025/204500250025_R08C01"
      idat8_tib <- prefixToIdat(prefix = idat8_prefix, verbose=opt$verbose)
      
      # Combine all tangos::
      idats_tib <- 
        dplyr::select(idat1_tib, Address) %>% dplyr::mutate(Src1=1) %>% dplyr::full_join(
          # dplyr::select(idat2_tib, Address) %>% dplyr::mutate(Src2=2), by="Address" ) %>% dplyr::full_join(
          #   dplyr::select(idat3_tib, Address) %>% dplyr::mutate(Src3=3), by="Address" ) %>% dplyr::full_join(
          #     dplyr::select(idat4_tib, Address) %>% dplyr::mutate(Src4=4), by="Address" ) %>% dplyr::full_join(
          #       dplyr::select(idat5_tib, Address) %>% dplyr::mutate(Src5=5), by="Address" ) %>% dplyr::full_join(
          #         dplyr::select(idat6_tib, Address) %>% dplyr::mutate(Src6=6), by="Address" ) %>% dplyr::full_join(
          #           dplyr::select(idat7_tib, Address) %>% dplyr::mutate(Src7=7), by="Address" ) %>% dplyr::full_join(
          dplyr::select(idat8_tib, Address) %>% dplyr::mutate(Src8=8), by="Address" ) 
      
      unq_idats_tib <- idats_tib %>% dplyr::filter(!is.na(Address)) %>% 
        dplyr::distinct(Address) %>% dplyr::arrange(Address)
      unq_cnt <- unq_idats_tib %>% base::nrow()
      idats_sum_tib <- idats_tib %>% dplyr::summarise(
        Tot_Cnt=n(),
        S1=sum(is.na(Src1)),
        # S2=sum(is.na(Src2)),
        # S3=sum(is.na(Src3)),
        # S4=sum(is.na(Src4)),
        # S5=sum(is.na(Src5)),
        # S6=sum(is.na(Src6)),
        # S7=sum(is.na(Src7)),
        S8=sum(is.na(Src8)),
      ) %>% dplyr::mutate(Unq_Cnt=unq_cnt, Dif_Cnt=Tot_Cnt-Unq_Cnt)
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} idat({ret_cnt})={RET}{RET}"))
      print(ret_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# End of file
