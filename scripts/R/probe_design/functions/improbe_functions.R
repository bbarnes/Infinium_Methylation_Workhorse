
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   improbe (Infinium Methylation) Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

suppressWarnings(suppressPackageStartupMessages(require("Biostrings")) )

COM <- ","
TAB <- "\t"
RET <- "\n"
BNG <- "|"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Docker improbe Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

addSeq48U = function(tib, field,
                     sidx=2, plen=50,
                     verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'addSeq48U'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    idx1 <- sidx
    len1 <- plen - 1
    idx2 <- sidx + 1
    len2 <- plen
    
    field_sym <- rlang::sym(field)

    ret_tib <- tib %>% dplyr::mutate(
      Seq_48U_1=stringr::str_sub(!!field_sym, idx1,len1) %>% 
        stringr::str_to_upper() %>% 
        stringr::str_replace_all('R', 'A') %>% 
        stringr::str_replace_all('Y', 'T'),
      
      Seq_48U_2=stringr::str_sub(!!field_sym, idx2,len2) %>% 
        stringr::str_to_upper() %>% 
        stringr::str_replace_all('R', 'A') %>% 
        stringr::str_replace_all('Y', 'T'),
    )
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

improbe_design_all = function(tib, ptype, outDir, gen,
                              image, shell,
                              seqKey="IUPAC_Sequence",strsSR="FR",reduce_imp=TRUE,
                              sidx=2, plen=50,
                              parse_din=FALSE, del='_',
                              parallel=TRUE,
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'improbe_design_all'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; ptype={ptype}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    idx1 <- sidx
    len1 <- plen - 1
    idx2 <- sidx + 1
    len2 <- plen

    if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
    
    imp_fwd_name <- paste(ptype,gen,"improbe_fwd-seq.tsv.gz", sep=del)
    imp_fwd_path <- file.path(outDir, imp_fwd_name)
    
    # These are NOT used...
    # imp_des_name <- paste(ptype,gen,"improbe-designOutput.tsv.gz", sep=del)
    # imp_des_path <- file.path(outDir, imp_des_name)
    
    imp_fwd_tib <- tib %>% 
      dplyr::select(Seq_ID, Sequence, Genome_Build, Chromosome, Coordinate, CpG_Island)
    readr::write_tsv(imp_fwd_tib, imp_fwd_path)
    
    imp_des_tib <- improbe_docker(
      dir=outDir, file=imp_fwd_name, 
      name=gen, image=image, shell=shell, 
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
      dplyr::mutate(
        Seq_48U_1=stringr::str_sub(UnMethyl_Probe_Sequence, idx1,len1) %>% 
          stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% 
          stringr::str_replace_all('Y', 'T'),
        
        Seq_48U_2=stringr::str_sub(UnMethyl_Probe_Sequence, idx2,len2) %>% 
          stringr::str_to_upper() %>% stringr::str_replace_all('R', 'A') %>% 
          stringr::str_replace_all('Y', 'T'),
        
        Chromosome=as.character(Chromosome),
        Strand_SR=Methyl_Allele_FR_Strand,
        Strand_TB=stringr::str_sub(Methyl_Allele_TB_Strand,1,1),
        Strand_CO=Methyl_Allele_CO_Strand,
        Probe_Score_Min=pmin(Methyl_Final_Score,UnMethyl_Final_Score),
        Underlying_CpG_Count=as.integer(Methyl_Underlying_CpG_Count),
        Underlying_CpG_Min_Dist=as.integer(Methyl_Underlying_CpG_Min_Dist)
        )
    
    if (parse_din) imp_des_tib <- imp_des_tib %>%
      tidyr::separate(Seq_ID, into=c("Seq_ID","DiNuc"), sep=del)
    
    #
    # TBD:: Add better reduced returned value...
    #
    if (reduce_imp) {
      imp_des_tib <- imp_des_tib %>%
        dplyr::select(
          Seq_ID,Forward_Sequence,Genome_Build,Chromosome,Coordinate,Top_Sequence,
          Strand_SR,Strand_TB,Strand_CO,Probe_Score_Min,
          Underlying_CpG_Count,Underlying_CpG_Min_Dist,
          Methyl_Probe_Sequence,UnMethyl_Probe_Sequence,Seq_48U_1,Seq_48U_2)
    }
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} tib={RET}"))
      print(tib)
    }
    
    # Build de-novo IUPAC designs
    #
    iup_des_tib <- desSeq_to_prbs(
      tib=tib,
      # idsKey="Seq_ID",seqKey="IUPAC_Sequence",prbKey="Probe_Type",
      idsKey="Seq_ID",seqKey=seqKey,prbKey="Probe_Type",
      strsSR=strsSR, parallel=parallel,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
      dplyr::mutate(
        
        Seq_48U_1=stringr::str_sub(PRB1_U_MAT, idx1,len1) %>% 
          stringr::str_to_upper() %>% 
          stringr::str_replace_all('R', 'A') %>% 
          stringr::str_replace_all('Y', 'T'),
        
        Seq_48U_2=stringr::str_sub(PRB1_U_MAT, idx2,len2) %>% 
          stringr::str_to_upper() %>% 
          stringr::str_replace_all('R', 'A') %>% 
          stringr::str_replace_all('Y', 'T'),
      )
    
    if (parse_din) iup_des_tib <- iup_des_tib %>%
      tidyr::separate(Seq_ID, into=c("Seq_ID","DiNuc"), sep=del)
    
    if (TRUE) {
      if (!reduce_imp) imp_des_tib <- imp_des_tib %>% 
          dplyr::rename(Probe_Type_IMP=Probe_Type)
      
      ret_tib <- dplyr::inner_join(
        imp_des_tib, iup_des_tib,
        by=c("Seq_ID", "Strand_SR", "Strand_CO",
             "Seq_48U_1", "Seq_48U_2"),
        suffix=c("_IMP", "_IUP")
      )
      ret_cnt <- ret_tib %>% base::nrow()
    } else {
      ret_tib$imp <- imp_des_tib
      ret_tib$iup <- iup_des_tib
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

improbe_docker = function(dir, file, name, image, shell, 
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'improbe_docker'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; name={name}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_log <- 
      file.path(dir, paste(name,'improbe-designOutput.log', sep='.'))
    ret_tsv <- 
      file.path(dir, paste(name,'improbe-designOutput.tsv.gz', sep='.'))
    
    if (file.exists(ret_log)) unlink(ret_log)
    if (file.exists(ret_tsv)) unlink(ret_tsv)
    
    system(glue::glue("touch {ret_log}"))
    system(glue::glue("touch {ret_tsv}"))
    
    # imp_doc_cmd <- glue::glue("docker run -il --rm ",
    
    imp_doc_cmd <- glue::glue("docker run -i --rm ",
                              "-v {dir}:/work -v {dir}:/output ",
                              "{image} {shell} {file} {name}")
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: Running improbe cmd='{imp_doc_cmd}'...{RET}"))
    system(imp_doc_cmd)
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: Done.{RET}{RET}"))
    
    ret_tib <- 
      suppressMessages(suppressWarnings( readr::read_tsv(ret_tsv) ))
    if (verbose>=vt+4) print(ret_tib)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Modern improbe IO Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

isTibIMP = function(tib,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'isTibIMP'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_val <- FALSE
  stime <- system.time({
    
    imp_col <- INIT_IMP_HEADER()
    imp_key <- imp_col$cols %>% names()
    sam_col <- tib %>% names()
    
    if (base::all.equal(imp_key,sam_col)) ret_val <- TRUE
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Matched Header={ret_val}, elapsed={etime}.{RET}{RET}"))
  
  ret_val
}

loadIMP = function(file, max=Inf, format=NULL, del='\t',
                   verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadIMP'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    imp_col <- INIT_IMP_HEADER()
    imp_key <- imp_col$cols %>% names()
    sam_col <- utils::read.table(file, header=TRUE, sep=del, nrows=1) %>% names()
    
    if (max==Inf) {
      if (base::all.equal(imp_key,sam_col)) {
        ret_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file, col_types = imp_col) ))
      } else {
        ret_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) ))
      }
    } else {
      if (!is.numeric(max)) {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: max must be an integer; max={max}!!!{RET}{RET}"))
        return(ret_tib)
      }
      max <- as.integer(max)
      
      if (base::all.equal(imp_key,sam_col)) {
        ret_tib <- suppressMessages(suppressWarnings( 
          utils::read.table(file, header=TRUE, sep=del, nrows=max, colClasses=imp_col) )) %>% 
          tibble::as_tibble() %>% dplyr::mutate(Chromosome=as.character(Chromosome))
      } else {
        ret_tib <- suppressMessages(suppressWarnings( 
          utils::read.table(file, header=TRUE, sep=del, nrows=max) )) %>% tibble::as_tibble()
      }
    }
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Size={ret_cnt}, elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Bowtie Alignment Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadProbeAlignBowtieInfI = function(sam, reduced=FALSE, filtered=FALSE, flipSeq=FALSE,
                                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadProbeAlignBowtieInfI'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} sam={sam}.{RET}"))
  
  col_vec <- c('QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL',
               'AS', 'XN', 'XM', 'XO', 'XG', 'NM', 'MD', 'YT')
  
  stime <- system.time({
    
    snp_raw_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file=sam, col_names=col_vec, comment='@') ))
    
    if (reduced)  snp_raw_tib <- snp_raw_tib %>% dplyr::select(QNAME:POS,SEQ,MD)
    if (filtered) snp_raw_tib <- snp_raw_tib %>% dplyr::filter(FLAG==0 | FLAG==16)
    if (flipSeq)  snp_raw_tib <- snp_raw_tib %>% dplyr::mutate(SEQ=case_when( FLAG==16 ~ revCmp(SEQ),TRUE ~ SEQ ))
    
    # Split by Infinium I Probe Design::
    sam_1A_tib <- snp_raw_tib %>% dplyr::filter(stringr::str_ends(QNAME,'_IA')) %>% dplyr::mutate(QNAME=stringr::str_remove(QNAME, '_IA$'))
    sam_1B_tib <- snp_raw_tib %>% dplyr::filter(stringr::str_ends(QNAME,'_IB')) %>% dplyr::mutate(QNAME=stringr::str_remove(QNAME, '_IB$'))
    sam_tib <- dplyr::inner_join(sam_1A_tib,sam_1B_tib, by=c("QNAME","FLAG","RNAME","POS"), suffix=c("_IA", "_IB")) # %>%

    MD_IA_cnt <- sam_tib %>% dplyr::filter(MD_IA %in% art_snp_md_vec) %>% base::nrow()
    MD_IB_cnt <- sam_tib %>% dplyr::filter(MD_IB %in% art_snp_md_vec) %>% base::nrow()
    if (verbose>=vt) cat(glue::glue("[{par$prgmTag}]: IA/IB = {MD_IA_cnt}, {MD_IB_cnt}{RET}"))
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  sam_tib
}

bowtieProbeAlign = function(exe, fas, gen, dir, lan=NULL, run=FALSE,
                            verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'bowtieProbeAlign'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} dir={dir}.{RET}"))
  
  stime <- system.time({
    fas_name <- fas %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
    gen_name <- gen %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
    gen_file <- gen %>% stringr::str_remove('.gz$')
    
    # Create Shell and Alignment directories::
    sh_dir <- file.path(dir, 'shells')
    al_dir <- file.path(dir, 'align')
    
    if (!dir.exists(sh_dir)) dir.create(sh_dir, recursive=TRUE)
    if (!dir.exists(al_dir)) dir.create(al_dir, recursive=TRUE)
    
    aln_ssh <- file.path(sh_dir, paste0('run_bow-',fas_name,'-',gen_name,'.sh') )
    aln_sam <- file.path(al_dir, paste0(fas_name,'-',gen_name,'.bowtie.sam.gz') )
    aln_cmd <- paste(exe, '-f -x',gen_file, '-U',fas, '| gzip -c ->',aln_sam, sep=' ')

    cat(glue::glue("[{funcTag}]:{TAB} Launching bowtie alignments: {fas_name} vs. {gen_name}...{RET}"))
    readr::write_lines(x=aln_cmd, file=aln_ssh, append=FALSE)
    Sys.chmod(paths=aln_ssh, mode="0777")
    if (!is.null(lan)) aln_ssh <- paste(lan,aln_ssh, sep=' ')
    if (run) base::system(aln_ssh)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))

  if (run) return(aln_sam)
  aln_ssh
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Bowtie Alignment Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

bsmapProbeAlign = function(exe, fas, gen, dir, 
                           opt="-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R", 
                           lan=NULL, run=FALSE,
                           verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'bsmapProbeAlign'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} dir={dir}.{RET}"))
  
  stime <- system.time({
    if (is.null(opt) || length(opt)==0)
      opt="-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
    
    fas_name <- fas %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
    gen_name <- gen %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
    
    # Create Shell and Alignment directories::
    sh_dir <- file.path(dir, 'shells')
    al_dir <- file.path(dir, 'align')
    
    if (!dir.exists(sh_dir)) dir.create(sh_dir, recursive=TRUE)
    if (!dir.exists(al_dir)) dir.create(al_dir, recursive=TRUE)

    aln_ssh <- file.path(sh_dir, paste0('run_bsp-',fas_name,'-',gen_name,'.sh') )
    aln_bsp <- file.path(al_dir, paste0(fas_name,'-',gen_name,'.bsmap.bsp') )
    ret_tsv <- file.path(al_dir, paste0(fas_name,'-',gen_name,'.bsmap.formatted.tsv.gz') )
    aln_cmd <- paste(exe, '-a',fas, '-d',gen, opt,'-o',aln_bsp, sep=' ')
    aln_cmd <- glue::glue("{aln_cmd}{RET}",
                          "gzip {aln_bsp}{RET}",
                          "gzip -dc {aln_bsp}.gz | cut -f 1,2,4-11 | perl -pe 's/_/\t/gi; s/:/\t/gi' | ",
                          "gzip -c - > {ret_tsv}{RET}")
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{TAB} Launching; CMD={aln_cmd}...{RET}"))
    
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{TAB} Launching bsmap alignments: {fas_name} vs. {gen_name}...{RET}"))
    readr::write_lines(x=aln_cmd, file=aln_ssh, append=FALSE)
    Sys.chmod(paths=aln_ssh, mode="0777")
    
    if (!is.null(lan)) aln_ssh <- paste(lan,aln_ssh, sep=' ')
    if (run) base::system(aln_ssh)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  if (run) return(ret_tsv)
  aln_ssh
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Infinium Methylation Probe toStrings::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tibToFas = function(tib, key,seq, prefix,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'tibToFas'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; key={key}, seq={seq}, prefix={prefix}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    dir <- base::dirname(prefix)
    if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    
    file_fas <- paste0(prefix, ".fa.gz")
    file_csv <- paste0(prefix, ".csv.gz")
    
    key_sym <- rlang::sym(key)
    seq_sym <- rlang::sym(seq)
    
    ret_tib <- tib %>% dplyr::mutate(
      # fas_line=paste0(">",!!key, "\n",!!seq)
      fas_line=paste0(">",!!key_sym, "\n",!!seq_sym)
    )  %>%
      dplyr::filter(!is.na(fas_line))
    
    fas_vec <- ret_tib %>%
      # dplyr::distinct(fas_line) %>%
      dplyr::pull(fas_line)
    
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Writing Data CSV={file_csv}...{RET}"))
    readr::write_csv(ret_tib, file_csv)
    
    if (opt$verbose>=1)
      cat(glue::glue("[{par$prgmTag}]: Writing Data FAS={file_fas}...{RET}"))
    readr::write_lines(x=fas_vec, file=file_fas)

    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib({ret_cnt})={RET}"))
      ret_tib %>% print()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

ordToFas = function(tib, dir, name, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'writeBedFas'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}, dir={dir}.{RET}"))
  
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  fas_file <- file.path(dir, paste0(name,'.cgn.fa.gz'))
  
  fas_tib <- tib %>% dplyr::mutate(
    line=dplyr::case_when(
      Normalization_Bin=='C' ~ paste0('>',AlleleA_Probe_Id,'\n',AlleleA_Probe_Sequence),
      Normalization_Bin!='C' ~ paste0('>',AlleleA_Probe_Id,'\n',AlleleA_Probe_Sequence,'\n',
                                      '>',AlleleB_Probe_Id,'\n',AlleleB_Probe_Sequence),
      TRUE ~ NA_character_ ) ) %>%
    dplyr::filter(!is.na(line)) %>% dplyr::pull(line)
  
  readr::write_lines(x=fas_tib, file=fas_file)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  fas_file
}

writeBedFas = function(tib, dir, name, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'writeBedFas'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}, dir={dir}.{RET}"))
  
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  fas_file <- file.path(dir, paste0(name,'.fa.gz'))
  bed_file <- file.path(dir, paste0(name,'.sorted.bed.gz'))
  
  tib <- tib %>% dplyr::arrange(-MinScore) %>% dplyr::distinct(PRB1_U, .keep_all=TRUE) %>%
    dplyr::arrange(Chromosome,Coordinate) %>% dplyr::rename(Chrom=Chromosome, Beg=Coordinate) %>%
    dplyr::mutate(End=Beg+stringr::str_length(PRB1_U), 
                  Full_ID=paste(Seq_ID,FR_Imp,TB_Imp,CO_Imp,Beg, sep='_') )

  fas_tib <- tib %>% dplyr::select(Full_ID, PRB1_U) %>% 
    dplyr::mutate(line=paste0('>',Full_ID,'\n',PRB1_U)) %>% dplyr::pull(line)
  bed_tib <- tib %>% dplyr::select(Chrom, Beg, End, Full_ID, MinScore, FR_Imp, TB_Imp, CO_Imp, everything())
  
  readr::write_lines(x=fas_tib, file=fas_file)
  readr::write_tsv(bed_tib, bed_file)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  bed_tib
}

printPrbs = function(tib, pr='cg', org=NULL, outDir, plotName, max=NULL,
                     verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'printPrbs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; pr={pr}, plotName={plotName}, outDir={outDir}.{RET}"))
  
  prb_mat_tibs <- tib %>% dplyr::filter(PRB_DES==pr) %>% dplyr::distinct()
  prb_mat_cnt  <- prb_mat_tibs %>% base::nrow()
  
  plot_ord_tib <- prb_mat_tibs %>% dplyr::distinct(Seq_ID, .keep_all=TRUE)
  plot_ord_cnt <- plot_ord_tib %>% base::nrow()
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} plot_ord_cnt={plot_ord_cnt}, prb_mat_cnt={prb_mat_cnt}.{RET}"))
  
  if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
  
  passFile <- file.path(outDir, paste(plotName,pr,'pass.txt', sep='_') )
  failFile <- file.path(outDir, paste(plotName,pr,'fail.txt', sep='_') )
  unlink(passFile)
  unlink(failFile)
  
  if (!is.null(org)) org <- org %>% dplyr::distinct()
  
  tibs <- NULL
  for (ii in seq(1,plot_ord_cnt)) {
    tag_tib <- plot_ord_tib[ii,]
    cur_SeqID <- tag_tib$Seq_ID[1]
    # print(tag_tib)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ii={ii}; cur_SeqID={cur_SeqID}.{RET}"))
    
    cur_org <- NULL
    org_str <- ''
    if (!is.null(org)) {
      cur_org <- org %>% dplyr::filter(Seq_ID==cur_SeqID) %>% head(n=1)
      
      if (FALSE) {
        for (jj in c(1:length(cur_org))) {
          org_str[jj] <- stringr::str_c(cur_org[jj], collapse=',')
          # org_str[jj] <- paste(cur_org[jj], collapse=', ')
          break
        }
      }
      org_str <- cur_org %>% paste(collapse=', ')
    }
    # org_str <- stringr::str_c(org_str, collapse=';')
    # org_str <- stringr::str_c(unique(org_str), collapse=';')
    # print(org_str)

    cur_tib <- prb_mat_tibs %>% dplyr::filter(Seq_ID==cur_SeqID)
    tibs <- cur_tib %>% srdsToBrac()
    # print(cur_tib)
    # return(tibs)
    # cat("\n\n\n\n\n")
    # cat("TIBS::\n")
    # print(tibs)
    # cat("\n\n\n\n\n")
    
    # Build Printable Strings for reach strand::
    strs <- NULL
    strs$F <- prbsToStr(tibs$F, pr=tag_tib$PRB_DES[1], verbose=verbose, vt=6)
    strs$R <- prbsToStr(tibs$R, pr=tag_tib$PRB_DES[1], verbose=verbose, vt=6)

    # Determine which file to write to::
    #   - Failed if both F/R strands have equal PRB1_U==PRB1_M
    #   - Passed otherwise
    fwd_cnt <- tibs$F %>% dplyr::filter(PRB1_U!=PRB1_M) %>% base::nrow()
    rev_cnt <- tibs$R %>% dplyr::filter(PRB1_U!=PRB1_M) %>% base::nrow()
    
    space_str <- paste(stringr::str_pad("# ", width=150, side="right", pad="-"), "\n", sep='')
    seq_id_str <- ''
    for (jj in c(1:length(org_str))) {
      # seq_id_str <- paste0("Probe_Type=",pr,"; Seq_ID: ", paste(unique(cur_tib$Seq_ID), collapse='\t'),"; Original: ",org_str,"\n")
      
      seq_id_str <- paste0("Probe_Type=",pr,"; Seq_ID: ", paste(unique(cur_tib$Seq_ID), collapse='\t'),"; Original: ",org_str[jj],"\n")
      # seq_id_str <- paste0("Probe_Type=",pr,"; Seq_ID: ", paste(unique(cur_tib$Seq_ID), collapse='\t'),"\n")
    }
    out_lines <- NULL
    out_lines[1] <- space_str
    out_lines[2] <- seq_id_str
    out_lines[3] <- stringr::str_c(strs$F)
    out_lines[4] <- stringr::str_c(strs$R)
    
    if (verbose>=vt) cat(out_lines)
    
    if (fwd_cnt>0 || rev_cnt>0) {
      readr::write_lines(out_lines, passFile, append=TRUE)
    } else {
      readr::write_lines(out_lines, failFile, append=TRUE)
    }
    
    # This was done for known comparisons::
    #
    if (FALSE) {
      if (FALSE) { # && cur_tib$Infinium_Design_Type[1]=='II') {
        cur_tib %>% dplyr::select(FR,CO, PRB2_D_IUP, PRB2_D_IMP,Man_MisMatch,Man_TarMatch,Bad_Design) %>% print()
      } else {
        prb1U_iup <- cur_tib %>% dplyr::select(FR,CO, PRB1_U_IUP, PRB1_U_IMP,Man_MisMatch,Man_TarMatch,Bad_Design)
        prb1M_iup <- cur_tib %>% dplyr::select(FR,CO, PRB1_M_IUP, PRB1_M_IMP,Man_MisMatch,Man_TarMatch,Bad_Design)
        # readr::write_file(str_c(prb1U_iup), passFile, append=TRUE)
        # readr::write_file(str_c(prb1M_iup), passFile, append=TRUE)
      }
    }
    
    if (!is.null(max) && ii>=max) break
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  tibs
}

srdsToBrac = function(tib, 
                      beg1=1, end1=60, mid1=61,
                      beg2=63,end2=122,mid2=62) {
  # TBD:: Calculate all data points based on sequence length
  
  tib %>% 
    dplyr::mutate(StrandFR=case_when(FR ~ 'F', !FR ~ 'R', TRUE ~ NA_character_),
                  StrandCO=case_when(CO ~ 'C', !CO ~ 'O', TRUE ~ NA_character_),
                  DesSeqN=paste0(stringr::str_sub(DesSeqN,beg1,end1),{BNG},
                                 stringr::str_sub(DesSeqN,mid1,mid2),{BNG},
                                 stringr::str_sub(DesSeqN,beg2,end2)),
                  
                  DesBscU=paste0(stringr::str_sub(DesBscU,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscU,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscU,beg2,end2)),
                  
                  DesBscM=paste0(stringr::str_sub(DesBscM,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscM,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscM,beg2,end2)),
                  
                  DesBscD=paste0(stringr::str_sub(DesBscD,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscD,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscD,beg2,end2))) %>%
    dplyr::arrange(StrandFR, StrandCO) %>% split(.$StrandFR)
  
}

prbsToStr = function(tib, pr, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'prbsToStr'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; pr={pr}.{RET}"))
  
  # fr1Key <- tib %>% dplyr::select(!!frKey) %>% head(n=1) %>% pull()
  # fr2Key <- tib$StrandFR[2]
  # co1Key <- tib$StrandCO[1]
  # co2Key <- tib$StrandCO[2]
  
  fr1Key <- tib$StrandFR[1]
  fr2Key <- tib$StrandFR[2]
  co1Key <- tib$StrandCO[1]
  co2Key <- tib$StrandCO[2]
  
  if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} fr1Key={fr1Key}, fr2Key={fr2Key}, co1Key={co1Key}, co2Key={co2Key}.{RET}"))
  
  mud <- list('U'='U','M'='M','D'='D')

  desSeq <- 'DesSeqN'
  bscKey <- lapply(mud, function(x) { paste('DesBsc',x,sep='')} )
  nxbKey <- lapply(mud, function(x) { paste('NXB',x,sep='_')} )
  cpnKey <- lapply(mud, function(x) { paste('CPN',x,sep='_')} )
  tarKey <- lapply(mud, function(x) { paste('TAR',x,sep='_')} )
  secKey <- lapply(mud, function(x) { paste('SEC',x,sep='_')} )
  bodKey <- lapply(mud, function(x) { paste('BOD',x,sep='_')} )
  endKey <- lapply(mud, function(x) { paste('END',x,sep='_')} )

  # Dertermine if the probes on this strand were designable::
  pass_cnt <- tib %>% dplyr::filter(PRB1_U!=PRB1_M) %>% base::nrow()
  pass_str <- 'PASS: '
  if (pass_cnt==0) pass_str <- 'FAIL: '
  
  # TBD:: Note on the Opposite Strand we should reverse all framents, but currently fragLen==1 are left alone for effiecntcy...
  # Sketch Output::
  #
  # F_C_N    DesSeqN[CG]DesSeqN
  #
  #                  D2 22222
  #                N M1 1111
  #                N U1 1111
  # F_C_U    DesBscU[tG]DesBscU
  # F_O_U    DesBscU[Ca]DesBscU
  #             1111 1U N
  #             1111 1M N
  #            22222 2D
  #
  if (pr=='rs'||pr=='ch') {
    if (fr1Key=='F' && fr2Key=='F') {
      bufC <- 0
      bufO <- 0
      str <- glue::glue(
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 61+bufC), collapse=''),"{tib[[tarKey$D]][1]}{tib[[secKey$D]][1]}{BNG}{tib[[bodKey$D]][1]}{tib[[endKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 59+bufC), collapse=''),"{tib[[nxbKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[secKey$M]][1]}{BNG}{tib[[bodKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 59+bufC), collapse=''),"{tib[[nxbKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[secKey$U]][1]}{BNG}{tib[[bodKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_U {TAB}{tib[[bscKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_M {TAB}{tib[[bscKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{tib[[bscKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{tib[[desSeq]][1]}{RET}",
        # "FwdSeq{TAB}{tib$Forward_Sequence[1]}{RET}",
        "{RET}",
        "{pass_str}{fr1Key}_{co2Key}_N {TAB}{cmpl(tib[[desSeq]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_M {TAB}{Biostrings::reverse(tib[[bscKey$M]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_U {TAB}{Biostrings::reverse(tib[[bscKey$U]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 11-bufO), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][2])}{tib[[secKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[nxbKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 11-bufO), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][2])}{tib[[secKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[nxbKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 10-bufO), collapse=''),"{tib[[endKey$D]][2]}{Biostrings::reverse(tib[[bodKey$D]][2])}{tib[[secKey$D]][2]}{BNG}{tib[[tarKey$D]][2]}{RET}",
        "{RET}")
    } else if (fr1Key=='R' && fr2Key=='R') {
      bufO <- 0
      bufC <- 0
      str <- glue::glue(
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 61+bufO), collapse=''),"{tib[[tarKey$D]][2]}{tib[[secKey$D]][2]}{BNG}{tib[[bodKey$D]][2]}{tib[[endKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 59+bufO), collapse=''),"{tib[[nxbKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[secKey$M]][2]}{BNG}{tib[[bodKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 59+bufO), collapse=''),"{tib[[nxbKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[secKey$U]][2]}{BNG}{tib[[bodKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{tib[[bscKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_N {TAB}{revCmp(tib[[desSeq]][2])}{RET}",
        # "{fr2Key}_{co2Key}_N {TAB}{tib[[desSeq]][2]}{RET}",
        "{RET}",
        # "{fr1Key}_{co1Key}_N {TAB}{cmpl(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{Biostrings::reverse(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 11-bufC), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][1])}{tib[[secKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[nxbKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 11-bufC), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][1])}{tib[[secKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[nxbKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 10-bufC), collapse=''),"{tib[[endKey$D]][1]}{Biostrings::reverse(tib[[bodKey$D]][1])}{tib[[secKey$D]][1]}{BNG}{tib[[tarKey$D]][1]}{RET}",
        "{RET}")
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: fr1Key={fr1Key}, fr2Key={fr2Key}, Allowed Values=[F,R]!{RET}{RET}"))
    }
  } else {
    if (fr1Key=='F' && fr2Key=='F') {
      buf <- 0
      str <- glue::glue(
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 61+buf), collapse=''),"{tib[[tarKey$D]][1]}{tib[[secKey$D]][1]}{BNG}{tib[[bodKey$D]][1]}{tib[[endKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 59+buf), collapse=''),"{tib[[nxbKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[secKey$M]][1]}{BNG}{tib[[bodKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 59+buf), collapse=''),"{tib[[nxbKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[secKey$U]][1]}{BNG}{tib[[bodKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{tib[[bscKey$D]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{tib[[desSeq]][1]}{RET}",
        "{RET}",
        "{pass_str}{fr1Key}_{co2Key}_N {TAB}{cmpl(tib[[desSeq]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][2])}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 11-buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][2])}{tib[[secKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[nxbKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 11-buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][2])}{tib[[secKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[nxbKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 10-buf), collapse=''),"{tib[[endKey$D]][2]}{Biostrings::reverse(tib[[bodKey$D]][2])}{tib[[secKey$D]][2]}{BNG}{tib[[tarKey$D]][2]}{RET}",
        "{RET}")
    } else if (fr1Key=='R' && fr2Key=='R') {
      buf <- 0
      buf <- 1
      str <- glue::glue(
        "{pass_str}{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 61+buf), collapse=''),"{tib[[tarKey$D]][2]}{BNG}{tib[[secKey$D]][2]}{tib[[bodKey$D]][2]}{tib[[endKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 60+buf), collapse=''),"{tib[[nxbKey$M]][2]}{tib[[tarKey$M]][2]}{BNG}{tib[[secKey$M]][2]}{tib[[bodKey$M]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 60+buf), collapse=''),"{tib[[nxbKey$U]][2]}{tib[[tarKey$U]][2]}{BNG}{tib[[secKey$U]][2]}{tib[[bodKey$U]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_D {TAB}{tib[[bscKey$D]][2]}{RET}",
        "{pass_str}{fr2Key}_{co2Key}_N {TAB}{revCmp(tib[[desSeq]][2])}{RET}",
        # "{fr2Key}_{co2Key}_N {TAB}{tib[[desSeq]][2]}{RET}",
        "{RET}",
        # "{fr1Key}_{co1Key}_N {TAB}{cmpl(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_N {TAB}{Biostrings::reverse(tib[[desSeq]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][1])}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 11+buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][1])}{BNG}{tib[[secKey$U]][1]}{tib[[tarKey$U]][1]}{BNG}{tib[[nxbKey$U]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 11+buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][1])}{BNG}{tib[[secKey$M]][1]}{tib[[tarKey$M]][1]}{BNG}{tib[[nxbKey$M]][1]}{RET}",
        "{pass_str}{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 10+buf), collapse=''),"{tib[[endKey$D]][1]}{Biostrings::reverse(tib[[bodKey$D]][1])}{BNG}{tib[[secKey$D]][1]}{tib[[tarKey$D]][1]}{RET}",
        "{RET}")
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: fr1Key={fr1Key}, fr2Key={fr2Key}, Allowed Values=[F,R]!{RET}{RET}"))
    }
  }
  if (verbose>=vt) cat(str)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  str
}

prbsToStrMUD = function(tib, mu='U', verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'prbsToStr'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  fr1Key <- tib$StrandFR[1]
  fr2Key <- tib$StrandFR[2]
  
  co1Key <- tib$StrandCO[1]
  co2Key <- tib$StrandCO[2]
  
  bscKey <- paste('DesBsc',mu, sep='')
  
  nxbKey <- paste('NXB',mu, sep='_')
  cpnKey <- paste('CPN',mu, sep='_')
  tarKey <- paste('TAR',mu, sep='_')
  secKey <- paste('SEC',mu, sep='_')
  bodKey <- paste('BOD',mu, sep='_')
  endKey <- paste('END',mu, sep='_')
  
  str <- glue::glue(
    paste0(rep(" ", 59+10), collapse=''),"{tib[[cpnKey]][1]}{tib[[secKey]][1]}{BNG}{tib[[bodKey]][1]}{tib[[endKey]][1]}{RET}",
    paste0(rep(" ", 59+8),  collapse=''),"{tib[[nxbKey]][1]}{BNG}{tib[[cpnKey]][1]}{tib[[secKey]][1]}{BNG}{tib[[bodKey]][1]}{RET}",
    "{fr1Key}_{co1Key}_{mu}{TAB}{tib[[bscKey]][1]}{RET}",
    "{fr2Key}_{co2Key}_{mu}{TAB}{Biostrings::reverse(tib[[bscKey]][2])}{RET}",
    paste0(rep(" ", 9+10), collapse=''),"{Biostrings::reverse(tib[[bodKey]][2])}{tib[[secKey]][2]}{BNG}{tib[[cpnKey]][2]}{tib[[nxbKey]][2]}{RET}",
    paste0(rep(" ", 9+9),  collapse=''),"{tib[[endKey]][2]}{Biostrings::reverse(tib[[bodKey]][2])}{tib[[secKey]][2]}{BNG}{tib[[cpnKey]][2]}{RET}",
    "{RET}")
  
  if (verbose>=vt) cat(str)
  
  str
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Infinium Methylation Probe Design Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Function to design all probes from a generic tibble
#  old name: tib2prbs
desSeq_to_prbs = function(tib, idsKey,seqKey,prbKey, strsSR='FR', strsCO='CO',
                          addMatSeq=TRUE, parallel=FALSE,del='_',max=0,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'desSeq_to_prbs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; SR={strsSR}, CO={strsCO}, ",
                                  "idsKey={idsKey}, prbKey={prbKey}, seqKey={seqKey}.{RET}"))
  
  # Unambigous Source Design Sequence Strand
  sr_vec <- stringr::str_split(strsSR, '', simplify=TRUE) %>% as.vector()
  co_vec <- stringr::str_split(strsCO, '', simplify=TRUE) %>% as.vector()
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    idsKey <- rlang::sym(idsKey)
    prbKey <- rlang::sym(prbKey)
    seqKey <- rlang::sym(seqKey)
    
    if (typeof(tib)=='character') {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading from file={tib}...{RET}"))
      if (stringr::str_ends(tib,'.tsv.gz') || stringr::str_ends(tib,'.tsv')) {
        tib <- suppressMessages(suppressWarnings( readr::read_tsv(tib) ))
      } else {
        tib <- suppressMessages(suppressWarnings( readr::read_csv(tib) ))
      }
    }
    
    if (max>0) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Will subset input to max={max}.{RET}"))
      tib <- tib %>% head(n=max)
    }
    
    src_man_tib <- tib %>% dplyr::select(!!idsKey, !!prbKey, !!seqKey) %>% 
      dplyr::mutate(Seq_ID:=!!idsKey, PRB_DES:=!!prbKey)
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} src_man_tib(original)={RET}"))
    if (verbose>=vt+4) print(src_man_tib)
    
    # Ensure we have 122 mer format 60[NN]60
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Validating design sequneces...{RET}"))
    src_man_tib <- src_man_tib %>%
      dplyr::mutate(!!seqKey := stringr::str_replace(!!seqKey, '\\[','_') %>% stringr::str_replace('\\]','_')
      ) %>%
      tidyr::separate(!!seqKey, into=c("PreSeqN", "MidSeqN", "PosSeqN"), sep='_') %>%
      dplyr::mutate(PreSeqN=stringr::str_sub(PreSeqN,   -60),
                    PosSeqN=stringr::str_sub(PosSeqN, 1, 60),
                    PreSeqN=stringr::str_pad(string=PreSeqN, width=60, side='left', pad='N'),
                    PosSeqN=stringr::str_pad(string=PosSeqN, width=60, side='right', pad='N'),
                    DesNucA=stringr::str_sub(MidSeqN, 1,1), DesNucB=stringr::str_sub(MidSeqN, 2,2),
                    !!seqKey :=paste0(PreSeqN,'[',MidSeqN,']',PosSeqN) )
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} src_man_tib(122 check)={RET}"))
    if (verbose>=vt+4) print(src_man_tib)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                      Calcluate Probes on All Strands::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building forward & reverse design sequneces...{RET}"))
    
    des_seq_F_C <- src_man_tib %>% 
      dplyr::mutate(SR=TRUE, CO=TRUE, DesSeqN=shearBrac(!!seqKey))
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} des_seq_F_C={RET}"))
    if (verbose>=vt+4) print(des_seq_F_C)
    
    des_seq_R_C <- des_seq_F_C %>% dplyr::mutate(
      SR=!SR,CO=CO, DesSeqN=revCmp(DesSeqN) )
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} des_seq_R_C={RET}"))
    if (verbose>=vt+4) print(des_seq_R_C)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Build All Design Strands::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building bisulfite converted design strands...{RET}"))
    
    bsc_tibs <- NULL
    
    # BSC-Forward-Converted::
    bsc_tibs$FC <- des_seq_F_C %>% dplyr::mutate(
      DesBscU = bscUs(DesSeqN),
      DesBscM = bscMs(DesSeqN),
      DesBscD = bscDs(DesSeqN) )
    if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Done.Building bisulfite(FC).{RET}{RET}"))
    if (verbose>=vt+4) bsc_tibs$FC %>% head(n=3) %>% print()
    
    # BSC-Foward-Opposite::
    bsc_tibs$FO <- bsc_tibs$FC %>% dplyr::mutate(
      SR=SR,CO=!CO,
      DesBscU=revCmp(DesBscU),
      DesBscM=revCmp(DesBscM),
      DesBscD=revCmp(DesBscD) )
    if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Done.Building bisulfite(FO).{RET}{RET}"))
    if (verbose>=vt+4) bsc_tibs$FO %>% head(n=3) %>% print()
    
    # BSC-Reverse-Converted::
    bsc_tibs$RC <- des_seq_R_C %>% dplyr::mutate(
      DesBscU = bscUs(DesSeqN),
      DesBscM = bscMs(DesSeqN),
      DesBscD = bscDs(DesSeqN) )
    if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Done.Building bisulfite(RC).{RET}{RET}"))
    if (verbose>=vt+4) bsc_tibs$RC %>% head(n=3) %>% print()
    
    # BSC-Reverse-Opposite::
    bsc_tibs$RO <- bsc_tibs$RC %>% dplyr::mutate(
      SR=SR,CO=!CO,
      DesBscU=revCmp(DesBscU),
      DesBscM=revCmp(DesBscM),
      DesBscD=revCmp(DesBscD) )
    if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Done.Building bisulfite(RO).{RET}{RET}"))
    if (verbose>=vt+4) bsc_tibs$RO %>% head(n=3) %>% print()
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.Building bisulfite converted design strands.{RET}{RET}"))
    # return(bsc_tibs)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Build all Probes on Each Strand::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    srd_names <- names(bsc_tibs)
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} Strand Names={RET}"))
    if (verbose>=vt+4) print(srd_names)
    
    if (parallel) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Build probes for each strand (Parallel)...{RET}"))
      ret_tib <- foreach (srd=srd_names, .combine=rbind) %dopar% {
        lapply(split(bsc_tibs[[srd]], bsc_tibs[[srd]]$PRB_DES), desAllPrbs, 
               verbose=verbose, vt=vt+1, tc=tc+1) %>% dplyr::bind_rows()
      }
    } else {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Build probes for each strand (Linear)...{RET}"))
      for (srd in srd_names) {
        ret_tib <- ret_tib %>% dplyr::bind_rows(
          lapply(split(bsc_tibs[[srd]], bsc_tibs[[srd]]$PRB_DES), desAllPrbs, 
                 verbose=verbose, vt=vt+1, tc=tc+1) %>% dplyr::bind_rows() )
      }
    }
    
    # Update Keys::
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Updating keys and strands.{RET}"))
    ret_tib <- ret_tib %>% dplyr::mutate(
      Strand_SR=case_when(SR ~ sr_vec[1], !SR ~ sr_vec[2], TRUE ~ NA_character_),
      Strand_CO=case_when(CO ~ co_vec[1], !CO ~ co_vec[2], TRUE ~ NA_character_),
      Seq_ID_Uniq=paste(Seq_ID,paste0(Strand_SR,Strand_CO), sep=del)
    ) %>% dplyr::arrange(Seq_ID_Uniq)
    
    # OLD Naming Scheme::
    #
    # ret_tib <- ret_tib %>% dplyr::mutate(
    #   SR_Str=case_when(SR ~ sr_vec[1], !SR ~ sr_vec[2], TRUE ~ NA_character_),
    #   CO_Str=case_when(CO ~ co_vec[1], !CO ~ co_vec[2], TRUE ~ NA_character_),
    #   Seq_ID_Uniq=paste(Seq_ID,paste0(SR_Str,CO_Str), sep=del)
    # ) %>% dplyr::arrange(Seq_ID_Uniq)
    
    # Add match probe sequences
    if (addMatSeq) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding match seqqunes.{RET}"))
      ret_tib <- ret_tib %>% 
        dplyr::mutate(PRB1_U_MAT=stringr::str_to_upper(PRB1_U),
                      PRB1_M_MAT=stringr::str_to_upper(PRB1_M),
                      PRB2_D_MAT=stringr::str_to_upper(PRB2_D) )
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tib({ret_cnt})::{RET}"))
    if (verbose>=vt+4) print(ret_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# Function to design all probes in a single call
desAllPrbs = function(tib, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'desAllPrbs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  ret_tib <- NULL
  stime <- system.time({
    
    fr <- tib %>% dplyr::distinct(SR) %>% base::as.logical()
    co <- tib %>% dplyr::distinct(CO) %>% base::as.logical()
    pr <- tib %>% dplyr::distinct(PRB_DES) %>% base::as.character()
    
    ret_tib <- dplyr::bind_rows(tib %>% 
                       des2prbs(fwd=fr, con=co, pr=pr, mu='U', desSeq='DesBscU', verbose=verbose, vt=vt, tc=tc+1) %>%
                       des2prbs(fwd=fr, con=co, pr=pr, mu='M', desSeq='DesBscM', verbose=verbose, vt=vt, tc=tc+1) %>%
                       des2prbs(fwd=fr, con=co, pr=pr, mu='D', desSeq='DesBscD', verbose=verbose, vt=vt, tc=tc+1) )
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tib({ret_cnt})::{RET}"))
    if (verbose>=vt+4) print(ret_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# Function to design probes from single design strand and orientations::
#   TBD:: Previous default was 'QC_CPN=TRUE' Not sure if that is needed...
#
des2prbs = function(tib, fwd, con, pr, mu, desSeq='DesSeqN', len=48, del='_',QC_CPN=FALSE,
                    verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'des2prbs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; fwd={fwd}, con={con}, pr={pr}, mu={mu}, desSeq={desSeq}.{RET}"))

  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({

    stopifnot(is.logical(fwd))
    stopifnot(is.logical(con))
    
    if (mu!='N' && mu!='U' && mu!='M' && mu!='D')
      stop(glue::glue("{RET}[{funcTag}]: ERROR: mu={mu} Only Supported=[N,U,M,D]!{RET}{RET}"))
    if (pr!='cg' && pr!='ch' & pr!='rs' & pr!='rp')
      stop(glue::glue("{RET}[{funcTag}]: ERROR: pr={pr} Only Supported=[cg,ch,rp,rs]!{RET}{RET}"))
    
    desSeq <- rlang::sym(desSeq)
    if (pr=='rs') {
      # $prb_F_C_I   = revCmp(substr($des_F_C, 60, 50));
      # $prb_R_C_I   = revCmp(substr($des_R_C, 61, 50));
      # $prb_F_C_II  = revCmp(substr($des_F_C, 60, 50));
      # $prb_R_C_II  = revCmp(substr($des_R_C, 61, 50));
      #
      # $prb_F_O_I   = revCmp(substr($des_F_O, 61, 50));
      # $prb_R_O_I   = revCmp(substr($des_R_O, 62, 50));
      # $prb_F_O_II  = revCmp(substr($des_F_O, 61, 50));
      # $prb_R_O_II  = revCmp(substr($des_R_O, 62, 50));
      if      ( fwd &&  con) nxb_pos <- 60
      else if (!fwd &&  con) nxb_pos <- 61
      else if ( fwd && !con) nxb_pos <- 61
      else if (!fwd && !con) nxb_pos <- 60
      # else if (!fwd && !con) nxb_pos <- 62  # Original 
      else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: unsupported combination fwd={fwd}, con={con}!!!{RET}{RET}"))
      }
    } else if (pr=='ch') {
      # Originally thsi was identical to rs format above, but for forward sequences needs to be shifted
      #  upstream for converted and downstream for opposite::
      if      ( fwd &&  con) nxb_pos <- 61 # Previously = 60
      else if (!fwd &&  con) nxb_pos <- 61
      else if ( fwd && !con) nxb_pos <- 60 # Previously = 61
      else if (!fwd && !con) nxb_pos <- 60
      else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: unsupported combination fwd={fwd}, con={con}!!!{RET}{RET}"))
      }
      
      #
      # NEW:: CpH Code::
      #
      if      ( fwd &&  con) nxb_pos <- 60 # Previously = 60
      else if (!fwd &&  con) nxb_pos <- 61
      else if ( fwd && !con) nxb_pos <- 61 # Previously = 61
      else if (!fwd && !con) nxb_pos <- 60
      else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: unsupported combination fwd={fwd}, con={con}!!!{RET}{RET}"))
      }

    } else if (pr=='cg' || pr=='rp' || pr=='mu' || stringr::str_starts(pr,'ct')) {
      # $prb_F_C_I  = revCmp(substr($des_F_C, 60, 50));
      # $prb_F_O_I  = revCmp(substr($des_F_O, 61, 50));
      # $prb_R_C_I  = revCmp(substr($des_R_C, 60, 50));
      # $prb_R_O_I  = revCmp(substr($des_R_O, 61, 50));
      # 
      # $prb_F_C_II  = revCmp(substr($des_F_C, 61, 50));
      # $prb_F_O_II  = revCmp(substr($des_F_O, 62, 50));
      # $prb_R_C_II  = revCmp(substr($des_R_C, 61, 50));
      # $prb_R_O_II  = revCmp(substr($des_R_O, 62, 50));
      nxb_pos <- 60
      if (!con) nxb_pos <- 61
      
    } else {
      stop(glue::glue("[{funcTag}]: ERROR: Probe_Type={pr} is currently not supported!!!{RET}{RET}"))
      
    }
    cpg_pos <- nxb_pos + 1
    sec_pos <- cpg_pos + 1
    bod_pos <- sec_pos + 1
    end_pos <- bod_pos + len
    
    # Special consideration is needed for U/M strands at the query site. 
    #  For CN (i.e. cg or ch) this is actually done naturally in U/M conversion
    #  However, for  non-CN probes (i.e. rs) this needs to be forced to U/M
    #
    # This is handled by the TAR (Target/Query Nucleotide). This should only change
    #  for U/M (QMAP_U/QMAP_M) for D its just itself.
    #
    ret_tib <- tib %>% dplyr::mutate(
      NXB=stringr::str_sub(!!desSeq, nxb_pos, nxb_pos),
      CPN=stringr::str_sub(!!desSeq, cpg_pos, cpg_pos),
      TAR=qmaps(CPN, mu=mu),
      SEC=stringr::str_sub(!!desSeq, sec_pos, sec_pos),
      BOD=stringr::str_sub(!!desSeq, bod_pos, end_pos-1),
      END=stringr::str_sub(!!desSeq, end_pos, end_pos)
    )
    
    #  QC TEST:: for CpN (cg or ch) verify that the probes are equal. Well call this
    #   PRB0 (CGN) and PRB1 (TAR). After testing remove PRB0
    #
    if (QC_CPN && (pr=='cg')) {
      ret_tib <- ret_tib %>%
        tidyr::unite(PRB0, CPN,SEC,BOD, sep='', remove=FALSE) %>%
        tidyr::unite(PRB1, TAR,SEC,BOD, sep='', remove=FALSE) %>%
        tidyr::unite(PRB2, SEC,BOD,END, sep='', remove=FALSE) %>%
        dplyr::mutate(PRB0=revCmp(PRB0), PRB1=revCmp(PRB1), PRB2=revCmp(PRB2))
      
      qc_tib <- ret_tib %>% filter(PRB0!=PRB1)
      qc_len <- qc_tib %>% base::nrow()
      if (qc_len != 0) {
        qc_tib %>% dplyr::select(1,PRB0,PRB1) %>% print()
        stop(glue::glue("{RET}[{funcTag}]: ERROR: pr={pr}, qc_len={qc_len} != 0!!!{RET}{RET}"))
      }
    } else {
      ret_tib <- ret_tib %>%
        tidyr::unite(PRB1, TAR,SEC,BOD, sep='', remove=FALSE) %>%
        tidyr::unite(PRB2, SEC,BOD,END, sep='', remove=FALSE) %>%
        dplyr::mutate(PRB1=revCmp(PRB1), PRB2=revCmp(PRB2))
    }
    
    # Add suffix to sequences for merging later
    ret_tib <- ret_tib %>%
      dplyr::select(PRB1,PRB2, NXB,CPN,TAR,BOD,END, everything()) %>%
      dplyr::rename(!!paste('PRB1',mu, sep=del):=PRB1,
                    !!paste('PRB2',mu, sep=del):=PRB2,
                    !!paste('NXB', mu, sep=del):=NXB,
                    !!paste('CPN', mu, sep=del):=CPN,
                    !!paste('TAR', mu, sep=del):=TAR,
                    !!paste('SEC', mu, sep=del):=SEC,
                    !!paste('BOD', mu, sep=del):=BOD,
                    !!paste('END', mu, sep=del):=END)
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ret_dat({ret_cnt})={RET}"))
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

des2prbsNOTES = function(srd, desSeq,
                         verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'des2prbs'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  # bscNewU <- bscNewU(desSeq)
  # bscNewM <- bscNewM(desSeq)
  # bscNewD <- bscNewD(desSeq)
  
  # my @desSetU = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewU);
  # my @desSetM = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewM);
  # my @desSetD = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewD);
  
  #    my @desSetU=unpack("A11"."A48"."AAAA"."A47"."A12",bscNewU($desSeq,$retUC));                                                                                                             
  #    my @desSetM=unpack("A11"."A48"."AAAA"."A47"."A12",bscNewM($desSeq,$retUC));                                                                                                             
  
  # my ($desNxbU, $desCpgU, $desSeqU, $desEndU);
  # my ($desNxbM, $desCpgM, $desSeqM, $desEndM);
  # my ($desNxbD, $desCpgD, $desSeqD, $desEndD);
  
  #         (  $desNxbU,    $desCpgU,               $desSeqU,                $desEndU) =                                                                                                       
  #   return( $desSetU[2], $desSetU[3], $desSetU[4].$desSetU[5].$desSetU[6], $desSetU[7],                                                                                                      
  #           $desSetM[2], $desSetM[3], $desSetM[4].$desSetM[5].$desSetM[6], $desSetM[7]) if ($desCO eq $C);                                                                                   
  #          ( $desNxbM,    $desCpgM,               $desSeqM,                $desEndM) =                                                                                                       
  
  #                  ( $desNxbU,            $desCpgU,                    $desSeqU,                  $desEndU) =                                                                                
  #   return( revCmpl($desSetU[4]), revCmpl($desSetU[3]), revCmpl($desSetU[1].$desSetU[2]), revCmpl($desSetU[0]),                                                                              
  #           revCmpl($desSetM[4]), revCmpl($desSetM[3]), revCmpl($desSetM[1].$desSetM[2]), revCmpl($desSetM[0])) if ($desCO eq $O);                                                           
  #                  ( $desNxbM,            $desCpgM,                    $desSeqM,                  $desEndM) =                                                                                
  
  
  # $$prbRef[$srd][$iU]=[ $desSetU[2], $desSetU[3], $desSetU[4].$desSetU[5].$desSetU[6], $desSetU[7] ];
  # $$prbRef[$srd][$iM]=[ $desSetM[2], $desSetM[3], $desSetM[4].$desSetM[5].$desSetM[6], $desSetM[7] ];
  # $$prbRef[$srd][$iD]=[ $desSetD[2], $desSetD[3], $desSetD[4].$desSetD[5].$desSetD[6], $desSetD[7] ];
  # 
  # $srd++;
  # $$prbRef[$srd][$iU]=[ revCmpl($desSetU[4]), revCmpl($desSetU[3]), revCmpl($desSetU[1].$desSetU[2]), revCmpl($desSetU[0]) ];
  # $$prbRef[$srd][$iM]=[ revCmpl($desSetM[4]), revCmpl($desSetM[3]), revCmpl($desSetM[1].$desSetM[2]), revCmpl($desSetM[0]) ];
  # $$prbRef[$srd][$iD]=[ revCmpl($desSetD[4]), revCmpl($desSetD[3]), revCmpl($desSetD[1].$desSetD[2]), revCmpl($desSetD[0]) ];
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  NULL
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Output improbe Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
writeImprobeInput = function(tib, name, dir, run=FALSE, 
                             exe=NULL, impTango=NULL, imp13mer=NULL, 
                             tbVar='BOTH', coVar='BOTH',
                             verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'writeImprobeInput'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}, dir={dir}.{RET}"))
  
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  
  imp_run_sh  <- file.path(dir, paste0(name,'.improbe-input.sh'))
  imp_inp_tsv <- file.path(dir, paste0(name,'.improbe-input.tsv'))
  imp_out_tsv <- file.path(dir, paste0(name,'.improbe-output.tsv'))
  imp_log_txt <- file.path(dir, paste0(name,'.improbe-output.log'))
  imp_fin_txt <- file.path(dir, paste0(name,'.improbe-output.fin.txt'))
  imp_cmd_str <- ''

  readr::write_tsv(tib, imp_inp_tsv)
  
  if (!is.null(exe) && ! is.null(impTango) && !is.null(imp13mer)) {
    # stopifnot(file.exists(exe))
    # stopifnot(file.exists(impTango))
    # stopifnot(file.exists(imp13mer))
    
    # ${CMD} -oASPE -tBOTH -cBoth -n${fn_13mer} -a${fn_tango} -V ${fn_in} >$fn_out 2>${fn_log}
    imp_cmd_str <- glue::glue("{exe} -oASPE -t{tbVar} -c{coVar} -n{imp13mer} -a{impTango} -V ",
                              "{imp_inp_tsv} >{imp_out_tsv} 2>{imp_log_txt}{RET}",
                              "touch {imp_fin_txt}{RET}")
    
    readr::write_file(imp_cmd_str, imp_run_sh)
    Sys.chmod(imp_run_sh, mode='0777')
    
    if (run) system(imp_run_sh, wait = TRUE)
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  imp_out_tsv
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Load/Format Input File Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
loadManifestRS = function(file, swap, revAllele=FALSE, 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(RS)={file}.{RET}"))
  
  snp_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(file) ))
  snp_swp_tib <- suppressMessages(suppressWarnings(readr::read_tsv(swap) ))
  
  # NOTE: Rediculous Swapping Needed for Infinium I M/U -> U/M from old manufacturing mistake::
  #  - This method should not be used...
  if (revAllele) {
    snp_man_tib <- snp_man_tib %>% 
      dplyr::mutate(
        TMP_ADD_A=AddressA_ID,
        TMP_ADD_B=AddressB_ID,
        TMP_SEQ_A=AlleleA_ProbeSeq,
        TMP_SEQ_B=AlleleB_ProbeSeq,
        IS_SWAP=case_when(
          Infinium_Design_Type=='I' & stringr::str_ends(AlleleA_ProbeSeq,'C') ~ TRUE,
          TRUE ~ FALSE
        ),
        AddressA_ID=case_when( IS_SWAP ~ TMP_ADD_B, TRUE ~ AddressA_ID),
        AddressB_ID=case_when( IS_SWAP ~ TMP_ADD_A, TRUE ~ AddressB_ID),
        AlleleA_ProbeSeq=case_when( IS_SWAP ~ TMP_SEQ_B, TRUE ~ AlleleA_ProbeSeq ),
        AlleleB_ProbeSeq=case_when( IS_SWAP ~ TMP_SEQ_A, TRUE ~ AlleleB_ProbeSeq )
      )
  }
  
  snp_man_tib <-  snp_man_tib %>% 
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$') ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq, 
                  Infinium_Design_Type, Next_Base, Color_Channel, CGN)
  
  snp_swp_tib <- snp_swp_tib %>% tidyr::separate(Seq_ID, into=c('IlmnID', 'diNUC'), sep='_') %>% 
    dplyr::rename(Forward_Sequence=Sequence, CHR=Chromosome, MAPINFO=Coordinate) %>% 
    dplyr::mutate(Forward_Sequence=stringr::str_replace(Forward_Sequence,'\\[CG\\]', paste0('[',diNUC,']') ) ) %>% 
    dplyr::select(-CpG_Island)
  
  snp_man_tib <- snp_man_tib %>% dplyr::inner_join(snp_swp_tib, by="IlmnID") %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  snp_man_tib
}                          

loadManifestCG = function(file, pr='cg', 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestCG'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CG)={file}.{RET}"))
  
  cpg_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(file, skip=7) )) %>%
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[([A-Za-z][A-Za-z])].*$', "\\$1"),'^\\\\') ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cpg_man_tib
}                          

loadManifestCH = function(file, pr='ch', ry='R', 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestCH'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CH)={file}.{RET}"))
  
  cph_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(cph_man_csv ))) %>% 
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(Strand=stringr::str_sub(IlmnID, -1), FR=Strand, CO='O',
                  TP=case_when( is.na(AlleleB_ProbeSeq) ~ 'II', TRUE ~ 'I' )) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=paste0(ry,stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[[A-Za-z]([A-Za-z])].*$', "\\$1"),'^\\\\')),
                  Forward_Sequence=stringr::str_replace(Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', paste0('[',diNUC,']')) ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cph_man_tib
}

loadManifestYH = function(file, pr='ch', 
                          verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadManifestCH'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CH)={file}.{RET}"))
  
  cph_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(cph_man_csv ))) %>% 
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(Strand=stringr::str_sub(IlmnID, -1), FR=Strand, CO='O',
                  TP=case_when( is.na(AlleleB_ProbeSeq) ~ 'II', TRUE ~ 'I' )) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=paste0('Y',stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[[A-Za-z]([A-Za-z])].*$', "\\$1"),'^\\\\')),
                  Forward_Sequence=stringr::str_replace(Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', paste0('[',diNUC,']')) ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cph_man_tib
}                          

#
# TBD:: This function should probably be deprecated::
#
loadImprobeDesign = function(file=NULL, src_des_tib=NULL, 
                             verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadImprobeDesign'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading improbe={file}.{RET}"))

  if (!is.null(file)) src_des_tib <- suppressMessages(suppressWarnings(readr::read_tsv(file)))
  stopifnot(!is.null(src_des_tib))
  
  src_des_tib <- src_des_tib %>%
    dplyr::rename(PRB1_U=UnMethyl_Probe_Sequence,
                  PRB1_M=Methyl_Probe_Sequence,
                  NXB_U=UnMethyl_Next_Base,
                  NXB_M=Methyl_Next_Base) %>%
    dplyr::mutate(Probe_ID=paste(Seq_ID, Methyl_Allele_FR_Strand, stringr::str_sub(Methyl_Allele_TB_Strand,1,1), 
                                 Methyl_Allele_CO_Strand, sep='_'),
                  FR=case_when(Methyl_Allele_FR_Strand=='F'   ~ TRUE, Methyl_Allele_FR_Strand=='R'   ~ FALSE, TRUE ~ NA),
                  TB=case_when(Methyl_Allele_TB_Strand=='TOP' ~ TRUE, Methyl_Allele_TB_Strand=='BOT' ~ FALSE, TRUE ~ NA),
                  CO=case_when(Methyl_Allele_CO_Strand=='C'   ~ TRUE, Methyl_Allele_CO_Strand=='O'   ~ FALSE, TRUE ~ NA),
                  diNUC='CG',
                  NXB_IMP=case_when(NXB_U==NXB_M ~ NXB_U,
                                    TRUE ~ NA_character_),
                  COL_U=case_when(NXB_U=='A'|NXB_U=='T' ~ FALSE, # Red == FALSE
                                  NXB_U=='C'|NXB_U=='G' ~ TRUE,  # Grn == TRUE
                                  TRUE ~ NA),
                  COL_M=case_when(NXB_M=='A'|NXB_M=='T' ~ FALSE, # Red == FALSE
                                  NXB_M=='C'|NXB_M=='G' ~ TRUE,  # Grn == TRUE
                                  TRUE ~ NA),
                  COL_IMP=case_when(COL_U==COL_M ~ COL_U,
                                    TRUE ~ NA)
                  
    ) %>%
    
    # Design Score Parameters::
    dplyr::rename(
      PRB_SCR_U=UnMethyl_Final_Score,
      PRB_SCR_M=Methyl_Final_Score,
      PRB_SCR_S=Probeset_Score,
      
      TM_RAW_M=Methyl_Tm,
      TM_SCR_M=Methyl_Tm_Score,
      GC_RAW_M=Methyl_GC_Percent,
      GC_SCR_M=Methyl_GC_Score,
      KM_RAW_M=Methyl_13mer_Count,
      KM_SCR_M=Methyl_13mer_Score,
      AD_RAW_M=Methyl_Address_Count,
      AD_SCR_M=Methyl_Address_Score,
      CM_RAW_M=Methyl_Self_Complementarity,
      CM_SCR_M=Methyl_Self_Complementarity_Score,
      MO_RAW_M=Methyl_Mono_Run,
      MO_SCR_M=Methyl_Mono_Run_Score,
      EP_RAW_M=Methyl_Ectopic_Count,
      EP_SCR_M=Methyl_Ectopic_Score,
      CG_RAW_M=Methyl_Underlying_CpG_Count,
      MD_RAW_M=Methyl_Underlying_CpG_Min_Dist,
      CG_SCR_M=Methyl_Underlying_CpG_Score,
      NB_SCR_M=Methyl_Next_Base_Score,
      
      TM_RAW_U=UnMethyl_Tm,
      TM_SCR_U=UnMethyl_Tm_Score,
      GC_RAW_U=UnMethyl_GC_Percent,
      GC_SCR_U=UnMethyl_GC_Score,
      KM_RAW_U=UnMethyl_13mer_Count,
      KM_SCR_U=UnMethyl_13mer_Score,
      AD_RAW_U=UnMethyl_Address_Count,
      AD_SCR_U=UnMethyl_Address_Score,
      CM_RAW_U=UnMethyl_Self_Complementarity,
      CM_SCR_U=UnMethyl_Self_Complementarity_Score,
      MO_RAW_U=UnMethyl_Mono_Run,
      MO_SCR_U=UnMethyl_Mono_Run_Score,
      EP_RAW_U=UnMethyl_Ectopic_Count,
      EP_SCR_U=UnMethyl_Ectopic_Score,
      CG_RAW_U=UnMethyl_Underlying_CpG_Count,
      MD_RAW_U=UnMethyl_Underlying_CpG_Min_Dist,
      CG_SCR_U=UnMethyl_Underlying_CpG_Score,
      NB_SCR_U=UnMethyl_Next_Base_Score
    ) %>%
    dplyr::select(Probe_ID, Seq_ID, FR, TB, CO, diNUC, NXB_IMP, COL_IMP, PRB1_U, PRB1_M, 
                  Forward_Sequence, 
                  Genome_Build, Chromosome, Coordinate,
                  dplyr::contains("_RAW_"), dplyr::contains("_SCR_"))
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  src_des_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Format Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
shearBrac = function(x) {
  x %>% stringr::str_remove('\\[') %>% stringr::str_remove('\\]')
}

addBrac = function(seq) {
  seq_len <- stringr::str_length(seq)
  mid_len <- 2
  snp_idx <- as.integer(seq_len / 2)
  
  pre_beg <- 1
  pre_end <- snp_idx-1
  
  mid_beg <- snp_idx
  mid_end <- mid_beg+mid_len-1
  
  pos_beg <- mid_end+1
  pos_end <- seq_len
  
  pre_seq <- seq %>% stringr::str_sub(pre_beg,pre_end)
  mid_seq <- seq %>% stringr::str_sub(mid_beg,mid_end)
  pos_seq <- seq %>% stringr::str_sub(pos_beg,pos_end)
  
  new_seq <- paste0(pre_seq,'[',mid_seq,']',pos_seq)
  
  new_seq
}

addDesignSeqCG = function(tib, seq, add, din=NULL,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'addDesignSeqCG'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    add_sym <- add %>% rlang::sym()
    seq_sym <- seq %>% rlang::sym()
    
    ret_tib <- tib %>% 
      dplyr::mutate(!!add:=stringr::str_replace(!!seq_sym, "\\[[a-zA-Z][a-zA-Z]\\]", "[CG]"))
    
    if (!is.null(din))
      ret_tib <- ret_tib %>% dplyr::mutate(
        !!din := !!seq_sym %>%
          stringr::str_remove("\\].*$") %>%
          stringr::str_remove("^.*\\[") )

    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

replaceDesignSeqCG = function(tib, seq, add, nuc,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'replaceDesignSeqCG'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    add_sym <- add %>% rlang::sym()
    seq_sym <- seq %>% rlang::sym()
    nuc_sym <- nuc %>% rlang::sym()
    
    ret_tib <- tib %>% 
      dplyr::mutate(PreSeq=stringr::str_remove(!!seq_sym, "\\[.*$"),
                    PosSeq=stringr::str_remove(!!seq_sym, "^.*\\]"),
                    !!add:=paste0(PreSeq,"[",!!nuc_sym,"]",PosSeq) ) %>%
      dplyr::select(-PreSeq, -PosSeq)

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
#                     Illumina Strand Methods:: TOP/BOT
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

setTopBot_tib = function(tib, seqKey, srdKey, topKey=NULL, max=0,
                         verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'setTopBot_tib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; select field={seqKey}, new field={srdKey}...{RET}"))
  
  seqKey <- rlang::sym(seqKey)
  srdKey <- rlang::sym(srdKey)
  if (!is.null(topKey)) topKey <- rlang::sym(topKey)
  
  # TBD:: Should remove any seqKey that is na or not the correct length...
  
  ret_tib <- NULL
  stime <- system.time({
    if (max!=0) tib <- tib %>% head(n=max)
    
    bit_tib <- tib %>% 
      dplyr::mutate(
        pre_seq = !!seqKey %>% stringr::str_remove('\\[.*$') %>% 
          Biostrings::reverse() %>%
          stringr::str_to_upper() %>% 
          stringr::str_replace_all('A', '1') %>% 
          stringr::str_replace_all('T', '1') %>% 
          stringr::str_replace_all('C', '2') %>% 
          stringr::str_replace_all('G', '2'),
        
        pos_seq = !!seqKey %>% stringr::str_remove('^.*\\]') %>%
          stringr::str_to_upper() %>% 
          stringr::str_replace_all('A', '1') %>% 
          stringr::str_replace_all('T', '1') %>% 
          stringr::str_replace_all('C', '2') %>% 
          stringr::str_replace_all('G', '2')
      ) %>% dplyr::select(pre_seq, pos_seq)
    
    pre_bit <- bit_tib$pre_seq %>% stringr::str_split('', simplify = TRUE)
    pos_bit <- bit_tib$pos_seq %>% stringr::str_split('', simplify = TRUE)
    
    pre_mat <- pre_bit %>% as.data.frame() %>% data.matrix()
    pos_mat <- pos_bit %>% as.data.frame() %>% data.matrix()
    dif_mat <- pre_mat-pos_mat
    dif_mat <- dif_mat %>% cbind(rep(2, dim(dif_mat)[1]))
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} dif_mat={RET}"))
    if (verbose>=vt+4) print(dif_mat)
    
    # This is just for a sanity check...
    #  dif_mat[1,1:60] = 0
    
    bit_vec <- apply(dif_mat,1, function(x) head(x[x!=0],1)) + 1
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} bit_vec={RET}"))
    if (verbose>=vt+4) print(bit_vec)
    
    # Convert numeric values back to characters (T=TOP,B=BOT,U=Unknown)
    tb_vec <- bit_vec %>% as.character() %>% 
      stringr::str_replace('0', 'T') %>% 
      stringr::str_replace('2', 'B') %>%
      stringr::str_replace('3', 'U')
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} tb_vec={RET}"))
    if (verbose>=vt+4) print(tb_vec)
    
    ret_tib <- tib %>% dplyr::mutate(!!srdKey := tb_vec)
    
    if (!is.null(topKey)) {
      revCmp_vec <- ret_tib %>% 
        dplyr::pull(!!seqKey) %>% shearBrac() %>% revCmp() %>% addBrac()
      ret_tib <- ret_tib %>% dplyr::mutate(TMP_REVCMP_SEQ=revCmp_vec)
      
      ret_tib <- ret_tib %>% dplyr::mutate(
        !!topKey:=dplyr::case_when(
          !!srdKey=='T' ~ !!seqKey,
          !!srdKey=='B' ~ TMP_REVCMP_SEQ,
          TRUE ~ NA_character_)
      ) %>% dplyr::select(-TMP_REVCMP_SEQ)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# Single variable method that is OUT-OF-DATE; Should use the method above...
isTopBot_single = function(seq, verbose=0, vt=4) {
  seq <- seq %>% stringr::str_to_upper() %>% 
    stringr::str_replace_all('T', 'A') %>% 
    stringr::str_replace_all('G', 'C')

  seq_vec <- stringr::str_split(seq,'\\[', simplify=TRUE) %>% as.vector()
  pre_seq <- seq_vec[1] %>% Biostrings::reverse()
  seq_vec <- stringr::str_split(seq,'\\]', simplify=TRUE) %>% as.vector()
  pos_seq <- seq_vec[2]
  
  pre_vec <- pre_seq %>% stringr::str_split('', simplify=TRUE) %>% as.vector()
  pos_vec <- pos_seq %>% stringr::str_split('', simplify=TRUE) %>% as.vector()
  
  if (verbose>=vt) {
    print(seq)
    print(pre_seq)
    print(pos_seq)
  }
  
  min_len <- min(stringr::str_length(pre_seq),stringr::str_length(pos_seq))
  for (ii in c(1:min_len)) {
    if (pre_vec[ii]==pos_vec[ii]) next
    if (pre_vec[ii]!=pos_vec[ii]) {
      if (pre_vec[ii]=='A') return('T')
      return('B')
    }
  }
  return('U')
}

isTopBots = function(x, verbose=0, vt=4) {
  x <- lapply(x,isTopBot_single)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Basic Bisulfite Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

MAPDi = function(x) {
  if (length(MAP_DI[[x]])==0) return(NA)
  MAP_DI[[x]]
}
mapDIs = function(x) {
  x <- lapply(x, MAPDi) %>% BiocGenerics::unlist()
}


bscU = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  if (uc) x <- tr(x, 'CYSMBHV', 'TTKWKWD')
  else    x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  x
}
bscUs = function(x, uc=FALSE) { bscU(x) }

MAPM = function(x) {
  if (length(MAP_M[[x]])==0) return(x)
  MAP_M[[x]]
}
bscM = function(x) { stringr::str_replace_all(x, '([CYSMBHV][GRSKBDV])', MAPM) }
bscMs = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  x <- lapply(x, bscM) %>% BiocGenerics::unlist()
  x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  if (uc) x <- stringr::str_to_upper(x)
  x
}

MAPD = function(x) {
  if (length(MAP_D[[x]])==0) return(x)
  MAP_D[[x]]
}
bscD = function(x) { stringr::str_replace_all(x, '([CYSMBHV][GRSKBDV])', MAPD) }
bscDs = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  x <- lapply(x, bscD) %>% BiocGenerics::unlist()
  x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  if (uc) x <- stringr::str_to_upper(x)
  x
}

QMAP = function(x, mu) {
  if (mu=='U') {
    return(QMAP_U[[x]])
  } else if (mu=='M') {
    return(QMAP_M[[x]])
  }
  x
}
qmaps = function(x, mu) {
  x <- lapply(x, QMAP, mu) %>% BiocGenerics::unlist()
}

cmpIUPAC = function(x) {
  if (base::is.element(x, names(IUPAC_EQ) )) return(IUPAC_EQ[[x]])
  # if (is.null(IUPAC_EQ[[x]])) return(FALSE)
  # if (length(IUPAC_EQ[[x]])==0) return(FALSE)
  FALSE
}

cmpIUPACs = function(x) {
  x <- lapply(x, cmpIUPAC) %>% BiocGenerics::unlist()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      MisMatch Probe Comparison Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cmpInfIMU_MisMatch = function(tib, fieldA, fieldB, mu, del='_',
                              verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'cmpInfIMU_MisMatch'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldA={fieldA}, fieldB={fieldB} mu={mu}{RET}"))
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(
      BOD_NumMM=mapply(
        adist,
        stringr::str_sub(stringr::str_to_upper(!!fieldA),1,stringr::str_length(!!fieldA)-1),
        stringr::str_sub(stringr::str_to_upper(!!fieldB),1,stringr::str_length(!!fieldB)-1) ),
      
      DI_NUC_AB=paste0(
        stringr::str_to_upper(stringr::str_sub(!!fieldA,stringr::str_length(!!fieldA),stringr::str_length(!!fieldA)) ),
        stringr::str_to_upper(stringr::str_sub(!!fieldB,stringr::str_length(!!fieldB),stringr::str_length(!!fieldB)) )
      ),
      TAR_EQU=cmpIUPACs(DI_NUC_AB)
    ) %>%
    dplyr::rename(!!paste('BOD_NumMM',mu, sep=del):=BOD_NumMM,
                  !!paste('TAR_EQU',  mu, sep=del):=TAR_EQU)
  
  tib
}

cmpInfI_MisMatch = function(tib, fieldAU, fieldBU, fieldAM, fieldBM, del='_',
                            verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'cmpInfI_MisMatch'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldAU={fieldAU}, fieldBU={fieldBU}{RET}"))
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldAM={fieldAM}, fieldBM={fieldBM}{RET}"))

  tib <- tib %>% cmpInfIMU_MisMatch(fieldAU, fieldBU, mu='U', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>% cmpInfIMU_MisMatch(fieldAM, fieldBM, mu='M', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>% dplyr::mutate(
    Man_MisMatch=(BOD_NumMM_U+BOD_NumMM_M)/2, #, na.rm=TRUE),
    Man_TarMatch=case_when(TAR_EQU_U & TAR_EQU_M ~ TRUE, TRUE ~ FALSE) ) %>%
    dplyr::select(-c(BOD_NumMM_U,BOD_NumMM_M,TAR_EQU_U,TAR_EQU_M))
  
  tib
}

cmpInfII_MisMatch = function(tib, fieldA, fieldB, mu='D', del='_',
                             verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'cmpInfI'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% cmpInfIMU_MisMatch(fieldA, fieldB, mu='D', del=del,verbose=verbose, vt=vt+1,tc=tc+1) %>%
    dplyr::rename(
      Man_MisMatch=BOD_NumMM_D,
      Man_TarMatch=TAR_EQU_D)
    # dplyr::select(-c(BOD_NumMM_U,BOD_NumMM_M,TAR_EQU_U,TAR_EQU_M))

  tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Exact Probe Comparison Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cmpInfIMU= function(tib, fieldA, fieldB, mu, del='_',
                    verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'cmpInfMU'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldA={fieldA}, fieldB={fieldB} mu={mu}{RET}"))
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(
      SUB_SEQ_A=stringr::str_to_upper(stringr::str_sub(!!fieldA,1,stringr::str_length(!!fieldA)-1)),
      SUB_SEQ_B=stringr::str_to_upper(stringr::str_sub(!!fieldB,1,stringr::str_length(!!fieldB)-1)),
      DI_NUC_AB=paste0(
        stringr::str_to_upper(stringr::str_sub(!!fieldA,stringr::str_length(!!fieldA),stringr::str_length(!!fieldA)) ),
        stringr::str_to_upper(stringr::str_sub(!!fieldB,stringr::str_length(!!fieldB),stringr::str_length(!!fieldB)) )
      ),
      
      BOD_EQU=case_when(SUB_SEQ_A==SUB_SEQ_B ~ TRUE, TRUE ~ FALSE),
      TAR_EQU=cmpIUPACs(DI_NUC_AB),
      Inf1_Match=case_when(BOD_EQU & BOD_EQU==TAR_EQU ~ TRUE, TRUE ~ FALSE)
    ) %>%
    dplyr::select(-c(SUB_SEQ_A,SUB_SEQ_B,BOD_EQU,DI_NUC_AB,TAR_EQU)) %>%
    dplyr::rename(!!paste('Inf1_Match',mu, sep=del):=Inf1_Match)
  
  if (verbose>=vt+1) print(tib)
  
  tib
}

cmpInfI = function(tib, fieldAU, fieldBU, fieldAM, fieldBM, del='_',
                   verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'cmpInfI'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldAU={fieldAU}, fieldBU={fieldBU}{RET}"))
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldAM={fieldAM}, fieldBM={fieldBM}{RET}"))
  
  tib <- tib %>% cmpInfIMU(fieldAU, fieldBU, mu='U', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>% cmpInfIMU(fieldAM, fieldBM, mu='M', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>%
    dplyr::mutate(Man_Match=case_when(Inf1_Match_U & Inf1_Match_M ~ TRUE, TRUE ~ FALSE) ) %>%
    dplyr::select(-c(Inf1_Match_U,Inf1_Match_M))
  
  tib
}

cmpInfII = function(tib, fieldA, fieldB, mu='D', del='_',
                    verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'cmpInfI'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(Man_Match=stringr::str_to_upper(!!fieldA)==stringr::str_to_upper(!!fieldB) )
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Basic Reverse/Complement Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
revCmp = function(x) {
  Biostrings::reverse(x) %>% cmpl()
}

revSeq = function(x) {
  Biostrings::reverse(x)
}

cmpl = function(x) {
  tr(x, 'ACTGRYSWKMBDHVactgryswkmbdhv[]','TGACYRSWMKVHDBtgacyrswmkvhdb][')
  # x <- tr(x, 'ACTGRYSWKMBDHV','TGACYRSWMKVHDB')
  # x <- tr(x, 'actgryswkmbdhv','tgacyrswmkvhdb')
  # x
}

tr = function(x, old, new) {
  Biostrings::chartr(old, new, x)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Improbe Header Mapping Structure::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

INIT_IMP_HEADER = function() {
  IMP_TYPES = cols(
    Seq_ID = col_character(),
    Forward_Sequence = col_character(),
    Genome_Build = col_character(),
    Chromosome = col_character(),
    Coordinate = col_double(),
    Design_State = col_character(),
    Seq_Length = col_double(),
    Forward_CpG_Coord = col_double(),
    TB_Strand = col_character(),
    Top_Sequence = col_character(),
    Top_CpG_Coord = col_double(),
    Probe_Type = col_character(),
    Probeset_ID = col_character(),
    Probeset_Score = col_double(),
    Methyl_Probe_ID = col_character(),
    Methyl_Probe_Sequence = col_character(),
    Methyl_Probe_Length = col_double(),
    Methyl_Start_Coord = col_double(),
    Methyl_End_Coord = col_double(),
    Methyl_Probe_Covered_Top_Sequence = col_character(),
    Methyl_Allele_FR_Strand = col_character(),
    Methyl_Allele_TB_Strand = col_character(),
    Methyl_Allele_CO_Strand = col_character(),
    Methyl_Allele_Type = col_character(),
    Methyl_Final_Score = col_double(),
    Methyl_Tm = col_double(),
    Methyl_Tm_Score = col_double(),
    Methyl_GC_Percent = col_double(),
    Methyl_GC_Score = col_double(),
    Methyl_13mer_Count = col_double(),
    Methyl_13mer_Score = col_double(),
    Methyl_Address_Count = col_double(),
    Methyl_Address_Score = col_double(),
    Methyl_Self_Complementarity = col_double(),
    Methyl_Self_Complementarity_Score = col_double(),
    Methyl_Mono_Run = col_double(),
    Methyl_Mono_Run_Score = col_double(),
    Methyl_Ectopic_Count = col_double(),
    Methyl_Ectopic_Score = col_double(),
    Methyl_Underlying_CpG_Count = col_double(),
    Methyl_Underlying_CpG_Min_Dist = col_double(),
    Methyl_Underlying_CpG_Score = col_double(),
    Methyl_In_CpG_Island_Relaxed = col_logical(),
    Methyl_CpG_Island_Score = col_double(),
    Methyl_Next_Base = col_character(),
    Methyl_Next_Base_Score = col_double(),
    UnMethyl_Probe_ID = col_character(),
    UnMethyl_Probe_Sequence = col_character(),
    UnMethyl_Probe_Length = col_double(),
    UnMethyl_Start_Coord = col_double(),
    UnMethyl_End_Coord = col_double(),
    UnMethyl_Probe_Covered_Top_Sequence = col_character(),
    UnMethyl_Allele_FR_Strand = col_character(),
    UnMethyl_Allele_TB_Strand = col_character(),
    UnMethyl_Allele_CO_Strand = col_character(),
    UnMethyl_Allele_Type = col_character(),
    UnMethyl_Final_Score = col_double(),
    UnMethyl_Tm = col_double(),
    UnMethyl_Tm_Score = col_double(),
    UnMethyl_GC_Percent = col_double(),
    UnMethyl_GC_Score = col_double(),
    UnMethyl_13mer_Count = col_double(),
    UnMethyl_13mer_Score = col_double(),
    UnMethyl_Address_Count = col_double(),
    UnMethyl_Address_Score = col_double(),
    UnMethyl_Self_Complementarity = col_double(),
    UnMethyl_Self_Complementarity_Score = col_double(),
    UnMethyl_Mono_Run = col_double(),
    UnMethyl_Mono_Run_Score = col_double(),
    UnMethyl_Ectopic_Count = col_double(),
    UnMethyl_Ectopic_Score = col_double(),
    UnMethyl_Underlying_CpG_Count = col_double(),
    UnMethyl_Underlying_CpG_Min_Dist = col_double(),
    UnMethyl_Underlying_CpG_Score = col_double(),
    UnMethyl_In_CpG_Island_Relaxed = col_logical(),
    UnMethyl_CpG_Island_Score = col_double(),
    UnMethyl_Next_Base = col_character(),
    UnMethyl_Next_Base_Score = col_double()
  )
  
  IMP_TYPES
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Nucelotide Mapping Structures::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# mapA = function(x) {
#   if (length(MAP_A[[x]])==0) return(x)
#   MAP_A[[x]]
# }
# INIT_MAP_A = function() {
#   MAP_A <- NULL
#   MAP_A[['AA']] <- 'aa'
#   MAP_A[['Aa']] <- 'aa'
#   MAP_A[['aA']] <- 'aa'
#   MAP_A[['aa']] <- 'aa'
#   
#   MAP_A
# }
# MAP_A <- INIT_MAP_A()

# Di-Nuc to IUPAC::
INIT_MAP_DI = function() {
  MAP <- NULL
  MAP[['AC']] <- 'M'
  MAP[['AG']] <- 'R'
  MAP[['AT']] <- 'W'
  MAP[['AA']] <- 'A'
  
  MAP[['CA']] <- 'M'
  MAP[['CT']] <- 'Y'
  MAP[['CG']] <- 'S'
  MAP[['CC']] <- 'C'
  
  MAP[['GA']] <- 'R'
  MAP[['GT']] <- 'K'
  MAP[['GC']] <- 'S'
  MAP[['GG']] <- 'G'
  
  MAP[['TC']] <- 'Y'
  MAP[['TG']] <- 'K'
  MAP[['TA']] <- 'W'
  MAP[['TT']] <- 'T'
  
  MAP
}
MAP_DI <- INIT_MAP_DI()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Initialize M/D Maps for Bisulfite Conversion::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Data generated with:: /Users/bbarnes/Documents/Projects/scripts/mapMD.pl
INIT_MAP_M = function() {
  MAP <- NULL
  MAP[['CG']] <- 'cG'
  MAP[['CR']] <- 'cR'
  MAP[['CS']] <- 'cS'
  MAP[['CK']] <- 'cK'
  MAP[['CB']] <- 'cB'
  MAP[['CD']] <- 'cD'
  MAP[['CV']] <- 'cV'
  MAP[['YG']] <- 'yG'
  MAP[['YR']] <- 'yR'
  MAP[['YS']] <- 'yS'
  MAP[['YK']] <- 'yK'
  MAP[['YB']] <- 'yB'
  MAP[['YD']] <- 'yD'
  MAP[['YV']] <- 'yV'
  MAP[['SG']] <- 'sG'
  MAP[['SR']] <- 'sR'
  MAP[['SS']] <- 'sS'
  MAP[['SK']] <- 'sK'
  MAP[['SB']] <- 'sB'
  MAP[['SD']] <- 'sD'
  MAP[['SV']] <- 'sV'
  MAP[['MG']] <- 'mG'
  MAP[['MR']] <- 'mR'
  MAP[['MS']] <- 'mS'
  MAP[['MK']] <- 'mK'
  MAP[['MB']] <- 'mB'
  MAP[['MD']] <- 'mD'
  MAP[['MV']] <- 'mV'
  MAP[['BG']] <- 'bG'
  MAP[['BR']] <- 'bR'
  MAP[['BS']] <- 'bS'
  MAP[['BK']] <- 'bK'
  MAP[['BB']] <- 'bB'
  MAP[['BD']] <- 'bD'
  MAP[['BV']] <- 'bV'
  MAP[['HG']] <- 'hG'
  MAP[['HR']] <- 'hR'
  MAP[['HS']] <- 'hS'
  MAP[['HK']] <- 'hK'
  MAP[['HB']] <- 'hB'
  MAP[['HD']] <- 'hD'
  MAP[['HV']] <- 'hV'
  MAP[['VG']] <- 'vG'
  MAP[['VR']] <- 'vR'
  MAP[['VS']] <- 'vS'
  MAP[['VK']] <- 'vK'
  MAP[['VB']] <- 'vB'
  MAP[['VD']] <- 'vD'
  MAP[['VV']] <- 'vV'
  
  MAP
}
MAP_M <- INIT_MAP_M()

INIT_MAP_D = function() {
  MAP <- NULL
  
  MAP[['CG']] <- 'yG'
  MAP[['CR']] <- 'yR'
  MAP[['CS']] <- 'yS'
  MAP[['CK']] <- 'yK'
  MAP[['CB']] <- 'yB'
  MAP[['CD']] <- 'yD'
  MAP[['CV']] <- 'yV'
  MAP[['YG']] <- 'yG'
  MAP[['YR']] <- 'yR'
  MAP[['YS']] <- 'yS'
  MAP[['YK']] <- 'yK'
  MAP[['YB']] <- 'yB'
  MAP[['YD']] <- 'yD'
  MAP[['YV']] <- 'yV'
  MAP[['SG']] <- 'bG'
  MAP[['SR']] <- 'bR'
  MAP[['SS']] <- 'bS'
  MAP[['SK']] <- 'bK'
  MAP[['SB']] <- 'bB'
  MAP[['SD']] <- 'bD'
  MAP[['SV']] <- 'bV'
  MAP[['MG']] <- 'hG'
  MAP[['MR']] <- 'hR'
  MAP[['MS']] <- 'hS'
  MAP[['MK']] <- 'hK'
  MAP[['MB']] <- 'hB'
  MAP[['MD']] <- 'hD'
  MAP[['MV']] <- 'hV'
  MAP[['BG']] <- 'bG'
  MAP[['BR']] <- 'bR'
  MAP[['BS']] <- 'bS'
  MAP[['BK']] <- 'bK'
  MAP[['BB']] <- 'bB'
  MAP[['BD']] <- 'bD'
  MAP[['BV']] <- 'bV'
  MAP[['HG']] <- 'hG'
  MAP[['HR']] <- 'hR'
  MAP[['HS']] <- 'hS'
  MAP[['HK']] <- 'hK'
  MAP[['HB']] <- 'hB'
  MAP[['HD']] <- 'hD'
  MAP[['HV']] <- 'hV'
  MAP[['VG']] <- 'nG'
  MAP[['VR']] <- 'nR'
  MAP[['VS']] <- 'nS'
  MAP[['VK']] <- 'nK'
  MAP[['VB']] <- 'nB'
  MAP[['VD']] <- 'nD'
  MAP[['VV']] <- 'nV'
  
  MAP
}
MAP_D <- INIT_MAP_D()



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Initialize Query Maps::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# These are for mapping Query positions to their possible outcomes based on
#  the assumption of unmethylated and methylated
#
# OLD_INIT_QMAP_* are previous versions that did not account for Infinium I
#  probes with some degenerate bases like W={A,T} that would still work
#

# TBD:: Address the tri-nucelotide version!!!

# Unmethylated will always use the letter earliest in the alphabet
#  e.g. 
#   A over T, 
#   A over C, 
#   A over G, 
#   C over G, 
#   C over T, 
#   G over T,
#
INIT_QMAP_U = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'
  MAP[['C']] <- 'C'  # '-'
  MAP[['G']] <- 'G'  # '-'
  MAP[['T']] <- 'T'
  
  MAP[['R']] <- 'A'
  MAP[['Y']] <- 'T'
  MAP[['S']] <- 'C'  # C over G 
  MAP[['W']] <- 'A'  # A over T
  # MAP[['S']] <- 'S'  # '-' Old Version
  # MAP[['W']] <- 'W'  # '-' Old Version
  
  MAP[['K']] <- 'T'
  MAP[['M']] <- 'A'
  MAP[['B']] <- 'T'
  MAP[['D']] <- 'W'
  MAP[['H']] <- 'W'
  MAP[['V']] <- 'A'
  
  MAP[['N']] <- 'W'
  
  # Lower Case::
  MAP[['a']] <- 'a'
  MAP[['c']] <- 'c'  # '-'
  MAP[['g']] <- 'g'  # '-'
  MAP[['t']] <- 't'
  
  MAP[['r']] <- 'a'
  MAP[['y']] <- 't'
  MAP[['s']] <- 'c'  # c over g
  MAP[['w']] <- 'a'  # a over t
  # MAP[['s']] <- 's'  # '-' Old Version
  # MAP[['w']] <- 'w'  # '-' Old Version
  
  MAP[['k']] <- 't'
  MAP[['m']] <- 'a'
  MAP[['b']] <- 't'
  MAP[['d']] <- 'w'
  MAP[['h']] <- 'w'
  MAP[['v']] <- 'a'
  
  MAP[['n']] <- 'w'
  
  MAP
}
QMAP_U <- INIT_QMAP_U()

# Methylated will always use the letter latest in the alphabet
#  e.g. 
#   T over A, 
#   T over C, 
#   T over G,
#   G over A, 
#   G over C, 
#   C over A, 
#
INIT_QMAP_M = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'  # '-'
  MAP[['C']] <- 'C'
  MAP[['G']] <- 'G'
  MAP[['T']] <- 'T'  # '-'
  
  MAP[['R']] <- 'G'
  MAP[['Y']] <- 'C'
  MAP[['S']] <- 'G'  # G over C
  MAP[['W']] <- 'T'  # T over A
  # MAP[['S']] <- 'S'  # '-' Old Version
  # MAP[['W']] <- 'W'  # '-' Old Version
  
  MAP[['K']] <- 'G'
  MAP[['M']] <- 'C'
  MAP[['B']] <- 'S'
  MAP[['D']] <- 'G'
  MAP[['H']] <- 'C'
  MAP[['V']] <- 'S'
  
  MAP[['N']] <- 'S'
  
  # Lower Case::
  MAP[['a']] <- 'a'  # '-'
  MAP[['c']] <- 'c'
  MAP[['g']] <- 'g'
  MAP[['t']] <- 't'  # '-'
  
  MAP[['r']] <- 'g'
  MAP[['y']] <- 'c'
  MAP[['s']] <- 'g'  # g over c
  MAP[['w']] <- 't'  # t over a
  # MAP[['s']] <- 's'  # '-' Old Version
  # MAP[['w']] <- 'w'  # '-' Old Version

  MAP[['k']] <- 'g'
  MAP[['m']] <- 'c'
  MAP[['b']] <- 's'
  MAP[['d']] <- 'g'
  MAP[['h']] <- 'c'
  MAP[['v']] <- 's'
  
  MAP[['n']] <- 's'
  
  MAP
}
QMAP_M <- INIT_QMAP_M()

# OLD Versions to be used now::
#
INIT_QMAP_U = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'
  MAP[['C']] <- 'C'  # '-'
  MAP[['G']] <- 'G'  # '-'
  MAP[['T']] <- 'T'
  
  MAP[['R']] <- 'A'
  MAP[['Y']] <- 'T'
  MAP[['S']] <- 'S'  # '-'
  MAP[['W']] <- 'W'  # '-'
  
  MAP[['K']] <- 'T'
  MAP[['M']] <- 'A'
  MAP[['B']] <- 'T'
  MAP[['D']] <- 'W'
  MAP[['H']] <- 'W'
  MAP[['V']] <- 'A'
  
  MAP[['N']] <- 'W'
  
  # Lower Case::
  MAP[['a']] <- 'a'
  MAP[['c']] <- 'c'  # '-'
  MAP[['g']] <- 'g'  # '-'
  MAP[['t']] <- 't'
  
  MAP[['r']] <- 'a'
  MAP[['y']] <- 't'
  MAP[['s']] <- 's'  # '-'
  MAP[['w']] <- 'w'  # '-'
  
  MAP[['k']] <- 't'
  MAP[['m']] <- 'a'
  MAP[['b']] <- 't'
  MAP[['d']] <- 'w'
  MAP[['h']] <- 'w'
  MAP[['v']] <- 'a'
  
  MAP[['n']] <- 'w'
  
  MAP
}
# OLD_QMAP_U <- OLD_INIT_QMAP_U()

OLD_INIT_QMAP_M = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'  # '-'
  MAP[['C']] <- 'C'
  MAP[['G']] <- 'G'
  MAP[['T']] <- 'T'  # '-'
  
  MAP[['R']] <- 'G'
  MAP[['Y']] <- 'C'
  MAP[['S']] <- 'S'  # '-'
  MAP[['W']] <- 'W'  # '-'
  
  MAP[['K']] <- 'G'
  MAP[['M']] <- 'C'
  MAP[['B']] <- 'S'
  MAP[['D']] <- 'G'
  MAP[['H']] <- 'C'
  MAP[['V']] <- 'S'
  
  MAP[['N']] <- 'S'
  
  # Lower Case::
  MAP[['a']] <- 'a'  # '-'
  MAP[['c']] <- 'c'
  MAP[['g']] <- 'g'
  MAP[['t']] <- 't'  # '-'
  
  MAP[['r']] <- 'g'
  MAP[['y']] <- 'c'
  MAP[['s']] <- 's'  # '-'
  MAP[['w']] <- 'w'  # '-'
  
  MAP[['k']] <- 'g'
  MAP[['m']] <- 'c'
  MAP[['b']] <- 's'
  MAP[['d']] <- 'g'
  MAP[['h']] <- 'c'
  MAP[['v']] <- 's'
  
  MAP[['n']] <- 's'
  
  MAP
}
# QMAP_M <- OLD_INIT_QMAP_M()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Initialize IUPAC Equivelency Tables::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# This is really for convience...
#
INIT_IUPAC_EQ = function() {
  MAP <- NULL
  
  MAP[['AA']] <- TRUE
  MAP[['AD']] <- TRUE
  MAP[['AH']] <- TRUE
  MAP[['AM']] <- TRUE
  MAP[['AN']] <- TRUE
  MAP[['AR']] <- TRUE
  MAP[['AV']] <- TRUE
  MAP[['AW']] <- TRUE
  MAP[['BB']] <- TRUE
  MAP[['BC']] <- TRUE
  MAP[['BG']] <- TRUE
  MAP[['BK']] <- TRUE
  MAP[['BN']] <- TRUE
  MAP[['BS']] <- TRUE
  MAP[['BT']] <- TRUE
  MAP[['BY']] <- TRUE
  MAP[['CB']] <- TRUE
  MAP[['CC']] <- TRUE
  MAP[['CH']] <- TRUE
  MAP[['CM']] <- TRUE
  MAP[['CN']] <- TRUE
  MAP[['CS']] <- TRUE
  MAP[['CV']] <- TRUE
  MAP[['CY']] <- TRUE
  MAP[['DA']] <- TRUE
  MAP[['DD']] <- TRUE
  MAP[['DG']] <- TRUE
  MAP[['DK']] <- TRUE
  MAP[['DN']] <- TRUE
  MAP[['DR']] <- TRUE
  MAP[['DT']] <- TRUE
  MAP[['DW']] <- TRUE
  MAP[['GB']] <- TRUE
  MAP[['GD']] <- TRUE
  MAP[['GG']] <- TRUE
  MAP[['GK']] <- TRUE
  MAP[['GN']] <- TRUE
  MAP[['GR']] <- TRUE
  MAP[['GS']] <- TRUE
  MAP[['GV']] <- TRUE
  MAP[['HA']] <- TRUE
  MAP[['HC']] <- TRUE
  MAP[['HH']] <- TRUE
  MAP[['HM']] <- TRUE
  MAP[['HN']] <- TRUE
  MAP[['HT']] <- TRUE
  MAP[['HW']] <- TRUE
  MAP[['HY']] <- TRUE
  MAP[['KB']] <- TRUE
  MAP[['KD']] <- TRUE
  MAP[['KG']] <- TRUE
  MAP[['KK']] <- TRUE
  MAP[['KN']] <- TRUE
  MAP[['KT']] <- TRUE
  MAP[['MA']] <- TRUE
  MAP[['MC']] <- TRUE
  MAP[['MH']] <- TRUE
  MAP[['MM']] <- TRUE
  MAP[['MN']] <- TRUE
  MAP[['MV']] <- TRUE
  MAP[['NA']] <- TRUE
  MAP[['NB']] <- TRUE
  MAP[['NC']] <- TRUE
  MAP[['ND']] <- TRUE
  MAP[['NG']] <- TRUE
  MAP[['NH']] <- TRUE
  MAP[['NK']] <- TRUE
  MAP[['NM']] <- TRUE
  MAP[['NN']] <- TRUE
  MAP[['NR']] <- TRUE
  MAP[['NS']] <- TRUE
  MAP[['NT']] <- TRUE
  MAP[['NV']] <- TRUE
  MAP[['NW']] <- TRUE
  MAP[['NY']] <- TRUE
  MAP[['RA']] <- TRUE
  MAP[['RD']] <- TRUE
  MAP[['RG']] <- TRUE
  MAP[['RN']] <- TRUE
  MAP[['RR']] <- TRUE
  MAP[['RV']] <- TRUE
  MAP[['SB']] <- TRUE
  MAP[['SC']] <- TRUE
  MAP[['SG']] <- TRUE
  MAP[['SN']] <- TRUE
  MAP[['SS']] <- TRUE
  MAP[['SV']] <- TRUE
  MAP[['TB']] <- TRUE
  MAP[['TD']] <- TRUE
  MAP[['TH']] <- TRUE
  MAP[['TK']] <- TRUE
  MAP[['TN']] <- TRUE
  MAP[['TT']] <- TRUE
  MAP[['TW']] <- TRUE
  MAP[['TY']] <- TRUE
  MAP[['VA']] <- TRUE
  MAP[['VC']] <- TRUE
  MAP[['VG']] <- TRUE
  MAP[['VM']] <- TRUE
  MAP[['VN']] <- TRUE
  MAP[['VR']] <- TRUE
  MAP[['VS']] <- TRUE
  MAP[['VV']] <- TRUE
  MAP[['WA']] <- TRUE
  MAP[['WD']] <- TRUE
  MAP[['WH']] <- TRUE
  MAP[['WN']] <- TRUE
  MAP[['WT']] <- TRUE
  MAP[['WW']] <- TRUE
  MAP[['YB']] <- TRUE
  MAP[['YC']] <- TRUE
  MAP[['YH']] <- TRUE
  MAP[['YN']] <- TRUE
  MAP[['YT']] <- TRUE
  MAP[['YY']] <- TRUE
  
  MAP
}
IUPAC_EQ <- INIT_IUPAC_EQ()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Temp Code to be deleted::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  ret_tib <- readr::read_tsv(file, # guess_max=1000000)
                             col_types = cols(
                               Seq_ID = col_character(),
                               Forward_Sequence = col_character(),
                               Genome_Build = col_character(),
                               Chromosome = col_character(),
                               Coordinate = col_double(),
                               Design_State = col_character(),
                               Seq_Length = col_double(),
                               Forward_CpG_Coord = col_double(),
                               TB_Strand = col_character(),
                               Top_Sequence = col_character(),
                               Top_CpG_Coord = col_double(),
                               Probe_Type = col_character(),
                               Probeset_ID = col_character(),
                               Probeset_Score = col_double(),
                               Methyl_Probe_ID = col_character(),
                               Methyl_Probe_Sequence = col_character(),
                               Methyl_Probe_Length = col_double(),
                               Methyl_Start_Coord = col_double(),
                               Methyl_End_Coord = col_double(),
                               Methyl_Probe_Covered_Top_Sequence = col_character(),
                               Methyl_Allele_FR_Strand = col_character(),
                               Methyl_Allele_TB_Strand = col_character(),
                               Methyl_Allele_CO_Strand = col_character(),
                               Methyl_Allele_Type = col_character(),
                               Methyl_Final_Score = col_double(),
                               Methyl_Tm = col_double(),
                               Methyl_Tm_Score = col_double(),
                               Methyl_GC_Percent = col_double(),
                               Methyl_GC_Score = col_double(),
                               Methyl_13mer_Count = col_double(),
                               Methyl_13mer_Score = col_double(),
                               Methyl_Address_Count = col_double(),
                               Methyl_Address_Score = col_double(),
                               Methyl_Self_Complementarity = col_double(),
                               Methyl_Self_Complementarity_Score = col_double(),
                               Methyl_Mono_Run = col_double(),
                               Methyl_Mono_Run_Score = col_double(),
                               Methyl_Ectopic_Count = col_double(),
                               Methyl_Ectopic_Score = col_double(),
                               Methyl_Underlying_CpG_Count = col_double(),
                               Methyl_Underlying_CpG_Min_Dist = col_double(),
                               Methyl_Underlying_CpG_Score = col_double(),
                               Methyl_In_CpG_Island_Relaxed = col_logical(),
                               Methyl_CpG_Island_Score = col_double(),
                               Methyl_Next_Base = col_character(),
                               Methyl_Next_Base_Score = col_double(),
                               UnMethyl_Probe_ID = col_character(),
                               UnMethyl_Probe_Sequence = col_character(),
                               UnMethyl_Probe_Length = col_double(),
                               UnMethyl_Start_Coord = col_double(),
                               UnMethyl_End_Coord = col_double(),
                               UnMethyl_Probe_Covered_Top_Sequence = col_character(),
                               UnMethyl_Allele_FR_Strand = col_character(),
                               UnMethyl_Allele_TB_Strand = col_character(),
                               UnMethyl_Allele_CO_Strand = col_character(),
                               UnMethyl_Allele_Type = col_character(),
                               UnMethyl_Final_Score = col_double(),
                               UnMethyl_Tm = col_double(),
                               UnMethyl_Tm_Score = col_double(),
                               UnMethyl_GC_Percent = col_double(),
                               UnMethyl_GC_Score = col_double(),
                               UnMethyl_13mer_Count = col_double(),
                               UnMethyl_13mer_Score = col_double(),
                               UnMethyl_Address_Count = col_double(),
                               UnMethyl_Address_Score = col_double(),
                               UnMethyl_Self_Complementarity = col_double(),
                               UnMethyl_Self_Complementarity_Score = col_double(),
                               UnMethyl_Mono_Run = col_double(),
                               UnMethyl_Mono_Run_Score = col_double(),
                               UnMethyl_Ectopic_Count = col_double(),
                               UnMethyl_Ectopic_Score = col_double(),
                               UnMethyl_Underlying_CpG_Count = col_double(),
                               UnMethyl_Underlying_CpG_Min_Dist = col_double(),
                               UnMethyl_Underlying_CpG_Score = col_double(),
                               UnMethyl_In_CpG_Island_Relaxed = col_logical(),
                               UnMethyl_CpG_Island_Score = col_double(),
                               UnMethyl_Next_Base = col_character(),
                               UnMethyl_Next_Base_Score = col_double()
                             )
  )
}


# End of file
