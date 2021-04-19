
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
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'template_func'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # ret_cnt <- ret_tib %>% base::nrow()
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
#                           BED/Tabix Methods::
#
# TBD:: Move this somewhere else
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

write_tabix_bed = function(tib, dir, name, 
                           head_char=NULL, s_idx=1, b_idx=2, e_idx=3,
                           header=FALSE,
                           ref=NULL, pipe=NULL,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'write_tabix_bed'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    
    ret_tib  <- tib
    bed_file <- file.path(dir, paste(name,'sorted.bed', sep='.'))
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing='{bed_file}'...{RET}"))
    if (header) {
      readr::write_tsv(ret_tib, bed_file) # , col_names = FALSE)
    } else {
      readr::write_tsv(ret_tib, bed_file, col_names = FALSE)
    }
    cmd <- glue::glue("bgzip -f {bed_file}")
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Running='{cmd}'...{RET}"))
    system(cmd)
    bed_file <- file.path(dir, paste(name,'sorted.bed.gz', sep='.'))
    
    cmd <- glue::glue("tabix")
    if (!is.null(head_char)) cmd <- glue::glue("{cmd} -c {head_char}")
    cmd <- glue::glue("{cmd} -s {s_idx} -b {b_idx} -e {e_idx} -f {bed_file}")
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Running='{cmd}'...{RET}"))
    system(cmd)
    
    if (!is.null(ref) && file.exists(ref)) {
      ref_name <- ref %>% base::basename() %>%
        stringr::str_remove(".gz$") %>% 
        stringr::str_remove(".bed$") %>%
        stringr::str_remove(".vcf$") %>%
        stringr::str_remove(".csv$") %>%
        stringr::str_remove(".tsv$")
        
      int_file <- file.path(dir, paste(name,ref_name,"intersect.tsv.gz", sep='.'))
      cmd <- glue::glue("tabix -h -R {bed_file} {ref}")
      if (!is.null(pipe)) cmd <- glue::glue("{cmd} | {pipe}")
      cmd <- glue::glue("{cmd} | gzip -c -> {int_file}")
      
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Running='{cmd}'...{RET}"))
      system(cmd)
      
      ret_tib <- load_dbSNP_vcf(int_file, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    } else {
      ret_tib <- bed_file
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
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

seq_to_prbs = function(tib, seq,
                       cgn="Seq_ID",chr="Chromosome", pos="Coordinate",
                       srd="F",
                       ups_len=60, seq_len=122,
                       iupac=NULL, add_flank=FALSE,
                       del="_", 
                       verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'seq_to_prbs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    cgn_sym <- rlang::sym(cgn)
    chr_sym <- rlang::sym(chr)
    pos_sym <- rlang::sym(pos)

    cur_begs <- NULL
    cur_ends <- NULL
    cur_begs <- tib %>% dplyr::pull(!!pos_sym) - ups_len - 2
    cur_ends <- cur_begs + seq_len - 1 + 4
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_begs={RET}"))
      cur_begs %>% head() %>% print()
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_ends={RET}"))
      cur_ends %>% head() %>% print()
    }
    
    # Sub-string ref sequences::
    #
    ref_seqs <- NULL
    ref_seqs <- 
      stringr::str_sub( as.character(dna_seqs[[chr_idx]]), cur_begs, cur_ends) %>%
      stringr::str_to_upper()
    
    if (srd=="R") ref_seqs <- ref_seqs %>% cmpl()
    
    # Build individual parts of the template sequence::
    # up01.up02...up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn02.dn01
    #                             [  C   G  ]
    #
    ref_up01 <- stringr::str_sub(ref_seqs, 1,1)
    ref_up02 <- stringr::str_sub(ref_seqs, 2,2)
    ref_up58 <- stringr::str_sub(ref_seqs, 3,ups_len-2+2)
    
    ref_up59 <- stringr::str_sub(ref_seqs, ups_len-2+3,ups_len-2+3)
    ref_up60 <- stringr::str_sub(ref_seqs, ups_len-2+4,ups_len-2+4)
    ref_up61 <- stringr::str_sub(ref_seqs, ups_len-2+5,ups_len-2+5)
    
    iupac_vec <- ref_up61
    if (!is.null(iupac)) iupac_vec <- cur_tib %>% dplyr::pull(!!iupac)
    
    ref_dn61 <- stringr::str_sub(ref_seqs, ups_len-2+6,ups_len-2+6)
    ref_dn60 <- stringr::str_sub(ref_seqs, ups_len-2+7,ups_len-2+7)
    ref_dn59 <- stringr::str_sub(ref_seqs, ups_len-2+8,ups_len-2+8)
    
    ref_dn58 <- stringr::str_sub(ref_seqs, ups_len-2+9,ups_len-2+ups_len-2+8)
    ref_dn02 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+8,ups_len-2+ups_len-2+8)
    ref_dn01 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+9,ups_len-2+ups_len-2+9)
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ref_seqs({chr_str})={RET}"))
      ref_seqs %>% head() %>% print()
    }
    
    # Add all fields to current tibble::
    #
    cur_tib <- cur_tib %>%
      dplyr::mutate(
        !!ext_sym := ref_seqs %>% addBrac(),
        
        up01=ref_up01,
        up02=ref_up02,
        up58=ref_up58,
        
        up59=ref_up59,
        up60=ref_up60,
        up61=ref_up61,
        
        iupac=iupac_vec,
        
        dn61=ref_dn61,
        dn60=ref_dn60,
        dn59=ref_dn59,
        
        dn58=ref_dn58,
        dn02=ref_dn02,
        dn01=ref_dn01,
        
        # Now we can assemble to optimal template::
        #
        !!des_sym := paste0(up58,
                            up59,up60,
                            iupac,dn61,
                            dn60,dn59,
                            dn58) %>%
          addBrac(),
        !!imp_sym := paste0(up58,
                            up59,up60,
                            "CG",
                            dn60,dn59,
                            dn58) %>%
          addBrac(),
        
        # TBD:: Pretty sure this right, but really needs more testing...
        #
        Din_Str=stringr::str_to_lower(paste0(iupac,dn61)),
        Des_Din=dplyr::case_when(
          up61=='C' & dn61=='G' ~ "cg",
          up61=='C' & dn61!='G' ~ "ch",
          up61!='C' ~ "rs",
          TRUE ~ NA_character_
        )
      )
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_tib({chr_str})={RET}"))
      cur_tib %>% print()
    }
    
    #  
    # TBD:: We can also break out tri-fecta sites if probe type == rs
    #
    ups_tib <- NULL
    dns_tib <- NULL
    
    if (add_flank) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Upstream CG Flank Seuqnces({chr_str})...{RET}"))
      
      ups_tib <- cur_tib %>% 
        dplyr::filter(up59=="C" & up60=="G") %>%
        dplyr::mutate(
          !!des_sym := paste0(up01,up02,up58,up59,up60,iupac,dn61,dn60,dn59,dn58) %>%
            stringr::str_sub(1,seq_len) %>%
            addBrac(),
          Din_Str=stringr::str_to_lower(paste0(up59,up60)),
          Des_Din="cg",
          Name=paste(!!name_sym,"CG_UP", sep=del),
          Coordinate=as.integer(Coordinate-2)
        )
      
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Downstream CG Flank Seuqnces({chr_str})...{RET}"))
      
      dns_tib <- cur_tib %>% 
        dplyr::filter(dn61=="C" & dn60=="G") %>%
        dplyr::mutate(
          !!des_sym := paste0(up58,up59,up60,iupac,dn61,dn60,dn59,dn58,dn02) %>%
            stringr::str_sub(2) %>%
            addBrac(),
          Din_Str=stringr::str_to_lower(paste0(dn61,dn60)),
          Des_Din="cg",
          Name=paste(!!name_sym,"CG_DN", sep=del),
          Coordinate=as.integer(Coordinate+2)
        )
    }
    
    # Add New Tri-fecta Probes::
    #  - NOTE:: These will have their names changed later after alignment...
    #
    if (verbose>=vt+1)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joining and sorting...{RET}"))
    
    cur_tib <- dplyr::bind_rows(cur_tib, ups_tib, dns_tib) %>%
      dplyr::arrange(!!name_sym)
    
    ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
    
    if (verbose>=vt+1)
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Substring chr_str={chr_str}.{RET}{RET}"))
    
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

bed_to_prbs = function(tib, fas,
                       file=NULL,
                       din="Probe_Type",gen="na",
                       cgn="Seq_ID",chr="Chromosome", pos="Coordinate",
                       fwd_seq="Fwd_Temp_Seq",
                       iup_seq="Iup_Temp_Seq",
                       imp_seq="Imp_Temp_Seq",
                       ups_len=60, seq_len=122, nrec=0,
                       iupac=NULL, add_flank=FALSE,
                       del="_", 
                       verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'bed_to_prbs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    cgn_sym <- rlang::sym(cgn)
    chr_sym <- rlang::sym(chr)
    pos_sym <- rlang::sym(pos)
    
    # Get List of BSC Genomes::
    # fas_dir <- fas %>% base::dirname()
    fas_pre <- fas %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$")
    if (is.null(gen)) gen <- base::basename(fas_pre)
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Prefix({gen})={fas_pre}.{RET}"))
    
    # There isn't an opposite genome, need to build that from the converted::
    # srd_cos <- c("C","O")
    
    srd_cos <- c("C")
    srd_frs <- c("F","R")
    srd_mus <- c("U","M","D")

    for (fr in srd_frs) {
      # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} fr={fr}...{RET}"))
      
      for (co in srd_cos) {
        # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} co={co}...{RET}"))
        
        for (mu in srd_mus) {
          if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} fr={fr}, co={co}, mu={mu}...{RET}"))
          gen_fas <- paste0(fas_pre,".", fr,co,mu, ".fa.gz")
          
          if (!file.exists(gen_fas)) {
            if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Skipping; gen_fas=({fr},{co},{mu})",
                                            "={gen_fas} Does NOT Exist!{RET}"))
            next
          } else {
            if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading; gen_fas({fr},{co},{mu})",
                                            "={gen_fas}...{RET}"))
          }
          
          cur_tib <- fas_to_seq(
            tib=tib, fas=gen_fas,
            file=file,
            name=cgn,din=din,gen=gen,
            chr1=chr,pos=pos,
            srd=fr,
            fwd_seq=fwd_seq,iup_seq=iup_seq,imp_seq=imp_seq,
            
            # See Above:: Will Pass these variables in later...
            # fwd_seq="Fwd_Temp_Seq",des_seq="Iup_Temp_Seq",imp_seq="Imp_Temp_Seq",
            
            iupac=iupac,
            nrec=nrec,
            add_flank=FALSE,
            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt
          ) %>%
            dplyr::mutate(
              Srd_FR=fr,
              Srd_CO=co,
              Prb_Des=mu
            )
          
          if (verbose>=vt+4) {
            cat(glue::glue("[{funcTag}]:{tabsStr} Current Tib Pre Join={RET}"))
            print(cur_tib)
          }
          
          # TBD:: Add the current FR/CO/TB fields to cur_tib
          # TBD:: Add the opposite designs
          # TBD:: Remove previous information from the input tibble
          
          prb_tib <- NULL
          if (fr=="F") {
            prb_tib <- tibble::tibble(
              # TBD:: This should be pulled by symbolic name
              Seq_ID = cur_tib$Seq_ID,
              
              # Forward Converted Designs:: Infinium I/II
              Prb_SeqC1 = paste0(
                cur_tib$up61,
                cur_tib$dn61,
                cur_tib$dn60,
                cur_tib$dn59,
                cur_tib$dn58 %>% stringr::str_sub(1,46)
              ) %>% revCmp(),
              Prb_SeqC2 = paste0(
                # cur_tib$up61,
                cur_tib$dn61,
                cur_tib$dn60,
                cur_tib$dn59,
                cur_tib$dn58 %>% stringr::str_sub(1,47)
              ) %>% revCmp(),
              
              # Forward Opposite Designs:: Infinium I/II
              #   Just need to back up by one base::
              Prb_SeqO1 = paste0(
                cur_tib$up58 %>% stringr::str_sub(12,58), # 13 -> 12
                cur_tib$up59,
                cur_tib$up60,
                cur_tib$up61
                # cur_tib$dn61 # Remove dn61
              ), # No complement
              Prb_SeqO2 = paste0(
                cur_tib$up58 %>% stringr::str_sub(11,58), # 12 -> 11
                cur_tib$up59,
                cur_tib$up60
                # cur_tib$up61 # Remove dn61
                # cur_tib$dn61,
                # cur_tib$dn60
              ) # No Complement
              
            )
          } else if (fr=="R") {
            #
            # NOTE: From the ch/rs perspective the reverse strand 
            #   coordinates frame should be moved back by 1
            #
            # if ()
            prb_tib <- tibble::tibble(
              Seq_ID = cur_tib$Seq_ID,

              # Reverse Converted Designs:: Infinium I/II
              Prb_SeqC1 = paste0(
                cur_tib$up58 %>% stringr::str_sub(13,58),
                cur_tib$up59,
                cur_tib$up60,
                cur_tib$up61,
                cur_tib$dn61
              ) %>% cmpl(),
              Prb_SeqC2 = paste0(
                cur_tib$up58 %>% stringr::str_sub(12,58),
                cur_tib$up59,
                cur_tib$up60,
                cur_tib$up61
                # cur_tib$dn61,
                # cur_tib$dn60
              ) %>% cmpl(),
              
              # Reverse Opposite Designs:: Infinium I/II
              #   Just need to move forward by 1::
              Prb_SeqO1 = paste0(
                # cur_tib$up61, # Remove up61
                cur_tib$dn61,
                cur_tib$dn60,
                cur_tib$dn59,
                cur_tib$dn58 %>% stringr::str_sub(1,47) # 46 -> 47
              ), # %>% revCmp(), # No revCmp
              Prb_SeqO2 = paste0(
                # cur_tib$up61,
                # cur_tib$dn61, # Remove dn61
                cur_tib$dn60,
                cur_tib$dn59,
                cur_tib$dn58 %>% stringr::str_sub(1,48) # 47 -> 48
              ) # %>% revCmp(), # No revCmp

            )
          } else {
            stop(glue::glue("{RET}[{par$prgmTag}]: Unsupported fr={fr}...{RET}"))
            return(NULL)
          }
          cur_tib <- cur_tib %>% dplyr::left_join(prb_tib, by="Seq_ID")
          ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
        }
      }
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

fas_to_seq = function(tib, fas, 
                      file=NULL,
                      name,din,gen="na",
                      chr1="Chromosome", pos="Coordinate",
                      chr2="Chrom_Char", srd="F",
                      fwd_seq="Fwd_Temp_Seq",
                      iup_seq="Iup_Temp_Seq",
                      imp_seq="Imp_Temp_Seq",
                      ups_len=60, seq_len=122, nrec=0,
                      iupac=NULL,
                      ref_col="Ref",alt_col="Alt",iup_col="Iupac",
                      add_flank=FALSE,del="_", 
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'fas_to_seq'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (nrec==0) {
      dna_seqs <- 
        Biostrings::readDNAStringSet(filepath=fas, format="fasta")
    } else {
      dna_seqs <- 
        Biostrings::readDNAStringSet(filepath=fas, format="fasta", nrec=nrec)
    }
    
    dna_maps <- dna_seqs %>% 
      names() %>% 
      stringr::str_remove(" .*$") %>% 
      stringr::str_remove("^chr") %>%
      tibble::tibble() %>% 
      purrr::set_names("Chrom_Char") %>% 
      dplyr::mutate(Idx=dplyr::row_number(),
                    Chrom_Char=paste0("chr",Chrom_Char)
      )
    
    pos_sym  <- rlang::sym(pos)
    name_sym <- rlang::sym(name)
    chr1_sym <- rlang::sym(chr1)
    chr2_sym <- rlang::sym(chr2)
    fwd_sym  <- rlang::sym(fwd_seq)
    iup_sym  <- rlang::sym(iup_seq)
    imp_sym  <- rlang::sym(imp_seq)
    din_sym  <- rlang::sym(din)
    
    ref_col_sym  <- rlang::sym(ref_col)
    alt_col_sym  <- rlang::sym(alt_col)
    iup_col_sym  <- rlang::sym(iup_col)
    
    chr_list <- tib %>% 
      dplyr::arrange(!!chr1_sym,!!pos_sym) %>%
      split(.[[chr1]])
    
    chr_maps <- dna_maps %>%
      split(.[[chr2]])
    
    chr_names <- chr_list %>% names()
    for (chr_str in chr_names) {
      cur_tib <- chr_list[[chr_str]]
      
      if (is.null(chr_maps[[chr_str]])) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Failed to find chr_str={chr_str} ",
                         "in map; Skipping...{RET}"))
        next
      }
      
      chr_idx <- chr_maps[[chr_str]] %>% 
        dplyr::filter(!!chr2_sym==chr_str) %>% 
        head(n=1) %>% pull(Idx) %>% as.integer()
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Substring chr(str/idx)=({chr_str}/{chr_idx})={RET}"))
        print(cur_tib)
      }
      
      # Subtract two for cg upstream 122mer formation and add two+two as well
      #
      cur_begs <- NULL
      cur_ends <- NULL
      cur_begs <- cur_tib %>% dplyr::pull(!!pos_sym) - ups_len - 2
      cur_ends <- cur_begs + seq_len - 1 + 4
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_begs={RET}"))
        cur_begs %>% head() %>% print()
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_ends={RET}"))
        cur_ends %>% head() %>% print()
      }
      
      # Sub-string ref sequences::
      #
      ref_seqs <- NULL
      ref_seqs <- 
        stringr::str_sub( as.character(dna_seqs[[chr_idx]]), cur_begs, cur_ends) %>%
        stringr::str_to_upper()
      
      if (srd=="R") ref_seqs <- ref_seqs %>% cmpl()
      
      # Build individual parts of the template sequence::
      # up01.up02...up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn02.dn01
      #                             [  C   G  ]
      #
      ref_up01 <- stringr::str_sub(ref_seqs, 1,1)
      ref_up02 <- stringr::str_sub(ref_seqs, 2,2)
      ref_up58 <- stringr::str_sub(ref_seqs, 3,ups_len-2+2)
      
      ref_up59 <- stringr::str_sub(ref_seqs, ups_len-2+3,ups_len-2+3)
      ref_up60 <- stringr::str_sub(ref_seqs, ups_len-2+4,ups_len-2+4)
      ref_up61 <- stringr::str_sub(ref_seqs, ups_len-2+5,ups_len-2+5)
      
      iupac_vec <- ref_up61
      if (!is.null(iupac)) iupac_vec <- cur_tib %>% dplyr::pull(!!iupac)
      
      ref_dn61 <- stringr::str_sub(ref_seqs, ups_len-2+6,ups_len-2+6)
      ref_dn60 <- stringr::str_sub(ref_seqs, ups_len-2+7,ups_len-2+7)
      ref_dn59 <- stringr::str_sub(ref_seqs, ups_len-2+8,ups_len-2+8)
      
      ref_dn58 <- stringr::str_sub(ref_seqs, ups_len-2+9,ups_len-2+ups_len-2+8)
      ref_dn02 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+8,ups_len-2+ups_len-2+8)
      ref_dn01 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+9,ups_len-2+ups_len-2+9)
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ref_seqs({chr_str})={RET}"))
        ref_seqs %>% head() %>% print()
      }
      
      # Add all fields to current tibble::
      #
      cur_tib <- cur_tib %>%
        dplyr::mutate(
          !!fwd_sym := ref_seqs %>% addBrac(),
          
          up01=ref_up01,
          up02=ref_up02,
          up58=ref_up58,
          
          up59=ref_up59,
          up60=ref_up60,
          up61=ref_up61,
          
          iupac=iupac_vec,
          
          dn61=ref_dn61,
          dn60=ref_dn60,
          dn59=ref_dn59,
          
          dn58=ref_dn58,
          dn02=ref_dn02,
          dn01=ref_dn01,
          
          # Now we can assemble to optimal template::
          #
          !!iup_sym := paste0(up58,
                              up59,up60,
                              iupac,dn61,
                              dn60,dn59,
                              dn58) %>%
            addBrac(),
          !!imp_sym := paste0(up58,
                              up59,up60,
                              "CG",
                              dn60,dn59,
                              dn58) %>%
            addBrac(),
          
          # TBD:: Pretty sure this right, but really needs more testing...
          #
          Din_Str=stringr::str_to_lower(paste0(iupac,dn61)),
          Des_Din=dplyr::case_when(
            up61=='C' & dn61=='G' ~ "cg",
            up61=='C' & dn61!='G' ~ "ch",
            up61!='C' ~ "rs",
            TRUE ~ NA_character_
          )
        )
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_tib({chr_str})={RET}"))
        cur_tib %>% print()
      }
      
      #  
      # TBD:: We can also break out tri-fecta sites if probe type == rs
      #
      ups_tib <- NULL
      dns_tib <- NULL
      
      if (add_flank) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Upstream CG Flank Seuqnces({chr_str})...{RET}"))
        
        ups_tib <- cur_tib %>% 
          dplyr::filter(!!din_sym == "rs") %>%
          dplyr::filter(up59=="C" & up60=="G") %>%
          dplyr::mutate(
            !!din_sym := "cg",
            !!fwd_sym := paste0(up01,up02,up58,up59,up60,up61,dn61,dn60,dn59,dn58) %>%
              stringr::str_sub(1,seq_len) %>%
              addBrac(),
            !!iup_sym := paste0(up01,up02,up58,up59,up60,iupac,dn61,dn60,dn59,dn58) %>%
              stringr::str_sub(1,seq_len) %>%
              addBrac(),
            !!imp_sym := paste0(up01,up02,up58,"CG",up61,dn61,dn60,dn59,dn58) %>%
              stringr::str_sub(1,seq_len) %>%
              addBrac(),
            !!ref_col_sym := "C",
            !!alt_col_sym := "T",
            !!iup_col_sym := "Y",
            Din_Str=stringr::str_to_lower(paste0(up59,up60)),
            Des_Din="cg",
            !!name_sym:=paste(!!name_sym,"CG-UP", sep='-'),
            Coordinate=as.integer(Coordinate-2)
          )
        
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Downstream CG Flank Seuqnces({chr_str})...{RET}"))
        
        dns_tib <- cur_tib %>% 
          dplyr::filter(!!din_sym == "rs") %>%
          dplyr::filter(dn61=="C" & dn60=="G") %>%
          dplyr::mutate(
            !!din_sym := "cg",
            !!fwd_sym := paste0(up58,up59,up60,up61,dn61,dn60,dn59,dn58,dn02) %>%
              stringr::str_sub(2) %>%
              addBrac(),
            !!iup_sym := paste0(up58,up59,up60,iupac,dn61,dn60,dn59,dn58,dn02) %>%
              stringr::str_sub(2) %>%
              addBrac(),
            !!imp_sym := paste0(up58,up59,up60,up61,"CG",dn59,dn58,dn02) %>%
              stringr::str_sub(2) %>%
              addBrac(),
            !!ref_col_sym := "C",
            !!alt_col_sym := "T",
            !!iup_col_sym := "Y",
            Din_Str=stringr::str_to_lower(paste0(dn61,dn60)),
            Des_Din="cg",
            !!name_sym:=paste(!!name_sym,"CG-DN", sep='-'),
            Coordinate=as.integer(Coordinate+2)
          )
      }
      
      # Add New Tri-fecta Probes::
      #  - NOTE:: These will have their names changed later after alignment...
      #
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joining and sorting...{RET}"))
      
      cur_tib <- dplyr::bind_rows(cur_tib, ups_tib, dns_tib) %>%
        dplyr::arrange(!!name_sym)
      
      ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
      
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Substring chr_str={chr_str}.{RET}{RET}"))
    }

    if (!is.null(file)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing improbe input TSV={file}...{RET}"))
      
      imp_col <- c("Seq_ID","Sequence","Genome_Build","Chromosome","Coordinate","CpG_Island")
      imp_tib <- ret_tib %>% 
        dplyr::mutate(Genome_Build=!!gen, CpG_Island="FALSE") %>%
        dplyr::select(dplyr::all_of(c(!!name, !!imp_seq, "Genome_Build", !!chr1, !!pos, "CpG_Island")) ) %>%
        purrr::set_names(imp_col)
      
      safe_write(x=imp_tib,type="tsv",file=file,funcTag=funcTag, 
                 verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} DIN Summary={RET}"))
      
      sum_tib <- ret_tib %>% 
        dplyr::group_by(up61,dn61,Des_Din) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      sum_tib %>% print(n=base::nrow(sum_tib))
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

dna_to_template1 = function(tib, fas, # map, 
                            file=NULL,
                            name="Name",gen="na",
                            chr1="Chromosome", pos="Coordinate",
                            chr2="Chrom_Char",
                            ext_seq="Fwd_Seq",
                            des_seq="Des_Seq",
                            imp_seq="Sequence",
                            ups_len=60, seq_len=122, 
                            iupac=NULL, add_flank=FALSE,
                            del="_", 
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'dna_to_template1'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    dna_seqs <- 
      Biostrings::readDNAStringSet(filepath = fas, format = "fasta", nrec = 2)
    
    dna_maps <- dna_seqs %>% 
      names() %>% 
      stringr::str_remove(" .*$") %>% 
      stringr::str_remove("^chr") %>%
      tibble::tibble() %>% 
      purrr::set_names("Chrom_Char") %>% 
      dplyr::mutate(Idx=dplyr::row_number(),
                    Chrom_Char=paste0("chr",Chrom_Char)
      )
    
    pos_sym  <- rlang::sym(pos)
    name_sym <- rlang::sym(name)
    chr1_sym <- rlang::sym(chr1)
    chr2_sym <- rlang::sym(chr2)
    ext_sym  <- rlang::sym(ext_seq)
    des_sym  <- rlang::sym(des_seq)
    imp_sym  <- rlang::sym(imp_seq)
    
    chr_list <- tib %>% 
      dplyr::arrange(!!chr1_sym,!!pos_sym) %>%
      split(.[[chr1]])
    
    chr_maps <- dna_maps %>%
      split(.[[chr2]])
    
    chr_names <- chr_list %>% names()
    for (chr_str in chr_names) {
      cur_tib <- chr_list[[chr_str]]
      
      if (is.null(chr_maps[[chr_str]])) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Failed to find chr_str={chr_str} ",
                         "in map; Skipping...{RET}"))
        next
      }
      
      chr_idx <- chr_maps[[chr_str]] %>% 
        dplyr::filter(!!chr2_sym==chr_str) %>% 
        head(n=1) %>% pull(Idx) %>% as.integer()
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Substring chr(str/idx)=({chr_str}/{chr_idx})={RET}"))
        print(cur_tib)
      }
      
      # Subtract two for cg upstream 122mer formation and add two+two as well
      #
      cur_begs <- NULL
      cur_ends <- NULL
      cur_begs <- cur_tib %>% dplyr::pull(!!pos_sym) - ups_len - 2
      cur_ends <- cur_begs + seq_len - 1 + 4
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_begs={RET}"))
        cur_begs %>% head() %>% print()
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_ends={RET}"))
        cur_ends %>% head() %>% print()
      }
      
      # Substring ref sequences::
      #
      ref_seqs <- NULL
      ref_seqs <- stringr::str_sub( as.character(dna_seqs[[chr_idx]]), cur_begs, cur_ends) %>%
        stringr::str_to_upper()
      
      # Build individual parts of the template sequence::
      # up01.up02...up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn02.dn01
      #                             [  C   G  ]
      #
      ref_up01 <- stringr::str_sub(ref_seqs, 1,1)
      ref_up02 <- stringr::str_sub(ref_seqs, 2,2)
      ref_up58 <- stringr::str_sub(ref_seqs, 3,ups_len-2+2)
      
      ref_up59 <- stringr::str_sub(ref_seqs, ups_len-2+3,ups_len-2+3)
      ref_up60 <- stringr::str_sub(ref_seqs, ups_len-2+4,ups_len-2+4)
      ref_up61 <- stringr::str_sub(ref_seqs, ups_len-2+5,ups_len-2+5)
      
      iupac_vec <- ref_up61
      if (!is.null(iupac)) iupac_vec <- cur_tib %>% dplyr::pull(!!iupac)
      
      ref_dn61 <- stringr::str_sub(ref_seqs, ups_len-2+6,ups_len-2+6)
      ref_dn60 <- stringr::str_sub(ref_seqs, ups_len-2+7,ups_len-2+7)
      ref_dn59 <- stringr::str_sub(ref_seqs, ups_len-2+8,ups_len-2+8)
      
      ref_dn58 <- stringr::str_sub(ref_seqs, ups_len-2+9,ups_len-2+ups_len-2+8)
      ref_dn02 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+8,ups_len-2+ups_len-2+8)
      ref_dn01 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+9,ups_len-2+ups_len-2+9)
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ref_seqs({chr_str})={RET}"))
        ref_seqs %>% head() %>% print()
      }
      
      # Add all fields to current tibble::
      #
      cur_tib <- cur_tib %>%
        dplyr::mutate(
          !!ext_sym := ref_seqs %>% addBrac(),
          
          up01=ref_up01,
          up02=ref_up02,
          up58=ref_up58,
          
          up59=ref_up59,
          up60=ref_up60,
          up61=ref_up61,
          
          iupac=iupac_vec,
          
          dn61=ref_dn61,
          dn60=ref_dn60,
          dn59=ref_dn59,
          
          dn58=ref_dn58,
          dn02=ref_dn02,
          dn01=ref_dn01,
          
          # Now we can assemble to optimal template::
          #
          !!des_sym := paste0(up58,
                              up59,up60,
                              iupac,dn61,
                              dn60,dn59,
                              dn58) %>%
            addBrac(),
          !!imp_sym := paste0(up58,
                              up59,up60,
                              "CG",
                              dn60,dn59,
                              dn58) %>%
            addBrac(),
          
          # TBD:: Pretty sure this right, but really needs more testing...
          #
          Din_Str=stringr::str_to_lower(paste0(iupac,dn61)),
          Des_Din=dplyr::case_when(
            up61=='C' & dn61=='G' ~ "cg",
            up61=='C' & dn61!='G' ~ "ch",
            up61!='C' ~ "rs",
            TRUE ~ NA_character_
          )
        )
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_tib({chr_str})={RET}"))
        cur_tib %>% print()
      }
      
      #  
      # TBD:: We can also break out trifecta sites if probe type == rs
      #
      ups_tib <- NULL
      dns_tib <- NULL
      
      if (add_flank) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Upstream CG Flank Seuqnces({chr_str})...{RET}"))
        
        ups_tib <- cur_tib %>% 
          dplyr::filter(up59=="C" & up60=="G") %>%
          dplyr::mutate(
            !!des_sym := paste0(up01,up02,up58,up59,up60,iupac,dn61,dn60,dn59,dn58) %>%
              stringr::str_sub(1,seq_len) %>%
              addBrac(),
            Din_Str=stringr::str_to_lower(paste0(up59,up60)),
            Des_Din="cg",
            Name=paste(!!name_sym,"CG_UP", sep=del),
            Coordinate=as.integer(Coordinate-2)
          )
        
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Downstream CG Flank Seuqnces({chr_str})...{RET}"))
        
        dns_tib <- cur_tib %>% 
          dplyr::filter(dn61=="C" & dn60=="G") %>%
          dplyr::mutate(
            !!des_sym := paste0(up58,up59,up60,iupac,dn61,dn60,dn59,dn58,dn02) %>%
              stringr::str_sub(2) %>%
              addBrac(),
            Din_Str=stringr::str_to_lower(paste0(dn61,dn60)),
            Des_Din="cg",
            Name=paste(!!name_sym,"CG_DN", sep=del),
            Coordinate=as.integer(Coordinate+2)
          )
      }
      
      # Add New Tri-fecta Probes::
      #  - NOTE:: These will have their names changed later after alignment...
      #
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joining and sorting...{RET}"))
      
      cur_tib <- dplyr::bind_rows(cur_tib, ups_tib, dns_tib) %>%
        dplyr::arrange(!!name_sym)
      
      ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
      
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Substring chr_str={chr_str}.{RET}{RET}"))
    }
    
    if (!is.null(file)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing improbe input TSV={file}...{RET}"))
      
      imp_col <- c("Seq_ID","Sequence","Genome_Build","Chromosome","Coordinate","CpG_Island")
      imp_tib <- ret_tib %>% 
        dplyr::mutate(Genome_Build=!!gen, CpG_Island="FALSE") %>%
        dplyr::select(dplyr::all_of(c(!!name, !!imp_seq, "Genome_Build", !!chr1, !!pos, "CpG_Island")) ) %>%
        purrr::set_names(imp_col)
      
      safe_write(x=imp_tib,type="tsv",file=file,funcTag=funcTag, 
                 verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} DIN Summary={RET}"))
      
      sum_tib <- ret_tib %>% 
        dplyr::group_by(up61,dn61,Des_Din) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      sum_tib %>% print(n=base::nrow(sum_tib))
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

dna_to_template2 = function(tib, dna, map, file=NULL,
                            name="Name",gen="na",
                            chr1="Chromosome", pos="Coordinate",
                            chr2="Chrom_Char",
                            ext_seq="Fwd_Seq",
                            des_seq="Des_Seq",
                            imp_seq="Sequence",
                            ups_len=60, seq_len=122, 
                            iupac=NULL, add_flank=FALSE,
                            del="_", 
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'dna_to_template2'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    pos_sym  <- rlang::sym(pos)
    name_sym <- rlang::sym(name)
    chr1_sym <- rlang::sym(chr1)
    chr2_sym <- rlang::sym(chr2)
    ext_sym  <- rlang::sym(ext_seq)
    des_sym  <- rlang::sym(des_seq)
    imp_sym  <- rlang::sym(imp_seq)
    
    chr_list <- tib %>% 
      dplyr::arrange(!!chr1_sym,!!pos_sym) %>%
      split(.[[chr1]])
    
    chr_maps <- map %>%
      split(.[[chr2]])
    
    chr_names <- chr_list %>% names()
    for (chr_str in chr_names) {
      cur_tib <- chr_list[[chr_str]]
      
      if (is.null(chr_maps[[chr_str]])) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Failed to find chr_str={chr_str} ",
                         "in map; Skipping...{RET}"))
        next
      }
      
      chr_idx <- chr_maps[[chr_str]] %>% 
        dplyr::filter(!!chr2_sym==chr_str) %>% 
        head(n=1) %>% pull(Idx) %>% as.integer()
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Substring chr(str/idx)=({chr_str}/{chr_idx})={RET}"))
        print(cur_tib)
      }
      
      # Subtract two for cg upstream 122mer formation and add two+two as well
      #
      cur_begs <- NULL
      cur_ends <- NULL
      cur_begs <- cur_tib %>% dplyr::pull(!!pos_sym) - ups_len - 2
      cur_ends <- cur_begs + seq_len - 1 + 4
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_begs={RET}"))
        cur_begs %>% head() %>% print()
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_ends={RET}"))
        cur_ends %>% head() %>% print()
      }
      
      # Substring ref sequences::
      #
      ref_seqs <- NULL
      ref_seqs <- stringr::str_sub( as.character(dna[[chr_idx]]), cur_begs, cur_ends) %>%
        stringr::str_to_upper()
      
      # Build individual parts of the template sequence::
      # up01.up02...up58...up59.up60.up61.dn61.dn60.dn59...dn58...dn02.dn01
      #                             [  C   G  ]
      #
      ref_up01 <- stringr::str_sub(ref_seqs, 1,1)
      ref_up02 <- stringr::str_sub(ref_seqs, 2,2)
      ref_up58 <- stringr::str_sub(ref_seqs, 3,ups_len-2+2)
      
      ref_up59 <- stringr::str_sub(ref_seqs, ups_len-2+3,ups_len-2+3)
      ref_up60 <- stringr::str_sub(ref_seqs, ups_len-2+4,ups_len-2+4)
      ref_up61 <- stringr::str_sub(ref_seqs, ups_len-2+5,ups_len-2+5)
      
      iupac_vec <- ref_up61
      if (!is.null(iupac)) iupac_vec <- cur_tib %>% dplyr::pull(!!iupac)
      
      ref_dn61 <- stringr::str_sub(ref_seqs, ups_len-2+6,ups_len-2+6)
      ref_dn60 <- stringr::str_sub(ref_seqs, ups_len-2+7,ups_len-2+7)
      ref_dn59 <- stringr::str_sub(ref_seqs, ups_len-2+8,ups_len-2+8)
      
      ref_dn58 <- stringr::str_sub(ref_seqs, ups_len-2+9,ups_len-2+ups_len-2+8)
      ref_dn02 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+8,ups_len-2+ups_len-2+8)
      ref_dn01 <- stringr::str_sub(ref_seqs, ups_len-2+ups_len-2+9,ups_len-2+ups_len-2+9)
      
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ref_seqs({chr_str})={RET}"))
        ref_seqs %>% head() %>% print()
      }
      
      # Add all fields to current tibble::
      #
      cur_tib <- cur_tib %>%
        dplyr::mutate(
          !!ext_sym := ref_seqs %>% addBrac(),
          
          up01=ref_up01,
          up02=ref_up02,
          up58=ref_up58,
          
          up59=ref_up59,
          up60=ref_up60,
          up61=ref_up61,
          
          iupac=iupac_vec,
          
          dn61=ref_dn61,
          dn60=ref_dn60,
          dn59=ref_dn59,
          
          dn58=ref_dn58,
          dn02=ref_dn02,
          dn01=ref_dn01,
          
          # Now we can assemble to optimal template::
          #
          !!des_sym := paste0(up58,
                              up59,up60,
                              iupac,dn61,
                              dn60,dn59,
                              dn58) %>%
            addBrac(),
          !!imp_sym := paste0(up58,
                              up59,up60,
                              "CG",
                              dn60,dn59,
                              dn58) %>%
            addBrac(),
          
          # TBD:: Pretty sure this right, but really needs more testing...
          #
          Din_Str=stringr::str_to_lower(paste0(iupac,dn61)),
          Des_Din=dplyr::case_when(
            up61=='C' & dn61=='G' ~ "cg",
            up61=='C' & dn61!='G' ~ "ch",
            up61!='C' ~ "rs",
            TRUE ~ NA_character_
          )
        )
      if (verbose>=vt+4) {
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} cur_tib({chr_str})={RET}"))
        cur_tib %>% print()
      }
      
      #  
      # TBD:: We can also break out trifecta sites if probe type == rs
      #
      ups_tib <- NULL
      dns_tib <- NULL
      
      if (add_flank) {
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Upstream CG Flank Seuqnces({chr_str})...{RET}"))
        
        ups_tib <- cur_tib %>% 
          dplyr::filter(up59=="C" & up60=="G") %>%
          dplyr::mutate(
            !!des_sym := paste0(up01,up02,up58,up59,up60,iupac,dn61,dn60,dn59,dn58) %>%
              stringr::str_sub(1,seq_len) %>%
              addBrac(),
            Din_Str=stringr::str_to_lower(paste0(up59,up60)),
            Des_Din="cg",
            Name=paste(!!name_sym,"CG_UP", sep=del),
            Coordinate=as.integer(Coordinate-2)
          )
        
        if (verbose>=vt+1)
          cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Downstream CG Flank Seuqnces({chr_str})...{RET}"))
        
        dns_tib <- cur_tib %>% 
          dplyr::filter(dn61=="C" & dn60=="G") %>%
          dplyr::mutate(
            !!des_sym := paste0(up58,up59,up60,iupac,dn61,dn60,dn59,dn58,dn02) %>%
              stringr::str_sub(2) %>%
              addBrac(),
            Din_Str=stringr::str_to_lower(paste0(dn61,dn60)),
            Des_Din="cg",
            Name=paste(!!name_sym,"CG_DN", sep=del),
            Coordinate=as.integer(Coordinate+2)
          )
      }
      
      # Add New Tri-fecta Probes::
      #  - NOTE:: These will have their names changed later after alignment...
      #
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joining and sorting...{RET}"))
      
      cur_tib <- dplyr::bind_rows(cur_tib, ups_tib, dns_tib) %>%
        dplyr::arrange(!!name_sym)
      
      ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
      
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Substring chr_str={chr_str}.{RET}{RET}"))
    }
    
    if (!is.null(file)) {
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing improbe input TSV={file}...{RET}"))
      
      imp_col <- c("Seq_ID","Sequence","Genome_Build","Chromosome","Coordinate","CpG_Island")
      imp_tib <- ret_tib %>% 
        dplyr::mutate(Genome_Build=!!gen, CpG_Island="FALSE") %>%
        dplyr::select(dplyr::all_of(c(!!name, !!imp_seq, "Genome_Build", !!chr1, !!pos, "CpG_Island")) ) %>%
        purrr::set_names(imp_col)
      
      safe_write(x=imp_tib,type="tsv",file=file,funcTag=funcTag, 
                 verbose=verbose,vt=vt,tc=tc,append=FALSE)
    }
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} DIN Summary={RET}"))
      
      sum_tib <- ret_tib %>% 
        dplyr::group_by(up61,dn61,Des_Din) %>% 
        dplyr::summarise(Count=n(), .groups="drop")
      sum_tib %>% print(n=base::nrow(sum_tib))
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
#                          Genome Studio IO Basics::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

genome_studio_header = function(name="MethylationEPIC_v-1-0_B4", format="Infinium 2",
                                count=866895,date=NULL,
                                verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'genome_studio_header'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    if (is.null(date))
      date <- Sys.Date() %>% stringr::str_replace_all("-","/")
    
    ret_str <- glue::glue(
      "Illumina, Inc.,,,,,,,,,,,,,,,,,,,,,,,,,{RET}",
      "[Heading],,,,,,,,,,,,,,,,,,,,,,,,,,{RET}",
      "Descriptor File Name,{name}.csv,,,,,,,,,,,,,,,,,,,,,,,,,{RET}",
      "Assay Format,{format},,,,,,,,,,,,,,,,,,,,,,,,,{RET}",
      "Date Manufactured,{date},,,,,,,,,,,,,,,,,,,,,,,,,{RET}",
      "Loci Count ,{count},,,,,,,,,,,,,,,,,,,,,,,,,{RET}",
      "[Assay],,,,,,,,,,,,,,,,,,,,,,,,,,"
    )
    
    # ret_cnt <- ret_tib %>% base::nrow()
    # ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_str
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Genomic Range Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

intersect_seq = function(ref, can, out, idxA=1, idxB=1,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'intersect_seq'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  int_seq_cols <-
    cols(
      Imp_Seq  = col_character(),
      Imp_Nuc  = col_character(),
      
      Imp_SrdI = col_integer(),
      Imp_Srd3 = col_character(),
      
      Imp_Key  = col_character(),
      Imp_Scr  = col_character(),
      
      Imp_Cnt  = col_integer(),
      aln_key  = col_character()
    )
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    cmd_str = glue::glue("gzip -dc {ref} | join -t $'\t' -1{idxA} -2{idxB} - {can} | gzip -c - > {out}")
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]: Running cmd={cmd_str}...{RET}"))
    cmd_ret <- system(cmd_str)
    
    if (opt$verbose>=1)
      cat(glue::glue("[{funcTag}]: Loading intersection output={out}...{RET}"))
    
    ret_tib <- suppressMessages(suppressWarnings( 
      readr::read_tsv(out, col_names=names(int_seq_cols$cols), col_types=int_seq_cols) )) %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor), as.character) )
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

intersect_GRS = function(can,ref, can_key="unq_can_key", ref_prefix=NULL, # ref_prefix="imp",
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'intersect_GRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    map_tib <- 
      GenomicRanges::findOverlaps(can,ref, ignore.strand=TRUE) %>%
      as.data.frame() %>% tibble::as_tibble()
    map_cnt <- print_tib(map_tib,funcTag, verbose,vt+4,tc, n="map")
    
    can_tib <- can %>% as.data.frame() %>%
      rownames_to_column(var=can_key) %>% 
      tibble::as_tibble() %>%
      dplyr::select(dplyr::all_of(can_key))
    can_cnt <- print_tib(can_tib,funcTag, verbose,vt+4,tc, n="can")
    
    ref_tib <- ref %>% as.data.frame() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(seqnames=as.character(seqnames),
                    strand=as.character(strand)) %>%
      dplyr::rename(chr=seqnames,
                    pos=start,
                    top_srd=strand)
    
    if (!is.null(ref_prefix)) ref_tib <- ref_tib %>% 
      purrr::set_names(paste(ref_prefix,names(.), sep="_"))
    
    ref_cnt <- print_tib(ref_tib,funcTag, verbose,vt+4,tc, n="ref")
    
    # Column Bind Data Sets::
    ret_tib <- dplyr::bind_cols(
      can_tib[map_tib$queryHits, ],
      ref_tib[map_tib$subjectHits,]
    )
    
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
#                       Address To Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

mutate_probe_id = function(tib, 
                           pid="Probe_ID", cgn="Imp_Cgn_Seq",
                           des="Ord_Des",  din="Ord_Din",
                           tb="Imp_TB_Seq", co="Imp_CO_Seq",
                           inf="Infinium_Design",pad=8,
                           verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutate_probe_id'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    pid_sym <- rlang::sym(pid)
    cgn_sym <- rlang::sym(cgn)
    des_sym <- rlang::sym(des)
    din_sym <- rlang::sym(din)
    tb_sym  <- rlang::sym(tb)
    co_sym  <- rlang::sym(co)
    inf_sym <- rlang::sym(inf)
    
    ret_tib <- tib %>% 
      dplyr::mutate(
        !!inf_sym:=dplyr::case_when(
          !!des_sym=="U" | !!des_sym=="M" ~ 1, 
          !!des_sym=="2" ~ 2, 
          TRUE ~ NA_real_) %>% as.integer(), 
        !!pid_sym:=paste(paste0(!!din_sym,stringr::str_pad(!!cgn_sym,width=pad, side="left", pad="0")), 
                         paste0(!!tb_sym,!!co_sym,!!inf_sym), sep="_")
      )
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

add_comb = function(tibA, tibB, field,
                    join,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'add_comb'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- 
      dplyr::inner_join(
        tibA, tibB, 
        by=dplyr::all_of(join),
        suffix=c("_U","_M")
      ) # %>%
    # dplyr::select(dplyr::starts_with(field)) %>%
    # dplyr::distinct() %>%
    # purrr::set_names(c("U","M"))
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

man_to_add = function(tib, 
                      pid="IlmnID", des="Man_Des", din="Man_Din", inf="Man_Inf",
                      addA="AddressA_ID", prbA="AlleleA_ProbeSeq", 
                      addB="AddressB_ID", prbB="AlleleB_ProbeSeq",
                      del="_",
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'man_to_add'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    pid_sym  = rlang::sym(pid)
    des_sym  = rlang::sym(des)
    din_sym  = rlang::sym(din)
    inf_sym  = rlang::sym(inf)
    
    addA_sym = rlang::sym(addA)
    prbA_sym = rlang::sym(prbA)
    
    addB_sym = rlang::sym(addB)
    prbB_sym = rlang::sym(prbB)
    
    inf_list <- tib %>%
      dplyr::mutate(
        !!din_sym:=stringr::str_sub(!!pid_sym, 1,2),
        !!inf_sym:=dplyr::case_when(
          is.na(!!addA_sym) | is.na(!!prbA_sym) ~ NA_real_,
          is.na(!!addB_sym) | is.na(!!prbB_sym) ~ 2,
          !is.na(!!addA_sym) & ! is.na(!!prbA_sym) &
            !is.na(!!addB_sym) & ! is.na(!!prbB_sym) ~ 1,
          TRUE ~ NA_real_) %>% as.integer()
      ) %>% split(.[[inf]])
    
    ret_tib <- dplyr::bind_rows(
      inf_list[["1"]] %>% 
        dplyr::mutate(!!des_sym:=as.character("U")) %>%
        dplyr::rename(Man_Add=!!addA_sym,
                      Man_Prb=!!prbA_sym,
                      Man_PID=!!pid_sym) %>%
        dplyr::select(Man_Add,!!des,Man_Din,!!inf,Man_Prb,Man_PID),
      
      inf_list[["1"]] %>% 
        dplyr::mutate(!!des_sym:=as.character("M")) %>%
        dplyr::rename(Man_Add=!!addB_sym,
                      Man_Prb=!!prbB_sym,
                      Man_PID=!!pid_sym) %>%
        dplyr::select(Man_Add,!!des,Man_Din,!!inf,Man_Prb,Man_PID),
      
      inf_list[["2"]] %>% 
        dplyr::mutate(!!des_sym:=as.character("2")) %>%
        dplyr::rename(Man_Add=!!addA_sym,
                      Man_Prb=!!prbA_sym,
                      Man_PID=!!pid_sym) %>%
        dplyr::select(Man_Add,!!des,Man_Din,!!inf,Man_Prb,Man_PID)
    ) %>% dplyr::arrange(Man_Add)
    
    ret_tib <- ret_tib %>%
      clean_tibble(verbose=verbose,tt=tt)
    
    # ret_cnt <- ret_tib %>% base::nrow()
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# col_key='Ord_Col'
# nxb_key='Imp_Nxb_Bsp_U'
add_to_man = function(tib, join, runName,
                      des_key="Ord_Des", pid_key="Ord_Key",
                      rep_key=NULL, rep_val=NULL,
                      col_key=NULL, nxb_key=NULL,
                      csv=NULL, validate=TRUE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'add_to_man'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    des_list <- NULL
    des_list <- tib %>% split(.[[des_key]])
    des_cnt  <- des_list %>% names() %>% length()
    
    # Build Infinium I::
    ret1_tib <- NULL
    if (!is.null(des_list[["U"]]) && !is.null(des_list[["M"]])) {
      ret1_tib <- dplyr::full_join(
        dplyr::select(des_list[["U"]], -dplyr::all_of(des_key)), 
        dplyr::select(des_list[["M"]], -dplyr::all_of(des_key)), 
        by=dplyr::all_of(join),
        suffix=c("_U","_M")
      ) %>% 
        # TBD:: We should allow these "singletons" to pass, but under
        #  a different classification...
        #
        dplyr::filter(!is.na(Address_U) & !is.na(Address_M))
      
      ret1_cnt <- print_tib(ret1_tib,funcTag, verbose,vt+4,tc, n="InfI")
    }
    
    # Build Infinium II::
    ret2_tib <- NULL
    if (!is.null(des_list[["2"]])) {
      ret2_tib <- dplyr::bind_cols(
        dplyr::select(des_list[["2"]],  dplyr::all_of(join)),
        dplyr::select(des_list[["2"]], -dplyr::all_of(join)) %>% 
          purrr::set_names(paste(names(.),"U", sep="_"))
      )
      ret2_cnt <- print_tib(ret2_tib,funcTag, verbose,vt+4,tc, n="InfII")
    }
    
    # Bind Infinium I/II into single manifest::
    ret_tib <- dplyr::bind_rows(ret1_tib, ret2_tib) %>%
      dplyr::mutate(
        Infinium_Design=dplyr::case_when(
          is.na(Address_M) ~ 2,
          !is.na(Address_M) ~ 1,
          TRUE ~ NA_real_
        ) %>% as.integer(),
        Infinium_Design_Type=dplyr::case_when(
          Infinium_Design==1 ~ 'I',
          Infinium_Design==2 ~ 'II',
          TRUE ~ NA_character_)
      )
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="InfI+II.0")
    
    if (!is.null(col_key)) {
      col_sym <- rlang::sym(col_key)
      ret_tib <- ret_tib %>%
        dplyr::mutate(
          col=dplyr::case_when(
            Infinium_Design==1 ~ !!col_sym,
            TRUE ~ NA_character_),
          Color_Channel=dplyr::case_when(
            Infinium_Design==2 ~ 'Both',
            col=="R" ~ 'Red',
            col=='G' ~ 'Grn',
            TRUE ~ NA_character_),
          Next_Base=dplyr::case_when(
            col=="R" ~ "A",
            col=="G" ~ "C",
            TRUE ~ NA_character_
          ),
          Probe_Source=runName
        )
    } else if (!is.null(nxb_key)) {
      nxb_sym <- rlang::sym(nxb_key)
      
      ret_tib <- ret_tib %>%
        dplyr::mutate(
          Color_Channel=dplyr::case_when(
            Infinium_Design==2 ~ 'Both',
            Infinium_Design==1 & (!!nxb_sym=='A' | !!nxb_sym=='T') ~ 'Red',
            Infinium_Design==1 & (!!nxb_sym=='C' | !!nxb_sym=='G') ~ 'Grn',
            TRUE ~ NA_character_),
          col=dplyr::case_when(
            Infinium_Design==2 ~ NA_character_,
            Infinium_Design==1 ~ stringr::str_sub(Color_Channel, 1,1),
            TRUE ~ NA_character_),
          Next_Base=dplyr::case_when(
            col=="R" ~ !!nxb_sym,
            col=="G" ~ !!nxb_sym,
            TRUE ~ NA_character_
          ),
          # Next_Base=!!nxb_sym,
          Probe_Source=runName
        )
    } else {
      # Do nothing for now, but this should not be a standard use case
      return(NULL)
    }
    
    if (!is.null(pid_key)) {
      pid_key_sym = rlang::sym(pid_key)
      ret_tib <- ret_tib %>% dplyr::arrange(pid_key)
    }
    
    if (!is.null(pid_key) && !is.null(rep_key) && !is.null(rep_val)) {
      rep_key_sym = rlang::sym(rep_key)
      rep_val_sym = rlang::sym(rep_val)
      
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Adding Probe Replicate({rep_key}/{rep_val})...{RET}"))
      ret_tib <- ret_tib %>% dplyr::group_by(!!rep_key_sym) %>%
        dplyr::mutate(!!rep_val_sym := dplyr::row_number(),
                      !!pid_key_sym := paste0(!!rep_key_sym,!!rep_val_sym)) %>%
        dplyr::ungroup()
    }
    
    if (!is.null(csv)) {
      outDir <- base::dirname(csv)
      if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
      safe_write(ret_tib,"csv",csv,
                 funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="InfI+II")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             FASTA File Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

add_to_fas = function(tib, prb_key="Prb_Seq", 
                      add_key="Address", des_key="Prb_Des", type_key="prb_type",
                      prb_fas=NULL,dat_csv=NULL,
                      u49_tsv=NULL,m49_tsv=NULL,
                      del="_",
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'add_to_fas'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; key={prb_key}...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    prb_sym <- rlang::sym(prb_key)
    des_sym <- rlang::sym(des_key)
    aln_vec <- c(add_key,des_key,type_key)
    
    # Build Alignment Keys/Seqs
    ret_tib <- tib %>% # tail(n=200) %>%
      tidyr::unite(Aln_Key, dplyr::all_of(aln_vec),sep=del,remove=FALSE) %>%
      dplyr::mutate(
        Aln_Prb=deMs(!!prb_sym, uc=TRUE),
        Aln_Rev=revCmp(Aln_Prb),
        Aln_P49=dplyr::case_when(
          !!des_sym == '2' ~ stringr::str_sub(Aln_Prb, 2),
          !!des_sym == 'U' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          !!des_sym == 'M' ~ stringr::str_remove(Aln_Prb, '[A-Z]$'),
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::distinct(Aln_Key,Aln_Prb, .keep_all=TRUE) %>%
      clean_tibble()
    # utils::type.convert() %>% 
    # dplyr::mutate(across(where(is.factor), as.character) )
    
    # Generate Remaining Data::
    if (!is.null(prb_fas))
      fas_vec <- ret_tib %>%
      dplyr::mutate(fas_line=paste0(">",Aln_Key,"\n",Aln_Prb) ) %>%
      dplyr::pull(fas_line)
    
    u49_tib <- NULL
    if (!is.null(u49_tsv))
      u49_tib <- ret_tib %>% 
      dplyr::filter(!!des_sym == '2' | !!des_sym=='U') %>%
      dplyr::filter(!is.na(Aln_P49)) %>%
      dplyr::select(Aln_P49,Aln_Key) %>%
      dplyr::arrange(Aln_P49)
    
    m49_tib <- NULL
    if (!is.null(m49_tsv))
      m49_tib <- ret_tib %>% 
      dplyr::filter(!!des_sym=='M') %>%
      dplyr::filter(!is.na(Aln_P49)) %>%
      dplyr::select(Aln_P49,Aln_Key) %>%
      dplyr::arrange(Aln_P49)
    
    # Safe Write Outputs::
    safe_write(ret_tib,"csv",dat_csv,  
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    safe_write(fas_vec,"line",prb_fas, 
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    safe_write(u49_tib,"tsv",u49_tsv,cols=FALSE, 
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    safe_write(m49_tib,"tsv",m49_tsv,cols=FALSE, 
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
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
#                        Bowtie Alignment Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_bsmap = function(exe, fas, gen, bsp, 
                     opt=NULL, lan=NULL, run=FALSE,
                     verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'run_bsmap'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; bsp={bsp}.{RET}"))
  
  aln_bsp <- bsp
  
  aln_ssh <- NULL
  ret_tsv <- NULL
  sys_ret <- "Not_Run"
  stime <- system.time({
    if (is.null(opt) || length(opt)==0)
      opt="-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R"
    
    build_file_dir(aln_bsp)
    aln_ssh <- paste(aln_bsp, 'run.sh', sep='.')
    ret_tsv <- paste(aln_bsp, 'tsv.gz', sep='.')
    
    aln_cmd <- glue::glue("{exe} -a {fas} -d {gen} {opt} -o {aln_bsp}{RET}",
                          "cat {aln_bsp} | ",
                          # "cut -f 1,2,4-11 | ",
                          "gzip -c - > {ret_tsv}{RET}",
                          "rm -f {aln_bsp}{RET}")
    
    safe_write(aln_cmd,"line",aln_ssh,
               funcTag=funcTag,verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    Sys.chmod(paths=aln_ssh, mode="0777")
    
    if (run) {
      fas_name <- fas %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
      gen_name <- gen %>% base::basename() %>% stringr::str_remove('.[A-Za-z]+$') %>% stringr::str_remove('.[A-Za-z]+$')
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Will Execute bsmap: {fas_name} vs. {gen_name}...{RET}"))
      
      exe_method="Running"
      if (!is.null(lan)) {
        exe_method="Launching"
        aln_ssh <- paste(lan,aln_ssh, sep=' ')
      }
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} {exe_method}; CMD={aln_cmd}...{RET}"))
      
      sys_ret <- base::system(aln_ssh)
      
      if (sys_ret!=0) {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: sys_ret({sys_ret}) != 0!!!{RET}{RET}"))
        return(NULL)
      }
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; sys={sys_ret}; elapsed={etime}.{RET}{RET}"))
  
  if (run) return(ret_tsv)
  aln_ssh
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             BSMAP Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_bsmap = function(bsp, sort=FALSE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_bsmap'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    bsp_cols <- cols(
      Bsp_Key  = col_character(),
      Bsp_Seq  = col_character(),
      Bsp_Qual = col_character(),
      Bsp_Tag  = col_character(),
      
      Bsp_Chr  = col_character(),
      Bsp_Beg  = col_integer(),
      Bsp_Srd  = col_character(),
      Bsp_Mis  = col_integer(),
      
      Bsp_Ref  = col_character(),
      Bsp_Gap  = col_integer(),
      
      Bsp_Str  = col_character()
    )
    
    # Load BSP
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading BSP={bsp}...{RET}"))
    
    ret_tib <- 
      readr::read_tsv(bsp,col_names=names(bsp_cols$cols),col_types=bsp_cols) %>%
      dplyr::select(-Bsp_Qual) %>%
      dplyr::mutate(Bsp_Chr=paste0('chr',stringr::str_remove(Bsp_Chr,'^chr')))
    
    # Sort by genomic position::
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(Bsp_Chr, Bsp_Beg)
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

join_seq_intersect = function(u49,m49,bed=NULL,org=NULL,
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'join_seq_intersect'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    imp_col_vec <- c("Address","Ord_Des","Ord_Din",
                     "Imp_Chr","Imp_Pos","Imp_Cgn",
                     "Imp_FR","Imp_TB","Imp_CO","Imp_Nxb",
                     "Bsp_Din_Ref","Bsp_Din_Scr","Aln_Prb")
    
    ret_tib <- 
      dplyr::bind_rows(u49,m49) %>%
      dplyr::select(-Imp_SrdI,-Imp_Scr) %>% 
      dplyr::mutate(Imp_Key=stringr::str_split(Imp_Key, pattern=",") ) %>% 
      tidyr::unnest(Imp_Key) %>%
      tidyr::separate(
        Imp_Key, 
        into=c("Imp_Cgn","Imp_Hit_hg38","Imp_Hit_hg37", "Imp_Hit_hg36", "Imp_Hit_mm10"), 
        sep="_", remove=TRUE) %>%
      tidyr::separate(
        aln_key, 
        into=c("Address", "Ord_Des", "Ord_Din"), 
        sep="_") %>%
      dplyr::rename(Aln_Prb=Imp_Seq, Aln_Nuc=Imp_Nuc) %>%
      tidyr::separate(Imp_Srd3, into=c("Imp_TB","Imp_CO", "Imp_Nxb"), 
                      sep=c(1,2)) %>%
      dplyr::select(Address, Ord_Des, Ord_Din, Imp_Cgn, 
                    Imp_TB, Imp_CO, Imp_Nxb, Aln_Prb, Aln_Nuc, 
                    dplyr::everything()) %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor), as.character) )
    
    if (!is.null(bed)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Adding cgn bed...{RET}"))
      
      ret_tib <- bed %>%
        dplyr::right_join(ret_tib,by=c("Imp_Cgn")) %>%
        dplyr::mutate(
          Imp_FR=dplyr::case_when(
            Imp_Top_Srd=="+" & Imp_TB=="T" ~ "F",
            Imp_Top_Srd=="-" & Imp_TB=="T" ~ "R",
            Imp_Top_Srd=="+" & Imp_TB=="B" ~ "R",
            Imp_Top_Srd=="-" & Imp_TB=="B" ~ "F",
            TRUE ~ NA_character_
          )
        )
    }
    
    if (!is.null(org)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Adding org bed...{RET}"))
      
      ret_tib <- ret_tib %>%
        dplyr::left_join(
          org, 
          by=c("Aln_Prb"="Org_49P",
               "Ord_Des"="Org_Des",
               "Imp_Chr"="Org_Chr",
               "Imp_Pos"="Org_Pos",
               "Imp_FR"="Org_FR",
               "Imp_TB"="Org_TB",
               "Imp_CO"="Org_CO")
        )
    }
    
    ret_tib <- ret_tib  %>%
      dplyr::select(dplyr::any_of(imp_col_vec),dplyr::everything()) %>%
      clean_tibble()
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

join_bsmap = function(add, bsp=NULL, bed=NULL, org=NULL, file=NULL, 
                      join_key, join_type="inner", 
                      prb_des_key="Ord_Des", prb_din_key="Ord_Din",
                      fields=NULL, sort=FALSE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'join_bsmap'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}",
                   "[{funcTag}]:{tabsStr}  join_key={join_key}{RET}",
                   "[{funcTag}]:{tabsStr} join_type={join_type}{RET}")
    )
  
  if (is.null(bsp) && is.null(file)) {
    stop(glue::glue("[{funcTag}]:{tabsStr} ERROR: bsp AND file are null!!!{RET}{RET}"))
    return(NULL)
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    imp_col_vec <- c("Address","Ord_Des","Ord_Din",
                     "Imp_Chr","Imp_Pos","Imp_Cgn",
                     "Imp_FR","Imp_TB","Imp_CO","Imp_Nxb",
                     "Bsp_Din_Ref","Bsp_Din_Scr","Aln_Prb")
    
    # Load from file if provided::
    if (!is.null(file))
      bsp <- load_bsmap(bsp=file, sort=sort,
                        verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    if (is.null(bsp)) {
      stop(glue::glue("[{funcTag}]:{tabsStr} ERROR: bsp is null!!!{RET}{RET}"))
      return(ret_tib)
    }
    print_tib(add,funcTag, verbose,vt+4,tc, n="address_tib")
    
    prb_des_sym  <- rlang::sym(prb_des_key)
    prb_din_sym  <- rlang::sym(prb_din_key)
    bsp_join_key <- names(bsp)[1]
    join_key_sym <- rlang::sym(join_key)
    bsp_join_sym <- rlang::sym( bsp_join_key )
    
    # Rename Bsp_Key to join key and join::
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Joining source data...{RET}",
                     "[{funcTag}]:{tabsStr} Join Keys::{RET}",
                     "[{funcTag}]:{tabsStr} add_join_key={join_key}{RET}",
                     "[{funcTag}]:{tabsStr} bsp_join_key={bsp_join_key}{RET}{RET}" )
      )
    
    if (join_type=="inner") {
      ret_tib <- bsp %>%
        dplyr::rename(!!join_key_sym := !!bsp_join_sym) %>%
        dplyr::inner_join(add, by=join_key) # %>% head(n=100)
    } else if (join_type=="left") {
      ret_tib <- bsp %>%
        dplyr::rename(!!join_key_sym := !!bsp_join_sym) %>%
        dplyr::left_join(add, by=join_key) # %>% head(n=100)
    } else if (join_type=="right") {
      ret_tib <- bsp %>%
        dplyr::rename(!!join_key_sym := !!bsp_join_sym) %>%
        dplyr::right_join(add, by=join_key) # %>% head(n=100)
    } else if (join_type=="full") {
      ret_tib <- bsp %>%
        dplyr::rename(!!join_key_sym := !!bsp_join_sym) %>%
        dplyr::full_join(add, by=join_key) # %>% head(n=100)
    } else {
      stop(glue::glue("[{funcTag}]:{tabsStr} ERROR: Unsupported join_type={join_type}!!!{RET}{RET}"))
      return(ret_tib)
    }
    
    #
    # Calculate new fields from join::
    #
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Calculating new fields from joined data...{RET}"))
    
    ret_tib <- ret_tib %>%
      dplyr::mutate(
        # Generate bisulfite converted Reference Sequences::
        Bsp_RefU=bscUs(Bsp_Ref,uc=TRUE),
        Bsp_RefM=bscMs(Bsp_Ref,uc=TRUE),
        Bsp_RefD=bscDs(Bsp_Ref,uc=TRUE),
        
        # Add Alignment Orientation wrt Probe Sequence/BSP Seq::
        #  This is used to determine NxB, CpG matching, etc.
        Bsp_Prb_Dir=dplyr::case_when(
          Aln_Prb==Bsp_Seq ~ "f",
          Aln_Rev==Bsp_Seq ~ "r",
          TRUE ~ NA_character_),
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                      Di-Nucleotide & Next Base:: Reference
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        # Extract Potential Next Base from Reference Sequence::
        NB_F1=stringr::str_sub(Bsp_Ref,-2,-2),
        NB_R1=stringr::str_sub(Bsp_Ref, 2, 2),
        NB_F2=stringr::str_sub(Bsp_Ref,-1,-1),
        NB_R2=stringr::str_sub(Bsp_Ref, 1, 1),
        
        # Extract Potential CG di-nucleotide from Reference Sequence::
        CG_F1=stringr::str_sub(Bsp_Ref,-4,-3),
        CG_R1=stringr::str_sub(Bsp_Ref, 3, 4),
        CG_F2=stringr::str_sub(Bsp_Ref,-3,-2),
        CG_R2=stringr::str_sub(Bsp_Ref, 2, 3),
        
        # Set Expected target Next Base from Alignment Orientation::
        Bsp_Nxb_Ref=dplyr::case_when(
          !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F1,
          !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R1,
          
          !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F1,
          !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R1,
          
          !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F2,
          !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R2,
          TRUE ~ NA_character_
        ), # %>% stringr::str_to_upper(),
        
        # Set Expected target CG di-nucleotide from Alignment Orientation::
        Bsp_Din_Ref=dplyr::case_when(
          !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F1,
          !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R1,
          
          !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F1,
          !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R1,
          
          !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F2,
          !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R2,
          TRUE ~ NA_character_
        ) %>% stringr::str_to_upper(),
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #             Di-Nucleotide & Next Base:: Bisulfite Converted
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        # Extract Potential Next Base from Reference Sequence::
        NB_F1M=stringr::str_sub(Bsp_RefM,-2,-2),
        NB_R1M=stringr::str_sub(Bsp_RefM, 2, 2),
        
        NB_F1U=stringr::str_sub(Bsp_RefU,-2,-2),
        NB_R1U=stringr::str_sub(Bsp_RefU, 2, 2),
        
        NB_F2D=stringr::str_sub(Bsp_RefD,-1,-1),
        NB_R2D=stringr::str_sub(Bsp_RefD, 1, 1),
        
        # Extract Potential CG di-nucleotide from Reference Sequence::
        CG_F1M=stringr::str_sub(Bsp_RefM,-4,-3),
        CG_R1M=stringr::str_sub(Bsp_RefM, 3, 4),
        
        CG_F1U=stringr::str_sub(Bsp_RefU,-4,-3),
        CG_R1U=stringr::str_sub(Bsp_RefU, 3, 4),
        
        CG_F2D=stringr::str_sub(Bsp_RefD,-3,-2),
        CG_R2D=stringr::str_sub(Bsp_RefD, 2, 3),
        
        # Set Expected target Next Base from Alignment Orientation::
        Bsp_Nxb_Bsc=dplyr::case_when(
          !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F1M,
          !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R1M,
          
          !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F1U,
          !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R1U,
          
          !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ NB_F2D,
          !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ NB_R2D,
          TRUE ~ NA_character_
        ), # %>% stringr::str_to_upper(),
        
        # Set Expected target CG di-nucleotide from Alignment Orientation::
        Bsp_Din_Bsc=dplyr::case_when(
          !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F1M,
          !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R1M,
          
          !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F1U,
          !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R1U,
          
          !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ CG_F2D,
          !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ CG_R2D,
          TRUE ~ NA_character_
        ) %>% stringr::str_to_upper(),
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                      Correct CG# Alignment Position::
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        #
        # Update Correct Genomic CG# Location based on alignment orientation::
        #
        # TBD:: Need to correct for non-CpG sites::
        #       - rs:: [+-] => +0
        #              [--] => +1
        #
        Bsp_Pos=dplyr::case_when(
          # prb_din_sym=="cg":: Original Code for only CpG! TO BE DELETED...
          #
          # !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          # !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          # 
          # !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          # !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          # 
          # !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          # !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          #
          # prb_din_sym=="cg"::
          #
          !!prb_din_sym=="cg" & !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          !!prb_din_sym=="cg" & !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="cg" & !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +48,
          !!prb_din_sym=="cg" & !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="cg" & !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="cg" & !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          #
          # prb_din_sym=="ch":: for now assume ch==rs
          #
          !!prb_din_sym=="ch" & !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="ch" & !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="ch" & !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="ch" & !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="ch" & !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          !!prb_din_sym=="ch" & !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          #
          # prb_din_sym=="rs"::
          #
          # !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          # !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          # 
          # !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          # !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          # 
          # !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="--" | Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          # !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="+-" | Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          #
          # prb_din_sym=="rs":: --/+-
          #
          !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="--") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="+-") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="--") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="+-") & Bsp_Prb_Dir=="r" ~ Bsp_Beg + 0,
          
          !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="--") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +50,
          !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="+-") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          # 
          # prb_din_sym=="rs":: --/+-
          #
          !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="rs" & !!prb_des_sym=="M" & (Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="rs" & !!prb_des_sym=="U" & (Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 1,
          
          !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="++") & Bsp_Prb_Dir=="f" ~ Bsp_Beg +49,
          !!prb_din_sym=="rs" & !!prb_des_sym=="2" & (Bsp_Srd=="-+") & Bsp_Prb_Dir=="r" ~ Bsp_Beg - 0,
          
          TRUE ~ NA_real_
        ) %>% as.integer(),
        Bsp_Din_Scr=dplyr::case_when(
          Bsp_Srd=="--" & Bsp_Din_Ref=="CG" ~ 0,
          Bsp_Srd=="+-" & Bsp_Din_Ref=="CG" ~ 1,
          Bsp_Srd=="-+" & stringr::str_starts(Bsp_Din_Ref,"G") ~ 2,
          Bsp_Srd=="++" & stringr::str_ends(Bsp_Din_Ref,"C") ~ 3,
          TRUE ~ 9) %>% as.integer()
      ) %>% 
      # Create Unique Aln Key for multiple hits::
      #
      dplyr::group_by(!!join_key_sym) %>%
      dplyr::mutate(Aln_Key_Unq=paste(!!join_key_sym, dplyr::row_number(), sep="_")) %>%
      dplyr::ungroup() %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    # Sort by genomic position::
    if (sort) ret_tib <- ret_tib %>% dplyr::arrange(Bsp_Chr, Bsp_Beg)
    
    if (!is.null(bed)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Adding cgn bed...{RET}"))
      
      ret_tib <- bed %>% 
        dplyr::right_join(ret_tib, by=c("Imp_Chr"="Bsp_Chr", "Imp_Pos"="Bsp_Pos")) %>%
        dplyr::mutate(
          Imp_FR=dplyr::case_when(
            Bsp_Srd=="+-" ~ "R",
            Bsp_Srd=="--" ~ "F",
            Bsp_Srd=="++" ~ "R",
            Bsp_Srd=="-+" ~ "F",
            TRUE ~ NA_character_
          ),
          Imp_CO=dplyr::case_when(
            Bsp_Srd=="+-" ~ "C",
            Bsp_Srd=="--" ~ "C",
            Bsp_Srd=="++" ~ "O",
            Bsp_Srd=="-+" ~ "O",
            TRUE ~ NA_character_
          ),
          Imp_TB=dplyr::case_when(
            Imp_Top_Srd=="+" & Imp_FR=="F" ~ "T",
            TRUE ~ "B"
          ),
          Imp_Nxb=stringr::str_to_upper(Bsp_Nxb_Ref)
        )
    }
    
    if (!is.null(org)) {
      if (verbose>=vt)
        cat(glue::glue("[{funcTag}]:{tabsStr} Adding org bed...{RET}"))
      
      ret_tib <- ret_tib %>%
        dplyr::left_join(
          org, 
          by=c("Aln_P49"="Org_49P",
               "Ord_Des"="Org_Des",
               "Imp_Chr"="Org_Chr",
               "Imp_Pos"="Org_Pos",
               "Imp_FR"="Org_FR",
               "Imp_TB"="Org_TB",
               "Imp_CO"="Org_CO")
        )
    }
    
    ret_tib <- ret_tib  %>%
      dplyr::select(dplyr::any_of(imp_col_vec),
                    dplyr::everything()) %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor), as.character) )
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

bsp_to_grs = function(bsp, rds=NULL,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'bsp_to_grs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (!is.null(rds)) {
      dir <- base::dirname(rds)
      if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    }
    
    #
    # TBD:: Temp Fix for Bsp_Pos == NA !!!
    #
    
    bsp <- bsp %>% dplyr::filter(!is.na(Bsp_Pos))
    
    ret_grs <- 
      GRanges(
        seqnames = Rle(paste0("chr",bsp$Bsp_Chr)),
        strand=Rle(stringr::str_sub( bsp$Bsp_Srd, 1,1 ) ),
        
        # Ord_Id  = bsp$Ord_Id,
        
        # Bsp_Add = bsp$Bsp_Add,
        # Bsp_Add = bsp$Address,
        
        # Bsp_Srd = bsp$Bsp_Srd,
        # Bsp_Src = bsp$Bsp_Src,
        
        Prb_Cgn = bsp$Prb_Cgn,
        Prb_Des = bsp$Prb_Des,
        
        Prb_Ord_Seq  = bsp$Prb_Seq,
        Prb_Aln_50U  = bsp$Prb_Seq %>%
          stringr::str_replace_all("R","A") %>% # A/G
          stringr::str_replace_all("Y","T") %>% # C/T
          
          stringr::str_replace_all("S","C") %>% # G/C
          stringr::str_replace_all("W","A") %>% # A/T
          stringr::str_replace_all("K","T") %>% # G/T
          stringr::str_replace_all("M","A") %>% # A/C
          
          stringr::str_replace_all("B","T") %>% # C/G/T
          stringr::str_replace_all("D","A") %>% # A/G/T
          stringr::str_replace_all("H","A") %>% # A/C/T
          stringr::str_replace_all("V","A") %>% # A/C/G
          
          stringr::str_replace_all("N","A"), # A/C/T/G
        
        Bsp_Ref_Seq = bsp$Bsp_Ref,
        Bsp_Din_Scr = bsp$Bsp_Din_Scr,
        Bsp_Din_Ref = bsp$Bsp_Din_Ref,
        Bsp_Nxb_Ref = bsp$Bsp_Nxb_Ref,
        Bsp_Din_Bsc = bsp$Bsp_Din_Bsc,
        Bsp_Nxb_Bsc = bsp$Bsp_Nxb_Bsc,
        
        Bsp_Tag = bsp$Bsp_Tag,
        Bsp_Srd = bsp$Bsp_Srd,
        Bsp_Mis = bsp$Bsp_Mis,
        Bsp_Gap = bsp$Bsp_Gap,
        
        Bsp_Hit0 = bsp$Bsp_Hit0,
        Bsp_Hit1 = bsp$Bsp_Hit1,
        Bsp_Hit2 = bsp$Bsp_Hit2,
        Bsp_Hit3 = bsp$Bsp_Hit3,
        Bsp_Hit4 = bsp$Bsp_Hit4,
        Bsp_Hit5 = bsp$Bsp_Hit5,
        
        IRanges(start=bsp$Bsp_Pos,
                end=bsp$Bsp_Pos+1,
                names=bsp$Unique_ID
        )
      )
    
    if (!is.null(rds)) {
      if (verbose>=1)
        cat(glue::glue("[{funcTag}]: Writing Probe GRS RDS={rds}...{RET}"))
      readr::write_rds(ret_grs, rds, compress="gz")
    }
    ret_cnt <- bsp %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_grs
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Manifest File Methods:: Formatting
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

format_ORD = function(tib, idx_key="IDX", uniq=TRUE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_ORD'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    idx_sym <- rlang::sym(idx_key)
    add_col  <- c(idx_key,"Ord_Id","Prb_Des","Prb_Seq","Prb_Par","Prb_Col")
    sel_colA <- c(idx_key,"Assay_Design_Id","DesA","PrbA","PrbB","PrbCol")
    sel_colB <- c(idx_key,"Assay_Design_Id","DesB","PrbB","PrbA","PrbCol")
    
    tmp_tib <- tib %>%
      dplyr::select(-AlleleA_Probe_Id, -AlleleB_Probe_Id) %>%
      dplyr::rename(PrbA=AlleleA_Probe_Sequence,
                    PrbB=AlleleB_Probe_Sequence) %>%
      dplyr::mutate(
        PrbA=stringr::str_to_upper(PrbA),
        PrbB=stringr::str_to_upper(PrbB),
        DesA=dplyr::case_when(
          !is.na(PrbA) & !is.na(PrbB) ~ "U",
          !is.na(PrbA) &  is.na(PrbB) ~ "2",
          TRUE ~ NA_character_
        ),
        DesB=dplyr::case_when(
          !is.na(PrbA) & !is.na(PrbB) ~ "M",
          TRUE ~ NA_character_
        ),
        PrbCol=dplyr::case_when(
          DesA=="2" ~ NA_character_,
          DesA=="U" & DesB=="M" & Normalization_Bin=="A" ~ "R",
          DesA=="U" & DesB=="M" & Normalization_Bin=="B" ~ "G",
          TRUE ~ NA_character_
        )
      )
    tmp_cnt <- print_tib(tmp_tib,funcTag,verbose,vt+6,tc, n="tmp")
    
    ret_tib <- dplyr::bind_rows(
      tmp_tib %>% dplyr::filter(DesB=="M") %>% 
        dplyr::select(dplyr::all_of(sel_colB)) %>%
        purrr::set_names(add_col),
      tmp_tib %>% dplyr::filter(DesA=="U") %>% 
        dplyr::select(dplyr::all_of(sel_colA)) %>%
        purrr::set_names(add_col),
      tmp_tib %>% dplyr::filter(DesA=="2") %>% 
        dplyr::select(dplyr::all_of(sel_colA)) %>%
        purrr::set_names(add_col)
    ) %>%
      dplyr::arrange(plyr::desc(!!idx_sym))
    
    ret1_cnt <- print_tib(ret_tib,funcTag,verbose,vt+6,tc, n="ret1")
    if (uniq) ret_tib <- ret_tib %>% dplyr::distinct(prb_seq, .keep_all=TRUE)
    unq1_cnt <- print_tib(ret_tib,funcTag,verbose,vt+6,tc, n="retUnq")
    
    ret_tib <- ret_tib %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    ret_cnt <- print_tib(ret_tib,funcTag,verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

format_MAT = function(tib, idx_key="IDX", uniq=TRUE, trim=TRUE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_MAT'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    idx_sym <- rlang::sym(idx_key)
    
    ret_tib <- tib
    if (c("address_names") %in% names(tib)) {
      
      if (trim) {
        trim_cnt <- ret_tib %>% 
          dplyr::filter(!stringr::str_starts(address_names, '1')) %>% base::nrow()
        if (trim_cnt==0) {
          ret_tib <- ret_tib %>% 
            dplyr::mutate(address_names=as.numeric(
              stringr::str_remove(stringr::str_remove(address_names, '^1'), '^0+')) )
        } else {
          stop(glue::glue("{RET}[{funcTag}]: ERROR: Attempting to trim new tango fomrat, ",
                          "but format doesn't match; trim_cnt={trim_cnt}!!!{RET}{RET}"))
          return(NULL)
        }
      }
      
      ret_tib <- ret_tib %>% 
        dplyr::rename(Address=address_names,
                      Sequence=bo_seq)
      
    } else if (c("Address") %in% names(tib)) {
      if (trim) {
        trim_cnt <- ret_tib %>% 
          dplyr::filter(!stringr::str_starts(Address, '1')) %>% base::nrow()
        if (trim_cnt==0) {
          ret_tib <- ret_tib %>% 
            dplyr::mutate(Address=as.numeric(
              stringr::str_remove(stringr::str_remove(Address, '^1'), '^0+')) )
        } else {
          stop(glue::glue("{RET}[{funcTag}]: ERROR: Attempting to trim new tango fomrat, ",
                          "but format doesn't match; trim_cnt={trim_cnt}!!!{RET}{RET}"))
          return(NULL)
        }
      }
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: Failed to match 'addres_names' or 'Address'!!!{RET}{RET}"))
      return(NULL)
    }
    
    ret_tib <- ret_tib %>% 
      dplyr::mutate(Sequence=stringr::str_to_upper(Sequence),
                    tan_seq=stringr::str_sub(Sequence,1,45),
                    Prb_Seq=stringr::str_sub(Sequence,46)
      ) %>%
      dplyr::filter(!is.na(Address)) %>%
      dplyr::select(Address,!!idx_sym,Prb_Seq, dplyr::everything()) %>% 
      dplyr::arrange(plyr::desc(!!idx_sym))
    
    
    if (uniq) ret_tib <- ret_tib %>% dplyr::distinct(Address, .keep_all=TRUE)
    
    ret_tib <- ret_tib %>% 
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib=({ret_cnt})={RET}"))
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

format_AQP = function(tib, idx_key="IDX", uniq=TRUE,filt=TRUE,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'format_AQP'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    idx_sym <- rlang::sym(idx_key)
    
    # Determine if its AQP & PQC::
    #
    ret_tib <- tib %>% 
      dplyr::mutate(Address=as.numeric(
        stringr::str_remove(stringr::str_remove(Address, '^1'), '^0+')) )
    
    if (c("Decode_Status") %in% names(tib)) {
      ret_tib <- ret_tib %>% 
        dplyr::select(Address,Decode_Status,!!idx_sym) %>%
        dplyr::rename(aqp_val=Decode_Status) %>%
        dplyr::mutate(aqp_src="AQP")
    } else if (c("Status") %in% names(tib)) {
      ret_tib <- ret_tib %>% 
        dplyr::select(Address,Status,!!idx_sym) %>%
        dplyr::rename(aqp_val=Status) %>%
        dplyr::mutate(aqp_src="PQC")
    } else {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Unable to match AQP/PQC!{RET}{RET}"))
      return(NULL)
    } %>% dplyr::arrange(plyr::desc(!!idx_sym))
    
    pre_cnt <- ret_tib %>% base::nrow()
    ret_tib <- ret_tib %>% 
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    if (filt) ret_tib <- ret_tib %>% dplyr::filter(aqp_val==0)
    if (uniq) ret_tib <- ret_tib %>% dplyr::distinct(Address, .keep_all=TRUE)
    
    ret_cnt <- ret_tib %>% base::nrow()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib=({ret_cnt}/{pre_cnt})={RET}"))
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
#                     Manifest File Methods:: Loading
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

guess_man_file = function(file, 
                          fields=c("Assay_Design_Id","Plate","address_names","Address"), 
                          cols=NULL,
                          n_max=100, # guess=100000,
                          verbose=0,vt=6,tc=1,tt=NULL) {
  funcTag <- 'guess_man_file'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  ret_val <- NULL
  stime <- system.time({
    
    file_del_str <- guess_file_del(file, n_max=n_max, verbose=verbose,vt=vt+4)
    if (is.null(file_del_str)) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: file_del_str=NULL!!!{RET}{RET}"))
      return(NULL)
    }
    
    dat_tib <- readr::read_lines(file, n_max=n_max) %>%
      tibble::as_tibble() %>% 
      dplyr::mutate(row_num=dplyr::row_number())
    
    ret_tib <- NULL
    for (field in fields) {
      if (verbose>=vt+6)
        cat(glue::glue("[{funcTag}]:{tabsStr} field={field}{RET}"))
      
      min_tib <- dat_tib %>% dplyr::filter(stringr::str_starts(value, field)) %>% tail(n=1)
      min_cnt <- base::nrow(min_tib)
      
      # Skip Empty Results::
      if (min_cnt==0) next
      
      # Extract Column Count Step::
      col_num <- min_tib$value[1] %>% 
        stringr::str_split(file_del_str, simplify=TRUE) %>% 
        as.vector() %>% length()
      
      # Combine Step::
      ret_tib <- ret_tib %>% dplyr::bind_rows(
        min_tib %>%
          dplyr::mutate(col_num=col_num,del_key=file_del_str,beg_key=field)
      )
      if (verbose>=vt+6)
        cat(glue::glue("[{funcTag}]:{tabsStr} NEXT...{RET}{RET}"))
    }
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

load_man_file = function(file, idx=NULL,
                         n_max=100, guess=1000,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_man_file'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  # Saftey Check for Empty Files::
  if (is.null(file) || length(file)==0)
    return(base::invisible(NULL))
  
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  val_cols <- list()
  val_cols$ord <- 
    cols(
      Assay_Design_Id        = col_character(),
      AlleleA_Probe_Id       = col_character(),
      AlleleA_Probe_Sequence = col_character(),
      AlleleB_Probe_Id       = col_character(),
      AlleleB_Probe_Sequence = col_character(),
      Normalization_Bin      = col_character()
    )
  
  val_cols$mat1 <- 
    cols(
      Plate    = col_character(),
      Row      = col_character(),
      Col      = col_integer(),
      Address  = col_integer(),
      Mod5     = col_character(),
      Sequence = col_character(),
      Mod3     = col_character(),
      Comments = col_character()
    )
  
  val_cols$mat2 <- 
    cols(
      address_names = col_integer(),
      probe_id      = col_character(),
      sequence      = col_character(),
      type_b        = col_character(),
      address_name  = col_integer(),
      bo_seq        = col_character()
    )
  
  val_cols$aqp <- 
    cols(
      Address           = col_integer(),
      Decode_Status     = col_integer(),
      Decode_Error_Code = col_integer(),
      Decode_Score      = col_integer(),
      Func_Status       = col_integer(),
      Func_Error_Code   = col_integer(),
      QC_Action         = col_integer()
    )
  
  val_cols$pqc <- 
    cols(
      Address      = col_integer(),
      Status       = col_integer(),
      Eval_Code    = col_integer(),
      Average_Rep  = col_integer(),
      Expected_Rep = col_integer()
    )
  
  sel_cols <- list()
  sel_cols$ord  <- c("Assay_Design_Id","AlleleA_Probe_Sequence","AlleleB_Probe_Sequence","Normalization_Bin")
  sel_cols$mat1 <- c("Address","Sequence")
  sel_cols$mat2 <- c("address_names","bo_seq")
  # sel_cols$mat2 <- c("address_names","sequence","probe_id")
  sel_cols$aqp  <- c("Address","Decode_Status")
  sel_cols$pqc  <- c("Address","Status")
  
  key_cols <- list()
  key_cols$ord  <- c("Ord_Key","Ord_PrbA","Ord_PrbB","Ord_Norm")
  key_cols$mat1 <- c("Address","Sequence")
  key_cols$mat2 <- c("Address","Sequence")
  key_cols$aqp  <- c("Address","Decode_Status")
  key_cols$pqc  <- c("Address","Decode_Status")
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    guess_tib <- guess_man_file(file, n_max=n_max,
                                verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    row_num=guess_tib$row_num[1]
    col_num=guess_tib$col_num[1]
    del_key=guess_tib$del_key[1]
    beg_key=guess_tib$beg_key[1]
    
    dat_key <- NULL
    val_col <- NULL
    sel_col <- NULL
    key_col <- NULL
    if (is.null(beg_key) || is.null(col_num)) {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Both beg_key AND col_num are NULL!!!{RET}{RET}"))
      return(ret_tib)
    } else if (beg_key==names(val_cols$ord$cols)[1] && col_num==length(val_cols$ord$cols)) {
      dat_key <- "ord"
    } else if (beg_key==names(val_cols$mat1$cols)[1] && col_num==length(val_cols$mat1$cols)) {
      dat_key <- "mat1"
    } else if (beg_key==names(val_cols$mat2$cols)[1] && col_num==length(val_cols$mat2$cols)) {
      dat_key <- "mat2"
    } else if (beg_key==names(val_cols$aqp$cols)[1] && col_num==length(val_cols$aqp$cols)) {
      dat_key <- "aqp"
    } else if (beg_key==names(val_cols$pqc$cols)[1] && col_num==length(val_cols$pqc$cols)) {
      dat_key <- "pqc"
    } else {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Failed to match beg_key({beg_key}) AND ",
                      "col_num({col_num}) to known formats!!!{RET}{RET}"))
      return(ret_tib)
    }
    
    # This sets all the proper valid col types, col selection and col renaming::
    val_col <- val_cols[[dat_key]]  # File Format Column Names/Data-Types
    sel_col <- sel_cols[[dat_key]]  # Selected Column Names
    key_col <- key_cols[[dat_key]]  # Selected Column New Names
    
    # Apply Delimiter Type and Column Names/Data-Types::
    if (del_key==COM) {
      ret_tib <- suppressMessages(suppressWarnings(
        readr::read_csv(file, skip=row_num, guess_max=guess, skip_empty_rows=FALSE,
                        col_names=names(val_col$cols), col_types=val_col)))
      
    } else if (del_key==TAB || del_key==" ") {
      ret_tib <- suppressMessages(suppressWarnings(
        readr::read_tsv(file, skip=row_num, guess_max=guess, skip_empty_rows=FALSE,
                        col_names=names(val_col$cols), col_types=val_col)))
      
    } else {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Unsupported del_key={del_key}!!!{RET}{RET}"))
      return(NULL)
    }
    ret1_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n="ret1")
    
    # Apply Selected Columns::
    ret_tib  <- ret_tib %>% dplyr::select(dplyr::all_of(sel_col))
    ret2_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n="ret2")
    
    # Apply Canonical-Renamed Columns::
    ret_tib  <- ret_tib %>% purrr::set_names(key_col)
    ret3_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n="ret3")
    
    # Apply Formatting Rules:: Order Data Type
    if (dat_key=="ord") {
      ord_tib <- ret_tib %>% 
        dplyr::mutate(
          Ord_PrbA=stringr::str_to_upper(Ord_PrbA),
          Ord_PrbB=stringr::str_to_upper(Ord_PrbB),
          
          Ord_DesA=dplyr::case_when(
            !is.na(Ord_PrbA) & !is.na(Ord_PrbB) ~ "U",
            !is.na(Ord_PrbA) &  is.na(Ord_PrbB) ~ "2",
            TRUE ~ NA_character_
          ),
          Ord_DesB=dplyr::case_when(
            !is.na(Ord_PrbA) & !is.na(Ord_PrbB) ~ "M",
            TRUE ~ NA_character_
          ),
          Ord_Col=dplyr::case_when(
            Ord_DesA=="2" ~ NA_character_,
            Ord_DesA=="U" & Ord_DesB=="M" & Ord_Norm=="A" ~ "R",
            Ord_DesA=="U" & Ord_DesB=="M" & Ord_Norm=="B" ~ "G",
            TRUE ~ NA_character_
          )
        )
      
      ord_cols <- c("Ord_Key","Ord_Des", "Ord_Prb", "Ord_Par", "Ord_Col")
      ord_colA <- c("Ord_Key","Ord_DesA","Ord_PrbA","Ord_PrbB","Ord_Col")
      ord_colB <- c("Ord_Key","Ord_DesB","Ord_PrbB","Ord_PrbA","Ord_Col")
      
      ret_tib <- dplyr::bind_rows(
        ord_tib %>% dplyr::filter(Ord_DesB=="M") %>%
          dplyr::select(dplyr::all_of(ord_colB)) %>%
          purrr::set_names(ord_cols),
        ord_tib %>% dplyr::filter(Ord_DesA=="U") %>%
          dplyr::select(dplyr::all_of(ord_colA)) %>%
          purrr::set_names(ord_cols),
        ord_tib %>% dplyr::filter(Ord_DesA=="2") %>%
          dplyr::select(dplyr::all_of(ord_colA)) %>%
          purrr::set_names(ord_cols)
      ) %>% 
        dplyr::mutate(Ord_Din=stringr::str_sub(Ord_Key, 1,2))
    }
    
    # Apply Formatting Rules:: Match/AQP/PQC Data Types:: Address
    if (dat_key=="mat1" || dat_key=="mat2" ||
        dat_key=="aqp"  || dat_key=="pqc") {
      trim_cnt <- ret_tib %>% 
        dplyr::filter(!stringr::str_starts(Address, '1')) %>% base::nrow()
      if (trim_cnt==0) {
        ret_tib <- ret_tib %>% 
          dplyr::mutate(Address=as.numeric(
            stringr::str_remove(stringr::str_remove(Address, '^1'), '^0+')) )
      } else {
        stop(glue::glue("{RET}[{funcTag}]: ERROR: Attempting to trim new tango fomrat, ",
                        "but format doesn't match; trim_cnt={trim_cnt}!!!{RET}{RET}"))
        return(NULL)
      }
    }
    
    # Apply Formatting Rules:: Match Files:: Sequence
    if (dat_key=="mat1" || dat_key=="mat2") {
      ret_tib <- ret_tib %>% 
        dplyr::mutate(Sequence=stringr::str_to_upper(Sequence),
                      Mat_Tan=stringr::str_sub(Sequence,1,45),
                      Mat_Prb=stringr::str_sub(Sequence,46)
        ) %>%
        dplyr::select(-Sequence) %>%
        dplyr::filter(!is.na(Address)) # Not sure if we need this...
    }
    
    # Add Dat_IDX if provided
    if (!is.null(idx))
      ret_tib <- ret_tib %>%
      dplyr::mutate(Dat_IDX=as.integer(idx))
    
    # Standard Column Data Type Clean-Up::
    ret_tib <- ret_tib %>%
      utils::type.convert() %>%
      dplyr::mutate(across(where(is.factor),  as.character) )
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_tib
}

load_manifestBuildFile = function(file, field="Assay_Design_Id", cols=NULL,
                                  # aqp=NULL, bpn=NULL,
                                  n_max=50, guess=100000,
                                  verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_manifestBuildFile'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    skip_cnt <- get_headerSkipCount(
      file=file, field=field, n_max=n_max,
      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} skip_cnt={skip_cnt}.{RET}"))
    
    file_type <- get_fileSuffix(file)
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} file_type={file_type}.{RET}"))
    
    if (file_type=="csv") {
      ret_tib <- suppressMessages(suppressWarnings(
        readr::read_csv(file, skip=skip_cnt, guess_max=guess, skip_empty_rows=FALSE)))
    } else if (file_type=="tsv" || file_type=="txt" || file_type=="match") {
      ret_tib <- suppressMessages(suppressWarnings(
        readr::read_tsv(file, skip=skip_cnt, guess_max=guess, skip_empty_rows=FALSE)))
    } else {
      stop(glue::glue("{RET}[{funcTag}]:ERROR: Unsupported fileType={file_type}!!!{RET}{RET}"))
      return(NULL)
    }
    ret1_cnt <- print_tib(ret_tib,funcTag, verbose,vt+6,tc, n="ret1")
    
    if (!is.null(cols)) {
      cur_len <- length(names(ret_tib))
      col_len <- length(cols)
      if (cur_len==col_len) {
        ret_tib <- ret_tib %>% purrr::set_names(cols)
      } else {
        stop(glue::glue("{RET}[{funcTag}]:ERROR: Unequal Column Lengths {cur_len} != {col_len} !!!{RET}{RET}"))
        return(NULL)
      }
    }
    
    ret_tib <- ret_tib %>% 
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
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
#                         File Format Predction Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

get_headerSkipCount = function(file, field="Assay_Design_Id", n_max=50,
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'get_headerSkipCount'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting; {RET}",
                   "file='{file}'{RET}",
                   "field='{field}'{RET}") )
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    ret_tib <- 
      readr::read_lines(file=file, n_max=n_max, skip_empty_rows=FALSE) %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(value=stringr::str_remove(value, "^\\s+") )
    
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="raw")
    
    ret_cnt <- ret_tib %>%
      dplyr::mutate(
        Row_Num=dplyr::row_number()) %>% 
      dplyr::filter(stringr::str_starts(value,!!field)) %>% 
      # head(n=1) %>% 
      tail(n=1) %>% 
      dplyr::pull(Row_Num) - 1
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_cnt
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Build Improbe Position Genomic Region:: RDS
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Example Build Commands::
#
# cur_pos_db_grs <- write_impTopGenomeGRS(genBuild="GRCm10", verbose=10, tt=pTracker)
# cur_pos_db_grs <- write_impTopGenomeGRS(genBuild="GRCh37", verbose=10, tt=pTracker)
# cur_pos_db_grs <- write_impTopGenomeGRS(genBuild="GRCh38", verbose=10, tt=pTracker)
write_impTopGenomeGRS = function(genBuild,
                                 impDir="/Users/bretbarnes/Documents/data/improbe/GRCh36-38.mm10.FWD_SEQ.21022021/designInput",
                                 verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'write_impTopGenomeGRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    imp_pos_col <- NULL
    imp_pos_col <- cols(
      imp_chr = col_integer(),
      imp_pos = col_integer(),
      imp_cgn = col_integer(),
      imp_srd = col_integer(),
      imp_top = col_character()
    )
    print(imp_pos_col)
    
    pos_tsv <- ""
    pos_rds <- ""
    
    pos_tsv <- file.path(impDir, paste0(genBuild,'.cgn-pos-fwd.base.tsv.gz') )
    pos_rds <- file.path(impDir, paste0(genBuild,'.cgn-pos-fwd.base.rds') )
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading TSV={pos_tsv}...{RET}"))
    
    pos_tib <- NULL
    pos_tib <- readr::read_tsv(pos_tsv, 
                               col_names=names(imp_pos_col$cols),
                               col_types=imp_pos_col)
    
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr} pos_tib={RET}"))
      print(pos_tib)
    }
    ret_cnt <- pos_tib %>% base::nrow()
    
    if (FALSE) {
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Building GRS..{RET}"))
      
      ret_grs <- NULL
      ret_grs <- 
        GRanges(
          seqnames = Rle(paste0("chr",pos_tib$imp_chr)),
          # strand=Rle(stringr::str_sub( bsp$Bsp_Srd, 1,1 ) ),
          
          imp_cg=pos_tib$imp_cgn,
          imp_fr=pos_tib$imp_fr,
          imp_tb=pos_tib$imp_tb,
          imp_co=pos_tib$imp_co,
          
          IRanges(start=pos_tib$imp_pos,
                  width=2,
                  # end=pos_tib$Coordinate+1,
                  names=paste(pos_tib$imp_cgn,
                              pos_tib$imp_chr,
                              pos_tib$imp_pos,
                              sep="_")
          )
        )
      
      if (verbose>=vt) {
        cat(glue::glue("[{funcTag}]:{tabsStr} ret_grs={RET}"))
        print(ret_grs)
      }
      
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Writing RDS={pos_rds}...{RET}"))
      readr::write_rds(ret_grs, pos_rds, compress="gz")
      
      ret_cnt <- ret_grs %>% length()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  # ret_grs
  pos_tib
}

# Example Build Commands::
#
# cur_pos_db_grs <- write_impGenomeGRS(genBuild="GRCm10", verbose=10, tt=pTracker)
# cur_pos_db_grs <- write_impGenomeGRS(genBuild="GRCh37", verbose=10, tt=pTracker)
# cur_pos_db_grs <- write_impGenomeGRS(genBuild="GRCh38", verbose=10, tt=pTracker)
write_impGenomeGRS = function(genBuild,
                              impDir="/Users/bretbarnes/Documents/data/improbe/designOutput_21092020",
                              verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'write_impGenomeGRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    imp_pos_col <- NULL
    imp_pos_col <- cols(
      imp_cgn = col_character(),
      imp_chr = col_character(),
      imp_pos = col_integer(),
      
      imp_fr  = col_character(),
      imp_tb  = col_character(),
      imp_co  = col_character()
    )
    
    pos_tsv <- ""
    pos_rds <- ""
    
    pos_tsv <- file.path(impDir, paste0(genBuild,'-21092020_improbe-designOutput.cgn-pos-srd.tsv.gz') )
    pos_rds <- file.path(impDir, paste0(genBuild,'-21092020_improbe-designOutput.cgn-pos-srd.rds') )
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Loading TSV={pos_tsv}...{RET}"))
    
    pos_tib <- NULL
    pos_tib <- readr::read_tsv(pos_tsv, 
                               col_names=names(imp_pos_col$cols),
                               col_types=imp_pos_col)
    
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr} pos_tib={RET}"))
      print(pos_tib)
    }
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Building GRS..{RET}"))
    
    
    ret_grs <- NULL
    ret_grs <- 
      GRanges(
        seqnames = Rle(paste0("chr",pos_tib$imp_chr)),
        # strand=Rle(stringr::str_sub( bsp$Bsp_Srd, 1,1 ) ),
        
        imp_cg=pos_tib$imp_cgn,
        imp_fr=pos_tib$imp_fr,
        imp_tb=pos_tib$imp_tb,
        imp_co=pos_tib$imp_co,
        
        IRanges(start=pos_tib$imp_pos,
                width=2,
                # end=pos_tib$Coordinate+1,
                names=paste(pos_tib$imp_cgn,
                            pos_tib$imp_chr,
                            pos_tib$imp_pos,
                            sep="_")
        )
      )
    
    if (verbose>=vt) {
      cat(glue::glue("[{funcTag}]:{tabsStr} ret_grs={RET}"))
      print(ret_grs)
    }
    
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr} Writing RDS={pos_rds}...{RET}"))
    readr::write_rds(ret_grs, pos_rds, compress="gz")
    
    ret_cnt <- ret_grs %>% length()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}",
                   "{tabsStr}# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #{RET}{RET}"))
  
  ret_grs
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Standard UCSC Loading Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

load_ncbi_gene = function(file,grs=FALSE,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_ncbi_gene'
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
                       source="NCBI",
                       source2=NA_character_),
        
        tibble::tibble(chr=dat_tib$chrom,
                       beg=dat_tib$tss_200_beg,
                       end=dat_tib$tss_200_end,
                       srd=dat_tib$strand,
                       name2=dat_tib$name2,
                       name=dat_tib$name,
                       class="TSS200",
                       source="NCBI",
                       source2=NA_character_),
        
        tibble::tibble(chr=dat_tib$chrom,
                       beg=dat_tib$tss_1500_beg,
                       end=dat_tib$tss_1500_end,
                       srd=dat_tib$strand,
                       name2=dat_tib$name2,
                       name=dat_tib$name,
                       class="TSS1500",
                       source="NCBI",
                       source2=NA_character_)
      ) %>% 
      dplyr::group_by(class) %>%
      dplyr::mutate(rank=dplyr::row_number(),
                    unq_key=paste(class,rank, sep="_")) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(chr,beg,end)
    
    if (grs) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
      ret_tib = 
        GRanges(
          seqnames=Rle(ret_tib$chr), 
          strand=Rle(ret_tib$srd),
          
          name=ret_tib$tran,
          name2=ret_tib$gene,
          class=ret_tib$class,
          source=ret_tib$source,
          
          IRanges(start=ret_tib$beg, 
                  end=ret_tib$end, 
                  names=ret_tib$unq_key)
        )
      ret_cnt <- ret_grs %>% names %>% length()
    } else {
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

load_ucsc_gene = function(file,grs=FALSE,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_ucsc_gene'
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
                       source="UCSC",
                       source2=NA_character_),
        
        tibble::tibble(chr=dat_tib$chrom,
                       beg=dat_tib$tss_200_beg,
                       end=dat_tib$tss_200_end,
                       srd=dat_tib$strand,
                       name2=dat_tib$proteinID,
                       name=dat_tib$name,
                       class="TSS200",
                       source="UCSC",
                       source2=NA_character_),
        
        tibble::tibble(chr=dat_tib$chrom,
                       beg=dat_tib$tss_1500_beg,
                       end=dat_tib$tss_1500_end,
                       srd=dat_tib$strand,
                       name2=dat_tib$proteinID,
                       name=dat_tib$name,
                       class="TSS1500",
                       source="UCSC",
                       source2=NA_character_)
      ) %>% 
      dplyr::group_by(class) %>%
      dplyr::mutate(rank=dplyr::row_number(),
                    unq_key=paste(class,rank, sep="_")) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(chr,beg,end)
    
    if (grs) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
      ret_tib = 
        GRanges(
          seqnames=Rle(ret_tib$chr), 
          strand=Rle(ret_tib$srd),
          
          name=ret_tib$name,
          name2=ret_tib$name2,
          Class=ret_tib$class,
          source=ret_tib$source,
          
          IRanges(start=ret_tib$beg, 
                  end=ret_tib$end, 
                  names=ret_tib$unq_key)
        )
      ret_cnt <- ret_grs %>% names %>% length()
    } else {
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    }
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

load_ucsc_cpgs = function(file, grs=FALSE,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'load_ucsc_cpgs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_grs <- NULL
  stime <- system.time({
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading Raw Data={file}...{RET}"))
    dat_tib <- suppressMessages(suppressWarnings( readr::read_tsv(file) )) %>% 
      dplyr::mutate(name=paste0(chrom,'-',chromStart,'-',chromEnd))
    
    ret_tib <- dplyr::bind_rows(
      
      tibble::tibble(chr=dat_tib$chrom,
                     beg=dat_tib$chromStart-4000,
                     end=dat_tib$chromStart-2000,
                     srd="+",
                     name2=NA_character_,
                     name=dat_tib$name,
                     class="NShelf",
                     source="UCSC",
                     source2=NA_character_),
      
      tibble::tibble(chr=dat_tib$chrom,
                     beg=dat_tib$chromStart-2000,
                     end=dat_tib$chromStart,
                     srd="+",
                     name2=NA_character_,
                     name=dat_tib$name,
                     class="NShore",
                     source="UCSC",
                     source2=NA_character_),
      
      tibble::tibble(chr=dat_tib$chrom,
                     beg=dat_tib$chromStart,
                     end=dat_tib$chromEnd,
                     srd="+",
                     name2=NA_character_,
                     name=dat_tib$name,
                     class="Island",
                     source="UCSC",
                     source2=NA_character_),
      
      tibble::tibble(chr=dat_tib$chrom,
                     beg=dat_tib$chromEnd,
                     end=dat_tib$chromEnd+2000,
                     srd="+",
                     name2=NA_character_,
                     name=dat_tib$name,
                     class="SShore",
                     source="UCSC",
                     source2=NA_character_),
      
      tibble::tibble(chr=dat_tib$chrom,
                     beg=dat_tib$chromEnd+2000,
                     end=dat_tib$chromEnd+4000,
                     srd="+",
                     name2=NA_character_,
                     name=dat_tib$name,
                     class="SShelf",
                     source="UCSC",
                     source2=NA_character_)
    ) %>% 
      dplyr::group_by(class) %>%
      dplyr::mutate(rank=dplyr::row_number(),
                    unq_key=paste(class,rank, sep="_")) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(chr,beg,end)
    
    if (grs) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building GRanges...{RET}"))
      ret_tib <-
        GRanges(
          seqnames=Rle(ret_tib$chr), 
          strand=Rle(ret_tib$srd),
          
          name=ret_tib$name,
          class=ret_tib$class,
          source=ret_tib$source,
          
          IRanges(start=ret_tib$beg, 
                  end=ret_tib$end, 
                  names=ret_tib$unq_key)
        )
      ret_cnt <- ret_grs %>% names %>% length()
    } else {
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    }
    
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# End of file
