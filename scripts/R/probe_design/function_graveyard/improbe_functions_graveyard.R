
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   improbe (Infinium Methylation) Methods::
#
#                           [FUNCTION GRAVEYARD]
#
#  The functions have been replaced. This file should be deleted, but is 
#           just a resting place for these functions for now :)
#
#              These functions should probably be deprecated
#
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
BRK <- paste0("# ",
              paste(rep("-----",6),collapse=" "),"|",
              paste(rep("-----",6),collapse=" ")," #")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Standard Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

template_func = function(tib,
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='template_func') {
  
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (verbose>=vt+2) {
    cat(glue::glue("{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr} Function Parameters::{RET}"))
    cat(glue::glue("[{funcTag}]:{tabsStr}   funcTag={funcTag}.{RET}"))
    cat(glue::glue("{RET}"))
  }
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- base::system.time({
    
    # verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    ret_key <- glue::glue("ret-FIN({funcTag})")
    ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n=ret_key)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    String Diff Method from the Internet::
#
#  NOTE:: Maybe useful, so we'll just keep it here for now. 
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

list.string.diff<-function(a="ATTCGA-",b="attTGTT",exclude=c("-","?"),ignore.case=TRUE)
{
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  if(ignore.case==TRUE)
  {
    a<-toupper(a)
    b<-toupper(b)
  }
  seq.a<-unlist(strsplit(a,split=""))
  seq.b<-unlist(strsplit(b,split=""))
  diff.d<-rbind(seq.a,seq.b)
  only.diff<-diff.d[,diff.d[1,]!=diff.d[2,]]
  pos<-which(diff.d[1,]!=diff.d[2,])
  only.diff<-rbind(pos,only.diff)
  for(ex.loop in 1:length(exclude))
  {
    only.diff<-only.diff[,!(only.diff["seq.a",]==exclude[ex.loop]|only.diff["seq.b",]==exclude[ex.loop])]
  }
  return(only.diff)
}

if (FALSE) {
  if (FALSE) {
    cur_tib <- chr_list[[chr_str]]
    chr_idx <- chr_maps[[chr_str]] %>% 
      dplyr::filter(!!chr_sym==chr_str) %>% 
      head(n=1) %>% pull(Idx) %>% as.integer()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 Parse Forward Template Sequence from Genome::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ref_seqs <- parse_genomic_seqs(tib=cur_tib,
                                   seq=as.character(dna_dat$seqs[[chr_idx]]), 
                                   
                                   srd_str=srd_str,
                                   pos_key=pos_key,
                                   chr_str=chr_str,
                                   chr_idx=chr_idx,
                                   
                                   ups_len=ups_len, 
                                   seq_len=seq_len,
                                   verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                Forward Template Sequence Generation:: s-improbe
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    cur_tib <- s_improbe_template(tib=cur_tib,
                                  seqs=ref_seqs,
                                  
                                  chr_str=chr_str, 
                                  chr_idx=chr_idx,
                                  
                                  ext_seq=ext_seq,
                                  iup_seq=iup_seq,
                                  imp_seq=imp_seq,
                                  
                                  iupac=iupac, 
                                  ups_len=ups_len,
                                  seq_len=seq_len,
                                  verbose=verbose, vt=vt+1,tc=tc+1,tt=tt)
  }
  
}
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Modern improbe IO Methods:: 
#                          Not so Modern Any More
#
#              These functions should probably be deprecated
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Common Conversion Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# NOTE:: Pretty sure these function is old and not needed. Its probably not
#   needed anymore and can be moved to the graveyard. Seems useful though...
#
if (FALSE) {
  
  seq_to_prbs = function(tib, seq, ids,
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
    
    stop(glue::glue("{RET}[{funcTag}]:{tabsStr} ERROR: This function isn't ready!!!{RET}{RET}"))
    return(ret_tib)
    
    stime <- base::system.time({
      
      cgn_sym <- rlang::sym(cgn)
      chr_sym <- rlang::sym(chr)
      pos_sym <- rlang::sym(pos)
      ids_sym <- rlang::sym(ids)
      
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
            Name=paste(!!ids_sym,"CG_UP", sep=del),
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
            Name=paste(!!ids_sym,"CG_DN", sep=del),
            Coordinate=as.integer(Coordinate+2)
          )
      }
      
      # Add New Tri-fecta Probes::
      #  - NOTE:: These will have their names changed later after alignment...
      #
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joining and sorting...{RET}"))
      
      cur_tib <- dplyr::bind_rows(cur_tib, ups_tib, dns_tib) %>%
        dplyr::arrange(!!ids_sym)
      
      ret_tib <- dplyr::bind_rows(ret_tib, cur_tib)
      
      if (verbose>=vt+1)
        cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Substring chr_str={chr_str}.{RET}{RET}"))
      
      
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    ret_tib
  }
  
  bed_to_prbs = function(tib, fas,
                         file=NULL,
                         din="Probe_Type",gen="na",
                         cgn="Seq_ID",chr="Chromosome", pos="Coordinate",
                         ext_seq="Fwd_Temp_Seq",
                         iup_seq="Iup_Temp_Seq",
                         imp_seq="Imp_Temp_Seq",
                         ups_len=60, seq_len=122, nrec=0,
                         iupac=NULL, add_flank=FALSE,
                         del="_", 
                         verbose=0,vt=3,tc=1,tt=NULL,
                         funcTag='bed_to_prbs') {
    
    tabsStr <- paste0(rep(TAB, tc), collapse='')
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
    
    ret_cnt <- 0
    ret_tib <- NULL
    stime <- base::system.time({
      
      cgn_sym <- rlang::sym(cgn)
      chr_sym <- rlang::sym(chr)
      pos_sym <- rlang::sym(pos)
      
      # Get List of BSC Genomes::
      # fas_dir <- fas %>% base::dirname()
      fas_pre <- fas %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".fa$")
      if (is.null(gen)) gen <- base::basename(fas_pre)
      if (verbose>=vt) 
        cat(glue::glue("[{funcTag}]:{tabsStr} Prefix({gen})={fas_pre}.{RET}"))
      
      # NOTE:: There isn't an opposite genome, build that from the converted::
      #  This is why we don't need 'srd_cos <- c("C","O")' below::
      #
      srd_cos <- c("C")
      srd_frs <- c("F","R")
      srd_mus <- c("U","M","D")
      
      for (fr in srd_frs) {
        # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} fr={fr}...{RET}"))
        
        for (co in srd_cos) {
          # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} co={co}...{RET}"))
          
          for (mu in srd_mus) {
            if (verbose>=vt) 
              cat(glue::glue("[{funcTag}]:{tabsStr} fr={fr}, co={co}, mu={mu}...{RET}"))
            gen_fas <- paste0(fas_pre,".", fr,co,mu, ".fa.gz")
            
            if (!file.exists(gen_fas)) {
              if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Skipping; gen_fas=({fr},{co},{mu})",
                                              "={gen_fas} Does NOT Exist!{RET}"))
              next
            } else {
              if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading; gen_fas({fr},{co},{mu})",
                                              "={gen_fas}...{RET}"))
            }
            
            cur_tib <- fas_to_seq(tib=tib, fas=gen_fas, file=file,
                                  name=cgn, din=din, gen=gen,
                                  chr1=chr, pos=pos, srd=fr,
                                  ext_seq=ext_seq,iup_seq=iup_seq,imp_seq=imp_seq,
                                  
                                  # See Above:: Will Pass these variables in later...
                                  # ext_seq="Fwd_Temp_Seq",des_seq="Iup_Temp_Seq",
                                  # imp_seq="Imp_Temp_Seq",
                                  
                                  iupac=iupac, mnrec=nrec, add_flank=FALSE,
                                  verbose=verbose,vt=vt+1,tc=tc+1,tt=tt
            ) %>%
              dplyr::mutate(Srd_FR=fr, Srd_CO=co, Prb_Des=mu)
            
            cur_cnt <- print_tib(cur_tib,funcTag, verbose,vt+4,tc, n="cur-pre-join")
            cgn_vec <- cur_tib %>% dplyr::pull(!!cgn_sym)
            
            # TBD:: Add the current FR/CO/TB fields to cur_tib
            # TBD:: Add the opposite designs
            # TBD:: Remove previous information from the input tibble
            
            prb_tib <- NULL
            if (fr=="F") {
              prb_tib <- tibble::tibble(
                !!cgn_sym := cgn_vec,
                
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
                # Seq_ID = cur_tib$Seq_ID,
                # !!cgn_sym := cur_tib[[!!cgn_sym]],
                !!cgn_sym := cgn_vec,
                
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
              stop(glue::glue("{RET}[{funcTag}]: Unsupported fr={fr}...{RET}"))
              return(NULL)
            }
            prb_cnt <- print_tib(prb_tib,funcTag, verbose,vt+4,tc, n="prb")
            
            cur_tib <- cur_tib %>% dplyr::left_join(prb_tib, by=cgn)
            ret_tib <- ret_tib %>% dplyr::bind_rows(cur_tib)
          }
        }
      }
      
      ret_cnt <- print_tib(ret_tib,funcTag, verbose,vt+4,tc, n="ret")
    })
    etime <- stime[3] %>% as.double() %>% round(2)
    if (!is.null(tt)) tt$addTime(stime,funcTag)
    if (verbose>=vt) cat(glue::glue(
      "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
      "{RET}{tabsStr}{BRK}{RET}{RET}"))
    
    ret_tib
  }
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Modern improbe IO Methods:: 
#                          Not so Modern Any More
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

addSeq48U = function(tib, field,
                     sidx=2, plen=50,
                     verbose=0,vt=3,tc=1,tt=NULL,
                     funcTag='addSeq48U') {
  
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
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

isTibIMP = function(tib,
                    verbose=0,vt=3,tc=1,tt=NULL,
                    funcTag='isTibIMP') {
  
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
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
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
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

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
      dplyr::filter(!is.na(fas_line)) %>%
      utils::type.convert() %>% 
      dplyr::mutate(across(where(is.factor),  as.character) )
    
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
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

ordToFas = function(tib, dir, name, verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'ordToFas'
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
  funcTag <- 'loadManifestYH'
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
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
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
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          BSMAP Alignment Functions::
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
  if (verbose>=vt) cat(glue::glue(
    "[{funcTag}]:{tabsStr} Done; Count={ret_cnt}; elapsed={etime}.{RET}",
    "{RET}{tabsStr}{BRK}{RET}{RET}"))
  
  if (run) return(ret_tsv)
  aln_ssh
}

# End of file
