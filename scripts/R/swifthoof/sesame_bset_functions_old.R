

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#            EPIC hg19 Infinium I Inferred SNP Table Construction::
#
#               Should move this scratch to somewhere else::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# 450k Genome Studio::
#
# hm45_man_csv <- '/Users/bretbarnes/Documents/data/manifests/HumanMethylation450_15017482_v.1.2.csv.gz'
# hm45_man_lst <- loadManifestGenomeStudio(file=hm45_man_csv, verbose=opt$verbose+4)
# hm45_man_tib <- hm45_man_lst$man %>% dplyr::select(IlmnID,Genome_Build,CHR,MAPINFO,Strand)
# hm45_ctl_tib <- hm45_man_lst$ctl %>% dplyr::mutate(Address=as.double(Address))

#
# SNP Table
#
epicI_ses_grs <- sesameDataPullVariantAnno_InfiniumI(platform = "EPIC")
epicR_ses_grs <- sesameDataPullVariantAnno_SNP()

epicI_ses_tib <- epicI_ses_grs %>% 
  as.data.frame() %>% tibble::rownames_to_column(var="IlmnID") %>% tibble::as_tibble()

# Five Examples by hand::
#  gs_ses_csv <- '/Users/bretbarnes/Documents/tools/notes/inferred-snps.GS-vs-Ses.csv'
#  gs_ses_tib <- readr::read_csv(gs_ses_csv)
#
#

# epic_man_csv <- file.path(par$datDir, 'manifest/base/EPIC-B4.manifest.sesame-base.cpg-sorted.csv.gz')
epic_man_csv <- '/Users/bretbarnes/Documents/data/manifests/MethylationEPIC_v-1-0_B4.csv.gz'
epic_man_lst <- loadManifestGenomeStudio(file=epic_man_csv, verbose=opt$verbose+4)
epic_ctl_tib <- epic_man_lst$ctl %>% dplyr::mutate(Address=as.double(Address))
epic_man_tib <- epic_man_lst$man %>% 
  dplyr::select(IlmnID,Genome_Build,CHR,MAPINFO,Strand,Next_Base,Infinium_Design_Type) %>%
  dplyr::filter(Infinium_Design_Type=='I') %>% 
  dplyr::mutate(CHR=paste0('chr',CHR)) %>% dplyr::arrange(CHR,MAPINFO)

#
# Validation Shown Below::
#
join_man_tib <- epic_man_tib %>% 
  dplyr::inner_join(epicI_ses_tib, by="IlmnID", suffix=c("_Man", "_Ses")) %>%
  dplyr::mutate(REF_CON=dplyr::case_when(
    strand=='+' & REF=='G' ~ 'A',
    strand=='-' & REF=='C' ~ 'T',
    TRUE ~ REF)
  )
join_sum_tib <- join_man_tib %>% dplyr::filter(CHR==seqnames) %>% 
  dplyr::mutate(Pos_Dif=MAPINFO-start) %>% 
  dplyr::group_by(Pos_Dif,Strand,strand,Next_Base,REF_CON) %>% 
  dplyr::summarise(Count=n(), .groups='drop')

#
# Better, but NOT working attempt below::
#
if (FALSE) {
  join_man_tib3 <- epic_man_tib %>% 
    dplyr::inner_join(epicI_ses_tib, by="IlmnID", suffix=c("_Man", "_Ses")) %>%
    dplyr::mutate(
      REF_RC=dplyr::case_when(
        strand=='-' ~ revCmp(REF),
        TRUE ~ REF),
      REF_CO=dplyr::case_when(
        strand=='+' & REF_RC=='G' ~ 'A',
        strand=='-' & REF_RC=='C' ~ 'T',
        TRUE ~ REF)
    )
  join_sum_tib3 <- join_man_tib3 %>% 
    dplyr::filter(CHR==seqnames) %>% 
    dplyr::mutate(Pos_Dif=MAPINFO-start) %>% 
    dplyr::group_by(Pos_Dif,Strand,strand,Next_Base,REF_CO) %>% 
    dplyr::summarise(Count=n(), .groups='drop')
}

#
# hub19 doesn't have the most up to date dbSNP
#
hub_hg19 <- subset(hub, (hub$species == "Homo sapiens") & (hub$genome == "hg19"))
length(hub_hg19)
hub_hg19$title[grep("SNP", hub_hg19$title)]

hub_snp137 <- subset(hub_hg19, title=='Common SNPs(137)')

hub_snp137_grs <- hub_snp137[['AH5105']]
hub_snp137_tib <- hub_snp137_grs %>% as.data.frame() %>% tibble::as_tibble()

#
# Check overlap of Sesame and dbSNP(137) list
#
epicI_ses_tib %>% dplyr::inner_join(hub_snp137_tib, by=c("rs"="name"), suffix=c("_sesI", "_db137"))



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#             Every thing below is not used and can be removed::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (FALSE) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Add Swapped Summary Percentages::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (!is.null(opt$skipSwap) && !opt$skipSwap) {
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Building raw_sset_tib...{RET}") )
    
    cur_sset_tib <- NULL
    raw_sset_tib <- NULL
    opt$writeSsetRaw <- FALSE
    raw_sset_tib <- sset2tib(sset=raw_sset, by="Probe_ID", des="Probe_Design",  
                             percision=opt$percisionSigs, sort=FALSE, 
                             save=opt$writeSsetRaw, csv=raw_sset_csv, 
                             verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    if (retData) ret$raw_sset_tib <- raw_sset_tib
  }
  
  if (!is.null(raw_sset_tib) && !is.null(cur_sset_tib)) {
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Adding Inferred Sample Swapped Stats.{RET}"))
    
    swap_sum_tib <- NULL
    swap_sum_tib <- dplyr::inner_join(joinSsetTibInfI(tib=raw_sset_tib), 
                                      joinSsetTibInfI(tib=cur_sset_tib), 
                                      by="Probe_ID", suffix=c("_Raw", "_Cur") ) %>%
      dplyr::mutate(isSwap=dplyr::case_when(
        Probe_Design_Inb_Raw==Probe_Design_Inb_Cur & Probe_Design_Oob_Raw==Probe_Design_Oob_Cur ~ 'Reference',
        Probe_Design_Inb_Raw!=Probe_Design_Inb_Cur & Probe_Design_Oob_Raw!=Probe_Design_Oob_Cur ~ 'Alternate',
        TRUE ~ NA_character_
      )) %>% dplyr::group_by(Probe_Design_Inb_Raw,isSwap) %>% 
      dplyr::summarise(Count=n(), .groups='drop') %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(Total=sum(Count)) %>% 
      tidyr::unite(Type, Probe_Design_Inb_Raw, isSwap, sep='_') %>% 
      dplyr::mutate(Perc=round(100*Count/Total, 3)) %>% 
      dplyr::select(Type, Perc) %>% tidyr::spread(Type, Perc) %>%
      purrr::set_names(paste(names(.),'Perc', sep='_') )
    
    ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(swap_sum_tib)
    
    ssheet_ncols <- ssheet_tib %>% base::ncol()
    if (verbose>=vt) 
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Done. Binding Sample Sheet (with channel-swap-stats) ssheet_ncols={ssheet_ncols}.{RET}{RET}"))
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   OLD Sesame SSET Manipulation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sesameStepAbbreviation = function(x, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sesameStepAbbreviation'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (x=='raw') return('R')
  if (x=='dyeBiasCorrTypeINorm') return('D')
  if (x=='detectionPnegEcdf') return('N')
  if (x=='pOOBAH') return('P')
  if (x=='noob') return('N')
  if (x=='noobsb') return('S')
  if (x=='inferTypeIChannel') return('I')
  stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported Sesame Abbreviation={x}!!!{RET}{RET}") )
  
  return('U')
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     SSET to BeadSET Conversion Methods::
#
#                                 NOT USED!
#                                 NOT USED!
#                                 NOT USED!
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sset2bset = function(sset, addSigs=TRUE, addNegs=TRUE, addPoob=TRUE, addBeta=TRUE,
                     
                     # getBeta Parameters::
                     quality.mask=FALSE, nondetection.mask=FALSE, 
                     mask.use.tcga=FALSE, pval.threshold=1, sum.TypeI=FALSE,
                     as.enframe=TRUE,
                     round_dat=TRUE, round_pval=6, round_beta=4,
                     verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sset2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  bset <- NULL
  stime <- system.time({
    sigs <- dplyr::bind_rows(
      sset@IG   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='IG'),
      sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='OG'),
      
      sset@IR   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='IR'),
      sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='OR'),
      
      sset@II   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col=NA)
    )
    if (round_dat) sigs <- sigs %>% dplyr::mutate_if(is.numeric, list(as.integer))
    
    if (addSigs) { bset <- sigs
    } else { bset <- dplyr::select(sigs, 'Probe_ID', 'col') }
    
    if (addNegs) {
      negs <- sset %>% sesame::detectionPnegEcdf() %>% sesame::pval() %>% tibble::enframe(name='Probe_ID', value='negs')
      if (round_dat && round_pval>0) negs <- negs %>% dplyr::mutate_if(is.numeric, list(round), round_pval)
      bset <- bset %>% dplyr::left_join(negs, by='Probe_ID')
    }
    
    if (addPoob) {
      poob <- sset %>% sesame::pOOBAH() %>% sesame::pval() %>% tibble::enframe(name='Probe_ID', value='poob')
      if (round_dat && round_pval>0) poob <- poob %>% dplyr::mutate_if(is.numeric, list(round), round_pval)
      bset <- bset %>% dplyr::left_join(poob, by='Probe_ID')
    }
    
    if (addBeta) {
      beta <- sesame::getBetas(sset=sset, quality.mask=quality.mask, 
                               nondetection.mask=nondetection.mask, 
                               mask.use.tcga=mask.use.tcga, 
                               pval.threshold=pval.threshold, 
                               sum.TypeI=sum.TypeI)
      if (as.enframe) beta <- beta %>% tibble::enframe(name='Probe_ID', value='beta')
      if (round_dat && round_beta>0) beta <- beta %>% dplyr::mutate_if(is.numeric, list(round), round_beta)
      bset <- bset %>% dplyr::left_join(beta, by='Probe_ID')
    }
    
    # bset <- bset %>% 
    #   dplyr::left_join(negs, by='Probe_ID') %>%
    #   dplyr::left_join(poob, by='Probe_ID') %>%
    #   dplyr::left_join(beta, by='Probe_ID')
  })
  bset_nrows <- bset %>% base::nrow()
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done bset_nrows={bset_nrows}, elapsed={etime}.{RET}{RET}"))
  
  bset
}

ssetBeta2bset = function(sset, bset, nkey, del='_', 
                         quality.mask=FALSE, nondetection.mask=FALSE,
                         mask.use.tcga=FALSE, pval.threshold=1, sum.TypeI=FALSE,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetBeta2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  datTag  <- 'beta'
  betaTag <- paste(nkey,datTag, sep=del) # %>% rlang::sym()
  
  stime <- system.time({
    dat <- NULL
    
    betas <- sesame::getBetas(sset=sset,
                              quality.mask=quality.mask, 
                              nondetection.mask=nondetection.mask, 
                              mask.use.tcga=mask.use.tcga, 
                              pval.threshold=pval.threshold, 
                              sum.TypeI=sum.TypeI) %>%
      tibble::enframe(name='Probe_ID', value=betaTag)
    
    dat <- add2bset(bset=bset, inf1=betas, keyA='Probe_ID', verbose=verbose,vt=1,tc=0,tt=tt)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

ssetPval2bset = function(sset, bset, nkey, pkey, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetPval2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  datTag  <- 'pval'
  pvalTag <- paste(nkey,pkey,datTag, sep=del) # %>% rlang::sym()
  
  stime <- system.time({
    dat <- NULL
    
    pvals <- sset@pval %>% tibble::enframe(name='Probe_ID', value=pvalTag)
    dat <- add2bset(bset=bset, inf1=pvals, keyA='Probe_ID', verbose=verbose,vt=1,tc=0,tt=tt)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

ssetSigs2bset = function(sset, bset, nkey, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetSig2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    dat <- NULL
    
    datTag <- 'sig'
    grnTag <- paste(nkey,'Grn',datTag, sep=del) %>% rlang::sym()
    redTag <- paste(nkey,'Red',datTag, sep=del) %>% rlang::sym()
    
    IG_tib <- sset@IG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    IR_tib <- sset@IR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    OG_tib <- sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    OR_tib <- sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    I1 <- dplyr::bind_rows(dplyr::inner_join(IG_tib,OR_tib, by=c('Probe_ID', 'Design_Type')),
                           dplyr::inner_join(OG_tib,IR_tib, by=c('Probe_ID', 'Design_Type')) )
    
    I2 <- sset@II %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      dplyr::rename(!!redTag :=U, !!grnTag :=M)
    
    dat <- add2bset(bset=bset, inf1=I1, inf2=I2, keyA='Probe_ID',keyB='Design_Type', 
                    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

add2bset = function(bset, inf1, inf2=NULL, keyA=NULL,keyB=NULL,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'add2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    dat <- NULL
    if (is.null(inf2)) inf2 <- inf1
    if (is.null(keyB)) {
      dat[['I1']] <- bset[[1]] %>% dplyr::inner_join(inf1, by=c(keyA))
    } else {
      dat[['I1']] <- bset[[1]] %>% dplyr::inner_join(inf1, by=c(keyA,keyB))
    }
    dat[['I2']] <- bset[[2]] %>% dplyr::inner_join(inf2, by=c(keyA))
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

ssetToBetaTib_old = function(sset, name, 
                             quality.mask = FALSE, nondetection.mask = FALSE, 
                             correct.switch = TRUE, mask.use.tcga = FALSE, 
                             pval.threshold = 1, 
                             pval.method = "pOOBAH", sum.TypeI = FALSE,
                             as.enframe=FALSE, percision=0,
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToBetaTib_old'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}.{RET}"))
  
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}      quality.mask={quality.mask}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} nondetection.mask={nondetection.mask}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}    correct.switch={correct.switch}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}     mask.use.tcga={mask.use.tcga}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}    pval.threshold={pval.threshold}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}         sum.TypeI={sum.TypeI}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}       pval.method={pval.method}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}        as.enframe={as.enframe}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}         percision={percision}.{RET}"))
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}              sset={RET}"))
  if (verbose>=vt+4) print(sset)
  if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{RET}{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    ret_tib <- sesame::getBetas(sset=sset,
                                quality.mask = quality.mask, 
                                nondetection.mask = nondetection.mask, 
                                correct.switch = correct.switch, 
                                mask.use.tcga = mask.use.tcga, 
                                pval.threshold = pval.threshold, 
                                pval.method = pval.method, 
                                sum.TypeI = sum.TypeI)
    
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} ret_tib={RET}"))
    if (verbose>=vt+4) ret_tib %>% head() %>% print()
    if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr}{RET}{RET}"))
    
    if (percision!=0) ret_tib <- round(ret_tib, percision)
    if (as.enframe) ret_tib <- ret_tib %>% 
        tibble::enframe(name='Probe_ID', value=name)
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}


#
# TBD:: ssetToPrbTib (Probe_ID,beta,pvals...)
#
# ssetToPrbTib(sset=rdat$raw_sset, verbose=10)
ssetToPrbTib_old = function(sset,
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToPrbTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # Validate p-value slots::
    #
    slot_names <- slotNames(sset)
    if (grep("pvals",slot_names) %>% length() == 0) {
      stop(glue::glue("{RET}[{funcTag}]: ERROR; Failed to find pval in sset!!!{RET}{RET}"))
      return(ret_tib)
    }
    
    # Process each p-value::
    #
    pval_names <- sset@pval %>% names
    for (pval_name in pval_names) {
      # list = ret_tib <- sset@pval[pval_name]
      # dobl = ret_tib <- sset@pval[[pval_name]]
      
      ret_tib <- sset@pval[[pval_name]]
    }
    
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; Return Count={ret_cnt}; elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

sset2calls_old = function(sset, workflow, 
                          quality.mask = FALSE, nondetection.mask = FALSE, 
                          correct.switch = TRUE, mask.use.tcga = FALSE, 
                          pval.threshold = 1, 
                          pval.method = "pOOBAH", sum.TypeI = FALSE,
                          
                          as.enframe=FALSE, percisionBeta=0, percisionPval=0,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sset2calls_old'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}.{RET}"))
  if (verbose>=vt+4) print(sset)
  
  ret_cnt <- 0
  ret_tib <- NULL
  stime <- system.time({
    
    # noob:: beta
    #
    name <- paste(workflow,'beta', sep='_')
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Mutating/Settting name={name}...{RET}"))
    
    #
    # ssetToBeta (provide return type=tib/dat)
    #
    beta <- ssetToBetaTib(sset=sset, name=name, 
                          as.enframe=as.enframe,
                          percision=percisionBeta, 
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    
    
    ret_tib <- tibble::enframe(beta, name='Probe_ID', value=name)
    if (verbose>=vt+4) head(ret_tib) %>% print()
    
    # PnegEcdf:: negs
    #
    name <- paste(workflow,'negs', sep='_')
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Mutating/Setting name={name}...{RET}"))
    
    ssetA <- mutateSesame(sset=sset, method='detectionPnegEcdf', 
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (verbose>=vt+4) print(ssetA)
    
    pvalA <- ssetToPvalTib(sset=ssetA, method='PnegEcdf', name=name, 
                           percision=percisionPval, 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    # if (verbose>=vt+4) head(pval) %>% print()
    ret_tib <- ret_tib %>% dplyr::left_join(pvalA, by="Probe_ID")
    
    # pOOBAH:: poob
    #
    name <- paste(workflow,'poob', sep='_')
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Mutating/Setting name={name}...{RET}"))
    
    ssetB <- mutateSesame(sset=sset, method='pOOBAH', 
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    if (verbose>=vt+4) print(ssetB)
    
    pvalB <- ssetToPvalTib(sset=sset, method='pOOBAH', name=name, 
                           percision=percisionPval, 
                           verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    # if (verbose>=vt+4) head(pval) %>% print()
    ret_tib <- ret_tib %>% dplyr::left_join(pvalB, by="Probe_ID")
    ret_cnt <- ret_tib %>% base::nrow()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done ret_cnt={ret_cnt}, elapsed={etime}.{RET}{RET}"))
  
  ret_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Extracted Sesame Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

inferSexKaryotypes_Copy = function (sset) 
{
  stopifnot(is(sset, "SigSet"))
  sex.info <- getSexInfo_Copy(sset)
  auto.median <- median(sex.info[paste0("chr", seq_len(22))], 
                        na.rm = TRUE)
  XdivAuto <- sex.info["medianX"]/auto.median
  YdivAuto <- sex.info["medianY"]/auto.median
  if (XdivAuto > 1.2) {
    if (sex.info["fracXlinked"] >= 0.5) 
      sexX <- "XaXi"
    else if (sex.info["fracXmeth"] > sex.info["fracXunmeth"]) 
      sexX <- "XiXi"
    else sexX <- "XaXa"
  }
  else {
    if (sex.info["fracXmeth"] > sex.info["fracXunmeth"]) 
      sexX <- "Xi"
    else sexX <- "Xa"
  }
  if ((sexX == "Xi" || sexX == "Xa") && XdivAuto >= 1 && sex.info["fracXlinked"] >= 
      0.5) 
    sexX <- "XaXi"
  if (YdivAuto > 0.3 || sex.info["medianY"] > 2000) 
    sexY <- "Y"
  else sexY <- ""
  karyotype <- paste0(sexX, sexY)
  karyotype
}

inferSex_Copy = function (sset) 
{
  stopifnot(is(sset, "SigSet"))
  sex.info <- getSexInfo_Copy(sset)[seq_len(3)]
  as.character(predict(sesameDataGet("sex.inference"), sex.info))
}

# 
# sesameData::sesameDataGet(paste0(rdat$raw_sset@platform, '.probeInfo'))$chrY.clean
# sex_info <- getSexInfo_Copy(sset=rdat$raw_sset)
# kar_info <- sesame::inferSexKaryotypes(sset=rdat$raw_sset)
getSexInfo_Copy = function (sset) 
{
  if (is(sset, "SigSetList")) 
    return(do.call(cbind, lapply(sset, getSexInfo)))
  stopifnot(is(sset, "SigSet"))
  cleanY <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$chrY.clean
  # cleanY %>% length() %>% print()
  
  xLinked <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$chrX.xlinked
  # xLinked %>% length() %>% print()
  
  probe2chr <- sesameDataGet(paste0(sset@platform, ".probeInfo"))$probe2chr.hg19
  # print(probe2chr)
  
  xLinkedBeta <- sesame::getBetas(sset=sesame::subsetSignal(sset, xLinked), 
                                  quality.mask = FALSE)
  intens <- sesame::totalIntensities(sset)
  probes <- intersect(names(intens), names(probe2chr))
  intens <- intens[probes]
  probe2chr <- probe2chr[probes]
  # print(probe2chr)
  
  # return( sesame::subsetSignal(sset, cleanY) )
  # return( median(sesame::totalIntensities(sesame::subsetSignal(sset, cleanY))) )
  
  c(medianY = median(sesame::totalIntensities(sesame::subsetSignal(sset, cleanY)), na.rm=TRUE), 
    medianX = median(sesame::totalIntensities(sesame::subsetSignal(sset, xLinked)), na.rm=TRUE), fracXlinked = 
      (sum(xLinkedBeta > 0.3 & xLinkedBeta < 0.7, na.rm = TRUE)/sum(!(is.na(xLinkedBeta))) ), 
    fracXmeth = (sum(xLinkedBeta > 0.7, na.rm = TRUE)/sum(!(is.na(xLinkedBeta)))), 
    fracXunmeth = (sum(xLinkedBeta < 0.3, na.rm = TRUE)/sum(!(is.na(xLinkedBeta)))), 
    tapply(intens, probe2chr, median, na.rm=TRUE))
}

# End of file
