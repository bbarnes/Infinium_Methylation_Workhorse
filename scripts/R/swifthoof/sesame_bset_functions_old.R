

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

# End of file
