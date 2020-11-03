
tru_epic_tsv <- '/Users/bretbarnes/Documents/data/CustomContent/TruDx/TruEPIC-Mikey.CGN-tangos.903487.tsv.gz'
tru_epic_tib <- suppressMessages(suppressWarnings( readr::read_tsv(tru_epic_tsv) )) %>% 
  dplyr::rename(IlmnID=TargetID, AddressA_ID=ProbeID_A, AddressB_ID=ProbeID_B)

b4_epic_csv <- '/Users/bretbarnes/Documents/data/manifests/MethylationEPIC_v-1-0_B2.csv.gz'
b2_450k_csv <- '/Users/bretbarnes/Documents/data/manifests/HumanMethylation450_15017482_v.1.2.csv.gz'

b4_epic_lst <- loadManifestGenomeStudio(b4_epic_csv, addSource=TRUE, normalize=TRUE, verbose=10)
b2_450k_lst <- loadManifestGenomeStudio(b2_450k_csv, addSource=TRUE, normalize=TRUE, verbose=10)

b4_epic_lst$man %>% dplyr::anti_join(tru_epic_tib, by="IlmnID") %>% base::nrow()
tru_epic_tib %>% dplyr::anti_join(b4_epic_lst$man, by="IlmnID") %>% base::nrow()

b2_450k_lst$man %>% dplyr::anti_join(tru_epic_tib, by="IlmnID") %>% base::nrow()
tru_epic_tib %>% dplyr::anti_join(b2_450k_lst$man, by="IlmnID") %>% base::nrow()

ses_cols <- c('Probe_ID','M','U','DESIGN','COLOR_CHANNEL','col','Probe_Type','Probe_Source','Next_Base')

b1_EPIC_tib <- dplyr::bind_rows(
  b2_450k_lst$man %>% dplyr::inner_join(tru_epic_tib, by="IlmnID") %>% dplyr::select(IlmnID,Next_Base,Color_Channel) %>% dplyr::arrange(IlmnID),
  b4_epic_lst$man %>% dplyr::inner_join(tru_epic_tib, by="IlmnID") %>% dplyr::select(IlmnID,Next_Base,Color_Channel) %>% dplyr::arrange(IlmnID) ) %>%
  dplyr::distinct(IlmnID, .keep_all=TRUE) %>%
  dplyr::rename(COLOR_CHANNEL=Color_Channel) %>%
  dplyr::mutate(col=stringr::str_sub(COLOR_CHANNEL,1,1) )

ctl_only_csv <- '/Users/bretbarnes/Documents/data/CustomContent/Genknowme/manifest/GKME_Sesame_EPIC-Controls-Only.csv'
ctl_only_tib <- readr::read_csv(ctl_only_csv, col_names=c('Probe_ID','M','U','DESIGN')) %>%
  dplyr::mutate(Probe_Class='ctl',Probe_Type=Probe_Class,IlmnID=Probe_ID,
                COLOR_CHANNEL='Both',col=NA_character_,Next_Base=NA_character_,
                Infinium_Design=2,Rep_Num=1)

fin_epic_tib <- 
  tru_epic_tib %>% dplyr::mutate(
    AddressB_ID=dplyr::case_when(AddressA_ID==AddressB_ID ~ NA_real_, TRUE ~ AddressB_ID), 
    DESIGN=dplyr::case_when(is.na(AddressB_ID) ~ 'II', TRUE ~ 'I')) %>%
  dplyr::rename(M=AddressB_ID,U=AddressA_ID) %>%
  dplyr::left_join(b1_EPIC_tib, by="IlmnID") %>%
  dplyr::rename(Probe_ID=IlmnID) %>%
  dplyr::bind_rows(ctl_only_tib)

fin_epic_csv <- file.path(par$datDir, paste('EPIC-B1.manifest.sesame-base.cpg-sorted.csv.gz') )
readr::write_csv(fin_epic_tib,fin_epic_csv)

#
# Scratch for COVID Direct
#

pos1_csv <- '/Users/bretbarnes/Documents/scratch/COVID_All/swifthoof_main/COVID-Direct-Set1/204756130015/24x1/UNK/204756130015_R01C01_1_COVID_C1_d-sset.csv.gz'
pos2_csv <- '/Users/bretbarnes/Documents/scratch/COVID_All/swifthoof_main/COVID-Direct-Set1/204756130015/24x1/UNK/204756130015_R03C01_1_COVID_C1_d-sset.csv.gz'
pos3_csv <- '/Users/bretbarnes/Documents/scratch/COVID_All/swifthoof_main/COVID-Direct-Set1/204756130015/24x1/UNK/204756130015_R05C01_1_COVID_C1_d-sset.csv.gz'

pos1_tib <- readr::read_csv(pos1_csv)
pos2_tib <- readr::read_csv(pos2_csv)
pos3_tib <- readr::read_csv(pos3_csv)


neg1_csv <- '/Users/bretbarnes/Documents/scratch/COVID_All/swifthoof_main/COVID-Direct-Set1/204756130015/24x1/UNK/204756130015_R08C02_1_COVID_C1_d-sset.csv.gz'
neg2_csv <- '/Users/bretbarnes/Documents/scratch/COVID_All/swifthoof_main/COVID-Direct-Set1/204756130015/24x1/UNK/204756130015_R10C02_1_COVID_C1_d-sset.csv.gz'
neg3_csv <- '/Users/bretbarnes/Documents/scratch/COVID_All/swifthoof_main/COVID-Direct-Set1/204756130015/24x1/UNK/204756130015_R12C02_1_COVID_C1_d-sset.csv.gz'

neg1_tib <- readr::read_csv(neg1_csv)
neg2_tib <- readr::read_csv(neg2_csv)
neg3_tib <- readr::read_csv(neg3_csv)


#
# Infinium II::
#
pos1_tib %>% dplyr::filter(Probe_Design=='II') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
pos2_tib %>% dplyr::filter(Probe_Design=='II') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
pos3_tib %>% dplyr::filter(Probe_Design=='II') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))

neg1_tib %>% dplyr::filter(Probe_Design=='II') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
neg2_tib %>% dplyr::filter(Probe_Design=='II') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
neg3_tib %>% dplyr::filter(Probe_Design=='II') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))


#
# Infinium I::
#

pos1_tib %>% dplyr::filter(Probe_Design=='IG') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
pos2_tib %>% dplyr::filter(Probe_Design=='IG') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
pos3_tib %>% dplyr::filter(Probe_Design=='IG') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))

neg1_tib %>% dplyr::filter(Probe_Design=='IG') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
neg2_tib %>% dplyr::filter(Probe_Design=='IG') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
neg3_tib %>% dplyr::filter(Probe_Design=='IG') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))

pos1_tib %>% dplyr::filter(Probe_Design=='IR') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
pos2_tib %>% dplyr::filter(Probe_Design=='IR') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
pos3_tib %>% dplyr::filter(Probe_Design=='IR') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))

neg1_tib %>% dplyr::filter(Probe_Design=='IR') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
neg2_tib %>% dplyr::filter(Probe_Design=='IR') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))
neg3_tib %>% dplyr::filter(Probe_Design=='IR') %>% dplyr::summarise(M_sum=sum(M,na.rm=TRUE), U_sum=sum(U,na.rm=TRUE))




# End of file
