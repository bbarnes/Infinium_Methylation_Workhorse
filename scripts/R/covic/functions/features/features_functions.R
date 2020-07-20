
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Machine Learning Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )
suppressWarnings(suppressPackageStartupMessages(require("grid")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Load Training Features IO::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadCgnFeatures = function(csv, id="Probe_ID",
                           verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'loadCgnFeatures'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  stopifnot(file.exists(csv))
  
  id <- rlang::sym(id)
  tib <- suppressMessages(suppressWarnings( readr::read_csv(csv) )) %>%
    dplyr::select(!!id) %>% dplyr::arrange(!!id)
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     DML Script Launching Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# launchDML(dir=opt$outDir, exe=par$dmlExe, rscript=opt$Rscript, verbose={opt$verbosity})
launchDML = function(dir, exe, rscript, 
                     verbose=0,vt=2,tc=1,tt=NULL) {
  funcTag <- 'launchDML'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  cmd_bool <- ''
  cmd_strs <- glue::glue(" --outDir={dir} --buildDir={dir}")
  cmd_full <- glue::glue("{rscript} {exe} {cmd_bool} {cmd_strs} --verbosity={verbose}{RET}")
  
  run_name <- paste(stringr::str_remove(base::basename(exe),'.R$'),'full-data', sep='-')
  runShell <- file.path(par$shellDir, paste0('run.',run_name,'.sh'))
  lanShell <- file.path(par$shellDir, paste0('lan.',run_name,'.sh'))
  
  # Change this to exePath -> exeLocPrepPath
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: Writing runShell={runShell}{RET}!"))
  readr::write_file(cmd_full, path=runShell)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: chmod 777 runShell={runShell}{RET}!"))
  Sys.chmod(runShell, mode="0777")
  
  # Add cluster execute if avialbel (i.e. linux)
  cmd_lan <- glue::glue("{runShell}{RET}")
  if (par$isLinux) cmd_lan <- glue::glue("{par$lan_exe} dmls-full {runShell}{RET}")
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: Writing lanShell={lanShell}{RET}!"))
  readr::write_file(cmd_lan, path=lanShell)
  Sys.chmod(lanShell, mode="0777")
  
  # Launch Script
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: Launching lanShell={lanShell}{RET}!"))
  if (opt$execute) system(lanShell)
  
  TRUE
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         dBL (deltaBeta Loci) Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

crossSampleLociDeltaBeta = function(tib, id='Probe_ID', max=NULL, retData=FALSE,
                                    cpp.verbose=0,verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'crossSampleLociDeltaBeta'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (!is.null(max)) tib <- tib %>% head(n=max)
  
  stime <- system.time({
    # Builds Samle Counts and Sample Names vectors::
    # sam_tib <- tib %>% dplyr::select(-id) %>% names() %>% stringr::str_remove('_.*$') %>% 
    sam_tib <- tib %>% dplyr::select(-all_of(id)) %>% names() %>% stringr::str_remove('_.*$') %>% 
      tibble::enframe() %>% dplyr::group_by(value) %>% summarise(Count=n())
    sam_vec <- sam_tib %>% dplyr::pull(1) %>% as.vector()
    cnt_vec <- sam_tib %>% dplyr::pull(2) %>% as.vector()
    
    # Builds Summary Matrix (c++) and converts to tibble::
    sam_mat <- tib %>% column_to_rownames(id) %>% as.matrix()
    sum_mat <- C_crossSampleLociRSquared(sam_mat, cnt_vec, sam_vec, verbose=cpp.verbose)
    
    # TBD:: Need to fix this::
    sum_tib <- sum_mat %>% tibble::as_tibble(rownames=id)
    # sum_tib <- sum_mat
    
    # cpg_vec <- tib %>% dplyr::pull(1) %>% as.vector()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  sum_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              DML Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

dmlsToTib = function(dml, verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'dmlsToTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  dml_tibs <- NULL
  stime <- system.time({
    
    for (type_key in names(dml)) {
      # dml_tib <- sesame::topLoci(dml[[type_key]]) %>% as.data.frame() %>% rownames_to_column('Probe_ID') %>%
      dml_tib <- topLoci_Local(dml[[type_key]]) %>% as.data.frame() %>% rownames_to_column('Probe_ID') %>%
        tibble::as_tibble() %>% purrr::set_names(stringr::str_replace_all(names(.), '\\.', '_')) %>%
        purrr::set_names(names(.) %>% stringr::str_replace_all("[^[A-Za-z,]]+", "_") %>% stringr::str_remove('_$') ) %>%
        dplyr::mutate_if(is.double, list(round), 6) %>%
        dplyr::mutate(Rank=row_number()) %>% dplyr::mutate(Type=type_key)
      
      dml_tibs <- dml_tibs %>% dplyr::bind_rows(dml_tib)
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  dml_tibs
}

betasToDMLs = function(betas, verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'betasToDMLs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  dml_ret <- NULL
  stime <- system.time({
    
    dml_dat <- NULL
    dml_dat$betas <- tibble::column_to_rownames(betas, "Probe_ID") %>% as.matrix()
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} dml_dat-betas::{RET}"))
      dml_dat$betas %>% head(n=5) %>% print()
    }
    dml_dat$sampleInfo <- dml_dat$betas %>% colnames() %>% tibble::enframe() %>% 
      dplyr::mutate(patient=value) %>%
      tidyr::separate(patient, into=c('type', 'rep'), sep='_', remove=FALSE) %>% 
      dplyr::select(-name) %>% tibble::column_to_rownames('value')
    
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} DML_PRE::betas{RET}"))
      dml_dat$betas %>% head(n=5) %>% print()
      cat(glue::glue("[{funcTag}]:{tabsStr} DML_PRE::sampleInfo{RET}"))
      dml_dat$sampleInfo %>% head(n=5) %>% print()
    }
    
    # DML::
    # dml_ret <- sesame::DML(dml_dat$betas, dml_dat$sampleInfo, ~type)
    dml_ret <- DML_Local(dml_dat$betas, dml_dat$sampleInfo, ~type)
    if (verbose>=vt+4) {
      cat(glue::glue("[{funcTag}]:{tabsStr} DML_RET::{RET}"))
      dml_ret %>% head(n=5) %>% print()
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  dml_ret
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Rebranded Copied Sesame Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# These methods address differenet versions of anaconda that may not have
#  support the most recent version of Sesame...

topLoci_Local <- function(cf1) {
  cf1[order(cf1[,'Pr(>|t|)']),]
}

DML_Local = function(betas, sample.data, formula, se.lb=0.06, balanced=FALSE, cf.test=NULL) {
  
  design <- model.matrix(formula, sample.data)
  ## convert to factor for faster processing
  design.fac <- data.frame(lapply(as.data.frame(design), as.factor))
  rdf0 <- unlist(
    lapply(design.fac, function(x) max(min(tabulate(x)),1)),
    recursive = FALSE, use.names=FALSE)
  
  n.cpg <- dim(betas)[1]
  n.cf <- dim(design)[2]
  
  ## cf.test specify the factors to be reported
  if (is.null(cf.test)) {
    cf.test <- colnames(design)
    cf.test <- cf.test[cf.test != '(Intercept)']
  } else if (cf.test == 'all') {
    cf.test <- colnames(design)
  }
  
  ## preprare output
  cf <- lapply(cf.test, function(cfi) matrix(
    data=NA, nrow=n.cpg, ncol=5,
    dimnames=list(rownames(betas),
                  c('Estimate', 'Std. Error', 't-stat', 'Pr(>|t|)', 'Effect size'))))
  
  names(cf) <- cf.test
  
  message('Testing differential methylation on each locus:')
  for (i in seq_len(n.cpg)) {
    
    if (i%%ceiling(n.cpg/80)==0) message('.', appendLF=FALSE);
    
    ## filter NA
    sample.is.na <- is.na(betas[i,])
    if (all(sample.is.na)) next;
    
    if (any(sample.is.na)) {      # this "if" improves performance
      design1 <- design[!sample.is.na,,drop=FALSE]
      design1.fac <- design.fac[!sample.is.na,]
      betas1 <- betas[i,!sample.is.na]
      
      if (sum(apply(
        design1[,2:n.cf,drop=FALSE], 2,
        function(x) length(unique(x)))==1) > 0) next;
      
      rdf <- unlist(lapply(design1.fac, function(x) max(min(
        tabulate(x)),1)), recursive = FALSE, use.names=FALSE)
      
    } else {
      design1 <- design
      design1.fac <- design.fac
      betas1 <- betas[i,]
      rdf <- rdf0
    }
    
    ## sigma is padded, boundary adjusted, and log2
    ## damped (should use beta distribution)
    ## wts <- rep(1, length(betas1))
    pad <- 0.01
    betas1.padded <- pmin(pmax(betas1, pad), 1-pad)
    wts <- 1/(betas1.padded*(1-betas1.padded)) # 1/var
    
    ## QR-solve weighted least square
    z <- .lm.fit(design1*wts, betas1*wts)
    names(z$coefficients) <- colnames(design1)
    p1 <- seq_len(z$rank)
    coefs <- z$coefficients[z$pivot[p1]]
    residuals <- z$residuals / wts
    
    if (balanced) {
      se <- sum(residuals^2)/(length(residuals)-1)
    } else {
      ## Welch-Satterthwaite se
      rs <- residuals^2
      se <- unlist(
        lapply(design1.fac, function(group1) {sqrt(sum(
          vapply(split(
            rs, group1), mean, numeric(1)), na.rm=TRUE))}),
        recursive = FALSE, use.names=FALSE)
    }
    
    ## lower bound coefficient se, this selects against
    ## high effect size differences
    se <- pmax(se, se.lb)
    
    ## t-statistics
    t.stat <- coefs / se
    pval <- 2*pt(abs(t.stat), rdf, lower.tail=FALSE)
    stopifnot(!any(is.na(pval)))
    
    ## output
    fitted.rg <- range(betas1 - residuals)
    eff <- fitted.rg[2] - fitted.rg[1]  # effect size
    ans <- cbind(coefs, se, t.stat, pval, eff)
    for (cfi in cf.test) {
      if (cfi %in% rownames(ans))
        cf[[cfi]][i,] <- ans[cfi,]
    }
    ## betas.fitted[i,!sample.is.na] <- betas1 - residuals
  }
  message('.\n', appendLF=FALSE)
  ## CpGs are correlated, should apply p-value adjustment on segments
  ## cf <- lapply(cf, function(cf1) cbind(cf1,
  ## P.adjusted=p.adjust(cf1[,'Pr(>|t|)'],method='BH')))
  
  message('Significant loci (p<0.05):')
  sigcnts <- lapply(cf, function(x) sum(x[,'Pr(>|t|)'] < 0.05, na.rm=TRUE))
  sigmsg <- lapply(seq_along(sigcnts), function(i) sprintf(
    ' - %s: %d significant loci.', names(sigcnts[i]), sigcnts[i][[1]]))
  message(do.call(paste0, list(sigmsg, collapse='\n')))
  
  cf
}

# End of file
