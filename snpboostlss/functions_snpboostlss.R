###---------- functions newly implemented for snpboostlss--------------------###
computeLoss <- function (phenotype, prediction.mu, prediction.sigma) {
  rho <- sum(prediction.sigma + (phenotype - prediction.mu)^2 / (2*exp(2*prediction.sigma)))
  return(rho)                       
}



computeResidualsLSS <- function(phenotype, prediction.mu, prediction.sigma){
  r_mu <- (phenotype-prediction.mu) / exp(2*prediction.sigma)
  r_sigma <- (phenotype-prediction.mu)^2 / exp(2*prediction.sigma) - 1
  list(r_mu=r_mu, r_sigma=r_sigma)
}


###---------- functions from snpboost but adapted for snpboostlss------------###

## setupConfigs: default value for standardize.variant = TRUE
setupConfigs <- function(configs, genotype.pfile, phenotype.file, phenotype, covariates.mu, covariates.sigma, split.col, status.col, mem){
  out.args <- as.list(environment())
  defaults <- list(
    missing.rate = 0.1,
    MAF.thresh = 0.001,
    nCores = 1,
    glmnet.thresh = 1e-07,
    num.snps.batch = 1000,
    vzs=FALSE, # no zipped file will be produced by PLINK
    increase.size = NULL,
    standardize.variant = TRUE,
    early.stopping = TRUE,
    stopping.lag = 2,
    niter = 50,
    keep = NULL,
    KKT.verbose = FALSE,
    save = FALSE,
    save.computeProduct = FALSE,
    prevIter = 0,
    results.dir = NULL,
    meta.dir = 'meta',
    save.dir = 'results',
    verbose = FALSE,
    KKT.check.aggressive.experimental = FALSE,
    gcount.basename.prefix = 'snpboost.train',
    gcount.full.prefix=NULL,
    endian="little",
    metric=NULL,
    plink2.path='plink2',
    zstdcat.path='zstdcat',
    zcat.path='zcat',
    rank = TRUE,
    excludeSNP = NULL,
    quantile = 0.5
  )
  out <- defaults
  
  # store additional params
  for (name in setdiff(names(out.args), "configs")) {
    out[[name]] <- out.args[[name]]
  }
  
  # update the defaults with the specified parameters and keep redundant parameters from configs
  for (name in names(configs)) {
    out[[name]] <- configs[[name]]
  }
  
  # update settings
  out[["early.stopping"]] <- ifelse(out[["early.stopping"]], out[['stopping.lag']], -1)
  if(is.null(out[['increase.size']]))  out[['increase.size']] <- out[['num.snps.batch']]/2
  
  # configure the temp file locations
  #   We will write some intermediate files to meta.dir and save.dir.
  #   those files will be deleted with snpnet::cleanUpIntermediateFiles() function.
  if (is.null(out[['results.dir']])) out[['results.dir']] <- tempdir(check = TRUE)
  dir.create(file.path(out[['results.dir']], out[["meta.dir"]]), showWarnings = FALSE, recursive = T)
  dir.create(file.path(out[['results.dir']], out[["save.dir"]]), showWarnings = FALSE, recursive = T)
  if(is.null(out[['gcount.full.prefix']])) out[['gcount.full.prefix']] <- file.path(
    out[['results.dir']], out[["meta.dir"]], out['gcount.basename.prefix']
  )
  
  out
}

# function for mu prediction
predict_snpboostlss <- function(fit, new_genotype_file, new_phenotype_file, phenotype, what=c("mu", "sigma"), subset=NULL,idx=NULL){
  if(is.null(idx)){
    if(!is.null(fit$metric.val)){
      idx <- which.min(fit$metric.val)
    }else{
      idx <- which.min(fit$metric.train)
    }
  }
  
  final_betas <- get_coefficients(fit$beta, idx, covariates=fit$configs[['covariates.mu']])
  final_gammas <- get_coefficients(fit$gamma, idx, covariates=fit$configs[['covariates.sigma']])
  
  chosen_SNPs_beta <- cbind(t(mapply(split_last, names(final_betas)[-(1:(length(fit$configs[['covariates.mu']])+1))], MoreArgs = list(pattern = "_"))),
                       final_betas[-(1:(length(fit$configs[['covariates.mu']])+1))])
  rownames(chosen_SNPs_beta)=NULL
  fwrite(chosen_SNPs_beta, paste0(fit$configs['results.dir'], "/chosen_SNPs_predict_beta.txt"),sep="\t")
  
  chosen_SNPs_gamma <- cbind(t(mapply(split_last, names(final_gammas)[-(1:(length(fit$configs[['covariates.sigma']])+1))], MoreArgs = list(pattern = "_"))),
                            final_gammas[-(1:(length(fit$configs[['covariates.sigma']])+1))])
  rownames(chosen_SNPs_gamma)=NULL
  fwrite(chosen_SNPs_gamma, paste0(fit$configs['results.dir'], "/chosen_SNPs_predict_gamma.txt"),sep="\t")
  
  phe <- fread(new_phenotype_file)
  if(!is.null(subset)){
    phe <- phe[phe[[fit$configs[['split.col']]]]==subset,]
  }
  
  if(("mu" %in% what) & nrow(chosen_SNPs_beta)>0){
    plink2_cmd <- paste(fit$configs['plink2.path'],"--pfile",new_genotype_file,"--score",
                        paste0(fit$configs['results.dir'],"/chosen_SNPs_predict_beta.txt"),1,2,3,"header cols=maybefid,denom,dosagesum,scoreavgs,scoresums","--out",paste0(fit$configs['results.dir'],"/PRS_beta"))
    system(plink2_cmd, intern=F, wait=T)
    
    PRS <- fread(paste0(fit$configs['results.dir'],"/PRS_beta.sscore")) %>% rename(PRS_mu = SCORE1_SUM)
    
    phe_PRS <- left_join(phe,PRS,by="IID")
  }else{
    phe_PRS <- phe %>% mutate(PRS_mu=0) 
  }
  pred.mu <- matrix(final_betas[1] + phe_PRS$PRS_mu + if(!is.null(fit$configs[['covariates.mu']])){as.matrix(phe_PRS%>%select(fit$configs[['covariates.mu']])) %*% final_betas[fit$configs[['covariates.mu']]]}else{0},ncol=1)
  rownames(pred.mu)=phe_PRS$IID
  
  if(("sigma" %in% what) & nrow(chosen_SNPs_gamma)>0){
    plink2_cmd <- paste(fit$configs['plink2.path'],"--pfile",new_genotype_file,"--score",
                        paste0(fit$configs['results.dir'],"/chosen_SNPs_predict_gamma.txt"),1,2,3,"header cols=maybefid,denom,dosagesum,scoreavgs,scoresums","--out",paste0(fit$configs['results.dir'],"/PRS_gamma"))
    system(plink2_cmd, intern=F, wait=T)
    
    PRS <- fread(paste0(fit$configs['results.dir'],"/PRS_gamma.sscore")) %>% rename(PRS_sigma = SCORE1_SUM)
    
    phe_PRS <- left_join(phe,PRS,by="IID")
  }else{
    phe_PRS <- phe %>% mutate(PRS_sigma=0) 
  }
  
  pred.sigma <- matrix(final_gammas[1] + phe_PRS$PRS_sigma + if(!is.null(fit$configs[['covariates.sigma']])){as.matrix(phe_PRS%>%select(fit$configs[['covariates.sigma']])) %*% final_gammas[fit$configs[['covariates.sigma']]]}else{0},ncol=1)
  rownames(pred.sigma)=phe_PRS$IID
  
  # if(fit$configs[['family']]=='cox'){
  #   surv <- survival::Surv(phe_PRS[[phenotype]], phe_PRS[[fit$configs[['status.col']]]])
  #   T_order <- order(surv[,1])
  #   IPCw <- IPCweights(surv)
  #   sum_help <- function(i){ifelse(i<length(surv[,1]),IPCw[T_order][i]^2*sum((surv[,1][T_order][(i+1):length(surv)]-surv[,1][T_order][i])!=0)+
  #                                    0.5*IPCw[T_order][i]^2*(sum((surv[,1][T_order]-surv[,1][T_order][i])==0)-1),0.5*IPCw[T_order][i]^2*(sum((surv[,1][T_order]-surv[,1][T_order][i])==0)-1))}
  #   
  #   w_denom <- sum(mapply(sum_help,1:length(surv)))
  #   wweights <- IPCw^2/w_denom
  #   
  #   residuals <- computeResiduals(surv, pred, fit$configs, IPCw, w_denom, T_order, wweights, fit$sigma)
  #   rownames(residuals) <- phe_PRS$IID
  #   
  #   metric <- computeMetric(residuals, surv, pred, fit$configs, surv, IPCw, fit$sigma, T_order)
  # }else{
  #   residuals <- computeResidualsLSS(phe_PRS[[phenotype]], pred, fit$configs)
  #   rownames(residuals) <- phe_PRS$IID
  #   
  #   metric <- computeLoss(residuals, phe_PRS[[phenotype]], pred, fit$configs, sigma = fit$sigma)
  # }
  
  # if(fit$configs[['family']]=='gaussian'){
  #   Rsquare <- cor(phe_PRS[[phenotype]],pred)^2
  # }else if(fit$configs[['family']]=='binomial'){
  #   AUC <- as.numeric(pROC::auc(phe_PRS[[phenotype]],exp(pred)/(1+exp(pred))))
  # }
  # out <- list(prediction = pred, residuals=residuals, metric = metric,
  #             Rsquare=if(fit$configs[['family']]=='gaussian'){Rsquare}else{NULL}, AUC=if(fit$configs[['family']]=='binomial'){AUC}else{NULL})

  out <- list(phenotype = phe_PRS[[phenotype]],
              pred_mu = if("mu" %in% what){pred.mu}else{NULL},
              pred_sigma = if("sigma" %in% what){pred.sigma}else{NULL},
              metric = computeLoss(phe_PRS[[phenotype]], pred.mu, pred.sigma))
  out
}

###---------- functions from snpnet but adapted for snpboostlss -------------###
readPheMaster <- function(phenotype.file, psam.ids, family, covariates.mu, covariates.sigma, phenotype, status, split.col, configs){
  
  covariates <- union(covariates.mu, covariates.sigma)    
  sort_order <- . <- ID <- NULL  # to deal with "no visible binding for global variable"
  
  if(!is.null(family) && family == 'cox'){
    selectCols <- c("FID", "IID", covariates, phenotype, status, split.col)
  } else{
    selectCols <- c("FID", "IID", covariates, phenotype, split.col)
  }
  
  phe.master.unsorted <- data.table::fread(
    cmd=if(.Platform$OS.type == "windows"){
      paste("powershell -command \"", cat_or_zcat(phenotype.file, configs), phenotype.file, ' | ', configs[['sed.path']], '-e s/^#//g \"')
      # sed is to edit file content.
      ## For Windows, need to specify the location of sed.
    }else{
      paste(cat_or_zcat(phenotype.file, configs), phenotype.file, ' | sed -e s/^#//g ') 
      ## For other OS, no need to specify the location of sed.
    },
    colClasses = c("FID" = "character", "IID" = "character"), select = selectCols
  )
  phe.master.unsorted$ID <- paste(phe.master.unsorted$FID, phe.master.unsorted$IID, sep = "_")
  
  # make sure the phe.master has the same individual ordering as in the genotype data
  # so that we don't have error when opening pgen file with sample subset option.
  phe.master <- phe.master.unsorted %>%
    dplyr::left_join(
      data.frame(ID = psam.ids, stringsAsFactors=F) %>%
        dplyr::mutate(sort_order = 1:n()),
      by='ID'
    ) %>%
    dplyr::arrange(sort_order) %>% dplyr::select(-sort_order) %>%
    data.table::as.data.table()
  rownames(phe.master) <- phe.master$ID
  
  
  for (name in c(covariates, phenotype)) {
    set(phe.master, i = which(phe.master[[name]] == -9), j = name, value = NA) # missing phenotypes are encoded with -9
  }
  
  # focus on individuals with complete covariates values
  if (is.null(covariates)) {
    phe.no.missing <- phe.master
  } else {
    phe.no.missing <- phe.master %>%
      dplyr::filter_at(dplyr::vars(covariates), dplyr::all_vars(!is.na(.)))
  }
  
  # focus on individuals with at least one observed phenotype values
  phe.no.missing <- phe.no.missing %>%
    dplyr::filter_at(dplyr::vars(phenotype), dplyr::any_vars(!is.na(.))) %>%
    dplyr::filter(ID %in% psam.ids) # check if we have genotype
  
  phe.no.missing.IDs <- phe.no.missing$ID
  
  if(!is.null(split.col)){
    # focus on individuals in training and validation set
    phe.no.missing.IDs <- intersect(
      phe.no.missing.IDs,
      phe.master$ID[ (phe.master[[split.col]] %in% c('train', 'val', 'test')) ]
    )
  }
  if(!is.null(configs[['keep']])){
    # focus on individuals in the specified keep file
    phe.no.missing.IDs <- intersect(phe.no.missing.IDs, readPlinkKeepFile(configs[['keep']]))
  }
  checkMissingPhenoWarning(phe.master, phe.no.missing.IDs)
  
  phe.master[ phe.master$ID %in% phe.no.missing.IDs, ]
}

