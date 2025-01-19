#### based on snpboost function, separate batch building for different parameters, 

snpboostlss <- function(data.dir, # data folder
                        genotype.pfile, # name of genotype file 
                        phenotype.file, # name of file containing the interested phenotype 
                        phenotype, # column name of interested phenotype
                        family = "gaussian", # currently only gaussian 
                        metric = "loss", # currently only loss function
                        covariates.mu = NULL, # name of covariates for mu
                        covariates.sigma = NULL, # name of covariates for sigma
                        split.col = NULL, # name of column indicating train/test/val grouping
                        mem = NULL, # memory for PLINK2.0 to harness for computation
                        give_residuals = FALSE, # TRUE if residuals are saved for each iteration. Attention: large files
                        configs = NULL, # other configurations
                        p_batch = 1000, # batch size
                        m_batch = 1000, # max no. of iterations per batch
                        b_max = 20000, # max no. of batches 
                        b_stop = 2, # stopping lag for outer stopping criteria
                        sl = NULL # if NULL, use adaptive step length; if a number specified, used as fixed step length for mu and sigma.
                        ){
  
  ####---------------------Preprocessing------------------------------------####
  ID <- ALT <- NULL
  genotype.pfile <- paste0(data.dir, genotype.pfile)
  phenotype.file <- paste0(data.dir, phenotype.file)
  
  time.start <- Sys.time()
  snpboostLogger('Start snpboostlss', log.time = time.start)
  
  snpboostLogger('Preprocessing start..')
  
  ### --- Read genotype IDs --- ###
  ids <- list(); phe <- list()
  ids[['psam']] <- readIDsFromPsam(paste0(genotype.pfile, '.psam'))
  
  ### --- combine the specified configs with the default values --- ###
  configs <- setupConfigs(configs, genotype.pfile, phenotype.file, phenotype, 
                          covariates.mu=covariates.mu, covariates.sigma=covariates.sigma, 
                          split.col=split.col)
  # setupConfigs can be updated as argument covariates seems to be useless
  
  ### --- Read phenotype file --- ###
  phe[['master']] <- suppressWarnings(readPheMaster(phenotype.file, ids[['psam']], family, covariates.mu, covariates.sigma, phenotype, status.col, split.col, configs))  
  rownames(phe[['master']]) <- phe[['master']]$ID
  # readPheMaster can be updated, some arguments are already in configs
  
  ### --- Define the set of individual IDs for training (and validation) set(s) --- ###
  if(is.null(split.col)){
    splits <- c('train')
    ids[['train']] <- phe[['master']]$ID
  }else{
    splits <- c('train', 'val') # can use levels in split.col instead of specifying manually
    for(s in splits){
      ids[[s]] <- phe[['master']]$ID[ phe[['master']][[split.col]] == s ]
    }
  }
  
  ### --- Prepare the feature matrix --- ###
  features.mu <- list()
  features.sigma <- list()
  for(s in splits){
    phe[[s]] <- phe[['master']][match(ids[[s]], phe[['master']]$ID), ]
    rownames(phe[[s]]) <- phe[[s]]$ID
    if (length(covariates.mu) > 0) {
      features.mu[[s]] <- phe[[s]][, covariates.mu, with = F]
    } else {
      features.mu[[s]] <- NULL
    }
    if (length(covariates.sigma) > 0) {
      features.sigma[[s]] <- phe[[s]][, covariates.sigma, with = F]
    } else {
      features.sigma[[s]] <- NULL
    }
    if(configs[['verbose']]) snpboostLogger(sprintf("The number of individuals in %s set: %d", s, dim(phe[[s]])[1]))
  }
  
  ### --- Prepare the response --- ###
  response <- list() ; pred <- list()
  for(s in splits){
    response[[s]] <- phe[[s]][[phenotype]]
  }
  
  ### --- Read genotypes --- ###
  vars <- dplyr::mutate(
    dplyr::rename(
      data.table::fread(file = paste0(genotype.pfile, '.pvar')),
      'CHROM'='#CHROM'),
    VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  configs[["excludeSNP"]] <- base::intersect(configs[["excludeSNP"]], vars)
  pvar <- pgenlibr::NewPvar(paste0(genotype.pfile, '.pvar'))
  pgen <- list()
  for(s in splits) pgen[[s]] <- pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), pvar=pvar, sample_subset=match(ids[[s]], ids[['psam']]))
  pgenlibr::ClosePvar(pvar)
  
  stats <- computeStats(genotype.pfile, phe[['train']]$ID, configs = configs)
  
  ### --- End --- ###
  snpboostLoggerTimeDiff("Preprocessing end.", time.start, indent=1)
  
  ####---------------------Initialization-----------------------------------####
  
  ### --- Iteration 0 --- ###
  time.init.start <- Sys.time()
  snpboostLogger("Iteration 0", log.time = time.init.start)
  
  mod <- gamlss(
    as.formula(paste(phenotype, " ~ ", paste(c(1, covariates.mu), collapse = " + "))),
    sigma.formula = as.formula(paste(" ~ ", paste(c(1, covariates.sigma), collapse = " + "))),
    family = NO(), data = phe[['train']])
  
  prediction.mu <- mod$mu.lp
  prediction.sigma <- mod$sigma.lp # linear predictor for log(sigma)
  residual <- computeResidualsLSS(phe[['train']][[phenotype]], prediction.mu, prediction.sigma)
  residual.mu <- matrix(residual$r_mu, ncol=1)
  residual.sigma <- matrix(residual$r_sigma, ncol=1)

  cov.mu_val <- model.matrix(as.formula(paste(phenotype, " ~ ", paste(c(1, covariates.mu), collapse = " + "))), data = phe[['val']])
  cov.sigma_val <- model.matrix(as.formula(paste(" ~ ", paste(c(1, covariates.sigma), collapse = " + "))), data = phe[['val']])
  prediction.mu_val <- Rfast::mat.mult(cov.mu_val, matrix(mod$mu.coefficients))
  prediction.sigma_val <- Rfast::mat.mult(cov.sigma_val, matrix(mod$sigma.coefficients))
  
  residual_val <- computeResidualsLSS(phe[['val']][[phenotype]], prediction.mu_val, prediction.sigma_val)
  residual.mu_val <- matrix(residual_val$r_mu, ncol=1)
  residual.sigma_val <- matrix(residual_val$r_sigma, ncol=1)
  
  rownames(residual.mu) <- rownames(phe[['train']])
  colnames(residual.mu) <- c('0')
  rownames(residual.sigma) <- rownames(phe[['train']])
  colnames(residual.sigma) <- c('0')
  rownames(residual.mu_val) <- rownames(phe[['val']])
  colnames(residual.mu_val) <- c('0')
  rownames(residual.sigma_val) <- rownames(phe[['val']])
  colnames(residual.sigma_val) <- c('0')
  
  snpboostLoggerTimeDiff("Iteration 0 end.", time.init.start, indent=1)
  
  ### --- create containers to store results --- ###
  
  # beta: store updated coefficients for mu in each iteration, 1st column is the updated SNP.
  beta <- data.frame(NA, t(mod$mu.coefficients), 0, row.names = NULL)
  names(beta) <- c("name", "intercept", covariates.mu, "variant")
  # gamma: store updated coefficients for sigma in each iteration, 1st column is the updated SNP.
  gamma <- data.frame(NA, t(mod$sigma.coefficients), 0, row.names = NULL)
  names(gamma) <- c("name", "intercept", covariates.sigma, "variant")
  
  # could be deleted. not used later.
  nobs <- nrow(phe[['train']])
  nvars <- length(vars)-length(stats[["excludeSNP"]])
  
  metric.train <- numeric()
  metric.val <- numeric()
  metric.train[1] <- computeLoss(phe[['train']][[phenotype]], prediction.mu, prediction.sigma)
  metric.val[1] <- computeLoss(phe[['val']][[phenotype]], prediction.mu_val, prediction.sigma_val)
  
  if(give_residuals){
    residual.track <- list(residual)
    residual.track.val <- list(residual_val)
  }
  
  if(is.null(sl)){
    adaptSL <- list()
  }
  
  m_batch_list <- numeric()
  
  snpboostLoggerTimeDiff("Initialization end.", time.init.start, indent=1)
  
  ####--------------------------Outer loop----------------------------------####
  
  for (iter in 1:b_max){
    
    time.iter.start <- Sys.time()
    snpboostLogger(paste0("Iteration ", iter), log.time=time.iter.start)
    
    ###----------Calculate correlations----------###
    if (configs[['verbose']]) snpboostLogger("  Start computing inner product ...")
    time.prod.start <- Sys.time()
    
    # compute correlations with residual.mu
    prod.full.mu <- computeProduct(residual.mu, genotype.pfile, vars, stats, configs, iter=iter) # excluded SNPs result in NA
    score.mu <- abs(prod.full.mu[, 1]) 
    rm(prod.full.mu)
    
    # compute correlations with residual.sigma
    prod.full.sigma <- computeProduct(residual.sigma, genotype.pfile, vars, stats, configs, iter=iter)
    score.sigma <- abs(prod.full.sigma[, 1]) # excluded SNPs result in NA
    rm(prod.full.sigma)
    
    if (configs[['verbose']]) snpboostLoggerTimeDiff("  End computing inner product.", time.prod.start)
    
    ###----------Create batches----------###
    if (configs[['verbose']]) snpboostLogger("Start updating feature matrix ...", indent=1)
    time.updateFeature.start <- Sys.time()
    
    sorted.score.mu <- sort(score.mu, decreasing = T, na.last = NA) # with na.last=NA, NA values in score.mu are removed
    c_stop.mu <- ifelse(length(sorted.score.mu)>p_batch, sorted.score.mu[p_batch+1], 0)
    sorted.score.sigma <- sort(score.sigma, decreasing = T, na.last = NA)
    c_stop.sigma <- ifelse(length(sorted.score.sigma)>p_batch, sorted.score.sigma[p_batch+1], 0)
    
    if (length(sorted.score.mu) > 0) {
      features.to.use.mu <- c(covariates.mu,names(sorted.score.mu)[1:min(p_batch, length(sorted.score.mu))])
      
      #features matrix set to NULL
      features.to.discard.mu <- setdiff(colnames(features.mu[['train']]),features.to.use.mu)
      
      if (length(features.to.discard.mu) > 0) {
        for(s in splits){features.mu[[s]][,features.to.discard.mu]= NULL}
      }
      
      features.to.add.mu <- setdiff(features.to.use.mu,colnames(features.mu[['train']]))
      
      if(length(features.to.add.mu)>0){
        for(s in splits){
          tmp.features.add.mu <- prepareFeatures(pgen[[s]], vars, features.to.add.mu, stats)
          if (!is.null(features.mu[[s]])) {
            features.mu[[s]][, colnames(tmp.features.add.mu) := tmp.features.add.mu]
          } else {
            features.mu[[s]] <- tmp.features.add.mu
          }
          rm(tmp.features.add.mu)
        }
      }
    } else {
      break
    }
    
    if (length(sorted.score.sigma) > 0) {
      features.to.use.sigma <- c(covariates.sigma,names(sorted.score.sigma)[1:min(p_batch, length(sorted.score.sigma))])
      
      #features matrix set to NULL
      features.to.discard.sigma <- setdiff(colnames(features.sigma[['train']]),features.to.use.sigma)
      
      if (length(features.to.discard.sigma) > 0) {
        for(s in splits){features.sigma[[s]][,features.to.discard.sigma]= NULL}
      }
      
      features.to.add.sigma <- setdiff(features.to.use.sigma,colnames(features.sigma[['train']]))
      
      if(length(features.to.add.sigma)>0){
        for(s in splits){
          tmp.features.add.sigma <- prepareFeatures(pgen[[s]], vars, features.to.add.sigma, stats)
          if (!is.null(features.sigma[[s]])) {
            features.sigma[[s]][, colnames(tmp.features.add.sigma) := tmp.features.add.sigma]
          } else {
            features.sigma[[s]] <- tmp.features.add.sigma
          }
          rm(tmp.features.add.sigma)
        }
      }
    } else {
      break
    }
    
    if (configs[['verbose']]) snpboostLoggerTimeDiff("End updating feature matrix.", time.updateFeature.start, indent=2)
    
    ###----------Inner loop: cyclical boosting----------###
    if (configs[['verbose']]) snpboostLogger("Start Boosting ...", indent=1)
    
    lmlss_boost_tmp <- lmlss_boost(r.mu = residual.mu,
                                   r.sigma = residual.sigma,
                                   x = features.mu[['train']],
                                   z = features.sigma[['train']],
                                   r.mu_val = residual.mu_val,
                                   r.sigma_val = residual.sigma_val, 
                                   x_val = features.mu[['val']],
                                   z_val = features.sigma[['val']],
                                   beta = beta, 
                                   gamma = gamma,
                                   m_batch = m_batch, 
                                   sl = sl, 
                                   covariates.mu = covariates.mu,
                                   covariates.sigma = covariates.sigma,
                                   give_residuals=give_residuals,
                                   c_stop.mu=c_stop.mu,
                                   c_stop.sigma=c_stop.sigma,
                                   configs=configs,
                                   phe = phe, 
                                   phenotype = phenotype, 
                                   prediction.mu = prediction.mu, 
                                   prediction.sigma = prediction.sigma, 
                                   prediction.mu_val = prediction.mu_val,
                                   prediction.sigma_val = prediction.sigma_val
                                   )
    
    ###----------Save results from inner loop----------###
    m_batch_list[iter] <- lmlss_boost_tmp$stop
    
    prediction.mu <- lmlss_boost_tmp$prediction.mu
    prediction.sigma <- lmlss_boost_tmp$prediction.sigma
    residual <- lmlss_boost_tmp$residual
    residual.mu <- matrix(residual$r_mu, ncol=1)
    residual.sigma <- matrix(residual$r_sigma, ncol=1)
    
    prediction.mu_val <- lmlss_boost_tmp$prediction.mu_val
    prediction.sigma_val <- lmlss_boost_tmp$prediction.sigma_val
    residual_val <- lmlss_boost_tmp$residual_val
    residual.mu_val <- matrix(residual_val$r_mu, ncol=1)
    residual.sigma_val <- matrix(residual_val$r_sigma, ncol=1) 
    
    rownames(residual.mu) <- rownames(phe[['train']])
    colnames(residual.mu) <- c(iter)
    rownames(residual.sigma) <- rownames(phe[['train']])
    colnames(residual.sigma) <- c(iter)
    rownames(residual.mu_val) <- rownames(phe[['val']])
    colnames(residual.mu_val) <- c(iter)
    rownames(residual.sigma_val) <- rownames(phe[['val']])
    colnames(residual.sigma_val) <- c(iter)
    
    metric.train <- c(metric.train, lmlss_boost_tmp$metric.train)
    metric.val <- c(metric.val, lmlss_boost_tmp$metric.val)
    
    beta <- lmlss_boost_tmp$beta
    gamma <- lmlss_boost_tmp$gamma
  
    if(give_residuals){
      residual.track <- c(residual.track, lmlss_boost_tmp$residual.track)
      residual.track.val <- c(residual.track.val, lmlss_boost_tmp$residual.track_val)
    }
    
    if(is.null(sl)){
      adaptSL[[iter]] <- list(adaptSL.mu=lmlss_boost_tmp$adaptSL.mu,
                              adaptSL.sigma=lmlss_boost_tmp$adaptSL.sigma)
    }
    
    rm(lmlss_boost_tmp)
    
    if (configs[['verbose']]) snpboostLogger("End Boosting.", indent=1)
    
    if (configs[['save']]) {
      save(metric.train, metric.val, beta, gamma, configs, residual, residual_val, m_batch_list, 
           prediction.mu, prediction.sigma, prediction.mu_val, prediction.sigma_val,features.mu,features.sigma,
           file = file.path(configs[['results.dir']], configs[["save.dir"]], paste0("output_iter_", iter, ".RData")))
    }
    
    time.iter.end <- Sys.time()
    snpboostLoggerTimeDiff(paste0("End iteration ", iter, '.'), time.iter.start, time.iter.end, indent=1)
    snpboostLoggerTimeDiff("The total time since start.", time.start, time.iter.end, indent=2)
    
    ###----------check stopping criteria for outer loop----------###
    if(checkEarlyStoppingBatches(metric.val,iter,b_stop,configs,m_batch_list)) break
    
  }
  
  
  ####--------------------------Determine mstop-----------------------------####
  # record number of batches gone through
  n_iter = iter
  # record the iteration number which minimize metric.val
  n_step = which.min(metric.val)[1] 
  
  snpboostLoggerTimeDiff("End snpboostlss.", time.start)
  if(! configs[['save']]) cleanUpIntermediateFiles(configs)
  if(configs[['verbose']]) print(gc())

  out <- list(metric.train = metric.train,
              metric.val = metric.val,
              beta = beta,
              gamma = gamma,
              configs = configs,
              stats = stats,
              residual.track = if(give_residuals){residual.track}else{NULL},
              residual.track.val = if(give_residuals){residual.track.val}else{NULL},
              adaptSL = if(is.null(sl)){adaptSL}else{NULL},
              n_iter= n_iter, 
              n_step = n_step, 
              m_batch_list=m_batch_list)
  out
  
}