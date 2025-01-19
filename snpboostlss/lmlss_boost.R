lmlss_boost <- function(r.mu, r.sigma, x, z,
                        r.mu_val=NULL,
                        r.sigma_val=NULL,
                        x_val=NULL,
                        z_val=NULL,
                        beta=data.frame(name=NA,intercept=0,variant=NA),
                        gamma=data.frame(name=NA,intercept=0,variant=NA),
                        m_batch=1,
                        sl=NULL, # if NULL, use adaptive step length; if a number specified, used as fixed step length for mu and sigma.
                        covariates.mu=NULL,
                        covariates.sigma=NULL,
                        give_residuals=FALSE,
                        c_stop.mu=0,
                        c_stop.sigma=0,
                        configs=NULL,
                        phe=NULL,
                        phenotype=NULL,
                        prediction.mu=NULL,
                        prediction.sigma=NULL,
                        prediction.mu_val=NULL,
                        prediction.sigma_val=NULL
                        ){
  
  
  ###----------Initialization----------###
  p.mu <- ncol(x)-length(covariates.mu)
  p.sigma <- ncol(z)-length(covariates.sigma)
  nobs <- length(r.mu)
  variants.mu <- colnames(x)[!colnames(x) %in% covariates.mu]
  variants.sigma <- colnames(z)[!colnames(z) %in% covariates.sigma]
  metric.train <- numeric(m_batch)
  pred.mu <- prediction.mu
  pred.sigma <- prediction.sigma
  
  previous.beta <- get_coefficients(coeff_path=beta,covariates=covariates.mu)
  previous.beta[variants.mu[!variants.mu %in% names(previous.beta)]]=0
  previous.beta <- previous.beta[c("(Intercept)",covariates.mu,variants.mu)]
  b <- beta %>% add_row(name=rep(NA,m_batch),intercept=rep(0,m_batch),variant=rep(0,m_batch))
  
  x_std <- as.matrix(dataPreparation::fast_scale(x%>%select(all_of(variants.mu)),verbose=F))
  X <- as.matrix(cbind(rep(1,nrow(x)), x))
  colnames(X)=c("(Intercept)",colnames(x))
  
  previous.gamma <- get_coefficients(coeff_path=gamma,covariates=covariates.sigma)
  previous.gamma[variants.sigma[!variants.sigma %in% names(previous.gamma)]]=0
  previous.gamma <- previous.gamma[c("(Intercept)",covariates.sigma,variants.sigma)]
  g <- gamma %>% add_row(name=rep(NA,m_batch),intercept=rep(0,m_batch),variant=rep(0,m_batch))
  
  z_std <- as.matrix(dataPreparation::fast_scale(z%>%select(all_of(variants.sigma)),verbose=F))
  Z <- as.matrix(cbind(rep(1,nrow(z)), z))
  colnames(Z)=c("(Intercept)",colnames(z))
  
  # for validation set
  metric.val <- numeric()
  pred.mu_val <- prediction.mu_val
  pred.sigma_val <- prediction.sigma_val
  X_val <- as.matrix(cbind(rep(1,nrow(x_val)),x_val))
  Z_val <- as.matrix(cbind(rep(1,nrow(z_val)),z_val))
  
  colnames(X_val)=c("(Intercept)",colnames(x_val))
  colnames(Z_val)=c("(Intercept)",colnames(z_val))
  
  if(give_residuals){
    residual.track <- vector(mode = "list", length = m_batch)
    residual.track_val <- vector(mode = "list", length = m_batch)
  }
  
  # if no input of fixed sl, create vectors to store adaptive sl for mu and sigma.
  if(is.null(sl)){
    adaptSL.mu <- numeric()
    adaptSL.sigma <- numeric()
  }
  
  # set early stopping flag for mu
  earlyStop.mu <- FALSE
  # set early stopping flag for sigma
  earlyStop.sigma <- FALSE
  
  
  ###----------Boosting----------###
  for(m in 1:m_batch){

    ##---Determine whether enter inner loop---##
    if(earlyStop.mu*earlyStop.sigma) break
    
    
    ##---Boosting for mu---###
    if(earlyStop.mu){
      # no need to update preivous.beta, pred.mu, pred.sigma, r.mu, r.sigma
      # save new beta as the same values of previous iteration
      b[nrow(beta)+m,] <- b[nrow(beta)+m-1,]
      # save new residuals as the same values of previous iteration
      if(give_residuals){
        residuals.track[[m]] <- list(r_mu = r.mu, r_sigma = r.sigma)
      }
    }else{
      # choose the variant with highest correlation with current residual
      if(p.mu==1){
        j_star <- 1
      }else{
        score.inbatch.mu <- abs(Rfast::Crossprod(x_std, matrix((r.mu-mean(r.mu))/sd(r.mu),ncol=1)))
        j_star <- which.max(score.inbatch.mu)
      }
      
      # Determine whether mu should be early stopped
      if(score.inbatch.mu[j_star]/(nobs-1) < c_stop.mu){
        # set early stopping flag for mu as TRUE
        earlyStop.mu <- TRUE
        # no need to update previous.beta, pred.mu, pred.sigma, r.mu, r.sigma
        # save new beta as the same values of previous iteration
        b[nrow(beta)+m,] <- b[nrow(beta)+m-1,]
        # save new residuals as the same values of previous iteration
        if(give_residuals){
          residuals.track[[m]] <- list(r_mu = r.mu, r_sigma = r.sigma)
        }
      }else{
        # fit linear model
        lm_tmp.mu <- .lm.fit(x=X[,c("(Intercept)",covariates.mu,variants.mu[j_star])],y=r.mu)
        # determine step length if adaptive sl is used
        if (is.null(sl)){
          # see Zhang et al. (2022) on adaptive length for boosting location-scale regression
          sl.mu <- 0.1 * (sum((as.vector(r.mu - lm_tmp.mu$residuals))^2) / sum((as.vector(r.mu - lm_tmp.mu$residuals))^2/exp(2*pred.sigma)))
          adaptSL.mu <- append(adaptSL.mu, sl.mu)
        }else{
          sl.mu <- sl
        }
        # update previous.beta
        previous.beta[c("(Intercept)",covariates.mu,variants.mu[j_star])] <- previous.beta[c("(Intercept)",covariates.mu,variants.mu[j_star])]+sl.mu*coef(lm_tmp.mu)
        # save new intercept and fitted coefficient into b
        b[nrow(beta)+m,1] <- variants.mu[j_star]
        b[nrow(beta)+m,-1] <- c(previous.beta[c("(Intercept)",covariates.mu)],previous.beta[variants.mu[j_star]])
        # update prediction
        pred.mu <- pred.mu + sl.mu * as.vector(r.mu - lm_tmp.mu$residuals)
        # update residuals
        r <- computeResidualsLSS(phe[['train']][[phenotype]], pred.mu, pred.sigma)
        r.mu <- matrix(r$r_mu, ncol=1)
        r.sigma <- matrix(r$r_sigma, ncol=1)
        
        # update pred.mu_val
        pred.mu_val <- pred.mu_val + sl.mu * Rfast::mat.mult(X_val[,c("(Intercept)", covariates.mu, variants.mu[j_star])], matrix(coef(lm_tmp.mu)))
        
        if(give_residuals){
          residual.track[[m]] = r
        }
        # no need to update metric (loss) here, update after boosting gamma 
        # no need to update validation set here, update after boosting gamma
        
        rm(lm_tmp.mu)
      }
    }
    
    
    ##---Boosting for sigma---##
    if(earlyStop.sigma){
      # no need to update preivous.beta, pred.mu, pred.sigma, r.mu, r.sigma
      # save new gamma as the same values of previous iteration
      g[nrow(gamma)+m,] <- g[nrow(gamma)+m-1,]
      # no need to save residuals, since the m-th item in residuals.track has already been saved when boosting mu.
    }else{
      # choose the variant with highest correlation with current residual
      if(p.sigma==1){
        h_star <- 1
      }else{
        score.inbatch.sigma <- abs(Rfast::Crossprod(z_std,matrix((r.sigma-mean(r.sigma))/sd(r.sigma),ncol=1)))
        h_star <- which.max(score.inbatch.sigma)
      }
      
      # Determine whether sigma should be early stopped
      if(score.inbatch.sigma[h_star]/(nobs-1) < c_stop.sigma){
        # set early stopping flag for sigma as TRUE
        earlyStop.sigma <- TRUE
        # no need to update previous.beta, pred.mu, pred.sigma, pred_tmp.mu, pred_tmp.sigma, r.mu, r.sigma
        # save new gamma as the same values of previous iteration
        g[nrow(gamma)+m,] <- g[nrow(gamma)+m-1,]
        # no need to save residuals, since the m-th item in residuals.track has already been saved when boosting mu.
      }else{
        # fit linear model
        lm_tmp.sigma <- .lm.fit(x=Z[,c("(Intercept)",covariates.sigma,variants.sigma[h_star])],y=r.sigma)
        # determine step length if adaptive sl is used
        if (is.null(sl)){
          # see Zhang et al. (2022) on adaptive length for boosting location-scale regression
          sl.sigma <- 0.05
          adaptSL.sigma <- append(adaptSL.sigma, sl.sigma)
        }else{
          sl.sigma <- sl
        }
        # update previous.gamma
        previous.gamma[c("(Intercept)",covariates.sigma,variants.sigma[h_star])] <- previous.gamma[c("(Intercept)",covariates.sigma,variants.sigma[h_star])]+sl.sigma*coef(lm_tmp.sigma)
        # save new intercept and fitted coefficient into g
        g[nrow(gamma)+m,1] <- variants.sigma[h_star]
        g[nrow(gamma)+m,-1] <- c(previous.gamma[c("(Intercept)",covariates.sigma)],previous.gamma[variants.sigma[h_star]])
        # update prediction
        pred.sigma <- pred.sigma + sl.sigma * as.vector(r.sigma - lm_tmp.sigma$residuals)
        # update residuals
        r <- computeResidualsLSS(phe[['train']][[phenotype]], pred.mu, pred.sigma)
        r.mu <- matrix(r$r_mu, ncol=1)
        r.sigma <- matrix(r$r_sigma, ncol=1)
        
        # update pred.sigma_val
        pred.sigma_val <- pred.sigma_val + sl.sigma * Rfast::mat.mult(Z_val[,c("(Intercept)", covariates.sigma, variants.sigma[h_star])], matrix(coef(lm_tmp.sigma)))
        
        if(give_residuals){
          residual.track[[m]] = r
        }
        
        rm(lm_tmp.sigma)
        
      }
    }
    
    ##--- update metrics ---##
    metric.train[m] <- computeLoss(phe[['train']][[phenotype]], pred.mu, pred.sigma)
    metric.val[m] <- computeLoss(phe[['val']][[phenotype]], pred.mu_val, pred.sigma_val)
    
    ##--- update the variables on validation set ---##
    r_val <- computeResidualsLSS(phe[['val']][[phenotype]], pred.mu_val, pred.sigma_val)
    if(give_residuals){
      residual.track_val[[m]] = r_val
    }
    
  }
  
  ###----------Outputs----------###
  stop= length(metric.val)
  b <- b[1:(nrow(beta)+stop),]
  g <- g[1:(nrow(gamma)+stop),]
  metric.train <- metric.train[1:stop]
  if(give_residuals){
    residual.track <- residual.track[1:stop]
    residual.track_val <- residual.track_val[1:stop]
  }
  
  out <- list(beta=b, 
              gamma=g, 
              residual=r,
              residual_val=r_val,
              metric.train=metric.train,
              metric.val=metric.val,
              residual.track=if(give_residuals){residual.track}else{NULL}, 
              residual.track_val=if(give_residuals){residual.track_val}else{NULL},
              adaptSL.mu=if(is.null(sl)){adaptSL.mu}else{NULL},
              adaptSL.sigma=if(is.null(sl)){adaptSL.sigma}else{NULL},
              stop=stop, 
              prediction.mu = pred.mu,
              prediction.sigma = pred.sigma,
              prediction.mu_val = pred.mu_val,
              prediction.sigma_val = pred.sigma_val,
              earlyStop.mu = earlyStop.mu,
              earlyStop.sigma = earlyStop.sigma
              )
  
  
  return(out)
  
}