rm(list=ls())

library(data.table)
library(tidyverse)
library(parallel)
library(bigsnpr)
library(cowplot)
library(readr)
library(dataPreparation)
library(pgenlibr)
library(Rfast)
library(gamlss)

### load snpboostlss functions 
code.dir <- "C:/RESEARCH/project1/snpboostlss/snpboostlss/"
source(paste0(code.dir, "snpboostlss.R"))
source(paste0(code.dir, "lmlss_boost.R"))
source(paste0(code.dir, "functions_snpboostlss.R"))
source(paste0(code.dir, "functions_snpboost.R"))
source(paste0(code.dir, "functions_snpnet.R"))

data.dir <- "C:/RESEARCH/project1/data/example_data/"
genotype.pfile <- "genotypes_Qiong" ### genotype file, file formats: .pgen .pvar .psam
phenotype.file <- "phenotypes_Qiong.phe" ### phenotype file, .phe, tabular data
phenotype <- "LDL_adj" ### name of phenotype column in phenotype file
covariates.mu <- c("age_at_2015", "sex") ### name of covariates for mu
covariates.sigma <- c("age_at_2015") ### name of covariates for sigma

configs <- list(plink2.path = "C:/RESEARCH/project1/prep/plink2/plink2", ## path to plink2
                zcat.path = "cat", ## cat is already in global environment
                sed.path = "C:/RESEARCH/project1/prep/GnuWin32/bin/sed", 
                ## If using OS other than windows, sed is automatically installed, no need to specify this path.
                ## For windows, need to install sed somewhere and spedify its path here.
                results.dir = paste0("C:/RESEARCH/project1/analysis/snpboostlss_trial_run/test_",phenotype), ## results folder
                missing.rate = 0.1, # QC criteria: genotyping rate per subject should be larger than 0.9.
                MAF.thresh = 0.001,  # QC criteria: MAF per SNP should be larger than 0.001.
                verbose = TRUE, # Default=FALSE. If TRUE, progress messages will be printed.
                mem=16000, # memory for PLINK2.0 to harness for computation
                nCores=16 # number of CPUs that you want to use
)

###fit fit_snpboost
time_snpboost_start <- Sys.time()

fit_snpboostlss <- snpboostlss(
  data.dir = data.dir,
  genotype.pfile = genotype.pfile,
  phenotype.file = phenotype.file,
  phenotype = phenotype,
  covariates.mu = covariates.mu,
  covariates.sigma = covariates.sigma,
  configs = configs,
  split.col = "split", ### name of column in phenotype file indicating split
  p_batch = 1000, # batch size
  m_batch = 1000, # maximum number of boosting steps per batch
  b_max = 20000, # maximum number of batches
  b_stop = 2, # stop if performance not improving for more than 2 batches
  sl= 0.1, # learning rate (default=NULL)
  give_residuals = FALSE, # TRUE if residuals should be saved for each boosting step, attention: large files!
  family = "gaussian", # family argument: currently only gaussian
  metric = "loss" # loss function/metric: currently only negative log-likelihood for gaussian
) 

time_snpboost_end <- Sys.time()
time_snpboost_end - time_snpboost_start


# extract coefficients
mstop <- fit_snpboostlss$n_step
betas <- get_coefficients(fit_snpboostlss$beta,
                          mstop,
                          covariates=fit_snpboostlss$configs[['covariates']])
gammas <- get_coefficients(fit_snpboostlss$gamma,
                          mstop,
                          covariates=fit_snpboostlss$configs[['covariates']])

# prediction on test set
pred <- predict_snpboostlss(
  fit_snpboostlss,
  new_genotype_file=paste0(data.dir,genotype.pfile),
  new_phenotype_file=paste0(data.dir,phenotype.file),
  phenotype=phenotype,
  subset="test")
