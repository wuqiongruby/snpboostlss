### functions from snpnet package that are used for snpboost and are not adapted

cat_or_zcat <- function(filename, configs=list(zstdcat.path='zstdcat', zcat.path='zcat')){
  if(stringr::str_ends(basename(filename), '.zst')){
    return(configs[['zstdcat.path']])
  }else if(stringr::str_ends(basename(filename), '.gz')){
    return(configs[['zcat.path']])
  }else{
    return('cat')
  }
}

readIDsFromPsam <- function(psam){
  FID <- IID <- NULL  # to deal with "no visible binding for global variable"
  df <- data.table::fread(psam) %>%
    dplyr::rename('FID' = '#FID') %>%
    dplyr::mutate(ID = paste(FID, IID, sep='_'))
  df$ID
}


####try to import from snpnet_compact

checkMissingPhenoWarning <- function(phe.master, phe.no.missing.IDs){
  # Show warning message if there are individuals (in phe file)
  # that have (genotype or phenotype) missing values.
  phe.missing.IDs <- phe.master$ID[ ! phe.master$ID %in% phe.no.missing.IDs ]
  if(length(phe.missing.IDs) > 0){
    warning(sprintf(
      'We detected missing values for %d individuals (%s ...).\n',
      length(phe.missing.IDs),
      paste(utils::head(phe.missing.IDs, 5), collapse=", ")
    ))
  }
}

readPlinkKeepFile <- function(keep_file){
  ID <- NULL  # to deal with "no visible binding for global variable"
  keep_df <- data.table::fread(keep_file, colClasses='character', stringsAsFactors=F)
  keep_df$ID <- paste(keep_df$V1, keep_df$V2, sep = "_")
  keep_df %>% dplyr::pull(ID)
}

inferFamily <- function(phe, phenotype, status){
  if (all(unique(phe[[phenotype]] %in% c(0, 1, 2, -9)))) {
    family <- "binomial"
  } else if(!is.null(status) && (status %in% colnames(phe))) {
    family <- "cox"
  } else {
    family <- "gaussian"
  }
  family
}

updateConfigsWithFamily <- function(configs, family){
  out <- configs
  out[['family']] <- family
  if (is.null(out[['metric']])) out[['metric']] <- setDefaultMetric(family)
  out
}

readBinMat <- function(fhead, configs){
  # This is a helper function to read binary matrix file (from plink2 --variant-score bin)
  rows <- data.table::fread(file = paste0(fhead, '.vars'), head=F)$V1
  cols <- data.table::fread(paste0(fhead, '.cols'), head=F)$V1
  bin.reader <- file(paste0(fhead, '.bin'), 'rb')
  M = matrix(
    readBin(bin.reader, 'double', n=length(rows)*length(cols), endian = configs[['endian']]),
    nrow=length(rows), ncol=length(cols), byrow = T
  )
  close(bin.reader)
  colnames(M) <- cols
  rownames(M) <- rows
  if (! configs[['save.computeProduct']]) system(paste(
    'rm', paste0(fhead, '.cols'), paste0(fhead, '.vars.zst'),
    paste0(fhead, '.bin'), sep=' '
  ), intern=F, wait=T)
  
  M
}

cleanUpIntermediateFiles <- function(configs){
  for(subdir in c(configs[["save.dir"]], configs[["meta.dir"]])){
    system(paste(
      'rm', '-rf', file.path(configs[['results.dir']], subdir), sep=' '
    ), intern=F, wait=T)
  }
}