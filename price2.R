#' price2 : Alternative decomposition of Price equation
#' 
#' Last updated: 09 August 2018
#' 
#' @param mat a species-by-community biomass matrix
#' @param standardize whether output should be scaled by maximum observed value
#' @baseline whether the best site should be used as the baseline, or a user-supplied rownumber
#' 
#' @return a data.frame with components
#' 
price2 <- function(mat, baseline = "best", standardize = TRUE) {
  
  mat <- as.matrix(mat)
  
  if(is.null(rownames(mat))) rownames(mat) <- 1:nrow(mat)
  
  # remove colSums == 0
  mat <- mat[, colSums(mat) != 0]
  
  # identify the baseline community
  if(baseline == "best") b <- which.max(rowSums(mat)) else b <- as.numeric(baseline)
  
  # get vectors of baseline biomass
  base <- mat[b, ]
  
  # loop over comparisons(
  ret <- do.call(rbind, lapply((1:nrow(mat))[-b], function(j) {
    
    # get vectors of comparison biomass
    comp <- mat[j, ]
    
    # get vector of shared species
    shared <- apply(rbind(base, comp), 2, function(x) ifelse(all(x > 0), 1, 0) )
    
    # number of species in common
    sc <- sum(shared)
    
    # number of species unique to baseline
    suB <- sum(base[!shared] > 0)
    
    # number of species unique to comparison
    suD <- sum(comp[!shared] > 0)
    
    # average value across species unique to the baseline
    zbaruB <- sum(base[!shared])/suB
    
    if(is.na(zbaruB)) zbaruB <- 0
    
    # average value across species unique to the comparison
    zbaruD <- sum(comp[!shared])/suD    
    
    if(is.na(zbaruD)) zbaruD <- 0
    
    # average value across shared species in the baseline
    zbarcB <- mean(base[shared > 0])
    
    if(is.na(zbarcB)) zbarcB <- 0
    
    # average value across shared species in the comparison
    zbarcD <- mean(comp[shared > 0])
    
    if(is.na(zbarcD)) zbarcD <- 0
    
    # return data.frame
    ret <- data.frame(
      baseline = rownames(mat)[b],
      comparison = rownames(mat)[j],
      baselineS = sum(base > 0),
      comparisonS = sum(comp > 0),
      sharedS = sc,
      deltaFunc = sum(comp) - sum(base),
      RICH_L = -suB * zbarcB,
      COMP_L = -suB * (zbaruB - zbarcB),
      RICH_G = suD * zbarcD,
      COMP_G = suD * (zbaruD - zbarcD),
      CDE = sc * (zbarcD - zbarcB)
    )
    
  } ) )
  
  # get total diversity effect
  ret$TOTAL_DIV <- rowSums(ret[, 7:10])
  
  if(standardize == TRUE) ret[, 7:12] <- ret[, 7:12] / rowSums(abs(ret[, 7:11]), na.rm = TRUE)
  
  return(ret)
  
}

#' randomizeMat: a function to generate random permutations of a species-biomass matrix
#' 
#' @param amat a community (rows)-by-species (columns) abundance matrix
#' @param bmat a community (rows)-by-species (columns) biomass matrix
#' 
#' @return a matrix object with the same dimensions as mat
#' 
randomizeMat <- function(amat, bmat) {
  
  # Get per capita biomass matrix
  pcmat <- bmat / amat
  
  pcmat[is.na(pcmat)] <- 0
  
  # Get average per capita contribution by species
  species.means <- apply(pcmat, 2, function(x) mean(x[x > 0]))
  
  # Create random communities
  newmat <- do.call(rbind, lapply(1:nrow(amat), function(k) {
    
    # Get vector of new abundances
    newabund <- sapply(1:ncol(amat), function(l) sample(min(amat[, l]):max(amat[, l]), 1) )
    
    # Get vector of new biomasses
    newbiom <- newabund * species.means
    
    # Get number of species absent from the original community
    s <- sum(amat[k, ] == 0)
    
    # Randomly set s species' biomasses to zero
    newbiom[sample(1:length(newbiom), s)] <- 0
    
    return(newbiom)
    
  } ) )
  
}