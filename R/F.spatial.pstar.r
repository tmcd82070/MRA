#' @title Compute individual's probability of capture during primaries 
#' 
#' @description Given estimates of detection parameters, trap locations, and estimated 
#'    activity center locations, this routine computes 
#'    the probability of each individual being captured at least once during the secondaries of each 
#'    primary occasion.  
#'    
#' @param g0 Vector of distance function height parameters (untransformed, in probability space).  Length is 
#'    number of primaries. 
#'    
#' @param sigma Vector of distance function width parameters (untransformed, in probability space). 
#'    Length is number of primaries.
#'    
#' @param traps Either a \code{SpatialPoints} object or a list of lists of trap locations.  
#'    See \code{F.spat.robust.loglik.X} for description.
#' 
#' @param aclocs A nan X nprimaries X 2 array of estimated activity center locations for every animal 
#'    every primary occasion. 
#'    
#' @return A nan X nprimaries matrix containing estimated probabilities of detection during each 
#'    primary.  This is the "p star" of Kendall's papers, which is the "p dot" of Borchers and Efford. 
#'       
F.spatial.pstar <- function(g0, sigma, traps, aclocs){


  # Loop over primaries -------------
  nan <- nrow(aclocs)
  nprim <- length(traps)
  
  pstar <- matrix(NA, nan, nprim)
  for( j in 1:nprim ){
    
    # Each loop, replicate the p. calculations in F.spat.loglik.X

    pstar[,j] <- F.spat.capProbs(c(g0=g0[j], sigma=sigma[j]), traps[[j]], aclocs[,j,] )
  }
  
  pstar
}
    
