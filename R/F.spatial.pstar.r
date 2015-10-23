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
#'  @param hab.mask A SpatialPixels or SpatialGrid object (or their *DataFrame analogs) containing allowable locations for activity centers. 
#'    Exactly like masks in the SECR package, \code{hab.mask} defines the outer limit of integration, defines 
#'    sites that are habitat and thus can be occupied, descritizes habitat covariates for use in models.
#'    The current implementation does not allow the habitat mask to change between primaries or secondaries. 
#'    The mask is constant throughout the study.  \code{hab.mask} is needed here because 
#'    use it to compute p.star when an animal was not captured during a primary. 
#'         
#'    
#' @return A nan X nprimaries matrix containing estimated probabilities of detection during each 
#'    primary.  This is the "p star" of Kendall's papers, which is the "p dot" of Borchers and Efford. 
#'       
F.spatial.pstar <- function(g0, sigma, traps, aclocs, hab.mask){


  # Loop over primaries -------------
  nan <- nrow(aclocs)
  nprim <- length(traps)
  
  pstar <- matrix(NA, nan, nprim)
  for( j in 1:nprim ){
    
    # Each loop, replicate the p. calculations in F.spat.loglik.X
    p.star <- F.spat.capProbs(c(g0[j], sigma[j]), traps[[j]], aclocs[,j,] )
    
     pstar[,j] <- p.star$pdot
     
     # Handle the animals that were not seen -------------------
     uncaught.ind <- is.na(aclocs[,j,1])
     if( any( uncaught.ind ) ){
       #  At least one animal was not captured this primary
       hab.pts <- coordinates(hab.mask)
       p.star <- F.spat.capProbs(c(g0[j], sigma[j]), traps[[j]], hab.pts )
       pstar[uncaught.ind,j] <- mean(p.star$pdot)  # or sum?  don't think it matters, but maybe
     }
     
#      cat("in F.spatial.pstar-----")
#      print(pstar[,j])
#      cat("-------\n")
     
     
  }
  

  pstar
}
    
