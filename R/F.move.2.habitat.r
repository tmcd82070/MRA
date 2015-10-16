#'  @title Check that a point is in habitat, and if not, move it to habitat
#'  
#'  @description Check that a point is containing in a habitat mask.  If so, snap it to 
#'    the center of the containing cell.  If not, move x to the center of the cell with 
#'    the highest likelihood. 
#'    
#'  @param x Vector of length 2 equal to the (x,y) coordinates of the point to test.
#'  
#'  @param hab Either a SpatialPixels or SpatialGrid object describing the habitat mask. 
#'  
#'  @param h Vector of length nsecondary occasions containing trap numbers that caught 
#'    a particular individual.  This is used in evaluation of the likelihood if x has to be moved. 
#'    
#'  @param tr A list of length nsecondary occasions.  Each element in tr is a T X 2 matrix of 
#'    trap locations that were "on" during the occasion. Used in evaluation of likelihhod if x has 
#'    to be moved. 
#'    
#'  @param g0 the sightability function height (or intercept) parameter.   Used in evaluation of likelihhod if x has 
#'    to be moved.
#'    
#'  @param sigma the sightability function width parameter. Used in evaluation of likelihhod if x has 
#'    to be moved.
#'    
#'  @param num.nneighbors Number of nearest neighbors to evaluate as potential places for x.  See Details.
#'  
#'  @details If \code{x} is not in habitat, the likelihood for containing the activity 
#'    center (i.e., \code{F.X.loglik}) is evaluated at each of the \code{num.nneighbors} pixels 
#'    in \code{hab} are nearest to \code{x}.  \code{x} is then moved to the center of the cell with 
#'    highest likelihood.  Note that the entire habitat grid is not searched for a maximum.  Only 
#'    those cells in the vicinity of the imput x, where vicinity is defined \code{num.nneighbors} 
#'    nearest neighbors. 
#'    
#'  @return A vector of length 2 containing the (x,y) coordinates of center of the cell containing x, 
#'    or of the cell with highest likelihood of containing an activity center. 
#'    
F.move.2.habitat <- function(x, hab, h, tr, g0, sigma, num.nneighbors=10){
 
  x <- SpatialPoints(data.frame(x=x[1], y=x[2]), proj4string = CRS(proj4string(hab)))
  x.in.habitat <- over(x,hab)  # because x first here, result is length 1 = index of pixel containing x
  
  if( !is.na(x.in.habitat) ){
    # We are good.  X is in habitat. Snap x to center of containing cell. 
    x <- SpatialPoints(hab[x.in.habitat,])
    
  } else {
    # We are not good.  X is outside habitat. 
    hab.pts <- as(hab, "SpatialPoints")  # centers of habitat pixels (or grid cells)
    hab.dists <- gDistance(x, hab.pts, byid=TRUE)  # Distance from point to all habit locations
    K <- min(num.nneighbors, length(hab.dists))
    d <- sort(hab.dists)[K]  # a distance that should catch num.nneighbors points
    
    # Find all habitat points within d distance of point
    x <- gBuffer( x, width = d)
    inbuff <- over(hab.pts, x)
    inbuff <- inbuff[!is.na(inbuff)]
    inbuff <- names(inbuff)
    
    funval <- rep(NA,length(inbuff))
    for( i in 1:length(inbuff) ){
      funval[i] <- F.X.loglik(coordinates(hab.pts[inbuff[i],]), h, traps, g0, sigma)
    }
    x <- hab.pts[inbuff[which.min(funval)],]
  }

  coordinates(x)

} 