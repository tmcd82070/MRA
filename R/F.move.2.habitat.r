F.move.2.habitat <- function(x, hab, h, tr, g0, sigma, num.nneighbors=10){
  # 
# x = vector of length 2 = (x,y) coordinate
# hab = SpatialPixels or SpatialGrid 
# h = nsecondary X 1 vector of capture history
# tr = list length nsecondary of "on" trap locations
# g0 = vector length 1 of g0 parameter (no variation across secondaries currently allowed. 
# sigma = vector length 1 of sigma parameter
# num.nneighbors = number of nearest neighbors to inspect when x is outside habitat.

  #Move x to the nearest point in hab.  Check neighborhood to make sure 
  # log like is still highest there. 

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