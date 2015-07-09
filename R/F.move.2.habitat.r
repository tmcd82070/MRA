F.move.2.habitat <- function(x, hab, searchdist, h, tr, g0, sigma){
  # Move x to the nearest point in hab.  Check neighborhood to make sure 
  # log like is still highest there. 
  #   hab = a SpatialPoints object
  
  x <- SpatialPoints(x)
  x.buffer <- gBuffer( x, width=searchdist )
  
  inbuff <- over( hab, x.buffer )
  inbuff <- inbuff[!is.na(inbuff)]
  
  if(length(inbuff)>0){  
    funval <- rep(NA,length(inbuff))
    for( i in 1:length(inbuff)){
      funval[i] <- F.X.loglik(hab[inbuff[i]], h, traps, g0, sigma)
    }
    x <- hab[inbuff[which.min(funval)]]
    
  } 
  
  coordinates(x)
  
}