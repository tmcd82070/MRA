F.spatial.pstar2 <- function(g0, sigma, traps, ch, buffer){
  #
  # Purpose: to estimate a p.star (or p. in spatial notation) for every captured individual every primary. 
  # do this by first estimating a location via Max Like. 
  
  nan <- dim(ch)[1]
  nprim <- dim(ch)[2]
  nsec <- dim(ch)[3]
  
  # for pdot to work below, traps has to be a "traps" object.
  tr <- as.data.frame(traps)
  rownames(tr)<-1:nrow(tr)
  attr(tr,"detector") <- "multi"  
  class(tr) <- c("traps","data.frame")
  
  p.star <- matrix(NA, nan, nprim)
  
  # First, loop over occasions when individual was caught. 
  for(j in 1:nprim){
    # for some reason, loop over primaries first.  Could do other way around. 
    
    for( i in 1:nan){
      if( sum(ch[i,j,]) > 0 ){
        h <- ch[i,j,]
        
        # Starting value for ML estimator
        x.cent <-  mean(tr[h,1]) # centroid x
        y.cent <-  mean(tr[h,2]) # centroid y
        
        # Find ML estimate of location
        X.ml <- optim( c(x.cent,y.cent), F.X.loglik, h=h, traps=tr, g0=g0[j], sigma=sigma[j])
        
        # Compute p. = p.star for ML location
        p. <- pdot(X.ml$par, tr, 
                   detectpar=list(g0=g0[j], sigma=sigma[j]), noccasions=nsec)
        
        # Store p. for output
        p.star[i,j] <- p.
        
      }
      
    }
  }
  
  # Now loop over every location in the mask, and average
  habmask <- SpatialPoints(traps,proj4string=CRS("+init=epsg:26915"))
  coordnames(habmask) <- c("x","y")
  sa <- gBuffer( habmask, width=buffer )
  sa <- spsample( sa, 1000, "regular" )
  sa <- coordinates(sa)
  P <- rep(0,nprim)
  for(j in 1:nprim){
    for( i in 1:nrow(sa)){
      p. <- pdot(sa[i,], tr, 
                 detectpar=list(g0=g0[j], sigma=sigma[j]), noccasions=nsec)
      P[j] <- P[j] + p.
    }
  }
  P <- P / nrow(sa)
  
  
  p.star <- sapply(1:nprim, function(j,ps,P){
    ps.j <- ps[,j]
    ps.j[ is.na(ps.j) ] <- P[j]
    ps.j
    
  }, ps=p.star, P=P)
  
  #   # Assign average for animal when it is missed. 
  #   p.star <- apply(p.star, 1, function(x){
  #     x[is.na(x)] <- mean(x,na.rm=TRUE)  
  #     x
  #   })
  #   p.star <- t(p.star)
  
  # Assign average of animals that were captured. this is not 
  # Average over entire habitat mask, but try it. 
#   p.star[is.na(p.star)] <- mean(p.star,na.rm=TRUE)
  
  p.star
}
