F.spatial.pstar2 <- function(g0, sigma, traps, ch, sa, pt.searchdist=NULL){
  #
  # Purpose: to estimate a p.star (or p. in spatial notation) for every captured individual every primary. 
  # do this by first estimating a location via Max Like.
  #
  # g0, sigma = parameters of distance function
  # traps = R X 2 matrix of trap coordinates
  # ch = nan X nprim X nsec matrix of trap-capture histories
  # sa = S X 2 matrix of possible locations in study area. Used for integration.  
  
  nan <- dim(ch)[1]
  nprim <- dim(ch)[2]
  nsec <- dim(ch)[3]
  
  # for pdot to work below, traps has to be a "traps" object.
  tr <- as.data.frame(traps)
  rownames(tr)<-1:nrow(tr)
  attr(tr,"detector") <- "multi"  
  class(tr) <- c("traps","data.frame")
  
  p.star <- matrix(NA, nan, nprim)
  
  # If study area is not a SpatialPOints object, convert it to one.  Later, we convert back to matrix
  if( length(grep("SpatialPoints", class(sa))) == 0 ){
    sa <- SpatialPoints(sa)
  }
  
  # Set distance to search for habitat points from ML estimat of point.  For gridded habitat, 
  # this should be sqrt(2)*grid spacing
  if(is.null(pt.searchdist)){
    pt.searchdist <- F.pt.searchdist(sa)
  }
  
  # First, loop over occasions when individual was caught. 
  for(j in 1:nprim){
    # for some reason, loop over primaries first.  Could do other way around. 
    
    for( i in 1:nan){
      if( sum(ch[i,j,]) > 0 ){
        h <- ch[i,j,]
        
        # Starting value for ML estimator
        x.cent <-  mean(tr[h,1]) # centroid x
        y.cent <-  mean(tr[h,2]) # centroid y
        
        if( any(is.na(c(x.cent,y.cent))) ){
          cat("Missing initial values:\n")
          print(c(x.cent,y.cent))
          x.cent <- mean(tr[,1])
          y.cent <- mean(tr[,2])
        }
        # Find ML estimate of location
        X.ml <- optim( c(x.cent,y.cent), F.X.loglik, h=h, traps=tr, g0=g0[j], sigma=sigma[j])
        
        # X.ml may not be in habitat. Check and move it to habitat point with highest loglik. 
        X.ml <- F.move.2.habitat(X.ml, sa, pt.searchdist, h, tr, g0[j], sigma[j])
        
        # Compute p. = p.star for location in habitat at least close to ML
        p. <- pdot(X.ml, tr, 
                   detectpar=list(g0=g0[j], sigma=sigma[j]), noccasions=nsec)
        
        # Store p. for output
        p.star[i,j] <- p.
        
      }
      
    }
  }
  
 
  
  # Now loop over every location in the mask, and average
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

# tmp.1 <- c(0.4049685, 0.4187768, 0.3477709, 0.3831983)
# tmp.1 <- 1/(1+exp(-tmp.1))
# tmp.2 <- c(2.541767, 2.548663, 2.54559, 2.547148)
# tmp.2 <- exp(tmp.2)
# 
# tmp <- F.spatial.pstar2( tmp.1, tmp.2, coordinates(X), H, 100)
#             
