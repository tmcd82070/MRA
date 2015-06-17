F.spatial.pstar <- function(g0, sigma, traps, ch){
  #
  # Purpose: to estimate a p.star (or p. in spatial notation) for every captured individual every primary. 
  # do this by first estimating a location via Max Like. 
  
  nan <- dim(ch)[1]
  nprim <- dim(ch)[2]
  nsec <- dim(ch)[3]

  p.star <- matrix(NA, nan, nprim)
  
  for(j in 1:nprim){
    # for some reason, loop over primaries first.  Could do other way around. 
    for( i in 1:nan){
      if( sum(ch[i,j,]) > 0 ){
        h <- ch[i,j,]
        
        # Starting value for ML estimator
        x.cent <-  mean(traps[h,1]) # centroid x
        y.cent <-  mean(traps[h,2]) # centroid y
        
        # Find ML estimate of location
        X.ml <- optim( c(x.cent,y.cent), F.X.loglik, h=h, traps=traps, g0=g0[j], sigma=sigma[j])
        
        # Compute p. = p.star for ML location
        p. <- pdot(X.ml$par, traps, 
                   detectpar=list(g0=g0[j], sigma=sigma[j]), noccasions=nsec)
        
        # Store p. for output
        p.star[i,j] <- p.
      
      }
    }
    # Assign average for animal when it is missed. 
    p.star[i,is.na(p.star[i,])] <- mean(p.star[i,!is.na(p.star[i,])])
  }


  p.star
}
