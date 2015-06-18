# Fits a Spatial Robust design model via maximum likelihood
# 
# Objective function -----
F.spat.robust.loglik <- function( beta, ch, traps, buffer ){
  #
  # inputs:
  #   beta = coefficients
  #   ch = 3-d array, rows are individuals, columns are secondary occasions, pages are primary occasions. 
  #        columns that are all NA are "not there", meaning variable number of secondary occasions are specified 
  #        by missing columns in one of the pages (e.g., ch[,3,2] could be all NA, but ch[,3,1] is there.  
  #        this means 3 secondary occasions during primary 1 but only 2 during primary 2). 
  #        Value in each cell is the number (row in trap location matrix) of trap that caught the individual.
  #        "0" means un-captured.
  #   traps = K x 2 matrix of trap locations
  #   buffer = distance to buffer trap locations
  #
  #
  
  if( !inherits(ch,"array") ) stop("ch must be an array")
  if( length(dim(ch)) != 3) stop(paste("ch must have 3 dimensions.", length(dim(ch)), "found."))
  
  # Find dimensions                                 
  d <- dim(ch)
  nan <- d[1]
  nprimary <- d[2]
  nsecondary <- apply(ch,2,function(x){
    xx<-apply(x,2,function(y){all(is.na(y))})
    sum(!xx)
  })
  
  # Evaluate the model to get real parameters
  f.real.model<-function(beta, nprim){
    # Normally, you would have a X matrix here.
    # For now, use names of beta to evaluate model.
    #
    # Current implementation allows small t or dot model for  s
    # This returns linear predictors.  Link functions are applied later.
    
    s.parms <- grep("^s",names(beta))
    g.parms <- grep("^gamma",names(beta))

    # For now, same SECR parameters in each primary session
    g0.parms <- grep("^g0",names(beta))
    sigma.parms <- grep("^sigma",names(beta))
    
    # s parameters
    s <- beta[s.parms]
    if( length(s) != (nprim-1) ){
      # Just use first 
      s <- rep(s[1],nprim-1)
    }
    
    # g parameters
    # This is the random emigration model, constant across time
    g <- beta[g.parms]
    if( length(g) == (nprim-2)){
      # constrain last g to second to last g.
      g <- c(NA,g,g[nprim-2])
    } else if(length(g) == (nprim-1)){
      g <- c(NA,g)  # use all gammas
    } else if( length(g) != (nprim-1)){
      g <- c(NA,rep(g[1],(nprim-1))) # first gamma not possible.
    } 
    # For rest of routine to work g must have length nprim.  g[1]=NA
    
    # g0 parameters
    g0 <- beta[g0.parms]
    if( length(g0) != nprim ){
      # Just use first 
      g0 <- rep(g0[1],nprim)
    }    
    
    # sigma parameters
    sigma <- beta[sigma.parms]
    if( length(sigma) != nprim ){
      # Just use first 
      sigma <- rep(sigma[1],nprim)
    } 
    
    list(s.eta=s, g.eta=g, g0.eta=g0, sigma.eta=sigma)
  }
  parms <- f.real.model(beta,nprimary)
  
  print(parms)
  print(nsecondary)
  # Compute SECR likelihood for each occasion ========================
  # For now, SECR parameters are constant accross (primary) sessions 
  closedLL <- sapply(1:nprimary,function(i,c.hist,b,ns,trps,buff){
    ch1 <- c.hist[,i,1:ns[i]]  # remove NA's here
    ch1 <- ch1[rowSums(ch1>0)>0,]    # remove all 0 lines here
    
    F.spat.loglik2(b,ch1,trps,buff)
  },
  c.hist=ch, 
  b=c(parms$g0.eta, parms$sigma.eta),
  ns=nsecondary, 
  trps=traps,
  buff=buffer
  )
  closedLL <- sum(closedLL)
  
  cat(paste("SECR part:", closedLL, "\n"))
  
  
  
  # Compute Open part of robust design likelihood =========================
  # Take links
  s <- 1/(1+exp(-parms$s.eta))
  gamma <- 1/(1+exp(-parms$gamma.eta))
  g0 <- 1/(1+exp(-parms$g0.eta))
  sigma <- exp(parms$sigma.eta)   # Note log link here, rather than logit
  
  p.star <- F.spatial.pstar(g0, sigma, traps, ch)   # returns nan X nprimary matrix or p.dots
  
  openLL <- F.robust.open.part(ch,p.star,s,gamma)
  
  
  # Done ======================================================
  ll <- openLL - closedLL

  print(ll)

  ll
}

