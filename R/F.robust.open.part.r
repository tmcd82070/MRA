F.robust.open.part <- function(ch, p.star, s, g){
  # Purpose: 
  # Evaluate the "open" part of the robust design likelihood.
  #
  # ch = 3D capture indicator matrix
  # p.star = either vector of length (nprimary) or matrix of size  dim(ch)[1] X dim(ch)[2] = 
  #   (num animals) X (num primaries)  containing Pr(capture for individual i during primary occasion j)
  
  nan <- dim(ch)[1]
  nprimary <- dim(ch)[2]
  
  # collapse secondaries to just indicators for primaries.  no trap info, just 0-1. 
  ch.reduced <- F.collapse.secondaries(ch)
  
  #   print(ch.reduced[1:5,])
  
  # First and last primary occasions of encounter
  f <- apply(col(ch.reduced)*(ch.reduced >= 1), 1, function(x){ min(x[x>0])})
  l <- apply(col(ch.reduced)*(ch.reduced >= 1), 1, function(x){ max(x[x>0])})
  
  #   print(cbind(f,l)[1:5,])
  
  # SURVIVAL PART =========================================================
  # Expand survival so can multiply later
  s.mat <- matrix( rep(s, nan), nan, nprimary-1, byrow=T )
  
  # "Live" cells in ch matrix
  alive <- ch.reduced
  alive[ f <= col(alive) & col(alive) <= l ] <- 1
  
  # Survival part of likelihood
  s.ints <- alive[,-1]*alive[,-ncol(alive)]   # intervals known to be alive
  #  cat("S part ====== \n")
  #  print(s.ints)
  s.part <- s.ints * log(s.mat)
  #  print(s.part)
  s.part <- sum( s.part )
#      cat(paste("s part:", s.part, "\n"))
  
  # IMMIGRATION PART =========================================================
  # Gamma = immigration part of likelihood
  f.gamma.recurse<-function(h,g,p){
    
    #     cat("h=")
    #     cat(h)
    #     cat("\n")
    #     cat("g=")
    #     cat(g)
    #     cat("\n")
    #     cat("p=")
    #     cat(p)
    #     cat("\n")
    
    if( length(h)<=1 ) return(1)
    
    if(h[1]==1 & h[2]==1){
      # These would be gamma prime primes if Markovian movement
      ans <- (1-g[2])*(p[2])*f.gamma.recurse(h[-1],g[-1],p[-1])
      
    } else if(h[1]==1 & h[2]==0){
      # These would be gamma prime primes if Markovian movement
      ans <- (1-g[2])*(1-p[2])*f.gamma.recurse(h[-1],g[-1],p[-1]) +
        (g[2])*f.gamma.recurse(h[-1],g[-1],p[-1]) 
      
    } else if(h[1]==0 & h[2]==1){
      # These would be gamma primes if Markovian movement
      ans <- (1-g[2])*(p[2])*f.gamma.recurse(h[-1],g[-1],p[-1]) 
      
    } else if(h[1]==0 & h[2]==0){
      # These would be gamma primes if Markovian movement      
      ans <- (1-g[2])*(1-p[2])*f.gamma.recurse(h[-1],g[-1],p[-1]) +
        (g[2])*f.gamma.recurse(h[-1],g[-1],p[-1]) 
    }
    ans
  }
  
  f.gamma.indiv <- function(i,ch,g,p,f,l){
    pos <- seq(along=ch[i,])
    ind <- f[i] <= pos & pos <= l[i]
    ch.i <- ch[i,ind]
    g.i <- g[ind]   # note, g and p are one element shorter than chi.i
    p.i <- p[i,ind]   # this is because condition on first capture.
    # computations in recursive part are based on the second of 
    # capture indicator pairs.
    #     cat("-----\n")
    ans <- f.gamma.recurse(ch.i, g.i, p.i)
  }
  
  # Make sure p.star is a matrix.
  if( length(p.star) == nprimary){
    # p.star is vector, rep across individuals
    #cat("p.star was not a matrix.  Is this okay?")
    p.star <- matrix( p.star, nan, nprimary, byrow=T )
#    print(p.star)
  }
  
  g.part <- sapply(1:nan, f.gamma.indiv, ch=ch.reduced, g=g, p=p.star, f=f, l=l)
#     print(cbind(ch.reduced,g.part)[1:10,])
  
  g.part <- sum(log(g.part))
#      cat(paste("g part:", g.part, "\n"))
  
  # CHI PART =========================================================
  f.chi.recurse <- function(i,g,p,s,K){
    if( i==K ) return(1)
    
    #     cat(paste("i=",i,"\n"))
    #     cat("g=")
    #     cat(g)
    #     cat("\n")
    #     cat("p=")
    #     cat(p)
    #     cat("\n")
    #     cat("s=")
    #     cat(s)
    #     cat("\n")
    
    
    
    # From Kendal et al 1997:
    
    ans <- 1 - s[i]*(1 - (1 - (1 - g[i+1])*p[i+1])*f.chi.recurse(i+1,g,p,s,K))
    ans
  }
  
  f.chi.indiv <- function(i,g,p,s,l,K){
    ans <- f.chi.recurse(l[i], g, p[i,], s, K)
  }
  
  chi.part <- sapply(1:nan, f.chi.indiv, g=g, p=p.star, l=l, s=s, K=nprimary)
#     print(cbind(ch.reduced,chi.part)[1:10,])
  
  chi.part <- sum(log(chi.part))
#      cat(paste("Chi part:", chi.part, "\n"))
  
  # Return the "real" log likelihood, not the negative of it.  
  # Assumably it will be negated in another routine. 

  openLL = s.part + g.part + chi.part
  openLL
  
}