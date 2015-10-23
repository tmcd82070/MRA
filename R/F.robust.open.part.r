#' @title Open part of robust design likelihood
#' 
#' @description Computes open portion of the robust design likelihood 
#' 
#' @param ch A 3-D array of capture histories.  Size is nan X nprim X nsecondaries.  Missing secondaries 
#' are specified with columns of all NA's (i.e., if the 4th secondary of the ith primary was not done, 
#' \code{all(is.na(ch[,i,4])) == TRUE}).  Cells in the array contain trap number that caught the individual.
#' Trap numbers are rows in a matrix embedded in the \code{traps} object.  \code{ch[i,j,k]} = 0 if 
#' animal i was uncaptured during 
#' secondary occasion k of primary occasion j. \code{ch[i,j,k]} = x (where x integer > 0) means 
#' animal i was captured 
#' during secondary occasion k of primary occasion j in trap x which has coordinates
#' \code{traps[[j]][[k]][x,]}.
#' 
#' @param p.star A nan X nprimaries matrix containing estimated probabilities of detection during each 
#'    primary.  \code{p.star[i,j]} = probability of capturing individual i at least once during primary 
#'    occasion j (over all secondaries contained in j).  
#'    This is the "p star" of Kendall's papers, which is the "p dot" of Borchers and Efford.
#'    
#' @param s Survival probabilities.  A vector of length (nprimary-1) corresponding to intervals between primary occasions. 
#'  logit link is assumed.
#'  
#' @param g Immigration and emigration parameters. A vector of length  (nprimary - 1) 
#' elements corresponding to intervals between primaries.  Only the random emigration model 
#' is implemented here.  In Kendall's (1997) notation, gamma prime = gamma prime prime. 
#'  logit link is assumed.
#'  
#'  @return A scalar.  The log likelihood for the open part of a robust design model.  Usually,
#'  one will want to multiply this number by -1 to minimize. 
#' 
F.robust.open.part <- function(ch, p.star, s, g){

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
  
  cat(paste("s part:", s.part, "\n"))
  
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
  
  print(p.star)
  
  # Make sure p.star is a matrix.
  if( length(p.star) == nprimary){
    # p.star is vector, rep across individuals
    #cat("p.star was not a matrix.  Is this okay?")
    p.star <- matrix( p.star, nan, nprimary, byrow=T )
#    print(p.star)
  }

  cat("in F.robust.open.part:")
  print(p.star)
  
assign("tmp.pstar2",p.star, pos=.GlobalEnv)

  g.part <- sapply(1:nan, f.gamma.indiv, ch=ch.reduced, g=g, p=p.star, f=f, l=l)
#     print(cbind(ch.reduced,g.part)[1:10,])
  
  g.part <- sum(log(g.part))

  cat(paste("g part:", g.part, "\n"))
  
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
     cat(paste("i=",i,"l[i]=",l[i],"K=",K,"ans=",ans,"\n"))   
    ans
  }
  
  chi.part <- sapply(1:nan, f.chi.indiv, g=g, p=p.star, l=l, s=s, K=nprimary)

#   print(cbind(ch.reduced,chi.part)[1:10,])
#   print(class(chi.part))

  chi.part <- sum(log(chi.part))
  
  cat(paste("Chi part:", chi.part, "\n"))
  
  # Return the "real" log likelihood, not the negative of it.  
  # Assumably it will be negated in another routine. 

  openLL = s.part + g.part + chi.part
  openLL
  
}