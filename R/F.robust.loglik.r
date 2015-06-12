# Fit Robust design model via maximum likelihood
# 
# Objective function -----
F.robust.loglik <- function( beta, ch ){
  #
  # inputs:
  #   ch = 3-d array, rows are individuals, columns are secondary occasions, pages are primary occasions. 
  #        columns that are all NA are "not there", meaning variable number of secondary occasions are specified 
  #        by missing columns in one of the pages (e.g., ch[,3,2] could be all NA, but ch[,3,1] is there.  
  #        this means 3 secondary occasions during primary 1 but only 2 during primary 2). 
  #
  #
  
  if( !inherits(ch,"array") ) stop("ch must be an array")
  if( length(dim(ch)) != 3) stop(paste("ch must have 3 dimensions.", length(dim(ch)), "found."))
             
  # Find dimensions                                 
  d <- dim(ch)
  nan <- d[1]
  nprimary <- d[3]
  nsecondary <- apply(ch,3,function(x){
      xx<-apply(x,2,function(y){all(is.na(y))})
      sum(!xx)
    })
  
  # Evaluate the model to get real parameters
  f.real.model<-function(beta, nprim){
    # Normally, you would have a X matrix here.
    # For now, use names of beta to evaluate model.
    #
    # Current implementation allows small t or dot model for p and s
    # This returns linear predictors.  Link functions are applied later.
    
    p.parms <- grep("^p",names(beta))
    s.parms <- grep("^s",names(beta))
    g.parms <- grep("^g",names(beta))
    
    # p parameters
    p <- beta[p.parms]
    if( length(p) != nprim ){
      # Just use first p
      p <- rep(p[1],nprim)
    }
    
    # s parameters
    s <- beta[s.parms]
    if( length(s) != (nprim-1) ){
      # Just use first 
      s <- rep(s[1],nprim-1)
    }
    
    # g parameters
    # This is the random emigration model, constant across time
    g <- beta[g.parms][1]
    if( length(g) != (nprim-1)){
      g <- c(NA,rep(g[1],(nprim-1))) # first gamma not possible.
    }
    
    list(p.eta=p, s.eta=s, g.eta=g)
  }
  parms <- f.real.model(beta,nprimary)
  
  print(parms)
  
  # For now, implement Huggins model inside each primary occasion.  No heterogeneity. 
  closedLL <- sapply(1:nprimary,function(i,c.hist,b,ns){
        F.hug.loglik(b[i],c.hist[,1:ns[i],i])
      },
      c.hist=ch, 
      b=parms$p.eta,
      ns=nsecondary
  )
  closedLL <- sum(closedLL)
  
  cat(paste("Closed part:", closedLL, "\n"))
  
  # Take links
  p <- 1/(1+exp(-parms$p.eta))
  s <- 1/(1+exp(-parms$s.eta))
  g <- 1/(1+exp(-parms$g.eta))
  
  # Evaluate robust design part of likelihood
  p.star <- 1 - (1-p)^nsecondary  # constant within secondary.  This is prob of catch each primary
  ch.reduced <- apply(ch,2,function(c.hist){
    as.numeric(rowSums(c.hist>0,na.rm=T)>0)
  })  # collapse secondary periods to captured or not.  This is nan X nprimary
  
#   print(ch.reduced[1:5,])
  
  # First and last primary occasions of encounter
  f <- apply(col(ch.reduced)*(ch.reduced >= 1), 1, function(x){ min(x[x>0])})
  l <- apply(col(ch.reduced)*(ch.reduced >= 1), 1, function(x){ max(x[x>0])})
  
#   print(cbind(f,l)[1:5,])
  
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
  s.part <- -sum( s.part )
  cat(paste("s part:", s.part, "\n"))
  
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
    p.i <- p[ind]   # this is because condition on first capture.
                    # computations in recursive part are based on the second of 
                    # capture indicator pairs.
#     cat("-----\n")
    ans <- f.gamma.recurse(ch.i, g.i, p.i)
  }
  
  g.part <- sapply(1:nan, f.gamma.indiv, ch=ch.reduced, g=g, p=p.star, f=f, l=l)
#   print(cbind(ch.reduced,g.part)[1:10,])
  
  g.part <- -sum(log(g.part))
  cat(paste("g part:", g.part, "\n"))
  
  # Chi part of likelihood
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
    ans <- f.chi.recurse(l[i], g, p, s, K)
  }
  
  chi.part <- sapply(1:nan, f.chi.indiv, g=g, p=p.star, l=l, s=s, K=nprimary)
#   print(cbind(ch.reduced,chi.part)[1:10,])

  chi.part <- -sum(log(chi.part))
  cat(paste("Chi part:", chi.part, "\n"))

  ll <- s.part + g.part + chi.part + closedLL
  ll
}

# =================================================================

source("./r/F.hug.loglik.r")

fn <- "c:/Program Files (x86)/MARK/Examples/Robust Design Huggins.inp"
ch <- read.table( fn, skip=3, colClasses=c("character","numeric",NULL), col.names=c("h","freq","sc"))
chh <- NULL
for(i in 1:nrow(ch)){
  chh <- c(chh, rep(ch$h[i], ch$freq[i]))
}

ch.expand <-function(ch){
  nc <- nchar(ch[1])
  matrix(as.numeric(t(sapply(ch,substring,first=1:nc,last=1:nc))),length(ch))
}
CH <- array(NA,c(length(chh),5,5))
CH[,1:2,1] <- ch.expand(substring(chh,1,2))
CH[,1:2,2] <- ch.expand(substring(chh,3,4))
CH[,1:4,3] <- ch.expand(substring(chh,5,8))
CH[,1:5,4] <- ch.expand(substring(chh,9,13))
CH[,1:2,5] <- ch.expand(substring(chh,14,15))

beta <- c(0,0,0)
names(beta) <- c("p","s","g")
tmp <- F.robust.loglik( beta, CH)
print(tmp)


beta <- c(1.7140078, -1.9046820, 0.6841988, 0.4789440, 0.4296271, 0.4280551, -0.2745019)
names(beta) <- c("s","g",rep("p",5))
tmp <- F.robust.loglik( beta, CH)
print(tmp)


fit1 <- optim( beta, F.robust.loglik, ch=CH,  hessian=T, method="Nelder-Mead", control=list(reltol=1e-10, maxit=500))



# ======== some checking =======

# This gives the correct huggins estimates for separate primaries. 
# ch1 <- CH[,1:5,4]
# ch1 <- ch1[rowSums(ch1)>=1,]
# fit1 <- optim( c(0), F.hug.loglik, ch=ch1, hessian=T, method="Brent", lower=-4, upper=4,  control=list(reltol=1e-10, maxit=500))

# Perspective plot =======

# library(rgl)
# pp <- seq(0.01, .99, length=50)
# oo <- seq(0.01, .99, length=50)
# pp <- log(pp/(1-pp))
# oo <- log(oo/(1-oo))
# df <- expand.grid(oo=oo, pp=pp)
# zz <- apply(df, 1, F.cjs.loglik, ch=ch)
# zz <- matrix(zz, length(oo), length(pp))
# ncols <- 20
# col.levs <- cut( zz, ncols, labels=FALSE)
# mycols <- heat.colors(ncols)[col.levs]
# persp3d(oo,pp,zz, col=mycols)
