# Fit Robust design model via maximum likelihood
# 
# Objective function -----
F.robust.loglik2 <- function( beta, ch ){
  #
  # inputs:
  #   ch = 3-d array, rows are individuals, columns are primary occasions, pages are secondary occasions within
  #        primary occasions. 
  #        pages that are all NA are "not there", meaning variable number of secondary occasions are specified 
  #        by missing pages in one of the columns (e.g., ch[,2,3] could be all NA, but ch[,1,3] is there.  
  #        this means 3 secondary occasions during primary 1 but only 2 during primary 2). 
  #
  # Note: this routine differs from F.robust.loglik in that this routine uses the m.array to compute the 
  # open part of the model.  F.robust.loglik uses recursive routines, and I could not get them to work. 
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
  
#   print(nan);print(nprimary);print(nsecondary)
  
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
    
    
    list(p.eta=p, s.eta=s, g.eta=g)
  }
  parms <- f.real.model(beta,nprimary)
  
#   print(parms)
  
  # Evaluate closed part of model ===================================================================
  # For now, implement Huggins model inside each primary occasion.  No heterogeneity. 
  closedLL <- sapply(1:nprimary,function(i,c.hist,b,ns){
#    cat("---------\n");print(i);print(b[i]);print(ns[i])
    ch1 <- c.hist[,i,1:ns[i]]  # remove NA's here
    ch1 <- ch1[rowSums(ch1>0)>0,]    # remove all 0 lines here
    ans<-F.hug.loglik(b[i],ch1)
#     print(ans)
    ans
  },
  c.hist=ch, 
  b=parms$p.eta,
  ns=nsecondary
  )
# print(closedLL)
  closedLL <- -sum(closedLL)  # F.hug.loglik returns -LL. re-negative it here so can add later.
  
  cat(paste("Closed part:", closedLL, "\n"))
  
  
  # Evaluate open part of model =====================================================================
  # Take links
  p <- 1/(1+exp(-parms$p.eta))
  s <- 1/(1+exp(-parms$s.eta))
  g <- 1/(1+exp(-parms$g.eta))
  
  # Evaluate robust design part of likelihood
  p.star <- 1 - (1-p)^nsecondary  # constant within secondary.  This is prob of catch each primary
  
  # Collapse secondaries
# print(dim(ch))
  ch.reduced <- F.collapse.secondaries(ch)
#   print(ch.reduced[1:50,])
  
  # Compute m-array
  m.array <- F.m.array(ch.reduced)
  m <- m.array$m
  R <- m.array$R

#   print(m)

  # Follow equation (2) in Kendall et al 1997 to evaluate the likelihood. 
  ns <- nprimary  # to save typing
  
#   print(p.star)
#   print(g)
  
  pp <- cbind(NA,matrix((1-g[2:ns])*p.star[2:ns],ns-1,ns-1,byrow=T))
  ss <- cbind(NA,matrix(s,ns-1,ns-1,byrow=T))

  ss[row(ss)>=col(ss)]<-NA
  pp[row(pp)>=col(pp)]<-NA

  pp.conj <- cbind(NA,1-pp[,1:(ns-1)])
  
  f.cumprod <-function(x){
    xx<-cumprod(x[!is.na(x)])
    x[!is.na(x)]<-xx
    x
  }
  ss.cprod <- t(apply(ss,1,f.cumprod))
  pp.conj.cprod <- t(apply(pp.conj,1,f.cumprod))
  pp.conj.cprod[col(pp.conj.cprod)-row(pp.conj.cprod) == 1] <- 1
  
#   print( log(ss.cprod))
#   print( log(pp.conj.cprod))
#   print( log(pp))
  
  cellp <- m * (log(ss.cprod) + log(pp.conj.cprod) + log(pp))
  
#  print(cellp)
  
  ll <- sum(cellp,na.rm=T)
  
  
  # Chi part of likelihood ===============================================================
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
  
# print(g)
# print(p.star)
# print(s)

  chi <- sapply(1:(ns-1), f.chi.recurse, g=g, p=p.star, s=s, K=ns)
#   print(cbind(1:(ns-1),chi))
  
  never.seen.again <- R[1:(ns-1)] - rowSums(m,na.rm=T)
#  print(never.seen.again)
  
  chi.part <- sum(never.seen.again*log(chi))

 cat(paste("ll part:", ll, "\n"))
 cat(paste("Chi part:", chi.part, "\n"))
 cat(paste("closed part:", closedLL, "\n"))

  ll <- ll + chi.part + closedLL
  print(-ll)
  -ll
}

# =================================================================

# source("./../r/F.hug.loglik.r")
# source("./../r/F.collapse.secondaries.r")
# source("./../r/F.m.array.r")
# 
# # An example from MARK ===================================================================================
# fn <- "c:/Program Files (x86)/MARK/Examples/Robust Design Huggins.inp"
# ch <- read.table( fn, skip=3, colClasses=c("character","numeric",NULL), col.names=c("h","freq","sc"))
# chh <- NULL
# for(i in 1:nrow(ch)){
#   chh <- c(chh, rep(ch$h[i], ch$freq[i]))
# }
# 
# ch.expand <-function(ch){
#   nc <- nchar(ch[1])
#   matrix(as.numeric(t(sapply(ch,substring,first=1:nc,last=1:nc))),length(ch))
# }
# CH <- array(NA,c(length(chh),5,5))
# CH[,1,1:2] <- ch.expand(substring(chh,1,2))
# CH[,2,1:2] <- ch.expand(substring(chh,3,4))
# CH[,3,1:4] <- ch.expand(substring(chh,5,8))
# CH[,4,1:5] <- ch.expand(substring(chh,9,13))
# CH[,5,1:2] <- ch.expand(substring(chh,14,15))
# 
# cat("M(t+1):\n")
# print(apply(CH, 2, function(x){sum(rowSums(x,na.rm=T)>0)}))
# 
# # beta <- c(0,0,0)
# # names(beta) <- c("p","s","g")
# # tmp <- F.robust.loglik( beta, CH)
# # print(tmp)
# 
# 
# # THis is the most general model S(t); g(t); p(t)
# beta <-c(
# 2.0475783,
# 1.9063523,
# 1.2076639,
# 1.3484628,
# -2.3646612,
# -2.0140449,
# -1.7334309,
# 0.6841976,
# 0.3841437,
# 0.4176103,
# 0.4364079,
# -0.0943085)
# #beta <- rep(0,4+3+5)
# names(beta)<-c(rep("s",4),rep("g",3),rep("p",5))
# 
# #beta <- fit1$par
# 
#  tmp <- F.robust.loglik2( beta, CH)
# print(tmp)


#fit1 <- optim( beta, F.robust.loglik2, ch=CH,  hessian=T, method="BFGS", control=list(reltol=1e-10, maxit=1000))

## THIS WORKS!!! I compared several models with same in MARK and everything matches.  
##  fit1$par are the coefficient estimates
##  fit1$convergence should equal 0
##  Standard errors are sqrt(diag(solve(fit1$hessian))), and match MARK's perfectly.
##  2*fit1$value equals MARK's -2loglik at solution.


# This is S(.),g(.),p(t)

# beta <- c(1.7, -1.9, 0.68, 0.47, 0.43, 0.43, -0.27)
# names(beta) <- c("s","g",rep("p",5))
# 
# fit3 <- optim( beta, F.robust.loglik2, ch=CH,  hessian=T, method="BFGS", control=list(reltol=1e-10, maxit=1000))

## THIS CONSTRAINED MODEL ALSO WORKS. Same answers as MARK.



# # An example with perfect data from Bill Kendall ===========================================================
# fn <- "genmssv4_trent.inp"
# ch <- read.table( fn,  colClasses=c("character","character"), col.names=c("h","freq"))
# ch$freq <- as.numeric(gsub(";","",ch$freq))
# chh <- NULL
# for(i in 1:nrow(ch)){
#   chh <- c(chh, rep(ch$h[i], round(ch$freq[i]*10)))
# }
# 
# ch.expand <-function(ch){
#   nc <- nchar(ch[1])
#   matrix(as.numeric(t(sapply(ch,substring,first=1:nc,last=1:nc))),length(ch))
# }
# CH <- array(NA,c(length(chh),4,2))
# CH[,1,1:2] <- ch.expand(substring(chh,1,2))
# CH[,2,1:2] <- ch.expand(substring(chh,3,4))
# CH[,3,1:2] <- ch.expand(substring(chh,5,6))
# CH[,4,1:2] <- ch.expand(substring(chh,7,8))
# 
# cat("M(t+1):\n")
# print(apply(CH, 2, function(x){sum(rowSums(x,na.rm=T)>0)}))
# 
# beta <- c(0.8038024, 0.8038024, 0.8038024, -0.9272944,-0.9272944,-0.9272944, -0.6435012, -0.6435012,-0.4115169,-0.4115169)
# beta <- c(rep(.86,3),rep(.1,3), rep(.2,2),rep(.3,2))
# beta <- log(beta/(1-beta))
# names(beta) <- c(rep("s",3),rep("g",3),rep("p",4))
# 
# tmp <- F.robust.loglik2( beta, CH)
#  print(tmp)
# 
# #fit1 <- optim( beta, F.robust.loglik2, ch=CH,  hessian=T, method="Nelder-Mead", control=list(reltol=1e-10, maxit=500))

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
