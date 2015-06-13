# Fit Huggins closed model via maximum likelihood
#
# Objective function ----
F.hug.loglik <- function(beta, ch){
  
  ns <- ncol(ch)
  nan <- nrow(ch)
  
  # Evaluate p model
  p <- exp(beta[1])
  p <- p / (1 + p)
  p <- matrix( rep(p, ns), nan, ns, byrow=T )
#  print(p)
  
  # Evaluate p part of likelihood ----
  p.part <- ch*log(p) + (1-ch)*log(1-p)
#  print(p.part)
  p.part <- sum( p.part ) 
#  print(p.part)
  
  # Evaluate prob of overal capturing each animal  ----
  # this is probability of observing a capture history
  row.capp <- 1 - apply( 1-p, 1, prod )
  row.capp <- log( row.capp )
  row.capp <- sum( row.capp )
#  print(row.capp)  

  # Full likelihood ----
  ll <- p.part - row.capp
  -ll
  
}

# Function calls ================================================

#F.closed.loglik( c(-1,-1), ch, 100)
# fit1 <- optim( c(0), F.hug.loglik, ch=ch, hessian=T, method="Brent", lower=-4, upper=4,  control=list(reltol=1e-10, maxit=500))
# 
# library(mra)
# fit2 <- F.huggins.estim( ~1, histories=ch)
# 
# print( rbind( c(fit1$par, sqrt(1/fit1$hessian), fit1$value), c(fit2$capcoef, fit2$se.capcoef, fit2$loglik)))

# All this checks out and verifies that F.hug.loglik is correct. 
# when you do the following, with a robust design data set in CH, you get 
# the same p's when you fit a robust design with no constraints (s(t),g'(t)=g''(t),p(t)=c(t)).  that is, 
# when the robust design is unconstrained, and no information passes between primaries, 
# you can estimate p's within primaries by separate Huggins models. 
# 
# library(mra)
# 
# for(j in 1:5){
#   ch1 <- CH[,j,]
#   ch1 <- ch1[,!is.na(colSums(ch1))]
#   ch1 <- ch1[rowSums(ch1>0)>0,]
#   print(dim(ch1))
#   fit1 <- optim( c(0), F.hug.loglik, ch=ch1, hessian=T, method="Brent", lower=-4, upper=4,  control=list(reltol=1e-10, maxit=500))
#   fit2 <- F.huggins.estim( ~1, histories=ch1)
# 
#   print(c(fit1$par,fit1$value))
#   print(fit2)
#   cat("-------------------------------------------------------\n")
# }
# 
# # Write out a single .inp to test in MARK
# j = 4
# ch1 <- CH[,j,]
# ch1 <- ch1[,!is.na(colSums(ch1))]
# ch1 <- ch1[rowSums(ch1>0)>0,]
# tmp <- paste(paste(ch1[,1],ch1[,2],ch1[,3],ch1[,4],ch1[,5],sep=""), "1;")
# write.table(tmp, "RobOcc4.inp", sep="",row.names=F, col.names=F, quote=F)
# fit1 <- optim( c(0), F.hug.loglik, ch=ch1, hessian=T, method="Brent", lower=-4, upper=4,  control=list(reltol=1e-10, maxit=500))
# 
# # When I fitted a Huggins model to data in RobOcc4.inp (p(.)), I get p.eta = .4364 and -2LL=3347.5956.  This 
# # matches what I get in R (i.e., fit1$par = .4362 and 2*fit1$val = 3347.5956)