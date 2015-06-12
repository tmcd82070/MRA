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

