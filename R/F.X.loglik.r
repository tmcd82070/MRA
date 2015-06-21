# Purpose ===========
#   An objective function for estimating activity center given captures
#

F.X.loglik <- function(x, h, traps, g0, sigma){
  # x = location of AC (2 X 1)
  # traps = T X 2 list of trap locations
  # g0 = distance function intercept
  # sigma = distance function std error
  # 

  xx <- x[1]
  xy <- x[2]
  
  d <- sqrt( (xx-traps[,1])^2 + (xy-traps[,2])^2 )  
  
  
  # Detection function
  g <- g0*exp(-0.5*(d/sigma)^2)
  
  h.k <- -log(1 - g)  # trap hazards, 
  h.  <- sum(h.k)  # sum over K = traps.  
  p.s  <- 1-exp(-h.)  # this is P(capture)
  
  PI <- p.s * h.k / (h. + 0.000000001)  # prob of capture by each trap. 
  PI <- c(1-p.s, PI)  # cell probabilities, with P(not cap) = first element
  
  PI.cells <- PI[h+1]
  PI.cells[PI.cells<.Machine$double.eps] <- .Machine$double.eps
  
  loglik <- log(PI.cells)

  #print(c(x,-sum(loglik)))
  
  -sum(loglik)
}


# ========================

# h <- popn$H[50,1,]
# 
# x.init <- colMeans(coordinates(X)[h,])
# 
# fit <- optim( x.init, F.X.loglik, h=h, traps=coordinates(X), g0=true.g0[1], sigma=true.sig[1])