#'  @title Likelihood location contains activity center. 
#'  
#'  @description Computes the likelihood that a particular location is 
#'  a (latent, unobserved) activity center location.  This is 
#'   the objective function for estimating activity center given captures. 
#'
#'  @param x A 2 X 1 vector containing (x,y) coordinates of the location (i.e., the potential 
#'    activity center location). 
#'    
#'  @param h A ns X 1 vector of trap capture history.  Element i is either 0 if individual went
#'    uncaptured, or the row number in \code{traps} of the trap that caught the individual.
#'    
#'  @param traps A list of trap location matricies.  Length of \code{traps} is ns = number of occasions. 
#'    Each element of \code{traps} is a T X 2 matrix of trap locations.  Trap location that caught 
#'    the individual associated with history \code{h} during occasion j is \code{traps[[j]][h[j],]}.
#'    
#'  @param g0 Vector of distance function intercept values
#'  
#'  @param sigma Vector of distance function width values
#'  
#'  @return The negative log likelihood that location \code{x} is the activity center 
#'    of the animal with history \code{h} in 
#'    traps defined by \code{traps},  given distance function parameters \code{g0} and \code{sigma}. 
#'    

F.X.loglik <- function(x, h, traps, g0, sigma){
  # x = location of AC (2 X 1)
  # traps = a nsecondary length list of T X 2 trap locations
  # g0 = distance function intercept
  # sigma = distance function std error
  # 

  # Get number of secondary occasions
  # Note: length(traps) should equal sum(!is.na(h)), but I'm not checking this
  nsecondary <- length(traps)
    
  # This would be the place to fix up g0 and sigma parameter vectors; however, 
  # I am currently allowing only time variation, so no fix up is needed.  g0 and sigma are length one. 
  

  # Loop over secondaries
  PI.cells <- rep(NA, nsecondary)
  for( j in 1:nsecondary){

    # each element of traps is TX2 trap coordinate matrix
    tr <- coordinates(traps[[j]])
    
    # Compute distances
    d <- sqrt( (x[1]-tr[,1])^2 + (x[2]-tr[,2])^2 )  

    # Compute Detection function
    # NOTE: half normal is hard-wired here
    g <- g0*exp(-0.5*(d/sigma)^2)
    
    h.k <- -log(1 - g)  # trap hazards, 
    h.  <- sum(h.k)  # sum over K = traps.  
    p.s  <- 1-exp(-h.)  # this is P(capture)
    
    PI <- p.s * h.k / (h. + 0.000000001)  # prob of capture by each trap. 
    PI <- c(1-p.s, PI)  # cell probabilities, with P(not cap) = first element
    
    PI.cells[j] <- PI[h[j] + 1]
    
  }

  PI.cells[PI.cells<.Machine$double.eps] <- .Machine$double.eps
  
  loglik <- log(PI.cells)

  #print(c(x,-sum(loglik)))
  
  -sum(loglik)
}


# ========================
# 
# THIS TESTING CODE DOES NOT WORK ANYMORE BECAUSE i CHANGED THE STRUCTURE OF TRAPS
#
# h <- H[1,3,]
# 
# x.init <- colMeans(coordinates(X)[h,])
# 
# h <- H[1,1,]
# 
# # (F.X.loglik(x.init, h, coordinates(X), .8, 15))
# # (F.X.loglik(c(116,45), h, coordinates(X), .8, 15))
# # (F.X.loglik(c(130,45), h, coordinates(X), .8, 15))
# 
# fit <- optim( c(50,70), F.X.loglik, h=h, traps=coordinates(X), g0=true.g0[1], sigma=true.sig[1])
# 
# ll <- NULL
# xx <- 50:300
# for(xxx in xx){
#   ll <- c(ll, F.X.loglik(c(xxx,xxx), h, coordinates(X), .8, 15))
# }
# plot(xx,ll, type="l")
# 
# # Conclusion----
# #  When h=c(0,0,0) (i.e., uncaught during primary), likelihood for X is zero (uncatchable) at ~3*sigma off the 
# #  grid. 