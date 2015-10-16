#'  @title Compute probability of capture at primary and secondary levels. 
#'  
#'  @description Given estimated activity center locations, computes probability 
#'  of capture at least once during the primary occasion (i.e., computes pdot).  Optionally, 
#'  computes probability of capture by each trap each occasion (i.e., computes p_ks) and 
#'  probability of capture each secondary (i.e., computes p.s). 
#'  
#'  @param beta  A named vector of coefficient estimates in transformed (linear) space.   
#'    See \code{F.spat.loglik.X} for 
#'    description of structure.  Should have at least named elements "g0" and "sigma".
#'    
#'  @param traps A list of lists. It is assumed that the locations 
#'    of "on" traps varied amoungst occasions. In this case, each element of  the list is 
#'    matrix of "on" trap locations.  Each element is either a matrix of size K X 2 or 
#'    a \code{SpatialPoints} object of "on" trap locations. In other words,
#'    \code{length(traps)} is the number of occasions. 
#'    \code{traps[[j]]} is either a T(j) X 2 matrix or a SpatialPoints 
#'    object containing T(j) trap locations that were "on" during the j-th occasion. 
#'  
#'  @param aclocs A nan X 2 matrix of (estimated) activity center locations to use in evaluating the 
#'  likelihood. It is assumed that all locations are in valid habitat, so no checking against the 
#'  habitat mask is performed here.
#'  
#'  @param type The type of trap, as defined in SECR package
#'  
#'  @param return.occasionp Whether to return occasion specific probabilities (see value).  Probability
#'    of capture at least once during all occasions is always returned. 
#'  
#'  @return A list with either one or three components, depending on whether 
#'    \code{return.cell} is TRUE or FALSE.  Components of the 
#'    list are: 
#'    \item{pdot}{A vector of size \code{nan == nrow(aclocs)} containing probability of capture 
#'    for an animal at the corresponding location.  i.e., estimated probability of capture for 
#'    animal at location \code{aclocs[i,]} at least 
#'    once in the trap grid is \code{pdot[i]}}.  \code{pdot} is useful 
#'    for computing activity center locations, p.stars in robust designs, etc. }
#'    
#'    \item{p.s}{If \code{return.cellp} == TRUE, a matrix of size 
#'    \code{nrow(aclocs)} X \code{length(traps)}.  Cell (i,j) contains
#'    probability that individual at location i is captured during secondary occasion j.}
#'    
#'    \item{p_ks}{If \code{return.cellp} == TRUE, a matrix of size \code{nrow(aclocs)} X 
#'    \code{nrow(traps)} X ns.  Cell (i,j,k) contains probability of animal i being 
#'    caught in trap j during secondary occasion k. }

F.spat.capProbs <- function(beta, traps, aclocs, type="multi", return.occasionp=FALSE ){
  
  
  # Basic sizes
  T <- nan <- nrow(aclocs)
  ns <- length(traps)
  
  K.max <- max(unlist(lapply(traps,nrow)))  # maximum number of traps used during an occasion
  
  # Pull parameters from beta. # COULD EASILY MAKE THIS SECONDARY DEPENDENT BY PUTTING THIS IN K LOOP BELOW.
  # Pull using names so beta could be in any order
  g0 <- exp(beta["g0"])/(1+exp(beta["g0"]))
  sigma <- exp(beta["sigma"])
  
  
  # Compute occasion specific g = trap distance functions. =========
  g <- array(NA, c(T,K.max,ns))
  
  for(k in 1:ns){    
    K.k <- nrow(traps[[k]])  # number of traps used during kth occasion
    
    # Compute distance from AC locations to  every trap =====
    # Pull activity center locations and coordinates.  Can't do this outside k loop because of K.k
    ACx <- matrix(aclocs[,1], T, K.k)
    ACy <- matrix(aclocs[,2], T, K.k)
    
    Tx <- matrix(traps[[k]][,1], T, K.k, byrow=T)
    Ty <- matrix(traps[[k]][,2], T, K.k, byrow=T)
    
    d <- sqrt( (ACx-Tx)^2 + (ACy-Ty)^2 )  # rows=AC location = individual; cols=Trap location
    
    # Put distance function  3D array. Dimensions are AC.centers (T) X max traps used (K.max) X Session (ns)
    # keep in mind that if the number of traps used varies, there are NA columns in some pages.
    # specifically those past K.jk.
    g[1:T,1:K.k,k] <- g0*exp(-d^2/(2*(sigma^2)))
  }

    
  # Compute trap hazard functions =============================================
  # This is trap hazard for animals at each HR center at each trap on occasion s. 
  if( type == "proximity"){
    # Do nothing: THIS NOT CORRECT. 
    stop("Proximity trap types not implemented yet.")
  } else if( type == "multi"){
    # Compute competing risks hazard-rate.  Keep in mind potential for NA's in g
    h.k <- -log(1 - g)  # trap hazards 
    h.  <- apply(h.k,c(1,3),sum, na.rm=T)  # sum over K = traps.  this is T x ns
    T_s <- 1  # Not needed, but just to remember could make occasion specific changes to p.s here.
    p.s  <- 1-exp(-T_s*h.)  # this is T x ns
    
    #     # debugging
    #     par(pty="s")
    #     for( j in 1:ns ){
    #       image(Xxx,Xyy, matrix(p.s[,j],length(Xxx)), main=j)
    #       points(traps[,1],traps[,2])
    #     }
    
    if( return.occasionp){
      # I think should be p.s/h. rather than (1-exp(-h.))/h.
      # Plus, make matrix same size as individual hazard mat, so can multiply next
      p_ks <- array( apply(p.s/h., 2, rep, times=K.max), c(T,K.max,ns)) # See note in paper.  I think T_s should be here. 
      p_ks <- p_ks*h.k  # this is p.s multiplied by proportion of trap hazard, this is T X K X ns
    }
    
    #     image(Xxx,Xyy, matrix(p_ks[,1,2],length(Xxx)), col=topo.colors(20))
    #     points(traps[,1],traps[,2], pch=16, col=0)
    #     image(Xxx,Xyy, matrix(h.k[,1,2],length(Xxx)), col=terrain.colors(20))
    #     points(traps[,1],traps[,2], pch=16, col=0)
    #     
    #     image(Xxx,Xyy, matrix(p_ks[,55,2],length(Xxx)), col=topo.colors(20))
    #     points(traps[,1],traps[,2], pch=16, col=0)
    #     image(Xxx,Xyy, matrix(p_ks[,100,2],length(Xxx)), col=topo.colors(20))
    #     points(traps[,1],traps[,2], pch=16, col=0)
  }
    
  # compute p. ===============================================================
  # Probability of being caught at least once over ns occasions for all T locations
  # CHECK: This p. agrees with secr's pdot() function.
  p. <- 1 - apply(1-p.s,1,prod)  # this is T x 1  
  
  ans <- list(pdot=p.)
  
  if(return.occasionp){
    ans <- c(ans, list(p.s=p.s, p_ks=p_ks))
  }
  
  ans
}

# EXTRA CODE===========================
# YOU CAN EVENTUALLY ERASE THE REST OF THIS.  tHIS WAS THE CODE THAT WORKED. AS SOON AS THE ABOVE WORKS, 
# ERASE THIS. 
# !!!!!!!  
#   # Basic sizes
#   ns <- ncol(ch)
#   nan <- nrow(ch)
#   K <- nrow(traps)  # number of traps
#   
#   # Pull parameters from beta
#   D <- exp(beta[1])
#   g0 <- exp(beta[2])/(1+exp(beta[2]))
#   sigma <- exp(beta[3])
#   
#   
#   #   # Set locations where home ranges can be. ===================================
#   #   Xxx <- seq(min(traps[,1])-buffer, max(traps[,1])+buffer, length=50 )
#   #   Xyy <- seq(min(traps[,2])-buffer, max(traps[,2])+buffer, length=50 )
#   #   X <- expand.grid( x=Xxx, y=Xyy )
#   #   T <- nrow(X)   # number of HR center locations
#   
#   
#   # Compute distance from AC location to  every trap =====
#   T <- nan   # number of AC center locations = number of animals
#   ACx <- matrix(aclocs[,1], T, K)
#   ACy <- matrix(aclocs[,2], T, K)
#   Tx <- matrix(traps[,1], T, K, byrow=T)
#   Ty <- matrix(traps[,2], T, K, byrow=T)
#   
#   
#   d <- sqrt( (ACx-Tx)^2 + (ACy-Ty)^2 )  # rows=AC location = individual; cols=Trap location
#   
#   # Apply distance function ===================================================
#   # Eventually, make this a call to a function that is passed in. i.e., add d.fund= parameter to this function call
#   g <- g0*exp(-d^2/(2*(sigma^2)))
#   
#   #   # This is how you plot the distance function for one trap
#   #     trp <- 1
#   #     par(pty="s")
#   #     image(Xxx,Xyy, matrix(g[,trp],length(Xxx)), main=paste("g for trap", trp))
#   #     points(traps[,1],traps[,2])
#   
#   # Make p_ks into a 3D array to account for occasions. dimensions are HR.centers (T) X Traps (K) X Session (ns)
#   g <- array( g, c(T,K,ns) )
#   
#   # Compute trap hazard functions =============================================
#   # This is trap hazard for animals at each HR center at each trap on occasion s. 
#   if( type == "proximity"){
#     # Do nothing: THIS NOT CORRECT. 
#   } else if( type == "multi"){
#     # Compute competing risks hazard-rate
#     h.k <- -log(1 - g)  # trap hazards 
#     h.  <- apply(h.k,c(1,3),sum)  # sum over K = traps.  this is T x ns
#     T_s <- 1  # Not needed, but just to remember could make occasion specific changes to p.s here.
#     p.s  <- 1-exp(-T_s*h.)  # this is T x ns
#     
#     #     # debugging
#     #     par(pty="s")
#     #     for( j in 1:ns ){
#     #       image(Xxx,Xyy, matrix(p.s[,j],length(Xxx)), main=j)
#     #       points(traps[,1],traps[,2])
#     #     }
#     
#     p_ks <- array( apply(p.s/h., 2, rep, times=K), c(T,K,ns)) # See note in paper.  I think T_s should be here. 
#     # I think should be p.s/h. rather than (1-exp(-h.))/h.
#     # Plus, make matrix same size as individual hazard mat, so can multiply next
#     p_ks <- p_ks*h.k  # this is p.s multiplied by proportion of trap hazard, this is T X K X ns
#     
#     #     image(Xxx,Xyy, matrix(p_ks[,1,2],length(Xxx)), col=topo.colors(20))
#     #     points(traps[,1],traps[,2], pch=16, col=0)
#     #     image(Xxx,Xyy, matrix(h.k[,1,2],length(Xxx)), col=terrain.colors(20))
#     #     points(traps[,1],traps[,2], pch=16, col=0)
#     #     
#     #     image(Xxx,Xyy, matrix(p_ks[,55,2],length(Xxx)), col=topo.colors(20))
#     #     points(traps[,1],traps[,2], pch=16, col=0)
#     #     image(Xxx,Xyy, matrix(p_ks[,100,2],length(Xxx)), col=topo.colors(20))
#     #     points(traps[,1],traps[,2], pch=16, col=0)
#   }
#   
#   # compute p. ===============================================================
#   # Probability of being caught at least once over ns occasions for all T locations
#   # CHECK: This p. agrees with secr's pdot() function.
#   p. <- 1 - apply(1-p.s,1,prod)  # this is T x 1
#   
#   #  assign("tmp.mypdot", p., pos=.GlobalEnv)
#   # p. has been checked many times and agrees with SECR. 
#   #
#   #   print(p.[50*9+9])
#   #   points(X[50*9+9,],col="blue",pch=16)
#   # 
#   #    par(pty="s")
#   #    #image(Xxx,Xyy, matrix(p.,length(Xxx)), main="p.", xlim=c(300,400),ylim=c(300,400))
#   #    contour(Xxx,Xyy, matrix(p.,length(Xxx)), main="p.", levels=(.9), xlim=c(300,400),ylim=c(300,400))
#   #    points(traps[,1],traps[,2])  
#   
#   
#   
#   
#   
# }