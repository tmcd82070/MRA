#' @title Compute spatially expliciate capture-recapture log likelihood
#'
#'  @description Given a set of activity center locations (mask), this routine computes the homogenous possion process liklihood. No covariates can be used at this time.
#'
#'  @param beta Coefficients to use when computing log likelihood.  This should be of length 3.  First element is
#'     log transform of density, second element is logit transform of g0, and third element is log transform of sigma.
#'
#'  @param ch  Data.frame of capture histories, size number of unique animals (nan) by number of secondaries (ns). The columns of ch should correpsond to each secondary occation, assumed to be in order. The rows of ch should correspond to a single animals capture history.  ch[i,j] = 0 if animal i was uncaptured during  occasion j. ch[i,j] = x (where x integer > 0) means animal i was captured at occasion j in trap x which has coordinates traps[x,].
#'
#'
#'
#'  @param traps  Data.frame of trap coordinates, size is K X 2. names(traps) should be 'x' and 'y'.
#'
#'
#'
#'  @param mask.pixel.area Size of every pixel in the habitat mask.  To calculate density, size of the
#'    habitat mask is needed.  Size of habitat mask is calculated by assuming all habitat points
#'    are centered in habitat pixels, and all habitat pixels are the same size. Size of habitat mask
#'    is calculated as number of pixels times \code{mask.pixel.area}.  Everything is unitless except this
#'    parameter.  density is per 1 unit of this parameter.
#'  @param type Type of trap, as defined by SECR. Note, currently only works for type='multi' the default.
#'
#'  @return The spatial (SECR) log likelihood. Note, at very last, log likelihood is
#'  multiply by -1 so this routine actually returns the negative of the log likelihood. Note, constants, such as the multinomial coefficient, are not included in the likelihood.
#'  This is done so that optim() can minimize and hessian estimation is correct.
#'
#'  @details   This routine using the notation of Borchera and Efford (2008) "Spatially
#'  Explicit Maximum Likelihood Methods for Capture–Recapture Studies", Biometrics, v64, 377-385.
#'
#' Another paper, with some ommisions I think, is "Density Estimation by spatially explicit capture-recapture: likelihood-based methods", by Efford, Borchers, and Byrom. This is unpublished, I think.
#' Either I don't fully comprehend, or this paper has some ommisions and is therefore less helpful.
#' This paper is for the homogenous Poisson case only.




## F.spat.loglik -----

F.spat.loglik <- function( beta, ch, traps, mask,mask.pixel.area, type="multi"){
  #
  # Compute closed population SECR likelihood for spatial capture models
  #
  # Input:
  #   beta = coefficients to minimize over.
  #   ch = matrix of capture histories.  ch[i,j] = 0 if animal i was uncaptured during
  #     occasion j. ch[i,j] = x (where x integer > 0) means animal i was captured
  #     at occasion j in trap x which has coordinates X[x,].  Eventually, this should be
  #     a SpatialPointDataFrame where points are coordinates of traps where each
  #     animal was caught and data frame contains animalID and occasion of capture. From
  #     that info, we could construct this ch matrix.
  #   traps = matrix of trap coordinates.  size is K X 2.  Eventually, this could be
  #     a SpatialPoints object.
  #   buffer = distance to buffer bounding box of trap locations to define polygon of
  #     possible home range centers.
  #
  # Value:
  # The spatial (SECR) log likelihood, Note, at very last, I multiply by -1 so this
  # routine actually returns the negative of the log likelihood.  This is done so
  # that optim() can minimize and hessian estimationis correct.
  #
  # I am using the notation of Borchera and Efford (2008) "Spatially Explicit Maximum Likelihood Methods
  # for Capture–Recapture Studies", Biometrics, v64, 377-385.
  #
  # Another paper, with some ommisions I think, is "Density Estimation by spatially explicit
  # capture-recapture: likelihood-based methods", by Efford, Borchers, and Byrom. This is unpublished, I think.
  # Either I don't fully comprehend, or this paper has some ommisions and is therefore less helpful.
  # This paper is for the homogenous Poisson case only.


  # Disallow 0 capture histories
  if( any(rowSums(ch>0) == 0) ) stop("Cannot have 0 capture histories.")

  # Basic sizes
  ns <- ncol(ch) # number of secondaries
  nan <- nrow(ch) # number of animals
  K <- nrow(traps)  # number of traps

  ## Pull parameters from beta
  D <- exp(beta[1])
  g0 <- exp(beta[2])/(1+exp(beta[2]))
  sigma <- exp(beta[3])



  ## mask should be data.frame with column names x and y
  ## these are the locations of the points in the habitat mask
  X <- mask

  T <- nrow(X)   # number of HR center locations

  # Compute distance from every potential home range center to every trap =====
  Xx <- matrix(X$x, T, K)
  Xy <- matrix(X$y, T, K)
  Kx <- matrix(traps$x, T, K, byrow=T)
  Ky <- matrix(traps$y, T, K, byrow=T)

  d <- sqrt( (Xx-Kx)^2 + (Xy-Ky)^2 )  # rows=HRlocation; cols=Trap location

  # Apply distance function ===================================================
  # Eventually, make this a call to a function that is passed in. i.e., add d.fund= parameter to this function call
  g <- g0*exp(-d^2/(2*(sigma^2)))

#   # This is how you plot the distance function for one trap
#     trp <- 1
#     par(pty="s")
#     image(Xxx,Xyy, matrix(g[,trp],length(Xxx)), main=paste("g for trap", trp))
#     points(traps[,1],traps[,2])

  # Make p_ks into a 3D array to account for occasions. dimensions are HR.centers (T) X Traps (K) X Session (ns)
  g <- array( g, c(T,K,ns) )

  # Compute trap hazard functions =============================================
  # This is trap hazard for animals at each HR center at each trap on occasion s.
  if( type == "proximity"){
    # Do nothing: THIS NOT CORRECT.
  } else if( type == "multi"){
    # Compute competing risks hazard-rate
    h.k <- -log(1 - g)  # trap hazards
    h.  <- apply(h.k,c(1,3),sum)  # sum over K = traps.  this is T x ns
    T_s <- 1  # Not needed, but just to remember could make occasion specific changes to p.s here.
    p.s  <- 1-exp(-T_s*h.)  # this is T x ns

#     # debugging
#     par(pty="s")
#     for( j in 1:ns ){
#       image(Xxx,Xyy, matrix(p.s[,j],length(Xxx)), main=j)
#       points(traps[,1],traps[,2])
#     }

    p_ks <- array( apply(p.s/h., 2, rep, times=K), c(T,K,ns)) # See note in paper.  I think T_s should be here.
                                                              # I think should be p.s/h. rather than (1-exp(-h.))/h.
                                                              # Plus, make matrix same size as individual hazard mat, so can multiply next
    p_ks <- p_ks*h.k  # this is p.s multiplied by proportion of trap hazard, this is T X K X ns

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

#
#   print(p.[50*9+9])
#   points(X[50*9+9,],col="blue",pch=16)
#
#    par(pty="s")
#    #image(Xxx,Xyy, matrix(p.,length(Xxx)), main="p.", xlim=c(300,400),ylim=c(300,400))
#    contour(Xxx,Xyy, matrix(p.,length(Xxx)), main="p.", levels=(.9), xlim=c(300,400),ylim=c(300,400))
#    points(traps[,1],traps[,2])

  # Compute Probability of capture histories =================================
  # Need p_ks and p.s here.
  tmp <- TRUE   # debugging only

  p.omega <- matrix(NA, T, nan)  # rows = HR centers, cols = animals
  for( i in 1:nan ){
    omega.i <- ch[i,]  # 0 histories checked above, so omega.i>0 for at least one element

    traps.i <- omega.i[omega.i > 0]
    occ.i <- which(omega.i > 0)
    p_ks.delta <- matrix(NA, T, length(traps.i))
    for(j in 1:length(traps.i)){
      p_ks.delta[,j] <- p_ks[,traps.i[j],occ.i[j]]
    }
    #ind <- c(outer(1:T, T*(omega.i[omega.i>0]-1) + T*ns*(which(omega.i>0)-1), "+")) # do this because don't want "cross indexes"
    #p_ks.delta <- matrix(p_ks[ind],T)   # this is T X (number of times i caught)
    p_ks.delta <- apply(p_ks.delta,1,prod)

#     if(tmp){
#       image(Xxx,Xyy, matrix(p_ks.delta,length(Xxx)), main=paste(omega.i,collapse="."))
#       points(traps[,1],traps[,2], pch=16)
#       cat("hit return...(0=don't plot any more, ie. finish)")
#       tmp <- readline() != 0
#     }

    if( any(omega.i==0) ){
      pUncap <- (1 - p.s)[,which(omega.i==0),drop=FALSE]# this is T X (num times i not captured)
      pUncap <- apply(pUncap,1,prod)
    } else {
      pUncap <- rep(1,T)
    }

#     if(tmp){
#       image(Xxx,Xyy, matrix(pUncap,length(Xxx)), main=paste(omega.i,collapse="."))
#       points(traps[,1],traps[,2], pch=16)
#       cat("That's pUncap.  Hit return..(0=don't plot any more, ie. finish)")
#       tmp <- readline() != 0
#     }


#    p.omega[,i] <- p_ks.delta*pUncap/p.   # this is T X 1; p.omega is T X nan
      p.omega[,i] <- p_ks.delta*pUncap   # this is T X 1; p.omega is T X nan
  }

#  debugging
#   for( i in 1:nan){
#     image(Xxx,Xyy, matrix(p.omega[,i],length(Xxx)), main=paste(ch[i,],collapse="."))
#     points(traps[,1],traps[,2], pch=16)
#     cat("hit return...(0=exit)")
#     if(readline()==0) break
#   }

  # Compute frequencies of capture histories.  Not sure this is needed =======
  ## This is only needed for the multinomial coefficient of the likelihood
  ## md <- ceiling(log10(nan))   # number of digits needed for nan
  ## n.freq <- table(apply(ch, 1, function(x,md){
  ##   x <- formatC(x,width=md,flag=0);
  ##   paste(x, collapse=".")},
  ##   md=2))

#   # The following computes capture frequencies without regard to traps, just captures
#   n.freq <- matrix(as.numeric(ch>0),nan)
#   n.freq <- table(apply(n.freq,1,paste,collapse=""))

#   cat("Capture frequencies: ========================\n")
#   print(n.freq)

  # In-Homogeneous point process case ========================================
  # More general than homogeneous. In future, allow D to vary by X. That is, D is T X 1 vector.
  # That is, D would be a function of a SpatialPointsDataFrame
#    Dp. <- D * p.
#    lambda <- sum(Dp.)  # same as a when D constant
#    f <- Dp. / lambda   # This is T X 1
#   ##L <- (factorial(nan) / prod(factorial(n.freq))) * prod(colSums(p.omega*f))   # Straight likelihood
#   logL <- lfactorial(nan) - sum(lfactorial(n.freq)) + sum(log(colSums(p.omega*f)))  # log Likelihood


# image(Xxx,Xyy, matrix(f,length(Xxx)), main="f")
# points(traps[,1],traps[,2], pch=16)
#
# plot(1:nan, colSums(p.omega*f), xlab="Animal ID/index")
# plot(1:nan, log(colSums(p.omega*f)), xlab="Animal ID/index")

  # Homogeneous point process case ===========================================
  # Compute a = integral of p.(X) over possible HR centers
  # total probability of capture (integeral of capture function
  # over all HR center locations)

  ## area give units to a
  ## currently assumes that area is the same for each point in the mask
  a <- sum(p.)*mask.pixel.area
  Da <- D*a
  #print(Da)
  # straight likelihood: L <- ((Da)^nan * exp(-Da) / factorial(nan)) * factorial(nan)/prod(factorial(n.freq)) * prod(colSums(p.omega*p.)) / a  # factorial(nan)'s cancel
  ##logL <- nan*log(Da) - Da  - sum(lfactorial(n.freq)) + sum(log(colSums(p.omega*p.))) - log(a)
  logL <- nan*log(Da) - Da  - lfactorial(nan) + sum(log(colSums(p.omega))) - nan*log(a)


  ##print(c(beta,area,-logL))

  return(-logL)

}



# ============================================================================
# Function calls

#library(secr)
#setwd(system.file('extdata', package='secr'))
#myCH <- read.capthist('capt.txt','trap.txt', fmt = 'XY')
#setwd("~/Programs/MRA/TestingVersion")

#secr0 <- secr.fit(myCH, model = g0~1, buffer = 100, trace = FALSE)

# # Compute SECR LL for one observation ----------
# obsn <- 1
# ll <- rep(NA,nrow(myCH))
# for( i in 1:1){
#   ch01 <- myCH[i,]
#   dim(ch01) <- c(length(obsn),ncol(myCH))
# #  ch01[1,2]<- 1
#   attr(ch01, "traps") <- attr(myCH,"traps")
#   attr(ch01, "session") <- attr(myCH,"session")
#   attr(ch01, "inject.time") <- rep(0, sum(ch01>0))
#   class(ch01)<- class(myCH)
#
#   tmp1 <- secr.fit(ch01, model = g0~1, start=secr0$fit$par, buffer = 100, details=list(LLonly=T))
#   cat(paste("SECR LL for observation", i, ":", tmp1, "\n"))
#   ll[i] <- tmp1
# }
# cat(paste("\nSum of SECR LL over all observations:", sum(ll,na.rm=T), "\n"))
#
#
# # Compute my LL for one observation ----------
# tmp <- F.spat.loglik( secr0$fit$par, ch01, attr(ch01,"traps"), 100 )
# cat(paste("My LL for observation", obsn, ":\n"))
# print(tmp)


# Prep and call secr likelihood
# tmp <- F.spat.loglik( secr0$fit$par, myCH, attr(myCH,"traps"), 100 )
# cat("Success...\n")
# print(tmp)

#fit1 <- nlminb(secr0$fit$par+ rnorm(3,0,.1), F.spat.loglik, ch=myCH, traps=attr(myCH,"traps"), buffer=100)
#fit2 <- optim(secr0$fit$par+ rnorm(3,0,.1), F.spat.loglik, ch=myCH, traps=attr(myCH,"traps"), buffer=100)
#print(fit2)
#
# tmp2 <- secr.fit(myCH, model = g0~1, start=fit1$par, buffer = 100, details=list(LLonly=T))
# cat("SECR LL at my params:\n")
# print(tmp2)
