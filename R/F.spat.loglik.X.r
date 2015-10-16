#' @title Compute spatial capture-recapture log likelihood given estimated activity center locations. Use this 
#' during the  "M" step of the EM algorithm
#' 
#'  @description Given a set of activity center locations, this routine computes the 
#'  liklihood for a set of coefficients.  
#'  
#'  @param beta A named vector of coefficients in the transformed (linear) space 
#'  to use when computing log likelihood.  This should be of length 3.  
#'  First element is log transform of density (named "D"), second element is logit transform of 
#'  g0 (named "g0"), 
#'  and third element is log transform of sigma (named "sigma").
#'  
#'  @param ch  Matrix of capture histories, size nan X ns.  ch[i,j] = 0 if animal i was uncaptured during 
#'     occasion j. ch[i,j] = x (where x integer > 0) means animal i was captured 
#'     at occasion j in trap x which has coordinates traps[x,].  Eventually, this should be 
#'     a SpatialPointDataFrame where points are coordinates of traps where each 
#'     animal was caught and data frame contains animalID and occasion of capture. From
#'     that info, we could construct this ch matrix.
#'     
#'  @param traps  Matrix of trap coordinates, size is K X 2.  
#'  
#'  @param aclocs A nan X 2 matrix of (estimated) activity center locations to use in evaluating the 
#'  likelihood. It is assumed that all locations are in valid habitat, so no checking against the 
#'  habitat mask is performed here.
#'  
#'  @param mask.pixel.area Size of every pixel in the habitat mask.  To calculate density, size of the
#'    habitat mask is needed.  Size of habitat mask is calculated by assuming all habitat points 
#'    are centered in habitat pixels, and all habitat pixels are the same size. Size of habitat mask 
#'    is calculated as number of pixels times \code{mask.pixel.area}.  Everything is unitless except this
#'    parameter.  Density is per 1 unit of this parameter.  For example, if \code{mask.pixel.area} is 
#'    in hectares, final density parameters are number of individuals per hectare.  It is incumbent 
#'    on the user to understand this and get it right.  
#'        
#'  @param type Type of trap, as defined by SECR. 
#'     
#'  @return The negative of the spatial (SECR) log likelihood. Negative is returned so 
#'    so  minimizers like \code{optim}  can be applied directly and (more importantly) 
#'    the hessian estimation is correct.
#'    
#'  @details   This routine using the notation of Borchera and Efford (2008) "Spatially 
#'  Explicit Maximum Likelihood Methods for Capture–Recapture Studies", Biometrics, v64, 377-385.
#'
#' Another paper, with some ommisions I think, is "Density Estimation by spatially explicit 
#' capture-recapture: likelihood-based methods", by Efford, Borchers, and Byrom. This is unpublished, I think.
#' Either I don't fully comprehend, or this paper has some ommisions and is therefore less helpful. 
#' This paper is for the homogenous Poisson case only. 
#' 
#' @examples 


F.spat.loglik.X <- function( beta, ch, traps, aclocs, mask.pixel.area=1, type="multi" ){

  # Disallow 0 capture histories
  if( any(rowSums(ch>0) == 0) ) stop("Cannot have 0 capture histories.")

  # Check for size restraints
  if( nrow(ch) != nrow(aclocs))  stop("Number of rows in ch must equal number of rows in aclocs.")
  
  # Compute probabilties given beta and trap configuration at all locations ============
  capProbs <- F.spat.capProbs(beta, traps, aclocs, type, return.cellp=TRUE)
  p. <- capProbs$pdot
  p.s <- capProbs$p.s
  p_ks <- capProbs$p_ks
  
 
  # Compute Probability of capture histories ===========================================
  # Need p_ks and p.s and ch here.
  # tmp <- TRUE   # debugging only
  
  p.omega <- rep(NA,nan)  # length = # animals = one p.omega for each animal's location
  for( i in 1:nan ){
    omega.i <- ch[i,]  # 0 histories checked above, so omega.i>0 for at least one element
    
    traps.i <- omega.i[omega.i > 0]
    occ.i <- which(omega.i > 0)
    p_ks.delta <- rep(NA, length(traps.i))
    for(j in 1:length(traps.i)){
      p_ks.delta[j] <- p_ks[i,traps.i[j],occ.i[j]]  
    }
    #ind <- c(outer(1:T, T*(omega.i[omega.i>0]-1) + T*ns*(which(omega.i>0)-1), "+")) # do this because don't want "cross indexes"
    #p_ks.delta <- matrix(p_ks[ind],T)   # this is T X (number of times i caught)
    p_ks.delta <- prod(p_ks.delta)
  
#     if(tmp){
#       image(Xxx,Xyy, matrix(p_ks.delta,length(Xxx)), main=paste(omega.i,collapse="."))
#       points(traps[,1],traps[,2], pch=16)
#       cat("hit return...(0=don't plot any more, ie. finish)")
#       tmp <- readline() != 0
#     }
    
    if( any(omega.i==0) ){
      pUncap <- (1 - p.s[i,])[which(omega.i==0),drop=FALSE]        # this has length (num times i not captured)
      pUncap <- prod(pUncap)
    } else {
      pUncap <- 1
    }

#     if(tmp){
#       image(Xxx,Xyy, matrix(pUncap,length(Xxx)), main=paste(omega.i,collapse="."))
#       points(traps[,1],traps[,2], pch=16)
#       cat("That's pUncap.  Hit return..(0=don't plot any more, ie. finish)")
#       tmp <- readline() != 0
#     }
    
    
#    p.omega[,i] <- p_ks.delta*pUncap/p.   # p.omega is T X nan
      p.omega[i] <- p_ks.delta*pUncap   # p.omega is length nan
  }  
  
#  debugging
#   for( i in 1:nan){
#     image(Xxx,Xyy, matrix(p.omega[,i],length(Xxx)), main=paste(ch[i,],collapse="."))
#     points(traps[,1],traps[,2], pch=16)
#     cat("hit return...(0=exit)")
#     if(readline()==0) break
#   }
  
  # Compute frequencies of capture histories.  Not sure this is needed =======
#   md <- ceiling(log10(nan))   # number of digits needed for nan
#   n.freq <- table(apply(ch, 1, function(x,md){
#     x <- formatC(x,width=md,flag=0); 
#     paste(x, collapse=".")}, 
#     md=2))

#   # The following computes capture frequencies without regard to traps, just captures
#   n.freq <- matrix(as.numeric(ch>0),nan)
#   n.freq <- table(apply(n.freq,1,paste,collapse=""))

#   cat("Capture frequencies: ========================\n")
#   print(n.freq)
  
  # In-Homogeneous point process case ========================================
  # More general than homogeneous. In future, allow D to vary by X. That is, D is T X 1 vector.
  # That is, D would be a function of a SpatialPointsDataFrame
#    Dp. <- D * p. * mask.pixel.area / nan
#    lambda <- sum(Dp.)  # same as a when D constant
#    f <- Dp. / lambda   # This is T X 1
#   ##L <- (factorial(nan) / prod(factorial(n.freq))) * prod(colSums(p.omega*f))   # Straight likelihood
#   logL <-  nan*log(lambda) - lambda + sum(log(p.omega)) + sum(log(f)) - sum(log(p.))  # log Likelihood


# image(Xxx,Xyy, matrix(f,length(Xxx)), main="f")
# points(traps[,1],traps[,2], pch=16)
# 
# plot(1:nan, colSums(p.omega*f), xlab="Animal ID/index")
# plot(1:nan, log(colSums(p.omega*f)), xlab="Animal ID/index")

  # Homogeneous point process case ===========================================
  # Compute a = integral of p.(X) over possible HR centers
  # total probability of capture (integeral of capture function
  # over all HR center locations)
  
  D <- exp(beta["D"])
  
   a <- sum(p.)  # this assumes all pixels in the same size
# 
   Da <- D*a
#   #print(Da)
#   # straight likelihood: L <- ((Da)^nan * exp(-Da) / factorial(nan)) * factorial(nan)/prod(factorial(n.freq)) * prod(colSums(p.omega*p.)) / a  # factorial(nan)'s cancel
#   ##logL <- nan*log(Da) - Da  - sum(lfactorial(n.freq)) + sum(log(colSums(p.omega*p.))) - log(a)                                                                                              
#   logL <- nan*log(Da) - Da  + sum(log(p.omega)) - nan*log(a)                                                                                              
  logL <- nan*log(D*mask.pixel.area) - (D*mask.pixel.area*a/nan)  + sum(log(p.omega))                                                                                               
  
#   print(c(beta,-logL ))

  # Return negative log likelihood
  -logL
  
}

  
  
# ============================================================================
# Function calls 

 library(secr)
# setwd(system.file('extdata', package='secr'))
# myCH <- read.capthist('capt.txt','trap.txt', fmt = 'XY')
# setwd("~/Programs/MRA/TestingVersion")

# m <- make.grid(nx=1, ny=1, spacing = 30, originxy = c(500,500))
# m <- make.mask(m, nx=1,ny=1)
# 
# secr0 <- secr.fit(myCH, model = g0~1, mask=m) 


# # Compare SECR and my likelihood observation by observation ----------
# m <- make.grid(nx=1, ny=1, spacing = 30, originxy = c(500,500))
# m <- make.mask(m, nx=1,ny=1)
# 
# obsn <- 2
# ll <- matrix(NA, nrow(myCH), 2)
# for(i in 1:nrow(myCH)){
#   ch01 <- myCH[i,]
#   dim(ch01) <- c(1,ncol(myCH))
#   attr(ch01, "traps") <- attr(myCH,"traps")
#   attr(ch01, "session") <- attr(myCH,"session")
#   attr(ch01, "inject.time") <- rep(0, sum(ch01>0))
#   class(ch01)<- class(myCH)
#   
#   tmp1 <- secr.fit(ch01, model = g0~1, start=secr0$fit$par, mask=m, details=list(LLonly=T))  
#   cat(paste("SECR LL for observation", i, ":", tmp1, "\n"))
#   ll[i,1] <- tmp1
# 
#   # Compute my LL for one observation ----------
#   tmp <- F.spat.loglik.X( secr0$fit$par, ch01, attr(ch01,"traps"), m )
#   cat(paste("My LL for observation", i, ":", tmp, "\n"))
#   ll[i,2] <- tmp
# }
# 
# cat(paste("\nSum of SECR LL over all observations:", colSums(ll,na.rm=T), "\n"))
# plot(ll[,1],ll[,2], xlab="SECR log like", ylab="My log like", pch=16)
# 
# ## Conclusion: when mask has one point in it, SECR and F.spat.loglik.X differ by a constant 
# ## for every observation.  Plot above shows a straight line, and apply(ll,1,function(x){diff(x)}) is constant.
# 
# 
# # Compare SECR and my likelihood for different values of parameters ----------
# m <- make.grid(nx=1, ny=1, spacing = 30, originxy = c(500,500))
# m <- make.mask(m, nx=1,ny=1)
# 
# 
# ch01 <- myCH[2,]
# ch01 <- myCH
# dim(ch01) <- c(nrow(ch01),ncol(myCH))
# attr(ch01, "traps") <- attr(myCH,"traps")
# attr(ch01, "session") <- attr(myCH,"session")
# attr(ch01, "inject.time") <- rep(0, sum(ch01>0))
# class(ch01)<- class(myCH)
# 
# secr0 <- secr.fit(ch01, model = g0~1, mask=m, start=c(2.953150, -4.652265, 16.590712)) 
# 
# acs <- matrix( as.numeric(m), nrow(ch01), 2, byrow = TRUE)
# 
# fit0 <- optim( c(2.953150, -4.652265, 16.590712), F.spat.loglik.X, ch=ch01, aclocs=acs, traps=traps(ch01), mask.pixel.area=attr(m,"area"), method="L-BFGS-B", hessian = TRUE, control=list(factr=5e9, pgtol=1e-8, maxit=1000))
# 
# par.vec <- secr0$fit$par
# D.vec <- seq(0.5*par.vec[1], 1.5*par.vec[1], length=40)
# ll <- matrix(NA, length(D.vec), 2)
# for(i in 1:length(D.vec)){
# 
#   tmp1 <- secr.fit(ch01, model = g0~1, start=c(D.vec[i],par.vec[2:3]), mask=m, details=list(LLonly=T))  
#   cat(paste("SECR LL for D=", D.vec[i], ":", tmp1, "\n"))
#   ll[i,1] <- tmp1
#   
#   tmp <- F.spat.loglik.X( c(D.vec[i],par.vec[2:3]), ch01, traps(ch01), acs )
#   cat(paste("  My LL for D=", D.vec[i], ":", tmp, "\n"))
#   ll[i,2] <- tmp
# }
# 
# plot(range(D.vec), range(ll), xlab="log(Density)", ylab="log like", type="n")
# lines(D.vec, ll[,1], col="black", pch=16, type="b")
# lines(D.vec, ll[,2], col="red", pch=16, type="b")
# abline(v=par.vec[1])

## Conclusion: At time of writing (14 oct 15), SECR and our likelihood do not agree in D dimension.  
## When minimized, SECR and our likelihood produce same detection parameters, but different densities. 
## From above code, secr0$fit$par[2:3] == fit0$par[2:3], but secr0$fit$par[1] != fit0$par[1]. 
## Jared found same thing by different methods. 

# # Check p. ---------------------------------
# get.pdot <- function(obj){
#   nocc <- ncol(obj$capthist)
#   trap <- attr(obj$capthist,'traps')
#   mask <- obj$mask
#   sigma <- exp(coef(obj)$beta[3])
#   g0 <- exp(coef(obj)$beta[2])/(1+exp(coef(obj)$beta[2]))
#   
#   P. <- pdot(mask,trap,detectpar=list(g0=g0,sigma=sigma),noccasions=nocc)
#   D <- nrow(obj$capthist)/sum(P.)
#   print(paste('n/a = ',D))
#   print(paste('log(n/a) = ',log(D)))
#   print(paste('SECR D=',exp(coef(obj)$beta[1])))
#   print(paste('SECR log(D)=',coef(obj)$beta[1]))
#         
#   P.
# }
# 
# m <- make.grid(nx=1, ny=1, spacing = 30, originxy = c(500,500))
# m <- make.mask(m, nx=1,ny=1)
# 
# secr0 <- secr.fit(myCH, model = g0~1, mask=m, trace = TRUE)
# tmp.pdot <- get.pdot(secr0)
# 
# acs <- matrix( as.numeric(m), nrow(myCH), 2, byrow = TRUE)
# tmp <- F.spat.loglik.X( secr0$fit$par, myCH, traps(myCH), acs )
# 
# 
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