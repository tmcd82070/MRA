#' @title Compute spatial robust-design log likelihood given estimated activity center locations. Use this 
#' during the  "M" step of the EM algorithm
#' 
#'  @description Given a set of activity center locations, this routine computes the 
#'  spatial robust design likelihood for a set of coefficients.  
#'  
#' @param beta A named vector of coefficient estimates.  
#' Names and number of elements should be as follows:  
#' . "surv?" (survival) = (nprimary-1) elements corresponding to intervals between primary occasions
#' . "gp?" (gamma prime, emigration) = (nprimary-1) elements corresponding to intervals between primary occasions.
#' . (Optional) "gpp?" (gamma prime prime, immigration) = (nprimary - 2) elements corresponding to intervals between 
#'    2nd primary and the last.  Recall there is no immigration parameter for the interval between primary 1 and 2. 
#' . "D?" (density) = nprimary elements. 
#' . "g0?" (capture probability at traps) = nprimary elements
#' . "sigma?" (capture probability width) = nprimary elements.
#' 
#' Current implementation allows constant or complete time variation (small t).  For constant model, 
#' beta has length 5 or 6 (6 if Markov emigration). For small t model, length of beta is either (5*nprimary - 2) or
#' (6*nprimary - 4).  
#' For example, in a problem with 4 primaries and 5 secondaries in each primary, and no gpp elements, beta 
#' should have length 18 = 3 + 3 + 4 + 4 + 4 and be named \code{c("surv1", "surv2", surv3", "gp1", "gp2",     
#' "gp3", "D1", "D2", "D3", "D4", "g01", "g02", "g03", "g04", "sigma1", "sigma2", "sigma3", "sigma4")}.  Note that 
#' order does not matter because elements are all named. 
#' If length of \code{beta} is 5*nprimary - 2, this routine assume that the 'random emigration' 
#' model is being fitted wherein emigration probability (gp) and immigration 
#' probability (gpp) are equal. If  \code{length(beta)} equals 6*nprimary - 4, this routine assumes the 
#' 'Markov emigration' model is being fitted 
#' wherein emigration (gp) and immigration (gpp) depend upon whether an individual is on or off the 
#' study area (see robust design literature for more details). 
#'
#' @parm ac.locs A nan X nprimary X 2 array of estimated activity center locations.  
#' \code{ac.locs[i,j,]} = c(x,y) = coordinates of activity center for animal i
#' during primary occasion j.  It is assumed that all locations are in valid habitat, so no checking against the 
#'  habitat mask is performed here.
#'  
#' @param ch A 3-D array of capture histories.  Size is nan X nprim X nsecondaries.  Missing secondaries 
#' are specified with columns of all NA's (i.e., if the 4th secondary of the ith primary was not done, 
#' \code{all(is.na(ch[,i,4])) == TRUE}).  Cells in the array contain trap number that caught the individual.
#' Trap numbers are rows in the \code{traps} matrix. ch[i,j,k] = 0 if animal i was uncaptured during 
#'     secondary occasion k of primary occasion j. ch[i,j,k] = x (where x integer > 0) means animal i was captured 
#'     during secondary occasion k of primary occasion j in trap x which has coordinates \code{traps[[j]][[k]][x,]}.
#' 
#' @param traps Either a \code{SpatialPoints} object, or a list of lists containing trap locations for each secondary. 
#' If \code{traps} is a \code{SpatialPoints} object, it is assumed that all traps were "on" during all secondary occasions (recall secondary
#' occasions are embedded in primary occasions).  If \code{traps} is a list of lists, it is assumed that the locations 
#' of "on" traps varied amoungst secondary occasions. In this case, the each element in the first layer of 
#' lists is a sub-lists of "on" trap locations 
#' for each primary occasion.  Each element in the second layer of lists is either a matrix or a \code{SpatialPoints} object of 
#' trap locations "on" during a secondary of a 
#' particular primary. In other words,
#' \code{length(traps)} is the number of primary occasions. 
#' \code{length(traps[[i]])} is number of secondaries during the i-th primary.  \code{traps[[i]][[j]]} is either a T(i,j) X 2 matrix
#' or a SpatialPoints object containing T(i,j) trap locations that were "on" during the j-th secondary occasion of the i-th primary. The trap location that caught 
#' animal m during the j-th secondary of the i-th primary is \code{traps[[i]][[j]][ch[m,i,j],]}.
#' 
#'  @param hab.mask A SpatialPixels or SpatialGrid object (or their *DataFrame analogs) containing allowable locations for activity centers. 
#'    Exactly like masks in the SECR package, \code{hab.mask} defines the outer limit of integration, defines 
#'    sites that are habitat and thus can be occupied, descritizes habitat covariates for use in models.
#'    The current implementation does not allow the habitat mask to change between primaries or secondaries. 
#'    The mask is constant throughout the study.  
#'    Note that the area of each pixel (grid cell) is computed
#'    as \code{prod(hab.mask@grid@cellsize)}, so \code{hab.mask} must be projected (e.g., in UTM's) 
#'    and coordinates in \code{traps} must be in the same units.  Resulting density estimates 
#'    are number of animals per squared unit of this system.  E.g., if using UTM coordinates in 
#'    units of kilometers, density comes out as inidividuals per square kilometer. 
#'     
#'  @return The spatial (SECR) log likelihood. Note, at very last, log likelihood is 
#'  multiply by -1 so this routine actually returns the negative of the log likelihood.  
#'  This is done so that optim() can minimize and hessian estimation is correct.
#'  
#'  @details   This routine using the notation of Borchera and Efford (2008) "Spatially 
#'  Explicit Maximum Likelihood Methods for Capture–Recapture Studies", Biometrics, v64, 377-385.
#'
#' Another paper, with some ommisions I think, is "Density Estimation by spatially explicit 
#' capture-recapture: likelihood-based methods", by Efford, Borchers, and Byrom. This is unpublished, I think.
#' Either I don't fully comprehend, or this paper has some ommisions and is therefore less helpful. 
#' This paper is for the homogenous Poisson case only.
#'  
#' Eventually, the \code{ch} array should be 
#'     a SpatialPointDataFrame where points are coordinates of traps where each 
#'     animal was caught and data frame contains animalID and occasion of capture. From
#'     that info, we could construct this ch matrix.
# 
# Objective function -----
F.spat.robust.loglik.X <- function( beta, ch, ac.locs, traps, hab.mask ){


  if( !inherits(ch,"array") ) stop("ch must be an array")
  if( length(dim(ch)) != 3) stop(paste("ch must have 3 dimensions.", length(dim(ch)), "found."))
  if( !inherits(hab.mask, "SpatialPixels") & !inherits(hab.mask, "SpatialGrid")) stop("hab.mask must be either a SpatialPixels or SpatialGrid object.")
  
  # Find dimensions                                 
  d <- dim(ch)
  nan <- d[1]
  nprimary <- d[2]
  nsecondary <- apply(ch,2,function(x){
    xx<-apply(x,2,function(y){all(is.na(y))})
    sum(!xx)
  })
  
  # Evaluate the model to get real parameters
  f.real.model<-function(beta, nprim){
    # Normally, you would have a X matrix here.
    # For now, use names of beta to evaluate model.
    #
    # Current implementation allows only time variation (small t or dot model are subsets).
    # This returns linear predictors.  Link functions are applied later.
    
    s.parms <- grep("^S",names(beta))
    gp.parms <- grep("^gp",names(beta))
    gpp.parms <- grep("^gpp",names(beta))
    
    
    # For now, same SECR parameters in each primary session
    D.parms <- grep("^D",names(beta))
    g0.parms <- grep("^g0",names(beta))
    sigma.parms <- grep("^sigma",names(beta))
    
    # s parameters
    s <- beta[s.parms]
    if( length(s) != (nprim-1) ){
      # Just use first 
      s <- rep(s[1],nprim-1)
    }
    
    # gp parameters
    gp <- beta[gp.parms]
    if( length(gp) == (nprim-2)){
      # constrain last g to second to last g.
      gp <- c(NA,gp,gp[nprim-2])
    } else if(length(gp) == (nprim-1)){
      gp <- c(NA,gp)  # use all gammas
    } else if( length(gp) == 0 ){
      # Constrain gamma to zero
      gp <- c(NA,rep(0,nprim-1))
    } else if( length(gp) != (nprim-1)){
      # Assume constant model
      gp <- c(NA,rep(gp[1],(nprim-1))) # first gamma not possible.
    } 
    # For rest of routine to work g must have length nprim.  gp[1]=NA
    
    # gpp parameters
    gpp <- beta[gpp.parms]
    if(length(gpp) == (nprim-2)){
      gpp <- c(NA,NA,gpp)  # use all gammas
    } else if( length(gpp) == 0 ){
      # Constrain gpp to equal gp
      gpp <- gp
      gpp[2] <- NA  # except gpp[2] is not possible to estimate
    } else {
      # Assume constant model
      gpp <- c(NA,NA,rep(gpp[1],(nprim-2))) # first two gpp's not possible.
    } 
    # For rest of routine to work gpp must have length nprim.  gpp[1:2]=NA
    
    # g0 parameters
    g0 <- beta[g0.parms]
    if( length(g0) != nprim ){
      # Just use first 
      g0 <- rep(g0[1],nprim)
    }    
    
    # sigma parameters
    sigma <- beta[sigma.parms]
    if( length(sigma) != nprim ){
      # Just use first 
      sigma <- rep(sigma[1],nprim)
    } 

    # Density parameters
    D <- beta[D.parms]
    if( length(D) != nprim ){
      # Just use first 
      D <- rep(D[1],nprim)
    } 
    
    list(s.eta=s, gp.eta=gp, gpp.eta=gpp, D.eta=D, g0.eta=g0, sigma.eta=sigma)
  }
  parms <- f.real.model(beta,nprimary)
  
    print(parms)
#   print(nsecondary)

  # Compute SECR likelihood for each occasion ========================
  hab.pixel.area <- prod(hab.mask@grid@cellsize)
    
  secrLL <- function(i,c.hist,b,ns,trps,aclocs,pix.area){

    ch1 <- c.hist[,i,1:ns[i]]  # histories from ith primary
    caught.this.primary <- rowSums(ch1>0)>0  # animals caught this primary
    ch1 <- ch1[caught.this.primary,]    # remove all 0 lines here
    aclocs <- aclocs[caught.this.primary,i,]  # ac locations of animals caught this primary
    
    # Note, one cannot compute and return pdots here for later use in computing 
    # pstar because uncaptured animals are dropped here. We need pstars for even uncaptured animals.
    # Hence, the call to F.spatial.pstar later.
    
    D.i.pos <- grep("^D",names(b))[i]
    g0.i.pos <- grep("^g0",names(b))[i]
    sigma.i.pos <- grep("^sigma",names(b))[i]
    trps <- trps[[i]]
    
    F.spat.loglik.X(b[c(D.i.pos, g0.i.pos, sigma.i.pos)],ch1,trps,aclocs,pix.area)
  }
  
  closedLL <- sapply(1:nprimary,secrLL,
    c.hist=ch, 
    b=c(parms$D.eta, parms$g0.eta, parms$sigma.eta),
    ns=nsecondary, 
    trps=traps,
    aclocs=ac.locs, 
    pix.area=hab.pixel.area
  )  
  closedLL <- sum(closedLL)

   cat(paste("SECR part:", closedLL, "\n"))
  
  
  
  # Compute Open part of robust design likelihood =========================
  # Take links
  s <- 1/(1+exp(-parms$s.eta))
  gamma <- 1/(1+exp(-parms$gamma.eta))
  g0 <- 1/(1+exp(-parms$g0.eta))
  sigma <- exp(parms$sigma.eta)   # Note log link here, rather than logit
  
  p.star <- F.spatial.pstar(g0, sigma, traps, ch)
  #print(p.star)

  # This returns the "real" log likelihood, not the negative
  openLL <- F.robust.open.part(ch,p.star,s,gamma)
  
  
  # Done ======================================================
  ll <- openLL + closedLL

#   cat(c(-ll))
#   cat(", ")

  cat(c(-ll,beta))
  cat("\n")

  -ll
}

