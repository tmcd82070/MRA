#' @title Estimate Activity Center locations. 
#' 
#' @description Given estimates of beta, compute estimates of latent activity center locations.
#' 
#' @param beta A named list of coefficient estimates.  Length of \code{beta} must be either 5 or 6. If  \code{length(beta)}
#' equals 5, assume that the 'random emigration' model is being fitted wherein emigration probability (gp) and immigration 
#' probability (gpp) are equal. If  \code{length(beta)} equals 6, assume 'Markov emigration' model is being fitted 
#' wherein emigration (gp) and immigration (gpp) depend upon whether an individual is on or off the 
#' study area (see robust design literature for more details). There are 2 or 3 elements for open part of likelihood
#' (Suvival, immigration, emigration) and 3 for closed part of model (density, g0 in sightability, sigma in sightability).   
#' Names and sizes of elements is as follows:  
#' . "surv" (survival) = a vector of length (nprimary-1) corresponding to intervals between primary occasions
#' . "gp" (gamma prime, emigration) = a vector of length (nprimary-1) corresponding to intervals between primary occasions.
#' . (Optional) "gpp" (gamma prime prime, immigration) = a vector of length (nprimary - 2) corresponding to intervals between 
#'    2nd primary and the last.  Recall there is no immigration parameter for the interval between primary 1 and 2. 
#' . "D" (density) = a vector of length nprimary. 
#' . "g0" (capture probability at traps) = a vector of length nprimary
#' . "sigma" (capture probability width) = vector or length nprimary
#' Here, the only the sightability coefficients (i.e., g0 and sigma) are used.
#' 
#' Note this structure for beta facilitates pure time varying models on all parameters, and subsets (e.g., constant).  This 
#' structure will not facilitate covariates or individual heterogeneity.  I don't know what you'll do if you need those. 
#' 
#' @param ch A 3-D array of capture histories.  Size is nan X nprim X nsecondaries.  Missing secondaries 
#' are specified with columns of all NA's (i.e., if the 4th secondary of the ith primary was not done, 
#' \code{all(is.na(ch[,i,4])) == TRUE}).  Cells in the array contain trap number that caught the individual.
#' Trap numbers are rows in the \code{traps} matrix.
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
#'  @param hab.mask A SpatialPixels or SpatialGrid object (or their *DataFrame analogs) containing allowable 
#'    locations for activity centers. 
#'    Exactly like masks in the SECR package, \code{hab.mask} defines the outer limit of integration, defines 
#'    sites that are habitat and thus can be occupied, descritizes habitat covariates for use in models.
#'    The current implementation does not allow the habitat mask to change between primaries or secondaries. 
#'    The mask is constant throughout the study.  Here, \code{hab.mask} is used to test whether 
#'    an estimated activity center is in valid habitat, and if not, moved to a location in valid habitat. 
#'    Because centers of the pixels or cells are used for this purpose, size (area) of the mask is not 
#'    needed, but the mask still must be projected in the same units as \code{traps}. 
#'
#' @details  This function estimates activity center locations given trapped locations using maximum likelihood. 
#' It finds the (X,y) location which maximized the probability of obsevering the capture history. 
#' 
#' @return A nan X nprimary X 2 array of estimated activity center locations.  return[i,j,] = c(x,y) = coordinates of 
#' activity center for animal i during primary occasion j.
#'  
#

F.AC.estim <- function( beta, ch, traps, hab.mask){
  
  if( !inherits(ch,"array") ) stop("ch must be an array")
  if( length(dim(ch)) != 3) stop(paste("ch must have 3 dimensions.", length(dim(ch)), "found."))
  
  
  nan <- dim(ch)[1]
  nprimary <- dim(ch)[2]
  nsecondary <- apply(ch,2,function(x){
    xx<-apply(x,2,function(y){all(is.na(y))})
    sum(!xx)
  })  
  
  # Internal functions -----------------------------------------
  f.get.trap.locs <- function(tr, h){
    # return trap locations where caught
    ans <- NULL
    for( i in 1:length(tr)){
      if( h[i] != 0 ){
        ans <- rbind( ans, coordinates(tr[[i]])[h[i],])
      }
    }
    ans
  }
  
  # Main code -----------------------------------------

  # Fix up beta vector --------------------------------
  # Extract coefficients from beta vector. Note, only time varying models here. No covariates or individual variation.
  g0 <- beta[["g0"]]  # vector length nprimary
  sigma <- beta[["sigma"]]  # vector length nprimary
  
    
  # Replicate traps object if necessary ---------------
  if( "SpatialPoints" %in% class(traps)){
    # All traps on all the time
    trps <- vector("list",nprimary)
    for(j in 1:nprimary){
      trps[[j]] <- rep(list(traps), nsecondary[j])
    }
    traps<-trps
  } else { 
    # Check length of traps list
    for( j in 1:nprimary){
      if(length(traps[[j]]) != nsecondary[j]) {
        stop("Length of trap list for primary occasion ", j, 
             " does not equal number of secondary occasions (", nsecondary[j], " secondaries found)")
      }
    }
  }
  
#   print(length(traps))
#   print(length(traps[[1]]))
    
  # Loop over animals first, then primaries ------------
  ac.locs <- array( NA, c(nan,nprimary,2) )
  for( i in 1:nan){
    for( j in 1:nprimary){
      # History for animal i during primary j
      h <- ch[i,j,]
      
      # Get list of traps "on" during secondaries of primary j
      trapsj <- traps[[j]]

      # debugging
#       cat(paste("------- animial", i, ", primary", j, "-------\n" ))
#       print(h)
      
      # Initial AC location
      if( !all(h==0) ){
        x.init <- f.get.trap.locs(trapsj, h)  # returns matrix of X,Y coordinates where captured
        x.init <- colMeans( x.init )
      } else {
        # Unknown what is best to do here.  Need a point way off grid.  Start at 
        # and extreme trap. Could pick random distance and direction from centroid.  
        # Do not start at centroid! optim bombs in some cases when you start at the centroid.
        x.init <- lapply( trapsj, function(x){apply(coordinates(x),2,max,na.rm=TRUE)} )
        x.init <- rowMeans(matrix(unlist(x.init),2))
      }

      # Debugging
#       cat("x.init=")
#       cat(x.init)
#       cat("\n")
#       
#       # Debugging
#       cat("b=")
#       cat(c(g0,sigma))
#       cat("\n")
      
      
      
      # Maximize the likelihood of h over all X,Y locations 
      ac.fit <-  optim( x.init, F.X.loglik, h=h, traps=trapsj, g0=g0[j], sigma=sigma[j])
      
      # print(ac.fit)
      
      if( ac.fit$convergence == 0){
        if( !all(h==0) ){
          # Animal was seen.  ac.fit may not be in habitat. Check and move it to habitat point with highest loglik. 
          ac.locs[i,j,] <- F.move.2.habitat(ac.fit$par, hab.mask, h, trapsj, g0[j], sigma[j])
        } else {
          # un seen animals AC is not required to be in habitat
          ac.locs[i,j,] <- ac.fit$par
        }
      } else {
        stop("No AC location for animal ",i,"during primary occasion",j,"\nConvergence code=",ac.fit$convergence)
      }

#       tmp <- readline()
#       if(nchar(tmp)>0) stop("okay, stop.")
      
    }
  }
  
  ac.locs
  
}

# # Testing ----------------------
# require(sp)
# 
# source("C:/Users/tmcdonald/Google Drive/Documents/Programs/MRA/R/F.X.loglik.r")
# 
# b <- list(surv = rep(0.9,3), 
#           gp = rep(.1,3),  
#           D = rep(11.724602,4),  
#           g0= rep(0.5974225,4), 
#           sigma=rep(12.764144,4)
#           )
# 
# tmp <- F.AC.estim(b, H, X)
# 
# 
# for( i in 1:nrow(H)){
#   plot(sa, main=paste("Animal",i),)
#   plot(X,  add=T)
#   for( j in 1:ncol(H)){
#     if( !all(H[i,j,]==0) ){
#       points(X[H[i,j,],], col=j, pch=16, cex=1.5)
#       points(tmp[i,j,1], tmp[i,j,2], pch=15, col=j, cex=1.5)
#     }
#   }
#   cat("Continue? (any character quits):")
#   tmp.2 <- readline()
#   if(nchar(tmp.2)>0) stop("okay, stop.")
# }