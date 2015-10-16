#'  @title Estimate parameters given activity center locations
#'  
#'  @description Given a set of activity center locations, maximize the spatial robust likelihood 
#'  to estimate parameters.  This is the "M" step of the EM algorithm. 
#'  
#'
#' @param ac.locs A nan X nprimary X 2 array of estimated activity center locations.  ac.locs[i,j,] = c(x,y) = coordinates of 
#' activity center for animal i during primary occasion j.
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
#  
#' @param beta.init A named list of initial coefficient estimates.  Length of \code{beta.init} is either 5 or 6. If  \code{length(beta.init)}
#' equals 5, assume that the 'random emigration' model is being fitted wherein emigration probability (gp) and immigration 
#' probability (gpp) are equal. If  \code{length(beta.init)} equals 6, assume 'Markov emigration' model is being fitted 
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
#' 
#' Note this structure for beta facilitates pure time varying models on all parameters, and subsets (e.g., constant).  This 
#' structure will not facilitate covariates or individual heterogeneity.  I don't know what you'll do if you need those. 
#' 
#' @return A named list of coefficient estimates, structured exactly like \code{beta.init}. 
#' 
#' @details Combining this routine with \code{F.AC.estim} one can implement the EM algorithm for estimating a spatial robust design. 
#' 
F.parm.estim <- function(ac.locs, ch, traps, beta.init){

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
  
  
  # Fix up beta.  I.e., take it from a list to a vector so optim works
  b.init <- unlist(beta.init)  # this concatinates names and element numbers. E.g., second element of "surv" vector gets named "surv2"
  
  # Do the maximization
  fit <- optim( b.init, F.spat.robust.loglik.X, ch=ch, ac.locs=ac.locs, traps=traps, method="L-BFGS-B", lower=llimit, upper=hlimit, 
                hessian = TRUE, control=list(factr=5e9, pgtol=1e-8, maxit=1000))
  
  
  fit$par
}  