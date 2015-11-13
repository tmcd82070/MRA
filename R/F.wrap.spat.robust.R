######################################################
## Jared Studyvin
## 11 Nov 2015
## Wrapper function for F.spat.robust.loglik
######################################################

#' inputs:
#' @param captureHistory = array animal X primary X secondary
#' @param maskSP spatial points data frame
#' @param trapSP at least a spatial points object. I'm not sure why data would be needed.
#' @param trap.usage array of 1, indicate trap is on, or 0, indicate trap is off. The array dim should be number of traps by number or primaries by number of secondaries.
#' @param D.model Default is intercept only. Can take covariates from the data frame in maskSP
#' @param g0.model Default is intercept only. Currently ignored.
#' @param sigma.model Default is intercept only. Currently ignored.
#' @param survival.model Default is intercept only.
#' @param gammap.model Default is intercept only. Currently ignored.
#' @param gammapp.model Default is setting gammapp = gammap. Currently ignored.
#'
#'
#' @param animalCovariate dataframe. Currently ignored
#'
#'
#'





code.dir <- '~/GoogleDrive/MRA/R/'
lapply(paste0(code.dir,c('F.collapse.secondaries.r','F.robust.open.part.r','F.spat.loglik.r','F.robust.loglik.r')),source)

## area from secr
0.007488869


formula <- list(D=~strata,g0=~1,sigma=~1)

D.model=~strata


F.wrap.spat.robust <- function(captureHistory,maskSP,trapSP,trap.usage,D.model=~1,g0.model=~1,sigma.model=~1,survival.model=~1,gammap.model=~1,gammapp.model=NULL,animalCovariate=NULL){

    ##needed packages
    require(sp) ##


    ## check the depend functions are loaded
    ## this may not be need for the finished package
    neededFunctions <- c('F.collapse.secondaries','F.robust.open.part','F.spat.loglik','F.robust.loglik')
    checkNeedFunction <- neededFunctions%in%ls()
    if(sum(checkNeedFunction)<length(neededFunctions)){
        stop(paste('The following functions are not found:\n',
                   paste(neededFunctions[!checkNeedFunction],collapse=' ')))
    }


    if( !inherits(captureHistory,"array") ) stop("ch must be an array")
    if( length(dim(captureHistory)) != 3) stop(paste("ch must have 3 dimensions.", length(dim(ch)), "found."))


#########################################################
### Other checks of the data need to go here


    ch <- captureHistory
    mask <- maskSP
    trap <- coordinates(trapSP)
    D.mod <- D.model
    g0.mod <- g0.model #ignored
    sig.mod <- sigma.model #ignored
    surv.mod <- survival.model
    gamp.mod <- gammap.model #ignored
    gampp <- gammapp.model #ignored
    animalCov <- animalCovariate #ignored

#########################################################





    ## Find dimensions
    d <- dim(ch)
    nan <- d[1]
    nprimary <- d[2]
    nsecondary <- apply(ch,2,function(x){
                            xx<-apply(x,2,function(y){all(is.na(y))})
                            sum(!xx)
                        })


    ## this the area (if utm it should be square meters) assocaited with one mask point
    maskArea <- prod(gridparameters(as(mask,'SpatialPixels'))[,'cellsize'])
    maskArea <-  maskArea/10000 ## convert to hectares. This matches secr, but I think it should be removed



    D.model.matrix <- model.matrix(D.mod,data=mask@data)
    B.beta <- rep(NA,ncol(D.model.matrix))
    head(maskSP@data$strata)



  f.real.model<-function(beta, nprim){
    # Normally, you would have a X matrix here.
    # For now, use names of beta to evaluate model.
    #
    # Current implementation allows small t or dot model for  s
    # This returns linear predictors.  Link functions are applied later.

    s.parms <- grep("^S",names(beta))
    g.parms <- grep("^gamma",names(beta))


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

    # g parameters
    # This is the random emigration model, constant across time
    g <- beta[g.parms]
    if( length(g) == (nprim-2)){
      # constrain last g to second to last g.
      g <- c(NA,g,g[nprim-2])
    } else if(length(g) == (nprim-1)){
      g <- c(NA,g)  # use all gammas
    } else if( length(g) == 0 ){
      # Constrain gamma to zero
      g <- c(NA,rep(0,nprim-1))
    } else if( length(g) != (nprim-1)){
      g <- c(NA,rep(g[1],(nprim-1))) # first gamma not possible.
    }
    # For rest of routine to work g must have length nprim.  g[1]=NA

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

    list(s.eta=s, gamma.eta=g, D.eta=D, g0.eta=g0, sigma.eta=sigma)
  }
  parms <- f.real.model(beta,nprimary)




# find MLE
F.spat.robust.loglik




} # end function
