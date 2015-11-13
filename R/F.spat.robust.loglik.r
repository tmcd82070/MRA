# Fits a Spatial Robust design model via maximum likelihood
#
# Objective function -----
F.spat.robust.loglik <- function( beta, ch, traps,hab.mask,pixel.area,subsetp.=1.96){

    require(plyr)
## this is a test of gitHub

  #
  # inputs:
  #   beta = coefficients
  #   ch = 3-d array, rows are individuals, columns are secondary occasions, pages are primary occasions.
  #        columns that are all NA are "not there", meaning variable number of secondary occasions are specified
  #        by missing columns in one of the pages (e.g., ch[,3,2] could be all NA, but ch[,3,1] is there.
  #        this means 3 secondary occasions during primary 1 but only 2 during primary 2).
  #        Value in each cell is the number (row in trap location matrix) of trap that caught the individual.
  #        "0" means un-captured.
  #   traps = K x 2 matrix of trap locations
  #   buffer = distance to buffer trap locations
  #
  #

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

  # Evaluate the model to get real parameters
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

    print(parms)
#   print(nsecondary)

  ## Compute SECR likelihood for each occasion ========================
  ## For now, SECR parameters are constant accross (primary) sessions

  getClosedLL <- function(i,c.hist,b,ns,trps,mask,area){
      ch1 <- c.hist[,i,1:ns[i]]  # remove NA's here
      ch1 <- ch1[rowSums(ch1>0)>0,]    # remove all 0 lines here

      D.i.pos <- grep("^D",names(b))[i]
      g0.i.pos <- grep("^g0",names(b))[i]
      sigma.i.pos <- grep("^sigma",names(b))[i]

      F.spat.loglik(beta=b[c(D.i.pos, g0.i.pos, sigma.i.pos)],ch=ch1,traps=trps,mask=mask,mask.pixel.area=area,only.loglik=FALSE)
  }

## for testing
##  getClosedLL(1,c.hist=ch,b=c(parms$D.eta, parms$g0.eta, parms$sigma.eta),ns=nsecondary,trps=traps,mask=hab.mask,area=pixel.area)

  ## closedLLresult <- sapply(1:nprimary,getClosedLL,
  ##                    c.hist=ch,
  ##                    b=c(parms$D.eta, parms$g0.eta, parms$sigma.eta),
  ##                    ns=nsecondary,
  ##                    trps=traps,
  ##                    mask=hab.mask,
  ##                    area=pixel.area
  ##                    )

    closedLLresult <- llply(1:nprimary,getClosedLL,
                            c.hist=ch,
                            b=c(parms$D.eta, parms$g0.eta, parms$sigma.eta),
                            ns=nsecondary,
                            trps=traps,
                            mask=hab.mask,
                            area=pixel.area
                            )
    ## this is the "real" log likelihood, no negative one
    closedLL <- sum(laply(closedLLresult,.fun=function(x){return(x$logL)}))



    ## the pdots from the closed portion of the likelihood
    ## each column of p.prim refers to each of the primaries
    p.prim <- laply(closedLLresult,.fun=function(x){return(x$p.)})
    p.prim <- t(p.prim)


#   cat(paste("SECR part:", closedLL, "\n"))

  # Compute Open part of robust design likelihood =========================
  # Take links
  s <- 1/(1+exp(-parms$s.eta))
  gamma <- 1/(1+exp(-parms$gamma.eta))
  g0 <- 1/(1+exp(-parms$g0.eta))
  sigma <- exp(parms$sigma.eta)   # Note log link here, rather than logit

    ## this creates the subset of pdots to create the pstars
    ## this code assumes there are no covariates for g0 or sigma
    subG <- dnorm(subsetp.*sigma,sd=sigma)
    subH <- -log(1-g0*subG)
    p.summary <- rep(NA,nprimary)
    for(i in 1:nprimary){
        hc <- rep(subH[i],nsecondary[i])
        p.sc <- 1-prod(1-hc)
        p.summary[i] <- mean(p.prim[p.prim[,i]>p.sc,i])
    }


## p.star <- F.spatial.pstar(g0, sigma, traps, ch)   # returns nan X nprimary matrix or p.dots


 ## p.star <- F.spatial.pstar2(g0, sigma, traps, ch, hab.mask)

  #print(p.star)
    ##colMeans(p.prim)
    p.star <- p.summary

  # This returns the "real" log likelihood, not the negative
  openLL <- F.robust.open.part(ch,p.star,s,gamma)


  # Done ======================================================
  ll <- openLL + closedLL

#   cat(c(-ll))
#   cat(", ")

#  cat(c(-ll,beta))
#  cat("\n")

  -ll
}

