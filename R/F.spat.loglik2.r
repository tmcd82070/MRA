# F.spat.loglik -----

F.spat.loglik2 <- function( beta, ch, traps, buffer, type="multi" ){
  # 
  # Compute closed SECR  likelihood with spatial locations using SECR routine
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
  # The spatial Huggins log likelihood, 
  
  require(secr)
  
  # Disallow 0 capture histories
  if( any(rowSums(ch>0) == 0) ) stop(paste("Cannot have capture histories of all zeros. check individual", 
                                           which(rowSums(ch>0)==0)[1]))
  traps(ch) <- traps
  logL <- secr.fit(ch, model = g0~1, start=beta, buffer = buffer, details=list(LLonly=T))

  print(c(beta,-logL))
  
  -logL
  
}



# ============================================================================
# Function calls 

#library(secr)
#setwd(system.file('extdata', package='secr'))
#myCH <- read.capthist('capt.txt','trap.txt', fmt = 'XY')
setwd("~/Programs/MRA/TestingVersion")

#secr0 <- secr.fit(myCH, model = g0~1, buffer = 100, trace = FALSE) 


# Compute my LL for one observation ----------
# tmp <- F.spat.loglik( secr0$fit$par, ch01, attr(ch01,"traps"), 100 )
# cat(paste("My LL for observation", obsn, ":\n"))
# print(tmp)


# Prep and call secr likelihood
# tmp <- F.spat.loglik( secr0$fit$par, myCH, attr(myCH,"traps"), 100 )
# cat("Success...\n")
# print(tmp)

#fit1 <- nlminb(secr0$fit$par+ rnorm(3,0,.1), F.spat.loglik, ch=myCH, traps=attr(myCH,"traps"), buffer=100)
fit2 <- optim(secr0$fit$par+ rnorm(3,0,.1), F.spat.loglik2, ch=myCH, traps=attr(myCH,"traps"), buffer=100)
print(fit2)
# 
# tmp2 <- secr.fit(myCH, model = g0~1, start=fit1$par, buffer = 100, details=list(LLonly=T))  
# cat("SECR LL at my params:\n")
# print(tmp2)