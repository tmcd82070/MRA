#########################################
## Jared Studyvin
## 29  Sept 2015
## Test SECR and Trent's likelihood for SECR
#########################################
rm(list=ls())
#library(mra) # Trent's package


library(secr)
setwd(system.file('extdata', package='secr'))
myCH <- read.capthist('capt.txt','trap.txt', fmt = 'XY')
setwd("~/GoogleDrive/MRA/R")


## create the habitat mask
mask <- make.mask(attr(myCH,'traps'),buffer=100,nx=50,ny=50,spacing=10)

plot(mask)
points(attr(myCH,'traps'))

## fit the secr model
secr0 <- secr.fit(myCH, model = g0~1,mask=mask, trace = FALSE)

## area associated with each mask location
area <- attr(mask,'area')

## check the n/a gives the density parameter
sigma <- exp(coef(secr0)$beta[3])
g0 <- exp(coef(secr0)$beta[2])/(1+exp(coef(secr0)$beta[2]))
P. <- pdot(mask,attr(myCH,'traps'),detectpar=list(g0=g0,sigma=sigma),noccasions=5)
D <- nrow(myCH)/(sum(P.)*area)
D;log(D)
coef(secr0)





source('F.spat.loglik.r')

CH <- as.data.frame(myCH[1:nrow(myCH),1:ncol(myCH)])
CH

tp <- data.frame(attr(myCH,'traps'))
names(tp)

X <- data.frame(x=mask$x,y=mask$y)

fitfull <- nlminb(secr0$fit$par+ rnorm(3,0,.1), F.spat.loglik, ch=CH, traps=tp,mask=X,area=area)

fitfull$par

secr0$fit$par





