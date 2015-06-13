F.JS.estim<-function( ch ){
#
#	Compute Cormack-Jolly-Seber open population capture-recapture
#	estimates.  This follows Jolly (1965).
#
#	Input:
#   ch = capture history matrix (nan X ns).  1 = capture and release, 2 = capture and not released.
#         0 = not captured.
#
# Output: 
# Jolly-Seber estimates, no standard errors for now. 
#
# Note, this routine relies on the m.array, so F.m.array must be available. 
#

  dat.list <- F.m.array(ch)
  
	tt <- dat.list$ns
	ii <- 1:(tt-1)
	ii <- c( t( cbind( ii, rev(ii) )))
	ii <- matrix(rep( rep(c(0,1), length(ii)/2), ii), byrow=T, ncol=tt)
#	print("ii")
#	print(ii)
	
	marr <- dat.list$m * ii    # mask out lower triangle

	r <- c( apply( marr, 1, sum, na.rm=T), NA)
	m <- apply( marr, 2, sum, na.rm=T)
	s <- dat.list$R
	n <- dat.list$n

	z <- 0
	for(i in (2:(tt-1))){
		j <- rep( c(0,1), c(i, (tt-i)) )
		j <- matrix( rep( j, i-1 ), byrow=T, ncol=tt)
		j <- rbind( j, (diag(tt)*0)[1:(tt-i),] )
#		print(i)
#		print(j)

		z <- c(z, sum( marr*j, na.rm=T ))
	}
	z <- c(z, NA)
#	print(r)
#	print(m)
#	print(z)
#	print(s)

	mhat <- (s*z)/r + m
	alpha <- m/n

#	print("M hat")
#	print(mhat)
#	print("alpha hat")
#	print(alpha)

	n.hat <- mhat/alpha

	surv.hat <- c(mhat[2:tt],NA) / (mhat - m + s)

	b.hat <- c(n.hat[2:tt],NA) - surv.hat*(n.hat - n + s)

	p.hat <- n/n.hat

	return(list(n.hat=n.hat, s.hat=surv.hat, p.hat=p.hat, b.hat=b.hat))
}

