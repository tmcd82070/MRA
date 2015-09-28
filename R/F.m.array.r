F.m.array<-function( hists ){
#
#	compute capture-recapture "m array" from capture histories.
#
#	Input:
#	hists = NaN X Ns matrix of 0's, 1's, and 2's representing capture histories
#		(NaN = number of animals; Ns = number of trap occations)
#	  Here 	0 = no capture, 1 = capture released alive, 
#		2 = capture but died before release
#		NA's are not allowed (use 0)
#
#	Output:
#	a list containing the m array and other info suitable for input to 
#	other routines (like F.JS.estim)
#	components are 
#			ns = number of trap occations
#			m = (t-1)Xt matrix of the n(ij), ie. the "m array"
#				  used by Burnham, White, and others. 
#				  Only the upper off-diagonal triangle is used.
#			R = vector of size t containing the number 
#				  released alive from the i-th sample after marking.
#			n = vector of size t containing number
#				  captured at each occation
  

	n <- apply( (hists>=1), 2, sum )
	s <- apply( (hists==2), 2, sum )
	s <- n - s
	tt <- ncol(hists)

	marr <- matrix( rep(NA,tt*(tt-1)), nrow=tt-1) 
  
	for( i in (1:(tt-1))){
		h <- hists[ hists[,i]>=1, ]
		for(j in ((i+1):tt)){
			k1 <- rep(0,tt)
			k2 <- rep(0,tt)
			k1[c(i,j)]<-1
			k2[i:j] <- 1

			k3 <- apply( (t(h) == k1)*k2, 2, sum)
			marr[i,j] <- sum( k3 == sum(k2) )
		}
	}

	return( list(m=marr, R=s, n=n, ns=tt))
}
			
# ====================================

#ch.reduced <- F.collapse.secondaries(CH)
#tmpm <- F.m.array(ch.reduced)


