F.pt.searchdist<-function(sa){
  
  if( length(grep("SpatialGrid", class(sa))) > 0 ){
    ans <- getGridTopology(sa)
    ans <- ans@cellsize 
    ans <- sqrt( sum(ans^2))
  } else if( length(grep("SpatialPoints", class(sa))) > 0 ){
    ans <- gDistance(sa,byid=TRUE)
    ans <- apply(ans, 1, function(x){sort(unique(x))[4]})
    ans <- mean(ans) + 2*sd(ans)  # mean 2nd nearest neighbor distance + 2 SE. Assuming normal, this is 95% of them
  }
  
  ans
}