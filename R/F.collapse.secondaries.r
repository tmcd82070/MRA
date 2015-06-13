F.collapse.secondaries <- function(ch){
  # collapse secondary periods of a Robust design to just primaries, 
  # to captured or not.  result is nan X nprimary
  #
  # ch = a 3D where rows= animals, columns = primaries, pages = secondaries within primaries.

    ch.reduced <- apply(ch,2,function(c.hist){
      as.numeric(rowSums(c.hist>0,na.rm=T)>0)
    })  
    
    ch.reduced

}
