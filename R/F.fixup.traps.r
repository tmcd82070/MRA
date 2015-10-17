#'  @title Fix up trap locations
#'  
#'  @description Check and fix up the traps location object for spatial robust designs.  
#'  
#'  @details If \code{trps} is a SpatialPoints object, replicate the locations in trps the
#'    appropriate number of times to produce the list of lists needed elsewhere in the estimation 
#'    process. 
#'    
#'  @param traps Either a SpatialPoints object, or list of lists containing trap location matricies.
#'  
#'   @param nprimary number of primary sessions
#'   
#'  @param nsecondary vector of length nprimary containing number of secondary sessions in each primary.
#'  
#'  @return A list of lists containing trap locations for each secondary. 
#'    Each element in the first layer of lists is a sub-lists of "on" trap locations 
#'    for each primary occasion.  Each element in the second layer of lists is a matrix 
#'    of trap location coordinates for "on" traps during each secondary occasion of the  
#'    primary. 
F.fixup.traps <- function( traps, nprimary, nsecondary ){

    
    # Replicate traps object if necessary ---------------
    if( inherits(traps, "SpatialPoints") ){
      # All traps on all the time
      trps <- vector("list",nprimary)
      for(j in 1:nprimary){
        trps[[j]] <- rep(list(coordinates(traps)), nsecondary[j])
      }
      traps <- trps
    } else { 
      # Check length of traps list
      for( j in 1:nprimary){
        if(length(traps[[j]]) != nsecondary[j]) {
          stop("Length of trap list for primary occasion ", j, 
               " does not equal number of secondary occasions (", nsecondary[j], " secondaries found)")
        }
      }
    }
  
    traps

}