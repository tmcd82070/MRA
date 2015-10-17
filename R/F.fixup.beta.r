#'  @title Fix up coefficient vector
#'  
#'  @description Apply the link or inverse link to beta coefficients. 
#'  
#'  @param beta A named list of coefficients.
#'  
#'  @param to.real If TRUE, assume \code{beta} is on the linear scale, and 
#'    transform output will be transformed to the response scale using the inverse-link.  Inverse-link
#'    scale is scale of the "real" parameters (probabilities, variances, etc.).  If 
#'    FALSE, transform beta from response scale to linear scale. 
#'    
#'  @return A named list, exactly like \code{beta}, except with parameters transformed 
#'    the correct direction.
#'     
F.fixup.beta <- function( b, to.real=TRUE ){
  if( !to.real ){
    b$surv <- log(b$surv/(1-b$surv))
    b$gp <- log(b$gp/(1-b$gp))
    b$g0 <- log(b$g0/(1-b$g0))
    b$D <- log(b$D)
    b$sigma <- log(b$sigma)
    if("gdp" %in% names(b)){
      b$gdp <- log(b$gdp/(1-b$gdp))
    }
  } else {
    b$surv <- 1 / (1+exp(-b$surv))
    b$gp <- 1 / (1+exp(-b$gp))
    b$g0 <- 1 / (1+exp(-b$g0))
    b$D <- exp(b$D)
    b$sigma <- exp(b$sigma)
    if("gdp" %in% names(b)){
      b$gdp <- 1 / (1+exp(-b$gdp))
    }
  }
  
  b
}