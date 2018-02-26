#' Data Generating process for Monte Carlo study
#' Alberto DGP
#' 
#' CREATED: 27 juillet 2017
#' 
#' @param delta is the order of the legendre polynomial (should be an integer)
#' @param n1 number of treated
#' @param n0 number of non treated
#' @param sigma is the standard deviation of the error terms
#' 
#' @autor Jeremy LHour

SimpleAlbertoDGP <- function(delta=1,n1=25,n0=50,sigma=1){
  d = c(rep(1,n1),rep(0,n0))
  X = c(runif(n1,min=-1,max=1),2*sqrt(runif(n0))-1)
  # Treated ~ uniform density [-1;1]
  # Controls ~ triangular density [-1;1]

  ### Create outcome function
  poly = legendre.polynomials(delta)
  polypart = vector(length=n1+n0)
  for(k in 1:(delta+1)){
    polypart = polypart + unlist(polynomial.values(poly[k],X))
  }
  polypart = polypart / sd(polypart) # normalize the part in X
  
  y = polypart + rnorm(n1+n0,sd=sigma)
  
  return(list(X=X,
              y=y,
              d=d))
}