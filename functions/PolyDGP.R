#' DGP for Exp 5 
#' See e-mail from 9/3/2018 by Alberto
#' 
#' CREATED: 10/3/2018
#' 
#' @param n1 number of treated
#' @param n0 number of non treated
#' @param p is the number of covariates
#' @param delta is the order of the polynomial (should be an integer)
#' 
#' @author Jeremy LHour

PolyDGP <- function(n1=25,n0=50,p=3,delta=2){
  d = c(rep(1,n1),rep(0,n0))
  X = rbind(matrix(runif(n1*p), ncol=p, nrow=n1),
            matrix(rnorm(n0*p), ncol=p, nrow=n0))
  # Treated ~ uniform density [0;1]
  # Controls ~ N(0,1) density
  # all iid (across i and j)
  
  ### Create outcome function
  stdev = sqrt(p/(2*delta+1) - p/(1+delta)^2)
  beta = rep(1,p)/stdev
  
  y = (X^delta)%*%beta + rnorm(n1+n0)
  
  return(list(X=X,y=y,d=d))
}