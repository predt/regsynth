#' Data Generating process for Monte Carlo study
#' DGP in paper
#' 
#' CREATED: 21/12/2017
#' 
#' @param n1 number of treated
#' @param n0 number of untreated
#' @param p is the number of covariates
#' @param a lower bound of support for treated
#' @param b upper bound of support for treated
#' @param h discrepency between treated and untreated support. Should be positive.
#' 
#' @author Jeremy LHour

UniformDGP <- function(n1=10,n0=100,p=1,a=.1,b=.9,h=.1){
  library("MASS")
  
  d = c(rep(1,n1),rep(0,n0))
  
  ### Draw covariates
  X = matrix(nrow=n1+n0,ncol=p)
  X[d==1,] = matrix(runif(p*sum(d),min=a,max=b), ncol=p)
  X[d==0,] = matrix(sqrt(runif(p*sum(1-d),min=a-h,max=b+h)), ncol=p)

  return(list(X=X,
              d=d))
}