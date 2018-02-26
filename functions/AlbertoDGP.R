#' Data Generating process for Monte Carlo study
#' Alberto DGP, multivariate
#' 
#' CREATED: 27 juillet 2017
#' EDITED: 10 aout 2017
#' 
#' @param n number of observations
#' @param p is the number of covariates
#' @param delta is the order of the legendre polynomial (should be an integer)
#' @param prop gives the expected proportion of treated
#' @param rho parametrizes the covariance between the covariates
#' @param a is the treatment effect value
#' @param sigma is the standard deviation of the error terms
#' 
#' @autor Jeremy LHour

AlbertoDGP <- function(n=20,p=1,delta=1,prop=.1,rho=.5,a=0,sigma=1){
  library("MASS")
  library("orthopolynom")
  
  ### Covariate variance matrix
  SigmaMat = matrix(0,nrow=p, ncol=p)
  for(k in 1:p){
    for(j in 1:p){
      SigmaMat[k,j] = rho^abs(k-j)
    }
  }
  
  L = (3/sqrt(2))*t(chol(SigmaMat)) # Cholesky transform, scaled
  
  ### Outcome equation coefficients
  b = rep(0,p)
  for(j in 1:p){
    b[j] = (-1)^(j+1) / (p-j+1)^2
  }
  
  ### Draw labels
  uu=0
  while(uu<2){
    d = runif(n) < prop
    uu = sum(d)
  }
  
  ### Draw covariates
  X = matrix(nrow=n,ncol=p)
  X[d==0,] = matrix(2*sqrt(runif(p*sum(1-d)))-1, ncol=p) # Distributed as triangular density on[-1;1] for untreated
  X[d==1,] = matrix(runif(p*sum(d),min=-1,max=1), ncol=p) # Distributed as uniform density on [-1;1] for untreated
  X = X %*% L
  
  ### Create outcome function
  poly = legendre.polynomials(delta)
  polypart = unlist(polynomial.values(poly[delta+1],X%*%b))
  polypart = polypart / sd(polypart) # normalize the part in X
  
  y = polypart + a*d + rnorm(n,sd=sigma)
  
  
  return(list(X=X,
              y=y,
              d=d,
              b=b))
}