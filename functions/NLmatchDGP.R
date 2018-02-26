#' Data Generating process for Monte Carlo study
#' Possible to module the non-linearity
#' 
#' CREATED: 19 avril 2017
#' EDITED: 24 juillet 2017
#' 
#' @param n number of observations wanted
#' @param p is the number of covariates
#' @param delta is the order of the legendre polynomial (should be an integer)
#' @param Ry gives the R-square wanted in outcome equation
#' @param Rd gives the R-square wanted in treatment equation
#' @param rho parametrizes the covariance between the covariates
#' @param a is the treatment effect value
#' 
#' @autor Jeremy LHour

NLmatchDGP <- function(n=20,p=2,delta=1,Ry=.5,Rd=.2,rho=.5,a=0){
  library("MASS")
  library("orthopolynom")
  
  ### Covariate variance matrix
  Sigma = matrix(0,nrow=p, ncol=p)
  for(k in 1:p){
    for(j in 1:p){
      Sigma[k,j] = rho^abs(k-j)
    }
  }
  
  ### Treatment variable coefficient
  gamma = rep(0,p)
  for(j in 1:p){
    gamma[j] = 1*(-1)^(j) / j^2
  }
  
  ### Outcome equation coefficients
  b = gamma
  for(j in 1:p){
    b[j] = (-1)^(j+1) / (p-j+1)^2
  }
  
  ### Adjustment to match R.squared
  c = sqrt((1/t(gamma)%*%Sigma%*%gamma)*(Rd/(1-Rd)))
  gamma = c*gamma
  
  c = sqrt((1/t(b)%*%Sigma%*%b)*(Ry/(1-Ry)))
  b = c*b
  
  X = mvrnorm(n = n, mu=rep(0,p), Sigma)
  d = as.numeric(runif(n) < pnorm(X%*%gamma))
  
  ### Create outcome function
  poly = legendre.polynomials(delta)
  polypart = unlist(polynomial.values(poly[delta+1],X%*%b))
  
  y = ifelse(abs(X%*%b)<1,(X%*%b+1)/2,polypart) + a*d + rnorm(n)
  
  return(list(X=X,
              y=y,
              d=d,
              b=b,
              g=gamma))
}