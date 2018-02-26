#' Data Generating process for the example
#' 
#' @param n number of observations wanted
#' @param p is the number of covariates
#' @param Ry gives the R-square wanted in outcome equation
#' @param Rd gives the R-square wanted in treatment equation
#' @param rho parametrizes the covariance between the covariates
#' @param a is the treatment effect value
#' 
#' @autor Jeremy LHour

matchDGP <- function(n=20,p=2,Ry=.5,Rd=.2,rho=.5,a=0){
  
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
  
  y = a*d + exp(-X%*%b) + rnorm(n)
  
  return(list(X=X,
              y=y,
              d=d,
              b=b,
              g=gamma,
              ATT=sum(d*a)/sum(d)))
}