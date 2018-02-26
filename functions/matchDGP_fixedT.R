#' Data Generating Process with fixed treated/control
#' 
#' @param n1 number of treated observations
#' @param n0 number of control observations
#' @param p is the number of covariates
#' @param Ry gives the R-square wanted in outcome equation
#' @param Rd gives the R-square wanted in treatment equation
#' @param rho parametrizes the covariance between the covariates
#' @param a is the treatment effect value
#' 
#' @autor Jeremy LHour

matchDGP_fixedT <- function(n1=20,n0=20,p=2,Ry=.5,Rd=.2,rho=.5,a=0){
  
  ### Covariate variance matrix
  Sigma = matrix(0,nrow=p, ncol=p)
  for(k in 1:p){
    for(j in 1:p){
      Sigma[k,j] = rho^abs(k-j)
    }
  }
  
  ### Outcome equation coefficients
  b = rep(0,p)
  for(j in 1:p){
    b[j] = (-1)^(j+1) / (p-j+1)^2
  }
  c = sqrt((1/t(b)%*%Sigma%*%b)*(Ry/(1-Ry)))
  b = c*b
  
  ### Mean of X for treated
  gamma = rep(0,p)
  for(j in 1:p){
    gamma[j] = 1*(-1)^(j) / j^2
  }
  c = sqrt((1/t(gamma)%*%Sigma%*%gamma)*(Rd/(1-Rd)))
  gamma = c*gamma
  
  d = c(rep(1,n1),rep(0,n0))
  
  X = mvrnorm(n = n1+n0, mu=rep(0,p), Sigma) + (d%o%gamma)
  y = a*d + exp(-X%*%b) + rnorm(n1+n0)
  
  return(list(X=X,
              y=y,
              d=d))
}