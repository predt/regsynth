#' Data Generating process for the example
#' 
#' @param n number of observations
#' 
#' 
#' 

matchDGP <- function(n=20,p=2,Ry=.5,Rd=.2,rho=.5){
  library("MASS")
  
  ### Covariate variance matrix
  Sigma <- matrix(0,nrow=p, ncol=p)
  
  for(k in 1:p){
    for(j in 1:p){
      Sigma[k,j] <- rho^abs(k-j)
    }
  }
  
  ### Treatment variable coefficient
  gamma <- rep(0,p)
  
  for(j in 1:p){
    gamma[j] <- 1*(-1)^(j) / j^2
  }
  
  ### Outcome equation coefficients
  b <- gamma
  
  for(j in 1:p){
    b[j] <- (-1)^(j+1) / (p-j+1)^2
  }
  
  ### Adjustment to match R.squared
  c <- sqrt((1/t(gamma)%*%Sigma%*%gamma)*(Rd/(1-Rd)))
  gamma <- c*gamma
  
  c <- sqrt((1/t(b)%*%Sigma%*%b)*(Ry/(1-Ry)))
  b <- c*b
  
  X <- mvrnorm(n = n, mu=rep(0,p), Sigma)
  d <- as.numeric(runif(n) < pnorm(X%*%gamma))
  
  ### Treatment effect
  a <- 0
  
  y <- a*d + exp(-X%*%b) + rnorm(n)
  
  return(list(X=X,
              y=y,
              d=d,
              b=b,
              g=gamma,
              ATT=sum(d*a)/sum(d)))
  
  
  
  n <- 20
  x <- matrix(rexp(2*n), ncol=2)
  g <- c(.5,-1)
  logit <- function(x) 1/(1+exp(-x))
  d <- as.numeric(runif(n) < logit(x %*%g))
  y <- 10*exp(-diag(x%*%t(x)))
  data <- data.frame(X=x,
                     d=d)
}