#' Permutation Test for regularized SYnthetic Control
#' 
#' Created: 9 septembre 2016
#' Returns R ATET estimates on reshuffled samples
#' and compute the p-value.
#' 
#' @param d is a vector of dimension n (treatment indicator)
#' @param X is a matrix of dimension n x p
#' @param y is a vector of dimension n (outcome)
#' @param V is a p x p matrix of weights
#' @param lambda is a positive penalty level
#' @param R is the number of replications
#' 
#' @autor Jeremy LHour

perm.test <- function(d,y,X,V,lambda,R=1000){
  
  # Compute ATET on original sample
  X0 = t(X[d==0,]); X1 = t(X[d==1,]);
  Y0 = y[d==0]; Y1 = y[d==1];
  sol1 = regsynth(X0,X1,y[d==0],y[d==1],V,lambda)
  theta.hat = sol1$ATT
  
  # Compute ATET on as many reshuffled samples
  theta.reshuffled = replicate(R, permutation.iter(d,y,X,V,lambda), simplify="vector")
  
  # Compute p-value
  p.val = (sum(theta.hat < theta.reshuffled)+1)/(R+1)
  
  return(list(p.val=p.val,
              theta.hat=theta.hat,
              theta.reshuffled=theta.reshuffled))
}


### Auxiliary function
permutation.iter = function(d,y,X,V,lambda){
  dstar = sample(d)
  X0 = t(X[dstar==0,]); X1 = t(X[dstar==1,]);
  Y0 = y[dstar==0]; Y1 = y[dstar==1];
  solstar = regsynth(X0,X1,Y0,Y1,V,lambda)
  return(solstar$ATT)
}