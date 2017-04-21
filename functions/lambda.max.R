#' Compute lambda.max
#' 
#' Function to compute the lambda.max, i.e. the minimal value of lambda
#' that gives the Nearest-neighbor solution.
#' 
#' WARNING: KKT are necessary but not sufficient in that case...
#' WARNING: Does not work.
#' 
#' @param X0 is a p x n0 matrix
#' @param X1 is a p x 1 vector
#' @param V is a p x p matrix of weights
#' 
#' @autor Jeremy LHour

lambda.max <- function(X0,X1,V){
  n = ncol(X0)
  Delta = matrix(t(X1)%*%V%*%X1, nrow=n, ncol=1) - 2*t(X0)%*%V%*%X1 + diag(t(X0)%*%V%*%X0)
  r = rank(Delta)
  i_star = which(r==1)
  
  if(length(i_star)>1){
    print("WARNING: Nearest-Neighbor is not unique. Smallest index taken.")
    i_star = min(i_star)
  } 
  
  lambda.pot = t(X0 - X0[,i_star]%*%matrix(1, ncol=n))%*%V%*%(X1-X0[,i_star]) / (Delta - Delta[i_star])
  lambda.pot[i_star] = 0
  lambda.max = max(lambda.pot)
  
  return(lambda.max)
}