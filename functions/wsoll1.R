#' Synthetic Control with l1-penalty
#' 
#' Function to compute synthetic control weights given a V matrix
#' 
#' @param X0 is a p x n0 matrix
#' @param X1 is a p x 1 vector
#' @param V is a p x p matrix of weights
#' @param pen is l1 penalty level
#' 
#' @autor Jeremy LHour

wsoll1 <- function(X0,X1,V,pen=0.0){
  n = ncol(X0)
  Delta = matrix(t(X1)%*%V%*%X1, nrow=n, ncol=1) - 2*t(X0)%*%V%*%X1 + diag(t(X0)%*%V%*%X0)
  
  P = 2*t(X0)%*%V%*%X0
  q = t(-2*t(X0)%*%V%*%X1 + pen*Delta)
  
  sol = LowRankQP(Vmat=P,dvec=q,Amat=matrix(1, ncol=n),bvec=1,uvec=rep(1,n), method="LU")
  return(sol$alpha)
}