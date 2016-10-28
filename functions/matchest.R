#' Matching estimator for a whole sample
#' 
#' From a sample of treated and control units, computes the weights
#' for each counterfactual, the individual treatment effects and counterfactuals,
#' and the Average Treatment on the Treated (ATT). Based on Neirest-Neighbor matching.
#' 
#' Edited: 28 octobre 2016
#' 
#' @param X0 is a p x n0 matrix
#' @param X1 is a p x n1 matrix
#' @param Y0 is a n0 x 1 vector
#' @param Y1 is a n1 x 1 vector
#' @param V is a p x p matrix of weights
#' @param m number of neighbors to find
#' 
#' @return ATT is the Average Treatment Effect on the Treated over the sample.
#' @return CATT is the individual treatment effect.
#' @return Wsol is the n1 x n0 matrix of weights to compute counterfactual.
#' @return y0_hat is the individual counterfactual for each treated unit.
#' 
#' @seealso \code{\link{matching}} for the computation of the matching weights.
#' 
#' @author Jeremy LHour

matchest <- function(X0,X1,Y0,Y1,V,m=3){ 
  n1 = ncol(X1); n0 = ncol(X0);
  
  Wsol = matrix(nrow=n1,ncol=n0)
  
  for(i in 1:n1){
    sol = matching(X0,X1[,i],V,m=m)
    Wsol[i,] = sol
  }
  
  y0_hat = Wsol%*%Y0
  tau = Y1 - y0_hat
  
  return(list(ATT = mean(tau),
              CATT = tau,
              Wsol = Wsol,
              y0_hat = y0_hat))
}
