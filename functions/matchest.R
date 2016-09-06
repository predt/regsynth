#' Matching estimator
#' 
#' Returns the weights for each control unit and the counterfactual
#' 
#' @param d is a vector of dimension n (treatment indicator)
#' @param X is a n x p matrix
#' @param y is a vector of dimension n (outcome)
#' @param V is a p x p matrix of weights
#' @param m number of neighbors to find

matchest <- function(d,X,y,V,m=3){
  
  X0 = t(X[d==0,]); X1 = t(X[d==1,]);
  Y0 = y[d==0]; Y1 = y[d==1]; 
  n1 = sum(d); n0 = sum(1-d);
  
  tau = vector(length=n1)
  y0_hat = vector(length=n1)
  Wsol = matrix(nrow=n1,ncol=n0)
  
  for(i in 1:n1){
    sol = matching(X0,X1[,i],V,m=m)
    Wsol[i,] = sol
    y0_hat[i] = t(Y0)%*%sol
    tau[i] = Y1[i] - t(Y0)%*%sol
  }
  
  return(list(ATT = mean(tau),
              CATT = tau,
              Wsol = Wsol,
              y0_hat = y0_hat))
}
