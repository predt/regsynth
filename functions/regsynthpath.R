#' Regularized Synthetic Control path
#' 
#' Created: 20 juillet 2016
#' Loops over regsynth to get regularized synthetic control solutions
#' For all given values of lambda
#' 
#' @param X0 is a p x n0 matrix
#' @param X1 is a p x n1 vector
#' @param Y0 is a n0 x 1 vector
#' @param Y1 is a n1 x 1 vector
#' @param V is a p x p matrix of weights
#' @param lambda is a vector of l1 penalty levels
#' @param tol gives the threshold for considering true zeros
#' 
#' @seealso \code{\link{regsynth}}
#' 
#' @autor Jeremy LHour

regsynthpath <- function(X0,X1,Y0,Y1,V,lambda,tol=1e-6){
  K = length(lambda); n0 = ncol(X0); n1 = ncol(X1)
  ATT = vector(length = K)
  tau = matrix(nrow=K,ncol=n1)
  Wsol = array(dim=c(K,n1,n0))
  
  t_start <- Sys.time()
  pb <- txtProgressBar(style = 3)
  for(k in 1:K){
    sol = regsynth(X0,X1,Y0,Y1,V,lambda[k],tol=1e-6)
    ATT[k] = sol$ATT
    tau[k,] = sol$CATT
    Wsol[k,,] = sol$Wsol
    setTxtProgressBar(pb, k/K)
  }
  close(pb)
  print(Sys.time()-t_start)
  
  return(list(ATT=ATT,
              CATT=tau,
              Wsol=Wsol))
}
