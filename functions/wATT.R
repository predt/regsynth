#' Compute weighted estimate for Treatment on the treated
#' 
#' @param d is a vector of dimension n (treatment indicator)
#' @param y is a vector of dimension n (outcome)
#' @param w is a vector of dimension n0 (weights)

wATT <- function(y,d,w){
  y = as.vector(y); d = as.vector(d);
  theta <- mean(y[d==1]) - t(y[d==0])%*%w
  return(c(theta))
}