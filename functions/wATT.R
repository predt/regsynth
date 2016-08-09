#' Compute weighted estimate for Treatment on the treated
#' 


wATT <- function(y,d,w){
  y <- as.vector(y)
  d <- as.vector(d)
  theta <- mean(y[d==1]) - t(y[d==0])%*%w
  return(c(theta))
}