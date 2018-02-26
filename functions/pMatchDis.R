#' Matching discrepancy density function
#' 
#' @param v is the point at which density is computed (should be a vector)
#' @param z is the discrepency point (should be a vector)
#' @param k is the number of continuous variables
#' @param f is the density of the vector

dMatchDis <- function(v,z,m=1,f=dnorm){
  k = length(v)
  C0 = (2*f(z)*pi^(k/2))/(k*gamma(k/2))
  y = (f(z) / factorial(m-1)) * (C0*enorm(v)^k)^(m-1) * exp(-C0*enorm(v)^k)
  return(y)
}

enorm = function(x) sqrt(sum(x^2))


v = seq(-10,10,.01)
z = 2
u = NULL
for(i in v){
  u = c(u,dMatchDis(i,z,m=1,dnorm))
}

plot(v,u,type="l")

### MC experiment to confirm
R = 10000
n = 1000
uh = vector(length=R)
for(r in 1:R){
  x = rnorm(n)
  rk = rank((x-2)^2)
  uh[r] = (x[which(rk==1)]-2)*n
}

library("ggplot2")
denMatch = function(x) mapply(function(y) dMatchDis(y,z,m=1,dnorm), x)
ggplot(data.frame(uh), aes(x=uh)) + 
    geom_histogram(bins=200,alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
    scale_x_continuous(name="U") +
    ggtitle("Discrepancy density for 3rd largest") + 
    stat_function(fun = denMatch, colour="darkorchid3", size=1) +
    theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")