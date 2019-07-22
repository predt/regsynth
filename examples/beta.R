m.beta <- function(a,b,k=1){
  return(beta(a+k,b)/beta(a,b))
}

cov.sum <- function(x){
  a = x[1]
  b = x[2]
  y = m.beta(a,b,2)*(1-m.beta(b,a,2))-2*m.beta(a,b,3)+m.beta(a,b,4)
  return(y)
}

optim(c(1,1),cov.sum)

for(a in c(.05,0.5,1:10)){
  for(b in c(.05,0.5,1:10)){
    print(cov.sum(a,b))
  }
}