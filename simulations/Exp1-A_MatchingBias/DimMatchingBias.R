### Dimension and matching bias: a simulation
### Jeremy L Hour
### 14/02/2018

setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")
rm(list=ls())

### 0. Settings

### Load packages
library("MASS")

### Load user functions
source("functions/matching.R")

### DGP
expDGP <- function(n=20,p=2,rho=.5){
  ### Covariate variance matrix
  Sigma = matrix(0,nrow=p, ncol=p)
  for(k in 1:p){
    for(j in 1:p){
      Sigma[k,j] = rho^abs(k-j)
    }
  }
  
  b = rep(0,p)
  for(j in 1:p){
    b[j] = 2*(-1)^(j+1) / (p-j+1)^2
  }
  
  X = mvrnorm(n = n, mu=rep(0,p), Sigma)
  y = exp(-X%*%b) + rnorm(n)
  
  return(list(X=X,
              y=y))
}


### Simulations: EXP A
R = 100; n=100; p=10

Results = matrix(ncol=1, nrow=R)
t_start = Sys.time()
pb = txtProgressBar(style = 3)
for(r in 1:R){
  ### 0. Generate data
  data = expDGP(n=n,p=p)
  X = data$X; y = data$y
  
  ### 1. Compute 1NN and differences in outcomes
  tau = vector(length=n)
  for(i in 1:n){
    X0 = t(X[-i,]); Y0 = y[-i]
    X1 = X[i,]; Y1 = y[i]
    W = matching(X0,X1,diag(p),m=1)
    tau[i] = Y1 - W%*%Y0
  }
  
  ### 2. Compute Statistics
  Results[r] = mean(tau)
  setTxtProgressBar(pb, r/R)
}

close(pb)
print(Sys.time()-t_start)

### Print results
print(abs(mean(Results)))

### Simulations: EXP B
R = 100; n=100; p=10

Results = matrix(ncol=1, nrow=R)
t_start = Sys.time()
pb = txtProgressBar(style = 3)
for(r in 1:R){
  ### 0. Generate data
  data = expDGP(n=n,p=p)
  X = data$X; y = data$y
  
  ### 1. Compute 1NN and discrepancy
  X0 = t(X[-1,]);X1 = X[1,]
  W = matching(X0,X1,diag(p),m=1)
  Results[r] = sqrt(sum((X1 - X0%*%W)^2))
  setTxtProgressBar(pb, r/R)
}

close(pb)
print(Sys.time()-t_start)

### Print results
print(abs(mean(Results)))

### Simulations: EXP C
R = 100; n=100; p=5

Results = matrix(ncol=1, nrow=R)
t_start = Sys.time()
pb = txtProgressBar(style = 3)
for(r in 1:R){
  ### 0. Generate data
  data = expDGP(n=n,p=p)
  X = data$X; y = data$y
  
  ### 1. Compute 1NN and discrepancy
  tau = vector(length=n)
  for(i in 1:n){
    X0 = t(X[-i,]); Y0 = y[-i]
    X1 = X[i,]; Y1 = y[i]
    W = matching(X0,X1,diag(p),m=1)
    tau[i] = sqrt(sum((X1 - X0%*%W)^2))
  }
  
  Results[r] = mean(tau)
  setTxtProgressBar(pb, r/R)
}

close(pb)
print(Sys.time()-t_start)

### Print results
print(abs(mean(Results)))