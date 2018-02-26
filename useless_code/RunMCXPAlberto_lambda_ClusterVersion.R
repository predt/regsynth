### Run MC_XP, with different Legendre polynomial order
### Goal: see how lambda varies
### Jeremy L Hour
### 24 juillet 2017

setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")
rm(list=ls())

### 0. Settings

### Load packages
library("MASS")
library("ggplot2")
library("gtable")
library("grid")
library("reshape2")
library("LowRankQP")
library("xtable")
library("orthopolynom")

### Load user functions
source("functions/wsol.R")
source("functions/wsoll1.R")
source("functions/wATT.R")
source("functions/matching.R")
source("functions/matchest.R")
source("functions/OBest.R")
source("functions/regsynth.R")
source("functions/regsynthpath.R")
source("functions/TZero.R")
source("functions/synthObj.R")
source("functions/SimpleAlbertoDGP.R")
source("functions/AlbertoDGP.R")


### MC XP: how does lambda moves with the degree of linearity?
set.seed(12071990)
lambda = seq(.01,5,.1) # set of lambda to be considered for optim, set it high enough at first
K = 5 # number of folds for optimal penalty level


### Settup of the MC experiment
AlbertoMCXP_setup <- function(R=1000,n=100,delta=1,prop=.2,a=.5,sigma=1,unifcase=F){
  lambda.opt.record = vector(length=R)
  t_start = Sys.time()
  pb = txtProgressBar(style = 3)
  
  for(r in 1:R){
    ### 0. Generate data
    data = SimpleAlbertoDGP(n,delta,prop,a,sigma)
    X = data$X; y = data$y; d = data$d
    
    X0 = t(as.matrix(X)[d==0,]); X1 = t(as.matrix(X)[d==1,]); V = diag(ncol(as.matrix(X)))
    Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)
    
    ### Splitting the sample for cross-validation
    uu=0 # reshuffle groups until no empty group
    while(uu==0){
      allocation = sample(1:K,n0,replace=T)
      uu=min(mapply(function(x) sum(allocation==x),1:K))
    }
    
    ### COmputing sevral lambda
    keeptau = matrix(nrow=length(lambda), ncol=length(Y0))
    for(k in 1:K){
      X1k = t(as.matrix(X0[,allocation==k]))
      X0k = t(as.matrix(X0[,allocation!=k]))
      Y1k = Y0[allocation==k]
      Y0k = Y0[allocation!=k]
      solpath = regsynthpath(X0k,X1k,Y0k,Y1k,V,lambda)
      keeptau[,allocation==k] = solpath$CATT
    }
    
    # The one that optimizes RMSE
    curve.RMSE = apply(keeptau^2,1,sum)/n0
    lambda.opt.RMSE = min(lambda[which(curve.RMSE==min(curve.RMSE))])
    lambda.opt.record[r] = lambda.opt.RMSE
    
    
    print("*** PROGRESS ***")
    print(100*r/R)
    setTxtProgressBar(pb, r/R)
  }
  
  close(pb)
  print(Sys.time()-t_start)
  
  fileN = paste("simulations/cluster/AlbertoXPoutput_delta",delta,"unifCASE.txt",sep="")
  
  write(c("",
          paste("Statistics on lambda CV")), fileN,append=TRUE)
  
  write.table(c(summary(lambda.opt.record)),sep="  ",col.names=F,file=fileN,append=TRUE)
  
  write(c("",
          paste("Nb. observations:",n),
          paste("Nb. replications:",R),
          paste("Order of Legendre polynomial:",delta),
          paste("Noise standard deviation:",sigma),
          "",
          paste(Sys.time())), fileN, append=TRUE)
}


### Launch simulations
for(deg_lin in c(7,5,4,3,2,1)){
  AlbertoMCXP_setup(R=1000,n=125,delta=deg_lin,prop=.2,a=0,sigma=.001,unifcase=T)
}