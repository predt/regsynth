### Behavior of the estimator as lambda varies
### Jeremy L Hour
### 25 juillet 2017

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

### Load user functions
source("functions/wsol.R")
source("functions/wsoll1.R")
source("functions/matchDGP.R")
source("functions/wATT.R")
source("functions/matching.R")
source("functions/matchest.R")
source("functions/OBest.R")
source("functions/regsynth.R")
source("functions/regsynthpath.R")
source("functions/TZero.R")
source("functions/synthObj.R")
source("simulations/MCXP_setup.R")


### 0. Settings
set.seed(12071990)
lambda = seq(0,2,.01) # set of lambda to be considered for optim
K = 5 # number of folds for optimal penalty level
R = 100
n = 20
p = 5

### 1. Loop
Results <- matrix(ncol=length(lambda)+3, nrow=R)
t_start <- Sys.time()
pb <- txtProgressBar(style = 3)
  
for(r in 1:R){
    ### A. Generate data
    data = matchDGP(n=n,p=p,Ry=.5,Rd=.2)
    X = data$X; y = data$y; d = data$d
    
    X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))
    Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)
    
    ### B. Compute solution for all lambdas
    solpath = regsynthpath(X0,X1,Y0,Y1,V,lambda)
    Results[r,1:length(lambda)] = solpath$ATT
    
    ### C. Splitting the sample for cross-validation
    uu=0 # reshuffle groups until no empty group
    while(uu==0){
      allocation = sample(1:K,n0,replace=T)
      uu=min(mapply(function(x) sum(allocation==x),1:K))
    }
    
    keeptau = matrix(nrow=length(lambda), ncol=length(Y0))
    for(k in 1:K){
      X1k = as.matrix(X0[,allocation==k])
      X0k = as.matrix(X0[,allocation!=k])
      Y1k = Y0[allocation==k]
      Y0k = Y0[allocation!=k]
      solpath = regsynthpath(X0k,X1k,Y0k,Y1k,V,lambda)
      keeptau[,allocation==k] = solpath$CATT
    }
    
    # The one that optimizes RMSE
    curve.RMSE = apply(keeptau^2,1,sum)/n0
    lambda.opt.RMSE = min(lambda[which(curve.RMSE==min(curve.RMSE))])
    Results[r,length(lambda)+1] = solpath$ATT[which(lambda==lambda.opt.RMSE)]
    
    # The one that optimizes bias
    curve.bias = abs(apply(keeptau,1,sum)/n0)
    lambda.opt.bias = min(lambda[which(curve.bias==min(curve.bias))])
    Results[r,length(lambda)+2] = solpath$ATT[which(lambda==lambda.opt.bias)]
    
    # The one that optimizes bias + variance
    curve.crit = curve.bias + apply(keeptau,1,sd)
    lambda.opt.crit = min(lambda[which(curve.crit==min(curve.crit))])
    Results[r,length(lambda)+3] = solpath$ATT[which(lambda==lambda.opt.crit)]
    
    print("*** PROGRESS ***")
    print(100*r/R)
    setTxtProgressBar(pb, r/R)
}
  
close(pb)
print(Sys.time()-t_start)
  
### 2. Compute bias and RMSE
StatDisplay <- data.frame()
StatDisplay[1:(length(lambda)+3),"bias"] <- abs(apply(Results,2,mean))
StatDisplay[1:(length(lambda)+3),"RMSE"]  <- sqrt(apply(Results^2,2,mean))
StatDisplay[1:(length(lambda)+3),"Variance"]  <- apply(Results,2,sd)

matplot(lambda, StatDisplay[1:(length(lambda)),"bias"], type="b", pch=20, lwd=1,
        main=expression("Average bias as a function of "*lambda*"."), col="steelblue",
        xlab=expression(lambda), ylab="bias")

matplot(lambda, StatDisplay[1:(length(lambda)),"Variance"], type="b", pch=20, lwd=1,
        main=expression("Average Variance as a function of "*lambda*"."), col="steelblue",
        xlab=expression(lambda), ylab="Variance")

matplot(lambda, StatDisplay[1:(length(lambda)),"RMSE"], type="b", pch=20, lwd=1,
        main=expression("Average MSE as a function of "*lambda*"."), col="steelblue",
        xlab=expression(lambda), ylab="MSE")
