### Exp 5: test
### See what is going wrong


setwd("/Users/jeremylhour/Documents/R/regsynth-master")
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
source("functions/PolyDGP.R")
source("functions/wATT.R")
source("functions/matching.R")
source("functions/matchest.R")
source("functions/OBest.R")
source("functions/regsynth.R")
source("functions/regsynthpath.R")
source("functions/TZero.R")
source("functions/synthObj.R")
source("simulations/Exp5_PolyDGP/Exp5_Poly_setup.R")


### MC XP
set.seed(12071990)
lambda = seq(0,2,.01) # set of lambda to be considered for optim
K = 5 # number of folds for optimal penalty level
n1 = 1000
n0 = 100
p = 2
delta = 2

### 0. Generate data
data = PolyDGP(n1,n0,p,delta)
X = data$X; y = data$y; d = data$d

summary(lm(y[d==1]  ~ X[d==1]))
  
X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))
Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)
    
    
# 1NN Matching
NN1 = matchest(X0,X1,Y0,Y1,V,m=1)
print(NN1$ATT)
solvec = vector(length=30)
for(i in 1:30){
  NN = matchest(X0,X1,Y0,Y1,V,m=i)
  solvec[i] = NN$ATT
}

plot(1:30,solvec)

# Large lambda penSynth
sol = regsynth(X0,X1,Y0,Y1,V,1e6)
print(sol$ATT)

lambda = seq(.01,5,.01)
solpath = regsynthpath(X0,X1,Y0,Y1,diag(p),lambda)

plot(lambda,solpath$ATT)

    ## B. K = 5
    NN5 = matchest(X0,X1,Y0,Y1,V,m=5)
    ## C. K = Kopt
    keeptauNN = matrix(nrow=10, ncol=length(Y0))
    for(k in 1:K){
      X1k = matrix(X0[,allocation==k], nrow=p)
      X0k = matrix(X0[,allocation!=k], nrow=p)
      Y1k = Y0[allocation==k]
      Y0k = Y0[allocation!=k]
      for(i in 1:10){
        soli = matchest(X0k,X1k,Y0k,Y1k,V,m=i)
        keeptauNN[i,allocation==k] = soli$CATT
      }
    }
    
    # The one that optimizes RMSE
    curve.RMSE = apply(keeptauNN^2,1,sum)/n0
    m.opt.RMSE = min(which(curve.RMSE==min(curve.RMSE)))
    sol = matchest(X0,X1,Y0,Y1,V,m=m.opt.RMSE)
    NN.opt.RMSE = sol$ATT
    
    # The one that optimizes bias
    curve.bias = abs(apply(keeptauNN,1,sum)/n0)
    m.opt.bias = min(which(curve.bias==min(curve.bias)))
    sol = matchest(X0,X1,Y0,Y1,V,m=m.opt.bias)
    NN.opt.bias = sol$ATT
    
    # The one that optimizes bias + variance
    curve.crit = curve.bias + apply(keeptauNN,1,sd)
    m.opt.crit = min(which(curve.crit==min(curve.crit)))
    sol = matchest(X0,X1,Y0,Y1,V,m=m.opt.crit)
    NN.opt.crit = sol$ATT
    
    # The one that optimizes MAE
    curve.mae = apply(abs(keeptauNN),1,sum)/n0
    m.opt.mae = min(which(curve.mae==min(curve.mae)))
    sol = matchest(X0,X1,Y0,Y1,V,m=m.opt.mae)
    NN.opt.mae = sol$ATT
    
    
    ### 3. Regularized Synthetic Control
    # A. lambda = .1
    sol = regsynth(X0,X1,Y0,Y1,V,.1)
    RSC.fixed = sol$ATT
    
    # B. lambda = lambdaopt
    keeptau = matrix(nrow=length(lambda), ncol=length(Y0))
    for(k in 1:K){
      X1k = matrix(X0[,allocation==k], nrow=p)
      X0k = matrix(X0[,allocation!=k], nrow=p)
      Y1k = Y0[allocation==k]
      Y0k = Y0[allocation!=k]
      solpath = regsynthpath(X0k,X1k,Y0k,Y1k,V,lambda)
      keeptau[,allocation==k] = solpath$CATT
    }
    
    # The one that optimizes RMSE
    curve.RMSE = apply(keeptau^2,1,sum)/n0
    lambda.opt.RMSE = min(lambda[which(curve.RMSE==min(curve.RMSE))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.RMSE)
    RSC.opt.RMSE = sol$ATT
    
    # The one that optimizes bias
    curve.bias = abs(apply(keeptau,1,sum)/n0)
    lambda.opt.bias = min(lambda[which(curve.bias==min(curve.bias))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.bias)
    RSC.opt.bias = sol$ATT
    
    # The one that optimizes bias + variance
    curve.crit = curve.bias + apply(keeptau,1,sd)
    lambda.opt.crit = min(lambda[which(curve.crit==min(curve.crit))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.crit)
    RSC.opt.crit = sol$ATT
    
    # The one that optimizes MAE
    curve.mae = apply(abs(keeptau),1,sum)/n0
    lambda.opt.mae = min(lambda[which(curve.mae==min(curve.mae))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.mae)
    RSC.opt.mae = sol$ATT
  