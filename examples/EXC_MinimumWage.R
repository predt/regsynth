### EXAMPLE 3: Minimum wage, Card and Krueger
### Sparse Synthetic Control
### Jeremy L Hour
### 13 decembre 2016

setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")

rm(list=ls())
set.seed(3101990)

### Load packages
library("MASS")
library("ggplot2")
library("gtable")
library("grid")
library("reshape2")
library("LowRankQP")
library("xtable")
library("foreign")

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

### Load data
data = read.dta("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth/data/fastfood.dta")

data[,"BK"] = data[,"chain"] == 1
data[,"KFC"] = data[,"chain"] == 2
data[,"RR"] = data[,"chain"] == 3


X = data[,c("BK","KFC","RR","co_owned",
            "empft","emppt","nmgrs","inctime")]
y = data[,"empft2"]
uu = complete.cases(cbind(X,y))
d = data[uu,"state"]==1
X = X[uu,]
y = y[uu]


### Use penalized synthetic control
X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))
Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)

# A. lambda = .1
V[5,5] = 100
sol = regsynth(X0,X1,Y0,Y1,V,.1)
RSC.fixed = sol$ATT

apply(X1 - X0%*%t(sol$Wsol) ,1,mean)

lambda = seq(0,2,.01) # set of lambda to be considered for optim
K = 5 # number of folds for optimal penalty level

### Splitting the sample for cross-validation
uu=0 # reshuffle groups until no empty group
while(uu==0){
  allocation = sample(1:K,n0,replace=T)
  uu=min(mapply(function(x) sum(allocation==x),1:K))
}


# B. lambda = lambdaopt
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

c(RSC.fixed,RSC.opt.RMSE,RSC.opt.bias,RSC.opt.crit)