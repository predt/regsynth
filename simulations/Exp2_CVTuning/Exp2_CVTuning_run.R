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
set.seed(2121988)
lambda = seq(.01,8,.1) # set of lambda to be considered for optim, set it high enough at first
K = 5 # number of folds for optimal penalty level


### Setup of the MC experiment
AlbertoMCXP_setup <- function(R=1000,delta=1,n1=1,n0=50,sigma=1){
  lambda.opt.record = vector(length=R)
  t_start = Sys.time()
  pb = txtProgressBar(style = 3)
  
  for(r in 1:R){
    ### 0. Generate data
    data = SimpleAlbertoDGP(delta,n1,n0,sigma)
    X = data$X; y = data$y; d = data$d
    
    X0 = t(as.matrix(X)[d==0,]); V = diag(1)
    Y0 = y[d==0]
    
    ### Splitting the sample for cross-validation
    uu=0 # reshuffle groups until no empty group
    while(uu==0){
      allocation = sample(1:K,n0,replace=T)
      uu=min(mapply(function(x) sum(allocation==x),1:K))
    }
    
    ### Computing several lambda
    keeptau = matrix(nrow=length(lambda), ncol=length(Y0))
    for(k in 1:K){
      X1k = t(as.matrix(X0[,allocation==k]))
      X0k = t(as.matrix(X0[,allocation!=k]))
      Y1k = Y0[allocation==k]
      Y0k = Y0[allocation!=k]
      solpath = regsynthpath(X0k,X1k,Y0k,Y1k,V,lambda)
      keeptau[,allocation==k] = solpath$CATT
    }
    
    ### Optimizes RMSE
    curve.RMSE = apply(keeptau^2,1,sum)/n0
    lambda.opt.RMSE = min(lambda[which(curve.RMSE==min(curve.RMSE))])
    lambda.opt.record[r] = lambda.opt.RMSE
  
    print("*** PROGRESS ***")
    print(100*r/R)
    setTxtProgressBar(pb, r/R)
  }
  
  close(pb)
  print(Sys.time()-t_start)

  fileN = paste("simulations/cluster/EXP2SumPoly_delta",delta,".txt",sep="")
  
  write(c("",
          paste("Statistics on lambda CV")), fileN,append=TRUE)
  
  write.table(c(summary(lambda.opt.record)),sep="  ",col.names=F,file=fileN,append=TRUE)
  
  write(c("",
          paste("Nb. treated:",n1),
          paste("Nb. controls:",n0),
          paste("Nb. replications:",R),
          paste("Order of Legendre polynomial:",delta),
          paste("Noise standard deviation:",sigma),
          "",
          paste(Sys.time())), fileN, append=TRUE)
}


### Launch simulations
kappa_set = c(1:12,15,20,25)
for(deg_lin in kappa_set){
  AlbertoMCXP_setup(R=300,n1=1,n0=50,delta=deg_lin,sigma=.1)
}


### Collect results and generate charts
substrRight = function(x, n) substr(x, nchar(x)-n+1, nchar(x))


lambdaCV = vector(length=length(kappa_set))
i = 0
for(kappa in kappa_set){
  i = i+1
  fileN = paste("simulations/cluster/EXP2SumPoly_delta",kappa,".txt",sep="")
  con = file(fileN, "r")
  uu = readLines(con, n = 6)
  close(con)
  lambdaCV[i] = as.numeric(substrRight(uu[6],6)) 
}



pdf("simulations/NonLinearity.pdf", width=6, height=6)
matplot(kappa_set, lambdaCV, type="b", pch=20, lwd=1,
        main=expression("Selected penalty level, "*lambda^{CV}*", as a function of Legendre polynomial order"), col="steelblue",
        xlab=expression(kappa), ylab=expression(lambda^{CV}))
dev.off()
