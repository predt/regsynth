### Run MC_XP_ Sub
### Why is the bias decreasing?
### Keep only NN estimators for speed
### Jeremy L Hour
### 19/02/2018

setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")
rm(list=ls())

### 0. Settings

### Load packages
library("MASS")
library("gtable")
library("grid")
library("reshape2")
library("xtable")

### Load user functions
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


### 0. Setup
MCXP_NN <- function(R=1000,n=100,p=50){
  Results = matrix(ncol=6, nrow=R)
  t_start = Sys.time()
  pb = txtProgressBar(style = 3)
  
  for(r in 1:R){
    ### 0. Generate data
    data = matchDGP(n=n,p=p)
    X = data$X; y = data$y; d = data$d
    
    X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))
    Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)
    
    ### Splitting the sample for cross-validation
    uu=0 # reshuffle groups until no empty group
    while(uu==0){
      allocation = sample(1:K,n0,replace=T)
      uu=min(mapply(function(x) sum(allocation==x),1:K))
    }
    
    ### 1. NN Matching
    ## A. K = 1
    NN1 = matchest(X0,X1,Y0,Y1,V,m=1)
    ## B. K = 5
    NN5 = matchest(X0,X1,Y0,Y1,V,m=5)
    ## C. K = Kopt
    keeptauNN = matrix(nrow=10, ncol=length(Y0))
    for(k in 1:K){
      X1k = as.matrix(X0[,allocation==k])
      X0k = as.matrix(X0[,allocation!=k])
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
    
    ### 4. ATT estimation
    Results[r,] <- c(NN1$ATT,NN5$ATT,NN.opt.RMSE,NN.opt.bias,NN.opt.crit,NN.opt.mae)
    setTxtProgressBar(pb, r/R)
  }
  
  close(pb)
  print(Sys.time()-t_start)
  
  ### Compute bias and RMSE
  StatDisplay <- data.frame()
  StatDisplay[1:6,"bias"] = abs(apply(Results,2,mean))
  StatDisplay[1:6,"RMSE"] = sqrt(apply(Results^2,2,mean))
  StatDisplay[1:6,"MAE"] = apply(abs(Results),2,mean)
  row.names(StatDisplay) = c("1NN Matching","5NN Matching",
                             "NN RMSE opt","NN bias opt","NN crit opt", "NN MAE opt")
  print(StatDisplay)
  
  fileN = paste("simulations/SquareOfunc_BiasDecreaseoutput_n",n,",p",p,".txt",sep="")
  
  print.xtable(xtable(StatDisplay, digits=3),type="latex",file=fileN)
  write(c(paste("Nb. observations:",n),
          paste("Nb. covariates:",p),
          paste("Nb. replications:",R),
          paste(Sys.time())), fileN, append=TRUE)
  
}

### MC XP
set.seed(2121988)
lambda = seq(0,2,.01) # set of lambda to be considered for optim
K = 5 # number of folds for optimal penalty level

for(n_xp in c(30,50,100)){
  for(p_xp in c(3,4,5,6,7,10,20,30,40)){
    print(paste("n=",n_xp))
    print(paste("p=",p_xp))
    MCXP_NN(R=1000,n=n_xp,p=p_xp)
  }
}

### Collect results and generate charts
substrRight = function(x, n) substr(x, nchar(x)-n+1, nchar(x))

p_set = c(3,4,5,6,7,10,20,30,40)

bias = matrix(nrow=length(p_set), ncol=3)
j = 0
for(n in c(30,50,100)){
  j=j+1; i=0
  for(p in p_set){
    i = i+1
    fileN = paste("simulations/BiasDecreaseoutput_n",n,",p",p,".txt",sep="")
    con = file(fileN, "r")
    uu = readLines(con, n = 9)
    close(con)
    bias[i,j] = as.numeric(substr(substrRight(uu[9],25),1,5)) 
  }
}



pdf("simulations/EXP_Bias_withdimension.pdf", width=6, height=6)
matplot(p_set, bias, type="b", pch=20, lwd=1,
        main=expression("Bias for 1NN as dimension p increases"), col="steelblue",
        xlab=expression(p), ylab="bias")
dev.off()

