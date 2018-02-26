### Run MC_XP, with different Legendre polynomial order
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


### I. MC XP: simple Alberto DGP
set.seed(12071990)
lambda = seq(.0001,2,.01) # set of lambda to be considered for optim
K = 5 # number of folds for optimal penalty level


### Settup of the MC experiment
AlbertoMCXP_setup <- function(R=1000,n=100,delta=1,prop=.1,a=.5,sigma=1){
  Results = matrix(ncol=10, nrow=R)
  lambda.opt.record = vector(length=R)
  t_start = Sys.time()
  pb = txtProgressBar(style = 3)
  
  for(r in 1:R){
    ### 0. Generate data
    data = SimpleAlbertoDGP(n,delta,prop,a,sigma)
    X = data$X; y = data$y; d = data$d
    
    X0 = t(as.matrix(X)[d==0,]); X1 = t(as.matrix(X)[d==1,]); V = diag(ncol(as.matrix(X)))
    Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)
    
    ### 1. Synthetic Control on mean of treated
    M = matrix(apply(X1,1,mean), ncol=1)
    AggSC = wATT(y,d,wsol(X0,M,V))
    
    ### Splitting the sample for cross-validation
    uu=0 # reshuffle groups until no empty group
    while(uu==0){
      allocation = sample(1:K,n0,replace=T)
      uu=min(mapply(function(x) sum(allocation==x),1:K))
    }
    
    ### 2. NN Matching
    ## A. K = 1
    NN1 = matchest(X0,X1,Y0,Y1,V,m=1)
    ## B. K = 5
    NN5 = matchest(X0,X1,Y0,Y1,V,m=5)
    ## C. K = Kopt
    keeptauNN = matrix(nrow=10, ncol=length(Y0))
    for(k in 1:K){
      X1k = t(as.matrix(X0[,allocation==k]))
      X0k = t(as.matrix(X0[,allocation!=k]))
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
    sol = matchest(X0k,X1k,Y0k,Y1k,V,m=m.opt.RMSE)
    NN.opt.RMSE = sol$ATT
    
    # The one that optimizes bias
    curve.bias = abs(apply(keeptauNN,1,sum)/n0)
    m.opt.bias = min(which(curve.bias==min(curve.bias)))
    sol = matchest(X0k,X1k,Y0k,Y1k,V,m=m.opt.bias)
    NN.opt.bias = sol$ATT
    
    # The one that optimizes bias + variance
    curve.crit = curve.bias + apply(keeptauNN,1,sd)
    m.opt.crit = min(which(curve.crit==min(curve.crit)))
    sol = matchest(X0k,X1k,Y0k,Y1k,V,m=m.opt.crit)
    NN.opt.crit = sol$ATT
    
    
    ### 3. Regularized Synthetic Control
    # A. lambda = .1
    sol = regsynth(X0,X1,Y0,Y1,V,.1)
    RSC.fixed = sol$ATT
    
    # B. lambda = lambdaopt
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
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.RMSE)
    RSC.opt.RMSE = sol$ATT
    lambda.opt.record[r] = lambda.opt.RMSE
    
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
    
    print("*** PROGRESS ***")
    print(100*r/R)
    
    ### 4. ATT estimation
    Results[r,] <- c(AggSC,NN1$ATT,NN5$ATT,
                     NN.opt.RMSE,NN.opt.bias,NN.opt.crit,
                     RSC.fixed,
                     RSC.opt.RMSE,RSC.opt.bias,RSC.opt.crit)
    setTxtProgressBar(pb, r/R)
  }
  
  close(pb)
  print(Sys.time()-t_start)
  
  ### Compute bias and RMSE
  StatDisplay <- data.frame()
  StatDisplay[1:10,"bias"] <- abs(apply(Results-a,2,mean))
  StatDisplay[1:10,"RMSE"]  <- sqrt(apply((Results-a)^2,2,mean))
  row.names(StatDisplay) <- c("Aggregate Synth","1NN Matching","5NN Matching",
                              "NN RMSE opt","NN bias opt","NN crit opt",
                              "Penalized Synth fixed",
                              "Penalized Synth RMSE opt","Penalized Synth bias opt","Penalized Synth crit opt")
  print(StatDisplay)
  
  fileN = paste("simulations/AlbertoXPoutput_delta",delta,".txt",sep="")
  
  print.xtable(xtable(StatDisplay, digits=3),type="latex",file=fileN)
  
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
for(deg_lin in seq(1,6,by=1)){
  AlbertoMCXP_setup(R=1000,n=50,delta=deg_lin,prop=.1,a=0)
}



### Compute the bias term
set.seed(12071990)
lambda = seq(.00001,2,.01)
delta = 4
data = SimpleAlbertoDGP(n=500,delta=delta,Ry=.5,prop=.05,a=0)
X = data$X; y = data$y; d = data$d; beta = data$b

X0 = t(as.matrix(X)[d==0,]); X1 = t(as.matrix(X)[d==1,]); V = diag(ncol(as.matrix(X)))
Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)

# example with one solution: how does bias term evolves with lambda?
poly = legendre.polynomials(delta)
BT = vector(length=length(lambda))
VAR = vector(length=length(lambda))
uu = 0
for(ll in lambda){
  uu = uu+1
  sol = regsynth(X0,X1,Y0,Y1,V,ll)
  polypart1 = unlist(polynomial.values(poly[delta+1],X1*beta))
  polypart0 = unlist(polynomial.values(poly[delta+1],X0*beta))
  
  BT[uu] = mean(polypart1 - sol$Wsol %*% polypart0) 
  Sj2 = apply(sol$Wsol,2,sum)^2
  VAR[uu] = mean((Y1-polypart1)^2) + (sum(1-d)/sum(d))*mean(Sj2*(Y0-polypart0)^2)
}

# plot(lambda,BT,type="line",ylim=c(-.1,.1))

plot(lambda,VAR,type="line")

