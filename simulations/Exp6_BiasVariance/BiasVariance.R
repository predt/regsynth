### Bias and Variance: matching vs. penalized synthetic control
### Jeremy L'Hour
### 21/12/2018

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
source("functions/UniformDGP.R")


### I. MC XP: simple Alberto DGP
set.seed(2363978)
lambda = seq(.0001,3,by=.25) # set of lambda to be considered for optim
M = 10 # maximum number of neighbors
n1 = 12; n0 = 100;
a=.1; b=.9; h=.1

data = UniformDGP(n1,n0,p=5,a=.1,b=.9,h=.1)
X = data$X; d = data$d
X0 = t(as.matrix(X)[d==0,]); X1 = t(as.matrix(X)[d==1,]); V = diag(ncol(as.matrix(X)))
n0 = sum(1-d);

rset = seq(1,4,by=.2)
Mresults = matrix(nrow=length(rset),ncol=2)
Sresults = matrix(nrow=length(rset),ncol=2)

for(i in 1:length(rset)){
  # Regression function
  mu0 <- function(x,r=1,a,b){
    C = (b^(2*r+1) - a^(2*r+1))/((b-a)*(2*r+1)) - ((b^(r+1)-a^(r+1))/((b-a)*(r+1)))^2
    return(sum(x^r) / sqrt(length(x)*C))
  }
  
  # Compute Expectation
  mu = apply(X,1,mu0,r=rset[i],a=a,b=b)
  Y0 = mu[d==0]; Y1 = mu[d==1]; 
  
  # Matching
  NN = matchest(X0,X1,Y0,Y1,V,1)
  Mresults[i,1] = abs(mean(Y1 - NN$Wsol%*%Y0)) # bias
  Mresults[i,2] = mean(apply(NN$Wsol,2,sum)^2) # variance
  
  Ssol = regsynth(X0,X1,Y0,Y1,V,0.01)
  Sresults[i,1] = abs(Ssol$ATT) # Bias
  Sresults[i,2] = mean(apply(Ssol$Wsol,2,sum)^2) # Variance
}

plot(Sresults[,2],type="l")
lines(Mresults[,2])




#####

# Regression function
mu0 <- function(x,r=1,a,b){
  C = (b^(2*r+1) - a^(2*r+1))/((b-a)*(2*r+1)) - ((b^(r+1)-a^(r+1))/((b-a)*(r+1)))^2
  return(sum(x^r) / sqrt(length(x)*C))
}

# Compute Expectation
mu = apply(X,1,mu0,r=2,a=a,b=b)
Y0 = mu[d==0]; Y1 = mu[d==1]; 

# Matching
Mresults = matrix(nrow=M,ncol=2)
mbias = vector(length=M)
mvar = vector(length=M)
for(m in 1:M){
  NN = matchest(X0,X1,Y0,Y1,V,m)
  mbias[m] = mean(Y1 - NN$Wsol%*%Y0)
  mvar[m] = mean(apply(NN$Wsol,2,sum)^2)
}

Mresults[,1] = abs(mbias) # bias
Mresults[,2] = mvar # variance

# Penalized Synthetic Control
Sresults = matrix(nrow=length(lambda),ncol=2)
Ssol = regsynthpath(X0,X1,Y0,Y1,V,lambda)
Sresults[,1] = abs(Ssol$ATT) # Bias
Sresults[,2] = apply(apply(Ssol$Wsol,c(1,3),sum)^2,1,mean) # Variance

### bias plot
plot(Mresults[,1], type="l")
lines(rev(Sresults[,1]))

### variance plot
plot(Mresults[,1]^2+Mresults[,2], type="l", ylim=c(0,0.6))
lines(rev(Sresults[,1]^2+Sresults[,2]))




### Settup of the MC experiment
BiasVariance_XP <- function(R=1000,n1=10,n0=100,p=1,a=.1,b=.9,h=.1){
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
    sol = regsynthpath(X0,X1,Y0,Y1,V,lambda)
    
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

