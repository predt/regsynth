### Lalonde dataset example
### Jeremy L Hour

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

### min-max scale
mMscale <- function(X){
  X = as.matrix(X)
  mins = apply(X,2,min)
  maxs = apply(X,2,max)
  return(scale(X, center=mins, scale=maxs-mins))
}

# Load data
library("causalsens")
data(lalonde.psid)

d = lalonde.psid[,"treat"]
y = lalonde.psid[,"re78"]

X = data.frame(lalonde.psid[,c("age","education","married","black","hispanic","re74","re75","nodegree")],
                     "NoIncome74"=as.numeric(lalonde.psid[,"re74"]==0),
                     "NoIncome75"=as.numeric(lalonde.psid[,"re75"]==0)
)
X[,c("age","education","re74","re75")] = mMscale(X[,c("age","education","re74","re75")])
X = as.matrix(X)


### Run Synthetic Control
X0 = t(X[d==0,])
X1 = t(X[d==1,])
V = diag(ncol(X))

set.seed(12071990)
lambda = seq(0,2,.01) # set of lambda to be considered for optim
K = 5 # number of folds for optimal penalty level

# A. lambda = .1
sol = regsynth(X0,X1,Y0,Y1,V,.1)
RSC.fixed = sol$ATT

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


### The dataset could be used to do some kind of Monte Carlo
R = 1000
n = 100
Results = matrix(nrow=R,ncol=6)
t_start <- Sys.time()
pb <- txtProgressBar(style = 3)

for(r in 1:R){
  draw = sample(1:nrow(X),n)
  Xr = X[draw,]
  dr = d[draw]
  yr = y[draw]
  X0 = t(Xr[dr==0,])
  X1 = t(Xr[dr==1,])
  V = diag(ncol(X))
  
  ### 2. SC on the mean
  M = matrix(apply(X1,1,mean), ncol=1)
  sol_mSC = wsol(X0,M,V)
  
  ### 2. SC on each unit
  cf = vector(length = sum(dr))
  for(i in 1:sum(dr)){
    sol = wsol(X0,X1[,i],V)
    cf[i] = t(yr[dr==0])%*%sol
  }
  attSCu <- mean(yr[dr==1] - cf)
  
  ### 3. 1NN matching 
  cf = vector(length = sum(dr))
  for(i in 1:sum(dr)){
    sol = matching(X0,X1[,i],V,m=1)
    cf[i] = t(yr[dr==0])%*%sol
  }
  attm1 <- mean(yr[dr==1] - cf)
  
  ### 4. 5NN matching
  cf = vector(length = sum(dr))
  for(i in 1:sum(dr)){
    sol = matching(X0,X1[,i],V,m=5)
    cf[i] = t(yr[dr==0])%*%sol
  }
  attm5 <- mean(yr[dr==1] - cf)
  
  ### 5. Sparse Synthetic Control, on mean
  sol_SparsemSC = wsoll1(X0,M,V,.1) # Penalty originally set to .1
  
  ### 6. Sparse Synthetic Control, on each unit
  cf = vector(length = sum(dr))
  for(i in 1:sum(dr)){
    sol = wsoll1(X0,X1[,i],V,.1)
    cf[i] = t(yr[dr==0])%*%sol
  }
  attSparseSCu <- mean(yr[dr==1] - cf)
  
  
  ### 6. Third step: ATT estimation
  Results[r,] <- c(wATT(yr,dr,sol_mSC),
                   attSCu,
                   attm1,
                   attm5,
                   wATT(yr,dr,sol_SparsemSC),
                   attSparseSCu)
  
  setTxtProgressBar(pb, r/R)
}

close(pb)
print(Sys.time()-t_start)

# Post-simulation treatment

# Draw the charts
id <- c(mapply(function(x) rep(x,R),1:6))
val <- c(Results)
data_res <- data.frame(val = val, model = id)

M <- max(abs(quantile(Results,.01)),abs(quantile(Results,.99)))
lb <- -1.1*M
ub <- 1.1*M
msd <- max(mapply(function(x)  sd(subset(data_res,model==x)[,1]),1:6))


### Function for plot
get.plot <- function(data,modelS,title="A Title",sdBCH){
  plot_res <- ggplot(subset(data, (model==modelS)), aes(x=val)) + 
    geom_histogram(binwidth = .02, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
    scale_x_continuous(limits=c(lb,ub), name="Treatment effect") +
    ggtitle(title) + 
    stat_function(fun = dnorm, args=list(mean=0, sd=sdBCH), colour="darkorchid3", size=1) +
    theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
  
  return(plot_res)
}

### Plots
get.plot(data_res,1,"Aggregate Synthetic Control", msd)
get.plot(data_res,2,"Individual Synthetic Control", msd)
get.plot(data_res,3,"1NN Matching", msd)
get.plot(data_res,4,"5NN Matching", msd)
get.plot(data_res,5,"Sparse Synthetic Control", msd)
get.plot(data_res,6,"Individual Sparse Synthetic Control", msd)

### Compute bias and RMSE
StatDisplay <- data.frame()
StatDisplay[1:6,"bias"] <- apply(Results,2,mean)
StatDisplay[1:6,"RMSE"]  <- sqrt(apply(Results^2,2,mean))
StatDisplay[1:6,"AsySD"]  <- apply(AsySD,2,mean)
StatDisplay[1:6,"ShapiroTest"]  <- apply(Results,2, function(x) shapiro.test(x)$p.value)
row.names(StatDisplay) <- c("AggregateSC","IndivSC","1nnMatching","5nnMatching","SparseAggregate","IndividualSparse")
print(StatDisplay)