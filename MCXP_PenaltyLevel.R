### Monte Carlo with selected penalty level
### Jeremy L Hour
### 29 aout 2016

setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")

rm(list=ls())
set.seed(12071990)

### 0. Settings

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
source("functions/regsynth.R")
source("functions/regsynthpath.R")
source("functions/TZero.R")
source("functions/synthObj.R")


### MC XP
lambda = seq(0,2,.01)
K = 5 # number of folds for optimal penalty level
R <- 100
Results <- matrix(ncol=5, nrow=R)
t_start <- Sys.time()
pb <- txtProgressBar(style = 3)

for(r in 1:R){
  ### 1. Generate data
  data = matchDGP(n=50,p=5,Ry=.5,Rd=.2)
  X = data$X
  y = data$y
  d = data$d
  
  X0 = t(X[d==0,])
  X1 = t(X[d==1,])
  V = diag(ncol(X))
  Y0 = y[d==0]
  n0 = sum(1-d)
  
  ### 2. SC on the mean
  M = matrix(apply(X1,1,mean), ncol=1)
  sol_mSC = wsol(X0,M,V)
  
  ### 3. 1NN matching 
  cf = vector(length = sum(d))
  for(i in 1:sum(d)){
    sol = matching(X0,X1[,i],V,m=1)
    cf[i] = t(y[d==0])%*%sol
  }
  attm1 <- mean(y[d==1] - cf)
  
  ### 4. 5NN matching
  cf = vector(length = sum(d))
  for(i in 1:sum(d)){
    sol = matching(X0,X1[,i],V,m=5)
    cf[i] = t(y[d==0])%*%sol
  }
  attm5 <- mean(y[d==1] - cf)
  
  ### 5. Regularized SC, fixed lambda
  cf = vector(length = sum(d))
  for(i in 1:sum(d)){
    sol = wsoll1(X0,X1[,i],V,.1)
    cf[i] = t(y[d==0])%*%sol
  }
  RegSC_fixed <- mean(y[d==1] - cf)
  
  ### 6. Regularized SC, optimized lambda
  uu=0
  while(uu==0){
    allocation = sample(1:K,n0,replace=T)
    uu=min(mapply(function(x) sum(allocation==x),1:K))
  }
  
  ATTcv = matrix(nrow=K, ncol=length(lambda))
  for(k in 1:K){
    X1k = as.matrix(X0[,allocation==k])
    X0k = as.matrix(X0[,allocation!=k])
    Y1k = Y0[allocation==k]
    Y0k = Y0[allocation!=k]
    solpath = regsynthpath(X0k,X1k,Y0k,Y1k,V,lambda)
    ATTcv[k,] = apply(solpath$CATT^2,1,sum)
  }
  
  curve = apply(ATTcv,2,sum)/n0
  lambda.opt = lambda[which(curve==min(curve))]
  
  cf = vector(length = sum(d))
  for(i in 1:sum(d)){
    sol = wsoll1(X0,X1[,i],V,lambda.opt)
    cf[i] = t(y[d==0])%*%sol
  }
  RegSC_opt <- mean(y[d==1] - cf)
  
  
  ### 7. Third step: ATT estimation
  Results[r,] <- c(wATT(y,d,sol_mSC),
                   attm1,
                   attm5,
                   RegSC_fixed,
                   RegSC_opt)
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
get.plot(data_res,2,"1NN Matching", msd)
get.plot(data_res,3,"5NN Matching", msd)
get.plot(data_res,4,"Reg SC fixed lambda", msd)
get.plot(data_res,5,"Reg SC opt lambda", msd)

### Compute bias and RMSE
StatDisplay <- data.frame()
StatDisplay[1:5,"bias"] <- apply(Results,2,mean)
StatDisplay[1:5,"RMSE"]  <- sqrt(apply(Results^2,2,mean))
StatDisplay[1:5,"ShapiroTest"]  <- apply(Results,2, function(x) shapiro.test(x)$p.value)
row.names(StatDisplay) <- c("AggregateSC","1nnMatching","5nnMatching","RegSCfixed","RegSCopt")
print(StatDisplay)