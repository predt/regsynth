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
source("functions/matchest.R")
source("functions/regsynth.R")
source("functions/regsynthpath.R")
source("functions/TZero.R")
source("functions/synthObj.R")


### MC XP
lambda = seq(0,2,.01)
K = 2 # number of folds for optimal penalty level
R = 1000
Results <- matrix(ncol=6, nrow=R)
t_start <- Sys.time()
pb <- txtProgressBar(style = 3)

for(r in 1:R){
  ### 0. Generate data
  data = matchDGP(n=50,p=5,Ry=.5,Rd=.2)
  X = data$X; y = data$y; d = data$d
  
  X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))
  Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)
  
  ### 1. Synthetic Control on mean of treated
  M = matrix(apply(X1,1,mean), ncol=1)
  AggSC = wATT(y,d,wsol(X0,M,V))
  
  ### 2. 1NN matching 
  NN1 = matchest(d,X,y,V,1)
  
  ### 3. 5NN matching
  NN5 = matchest(d,X,y,V,5)
  
  ### 4. Regularized Synthetic Control, fixed lambda
  sol = regsynth(X0,X1,Y0,Y1,V,.1)
  RSC.fixed = sol$ATT
  
  ### 5. Regularized SC, optimized lambda
  uu=0 # reshuffle groups until no empty group
  while(uu==0){
    allocation = sample(1:K,n0,replace=T)
    uu=min(mapply(function(x) sum(allocation==x),1:K))
  }
  
  print("*** PROGRESS ***")
  print(100*r/R)
  
  RMSEcv = matrix(nrow=K, ncol=length(lambda))
  biascv = matrix(nrow=K, ncol=length(lambda))
  for(k in 1:K){
    X1k = as.matrix(X0[,allocation==k])
    X0k = as.matrix(X0[,allocation!=k])
    Y1k = Y0[allocation==k]
    Y0k = Y0[allocation!=k]
    solpath = regsynthpath(X0k,X1k,Y0k,Y1k,V,lambda)
    RMSEcv[k,] = apply(solpath$CATT^2,1,sum)
    biascv[k,] = apply(solpath$CATT,1,sum)
  }
  
  # The one that optimizes RMSE
  curve.RMSE = apply(RMSEcv,2,sum)/n0
  lambda.opt.RMSE = lambda[which(curve.RMSE==min(curve.RMSE))]
  sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.RMSE)
  RSC.opt.RMSE = sol$ATT
  
  # The one that optimizes bias
  curve.bias = abs(apply(biascv,2,sum)/n0)
  lambda.opt.bias = lambda[which(curve.bias==min(curve.bias))]
  sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.bias)
  RSC.opt.bias = sol$ATT
  
  
  ### 6. Third step: ATT estimation
  Results[r,] <- c(AggSC,
                   NN1$ATT,
                   NN5$ATT,
                   RSC.fixed,
                   RSC.opt.RMSE,
                   RSC.opt.bias)
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
get.plot(data_res,5,"Reg SC opt lambda RMSE", msd)
get.plot(data_res,6,"Reg SC opt lambda bias", msd)


### Compute bias and RMSE
StatDisplay <- data.frame()
StatDisplay[1:6,"bias"] <- apply(Results,2,mean)
StatDisplay[1:6,"RMSE"]  <- sqrt(apply(Results^2,2,mean))
StatDisplay[1:6,"ShapiroTest"]  <- apply(Results,2, function(x) shapiro.test(x)$p.value)
row.names(StatDisplay) <- c("AggregateSC","1nnMatching","5nnMatching","RegSCfixed","RegSCopt","RegSCoptbias")
print(StatDisplay)