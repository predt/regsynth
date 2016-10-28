### Matching and synthetic control : Monte Carlo experiment
### Jeremy L Hour
### 11 juillet 2016

setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")

rm(list=ls())
set.seed(3101990)

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


### MC XP
R <- 10000
Results <- matrix(ncol=6, nrow=R)
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
  
  ### 2. SC on the mean
  M = matrix(apply(X1,1,mean), ncol=1)
  sol_mSC = wsol(X0,M,V)
  
  ### 2. SC on each unit
  cf = vector(length = sum(d))
  for(i in 1:sum(d)){
    sol = wsol(X0,X1[,i],V)
    cf[i] = t(y[d==0])%*%sol
  }
  attSCu <- mean(y[d==1] - cf)
  
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
  
  ### 5. Sparse Synthetic Control, on mean
  sol_SparsemSC = wsoll1(X0,M,V,.1) # Penlaty originally set to .1
  
  ### 6. Sparse Synthetic Control, on each unit
  cf = vector(length = sum(d))
  for(i in 1:sum(d)){
    sol = wsoll1(X0,X1[,i],V,.1)
    cf[i] = t(y[d==0])%*%sol
  }
  attSparseSCu <- mean(y[d==1] - cf)
  
  
  ### 6. Third step: ATT estimation
  Results[r,] <- c(wATT(y,d,sol_mSC),
                   attSCu,
                   attm1,
                   attm5,
                   wATT(y,d,sol_SparsemSC),
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