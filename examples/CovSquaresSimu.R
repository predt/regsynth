### Simulations on Cov(S^2,S^2)
### Jeremy L'Hour
### 09/07/2019


setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")
rm(list=ls())

### 0. Settings

### Load packages
library("MASS")
library("LowRankQP")
library("xtable")

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

set.seed(12071990)

### 0. Parameters
R = 10000 # nb sim

n1_set = c(1,4,5,10,20)
n0_set = c(5,10,15,20,50) # all elements must be larger than 2
p_set = c(2,5,8,10)


### 1. Simulations

Results = matrix(,nrow=length(n0_set)*length(p_set),ncol=length(n1_set))
Results_T = matrix(,nrow=length(n0_set)*length(p_set),ncol=length(n1_set))
t_start = Sys.time()
pb = txtProgressBar(style = 3)
i=0; c=0;

for(p in p_set){
  i=i+1; k=0;
  for(n0 in n0_set){
    j=0; k=k+1
    for(n1 in n1_set){
      j=j+1
      c=c+1
      simu = matrix(,nrow=R,ncol=n0)
      for(r in 1:R){
        # Simulate covariate
        X1 = matrix(runif(n1*p,min=.1,max=.9), nrow=p,ncol=n1)
        X0 = matrix(sqrt(runif(n0*p)), nrow=p,ncol=n0)
        Y0 = rep(0,n0); Y1 = rep(1,n1)
        
        # Compute synthetic weights
        sol = regsynth(X0,X1,Y0,Y1,V=diag(p),.1)
        simu[r,] = apply(sol$Wsol,2,sum)
      }
      covariances = cov(simu^2, method="pearson")
      
      # Covariance
      Results[(i-1)*length(n0_set)+k,j] = covariances[1,2]
      
      # T-stat from linear model
      reg = lm(simu[,1]^2 ~ simu[,2]^2)
      Results_T[(i-1)*length(n0_set)+k,j] = summary(reg)[["coefficients"]][2, "t value"] # rejet si plus grand que 1.64
      
      setTxtProgressBar(pb, c/(length(n1_set)*length(n0_set)*length(p_set)))
    }
  }
}

close(pb)
print(Sys.time()-t_start)

### 2. Mise en forme
j=0
bracket <- function(x, i, j){
  x[i,j] <- sprintf("(%s)", x[i,j])
  x
}

for(p in p_set){
  table = Results[(j*length(n0_set)+1):((j+1)*length(n0_set)),]
  table_T = Results_T[(j*length(n0_set)+1):((j+1)*length(n0_set)),]
  
  Sup_Table = matrix(ncol=ncol(table),nrow=2*nrow(table))
  
  for(k in 1:nrow(table)){
    Sup_Table[2*(k-1)+1,] = table[k,]
    Sup_Table[2*k,] = table_T[k,]
  }
  rownames(Sup_Table) = c(rbind(n0_set,n0_set))
  colnames(Sup_Table) <- n1_set
  Sup_Table = round(Sup_Table, digits=3)
  
  ind = which(row(Sup_Table) %% 2 == 0, arr.ind = TRUE)
  Sup_Table <- bracket(Sup_Table, ind[,1], ind[,2])

  TeX_Table = xtable(Sup_Table,digits=3, caption=paste("Covariance of sum of weights squared, p=",p,". T-statistics of nullity test between brackets. ", R," simulations.",sep=""))
  
  print(TeX_Table, include.rownames=T,file=paste("CovSquared_p",p,".txt",sep=""))
  j=j+1
}