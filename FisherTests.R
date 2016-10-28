### Penalized Synthetic Control
### Inference and permutation tests
### Jeremy L Hour
### 8 septembre 2016
### EDIT: 21 octobre 2016

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
source("functions/perm.test.R")

### 0. Generate data
data = matchDGP(n=100,p=10,Ry=.5,Rd=.2,a=1)
X = data$X; y = data$y; d = data$d

X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))
Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)


### 2. Compute ATET on the original sample
sol1 = regsynth(X0,X1,y[d==0],y[d==1],V,.1)
print(sol1$ATT)

### 3. Inference based on Permutation Tests
# Reshuffle the treatment
ptest = perm.test(d,y,X,V,.1,R=1000)

titer = data.frame(val=ptest$theta.reshuffled)

ggplot(titer, aes(x=val)) + 
    geom_histogram(binwidth = .06, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
    scale_x_continuous(name="Estimated ATET") +
    ggtitle("Distribution of permutated ATETs") + 
    geom_vline(xintercept = ptest$theta.hat, colour="darkorchid3", size=1) +
    theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
  

### 4. A Monte Carlo experiment

# Set the function
InferenceMCXP <- function(R=1000,B=100000,n=100,p=50){
  Results = vector(length=R)
  t_start = Sys.time()
  pb = txtProgressBar(style = 3)
  
  for(r in 1:R){
    ### 0. Generate data
    data = matchDGP(n=n,p=p,Ry=.5,Rd=.2)
    X = data$X; y = data$y; d = data$d; V = diag(ncol(X))
    
    ### Perform test
    ptest = perm.test(d,y,X,V,.1,R=B)
    
    Results[r] = c(ptest$p.val)
    setTxtProgressBar(pb, r/R)
  }
  
  close(pb)
  print(Sys.time()-t_start)

  # Screen print
  print(paste("Rejection rate at .1 pct:",mean(Results < .1)))
  print(paste("Rejection rate at .05 pct:",mean(Results < .05)))
  print(paste("Rejection rate at .01 pct:",mean(Results < .01)))
  
  fileN = paste("simulations/TestOutput_n",n,",p",p,".txt",sep="")
  
  write(c(paste("Nb. observations:",n),
          paste("Nb. covariates:",p),
          paste("Nb. replications:",R),
          paste("Nb. draws:",B),
          paste("Rejection rate .1:",mean(Results < .1)),
          paste("Rejection rate .05:",mean(Results < .05)),
          paste("Rejection rate .01:",mean(Results < .01)),
          paste(Sys.time())), fileN, append=TRUE)
}

set.seed(12071990)
for(n_xp in c(30,50,100)){
  for(p_xp in c(3,10,20)){
    InferenceMCXP(R=1000,B=1000,n=n_xp,p=p_xp)
  }
}
