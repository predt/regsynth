### Penalized Synthetic Control
### Inference and permutation tests
### Jeremy L Hour
### 8 septembre 2016

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

p.val = sum(ptest$theta.hat < ptest$theta.reshuffled)/1000
titer = data.frame(val=ptest$theta.reshuffled)

ggplot(titer, aes(x=val)) + 
    geom_histogram(binwidth = .06, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
    scale_x_continuous(name="Estimated ATET") +
    ggtitle("Distribution of permutated ATETs") + 
    geom_vline(xintercept = ptest$theta.hat, colour="darkorchid3", size=1) +
    theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
  

### 4. A Monte Carlo experiment
set.seed(12071990)
R = 5000
Results <- matrix(ncol=1, nrow=R)
t_start <- Sys.time()
pb <- txtProgressBar(style = 3)

for(r in 1:R){
  ### 0. Generate data
  data = matchDGP(n=50,p=5,Ry=.5,Rd=.2)
  X = data$X; y = data$y; d = data$d; V = diag(ncol(X))
  
  ### Perform test
  ptest = perm.test(d,y,X,V,.1,R=1000)
  
  ### 6. Third step: ATT estimation
  Results[r,] <- c(ptest$p.val)
  setTxtProgressBar(pb, r/R)
}

close(pb)
print(Sys.time()-t_start)

sum(Results < .05)/R