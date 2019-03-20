### Distribution of S_j
### Jeremy L Hour
### 23/08/2017

setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")
rm(list=ls())

### 0. Settings

### Load packages
library("MASS")
library("gtools")
library("ggplot2")
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


### 0. Parameters
R = 1000 # nb simulations
n1 = 10 # nb treated
n0 = 20 # nb controls
p = 20

### DGP Setup
### Covariate variance matrix
Sigma = matrix(0,nrow=p, ncol=p)
for(k in 1:p){
  for(j in 1:p){
    Sigma[k,j] = .5^abs(k-j)
  }
}

b = rep(0,p)
for(j in 1:p){
  b[j] = 1*(-1)^(j) / j^2
}

set.seed(12071990)

### 1. Simulations

Results = matrix(,nrow=R,ncol=n0)
t_start = Sys.time()
pb = txtProgressBar(style = 3)

for(r in 1:R){
  # Simulate covariate
  
  # SETUP 1 
  #X0 = matrix(2*sqrt(runif(n0))-1, nrow=1)
  #X1 = matrix(runif(n1,min=-1,max=1), nrow=1)
  #Y0 = rep(0,n0); Y1 = rep(1,n1)
  
  # SETUP 2
  #X = mvrnorm(n = n0+n1, mu=rep(0,p), Sigma)
  #X1 = t(X[1:n1,]); X0 = t(X[(n1+1):(n1+n0),]);
  #y = exp(-X%*%b) + rnorm(n0+n1)
  #Y1 = y[1:n1]; Y0 = y[(n1+1):(n1+n0)]
  
  # SETUP 3: Case of small common support
  #X0 = matrix(rbeta(n0,2,5), nrow=1)
  #X1 = matrix(rbeta(n1,5,1), nrow=1)
  #Y0 = rep(0,n0); Y1 = rep(1,n1)
  
  # SETUP 3: Case of small common support
  X0 = matrix(rbeta(n0,2,5), nrow=1)
  X1 = matrix(rbeta(n1,2,5), nrow=1)
  Y0 = rep(0,n0); Y1 = rep(1,n1)
  
  # Compute synthetic weights
  sol = regsynth(X0,X1,Y0,Y1,V=1,.1)
  Results[r,] = apply(sol$Wsol,2,sum)
  
  # Progress bar
  setTxtProgressBar(pb, r/R)
}

close(pb)
print(Sys.time()-t_start)


### 2. Features of the distribution

### Second moment
print("Empirical Second moments:")
print(apply(Results^2,2,mean))

print("Upper bound:")
print(n1^2/n0)

### Variance
print("Empirical variance:")
print(apply(Results,2,var))
print("Upper bound:")
print( (n1^2/n0)*((n0-1)/n0) )

### Variance of the square
print("Empirical variance of the square:")
print(apply(Results^2,2,var))

### Correlation of the squares
print("Empirical correlations of squares:")
cor(Results^2)

### Probability that S_1 is an integer
cumuu = vector(length=n0)
for(i in 0:n1){
  cumuu = cumuu + mapply(function(k) mean(Results[,k]==i),1:n0)
  print(paste("Probability that S_1 =",i))
  print(mapply(function(k) mean(Results[,k]==i),1:n0))
}

print(paste("Probability that S_1 is an integer"))
print(cumuu)

### Distribution of the max
data = data.frame("max"=apply(Results,1,max))
summary(data)

ggplot(data, aes(x=max)) + 
  geom_histogram(binwidth = .2, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(name="Max of S_j") +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")

print(quantile(data[,"max"], probs = c(.95,.975,.99)))

### Distribution of an S_j
data = data.frame("S_1"=Results[,6])
summary(data)

pdf("plot/DistributionSj.pdf", width=6, height=6)
ggplot(data, aes(x=S_1)) + 
  geom_histogram(binwidth = .2, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(name="Distribution of S_1", limits=c(0,n1)) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
dev.off()

print(quantile(data[,"S_1"], probs = c(.95,.975,.99)))


### Binomial distribution: moments
n1 = 20
n0 = 500
simBinom = rbinom(R,n1,1/n0)
  
m1 = n1/n0
print(paste("Theory:",m1," - Empirical value:", mean(simBinom)))
  
m2 = (n1/n0)*((n0+n1-1)/n0)
print(paste("Theory:",m2," - Empirical value:", mean(simBinom^2)))
  
m3 = n1*((n1-1)*(3*n0+n1-2)+n0^2)/(n0^3)
print(paste("Theory:",m3," - Empirical value:", mean(simBinom^3)))
  
m4 = (n1/n0)*(1 + ((n1-1)/n0) * ( ((n1-2)/n0)*(6+(n1-3)/n0) +7) )
print(paste("Theory:",m4," - Empirical value:", mean(simBinom^4)))
  
print(paste("Theory:",m4-m2^2," - Empirical value:", var(simBinom^2)))
