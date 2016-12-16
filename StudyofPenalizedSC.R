### Sparse synthetic control: studying the penalty
### Jeremy L Hour
### 11 juillet 2016
### Edited: 22 juillet 2016

setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")


rm(list=ls())
set.seed(14021989)

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
source("functions/TZero.R")
source("functions/synthObj.R")
source("functions/regsynth.R")
source("functions/regsynthpath.R")
source("functions/lambda.max.R")

### Generate data
data = matchDGP(n=200,p=5,Ry=.5,Rd=.2)
X = data$X
y = data$y
d = data$d
  
X0 = t(X[d==0,])
X1 = t(X[d==1,])
V = diag(ncol(X))
Y1 = y[d==1]

### Compute solution
lambda = seq(0,2,.01)
solpath = regsynthpath(X0,X1,y[d==0],y[d==1],V,lambda)

plot(lambda,solpath$ATT,col="steelblue",pch=19)
abline(h=0, lty=6)

### How to select lambda ?
# Let's take a OB estimate
dataOB = data.frame(y=y,X=X)
OBreg = nls(y ~ exp(X.1*b1 + X.2*b2 + X.3*b3 + X.4*b4 + X.5*b5),dataOB,start=list(b1=0,b2=0,b3=0,b4=0,b5=0),subset=(d==0))
#OBreg = nls(y ~ exp(X.1*b1 + X.2*b2 + X.3*b3 + X.4*b4 + X.5*b5 + X.6*b6 + X.7*b7 + X.8*b8 + X.9*b9 + X.10*b10),
#            dataOB,start=list(b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,b7=0,b8=0,b9=0,b10=0),subset=(d==0))
summary(OBreg)
tau_hat2 = Y1 - exp(X[d==1,]%*%coef(OBreg))

# Measure fit
RSS = mapply(function(x) mean((solpath$CATT[x,]-tau_hat2)^2),x=1:nrow(solpath$CATT))
plot(lambda,RSS,col="steelblue",pch=19)

### See how the bias evolves
# True bias
bias = vector(length=length(lambda))
for(k in 1:length(lambda)){
  bias[k] = mean(exp(-X[d==1,]%*%data$b) - solpath$Wsol[k,,]%*%exp(-X[d==0,]%*%data$b))
}

# Estimated bias
bias_hat = vector(length=length(lambda))
for(k in 1:length(lambda)){
  bias_hat[k] = mean(exp(X[d==1,]%*%coef(OBreg)) - solpath$Wsol[k,,]%*%exp(X[d==0,]%*%coef(OBreg)))
}

matplot(lambda,cbind(solpath$ATT,abs(bias),abs(bias_hat)),type="l")

# Estimated variance
r = matrix(nrow=length(lambda), ncol=nrow(X))
for(k in 1:length(lambda)){
  r[k,d==1] = -1
  r[k,d==0] = apply(solpath$Wsol[k,,],2,sum)^2
}

# now need to estimate conditional variance
# do this using NN matching
cvar = vector(length=nrow(X))
cvar[d==0] = (y[d==0] - exp(X[d==0,]%*%coef(OBreg)))^2
jNN = vector(length=ncol(X1))
for(k in 1:ncol(X1)){
  n0 = ncol(X0)
  M = X1[,k]
  Delta = matrix(t(M)%*%V%*%M, nrow=n0, ncol=1) - 2*t(X0)%*%V%*%M + diag(t(X0)%*%V%*%M)
  jNN[k] = which(Delta == min(Delta))
}
cvar[d==1] = (cvar[d==0])[jNN]
plot(lambda,r%*%cvar/sum(d),col="steelblue",pch=19)

## Final plot
matplot(lambda,cbind(solpath$ATT,abs(bias)^2,abs(bias_hat)^2,r%*%cvar/sum(d)),t="o",pch=20)

### Get optimal lambda value
lambda[which(abs(bias)==min(abs(bias)))]
lambda[which(abs(bias_hat)==min(abs(bias_hat)))]
lambda[which(r%*%cvar/sum(d)==min(r%*%cvar/sum(d)))]

### Number of neighbors used
index = solpath$Wsol!=0
nonz = apply(index,c(1,2),sum)
matplot(lambda,cbind(apply(nonz,1,mean),apply(nonz,1,min),apply(nonz,1,max)), type="o",pch=20)

###############################
###############################
### ANALYSIS OF THE PROGRAM ###
###############################
###############################

### Trade-off for a point in the convex hull
M = apply(X1,1,mean)
lambda = seq(0,4,.01)
Wsol = matrix(nrow=length(lambda), ncol=sum(1-d))
l1norm_hat = vector(length = length(lambda))
loss_hat = vector(length = length(lambda))
att = vector(length = length(lambda))

for(i in 1:length(lambda)){
  sol = wsoll1(X0,M,V,lambda[i])
  sol = TZero(sol)
  Wsol[i,] = sol
  att[i] <- mean(y[d==1]) - y[d==0]%*%sol
  func <- synthObj(sol,X0,M,V)
  l1norm_hat[i] = func$l1norm
  loss_hat[i] = func$loss
}

# ATT as function of penalty level
plot(lambda,att,col="steelblue",pch=19)
abline(h=0, lty=6)

# Weight as function of penalty level
matplot(lambda,Wsol, type="l",
        main="Sparse Synthetic Control Regularization path",
        xlab="Penalty level", ylab="weight", ylim=c(0,1))

# Trade-off illustrated
plot(l1norm_hat,loss_hat,
     main="Fit-Sparsity Trade-off",
     xlab="l1-norm", ylab="loss",
     col="steelblue", pch=20)
# The value on the right is the pure synthetic control

### Bound on sum of distance
Delta = synthObj(sol,X0,M,V)$Delta
jnn = which(Delta == min(Delta))
sumofd_bound = Delta[jnn]*(1+lambda)/lambda
sumofd <- l1norm_hat

matplot(lambda,cbind(sumofd,sumofd_bound,Delta[jnn]), type="l",
        main="Bound on sum of distances",
        xlab="Penalty level", ylab="Sum of distances and bounds",
        xlim=c(0,.5), ylim=c(0,5))

### bound for the overall function
f = loss_hat + lambda*l1norm_hat
matplot(lambda,cbind(f,Delta[jnn]*(1+lambda),lambda^.9), type="l",
        main="Bound on sum of obj func",
        xlab="Penalty level", ylab="func",
        xlim=c(0,1),ylim=c(0,.5)) 

matplot(lambda,cbind(f,Delta[jnn]*(1+lambda),lambda^.9), type="l",
        main="Bound on sum of obj func",
        xlab="Penalty level", ylab="func")

### To be modified:
### Bound for the treated-synthetic discrepency
dp2 = loss_hat[200] + lambda*(l1norm_hat[200] - Delta[jnn])
matplot(lambda,cbind(loss_hat,dp2), type="l",
        main="Bound on distance to synthetic unit",
        xlab="Penalty level", ylab="Actual distance and bound")

#######################
#######################
### FIND lambda.max ###
#######################
#######################

### Last edited = 3/11/2016

### Generate data
set.seed(21101980)
data = matchDGP(n=50,p=5,Ry=.5,Rd=.2)
X = data$X; y = data$y; d = data$d

X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))

###
lambda = seq(0,3,.01)
solpath = regsynthpath(X0,t(X[d==1,]),y[d==0],y[d==1],V,lambda)

# Weight as function of penalty level
matplot(lambda,solpath$Wsol[,1,], type="l", lwd=2,
        main="Sparse Synthetic Control Regularization path",
        xlab="Penalty level", ylab="weight", ylim=c(0,1))
solpath$Wsol[which(lambda==2.50),1,]
# lambda.max seems before 2.5

M = X1[,1]
l.max = lambda.max(X0,M,V)
# but here given by 1.16
solpath$Wsol[which(lambda==1.16),1,]

### Which solution gives the best objective function with lambda.max?
Delta = matrix(t(M)%*%V%*%M, nrow=ncol(X0), ncol=1) - 2*t(X0)%*%V%*%M + diag(t(X0)%*%V%*%X0)
P = 2*t(X0)%*%V%*%X0
q = t(-2*t(X0)%*%V%*%M + l.max*Delta)

# Case 1
w1 = solpath$Wsol[which(lambda==1.16),1,]
print(t(w1)%*%P%*%w1/2 + q%*%w1)

# Case 2
w2 = rep(0, ncol(X0))
w2[which(Delta==min(Delta))] = 1
print(t(w2)%*%P%*%w2/2 + q%*%w2)


### Other example
data = matchDGP(n=10,p=50,Ry=.5,Rd=.2)
X = data$X; y = data$y; d = data$d

X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))

###
lambda = seq(0,15,.1)
solpath = regsynthpath(X0,t(X[d==1,]),y[d==0],y[d==1],V,lambda)

# Weight as function of penalty level
matplot(lambda,solpath$Wsol[,1,], type="l", lwd=2,
        main="Sparse Synthetic Control Regularization path",
        xlab="Penalty level", ylab="weight", ylim=c(0,1))
solpath$Wsol[which(lambda==11),1,]
# lambda.max seems before 11

M = X1[,1]
l.max = lambda.max(X0,M,V)
solpath$Wsol[which(lambda==5),1,]

### Which solution gives the best objective function with lambda.max?
Delta = matrix(t(M)%*%V%*%M, nrow=ncol(X0), ncol=1) - 2*t(X0)%*%V%*%M + diag(t(X0)%*%V%*%X0)
P = 2*t(X0)%*%V%*%X0
q = t(-2*t(X0)%*%V%*%M + l.max*Delta)

# Case 1
w1 = solpath$Wsol[which(lambda==5),1,]
print(t(w1)%*%P%*%w1/2 + q%*%w1)

# Case 2
w2 = rep(0, ncol(X0))
w2[which(Delta==min(Delta))] = 1
print(t(w2)%*%P%*%w2/2 + q%*%w2)
