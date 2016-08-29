### Synthetic control, matching and bias
### VERY INTUITIVE EXAMPLE
### Jeremy L Hour
### 8 july 2016

# setwd("/Volumes/USB_KEY/New_project/SyntheticControlMatching/code") # mac
setwd("E:/New_project/SyntheticControlMatching/code")

rm(list=ls())
set.seed(1214021989)

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

### 1. A very simple example
n = 20
x = matrix(rexp(2*n), ncol=2)
g = c(.5,-1)
logit = function(x) 1/(1+exp(-x))
d = as.numeric(runif(n) < logit(x %*%g))
d = as.factor(d)
levels(d)[levels(d)=="1"] = "Treated"
levels(d)[levels(d)=="0"] ="Control"
y = 10*exp(-diag(x%*%t(x)))
data = data.frame(X=x,d=d)
ATT = data.frame()


### 0. Display the dots
X0 = t(x[d=="Control",])
X1 = t(x[d=="Treated",])
X1 = matrix(apply(X1,1,mean), ncol=1)
V = diag(2)

### plot the average of the treated
data2 <- data
data2$d <- as.character(data2$d)
data2 <- rbind(data2,c(X1,"Aggregate Mean"))
data2$d <- factor(data2$d, levels=c("Treated", "Control", "Aggregate Mean"))
data2$X.1 <- as.numeric(data2$X.1)
data2$X.2 <- as.numeric(data2$X.2)

# The blue point is the average of the treated
pdf("plot/datavis.pdf")
ggplot(data=data2, aes(x=X.1, y=X.2, color=d)) +
  geom_point(size=6) +
  scale_x_continuous(name="X1") +
  scale_y_continuous(name="X2") +
  ggtitle("X distribution for Treated / Control group") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  scale_fill_discrete(name="Group",
                      labels=c("COntrol","Treated","Aggregate Mean")) +
  theme(legend.position="bottom") 
dev.off()

# Naive estimator
ATT[1,"TE"] = mean(y[d=="Treated"]) - mean(y[d=="Control"])

### 1. Run Synthetic Control on Aggregate unit
sol = wsol(X0,X1,diag(2),0)

data2$d <- as.character(data2$d)
data3 <- rbind(data2,c(X0%*%sol,"Synthetic Unit"))
data3$d <- factor(data3$d, levels=c("Treated", "Control", "Aggregate Mean","Synthetic Unit"))
data3$X.1 <- as.numeric(data3$X.1)
data3$X.2 <- as.numeric(data3$X.2)

# The purple point is the synthetic control,
# It is exactly matched
pdf("plot/datavis2.pdf")
ggplot(data=data3, aes(x=X.1, y=X.2, color=d)) +
  geom_point(size=6) +
  scale_x_continuous(name="X1") +
  scale_y_continuous(name="X2") +
  ggtitle("X distribution for Treated / Control group") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  scale_fill_discrete(name="Group",
                      labels=c("COntrol","Treated","Aggregate Mean")) +
  theme(legend.position="bottom") 
dev.off()

### My guess here:
# the outcome function of the average is not the same as the average outcome function
# so there might be a big bias by aggregating if the outcome equation is not linear

# Let's see:
ATT[2,"TE"] = mean(y[d=="Treated"]) - t(y[d=="Control"])%*%sol

### 2. Synthetic Control on each treated unit
X1 = t(x[d=="Treated",])
synthetic = matrix(,nrow=sum(d=="Treated"),ncol=3)
cf = vector(length = sum(d=="Treated"))
for(i in 1:sum(d=="Treated")){
  sol = wsol(X0,X1[,i],diag(2),0)
  cf[i] = t(y[d=="Control"])%*%sol
  synthetic[i,] = c(X0%*%sol,"Synthetic Unit")
}

### Now let's match all the synthetic control units
matdata = as.matrix(data2)
matdata = rbind(matdata,synthetic)
data3 = data.frame(X.1=as.numeric(matdata[,1]),
                   X.2=as.numeric(matdata[,2]),
                   d=factor(matdata[,3], levels=c("Treated", "Control", "Aggregate Mean","Synthetic Unit")))


pdf("plot/datavis3.pdf")
ggplot(data=data3, aes(x=X.1, y=X.2, color=d)) +
  geom_point(size=6) +
  scale_x_continuous(name="X1") +
  scale_y_continuous(name="X2") +
  ggtitle("X distribution for Treated / Control group") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  scale_fill_discrete(name="Group",
                      labels=c("COntrol","Treated","Aggregate Mean")) +
  theme(legend.position="bottom")
dev.off()
### The  points outside the convex hull are reproduced with neighbors
### The points in the convex hull are reproduced with anything that goes !!!
### Big bias to be expected ?
ATT[3,"TE"] = mean(y[d=="Treated"] - cf)
# -.40


### 3. Let's see if the penalty can reduce the problem
synthetic = matrix(,nrow=sum(d=="Treated"),ncol=3)
cf = vector(length = sum(d=="Treated"))
for(i in 1:sum(d=="Treated")){
  sol = wsoll1(X0,X1[,i],diag(2),4)
  cf[i] = t(y[d=="Control"])%*%sol
  synthetic[i,] = c(X0%*%sol,"Synthetic Unit")
}

matdata = as.matrix(data2)
matdata = rbind(matdata,synthetic)
data3 = data.frame(X.1=as.numeric(matdata[,1]),
                   X.2=as.numeric(matdata[,2]),
                   d=factor(matdata[,3], levels=c("Treated", "Control", "Aggregate Mean","Synthetic Unit")))

pdf("plot/datavis4.pdf")
ggplot(data=data3, aes(x=X.1, y=X.2, color=d)) +
  geom_point(size=6) +
  scale_x_continuous(name="X1") +
  scale_y_continuous(name="X2") +
  ggtitle("X distribution for Treated / Control group") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  scale_fill_discrete(name="Group",
                      labels=c("COntrol","Treated","Aggregate Mean")) +
  theme(legend.position="bottom")
dev.off()

# Bias is reduced 
ATT[4,"TE"] = mean(y[d=="Treated"] - cf)

### 4. 1-NN matching
matched = matrix(,nrow=sum(d=="Treated"),ncol=3)
cf = vector(length = sum(d=="Treated"))
for(i in 1:sum(d=="Treated")){
  M = X1[,i]
  dis = mapply(function(i) t(M-X0[,i]) %*%(M-X0[,i]),1:sum(d=="Control"))
  n <- length(dis)
  NN1 = which(dis==min(dis))
  Y0 = y[d=="Control"]
  cf[i] = mean(Y0[NN1])
  matched[i,] = c(X0[,NN1],"Matched Unit")
}

matdata = as.matrix(data2)
matdata = rbind(matdata,matched)
data3 = data.frame(X.1=as.numeric(matdata[,1]),
                   X.2=as.numeric(matdata[,2]),
                   d=factor(matdata[,3], levels=c("Treated", "Control", "Aggregate Mean","Matched Unit")))

pdf("plot/datavis5.pdf")
ggplot(data=data3, aes(x=X.1, y=X.2, color=d)) +
  geom_point(size=6) +
  scale_x_continuous(name="X1") +
  scale_y_continuous(name="X2") +
  ggtitle("X distribution for Treated / Control group") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  scale_fill_discrete(name="Group",
                      labels=c("COntrol","Treated","Aggregate Mean")) +
  theme(legend.position="bottom")
dev.off()

## Matching effect
ATT[5,"TE"] = mean(y[d=="Treated"] - cf)
# .20

### 5. How does 2-NN matching do ??
matched = matrix(,nrow=sum(d=="Treated"),ncol=3)
cf = vector(length = sum(d=="Treated"))
for(i in 1:sum(d=="Treated")){
  M = X1[,i]
  dis = mapply(function(i) t(M-X0[,i]) %*%(M-X0[,i]),1:sum(d=="Control"))
  n <- length(dis)
  NN1 = which(dis==min(dis))
  NN2 = which(dis==sort(dis,partial=2)[2])
  Y0 = y[d=="Control"]
  cf[i] = mean(Y0[c(NN1,NN2)])
  matched[i,] = c(apply(X0[,c(NN1,NN2)],1,mean),"Matched Unit")
}

matdata = as.matrix(data2)
matdata = rbind(matdata,matched)
data3 = data.frame(X.1=as.numeric(matdata[,1]),
                   X.2=as.numeric(matdata[,2]),
                   d=factor(matdata[,3], levels=c("Treated", "Control", "Aggregate Mean","Matched Unit")))

pdf("plot/datavis6.pdf")
ggplot(data=data3, aes(x=X.1, y=X.2, color=d)) +
  geom_point(size=6) +
  scale_x_continuous(name="X1") +
  scale_y_continuous(name="X2") +
  ggtitle("X distribution for Treated / Control group") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  scale_fill_discrete(name="Group",
                      labels=c("COntrol","Treated","Aggregate Mean")) +
  theme(legend.position="bottom")
dev.off()

## Matching effect
ATT[6,"TE"] = mean(y[d=="Treated"] - cf)
# .26



row.names(ATT) = c("Naive", "Aggregate SC","Individual SC", "Penalized SC","1NN Matching","2NN Matching")
print(ATT)