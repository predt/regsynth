### Fisher Exact Testing
### Jeremy L Hour
### 10 octobre 2016

### Very simple example of Fisher Tests for regression
### NB: not directly related to the project


setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")
rm(list=ls())

### 0. Settings

### Load packages
library("MASS")
library("ggplot2")


### Monte Carlo Experiment

### 1. DGP
n = 500
p = 5
X = matrix(rnorm(2*n*p), ncol=p)
d = c(rep(0,n),rep(1,n))
mu = 1
beta = c(1,-.8,.6,-.4,.2)
y = d*mu + X%*%beta +  rnorm(2*n, sd=1.5)

### 2. Original sample stats
lm.orig = lm(y ~ d + X)
summary(lm.orig)
theta.hat = coef(lm.orig)['d']
print(theta.hat)

### 3. Permutation function
permutation.iter = function(d,y,X){
  dstar = sample(d)
  lm.b = lm(y ~ dstar + X)
  return(coef(lm.b)['dstar'])
}

permutation.iter(d,y,X)

### 4. Fisher Testing
B = 1000
theta.reshuffled = replicate(B, permutation.iter(d,y,X), simplify="vector")
p.val = sum(theta.hat < theta.reshuffled)/B

titer = data.frame(val=theta.reshuffled)

ggplot(titer, aes(x=val)) + 
  geom_histogram(binwidth = .06, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(name="Estimated ATET") +
  ggtitle("Distribution of permutated ATETs") + 
  geom_vline(xintercept = theta.hat, colour="darkorchid3", size=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")