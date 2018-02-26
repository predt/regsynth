### Run MC_XP, with different Legendre polynomial order
### Jeremy L Hour
### 24 juillet 2017

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
library("xtable")
library("orthopolynom")

### Load user functions
source("functions/wsol.R")
source("functions/wsoll1.R")
source("functions/wATT.R")
source("functions/matching.R")
source("functions/matchest.R")
source("functions/OBest.R")
source("functions/regsynth.R")
source("functions/regsynthpath.R")
source("functions/TZero.R")
source("functions/synthObj.R")
source("functions/NLmatchDGP.R")
source("simulations/NonLinearMCXP_setup.R")


### MC XP
set.seed(12071990)
lambda = seq(0,5,.01) # set of lambda to be considered for optim
K = 5 # number of folds for optimal penalty level
n_xp = 50 # number of individuals
p_xp = 10 # number of covariates

for(deg_lin in seq(0,6,by=1)){
    NonLinearMCXP_setup(R=1000,n=n_xp,p=p_xp,delta=deg_lin)
}