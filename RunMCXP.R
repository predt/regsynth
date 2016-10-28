### Run MC_XP
### Jeremy L Hour
### 12 octobre 2016

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
source("simulations/MCXP_setup.R")


### MC XP
set.seed(12071990)
lambda = seq(0,2,.01) # set of lambda to be considered for optim
K = 5 # number of folds for optimal penalty level

for(n_xp in c(30,50,100)){
  for(p_xp in c(3,10,20)){
    MCXP_setup(R=10000,n=n_xp,p=p_xp)
  }
}
