### Exp 5 : performance of different procedures
### comparision with matching estimators
### NEW DGP
### Jeremy L Hour
### 12 octobre 2016

setwd("/Users/jeremylhour/Documents/R/regsynth-master")
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
source("functions/PolyDGP.R")
source("functions/wATT.R")
source("functions/matching.R")
source("functions/matchest.R")
source("functions/regsynth.R")
source("functions/regsynthpath.R")
source("functions/TZero.R")
source("functions/synthObj.R")
source("simulations/Exp5_PolyDGP/Exp5_Poly_setup.R")


### MC XP
set.seed(12071990)
lambda = seq(0,5,.005) # set of lambda to be considered for optim
K = 5 # number of folds for optimal penalty level

Exp5_Poly_setup(R=10,n1=10,n0=100,p=2,delta=2)

for(n_xp in c(2,25,50,100)){
  for(p_xp in c(1,4,8)){
    for(delta_xp in c(1,2,4)){
      Exp5_Poly_setup(R=10,n1=n_xp,n0=50,p=p_xp,delta=delta_xp)
    }
  }
}


# Problemes: (delta=2)
# bias toujours positif (MAE =  bias)
# performance du PenSynth bien inférieure
# Soit un probleme de support commun ? Pb du grid search (pass assez de lambda?) ?
# R2 is .5 for the treated but not for controls. Do the think to see where lambda CV converges
# PenSynth works when not so many controls