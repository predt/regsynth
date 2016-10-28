### Fisher Testing and Confidence Intervals

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

### DGP

FisherDGP <- function(n=100,tau=0){
  d = runif(n) > .5
  y = rnorm(n) + tau*d
  
  return(list(d=d,y=y))
}

data = FisherDGP(1000,tau=1)