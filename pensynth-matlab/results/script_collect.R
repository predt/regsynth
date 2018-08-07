### Script to combine results
### 23/07/2018
### Jérémy L'Hour

setwd("/Users/jeremylhour/Downloads/pensynth-matlab/results/")
rm(list=ls())

### 0. Settings

### Load packages
library("MASS")
library("xtable")

k_rng = c(2,5,10,15,30,50)
r_rng = c(1,1.2,1.5,1.8,2,2.2,2.5)

AA = read.table("n1_10_n0_100_k_10_r_1_1090_10_T_1000.txt", header = TRUE, sep = "", dec = ".")