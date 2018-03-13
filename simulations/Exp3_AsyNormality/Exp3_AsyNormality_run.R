### Exp3: Asymproric Normality, Run file
### Jeremy L Hour
### 21/02/2018

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
library("gridExtra")

### Load user functions
source("functions/wsol.R")
source("functions/wsoll1.R")
source("functions/matchDGP_fixedT.R")
source("functions/wATT.R")
source("functions/matching.R")
source("functions/matchest.R")
source("functions/OBest.R")
source("functions/regsynth.R")
source("functions/regsynthpath.R")
source("functions/TZero.R")
source("functions/synthObj.R")
source("simulations/Exp3_AsyNormality/Exp3_AsyNormality_setup.R")


### MC XP
set.seed(2121988)
lambda = seq(0,2,.01) # set of lambda to be considered for optim

xp = Exp3_setup(R=5000,n1=20,n0=70,p=20,K=5)

Results = xp
R = nrow(Results)

# Draw the charts
id = c(mapply(function(x) rep(x,R),1:5))
val = c(Results)
data_res = data.frame(val = val, model = id)

M = max(abs(quantile(Results,.01)),abs(quantile(Results,.99)))
lb = -1.1*M; ub = 1.1*M

sdBCH = sd(Results[,1])

### Function for plot
get.plot <- function(data,modelS,title="A Title",sdBCH){
  plot_res = ggplot(subset(data, (model==modelS)), aes(x=val)) + 
    geom_histogram(binwidth = .1, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
    scale_x_continuous(limits=c(lb,ub), name="Treatment effect") +
    ggtitle(title) + 
    stat_function(fun = dnorm, args=list(mean=0, sd=sdBCH), colour="darkorchid3", size=1) +
    theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
  
  return(plot_res)
}


grid.arrange(get.plot(data_res,1,"Fixed lambda", sdBCH), get.plot(data_res,2,"RMSE opt", sdBCH), ncol=2)


