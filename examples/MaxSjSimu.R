### Distribution of maximum of S_j
### when w_i follows a Dirichlet
### Jeremy L Hour
### 15/06/2017

setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")
rm(list=ls())

### 0. Settings

### Load packages
library("MASS")
library("gtools")
library("ggplot2")

### Parameters
R = 50000
n1 = 12 # number of treated
n0 = 200 # number of controls

### 1. A first look
set.seed(18041991)

alpha = rep(.5, n0) # It gives a lot of positive and negative weights, which is similar to what we have with our penalization type
data = vector(length=R)

for(r in 1:R){
  draw = rdirichlet(n1, alpha)
  Sj = apply(draw,2,sum)
  data[r] = max(Sj)
}

data = data.frame("max"=data)

ggplot(data, aes(x=max)) + 
  geom_histogram(binwidth = .002, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(name="Max of S_j") +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")

print(quantile(data[,"max"], probs = c(.95,.975,.99)))

### Further study
# Here we move n1 and n0
set.seed(18041991)
R = 10000
N1 = c(10,20,50,80,100,150,200)
N0 = c(10,20,50,80,100,150,200)

sim = array(dim=c(length(N1),length(N0),R))
QTable = data.frame("N1"=integer(),
                    "N0"=integer(),
                    "Q95"=double(),
                    "Q975"=double(),
                    "Q99"=double())

for(i in 1:length(N0)){
  alpha = rep(.5, N0[i]) 
  for(j in 1:length(N1)){
    
    for(r in 1:R){
      draw = rdirichlet(N1[j], alpha)
      Sj = apply(draw,2,sum)
      sim[j,i,r] = max(Sj)/N1[j] #max normalized by treated group size
    }
    QTable = rbind(QTable,
                   c(N0[i],N1[j],quantile(sim[j,i,], probs = c(.95,.975,.99)))
                   )
    print(paste("N1 count:",100*j/length(N1)))
  }
  print(paste("N0 count:",100*i/length(N0)))
}

colnames(QTable)=c("N1","N0","Q95","Q975","Q99")

### Final Table
FTable = matrix(nrow=length(N1), ncol=length(N0))
for(i in 1:length(N0)){
  for(j in 1:length(N1)){
    FTable[j,i] = QTable[which(QTable[,"N1"]==N1[j] & QTable[,"N0"]==N0[i]),"Q95"]
  }
}

heatmap(FTable)