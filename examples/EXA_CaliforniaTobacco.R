### EXAMPLE 1: California Tobacco COnsumption
### Sparse Synthetic Control
### Jeremy L Hour
### 11 Juillet 2016

setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")




rm(list=ls())
set.seed(3101990)

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

### Loading data
data <- data.frame(t(read.table("R:/Simulations/BEAST/datasets/MLAB_data.txt")))

Names <- c("State_ID","Income","RetailPrice", "Young", "BeerCons","Smoking1988", "Smoking1980","Smoking1975",
           mapply(function(x) paste("SmokingCons",x,sep=""),1970:2000))
colnames(data) <- Names
States <- c("Alabama", "Arkansas","Colorado","Connecticut","Delaware",
                    'Georgia',  'Idaho',  'Illinois',  'Indiana', 'Iowa', 'Kansas',
                    'Kentucky', 'Louisiana', 'Maine', 'Minnesota', 'Mississippi',
                    'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire',
                    'New Mexico', 'North Carolina', 'North Dakota', 'Ohio', 'Oklahoma',
                    'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota',
                    'Tennessee', 'Texas', 'Utah', 'Vermont', 'Virginia' , 'West Virginia',
                    'Wisconsin', 'Wyoming', 'California')
rownames(data) <- States
data[,"Treated"] <- as.numeric(data[,"State_ID"]==3) #California is state with ID=3

CaliSmoke <- ts(unlist(data[data[,"Treated"]==1, mapply(function(x) paste("SmokingCons",x,sep=""),1970:2000)]),
                start=c(1970), freq=1)

# Original version
X <- cbind(data[,c("Income","RetailPrice", "Young", "BeerCons"
                   , "SmokingCons1970", "SmokingCons1971", "SmokingCons1972", "SmokingCons1973", "SmokingCons1974", "SmokingCons1975"
                   , "SmokingCons1980", "SmokingCons1988"
           )])

d <- data[,"Treated"]

# Choice of V
#1. standard deviation
S = sqrt(diag(cov(X)))
V = diag(S)

#2. absolute value of the regression coefficient in 1988
reg = lm(SmokingCons1988 ~ ., data=X)
V = diag(c(abs(coef(reg)[-1]),1))
### Running Sparse Synthetic Control
X0 = t(X[d==0,])
X1 = t(X[d==1,])


### Regularization path
lambda = seq(0,3.5,.001)
Wsol = matrix(nrow=length(lambda), ncol=sum(1-d))
att = vector(length = length(lambda))

for(i in 1:length(lambda)){
  sol = wsoll1(X0,X1,V,lambda[i])
  sol = TZero(sol)
  Wsol[i,] = sol
}

colnames(Wsol) <- States[States!="California"]
Wsol[1,Wsol[1,]!=0]


# Weight as function of penalty level
pdf("plot/RegPath.pdf", width=10, height=6)
matplot(lambda,Wsol, type="l", lwd=2,
        main="Sparse Synthetic Control Regularization path",
        xlab="Penalty level", ylab="weight", ylim=c(0,1))
dev.off()

# All possible counterfactuals
varname = 1970:2000
i=0
for(t in 1970:2000){
  i=i+1
  varname[i] = paste("SmokingCons",t,sep="")
}
  
y <- data[rownames(data)!="California",varname]  
SyntheticControl <- Wsol %*% as.matrix(y)

pdf("plot/All.pdf", width=10, height=10)
matplot(1970:2000,t(SyntheticControl), type="l",
        main="All possible counterfactuals",
        xlab="Year", ylab="weight", ylim=c(35,150))
lim <- par("usr")
rect(1988, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
dev.off()

### Synthetic control plot with simple solution + sparse
plotdata <- ts(cbind(t(data[rownames(data)=="California",varname]), SyntheticControl[1,], SyntheticControl[2,]),start=c(1970), freq=1)

pdf("plot/SCvsSSC.pdf", width=10, height=10)
plot(plotdata, plot.type="single",
     col=c("steelblue","firebrick","forestgreen"), lwd=2,
     lty=c(1,6,6),xlab="", ylab="Cigarette consumption (Packs per capita)",
     ylim=c(35,150))
lim <- par("usr")
rect(1988, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(1971,80,
       legend=c("Real California", "Synthetic Control", "Sparse Synthetic Control"),
       col=c("steelblue","firebrick","forestgreen"), lwd=2,
       lty=c(1,6,6))
dev.off()


### MACHINE LEARNING: Select penalty level from
### training and test sample
varname = 1970:1980
i=0
for(t in 1970:1980){
  i=i+1
  varname[i] = paste("SmokingCons",t,sep="")
}
Xtrain <- cbind(data[,c("Income","RetailPrice", "Young", "BeerCons", varname)])

S = sqrt(diag(cov(Xtrain)))
V = diag(S)
X0 = t(Xtrain[d==0,])
X1 = t(Xtrain[d==1,])

### Regularization path
lambda = seq(0,3.5,.001)
Wsol = matrix(nrow=length(lambda), ncol=sum(1-d))
att = vector(length = length(lambda))

for(i in 1:length(lambda)){
  sol = wsoll1(X0,X1,V,lambda[i])
  sol = TZero(sol)
  Wsol[i,] = sol
}

colnames(Wsol) <- States[States!="California"]

### See performance on test sample
varname = 1981:1988
i=0
for(t in 1981:1988){
  i=i+1
  varname[i] = paste("SmokingCons",t,sep="")
}
Xtest = data[rownames(data)!="California",varname]
SyntheticControl = Wsol %*% as.matrix(Xtest)
MSPE = (SyntheticControl - kronecker(matrix(1,nrow=length(lambda)),as.matrix(data[rownames(data)=="California",varname])))^2
MSPE = apply(MSPE,1,mean)

pdf("plot/MSPE.pdf", width=10, height=10)
matplot(lambda,MSPE, type="o", pch=20,
        main="MSPE", col="steelblue",
        xlab=expression(lambda), ylab="MSPE")
dev.off()