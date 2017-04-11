### EXAMPLE 4: Geithner connections
### Jeremy L Hour
### 11 avril 2017

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
library("R.matlab")

### Load user functions
source("functions/wsol.R")
source("functions/wsoll1.R")
source("functions/matchDGP.R")
source("functions/wATT.R")
source("functions/matching.R")
source("functions/matchest.R")
source("functions/regsynth.R")
source("functions/regsynthpath.R")
source("functions/TZero.R")
source("functions/synthObj.R")

### 0. Loading data
data = readMat("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth/data/GeithnerConnexions/Matlab Files/Data.mat")

### 1. Data cleaning and setting
ticker = data$ticker # firms tickers

X = data.frame(data$num) # firms characteristics
names(X) = c(unlist(data$VarNames)) # setting variable names
for(i in 1:603){
  if(length(unlist(ticker[i])) > 0){
    row.names(X)[i] = unlist(ticker[i])
  } else {
    row.names(X)[i] = "Uknown"
  }
}

ind = is.na(X[,8]) | is.na(X[,9]) | is.na(X[,10])  # eliminating firms with no data for 'ta2008_log','roe2008','tdtc2008'
X = X[!ind,]
y = as.matrix(data$Re) # Returns
y = y[,!ind]

y[is.na(y)] = 0 # replacing missing returns with zero

# Identification of the event
ConnMeasure = 3 # 1: Shared Board 2: NY Connection 3: Geithner Schedule 4: Geithner Schedule 2007, position in data frame
GeiNomDate = 355 # Geithner nomination date
EventDate = GeiNomDate-1 
PreTreatPeriod = (GeiNomDate-281):(GeiNomDate-32) # Window of 250 days ending 30 days prior to Geithner nomination
FalsifTest = c(340:353,357:389,395:447) # Window for falsification test

# Correlation with Citi and BoA on Pre-treatment period
# We want to exclude the effect of the CitiGroup bailout
# Not sure whether we need to exclude BoA correlated firms (check page 30 of Acemoglu paper)
Citi = which(X[,5]==140)  # Citi Group 
BoA = which(X[,5]==56)    # Bank of America

corrCiti = cor(y[PreTreatPeriod,Citi], y[PreTreatPeriod,])
corrBoA = cor(y[PreTreatPeriod,BoA], y[PreTreatPeriod,])

# Compute Q10 for correlation distributions
corrCitiTr = sort(corrCiti,decreasing=T)[58]
corrBoATr = sort(corrBoA,decreasing=T)[58]

# Treatment variable
d = X[,ConnMeasure] > 0

# Control variable other than pre-treatment outcomes
# Include:
# - ta2008_log : firm size
# - roe2008 : profitability
# - tdtc2009 : leverage
Z = cbind(X[,c(8,9,10)], X[,c(8,9,10)]^2, X[,c(8,9,10)]^3)


### 2. Some descriptive statistics
ConnReturns = ts(apply(y[PreTreatPeriod,d==1],1,mean),start=c(1), freq=365)
NConnReturns = ts(apply(y[PreTreatPeriod,d==0],1,mean),start=c(1), freq=365)

# TO DO: charts and balance checks...

### 2. Estimation and Synthetic Control
X0 = y[PreTreatPeriod,d==0] # For now only pre-treatment outcomes
X1 = y[PreTreatPeriod,d==1]

Y0 = y[GeiNomDate+1,d==0]
Y1 = y[GeiNomDate+1,d==1]

V = cor(t(X0))

solution = regsynth(X0,X1,Y0,Y1,V,1,tol=1e-6) # setting lambda = 1 gives a non-NN solution


# Abnormal returns are defined as returns minus synthetic control returns
AR = y[,d==1] - y[,d==0] %*% t(solution$Wsol)

# Draw cumulated returns chart
# select a windo with few dates to better see the chart
TV0 = GeiNomDate - 5
TVF = GeiNomDate + 30
cumy = apply(y,2,cumsum)

plotdata = ts(cbind(apply(cumy[TV0:TVF,d==1],1,mean), apply(cumy[TV0:TVF,d==0] %*% t(solution$Wsol),1,mean)),start=c(1), freq=1)


plot(plotdata, plot.type="single",
     col=c("steelblue","firebrick"), lwd=2,
     lty=c(1,6),xlab="", ylab="Cumulated returns",
     ylim=c(-0.2,0.2))
lim <- par("usr")
rect(5, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(1971,80,
       legend=c("Real Connected Firms", "Synthetic Firms"),
       col=c("steelblue","firebrick","forestgreen"), lwd=2,
       lty=c(1,6,6))






data = data.frame(t(read.table("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth/data/MLAB_data.txt")))

Names = c("State_ID","Income","RetailPrice", "Young", "BeerCons","Smoking1988", "Smoking1980","Smoking1975",
          mapply(function(x) paste("SmokingCons",x,sep=""),1970:2000))
colnames(data) = Names
States = c("Alabama", "Arkansas","Colorado","Connecticut","Delaware",
           'Georgia',  'Idaho',  'Illinois',  'Indiana', 'Iowa', 'Kansas',
           'Kentucky', 'Louisiana', 'Maine', 'Minnesota', 'Mississippi',
           'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire',
           'New Mexico', 'North Carolina', 'North Dakota', 'Ohio', 'Oklahoma',
           'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota',
           'Tennessee', 'Texas', 'Utah', 'Vermont', 'Virginia' , 'West Virginia',
           'Wisconsin', 'Wyoming', 'California')
rownames(data) = States
data[,"Treated"]= as.numeric(data[,"State_ID"]==3) #California is state with ID=3

CASmoke = ts(unlist(data[data[,"Treated"]==1, mapply(function(x) paste("SmokingCons",x,sep=""),1970:2000)]),
             start=c(1970), freq=1)

### 1. Arbitrary tuning parameter
X = data[,c("Income","RetailPrice", "Young", "BeerCons",
            mapply(function(x) paste("SmokingCons",x,sep=""),1970:1988))]
d = data[,"Treated"]
V = diag(ncol(X))

X0 = t(X[d==0,])
X1 = t(X[d==1,])

### Regularization path
lambda = seq(0,1,.001)
Wsol = matrix(nrow=length(lambda), ncol=sum(1-d))
att = vector(length = length(lambda))

for(i in 1:length(lambda)){
  sol = wsoll1(X0,X1,V,lambda[i])
  sol = TZero(sol)
  Wsol[i,] = sol
}

colnames(Wsol) = States[States!="California"]
print("Non-zero weights for the pure synthetic control:")
Wsol[1,Wsol[1,]!=0]

# Weight as function of penalty level
matplot(lambda,Wsol, type="l", lwd=2,
        main="Regularization path",
        xlab="Penalty level", ylab="weight", ylim=c(0,1))


# All possible counterfactuals
y = data[rownames(data)!="California",mapply(function(x) paste("SmokingCons",x,sep=""),1970:2000)]  
SyntheticControl = Wsol %*% as.matrix(y)

matplot(1970:2000,t(SyntheticControl), type="l",
        main="All possible counterfactuals",
        xlab="Year", ylab="weight", ylim=c(35,150))
lim <- par("usr")
rect(1988, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 


### 2. Cross-validation to select tuning parameter (time-based)
varname = mapply(function(x) paste("SmokingCons",x,sep=""),1970:1980)
Xtrain = cbind(data[,c("Income","RetailPrice", "Young", "BeerCons", varname)])

V = diag(ncol(Xtrain))
X0 = t(Xtrain[d==0,])
X1 = t(Xtrain[d==1,])

### Computing with all lambda's on training sample
lambda = seq(0,3,.01)
Wsol = matrix(nrow=length(lambda), ncol=sum(1-d))
att = vector(length = length(lambda))

for(i in 1:length(lambda)){
  sol = wsoll1(X0,X1,V,lambda[i])
  sol = TZero(sol)
  Wsol[i,] = sol
}

colnames(Wsol) = States[States!="California"]

### See performance on test sample
varname = mapply(function(x) paste("SmokingCons",x,sep=""),1983:1988)
Xtest = data[rownames(data)!="California",varname]
SyntheticControl = Wsol %*% as.matrix(Xtest)
MSPE = (SyntheticControl - kronecker(matrix(1,nrow=length(lambda)),as.matrix(data[rownames(data)=="California",varname])))^2
MSPE = apply(MSPE,1,mean)


matplot(lambda,MSPE, type="o", pch=20,
        main="MSPE", col="steelblue",
        xlab=expression(lambda), ylab="MSPE")

lambda.opt.MSPE = min(lambda[which(MSPE==min(MSPE))])

varname = mapply(function(x) paste("SmokingCons",x,sep=""),1970:2000)
y = data[rownames(data)!="California",varname]  
SyntheticControl = t(Wsol[which(MSPE==min(MSPE)),] %*% as.matrix(y))
OriginalSC = t(Wsol[1,] %*% as.matrix(y))

plotdata = ts(cbind(t(data[rownames(data)=="California",varname]), SyntheticControl, OriginalSC),start=c(1970), freq=1)


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
       legend=c("Real California", "Synthetic Control, opt lambda", "Original Synthetic Control"),
       col=c("steelblue","firebrick","forestgreen"), lwd=2,
       lty=c(1,6,6))


### 3. Cross-validation (unit based)
varname = mapply(function(x) paste("SmokingCons",x,sep=""),1970:2000)
Xtrain = cbind(data[d==0,c("Income","RetailPrice", "Young", "BeerCons", varname)])
testvarname = mapply(function(x) paste("SmokingCons",x,sep=""),1989:2000)
Xtest = cbind(data[d==0,testvarname])

V = diag(ncol(Xtrain))
dstar = rep(0,sum(d==0))
lambda = seq(0,3.5,.01)
MSPE = matrix(nrow=sum(d==0), ncol=length(lambda))

for(j in 1:sum(d==0)){
  dstar[j] = 1 
  X0 = t(Xtrain[dstar==0,])
  X1 = t(Xtrain[dstar==1,])
  
  Wsol = matrix(nrow=length(lambda), ncol=sum(1-dstar))
  
  for(i in 1:length(lambda)){
    sol = wsoll1(X0,X1,V,lambda[i])
    sol = TZero(sol)
    Wsol[i,] = sol
  }
  
  SyntheticControl = Wsol %*% as.matrix(Xtest[dstar==0,])
  MSPE_i = (SyntheticControl - kronecker(matrix(1,nrow=length(lambda)),as.matrix(Xtest[dstar==1,])))^2
  MSPE[j,] = apply(MSPE_i,1,mean)
  
  dstar[j] = 0 
}

MSPETot = apply(MSPE,2,mean)

matplot(lambda,MSPETot, type="o", pch=20,
        main="MSPE", col="steelblue",
        xlab=expression(lambda), ylab="MSPE")


### Computing with all lambda's on training sample
lambda = seq(0,3.5,.001)
Wsol = matrix(nrow=length(lambda), ncol=sum(1-d))

for(i in 1:length(lambda)){
  sol = wsoll1(X0,X1,V,lambda[i])
  sol = TZero(sol)
  Wsol[i,] = sol
}

colnames(Wsol) = States[States!="California"]

### See performance on test sample
varname = mapply(function(x) paste("SmokingCons",x,sep=""),1983:1988)
Xtest = data[rownames(data)!="California",varname]
SyntheticControl = Wsol %*% as.matrix(Xtest)
MSPE = (SyntheticControl - kronecker(matrix(1,nrow=length(lambda)),as.matrix(data[rownames(data)=="California",varname])))^2
MSPE = apply(MSPE,1,mean)


matplot(lambda,MSPE, type="o", pch=20,
        main="MSPE", col="steelblue",
        xlab=expression(lambda), ylab="MSPE")

lambda.opt.MSPE = min(lambda[which(MSPE==min(MSPE))])