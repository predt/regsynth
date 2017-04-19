### EXAMPLE 4: Geithner connections
### Jeremy L Hour
### 11 avril 2017

setwd("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth")
# Maison:
# setwd("/Users/jeremylhour/Documents/R/regsynth")

rm(list=ls())

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
source("functions/perm.test.R")
source("functions/conf.interval.R")
source("functions/conf.interval.Geithner.R")

### 0. Loading data
data = readMat("//ulysse/users/JL.HOUR/1A_These/A. Research/RegSynthProject/regsynth/data/GeithnerConnexions/Matlab Files/Data.mat")

# Maison:
#data = readMat("data/GeithnerConnexions/Matlab Files/Data.mat")


### 1. Data cleaning and setting
ticker = data$ticker # firms tickers

# collect names and tickers
FirmID = data.frame()
for(i in 1:603){
  if(length(unlist(ticker[i])) > 0 & length(unlist(ticker[603+i])) > 0){
    FirmID[i,"Ticker"] = unlist(ticker[i])
    FirmID[i,"Name"] = unlist(ticker[603+i])
  } else {
    FirmID[i,"Ticker"] = "Unknown"
    FirmID[i,"Name"] = "Unknown"
  }
}

X = data.frame(data$num) # firms characteristics
names(X) = c(unlist(data$VarNames)) # setting variable names
row.names(X) = FirmID[,"Ticker"] # setting firms name

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
# No need to exclude BoA (tax problem after nomination, check page 29)
Citi = which(X[,5]==140)  # Citi Group 
corrCiti = cor(y[PreTreatPeriod,Citi], y[PreTreatPeriod,])
# Compute Q10 for correlation distributions
corrCitiTr = sort(corrCiti,decreasing=T)[58]

X = X[corrCiti<corrCitiTr,]
y = y[,corrCiti<corrCitiTr]

# Treatment variable
d = X[,ConnMeasure] > 0 # one or more meetings in 2007-08

# What are the treated?
print("Treated frim and number of meetings in 2007-08")
cbind(FirmID[match(rownames(X[d==1,]), FirmID[,1]),2], X[d==1,ConnMeasure])

# Control variable other than pre-treatment outcomes, useless for now
# Include:
# - ta2008_log : firm size
# - roe2008 : profitability
# - tdtc2009 : leverage
Z = cbind(X[,c(8,9,10)], X[,c(8,9,10)]^2, X[,c(8,9,10)]^3)


### 2. Some descriptive statistics (TO BE CONTINUED ?)
ConnReturns = ts(apply(y[PreTreatPeriod,d==1],1,mean),start=c(1), freq=365)
NConnReturns = ts(apply(y[PreTreatPeriod,d==0],1,mean),start=c(1), freq=365)

# Balance check
# Treated
apply(X[d==1,c(8,9,10)],2,summary)
# Control
apply(X[d==0,c(8,9,10)],2,summary)

# Nested support seems to hold

### 3. CV for selecting optimal lambda
X0 = y[PreTreatPeriod,d==0]; X1 = y[PreTreatPeriod,d==1]
Y0 = y[GeiNomDate+1,d==0]; Y1 = y[GeiNomDate+1,d==1]

V = diag(nrow(X0)) # Not sure we need to reweight as all X are currently the same scale

lambda = seq(0,1.5,.1) # sequence of lambdas to test
estval = regsynthpath(X0,X1,Y0,Y1,V,lambda,tol=1e-6)
MSPE = vector(length=length(lambda))

for(k in 1:length(lambda)){
  MSPE[k] = mean(apply((y[324:354,d==1] - y[324:354,d==0]%*%t(estval$Wsol[k,,]))^2,2,mean))
}
lambda.opt.MSPE = min(lambda[which(MSPE==min(MSPE))]) # Optimal lambda is .2

### Figure 1: MSPE
pdf("plot/Geithner_MSPE.pdf", width=6, height=6)
matplot(lambda, MSPE, type="b", pch=20, lwd=1,
        main=expression("MSPE, "*lambda^{opt}*"= .2"), col="steelblue",
        xlab=expression(lambda), ylab="MSPE")
abline(v=lambda.opt.MSPE,lty=2,lwd=2,col="grey")
dev.off()



### 4. Estimation
Wsol = estval$Wsol[which(MSPE==min(MSPE)),,]
colnames(Wsol) = rownames(X[d==0,])

# Number of active controls
apply(Wsol>0,1,sum)

# Balance check on Z
apply(X[d==1,c(8,9,10)],2,mean)
apply(Wsol%*%as.matrix(X[d==0,c(8,9,10)]),2,mean)

# Compute the statistics (see paper)
TestPeriod = (GeiNomDate-15):(GeiNomDate+30)
sigma = sqrt(apply((X1 - X0%*%t(Wsol))^2,2,mean)) # Goodness of fit for each treated over pre-treatment period, used in the original paper
# omega = 1/(sigma*sum(1/sigma))
# phi = (y[TestPeriod,d==1] - y[TestPeriod,d==0]%*%t(Wsol))%*%omega
phi = apply((y[TestPeriod,d==1] - y[TestPeriod,d==0]%*%t(Wsol)),1,mean)
sigma_cutoff = mean(sigma) # for later use: correction during Fisher test
# Shows the path of abnormal returns
plot(phi,type="l")


### 5. Fisher Test of No-Effect Assumption (C=0)
set.seed(3101990)
R = 5000 # Number of replications
alpha = sqrt(3) # correction cut-off (see paper)
Result = matrix(nrow=R, ncol=length(TestPeriod))
Result_C = matrix(nrow=R, ncol=length(TestPeriod))
t_start = Sys.time()
pb = txtProgressBar(style = 3)
for(r in 1:R){
  dstar = sample(d)
  X0star = y[PreTreatPeriod,dstar==0]; X1star = y[PreTreatPeriod,dstar==1]
  solstar = regsynth(X0star,X1star,Y0,Y1,V,lambda)
  
  # Not corrected
  Result[r,] = apply((y[TestPeriod,dstar==1] - y[TestPeriod,dstar==0]%*%t(solstar$Wsol)),1,mean)
  
  # Corrected
  omegastar_C = rep(1,sum(d))
  sigmastar = sqrt(apply((X1star - X0star%*%t(solstar$Wsol))^2,2,mean))
  omegastar_C[sigmastar>alpha*sigma_cutoff] = 0
  omegastar_C = omegastar_C/sum(omegastar_C)
  Result_C[r,] = (y[TestPeriod,dstar==1] - y[TestPeriod,dstar==0]%*%t(solstar$Wsol))%*%omegastar_C
  
  setTxtProgressBar(pb, r/R)
}
close(pb)
print(Sys.time()-t_start)

### A. Not corrected

# Compute .025 and .975 quantiles of CAR for each date
phi_q005 = mapply(function(t) quantile(Result[,t], probs = .005), 1:length(TestPeriod))
phi_q025 = mapply(function(t) quantile(Result[,t], probs = .025), 1:length(TestPeriod))
phi_q975 = mapply(function(t) quantile(Result[,t], probs = .975), 1:length(TestPeriod))
phi_q995 = mapply(function(t) quantile(Result[,t], probs = .995), 1:length(TestPeriod))

ATTdata = ts(cbind(phi_q005,phi_q025,phi,phi_q975,phi_q995),start=c(-15), freq=1)

### Figure 2: Geithner connected firms effect vs. random permutations
pdf("plot/GeithnerAR_FisherTest.pdf", width=10, height=6)
plot(ATTdata, plot.type="single",
     col=c("firebrick","firebrick","firebrick","firebrick","firebrick"), lwd=c(1,1,2,1,1),
     lty=c(3,4,1,4,3),xlab="Day", ylab="AR, in pp",
     ylim=c(-.1,.1),
     main="Abnormal Returns (AR) for True Treatment vs. Random Permutations")
abline(h=0,
       lty=2,col="grey")
lim <- par("usr")
rect(0, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(-15,-.075,
       legend=c("Estimate", ".95 confidence bands of Fisher distrib.",".99 confidence bands of Fisher distrib."),
       col=c("firebrick","firebrick"), lwd=c(2,1,1),
       lty=c(1,4,3))
dev.off()

### B. Corrected

phi_q005 = mapply(function(t) quantile(Result_C[,t], probs = .005), 1:length(TestPeriod))
phi_q025 = mapply(function(t) quantile(Result_C[,t], probs = .025), 1:length(TestPeriod))
phi_q975 = mapply(function(t) quantile(Result_C[,t], probs = .975), 1:length(TestPeriod))
phi_q995 = mapply(function(t) quantile(Result_C[,t], probs = .995), 1:length(TestPeriod))

ATTdata = ts(cbind(phi_q005,phi_q025,phi,phi_q975,phi_q995),start=c(-15), freq=1)

### Figure 3: Geithner connected firms effect vs. random permutations, corrected
pdf("plot/GeithnerAR_FisherTestCorrected.pdf", width=10, height=7)
plot(ATTdata, plot.type="single",
     col=c("firebrick","firebrick","firebrick","firebrick","firebrick"), lwd=c(1,1,2,1,1),
     lty=c(3,4,1,4,3),xlab="Day", ylab="AR, in pp",
     ylim=c(-.1,.1),
     main="")
abline(h=0,
       lty=2,col="grey")
lim <- par("usr")
rect(0, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(0.5,-.065,
       legend=c("Estimate", ".95 confidence bands of Fisher distrib.",".99 confidence bands of Fisher distrib."),
       col=c("firebrick","firebrick"), lwd=c(2,1,1),
       lty=c(1,4,3))
dev.off()


### CAR[0,1] and CAR[0,10] Table, non corrected version
cumphi = cumsum(phi[16:length(phi)])

cumResult = t(apply(Result[,16:length(phi)],1,cumsum))
cumphi_q005 = mapply(function(t) quantile(cumResult[,t], probs = .005), 1:ncol(cumResult))
cumphi_q025 = mapply(function(t) quantile(cumResult[,t], probs = .025), 1:ncol(cumResult))
cumphi_q975 = mapply(function(t) quantile(cumResult[,t], probs = .975), 1:ncol(cumResult))
cumphi_q995 = mapply(function(t) quantile(cumResult[,t], probs = .995), 1:ncol(cumResult))

Table5 = data.frame("Estimate"=cumphi,"Q0.05"=cumphi_q005,"Q0.25"=cumphi_q025,
                    "Q97.5"=cumphi_q975,"Q99.5"=cumphi_q995)

print("Event day 0")
print(Table5[1,])

print("Event day 10")
print(Table5[11,])

ATTdata = ts(cbind(cumphi_q005,cumphi_q025,cumphi,cumphi_q975,cumphi_q995),start=c(1), freq=1)
### C. Figure 3: Geithner connected firms effect vs. random permutations
pdf("plot/GeithnerCAR_FisherTest.pdf", width=10, height=6)
plot(ATTdata, plot.type="single",
     col=c("firebrick","firebrick","firebrick","firebrick","firebrick"), lwd=c(1,1,2,1,1),
     lty=c(3,4,1,4,3),xlab="Day", ylab="CAR, in pp",
     ylim=c(-.25,.25),
     main="Cumulative Abnormal Returns (CAR) for True Treatment vs. Random Permutations")
abline(h=0,
       lty=2,col="grey")
legend(1,-.15,
       legend=c("Estimate", ".95 confidence bands of Fisher distrib.",".99 confidence bands of Fisher distrib."),
       col=c("firebrick","firebrick"), lwd=c(2,1,1),
       lty=c(1,4,3))
dev.off()

### CAR[0,1] and CAR[0,10] Table, corrected version
cumResult_C = t(apply(Result_C[,16:length(phi)],1,cumsum))
cumphi_q005 = mapply(function(t) quantile(cumResult_C[,t], probs = .005), 1:ncol(cumResult))
cumphi_q025 = mapply(function(t) quantile(cumResult_C[,t], probs = .025), 1:ncol(cumResult))
cumphi_q975 = mapply(function(t) quantile(cumResult_C[,t], probs = .975), 1:ncol(cumResult))
cumphi_q995 = mapply(function(t) quantile(cumResult_C[,t], probs = .995), 1:ncol(cumResult))

Table5_Corrected = data.frame("Estimate"=cumphi,"Q0.05"=cumphi_q005,"Q0.25"=cumphi_q025,
                    "Q97.5"=cumphi_q975,"Q99.5"=cumphi_q995)

print("Event day 0")
print(Table5_Corrected[1,])

print("Event day 10")
print(Table5_Corrected[11,])

ATTdata = ts(cbind(cumphi_q005,cumphi_q025,cumphi,cumphi_q975,cumphi_q995),start=c(1), freq=1)
### D. Figure 4: Geithner connected firms effect vs. random permutations
pdf("plot/GeithnerCAR_FisherTestCorrected.pdf", width=10, height=6)
plot(ATTdata, plot.type="single",
     col=c("firebrick","firebrick","firebrick","firebrick","firebrick"), lwd=c(1,1,2,1,1),
     lty=c(3,4,1,4,3),xlab="Day", ylab="CAR, in pp",
     ylim=c(-.25,.25),
     main="Cumulative Abnormal Returns (CAR) for True Treatment vs. Random Permutations, \n Corrected Inference")
abline(h=0,
       lty=2,col="grey")
legend(1,-.15,
       legend=c("Estimate", ".95 confidence bands of Fisher distrib.",".99 confidence bands of Fisher distrib."),
       col=c("firebrick","firebrick"), lwd=c(2,1,1),
       lty=c(1,4,3))
dev.off()

library(stargazer)
ToPrint = t(rbind(Table5[1,],Table5[11,],Table5_Corrected[1,],Table5_Corrected[11,]))
stargazer(t(ToPrint))


### .95 Confidence interval for CAR[0] and CAR[10]
GeiCI0 = conf.interval.Geithner(d,as.matrix(y[GeiNomDate,]),t(y[PreTreatPeriod,]),V,lambda=.1,B=5000,alpha=.05)
fileConn = file("plot/outputCI0.txt")
writeLines(paste(GeiCI0$c.int), fileConn)
close(fileConn)

GeiCI10 = conf.interval.Geithner(d,t(y[GeiNomDate:(GeiNomDate+10),]),t(y[PreTreatPeriod,]),V,lambda=.1,B=5000,alpha=.05)
fileConn = file("plot/outputCI10.txt")
writeLines(paste(GeiCI10$c.int), fileConn)
close(fileConn)