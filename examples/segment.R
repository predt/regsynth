### Is it on the segment?
### Jeremy L Hour
### 11 juillet 2018

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
library("deldir")
library("plotrix")

### Load user functions
source("functions/wsol.R")
source("functions/wsoll1.R")
source("functions/PolyDGP.R")
source("functions/wATT.R")
source("functions/matching.R")
source("functions/matchest.R")
source("functions/OBest.R")
source("functions/regsynth.R")
source("functions/regsynthpath.R")
source("functions/TZero.R")
source("functions/synthObj.R")
source("functions/which.tile.R")

### MC XP
set.seed(2121988)
lambda = seq(.001,30,.1) # set of lambda to be considered for optim
n1 = 1
n0 = 100
p = 2
delta = 2

# Setting up data
data = PolyDGP(n1,n0,p,delta)
X = data$X; y = data$y; d = data$d

X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))
Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)

# Synthetic control for each lambda
solpath = regsynthpath(X0,t(X1),Y0,Y1,diag(p),lambda)

W0 = drop(solpath$Wsol)
SC = W0%*%t(X0)

# Find the nearest neighbor
NN = matching(X0,t(X1),diag(p),m=1)
Xnn = X0%*%NN

# Find Delaunay tesselation (control only and control+treatd)
DTco = deldir(X0[1,], X0[2,])
DTct = deldir(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]))

# Connections in DT...

# in augmented DT
# treated unit is index n0+1
union(DTct$delsgs[DTct$delsgs[,"ind2"]==n0+1,"ind1"],
      DTct$delsgs[DTct$delsgs[,"ind1"]==n0+1,"ind2"])

# in control DT
tril = triang.list(DTco)
which.tile(X1[,1], X1[,2], tril)

# which Voronoi cell?
tl = tile.list(deldir(X0[1,], X0[2,]))
which.tile(X1[,1], X1[,2], tl)

# used in synthetic controls (across all lambdas)
which(apply(W0,2,sum)>0)

# Plotting

# PLOT 1
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,ylim=c(.2,.5),xlim=c(.2,.5))
points(X0[1,], X0[2,], pch=20, col="orange", cex=0.5)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=.7)
plot(DTco, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=1,col="grey")
points(SC,type="l",col="red", lwd=2)
title(main="Synthetic control in DT of control units")

# PLOT 2
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,ylim=c(.2,.5),xlim=c(.2,.5))
points(X0[1,], X0[2,], pch=20, col="orange", cex=0.5)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=.7)
plot(DTct, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
points(SC,type="l",col="red", lwd=2)
title(main="Synthetic control in augmented DT")


### Try with another X1
X1 = matrix(c(-.5,-1.2),ncol=2)

# Synthetic control for each lambda
solpath = regsynthpath(X0,t(X1),Y0,Y1,diag(p),lambda)

W0 = drop(solpath$Wsol)
SC = W0%*%t(X0)

# Find the nearest neighbor
NN = matching(X0,t(X1),diag(p),m=1)
Xnn = X0%*%NN


# Find Delaunay tesselation (control only and control+treatd)
DTct = deldir(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]))

# Connections in DT...

# in augmented DT
# treated unit is index n0+1
union(DTct$delsgs[DTct$delsgs[,"ind2"]==n0+1,"ind1"],
      DTct$delsgs[DTct$delsgs[,"ind1"]==n0+1,"ind2"])

# which Voronoi cell?
tl = tile.list(deldir(X0[1,], X0[2,]))
which.tile(X1[,1], X1[,2], tl)

# used in synthetic controls (across all lambdas)
active = which(apply(W0,2,sum)>0)


# PLOT 1
#pdf("DelaunaySynthetic.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(-1,0),ylim=c(-2,-.8),)
points(X0[1,], X0[2,], pch=20, col="orange", cex=0.5)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=.7)
plot(DTco, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
points(SC,type="line",col="red", lwd=2)
points(t(X0[,active]), pch=17, col="purple", cex=1.5)
points(t(Xnn), pch=17, col="yellow", cex=1.5)
title(main="Synthetic control in DT of control units")
#dev.off()


# PLOT 2: circle
rad = sqrt(sum((t(X1)-Xnn)^2))
#pdf("DelaunaySynthetic.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(-1,0),ylim=c(-2,-.8),)
points(X0[1,], X0[2,], pch=20, col="orange", cex=0.5)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=.7)
plot(DTct, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
draw.circle(X1[,1], X1[,2],radius=rad,nv=100,border=NULL,lty=1)
points(t(X0[,active]), pch=17, col="purple", cex=1.5)
points(t(Xnn), pch=17, col="yellow", cex=1.5)
title(main="Synthetic control in DT of control units")
#dev.off()







####################
### Try with another X1
X1 = matrix(c(-.5,-1.4),ncol=2)

# Synthetic control for each lambda
solpath = regsynthpath(X0,t(X1),Y0,Y1,diag(p),lambda)

W0 = drop(solpath$Wsol)
SC = W0%*%t(X0)

# Find the nearest neighbor
NN = matching(X0,t(X1),diag(p),m=1)
Xnn = X0%*%NN


# Find Delaunay tesselation (control only and control+treatd)
DTct = deldir(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]))

# Connections in DT...

# in augmented DT
# treated unit is index n0+1
union(DTct$delsgs[DTct$delsgs[,"ind2"]==n0+1,"ind1"],
      DTct$delsgs[DTct$delsgs[,"ind1"]==n0+1,"ind2"])

# which Voronoi cell?
tl = tile.list(deldir(X0[1,], X0[2,]))
which.tile(X1[,1], X1[,2], tl)

# used in synthetic controls (across all lambdas)
active = which(apply(W0,2,sum)>0)


# PLOT 1
#pdf("DelaunaySynthetic.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(-1,0),ylim=c(-2,-.8),)
points(X0[1,], X0[2,], pch=20, col="orange", cex=0.5)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=.7)
plot(DTco, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
points(SC,type="l",col="red", lwd=2)
points(t(X0[,active]), pch=17, col="purple", cex=1.5)
points(t(Xnn), pch=17, col="yellow", cex=1.5)
title(main="Synthetic control in DT of control units")
#dev.off()


# PLOT 2: circle
rad = sqrt(sum((t(X1)-Xnn)^2))
#pdf("DelaunaySynthetic.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(-1,0),ylim=c(-2,-.8),)
points(X0[1,], X0[2,], pch=20, col="orange", cex=0.5)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=.7)
plot(DTct, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
draw.circle(X1[,1], X1[,2],radius=rad,nv=100,border=NULL,lty=1)
points(t(X0[,active]), pch=17, col="purple", cex=1.5)
points(t(Xnn), pch=17, col="yellow", cex=1.5)
points(SC,type="l",col="red", lwd=2)
draw.circle(X1[,1], X1[,2],radius=.8*rad,nv=100,border=NULL,lty=1)
title(main="Synthetic control in DT of control and treated units")
#dev.off()


#############################
### X1 outside convex hull
X1 = matrix(c(2,-2),ncol=2)

# Synthetic control for each lambda
solpath = regsynthpath(X0,t(X1),Y0,Y1,diag(p),lambda)

W0 = drop(solpath$Wsol)
SC = W0%*%t(X0)

# Find the nearest neighbor
NN = matching(X0,t(X1),diag(p),m=1)
Xnn = X0%*%NN


# Find Delaunay tesselation (control only and control+treatd)
DTct = deldir(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]))

# Connections in DT...

# in augmented DT
# treated unit is index n0+1
union(DTct$delsgs[DTct$delsgs[,"ind2"]==n0+1,"ind1"],
      DTct$delsgs[DTct$delsgs[,"ind1"]==n0+1,"ind2"])

# which Voronoi cell?
tl = tile.list(deldir(X0[1,], X0[2,]))
which.tile(X1[,1], X1[,2], tl)

# used in synthetic controls (across all lambdas)
active = which(apply(W0,2,sum)>0)


# PLOT 1
#pdf("DelaunaySynthetic.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1)
points(X0[1,], X0[2,], pch=20, col="orange", cex=0.5)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=.7)
plot(DTco, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
points(SC,type="l",col="red", lwd=2)
points(t(X0[,active]), pch=17, col="purple", cex=1.5)
points(t(Xnn), pch=17, col="yellow", cex=1.5)
title(main="Synthetic control in DT of control units")
#dev.off()


# PLOT 2: circle
rad = sqrt(sum((t(X1)-Xnn)^2))
#pdf("DelaunaySynthetic.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(1,2.5),ylim=c(-2.5,.5),)
points(X0[1,], X0[2,], pch=20, col="orange", cex=0.5)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=.7)
plot(DTct, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
draw.circle(X1[,1], X1[,2],radius=rad,nv=100,border=NULL,lty=1)
points(t(X0[,active]), pch=17, col="purple", cex=1.5)
points(t(Xnn), pch=17, col="yellow", cex=1.5)
title(main="Synthetic control in DT of control and treated units")
#dev.off()