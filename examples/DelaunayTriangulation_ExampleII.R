### Is it on the segment?
### Jeremy L Hour
### 11 juillet 2018

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
SyntheticUnit = W0%*%t(X0)

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
points(SyntheticUnit,type="l",col="red", lwd=2)
title(main="Synthetic control in DT of control units")

# PLOT 2
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,ylim=c(.2,.5),xlim=c(.2,.5))
points(X0[1,], X0[2,], pch=20, col="orange", cex=0.5)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=.7)
plot(DTct, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
points(SyntheticUnit,type="l",col="red", lwd=2)
title(main="Synthetic control in augmented DT")

############################
############################
############################
### BEGINNING OF EXAMPLE ###
############################
############################
############################

# Take 100 control units distributed uniformly in [0,1] x [0,1]
# An 1 treated unit with same distribution.
X1 = matrix(c(.15,.385),ncol=2)

### 1. Here is the plot
pdf("plot/DelaunayTri/1_Points.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="")
points(X0[1,], X0[2,], pch=20, col="orange", cex=1,lwd=7)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=1, lwd=7)
dev.off()

# Take the Voronoi Tesselation of the control units
# It separates the space in regions of nearest-neighbors.

# which Voronoi cell?
tl = tile.list(deldir(X0[1,], X0[2,]))
which.tile(X1[,1], X1[,2], tl)

### 2. Same graph with Voronoi regions 
pdf("plot/DelaunayTri/2_Voronoi.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="")
points(X0[1,], X0[2,], pch=20, col="orange", cex=1, lwd=7)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=1, lwd=7)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
dev.off()

# Now the Dealaunay Triangulation is the dual of the Voronoi.
# This is the object of interest for us.

### 3. Same graph with Delaunay triangulation instead
pdf("plot/DelaunayTri/3_Delaunay.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="")
points(X0[1,], X0[2,], pch=20, col="orange", cex=1, lwd=7)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=1, lwd=7)
plot(DTco, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
dev.off()

# The union of all Delaunay triangles create a partition of the convex hull.

### 4. Now represent both.
pdf("plot/DelaunayTri/4_Voronoi+Delaunay.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="")
points(X0[1,], X0[2,], pch=20, col="orange", cex=1, lwd=7)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=1, lwd=7)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
plot(DTco, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
dev.off()

### 5. Let's zoom in
pdf("plot/DelaunayTri/5_Zoom.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(0,.3),ylim=c(.2,.6),
     xlab="",ylab="")
points(X0[1,], X0[2,], pch=20, col="orange", cex=1, lwd=7)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=1, lwd=7)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
plot(DTco, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
dev.off()

# Let's represent the synthetic control units for a wide range of lambda's.
# Synthetic control for each lambda
solpath = regsynthpath(X0,t(X1),Y0,Y1,diag(p),lambda)

W0 = drop(solpath$Wsol)
SyntheticUnit = W0%*%t(X0)

# Find the nearest neighbor
NN = matching(X0,t(X1),diag(p),m=1)
Xnn = X0%*%NN

### 6. Synthetic Unit for lambda's
pdf("plot/DelaunayTri/6_SyntheticUnit.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(0,.3),ylim=c(.2,.6),
     xlab="",ylab="")
points(X0[1,], X0[2,], pch=20, col="orange", cex=1, lwd=7)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=1, lwd=7)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
plot(DTco, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
points(SyntheticUnit,type="line",col="red", lwd=6)
dev.off()

# What are the active control units?
# used in synthetic controls (across all lambdas)
active = which(apply(W0,2,sum)>0)

### 7. 6+active units (Theorem 1)
pdf("plot/DelaunayTri/7_ActiveUnits.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(0,.3),ylim=c(.2,.6),
     xlab="",ylab="")
points(X0[1,], X0[2,], pch=20, col="orange", cex=1, lwd=7)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=1, lwd=7)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
plot(DTco, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
points(SyntheticUnit,type="line",col="red", lwd=6)
points(t(X0[,active]), pch=20, col="deepskyblue", cex=1, lwd=9)
points(t(Xnn), pch=20, col="purple", cex=1, lwd=9)
dev.off()

# If we create the Augmented Delaunay Triangulation (adding the treated unit),
# we can determine who might be active or not.
# Find Delaunay tesselation (control only and control+treatd)
DTct = deldir(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]))

# Connections in augmented DT
# treated unit is index n0+1
union(DTct$delsgs[DTct$delsgs[,"ind2"]==n0+1,"ind1"],
      DTct$delsgs[DTct$delsgs[,"ind1"]==n0+1,"ind2"])

# 8. Augmented Delaunay Triangulation (Lemma A.6)
pdf("plot/DelaunayTri/8_AugmentedDT.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(0,.3),ylim=c(.2,.6),
     xlab="",ylab="")
points(X0[1,], X0[2,], pch=20, col="orange", cex=1, lwd=6)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=1, lwd=6)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
plot(DTct, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
points(SyntheticUnit,type="line",col="red", lwd=6)
points(t(X0[,active]), pch=20, col="deepskyblue", cex=1, lwd=9)
points(t(Xnn), pch=20, col="purple", cex=1, lwd=9)
dev.off()

# Represent NN circle
rad = sqrt(sum((t(X1)-Xnn)^2))

# 9. Augmented Delaunay Triangulation + NN Circle (Lemma 1)
pdf("plot/DelaunayTri/9_AugmentedDT+Circle.pdf",width=10,height=10)
plot(c(X0[1,],X1[,1]), c(X0[2,],X1[,2]), type="n", asp=1,xlim=c(0,.3),ylim=c(.2,.6),
     xlab="",ylab="")
points(X0[1,], X0[2,], pch=20, col="orange", cex=1, lwd=6)
points(X1[,1], X1[,2], pch=19, col="forestgreen", cex=1, lwd=6)
plot(DTco, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=2,col="grey")
plot(DTct, wlines="triang", wpoints="none", number=FALSE, add=TRUE, lty=1)
draw.circle(X1[,1], X1[,2],radius=rad,nv=100,border="blue",lty=1,lwd=2)
points(SyntheticUnit,type="line",col="red", lwd=6)
points(t(X0[,active]), pch=20, col="deepskyblue", cex=1, lwd=9)
points(t(Xnn), pch=20, col="purple", cex=1, lwd=9)
dev.off()