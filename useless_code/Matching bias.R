### Bias in matching
### Jeremy L'Hour
### 7 july 2016

setwd("/Volumes/USB_KEY/New_project/SyntheticControlMatching/code")

rm(list=ls())
set.seed(30031987)


### 0. Settings

### Load packages
library("MASS")
library("ggplot2")
library("gtable")
library("grid")

### Load user functions


### 1. A very simple example
n <- 500
x <- rexp(n)
#d <- as.numeric(runif(n) < 0*pnorm(x))
d <- as.numeric(runif(n) < .5)

data <- data.frame(X=x,
                   d=d)

### Control group is concentrated around small X values
p <- ggplot(data, aes(x=X, fill=as.factor(d))) + 
        geom_density(alpha=.3) +
        scale_x_continuous(name="X", limits=c(0,7)) +
        ggtitle("X distribution for Treated / Control group") + 
        theme(plot.title = element_text(lineheight=.8, face="bold")) +
        scale_fill_discrete(name="Group",
                      breaks=c("1", "0"),
                      labels=c("Treatment","Control")) +
        theme(legend.position="bottom") 

y <- .5 + x + rnorm(n, sd=.5)
data[,"y"] <- y

### Plotting on second axis
p2 = ggplot(data, aes(x=X, fill=as.factor(d))) +
  geom_point(data=data, aes(x=X, y=y), colour="blue")

### Matching
X1 = x[d==1]
X0 = x[d==0]
matched <- vector(length=sum(d))
for(i in 1:sum(d)){
  val <- X1[i]
  dist <- (X0 - val)^2
  matched [i] <- which.min(dist)
}

mean(X1 - X0[matched])
y0 <- y[d==0]
# There's a big bias because the support is not the same
# and matching is too crude
ATT <- mean(y[d==1] - y0[matched])
print(ATT)

# regression
reg <- lm(y ~ x, data=data, subset=(d==0))