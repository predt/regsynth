### Loi de la plus proche à droite
### 21/06/2017
### Jeremy L Hour

### Loading libraries
library("ggplot2")


### 1. Loi conditionelle sachant X_j
R = 10000
n0 = 15
X_j = .3

data = vector(length=R)

for(r in 1:R){
	X = runif(n0-1)
	if(X_j > max(X)){
		  data[r] = 1
		} else {
		  data[r] = min(X[X>X_j])
		}
}


data = data.frame("draw"=data)

ggplot(data, aes(x=draw)) + 
  geom_histogram(binwidth = .002, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(limits=c(X_j,1.01), name="u") +
  ggtitle("Loi du plus proche à droite") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")


### 2 bis. fonction de repartition
Frep <- function(x,b,n0){
  if(x==1){
    return(1)
  } else {
    return(1-(1-x+b)^(n0-1))
  }
}

Frep(X_j,b=X_j,n0)

# Empirical
ggplot(data, aes(x=draw)) +
  stat_ecdf(geom = "step") +
  scale_x_continuous(limits=c(X_j,1.01), name="u") +
  ggtitle("Loi du plus proche à droite") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
  
# Theoretical  
ggplot(data.frame(x=c(X_j, 1)), aes(x)) +
  stat_function(fun=Frep,  args=list(b=X_j,n0=n0)) +
  scale_x_continuous(limits=c(X_j,1.01), name="u") +
  ggtitle("Loi du plus proche à droite") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")


Frep(.5,b=X_j,n0) - mean(data[data<.5])

### 2. Loi inconditionelle
R = 100000
n0 = 10

data = vector(length=R)

for(r in 1:R){
	X = runif(n0)
	X_j = X[1]
	if(X_j == max(X)){
		  data[r] = 1
		} else {
		  data[r] = min(X[X>X_j])
		}
}

### Test sur la probabilité que ce soit le max
sqrt(R) * (mean(data==1) - 1/(n0+1))/sqrt(n0/(n0+1)^2)

data = data.frame("draw"=data)

ggplot(data, aes(x=draw)) + 
  geom_histogram(binwidth = .005, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(limits=c(0,1.01), name="u") +
  ggtitle("Loi du plus proche à droite") + 
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")

