
f <- function(x,l){
  y = (-2+3*x)^2 + l*(1+8*(1-x))
  return(y)
}

xset = seq(0,1,.00001)

yset = f(xset,0)

plot(xset, yset)

i = which(yset==min(yset))

print(xset[i])

X1 = as.matrix(2)
X0 = matrix(c(1,4,5),nrow=1)

lambda = .0005
sol = wsoll1(X0,X1,as.matrix(1),lambda)
sol = TZero(sol,1e-6)
print(sol)

# sol_th = (2+(4/3)*lambda)/3
sol_th = (2+lambda/2)/3
print(sol_th)