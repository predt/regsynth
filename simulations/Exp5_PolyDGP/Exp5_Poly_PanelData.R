### Setup of Exp 5, Panel Data style
### 29/03/2018

Exp5_Poly_PanelData <- function(R=1000,n1=100,n0,p,delta,a=.1,b=.9){
  Results = matrix(ncol=13, nrow=R)
  t_start = Sys.time()
  pb = txtProgressBar(style = 3)
  
  for(r in 1:R){
    ### 0. Generate data
    data = PanelPolyDGP(n1,n0,p,delta,a,b)
    X = data$X; y = data$y; d = data$d
    
    X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))
    Y0 = c(y[d==0,2]); Y1 = c(y[d==1,2]); # True outcomes are at period 2 (col 2)
    Y0t = c(y[d==0,1]); Y1t = c(y[d==1,1]);
    n0 = sum(1-d) 
    
    ### 1. Synthetic Control on mean of treated
    M = matrix(apply(X1,1,mean), ncol=1)
    AggSynth = wATT(y[,2],d,wsoll1(X0,M,V))
    
    ### 2. NN Matching
    ## A. K = 1
    NN1 = matchest(X0,X1,Y0,Y1,V,m=1)
    ## B. K = 5
    NN5 = matchest(X0,X1,Y0,Y1,V,m=5)
    ## C. K = Kopt
    keeptauNN = matrix(nrow=10, ncol=n0)
      for(i in 1:10){
        sol = matchest(X0,X1,Y0t,Y1t,V,m=i)
        keeptauNN[i,] = sol$CATT
      }
    
    # The one that optimizes MSE
    curve.MSE = apply(keeptauNN^2,1,sum)/n0
    m.opt.MSE = min(which(curve.MSE==min(curve.MSE)))
    sol = matchest(X0,X1,Y0,Y1,V,m=m.opt.MSE)
    NN.opt.MSE = sol$ATT
    
    # The one that optimizes bias
    curve.bias = abs(apply(keeptauNN,1,sum)/n0)
    m.opt.bias = min(which(curve.bias==min(curve.bias)))
    sol = matchest(X0,X1,Y0,Y1,V,m=m.opt.bias)
    NN.opt.bias = sol$ATT
    
    # The one that optimizes bias + variance
    curve.crit = curve.bias + apply(keeptauNN,1,var)
    m.opt.crit = min(which(curve.crit==min(curve.crit)))
    sol = matchest(X0,X1,Y0,Y1,V,m=m.opt.crit)
    NN.opt.crit = sol$ATT
    
    # The one that optimizes MAE
    curve.MAE = apply(abs(keeptauNN),1,sum)/n0
    m.opt.MAE = min(which(curve.MAE==min(curve.MAE)))
    sol = matchest(X0,X1,Y0,Y1,V,m=m.opt.MAE)
    NN.opt.MAE = sol$ATT
    
    
    ### 3. Penalized Synthetic Control
    # A. lambda = 0
    sol = regsynth(X0,X1,Y0,Y1,V,0)
    Synth.NoPen = sol$ATT
    
    # B. lambda = .01
    sol = regsynth(X0,X1,Y0,Y1,V,.01)
    Synth.fixed = sol$ATT
    
    # C. lambda = lambdaopt
    solpath = regsynthpath(X0,X1,Y0t,Y1t,V,lambda,bar=T)
    keeptau = solpath$CATT

    
    # The one that optimizes MSE
    curve.MSE = apply(keeptau^2,1,sum)/n0
    lambda.opt.MSE = min(lambda[which(curve.MSE==min(curve.MSE))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.MSE)
    Synth.opt.MSE = sol$ATT
    
    # The one that optimizes bias
    curve.bias = abs(apply(keeptau,1,sum)/n0)
    lambda.opt.bias = min(lambda[which(curve.bias==min(curve.bias))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.bias)
    Synth.opt.bias = sol$ATT
    
    # The one that optimizes bias + variance
    curve.crit = curve.bias + apply(keeptau,1,var)
    lambda.opt.crit = min(lambda[which(curve.crit==min(curve.crit))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.crit)
    Synth.opt.crit = sol$ATT
    
    # The one that optimizes MAE
    curve.MAE = apply(abs(keeptau),1,sum)/n0
    lambda.opt.MAE = min(lambda[which(curve.MAE==min(curve.MAE))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.MAE)
    Synth.opt.MAE = sol$ATT
    
    print("*** PROGRESS ***")
    print(100*r/R)
    
    ### 4. ATT estimation
    Results[r,] <- c(AggSynth,NN1$ATT,NN5$ATT,
                     NN.opt.MSE,NN.opt.bias,NN.opt.crit,NN.opt.MAE,
                     Synth.NoPen,Synth.fixed,
                     Synth.opt.MSE,Synth.opt.bias,Synth.opt.crit,Synth.opt.MAE)
    setTxtProgressBar(pb, r/R)
  }
  
  close(pb)
  print(Sys.time()-t_start)
  
  ### Compute bias and RMSE
  StatDisplay <- data.frame()
  StatDisplay[1:13,"bias"] = abs(apply(Results,2,mean))
  StatDisplay[1:13,"MSE"] = apply(Results^2,2,mean)
  StatDisplay[1:13,"MAE"] = apply(abs(Results),2,mean)
  row.names(StatDisplay) = c("Aggregate Synth","1NN Matching","5NN Matching",
                             "NN MSE opt","NN bias opt","NN crit opt", "NN MAE opt",
                             "Synth No Pen","Synth fixed",
                             "Synth MSE opt","Synth bias opt","Synth crit opt","Synth MAE opt")
  print(StatDisplay)
  
  fileN = paste("simulations/Exp5_PolyDGP/Paneloutput_n",n1,",p",p,",delta",delta,".txt",sep="")
  
  print.xtable(xtable(StatDisplay, digits=3),type="latex",file=fileN)
  write(c(paste("Nb. treated:",n1),
          paste("Nb. controls:",n0),
          paste("Nb. covariates:",p),
          paste("Degree poly:",delta),
          paste("Nb. replications:",R),
          paste(Sys.time())), fileN, append=TRUE)
  
}