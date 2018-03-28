### Setup of Exp 5
### Modified: 13/03/2018

Exp5_Poly_setup <- function(R=1000,n1=100,n0,p,delta){
  Results = matrix(ncol=12, nrow=R)
  t_start = Sys.time()
  pb = txtProgressBar(style = 3)
  
  for(r in 1:R){
    ### 0. Generate data
    data = PolyDGP(n1,n0,p,delta)
    X = data$X; y = data$y; d = data$d
  
    X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))
    Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)
    
    ### 1. Synthetic Control on mean of treated
    M = matrix(apply(X1,1,mean), ncol=1)
    AggSC = wATT(y,d,wsol(X0,M,V))
    
    ### Splitting the sample for cross-validation
    split = runif(n0)
    allocation = as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = T))  
    
    ### 2. NN Matching
    ## A. K = 1
    NN1 = matchest(X0,X1,Y0,Y1,V,m=1)
    ## B. K = 5
    NN5 = matchest(X0,X1,Y0,Y1,V,m=5)
    ## C. K = Kopt
    keeptauNN = matrix(nrow=10, ncol=n0)
    for(k in 1:K){
      X1k = matrix(X0[,allocation==k], nrow=p)
      X0k = matrix(X0[,allocation!=k], nrow=p)
      Y1k = Y0[allocation==k]
      Y0k = Y0[allocation!=k]
      for(i in 1:10){
        soli = matchest(X0k,X1k,Y0k,Y1k,V,m=i)
        keeptauNN[i,allocation==k] = soli$CATT
      }
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
    
    
    ### 3. Regularized Synthetic Control
    # A. lambda = .1
    sol = regsynth(X0,X1,Y0,Y1,V,.1)
    RSC.fixed = sol$ATT
    
    # B. lambda = lambdaopt
    keeptau = matrix(nrow=length(lambda), ncol=n0)
    for(k in 1:K){
      X1k = matrix(X0[,allocation==k], nrow=p)
      X0k = matrix(X0[,allocation!=k], nrow=p)
      Y1k = Y0[allocation==k]
      Y0k = Y0[allocation!=k]
      solpath = regsynthpath(X0k,X1k,Y0k,Y1k,V,lambda,bar=T)
      keeptau[,allocation==k] = solpath$CATT
    }
    
    # The one that optimizes MSE
    curve.MSE = apply(keeptau^2,1,sum)/n0
    lambda.opt.MSE = min(lambda[which(curve.MSE==min(curve.MSE))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.MSE)
    RSC.opt.MSE = sol$ATT
    
    # The one that optimizes bias
    curve.bias = abs(apply(keeptau,1,sum)/n0)
    lambda.opt.bias = min(lambda[which(curve.bias==min(curve.bias))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.bias)
    RSC.opt.bias = sol$ATT
    
    # The one that optimizes bias + variance
    curve.crit = curve.bias + apply(keeptau,1,var)
    lambda.opt.crit = min(lambda[which(curve.crit==min(curve.crit))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.crit)
    RSC.opt.crit = sol$ATT
    
    # The one that optimizes MAE
    curve.MAE = apply(abs(keeptau),1,sum)/n0
    lambda.opt.MAE = min(lambda[which(curve.MAE==min(curve.MAE))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.MAE)
    RSC.opt.MAE = sol$ATT
    
    print("*** PROGRESS ***")
    print(100*r/R)
    
    ### 4. ATT estimation
    Results[r,] <- c(AggSC,NN1$ATT,NN5$ATT,
                     NN.opt.MSE,NN.opt.bias,NN.opt.crit,NN.opt.MAE,
                     RSC.fixed,
                     RSC.opt.MSE,RSC.opt.bias,RSC.opt.crit,RSC.opt.MAE)
    setTxtProgressBar(pb, r/R)
   }
  
  close(pb)
  print(Sys.time()-t_start)
  
  ### Compute bias and RMSE
  StatDisplay <- data.frame()
  StatDisplay[1:12,"bias"] = abs(apply(Results,2,mean))
  StatDisplay[1:12,"MSE"] = apply(Results^2,2,mean)
  StatDisplay[1:12,"MAE"] = apply(abs(Results),2,mean)
  row.names(StatDisplay) = c("Aggregate Synth","1NN Matching","5NN Matching",
                              "NN MSE opt","NN bias opt","NN crit opt", "NN MAE opt",
                              "PenSynth fixed",
                              "PenSynth MSE opt","PenSynth bias opt","PenSynth crit opt","PenSynth MAE opt")
  print(StatDisplay)
  
  fileN = paste("simulations/Exp5_PolyDGP/output_n",n1,",p",p,",delta",delta,".txt",sep="")
  
  print.xtable(xtable(StatDisplay, digits=3),type="latex",file=fileN)
  write(c(paste("Nb. treated:",n1),
          paste("Nb. controls:",n0),
          paste("Nb. covariates:",p),
          paste("Degree poly:",delta),
          paste("Nb. replications:",R),
          paste(Sys.time())), fileN, append=TRUE)

}