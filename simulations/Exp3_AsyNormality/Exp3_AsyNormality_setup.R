### Setup of the Exp3
### Goal I: show that the estimator is asymptotically normal
### Goal II: show that selecting lambda in a specific way can help reduce the bias
### Started: 21/02/2018
### Idea is from Xavier

Exp3_setup <- function(R=1000,n1=30,n0=30,p=50,K=5){
  Results = matrix(ncol=5, nrow=R)
  t_start = Sys.time()
  pb = txtProgressBar(style = 3)
  
  for(r in 1:R){
    ### 0. Generate data
    data = matchDGP_fixedT(n1=n1,n0=n0,p=p)
    X = data$X; y = data$y; d = data$d
  
    X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))
    Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)
    
    ### Splitting the sample for cross-validation
    uu=0 # reshuffle groups until no empty group
    while(uu==0){
      allocation = sample(1:K,n0,replace=T)
      uu=min(mapply(function(x) sum(allocation==x),1:K))
    }
    
    ### 1. COmpute different estimators
    # A. lambda = .1
    sol = regsynth(X0,X1,Y0,Y1,V,.1)
    RSC.fixed = sol$ATT
    
    # B. lambda = lambdaopt
    keeptau = matrix(nrow=length(lambda), ncol=length(Y0))
    for(k in 1:K){
      X1k = as.matrix(X0[,allocation==k])
      X0k = as.matrix(X0[,allocation!=k])
      Y1k = Y0[allocation==k]
      Y0k = Y0[allocation!=k]
      solpath = regsynthpath(X0k,X1k,Y0k,Y1k,V,lambda)
      keeptau[,allocation==k] = solpath$CATT
    }
    
    # The one that optimizes RMSE
    curve.RMSE = apply(keeptau^2,1,sum)/n0
    lambda.opt.RMSE = min(lambda[which(curve.RMSE==min(curve.RMSE))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.RMSE)
    RSC.opt.RMSE = sol$ATT
    
    # The one that optimizes bias
    curve.bias = abs(apply(keeptau,1,sum)/n0)
    lambda.opt.bias = min(lambda[which(curve.bias==min(curve.bias))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.bias)
    RSC.opt.bias = sol$ATT
    
    # The one that optimizes bias + variance
    curve.crit = curve.bias + apply(keeptau,1,sd)
    lambda.opt.crit = min(lambda[which(curve.crit==min(curve.crit))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.crit)
    RSC.opt.crit = sol$ATT
    
    # The one that optimizes MAE
    curve.mae = apply(abs(keeptau),1,sum)/n0
    lambda.opt.mae = min(lambda[which(curve.mae==min(curve.mae))])
    sol = regsynth(X0,X1,Y0,Y1,V,lambda.opt.mae)
    RSC.opt.mae = sol$ATT
    
    print("*** PROGRESS ***")
    print(100*r/R)
    
    ### 2. Collect ATT estimate
    Results[r,] <- c(RSC.fixed,RSC.opt.RMSE,RSC.opt.bias,RSC.opt.crit,RSC.opt.mae)
    setTxtProgressBar(pb, r/R)
   }
  
  close(pb)
  print(Sys.time()-t_start)
  
  ### Compute bias and RMSE
  StatDisplay <- data.frame()
  StatDisplay[1:5,"bias"] = abs(apply(Results,2,mean))
  StatDisplay[1:5,"RMSE"] = sqrt(apply(Results^2,2,mean))
  StatDisplay[1:5,"MAE"] = apply(abs(Results),2,mean)
  row.names(StatDisplay) = c("Penalized Synth fixed","Penalized Synth RMSE opt","Penalized Synth bias opt","Penalized Synth crit opt","Penalized Synth MAE opt")
  print(StatDisplay)
  
  fileN = paste("simulations/Exp3_AsyNormality/output_n",n1+n0,",p",p,".txt",sep="")
  
  print.xtable(xtable(StatDisplay, digits=3),type="latex",file=fileN)
  write(c(paste("Nb. treated:",n1),
          paste("Nb. controls:",n0),
          paste("Nb. covariates:",p),
          paste("Nb. replications:",R),
          paste(Sys.time())), fileN, append=TRUE)
  
  return(Results)
}