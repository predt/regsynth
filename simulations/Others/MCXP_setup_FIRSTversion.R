### Setup of the MC experiment

MCXP_setup <- function(R=1000,n=100,p=50){
  Results <- matrix(ncol=10, nrow=R)
  t_start <- Sys.time()
  pb <- txtProgressBar(style = 3)
  
  for(r in 1:R){
    ### 0. Generate data
    data = matchDGP(n=n,p=p,Ry=.5,Rd=.2)
    X = data$X; y = data$y; d = data$d
  
    X0 = t(X[d==0,]); X1 = t(X[d==1,]); V = diag(ncol(X))
    Y0 = y[d==0]; Y1 = y[d==1]; n0 = sum(1-d)
    
    ### 1. Synthetic Control on mean of treated
    M = matrix(apply(X1,1,mean), ncol=1)
    AggSC = wATT(y,d,wsol(X0,M,V))
    
    ### Splitting the sample for cross-validation
    uu=0 # reshuffle groups until no empty group
    while(uu==0){
      allocation = sample(1:K,n0,replace=T)
      uu=min(mapply(function(x) sum(allocation==x),1:K))
    }
    
    ### 2. NN Matching
    ## A. K = 1
    NN1 = matchest(X0,X1,Y0,Y1,V,m=1)
    ## B. K = 5
    NN5 = matchest(X0,X1,Y0,Y1,V,m=5)
    ## C. K = Kopt
    keeptauNN = matrix(nrow=10, ncol=length(Y0))
    for(k in 1:K){
      X1k = as.matrix(X0[,allocation==k])
      X0k = as.matrix(X0[,allocation!=k])
      Y1k = Y0[allocation==k]
      Y0k = Y0[allocation!=k]
      for(i in 1:10){
        soli = matchest(X0k,X1k,Y0k,Y1k,V,m=i)
        keeptauNN[i,allocation==k] = soli$CATT
      }
    }
    
    # The one that optimizes RMSE
    curve.RMSE = apply(keeptauNN^2,1,sum)/n0
    m.opt.RMSE = min(which(curve.RMSE==min(curve.RMSE)))
    sol = matchest(X0k,X1k,Y0k,Y1k,V,m=m.opt.RMSE)
    NN.opt.RMSE = sol$ATT
    
    # The one that optimizes bias
    curve.bias = abs(apply(keeptauNN,1,sum)/n0)
    m.opt.bias = min(which(curve.bias==min(curve.bias)))
    sol = matchest(X0k,X1k,Y0k,Y1k,V,m=m.opt.bias)
    NN.opt.bias = sol$ATT
    
    # The one that optimizes bias + variance
    curve.crit = curve.bias + apply(keeptauNN,1,sd)
    m.opt.crit = min(which(curve.crit==min(curve.crit)))
    sol = matchest(X0k,X1k,Y0k,Y1k,V,m=m.opt.crit)
    NN.opt.crit = sol$ATT
    
    
    ### 3. Regularized Synthetic Control
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
    
    print("*** PROGRESS ***")
    print(100*r/R)
    
    ### 4. ATT estimation
    Results[r,] <- c(AggSC,NN1$ATT,NN5$ATT,
                     NN.opt.RMSE,NN.opt.bias,NN.opt.crit,
                     RSC.fixed,
                     RSC.opt.RMSE,RSC.opt.bias,RSC.opt.crit)
    setTxtProgressBar(pb, r/R)
   }
  
  close(pb)
  print(Sys.time()-t_start)
  
  ### Compute bias and RMSE
  StatDisplay <- data.frame()
  StatDisplay[1:10,"bias"] <- abs(apply(Results,2,mean))
  StatDisplay[1:10,"RMSE"]  <- sqrt(apply(Results^2,2,mean))
  if(R<5000) StatDisplay[1:10,"ShapiroTest"]  <- apply(Results,2, function(x) shapiro.test(x)$p.value)
  row.names(StatDisplay) <- c("Aggregate Synth","1NN Matching","5NN Matching",
                              "NN RMSE opt","NN bias opt","NN crit opt",
                              "Penalized Synth fixed",
                              "Penalized Synth RMSE opt","Penalized Synth bias opt","Penalized Synth crit opt")
  print(StatDisplay)
  
  fileN = paste("simulations/output_n",n,",p",p,".txt",sep="")
  
  print.xtable(xtable(StatDisplay, digits=3),type="latex",file=fileN)
  write(c(paste("Nb. observations:",n),
          paste("Nb. covariates:",p),
          paste("Nb. replications:",R),
          paste(Sys.time())), fileN, append=TRUE)

}