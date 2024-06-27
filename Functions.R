PIP_stratified_bs_robust_models = function(data,B,formula0,formula1,strat,m){
  n = nrow(data)
  pip = vector(length=B)
  mod_largest = rlm(formula1,data,method="MM",maxit=500)
  r=mod_largest$res
  rank=rank(abs(r))
  
  if(strat==8){  
    s1=which(rank<=n/8)
    s2=which(rank>n/8 & rank<=2*n/8)
    s3=which(rank>2*n/8 & rank<=3*n/8)
    s4=which(rank>3*n/8 & rank<=4*n/8)
    s5=which(rank>4*n/8 & rank<=5*n/8)
    s6=which(rank>5*n/8 & rank<=6*n/8)
    s7=which(rank>6*n/8 & rank<=7*n/8)
    s8=which(rank>7*n/8 & rank<=n)
  }
  if(strat==4){  
    s1=which(rank<=n/4)
    s2=which(rank>n/4 & rank<=2*n/4)
    s3=which(rank>2*n/4 & rank<=3*n/4)
    s4=which(rank>3*n/4 & rank<=n)
  }
  if(strat==2){  
    s1=which(rank<=n/2)
    s2=which(rank>n/2 & rank<=2*n/2)
  }
  
  i = 0
  while(i<B){
    bs.length = 0
    while (bs.length < p+1){
      if(strat==8){
        bs=c(sample(s1,m/8,replace=TRUE),sample(s2,m/8,replace=TRUE),sample(s3,m/8,replace=TRUE),sample(s4,m/8,replace=TRUE),sample(s5,m/8,replace=TRUE),sample(s6,m/8,replace=TRUE),sample(s7,m/8,replace=TRUE),sample(s8,m/8,replace=TRUE)) 
        bs.length=dim(table(bs))
        y.st=y[bs]
      }
      if(strat==4){
        bs=c(sample(s1,m/4,replace=TRUE),sample(s2,m/4,replace=TRUE),sample(s3,m/4,replace=TRUE),sample(s4,m/4,replace=TRUE)) 
        bs.length=dim(table(bs))
        y.st=y[bs]
      }   
      if(strat==2){
        bs=c(sample(s1,m/2,replace=TRUE),sample(s2,m/2,replace=TRUE)) 
        bs.length=dim(table(bs))
        y.st=y[bs]
      }   
    }
    
    indices_in = bs
    indices_out = (1:n)[which(!(1:n) %in%indices_in)]
    
    mod0 = tryCatch({rlm(formula0,data=data[indices_in,],method="MM",maxit=500)}, error=function(e){})
    mod1 = tryCatch({rlm(formula1,data=data[indices_in,],method="MM",maxit=500)}, error=function(e){})
    
    if(!is.null(mod0) & !is.null(mod1) & !any(is.na(coef(mod0)))& !any(is.na(coef(mod1)))){ 
      pred0 =  predict(mod0,newdata=data[indices_out,])
      pred1 =  predict(mod1,newdata=data[indices_out,])
      
      i=i+1
      pip[i] = mean(((pred1-data[indices_out,"y"])^2) < ((pred0-data[indices_out,"y"])^2))+0.5*mean(((pred1-data[indices_out,"y"])^2) == ((pred0-data[indices_out,"y"])^2)) 
    }
  }
  return(list("pip_bs" = pip))
}



f_measures = function(data,K,reps,seed=1988,alpha=0.05){
  set.seed(seed)
  
  PIP_cv = c()
  mse0_CV = c()
  mse1_CV = c()
  mae0_CV = c()
  mae1_CV = c()
  
  for(rep in 1:reps){
    
    #Perform K fold cross validation
    
    data <-data[sample(nrow(data)),]
    cvIndex <- createFolds(factor(data$x), K, returnTrain = T)
    
    
    pip_cv = c()
    mse0_CV_sub=c()
    mse1_CV_sub=c()
    mae0_CV_sub=c()
    mae1_CV_sub=c()
    
    for(j in names(cvIndex)){
      trainData = data[cvIndex[[j]],]
      testData = data[-cvIndex[[j]],]
      
      mod0 = lm(y~x,data=trainData)
      mod1 = rlm(y~x,method="MM",data=trainData)
      
      pred0 = predict(mod0,testData)
      pred1 = predict(mod1,testData)
      
      pip_cv = c(pip_cv,mean(((pred1-testData$y)^2) < ((pred0-testData$y)^2)))
      mse0_CV_sub = c(mse0_CV_sub,mean((pred0-testData$y)^2))
      mse1_CV_sub = c(mse1_CV_sub,mean((pred1-testData$y)^2))
      
      mae0_CV_sub=c(mae0_CV_sub,median(abs(pred0-testData$y)))
      mae1_CV_sub=c(mae1_CV_sub,median(abs(pred1-testData$y)))
    }
    PIP_cv = c(PIP_cv,mean(pip_cv))
    
    mse0_CV = c(mse0_CV,mean(mse0_CV_sub))
    mse1_CV = c(mse1_CV,mean(mse1_CV_sub))
    
    mae0_CV = c(mae0_CV,mean(mae0_CV_sub))
    mae1_CV = c(mae1_CV,mean(mae1_CV_sub))
  }
  return(list("PIP_cv"=mean(PIP_cv),"PIP_cv_lower"=quantile(PIP_cv,alpha),"PIP_cv_upper" = quantile(PIP_cv,1-alpha),"mse0"=mean(mse0_CV),"mse1"=mean(mse1_CV),"mae0"=mean(mae0_CV),"mae1"=mean(mae1_CV)))
}


PIP_stratified_bs_gbm = function(data,B,formula0,formula1,strat,m){
  n = nrow(data)
  pip = vector(length=B)
  mod_largest = gbm(formula1,data=data,distribution = "laplace")
  r= data$y - predict(mod_largest,data,n.trees=mod_largest$n.trees)
  rank=rank(abs(r))
  
  if(strat==8){  
    s1=which(rank<=n/8)
    s2=which(rank>n/8 & rank<=2*n/8)
    s3=which(rank>2*n/8 & rank<=3*n/8)
    s4=which(rank>3*n/8 & rank<=4*n/8)
    s5=which(rank>4*n/8 & rank<=5*n/8)
    s6=which(rank>5*n/8 & rank<=6*n/8)
    s7=which(rank>6*n/8 & rank<=7*n/8)
    s8=which(rank>7*n/8 & rank<=n)
  }
  if(strat==4){  
    s1=which(rank<=n/4)
    s2=which(rank>n/4 & rank<=2*n/4)
    s3=which(rank>2*n/4 & rank<=3*n/4)
    s4=which(rank>3*n/4 & rank<=n)
  }
  
  for(i in 1:B){
    bs.length = 0
    while (bs.length < p+1){
      if(strat==8){
        bs=c(sample(s1,m/8,replace=TRUE),sample(s2,m/8,replace=TRUE),sample(s3,m/8,replace=TRUE),sample(s4,m/8,replace=TRUE),sample(s5,m/8,replace=TRUE),sample(s6,m/8,replace=TRUE),sample(s7,m/8,replace=TRUE),sample(s8,m/8,replace=TRUE)) 
        bs.length=dim(table(bs))
        y.st=y[bs]
      }
      if(strat==4){
        bs=c(sample(s1,m/4,replace=TRUE),sample(s2,m/4,replace=TRUE),sample(s3,m/4,replace=TRUE),sample(s4,m/4,replace=TRUE)) 
        bs.length=dim(table(bs))
        y.st=y[bs]
      }   
    }
    
    indices_in = bs
    indices_out = (1:n)[which(!(1:n) %in%indices_in)]
    
    if(formula0 == "y ~ 1"){pred0 = rep(mean(data[indices_in,"y"]),length=length(indices_out))}
    if(formula0 != "y ~ 1"){mod0 = gbm(formula0,data=data[indices_in,],distribution = "laplace",interaction.depth = 2,n.trees=500)
    pred0 =  predict(mod0,newdata=data[indices_out,],n.trees=mod0$n.trees)
    }
    
    mod1 = gbm(formula1,data=data[indices_in,],distribution = "laplace",interaction.depth = 2,n.trees=500)
    pred1 =  predict(mod1,newdata=data[indices_out,],n.trees=mod1$n.trees)
    
    pip[i] = mean(((pred1-data[indices_out,"y"])^2) < ((pred0-data[indices_out,"y"])^2))+0.5*mean(((pred1-data[indices_out,"y"])^2) == ((pred0-data[indices_out,"y"])^2)) 
  }
  return(list("pip_bs" = pip))
}


PIP_stratified_bs_robust_vs_gbm = function(data,B,formula0,formula1,strat,m,cap=50){
  n = nrow(data)
  pip = vector(length=B)
  mod_largest = rlm(formula0,data,method="MM",maxit=500)
  r=mod_largest$res
  rank=rank(abs(r))
  
  mse0 = vector(length=B)
  mse1 = vector(length=B)
  
  tmse0 = vector(length=B)
  tmse1 = vector(length=B)
  
  if(strat==8){  
    s1=which(rank<=n/8)
    s2=which(rank>n/8 & rank<=2*n/8)
    s3=which(rank>2*n/8 & rank<=3*n/8)
    s4=which(rank>3*n/8 & rank<=4*n/8)
    s5=which(rank>4*n/8 & rank<=5*n/8)
    s6=which(rank>5*n/8 & rank<=6*n/8)
    s7=which(rank>6*n/8 & rank<=7*n/8)
    s8=which(rank>7*n/8 & rank<=n)
  }
  if(strat==4){  
    s1=which(rank<=n/4)
    s2=which(rank>n/4 & rank<=2*n/4)
    s3=which(rank>2*n/4 & rank<=3*n/4)
    s4=which(rank>3*n/4 & rank<=n)
  }
  if(strat==2){  
    s1=which(rank<=n/2)
    s2=which(rank>n/2 & rank<=2*n/2)
  }
  
  i = 0
  while(i<B){
    bs.length = 0
    while (bs.length < p+1){
      if(strat==8){
        bs=c(sample(s1,m/8,replace=TRUE),sample(s2,m/8,replace=TRUE),sample(s3,m/8,replace=TRUE),sample(s4,m/8,replace=TRUE),sample(s5,m/8,replace=TRUE),sample(s6,m/8,replace=TRUE),sample(s7,m/8,replace=TRUE),sample(s8,m/8,replace=TRUE)) 
        bs.length=dim(table(bs))
        y.st=y[bs]
      }
      if(strat==4){
        bs=c(sample(s1,m/4,replace=TRUE),sample(s2,m/4,replace=TRUE),sample(s3,m/4,replace=TRUE),sample(s4,m/4,replace=TRUE)) 
        bs.length=dim(table(bs))
        y.st=y[bs]
      }   
      if(strat==2){
        bs=c(sample(s1,m/2,replace=TRUE),sample(s2,m/2,replace=TRUE)) 
        bs.length=dim(table(bs))
        y.st=y[bs]
      }   
    }
    
    indices_in = bs
    indices_out = (1:n)[which(!(1:n) %in%indices_in)]
    
    mod0 = tryCatch({rlm(formula0,data=data[indices_in,],method="MM",maxit=500)}, error=function(e){})
    mod1 = tryCatch({gbm(formula1,data=data[indices_in,],distribution = "laplace",interaction.depth = 2,n.trees=500)}, error=function(e){})
    
    if(!is.null(mod0) & !is.null(mod1) & !any(is.na(coef(mod0)))& !any(is.na(coef(mod1)))){ 
      pred0 =  predict(mod0,newdata=data[indices_out,])
      pred1 =  predict(mod1,newdata=data[indices_out,],n.trees=mod1$n.trees)
      
      i=i+1
      pip[i] = mean(((pred1-data[indices_out,"y"])^2) < ((pred0-data[indices_out,"y"])^2))+0.5*mean(((pred1-data[indices_out,"y"])^2) == ((pred0-data[indices_out,"y"])^2)) 
      
      squares0 = (pred0-data[indices_out,"y"])^2
      squares1 = (pred1-data[indices_out,"y"])^2
      
      mse0[i] = mean(squares0)
      mse1[i] = mean(squares1)
      
      tsquares0 = squares0
      tsquares0[tsquares0>cap] = cap
      
      tsquares1 = squares1
      tsquares1[tsquares1>cap] = cap
      
      tmse0[i] = mean(tsquares0)
      tmse1[i] = mean(tsquares1)
      
    }
  }
  return(list("pip_bs" = pip,"mse0"=mse0,"mse1"=mse1,"tmse0"=tmse0,"tmse1"=tmse1))
}
