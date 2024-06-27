logit = function(x){
  return(log(x/(1-x)))
}
expit = function(x){
  return((exp(x)/(1+exp(x))))
}

PIP_K_cv = function(data,K,type,alpha=0.05){
  if(K>nrow(data)){print("error: K should be less than or equal to n")}
  #Randomly shuffle the data
  yourData<-data[sample(nrow(data)),]
  
  #Create K equally size folds
  folds <- cut(seq(1,nrow(yourData)),breaks=K,labels=FALSE)
  
  #Perform 10 fold cross validation
  pip_cv = c()
  mse0_LOO_sub=c()
  mse1_LOO_sub=c()
  for(i in 1:K){
    #Segement your data by fold using the which() function
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- yourData[testIndexes, ]
    trainData <- yourData[-testIndexes, ]
    #Use the test and train data partitions however you desire...
    
    if(type=="gaussian"){
      mod1 = glm(y ~ x, family="gaussian", data=trainData)
      mod0 = glm(y ~ 1, family="gaussian", data=trainData)
    }
    
    if(type=="poisson"){
      mod1 = glm(y ~ x, family="poisson", data=trainData, maxit = 5000)
      mod0 = glm(y ~ 1, family="poisson", data=trainData, maxit = 5000)
    }
    
    if(type=="gamma"){
      mod1 = glm(y ~ x, Gamma(link = "log"), data=trainData, maxit = 5000,start=c(beta0+0.05,beta1-0.5))
      mod0 = glm(y ~ 1, Gamma(link = "log"), data=trainData, maxit = 5000)
    }
    
    if(type=="binomial"){
      mod1 = glm(y ~ x, family="binomial", data=trainData, maxit = 5000)
      mod0 = glm(y ~ 1, family="binomial", data=trainData, maxit = 5000)
    }
    
    pred0 = predict(mod0,testData,type="response")
    pred1 = predict(mod1,testData,type="response")
    
    pip_cv = c(pip_cv,mean((pred1-testData$y)^2 < (pred0-testData$y)^2))
    mse0_LOO_sub = c(mse0_LOO_sub,mean((pred0-testData$y)^2))
    mse1_LOO_sub = c(mse1_LOO_sub,mean((pred1-testData$y)^2))
  }
  PIP_cv = mean(pip_cv)
  PIP_cv_lower = PIP_cv - qnorm(1-alpha/2)*sqrt(PIP_cv*(1-PIP_cv)/K)
  PIP_cv_upper = PIP_cv + qnorm(1-alpha/2)*sqrt(PIP_cv*(1-PIP_cv)/K)
  
  PIP_cv_mean = mean(pip_cv)
  PIP_cv_var = var(pip_cv)
  
  # PIP_cv_lower2 = expit(logit( PIP_cv_mean) - qnorm(1-alpha/2)*sqrt(PIP_cv_var/(PIP_cv_mean*(PIP_cv_mean))))
  # PIP_cv_upper2 = expit(logit( PIP_cv_mean) + qnorm(1-alpha/2)*sqrt(PIP_cv_var/(PIP_cv_mean*(PIP_cv_mean))))
  
  PIP_cv_lower2 = expit(logit( PIP_cv_mean) - qt(1-alpha/2,K-1)*sqrt(PIP_cv_var/(PIP_cv_mean*(PIP_cv_mean))))
  PIP_cv_upper2 = expit(logit( PIP_cv_mean) + qt(1-alpha/2,K-1)*sqrt(PIP_cv_var/(PIP_cv_mean*(PIP_cv_mean))))
  
  
  mse0_LOO = mean(mse0_LOO_sub)
  mse1_LOO = mean(mse1_LOO_sub)
  
  if(K==nrow(data)){
    return(list("PIP_cv"=PIP_cv,"mse0"=mse0_LOO,"mse1"=mse1_LOO,"PIP_cv_lower"=PIP_cv_lower,"PIP_cv_upper"=PIP_cv_upper))}
  if(K<nrow(data)){
    return(list("PIP_cv"=PIP_cv,"mse0"=mse0_LOO,"mse1"=mse1_LOO,"PIP_cv_lower"=PIP_cv_lower2,"PIP_cv_upper"=PIP_cv_upper2))}
}
