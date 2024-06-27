##------------------
## Data application
##------------------

# Load packages

library(MASS)
library(tidyverse)
library(robustbase)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(gbm)
library(caret)
library(mlbench)

source("R/Functions_paper2.R")

# Example 1: Stack Loss data

data(stackloss)
head(stackloss)
stackloss_use = stackloss%>%mutate(y=stack.loss,X1=1,X2 =Air.Flow,X3=Water.Temp,X4=Acid.Conc.)
stackloss_use = data.frame(stackloss_use %>% dplyr::select(c("y","X1","X2","X3","X4")))

y = stackloss_use$y
X = as.matrix(stackloss_use[,2:5])
dat = data.frame(y,X)

# PIP model selection (cfr. Table 1 in paper The Probability of Improved Prediction as a New Concept for Model Selection in the Presence of Outliers)

formula0 = as.formula(paste("y~",1))
formula1a = as.formula(paste("y~","X2"))
formula1b = as.formula(paste("y~","X3"))
formula1c = as.formula(paste("y~","X4"))
formula2a = as.formula(paste("y~","X2+X3"))
formula2b = as.formula(paste("y~","X2+X4"))
formula2c = as.formula(paste("y~","X3+X4"))
formula3a = as.formula(paste("y~","X2+X3+X4"))

models_stackloss =  c(formula1a,formula1b,formula1c,formula2a,formula2b,formula2c,formula3a)

p=3; m=12
K=500
models = models_stackloss
formula0 = as.formula(paste("y~",1))
pips = c()

set.seed(1988)
for(mod in models[1:3]){
  pips = c(pips,mean(PIP_stratified_bs_robust_models(dat,K,formula0,mod,4,m)$pip_bs))
}
print(pips)

if(max(pips,na.rm=TRUE)>0.5){
  new_mod0 = models[1:3][which.max(pips)][[1]]

  comparators = models[4:6][unlist(lapply(lapply(models[4:6], function(x) attr(terms(x), "term.labels")),function(x) attr(terms(new_mod0), "term.labels") %in%x))]
  pips = c()
  for(mod in comparators){
    pips = c(pips,mean(PIP_stratified_bs_robust_models(dat,K,new_mod0,mod,4,m)$pip_bs))
  }
  print(pips)

  if(max(pips,na.rm=TRUE)>0.5){
    new_mod0 = comparators[which.max(pips)][[1]]

    comparator = models[7][[1]]
    pips = mean(PIP_stratified_bs_robust_models(dat,K,new_mod0,comparator,4,m)$pip_bs)
    if(pips>0.5) new_mod0 = comparator
  }}

print(pips)

# The final selected model y ~ X2 + X3
selected = new_mod0
print(selected)

# Application 2: Los Angeles Ozone Pollution Data

data("Ozone")
ozone_use = Ozone%>%mutate(y=V4,X1=V8,X2 =V10,X3=V11,X4=V13,
                           X5=V5 ,X6=V7, X7= V12,X8=V6)%>%
  dplyr::select(all_of(c("y","X1","X2","X3","X4","X5","X6","X7","X8")))%>%drop_na()

# 1	Month: 1 = January, ..., 12 = December
# 2	Day of month
# 3	Day of week: 1 = Monday, ..., 7 = Sunday
# 4	Daily maximum one-hour-average ozone reading
# 5	500 millibar pressure height (m) measured at Vandenberg AFB
# 6	Wind speed (mph) at Los Angeles International Airport (LAX)
# 7	Humidity (%) at LAX
# 8	Temperature (degrees F) measured at Sandburg, CA
# 9	Temperature (degrees F) measured at El Monte, CA
# 10	Inversion base height (feet) at LAX
# 11	Pressure gradient (mm Hg) from LAX to Daggett, CA
# 12	Inversion base temperature (degrees F) at LAX
# 13	Visibility (miles) measured at LAX

# Selected variables CSDA Salibian-Barrera & Van Aelst (2008) shown in design below

# 7 predictoren: 1 12 15 27 29 33 43
# 10 predictoren: 1  8 12 14 17 23 34 40 41 43
# 23 predictoren:  1 4  5  6  7  8  9 11 12 13 14 15 17 24 29 31 32 33 39 40 41 42 43

design = read.table("R/Ozone_design.txt",sep='\t',header=TRUE)
vars = names(design)[2:46]

mod7_vars = vars[c(12, 15, 27, 29, 33, 43)]
mod10_vars = vars[c(8 ,12, 14, 17 ,23, 34, 40, 41 ,43)]
mod23_vars = vars[c(4,  5,  6,  7,  8,  9, 11, 12, 13, 14, 15, 17, 24, 29, 31, 32, 33, 39, 40, 41, 42, 43)]

design$y=ozone_use$y

# Models and predictions from Salibian-Barrera & Van Aelst (2008)
mod7_ozone = rlm(as.formula(paste0("y~",paste(mod7_vars,collapse="+"))),data=design ,method="MM",maxit=1000)
mod10_ozone = rlm(as.formula(paste0("y~",paste(mod10_vars,collapse="+"))),data=design ,method="MM",maxit=1000)
mod23_ozone = rlm(as.formula(paste0("y~",paste(mod23_vars,collapse="+"))),data=design ,method="MM",maxit=1000)

preds_mod7 = predict(mod7_ozone,design)
preds_mod10 = predict(mod10_ozone,design)
preds_mod23 = predict(mod23_ozone,design)


par(mfrow=c(2,2))
plot(preds_mod7,design$y)
title(paste0("Pearson cor: ",cor(preds_mod7,design$y)))
plot(preds_mod10,design$y)
title(paste0("Pearson cor: ",cor(preds_mod10,design$y)))
plot(preds_mod23,design$y)
title(paste0("Pearson cor: ",cor(preds_mod23,design$y)))


# Full linear model
par(mfrow=c(1,2))

ozone_use_lm = ozone_use%>%mutate(X1sq = X1^2,X2sq = X2^2,X3sq = X3^2,X4sq = X4^2,X5sq = X5^2,X6sq = X6^2,
                                  X7sq = X7^2,X8sq = X8^2)
full_linear_mod_ozone = rlm(y~(X1+X2+X3+X4+X5+X6+X7+X8)^2 + X1sq+ X2sq+ X3sq+ X4sq+ X5sq+ X6sq+ X7sq+ X8sq,data=ozone_use_lm ,method="MM")
preds_lm = predict(full_linear_mod_ozone,ozone_use_lm,n.trees=mod_ozone$n.trees)
plot(preds_lm,ozone_use$y)
title(paste0("Pearson cor: ",cor(preds_lm,ozone_use$y)))
scale = mad(ozone_use_lm$y - predict(full_linear_mod_ozone ))

mod_ozone = gbm(y~.,data=ozone_use,distribution = "laplace",interaction.depth = 2,n.trees=5000)
preds_gbm = predict(mod_ozone,ozone_use,n.trees=mod_ozone$n.trees)
plot(preds_gbm,ozone_use$y,col="red")
title(paste0("Pearson cor: ",cor(preds_gbm,ozone_use$y)))


# PIP model selection (cfr. Table 2 in paper The Probability of Improved Prediction as a New Concept for Model Selection in the Presence of Outliers)

# Stage 1 models

formula0 = as.formula(paste("y~",1))
formula1a = as.formula(paste("y~","X1"))
formula1b = as.formula(paste("y~","X2"))
formula1c = as.formula(paste("y~","X3"))
formula1d = as.formula(paste("y~","X4"))
formula1e = as.formula(paste("y~","X5"))
formula1f = as.formula(paste("y~","X6"))
formula1g = as.formula(paste("y~","X7"))
formula1h = as.formula(paste("y~","X8"))

K=500
pips = c()

set.seed(1988)
for(mod in c(formula1a,formula1b,formula1c,formula1d,formula1e,formula1f,formula1g,formula1h)){
  pips = c(pips,mean(PIP_stratified_bs_gbm(ozone_use,K,formula0,mod,4,200)$pip_bs))
}
print(pips) # X1 added


# Stage 2 models

formula2b = as.formula(paste("y~","X1+X2"))
formula2c = as.formula(paste("y~","X1+X3"))
formula2d = as.formula(paste("y~","X1+X4"))
formula2e = as.formula(paste("y~","X1+X5"))
formula2f = as.formula(paste("y~","X1+X6"))
formula2g = as.formula(paste("y~","X1+X7"))
formula2h = as.formula(paste("y~","X1+X8"))

K=500
pips = c()

#set.seed(1988)
for(mod in c(formula2b,formula2c,formula2d,formula2e,formula2f,formula2g,formula2h)){
  pips = c(pips,mean(PIP_stratified_bs_gbm(ozone_use,K,formula1a,mod,4,200)$pip_bs))
}
print(pips) # X6 added

# Stage 3 models
formula3b = as.formula(paste("y~","X1+X6+X2"))
formula3c = as.formula(paste("y~","X1+X6+X3"))
formula3d = as.formula(paste("y~","X1+X6+X4"))
formula3e = as.formula(paste("y~","X1+X6+X5"))
formula3g = as.formula(paste("y~","X1+X6+X7"))
formula3h = as.formula(paste("y~","X1+X6+X8"))


K=500
pips = c()

#set.seed(1988)
for(mod in c(formula3b,formula3c,formula3d,formula3e,formula3g,formula3h)){
  pips = c(pips,mean(PIP_stratified_bs_gbm(ozone_use,K,formula2f,mod,4,200)$pip_bs))
}
print(pips) # X7 added


# Stage 4 models
formula4b = as.formula(paste("y~","X1+X6+X7+X2"))
formula4c = as.formula(paste("y~","X1+X6+X7+X3"))
formula4d = as.formula(paste("y~","X1+X6+X7+X4"))
formula4e = as.formula(paste("y~","X1+X6+X7+X5"))
formula4h = as.formula(paste("y~","X1+X6+X7+X8"))


K=500
pips = c()

#set.seed(1988)
for(mod in c(formula4b,formula4c,formula4d,formula4e,formula4h)){
  pips = c(pips,mean(PIP_stratified_bs_gbm(ozone_use,K,formula3g,mod,4,200)$pip_bs))
}
print(pips) # X2 added


# Stage 5 models
formula5c = as.formula(paste("y~","X1+X6+X7+X2+X3"))
formula5d = as.formula(paste("y~","X1+X6+X7+X2+X4"))
formula5e = as.formula(paste("y~","X1+X6+X7+X2+X5"))
formula5h = as.formula(paste("y~","X1+X6+X7+X2+X8"))

K=500
pips = c()

#set.seed(1988)
for(mod in c(formula5c,formula5d,formula5e,formula5h)){
  pips = c(pips,mean(PIP_stratified_bs_gbm(ozone_use,K,formula4b,mod,4,200)$pip_bs))
}
print(pips) # X4 added

# Stage 6 models
formula6c = as.formula(paste("y~","X1+X6+X7+X2+X4+X3"))
formula6e = as.formula(paste("y~","X1+X6+X7+X2+X4+X5"))
formula6h = as.formula(paste("y~","X1+X6+X7+X2+X4+X8"))

K=500
pips = c()

#set.seed(1988)
for(mod in c(formula6c,formula6e,formula6h)){
  pips = c(pips,mean(PIP_stratified_bs_gbm(ozone_use,K,formula5d,mod,4,200)$pip_bs))
}
print(pips) # X3 added

# Stage 7 models
formula7e = as.formula(paste("y~","X1+X6+X7+X2+X4+X3+X5"))
formula7h = as.formula(paste("y~","X1+X6+X7+X2+X4+X3+X8"))

K=500
pips = c()

#set.seed(1988)
for(mod in c(formula7e,formula7h)){
  pips = c(pips,mean(PIP_stratified_bs_gbm(ozone_use,K,formula6c,mod,4,200)$pip_bs))
}
print(pips) # no further additions

## Selected model
par(mfrow=c(1,2))
set.seed(1988)
mod_ozone = gbm(y~X1+X2+X3+X4+X6+X7,data=ozone_use,distribution = "laplace",interaction.depth = 2,n.trees=5000)
preds_gbm = predict(mod_ozone,ozone_use,n.trees=mod_ozone$n.trees)
plot(preds_gbm,ozone_use$y,col="red")
title(paste0("Pearson cor: ",cor(preds_gbm,ozone_use$y)))
#summary.gbm(mod_ozone)

set.seed(1988)
mod_ozone_compare = gbm(as.formula(paste0("X~",paste(vars[c(2,7,8,3,5,4)],collapse = "+"))),data=design,distribution = "laplace",interaction.depth = 2,n.trees=5000)
preds_gbm_compare = predict(mod_ozone_compare,design,n.trees=mod_ozone_compare$n.trees)
plot(preds_gbm_compare,design$x,col="red")
title(paste0("Pearson cor: ",cor(preds_gbm_compare,design$X)))


# Compare full linear model with selected GBM

set.seed(1988)
comp_linear_gbm = PIP_stratified_bs_robust_vs_gbm(ozone_use_lm,K,as.formula("y~(X1+X2+X3+X4+X5+X6+X7+X8)^2 + X1sq+ X2sq+ X3sq+ X4sq+ X5sq+ X6sq+ X7sq+ X8sq"),as.formula("y ~ X1 + X6 + X7 + X2 + X4 + X3"),4,200,cap=20)
mean(comp_linear_gbm$pip_bs) # GBM model is better in 57% of made predictions
mean(comp_linear_gbm$mse0)
mean(comp_linear_gbm$mse1)
mean(comp_linear_gbm$tmse0)
mean(comp_linear_gbm$tmse1)


