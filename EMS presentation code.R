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

source("Functions.R")

# Example 1: Stack Loss data

data(stackloss)
head(stackloss)
stackloss_use = stackloss%>%mutate(y=stack.loss,X1=1,X2 =Air.Flow,X3=Water.Temp,X4=Acid.Conc.)
stackloss_use = data.frame(stackloss_use %>% dplyr::select(c("y","X1","X2","X3","X4")))

y = stackloss_use$y
X = as.matrix(stackloss_use[,2:5])
dat = data.frame(y,X)

# Visualisations and summary measures presentation EMS 

names(dat)[2:5] = paste0('X',c(0,1,2,3))
dat$id = rownames(dat)

dat  = dat %>% mutate(Outlier = dplyr::case_when(!(dat$id %in% c(1,3,4,21)) ~ "No", 
                                                 dat$id %in% c(1,3,4,21) ~ "Yes"))

p1b = ggplot(data=dat, aes(x= X1, y= y))+
  geom_point(aes(shape=Outlier),size=3)+
  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                               data,
                                                               method="MM"),
              fullrange=TRUE, aes(colour="Robust MM")) +
  stat_smooth(method="lm", aes(colour="Least Squares")) +
  scale_colour_manual(name="Estimators", values=c("blue", "red"))+ theme(legend.title=element_text(size=12),legend.text=element_text(size=12))


p2b =ggplot(data=dat, aes(x= X2, y= y))+
  geom_point(aes(shape=Outlier),size=3)+
  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                               data,
                                                               method="MM"),
              fullrange=TRUE, aes(colour="Robust MM")) +
  stat_smooth(method="lm", aes(colour="Least Squares")) +
  scale_colour_manual(name="Estimators", values=c("blue", "red"))+ theme(legend.title=element_text(size=12),legend.text=element_text(size=12))



p3b =ggplot(data=dat, aes(x= X3, y= y))+
  geom_point(aes(shape=Outlier),size=3)+
  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                               data,
                                                               method="MM"),
              fullrange=TRUE, aes(colour="Robust MM")) +
  stat_smooth(method="lm", aes(colour="Least Squares")) +
  scale_colour_manual(name="Estimators", values=c("blue", "red"))+ theme(legend.title=element_text(size=12),legend.text=element_text(size=12))


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1b)

grid.arrange(arrangeGrob(p1b + theme(legend.position="none"),
                         p2b + theme(legend.position="none"),
                         p3b + theme(legend.position="none"),
                         mylegend,
                         nrow=2))



mod_robust = rlm(y~X3,data=dat,method="MM",maxit=500)
mod_ls = lm(y~X3,data=dat)

mean((dat$y-predict(mod_robust))^2)
mean((dat$y-predict(mod_ls))^2)

median(abs(dat$y-predict(mod_robust)))
median(abs(dat$y-predict(mod_ls)))

## Comparing Least squares and robust MM modes with PIP, MSE and MAE

dat_use = dat[,c("y","X3")]
names(dat_use)[2] = "x"
f_measures(dat_use,5,50)
