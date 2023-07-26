###############################
# Function for producing rwa results with confidence intervals for logistic regression
###############################
library(boot)
library(purrr)
suppressPackageStartupMessages(library(tidyverse))
library(rwa)
library(broom)
library(knitr)

#Function to get relative weights on the data
log.rwa<-function(df, outcome = outcome, predictors = predictors){
  #singular value decompostion of predictors
  X<- df %>% dplyr::select(all_of(predictors)) %>% 
    scale() 
  X.svd <-  svd(X)
  Q<-X.svd$v
  P<-X.svd$u
  Z<-P%*%t(Q)
  Z.stand<-scale(Z)
  Lambda<-solve(t(Z.stand)%*%Z.stand)%*%t(Z.stand)%*%X #Obtaining Lambda from equation 7 from Johnson (2000) pg 8
  logrfit<-glm(unlist(df[outcome])~Z.stand,family=binomial)
  b<-coef(logrfit)[2:length(coef(logrfit))] #grab coeff/drop intercept
  LpredY<-predict(logrfit,newdata=df,type="response") #get predicted values
  lYhat<-log(LpredY/(1-LpredY))#Creating logit-Y-hat
  stdlYhat<-sd(lYhat)#Getting stdev of logit-Y-hat
  getting.Rsq<-lm(LpredY~unlist(df[outcome]))#Getting R-sq
  Rsq<-summary(getting.Rsq)$r.squared
  beta<-b*((sqrt(Rsq))/stdlYhat)#Computing standardized logistic regression coefficients
  epsilon<-Lambda^2%*%beta^2
  R.sq <-sum(epsilon)
  PropWeights<-(epsilon/R.sq)
  result <-data.frame(predictors, Raw.RelWeight=epsilon, Rescaled.RelWeight=PropWeights)
}

#function to get relative weights in presence of a random variable for tests of significance
log.rwa.rand <- function (df, outcome, predictors) 
{
  X <- df %>% dplyr::select(all_of(predictors)) %>% 
    ##Add Random variable to data set +pipes----
    dplyr::mutate(rand = rnorm(n(),0,1)) %>% 
    scale()  
  X.svd <-  svd(X)
  Q<-X.svd$v
  P<-X.svd$u
  Z<-P%*%t(Q)
  Z.stand<-scale(Z)
  Lambda<-solve(t(Z.stand)%*%Z.stand)%*%t(Z.stand)%*%X #Obtaining Lambda from equation 7 from Johnson (2000) pg 8
  logrfit<-glm(unlist(df[outcome])~Z.stand,family=binomial)
  b<-coef(logrfit)[2:length(coef(logrfit))] #grab coeff/drop intercept
  LpredY<-predict(logrfit,newdata=df,type="response") #get predicted values
  lYhat<-log(LpredY/(1-LpredY))#Creating logit-Y-hat
  stdlYhat<-sd(lYhat)#Getting stdev of logit-Y-hat
  getting.Rsq<-lm(LpredY~unlist(df[outcome]))#Getting R-sq
  Rsq<-summary(getting.Rsq)$r.squared
  beta<-b*((sqrt(Rsq))/stdlYhat)#Computing standardized logistic regression coefficients
  epsilon<-Lambda^2%*%beta^2
  ##reformat weights to be each value minus the last one produced by the random variable ---
  RawWgt <- as.vector(epsilon) #turn into vector for later
  RawWgt <- RawWgt - tail(RawWgt, n=1) #calculate diff
  head(RawWgt, -1) #drop last value
}

#Function to calculate difference between relative weights of two predictors
log.rwa.comp <- function (df, outcome, predictors, focal) 
{ 
  #Create list of comparisons
  comparisons <<- paste(predictors[-(grep(focal, predictors))], focal, sep = " - ") 
  ###
  X <- df %>% dplyr::select(all_of(predictors)) %>% 
    ##puts focal variable in last column
    relocate(all_of(focal), .after = last_col()) %>% 
    scale()  
  X.svd <-  svd(X)
  Q<-X.svd$v
  P<-X.svd$u
  Z<-P%*%t(Q)
  Z.stand<-scale(Z)
  Lambda<-solve(t(Z.stand)%*%Z.stand)%*%t(Z.stand)%*%X #Obtaining Lambda from equation 7 from Johnson (2000) pg 8
  logrfit<-glm(unlist(df[outcome])~Z.stand,family=binomial)
  b<-coef(logrfit)[2:length(coef(logrfit))] #grab coeff/drop intercept
  LpredY<-predict(logrfit,newdata=df,type="response") #get predicted values
  lYhat<-log(LpredY/(1-LpredY))#Creating logit-Y-hat
  stdlYhat<-sd(lYhat)#Getting stdev of logit-Y-hat
  getting.Rsq<-lm(LpredY~unlist(df[outcome]))#Getting R-sq
  Rsq<-summary(getting.Rsq)$r.squared
  beta<-b*((sqrt(Rsq))/stdlYhat)#Computing standardized logistic regression coefficients
  epsilon<-Lambda^2%*%beta^2
  ##reformat weights to be each value minus the last one produced by the focal variable ---
  RawWgt <- as.vector(epsilon) #turn into vector for later
  RawWgt <- RawWgt - tail(RawWgt, n=1) #calculate diff
  head(RawWgt, -1) #drop last value
}

#function for bootsrapping
#build a function to calculate 3 things that I can bootstrap later
#1 raw weights
#2 diff between raw Weight  and random variable
#3 diff between raw Weight for focal variable and random variable
my.boot <- function(data, indices, outcome, predictors, focal = NULL){
  if (compare == "Yes") { 
    dt <- data[indices,]
    c(
      log.rwa(dt,  outcome, predictors)$Raw.RelWeight,
      log.rwa.rand(dt, outcome, predictors),
      log.rwa.comp(dt,  outcome, predictors, focal)
    )
  }
  else if (compare == "No") {
    dt <- data[indices,]
    c(
      log.rwa(dt,  outcome, predictors)$Raw.RelWeight,
      log.rwa.rand(dt, outcome, predictors)
    )
  }
}


#Run rwa on the data
myrwa <- log.rwa(df, outcome = outcome, predictors = predictors)

#Call the function above {my.boot} within the bootstrap function built in R
myBootstrap <- boot(df, my.boot, R=num.bootstraps, outcome = outcome, predictors = predictors, focal = focal) #run the bootstrap

#Get the bca CIs 
ci.results <- tidy(myBootstrap,conf.int=TRUE,conf.method="bca")

#Print Results----
num.predictors <- nrow(myrwa) #calculate the number of predictors for formatting
print(kable(myrwa, caption = "Raw and Rescaled Relative Weights"))
print(kable(cbind(myrwa[1], ci.results[1:num.predictors,c(4,5)]), caption = "CIs around the raw Relative Weights"))
print(kable(cbind(myrwa[1], ci.results[I(num.predictors+1):I(num.predictors*2),c(4,5)]), caption = "CIs around the difference between a Relative Weight for a substantive variable and a random variable. If zero is included in the interval that predictor is not significant"))
if (compare == "Yes") { 
  print(kable(cbind(comparisons, ci.results[I((num.predictors*2)+1):nrow(ci.results),c(4,5)]), caption = "CIs around the difference between two substantive variables. If zero is included in the interval, those two predictors are not significantly different from one another"))
}

