
library(dplyr)
library(rms)
library(tidyverse)
library(xtable)
library(survivalROC)
library(survival)
library(survcomp)
library(survAUC)
simulate_data <- function(seed) {
  library(truncdist)
  #Data<-read.csv("FULLDATA1.csv")
  Data<-read.csv("DataNONSIG.csv")
  split<-split(Data, f =Data$ID>=160) 
  
  train <-split$`FALSE`
  test <-split$`TRUE`
  
  
  dftrain <-train
  dftest <- test
  
  surv_train <- dftrain$Stop
  cen_train <- dftrain$Event
  
  surv_test <- dftest$Stop
  cen_test <- dftest$Event
  
  ## Cox Model
  s_extract <- function(summary_entry){
    separate(as_tibble(summary_entry),
             sep = ":",
             col = value, 
             remove = FALSE, 
             into = c("bad", "good"))[[3]] %>% 
      as.numeric() 
  }
  
  Deep1 = list()
  Deep2 = list()
  Deep3 = list()
  Deep4 = list()
  Deep5 = list()
  CoxM  = list()
  coxmodel<-cph(Surv(time = Start, time2 = Stop, event = Event) ~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+cluster(ID)
                , data =dftrain, x=T, y=T) 
  valecox<-summary(validate(coxmodel,method="boot", B=50,dxy=TRUE))
  valecox<-s_extract(valecox[4,3])
  ###
  qt<-c(2,20,30,50,60)
  pseudo <- pseudosurv(time=dftrain$Stop,event=dftrain$Event,tmax=qt)
  
  btrain <- NULL
  for(it in 1:length(pseudo$time)){
    btrain <- rbind(btrain,cbind(dftrain,pseudo = pseudo$pseudo[,it],
                                 tpseudo = pseudo$time[it],id=1:nrow(dftrain)))
  }
  dftrain <- btrain
  x_train <- dftrain[,c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]
  x_test <- dftest[,c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]
  
  smatrix=model.matrix(~as.factor(btrain$tpseudo)+0)
  
  #create input predictors 
  x_train.all <- cbind(x_train, smatrix)
  x_train.all<-as.matrix(x_train.all)
  # Isolate the Y
  y_train.all <- btrain$pseudo
  
  # Train the model
  
  model = pseudoDNN.train(x_train.all, y_train.all)
  
  ## predict
 
  x_test.all=do.call(rbind, replicate(length(qt), as.matrix(x_test), simplify=FALSE))
  s_test=rep(qt,each=nrow(x_test))
  s_test
  smatrix.test=model.matrix(~as.factor(s_test)+0)
  x_test.all=cbind(x_test.all,smatrix.test)
  ypred.con <- pseudoDNN.predict(model, x_test.all)
  ypred.con
  # obtain the marginal survival probability by multiple series of conditional probabilities
  ypred.con <- matrix(ypred.con, nrow=nrow(x_test))
  ypred <- lapply(1:length(qt), function(i) apply(ypred.con[,1:i, drop=FALSE], 1, prod))
  surv_prob <- Reduce(cbind, ypred)
  c_1<-concordance.index(x=1-surv_prob[,1], surv.time=surv_test, surv.event=cen_test, method="noether")$c.index
  c_2<-concordance.index(x=1-surv_prob[,2], surv.time=surv_test, surv.event=cen_test, method="noether")$c.index
  c_3<-concordance.index(x=1-surv_prob[,3], surv.time=surv_test, surv.event=cen_test, method="noether")$c.index
  c_4<-concordance.index(x=1-surv_prob[,4], surv.time=surv_test, surv.event=cen_test, method="noether")$c.index
  c_5<-concordance.index(x=1-surv_prob[,5], surv.time=surv_test, surv.event=cen_test, method="noether")$c.index
  valcox <-(valecox)/2+0.5 
  Deep1<-c(Deep1, c_1)
  Deep2<-c(Deep2, c_2)
  Deep3<-c(Deep3, c_3)
  Deep4<-c(Deep4, c_4)
  Deep5<-c(Deep5, c_5)
  CoxM<-c(CoxM,valcox)
  c(Deep1,Deep2,Deep3,Deep4,Deep5, CoxM,seed)
}

simulate_trial <- function(n_sims, seed) {
  set.seed(seed)
  results <- replicate(n_sims, simulate_data(seed))
  data.frame(t(results))
}
build_design_matrix <- function(initial_seed, num_seeds) {
  set.seed(initial_seed)
  seeds <- sample.int(100000, num_seeds)
  design <- expand.grid(
    seed = seeds
  )
}
n_sims <-100 # At each sample size we repeat n_sims times
setup_seed <-42  # 
n_seeds <- 5# trying the experiment at several number of seeds
design <- build_design_matrix(setup_seed, n_seeds)
results <- design %>%
  rowwise() %>%do(simulate_trial(n_sims,.$seed))
Results<-cbind(results)
colnames(Results)<-c("Deep1","Deep2","Deep3","Deep4","Deep5","CoxM","Seed")
write_csv(Results, file="New_data_setSimu200NonSIG.csv")

coxmodel<-coxph(Surv(time = Start, time2 = Stop, event = Event) ~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+cluster(ID)
              , data =dftrain) 
