
#################################################
# Code for replication of analysis in Chapter 3 #
#################################################

#########################
## CASE STUDY ANALYSIS ##
#########################

setwd("/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github")

library(Surrogate)    # load the Surrogate library
library(tidyr)        # load the tidyr library (function: gather)
library(nlme)         # load the nlme library (function: lme)
library(lme4)         # load the lme4 library (function: lmer)
library("openxlsx")   # load the openxlsx library (function: write.xlsx)

data("Schizo")        # load the schizophrenia data set
head(Schizo)          # show the first rows of the data set

schizo_panss_bprs <- subset(Schizo, select = -c(CGI,PANSS_Bin,
                                                BPRS_Bin,CGI_Bin))
schizo_pb_long <- gather(schizo_panss_bprs, endpoint, outcome,
                         BPRS:PANSS, factor_key = TRUE)     # transpose the data to long format
head(schizo_pb_long)

#########################
### Analysis with lme ###
#########################

lme1 <- lme(outcome~-1+as.factor(endpoint)+as.factor(endpoint):Treat, 
            random=~-1+as.factor(endpoint)+as.factor(endpoint):Treat|as.factor(InvestId),
            correlation=nlme::corSymm(form=~1|as.factor(InvestId)/as.factor(Id)), 
            data=schizo_pb_long,
            weights=nlme::varIdent(form=~1|endpoint),
            na.action=na.exclude,
            control=list(maxIter=100, msMaxIter=500))
D <- lme1$modelStruct$reStruct$`as.factor(InvestId)`
A <- matrix(c(D[4,1], D[4,3]), nrow=2, ncol=1)
B <- matrix(c(D[1,1], D[3,1], D[3,1], D[3,3]), nrow=2, ncol=2)
C <- as.matrix(c(D[4,4]))
R2trial <- (t(A) %*% solve(B) %*% A)/C
ME <- min(eigen(D)$values)
CN <- max(svd(D)$d)/min(svd(D)$d)

##########################
### Analysis with lmer ###
##########################

lmer <- lmer(outcome~-1+as.factor(endpoint)+as.factor(endpoint):Treat + 
               (-1+as.factor(endpoint)+as.factor(endpoint):Treat|InvestId), 
             data=schizo_pb_long)
D <- matrix(summary(lmer)$varcor$InvestId[1:16], ncol=4)
A <- matrix(c(D[4,1], D[4,3]), nrow=2, ncol=1)
B <- matrix(c(D[1,1], D[3,1],
                D[3,1], D[3,3]), nrow=2, ncol=2)
C <- as.matrix(D[4,4])
R2trial <- (t(A) %*% solve(B) %*% A)/C
ME <- min(eigen(D)$values)
CN <- max(svd(D)$d)/min(svd(D)$d)

#####################################
### Analysis with BifixedContCont ###
### Two-stage approach            ###
#####################################

two_stage <- BifixedContCont(Dataset=Schizo, Surr=BPRS, True=PANSS,
                             Treat=Treat, Trial.ID=InvestId, Pat.ID=Id,
                             Weighted=TRUE)
summary(two_stage)

###########################################
### Analysis with PROC MIXED (from SAS) ###
###########################################

data <- read.table("C:/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github/Schizo_ori_pb UN.txt",header = TRUE)
data <- read.table("C:/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github/Schizo_ori_pb FA.txt",header = TRUE)

D <- matrix(c(data$D_11,data$D_21,data$D_31,data$D_41,
              data$D_21,data$D_22,data$D_32,data$D_42,
              data$D_31,data$D_32,data$D_33,data$D_43,
              data$D_41,data$D_42,data$D_43,data$D_44), 
            nrow = 4, ncol = 4)
A <- matrix(c(D[4,1], D[4,3]), nrow=2, ncol=1)
B <- matrix(c(D[1,1], D[3,1], D[3,1], D[3,3]), nrow=2, ncol=2)
C <- as.matrix(c(D[4,4]))
R2trial <- (t(A) %*% solve(B) %*% A)/C
ME <- min(eigen(D)$values)
CN <- max(svd(D)$d)/min(svd(D)$d)


######################
## SIMULATION STUDY ##
######################

##########################################
### Data generation                    ###
### Code for simulation setting 1 only ###
##########################################

beta <- c(45, 50, 30, 50)*10     # mu_S, mu_T, alpha, beta
gamma <- c(.1,1)                 # gamma (ratio of between-cluster to within-cluster variability)
N.trial <- c(5,10,20)            # number of clusters 
ni.trial <- 20                   # cluster size (change when necessary)
runs <- 500                      # number of data sets
all_results <- NULL

for (c in 1:length(gamma)){
  for (b in 1:length(ni.trial)){
    for (a in 1:length(N.trial)){
      
      gamma_value <- gamma[c]  
      N.trial_value <- N.trial[a]  
      ni.trial_value <- ni.trial[b]
      
      setting <- paste("N.trial = ", N.trial_value, ", ni.trial = ", ni.trial_value, ", gamma = ", gamma_value, sep="")
      setwd("C:/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github")
      dir.create(setting)
      setwd(setting)
      
      for (i in 1:runs){
        
        Sigma_SS.Target <- c(3) 
        Sigma_ST.Target <- c(sqrt(.5*9)) 
        Sigma_TT.Target <- c(3) 
        Sigma <- matrix(c(Sigma_SS.Target, Sigma_ST.Target, Sigma_ST.Target, Sigma_TT.Target), ncol=2)
        Sigma <- Sigma * 100
        Sigma_SS.Target <- Sigma[1, 1]
        Sigma_ST.Target <- Sigma[1, 2]
        Sigma_TT.Target <- Sigma[2, 2]
        R2.Indiv.Target <- Sigma_ST.Target**2 / (Sigma_SS.Target * Sigma_TT.Target)
        
        D <- matrix(rep(0, times=16), ncol=4)
        diag(D) <- c(1, 1, 1, 1)
        D[1,2] <- D[2,1] <- c(.4)
        D[3,4] <- D[4,3] <- c(sqrt(.5))
        D <- D * gamma_value
        D <- D * 1000
        
        A <- matrix(c(D[4,1], D[4,3]), nrow=2, ncol=1)
        B <- matrix(c(D[1,1], D[3,1], D[3,1], D[3,3]), nrow=2, ncol=2)
        C <- as.matrix(D[4,4])
        R2.Trial.Target <- (t(A) %*% solve(B) %*% A)/C
        D_SS.Target <- D[1,1]
        D_TT.Target <- D[2,2]
        D_aa.Target <- D[3,3]
        D_bb.Target <- D[4,4]
        D_TS.Target <- D[2,1]
        D_aS.Target <- D[3,1]
        D_bS.Target <- D[4,1]
        D_aT.Target <- D[3,2]
        D_bT.Target <- D[4,2]
        D_ba.Target <- D[4,3]
        Min.Eigen.D <- min(eigen(D)$values)
        singular  <- svd(D)$d
        Cond.Number.D.Matrix <- max(singular)/min(singular)
        
        all.results <- as.vector(NULL)
        set.seed(i)
        btw.err <- mvrnorm(N.trial_value, c(0, 0, 0, 0), D)   
        wth.err <- mvrnorm(N.trial_value*ni.trial_value, c(0, 0), Sigma)
        
        for (n in 1:N.trial_value){
          Trial <- n
          Z <- sample(x=rep(c(-1, 1), each=(ceiling(ni.trial_value/2))), ni.trial_value, replace=FALSE)  
          Surr <- beta[1] + btw.err[n,1] + (beta[3] + btw.err[n,3])*Z 
          True <- beta[2] + btw.err[n,2] + (beta[4] + btw.err[n,4])*Z  
          results <- cbind(Trial, Z, Surr, True)
          all.results <- rbind(all.results, results)
        }
        
        all.results[,3] <- all.results[,3] +  wth.err[,1]
        all.results[,4] <- all.results[,4] +  wth.err[,2]
        
        Pat.ID <- c(1:c(N.trial_value*ni.trial_value))
        all.results <- cbind(Pat.ID, all.results)
        dataset <- data.frame(all.results)
        names(dataset) <- c("Pat.ID", "Trial.ID", "Treat", "Surr", "True")
        
        write.table(dataset, file=paste("Dataset ", i, ".txt", sep=""), row.names = FALSE)
      }
    }
  }
}

#########################################################
### Analysis with lme                                 ###
### (with description of the within-group correlation ### 
### and heteroscedasticity)                           ###
#########################################################

runs <- 500
all_results <- NULL

for (i in 1:runs){
  dataset <- read.csv(paste("C:/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github/N.trial = 5, ni.trial = 20, gamma = 0.1/Dataset ",i,".txt",sep=""), sep="")
  dataset_long <- gather(dataset, endpoint, outcome, Surr:True)
  
  fit.lme.1 <- try(lme(outcome~ -1 + as.factor(endpoint) + as.factor(endpoint):Treat, 
                       random=~ -1 + as.factor(endpoint) + as.factor(endpoint):Treat|as.factor(Trial.ID),
                       data=dataset_long,
                       correlation=nlme::corSymm(form=~1|as.factor(Trial.ID)/as.factor(Pat.ID)),
                       weights=nlme::varIdent(form=~1|endpoint)), 
                   silent=TRUE)
  
  if(class(fit.lme.1)=="try-error"){
    result.lme.1 <- data.frame(NA, NA, NA, NA)
    names(result.lme.1) <- c("R2.trial.lme.1", "Cond.num.lme.1", "Min.eigen.lme.1", "PD.lme.1")
  }
  if (class(fit.lme.1)!="try-error"){
    D.lme.1 <- fit.lme.1$modelStruct$reStruct$`as.factor(Trial.ID)`
    A.lme.1 <- matrix(c(D.lme.1[4,1], D.lme.1[4,3]), nrow=2, ncol=1)
    B.lme.1 <- matrix(c(D.lme.1[1,1], D.lme.1[3,1], D.lme.1[3,1], D.lme.1[3,3]), nrow=2, ncol=2)
    C.lme.1 <- as.matrix(c(D.lme.1[4,4]))
    R2.trial.lme.1 <- (t(A.lme.1) %*% solve(B.lme.1) %*% A.lme.1)/C.lme.1
    Min.eigen.D.lme.1 <- min(eigen(D.lme.1)$values)
    if (Min.eigen.D.lme.1 <= 0){
      PD.lme.1 <- "No"}
    else{
      PD.lme.1 <- "Yes"}
    singular.lme.1  <- svd(D.lme.1)$d
    Cond.num.lme.1 <- max(singular.lme.1)/min(singular.lme.1)
    result.lme.1 <- data.frame(R2.trial.lme.1, Cond.num.lme.1, Min.eigen.D.lme.1, PD.lme.1)
    names(result.lme.1) <- c("R2.trial.lme.1", "Cond.num.lme.1", "Min.eigen.lme.1", "PD.lme.1")
  }
  result <- data.frame(i, result.lme.1)
  all_results <- rbind(all_results, result)
}
names(all_results)[1] <- c("Dataset")
write.xlsx(all_results, file="Sim lme1 setting1_5-0.1.xlsx")

############################################################
### Analysis with lme                                    ###
### (without description of the within-group correlation ### 
### and heteroscedasticity)                              ###
############################################################

runs <- 500
all_results <- NULL

for (i in 1:runs){
  dataset <- read.csv(paste("C:/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github/N.trial = 5, ni.trial = 20, gamma = 0.1/Dataset ",i,".txt",sep=""), sep="")
  dataset_long <- gather(dataset, endpoint, outcome, Surr:True)
  
  fit.lme.2 <- try(lme(outcome~ -1 + as.factor(endpoint) + as.factor(endpoint):Treat, 
                       random=~ -1 + as.factor(endpoint) + as.factor(endpoint):Treat|as.factor(Trial.ID),
                       data=dataset_long), 
                   silent=TRUE)
  
  if(class(fit.lme.2)=="try-error"){
    result.lme.2 <- data.frame(NA, NA, NA, NA)
    names(result.lme.2) <- c("R2.trial.lme.2", "Cond.num.lme.2", "Min.eigen.lme.2", "PD.lme.2")
  }
  if (class(fit.lme.2)!="try-error"){
    D.lme.2 <- fit.lme.2$modelStruct$reStruct$`as.factor(Trial.ID)`
    A.lme.2 <- matrix(c(D.lme.2[4,1], D.lme.2[4,3]), nrow=2, ncol=1)
    B.lme.2 <- matrix(c(D.lme.2[1,1], D.lme.2[3,1], D.lme.2[3,1], D.lme.2[3,3]), nrow=2, ncol=2)
    C.lme.2 <- as.matrix(c(D.lme.2[4,4]))
    R2.trial.lme.2 <- (t(A.lme.2) %*% solve(B.lme.2) %*% A.lme.2)/C.lme.2
    Min.eigen.D.lme.2 <- min(eigen(D.lme.2)$values)
    if(Min.eigen.D.lme.2 <= 0){
      PD.lme.2 <- "No"}
    else{
      PD.lme.2 <- "Yes"}
    singular.lme.2  <- svd(D.lme.2)$d
    Cond.num.lme.2 <- max(singular.lme.2)/min(singular.lme.2)
    result.lme.2 <- data.frame(R2.trial.lme.2, Cond.num.lme.2, Min.eigen.D.lme.2, PD.lme.2)
    names(result.lme.2) <- c("R2.trial.lme.2", "Cond.num.lme.2", "Min.eigen.lme.2", "PD.lme.2")
  }
  result <- data.frame(i, result.lme.2)
  all_results <- rbind(all_results, result)
}
names(all_results)[1] <- c("Dataset")
write.xlsx(all_results, file="Sim lme2 setting1_5-0.1.xlsx")

##########################
### Analysis with lmer ###
##########################

runs <- 500
all_results <- NULL

for (i in 1:runs){
  dataset <- read.csv(paste("C:/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github/N.trial = 5, ni.trial = 20, gamma = 0.1/Dataset ",i,".txt",sep=""), sep="")
  dataset_long <- gather(dataset, endpoint, outcome, Surr:True)
  
  fit.lmer <- try(lmer(outcome ~ -1 + as.factor(endpoint) + as.factor(endpoint):Treat + 
                         (-1 + as.factor(endpoint) + as.factor(endpoint):Treat|Trial.ID), 
                       data=dataset_long), 
                  silent=TRUE)
  
  if (is.null(summary(fit.lmer)$optinfo$conv$lme4$code)==FALSE){  
    result.lmer <- data.frame(NA, NA, NA, NA)
    names(result.lmer) <- c("R2.trial.lmer", "Cond.num.lmer", "Min.eigen.lmer", "PD.lmer")
  } 
  if (is.null(summary(fit.lmer)$optinfo$conv$lme4$code)==TRUE){  
    D.lmer <- matrix(summary(fit.lmer)$varcor$Trial.ID[1:16], ncol=4)
    A.lmer <- matrix(c(D.lmer[4,1], D.lmer[4,3]), nrow=2, ncol=1)
    B.lmer <- matrix(c(D.lmer[1,1], D.lmer[3,1], D.lmer[3,1], D.lmer[3,3]), nrow=2, ncol=2)
    C.lmer <- as.matrix(D.lmer[4,4])
    R2.trial.lmer <- try((t(A.lmer) %*% solve(B.lmer) %*% A.lmer)/C.lmer, silent=TRUE)
    if(class(R2.trial.lmer)=="try-error"){
      R2.trial.lmer <- NA}
    Min.eigen.D.lmer <- min(eigen(D.lmer)$values)
    if (Min.eigen.D.lmer <= 0){
      PD.lmer <- "No"}
    else{
      PD.lmer <- "Yes"}
    singular.lmer <- svd(D.lmer)$d
    Cond.num.lmer <- max(singular.lmer)/min(singular.lmer)
    result.lmer <- data.frame(R2.trial.lmer, Cond.num.lmer, Min.eigen.D.lmer, PD.lmer)
    names(result.lmer) <- c("R2.trial.lmer", "Cond.num.lmer", "Min.eigen.lmer", "PD.lmer")
  }
  result <- data.frame(i, result.lmer)
  all_results <- rbind(all_results, result)
}
names(all_results)[1] <- c("Dataset")
write.xlsx(all_results, file="Sim lmer setting1_5-0.1.xlsx")

#####################################
### Analysis with BifixedContCont ###
### Two-stage approach            ###
#####################################

runs <- 500
all_results <- NULL

for (i in 1:runs){
  dataset <- read.csv(paste("C:/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github/N.trial = 5, ni.trial = 20, gamma = 0.1/Dataset ",i,".txt",sep=""), sep="")
  
  fit.bifix <- try(BifixedContCont(Dataset = dataset, 
                                   Surr = Surr, 
                                   True = True, 
                                   Treat = Treat, 
                                   Trial.ID = Trial.ID, 
                                   Pat.ID = Pat.ID, 
                                   Model = "Full", 
                                   Weighted = TRUE), 
                   silent=TRUE)
  
  if(class(fit.bifix)=="try-error"){
    result.bifix <- data.frame(NA, NA, NA, NA)
    names(result.bifix) <- c("R2.trial.bifix", "Cond.num.bifix", "Min.eigen.bifix", "PD.bifix")
  }
  if (class(fit.bifix)!="try-error"){
    R2.trial.bifix <- fit.bifix$Trial.R2[1]
    D.bifix <- fit.bifix$D.Equiv
    Min.eigen.D.bifix <- min(eigen(D.bifix)$values)
    if(Min.eigen.D.bifix <= 0){
      PD.bifix <- "No"}
    else{
      PD.bifix <- "Yes"}
    singular.bifix  <- svd(D.bifix)$d
    Cond.num.bifix <- max(singular.bifix)/min(singular.bifix)
    result.bifix <- data.frame(R2.trial.bifix, Cond.num.bifix, Min.eigen.D.bifix, PD.bifix)
    names(result.bifix) <- c("R2.trial.bifix", "Cond.num.bifix", "Min.eigen.bifix", "PD.bifix")
  }
  result <- data.frame(i, result.bifix)
  all_results <- rbind(all_results, result)
}
names(all_results)[1] <- c("Dataset")
write.xlsx(all_results, file="Sim bifix setting1_5-0.1.xlsx")

###########################################
### Analysis with PROC MIXED (from SAS) ###
###########################################

############
#### UN ####
############

runs <- 500
all_results <- NULL

for (i in 1:runs){
  data <- read.csv(paste("C:/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github/N.trial = 5, ni.trial = 20, gamma = 0.1/Output ",i," SAS.txt", sep=""), sep="")
  
  SAS_converged <- data$SAS_converged
  SAS_pd_G <- data$SAS_posdef_G
  SAS_pd_H <- data$SAS_posdef_H
  
  D.sas.un <- matrix(c(data$D_11,data$D_21,data$D_31,data$D_41,
                       data$D_21,data$D_22,data$D_32,data$D_42,
                       data$D_31,data$D_32,data$D_33,data$D_43,
                       data$D_41,data$D_42,data$D_43,data$D_44), 
                     nrow = 4, ncol = 4)
  
  A.sas.un <- matrix(c(D.sas.un[4,1], D.sas.un[4,3]), nrow=2, ncol=1)
  B.sas.un <- matrix(c(D.sas.un[1,1], D.sas.un[3,1], D.sas.un[3,1], D.sas.un[3,3]), nrow=2, ncol=2)
  C.sas.un <- as.matrix(c(D.sas.un[4,4]))
  R2.trial.sas.un <- (t(A.sas.un) %*% solve(B.sas.un) %*% A.sas.un)/C.sas.un
  
  Min.eigen.D.sas.un <- min(eigen(D.sas.un)$values)
  if(Min.eigen.D.sas.un <= 0){
    PD.sas.un <- "No"}
  else{
    PD.sas.un <- "Yes"}
  singular.sas.un  <- svd(D.sas.un)$d
  Cond.num.sas.un <- max(singular.sas.un)/min(singular.sas.un)
  result.sas.un <- data.frame(R2.trial.sas.un, Cond.num.sas.un, Min.eigen.D.sas.un, PD.sas.un)
  names(result.sas.un) <- c("R2.trial.sas.un", "Cond.num.sas.un", "Min.eigen.sas.un", "PD.sas.un")
  
  result <- data.frame(i, SAS_converged, SAS_pd_G, SAS_pd_H, result.sas.un)
  all_results <- rbind(all_results, result)
}
names(all_results)[1] <- c("Dataset")
write.xlsx(all_results, file="Sim sas-un setting1_5-0.1.xlsx") 

