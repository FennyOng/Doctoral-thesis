
#######################################################
# Code for replication of analysis in Chapter 4 and 5 #
#######################################################

######################
## SIMULATION STUDY ##
######################

setwd("/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github")

library(Surrogate)     # load the Surrogate library
library(tidyverse)     # load the tidyverse library
library(Hmisc)         # load the Hmisc library (function: rMultinom)
library(MASS)          # load the MASS library (function: mvrnorm)
library("openxlsx")    # load the openxlsx library (function: write.xlsx)

############################################################
### Data generation                                      ###
### Code only for simulation when the true ICA is 0.2313 ###
### For other settings, change the values of parameters  ###
### as summarized in Appendix B.2.1 and B.2.2            ###
############################################################

N_total <- 2000
runs <- 300

setting <- paste("Setting 1 R2H = 0.2313", sep="")
setwd("/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github")
dir.create(setting)
setwd(setting)
set.seed(87780)

for(i in 1:runs){
  pi_00 <- 0.2156
  pi_01 <- 0.0906
  pi_10 <- 0.0201
  pi_11 <- 0.6737
  pi_T_all <- matrix(c(pi_00, pi_01, pi_10, pi_11), nrow=1)
  
  mean_S0_00 <- 1.8217
  mean_S0_01 <- 0.7258
  mean_S0_10 <- 2.3944
  mean_S0_11 <- 0.7258
  
  mean_S1_00 <- 1.3282
  mean_S1_01 <- 1.3282
  mean_S1_10 <- 2.4556
  mean_S1_11 <- 2.3304
  
  sigma_S00_00 <- 0.1184
  sigma_S00_01 <- 0.1184
  sigma_S00_10 <- 0.1184
  sigma_S00_11 <- 0.1184
  
  sigma_S11_00 <- 0.3121
  sigma_S11_01 <- 0.3121
  sigma_S11_10 <- 0.3121
  sigma_S11_11 <- 0.3121
  
  rho_S01_00 <- 0.3987
  rho_S01_01 <- 0.1920
  rho_S01_10 <- 0.5669
  rho_S01_11 <- 0.5169
  sigma_S01 <- c((sqrt(sigma_S00_00*sigma_S11_00)*rho_S01_00), 
                 (sqrt(sigma_S00_01*sigma_S11_01)*rho_S01_01), 
                 (sqrt(sigma_S00_10*sigma_S11_10)*rho_S01_10),
                 (sqrt(sigma_S00_11*sigma_S11_11)*rho_S01_11))
  
  Sigma_S00 <- matrix(c(sigma_S00_00, sigma_S01[1], sigma_S01[1], sigma_S11_00), ncol=2, byrow=T)
  Sigma_S01 <- matrix(c(sigma_S00_01, sigma_S01[2], sigma_S01[2], sigma_S11_01), ncol=2, byrow=T)
  Sigma_S10 <- matrix(c(sigma_S00_10, sigma_S01[3], sigma_S01[3], sigma_S11_10), ncol=2, byrow=T)
  Sigma_S11 <- matrix(c(sigma_S00_11, sigma_S01[4], sigma_S01[4], sigma_S11_11), ncol=2, byrow=T)
  
  min_eigen_S00 <- try(min(eigen(Sigma_S00)$values), TRUE)
  min_eigen_S01 <- try(min(eigen(Sigma_S01)$values), TRUE)
  min_eigen_S10 <- try(min(eigen(Sigma_S10)$values), TRUE)
  min_eigen_S11 <- try(min(eigen(Sigma_S11)$values), TRUE)
  
  singular_S00 <- try(svd(Sigma_S00)$d, TRUE)
  CN_S00 <- try(max(singular_S00)/min(singular_S00), TRUE)
  singular_S01 <- try(svd(Sigma_S01)$d, TRUE)
  CN_S01 <- try(max(singular_S01)/min(singular_S01), TRUE)
  singular_S10 <- try(svd(Sigma_S10)$d, TRUE)
  CN_S10 <- try(max(singular_S10)/min(singular_S10), TRUE)
  singular_S11 <- try(svd(Sigma_S11)$d, TRUE)
  CN_S11 <- try(max(singular_S11)/min(singular_S11), TRUE)
  
  if(min_eigen_S00 > 0 && min_eigen_S01 > 0 && min_eigen_S10 > 0 && min_eigen_S11 > 0 && 
     CN_S00 < 50 && CN_S01 < 50 && CN_S10 < 50 && CN_S11 < 50){
    
    data <- data.frame(t(rMultinom(pi_T_all, N_total)))
    colnames(data) <- c("Category")
    data$ID <- c(1:dim(data)[1])
    data$T0 <- ifelse(data$Category==1|data$Category==2, 0, 1)
    data$T1 <- ifelse(data$Category==1|data$Category==3, 0, 1)
    
    data$Treat <- sample(x=rep(c(-1,1), each=ceiling(dim(data)[1]/2)), dim(data)[1], replace=F)
    
    data <- data[order(data$Category),]
    
    S_00 <- try(mvrnorm(n=sum(data$Category==1), mu=c(mean_S0_00, mean_S1_00), Sigma=Sigma_S00), silent=T)
    S_01 <- try(mvrnorm(n=sum(data$Category==2), mu=c(mean_S0_01, mean_S1_01), Sigma=Sigma_S01), silent=T)
    S_10 <- try(mvrnorm(n=sum(data$Category==3), mu=c(mean_S0_10, mean_S1_10), Sigma=Sigma_S10), silent=T)
    S_11 <- try(mvrnorm(n=sum(data$Category==4), mu=c(mean_S0_11, mean_S1_11), Sigma=Sigma_S11), silent=T)
    
    data_S <- data.frame(rbind(S_00, S_01, S_10, S_11), stringsAsFactors=T)
    colnames(data_S) <- c("S0", "S1")
    
    data <- cbind(data, data_S)
    data <- data[order(data$ID),]
    
    data_obs <- data.frame(matrix(NA, nrow=dim(data)[1], ncol=4), stringsAsFactors=T)
    colnames(data_obs) <- c("ID", "True", "Surr", "Treat")
    for(m in 1:dim(data)[1]){
      if(data$Treat[m]==-1){
        data_obs$ID[m] <- data$ID[m]
        data_obs$True[m] <- data$T0[m]
        data_obs$Surr[m] <- data$S0[m]
        data_obs$Treat[m] <- data$Treat[m]
      }
      if(data$Treat[m]==1){
        data_obs$ID[m] <- data$ID[m] 
        data_obs$True[m] <- data$T1[m]
        data_obs$Surr[m] <- data$S1[m]
        data_obs$Treat[m] <- data$Treat[m]
      }
    }
  }
  data_complete <- data.frame(cbind(data$ID, data$T0, data$T1, data$S0, data$S1, data$Treat))
  names(data_complete) <- c("ID", "T0", "T1", "S0", "S1", "Treat")
  write.table(data_complete, file=paste("Dataset counterfactual ", i, ".txt", sep=""), row.names=F)
  write.table(data_obs, file=paste("Dataset TSZ ", i, ".txt", sep=""), row.names=F)
}

########################################
### Analysis with ICA.BinCont        ###
### General analysis / no assumption ###
### without bootstrap                ###
########################################

runs <- 300
all_results <- NULL

for (i in 1:runs){
  dataset <- read.csv(paste("C:/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github/Setting 1 R2H = 0.2313/Dataset TSZ ",i,".txt",sep=""), sep="")
  
  fit.ica.bincont <- try(ICA.BinCont(Dataset = dataset, 
                                     Surr = Surr, 
                                     True = True, 
                                     Treat = Treat,
                                     BS = FALSE,
                                     G_pi_10 = c(0,1),
                                     G_rho_01_00=c(-1,1), 
                                     G_rho_01_01=c(-1,1), 
                                     G_rho_01_10=c(-1,1), 
                                     G_rho_01_11=c(-1,1),
                                     Theta.S_0=c(1.8, 0.7, 2.4, 0.7, 0.3, 0.3, 0.3, 0.3), 
                                     Theta.S_1=c(1.3, 1.3, 2.5, 2.3, 0.6, 0.6, 0.6, 0.6),
                                     M = 1000, Seed = 123,
                                     Monotonicity = FALSE,
                                     Independence = FALSE,
                                     HAA = FALSE,
                                     Cond_ind = FALSE,
                                     Plots = FALSE), 
                         silent=T)
  
  if(class(fit.ica.bincont)=="try-error"){
    result.ica <- data.frame(NA, NA, NA, NA, NA, NA, NA, NA)
    names(result.ica) <- c("R2H.mean", "R2H.median", "R2H.min", "R2H.max", "R2H.q10", "R2H.q90", "R2H.q2.5", "R2H.q97.5")
  }
  if (class(fit.ica.bincont)!="try-error"){
    fit.ica.bincont$R2_H <- na.exclude(fit.ica.bincont$R2_H)
    R2H.mean <- mean(fit.ica.bincont$R2_H)
    R2H.median <- median(fit.ica.bincont$R2_H)
    R2H.min <- min(fit.ica.bincont$R2_H)
    R2H.max <- max(fit.ica.bincont$R2_H)
    R2H.q10 <- quantile(fit.ica.bincont$R2_H, probs=c(0.10))
    R2H.q90 <- quantile(fit.ica.bincont$R2_H, probs=c(0.90))
    R2H.q2.5 <- quantile(fit.ica.bincont$R2_H, probs=c(0.025))
    R2H.q97.5 <- quantile(fit.ica.bincont$R2_H, probs=c(0.975))
    
    result.ica <- data.frame(R2H.mean, R2H.median, R2H.min, R2H.max, R2H.q10, R2H.q90, R2H.q2.5, R2H.q97.5)
    names(result.ica) <- c("R2H.mean", "R2H.median", "R2H.min", "R2H.max", "R2H.q10", "R2H.q90", "R2H.q2.5", "R2H.q97.5")
  }
  results <- data.frame(i, result.ica)
  all_results <- rbind(all_results, results)
  save(fit.ica.bincont, file = paste("Sim_setting1_general", i))
}
names(all_results)[1] <- c("Dataset")
write.xlsx(all_results, file="Sim_setting1_general.xlsx")

######################################
### Adding information of coverage ###
######################################

data <- read.xlsx("/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github/Setting 1 R2H = 0.2313/Sim_setting1_general.xlsx")
data <- data %>% 
  mutate(data, coverage = ifelse(0.2313>=data$R2H.min & 0.2313<=data$R2H.max,1,0)) %>%
  mutate(data, coverage0.95 = ifelse(0.2313>=data$R2H.q2.5 & 0.2313<=data$R2H.q97.5,1,0)) %>%
  mutate(data, coverage0.80 = ifelse(0.2313>=data$R2H.q10 & 0.2313<=data$R2H.q90,1,0))
write.xlsx(data, file="Sim_setting1_general.xlsx")

#########################################
### Frequency distribution of the ICA ###
#########################################

runs <- 300
all_result <- NULL
for (i in 1:runs){
  load(file=paste("/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github/Setting 1 R2H = 0.2313/Sim_setting1_general ",i,"", sep=""))
  R2H <- fit.ica.bincont$R2_H
  results <- data.frame(i, R2H)
  names(results)[1] <- c("Dataset")
  all_result <- rbind(all_result, results)
}
write.xlsx(all_result, file="All ICA sim setting1 general.xlsx")

data <- read.xlsx("/Users/lucp10894/Documents/Hasselt_2020/Thesis/Github/Setting 1 R2H = 0.2313/All ICA sim setting1 general.xlsx")

pdf('Histogram ICA small.pdf')
hist(data$R2H, breaks = 50,
     main = "",
     xlim = c(0,1), ylim = c(0,30000),
     xlab = "True ICA = 0.2313")
abline(v=0.2313, col = "red", lty = 2, lwd = 2)
dev.off()


#########################
## CASE STUDY ANALYSIS ##
#########################

############################################################################
### Since the data from the Influenza Vaccine study                      ###
### cannot be made public, the data from the Schizophrenia study         ###
### is used to illustrate the implementation of the proposed methodology ###
############################################################################

########################################
### Analysis with ICA.BinCont.BS     ###
### General analysis / no assumption ###
### with 300 bootstrap samples       ###
########################################

data("Schizo_BinCont")     # load the schizophrenia data set
Schizo_BinCont_recode <- Schizo_BinCont %>%
  mutate(CGI_Bin = replace(CGI_Bin, CGI_Bin == "0", -1)) %>%     # recode the variable CGI_Bin such that T1-T0 = 1 
  mutate(CGI_Bin = replace(CGI_Bin, CGI_Bin == "1", 0)) %>%      # corresponds to a positive result
  mutate(CGI_Bin = replace(CGI_Bin, CGI_Bin == "-1", 1)) %>%
  mutate(PANSS = PANSS * -1)                                     # recode the variable PANSS such that a more positive
                                                                 # value corresponds to a positive result
schizo_no_miss <- data.frame(na.exclude(cbind(Schizo_BinCont_recode$PANSS, Schizo_BinCont_recode$CGI_Bin, 
                                              Schizo_BinCont_recode$Treat)))
colnames(schizo_no_miss) <- c("PANSS", "CGI_Bin", "Treat")

fit.ica.bincont <- ICA.BinCont.BS(Dataset = schizo_no_miss, 
                                  Surr = PANSS, 
                                  True = CGI_Bin, 
                                  Treat = Treat,
                                  BS = TRUE, nb = 300,
                                  G_pi_10=c(0,1),
                                  G_rho_01_00=c(-1,1),
                                  G_rho_01_01=c(-1,1),
                                  G_rho_01_10=c(-1,1),
                                  G_rho_01_11=c(-1,1),
                                  Theta.S_0 = c(50, 15, -10, -30, 4, 4, 4, 4),
                                  Theta.S_1 = c(60, 25, 0, -25, 4, 4, 4, 4),
                                  M=1000, Seed=123,
                                  Monotonicity=FALSE,
                                  Independence=FALSE,
                                  HAA=FALSE,
                                  Cond_ind=FALSE,
                                  Plots=FALSE)                          
save(fit.ica.bincont, file = "ICA schizo")
summary(fit.ica.bincont)

pdf('Histogram ICA schizo.pdf')
hist(fit.ica.bincont$R2_H,
     main = "", breaks = 50,
     xlim = c(0,1), ylim = NULL,
     xlab = "ICA")
dev.off()

#################################
### Analysis with SPF.BinCont ###
#################################

spf.schizo <- SPF.BinCont(fit.ica.bincont, a=-15, b=15)
save(spf.schizo, file = "SPF schizo (-15,15)")

summary_func <- function(x){
  mean <- round(mean(x), 4)
  range <- round(range(x), 4)
  SDI80 <- round(quantile(x,c(0.1,0.9)), 4)
  SDI95 <- round(quantile(x,c(0.025,0.975)), 4)
  return(list(mean=mean, range=range, SDI80=SDI80, SDI95=SDI95))
}

summary_func(spf.schizo$r_min1_min1)     # get the summary of P[Delta T=-1|Delta S in (-inf,-15)]
summary_func(spf.schizo$r_0_min1)        # get the summary of P[Delta T=0|Delta S in (-inf,-15)]
summary_func(spf.schizo$r_1_min1)        # get the summary of P[Delta T=1|Delta S in (-inf,-15)]
summary_func(spf.schizo$r_min1_0)        # get the summary of P[Delta T=-1|Delta S in (-15,15)]
summary_func(spf.schizo$r_0_0)           # get the summary of P[Delta T=0|Delta S in (-15,15)]
summary_func(spf.schizo$r_1_0)           # get the summary of P[Delta T=1|Delta S in (-15,15)]
summary_func(spf.schizo$r_min1_1)        # get the summary of P[Delta T=-1|Delta S in (15,inf)]
summary_func(spf.schizo$r_0_1)           # get the summary of P[Delta T=0|Delta S in (15,inf)]
summary_func(spf.schizo$r_1_1)           # get the summary of P[Delta T=1|Delta S in (15,inf)]

pdf('Histogram SPF schizo (-15,15).pdf', width=6, height=8, paper='special')
par(mfrow=c(3,3))
hist(spf.schizo$r_min1_min1,
     main = "",
     xlim = c(0,1), ylim = NULL,
     xlab=expression(paste("P[", Delta, "T = -1|", Delta, "S ", epsilon, " (-", infinity, ",-15)]")))
hist(spf.schizo$r_min1_0,
     main = "",
     xlim = c(0,1), ylim = NULL,
     xlab=expression(paste("P[", Delta, "T = -1|", Delta, "S ", epsilon, " (-15, 15)]")))
hist(spf.schizo$r_min1_1,
     main = "",
     xlim = c(0,1), ylim = NULL,
     xlab=expression(paste("P[", Delta, "T = -1|", Delta, "S ", epsilon, " (15,", infinity, ")]")))
hist(spf.schizo$r_0_min1,
     main = "",
     xlim = c(0,1), ylim = NULL,
     xlab=expression(paste("P[", Delta, "T = 0|", Delta, "S ", epsilon, " (-", infinity, ",-15)]")))
hist(spf.schizo$r_0_0,
     main = "",
     xlim = c(0,1), ylim = NULL,
     xlab=expression(paste("P[", Delta, "T = 0|", Delta, "S ", epsilon, " (-15, 15)]")))
hist(spf.schizo$r_0_1,
     main = "",
     xlim = c(0,1), ylim = NULL,
     xlab=expression(paste("P[", Delta, "T = 0|", Delta, "S ", epsilon, " (15,", infinity, ")]")))
hist(spf.schizo$r_1_min1,
     main = "",
     xlim = c(0,1), ylim = NULL,
     xlab=expression(paste("P[", Delta, "T = 1|", Delta, "S ", epsilon, " (-", infinity, ",-15)]")))
hist(spf.schizo$r_1_0,
     main = "",
     xlim = c(0,1), ylim = NULL,
     xlab=expression(paste("P[", Delta, "T = 1|", Delta, "S ", epsilon, " (-15, 15)]")))
hist(spf.schizo$r_1_1,
     main = "",
     xlim = c(0,1), ylim = NULL,
     xlab=expression(paste("P[", Delta, "T = 1|", Delta, "S ", epsilon, " (15,", infinity, ")]")))
dev.off()

summary_func(spf.schizo$P_DT_0_DS_0)     # get the summary of P[Delta T=0|Delta S=0]

pdf('Histogram CN schizo.pdf')
hist(spf.schizo$P_DT_0_DS_0,
     breaks = 50,
     main = "",
     xlim = c(0,1), ylim = NULL,
     xlab=expression(paste("P(", Delta, "T = 0|", Delta, "S = 0)")))
dev.off()

summary_func(spf.schizo$P_DT_psi_DS_max)     # maximum P[Delta T = psi(Delta S)]

pdf('Histogram max DT_psi_DS schizo (-15,15).pdf')
hist(spf.schizo$P_DT_psi_DS_max, breaks = 30,
     main = "",
     xlim = c(0,1), ylim = NULL,
     xlab=expression(paste("P(", Delta, "T = ", psi, "(", Delta, "S))")))
dev.off()

pdf('Best prediction function schizo (-15,15).pdf', width=6, height=3, paper='special')
par(mfrow=c(1,3))
barplot(table(spf.schizo$best.pred.min1), 
        main = "", 
        xlab = expression(paste(psi["ab"],"(-", infinity, ", -15)")), 
        ylab = "Frequency")
barplot(table(spf.schizo$best.pred.0), 
        main = "", 
        xlab = expression(paste(psi["ab"],"(-15, 15)")), 
        ylab = "Frequency")
barplot(table(spf.schizo$best.pred.1), 
        main = "", 
        xlab = expression(paste(psi["ab"],"(15, ", infinity,")")), 
        ylab = "Frequency")
dev.off()

