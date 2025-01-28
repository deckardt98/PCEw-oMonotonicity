### This file contains simulation code for the article entitled "Quadruply robust 
### principal stratification analysis without monotonicity" by Tong et al..
### SACE

library(truncnorm)
library(SuperLearner)
library(randomForest)
library(geex)
library(parallel)
library(data.table)
library(boot)
library(caret)

setwd("~/Desktop/Fan Li Project/Hayden_Biometrics/Simulation/Code")
source("EST_FUN11.R")

#Function simulating full data vector
simu_full_data <- function(n, theta){
  #Input:
  #n: sample size
  #theta: odds ratio
  
  Z <- D <- c()
  #simulate covariates
  X1 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X2 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X3 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X4 <- rbinom(n, size = 1, prob = 0.5)
  if (theta>9999){
    X1 <- abs(X1)
    X2 <- abs(X2)
    X3 <- abs(X3)
    probZ <- expit(0.1*(X1+X2+X3)+0.5*X4)
    #max(abs(probZ-0.5))
    Z <- rbinom(n, size = 1, prob = probZ)
    #mean(Z)
    
    #simulate G=(D(0),D(1))
    probD1 <- expit(1.2*X4)
    probD0 <- expit(-0.4-0.2*X1-0.2*X2-0.2*X3-0.2*X4)
    prob11 <- probD0
    prob01 <- probD1 - probD0
    prob00 <- 1 - probD1
    prob_matrix <- cbind(prob00, prob01, prob11)
    # Simulate G from the categorical distribution
    G <- apply(prob_matrix, 1, function(p) sample(0:2, size = 1, prob = p))
    #simulate D(1)
    D1 <- as.numeric(I(G!=0))
    #simulate D(0)
    D0 <- as.numeric(I(G==2))
    #compute observed intermediate outcome
    D <- D1*Z+(1-Z)*D0
  } else {
  probZ <- expit(0.1*(X1+X2+X3)+0.5*X4)
  #max(abs(probZ-0.5))
  Z <- rbinom(n, size = 1, prob = probZ)
  #mean(Z)
  #simulate (D(0),D(1))
  probD1 <- expit(0.3*X1+0.4*X2+0.3*X3+0.5*X4)
  probD0 <- expit(0.4*X1+0.3*X2+0.4*X3+0.5*X4)
  deltaX <- (1+(theta-1)*(probD1+probD0))^2-4*theta*(theta-1)*probD1*probD0
  prob11 <- (1+(theta-1)*(probD1+probD0)-sqrt(deltaX))/2/(theta-1)
  prob10 <- probD0 - prob11
  prob01 <- probD1 - prob11
  prob00 <- 1-prob11-prob10-prob01
  prob_matrix <- cbind(prob10, prob00, prob01, prob11)
  jointD <- apply(prob_matrix, 1, function(p) sample(1:4, size = 1, prob = p))
  #simulate D(0)
  D0 <- as.numeric(I(jointD==1|jointD==4))
  #simulate D(1)
  D1 <- as.numeric(I(jointD==3|jointD==4))
  #compute observed intermediate outcome
  D <- D1*Z+(1-Z)*D0
  #max(max(abs(0.5-probD1)),max(abs(0.5-probD0)))
  #table(D1,D0)/n
  }
  #simulate potential final outcome
  meanY1 <- -1+D1+X1+3*X2+3*X3+3*X4
  meanY0 <- 3-D0-1.5*X1+2*X2+2*X3-2*X4
  Y1 <- rnorm(n = n, mean = meanY1, sd = 1)
  Y0 <- rnorm(n = n, mean = meanY0, sd = 1)
  Y <- Y1*Z+(1-Z)*Y0
  
  return(list(X1, X2, X3, X4, Z, D1, D0, D, Y1, Y0, Y))
}


#Approximate truth
cal_tr <-  function(n, theta){
  full_data <-  simu_full_data(n, theta=theta)
  X1 <- full_data[[1]]
  X2 <- full_data[[2]]
  X3 <- full_data[[3]]
  X4 <- full_data[[4]]
  Z <- full_data[[5]]
  D1 <- full_data[[6]]
  D0 <- full_data[[7]]
  D <- full_data[[8]]
  Y1 <- full_data[[9]]
  Y0 <- full_data[[10]]
  Y <- full_data[[11]]
  #approximate truth
  truepce <- mean(D0*D1*(Y1-Y0))/mean(D0*D1)
  return(truepce)
}


sim_bias <-  function(n, B, truepce, theta){
  bias <- matrix(NA, nrow = B, ncol = 41)
  aese <- matrix(NA, nrow = B, ncol = 41)
  cov <- matrix(NA, nrow = B, ncol = 41)
  for (j in 1:B){
    tryCatch({
      full_data <-  simu_full_data(n, theta)
      X1 <- full_data[[1]]
      X2 <- full_data[[2]]
      X3 <- full_data[[3]]
      X4 <- full_data[[4]]
      Z <- full_data[[5]]
      D <- full_data[[8]]
      Y <- full_data[[11]]
      df <- as.data.frame(cbind(X1,X2,X3,X4,Z,D,Y))
      start <- Sys.time()
      result_est <- cal_estimates(df, theta, truepce)
      print( Sys.time() - start )
      bias[j,] <- result_est[seq(1,82, by =2)]
      aese[j,] <- result_est[-seq(1,82, by =2)]
      #print(round(aese[j,],2))
      cov[j,] <- I(truepce>=(bias[j,]-aese[j,]*1.96))*I(truepce<=(bias[j,]+aese[j,]*1.96))
    }, 
    warning = function(w){
      message("Warning in iteration ", j, ": ", conditionMessage(w))
      return(NULL)
    },
    error = function(e) {
      message("Error in iteration ", j, ": ", conditionMessage(e))
      return(NULL)})
  }
  
  return(list(bias,aese,cov))
}


set.seed(920784642)
theta <- 999999
truepce <- cal_tr(10000000, theta)
n <- 500
B <- 1000
result <- sim_bias(n,B,truepce, theta)
copy_result <- result
mcsd <- round(apply(result[[1]], 2, sd, na.rm = TRUE),3)
mcsd
mbias <- round(colMeans(result[[1]], na.rm = TRUE)-truepce,3)
mbias
aese <- round(colMeans(result[[2]], na.rm = TRUE),3)
aese
cp <- round(colMeans(result[[3]], na.rm = TRUE),3)
cp

out_table <- as.data.frame(cbind(mbias,mcsd,aese,cp))
rownames(out_table) <- c("Weighting",
                         "PAR0.5_1", "PAR1_1", "PAR2_1", "PARinfty_1", 
                         "ML0.5_1", "ML1_1", "ML2_1", "MLinfty_1", 
                         "PAR0.5_2", "PAR1_2", "PAR2_2", "PARinfty_2", 
                         "ML0.5_2", "ML1_2", "ML2_2", "MLinfty_2", 
                         "PAR0.5_3", "PAR1_3", "PAR2_3", "PARinfty_3", 
                         "ML0.5_3", "ML1_3", "ML2_3", "MLinfty_3", 
                         "PAR0.5_4", "PAR1_4", "PAR2_4", "PARinfty_4", 
                         "ML0.5_4", "ML1_4", "ML2_4", "MLinfty_4",
                         "PAR0.5_5", "PAR1_5", "PAR2_5", "PARinfty_5", 
                         "ML0.5_5", "ML1_5", "ML2_5", "MLinfty_5")
out_table
setwd("/Users/deckard/Desktop/Fan Li Project/Hayden_Biometrics/Simulation/results")
write.csv(out_table, file= "11table_new_MONO.csv")

