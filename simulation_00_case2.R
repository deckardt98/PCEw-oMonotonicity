###This file contains simulation code for the article entitled "Quadruply robust 
###principal stratification analysis without monotonicity" by Tong et al..

library(truncnorm)
library(SuperLearner)
library(randomForest)
library(geex)
library(parallel)
library(data.table)
library(boot)
library(caret)

#Ensemble of learning algorithms using super learner
SLmethods <- c("SL.glm", "SL.rpart", "SL.nnet")

custom_jacobian <- function(func, x, ...) {
  numDeriv::jacobian(func, x, ...)
}

gene_Fi <- function(df,n_group){
  df$strata <- interaction(df$Z, df$D)
  folds <- createFolds(df$strata, k = n_group, list = TRUE, returnTrain = FALSE)
  df$Group <- NA
  for (i in seq_along(folds)) {
    df$Group[folds[[i]]] <- i
  }
  return(df$G)
}

#expit(x)=exp(x)/(1+exp(x))
expit <- function(x){return(exp(x)/(1+exp(x)))}

#Function simulating full data vector
simu_full_data <- function(n){
  #Input:
  #n:sample size
  
  Z <- D <- c()
  #simulate covariates
  X1 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X2 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X3 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X1ab <- abs(X1)
  X2ab <- abs(X2)
  X3ab <- abs(X3)
  X4 <- rbinom(n, size = 1, prob = 0.5)
  probZ <- expit(0.1*(X1+X2+X3)+0.5*X4)
  #max(abs(probZ-0.5))
  Z <- rbinom(n, size = 1, prob = probZ)
  #mean(Z)
  
  #simulate G=(D(0),D(1))
  probD1 <- expit(1.2*X4)
  probD0 <- expit(-0.4-0.2*X1ab-0.2*X2ab-0.2*X3ab-0.2*X4)
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
  #max(max(abs(0.5-probD1)),max(abs(0.5-probD0)))
  #table(D1,D0)/n
  
  #simulate potential final outcome
  meanY1 <- -1+D1+X1+3*X2+3*X3+3*X4
  meanY0 <- 3-D0-1.5*X1+2*X2+2*X3-2*X4
  Y1 <- rnorm(n = n, mean = meanY1, sd = 1)
  Y0 <- rnorm(n = n, mean = meanY0, sd = 1)
  Y <- Y1*Z+(1-Z)*Y0
  
  return(list(X1, X2, X3, X4, Z, D1, D0, D, Y1, Y0, Y))
}

#Function calculating parametric quadruply robust estimator
pqrc_point <- function(Z,D,Y,proscore,prinscore,om){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]
  
  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)
  
  
  #parametric quadruply robust estimator
  pqr1 <- mean((1-prinscore0)*Z*(1-D)*(Y-om10)/proscore+om10*((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))/mean(((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))
  pqr0 <- mean((1-prinscore1)*(1-Z)*(1-D)*(Y-om00)/(1-proscore)+om00*((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))/mean(((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))
  pqrc <- pqr1-pqr0
  denom_est <- mean(((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))
  num_est <- mean((1-prinscore0)*Z*(1-D)*(Y-om10)/proscore+om10*((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))-
    mean((1-prinscore1)*(1-Z)*(1-D)*(Y-om00)/(1-proscore)+om00*((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))
  
  return(c(pqr1,pqr0,pqrc,num_est,denom_est))
}

#Function calculating the variance estimates for the nonparametric machine learning estimator proposed by us
ml_qr_var <- function(Z,D,Y,proscore,prinscore,om){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]
  
  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)
  
  #Estimate the denominator E[p0(X)p1(X)]
  denom_est <- mean(((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))
  
  #Compute variance estimates
  pqr1 <- mean((1-prinscore0)*Z*(1-D)*(Y-om10)/proscore+om10*((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))/mean(((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))
  pqr0 <- mean((1-prinscore1)*(1-Z)*(1-D)*(Y-om00)/(1-proscore)+om00*((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))/mean(((1-prinscore0)*(1-psi_D1)+(-psi_D0+prinscore0)*(1-prinscore1)))
  var_single <- mean(((((1-psi_D0)-(1-prinscore0))*(1-prinscore1)*(om10-pqr1)+(1-prinscore0)*(psi_Y1_D1-pqr1*(1-psi_D1)))-(((1-psi_D1)-(1-prinscore1))*(1-prinscore0)*(om00-pqr0)+(1-prinscore1)*(psi_Y1_D0-pqr0*(1-psi_D0))))^2)/denom_est^2
  
  return(var_single)
}

#Function calculating parametric triply robust estimator in Jiang et al.
ptrc_point <- function(Z,D,Y,proscore,prinscore,om){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]
  
  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)
  
  #parametric quadruply robust estimator
  ptr1 <- mean(psi_Y1_D1)/mean(1-psi_D1)
  ptr0 <- mean((1-prinscore1)/(1-prinscore0)*psi_Y1_D0+om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0)))/mean(1-psi_D1)
  ptrc <- ptr1-ptr0
  denom_est <- mean(1-psi_D1)
  num_est <- mean(psi_Y1_D1)-
    mean((1-prinscore1)/(1-prinscore0)*psi_Y1_D0+om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0)))
  
  return(c(ptr1,ptr0,ptrc,num_est,denom_est))
}

#Function calculating the variance estimates for the machine learning estimator using Jiang et al.'s approach
ml_tr_var <- function(Z,D,Y,proscore,prinscore,om){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]
  
  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)
  
  #Estimate the denominator 
  denom_est <- mean(1-psi_D1)
  
  #Compute variance estimates
  ptr1 <- mean(psi_Y1_D1)/mean(1-psi_D1)
  ptr0 <- mean((1-prinscore1)/(1-prinscore0)*psi_Y1_D0+om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0)))/mean(1-psi_D1)
  var_single <- mean(((psi_Y1_D1-ptr1*(1-psi_D1))-((1-prinscore1)/(1-prinscore0)*psi_Y1_D0+om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0))-ptr0*(1-psi_D1)))^2)/denom_est^2
  
  return(var_single)
}

#Estimating equations for weighing estimator
eq_psw <- function(data, models){
  Z <- data$Z
  D <- data$D
  Y <- data$Y
  
  Xpi <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$pi_model))
  Xp0 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p0_model))
  Xp1 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p1_model))
  
  pi_pos <- 1:ncol(Xpi)
  p0_pos <- (max(pi_pos) + 1):(max(pi_pos) + ncol(Xp0))
  p1_pos <- (max(p0_pos) + 1):(max(p0_pos) + ncol(Xp1))
  
  pi_scores <- grab_psiFUN(models$pi_model, data)
  p0_scores <- grab_psiFUN(models$p0_model, data)
  p1_scores <- grab_psiFUN(models$p1_model, data)
  
  function(theta){
    proscore <- plogis(Xpi %*% theta[pi_pos])
    p0 <- plogis(Xp0 %*% theta[p0_pos])
    p1 <- plogis(Xp1 %*% theta[p1_pos])
    c(pi_scores(theta[pi_pos]),
      p0_scores(theta[p0_pos])*I(Z==0),
      p1_scores(theta[p1_pos])*I(Z==1),
      Z*(1-D)*(1-p0)/proscore*(Y-theta[max(p1_pos)+1]),
      (1-Z)*(1-D)*(1-p1)/(1-proscore)*(Y-theta[max(p1_pos)+2]),
      theta[max(p1_pos)+1]-theta[max(p1_pos)+2]-theta[max(p1_pos)+3]) 
  }
}

#Estimating equations for parametric quadruply robust estimator
eq_qr <- function(data, models){
  Z <- data$Z
  D <- data$D
  Y <- data$Y
  
  Xpi <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$pi_model))
  Xp0 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p0_model))
  Xp1 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p1_model))
  Xm00 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m00_model))
  Xm10 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m10_model))
  
  pi_pos <- 1:ncol(Xpi)
  p0_pos <- (max(pi_pos) + 1):(max(pi_pos) + ncol(Xp0))
  p1_pos <- (max(p0_pos) + 1):(max(p0_pos) + ncol(Xp1))
  m00_pos <- (max(p1_pos) + 1):(max(p1_pos) + ncol(Xm00))
  m10_pos <- (max(m00_pos) + 1):(max(m00_pos) + ncol(Xm10))
  
  pi_scores <- grab_psiFUN(models$pi_model, data)
  p0_scores <- grab_psiFUN(models$p0_model, data)
  p1_scores <- grab_psiFUN(models$p1_model, data)
  m00_scores <- grab_psiFUN(models$m00_model, data)
  m10_scores <- grab_psiFUN(models$m10_model, data)
  
  function(theta){
    proscore <- plogis(Xpi %*% theta[pi_pos])
    p0 <- plogis(Xp0 %*% theta[p0_pos])
    p1 <- plogis(Xp1 %*% theta[p1_pos])
    om00 <- Xm00 %*% theta[m00_pos]
    om10 <- Xm10 %*% theta[m10_pos]
    psi_D0 <- (1-Z)/(1-proscore)*(D-p0)+p0
    psi_D1 <- Z/proscore*(D-p1)+p1
    psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-p1))+om10*(1-p1)
    psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-p0))+om00*(1-p0)
    c(pi_scores(theta[pi_pos]),
      p0_scores(theta[p0_pos])*I(Z==0),
      p1_scores(theta[p1_pos])*I(Z==1),
      m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
      m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
      (p0-psi_D0)*(1-p1)*(om10-theta[max(m10_pos)+1])+(1-p0)*(psi_Y1_D1-theta[max(m10_pos)+1]*(1-psi_D1)),
      (p1-psi_D1)*(1-p0)*(om00-theta[max(m10_pos)+2])+(1-p1)*(psi_Y1_D0-theta[max(m10_pos)+2]*(1-psi_D0)),
      theta[max(m10_pos)+1]-theta[max(m10_pos)+2]-theta[max(m10_pos)+3]) 
  }
}

#Estimating equations for parametric triply robust estimator in Jiang et al.
eq_tr <- function(data, models){
  Z <- data$Z
  D <- data$D
  Y <- data$Y
  
  Xpi <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$pi_model))
  Xp0 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p0_model))
  Xp1 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p1_model))
  Xm00 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m00_model))
  Xm10 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m10_model))
  
  pi_pos <- 1:ncol(Xpi)
  p0_pos <- (max(pi_pos) + 1):(max(pi_pos) + ncol(Xp0))
  p1_pos <- (max(p0_pos) + 1):(max(p0_pos) + ncol(Xp1))
  m00_pos <- (max(p1_pos) + 1):(max(p1_pos) + ncol(Xm00))
  m10_pos <- (max(m00_pos) + 1):(max(m00_pos) + ncol(Xm10))
  
  pi_scores <- grab_psiFUN(models$pi_model, data)
  p0_scores <- grab_psiFUN(models$p0_model, data)
  p1_scores <- grab_psiFUN(models$p1_model, data)
  m00_scores <- grab_psiFUN(models$m00_model, data)
  m10_scores <- grab_psiFUN(models$m10_model, data)
  
  function(theta){
    proscore <- plogis(Xpi %*% theta[pi_pos])
    p0 <- plogis(Xp0 %*% theta[p0_pos])
    p1 <- plogis(Xp1 %*% theta[p1_pos])
    om00 <- Xm00 %*% theta[m00_pos]
    om10 <- Xm10 %*% theta[m10_pos]
    psi_D0 <- (1-Z)/(1-proscore)*(D-p0)+p0
    psi_D1 <- Z/proscore*(D-p1)+p1
    psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-p1))+om10*(1-p1)
    psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-p0))+om00*(1-p0)
    c(pi_scores(theta[pi_pos]),
      p0_scores(theta[p0_pos])*I(Z==0),
      p1_scores(theta[p1_pos])*I(Z==1),
      m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
      m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
      psi_Y1_D1-theta[max(m10_pos)+1]*(1-psi_D1),
      (1-p1)/(1-p0)*psi_Y1_D0+om00*(1-psi_D1-(1-p1)/(1-p0)*(1-psi_D0))-theta[max(m10_pos)+2]*(1-psi_D1),
      theta[max(m10_pos)+1]-theta[max(m10_pos)+2]-theta[max(m10_pos)+3]) 
  }
}


cal_estimates <- function(df, truepce){
  
  n <- nrow(df)
  Z <- df$Z
  D <- df$D
  Y <- df$Y
  X1 <- df$X1
  X2 <- df$X2
  X3 <- df$X3
  X4 <- df$X4
  df$mis_X1 <- log(atan(X1+X4)+1.570796)
  df$mis_X2 <- (X2/5)^3
  df$mis_X3 <- atan(X3+X4)
  df$mis_X4 <- X4
  df$mis_X1ab <- log(atan(abs(X1)+abs(X4))+1.570796)
  df$mis_X2ab <- (abs(X2)/5)^3
  df$mis_X3ab <- atan(abs(X3)+abs(X4))
  df$mis_X4ab <- X4
  
  # Combine the data into a single dataframe
  data <- data.frame(Z = Z, D = D, Y = Y, X1 = X1, X2 = X2, X3 = X3, X4 = X4)
  
  ################Estimate nuisance functions by parametric approaches#############
  #Estimate propensity score
  propen_m <- glm(Z ~ X1+X2+X3+X4, data = df, family = binomial)
  proscore <- predict(propen_m, newdata=df, type="response")
  
  #Estimate propensity score with incorrect model specifications
  mis_propen_m <- glm(Z ~ mis_X1+mis_X2+mis_X3+mis_X4, data = df, family = binomial)
  mis_proscore <- predict(mis_propen_m, newdata=df, type="response")
  
  #Estimate principal score
  ps_m1 <- glm(D ~ abs(X1)+abs(X2)+abs(X3)+X4, data = df, subset = (Z==1), family = binomial)
  ps_m0 <- glm(D ~ abs(X1)+abs(X2)+abs(X3)+X4, data = df, subset = (Z==0), family = binomial)
  prinscore1 <-  predict(ps_m1, newdata=df, type="response")
  prinscore0 <-  predict(ps_m0, newdata=df, type="response")
  prinscore <- list(prinscore1,prinscore0)
  
  #Estimate principal score with incorrect model specifications
  mis_ps_m1 <- glm(D ~ mis_X1ab+mis_X2ab+mis_X3ab+mis_X4ab, data = df, subset = (Z==1), family = binomial)
  mis_ps_m0 <- glm(D ~ mis_X1ab+mis_X2ab+mis_X3ab+mis_X4ab, data = df, subset = (Z==0), family = binomial)
  mis_prinscore1 <-  predict(mis_ps_m1, newdata=df, type="response")
  mis_prinscore0 <-  predict(mis_ps_m0, newdata=df, type="response")    
  mis_prinscore <- list(mis_prinscore1,mis_prinscore0)
  
  #Estimate outcome mean
  om_m11 <- lm(Y ~ X1+X2+X3+X4, data = df, subset = (Z==1)&(D==1))
  om_m01 <- lm(Y ~ X1+X2+X3+X4, data = df, subset = (Z==0)&(D==1))
  om_m10 <- lm(Y ~ X1+X2+X3+X4, data = df, subset = (Z==1)&(D==0))
  om_m00 <- lm(Y ~ X1+X2+X3+X4, data = df, subset = (Z==0)&(D==0))
  om11 <-  predict(om_m11, newdata=df, type="response")
  om01 <-  predict(om_m01, newdata=df, type="response") 
  om10 <-  predict(om_m10, newdata=df, type="response")
  om00 <-  predict(om_m00, newdata=df, type="response") 
  om <- list(om11,om01,om10,om00)
  
  #Estimate outcome mean with incorrect model specifications
  mis_om_m11 <- lm(Y ~ mis_X1+mis_X2+mis_X3+mis_X4, data = df, subset = (Z==1)&(D==1))
  mis_om_m01 <- lm(Y ~ mis_X1+mis_X2+mis_X3+mis_X4, data = df, subset = (Z==0)&(D==1))
  mis_om_m10 <- lm(Y ~ mis_X1+mis_X2+mis_X3+mis_X4, data = df, subset = (Z==1)&(D==0))
  mis_om_m00 <- lm(Y ~ mis_X1+mis_X2+mis_X3+mis_X4, data = df, subset = (Z==0)&(D==0))
  mis_om11 <-  predict(mis_om_m11, newdata=df, type="response")
  mis_om01 <-  predict(mis_om_m01, newdata=df, type="response") 
  mis_om10 <-  predict(mis_om_m10, newdata=df, type="response")
  mis_om00 <-  predict(mis_om_m00, newdata=df, type="response") 
  mis_om <- list(mis_om11,mis_om01,mis_om10,mis_om00)
  
  ##################################################Weighting###############################################################
  ####Point estimate
  est_weighting1 <- mean(Z*(1-D)*(1-prinscore0)/proscore*Y)/mean(Z*(1-D)*(1-prinscore0)/proscore)
  est_weighting0 <- mean((1-Z)*(1-D)*(1-prinscore1)/(1-proscore)*Y)/mean((1-Z)*(1-D)*(1-prinscore1)/(1-proscore))
  est_weightingc <- est_weighting1 - est_weighting0
  
  ####Variance estimate
  theta_ps <- c(as.vector(coef(propen_m)),
                as.vector(coef(ps_m0)),
                as.vector(coef(ps_m1)),
                est_weighting1,
                est_weighting0,
                est_weightingc)
  estimate_psw <- function(data){
    pi_model <- glm(Z~X1+X2+X3+X4, data=data, family = binomial)
    p0_model <- glm(D~abs(X1)+abs(X2)+abs(X3)+X4, subset = (Z==0), data=data, family = binomial)
    p1_model <- glm(D~abs(X1)+abs(X2)+abs(X3)+X4, subset = (Z==1), data=data, family = binomial) 
    models <- list(pi_model=pi_model, p0_model=p0_model, p1_model=p1_model)
    m_estimate(estFUN = eq_psw, data = data,
               outer_args = list(models = models),
               #root_control = setup_root_control(start = theta_ps+rnorm(length(theta_ps),mean=0,sd=0.5))
               compute_roots = FALSE,
               roots = theta_ps,
               deriv_control = setup_deriv_control(
                 FUN = custom_jacobian, 
                 method = "simple"
               )
    )
  }
  cov_matrix_psw <- vcov(estimate_psw(data=df))
  se_weighting <- sqrt(cov_matrix_psw[ncol(cov_matrix_psw),ncol(cov_matrix_psw)])
  
  
  
  
  
  ###################################################Parametric quadruply robust estimator##################################################  
  ####Point estimates
  #All models are correctly specified
  pqr1 <- pqrc_point(Z,D,Y,proscore,prinscore,om)
  pqr11 <- pqr1[1]
  pqr01 <- pqr1[2]
  pqrc1 <- pqr1[3]
  #Misspecification scenario 1
  pqr2 <- pqrc_point(Z,D,Y,proscore,prinscore,mis_om)
  pqr12 <- pqr2[1]
  pqr02 <- pqr2[2]
  pqrc2 <- pqr2[3]
  #Misspecification scenario 2
  pqr3 <- pqrc_point(Z,D,Y,mis_proscore,prinscore,om)
  pqr13 <- pqr3[1]
  pqr03 <- pqr3[2]
  pqrc3 <- pqr3[3]
  #Misspecification scenario 3
  pqr4 <- pqrc_point(Z,D,Y,proscore,list(mis_prinscore1,prinscore0),list(mis_om11,mis_om01,mis_om10,om00))
  pqr14 <- pqr4[1]
  pqr04 <- pqr4[2]
  pqrc4 <- pqr4[3]
  #Misspecification scenario 4
  pqr5 <- pqrc_point(Z,D,Y,proscore,list(prinscore1,mis_prinscore0),list(mis_om11,mis_om01,om10,mis_om00))
  pqr15 <- pqr5[1]
  pqr05 <- pqr5[2]
  pqrc5 <- pqr5[3]
  #All models are wrongly specified
  pqr6 <- pqrc_point(Z,D,Y,mis_proscore,mis_prinscore,mis_om)
  pqr16 <- pqr6[1]
  pqr06 <- pqr6[2]
  pqrc6 <- pqr6[3]
  
  ####Sandwich variance estimates
  create_models <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m10_formula) {
    list(
      pi_model = glm(pi_formula, data = data, family = binomial),
      p0_model = glm(p0_formula, subset = (Z == 0), data = data, family = binomial),
      p1_model = glm(p1_formula, subset = (Z == 1), data = data, family = binomial),
      m00_model = glm(m00_formula, subset = (Z == 0) & (D == 0), data = data, family = gaussian),
      m10_model = glm(m10_formula, subset = (Z == 1) & (D == 0), data = data, family = gaussian)
    )
  }
  
  estimate_qr <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m10_formula, roots) {
    models <- create_models(data, pi_formula, p0_formula, p1_formula, m00_formula, m10_formula)
    m_estimate(
      estFUN = eq_qr, 
      data = data,
      outer_args = list(models = models),
      #root_control = setup_root_control(start = roots+rnorm(length(roots),mean=0,sd=0.1)),
      compute_roots = FALSE,
      roots = roots,
      deriv_control = setup_deriv_control(
        FUN = custom_jacobian, 
        method = "simple"
      )
    )
  }
  
  #roots of the joint estimating equations
  theta_list <- list(
    theta_qr1 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(om_m00)), as.vector(coef(om_m10)), pqr11, pqr01, pqrc1),
    theta_qr2 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(mis_om_m00)), as.vector(coef(mis_om_m10)), pqr12, pqr02, pqrc2),
    theta_qr3 = c(as.vector(coef(mis_propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(om_m00)), as.vector(coef(om_m10)), pqr13, pqr03, pqrc3),
    theta_qr4 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(mis_ps_m1)), 
                  as.vector(coef(om_m00)), as.vector(coef(mis_om_m10)), pqr14, pqr04, pqrc4),
    theta_qr5 = c(as.vector(coef(propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(mis_om_m00)), as.vector(coef(om_m10)), pqr15, pqr05, pqrc5),
    theta_qr6 = c(as.vector(coef(mis_propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(mis_ps_m1)), 
                  as.vector(coef(mis_om_m00)), as.vector(coef(mis_om_m10)), pqr16, pqr06, pqrc6)
  )
  
  formula_list <- list(
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, p1 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, 
         m00 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, p1 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, 
         m00 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4, m10 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4),
    list(pi = Z ~ mis_X1+mis_X2+mis_X3+mis_X4, p0 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, p1 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, 
         m00 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, p1 = D ~ mis_X1ab+mis_X2ab+mis_X3ab+mis_X4ab, 
         m00 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ mis_X1ab+mis_X2ab+mis_X3ab+mis_X4ab, p1 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, 
         m00 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4, m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ mis_X1+mis_X2+mis_X3+mis_X4, p0 = D ~ mis_X1ab+mis_X2ab+mis_X3ab+mis_X4ab, p1 = D ~ mis_X1ab+mis_X2ab+mis_X3ab+mis_X4ab, 
         m00 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4, m10 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4)
  )
  
  #estimates <- lapply(seq_along(theta_list), function(i) {
  #  estimate_qr(df, formula_list[[i]]$pi, formula_list[[i]]$p0, formula_list[[i]]$p1, 
  #              formula_list[[i]]$m00, formula_list[[i]]$m11, theta_list[[i]])
  #})
  
  #  for (k in 1:6){
  #    print(round(lapply(estimates, coef)[[k]]-theta_list[[k]],2))
  #  }
  
  #Compute standard errors
  #cov_matrices <- lapply(estimates, vcov)
  #se_list <- lapply(cov_matrices, function(cov_matrix) {
  #  sqrt(cov_matrix[ncol(cov_matrix), ncol(cov_matrix)])
  #})
  # Compute standard errors directly
  
  se_list <- mclapply(seq_along(theta_list), function(i) {
    # Compute the estimate (but do not keep it)
    estimate <- estimate_qr(df, formula_list[[i]]$pi, formula_list[[i]]$p0, formula_list[[i]]$p1, 
                            formula_list[[i]]$m00, formula_list[[i]]$m10, theta_list[[i]])
    
    # Compute the standard error
    cov_matrix <- vcov(estimate)
    sqrt(cov_matrix[ncol(cov_matrix), ncol(cov_matrix)])
  }, mc.cores = 3)
  
  # Assign standard errors to individual variables
  se_qr1 <- se_list[[1]]
  se_qr2 <- se_list[[2]]
  se_qr3 <- se_list[[3]]
  se_qr4 <- se_list[[4]]
  se_qr5 <- se_list[[5]]
  se_qr6 <- se_list[[6]]
  
  
  ###################################################Parametric triply robust estimator parallel to the quadruply robust estimator regarding model specifications##################################################  
  #### Point estimates
  #All models are correctly specified
  ptr1 <- ptrc_point(Z,D,Y,proscore,prinscore,om)
  ptr11 <- ptr1[1]
  ptr01 <- ptr1[2]
  ptrc1 <- ptr1[3]
  #Misspecification scenario 1
  ptr2 <- ptrc_point(Z,D,Y,proscore,prinscore,mis_om)
  ptr12 <- ptr2[1]
  ptr02 <- ptr2[2]
  ptrc2 <- ptr2[3]
  #Misspecification scenario 2
  ptr3 <- ptrc_point(Z,D,Y,mis_proscore,prinscore,om)
  ptr13 <- ptr3[1]
  ptr03 <- ptr3[2]
  ptrc3 <- ptr3[3]
  #Misspecification scenario 3
  ptr4 <- ptrc_point(Z,D,Y,proscore,list(mis_prinscore1,prinscore0),list(mis_om11,mis_om01,mis_om10,om00))
  ptr14 <- ptr4[1]
  ptr04 <- ptr4[2]
  ptrc4 <- ptr4[3]
  #Misspecification scenario 4
  ptr5 <- ptrc_point(Z,D,Y,proscore,list(prinscore1,mis_prinscore0),list(mis_om11,mis_om01,om10,mis_om00))
  ptr15 <- ptr5[1]
  ptr05 <- ptr5[2]
  ptrc5 <- ptr5[3]
  #All models are wrongly specified
  ptr6 <- ptrc_point(Z,D,Y,mis_proscore,mis_prinscore,mis_om)
  ptr16 <- ptr6[1]
  ptr06 <- ptr6[2]
  ptrc6 <- ptr6[3]
  
  ### Sandwich variance estimator
  estimate_tr <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m10_formula, roots) {
    models <- create_models(data, pi_formula, p0_formula, p1_formula, m00_formula, m10_formula)
    m_estimate(
      estFUN = eq_tr, 
      data = data,
      outer_args = list(models = models),
      #root_control = setup_root_control(start = roots+rnorm(length(roots),mean=0,sd=0.1)),
      compute_roots = FALSE,
      roots = roots,
      deriv_control = setup_deriv_control(
        FUN = custom_jacobian, 
        method = "simple"
      )
    )
  }
  
  #roots of the joint estimating equations
  theta_list_tr <- list(
    theta_tr1 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(om_m00)), as.vector(coef(om_m10)), ptr11, ptr01, ptrc1),
    theta_tr2 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(mis_om_m00)), as.vector(coef(mis_om_m10)), ptr12, ptr02, ptrc2),
    theta_tr3 = c(as.vector(coef(mis_propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(om_m00)), as.vector(coef(om_m10)), ptr13, ptr03, ptrc3),
    theta_tr4 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(mis_ps_m1)), 
                  as.vector(coef(om_m00)), as.vector(coef(mis_om_m10)), ptr14, ptr04, ptrc4),
    theta_tr5 = c(as.vector(coef(propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(mis_om_m00)), as.vector(coef(om_m10)), ptr15, ptr05, ptrc5),
    theta_tr6 = c(as.vector(coef(mis_propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(mis_ps_m1)), 
                  as.vector(coef(mis_om_m00)), as.vector(coef(mis_om_m10)), ptr16, ptr06, ptrc6)
  )
  
  formula_list_tr <- list(
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, p1 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, 
         m00 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, p1 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, 
         m00 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4, m10 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4),
    list(pi = Z ~ mis_X1+mis_X2+mis_X3+mis_X4, p0 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, p1 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, 
         m00 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, p1 = D ~ mis_X1ab+mis_X2ab+mis_X3ab+mis_X4ab, 
         m00 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ mis_X1ab+mis_X2ab+mis_X3ab+mis_X4ab, p1 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, 
         m00 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4, m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ mis_X1+mis_X2+mis_X3+mis_X4, p0 = D ~ mis_X1ab+mis_X2ab+mis_X3ab+mis_X4ab, p1 = D ~ mis_X1ab+mis_X2ab+mis_X3ab+mis_X4ab, 
         m00 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4, m10 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4)
  )
  
  
  #estimates_tr <- lapply(seq_along(theta_list_tr), function(i) {
  #  estimate_tr(df, formula_list_tr[[i]]$pi, formula_list_tr[[i]]$p0, formula_list_tr[[i]]$p1, 
  #              formula_list_tr[[i]]$m00, formula_list_tr[[i]]$m11, theta_list_tr[[i]])
  #})
  #    for (k in 1:6){
  #      print(round(lapply(estimates_tr, coef)[[k]]-theta_list_tr[[k]],8))
  #    }
  
  # compute standard errors
  #cov_matrices_tr <- lapply(estimates_tr, vcov)
  #se_list_tr <- lapply(cov_matrices_tr, function(cov_matrix) {
  #  sqrt(cov_matrix[ncol(cov_matrix), ncol(cov_matrix)])
  #})
  
  se_list_tr <- mclapply(seq_along(theta_list_tr), function(i) {
    # Compute the estimate (but do not keep it)
    estimate <- estimate_tr(df, formula_list_tr[[i]]$pi, formula_list_tr[[i]]$p0, formula_list_tr[[i]]$p1, 
                            formula_list_tr[[i]]$m00, formula_list_tr[[i]]$m10, theta_list_tr[[i]])
    
    # Compute the standard error
    cov_matrix <- vcov(estimate)
    sqrt(cov_matrix[ncol(cov_matrix), ncol(cov_matrix)])
  }, mc.cores = 3)
  
  se_tr1 <- se_list_tr[[1]]
  se_tr2 <- se_list_tr[[2]]
  se_tr3 <- se_list_tr[[3]]
  se_tr4 <- se_list_tr[[4]]
  se_tr5 <- se_list_tr[[5]]
  se_tr6 <- se_list_tr[[6]]
  
  
  ################################Nonparametric quadruply machine learning estimator, 5-fold###############################
  # Define a function to apply transformations and fit SuperLearner
  fit_superlearner <- function(Y, X_list, indices, family, SLmethods, i, id) {
    lapply(X_list, function(X) {
      SuperLearner(Y[indices], X = X[indices, ], newX = X[which(id == i),], family = family, 
                   SL.library = SLmethods, cvControl = list(V = 5))$SL.predict
    })
  }
  
  # Transform the covariates
  X_list <- list(as.data.frame(cbind(X1,X2,X3,X4)),
                 as.data.frame(cbind(df$mis_X1,df$mis_X2,df$mis_X3,df$mis_X4)))
  X_list_abs <- list(as.data.frame(cbind(X1,X2,X3,X4)),
                     as.data.frame(cbind(df$mis_X1ab,df$mis_X2ab,df$mis_X3ab,df$mis_X4ab)))
  
  # Function to execute for each fold
  process_fold <- function(i, idc) {
    # Define indices for filtering
    idx_F <- which(idc != i)
    idx_Z1 <- which(idc != i & Z == 1)
    idx_Z0 <- which(idc != i & Z == 0)
    idx_Z1D1 <- which(idc != i & Z == 1 & D == 1)
    idx_Z0D1 <- which(idc != i & Z == 0 & D == 1)
    idx_Z1D0 <- which(idc != i & Z == 1 & D == 0)
    idx_Z0D0 <- which(idc != i & Z == 0 & D == 0)
    
    # Estimate nuisance functions using SuperLearner
    ml_proscore <- fit_superlearner(Z, X_list, idx_F, family = "binomial", SLmethods, i, idc)
    ml_ps_m1 <- fit_superlearner(D, X_list_abs, idx_Z1, family = "binomial", SLmethods, i, idc)
    ml_ps_m0 <- fit_superlearner(D, X_list_abs, idx_Z0, family = "binomial", SLmethods, i, idc)
    ml_om_m11 <- fit_superlearner(Y, X_list, idx_Z1D1, family = gaussian(), SLmethods, i, idc)
    ml_om_m01 <- fit_superlearner(Y, X_list, idx_Z0D1, family = gaussian(), SLmethods, i, idc)
    ml_om_m10 <- fit_superlearner(Y, X_list, idx_Z1D0, family = gaussian(), SLmethods, i, idc)
    ml_om_m00 <- fit_superlearner(Y, X_list, idx_Z0D0, family = gaussian(), SLmethods, i, idc)
    
    # MLQR
    mlqck1 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                         list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                         list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]))
    mlqck_num1 <- mlqck1[4]
    mlqck_denom1 <- mlqck1[5]
    
    mlqck2 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                         list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                         list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))
    mlqck_num2 <- mlqck2[4]
    mlqck_denom2 <- mlqck2[5]
    
    mlqck3 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                         list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                         list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]))
    mlqck_num3 <- mlqck3[4]
    mlqck_denom3 <- mlqck3[5]
    
    mlqck4 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                         list(ml_ps_m1[[2]], ml_ps_m0[[1]]), 
                         list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[1]]))
    mlqck_num4 <- mlqck4[4]
    mlqck_denom4 <- mlqck4[5]
    
    mlqck5 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                         list(ml_ps_m1[[1]], ml_ps_m0[[2]]), 
                         list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[1]], ml_om_m00[[2]]))
    mlqck_num5 <- mlqck5[4]
    mlqck_denom5 <- mlqck5[5]
    
    mlqck6 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                         list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                         list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))
    mlqck_num6 <- mlqck6[4]
    mlqck_denom6 <- mlqck6[5]
    
    
    mlqck_var1 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                            list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                            list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]))
    
    mlqck_var2 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                            list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                            list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))
    
    mlqck_var3 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                            list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                            list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]))
    
    mlqck_var4 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                            list(ml_ps_m1[[2]], ml_ps_m0[[1]]), 
                            list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[1]]))
    
    mlqck_var5 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                            list(ml_ps_m1[[1]], ml_ps_m0[[2]]), 
                            list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[1]], ml_om_m00[[2]]))
    
    mlqck_var6 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                            list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                            list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))
    
    
    # MLTR
    mltck1 <- ptrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                         list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                         list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]))
    mltck_num1 <- mltck1[4]
    mltck_denom1 <- mltck1[5]
    
    mltck2 <- ptrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                         list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                         list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))
    mltck_num2 <- mltck2[4]
    mltck_denom2 <- mltck2[5]
    
    mltck3 <- ptrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                         list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                         list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]))
    mltck_num3 <- mltck3[4]
    mltck_denom3 <- mltck3[5]
    
    mltck4 <- ptrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                         list(ml_ps_m1[[2]], ml_ps_m0[[1]]), 
                         list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[1]]))
    mltck_num4 <- mltck4[4]
    mltck_denom4 <- mltck4[5]
    
    mltck5 <- ptrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                         list(ml_ps_m1[[1]], ml_ps_m0[[2]]), 
                         list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[1]], ml_om_m00[[2]]))
    mltck_num5 <- mltck5[4]
    mltck_denom5 <- mltck5[5]
    
    mltck6 <- ptrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                         list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                         list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))
    mltck_num6 <- mltck6[4]
    mltck_denom6 <- mltck6[5]
    
    
    
    mltck_var1 <- ml_tr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                            list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                            list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]))
    
    mltck_var2 <- ml_tr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                            list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                            list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))
    
    mltck_var3 <- ml_tr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                            list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                            list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]))
    
    mltck_var4 <- ml_tr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                            list(ml_ps_m1[[2]], ml_ps_m0[[1]]), 
                            list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[1]]))
    
    mltck_var5 <- ml_tr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                            list(ml_ps_m1[[1]], ml_ps_m0[[2]]), 
                            list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[1]], ml_om_m00[[2]]))
    
    mltck_var6 <- ml_tr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                            list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                            list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))   
    
    # Return results as a vector
    c(mlqck_num1, mlqck_num2, mlqck_num3, mlqck_num4, mlqck_num5, mlqck_num6, 
      mlqck_denom1, mlqck_denom2, mlqck_denom3, mlqck_denom4, mlqck_denom5, mlqck_denom6,
      mlqck_var1, mlqck_var2, mlqck_var3, mlqck_var4, mlqck_var5, mlqck_var6,
      mltck_num1, mltck_num2, mltck_num3, mltck_num4, mltck_num5, mltck_num6,
      mltck_denom1, mltck_denom2, mltck_denom3, mltck_denom4, mltck_denom5, mltck_denom6,
      mltck_var1, mltck_var2, mltck_var3, mltck_var4, mltck_var5, mltck_var6)
  }
  ml_error <- TRUE
  n_fold <- 5
  while (ml_error) {
    idc <- gene_Fi(df,n_fold)
    tryCatch({
      results <- mclapply(1:n_fold, process_fold, mc.cores = 5, idc = idc)
      ml_error <- FALSE  
    }, 
    warning = function(w) {
      idc <- gene_Fi(df,n_fold)  
      ml_error <- TRUE  
    },
    error = function(e) {
      idc <- gene_Fi(df,n_fold)  
      ml_error <- TRUE  
    })
  }
  mlqck_num1 <- sapply(results, `[`, 1)
  mlqck_num2 <- sapply(results, `[`, 2)
  mlqck_num3 <- sapply(results, `[`, 3)
  mlqck_num4 <- sapply(results, `[`, 4)
  mlqck_num5 <- sapply(results, `[`, 5)
  mlqck_num6 <- sapply(results, `[`, 6)
  mlqck_denom1 <- sapply(results, `[`, 7)
  mlqck_denom2 <- sapply(results, `[`, 8)
  mlqck_denom3 <- sapply(results, `[`, 9)
  mlqck_denom4 <- sapply(results, `[`, 10)
  mlqck_denom5 <- sapply(results, `[`, 11)
  mlqck_denom6 <- sapply(results, `[`, 12)
  mlqck_var1 <- sapply(results, `[`, 13)
  mlqck_var2 <- sapply(results, `[`, 14)
  mlqck_var3 <- sapply(results, `[`, 15)
  mlqck_var4 <- sapply(results, `[`, 16)
  mlqck_var5 <- sapply(results, `[`, 17)
  mlqck_var6 <- sapply(results, `[`, 18)
  
  mlqc1 <- (sum(mlqck_num1 * as.numeric(table(idc))) / n)/(sum(mlqck_denom1 * as.numeric(table(idc))) / n)
  mlqc2 <- (sum(mlqck_num2 * as.numeric(table(idc))) / n)/(sum(mlqck_denom2 * as.numeric(table(idc))) / n)
  mlqc3 <- (sum(mlqck_num3 * as.numeric(table(idc))) / n)/(sum(mlqck_denom3 * as.numeric(table(idc))) / n)
  mlqc4 <- (sum(mlqck_num4 * as.numeric(table(idc))) / n)/(sum(mlqck_denom4 * as.numeric(table(idc))) / n)
  mlqc5 <- (sum(mlqck_num5 * as.numeric(table(idc))) / n)/(sum(mlqck_denom5 * as.numeric(table(idc))) / n)
  mlqc6 <- (sum(mlqck_num6 * as.numeric(table(idc))) / n)/(sum(mlqck_denom6 * as.numeric(table(idc))) / n)
  
  mlqr_var1 <- sum(mlqck_var1 * as.numeric(table(idc))) / n^2
  mlqr_var2 <- sum(mlqck_var2 * as.numeric(table(idc))) / n^2
  mlqr_var3 <- sum(mlqck_var3 * as.numeric(table(idc))) / n^2
  mlqr_var4 <- sum(mlqck_var4 * as.numeric(table(idc))) / n^2
  mlqr_var5 <- sum(mlqck_var5 * as.numeric(table(idc))) / n^2
  mlqr_var6 <- sum(mlqck_var6 * as.numeric(table(idc))) / n^2
  
  
  
  
  
  ################################Nonparametric triply machine learning estimator, 5-fold################################
  mltck_num1 <- sapply(results, `[`, 19)
  mltck_num2 <- sapply(results, `[`, 20)
  mltck_num3 <- sapply(results, `[`, 21)
  mltck_num4 <- sapply(results, `[`, 22)
  mltck_num5 <- sapply(results, `[`, 23)
  mltck_num6 <- sapply(results, `[`, 24)
  mltck_denom1 <- sapply(results, `[`, 25)
  mltck_denom2 <- sapply(results, `[`, 26)
  mltck_denom3 <- sapply(results, `[`, 27)
  mltck_denom4 <- sapply(results, `[`, 28)
  mltck_denom5 <- sapply(results, `[`, 29)
  mltck_denom6 <- sapply(results, `[`, 30)
  mltck_var1 <- sapply(results, `[`, 31)
  mltck_var2 <- sapply(results, `[`, 32)
  mltck_var3 <- sapply(results, `[`, 33)
  mltck_var4 <- sapply(results, `[`, 34)
  mltck_var5 <- sapply(results, `[`, 35)
  mltck_var6 <- sapply(results, `[`, 36)
  mltr1 <- (sum(mltck_num1 * as.numeric(table(idc))) / n)/(sum(mltck_denom1 * as.numeric(table(idc))) / n)
  mltr2 <- (sum(mltck_num2 * as.numeric(table(idc))) / n)/(sum(mltck_denom2 * as.numeric(table(idc))) / n)
  mltr3 <- (sum(mltck_num3 * as.numeric(table(idc))) / n)/(sum(mltck_denom3 * as.numeric(table(idc))) / n)
  mltr4 <- (sum(mltck_num4 * as.numeric(table(idc))) / n)/(sum(mltck_denom4 * as.numeric(table(idc))) / n)
  mltr5 <- (sum(mltck_num5 * as.numeric(table(idc))) / n)/(sum(mltck_denom5 * as.numeric(table(idc))) / n)
  mltr6 <- (sum(mltck_num6 * as.numeric(table(idc))) / n)/(sum(mltck_denom6 * as.numeric(table(idc))) / n)
  mltr_var1 <- sum(mltck_var1 * as.numeric(table(idc))) / n^2
  mltr_var2 <- sum(mltck_var2 * as.numeric(table(idc))) / n^2
  mltr_var3 <- sum(mltck_var3 * as.numeric(table(idc))) / n^2
  mltr_var4 <- sum(mltck_var4 * as.numeric(table(idc))) / n^2
  mltr_var5 <- sum(mltck_var5 * as.numeric(table(idc))) / n^2
  mltr_var6 <- sum(mltck_var6 * as.numeric(table(idc))) / n^2
  
  
  return(c(est_weightingc,se_weighting,
           pqrc1,se_qr1,
           ptrc1,se_tr1,
           mlqc1,sqrt(mlqr_var1),
           mltr1,sqrt(mltr_var1),
           pqrc2,se_qr2,
           ptrc2,se_tr2,
           mlqc2,sqrt(mlqr_var2),
           mltr2,sqrt(mltr_var2),
           pqrc3,se_qr3,
           ptrc3,se_tr3,
           mlqc3,sqrt(mlqr_var3),
           mltr3,sqrt(mltr_var3),
           pqrc4,se_qr4,
           ptrc4,se_tr4,
           mlqc4,sqrt(mlqr_var4),
           mltr4,sqrt(mltr_var4),
           pqrc5,se_qr5,
           ptrc5,se_tr5,
           mlqc5,sqrt(mlqr_var5),
           mltr5,sqrt(mltr_var5),
           pqrc6,se_qr6,
           ptrc6,se_tr6,
           mlqc6,sqrt(mlqr_var6),
           mltr6,sqrt(mltr_var6)))
}

#Approximate truth
cal_tr <-  function(n){
  full_data <-  simu_full_data(n)
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
  truepce <- mean((1-D0)*(1-D1)*(Y1-Y0))/mean((1-D0)*(1-D1))
  return(truepce)
}


sim_bias <-  function(n, B, truepce){
  bias <- matrix(NA, nrow = B, ncol = 25)
  aese <- matrix(NA, nrow = B, ncol = 25)
  cov <- matrix(NA, nrow = B, ncol = 25)
  for (j in 1:B){
    tryCatch({
      full_data <-  simu_full_data(n)
      X1 <- full_data[[1]]
      X2 <- full_data[[2]]
      X3 <- full_data[[3]]
      X4 <- full_data[[4]]
      Z <- full_data[[5]]
      D <- full_data[[8]]
      Y <- full_data[[11]]
      df <- as.data.frame(cbind(X1,X2,X3,X4,Z,D,Y))
      start <- Sys.time()
      result_est <- cal_estimates(df, truepce)
      print( Sys.time() - start )
      bias[j,] <- result_est[seq(1,50, by =2)]
      aese[j,] <- result_est[-seq(1,50, by =2)]
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
truepce <- cal_tr(10000000)
n <- 500
B <- 1050
result <- sim_bias(n,B,truepce)
copy_result <- result
mcsd <- round(apply(result[[1]], 2, sd, na.rm = TRUE),3)
mcsd
mbias <- round(colMeans(result[[1]], na.rm = TRUE)-truepce,3)
mbias
aese <- round(colMeans(result[[2]], na.rm = TRUE),3)
aese
cp <- round(colMeans(result[[3]], na.rm = TRUE),3)
cp

out_table <- as.data.frame(cbind(mbias,cp,mcsd,aese))
rownames(out_table) <- c("Weighting",
                         "QR1", "TR1", "MLQR1", "MLTR1",
                         "QR2", "TR2", "MLQR2", "MLTR2",
                         "QR3", "TR3", "MLQR3", "MLTR3",
                         "QR4", "TR4", "MLQR4", "MLTR4",
                         "QR5", "TR5", "MLQR5", "MLTR5",
                         "QR6", "TR6", "MLQR6", "MLTR6")
out_table
setwd("/Users/deckard/Desktop/Fan Li Project/Hayden_Biometrics/Simulation/results")
write.csv(out_table, file= "sim_table_00_case2.csv")

