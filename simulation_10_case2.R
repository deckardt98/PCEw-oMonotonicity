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
  X1 <- rtruncnorm(n, a=-20, b=20, mean = -2, sd = 1)
  X2 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X3 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X1ab <- abs(X1)
  X2ab <- abs(X2)
  X3ab <- abs(X3)
  X4 <- rbinom(n, size = 1, prob = 0.5)
  probZ <- expit(0.75+0.5*(X1+X2+X3+X4))
  #max(abs(probZ-0.5))
  Z <- rbinom(n, size = 1, prob = probZ)
  #mean(Z)
  
  #simulate G=(D(0),D(1))
  probD1 <- expit(-0.6+0.4*X1ab+0.2*X2ab+0.1*X3ab+0.1*X4)
  probD0 <- expit(-1+0.1*X1ab+0.1*X2ab-0.1*X3ab+0.1*X4)
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
  meanY1 <- 1+D1+X1+3*X2+3*X3+3*X4
  meanY0 <- -D0-1.5*X1+2*X2+2*X3-2*X4
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
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)
  
  
  #parametric quadruply robust estimator
  pqr1 <- mean(prinscore0*Z*(1-D)*(Y-om10)/proscore+om10*((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1))))/mean((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1)))
  pqr0 <- mean((1-prinscore1)*(1-Z)*D*(Y-om01)/(1-proscore)+om01*((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1))))/mean((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1)))
  pqrc <- pqr1-pqr0
  denom_est <- mean((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1)))
  num_est <- mean(prinscore0*Z*(1-D)*(Y-om10)/proscore+om10*((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1))))-
    mean((1-prinscore1)*(1-Z)*D*(Y-om01)/(1-proscore)+om01*((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1))))
  
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
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)
  
  #Estimate the denominator E[p0(X)(1-p1(X))]
  denom_est <- mean((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1)))
  
  #Compute variance estimates
  pqr1 <- mean(prinscore0*Z*(1-D)*(Y-om10)/proscore+om10*((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1))))/mean((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1)))
  pqr0 <- mean((1-prinscore1)*(1-Z)*D*(Y-om01)/(1-proscore)+om01*((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1))))/mean((prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1)))
  var_single <- mean((((psi_D0-prinscore0)*(1-prinscore1)*(om10-pqr1)+prinscore0*(psi_Y1_D1-pqr1*(1-psi_D1)))-((1-psi_D1-(1-prinscore1))*prinscore0*(om01-pqr0)+(1-prinscore1)*(psi_YD0-pqr0*psi_D0)))^2)/denom_est^2
  
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
      Z*(1-D)*p0/proscore*(Y-theta[max(p1_pos)+1]),
      (1-Z)*D*(1-p1)/(1-proscore)*(Y-theta[max(p1_pos)+2]),
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
  Xm01 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m01_model))
  Xm10 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m10_model))
  
  pi_pos <- 1:ncol(Xpi)
  p0_pos <- (max(pi_pos) + 1):(max(pi_pos) + ncol(Xp0))
  p1_pos <- (max(p0_pos) + 1):(max(p0_pos) + ncol(Xp1))
  m01_pos <- (max(p1_pos) + 1):(max(p1_pos) + ncol(Xm01))
  m10_pos <- (max(m01_pos) + 1):(max(m01_pos) + ncol(Xm10))
  
  pi_scores <- grab_psiFUN(models$pi_model, data)
  p0_scores <- grab_psiFUN(models$p0_model, data)
  p1_scores <- grab_psiFUN(models$p1_model, data)
  m01_scores <- grab_psiFUN(models$m01_model, data)
  m10_scores <- grab_psiFUN(models$m10_model, data)
  
  function(theta){
    proscore <- plogis(Xpi %*% theta[pi_pos])
    p0 <- plogis(Xp0 %*% theta[p0_pos])
    p1 <- plogis(Xp1 %*% theta[p1_pos])
    om01 <- Xm01 %*% theta[m01_pos]
    om10 <- Xm10 %*% theta[m10_pos]
    psi_D0 <- (1-Z)/(1-proscore)*(D-p0)+p0
    psi_D1 <- Z/proscore*(D-p1)+p1
    psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*p0)+om01*p0
    psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-p1))+om10*(1-p1)
    c(pi_scores(theta[pi_pos]),
      p0_scores(theta[p0_pos])*I(Z==0),
      p1_scores(theta[p1_pos])*I(Z==1),
      m01_scores(theta[m01_pos])*I(Z==0)*I(D==1),
      m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
      (psi_D0-p0)*(1-p1)*(om10-theta[max(m10_pos)+1])+p0*(psi_Y1_D1-theta[max(m10_pos)+1]*(1-psi_D1)),
      (p1-psi_D1)*p0*(om01-theta[max(m10_pos)+2])+(1-p1)*(psi_YD0-theta[max(m10_pos)+2]*psi_D0),
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
  # Combine the data into a single dataframe
  data <- data.frame(Z = Z, D = D, Y = Y, X1 = X1, X2 = X2, X3 = X3, X4 = X4)
  
  ################Estimate nuisance functions by parametric approaches#############
  #Estimate propensity score
  propen_m <- glm(Z ~ X1+X2+X3+X4, data = df, family = binomial)
  proscore <- predict(propen_m, newdata=df, type="response")
  
  #Estimate propensity score with incorrect model specifications
  mis_propen_m <- glm(Z ~ X1*cos(exp(X1/20)), data = df, family = binomial)
  mis_proscore <- predict(mis_propen_m, newdata=df, type="response")
  
  #Estimate principal score
  ps_m1 <- glm(D ~ abs(X1)+abs(X2)+abs(X3)+X4, data = df, subset = (Z==1), family = binomial)
  ps_m0 <- glm(D ~ abs(X1)+abs(X2)+abs(X3)+X4, data = df, subset = (Z==0), family = binomial)
  prinscore1 <-  predict(ps_m1, newdata=df, type="response")
  prinscore0 <-  predict(ps_m0, newdata=df, type="response")
  prinscore <- list(prinscore1,prinscore0)
  
  #Estimate principal score with incorrect model specifications
  mis_ps_m1 <- glm(D ~ X1*cos(exp(X1/20)), data = df, subset = (Z==1), family = binomial)
  mis_ps_m0 <- glm(D ~ X1*cos(exp(X1/20)), data = df, subset = (Z==0), family = binomial)
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
  mis_om_m11 <- lm(Y ~ X1*cos(exp(X1/20)), data = df, subset = (Z==1)&(D==1))
  mis_om_m01 <- lm(Y ~ X1*cos(exp(X1/20)), data = df, subset = (Z==0)&(D==1))
  mis_om_m10 <- lm(Y ~ X1*cos(exp(X1/20)), data = df, subset = (Z==1)&(D==0))
  mis_om_m00 <- lm(Y ~ X1*cos(exp(X1/20)), data = df, subset = (Z==0)&(D==0))
  mis_om11 <-  predict(mis_om_m11, newdata=df, type="response")
  mis_om01 <-  predict(mis_om_m01, newdata=df, type="response") 
  mis_om10 <-  predict(mis_om_m10, newdata=df, type="response")
  mis_om00 <-  predict(mis_om_m00, newdata=df, type="response") 
  mis_om <- list(mis_om11,mis_om01,mis_om10,mis_om00)
  
  ##################################################Weighting###############################################################
  ####Point estimate
  est_weighting1 <- mean(Z*(1-D)*prinscore0/proscore*Y)/mean(Z*(1-D)*prinscore0/proscore)
  est_weighting0 <- mean((1-Z)*D*(1-prinscore1)/(1-proscore)*Y)/mean((1-Z)*D*(1-prinscore1)/(1-proscore))
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
  pqr4 <- pqrc_point(Z,D,Y,proscore,list(mis_prinscore1,prinscore0),list(mis_om11,om01,mis_om10,mis_om00))
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
  create_models <- function(data, pi_formula, p0_formula, p1_formula, m01_formula, m10_formula) {
    list(
      pi_model = glm(pi_formula, data = data, family = binomial),
      p0_model = glm(p0_formula, subset = (Z == 0), data = data, family = binomial),
      p1_model = glm(p1_formula, subset = (Z == 1), data = data, family = binomial),
      m01_model = glm(m01_formula, subset = (Z == 0) & (D == 1), data = data, family = gaussian),
      m10_model = glm(m10_formula, subset = (Z == 1) & (D == 0), data = data, family = gaussian)
    )
  }
  
  estimate_qr <- function(data, pi_formula, p0_formula, p1_formula, m01_formula, m10_formula, roots) {
    models <- create_models(data, pi_formula, p0_formula, p1_formula, m01_formula, m10_formula)
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
                  as.vector(coef(om_m01)), as.vector(coef(om_m10)), pqr11, pqr01, pqrc1),
    theta_qr2 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(mis_om_m01)), as.vector(coef(mis_om_m10)), pqr12, pqr02, pqrc2),
    theta_qr3 = c(as.vector(coef(mis_propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(om_m01)), as.vector(coef(om_m10)), pqr13, pqr03, pqrc3),
    theta_qr4 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(mis_ps_m1)), 
                  as.vector(coef(om_m01)), as.vector(coef(mis_om_m10)), pqr14, pqr04, pqrc4),
    theta_qr5 = c(as.vector(coef(propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(mis_om_m01)), as.vector(coef(om_m10)), pqr15, pqr05, pqrc5),
    theta_qr6 = c(as.vector(coef(mis_propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(mis_ps_m1)), 
                  as.vector(coef(mis_om_m01)), as.vector(coef(mis_om_m10)), pqr16, pqr06, pqrc6)
  )
  
  formula_list <- list(
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, p1 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, 
         m01 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, p1 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, 
         m01 = Y ~ X1*cos(exp(X1/20)), m10 = Y ~ X1*cos(exp(X1/20))),
    list(pi = Z ~ X1*cos(exp(X1/20)), p0 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, p1 = D ~ abs(X1) + abs(X2) + abs(X3) + X4, 
         m01 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ X1 + X2 + X3 + X4, p1 = D ~ X1*cos(exp(X1/20)), 
         m01 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ X1*cos(exp(X1/20))),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ X1*cos(exp(X1/20)), p1 = D ~ X1 + X2 + X3 + X4, 
         m01 = Y ~ X1*cos(exp(X1/20)), m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1*cos(exp(X1/20)), p0 = D ~ X1*cos(exp(X1/20)), p1 = D ~ X1*cos(exp(X1/20)), 
         m01 = Y ~ X1*cos(exp(X1/20)), m10 = Y ~ X1*cos(exp(X1/20)))
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
                            formula_list[[i]]$m01, formula_list[[i]]$m10, theta_list[[i]])
    
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
  ptrc1 <- ptrc2 <- ptrc3 <- ptrc4 <- ptrc5 <- ptrc6 <- 0
  
  se_tr1 <- 0 
  se_tr2 <- 0
  se_tr3 <- 0
  se_tr4 <- 0
  se_tr5 <- 0
  se_tr6 <- 0
  
  
  
  ################################Nonparametric quadruply machine learning estimator, 5-fold###############################
  # Define a function to apply transformations and fit SuperLearner
  fit_superlearner <- function(Y, X_list, indices, family, SLmethods, i, id) {
    lapply(X_list, function(X) {
      SuperLearner(Y[indices], X = X[indices, ], newX = X[which(id == i),], family = family, 
                   SL.library = SLmethods, cvControl = list(V = 5))$SL.predict
    })
  }
  
  # List of transformations of covariates
  transformations <- list(
    function(X1, X2, X3, X4) as.data.frame(cbind(X1, X2, X3, X4)),
    function(X1, X2, X3, X4) as.data.frame(cbind(log(X1 + 40), log(X2 + 40), log(X3 + 40), X4)),
    function(X1, X2, X3, X4) as.data.frame(cbind(cos(X1), cos(X2), cos(X3), cos(as.numeric(X4) - 1)))
  )
  
  # Transform the covariates
  X_list <- lapply(transformations, function(f) f(X1, X2, X3, X4))
  
  
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
    ml_ps_m1 <- fit_superlearner(D, X_list, idx_Z1, family = "binomial", SLmethods, i, idc)
    ml_ps_m0 <- fit_superlearner(D, X_list, idx_Z0, family = "binomial", SLmethods, i, idc)
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
    mlqck2 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                         list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                         list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))
    mlqck_num2 <- mlqck2[4]
    mlqck_denom2 <- mlqck2[5]
    mlqck3 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[3]], 
                         list(ml_ps_m1[[3]], ml_ps_m0[[3]]), 
                         list(ml_om_m11[[3]], ml_om_m01[[3]], ml_om_m10[[3]], ml_om_m00[[3]]))
    mlqck_num3 <- mlqck3[4]
    mlqck_denom3 <- mlqck3[5]
    mlqck_var1 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                            list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                            list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]))
    mlqck_var2 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                            list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                            list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))
    mlqck_var3 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[3]], 
                            list(ml_ps_m1[[3]], ml_ps_m0[[3]]), 
                            list(ml_om_m11[[3]], ml_om_m01[[3]], ml_om_m10[[3]], ml_om_m00[[3]]))
    
    
    
    # Return results as a vector
    c(mlqck_num1, mlqck_num2, mlqck_num3, mlqck_denom1, mlqck_denom2, mlqck_denom3, mlqck_var1, mlqck_var2, mlqck_var3,
      0, 0, 0, 0, 0, 0, 0, 0, 0)
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
  mlqck_denom1 <- sapply(results, `[`, 4)
  mlqck_denom2 <- sapply(results, `[`, 5)
  mlqck_denom3 <- sapply(results, `[`, 6)
  mlqck_var1 <- sapply(results, `[`, 7)
  mlqck_var2 <- sapply(results, `[`, 8)
  mlqck_var3 <- sapply(results, `[`, 9)
  mlqc1 <- (sum(mlqck_num1 * as.numeric(table(idc))) / n)/(sum(mlqck_denom1 * as.numeric(table(idc))) / n)
  mlqc2 <- (sum(mlqck_num2 * as.numeric(table(idc))) / n)/(sum(mlqck_denom2 * as.numeric(table(idc))) / n)
  mlqc3 <- (sum(mlqck_num3 * as.numeric(table(idc))) / n)/(sum(mlqck_denom3 * as.numeric(table(idc))) / n)
  mlqr_var1 <- sum(mlqck_var1 * as.numeric(table(idc))) / n^2
  mlqr_var2 <- sum(mlqck_var2 * as.numeric(table(idc))) / n^2
  mlqr_var3 <- sum(mlqck_var3 * as.numeric(table(idc))) / n^2
  
  
  
  
  
  
  ################################Nonparametric triply machine learning estimator, 5-fold################################
  mltck_num1 <- sapply(results, `[`, 10)
  mltck_num2 <- sapply(results, `[`, 11)
  mltck_num3 <- sapply(results, `[`, 12)
  mltck_denom1 <- sapply(results, `[`, 13)
  mltck_denom2 <- sapply(results, `[`, 14)
  mltck_denom3 <- sapply(results, `[`, 15)
  mltck_var1 <- sapply(results, `[`, 16)
  mltck_var2 <- sapply(results, `[`, 17)
  mltck_var3 <- sapply(results, `[`, 18)
  mltr1 <- (sum(mltck_num1 * as.numeric(table(idc))) / n)/(sum(mltck_denom1 * as.numeric(table(idc))) / n)
  mltr2 <- (sum(mltck_num2 * as.numeric(table(idc))) / n)/(sum(mltck_denom2 * as.numeric(table(idc))) / n)
  mltr3 <- (sum(mltck_num3 * as.numeric(table(idc))) / n)/(sum(mltck_denom3 * as.numeric(table(idc))) / n)
  mltr_var1 <- sum(mltck_var1 * as.numeric(table(idc))) / n^2
  mltr_var2 <- sum(mltck_var2 * as.numeric(table(idc))) / n^2
  mltr_var3 <- sum(mltck_var3 * as.numeric(table(idc))) / n^2
  
  
  return(c(est_weightingc,se_weighting,
           pqrc1,se_qr1,
           pqrc2,se_qr2,
           pqrc3,se_qr3,
           pqrc4,se_qr4,
           pqrc5,se_qr5,
           pqrc6,se_qr6,
           ptrc1,se_tr1,
           ptrc2,se_tr2,
           ptrc3,se_tr3,
           ptrc4,se_tr4,
           ptrc5,se_tr5,
           ptrc6,se_tr6,
           mlqc1,sqrt(mlqr_var1),
           mlqc2,sqrt(mlqr_var2),
           mlqc3,sqrt(mlqr_var3),
           mltr1,sqrt(mltr_var1),
           mltr2,sqrt(mltr_var2),
           mltr3,sqrt(mltr_var3)))
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
  truepce <- mean(D0*(1-D1)*(Y1-Y0))/mean(D0*(1-D1))
  return(truepce)
}


sim_bias <-  function(n, B, truepce){
  bias <- matrix(NA, nrow = B, ncol = 19)
  aese <- matrix(NA, nrow = B, ncol = 19)
  cov <- matrix(NA, nrow = B, ncol = 19)
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
      bias[j,] <- result_est[seq(1,38, by =2)]
      aese[j,] <- result_est[-seq(1,38, by =2)]
      print(round(aese[j,],3))
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


set.seed(20240824)
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

#result[[2]][,c(14:16)]
sort(result[[2]][,19], decreasing = T)
index_badcase <- apply(result[[2]][, c(17:18)], 2, which.max)
for (i in 1:3) {
  result[[i]][index_badcase, ] <- NA
}



out_table <- as.data.frame(cbind(mbias,cp,mcsd,aese))
rownames(out_table) <- c("Weighting",
                         "QR1","QR2","QR3","QR4","QR5","QR6",
                         "TR1","TR2","TR3","TR4","TR5","TR6",
                         "MLQR_raw","MLQR_log","MLQR_cos",
                         "MLTR_raw","MLTR_log","MLTR_cos")
out_table
setwd("/Users/deckard/Desktop/Fan Li Project/Hayden_Biometrics/Simulation/results")
write.csv(out_table, file= "sim_table_10_case2.csv")

