### This file contains simulation code for the article entitled "Quadruply robust 
### principal stratification analysis without monotonicity" by Tong et al..
### CACE

library(truncnorm)
library(SuperLearner)
library(randomForest)
library(geex)
library(parallel)
library(data.table)
library(boot)
library(caret)

setwd("~/Desktop/Fan Li Project/Hayden_Biometrics/Simulation/Code")
source("EST_FUN01.R")

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

cal_estimates <- function(df, theta, truepce){
  
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
  
  # Combine the data into a single dataframe
  data <- data.frame(Z = Z, D = D, Y = Y, X1 = X1, X2 = X2, X3 = X3, X4 = X4)
  
  ################Estimate nuisance functions by parametric approaches#############
  # Estimate propensity score
  propen_m <- glm(Z ~ X1+X2+X3+X4, data = df, family = binomial)
  proscore <- predict(propen_m, newdata=df, type="response")
  
  #Estimate propensity score with incorrect model specifications
  mis_propen_m <- glm(Z ~ mis_X1+mis_X2+mis_X3+mis_X4, data = df, family = binomial)
  mis_proscore <- predict(mis_propen_m, newdata=df, type="response")
  
  # Estimate principal score
  ps_m1 <- glm(D ~ X1+X2+X3+X4, data = df, subset = (Z==1), family = binomial)
  ps_m0 <- glm(D ~ X1+X2+X3+X4, data = df, subset = (Z==0), family = binomial)
  prinscore1 <-  predict(ps_m1, newdata=df, type="response")
  prinscore0 <-  predict(ps_m0, newdata=df, type="response")
  prinscore <- list(prinscore1,prinscore0)
  
  # Estimate principal score with incorrect model specifications
  mis_ps_m1 <- glm(D ~ mis_X1+mis_X2+mis_X3+mis_X4, data = df, subset = (Z==1), family = binomial)
  mis_ps_m0 <- glm(D ~ mis_X1+mis_X2+mis_X3+mis_X4, data = df, subset = (Z==0), family = binomial)
  mis_prinscore1 <-  predict(mis_ps_m1, newdata=df, type="response")
  mis_prinscore0 <-  predict(mis_ps_m0, newdata=df, type="response")    
  mis_prinscore <- list(mis_prinscore1,mis_prinscore0)
  
  # Estimate outcome mean
  om_m11 <- lm(Y ~ X1+X2+X3+X4, data = df, subset = (Z==1)&(D==1))
  om_m01 <- lm(Y ~ X1+X2+X3+X4, data = df, subset = (Z==0)&(D==1))
  om_m10 <- lm(Y ~ X1+X2+X3+X4, data = df, subset = (Z==1)&(D==0))
  om_m00 <- lm(Y ~ X1+X2+X3+X4, data = df, subset = (Z==0)&(D==0))
  om11 <-  predict(om_m11, newdata=df, type="response")
  om01 <-  predict(om_m01, newdata=df, type="response") 
  om10 <-  predict(om_m10, newdata=df, type="response")
  om00 <-  predict(om_m00, newdata=df, type="response") 
  om <- list(om11,om01,om10,om00)
  
  # Estimate outcome mean with incorrect model specifications
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
  #### Point estimate
  # Compute delta(X) and lambda(X)
  if (theta>9999){
    est_weighting1 <- mean(Z*D*(prinscore1-prinscore0)/proscore/prinscore1*Y)/mean(prinscore1-prinscore0)
    est_weighting0 <- mean((1-Z)*(1-D)*(prinscore1-prinscore0)/(1-proscore)/(1-prinscore0)*Y)/mean(prinscore1-prinscore0)
    est_weightingc <- est_weighting1 - est_weighting0
  } else {
    delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
    mis_delta <- (1+(theta-1)*(mis_prinscore1+mis_prinscore0))^2-4*theta*(theta-1)*mis_prinscore1*mis_prinscore0
    
    lambda <- 1+(theta-1)*(prinscore0-prinscore1)-sqrt(delta)
    mis_lambda <- 1+(theta-1)*(mis_prinscore0-mis_prinscore1)-sqrt(mis_delta)
    
    est_weighting1 <- mean(lambda*Z*D*Y/proscore/prinscore1)/mean(lambda*Z*D/proscore/prinscore1)
    est_weighting0 <- mean(lambda*(1-Z)*(1-D)*Y/(1-proscore)/(1-prinscore0))/mean(lambda*(1-Z)*(1-D)/(1-proscore)/(1-prinscore0))
    est_weightingc <- est_weighting1 - est_weighting0
  }
  
  #### Variance estimate
  theta_ps <- c(as.vector(coef(propen_m)),
                as.vector(coef(ps_m0)),
                as.vector(coef(ps_m1)),
                est_weighting1,
                est_weighting0,
                est_weightingc)
  estimate_psw <- function(data, theta){
    pi_model <- glm(Z~X1+X2+X3+X4, data=data, family = binomial)
    p0_model <- glm(D~X1+X2+X3+X4, subset = (Z==0), data=data, family = binomial)
    p1_model <- glm(D~X1+X2+X3+X4, subset = (Z==1), data=data, family = binomial) 
    models <- list(pi_model=pi_model, p0_model=p0_model, p1_model=p1_model)
    m_estimate(estFUN = eq_psw, data = data,
               outer_args = list(models = models, oddsratio = theta),
               #root_control = setup_root_control(start = theta_ps+rnorm(length(theta_ps),mean=0,sd=0.5))
               compute_roots = FALSE,
               roots = theta_ps,
               deriv_control = setup_deriv_control(
                 FUN = custom_jacobian, 
                 method = "simple"
               )
    )
  }
  cov_matrix_psw <- vcov(estimate_psw(data=df, theta))
  se_weighting <- sqrt(cov_matrix_psw[ncol(cov_matrix_psw),ncol(cov_matrix_psw)])
  
  
  
  
  theta_com <- c(0.5,0.999999,2)
  ###################################################
  # Parametric multiply robust estimator
  ###################################################
  
  # Function to calculate point estimates for different model specifications
  calculate_pqrc <- function(proscore, prinscore, om, theta) {
    pqrc_point(Z, D, Y, proscore, prinscore, om, theta)
  }
  
  # Function to calculate and return the results for all working design matrices
  calculate_all_cases <- function(theta) {
    list(
      all_correct = pqrc_point(Z, D, Y, proscore, prinscore, om, theta),
      pi_p = pqrc_point(Z, D, Y, proscore, prinscore, mis_om, theta),
      p_om = pqrc_point(Z, D, Y, mis_proscore, prinscore, om, theta),
      pi_om = pqrc_point(Z, D, Y, proscore, mis_prinscore, om, theta),
      all_mis = pqrc_point(Z, D, Y, mis_proscore, mis_prinscore, mis_om, theta)
    )
  }
  
  results <- lapply(theta_com, calculate_all_cases)
  
  # Extract results for each case (for all theta_com values)
  # Initialize matrices to store the results for all cases
  pqr11 <- pqr01 <- pqrc1 <- numeric(length(theta_com))
  pqr12 <- pqr02 <- pqrc2 <- numeric(length(theta_com))
  pqr13 <- pqr03 <- pqrc3 <- numeric(length(theta_com))
  pqr14 <- pqr04 <- pqrc4 <- numeric(length(theta_com))
  pqr15 <- pqr05 <- pqrc5 <- numeric(length(theta_com))
  
  for (i in seq_along(theta_com)) {
    # All models correctly specified
    pqr11[i] <- results[[i]]$all_correct[1]
    pqr01[i] <- results[[i]]$all_correct[2]
    pqrc1[i] <- results[[i]]$all_correct[3]
    
    # pi+p (mis_om)
    pqr12[i] <- results[[i]]$pi_p[1]
    pqr02[i] <- results[[i]]$pi_p[2]
    pqrc2[i] <- results[[i]]$pi_p[3]
    
    # p+om (mis_proscore)
    pqr13[i] <- results[[i]]$p_om[1]
    pqr03[i] <- results[[i]]$p_om[2]
    pqrc3[i] <- results[[i]]$p_om[3]
    
    # pi+om (mis_prinscore)
    pqr14[i] <- results[[i]]$pi_om[1]
    pqr04[i] <- results[[i]]$pi_om[2]
    pqrc4[i] <- results[[i]]$pi_om[3]
    
    # All models wrongly specified (mis_proscore, mis_prinscore, mis_om)
    pqr15[i] <- results[[i]]$all_mis[1]
    pqr05[i] <- results[[i]]$all_mis[2]
    pqrc5[i] <- results[[i]]$all_mis[3]
  }
  
  ## Sandwich variance estimates
  create_models <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m11_formula) {
    list(
      pi_model = glm(pi_formula, data = data, family = binomial),
      p0_model = glm(p0_formula, subset = (Z == 0), data = data, family = binomial),
      p1_model = glm(p1_formula, subset = (Z == 1), data = data, family = binomial),
      m00_model = glm(m00_formula, subset = (Z == 0) & (D == 0), data = data, family = gaussian),
      m11_model = glm(m11_formula, subset = (Z == 1) & (D == 1), data = data, family = gaussian)
    )
  }
  
  estimate_qr <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m11_formula, roots, oddsratio) {
    models <- create_models(data, pi_formula, p0_formula, p1_formula, m00_formula, m11_formula)
    m_estimate(
      estFUN = eq_qr, 
      data = data,
      outer_args = list(models = models, oddsratio = oddsratio),
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
  theta_list <- lapply(seq_along(theta_com), function(i) {
    list(
      theta_qr1 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                    as.vector(coef(om_m00)), as.vector(coef(om_m11)), pqr11[i], pqr01[i], pqrc1[i]),
      theta_qr2 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                    as.vector(coef(mis_om_m00)), as.vector(coef(mis_om_m11)), pqr12[i], pqr02[i], pqrc2[i]),
      theta_qr3 = c(as.vector(coef(mis_propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                    as.vector(coef(om_m00)), as.vector(coef(om_m11)), pqr13[i], pqr03[i], pqrc3[i]),
      theta_qr4 = c(as.vector(coef(propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(mis_ps_m1)), 
                    as.vector(coef(om_m00)), as.vector(coef(om_m11)), pqr14[i], pqr04[i], pqrc4[i]),
      theta_qr5 = c(as.vector(coef(mis_propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(mis_ps_m1)), 
                    as.vector(coef(mis_om_m00)), as.vector(coef(mis_om_m11)), pqr15[i], pqr05[i], pqrc5[i])
    )
  })
  
  formula_list <- list(
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ X1 + X2 + X3 + X4, p1 = D ~ X1 + X2 + X3 + X4, 
         m00 = Y ~ X1 + X2 + X3 + X4, m11 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ X1 + X2 + X3 + X4, p1 = D ~ X1 + X2 + X3 + X4, 
         m00 = Y ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, m11 = Y ~ mis_X1 + mis_X2 + mis_X3 + mis_X4),
    list(pi = Z ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, p0 = D ~ X1 + X2 + X3 + X4, p1 = D ~ X1 + X2 + X3 + X4, 
         m00 = Y ~ X1 + X2 + X3 + X4, m11 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, p1 = D ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, 
         m00 = Y ~ X1 + X2 + X3 + X4, m11 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, p0 = D ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, p1 = D ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, 
         m00 = Y ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, m11 = Y ~ mis_X1 + mis_X2 + mis_X3 + mis_X4)
  )
  
  # Function to calculate estimates and standard errors for all theta_com values
  se_list <- mclapply(c(1:5), function(i) {
    lapply(c(1:3), function(j) {
      # Compute the estimate (but do not keep it)
      estimate <- estimate_qr(df, formula_list[[i]]$pi, formula_list[[i]]$p0, formula_list[[i]]$p1, 
                              formula_list[[i]]$m00, formula_list[[i]]$m11, theta_list[[j]][[i]], oddsratio = theta_com[j])
      
      # Compute the standard error
      cov_matrix <- vcov(estimate)
      sqrt(cov_matrix[ncol(cov_matrix), ncol(cov_matrix)])
    })
  }, mc.cores = 5)
  
  # Organize the results into separate standard errors for each model and theta_com
  se_qr_list <- lapply(1:5, function(i) {
    sapply(1:length(theta_com), function(j) se_list[[i]][[j]])
  })
  
  
  # Assign standard errors to individual variables for all theta_com values
  se_qr1 <- se_qr_list[[1]]
  se_qr2 <- se_qr_list[[2]]
  se_qr3 <- se_qr_list[[3]]
  se_qr4 <- se_qr_list[[4]]
  se_qr5 <- se_qr_list[[5]]
  
  
  ###################################################
  # theta = infty, under monotonicity
  ###################################################
  ## Point estimates
  # All models are correctly specified
  ptr1 <- ptrc_point(Z,D,Y,proscore,prinscore,om)
  ptr11 <- ptr1[1]
  ptr01 <- ptr1[2]
  ptrc1 <- ptr1[3]
  # pi+p
  ptr2 <- ptrc_point(Z,D,Y,proscore,prinscore,mis_om)
  ptr12 <- ptr2[1]
  ptr02 <- ptr2[2]
  ptrc2 <- ptr2[3]
  # p+om
  ptr3 <- ptrc_point(Z,D,Y,mis_proscore,prinscore,om)
  ptr13 <- ptr3[1]
  ptr03 <- ptr3[2]
  ptrc3 <- ptr3[3]
  # pi+om
  ptr4 <- ptrc_point(Z,D,Y,proscore,mis_proscore,om)
  ptr14 <- ptr4[1]
  ptr04 <- ptr4[2]
  ptrc4 <- ptr4[3]
  # All models are wrongly specified
  ptr5 <- ptrc_point(Z,D,Y,mis_proscore,mis_prinscore,mis_om)
  ptr15 <- ptr5[1]
  ptr05 <- ptr5[2]
  ptrc5 <- ptr5[3]
  
  ### Sandwich variance estimator
  estimate_tr <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m11_formula, roots) {
    models <- create_models(data, pi_formula, p0_formula, p1_formula, m00_formula, m11_formula)
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
                  as.vector(coef(om_m00)), as.vector(coef(om_m11)), ptr11, ptr01, ptrc1),
    theta_tr2 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(mis_om_m00)), as.vector(coef(mis_om_m11)), ptr12, ptr02, ptrc2),
    theta_tr3 = c(as.vector(coef(mis_propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                  as.vector(coef(om_m00)), as.vector(coef(om_m11)), ptr13, ptr03, ptrc3),
    theta_tr4 = c(as.vector(coef(propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(mis_ps_m1)), 
                  as.vector(coef(om_m00)), as.vector(coef(om_m11)), ptr14, ptr04, ptrc4),
    theta_tr5 = c(as.vector(coef(mis_propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(mis_ps_m1)), 
                  as.vector(coef(mis_om_m00)), as.vector(coef(mis_om_m11)), ptr15, ptr05, ptrc5)
  )
  
  formula_list_tr <- list(
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ X1 + X2 + X3 + X4, p1 = D ~ X1 + X2 + X3 + X4, 
         m00 = Y ~ X1 + X2 + X3 + X4, m11 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ X1 + X2 + X3 + X4, p1 = D ~ X1 + X2 + X3 + X4, 
         m00 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4, m11 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4),
    list(pi = Z ~ mis_X1+mis_X2+mis_X3+mis_X4, p0 = D ~ X1 + X2 + X3 + X4, p1 = D ~ X1 + X2 + X3 + X4, 
         m00 = Y ~ X1 + X2 + X3 + X4, m11 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ mis_X1+mis_X2+mis_X3+mis_X4, p1 = D ~ mis_X1+mis_X2+mis_X3+mis_X4, 
         m00 = Y ~ X1 + X2 + X3 + X4, m11 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ mis_X1+mis_X2+mis_X3+mis_X4, p0 = D ~ mis_X1+mis_X2+mis_X3+mis_X4, p1 = D ~ mis_X1+mis_X2+mis_X3+mis_X4, 
         m00 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4, m11 = Y ~ mis_X1+mis_X2+mis_X3+mis_X4)
  )
  
  
  se_list_tr <- mclapply(seq_along(theta_list_tr), function(i) {
    # Compute the estimate (but do not keep it)
    estimate <- estimate_tr(df, formula_list_tr[[i]]$pi, formula_list_tr[[i]]$p0, formula_list_tr[[i]]$p1, 
                            formula_list_tr[[i]]$m00, formula_list_tr[[i]]$m11, theta_list_tr[[i]])
    
    # Compute the standard error
    cov_matrix <- vcov(estimate)
    sqrt(cov_matrix[ncol(cov_matrix), ncol(cov_matrix)])
  }, mc.cores = 5)
  
  se_tr1 <- se_list_tr[[1]]
  se_tr2 <- se_list_tr[[2]]
  se_tr3 <- se_list_tr[[3]]
  se_tr4 <- se_list_tr[[4]]
  se_tr5 <- se_list_tr[[5]]
  
  
  
  ###################################################
  # Multiply machine learning estimator-5-fold
  ###################################################
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
  
  # Function to execute for each fold
  process_fold <- function(i, idc, theta_com) {
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
    
    # Initialize result vectors
    num_results <- numeric(length(theta_com) * 5)  # 5 numerators
    denom_results <- numeric(length(theta_com) * 5)  # 5 denominators
    var_results <- numeric(length(theta_com) * 5)  # 5 variances
    
    # Loop through each value of theta_com
    for (j in seq_along(theta_com)) {
      # ML0.5 for each theta_com
      mlqck1 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                           list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                           list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j])
      num_results[j] <- mlqck1[4]
      denom_results[j] <- mlqck1[5]
      
      mlqck2 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                           list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                           list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]), theta_com[j])
      num_results[length(theta_com) + j] <- mlqck2[4]
      denom_results[length(theta_com) + j] <- mlqck2[5]
      
      mlqck3 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                           list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                           list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j])
      num_results[2 * length(theta_com) + j] <- mlqck3[4]
      denom_results[2 * length(theta_com) + j] <- mlqck3[5]
      
      mlqck4 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                           list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                           list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j])
      num_results[3 * length(theta_com) + j] <- mlqck4[4]
      denom_results[3 * length(theta_com) + j] <- mlqck4[5]
      
      mlqck5 <- pqrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                           list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                           list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]), theta_com[j])
      num_results[4 * length(theta_com) + j] <- mlqck5[4]
      denom_results[4 * length(theta_com) + j] <- mlqck5[5]
      
      # Variance calculations for each theta_com
      mlqck_var1 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                              list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                              list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j])
      
      mlqck_var2 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                              list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                              list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]), theta_com[j])
      
      mlqck_var3 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                              list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                              list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j]);
      
      mlqck_var4 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                              list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                              list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j]);
      
      mlqck_var5 <- ml_qr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                              list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                              list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]), theta_com[j]);
      
      # Store variances
      var_results[j] <- mlqck_var1
      var_results[length(theta_com) + j] <- mlqck_var2
      var_results[2 * length(theta_com) + j] <- mlqck_var3
      var_results[3 * length(theta_com) + j] <- mlqck_var4
      var_results[4 * length(theta_com) + j] <- mlqck_var5
    }
    
    # MLinfty
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
                         list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                         list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]))
    mltck_num4 <- mltck4[4]
    mltck_denom4 <- mltck4[5]
    
    mltck5 <- ptrc_point(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                         list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                         list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))
    mltck_num5 <- mltck5[4]
    mltck_denom5 <- mltck5[5]
    
    
    
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
                            list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                            list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]))
    
    mltck_var5 <- ml_tr_var(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                            list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                            list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]))   
    
    # Return results as a vector
    c(num_results, denom_results, var_results,
      mltck_num1, mltck_num2, mltck_num3, mltck_num4, mltck_num5,
      mltck_denom1, mltck_denom2, mltck_denom3, mltck_denom4, mltck_denom5,
      mltck_var1, mltck_var2, mltck_var3, mltck_var4, mltck_var5)
  }
  ml_error <- TRUE
  n_fold <- 5
  while (ml_error) {
    idc <- gene_Fi(df,n_fold)
    tryCatch({
      results <- mclapply(1:n_fold, process_fold, mc.cores = 5, idc = idc, theta_com = theta_com)
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
  mlqck_num1 <- sapply(results, `[`, c(1:3))
  mlqck_num2 <- sapply(results, `[`, c(4:6))
  mlqck_num3 <- sapply(results, `[`, c(7:9))
  mlqck_num4 <- sapply(results, `[`, c(10:12))
  mlqck_num5 <- sapply(results, `[`, c(13:15))
  
  mlqck_denom1 <- sapply(results, `[`, c(16:18))
  mlqck_denom2 <- sapply(results, `[`, c(19:21))
  mlqck_denom3 <- sapply(results, `[`, c(22:24))
  mlqck_denom4 <- sapply(results, `[`, c(25:27))
  mlqck_denom5 <- sapply(results, `[`, c(28:30))
  
  mlqck_var1 <- sapply(results, `[`, c(31:33))
  mlqck_var2 <- sapply(results, `[`, c(34:36))
  mlqck_var3 <- sapply(results, `[`, c(37:39))
  mlqck_var4 <- sapply(results, `[`, c(40:42))
  mlqck_var5 <- sapply(results, `[`, c(43:45))
  
  mltck_num1 <- sapply(results, `[`, 46)
  mltck_num2 <- sapply(results, `[`, 47)
  mltck_num3 <- sapply(results, `[`, 48)
  mltck_num4 <- sapply(results, `[`, 49)
  mltck_num5 <- sapply(results, `[`, 50)
  
  mltck_denom1 <- sapply(results, `[`, 51)
  mltck_denom2 <- sapply(results, `[`, 52)
  mltck_denom3 <- sapply(results, `[`, 53)
  mltck_denom4 <- sapply(results, `[`, 54)
  mltck_denom5 <- sapply(results, `[`, 55)
  
  mltck_var1 <- sapply(results, `[`, 56)
  mltck_var2 <- sapply(results, `[`, 57)
  mltck_var3 <- sapply(results, `[`, 58)
  mltck_var4 <- sapply(results, `[`, 59)
  mltck_var5 <- sapply(results, `[`, 60)
  
  weig <- as.numeric(table(idc)) / n
  mlqc1 <- apply(mlqck_num1, 1, function(row) sum(row * weig))/apply(mlqck_denom1, 1, function(row) sum(row * weig))
  mlqc2 <- apply(mlqck_num2, 1, function(row) sum(row * weig))/apply(mlqck_denom2, 1, function(row) sum(row * weig))
  mlqc3 <- apply(mlqck_num3, 1, function(row) sum(row * weig))/apply(mlqck_denom3, 1, function(row) sum(row * weig))
  mlqc4 <- apply(mlqck_num4, 1, function(row) sum(row * weig))/apply(mlqck_denom4, 1, function(row) sum(row * weig))
  mlqc5 <- apply(mlqck_num5, 1, function(row) sum(row * weig))/apply(mlqck_denom5, 1, function(row) sum(row * weig))
  
  
  
  mlqr_var1 <- apply(mlqck_var1, 1, function(row) sum(row * weig))/n
  mlqr_var2 <- apply(mlqck_var2, 1, function(row) sum(row * weig))/n
  mlqr_var3 <- apply(mlqck_var3, 1, function(row) sum(row * weig))/n
  mlqr_var4 <- apply(mlqck_var4, 1, function(row) sum(row * weig))/n
  mlqr_var5 <- apply(mlqck_var5, 1, function(row) sum(row * weig))/n
  
  mltr1 <- (sum(mltck_num1 * as.numeric(table(idc))) / n)/(sum(mltck_denom1 * as.numeric(table(idc))) / n)
  mltr2 <- (sum(mltck_num2 * as.numeric(table(idc))) / n)/(sum(mltck_denom2 * as.numeric(table(idc))) / n)
  mltr3 <- (sum(mltck_num3 * as.numeric(table(idc))) / n)/(sum(mltck_denom3 * as.numeric(table(idc))) / n)
  mltr4 <- (sum(mltck_num4 * as.numeric(table(idc))) / n)/(sum(mltck_denom4 * as.numeric(table(idc))) / n)
  mltr5 <- (sum(mltck_num5 * as.numeric(table(idc))) / n)/(sum(mltck_denom5 * as.numeric(table(idc))) / n)
  
  mltr_var1 <- sum(mltck_var1 * as.numeric(table(idc))) / n^2
  mltr_var2 <- sum(mltck_var2 * as.numeric(table(idc))) / n^2
  mltr_var3 <- sum(mltck_var3 * as.numeric(table(idc))) / n^2
  mltr_var4 <- sum(mltck_var4 * as.numeric(table(idc))) / n^2
  mltr_var5 <- sum(mltck_var5 * as.numeric(table(idc))) / n^2
  
  
  
  
  return(c(est_weightingc,se_weighting,
           pqrc1[1],se_qr1[1],
           pqrc1[2],se_qr1[2],
           pqrc1[3],se_qr1[3],
           ptrc1,se_tr1,
           mlqc1[1],sqrt(mlqr_var1)[1],
           mlqc1[2],sqrt(mlqr_var1)[2],
           mlqc1[3],sqrt(mlqr_var1)[3],
           mltr1,sqrt(mltr_var1),
           pqrc2[1],se_qr2[1],
           pqrc2[2],se_qr2[2],
           pqrc2[3],se_qr2[3],
           ptrc2,se_tr2,
           mlqc2[1],sqrt(mlqr_var2)[1],
           mlqc2[2],sqrt(mlqr_var2)[2],
           mlqc2[3],sqrt(mlqr_var2)[3],
           mltr2,sqrt(mltr_var2),
           pqrc3[1],se_qr3[1],
           pqrc3[2],se_qr3[2],
           pqrc3[3],se_qr3[3],
           ptrc3,se_tr3,
           mlqc3[1],sqrt(mlqr_var3)[1],
           mlqc3[2],sqrt(mlqr_var3)[2],
           mlqc3[3],sqrt(mlqr_var3)[3],
           mltr3,sqrt(mltr_var3),
           pqrc4[1],se_qr4[1],
           pqrc4[2],se_qr4[2],
           pqrc4[3],se_qr4[3],
           ptrc4,se_tr4,
           mlqc4[1],sqrt(mlqr_var4)[1],
           mlqc4[2],sqrt(mlqr_var4)[2],
           mlqc4[3],sqrt(mlqr_var4)[3],
           mltr4,sqrt(mltr_var4),
           pqrc5[1],se_qr5[1],
           pqrc5[2],se_qr5[2],
           pqrc5[3],se_qr5[3],
           ptrc5,se_tr5,
           mlqc5[1],sqrt(mlqr_var5)[1],
           mlqc5[2],sqrt(mlqr_var5)[2],
           mlqc5[3],sqrt(mlqr_var5)[3],
           mltr5,sqrt(mltr_var5)))
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
  truepce <- mean((1-D0)*D1*(Y1-Y0))/mean((1-D0)*D1)
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
theta <- 0.999999
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
write.csv(out_table, file= "01table_new_1.csv")

