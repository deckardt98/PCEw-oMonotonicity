########################Estimating Equations for defiers###########################
#Function for calculating parametric multiply robust estimator for theta<infty
pqrc_point10 <- function(Z,D,Y,proscore,prinscore,om,theta){
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
  
  #tau and omega
  delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
  lambda <- 1+(theta-1)*(prinscore1-prinscore0)-sqrt(delta)
  tau <- lambda + (theta-1)/sqrt(delta)*((psi_D0-prinscore0)*(2*((theta-1)*(prinscore1-prinscore0)+prinscore1-sqrt(delta))-lambda)+(psi_D1-prinscore1)*(2*prinscore0-lambda))
  omega1 <- lambda/(1-prinscore1)*(psi_Y1_D1-om10*(1-psi_D1))+tau*om10
  omega0 <- lambda/prinscore0*(psi_YD0-om01*psi_D0)+tau*om01
  
  #parametric multiply robust estimator
  pqr1 <- mean(omega1)/mean(tau)
  pqr0 <- mean(omega0)/mean(tau)
  pqrc <- pqr1-pqr0
  denom_est <- mean(tau)
  num_est <- mean(omega1-omega0)
  
  return(c(pqr1,pqr0,pqrc,num_est,denom_est))
}

#Function for calculating the variance of the multiply machine learning estimator for theta<infty
ml_qr_var10 <- function(Z,D,Y,proscore,prinscore,om,theta){
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
  
  #tau and omega
  delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
  lambda <- 1+(theta-1)*(prinscore1-prinscore0)-sqrt(delta)
  tau <- lambda + (theta-1)/sqrt(delta)*((psi_D0-prinscore0)*(2*((theta-1)*(prinscore1-prinscore0)+prinscore1-sqrt(delta))-lambda)+(psi_D1-prinscore1)*(2*prinscore0-lambda))
  omega1 <- lambda/(1-prinscore1)*(psi_Y1_D1-om10*(1-psi_D1))+tau*om10
  omega0 <- lambda/prinscore0*(psi_YD0-om01*psi_D0)+tau*om01
  
  
  #Estimate the denominator 
  denom_est <- mean(tau)
  
  #Compute the variance estimates
  pqr1 <- mean(omega1)/mean(tau)
  pqr0 <- mean(omega0)/mean(tau)
  xi1 <- omega1-pqr1*tau
  xi0 <- omega0-pqr0*tau  
  var_single <- mean((xi1-xi0)^2)/denom_est^2
  
  return(var_single)
}

#Estimating equations for weighing estimator
eq_psw10 <- function(data, models, oddsratio){
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
      delta <- (1+(oddsratio-1)*(p1+p0))^2-4*oddsratio*(oddsratio-1)*p1*p0
      lambda <- 1+(oddsratio-1)*(p1-p0)-sqrt(delta)
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        Z*(1-D)*lambda/proscore/(1-p1)*(Y-theta[max(p1_pos)+1]),
        (1-Z)*D*lambda/(1-proscore)/p0*(Y-theta[max(p1_pos)+2]),
        theta[max(p1_pos)+1]-theta[max(p1_pos)+2]-theta[max(p1_pos)+3]) 
    }
}

#Estimating equations for parametric multiply robust estimator
eq_qr10 <- function(data, models, oddsratio){
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
    delta <- (1+(oddsratio-1)*(p1+p0))^2-4*oddsratio*(oddsratio-1)*p1*p0
    lambda <- 1+(oddsratio-1)*(p1-p0)-sqrt(delta)
    tau <- lambda + (oddsratio-1)/sqrt(delta)*((psi_D0-p0)*(2*((oddsratio-1)*(p1-p0)+p1-sqrt(delta))-lambda)+(psi_D1-p1)*(2*p0-lambda))
    omega1 <- lambda/(1-p1)*(psi_Y1_D1-om10*(1-psi_D1))+tau*om10
    omega0 <- lambda/p0*(psi_YD0-om01*psi_D0)+tau*om01
    c(pi_scores(theta[pi_pos]),
      p0_scores(theta[p0_pos])*I(Z==0),
      p1_scores(theta[p1_pos])*I(Z==1),
      m01_scores(theta[m01_pos])*I(Z==0)*I(D==1),
      m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
      omega1-theta[max(m10_pos)+1]*tau,
      omega0-theta[max(m10_pos)+2]*tau,
      theta[max(m10_pos)+1]-theta[max(m10_pos)+2]-theta[max(m10_pos)+3]) 
  }
}

cal_estimates10 <- function(df, theta, truepce){
  
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
    est_weighting1 <- 0
    est_weighting0 <- 0
    est_weightingc <- 0
    se_weighting <- 0
  } else {
    delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
    mis_delta <- (1+(theta-1)*(mis_prinscore1+mis_prinscore0))^2-4*theta*(theta-1)*mis_prinscore1*mis_prinscore0
    
    lambda <- 1+(theta-1)*(prinscore1-prinscore0)-sqrt(delta)
    mis_lambda <- 1+(theta-1)*(mis_prinscore1-mis_prinscore0)-sqrt(mis_delta)
    
    est_weighting1 <- mean(lambda*Z*(1-D)*Y/proscore/(1-prinscore1))/mean(lambda*Z*(1-D)/proscore/(1-prinscore1))
    est_weighting0 <- mean(lambda*(1-Z)*D*Y/(1-proscore)/prinscore0)/mean(lambda*(1-Z)*D/(1-proscore)/prinscore0)
    est_weightingc <- est_weighting1 - est_weighting0
    
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
      m_estimate(estFUN = eq_psw10, data = data,
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
  }
  
  
  
  
  
  theta_com <- c(0.5,0.999999,2)
  ###################################################
  # Parametric multiply robust estimator
  ###################################################
  
  # Function to calculate point estimates for different model specifications
  calculate_pqrc <- function(proscore, prinscore, om, theta) {
    pqrc_point10(Z, D, Y, proscore, prinscore, om, theta)
  }
  
  # Function to calculate and return the results for all working design matrices
  calculate_all_cases <- function(theta) {
    list(
      all_correct = pqrc_point10(Z, D, Y, proscore, prinscore, om, theta),
      pi_p = pqrc_point10(Z, D, Y, proscore, prinscore, mis_om, theta),
      p_om = pqrc_point10(Z, D, Y, mis_proscore, prinscore, om, theta),
      pi_om = pqrc_point10(Z, D, Y, proscore, mis_prinscore, om, theta),
      all_mis = pqrc_point10(Z, D, Y, mis_proscore, mis_prinscore, mis_om, theta)
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
  create_models <- function(data, pi_formula, p0_formula, p1_formula, m01_formula, m10_formula) {
    list(
      pi_model = glm(pi_formula, data = data, family = binomial),
      p0_model = glm(p0_formula, subset = (Z == 0), data = data, family = binomial),
      p1_model = glm(p1_formula, subset = (Z == 1), data = data, family = binomial),
      m01_model = glm(m01_formula, subset = (Z == 0) & (D == 1), data = data, family = gaussian),
      m10_model = glm(m10_formula, subset = (Z == 1) & (D == 0), data = data, family = gaussian)
    )
  }
  
  estimate_qr <- function(data, pi_formula, p0_formula, p1_formula, m01_formula, m10_formula, roots, oddsratio) {
    models <- create_models(data, pi_formula, p0_formula, p1_formula, m01_formula, m10_formula)
    m_estimate(
      estFUN = eq_qr10, 
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
                    as.vector(coef(om_m01)), as.vector(coef(om_m10)), pqr11[i], pqr01[i], pqrc1[i]),
      theta_qr2 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                    as.vector(coef(mis_om_m01)), as.vector(coef(mis_om_m10)), pqr12[i], pqr02[i], pqrc2[i]),
      theta_qr3 = c(as.vector(coef(mis_propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)), 
                    as.vector(coef(om_m01)), as.vector(coef(om_m10)), pqr13[i], pqr03[i], pqrc3[i]),
      theta_qr4 = c(as.vector(coef(propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(mis_ps_m1)), 
                    as.vector(coef(om_m01)), as.vector(coef(om_m10)), pqr14[i], pqr04[i], pqrc4[i]),
      theta_qr5 = c(as.vector(coef(mis_propen_m)), as.vector(coef(mis_ps_m0)), as.vector(coef(mis_ps_m1)), 
                    as.vector(coef(mis_om_m01)), as.vector(coef(mis_om_m10)), pqr15[i], pqr05[i], pqrc5[i])
    )
  })
  
  formula_list <- list(
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ X1 + X2 + X3 + X4, p1 = D ~ X1 + X2 + X3 + X4, 
         m01 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ X1 + X2 + X3 + X4, p1 = D ~ X1 + X2 + X3 + X4, 
         m01 = Y ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, m10 = Y ~ mis_X1 + mis_X2 + mis_X3 + mis_X4),
    list(pi = Z ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, p0 = D ~ X1 + X2 + X3 + X4, p1 = D ~ X1 + X2 + X3 + X4, 
         m01 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ X1 + X2 + X3 + X4, p0 = D ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, p1 = D ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, 
         m01 = Y ~ X1 + X2 + X3 + X4, m10 = Y ~ X1 + X2 + X3 + X4),
    list(pi = Z ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, p0 = D ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, p1 = D ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, 
         m01 = Y ~ mis_X1 + mis_X2 + mis_X3 + mis_X4, m10 = Y ~ mis_X1 + mis_X2 + mis_X3 + mis_X4)
  )
  
  # Function to calculate estimates and standard errors for all theta_com values
  se_list <- mclapply(c(1:5), function(i) {
    lapply(c(1:3), function(j) {
      # Compute the estimate (but do not keep it)
      estimate <- estimate_qr(df, formula_list[[i]]$pi, formula_list[[i]]$p0, formula_list[[i]]$p1, 
                              formula_list[[i]]$m01, formula_list[[i]]$m10, theta_list[[j]][[i]], oddsratio = theta_com[j])
      
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
      mlqck1 <- pqrc_point10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                           list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                           list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j])
      num_results[j] <- mlqck1[4]
      denom_results[j] <- mlqck1[5]
      
      mlqck2 <- pqrc_point10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                           list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                           list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]), theta_com[j])
      num_results[length(theta_com) + j] <- mlqck2[4]
      denom_results[length(theta_com) + j] <- mlqck2[5]
      
      mlqck3 <- pqrc_point10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                           list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                           list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j])
      num_results[2 * length(theta_com) + j] <- mlqck3[4]
      denom_results[2 * length(theta_com) + j] <- mlqck3[5]
      
      mlqck4 <- pqrc_point10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                           list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                           list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j])
      num_results[3 * length(theta_com) + j] <- mlqck4[4]
      denom_results[3 * length(theta_com) + j] <- mlqck4[5]
      
      mlqck5 <- pqrc_point10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                           list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                           list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]), theta_com[j])
      num_results[4 * length(theta_com) + j] <- mlqck5[4]
      denom_results[4 * length(theta_com) + j] <- mlqck5[5]
      
      # Variance calculations for each theta_com
      mlqck_var1 <- ml_qr_var10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                              list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                              list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j])
      
      mlqck_var2 <- ml_qr_var10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                              list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                              list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]), theta_com[j])
      
      mlqck_var3 <- ml_qr_var10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                              list(ml_ps_m1[[1]], ml_ps_m0[[1]]), 
                              list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j]);
      
      mlqck_var4 <- ml_qr_var10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[1]], 
                              list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                              list(ml_om_m11[[1]], ml_om_m01[[1]], ml_om_m10[[1]], ml_om_m00[[1]]), theta_com[j]);
      
      mlqck_var5 <- ml_qr_var10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore[[2]], 
                              list(ml_ps_m1[[2]], ml_ps_m0[[2]]), 
                              list(ml_om_m11[[2]], ml_om_m01[[2]], ml_om_m10[[2]], ml_om_m00[[2]]), theta_com[j]);
      
      # Store variances
      var_results[j] <- mlqck_var1
      var_results[length(theta_com) + j] <- mlqck_var2
      var_results[2 * length(theta_com) + j] <- mlqck_var3
      var_results[3 * length(theta_com) + j] <- mlqck_var4
      var_results[4 * length(theta_com) + j] <- mlqck_var5
    }
    
    # Return results as a vector
    c(num_results, denom_results, var_results)
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
  
  
  
  
  return(c(est_weightingc,se_weighting,
           pqrc1[1],se_qr1[1],
           pqrc1[2],se_qr1[2],
           pqrc1[3],se_qr1[3],
           0,0,
           mlqc1[1],sqrt(mlqr_var1)[1],
           mlqc1[2],sqrt(mlqr_var1)[2],
           mlqc1[3],sqrt(mlqr_var1)[3],
           0,0,
           pqrc2[1],se_qr2[1],
           pqrc2[2],se_qr2[2],
           pqrc2[3],se_qr2[3],
           0,0,
           mlqc2[1],sqrt(mlqr_var2)[1],
           mlqc2[2],sqrt(mlqr_var2)[2],
           mlqc2[3],sqrt(mlqr_var2)[3],
           0,0,
           pqrc3[1],se_qr3[1],
           pqrc3[2],se_qr3[2],
           pqrc3[3],se_qr3[3],
           0,0,
           mlqc3[1],sqrt(mlqr_var3)[1],
           mlqc3[2],sqrt(mlqr_var3)[2],
           mlqc3[3],sqrt(mlqr_var3)[3],
           0,0,
           pqrc4[1],se_qr4[1],
           pqrc4[2],se_qr4[2],
           pqrc4[3],se_qr4[3],
           0,0,
           mlqc4[1],sqrt(mlqr_var4)[1],
           mlqc4[2],sqrt(mlqr_var4)[2],
           mlqc4[3],sqrt(mlqr_var4)[3],
           0,0,
           pqrc5[1],se_qr5[1],
           pqrc5[2],se_qr5[2],
           pqrc5[3],se_qr5[3],
           0,0,
           mlqc5[1],sqrt(mlqr_var5)[1],
           mlqc5[2],sqrt(mlqr_var5)[2],
           mlqc5[3],sqrt(mlqr_var5)[3],
           0,0))
}

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

expit <- function(x){return(exp(x)/(1+exp(x)))}