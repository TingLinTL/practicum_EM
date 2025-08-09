#combine test
library(survival)
set.seed(2000)
n <- 500 #sample size 
tau <- 5.5 #maximum follow-up time
n_sim <- 9 #number of simulation
t_pred <- 2  #time point for estimates 
true_spce_aft <- 0.1901985 #true SPCE at time=2 from numeric iteration for weibull aft model
n_boot <-100

spce_estimates_em_aft <- spce_estimates_em_cox <- spce_estimates_Sto_aft <- spce_estimates_Sto_cox <- numeric(n_sim)
sd_em_aft <- sd_em_cox <- sd_Sto_aft <- sd_Sto_cox <- numeric(n_sim)


emU_spce_cox <- function(t, d, z, x, zetat, zetaz, theta, t0, iter = 20) {
  x <- data.matrix(x)
  n <- length(t)
  nx <- ncol(x)
  
  fn <- function(beta, x, z, p, zetaz) {
    -sum(z * (log(pnorm(x %*% beta + zetaz))*p + log(pnorm(x %*% beta))*(1-p)) +
           (1-z) * (log(1 - pnorm(x %*% beta + zetaz))*p + log(1 - pnorm(x %*% beta))*(1 - p)))
  }
  
  # === Initialization ===
  U.fit1 <- coxph(Surv(t, d) ~ z + x)
  t.coef1 <- c(U.fit1$coef, zetat)
  z.coef <- c(glm(z ~ x, family = binomial(link = "probit"))$coef, zetaz)
  
  p <- rep(0, n)
  
  # === EM Iteration ===
  for (j in 1:iter) {
    bh1 <- basehaz(U.fit1, centered = FALSE)
    index1 <- match(t, bh1$time)
    
    ptzu1 <- (1 - pnorm(cbind(1, x, 1) %*% z.coef))^(1 - z) * pnorm(cbind(1, x, 1) %*% z.coef)^z * theta *
      exp(cbind(z, x, 1) %*% t.coef1)^d *
      exp(-bh1[index1, 1] / exp(mean(log(exp(zetat) * p + (1 - p)))) * 
            exp(cbind(z, x, 1) %*% t.coef1))
    
    ptzu0 <- (1 - pnorm(cbind(1, x, 0) %*% z.coef))^(1 - z) * pnorm(cbind(1, x, 0) %*% z.coef)^z * (1 - theta) *
      exp(cbind(z, x, 0) %*% t.coef1)^d *
      exp(-bh1[index1, 1] / exp(mean(log(exp(zetat) * p + (1 - p)))) * 
            exp(cbind(z, x, 0) %*% t.coef1))
    
    p <- ptzu1 / (ptzu1 + ptzu0)
    p[is.na(p)] <- 0
    
    # Update model fits
    U.fit1 <- coxph(Surv(t, d) ~ z + x + offset(log(exp(zetat) * p + (1 - p))))
    t.coef1 <- c(U.fit1$coef, zetat)
    
    z.fit <- optim(par = z.coef[1:(nx + 1)], fn, x = cbind(1, x), z = z, p = p, zetaz = zetaz)
    z.coef <- c(z.fit$par, zetaz)
  }
  
  # === Compute SPCE ===
  linpred <- cbind(1, x, p) %*% matrix(z.coef, ncol = 1)
  ps <- pnorm(linpred)
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  pZ <- mean(z)
  w <- (pZ * z / ps) + ((1 - pZ) * (1 - z) / (1 - ps))
  w <- pmin(pmax(w, 0.1), 10)
  
  fit_w <- survfit(Surv(t, d) ~ z, weights = w)
  surv_summary <- summary(fit_w, times = t0)
  S_Z0 <- surv_summary$surv[surv_summary$strata == "z=0"]
  S_Z1 <- surv_summary$surv[surv_summary$strata == "z=1"]
  SPCE <- S_Z1 - S_Z0
  
  return(list(
    SPCE = SPCE,
    S_Z1 = S_Z1,
    S_Z0 = S_Z0,
    p = p,
    z.coef = z.coef,
    t.coef = t.coef1,
    weights = w,
    basehaz1 = basehaz(U.fit1, centered = FALSE)
  ))
}

emU_spce_aft <- function(t, d, z, x, gammat, gammaz, theta, t0, iter = 20) {
  fn_outcome <- function(params, t, d, x, z, p, gammat) {
    beta      <- params[1:(ncol(x)+1)]
    alpha     <- params[ncol(x)+2]
    log_sigma <- params[ncol(x)+3]
    sigma     <- exp(log_sigma)
    
    log_t <- log(pmax(t, 1e-8))
    xb <- cbind(1, x) %*% beta
    lp_meanU <- xb + alpha*z + gammat*p
    
    term1 <- d * (-log(sigma) + ((1/sigma - 1)*log_t) - (1/sigma)*lp_meanU)
    
    eta0 <- (log_t - xb - alpha*z - gammat*0)/sigma
    eta1 <- (log_t - xb - alpha*z - gammat*1)/sigma
    eta0 <- pmin(eta0, 50)
    eta1 <- pmin(eta1, 50)
    
    exp0 <- exp(eta0)
    exp1 <- exp(eta1)
    term2 <- p * exp1 + (1 - p) * exp0
    term2 <- pmax(term2, 1e-12)
    
    nll <- -sum(term1 - term2)
    if (!is.finite(nll)) return(1e12)
    return(nll)
  }
  
  fn_treat <- function(beta, x, z, p, gammaz) {
    lp1 <- cbind(1, x) %*% beta + gammaz
    lp0 <- cbind(1, x) %*% beta
    term1 <- z * (log(plogis(lp1)) * p + log(plogis(lp0)) * (1 - p))
    term2 <- (1 - z) * (log(1 - plogis(lp1)) * p + log(1 - plogis(lp0)) * (1 - p))
    return(-sum(term1 + term2))
  }
  
  x <- data.matrix(x)
  n <- length(t)
  nx <- ncol(x)
  
  # === Initialization ===
  U.fit1 <- survreg(Surv(t, d) ~ x + z, dist = "weibull")
  sigma <- U.fit1$scale
  t.coef1 <- c(coef(U.fit1), gammat)
  init_params <- c(coef(U.fit1), log(sigma))
  
  Z.fit1 <- glm(z ~ x, family = binomial(link = "logit"))
  z.coef <- c(coef(Z.fit1), gammaz)
  
  p <- rep(0.5, n)
  
  for (j in 1:iter) {
    ## E-step
    mu_u1 <- cbind(1, x, z, 1) %*% matrix(t.coef1, ncol = 1)
    mu_u0 <- cbind(1, x, z, 0) %*% matrix(t.coef1, ncol = 1)
    alpha_weibull <- 1 / sigma
    lambda_u1 <- exp(-mu_u1)
    lambda_u0 <- exp(-mu_u0)
    
    lik_u1 <- (alpha_weibull * lambda_u1^alpha_weibull * t^(alpha_weibull - 1))^d *
      exp(-(lambda_u1 * t)^alpha_weibull)
    lik_u0 <- (alpha_weibull * lambda_u0^alpha_weibull * t^(alpha_weibull - 1))^d *
      exp(-(lambda_u0 * t)^alpha_weibull)
    
    logit_zu1 <- (1 - plogis(cbind(1, x, 1) %*% z.coef))^(1 - z) *
      plogis(cbind(1, x, 1) %*% z.coef)^z
    logit_zu0 <- (1 - plogis(cbind(1, x, 0) %*% z.coef))^(1 - z) *
      plogis(cbind(1, x, 0) %*% z.coef)^z
    
    ptzu1 <- lik_u1 * logit_zu1 * theta
    ptzu0 <- lik_u0 * logit_zu0 * (1 - theta)
    
    p <- ptzu1 / (ptzu1 + ptzu0)
    p[ptzu1 == 0 & ptzu0 == 0] <- 0.5
    
    ## M-step: outcome model
    outcome_fit <- optim(par = init_params,
                         fn  = fn_outcome,
                         t   = t, d = d, x = x, z = z, p = p, gammat = gammat,
                         method = "BFGS", control = list(maxit = 1000))
    
    params  <- outcome_fit$par
    t.coef1 <- c(params[1:(ncol(x)+2)], gammat)
    sigma   <- exp(params[ncol(x)+3])
    init_params <- c(params[1:(ncol(x)+2)], log(sigma))
    
    ## M-step: treatment model
    z.fit <- optim(par = z.coef[1:(nx+1)],
                   fn  = fn_treat,
                   x   = x, z = z, p = p, gammaz = gammaz, method = "BFGS")
    z.coef <- c(z.fit$par, gammaz)
  }
  
  # === Compute SPCE ===
  linpred <- cbind(1, x, p) %*% matrix(z.coef, ncol = 1)
  ps <- plogis(linpred)
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  pZ <- mean(z)
  w <- (pZ * z / ps) + ((1 - pZ) * (1 - z) / (1 - ps))
  w <- pmin(pmax(w, 0.1), 10)
  
  fit_w <- survfit(Surv(t, d) ~ z, weights = w)
  surv_summary <- summary(fit_w, times = t0)
  
  S_Z0 <- surv_summary$surv[surv_summary$strata == "z=0"]
  S_Z1 <- surv_summary$surv[surv_summary$strata == "z=1"]
  SPCE <- S_Z1 - S_Z0
  
  return(list(
    SPCE     = SPCE,
    S_Z1     = S_Z1,
    S_Z0     = S_Z0,
    weights  = w,
    pU       = p,
    z.coef   = z.coef,
    t.coef   = t.coef1,
    sigma    = sigma
  ))
}

surv_stoEM_ipw_aft <- function(t, d, Z, X, zetat, zetaz, B, theta = 0.5){
  
  nx = 2
  n = length(t)
  
  
  #Record coefficients with simulated U
  spce = numeric(B)
  
  for (j in 1:B){
    Usim = SimulateU_surv_aft(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta)
    
    Z.fit = glm(Z ~ X, offset = zetaz * Usim$U, family=binomial(link="probit"))
    
    ps = Z.fit$fitted.values
    ipw = (sum(Z)/n) * Z/ps + (1-(sum(Z)/n))*(1-Z)/(1-ps)
    ipw = pmin(ipw, 10)
    ipw = pmax(ipw, 0.1)
    
    t1.ipw<- survreg(Surv(t,d) ~ Z, dist = "weibull", weights = ipw, robust = TRUE)
    
    scale <- t1.ipw$scale
    shape <- 1 / scale
    coef <- t1.ipw$coefficients
    
    new_Z1 <- data.frame(Z = rep(1, n))
    new_Z0 <- data.frame(Z = rep(0, n))
    # Linear predictors
    lp1 <- model.matrix(~ Z, data = new_Z1) %*% coef
    lp0 <- model.matrix(~ Z, data = new_Z0) %*% coef
    
    # Linear scale (i.e., time scale)
    scale1 <- exp(lp1)
    scale0 <- exp(lp0)
    
    # Predict survival probability at t = 2
    t0 <- 2
    S1 <- exp(- (t0 / scale1)^shape)
    S0 <- exp(- (t0 / scale0)^shape)
    
    # Compute SPCE
    spce[j] <- mean(S1) - mean(S0)
    
  }
  
  spce = mean(spce)
  
  
  return(list(spce = spce))
}




surv_stoEM_ipw_cox <- function(t, d, Z, X, zetat, zetaz, B, theta = 0.5){
  
  nx = 2
  n = length(t)
  
  
  #Record coefficients with simulated U
  spce = numeric(B)
  
  for (j in 1:B){
    Usim = SimulateU_surv_cox(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta)
    
    Z.fit = glm(Z ~ X, offset = zetaz * Usim$U, family=binomial(link="probit"))
    
    ps = Z.fit$fitted.values
    ipw = (sum(Z)/n) * Z/ps + (1-(sum(Z)/n))*(1-Z)/(1-ps)
    ipw = pmin(ipw, 10)
    ipw = pmax(ipw, 0.1)
    
    t1.ipw = coxph(Surv(t, d) ~ Z, weights = ipw, robust = TRUE)
    
    spce[j] = mean(predict(t1.ipw, newdata=data.frame(Z=c(rep(1,n))),type = "survival", times = 2))-
      mean(predict(t1.ipw, newdata=data.frame(Z=c(rep(0,n))),type = "survival", times = 2))
    
  }
  
  spce = mean(spce)
  
  
  return(list(spce = spce))
}

SimulateU_surv_cox <- function(t, d, z, x, zetat, zetaz, theta, iter = 20, weights = NULL){
  #t is a vector of n, time to event
  #d is a vector of n, indicator of event
  #z is a vector of n, treatment
  #x is a matrix, covariates
  #zetat is a scaler, zetaz is a scaler, sensitivity parameters
  
  n = length(t) #number of observations
  x = data.matrix(x)
  nx = dim(x)[2] #number of observed covariates
  
  p = theta
  U = rbinom(n,1,p)
  #Upath = U (useful for checking convergence)
  
  for(j in 1:iter){
    
    #fit time to event
    U.fit1 = coxph(Surv(t,d) ~ z + x + U, weights=weights)
    t.coef1 = U.fit1$coef
    t.coef1[length(t.coef1)]  = zetat
    #fit treatment
    z.coef = glm(z ~ x + U, family = binomial(link="probit"), control = glm.control(epsilon = 1e-6, maxit = 50))$coef
    z.coef[length(z.coef)] = zetaz
    
    
    t.coef1[is.na(t.coef1)] = 0
    z.coef[is.na(z.coef)] = 0
    
    bh1 = basehaz(U.fit1, centered=F) #cumulative baseline hazard for remission
    index1 = match(t,bh1$time)
    
    
    ptzu1 = (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta*
      exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1))^d * exp(-bh1[index1,1] * exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1)))
    
    ptzu0 = (1-pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1))^z*(1-theta)*
      exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1))^d * exp(-bh1[index1,1] * exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1)))
    
    
    p = ptzu1/(ptzu1 + ptzu0)
    p[ptzu1==0 & ptzu0==0] = 0
    
    U = rbinom(n,1,p)
    #Upath = cbind(Upath, U)
  }
  
  return(list(U = U, p = p))
}

SimulateU_surv_aft <- function(t, d, z, x, zetat, zetaz, theta, iter = 20){
  #t is a vector of n, time to event
  #d is a vector of n, indicator of event
  #z is a vector of n, treatment
  #x is a matrix, covariates
  #zetat is a scaler, zetaz is a scaler, sensitivity parameters
  
  n = length(t) #number of observations
  x = data.matrix(x)
  nx = dim(x)[2] #number of observed covariates
  
  p = theta
  U = rbinom(n,1,p)
  #Upath = U (useful for checking convergence)
  
  for(j in 1:iter){
    #fit time to event
    U.fit1 = survreg(Surv(t,d) ~ x + z + U, dist = "weibull")
    t.coef1 = U.fit1$coef
    t.coef1[length(t.coef1)]  = zetat
    sigma <- U.fit1$scale
    
    #fit treatment
    z.coef = glm(z ~ x + U, family = binomial(link = "logit"), control = glm.control(epsilon = 1e-6, maxit = 50))$coef
    z.coef[length(z.coef)] = zetaz
    
    
    t.coef1[is.na(t.coef1)] = 0
    z.coef[is.na(z.coef)] = 0
    
    # Logistic treatment model linear predictors
    logit_zu1 <- (1 - plogis(cbind(1, x, 1) %*% matrix(z.coef, ncol = 1)))^(1 - z) *
      plogis(cbind(1, x, 1) %*% matrix(z.coef, ncol = 1))^z
    
    logit_zu0 <- (1 - plogis(cbind(1, x, 0) %*% matrix(z.coef, ncol = 1)))^(1 - z) *
      plogis(cbind(1, x, 0) %*% matrix(z.coef, ncol = 1))^z
    
    # weibull aft linear prediction under U = 1 and U = 0
    mu_u1 <- cbind(1, x, z, 1) %*% matrix(t.coef1, ncol = 1)
    mu_u0 <- cbind(1, x, z, 0) %*% matrix(t.coef1, ncol = 1)
    alpha_weibull <- 1 / sigma
    
    # Lambda (rate-like scale)
    lambda_u1 <- exp(-mu_u1)
    lambda_u0 <- exp(-mu_u0)
    
    #weibull density and survival
    lik_u1 <- ( alpha_weibull * lambda_u1^alpha_weibull * t^(alpha_weibull - 1) )^d *
      exp( - (lambda_u1 * t)^alpha_weibull )
    
    lik_u0 <- ( alpha_weibull * lambda_u0^alpha_weibull * t^(alpha_weibull - 1) )^d *
      exp( - (lambda_u0 * t)^alpha_weibull )
    
    
    ptzu1 = logit_zu1*theta* lik_u1
    ptzu0 = logit_zu0*(1-theta)* lik_u0
    
    p = ptzu1/(ptzu1 + ptzu0)
    p[ptzu1==0 & ptzu0==0] = 0
    
    U = rbinom(n,1,p)
    #Upath = cbind(Upath, U)
  }
  
  return(list(U = U, p = p))
}


#================ simulation =============================#


for (i in 1:n_sim) {
  cat("Running simulation", i, "of", n_sim, "\n")
  
  X1 <- rbinom(n, size = 1, prob = 0.4)
  X2 <- rbinom(n, size = 1, prob = 0.6)
  U1 <- rbinom(n, size = 1, prob = 0.5)
  alpha_0  <- 0.1; alpha_x1 <- 0.3; alpha_x2 <- 0.5; alpha_u  <- -0.8
  linpred_A <- alpha_0 + alpha_x1 * X1 + alpha_x2 * X2 + alpha_u * U1
  pi_A      <- pnorm(linpred_A)
  A         <- rbinom(n, size = 1, prob = pi_A)
  eta_intercept_a0 <- 0.7; eta_intercept_a1 <- 0.2; eta_x1 <- -0.1; eta_x2 <- 0.4; eta_u  <- -0.8
  random_simu <- runif(n) #1-F(t) ~ Unif(0,1)
  v = 2 #weibull with shape = v
  lambda = 0.1 #base line rate, T~Weibull(v, rate) with rate for each individual is lambda*e^(X*beta)
  #T_i = (-log(U)/ lambda_i)^(1/v), lambda_i = lambda * e^(X*beta)
  #survival time D=(-log(random_simu) / (lambda * exp(eta_intercept + eta_x1 * X1 + eta_x2 * X2 + eta_u * U))) ^ (1 / v)
  lp_a1   <- eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U1
  lp_a0   <- eta_intercept_a0 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U1
  D_a0 = (-log(random_simu) / (lambda * exp(lp_a0))) ^ (1 / v)
  D_a1 = (-log(random_simu) / (lambda * exp(lp_a1))) ^ (1 / v)
  
  # sigma0 <- sigma1 <- 1/2; alpha0 <- alpha1 <- 1/sigma0
  # lp_a1   <- eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U1
  # lp_a0   <- eta_intercept_a0 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U1
  # lambda1 <- exp(-lp_a1); lambda0 <- exp(-lp_a0)
  # D_a0 <- rweibull(n, shape = alpha0, scale = 1/lambda0)
  # D_a1 <- rweibull(n, shape = alpha1, scale = 1/lambda1)
  
  # ---- Non-informative censoring ----
  C_a0 <- runif(n, 0.1, tau)
  C_a1 <- runif(n, 0.1, tau)
  
  # ---- Observed time & event indicator ----
  M <- ifelse(A == 1, pmin(tau, C_a1, D_a1),
              pmin(tau, C_a0, D_a0))
  delta <- ifelse(A == 1,
                  as.numeric(D_a1 <= pmin(tau, C_a1)),
                  as.numeric(D_a0 <= pmin(tau, C_a0)))
  
  ## ---- Observed dataset ----
  data_sim <- data.frame(
    M     = M,
    delta = delta,
    A     = A,
    X1    = X1,
    X2    = X2,
    U1    = U1
  )
  
  x_mat <- as.matrix(data_sim[, c("X1", "X2")]) # Prepare the covariate matrix
  
  EM_SPCE_aft <- emU_spce_aft(t = data_sim$M,
                              d = data_sim$delta,
                              z = data_sim$A, 
                              x = x_mat, gammat =-0.8, gammaz =-0.8, theta=0.5, t0=t_pred)
  
  
  EM_SPCE_cox <- emU_spce_cox(t = data_sim$M,
                              d = data_sim$delta,
                              z = data_sim$A, 
                              x = x_mat, zetat=-0.8, zetaz=-0.8, theta=0.5, t0=t_pred)
  
  
  stoEM_SPCE_aft <- surv_stoEM_ipw_aft(t= data_sim$M,
                                       d= data_sim$delta,
                                       Z= data_sim$A,
                                       X= x_mat, zetat=-0.8, zetaz=-0.8, B=10, theta = 0.5)
  
  
  stoEM_SPCE_cox <- surv_stoEM_ipw_cox(t= data_sim$M,
                                       d= data_sim$delta,
                                       Z= data_sim$A,
                                       X= x_mat, zetat=-0.8, zetaz=-0.8, B=10, theta = 0.5)
  
  spce_estimates_em_aft[i] <- EM_SPCE_aft$SPCE
  spce_estimates_em_cox[i] <- EM_SPCE_cox$SPCE
  spce_estimates_Sto_aft[i] <- stoEM_SPCE_aft$spce
  spce_estimates_Sto_cox[i] <- stoEM_SPCE_cox$spce
  
  ## --------------------------------------- ##
  ## === BOOTSTRAP inside the iteration  === ##
  ## --------------------------------------- ##
  B_boot <- n_boot
  boot_spce_em_aft  <- numeric(B_boot)
  boot_spce_em_cox  <- numeric(B_boot)
  boot_spce_Sto_aft <- numeric(B_boot)
  boot_spce_Sto_cox <- numeric(B_boot)
  
  for (b in 1:B_boot) {
    boot_idx <- sample(1:n, size = n, replace = TRUE)
    data_boot <- data_sim[boot_idx, ]
    x_boot <- as.matrix(data_boot[, c("X1", "X2")])
    
    
    ## ---------- EM AFT ----------
    boot_spce_em_aft[b] <- emU_spce_aft(t= data_boot$M,
                                        d = data_boot$delta,
                                        z = data_boot$A, 
                                        x = x_boot, gammat=-0.8, gammaz=-0.8, theta=0.5, t0=t_pred)$SPCE
    
    ## ---------- EM Cox ----------
    boot_spce_em_cox[b] <- emU_spce_cox(t= data_boot$M,
                                        d = data_boot$delta,
                                        z = data_boot$A, 
                                        x = x_boot, zetat=-0.8, zetaz=-0.8, theta=0.5, t0=t_pred)$SPCE
    
    ## ---------- StoEM AFT ----------
    boot_spce_Sto_aft[b] <- surv_stoEM_ipw_aft(t= data_boot$M,
                                               d = data_boot$delta,
                                               Z = data_boot$A, X = x_boot, zetat=-0.8, zetaz=-0.8, B=10, theta = 0.5)$spce
    
    
    ## ---------- StoEM Cox ----------
    boot_spce_Sto_cox[b] <- surv_stoEM_ipw_cox(t= data_boot$M,
                                               d = data_boot$delta,
                                               Z = data_boot$A, X = x_boot, zetat=-0.8, zetaz=-0.8, B=10, theta = 0.5)$spce
    
  }
  
  sd_em_aft[i] <- sd(boot_spce_em_aft)
  sd_em_cox[i] <- sd(boot_spce_em_cox)
  sd_Sto_aft[i] <- sd(boot_spce_Sto_aft)
  sd_Sto_cox[i] <- sd(boot_spce_Sto_cox)
  
}

spce_estimates_em_aft;spce_estimates_em_cox;spce_estimates_Sto_aft;spce_estimates_Sto_cox
sd_em_aft;sd_em_cox;sd_Sto_aft;sd_Sto_cox
