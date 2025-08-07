library(survival)
library(parallel)
compiler::enableJIT(3)

# ================== functions  ================== #
emU_surv_aft <- function(t, d, z, x, gammat, gammaz, theta, iter = 20){
  fn_outcome <- function(params, t, d, x, z, p, gammat) {
    # Unpack params
    beta      <- params[1:(ncol(x)+1)]
    alpha     <- params[ncol(x)+2]
    log_sigma <- params[ncol(x)+3]
    sigma     <- exp(log_sigma)  # always >0
    
    # Safe log(t)
    log_t <- log(pmax(t, 1e-8))
    
    # Linear predictor
    xb <- cbind(1, x) %*% beta
    lp_meanU <- xb + alpha*z + gammat*p
    
    # First term
    term1 <- d * (-log(sigma) + ((1/sigma - 1)*log_t) - (1/sigma)*lp_meanU)
    
    # Second term: prevent exp overflow
    eta0 <- (log_t - xb - alpha*z - gammat*0)/sigma
    eta1 <- (log_t - xb - alpha*z - gammat*1)/sigma
    eta0 <- pmin(eta0, 50)
    eta1 <- pmin(eta1, 50)
    
    exp0 <- exp(eta0)
    exp1 <- exp(eta1)
    term2 <- p * exp1 + (1 - p) * exp0
    term2 <- pmax(term2, 1e-12)  # avoid log(0)
    
    nll <- -sum(term1 - term2)
    
    # If nll is non-finite, return huge penalty
    if (!is.finite(nll)) return(1e12)
    
    return(nll)
  }
  
  
  fn_treat <- function(beta, x, z, p, gammaz) {
    # linear predictors
    lp1 <- cbind(1, x) %*% beta + gammaz  # when U_i = 1
    lp0 <- cbind(1, x) %*% beta           # when U_i = 0
    
    term1 <- z * (log(1 / (1 + exp(-lp1))) * p + log(1 / (1 + exp(-lp0))) * (1 - p))
    term2 <- (1 - z) * (log(1 - 1 / (1 + exp(-lp1))) * p + log(1 - 1 / (1 + exp(-lp0))) * (1 - p))
    
    #negative tretment log-likelihood
    return(-sum(term1 + term2))
  }
  
  
  x = data.matrix(x)
  n = length(t) #number of observations
  nx = dim(x)[2] #number of covariates ,or nx = ncol(x)
  
  #initialize parameters
  #fit time to event
  U.fit1 <- survreg(Surv(t, d) ~ x + z, dist = "weibull")
  # survreg fits AFT model: log(T) = Xβ + error
  # extract initial values for next iteration
  t.coef1 = c(coef(U.fit1), gammat) #beta0...beta_nx, alpha,gammat
  sigma <- U.fit1$scale
  init_params  <- c(coef(U.fit1), log(sigma)) 
  
  #fit treatment
  Z.fit1 <- glm(z ~ x, family = binomial(link = "logit"))
  z.coef <- c(coef(Z.fit1), gammaz) #beta_z1,beta_z2...beta_znx,gammaz
  
  t.coef1[is.na(t.coef1)] = 0
  z.coef[is.na(z.coef)] = 0
  
  p = rep(0,n)
  coef_trace <- matrix(NA, nrow = iter, ncol = length(init_params))  # track beta, alpha, log_sigma
  colnames(coef_trace) <- c(
    paste0("beta", 0:ncol(x)),   # intercept & covariates
    "alpha",
    "log_sigma"
  )
  
  for (j in 1:iter) {
    
    ## === E-step === ##
    mu_u1 <- cbind(1, x, z, 1) %*% matrix(t.coef1, ncol = 1)
    mu_u0 <- cbind(1, x, z, 0) %*% matrix(t.coef1, ncol = 1)
    alpha_weibull <- 1 / sigma
    
    lambda_u1 <- exp(-mu_u1)
    lambda_u0 <- exp(-mu_u0)
    
    lik_u1 <- (alpha_weibull * lambda_u1^alpha_weibull * t^(alpha_weibull - 1))^d *
      exp(-(lambda_u1 * t)^alpha_weibull)
    lik_u0 <- (alpha_weibull * lambda_u0^alpha_weibull * t^(alpha_weibull - 1))^d *
      exp(-(lambda_u0 * t)^alpha_weibull)
    
    logit_zu1 <- (1 - plogis(cbind(1, x, 1) %*% matrix(z.coef, ncol=1)))^(1 - z) *
      plogis(cbind(1, x, 1) %*% matrix(z.coef, ncol=1))^z
    logit_zu0 <- (1 - plogis(cbind(1, x, 0) %*% matrix(z.coef, ncol=1)))^(1 - z) *
      plogis(cbind(1, x, 0) %*% matrix(z.coef, ncol=1))^z
    
    ptzu1 <- lik_u1 * logit_zu1 * theta
    ptzu0 <- lik_u0 * logit_zu0 * (1 - theta)
    
    p <- ptzu1 / (ptzu1 + ptzu0)
    p[ptzu1 == 0 & ptzu0 == 0] <- 0
    
    ## === M-step: time-to-event === ##
    # use last params as init
    outcome_fit <- optim(par = init_params,
                         fn  = fn_outcome,
                         t   = t, d = d, x = x, z = z, p = p, gammat = gammat,
                         method = "BFGS",
                         control = list(maxit = 1000))
    
    params  <- outcome_fit$par
    t.coef1 <- c(params[1:(ncol(x)+2)], gammat)  # keep same design columns
    sigma   <- exp(params[ncol(x)+3])            # back-transform log_sigma
    init_params <- c(params[1:(ncol(x)+2)], log(sigma)) 
    
    # store the *updated* params for plotting
    coef_trace[j, ] <- init_params
    
    ## === M-step: treatment model === ##
    z.fit <- optim(par = z.coef[1:(nx+1)],
                   fn  = fn_treat,
                   x   = x, z = z, p = p, gammaz = gammaz,method = "BFGS")
    z.coef <- c(z.fit$par, gammaz)
    
    t.coef1[is.na(t.coef1)] <- 0
    z.coef[is.na(z.coef)]   <- 0
  }
  
  return(list(p = p, t.coef1 = t.coef1, z.coef = z.coef, sigma = sigma, trace = coef_trace))
  #all coefficients
}


emU_surv_cox <- function(t, d, z, x, zetat, zetaz, theta, iter = 20){
  #t is a vector of n, time to event
  #d is a vector of n, indicator of event
  #z is a vector of n, treatment
  #x is a matrix, covariates
  #zetat is a scaler, zetaz is a scaler, sensitivity parameters
  
  fn <- function(beta, x, z, p, zetaz) {
    -sum(z * (log(pnorm(x %*% beta + zetaz))*p + log(pnorm(x %*% beta))*(1-p)) +
           (1-z) * (log(1-pnorm(x %*% beta + zetaz))*p + log(1-pnorm(x %*% beta))*(1-p)))
  }
  
  x = data.matrix(x)
  n = length(t) #number of observations
  nx = dim(x)[2] #number of covariates
  
  #initialize parameters
  #fit time to event
  U.fit1 = coxph(Surv(t,d) ~ z + x)
  t.coef1 = c(U.fit1$coef, zetat)
  #fit treatment
  z.coef = c(glm(z ~ x, family=binomial(link="probit"), control = glm.control(epsilon = 1e-6, maxit = 50))$coef, zetaz)
  
  t.coef1[is.na(t.coef1)] = 0
  z.coef[is.na(z.coef)] = 0
  
  p = rep(0,n)
  
  for(j in 1:iter){
    bh1 = basehaz(U.fit1, centered=F) #cumulative baseline hazard for event with mean(offset)
    index1 = match(t,bh1$time)
    
    ptzu1 = (1-pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,1)%*%matrix(z.coef, ncol = 1))^z*theta*
      exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1))^d * exp(-bh1[index1,1]/exp(mean(log(exp(zetat)*p+(1-p)))) * exp(cbind(z,x,1)%*%matrix(t.coef1, ncol = 1)))
    
    ptzu0 = (1-pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1)))^(1-z)*
      pnorm(cbind(1,x,0)%*%matrix(z.coef, ncol = 1))^z*(1-theta)*
      exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1))^d * exp(-bh1[index1,1]/exp(mean(log(exp(zetat)*p+(1-p)))) * exp(cbind(z,x,0)%*%matrix(t.coef1, ncol = 1)))
    
    p = ptzu1/(ptzu1 + ptzu0) #posterior dist of U
    p[ptzu1==0 & ptzu0==0] = 0
    
    #fit time to event
    U.fit1 = coxph(Surv(t,d) ~ z + x + offset(log(exp(zetat)*p+(1-p))))
    t.coef1 = c(U.fit1$coef, zetat[1])
    #fit treatment
    z.fit = optim(par = z.coef[1:(nx+1)], fn, x=cbind(1,x), z=z, p=p, zetaz=zetaz)
    z.coef = c(z.fit$par, zetaz)
    
    t.coef1[is.na(t.coef1)] = 0
    z.coef[is.na(z.coef)] = 0
  }
  
  return(list(
    p        = p,
    t.coef1  = t.coef1,
    z.coef   = z.coef,
    U.fit1   = U.fit1,
    basehaz1 = basehaz(U.fit1, centered = FALSE)
  ))
}

compute_SPCE_EM_aft <- function(res_em, X, Z, M, delta, t0) {
  # unpack converged treatment model parameters
  z.coef_final <- res_em$z.coef    
  U_final <- res_em$p             
  
  # fitted propensity scores P(Z=1|X,U)
  linpred <- cbind(1, X, U_final) %*% matrix(z.coef_final, ncol = 1)
  ps <- plogis(linpred)
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  # stabilized weights
  pZ <- mean(Z) # marginal treatment probability
  w <- (pZ * Z / ps) + ((1 - pZ) * (1 - Z) / (1 - ps))
  
  # truncate extreme weights
  w <- pmin(w, 10)
  w <- pmax(w, 0.1)
  
  fit_w <- survfit(Surv(M, delta) ~ Z, weights = w)
  surv_summary <- summary(fit_w, times = t0)
  # Extract survival for Z=0 and Z=1
  S_Z0 <- surv_summary$surv[surv_summary$strata == "Z=0"]
  S_Z1 <- surv_summary$surv[surv_summary$strata == "Z=1"]
  
  SPCE <- S_Z1 - S_Z0
  
  return(list(
    SPCE = SPCE,
    S_Z1 = S_Z1,
    S_Z0 = S_Z0,
    weights = w
  ))
}


compute_SPCE_EM_cox <- function(res_em, X, Z, M, delta, t0) {
  
  z.coef_final <- res_em$z.coef   # intercept, X, and zetaz
  pU           <- res_em$p        # posterior mean of U
  
  # compute linear predictor for P(Z=1|X,U) using probit
  linpred <- cbind(1, X, pU) %*% matrix(z.coef_final, ncol = 1)
  ps <- pnorm(linpred)  # probit link for treatment model
  
  # compute stabilized weights
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)  # avoid extreme probabilities
  pZ <- mean(Z)  # marginal treatment probability
  w <- (pZ * Z / ps) + ((1 - pZ) * (1 - Z) / (1 - ps))
  
  # truncate extreme weights
  w <- pmin(w, 10)
  w <- pmax(w, 0.1)
  
  # fit weighted Cox model
  fit_w <- survfit(Surv(M, delta) ~ Z, weights = w)
  
  # survival probabilities at time t0
  surv_summary <- summary(fit_w, times = t0)
  S_Z0 <- surv_summary$surv[surv_summary$strata == "Z=0"]
  S_Z1 <- surv_summary$surv[surv_summary$strata == "Z=1"]
  
  # compute SPCE
  SPCE <- S_Z1 - S_Z0
  
  return(list(
    SPCE   = SPCE,
    S_Z1   = S_Z1,
    S_Z0   = S_Z0,
    weights = w
  ))
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


SimulateU_surv_cox <- function(t, d, z, x, zetat, zetaz, theta, iter = 20){
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
    U.fit1 = coxph(Surv(t,d) ~ z + x + U)
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

surv_stoEM_aft <- function(t, d, Z, X, zetat, zetaz, B, theta) {
  
  n <- length(t)
  
  final_beta <- NULL
  final_sigma <- NULL
  
  for (j in 1:B){
    Usim = SimulateU_surv_aft(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta) 
    
    
    t1.fit<- survreg(Surv(t,d) ~ X + Z + Usim$U, dist = "weibull")
    
    # extract parameters
    final_beta <- t1.fit$coef       # β0, β_X, β_Z, β_U
    final_sigma <- t1.fit$scale      # Weibull scale parameter
    #alpha <- 1 / sigma_hat        # Weibull shape 
    
  }
  # return only the *final* coefficients after B iterations
  return(list(
    beta = final_beta,    # includes β0, β_X, β_Z, β_U
    sigma = final_sigma   # Weibull scale
  ))
}

surv_stoEM_cox <- function(t, d, Z, X, zetat, zetaz, B, theta) {
  nx <- ncol(X)
  n <- length(t)
  
  final_beta <- NULL
  final_fit  <- NULL
  
  for (j in 1:B) {
    Usim <- SimulateU_surv_cox(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta)
    
    t1.fit <- coxph(Surv(t, d) ~ X + Z + Usim$U)
    
    final_beta <- t1.fit$coef
    final_fit  <- t1.fit  
  }
  
  return(list(
    beta    = final_beta,
    fit     = final_fit,
    basehaz = basehaz(final_fit, centered = FALSE)
  ))
}

surv_stoEM_ipw_aft  <- function(beta_final, sigma_final,
                                t, d, Z, X,
                                t0,
                                stabilize) {
  alpha <- 1 / sigma_final
  n <- nrow(X)
  
  # --- Estimate propensity scores ---
  ps_model <- glm(Z ~ ., data = data.frame(Z = Z, X), family = binomial)
  ps_hat <- predict(ps_model, type = "response")
  
  # --- Compute IP weights ---
  if (stabilize) {
    # stabilized weights
    pZ <- mean(Z)
    w <- ifelse(Z == 1,
                pZ / ps_hat,
                (1 - pZ) / (1 - ps_hat))
  } else {
    # unstabilized weights
    w <- ifelse(Z == 1, 1 / ps_hat, 1 / (1 - ps_hat))
  }
  
  # --- Compute survival under Z=1 & Z=0 ---
  # we marginalize U as 0.5 (prior)
  mu_Z1 <- cbind(1, X, 1, 0.5) %*% beta_final
  mu_Z0 <- cbind(1, X, 0, 0.5) %*% beta_final
  
  lambda_Z1 <- exp(-mu_Z1)
  lambda_Z0 <- exp(-mu_Z0)
  
  S_Z1_i <- exp(-(lambda_Z1 * t0)^alpha)
  S_Z0_i <- exp(-(lambda_Z0 * t0)^alpha)
  
  # --- weighted marginal survival ---
  S_Z1 <- sum(w * S_Z1_i) / sum(w)
  S_Z0 <- sum(w * S_Z0_i) / sum(w)
  
  # --- SPCE ---
  SPCE <- S_Z1 - S_Z0
  
  return(list(
    S_Z1 = S_Z1,
    S_Z0 = S_Z0,
    SPCE = SPCE,
    weights = w,
    ps_hat = ps_hat
  ))
}


surv_stoEM_ipw_cox <- function(beta_final,
                               basehaz,     
                               t0,
                               Z, X,
                               stabilize = TRUE) {
  n <- nrow(X)
  
  # --- Estimate propensity scores ---
  ps_model <- glm(Z ~ ., data = data.frame(Z = Z, X), family = binomial)
  ps_hat <- predict(ps_model, type = "response")
  
  # --- Compute IP weights ---
  if (stabilize) {
    pZ <- mean(Z)
    w <- ifelse(Z == 1, pZ / ps_hat, (1 - pZ) / (1 - ps_hat))
  } else {
    w <- ifelse(Z == 1, 1 / ps_hat, 1 / (1 - ps_hat))
  }
  
  # --- Extract Cox coefficients ---
  beta_x <- beta_final[1:ncol(X)]
  beta_z <- beta_final["Z"]
  beta_u <- beta_final[length(beta_final)]  
  
  # --- Linear predictors (marginalize over U = 0.5) ---
  lp_Z1 <- as.vector(X %*% beta_x + beta_z * 1 + beta_u * 0.5)
  lp_Z0 <- as.vector(X %*% beta_x + beta_z * 0 + beta_u * 0.5)
  
  # --- Baseline cumulative hazard at t0 ---
  H0_t0 <- approx(basehaz$time, basehaz$hazard, xout = t0, method = "linear", rule = 2)$y
  
  # --- Individual survival probabilities ---
  S_Z1_i <- exp(-H0_t0 * exp(lp_Z1))
  S_Z0_i <- exp(-H0_t0 * exp(lp_Z0))
  
  # --- Weighted marginal survival ---
  S_Z1 <- sum(w * S_Z1_i) / sum(w)
  S_Z0 <- sum(w * S_Z0_i) / sum(w)
  
  SPCE <- S_Z1 - S_Z0
  
  return(list(
    SPCE = SPCE,
    S_Z1 = S_Z1,
    S_Z0 = S_Z0,
    weights = w,
    ps_hat = ps_hat
  ))
}


one_sim_run <- function(i) {
  set.seed(2000 + i)
  n <- 500; tau <- 5.5; t_pred <- 2; n_boot <- 100  
  # 1. Data generation 
  X1 <- rbinom(n, size = 1, prob = 0.4)
  X2 <- rbinom(n, size = 1, prob = 0.6)
  U1 <- rbinom(n, size = 1, prob = 0.5)
  alpha_0  <- 0.1; alpha_x1 <- 0.3; alpha_x2 <- 0.5; alpha_u  <- -0.8
  linpred_A <- alpha_0 + alpha_x1 * X1 + alpha_x2 * X2 + alpha_u * U1
  pi_A      <- 1 / (1 + exp(-linpred_A))
  A         <- rbinom(n, size = 1, prob = pi_A)
  eta_intercept_a0 <- 0.2; eta_intercept_a1 <- 0.7; eta_x1 <- -0.1; eta_x2 <- 0.4; eta_u  <- -0.8
  sigma0 <- sigma1 <- 1/2; alpha0 <- alpha1 <- 1/sigma0
  lp_a1   <- eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U1
  lp_a0   <- eta_intercept_a0 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U1
  lambda1 <- exp(-lp_a1); lambda0 <- exp(-lp_a0)
  D_a0 <- rweibull(n, shape = alpha0, scale = 1/lambda0)
  D_a1 <- rweibull(n, shape = alpha1, scale = 1/lambda1)
  C_a0 <- runif(n, 0.1, tau); C_a1 <- runif(n, 0.1, tau)
  M <- ifelse(A == 1, pmin(tau, C_a1, D_a1), pmin(tau, C_a0, D_a0))
  delta <- ifelse(A == 1, as.numeric(D_a1 <= pmin(tau, C_a1)), as.numeric(D_a0 <= pmin(tau, C_a0)))
  data_sim <- data.frame(M = M, delta = delta, A = A, X1 = X1, X2 = X2, U1 = U1)
  x_mat <- as.matrix(data_sim[, c("X1", "X2")])
  
  # 2. Model fitting as in the loop
  ## ===  regular EM algorithm === ##
  
  ## =============== aft  =============== ##
  
  result_em_aft <- emU_surv_aft(
    t      = data_sim$M,
    d      = data_sim$delta,
    z      = data_sim$A,
    x      = x_mat,
    theta  = 0.5,
    gammat = -0.8,
    gammaz = -0.8
  )
  
  result_em_cox <- emU_surv_cox(
    t      = data_sim$M,
    d      = data_sim$delta,
    z      = data_sim$A,
    x      = x_mat,
    theta  = 0.5,
    zetat  = -0.8,
    zetaz  = -0.8
  )
  
  #compute spce by weighting regular EM algorithm
  res_spce_weight_aft <- compute_SPCE_EM_aft(
    res_em = result_em_aft,
    X      = x_mat,
    Z      = data_sim$A,
    M      = data_sim$M,
    delta  = data_sim$delta,
    t0     = t_pred
  )
  
  res_spce_weight_cox <- compute_SPCE_EM_cox(
    res_em = result_em_cox,
    X      = x_mat,
    Z      = data_sim$A,
    M      = data_sim$M,
    delta  = data_sim$delta,
    t0     = t_pred
  )
  
  ## ===  Stochastic EM algorithm === ##
  
  stoEM_aft <- surv_stoEM_aft(t= data_sim$M,
                              d= data_sim$delta,
                              Z= data_sim$A,
                              X= x_mat,
                              zetat=-0.8, zetaz=-0.8, B=10, theta=0.5)
  
  
  stoEM_cox <- surv_stoEM_cox(t= data_sim$M,
                              d= data_sim$delta,
                              Z= data_sim$A,
                              X= x_mat,
                              zetat=-0.8, zetaz=-0.8, B=10, theta=0.5)
  
  #stochastic EM for SPCE by weighting
  stoEM_SPCE_W_aft <- surv_stoEM_ipw_aft(beta_final=stoEM_aft$beta, 
                                         sigma_final=stoEM_aft$sigma, 
                                         t= data_sim$M,
                                         d= data_sim$delta,
                                         Z= data_sim$A,
                                         X= x_mat,
                                         t0=t_pred,
                                         stabilize=TRUE)
  
  stoEM_SPCE_W_cox <- surv_stoEM_ipw_cox(beta_final=stoEM_cox$beta, 
                                         basehaz=stoEM_cox$basehaz,
                                         t0=t_pred, Z= data_sim$A,X= x_mat,
                                         stabilize=TRUE)
  
  ## --------------------------------------- ##
  ## === BOOTSTRAP inside the iteration  === ##
  ## --------------------------------------- ##
  B_boot <- n_boot
  boot_spce_EM_aft  <- numeric(B_boot)
  boot_spce_EM_cox  <- numeric(B_boot)
  boot_spce_Sto_aft <- numeric(B_boot)
  boot_spce_Sto_cox <- numeric(B_boot)
  for (b in 1:B_boot) {
    boot_idx <- sample(1:n, size = n, replace = TRUE)
    data_boot <- data_sim[boot_idx, ]
    x_boot <- as.matrix(data_boot[, c("X1", "X2")])
    
    ## ---------- EM AFT ----------
    boot_em_aft <- emU_surv_aft(
      t      = data_boot$M,
      d      = data_boot$delta,
      z      = data_boot$A,
      x      = x_boot,
      theta  = 0.5,
      gammat = -0.8,
      gammaz = -0.8
    )
    
    boot_spce_EM_aft[b] <- compute_SPCE_EM_aft(
      res_em = boot_em_aft,
      X      = x_boot,
      Z      = data_boot$A,
      M      = data_boot$M,
      delta  = data_boot$delta,
      t0     = t_pred
    )$SPCE
    
    ## ---------- EM Cox ----------
    boot_em_cox <- emU_surv_cox(
      t      = data_boot$M,
      d      = data_boot$delta,
      z      = data_boot$A,
      x      = x_boot,
      theta  = 0.5,
      zetat  = -0.8,
      zetaz  = -0.8
    )
    
    boot_spce_EM_cox[b] <- compute_SPCE_EM_cox(
      res_em = boot_em_cox,
      X      = x_boot,
      Z      = data_boot$A,
      M      = data_boot$M,
      delta  = data_boot$delta,
      t0     = t_pred
    )$SPCE
    
    ## ---------- StoEM AFT ----------
    boot_sto_aft <- surv_stoEM_aft(
      t = data_boot$M,
      d = data_boot$delta,
      Z = data_boot$A,
      X = x_boot,
      zetat = -0.8, zetaz = -0.8, B = 10, theta = 0.5
    )
    
    boot_spce_Sto_aft[b] <- surv_stoEM_ipw_aft(
      beta_final = boot_sto_aft$beta,
      sigma_final = boot_sto_aft$sigma,
      t = data_boot$M,
      d = data_boot$delta,
      Z = data_boot$A,
      X = x_boot,
      t0 = t_pred,
      stabilize = TRUE
    )$SPCE
    
    ## ---------- StoEM Cox ----------
    boot_sto_cox <- surv_stoEM_cox(
      t = data_boot$M,
      d = data_boot$delta,
      Z = data_boot$A,
      X = x_boot,
      zetat = -0.8, zetaz = -0.8, B = 10, theta = 0.5
    )
    
    boot_spce_Sto_cox[b] <- surv_stoEM_ipw_cox(
      beta_final = boot_sto_cox$beta,
      basehaz    = boot_sto_cox$basehaz,
      Z = data_boot$A,
      X = x_boot,
      t0 = t_pred,
      stabilize = TRUE
    )$SPCE
  }
  
  c(
    EM_AFT = res_spce_weight_aft$SPCE,
    EM_Cox = res_spce_weight_cox$SPCE,
    StoEM_AFT = stoEM_SPCE_W_aft$SPCE,
    StoEM_Cox = stoEM_SPCE_W_cox$SPCE,
    EM_AFT_SD = sd(boot_spce_EM_aft),
    EM_Cox_SD = sd(boot_spce_EM_cox),
    StoEM_AFT_SD = sd(boot_spce_Sto_aft),
    StoEM_Cox_SD = sd(boot_spce_Sto_cox)
  )
}

n_sim <- 10
cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(cores)

# Export all custom functions and variables!
clusterExport(cl, varlist = ls())
clusterEvalQ(cl, { library(survival) })

# Parallel run:
res_list <- parLapplyLB(cl, 1:n_sim, one_sim_run)
stopCluster(cl)

results_df <- as.data.frame(do.call(rbind, res_list))
colnames(results_df) <- c("EM_AFT", "EM_Cox", "StoEM_AFT", "StoEM_Cox",
                          "EM_AFT_SD", "EM_Cox_SD", "StoEM_AFT_SD", "StoEM_Cox_SD")
str(results_df)  
write.csv(results_df, "sim_results.csv", row.names = FALSE)

