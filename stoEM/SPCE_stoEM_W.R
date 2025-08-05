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
