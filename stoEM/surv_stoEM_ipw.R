surv_stoEM_ipw  <- function(beta_final, sigma_final,
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
  