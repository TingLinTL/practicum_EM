surv_stoEM_ipw <- function(data, zetat, zetaz, B = 100, theta = 0.5, t0 = 5) {
  # t0 = SPCE time point
  
  t <- data$t
  d <- data$d
  Z <- data$Z
  X <- data.matrix(data$X)
  n <- length(t)
  
  SPCE_draws <- numeric(B)
  
  for (j in 1:B) {
    Usim <- SimulateU_surv(
      t, d, Z, X,
      zetat = zetat,
      zetaz = zetaz,
      theta = theta
    )
    
    #fit treatment model P(Z|X,U)
    fit_z <- glm(Z ~ X + Usim$U, family = binomial("logit"))
    ps <- fit_z$fitted.values  # predicted P(Z=1|X,U)
    
    # Stabilized weights
    pZ <- mean(Z)  # marginal probability of treatment
    w <- (pZ * Z / ps) + ((1 - pZ) * (1 - Z) / (1 - ps))
    
    #fit marginal Weibull AFT for survival ~ Z ONLY, weighted
    fit_w <- survreg(Surv(t, d) ~ Z, weights = w, dist = "weibull")
    beta_hat  <- fit_w$coef   # intercept + beta_Z
    sigma_hat <- fit_w$scale
    alpha     <- 1 / sigma_hat
    
    #predict marginal survival at t0 under Z=1 vs Z=0
    # marginal model,just use the intercept & beta_Z
    
    mu_Z1 <- beta_hat[1] + beta_hat["Z"] * 1
    mu_Z0 <- beta_hat[1] + beta_hat["Z"] * 0
    
    lambda_Z1 <- exp(-mu_Z1)
    lambda_Z0 <- exp(-mu_Z0)
    
    S_Z1 <- exp(-(lambda_Z1 * t0)^alpha)
    S_Z0 <- exp(-(lambda_Z0 * t0)^alpha)
    
    SPCE_draws[j] <- (S_Z1 - S_Z0)
  }
  
  # Monte Carlo mean & SE
  SPCE_mean <- mean(SPCE_draws)
  SPCE_se   <- sqrt(var(SPCE_draws) * (1 + 1/B))
  
  return(list(
    SPCE = SPCE_mean,
    SPCE_se = SPCE_se,
    all_draws = SPCE_draws
  ))
}
