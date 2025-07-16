surv_stoEM_SPCE <- function(data, zetat, zetaz, B = 100, theta = 0.5, t0 = 5) {
  #t0 is the SPCE time point.
    
  t = data$t
  d = data$d
  Z = data$Z
  X = data.matrix(data$X)
  #nx = dim(X)[2]
  n = length(t)
  
  SPCE_draws <- numeric(B)   # storage for SPCE per draw
  
  for (j in 1:B){
    Usim = SimulateU_surv(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta) 
    
    
    t1.fit<- survreg(Surv(t,d) ~ X + Z + Usim$U, dist = "weibull")
    beta_hat  <- t1.fit$coef      # β0, β_X, β_Z, β_U
    sigma_hat <- t1.fit$scale
    alpha     <- 1 / sigma_hat   # Weibull shape
    
    mu_Z1 <- cbind(1, X, 1, Usim$U) %*% beta_hat
    mu_Z0 <- cbind(1, X, 0, Usim$U) %*% beta_hat
    
    # Weibull *rate* λ = exp(-μ)
    lambda_Z1 <- exp(-mu_Z1)
    lambda_Z0 <- exp(-mu_Z0)
    
    # survival at time t0
    S_Z1 <- exp(-(lambda_Z1 * t0)^alpha)
    S_Z0 <- exp(-(lambda_Z0 * t0)^alpha)
    
    SPCE_draws[j] <- mean(S_Z1 - S_Z0)
    
  }
  ## Monte Carlo average & SE
  SPCE_mean <- mean(SPCE_draws)
  SPCE_se   <- sqrt(var(SPCE_draws) * (1 + 1/B))  # uncertainty comes only from between-draw variability
  
  return(list(
    SPCE = SPCE_mean,
    SPCE_se = SPCE_se,
    all_draws = SPCE_draws
  ))
}
