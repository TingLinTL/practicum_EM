compute_SPCE_G_aft <- function(beta_final, sigma_final, X, t0) {
  alpha_weibull <- 1 / sigma_final
  
  # --- Z=1 ---
  mu_1_U1 <- cbind(1, X, 1, 1) %*% beta_final
  mu_1_U0 <- cbind(1, X, 1, 0) %*% beta_final
  
  lam_1_U1 <- exp(-mu_1_U1)
  lam_1_U0 <- exp(-mu_1_U0)
  S_1_U1 <- exp(-(lam_1_U1 * t0)^alpha_weibull)
  S_1_U0 <- exp(-(lam_1_U0 * t0)^alpha_weibull)
  
  # --- Z=0 ---
  mu_0_U1 <- cbind(1, X, 0, 1) %*% beta_final
  mu_0_U0 <- cbind(1, X, 0, 0) %*% beta_final
  
  lam_0_U1 <- exp(-mu_0_U1)
  lam_0_U0 <- exp(-mu_0_U0)
  S_0_U1 <- exp(-(lam_0_U1 * t0)^alpha_weibull)
  S_0_U0 <- exp(-(lam_0_U0 * t0)^alpha_weibull)
  
  # --- Equal mixture 0.5 ---
  S_Z1_i <- 0.5 * S_1_U1 + 0.5 * S_1_U0
  S_Z0_i <- 0.5 * S_0_U1 + 0.5 * S_0_U0
  
  # --- Marginal survival ---
  S_Z1 <- mean(S_Z1_i)
  S_Z0 <- mean(S_Z0_i)
  
  # --- SPCE ---
  SPCE <- S_Z1 - S_Z0
  
  return(list(S_Z1 = S_Z1, S_Z0 = S_Z0, SPCE = SPCE))
}
