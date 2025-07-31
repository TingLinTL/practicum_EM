compute_SPCE_EM_aft <- function(res_em, X, t0) {
  beta_final  <- res_em$t.coef1
  sigma_final <- res_em$sigma
  p_final     <- res_em$p
  
  alpha_weibull <- 1 / sigma_final
  
  # Survival under Z=1 vs Z=0, for both U=1 and U=0 separately
  # Z=1
  mu_1_U1 <- cbind(1, X, 1, 1) %*% beta_final  # U=1
  mu_1_U0 <- cbind(1, X, 1, 0) %*% beta_final  # U=0
  lam_1_U1 <- exp(-mu_1_U1)
  lam_1_U0 <- exp(-mu_1_U0)
  S_1_U1 <- exp(-(lam_1_U1 * t0)^alpha_weibull)
  S_1_U0 <- exp(-(lam_1_U0 * t0)^alpha_weibull)
  
  # Z=0
  mu_0_U1 <- cbind(1, X, 0, 1) %*% beta_final
  mu_0_U0 <- cbind(1, X, 0, 0) %*% beta_final
  lam_0_U1 <- exp(-mu_0_U1)
  lam_0_U0 <- exp(-mu_0_U0)
  S_0_U1 <- exp(-(lam_0_U1 * t0)^alpha_weibull)
  S_0_U0 <- exp(-(lam_0_U0 * t0)^alpha_weibull)
  
  # Posterior mixture for each subject
  S_Z1_i <- p_final * S_1_U1 + (1 - p_final) * S_1_U0
  S_Z0_i <- p_final * S_0_U1 + (1 - p_final) * S_0_U0
  
  # Average over sample for marginal survival
  S_Z1 <- mean(S_Z1_i)
  S_Z0 <- mean(S_Z0_i)
  
  SPCE <- S_Z1 - S_Z0
  
  return(list(
    SPCE = SPCE,
    S_Z1 = S_Z1,
    S_Z0 = S_Z0
  ))
}

# the following is pluggin-in, not marginalize over U
# compute_SPCE_EM <- function(res_em, X, t0) {
#   # unpack stable params from EM result
#   beta_final <- res_em$t.coef1
#   sigma_final <- res_em$sigma
#   p_final <- res_em$p
#   
#   # predict counterfactual linear predictors
#   mu_Z1 <- cbind(1, X, 1, p_final) %*% matrix(beta_final, ncol = 1)  # Z=1 for all
#   mu_Z0 <- cbind(1, X, 0, p_final) %*% matrix(beta_final, ncol = 1)  # Z=0 for all
#   
#   lambda_Z1 <- exp(-mu_Z1)
#   lambda_Z0 <- exp(-mu_Z0)
#   
#   alpha_weibull <- 1 / sigma_final
#   
#   # survival at time t0
#   S_Z1 <- exp(-(lambda_Z1 * t0)^alpha_weibull)
#   S_Z0 <- exp(-(lambda_Z0 * t0)^alpha_weibull)
#   
#   SPCE <- mean(S_Z1 - S_Z0)
#   
#   return(list(
#     SPCE = SPCE,
#     S_Z1 = mean(S_Z1),
#     S_Z0 = mean(S_Z0)
#   ))
# }


compute_SPCE_EM_cox <- function(res_em, X, t0) {
  beta_all <- res_em$t.coef1
  pU       <- res_em$p
  baseH    <- res_em$basehaz1
  
  # Separate coefficients
  zetat <- tail(beta_all, 1)       
  beta  <- beta_all[-length(beta_all)]  # drop zetat
  beta_z <- beta[1]
  beta_x <- beta[-1]
  
  # linear predictors under Z=1 and Z=0
  lp_Z1 <- as.vector(X %*% beta_x + beta_z * 1 + zetat * pU)
  lp_Z0 <- as.vector(X %*% beta_x + beta_z * 0 + zetat * pU)
  
  # cumulative baseline hazard at t0
  H0_t0 <- baseH$hazard[which.min(abs(baseH$time - t0))]
  
  # compute survival probabilities
  S_Z1_i <- exp(-H0_t0 * exp(lp_Z1))
  S_Z0_i <- exp(-H0_t0 * exp(lp_Z0))
  
  # average over individuals
  S_Z1 <- mean(S_Z1_i)
  S_Z0 <- mean(S_Z0_i)
  SPCE <- S_Z1 - S_Z0
  
  return(list(
    SPCE = SPCE,
    S_Z1 = S_Z1,
    S_Z0 = S_Z0
  ))
}

