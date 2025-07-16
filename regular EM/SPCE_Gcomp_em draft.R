compute_SPCE_EM <- function(res_em, X, t0) {
  # unpack stable params from EM result
  beta_final <- res_em$t.coef1
  sigma_final <- res_em$sigma
  p_final <- res_em$p
  
  # predict counterfactual linear predictors
  mu_Z1 <- cbind(1, X, 1, p_final) %*% matrix(beta_final, ncol = 1)  # Z=1 for all
  mu_Z0 <- cbind(1, X, 0, p_final) %*% matrix(beta_final, ncol = 1)  # Z=0 for all
  
  lambda_Z1 <- exp(-mu_Z1)
  lambda_Z0 <- exp(-mu_Z0)
  
  alpha_weibull <- 1 / sigma_final
  
  # survival at time t0
  S_Z1 <- exp(-(lambda_Z1 * t0)^alpha_weibull)
  S_Z0 <- exp(-(lambda_Z0 * t0)^alpha_weibull)
  
  SPCE <- mean(S_Z1 - S_Z0)
  
  return(list(
    SPCE = SPCE,
    S_Z1 = mean(S_Z1),
    S_Z0 = mean(S_Z0)
  ))
}
