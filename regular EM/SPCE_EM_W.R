compute_IPW_SPCE_from_EM <- function(res_em, X, Z, M, delta, t0) {
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
