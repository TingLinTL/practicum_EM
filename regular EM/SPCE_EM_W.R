compute_SPCE_EM_aft_W <- function(res_em, X, Z, M, delta, t0) {
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


compute_SPCE_EM_cox_W <- function(res_em, X, Z, M, delta, t0) {

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

