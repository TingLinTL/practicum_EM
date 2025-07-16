compute_IPW_from_EM <- function(res_em, X, Z) {
  # unpack converged treatment model parameters
  z.coef_final <- res_em$z.coef    
  U_final <- res_em$p             
  
  # fitted propensity scores P(Z=1|X,U)
  linpred <- cbind(1, X, U_final) %*% matrix(z.coef_final, ncol = 1)
  ps <- plogis(linpred)
  
  # stabilized weights
  pZ <- mean(Z) # marginal treatment probability
  w <- (pZ * Z / ps) + ((1 - pZ) * (1 - Z) / (1 - ps))
  
  # truncate extreme weights
  w <- pmin(w, 10)
  w <- pmax(w, 0.1)
  
  return(w)
}
