surv_stoEM_cox <- function(t, d, Z, X, zetat, zetaz, B, theta) {
  nx <- ncol(X)
  n <- length(t)
  
  final_beta <- NULL
  final_fit  <- NULL
  
  for (j in 1:B) {
    Usim <- SimulateU_surv_cox(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta)
    
    t1.fit <- coxph(Surv(t, d) ~ X + Z + Usim$U)
    
    final_beta <- t1.fit$coef
    final_fit  <- t1.fit  
  }
  
  return(list(
    beta    = final_beta,
    fit     = final_fit,
    basehaz = basehaz(final_fit, centered = FALSE)
  ))
}
