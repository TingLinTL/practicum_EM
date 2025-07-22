surv_stoEM <- function(t, d, Z, X, zetat, zetaz, B, theta) {
  
  n <- length(t)
  
  final_beta <- NULL
  final_sigma <- NULL
  
  for (j in 1:B){
    Usim = SimulateU_surv(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta) 
    
    
    t1.fit<- survreg(Surv(t,d) ~ X + Z + Usim$U, dist = "weibull")
    
    # extract parameters
    final_beta <- t1.fit$coef       # β0, β_X, β_Z, β_U
    final_sigma <- t1.fit$scale      # Weibull scale parameter
    #alpha <- 1 / sigma_hat        # Weibull shape 
    
  }
  # return only the *final* coefficients after B iterations
  return(list(
    beta = final_beta,    # includes β0, β_X, β_Z, β_U
    sigma = final_sigma   # Weibull scale
  ))
}
