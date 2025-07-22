surv_stoEM_trace <- function(t, d, Z, X, zetat, zetaz, B, theta) {

  n <- length(t)
  
  ## ---- fit initial model (no U) ----
  init_fit <- survreg(Surv(t, d) ~ X + Z, dist = "weibull")
  init_beta <- init_fit$coef      # β0, β_X, β_Z
  init_sigma <- init_fit$scale    # Weibull scale
  init_Ucoef <- 1                 # manually set initial U coefficient
  
  ## ---- total columns = (init betas) + U + sigma ----
  n_coef <- length(init_beta) + 2  
  
  ## ---- storage for coefficient trajectory ----
  coef_trace <- matrix(NA, nrow = B + 1, ncol = n_coef)
  colnames(coef_trace) <- c(names(init_beta), "U_coeff", "sigma_hat")
  
  ## ---- store initial model ----
  coef_trace[1, ] <- c(init_beta, init_Ucoef, init_sigma)
  

  for (j in 1:B){
    Usim = SimulateU_surv(t, d, Z, X, zetat = zetat, zetaz = zetaz, theta = theta) 
    
    
    t1.fit<- survreg(Surv(t,d) ~ X + Z + Usim$U, dist = "weibull")
    
    # extract parameters
    beta_hat  <- t1.fit$coef        # β0, β_X, β_Z, β_U
    sigma_hat <- t1.fit$scale       # Weibull scale parameter
    #alpha <- 1 / sigma_hat        # Weibull shape 
    
    # store into trace matrix (shifted by +1)
    coef_trace[j + 1, ] <- c(beta_hat, sigma_hat)
  }
  return(coef_trace)
}
