# Function to compute true SPCE based on CoxPH data generation process
true_spce_population <- function(t0 = 2, 
                                 pX1 = 0.4, pX2 = 0.6, pU = 0.5,
                                 eta_A = -0.5, eta_x1 = -0.1, eta_x2 = 0.4, eta_u = -0.8) {
  
  vals <- 0
  
  # Iterate over all combinations of X1, X2, and U (since they are binary)
  for (x1 in 0:1) for (x2 in 0:1) for (u in 0:1) {
    
    # Compute the weight for this combination of covariates
    w <- (if (x1 == 1) pX1 else 1 - pX1) *
      (if (x2 == 1) pX2 else 1 - pX2) *
      (if (u == 1) pU else 1 - pU)
    
    # Calculate the linear predictor for covariates
    bx <- eta_x1 * x1 + eta_x2 * x2 + eta_u * u
    
    # Hazard for treatment group (Z = 1) (note: eta_A contributes to treatment effect)
    lambda1 <- exp(eta_A + bx)  # Hazard for treatment
    
    # Hazard for control group (Z = 0) (no treatment effect, so no eta_A term)
    lambda0 <- exp(bx)  # Hazard for control
    
    # Compute the survival probability for treatment group (Z=1) at t_0
    S1 <- exp(-lambda1 * t0)
    
    # Compute the survival probability for control group (Z=0) at t_0
    S0 <- exp(-lambda0 * t0)
    
    # Compute the contribution to the true SPCE from this covariate combination
    vals <- vals + w * (S1 - S0)
  }
  
  # Return the computed true SPCE
  return(vals)
}

# Example: Compute True SPCE at time t0 = 2
true_spce_population(t0 = 2, eta_A = -0.5, eta_x1 = -0.1, eta_x2 = 0.4, eta_u = -0.8)



#true 0.1573018