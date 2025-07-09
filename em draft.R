emU_surv <- function(t, d, z, x, gammat, gammaz, theta = 0.5, iter = 20){
  
  fn_outcome <- function(params, t, d, x, z, p) {
    # unpack parameters
    beta <- params[1:ncol(x)]        # regression coefficients
    alpha <- params[ncol(x) + 1]     # treatment effect
    gammat <- params[ncol(x) + 2]     # unmeasured confounder effect
    sigma <- params[ncol(x) + 3]     # scale parameter (> 0)
    
    # linear predictor
    xb <- as.vector(x %*% beta)
    
    # E[u_i] = p
    lp_meanU <- xb + alpha * z + gamma * p
    
    # first term
    term1 <- d * (-log(sigma) + ((1 / sigma - 1) * log(t)) - (1 / sigma) * lp_meanU )
    
    #second term
    exp0 <- exp((log(t) - xb - alpha * z - gamma * 0) / sigma)
    exp1 <- exp((log(t) - xb - alpha * z - gamma * 1) / sigma)
    term2<- p * exp1 + (1 - p) * exp0
    
    # Negative Q1 (for minimization)
    return(-sum(term1 - term2))
  }
  
    
  fn_treat <- function(beta, x, z, p, gammaz) {
    # linear predictors
    lp1 <- x %*% beta + gammaz  # when U_i = 1
    lp0 <- x %*% beta           # when U_i = 0
    #Q2
    term1 <- z * (log(1 / (1 + exp(-lp1))) * p + log(1 / (1 + exp(-lp0))) * (1 - p))
    term2 <- (1 - z) * (log(1 - 1 / (1 + exp(-lp1))) * p + log(1 - 1 / (1 + exp(-lp0))) * (1 - p))
      
    #negative log-likelihood
    return(-sum(term1 + term2))
  }
  
  
  #initialize parameters
  #fit time to event
  U.fit1 <- survreg(Surv(t, d) ~ z + x, dist = "weibull")
  # survreg fits AFT model: log(T) = Xβ + error
  # extract initial values for next iteration
  t.coef1 = c(coef(U.fit1), U.fit1$scale, gammat)
  
  #fit treatment
  Z.fit1 <- glm(z ~ x, family = binomial(link = "logit"))
  z.coef <- c(coef(Z.fit1), gammaz)
}