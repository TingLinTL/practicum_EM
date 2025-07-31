SimulateU_surv_aft <- function(t, d, z, x, zetat, zetaz, theta, iter = 20){
  #t is a vector of n, time to event
  #d is a vector of n, indicator of event
  #z is a vector of n, treatment
  #x is a matrix, covariates
  #zetat is a scaler, zetaz is a scaler, sensitivity parameters
  
  n = length(t) #number of observations
  x = data.matrix(x)
  nx = dim(x)[2] #number of observed covariates
  
  p = theta
  U = rbinom(n,1,p)
  #Upath = U (useful for checking convergence)
  
  for(j in 1:iter){
      #fit time to event
      U.fit1 = survreg(Surv(t,d) ~ x + z + U, dist = "weibull")
      t.coef1 = U.fit1$coef
      t.coef1[length(t.coef1)]  = zetat
      sigma <- U.fit1$scale
      
      #fit treatment
      z.coef = glm(z ~ x + U, family = binomial(link = "logit"), control = glm.control(epsilon = 1e-6, maxit = 50))$coef
      z.coef[length(z.coef)] = zetaz
    
    
    t.coef1[is.na(t.coef1)] = 0
    z.coef[is.na(z.coef)] = 0
    
    # Logistic treatment model linear predictors
    logit_zu1 <- (1 - plogis(cbind(1, x, 1) %*% matrix(z.coef, ncol = 1)))^(1 - z) *
      plogis(cbind(1, x, 1) %*% matrix(z.coef, ncol = 1))^z
    
    logit_zu0 <- (1 - plogis(cbind(1, x, 0) %*% matrix(z.coef, ncol = 1)))^(1 - z) *
      plogis(cbind(1, x, 0) %*% matrix(z.coef, ncol = 1))^z
    
    # weibull aft linear prediction under U = 1 and U = 0
    mu_u1 <- cbind(1, x, z, 1) %*% matrix(t.coef1, ncol = 1)
    mu_u0 <- cbind(1, x, z, 0) %*% matrix(t.coef1, ncol = 1)
    alpha_weibull <- 1 / sigma
    
    # Lambda (rate-like scale)
    lambda_u1 <- exp(-mu_u1)
    lambda_u0 <- exp(-mu_u0)
    
    #weibull density and survival
    lik_u1 <- ( alpha_weibull * lambda_u1^alpha_weibull * t^(alpha_weibull - 1) )^d *
      exp( - (lambda_u1 * t)^alpha_weibull )
    
    lik_u0 <- ( alpha_weibull * lambda_u0^alpha_weibull * t^(alpha_weibull - 1) )^d *
      exp( - (lambda_u0 * t)^alpha_weibull )
    
    
    ptzu1 = logit_zu1*theta* lik_u1
    ptzu0 = logit_zu0*(1-theta)* lik_u0
    
    p = ptzu1/(ptzu1 + ptzu0)
    p[ptzu1==0 & ptzu0==0] = 0
    
    U = rbinom(n,1,p)
    #Upath = cbind(Upath, U)
  }
  
  return(list(U = U, p = p))
}