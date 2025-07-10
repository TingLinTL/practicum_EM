emU_surv <- function(t, d, z, x, gammat, gammaz, theta = 0.5, iter = 20){
  fn_outcome <- function(params, t, d, x, z, p, gammat) {
    # unpack parameters
    beta   <- params[1:ncol(x)]        # regression coefficients
    alpha  <- params[ncol(x) + 1]      # treatment effect
    sigma  <- params[ncol(x) + 2]      # scale parameter (> 0)
    
    # linear predictor
    xb <- as.vector(x %*% beta)
    
    # E[u_i] = p
    lp_meanU <- xb + alpha * z + gammat * p
    
    # first term
    term1 <- d * (-log(sigma) + ((1 / sigma - 1) * log(t)) - (1 / sigma) * lp_meanU)
    
    # second term
    exp0 <- exp((log(t) - xb - alpha * z - gammat * 0) / sigma)
    exp1 <- exp((log(t) - xb - alpha * z - gammat * 1) / sigma)
    term2 <- p * exp1 + (1 - p) * exp0
    
    return(-sum(term1 - term2))  # negative log-likelihood
  }
  
  
  fn_treat <- function(beta, x, z, p, gammaz) {
    # linear predictors
    lp1 <- cbind(1, x) %*% beta + gammaz  # when U_i = 1
    lp0 <- cbind(1, x) %*% beta           # when U_i = 0
    
    term1 <- z * (log(1 / (1 + exp(-lp1))) * p + log(1 / (1 + exp(-lp0))) * (1 - p))
    term2 <- (1 - z) * (log(1 - 1 / (1 + exp(-lp1))) * p + log(1 - 1 / (1 + exp(-lp0))) * (1 - p))
    
    #negative tretment log-likelihood
    return(-sum(term1 + term2))
  }
  
  x = data.matrix(x)
  n = length(t) #number of observations
  nx = dim(x)[2] #number of covariates ,or nx = ncol(x)
  
  #initialize parameters
  #fit time to event
  U.fit1 <- survreg(Surv(t, d) ~ x + z, dist = "weibull")
  # survreg fits AFT model: log(T) = Xβ + error
  # extract initial values for next iteration
  t.coef1 = c(coef(U.fit1), gammat) #beta0...beta_nx, alpha,gammat
  sigma <- U.fit1$scale
  
  #fit treatment
  Z.fit1 <- glm(z ~ x, family = binomial(link = "logit"))
  z.coef <- c(coef(Z.fit1), gammaz) #beta_z1,beta_z2...beta_znx,gammaz
  
  t.coef1[is.na(t.coef1)] = 0
  z.coef[is.na(z.coef)] = 0
  
  p = rep(0,n)
  
  for(j in 1:iter){
    
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
    
    # Logistic treatment model linear predictors
    logit_zu1 <- (1 - plogis(cbind(1, x, 1) %*% matrix(z.coef, ncol = 1)))^(1 - z) *
      plogis(cbind(1, x, 1) %*% matrix(z.coef, ncol = 1))^z
    
    logit_zu0 <- (1 - plogis(cbind(1, x, 0) %*% matrix(z.coef, ncol = 1)))^(1 - z) *
      plogis(cbind(1, x, 0) %*% matrix(z.coef, ncol = 1))^z
    
    # Posterior numerator
    ptzu1 <- lik_u1 * logit_zu1 * theta
    ptzu0 <- lik_u0 * logit_zu0 * (1-theta)
    
    # Posterior of U = 1
    p <- ptzu1 / (ptzu1 + ptzu0) #vector of length n, different p_i for each individual
    p[ptzu1 == 0 & ptzu0 == 0] <- 0  
    
    #fit time to event
    init_params <- c(coef(U.fit1), sigma)
    outcome_fit <- optim(par = init_params,
                         fn = fn_outcome,
                         t = t, d = d, x = x, z = z, p = p, gammat=gammat,
                         method = "BFGS",
                         control = list(maxit = 1000))
    params <- outcome_fit$par
    t.coef1 <- c(params[1:(ncol(x) + 2)], gammat)
    sigma <- params[ncol(x) + 3]
      
    #fit treatment
    z.fit <- optim(par = z.coef[1:(nx+1)],
                   fn = fn_treat,
                   x = x, z = z, p = p, gammaz = gammaz)
    z.coef <- c(z.fit$par, gammaz)
    
    t.coef1[is.na(t.coef1)] = 0
    z.coef[is.na(z.coef)] = 0
    }
  return(list(p = p, t.coef1 = t.coef1, z.coef = z.coef, sigma = sigma))
  }
  
  
