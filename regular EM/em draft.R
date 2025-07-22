emU_surv <- function(t, d, z, x, gammat, gammaz, theta = 0.5, iter = 50){
  fn_outcome <- function(params, t, d, x, z, p, gammat) {
    # Unpack params
    beta      <- params[1:(ncol(x)+1)]
    alpha     <- params[ncol(x)+2]
    log_sigma <- params[ncol(x)+3]
    sigma     <- exp(log_sigma)  # always >0
    
    # Safe log(t)
    log_t <- log(pmax(t, 1e-8))
    
    # Linear predictor
    xb <- cbind(1, x) %*% beta
    lp_meanU <- xb + alpha*z + gammat*p
    
    # First term
    term1 <- d * (-log(sigma) + ((1/sigma - 1)*log_t) - (1/sigma)*lp_meanU)
    
    # Second term: prevent exp overflow
    eta0 <- (log_t - xb - alpha*z - gammat*0)/sigma
    eta1 <- (log_t - xb - alpha*z - gammat*1)/sigma
    eta0 <- pmin(eta0, 50)
    eta1 <- pmin(eta1, 50)
    
    exp0 <- exp(eta0)
    exp1 <- exp(eta1)
    term2 <- p * exp1 + (1 - p) * exp0
    term2 <- pmax(term2, 1e-12)  # avoid log(0)
    
    nll <- -sum(term1 - term2)
    
    # If nll is non-finite, return huge penalty
    if (!is.finite(nll)) return(1e12)
    
    return(nll)
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
  init_params  <- c(coef(U.fit1), log(sigma)) 
  
  #fit treatment
  Z.fit1 <- glm(z ~ x, family = binomial(link = "logit"))
  z.coef <- c(coef(Z.fit1), gammaz) #beta_z1,beta_z2...beta_znx,gammaz
  
  t.coef1[is.na(t.coef1)] = 0
  z.coef[is.na(z.coef)] = 0
  
  p = rep(0,n)
  coef_trace <- matrix(NA, nrow = iter, ncol = length(init_params))  # track beta, alpha, log_sigma
  colnames(coef_trace) <- c(
    paste0("beta", 0:ncol(x)),   # intercept & covariates
    "alpha",
    "log_sigma"
  )
  
  for (j in 1:iter) {
    
    ## === E-step === ##
    mu_u1 <- cbind(1, x, z, 1) %*% matrix(t.coef1, ncol = 1)
    mu_u0 <- cbind(1, x, z, 0) %*% matrix(t.coef1, ncol = 1)
    alpha_weibull <- 1 / sigma
    
    lambda_u1 <- exp(-mu_u1)
    lambda_u0 <- exp(-mu_u0)
    
    lik_u1 <- (alpha_weibull * lambda_u1^alpha_weibull * t^(alpha_weibull - 1))^d *
      exp(-(lambda_u1 * t)^alpha_weibull)
    lik_u0 <- (alpha_weibull * lambda_u0^alpha_weibull * t^(alpha_weibull - 1))^d *
      exp(-(lambda_u0 * t)^alpha_weibull)
    
    logit_zu1 <- (1 - plogis(cbind(1, x, 1) %*% matrix(z.coef, ncol=1)))^(1 - z) *
      plogis(cbind(1, x, 1) %*% matrix(z.coef, ncol=1))^z
    logit_zu0 <- (1 - plogis(cbind(1, x, 0) %*% matrix(z.coef, ncol=1)))^(1 - z) *
      plogis(cbind(1, x, 0) %*% matrix(z.coef, ncol=1))^z
    
    ptzu1 <- lik_u1 * logit_zu1 * theta
    ptzu0 <- lik_u0 * logit_zu0 * (1 - theta)
    
    p <- ptzu1 / (ptzu1 + ptzu0)
    p[ptzu1 == 0 & ptzu0 == 0] <- 0
    
    ## === M-step: time-to-event === ##
    # use last params as init
    outcome_fit <- optim(par = init_params,
                         fn  = fn_outcome,
                         t   = t, d = d, x = x, z = z, p = p, gammat = gammat,
                         method = "BFGS",
                         control = list(maxit = 1000))
    
    params  <- outcome_fit$par
    t.coef1 <- c(params[1:(ncol(x)+2)], gammat)  # keep same design columns
    sigma   <- exp(params[ncol(x)+3])            # back-transform log_sigma
    init_params <- c(params[1:(ncol(x)+2)], log(sigma)) 
    
    # store the *updated* params for plotting
    coef_trace[j, ] <- init_params
    
    ## === M-step: treatment model === ##
    z.fit <- optim(par = z.coef[1:(nx+1)],
                   fn  = fn_treat,
                   x   = x, z = z, p = p, gammaz = gammaz)
    z.coef <- c(z.fit$par, gammaz)
    
    t.coef1[is.na(t.coef1)] <- 0
    z.coef[is.na(z.coef)]   <- 0
  }
  
  return(list(p = p, t.coef1 = t.coef1, z.coef = z.coef, sigma = sigma, trace = coef_trace))
  #all coefficients
  }
  
  
