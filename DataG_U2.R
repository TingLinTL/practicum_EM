simulate_data_U2 <- function(n, tau) {
  
  ## ---- Generate baseline covariates & unmeasured confounder ----
  X1 <- rbinom(n, size = 1, prob = 0.4)   # X1 ~ Bern(0.4)
  X2 <- rbinom(n, size = 1, prob = 0.6)   # X2 ~ Bern(0.6)
  U2 <- rnorm(n, mean = 0, sd = 1)        # U2 ~ Normal(0,1)
  
  ## ---- Treatment assignment ----
  alpha_0  <- 0.1
  alpha_x1 <- 0.3
  alpha_x2 <- 0.5
  alpha_u  <- -0.8   # sensitivity parameter for U2
  
  # Treatment linear predictor now depends on continuous U2
  linpred_A <- alpha_0 + alpha_x1 * X1 + alpha_x2 * X2 + alpha_u * U2
  pi_A      <- 1 / (1 + exp(-linpred_A))
  A         <- rbinom(n, size = 1, prob = pi_A)
  
  ## ---- Potential outcome model params ----
  eta_intercept_a0 <- 0.2  # baseline for A=0
  eta_intercept_a1 <- 0.7  # baseline for A=1
  eta_x1 <- -0.1
  eta_x2 <- 0.4
  eta_u  <- -0.8            # sensitivity parameter for U2
  
  ## ---- Non-informative censoring ----
  C_a0 <- runif(n, 0.1, 5.5)
  C_a1 <- runif(n, 0.1, 5.5)
  
  ## ---- Weibull AFT potential event times ----
  sigma1 <- 1/2   # scale parameter for A=1
  sigma0 <- 1/2   # scale parameter for A=0
  alpha1 <- 1/sigma1  # Weibull shape for A=1
  alpha0 <- 1/sigma0  # Weibull shape for A=0
  
  # Linear predictors now include continuous U2
  lp_a1   <- eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U2
  lp_a0   <- eta_intercept_a0 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U2
  
  # Weibull scale parameters
  lambda1 <- exp(-lp_a1)
  lambda0 <- exp(-lp_a0)
  
  # Potential event times
  D_a1 <- rweibull(n, shape = alpha1, scale = 1/lambda1)
  D_a0 <- rweibull(n, shape = alpha0, scale = 1/lambda0)
  
  ## ---- Observed time & event indicator ----
  M <- ifelse(A == 1,
              pmin(tau, C_a1, D_a1),
              pmin(tau, C_a0, D_a0))
  
  delta <- ifelse(A == 1,
                  as.numeric(D_a1 <= pmin(tau, C_a1)),
                  as.numeric(D_a0 <= pmin(tau, C_a0)))
  
  ## ---- Observed dataset ----
  data_sim <<- data.frame(
    M     = M,
    delta = delta,
    A     = A,
    X1    = X1,
    X2    = X2,
    U2    = U2
  )
  
  ## ---- Full potential outcome dataset ----
  data_sim_com <<- data.frame(
    D_a0  = D_a0,
    D_a1  = D_a1,
    C_a0  = C_a0,
    C_a1  = C_a1,
    M     = M,
    delta = delta,
    A     = A,
    X1    = X1,
    X2    = X2,
    U2    = U2
  )
  
  invisible(NULL)  # no printing
}
