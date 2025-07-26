simulate_data <- function(n, tau, U_type = c("binary", "normal", "gamma", "binary+normal")) {
  
  # ---- Generate baseline covariates ----
  X1 <- rbinom(n, size = 1, prob = 0.4)
  X2 <- rbinom(n, size = 1, prob = 0.6)
  
  # ---- Generate unmeasured confounders ----
  U <- list()
  if (U_type == "binary") {
    U$U1 <- rbinom(n, size = 1, prob = 0.5)
  } else if (U_type == "normal") {
    U$U2 <- rnorm(n, mean = 0, sd = 1)
  } else if (U_type == "gamma") {
    U$U3 <- rgamma(n, shape = 0.5, scale = 0.5)
  } else if (U_type == "binary+normal") {
    U$U1 <- rbinom(n, size = 1, prob = 0.5)
    U$U2 <- rnorm(n, mean = 0, sd = 1)
  }
  
  # ---- Treatment assignment parameters ----
  alpha_0  <- 0.1
  alpha_x1 <- 0.3
  alpha_x2 <- 0.5
  alpha_u  <- -0.8
  
  # ---- Treatment linear predictor ----
  linpred_A <- alpha_0 + alpha_x1 * X1 + alpha_x2 * X2
  
  if (!is.null(U$U1)) linpred_A <- linpred_A + alpha_u * U$U1
  if (!is.null(U$U2)) linpred_A <- linpred_A + alpha_u * U$U2
  if (!is.null(U$U3)) linpred_A <- linpred_A + alpha_u * U$U3
  
  pi_A <- plogis(linpred_A)
  A    <- rbinom(n, 1, pi_A)
  
  # ---- Weibull AFT outcome params ----
  eta_intercept_a0 <- 0.2
  eta_intercept_a1 <- 0.7
  eta_x1 <- -0.1
  eta_x2 <- 0.4
  eta_u  <- -0.8
  
  sigma0 <- sigma1 <- 1/2
  alpha0 <- alpha1 <- 1/sigma0  # Weibull shape
  
  # ---- Linear predictors for potential outcomes ----
  lp_a0 <- eta_intercept_a0 + eta_x1*X1 + eta_x2*X2
  lp_a1 <- eta_intercept_a1 + eta_x1*X1 + eta_x2*X2
  
  if (!is.null(U$U1)) { lp_a0 <- lp_a0 + eta_u * U$U1; lp_a1 <- lp_a1 + eta_u * U$U1 }
  if (!is.null(U$U2)) { lp_a0 <- lp_a0 + eta_u * U$U2; lp_a1 <- lp_a1 + eta_u * U$U2 }
  if (!is.null(U$U3)) { lp_a0 <- lp_a0 + eta_u * U$U3; lp_a1 <- lp_a1 + eta_u * U$U3 }
  
  lambda0 <- exp(-lp_a0)
  lambda1 <- exp(-lp_a1)
  
  # ---- Potential event times ----
  D_a0 <- rweibull(n, shape = alpha0, scale = 1/lambda0)
  D_a1 <- rweibull(n, shape = alpha1, scale = 1/lambda1)
  
  # ---- Non-informative censoring ----
  C_a0 <- runif(n, 0.1, tau)
  C_a1 <- runif(n, 0.1, tau)
  
  # ---- Observed time & event indicator ----
  M <- ifelse(A == 1, pmin(tau, C_a1, D_a1),
              pmin(tau, C_a0, D_a0))
  delta <- ifelse(A == 1,
                  as.numeric(D_a1 <= pmin(tau, C_a1)),
                  as.numeric(D_a0 <= pmin(tau, C_a0)))
  
  # ---- Build the final data.frame robustly ----
  extra_cols <- if (length(U) > 0) U else NULL
  
  data_sim <<- as.data.frame(cbind(
    M = M, delta = delta, A = A, X1 = X1, X2 = X2,
    if (!is.null(extra_cols)) do.call(cbind, extra_cols)
  ))
  
  data_sim_com <<- as.data.frame(cbind(
    D_a0 = D_a0, D_a1 = D_a1, C_a0 = C_a0, C_a1 = C_a1,
    M = M, delta = delta, A = A, X1 = X1, X2 = X2,
    if (!is.null(extra_cols)) do.call(cbind, extra_cols)
  ))
  
  invisible(NULL) 
}
