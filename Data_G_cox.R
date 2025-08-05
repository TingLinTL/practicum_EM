simulate_data_cox_h0_1 <- function(n, tau, U_type = c("binary", "normal", "gamma", "binary+normal")) {
  
  # --- Covariates ---
  X1 <- rbinom(n, 1, 0.4)
  X2 <- rbinom(n, 1, 0.6)
  
  # --- Unmeasured confounders ---
  U <- list()
  if (U_type == "binary") {
    U$U1 <- rbinom(n, 1, 0.5)
  } else if (U_type == "normal") {
    U$U2 <- rnorm(n)
  } else if (U_type == "gamma") {
    U$U3 <- rgamma(n, shape = 0.5, scale = 0.5)
  } else if (U_type == "binary+normal") {
    U$U1 <- rbinom(n, 1, 0.5)
    U$U2 <- rnorm(n)
  }
  
  # --- Treatment assignment ---
  alpha_0 <- 0.1
  alpha_x1 <- 0.3
  alpha_x2 <- 0.5
  alpha_u <- -0.8
  
  linpred_A <- alpha_0 + alpha_x1 * X1 + alpha_x2 * X2
  if (!is.null(U$U1)) linpred_A <- linpred_A + alpha_u * U$U1
  if (!is.null(U$U2)) linpred_A <- linpred_A + alpha_u * U$U2
  if (!is.null(U$U3)) linpred_A <- linpred_A + alpha_u * U$U3
  
  pi_A <- pnorm(linpred_A)  #  PROBIT link
  A <- rbinom(n, 1, pi_A)
  
  # --- Log hazard linear predictor ---
  beta_A <- -0.5
  beta_x1 <- -0.1
  beta_x2 <- 0.4
  beta_u <- -0.8
  
  linpred_a1 <- beta_A * 1 + beta_x1 * X1 + beta_x2 * X2
  linpred_a0 <- beta_A * 0 + beta_x1 * X1 + beta_x2 * X2
  if (!is.null(U$U1)) {
    linpred_a1 <- linpred_a1 + beta_u * U$U1
    linpred_a0 <- linpred_a0 + beta_u * U$U1
  }
  if (!is.null(U$U2)) {
    linpred_a1 <- linpred_a1 + beta_u * U$U2
    linpred_a0 <- linpred_a0 + beta_u * U$U2
  }
  if (!is.null(U$U3)) {
    linpred_a1 <- linpred_a1 + beta_u * U$U3
    linpred_a0 <- linpred_a0 + beta_u * U$U3
  }
  
  # --- Potential event times under A = 1 and A = 0 ---
  D_a1 <- -log(runif(n)) / exp(linpred_a1)
  D_a0 <- -log(runif(n)) / exp(linpred_a0)
  
  # --- Censoring times under A = 1 and A = 0 ---
  C_a1 <- runif(n, 0.1, tau)
  C_a0 <- runif(n, 0.1, tau)
  
  # --- Observed time and event indicator ---
  M <- ifelse(A == 1, pmin(tau, D_a1, C_a1),
              pmin(tau, D_a0, C_a0))
  
  delta <- ifelse(A == 1,
                  as.numeric(D_a1 <= pmin(tau, C_a1)),
                  as.numeric(D_a0 <= pmin(tau, C_a0)))
  
  # --- Output ---
  extra_cols <- if (length(U) > 0) U else NULL
  
  data_sim <<- as.data.frame(cbind(
    M = M, delta = delta, A = A, X1 = X1, X2 = X2,
    if (!is.null(extra_cols)) do.call(cbind, extra_cols)
  ))
  
  data_sim_com <<- as.data.frame(cbind(
    D_a0 = D_a0, D_a1 = D_a1,
    C_a0 = C_a0, C_a1 = C_a1,
    M = M, delta = delta, A = A, X1 = X1, X2 = X2,
    if (!is.null(extra_cols)) do.call(cbind, extra_cols)
  ))
  
  invisible(NULL)
}
