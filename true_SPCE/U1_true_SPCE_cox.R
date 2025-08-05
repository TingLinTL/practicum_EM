beta_A <- -0.5
beta_x1 <- -0.1
beta_x2 <- 0.4
beta_u <- -0.8

t_val <- 2  # time point for SPCE

# S(t | A=1, x1, x2, u)
cox_survival_treat1 <- function(t, x1, x2, u) {
  linpred <- beta_A * 1 + beta_x1 * x1 + beta_x2 * x2 + beta_u * u
  S_t <- exp(-t * exp(linpred))
  return(S_t)
}

# S(t | A=0, x1, x2, u)
cox_survival_treat0 <- function(t, x1, x2, u) {
  linpred <- beta_A * 0 + beta_x1 * x1 + beta_x2 * x2 + beta_u * u
  S_t <- exp(-t * exp(linpred))
  return(S_t)
}


# Pr(X1=1) = 0.4, Pr(X2=1) = 0.6
integrate_X_given_U_treat1 <- function(t, u) {
  x1_vals <- c(0, 1)
  x2_vals <- c(0, 1)
  sum_val <- 0
  for (x1 in x1_vals) {
    for (x2 in x2_vals) {
      prob_x1 <- ifelse(x1 == 1, 0.4, 0.6)
      prob_x2 <- ifelse(x2 == 1, 0.6, 0.4)
      sum_val <- sum_val + cox_survival_treat1(t, x1, x2, u) * prob_x1 * prob_x2
    }
  }
  return(sum_val)
}

integrate_X_given_U_treat0 <- function(t, u) {
  x1_vals <- c(0, 1)
  x2_vals <- c(0, 1)
  sum_val <- 0
  for (x1 in x1_vals) {
    for (x2 in x2_vals) {
      prob_x1 <- ifelse(x1 == 1, 0.4, 0.6)
      prob_x2 <- ifelse(x2 == 1, 0.6, 0.4)
      sum_val <- sum_val + cox_survival_treat0(t, x1, x2, u) * prob_x1 * prob_x2
    }
  }
  return(sum_val)
}

marginal_survival_treat1 <- function(t) {
  0.5 * integrate_X_given_U_treat1(t, u = 0) +
    0.5 * integrate_X_given_U_treat1(t, u = 1)
}

marginal_survival_treat0 <- function(t) {
  0.5 * integrate_X_given_U_treat0(t, u = 0) +
    0.5 * integrate_X_given_U_treat0(t, u = 1)
}

true_SPCE_cox <- function(t) {
  marginal_survival_treat1(t) - marginal_survival_treat0(t)
}

marginal_survival_treat1(t_val)  # survival probability under treatment
marginal_survival_treat0(t_val)  # survival probability under control
true_SPCE_cox(t_val)             # true SPCE 0.1573018
