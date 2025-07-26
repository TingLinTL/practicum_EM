### --- Parameters ---
eta_intercept    <- 0.2
eta_intercept_a1 <- 0.7
eta_x1 <- -0.1
eta_x2 <- 0.4
eta_u  <- -0.8
sigma1 <- 0.5
sigma0 <- 0.5
t_val  <- 2

### --- Gamma(2,2) PDF ---
gamma22_pdf <- function(u) dgamma(u, shape = 0.5, scale = 0.5)

### --- Weibull AFT survival ---
weibull_aft_survival_treat1 <- function(t, x1, x2, u) {
  lp <- eta_intercept_a1 + eta_x1*x1 + eta_x2*x2 + eta_u*u
  exp(-(exp(-lp) * t)^(1/sigma1))
}

weibull_aft_survival_treat0 <- function(t, x1, x2, u) {
  lp <- eta_intercept    + eta_x1*x1 + eta_x2*x2 + eta_u*u
  exp(-(exp(-lp) * t)^(1/sigma0))
}

### --- Numeric quadrature over U using integrate() ---
integrate_over_U_treat1 <- function(t, x1, x2) {
  res <- integrate(function(u) weibull_aft_survival_treat1(t,x1,x2,u)*gamma22_pdf(u),
                   lower=0, upper=Inf, rel.tol=1e-8)
  res$value
}

integrate_over_U_treat0 <- function(t, x1, x2) {
  res <- integrate(function(u) weibull_aft_survival_treat0(t,x1,x2,u)*gamma22_pdf(u),
                   lower=0, upper=Inf, rel.tol=1e-8)
  res$value
}

### --- Marginalize over discrete X1,X2 ---
marginal_survival_treat1 <- function(t) {
  sum_val <- 0
  for (x1 in c(0,1)) {
    for (x2 in c(0,1)) {
      p_x1 <- ifelse(x1==1,0.4,0.6)
      p_x2 <- ifelse(x2==1,0.6,0.4)
      sum_val <- sum_val + p_x1*p_x2*integrate_over_U_treat1(t,x1,x2)
    }
  }
  sum_val
}

marginal_survival_treat0 <- function(t) {
  sum_val <- 0
  for (x1 in c(0,1)) {
    for (x2 in c(0,1)) {
      p_x1 <- ifelse(x1==1,0.4,0.6)
      p_x2 <- ifelse(x2==1,0.6,0.4)
      sum_val <- sum_val + p_x1*p_x2*integrate_over_U_treat0(t,x1,x2)
    }
  }
  sum_val
}

### --- True SPCE ---
true_SPCE <- function(t) marginal_survival_treat1(t) - marginal_survival_treat0(t)

S1_u3 <- marginal_survival_treat1(t_val)
S0_u3 <- marginal_survival_treat0(t_val)
SPCE_u3 <- true_SPCE(t_val)

S1_u3 #0.1908563 #gamma(0.5,0.5)  0.3902935
S0_u3 #0.04417639 #0.1185708
SPCE_u3 # 0.14668 #0.2717227