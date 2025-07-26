#Numeric Iteration for computing the true causal effect, SPCE, survival probability causal effect
#pre-specify the true values of parameters from data generating process
eta_intercept <- 0.2#intercept for A=0
eta_intercept_a1<- 0.7#intercept for A=1
eta_x1 <- -0.1
eta_x2 <- 0.4
eta_u1 <- -0.8
eta_u2<- -0.8
sigma1 <- 0.5
sigma0 <- 0.5

#true SPCE at time=2
t_val <- 2


# PDF of U2 ~ Normal(0,1)
normal_pdf <- function(u2) dnorm(u2, 0, 1)

# Survival function with U1 and U2
weibull_surv <- function(t, a, x1, x2, u1, u2) {
  lp <- (ifelse(a==1, eta_intercept_a1, eta_intercept)) +
    eta_x1*x1 + eta_x2*x2 +
    eta_u1*u1 + eta_u2*u2
  return(exp(-(exp(-lp)*t)^(1/(ifelse(a==1,sigma1,sigma0)))))
}

# First marginalize over U1 (binary)
surv_avg_over_U1 <- function(t, a, x1, x2, u2) {
  0.5 * weibull_surv(t, a, x1, x2, u1=0, u2) +
    0.5 * weibull_surv(t, a, x1, x2, u1=1, u2)
}

# Then integrate over U2
surv_avg_over_U2 <- function(t, a, x1, x2) {
  integrate(function(u2) surv_avg_over_U1(t, a, x1, x2, u2) * normal_pdf(u2),
            lower=-Inf, upper=Inf)$value
}

# Finally sum over X1,X2
marginal_survival <- function(t, a) {
  x1_vals <- c(0,1); x2_vals <- c(0,1)
  sum_val <- 0
  for (x1 in x1_vals) {
    for (x2 in x2_vals) {
      prob_x1 <- ifelse(x1==1,0.4,0.6)
      prob_x2 <- ifelse(x2==1,0.6,0.4)
      sum_val <- sum_val + prob_x1*prob_x2*surv_avg_over_U2(t,a,x1,x2)
    }
  }
  sum_val
}

# SPCE
true_SPCE <- function(t) marginal_survival(t,1) - marginal_survival(t,0)


marginal_survival(t_val,1) # marginal survival probability for treat=1 #0.3407907
marginal_survival(t_val,0) # marginal survival probability for treat=0  #0.1898964
true_SPCE(t_val) #0.1508942