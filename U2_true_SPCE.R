#Numeric Iteration for computing the true causal effect, SPCE, survival probability causal effect
#pre-specify the true values of parameters from data generating process
eta_intercept <- 0.2#intercept for A=0
eta_intercept_a1<- 0.7#intercept for A=1
eta_x1 <- -0.1
eta_x2 <- 0.4
eta_u <- -0.8
sigma1 <- 0.5
sigma0 <- 0.5

#true SPCE at time=2
t_val <- 2

# Define the PDF of Normal distribution for U (mean=0, sd=1)
normal_pdf <- function(u) {
  return(dnorm(u, mean = 0, sd = 1))  # PDF of N(0,1)
}

# Survival function for Weibull AFT model with continuous U
weibull_aft_survival_treat1 <- function(t, x1, x2, u){
  lp <- eta_intercept_a1 + eta_x1 * x1 + eta_x2 * x2 + eta_u * u
  S_t <- exp(-(exp(-lp) * t)^(1/sigma1))
  return(S_t)
}

weibull_aft_survival_treat0 <- function(t, x1, x2, u){
  lp <- eta_intercept + eta_x1 * x1 + eta_x2 * x2 + eta_u * u
  S_t <- exp(-(exp(-lp) * t)^(1/sigma0))
  return(S_t)
}

# Function to integrate over U for treatment = 1 using integrate()
integrate_X_given_U_treat1 <- function(t) {
  x1_vals <- c(0, 1)
  x2_vals <- c(0, 1)
  
  # Integrate over U using the integrate function for each x1 and x2
  sum_val <- 0
  for (x1 in x1_vals) {
    for (x2 in x2_vals) {
      prob_x1 <- ifelse(x1 == 1, 0.4, 0.6)
      prob_x2 <- ifelse(x2 == 1, 0.6, 0.4)
      
      # Numerically integrate over U with normal distribution (mean=0, sd=1)
      integral_value <- integrate(function(u) {
        weibull_aft_survival_treat1(t, x1, x2, u) * normal_pdf(u)
      }, lower = -Inf, upper = Inf)$value  # Integrating over all U (-∞ to ∞)
      
      sum_val <- sum_val + integral_value * prob_x1 * prob_x2
    }
  }
  
  return(sum_val)
}

# Function to integrate over U for treatment = 0 using integrate()
integrate_X_given_U_treat0 <- function(t) {
  x1_vals <- c(0, 1)
  x2_vals <- c(0, 1)
  
  # Integrate over U using the integrate function for each x1 and x2
  sum_val <- 0
  for (x1 in x1_vals) {
    for (x2 in x2_vals) {
      prob_x1 <- ifelse(x1 == 1, 0.4, 0.6)
      prob_x2 <- ifelse(x2 == 1, 0.6, 0.4)
      
      # Numerically integrate over U with normal distribution (mean=0, sd=1)
      integral_value <- integrate(function(u) {
        weibull_aft_survival_treat0(t, x1, x2, u) * normal_pdf(u)
      }, lower = -Inf, upper = Inf)$value  # Integrating over all U (-∞ to ∞)
      
      sum_val <- sum_val + integral_value * prob_x1 * prob_x2
    }
  }
  
  return(sum_val)
}

# Marginal survival for treatment = 1 (averaging over continuous U)
marginal_survival_treat1 <- function(t) {
  return(integrate_X_given_U_treat1(t))
}

# Marginal survival for treatment = 0 (averaging over continuous U)
marginal_survival_treat0 <- function(t) {
  return(integrate_X_given_U_treat0(t))
}

# True Survival Probability Causal Effect (SPCE)
true_SPCE <- function(t) {
  marginal_survival_treat1(t) - marginal_survival_treat0(t)
}


marginal_survival_treat1(t_val) # marginal survival probability for treat=1 #0.4833
marginal_survival_treat0(t_val) # marginal survival probability for treat=0  #0.2941053
true_SPCE(t_val) #0.1891947