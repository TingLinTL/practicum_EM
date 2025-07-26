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

#True SPCE
#library(cubature)
#log(T)=mu + sigma * error, error~Gumble (0,1), mu=linear prediction
#t_i* ~ weibull (alpha=1/sigma, lambda = e^-{X*beta}), lamda is called scale
#CDF F(t) = 1 -exp(-(lambda*t)^alpha)
#survival probability = 1-F(t) = exp(-(lambda*t)^alpha) = exp(-(exp(-X*beta)*t)^(1/sigma)) 

# Survival function for weibull aft model A = 1 
weibull_aft_survival_treat1 <- function(t, x1, x2, u){
  lp <- eta_intercept_a1 + eta_x1 * x1 + eta_x2 * x2 + eta_u * u
  # z<- (log(t)-lp)/sigma
  # S_t <- exp(-exp(z))
  S_t <- exp(-(exp(-lp) * t)^(1/sigma1))
  return(S_t)
}
# Survival function for weibull aft model A = 0
weibull_aft_survival_treat0 <- function(t, x1, x2, u){
  lp <- eta_intercept + eta_x1 * x1 + eta_x2 * x2 + eta_u * u
  # z<- (log(t)-lp)/sigma
  # S_t <- exp(-exp(z))
  S_t <- exp(-(exp(-lp) * t)^(1/sigma0))
  return(S_t)
}


# X1 ~ Bern(0.4), X2 ~ Bern(0.6)
#both X1 and X2 are binary covariate

# Integrate (sum) over X1 and X2 for given U, treatment = 1
integrate_X_given_U_treat1 <- function(t, u) {
  x1_vals <- c(0, 1)
  x2_vals <- c(0, 1)
  sum_val <- 0
  for (x1 in x1_vals) {
    for (x2 in x2_vals) {
      prob_x1 <- ifelse(x1 == 1, 0.4, 0.6)
      prob_x2 <- ifelse(x2 == 1, 0.6, 0.4)
      sum_val <- sum_val + weibull_aft_survival_treat1(t, x1, x2, u) * prob_x1 * prob_x2
    }
  }
  return(sum_val)
}

# Integrate (sum) over X1 and X2 for given U, treatment = 0
integrate_X_given_U_treat0 <- function(t, u) {
  x1_vals <- c(0, 1)
  x2_vals <- c(0, 1)
  sum_val <- 0
  for (x1 in x1_vals) {
    for (x2 in x2_vals) {
      prob_x1 <- ifelse(x1 == 1, 0.4, 0.6)
      prob_x2 <- ifelse(x2 == 1, 0.6, 0.4)
      sum_val <- sum_val + weibull_aft_survival_treat0(t, x1, x2, u) * prob_x1 * prob_x2
    }
  }
  return(sum_val)
}

# Marginal survival for treatment = 1 (averaging over U)
marginal_survival_treat1 <- function(t) {
  0.5 * integrate_X_given_U_treat1(t, u = 0) +
    0.5 * integrate_X_given_U_treat1(t, u = 1)
}

# Marginal survival for treatment = 0 (averaging over U)
marginal_survival_treat0 <- function(t) {
  0.5 * integrate_X_given_U_treat0(t, u = 0) +
    0.5 * integrate_X_given_U_treat0(t, u = 1)
}


# True Survival Probability Causal Effect (SPCE)
true_SPCE <- function(t) {
  marginal_survival_treat1(t) - marginal_survival_treat0(t)
}



marginal_survival_treat1(t_val) # marginal survival probability for treat=1 #0.2834901
marginal_survival_treat0(t_val) # marginal survival probability for treat=0  #0.09329167
true_SPCE(t_val) #0.1901985