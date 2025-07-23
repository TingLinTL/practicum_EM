#Simulation
source("DataG_U1.R")

source("regular EM/em draft.R")
source("regular EM/SPCE_EM_G.R")
source("regular EM/SPCE_EM_W.R")

source("stoEM/SimulateU_surv.R")
source("stoEM/stoEM trace.R")
source("stoEM/stoEM draft.R")
source("stoEM/SPCE_stoEM_G.R")
source("stoEM/SPCE_stoEM_W.R")

library(survival)
set.seed(2000)
n <- 1000 #sample size 
tau <- 5.5 #maximum follow-up time
n_sim <-200 #number of simulation
t_pred <- 2  #time point for estimates 
true_spce <- 0.1901985 #true SPCE at time=2 from numeric iteration for weibull aft model

spce_estimates_G <- numeric(n_sim)
spce_estimates_W <- numeric(n_sim)
spce_estimates_Sto_G <- numeric(n_sim)
spce_estimates_Sto_W <- numeric(n_sim)


for (i in 1:n_sim) {
  cat("Running simulation", i, "of", n_sim, "\n")
  
  simulate_data(n,tau) #use data_sim 
  x_mat <- as.matrix(data_sim[, c("X1", "X2")]) # Prepare the covariate matrix
  
  ## ===  regular EM algorithm === ##
  result_em <- emU_surv(
    t      = data_sim$M,
    d      = data_sim$delta,
    z      = data_sim$A,
    x      = x_mat,
    gammat = -0.8,
    gammaz = -0.8
  )
  
  #========================================================
  # plot convergence of EM algorithm coefficients 
  # trace <- result_em$trace
  # matplot(
  #   1:nrow(trace), trace,
  #   type = "l", lty = 1,
  #   xlab = "EM Iteration",
  #   ylab = "Coefficient Value",
  #   main = "Convergence of EM Algorithm Coefficients"
  # )
  # legend("topright", legend = colnames(trace), col = 1:ncol(trace), lty = 1)
  #============================================================================
  
  #compute spce by g-computaion for regular EM algorithm
  res_spce_G <- compute_SPCE_EM(
    res_em = result_em,
    X      = x_mat,
    t0     = t_pred
  )
  
  #compute spce by weighting regular EM algorithm
  res_spce_weight <- compute_IPW_SPCE_from_EM(
    res_em = result_em,
    X      = x_mat,
    Z      = data_sim$A,
    M      = data_sim$M,
    delta  = data_sim$delta,
    t0     = t_pred
  )
  
  
  ## ===  Stochastic EM algorithm === ##
  
  #============================================================================
  #Convergence plot
  # stoEM_trace <-surv_stoEM_trace(t= data_sim$M,
  #                       d= data_sim$delta,
  #                       Z= data_sim$A,
  #                       X= x_mat,
  #                       zetat=-0.8, zetaz=-0.8, B=100, theta=0.5)
  # 
  # matplot( stoEM_trace, type = "l", lty = 1,
  #         main = "Coefficient Trace (Initial + Updates)",
  #         xlab = "Iteration", ylab = "Coefficient")
  # legend("topright", legend = colnames( stoEM_trace), col = 1:ncol(stoEM_trace), lty = 1)
  # ====================================================================================
  
  stoEM <- surv_stoEM(t= data_sim$M,
                      d= data_sim$delta,
                      Z= data_sim$A,
                      X= x_mat,
                      zetat=-0.8, zetaz=-0.8, B=100, theta=0.5)
    
    
  #stochastic EM for SPCE by G-computation
  stoEM_SPCE_G <- compute_SPCE_G(beta_final=stoEM$beta, 
                                 sigma_final=stoEM$sigma, 
                                 X= x_mat,t0=t_pred)
  
  #stochastic EM for SPCE by weighting
  stoEM_SPCE_W <- surv_stoEM_ipw(beta_final=stoEM$beta, 
                                 sigma_final=stoEM$sigma, 
                                 t= data_sim$M,
                                 d= data_sim$delta,
                                 Z= data_sim$A,
                                 X= x_mat,
                                 t0=t_pred,
                                 stabilize=TRUE)

  
  
  
  
  
  #store values 
  
  spce_estimates_G[i] <- res_spce_G$SPCE
  spce_estimates_W[i] <- res_spce_weight$SPCE
  
  spce_estimates_Sto_G[i] <- stoEM_SPCE_G$SPCE
  spce_estimates_Sto_W[i]<- stoEM_SPCE_W$SPCE
  
}

mean_spce_G <- mean(spce_estimates_G);mean_spce_G #0.1908726
mean_spce_W <- mean(spce_estimates_W);mean_spce_W #0.193349
mean_spce_stoG <-  mean(spce_estimates_Sto_G); mean_spce_stoG #0.1895913
mean_spce_stoW <-  mean(spce_estimates_Sto_W); mean_spce_stoW #0.2080809

#500iterations