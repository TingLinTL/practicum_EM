#Simulation
#source("DataGenerating.R")
source("Data_G.R")

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
n_sim <-50 #number of simulation
t_pred <- 2  #time point for estimates 
true_spce <- 0.1901985 #true SPCE at time=2 from numeric iteration for weibull aft model
n_boot <-100

spce_estimates_G <- numeric(n_sim)
spce_estimates_W <- numeric(n_sim)
spce_estimates_Sto_G <- numeric(n_sim)
spce_estimates_Sto_W <- numeric(n_sim)


for (i in 1:n_sim) {
  cat("Running simulation", i, "of", n_sim, "\n")
  
  #simulate_data(n,tau) #use data_sim 
  #simulate_data_U2(n, tau) #U2
  #simulate_data_U1_U2(n, tau) #U1&U2
  #simulate_data_U3(n, tau) #U3
  simulate_data(n, tau, U_type = "gamma") #binary, normal, gamma, binary+normal
  
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
  
  ## --------------------------------------- ##
  ## === BOOTSTRAP inside this iteration === ##
  ## --------------------------------------- ##
  boot_spce_EM_G <- boot_spce_EM_W <- boot_spce_stoEM_G <- boot_spce_stoEM_W <-numeric(n_boot)
  boot_SE_G <- boot_SE_W <- boot_stoSE_G <- boot_stoSE_W <- numeric(n_sim)
  boot_cover_G <- boot_cover_W <- boot_cover_sto_G <- boot_cover_sto_W <- numeric(n_sim)  
  
  for (b in 1:n_boot) {
    
    cat("Running boot", b, "of", n_boot, "\n")
    
    boot_idx <- sample(1:nrow(data_sim), size = nrow(data_sim), replace = TRUE)
    data_boot <- data_sim[boot_idx, ]
    x_boot    <- as.matrix(data_boot[, c("X1", "X2")])
    
    boot_em <- emU_surv(
      t      = data_boot$M,
      d      = data_boot$delta,
      z      = data_boot$A,
      x      = x_boot,
      gammat = -0.8,
      gammaz = -0.8
    )
    
    # compute SPCE for bootstrap sample G-computation
    boot_res_spce_G <- compute_SPCE_EM(
      res_em = boot_em,
      X      = x_boot,
      t0     = t_pred
    )
    
    boot_spce_EM_G[b] <- boot_res_spce_G$SPCE
  }
  
  ## bootstrap SE for this iteration
  boot_SE_G[i] <- sd(boot_spce_EM_G)
  
  ## bootstrap CI (percentile)
  ci_boot <- quantile(boot_spce_EM_G, probs = c(0.025, 0.975))
  
  ## check coverage
  boot_cover_G[i] <- (true_spce >= ci_boot[1] & true_spce <= ci_boot[2])
  
  
  
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


mean_spce_G <- mean(spce_estimates_G);mean_spce_G #0.1891319 #u2 0.3571589 #u1&u2 0.3346073
mean_spce_W <- mean(spce_estimates_W);mean_spce_W #0.1909884 #u2 0.3265647 #u1&u2 0.3006539
mean_spce_stoG <- mean(spce_estimates_Sto_G); mean_spce_stoG #0.1895913 #u2 0.3482394 #u1&u20.3253651
mean_spce_stoW <- mean(spce_estimates_Sto_W); mean_spce_stoW #0.2080809 #u2 0.4058067 #u1&u2 0.3736487

#===============================================================
## After all 200 iterations:
# Empirical SE across simulations
ESE_G <- sd(spce_estimates_G)

# Average bootstrap SE
ASE_G <- mean(boot_SE_G)

# Coverage Probability
CP_G  <- mean(boot_cover_G)


# ## RMSD function (Root Mean Squared Deviation)
# compute_RMSD <- function(estimates, true_value) {
#   sqrt(mean((estimates - true_value)^2))
# }
# 
# ## ESE (empirical SD across iterations)
# compute_ESE <- function(estimates) {
#   sd(estimates)  # sample SD across iterations
# }
# 
# #RMSD and ESE for G-computation
# RMSD_G <- compute_RMSD(spce_estimates_G, true_spce)
# ESE_G  <- compute_ESE(spce_estimates_G)
# 
# #RMSD and ESE for weighting
# RMSD_W <- compute_RMSD(spce_estimates_W, true_spce)
# ESE_W  <- compute_ESE(spce_estimates_W)
# 
# 
# RMSD_G ;ESE_G ;RMSD_W;ESE_W  
