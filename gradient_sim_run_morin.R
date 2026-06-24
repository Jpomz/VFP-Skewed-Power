# undersampling gradient simulation

# source("gradient_sim_run_morin.R")

# this script samples body sizes from bounded power law with known lambda
# lambdas change across a hypothetical gradient
# biases data according to 4 mesh sizes from Morin et al. 2004
# estimates lambda with biased and censored data
# estimates relationship of change in lambda across gradient (beta)

library(parallel)
library(foreach)
library(doParallel)

# source the custom functions written for this simulation study
source("custom_functions.R")
source("master_variable_designation.R")

rep = n_iter

# beta_groups <- data.frame(
#   group = rep(LETTERS[1:4], each = 5),
#   known_beta = rep(c(0, -0.1, -0.25, -0.5),
#                    each = 5),
#   known_lambda = c(-1.5, -1.5, -1.5, -1.5, -1.5, 
#                    -1.4, -1.45, -1.5, -1.55, -1.6,
#                    -1.25, -1.375, -1.5, -1.625, -1.75,
#                    -1., -1.25, -1.5, -1.75, -2),
#   env_gradient = rep(c(-1, -0.5, 0, 0.5, 1), 4),
#   xmin = xmin, 
#   xmax = xmax,
#   vecDiff = vecDiff)

beta_groups
morin_scenarios

df <- tidyr::expand_grid(beta_groups,
                  morin_scenarios,
                  rep = 1:rep,
                  n = n,
                  LWa = 0.0064,
                  LWb = 2.788)
# 4 betas
# 5 lambdas
# 4 scenarios
# 1 n
# 500 reps
4*5*4*1*500
# 40 000 sims ~

# set up parallel processing
cores <- detectCores()-1 # when running on its own
#cores <- 7 # when running with other simulations

cluster <- makeCluster(cores)
registerDoParallel(cluster)


tictoc::tic()
results <- foreach(j = 1:nrow(df), 
                   .combine = rbind,
                   .errorhandling = "remove",
                   .packages = c("sizeSpectra",
                                 "tidyverse",
                                 "poweRlaw")) %dopar%{
                                   # sets a new seed for each parallel instance
                                   # but makes this reproducible for re-running simulation
                                   set.seed(j+1130)
                                   df_out <- morin_bias_lambda(
                                     n = df[j,]$n, 
                                     lambda = df[j,]$known_lambda,
                                     xmin = df[j,]$xmin, 
                                     xmax = df[j,]$xmax, 
                                     M = df[j,]$M,
                                     LWa = df[j,]$LWa,
                                     LWb = df[j,]$LWb,
                                     vecDiff = df[j,]$vecDiff)
                                   
                                   df_out$bias_level <- df[j,]$bias_level
                                   df_out$M <- df[j,]$M
                                   df_out$rep <- df[j,]$rep
                                   df_out$group <- df[j,]$group
                                   df_out$known_beta <- df[j,]$known_beta
                                   df_out$env_gradient <- df[j,]$env_gradient
                                   df_out$LWa = df[j,]$LWa
                                   df_out$LWb = df[j,]$LWb
                                   
                                   df_out
                                 }
#results
stopCluster(cl = cluster)
end <- tictoc::toc()
run <- end$callback_msg

saveRDS(run, paste0("simulation_results/gradient_run_morin_time_s_", Sys.Date(), ".rds"))

# results
saveRDS(results, "simulation_results/gradient_sim_run_morin.rds")

# remove Rplots ####
file.remove(list.files(pattern = "Rplot*"))
