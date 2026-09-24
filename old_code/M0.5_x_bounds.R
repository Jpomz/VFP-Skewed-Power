# Investigating if %CIs with true value for M = 0.5 improves when shifting x bounds to be "bigger"
# i.e., not so much of the data "mass" is in the area with low retention probability 

# source("M0.5_x_bounds.R")

# this script samples body sizes from bounded power law with known lambdas
# This is different from the main text in that only the M = 0.5 bias is used
# additionally, the range of body sizes is shifted to start and end at larger values
# Estimates are corrected using the inverse Morin correction i.e., weight = 1 / retention probability

library(parallel)
library(foreach)
library(doParallel)

# source the custom functions written for this simulation study
source("custom_functions.R")
set.seed(965)
vecDiff = 2
xmin = 0.01
xmax = 100000
n = 10000
n_iter <- 500

lambda <- data.frame(
  known_lambda = c(-1.9, -2, -2.1),
  xmin = xmin, 
  xmax = xmax,
  vecDiff = vecDiff)

rep = n_iter # rep = 2

morin_scenarios <- data.frame(
  M = c(0.5),
  bias_level = factor(c("strong"),
    levels = c("strong"))
)


df <- tidyr::expand_grid(
  lambda,
  morin_scenarios,
  rep = 1:rep,
  n = 10000,
  LWa = 0.0064,
  LWb = 2.788)

3 * #  lambdas
  1 * #  scenarios
  1 * #  n's
  500  #  reps
# 1500 simulations 
  
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
                                   df_out <- morin_correction(
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
saveRDS(run, paste0("simulation_results/M0.5_x_bounds_run_time_s_", Sys.Date(), ".rds"))

# results
saveRDS(results, "simulation_results/M0.5_x_bounds_sim_run.rds")

# remove Rplots ####
file.remove(list.files(pattern = "Rplot*"))
