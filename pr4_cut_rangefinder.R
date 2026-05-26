# Determining cutoff

# source("pr4_cut_rangefinder.R")

# this script samples body sizes from bounded power law with known lambdas
# lambdas chnage across a hypothetical gradient
# biases data according to extreme sampling probabilities
# tests 4 set cutoff values
# estimates lambda with biased and censored data
# estimates relationship of change in lambda across gradient (beta)

library(parallel)
library(foreach)
library(doParallel)

# source the custom functions written for this simulation study
source("custom_functions.R")
source("master_variable_designation.R")

rep = n_iter

beta_groups <- data.frame(
  group = rep(LETTERS[1:4], each = 5),
  known_beta = rep(c(0, -0.1, -0.25, -0.5),
                   each = 5),
  known_lambda = c(-2, -2, -2, -2, -2, 
                   -1.9, -1.95, -2, -2.05, -2.1,
                   -1.75, -1.875, -2, -2.125, -2.25,
                   -1.5, -1.75, -2, -2.25, -2.5 ),
  env_gradient = c(-1,-0.5,0,0.5,1),
  xmin = xmin, 
  xmax = xmax,
  vecDiff = vecDiff)

# beta_groups
pr_scenarios <- pr_scenarios[c(1, 4),]
pr_scenarios

df <- tidyr::expand_grid(
  beta_groups,
  pr_scenarios,
  rep = 1:rep,
  n = 5000,
  cutoff = c(0.001, 0.01, 0.1, 1))

5 * #  lambdas
  2 * #  scenarios
  1 * #  n's
  500 * #  reps
  4  #  cutoffs
# 10 000 simulations ~ 79 seconds
# adding scenarios, and more betas
  
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
                                   df_out <- fixed_cut_lambda(
                                     n = df[j,]$n, 
                                     lambda = df[j,]$known_lambda,
                                     xmin = df[j,]$xmin, 
                                     xmax = df[j,]$xmax, 
                                     h = df[j,]$h, 
                                     b = df[j,]$b,
                                     cutoff = df[j,]$cutoff,
                                     vecDiff = df[j,]$vecDiff)
                                   
                                   df_out$pr_scenario <- df[j,]$pr_scenario
                                   df_out$rep <- df[j,]$rep
                                   df_out$group <- df[j,]$group
                                   df_out$known_beta <- df[j,]$known_beta
                                   df_out$env_gradient <- df[j,]$env_gradient
                                   
                                   df_out
                                 }
#results
stopCluster(cl = cluster)

end <- tictoc::toc()
run <- end$callback_msg
saveRDS(run, paste0("simulation_results/pr1_4_range_run_time_s_", Sys.Date(), ".rds"))

# results
saveRDS(results, "simulation_results/pr1_4_range_sim_run.rds")
