# undersampling gradient simulation

# this script samples body sizes from bounded power law with known lambda
# lambdas chnage across a hypothetical gradient
# biases data according to 2 sampling probabilities (minimal and extreme)
# estimates lambda with biased and censored data
# estimates relationship of change in lambda across gradient (beta)

library(parallel)
library(foreach)
library(doParallel)

# source the custom functions written for this simulation study
source("custom_functions.R")

rep = 500

beta_groups <- data.frame(
  group = rep(LETTERS[1:4], each = 3),
  known_beta = rep(c(0, -0.1, -0.25, -0.5),
                   each = 3),
  known_lambda = c(-2, -2, -2, 
                   -1.9, -2, -2.1,
                   -1.75, -2, -2.25,
                   -1.5, -2, -2.5 ),
  env_gradient = rep(c(-1, 0, 1), 4),
  xmin = 0.001, 
  xmax = 100,
  vecDiff = 2)

# beta_groups

pr_scenarios <- data.frame(h = c(
  0.00001,
  0.01),
  b = c(2),
  pr_scenario = c("01", "04"))

df <- tidyr::expand_grid(beta_groups,
                  pr_scenarios,
                  rep = 1:rep,
                  n = c(500, 1000, 5000))
# 4 betas
# 3 lambdas
# 2 scenarios
# 3 n's
# 2 reps
4*3*2*3*2

# set up parallel processing
#cores <- detectCores()-1 # when running on its own
cores <- 6 # when running with other simulations

cluster <- makeCluster(cores)
registerDoParallel(cluster)


{tictoc::tic()
results <- foreach(j = 1:nrow(df), 
                   .combine = rbind,
                   .errorhandling = "remove",
                   .packages = c("sizeSpectra",
                                 "tidyverse",
                                 "poweRlaw")) %dopar%{
                                   # sets a new seed for each parallel instance
                                   # but makes this reproducible for re-running simulation
                                   set.seed(j+1130)
                                   df_out <- plot_sub_lambda(
                                     n = df[j,]$n, 
                                     lambda = df[j,]$known_lambda,
                                     xmin = df[j,]$xmin, 
                                     xmax = df[j,]$xmax, 
                                     h = df[j,]$h, 
                                     b = df[j,]$b,
                                     vecDiff = df[j,]$vecDiff,
                                     plot = FALSE)
                                   
                                   df_out$pr_scenario <- df[j,]$pr_scenario
                                   df_out$rep <- df[j,]$rep
                                   df_out$group <- df[j,]$group
                                   df_out$known_beta <- df[j,]$known_beta
                                   df_out$env_gradient <- df[j,]$env_gradient
                                   
                                   df_out
                                 }
#results
tictoc::toc()}
stopCluster(cl = cluster)
# 24 rows = 15.77 seconds
# 240 rows = 156 seconds
# 144 rows = 65
# results


saveRDS(results, "simulation_results/gradient_sim_run.rds")
