# Determining cutoff

# source("fixed_cut_morin_search_probability.R")

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

lambda_parameters <- data.frame(
  known_lambda = -1.9,
  xmin = xmin, 
  xmax = xmax,
  vecDiff = vecDiff,
  bias_level = "strong")

cutoff_5 <- data.frame(cutoff_mass = c(0.001726209,
                                    0.013878422,
                                    0.111821289,
                                    0.150001767,
                                    0.159636867,
                                    0.186210813,
                                    0.225830864,
                                    0.439920252, 
                                    1.088082505,
                                    9.068814540),
                    cutoff = c(0.625,# 10% 
                               1.32,# 50% 
                               2.79,# 90% 
                               3.1, #92.5%
                               3.17, #93%
                               3.35, #94%
                               3.59,# 95% 
                               4.56, # 97.5% 
                               6.31,# 99.0% 
                               13.5), # 99.9%
                    cutoff_probabilities = c(
                      "10%",
                      "50%",
                      "90%",
                      "92.5%",
                      "93%",
                      "94%",
                      "95%",
                      "97.5%",
                      "99.0%",
                      "99.9%"),
                    M = 0.5)

cutoff_25 <- data.frame(cutoff_mass = c(0.0002209791,
                                        0.0013881449,
                                        0.0047414852,
                                        0.0116584183,
                                        0.0127374171,
                                        0.0144726531,
                                        0.0166802617,
                                        0.0295020714,
                                        0.0676656537,
                                        0.4453207447),
                       cutoff = c(0.299,# 10% 
                                  0.578,# 50% 
                                  0.898,# 90% 
                                  1.24, #92.5%
                                  1.28, #93%
                                  1.34, #94%
                                  1.41,# 95% 
                                  1.73, # 97.5% 
                                  2.33,# 99.0% 
                                  4.58), # 99.9%
                       cutoff_probabilities = c(
                         "10%",
                         "50%",
                         "90%",
                         "92.5%",
                         "93%",
                         "94%",
                         "95%",
                         "97.5%",
                         "99.0%",
                         "99.9%"),
                       M = 0.25)

cutoff_df <- rbind(cutoff_5, cutoff_25)

# 10% = 0.625
# 50% = 1.32
# 90% = 2.79
# 95% = 3.59
# 97/5 = 4.56
# 99.0% = 6.31
# 99.9% = 13.5



df <- tidyr::expand_grid(
  lambda_parameters, 
  rep = 1:rep,
  n = 10000,
  # cutoffs are body lengths
  cutoff_df,
  LWa = 0.0064,
  LWb = 2.788)

1 * #  lambdas
  2 * #  scenarios
  1 * #  n's
  500 * #  reps
  10  #  cutoffs
# 10000 sims ~ 

  
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
                                   df_out <- morin_fixed_cut_lambda(
                                     n = df[j,]$n, 
                                     lambda = df[j,]$known_lambda,
                                     xmin = df[j,]$xmin, 
                                     xmax = df[j,]$xmax, 
                                     M = df[j,]$M,
                                     cutoff = df[j,]$cutoff,
                                     LWa = df[j,]$LWa,
                                     LWb = df[j,]$LWb,
                                     vecDiff = df[j,]$vecDiff)
                                   
                                   df_out$bias_level <- df[j,]$bias_level
                                   df_out$M <- df[j,]$M
                                   df_out$rep <- df[j,]$rep
                                   df_out$cutoff_mass <- df[j,]$cutoff_mass
                                   df_out$cutoff_probabilities <- df[j,]$cutoff_probabilities
                                   df_out$LWa = df[j,]$LWa
                                   df_out$LWb = df[j,]$LWb
                                   
                                   df_out
                                 }
#results
stopCluster(cl = cluster)

end <- tictoc::toc()
run <- end$callback_msg
saveRDS(run, paste0("simulation_results/fixed_cut_morin_search_run_time_s_", Sys.Date(), ".rds"))

# results
saveRDS(results, "simulation_results/fixed_cut_morin_search_sim_run.rds")

# remove Rplots ####
file.remove(list.files(pattern = "Rplot*"))
