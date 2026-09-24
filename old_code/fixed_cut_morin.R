# Morin sampling probabilities with fixed cutoffs. 
# originally tried 99% cutoffs
# after testing different cutoffs for M = 0.5, 95% sampling probability seems to be the sweet spot. 

# source("fixed_cut_morin.R")

# this script samples body sizes from bounded power law with known lambdas
# lambdas change across a hypothetical gradient
# biases data according to 4 Mesh sampling probabilities
# tests 1 set cutoff value for each Mesh size (~99.0% retention probability)
# estimates lambda with biased and censored data
# estimates relationship of change in lambda across gradient (beta)

library(parallel)
library(foreach)
library(doParallel)

# source the custom functions written for this simulation study
source("custom_functions.R")
source("master_variable_designation.R")

rep = n_iter # rep = 2

beta_groups 

# upon closer examination, 97.5% retention probability had better coverage in the 95% CIs
# lengths with ~99.0% retention probability
plogis(morin_ln_p(L = 0.705, M = 0.125))
plogis(morin_ln_p(L = 1.74, M = 0.25))
plogis(morin_ln_p(L = 4.58, M = 0.5))
plogis(morin_ln_p(L = 13.25, M = 1))

# masses with 99.0% RP
sizeSpectra::lengthToMass(c(0.705,
                            1.74,
                            4.58,
                            13.25),
                          LWa = 0.0064,
                          LWb = 2.788)


# 95% = 0.584, 1.41, 3.63, 10.1
# 97.5% = 0.0024, 0.0300, 0.4453, 8.6083
# 99% = 0.901, 2.30, 6.31, 19.15
cutoffs <- data.frame(
  cutoff = c(0.705,# mass = 0.00492
             1.74,# mass = 0.0653
             4.58,# mass = 1.083
             13.25)# mass = 25.281
)

cutoffs <- cbind(cutoffs, morin_scenarios)

df <- tidyr::expand_grid(
  beta_groups,
  cutoffs,
  rep = 1:rep,
  n = 10000,
  LWa = 0.0064,
  LWb = 2.788)

5 * #  lambdas
  4 * #  scenarios
  1 * #  n's
  4 * #  cutoffs
  500  #  reps
# 80 sets
# 40 000 simulations ~3 minutes
  
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
saveRDS(run, paste0("simulation_results/fixed_cut_morin_run_time_s_", Sys.Date(), ".rds"))

# results
saveRDS(results, "simulation_results/fixed_cut_morin_sim_run.rds")

# remove Rplots ####
file.remove(list.files(pattern = "Rplot*"))
