# Morin sampling probabilities with fixed cutoffs. 
# morin recommends a generic cutoff of L = 10X Mesh size
# testing that here. 

# source("fixed_cut_morin_10xMesh.R")

# this script samples body sizes from bounded power law with known lambdas
# lambdas change across a hypothetical gradient
# biases data according to 4 Mesh sampling probabilities
# tests 1 set cutoff value for each Mesh size (L = 10X mesh opening size)
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

# M sizes (mm) = (0.125, 0.25, 0.5, 1.0)
# cutoff lengths (mm) = (1.25, 2.5, 5, 10)

# masses of cutoff lengths
sizeSpectra::lengthToMass(c(1.25,
                            2.5,
                            5,
                            10),
                          LWa = 0.0064,
                          LWb = 2.788)

cutoffs <- data.frame(
  cutoff = c(1.25,# mass = 0.0120
             2.5,# mass =  0.082
             5,# mass =    0.569
             10)# mass =   3.928
)

cutoffs <- cbind(cutoffs, morin_scenarios)
cutoffs_out <- cutoffs
cutoffs_out$cutoff_mass <- sizeSpectra::lengthToMass(
  cutoffs_out$cutoff, 
  LWa, 
  LWb
)
cutoffs_out$cutoff_length <- cutoffs_out$cutoff
saveRDS(cutoffs_out[,2:5], "simulation_results/Morin_10x_cutoff.RDS")


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
saveRDS(run, paste0("simulation_results/fixed_cut_morin_10xM_run_time_s_", Sys.Date(), ".rds"))

# results
saveRDS(results, "simulation_results/fixed_cut_morin_10xM_sim_run.rds")

# remove Rplots ####
file.remove(list.files(pattern = "Rplot*"))
