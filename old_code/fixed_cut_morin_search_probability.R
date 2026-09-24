# Determining cutoff

# source("fixed_cut_morin_search_probability.R")

# this script samples body sizes from bounded power law with known lambdas
# biases data according to 4 retention probabilities
# tests 11 cutoff values

library(parallel)
library(foreach)
library(doParallel)

# source the custom functions written for this simulation study
source("custom_functions.R")
source("master_variable_designation.R")

rep = 1000
LWa = 0.0064
LWb = 2.788

lambda_parameters <- data.frame(
  known_lambda = -1.9,
  xmin = xmin, 
  xmax = xmax,
  vecDiff = vecDiff,
  bias_level = "strong")

cutoff125 = c(
  .1485,# 10% 
  0.266,# 50% 
  0.477,# 90% 
  0.518, #92.5%
  0.5282, #93%
  0.552, #94%
  0.581,# 95% 
  0.703, # 97.5% 
  0.91,# 99.0% 
  1.085, # 99.5%
  1.25) # 10X mesh size ~99.7 
plogis(morin_ln_p(L = cutoff125,
                  M = 0.125))
cutoff_125 <- data.frame(
  cutoff = cutoff125, 
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
    "99.5%",
    "10xM"),
  M = 0.125)
cutoff_125$cutoff_mass <- sizeSpectra::lengthToMass(
  cutoff_125$cutoff,
  LWa = LWa, 
  LWb = LWb
)

cutoff25 = c(0.305,# 10% 
           0.578,# 50% 
           0.885,# 90% 
           1.24, #92.5%
           1.28, #93%
           1.34, #94%
           1.41,# 95% 
           1.740, # 97.5% 
           2.3,# 99.0% 
           2.5, # 10 x Mesh ~99.2
           2.85) # 99.5%
plogis(morin_ln_p(L = cutoff25,
                  M = 0.25))


cutoff_25 <- data.frame(
  cutoff = cutoff25, 
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
    "10xM",
    "99.5%"),
  M = 0.25)
cutoff_25$cutoff_mass <- sizeSpectra::lengthToMass(
  cutoff_25$cutoff,
  LWa = LWa, 
  LWb = LWb
)

cutoff5 = c(0.625,# 10% 
            1.32,# 50% 
            2.79,# 90% 
            3.1, #92.5%
            3.18, #93%
            3.36, #94%
            3.59,# 95% 
            4.58, # 97.5% 
            5, # 10X mesh size ~98
            6.31,# 99.0% 
            8) # 99.5%
plogis(morin_ln_p(L = cutoff5,
                  M = 0.5))

cutoff_5 <- data.frame(
  cutoff = cutoff5, 
  cutoff_probabilities = c(
    "10%",
    "50%",
    "90%",
    "92.5%",
    "93%",
    "94%",
    "95%",
    "97.5%",
    "10xM",
    "99.0%",
    "99.5%"),
  M = 0.5)
cutoff_5$cutoff_mass <- sizeSpectra::lengthToMass(
  cutoff_5$cutoff,
  LWa = LWa, 
  LWb = LWb
)

cutoff1 = c(1.291,# 10% 
            3.09,# 50% 
            7.4,# 90% 
            8.37, #92.5%
            8.62, #93%
            9.2, #94%
            9.94,# 95% 
            10, # 10X mesh size ~95.07
            13.23, # 97.5% 
            19.15,# 99.0% 
            25.25) # 99.5%
plogis(morin_ln_p(L = cutoff1,
                  M = 1))

cutoff_1 <- data.frame(
  cutoff = cutoff1, 
  cutoff_probabilities = c(
    "10%",
    "50%",
    "90%",
    "92.5%",
    "93%",
    "94%",
    "95%",
    "10xM",
    "97.5%",
    "99.0%",
    "99.5%"),
  M = 1)
cutoff_1$cutoff_mass <- sizeSpectra::lengthToMass(
  cutoff_1$cutoff,
  LWa = LWa, 
  LWb = LWb
)


cutoff_df <- rbind(cutoff_125,
                   cutoff_25,
                   cutoff_5,
                   cutoff_1)

saveRDS(cutoff_df, 
        "simulation_results/morin_fixed_cutoff_search.rds")

df <- tidyr::expand_grid(
  lambda_parameters, 
  rep = 1:rep,
  n = 10000,
  # cutoffs are body lengths
  cutoff_df,
  LWa,
  LWb)

1 * #  lambdas
  2 * #  scenarios
  1 * #  n's
  10000 * #  reps
  11  #  cutoffs
# 22000 sims ~ 4 minutes

  
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
