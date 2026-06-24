# master script of all variables in simulations

set.seed(965)
vecDiff = 2
xmin = 0.001
xmax = 10000
# pr_scenarios <- data.frame(h = c(
#   0.00001,
#   0.0001,
#   0.001,
#   0.01),
#   b = c(1.5),
#   pr_scenario = c("01", "02", "03", "04"))
# h = c(0.0001, 0.00001)
# b = c(1.5, 2)
n = 10000
lambda = c(-1.5, -2, -2.5)
n_iter <- 500

morin_scenarios <- data.frame(
  M = c(0.125, 0.25, 0.5, 1),
  bias_level = factor(
    c("minimal", "moderate", "strong", "extreme"),
    levels = c("minimal", "moderate", "strong", "extreme"))
)

beta_groups <- data.frame(
  group = rep(LETTERS[1:4], each = 5),
  known_beta = rep(c(0, -0.05, -0.1, -0.2),
                   each = 5),
  known_lambda = c(-1.9, -1.9, -1.9, -1.9, -1.9, 
                   -1.825, -1.875, -1.9, -1.925, -1.95,
                   -1.8, -1.85, -1.9, -1.95, -2.00,
                   -1.7, -1.8, -1.9, -2, -2.1),
  env_gradient = rep(c(-1, -0.5, 0, 0.5, 1), 4),
  xmin = xmin, 
  xmax = xmax,
  vecDiff = vecDiff)
