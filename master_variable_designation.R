# master script of all variables in simulations

set.seed(965)
vecDiff = 2
xmin = 0.001
xmax = 100
pr_scenarios <- data.frame(h = c(
  0.00001,
  0.0001,
  0.001,
  0.01),
  b = c(2),
  pr_scenario = c("01", "02", "03", "04"))
# h = c(0.0001, 0.00001)
# b = c(1.5, 2)
n = c(500, 1000, 5000)
lambda = c(-1.5, -2, -2.5)

# pr_scenarios <- expand_grid(h = h,
#                             b = b) |>
#   mutate(pr_scenario = paste0("0", as.character(1:n())))

n_iter <- 500
