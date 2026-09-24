# sfs 2026
library(tidyverse)
library(tidybayes)
library(sizeSpectra)
library(poweRlaw)

# load data ####
# undersampling simulation ####
sim <- as_tibble(readRDS("simulation_results/undersampling_sim_run.rds"))

# sampling scenario plots ####
sample_pr <- function(x, h, b){
  pr = 1 / (1 + (h / x**b))
  return(pr)
}

# SI figures for sampling probability with simulation values

set.seed(2112)
pr_dat <- expand_grid(
  h = c(0.00001,
        0.0001,
        0.001,
        0.01),
  b = c(2),
  x = c(exp(seq(log(0.001), log(10), length.out = 1000)))) |>
  mutate(pr = sample_pr(x = x, h = h, b = b),
         scenario = case_when(
           h == 0.00001 & b == 2 ~ "minimal",
           h == 0.0001 & b == 2 ~ "moderate",
           h == 0.001 & b == 2 ~ "strong",
           h == 0.01 & b == 2 ~ "extreme",
         ))

## 4 scenarios on one plot
pr_dat |>
  mutate(scenario = factor(scenario, 
                    levels = c("minimal", "moderate", "strong", "extreme"))) |>
  ggplot(aes(x = x, 
             y = pr, 
             color = scenario)) +
  geom_line(linewidth = 1.5) +
  geom_vline(aes(xintercept = 0.01), linetype = "dashed") +
  geom_hline(aes(yintercept = 0.9), linetype = "dashed") +
  scale_x_log10() +
  theme_bw(base_size = 18) +
  labs(title = "Sampling probability as a function of body mass",
       x =expression(Log[10]~dry~mass),
       y = "Sampling probability") +
  scale_colour_manual(values = c("darkorchid1",
                                 "darkorchid2",
                                 "darkorchid3",
                                 "darkorchid4")) +
  NULL
ggsave("plots/sample_pr_4_scenarios.png",
       units = "in",
       height = 4,
       width = 8)

## 2 scenarios on one plot
pr_dat |>
  filter(scenario == "minimal" |
           scenario == "strong")|>
  mutate(scenario = factor(scenario, 
                           levels = c(
                             "minimal", 
                             "strong"))) |>
  ggplot(aes(x = x, 
             y = pr, 
             color = scenario)) +
  geom_line(linewidth = 1.5) +
  geom_vline(aes(xintercept = 0.01), linetype = "dashed") +
  geom_hline(aes(yintercept = 0.9), linetype = "dashed") +
  scale_x_log10() +
  theme_bw(base_size = 18) +
  labs(title = "Sampling probability as a function of body mass",
       x =expression(Log[10]~dry~mass),
       y = "Sampling probability") +
  scale_colour_manual(values = c("darkorchid1",
                                 #"darkorchid2",
                                 "darkorchid3"#,
                                 #"darkorchid4"
  )) +
  NULL
ggsave("plots/sample_pr_2_scenarios.png",
       units = "in",
       height = 4,
       width = 8)

pr_dat |>
  filter(scenario == "minimal")|>
  mutate(scenario = factor(scenario, 
                           levels = c(
                             "minimal"))) |>
  ggplot(aes(x = x, 
             y = pr, 
             color = scenario)) +
  geom_line(linewidth = 1.5) +
  geom_vline(aes(xintercept = 0.01), linetype = "dashed") +
  geom_hline(aes(yintercept = 0.9), linetype = "dashed") +
  scale_x_log10() +
  theme_bw(base_size = 18) +
  labs(title = "Sampling probability as a function of body mass",
       x =expression(Log[10]~dry~mass),
       y = "Sampling probability") +
  scale_colour_manual(values = c("darkorchid1")) +
  NULL
ggsave("plots/sample_pr_min_scenarios.png",
       units = "in",
       height = 4,
       width = 8)

pr_dat |>
  filter(scenario == "strong")|>
  mutate(scenario = factor(scenario, 
                           levels = c(
                             "strong"))) |>
  ggplot(aes(x = x, 
             y = pr, 
             color = scenario)) +
  geom_line(linewidth = 1.5) +
  geom_vline(aes(xintercept = 0.01), linetype = "dashed") +
  geom_hline(aes(yintercept = 0.9), linetype = "dashed") +
  scale_x_log10() +
  theme_bw(base_size = 18) +
  labs(title = "Sampling probability as a function of body mass",
       x =expression(Log[10]~dry~mass),
       y = "Sampling probability") +
  scale_colour_manual(values = c("darkorchid1")) +
  NULL
ggsave("plots/sample_pr_str_scenarios.png",
       units = "in",
       height = 4,
       width = 8)



##=======##
# simulate x from bounded power law
set.seed(2015)
xmin = 0.001
xmax = 100
original_dat <- tibble(
  x = rPLB(5000,
           b = -2,
           xmin = xmin,
           xmax = xmax))

# expand grid to include each variable combination with each body size
sample_dat <- expand_grid(
  original_dat, 
  h = c(0.00001, 0.001),
  b = 2) |>
  mutate(scenario = 
           case_when(
             h == 0.00001 & b == 2 ~ "minimal",
             h == 0.001 & b == 2 ~ "strong"))

# add sample probabilities as a function of x (body mass)
sampled_dat <- sample_dat |> 
  mutate(pr = sample_pr(
    x, 
    h = h, 
    b = b
  ), 
  sampled = rbinom(n(), 1, prob = pr),
  fill = case_when(sampled == 1 ~ "sampled", 
                   .default = "not sampled"))

# add "fill = original" to simulated data
original_x <- original_dat |>
  mutate(fill = "original")

# combine "original and sampled data
both_dats <- left_join(sampled_dat, original_x, by = "x") |>
  pivot_longer(fill.x:fill.y)


# plot matching new write up terms and colors
# estimate x_min for this scenario
sampled_x_min <- both_dats |>
  filter(value == "sampled",
         scenario == "minimal") |>
  pull(x)
x_power_min <- conpl$new(sampled_x_min)
x_xmin_min <- estimate_xmin(x_power_min)$xmin
x_xmin_min
# strong
sampled_x_str <- both_dats |>
  filter(value == "sampled",
         scenario == "strong") |>
  pull(x)
x_power_str <- conpl$new(sampled_x_str)
x_xmin_str <- estimate_xmin(x_power_str)$xmin
x_xmin_str

both_dats |>
  mutate(x_min = case_when(
    scenario == "minimal" ~ 0.0108,
    scenario == "strong" ~ 0.066
  )) |>
  filter(value!= "not sampled") |>
  mutate(value = case_when(
    value == "sampled" ~ "Biased", 
    .default = "Original"
  ),
  value = factor(value, 
                 levels = c("Original", "Biased"))) |>
  ggplot(aes(x = x, 
             fill = value)) +
  geom_histogram(binwidth = 0.05, 
                 position = "dodge",
                 alpha = 0.75) +
  scale_x_log10() +
  facet_wrap(~scenario) +
  geom_segment(aes(x = x_min,
                   y = 400,
                   xend = x_min,
                   yend = 200),
               arrow = arrow(length = unit(0.5, "cm")),
               linewidth = 2,
               color = "dodgerblue") +
  scale_fill_manual(values = c("black", "#FF1984")) + 
  theme_bw(base_size = 18) +
  guides(fill = guide_legend(title= "Data")) +
  labs(x =expression(Log[10]~dry~mass)) 
  
ggsave("plots/orig_bias_xmin_sfs.png",
       units = "in",
       height = 4,
       width = 8)

# minimal bias estimate
bias_min <- calcLike(
  negLL.fn = negLL.PLB,
  x = sampled_x_min,
  xmin = min(sampled_x_min), 
  xmax = max(sampled_x_min), 
  n = length(sampled_x_min), 
  sumlogx = sum(log(sampled_x_min)), 
  p = -1.5,
  suppress.warnings = TRUE,
  vecDiff = 2)

# strong bias estimate
bias_str <- calcLike(
  negLL.fn = negLL.PLB,
  x = sampled_x_str,
  xmin = min(sampled_x_str), 
  xmax = max(sampled_x_str), 
  n = length(sampled_x_str), 
  sumlogx = sum(log(sampled_x_str)), 
  p = -1.5,
  suppress.warnings = TRUE,
  vecDiff = 2)

# estimates with correction
x_min_correct <- sampled_x_min[sampled_x_min>=x_xmin_min]
lambda_min_corr <- calcLike(
  negLL.fn = negLL.PLB,
  x = x_min_correct,
  xmin = min(x_min_correct), 
  xmax = max(x_min_correct), 
  n = length(x_min_correct), 
  sumlogx = sum(log(x_min_correct)), 
  p = -1.5,
  suppress.warnings = TRUE,
  vecDiff = 2)

x_str_correct <- sampled_x_str[sampled_x_str>=x_xmin_str]
lambda_str_corr <- calcLike(
  negLL.fn = negLL.PLB,
  x = x_str_correct,
  xmin = min(x_str_correct), 
  xmax = max(x_str_correct), 
  n = length(x_str_correct), 
  sumlogx = sum(log(x_str_correct)), 
  p = -1.5,
  suppress.warnings = TRUE,
  vecDiff = 2)
lambda_str_corr


tibble(bias_est = c(bias_min$MLE,
                    bias_str$MLE),
       corr_est = c(lambda_min_corr$MLE,
                    lambda_str_corr$MLE))

# sim summary ####
# filter out a few examples
# minimal and extreme sampling
# only one N
distinct(sim, original_n)
distinct(sim, original_n, pr_scenario, known_lambda)
dat <- sim |>
  filter(original_n == 5000,
         pr_scenario == "01" | pr_scenario == "03")

plot_dat <- dat |>
  select(-lambda_under_lo, -lambda_under_hi) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed,
         `Bias level` = pr_scenario) |>
  mutate(`Bias level` = case_when(
    `Bias level` == "01" ~ "Minimal",
    `Bias level` == "03" ~ "Strong")) |>
  pivot_longer(`Biased`:`Censored`)
  
plot_dat |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.9,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#019AFF"))) +
  facet_grid(`Bias level`~known_lambda,
             scales = "free",
             labeller = labeller(
               `Bias level` = label_value,
               known_lambda = label_value)) +
  theme_bw(base_size = 18) +
  labs(x = "\u03bb estimate",
       y = "density",
       fill = "Data source"
  )
ggsave("plots/lambda_ests_sfs2.png",
       units = "in",
       height = 6,
       width = 10)
