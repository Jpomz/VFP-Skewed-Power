# Plots showing sampling probablity function and examples

# function ####
sample_pr <- function(x, h, b){
  pr = 1 / (1 + (h / x**b))
  return(pr)
}

# SI figures for sampling probability with simulation values
library(tidyverse)
library(sizeSpectra)
library(poweRlaw)

source("master_variable_designation.R")
set.seed(2112)
pr_scenarios
dat <- expand_grid(
  pr_scenarios,
  x = c(exp(seq(log(0.001), log(10), length.out = 1000)))) |>
mutate(pr = sample_pr(x = x, h = h, b = b))


## 4 scenarios on one plot
dat |>
  mutate(scenario = case_when(
    h == 0.00001 & b == 1.5 ~ "minimal",
    h == 0.0001 & b == 1.5 ~ "moderate",
    h == 0.001 & b == 1.5 ~ "strong",
    h == 0.01 & b == 1.5 ~ "extreme",
  ),
  scenario = factor(scenario, 
                       levels = c("minimal", "moderate", "strong", "extreme"))) |>
  ggplot(aes(x = x, 
             y = pr, 
             color = scenario)) +
  geom_line(linewidth = 2) +
  geom_vline(aes(xintercept = 0.01), linetype = "dashed") +
  geom_hline(aes(yintercept = 0.9), linetype = "dashed") +
  scale_x_log10() +
  theme_bw() +
  scale_colour_manual(values = c("darkorchid1",
                                 "darkorchid2",
                                 "darkorchid3",
                                 "darkorchid4"),
                      name = "Bias") +
  labs(title = "Sampling probability as a function of body mass",
       x =expression(Log[10]~dry~mass),
       y = "Sampling probability") +
  guides(x = "axis_logticks")

ggsave("plots/bias_scenarios_MS.png",
       units = "in",
       height = 4,
       width = 6.5)



# Example of undersampling  -----------------------------------------------


# make SI figures of "real" data and undersampling results

# simulate x from bounded power law
set.seed(205)
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
  h = c(0.00001, 0.01),
  b = 1.5) |>
  mutate(scenario = 
           case_when(
             h == 0.00001 & b == 1.5 ~ "minimal",
             h == 0.01 & b == 1.5 ~ "extreme"))

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
         scenario == "extreme") |>
  pull(x)
x_power_str <- conpl$new(sampled_x_str)
x_xmin_str <- estimate_xmin(x_power_str)$xmin
x_xmin_str

both_dats |>
  mutate(x_min = case_when(
    scenario == "minimal" ~ 0.00279,
    scenario == "extreme" ~ 0.05501
  )) |>
  filter(value!= "not sampled") |>
  mutate(value = case_when(
    value == "sampled" ~ "Biased", 
    .default = "Original"
  ),
  value = factor(value, 
                 levels = c("Original", "Biased")),
  scenario = factor(scenario,
                    levels = c("minimal", 
                               "extreme"))) |>
  ggplot(aes(x = x, 
             fill = value)) +
  geom_histogram(binwidth = 0.05, 
                 position = "dodge",
                 alpha = 0.75) +
  scale_x_log10() +
  facet_wrap(~scenario) +
  geom_segment(aes(x = x_min,
                   y = 550,
                   xend = x_min,
                   yend = 450),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 1.5,
               color = "dodgerblue") +
  scale_fill_manual(values = c("black", "#FF1984")) + 
  theme_bw(base_size = 18) +
  guides(fill = guide_legend(title= "Data")) +
  labs(x =expression(Log[10]~dry~mass)) +
  guides(x = "axis_logticks")

ggsave("plots/orig_bias_xmin_SI.png",
       units = "in",
       height = 4,
       width = 8)



### Move this to the summary script
x_pr_subsample <- function(
    df
    #N_orig,
    #N_sub,
    #lambda, 
    #xmin, 
    #xmax,
    #h,
    #b
    #x_pr
    ) {
  # x <- rPLB(n = N_orig,
  #           b = lambda,
  #           xmin = xmin, 
  #           xmax = xmax)
  
  # x_pr <- data.frame(x = x) |>
  #   mutate(pr = sample_pr(x, 
  #                         b = 2,
  #                         h = h))
  x_90 <- df |>
    filter(pr >=0.89,
           pr <=0.91) |>
    group_by(h, b) |>
    summarize(x_90 = mean(x, na.rm = TRUE)) |>
    ungroup() |>
    select(-h, -b)
  
  x_80 <- df |>
    filter(pr >=0.79,
           pr <=0.81) |>
    group_by(h, b) |>
    summarize(x_80 = mean(x, na.rm = TRUE))
  
  x_95 <- df |>
    filter(pr >=0.94,
           pr <=0.96) |>
    group_by(h, b) |>
    summarize(x_95 = mean(x, na.rm = TRUE))|>
    ungroup() |>
    select(-h, -b)
  
  x_85 <- df |>
    filter(pr >=0.84,
           pr <=0.86) |>
    group_by(h, b) |>
    summarize(x_85 = mean(x, na.rm = TRUE))|>
    ungroup() |>
    select(-h, -b)
  
  out <- bind_cols(x_80, 
                   x_85,
                   x_90,
                   x_95)
  #out$N_orig = N_orig
  #out$N_sub = N_sub
  return(out)
}

dat %>%
  #group_by(h, b) %>%
  reframe(x_pr_subsample(.)) |>
  arrange(b, h) |>
  mutate(across(starts_with("x_"), ~ round(.x, 3))) |>
  write_csv("simulation_summaries/body_size_sampling_probabillities.csv")

dat %>%
  #group_by(h, b) %>%
  reframe(x_pr_subsample(.)) |>
  pivot_longer(x_80:x_95) |>
  separate(name, into = c("x", "probability")) |>
  mutate(probability = as.numeric(probability)) |>
  ggplot(aes(x = value, 
             y = probability,
             color = interaction(h, b))) +
  geom_point()
  