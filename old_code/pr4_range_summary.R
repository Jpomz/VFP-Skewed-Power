#gradient sim summary

library(tidyverse)
library(broom)
library(tidybayes)

results <- readRDS("simulation_results/pr1_4_range_sim_run.rds")

results <- as_tibble(results)
results

# how many parameter sets?
distinct(results,
         original_n, group, known_beta, pr_scenario, cutoff)

# xmin for under and trimmed data
results |>
  group_by(pr_scenario, cutoff) |>
  summarize(under_xmin_mean = mean(xmin_under),
            under_xmin_sd = sd(xmin_under),
            trimmed_xmin_mean = mean(xmin_trimmed),
            trimmed_xmin_sd = sd(xmin_trimmed))
# xmin for trimmed data sets approx = fixed cutoff value
results |>
  ggplot(aes(x = xmin_trimmed, 
         fill = pr_scenario))+
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = cutoff)) +
  facet_grid(known_lambda~cutoff,
             scales = "free")

# lambda estimates --------------------------------------------------------
results |>
  filter(pr_scenario == "01") |>
  select(cutoff,
         known_lambda,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.6,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#019AFF"))) +
  facet_grid(cutoff ~ known_lambda,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "density",
       title = "Minimal Bias")

results |>
  filter(pr_scenario == "01",
         known_lambda == -2.5 |
           known_lambda == -2 |
           known_lambda == -1.5) |>
  select(cutoff,
         known_lambda,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.6,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#019AFF"))) +
  facet_grid(cutoff ~ known_lambda,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "density",
       title = "Minimal Bias")
ggsave("plots/lambda_fixed_cut_minimal.png",
       units = "in",
       height = 6,
       width = 10)


results |>
  filter(pr_scenario == "04",
         known_lambda == -2.5 |
           known_lambda == -2 |
           known_lambda == -1.5) |>
  select(cutoff,
         known_lambda,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.6,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#019AFF"))) +
  facet_grid(cutoff ~ known_lambda,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "density",
       title = "Extreme Bias")
ggsave("plots/lambda_fixed_cut_extreme.png",
       units = "in",
       height = 6,
       width = 10)

# gradient plots ----------------------------------------------------------
results |>
  filter(pr_scenario == "01") |>
  ggplot(aes(x = env_gradient,
             y = lambda_under,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.25,
              linewidth = 0.25,
              color = "#FF1984") +
  facet_grid(cutoff ~ known_beta,
             labeller = label_value) +
  theme_bw() +
  labs(title = "Minimally Biased Data")

results |>
  filter(pr_scenario == "04") |>
  ggplot(aes(x = env_gradient,
             y = lambda_under,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.1,
              linewidth = 0.1,
              color = "#FF1984") +
  facet_grid(cutoff ~ known_beta,
             labeller = label_value) +
  theme_bw() +
  labs(title = "Extreme Biased data")

results |>
  filter(pr_scenario == "01") |>
  ggplot(aes(x = env_gradient,
             y = lambda_trimmed,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.25,
              linewidth = 0.25,
              color = "#019AFF") +
  facet_grid(cutoff ~ known_beta,
             labeller = label_value) +
  theme_bw() +
  labs(title = "Censored data, Minimal Bias")

results |>
  filter(pr_scenario == "04") |>
  ggplot(aes(x = env_gradient,
             y = lambda_trimmed,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.25,
              linewidth = 0.25,
              color = "#019AFF") +
  facet_grid(cutoff ~ known_beta,
             labeller = label_value) +
  theme_bw() +
  labs(title = "Censored data, Extreme Bias")


# plot rep lines ----------------------------------------------------------

results |>
  filter(rep %in% 1:100,
         pr_scenario == "04") |>
  pivot_longer(c(lambda_under,
                 lambda_trimmed),
               names_to = "data_model",
               values_to = "lambda_est") |>
  ggplot(aes(x = env_gradient,
             y = lambda_est,
             color = data_model,
             group = interaction(rep, data_model))) +
  geom_point(size = 1,
             alpha = 0.25) +
  geom_line(
    stat = "smooth",
    method = "lm", 
    alpha = 0.1, 
    linewidth = 1,
    se = FALSE) + 
  facet_grid(cutoff~known_beta) +
  scale_color_manual(values = c(c("#019AFF",
                                  "#FF1984"))) +
  theme_bw() +
  geom_line(inherit.aes = FALSE,
            aes(x = env_gradient,
                y = known_lambda),
            color = "black",
            linetype = "dashed",
            linewidth = 1) +
  labs(title = "Extreme Bias")
ggsave("plots/rep_lines_extreme_MS.png",
       units = "in",
       height = 6,
       width = 10)

results |>
  filter(rep %in% 1:100,
         pr_scenario == "01") |>
  pivot_longer(c(lambda_under,
                 lambda_trimmed),
               names_to = "data_model",
               values_to = "lambda_est") |>
  ggplot(aes(x = env_gradient,
             y = lambda_est,
             color = data_model,
             group = interaction(rep, data_model))) +
  geom_point(size = 1,
             alpha = 0.25) +
  geom_line(
    stat = "smooth",
    method = "lm", 
    alpha = 0.1, 
    linewidth = 1,
    se = FALSE) +
  facet_grid(cutoff~known_beta) +
  scale_color_manual(values = c(c("#019AFF",
                                  "#FF1984"))) +
  theme_bw() +
  geom_line(inherit.aes = FALSE,
            aes(x = env_gradient,
                y = known_lambda),
            color = "black",
            linetype = "dashed",
            linewidth = 1)+
  labs(title = "Minimal Bias")
ggsave("plots/rep_lines_minimal_MS.png",
       units = "in",
       height = 6,
       width = 10)

# lm models ---------------------------------------------------------------

lm_models <- results |>
  select(group,
         lambda_under,
         lambda_trimmed,
         pr_scenario,
         rep,
         cutoff,
         original_n,
         known_beta,
         env_gradient) |>
  pivot_longer(c(lambda_under,
                 lambda_trimmed),
               names_to = "data_model", 
               values_to = "lambda_est") |>
  group_by(group, 
           original_n, 
           pr_scenario,
           cutoff,
           rep,
           known_beta,
           data_model) |>
  nest() |>
  mutate(lm_model = map(
    data, ~lm(
      lambda_est ~ env_gradient, .))) |>
  mutate(coefs = map(
    lm_model, tidy)) |>
  unnest(coefs) 

# how many models incorrectly assumed a relationship when there was not one?
lm_models |>
  filter(group == "A") |>
  group_by(data_model,
           pr_scenario,
           cutoff) |>
  mutate(tot_rep = n())|>
  filter(term == "env_gradient", 
         known_beta == 0,
         p.value <0.05) |>
  group_by(data_model, 
           pr_scenario,
           cutoff, 
           tot_rep) |>
  summarise(
    result_count = n()) |>
  mutate(rate = result_count / tot_rep)


# how many models failed to detect a relationship when there was one?
lm_models |>
  filter(group != "A") |>
  group_by(data_model, pr_scenario) |>
  mutate(tot_rep = n())|>
  filter(term == "env_gradient", 
         known_beta != 0,
         p.value >0.05)|>
  group_by(data_model, 
           pr_scenario,
           cutoff, 
           tot_rep) |>
  summarise(
    result_count = n()) |>
  mutate(rate = result_count / tot_rep)


# beta distributions ------------------------------------------------------

lm_models |>
  filter(original_n == 1000,
         term == "env_gradient",
         cutoff != 1,
         pr_scenario == "01") |>
  mutate(`Bias level` = case_when(
    pr_scenario == "01" ~ "Minimal",
    pr_scenario =="04" ~ "Extreme"
  ),
  `Bias level` = factor(`Bias level`,
                        levels = c("Minimal", 
                                   "Extreme"))) |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(cutoff~known_beta,
             scales = "free") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Minimal Bias",
       x = "beta estimate")
ggsave("plots/beta_fixed_minimal_MS.png",
       units = "in",
       height = 6,
       width = 10)


lm_models |>
  filter(original_n == 1000,
         term == "env_gradient",
         #cutoff != 0.001,
         pr_scenario == "04") |>
  mutate(`Bias level` = case_when(
    pr_scenario == "01" ~ "Minimal",
    pr_scenario =="04" ~ "Extreme"
  ),
  `Bias level` = factor(`Bias level`,
                        levels = c("Minimal", 
                                   "Extreme"))) |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(cutoff~known_beta,
             scales = "free") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Extreme Bias",
       x = "beta estimate")
ggsave("plots/beta_fixed_extreme_MS.png",
       units = "in",
       height = 6,
       width = 10)


# false positive beta distribution ----------------------------------------


lm_models |>
  filter(term == "env_gradient",
         p.value < 0.05,
         known_beta == 0,
         pr_scenario == "01") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(pr_scenario~cutoff) +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  labs(title = "False Positive Relationship, minimal bias")

lm_models |>
  filter(term == "env_gradient",
         p.value < 0.05,
         known_beta == 0,
         pr_scenario == "04") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(pr_scenario~cutoff) +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  labs(title = "False Positive Relationship, extreme bias")
  

# false negative beta distribution ----------------------------------------


lm_models |>
  filter(term == "env_gradient",
         p.value >= 0.05,
         known_beta !=0,
         pr_scenario == "01") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff) +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  labs(title = "False Negative Relationship, Minimal Bias")

lm_models |>
  filter(term == "env_gradient",
         p.value >= 0.05,
         known_beta !=0,
         pr_scenario == "04") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff) +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  labs(title = "False Negative Relationship, Extreme Bias")


# True positive -----------------------------------------------------------
lm_models |>
  filter(term == "env_gradient",
         p.value < 0.05,
         known_beta != 0,
         pr_scenario == "01") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff) +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  labs(title = "True Positive Relationship, Minimal Bias")

lm_models |>
  filter(term == "env_gradient",
         p.value < 0.05,
         known_beta != 0,
         pr_scenario == "04") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff) +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  labs(title = "True Positive Relationship, Extreme Bias")


# True negative -----------------------------------------------------------


lm_models |>
  filter(term == "env_gradient",
         p.value >= 0.05,
         known_beta == 0,
         pr_scenario== "01") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff) +
  labs(title = "True Negative Relationship, Minimal Bias") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984")))

lm_models |>
  filter(term == "env_gradient",
         p.value >= 0.05,
         known_beta == 0,
         pr_scenario== "04") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff) +
  labs(title = "True Negative Relationship, Extreme Bias") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984")))



# bias tables -------------------------------------------------------------

lm_table_dat <- lm_models |>
  ungroup() |>
  select(group,
         pr_scenario,
         rep,
         cutoff,
         known_beta,
         data_model,
         term, 
         estimate,
         std.error) |>
  filter(term == "env_gradient")
  
lm_table_dat |>
  distinct(known_beta,
           pr_scenario,
           cutoff,
           data_model)
# 2 scenarios, 4 groups, 4 cutoffs, 2 models
2*4*4 *2
# 64 groups

beta_bias <- lm_table_dat |> 
  mutate(target = known_beta,
         minCI = estimate - 1.96*std.error,
         maxCI = estimate + 1.96*std.error,
         target_name = "Regression Slope") |>
  mutate(conf_width = maxCI - minCI,
         diff = estimate - known_beta,
         abs_bias = abs(diff)) |>
  group_by(known_beta, 
           cutoff, 
           pr_scenario, 
           data_model) |>
  add_count() |>
  #group_by(target_name, name, n) %>%
  summarize(
    median_ci_range = median(conf_width,
                             na.rm = TRUE),
    median_abs_bias = median(abs_bias,
                             na.rm = TRUE),
    sd_abs_bias = sd(abs_bias,
                     na.rm = TRUE)) |>
  arrange(cutoff,
          known_beta,
          pr_scenario,
          data_model) 
write_csv(beta_bias, 
          "simulation_summaries/beta_bias_table_fixed_cutoff.csv")

lm_table_dat |>
  filter(known_beta == -0.5, 
         pr_scenario == "03", 
         original_n == 500)
