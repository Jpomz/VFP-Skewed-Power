#gradient sim summary

library(tidyverse)
library(broom)
library(tidybayes)

results <- readRDS("simulation_results/gradient_sim_run.rds")

results <- as_tibble(results)
results

# how many parameter sets?
distinct(results, original_n, group, pr_scenario)
# 48 sets
48 * #sets
  500 * # reps # 24000 
  5 # lambdas across gradient
# 120,000 rows


# plot rep lines ----------------------------------------------------------

results |>
  filter(rep %in% 1:100,
         pr_scenario == "01"|
         pr_scenario == "04") |>
  mutate(`Bias level` = case_when(
    pr_scenario == "01" ~ "Minimal",
    pr_scenario == "04" ~ "Extreme"
  ),
  `Bias level` = factor(`Bias level`,
                        levels = c("Minimal", 
                                   "Extreme"))) |>
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
  facet_grid(known_beta~`Bias level`) +
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
ggsave("plots/rep_lines_clauset_MS.png",
       units = "in",
       height = 6,
       width = 10)




lm_models <- results |>
  select(group,
         lambda_under,
         lambda_trimmed,
         pr_scenario,
         rep,
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
  group_by(data_model) |>
  count()
# 12000
lm_models |>
  filter(term == "env_gradient", 
         known_beta == 0,
         p.value <0.05) |>
  group_by(data_model) |>
  count()
311 / 12000
293 / 12000
# ~ 2.5% had false positive relationships

# how many models failed to detect a relationship when there was one?
lm_models |>
  filter(group != "A") |>
  group_by(data_model) |>
  count()
# 36000 sets per data model should be significant

lm_models |>
  filter(term == "env_gradient", 
         known_beta != 0,
         p.value >0.05)|>
  group_by(data_model) |>
  count()
12000*3 # total reps
9831 / 36000 # 2.7%
6791 / 36000 # 1.8%

lm_models |>
  filter(term == "env_gradient",
         p.value < 0.05,
         known_beta == 0,
         original_n == 5000) |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(pr_scenario~group) +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  labs(title = "False Positive Relationship")
  
lm_models |>
  filter(term == "env_gradient",
         p.value >= 0.05,
         known_beta != 0,
         original_n == 5000) |>
  #mutate(known_beta = known_beta *-1) |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(pr_scenario~known_beta,
             labeller = 
               labeller(pr_scenario = label_value,
                        known_beta = label_both)) +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  labs(title = "False Negative Relationship")

lm_models |>
  filter(term == "env_gradient",
         p.value < 0.05,
         known_beta != 0,
         original_n == 5000) |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(pr_scenario~known_beta,
             labeller = label_both) +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  labs(title = "True Positive Relationship")


lm_models |>
  filter(term == "env_gradient",
         p.value >= 0.05,
         known_beta == 0,
         original_n == 5000) |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(pr_scenario~known_beta,
             labeller = 
               labeller(pr_scenario = label_value,
                        known_beta = label_both)) +
  labs(title = "True Negative Relationship") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984")))

lm_models |>
  filter(original_n == 1000,
         term == "env_gradient",
         pr_scenario == "01" |
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
  facet_grid(known_beta~`Bias level`,
             scales = "free") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  theme_bw()
ggsave("plots/beta_clauset_MS.png",
       units = "in",
       height = 6,
       width = 10)



# lm summaries ------------------------------------------------------------
# lm_models |>
#   ungroup() |>
#   group_by(original_n, 
#            pr_scenario,
#            known_beta,
#            data_model) |>
#   mutate(rep_n = n(),
#          fasle_pos = case_when(
#            known_beta == 0 &
#            p.value < 0.05 ~ 1
#          ))




# Plot reps ---------------------------------------------------------------
results |>
  group_by(known_beta, 
           pr_scenario,
           original_n, 
           rep) |>
  count() |>
  group_by(n) |>
  count()

results |>
  filter(original_n == 5000) |>
  ggplot(aes(x = env_gradient,
             y = lambda_under,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.25,
              linewidth = 0.25,
              color = "#FF1984") +
  facet_grid(pr_scenario~known_beta) +
  theme_bw() +
  labs(title = "Biased data")

results |>
  filter(original_n == 5000) |>
  ggplot(aes(x = env_gradient,
             y = lambda_trimmed,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.25,
              linewidth = 0.25,
              color = "#019AFF") +
  facet_grid(pr_scenario~known_beta) +
  theme_bw() +
  labs(title = "Censored data")

results |>
  filter(original_n == 5000) |>
  pivot_longer(c(lambda_under,
                 lambda_trimmed),
               names_to = "data_model", 
               values_to = "lambda_est") |>
  group_by(known_beta, pr_scenario, original_n) |>
  sample_n(10) |>
  ggplot(aes(x = env_gradient,
             y = lambda_est,
             group = rep,
             color = data_model)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.5,
              linewidth = 0.5) +
  facet_grid(pr_scenario~known_beta) +
  theme_bw()


# bias tables -------------------------------------------------------------

lm_table_dat <- lm_models |>
  ungroup() |>
  select(group, pr_scenario, rep, original_n,
         known_beta, data_model, term, 
         estimate, std.error) |>
  filter(term == "env_gradient")
  
lm_table_dat |>
  distinct(original_n, known_beta, pr_scenario)

beta_bias <- lm_table_dat |> 
  mutate(target = known_beta,
         minCI = estimate - 1.96*std.error,
         maxCI = estimate + 1.96*std.error,
         target_name = "Regression Slope") |>
  mutate(conf_width = maxCI - minCI,
         diff = estimate - known_beta,
         abs_bias = abs(diff)) |>
  group_by(known_beta, 
           original_n, 
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
  arrange(original_n,
          known_beta,
          pr_scenario,
          data_model) 
write_csv(beta_bias, 
          "simulation_summaries/beta_bias_table.csv")

lm_table_dat |>
  filter(known_beta == -0.5, 
         pr_scenario == "03", 
         original_n == 500)
