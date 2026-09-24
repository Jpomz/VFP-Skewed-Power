#gradient sim summary

library(tidyverse)
library(broom)
library(tidybayes)

results <- readRDS("simulation_results/gradient_sim_run_morin_correction.rds")

results <- as_tibble(results) 
results


# lambda estimates --------------------------------------------------------
results |>
  ggplot(aes(x = lambda_trimmed)) +
  stat_halfeye(normalize = "groups",
               fill = "springgreen") +
  facet_grid(bias_level~known_lambda,
             scales = "free_x") +
  geom_vline(aes(xintercept = known_lambda))

# lambda bias ####
results |>
  filter(known_lambda == -1.9 |
           known_lambda == -2 |
           known_lambda == -2.1) |>
  ungroup() |>
  mutate(abs_delta = abs(known_lambda - lambda_trimmed)) |>
  ggplot(aes(x = abs_delta)) +
  stat_halfeye(normalize = "groups",
               fill = "springgreen") +
  facet_grid(known_lambda ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw()+
  labs(x = "\u0394 in \u03bb estimates",
       y = "") +
  scale_x_continuous(n.breaks = 4)
ggsave("plots/lambda_deltas_inverse.png",
       units = "in",
       height = 10,
       width = 10)

results |>
  filter(known_lambda == -1.9 |
           known_lambda == -2 |
           known_lambda == -2.1) |>
  mutate(
    abs_delta = abs(known_lambda - lambda_trimmed)) |>
  group_by(known_lambda, bias_level)|>
  summarise(median_abs_delta = median(abs_delta),
            sd_abs_delta = sd(abs_delta)) |>
  write_csv("simulation_summaries/inverse_lambda_deviation.csv")

# lambda CIs ####
corrected_lambda_ci <- results |>
  select(known_lambda:lambda_trimmed_hi, M) |>
  filter(known_lambda == -1.9|
           known_lambda == -2 |
           known_lambda == -2.1) |>
  mutate(in_fixed_ci = 
           known_lambda > lambda_trimmed_lo &
           known_lambda < lambda_trimmed_hi) |>
  select(known_lambda,
         M,
         lambda_trimmed,
         lambda_trimmed_lo,
         lambda_trimmed_hi,
         in_fixed_ci)

# fixed ci plot ####
corrected_lambda_ci |>
  group_by(M, known_lambda) |>
  arrange(lambda_trimmed_lo) |>
  mutate(y = (1:n())) |>
  ggplot(aes(x = lambda_trimmed, 
             xmin = lambda_trimmed_lo,
             xmax = lambda_trimmed_hi,
             y = y,
             color = in_fixed_ci)) +
  geom_pointrange(fatten = 0,
                  linewidth = 0.25,
                  alpha = 0.75) +
  geom_vline(aes(xintercept = known_lambda)) +
  facet_grid(known_lambda~M,
             scales = "free") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(title = "Corrected: 1/retention probability",
       y = "",
       x = "95% CI") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("plots/corrected_lambda_CIs.png",
       units = "in",
       height = 10,
       width = 10)

inverse_lambda_ci_percent <- corrected_lambda_ci |>
  mutate(
    ci_width = lambda_trimmed_hi - lambda_trimmed_lo) |>
  group_by(known_lambda, M) |>
  summarize(n = n(), 
            prop_in = sum(in_fixed_ci),
            mean_ci_width = mean(ci_width),
            sd_ci_width = sd(ci_width)) |>
  mutate(percent_ci = prop_in / n,
         data = "Corrected") |>
  select(known_lambda,
         M,
         data,
         percent_ci,
         mean_ci_width, 
         sd_ci_width) 
inverse_lambda_ci_percent |>
  mutate(bias_level = case_when(
    M == 0.125 ~ "minimal",
    M == 0.25 ~ "moderate",
    M == 0.5 ~ "strong",
    M == 1 ~ "extreme"
  )) |>
  arrange(-known_lambda, M, data) |>
  select(known_lambda, bias_level, data:sd_ci_width) |>
  write_csv("simulation_summaries/inverse_lambda_CI_percent.csv")

# rep lines ---------------------------------------------------------------

results |>
  filter(rep %in% 1:100) |>
  pivot_longer(c(lambda_under,
                 lambda_trimmed),
               names_to = "data_model",
               values_to = "lambda_est") |>
  select(data_model,
         rep, 
         lambda_est,
         known_beta, 
         bias_level,
         known_lambda,
         env_gradient 
  ) |>
  group_by(bias_level,
           rep,
           known_beta,
           data_model) |>
  add_count() |>
  filter(n ==5) |>
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
  facet_grid(known_beta~bias_level) +
  scale_color_manual(values = c(c("springgreen",
                                  "#FF1984"))) +
  theme_bw() +
  geom_line(inherit.aes = FALSE,
            aes(x = env_gradient,
                y = known_lambda),
            color = "black",
            linetype = "dashed",
            linewidth = 1) +
  labs(title = "Inverse Correction")
ggsave("plots/rep_lines_morin_inverse_SI.png",
       units = "in",
       height = 6,
       width = 10)

# lm models ---------------------------------------------------------------

lm_models <- results |>
  select(group,
         lambda_trimmed,
         bias_level,
         rep,
         original_n,
         known_beta,
         env_gradient) |>
  group_by(group, 
           original_n, 
           bias_level,
           rep,
           known_beta) |>
  add_count() |>
  filter(n ==5) |>
  nest() |>
  mutate(lm_model = map(
    data, ~lm(
      lambda_trimmed ~ env_gradient, .))) |>
  mutate(coefs = map(
    lm_model, tidy)) |>
  unnest(coefs) 

# beta distributions ------------------------------------------------------

lm_models |>
  filter(term == "env_gradient") |>
  ggplot(aes(x = estimate)) +
  stat_halfeye(normalize = "groups",
               fill = "springgreen") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~bias_level,
             scales = "free")  +
  theme_bw() +
  labs(title = "Inverse correction",
       x = "beta estimate")
ggsave("plots/beta_morin_inverse.png",
       units = "in",
       height = 6,
       width = 10)

lm_models |>
  ungroup() |>
  filter(term == "env_gradient") |>
  mutate(abs_delta = abs(known_beta - estimate)) |>
  group_by(known_beta, bias_level) |>
  summarise(median_abs_delta = median(abs_delta),
            sd_abs_delta = sd(abs_delta))|>
  write_csv("simulation_summaries/inverse_betas_deviation.csv")


lm_models |>
  ungroup() |>
  filter(term == "env_gradient") |>
  mutate(abs_delta = abs(known_beta - estimate)) |>
  ggplot(aes(x = abs_delta)) +
  stat_halfeye(normalize = "groups",
               fill = "springgreen") +
  facet_grid(known_beta ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw()+
  labs(x = "\u0394 in \u03b2 estimates",
       y = "") +
  scale_x_continuous(n.breaks = 4)
ggsave("plots/beta_deltas_inverse.png",
       units = "in",
       height = 10,
       width = 10)

# beta by bias_level ---------------------------------------------------------
beta_plot_dat <- lm_models |>
  mutate(keep = case_when(
    cutoff == 0.584 & bias_level == "minimal" ~ TRUE,
    cutoff == 1.410 & bias_level == "moderate" ~ TRUE,
    cutoff == 3.63 & bias_level == "strong" ~ TRUE,
    cutoff == 10.1 & bias_level == "extreme" ~ TRUE,
    .default = FALSE
  )) |> 
  filter(keep == TRUE)
saveRDS(beta_plot_dat,
        "simulation_summaries/beta_fixed_plot_data.RDS")

beta_fixed_00_panel <- beta_plot_dat |>
  filter(known_beta == 0,
         term == "env_gradient") |>
  select(known_beta, 
         bias_level,
         data_model,
         estimate) |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_beta),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(known_beta ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "")+
  scale_x_continuous(n.breaks = 5)

beta_fixed_01_panel <- beta_plot_dat |>
  filter(known_beta == -0.1,
         term == "env_gradient") |>
  select(known_beta, 
         bias_level,
         data_model,
         estimate) |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_beta),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(known_beta ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "")+
  scale_x_continuous(n.breaks = 5)

beta_fixed_25_panel <- beta_plot_dat |>
  filter(known_beta == -0.25,
         term == "env_gradient") |>
  select(known_beta, 
         bias_level,
         data_model,
         estimate) |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_beta),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(known_beta ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "")+
  scale_x_continuous(n.breaks = 5)

beta_fixed_05_panel <- beta_plot_dat |>
  filter(known_beta == -0.5,
         term == "env_gradient") |>
  select(known_beta, 
         bias_level,
         data_model,
         estimate) |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_beta),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(known_beta ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03b2 estimate",
       y = "")+
  scale_x_continuous(n.breaks = 5)

ggpubr::ggarrange(beta_fixed_00_panel,
                  beta_fixed_01_panel,
                  beta_fixed_25_panel,
                  beta_fixed_05_panel,
                  common.legend = TRUE,
                  legend = "right",
                  ncol = 1)
ggsave("plots/beta_fixed.png",
       units = "in",
       height = 6,
       width = 10)

# false positive beta distribution ----------------------------------------


lm_models |>
  filter(term == "env_gradient",
         p.value < 0.05,
         known_beta == 0,
         bias_level == "minimal") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(.~cutoff,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  labs(title = "False Positive Relationship, minimal bias")

lm_models |>
  filter(term == "env_gradient",
         p.value < 0.05,
         known_beta == 0,
         bias_level == "extreme") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(bias_level~cutoff,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  labs(title = "False Positive Relationship, extreme bias")
  

# false negative beta distribution ----------------------------------------


lm_models |>
  filter(term == "env_gradient",
         p.value >= 0.05,
         known_beta !=0,
         bias_level == "minimal") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  labs(title = "False Negative Relationship, Minimal Bias")

lm_models |>
  filter(term == "env_gradient",
         p.value >= 0.05,
         known_beta !=0,
         bias_level == "extreme") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  labs(title = "False Negative Relationship, Extreme Bias")


# True positive -----------------------------------------------------------
lm_models |>
  filter(term == "env_gradient",
         p.value < 0.05,
         known_beta != 0,
         bias_level == "minimal") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  labs(title = "True Positive Relationship, Minimal Bias")

lm_models |>
  filter(term == "env_gradient",
         p.value < 0.05,
         known_beta != 0,
         bias_level == "extreme") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  labs(title = "True Positive Relationship, Extreme Bias")


# True negative -----------------------------------------------------------


lm_models |>
  filter(term == "env_gradient",
         p.value >= 0.05,
         known_beta == 0,
         bias_level == "minimal") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff,
             scales = "free") +
  labs(title = "True Negative Relationship, Minimal Bias") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984")))

lm_models |>
  filter(term == "env_gradient",
         p.value >= 0.05,
         known_beta == 0,
         bias_level == "extreme") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff,
             scales = "free") +
  labs(title = "True Negative Relationship, Extreme Bias") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984")))



# bias tables -------------------------------------------------------------

lm_table_dat <- lm_models |>
  ungroup() |>
  select(group,
         bias_level,
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
           bias_level,
           cutoff,
           data_model)
# 4 bias levels, 4 groups, 4 cutoffs, 2 models
4*4*4 *2
# 128 groups

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
           bias_level, 
           data_model) |>
  #add_count() |>
  #group_by(target_name, name, n) %>%
  summarize(
    median_ci_range = median(conf_width,
                             na.rm = TRUE),
    median_abs_bias = median(abs_bias,
                             na.rm = TRUE),
    sd_abs_bias = sd(abs_bias,
                     na.rm = TRUE),
    n = n()) |>
  arrange(cutoff,
          known_beta,
          bias_level,
          data_model) 
write_csv(beta_bias, 
          "simulation_summaries/beta_bias_table_morin_fixed_cutoff.csv")

beta_bias |>
  filter(bias_level == "minimal",
         cutoff == 0.584)

beta_bias |>
  filter(bias_level == "strong",
         cutoff == 3.63)
beta_bias |>
  filter(bias_level == "extreme",
         cutoff == 10.1)
