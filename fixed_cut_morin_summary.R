#gradient sim summary

library(tidyverse)
library(broom)
library(tidybayes)

results <- readRDS("simulation_results/fixed_cut_morin_sim_run.rds")

results <- as_tibble(results) 
results

# length cutoffs = 0.584, 1.41, 3.63, 10.1
# masses ~ 0.00143, 0.0167, 0.233, 4.04
cut_mass <- tibble(cutoff = unique(results$cutoff),
       cutoff_mass = sizeSpectra::lengthToMass(cutoff, 
                                               LWa = 0.0064,
                                               LWb = 2.788))

results <- results |>
  left_join(cut_mass)
# how many parameter sets?
distinct(results,
         original_n, group, known_beta, M, cutoff)

# xmin for under and trimmed data
results |>
  group_by(cutoff_mass) |>
  summarize(under_xmin_median = median(xmin_under),
            trimmed_xmin_median = median(xmin_trimmed),
            trimmed_xmin_sd = sd(xmin_trimmed)) |>
  mutate(mass_bigger = trimmed_xmin_median > cutoff_mass)
# length cutoffs match min trimmed masses


# xmin for trimmed data sets approx = fixed cutoff value
results |>
  filter(known_lambda == -2 |
           known_lambda == -1.5 |
           known_lambda == -1) |>
  ggplot(aes(x = xmin_trimmed, 
         fill = bias_level))+
  stat_halfeye(alpha = 0.5,
               normalize = "panels") +
  geom_vline(aes(xintercept = cutoff_mass)) +
  facet_grid(known_lambda~cutoff_mass,
             scales = "free") +
  scale_x_log10(guide = "axis_logticks")

# lambda estimates --------------------------------------------------------
results |>
  filter(bias_level == "minimal",
         known_lambda == -1.25 |
           known_lambda == -1.5 |
           known_lambda == -2) |>
  select(cutoff,
         known_lambda,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(
               normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(cutoff ~ known_lambda,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "density",
       title = "Minimal Bias",
       subtitle = "Fixed cutoffs")

ggsave("plots/lambda_morin_fixed_cut_minimal.png",
       units = "in",
       height = 6,
       width = 10)


results |>
  filter(bias_level == "moderate",
         known_lambda == -1 |
           known_lambda == -1.5 |
           known_lambda == -2) |>
  select(cutoff,
         known_lambda,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value, 
             y = name,
             fill = name)) +
  stat_halfeye(
               normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(cutoff ~ known_lambda,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "density",
       title = "Moderate Bias",
       subtitle = "Fixed cutoffs")
ggsave("plots/lambda_morin_fixed_cut_moderate.png",
       units = "in",
       height = 6,
       width = 10)


results |>
  filter(bias_level == "strong",
         known_lambda == -1 |
           known_lambda == -1.5 |
           known_lambda == -2) |>
  select(cutoff,
         known_lambda,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value, 
             y = name,
             fill = name)) +
  stat_halfeye(
               normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(cutoff ~ known_lambda,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "density",
       title = "Strong Bias",
       subtitle = "Fixed cutoffs")
ggsave("plots/lambda_morin_fixed_cut_strong.png",
       units = "in",
       height = 6,
       width = 10)

results |>
  filter(bias_level == "extreme",
         known_lambda == -1 |
           known_lambda == -1.5 |
           known_lambda == -2) |>
  select(cutoff,
         known_lambda,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value, 
             y = name,
             fill = name)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(cutoff ~ known_lambda,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "density",
       title = "Extreme Bias",
       subtitle = "Fixed cutoffs")
ggsave("plots/lambda_morin_fixed_cut_extreme.png",
       units = "in",
       height = 6,
       width = 10)


# Lambda with 95% retention cutoffs ---------------------------------------
# lambda_fixed_minimal_panel <- results |>
#   filter(bias_level == "minimal",
#          cutoff == 0.584,
#          known_lambda == -1.25 |
#            known_lambda == -1.5 |
#            known_lambda == -1.75) |>
#   select(cutoff,
#          bias_level,
#          known_lambda,
#          lambda_under,
#          lambda_trimmed) |>
#   rename(`Biased` = lambda_under,
#          `Censored` = lambda_trimmed) |>
#   pivot_longer(`Biased`:`Censored`) |>
#   ggplot(aes(x = value,
#              y = name,
#              fill = name)) +
#   stat_halfeye(alpha = 0.6,
#                normalize = "groups") +
#   geom_vline(aes(xintercept = known_lambda),
#              linetype = "dashed") +
#   scale_fill_manual(values = c(c("#FF1984",
#                                  "#FF914A"))) +
#   facet_grid(known_lambda~bias_level,
#              labeller = label_value,
#              scales = "free") +
#   theme_bw() +
#   labs(x = "",
#        y = "") +
#   scale_x_continuous(n.breaks = 3)
# 
# lambda_fixed_moderate_panel <- results |>
#   filter(bias_level == "moderate",
#          cutoff == 1.41,
#          known_lambda == -1.25 |
#            known_lambda == -1.5 |
#            known_lambda == -1.75) |>
#   select(cutoff,
#          bias_level,
#          known_lambda,
#          lambda_under,
#          lambda_trimmed) |>
#   rename(`Biased` = lambda_under,
#          `Censored` = lambda_trimmed) |>
#   pivot_longer(`Biased`:`Censored`) |>
#   ggplot(aes(x = value,
#              y = name,
#              fill = name)) +
#   stat_halfeye(alpha = 0.6,
#                normalize = "groups") +
#   geom_vline(aes(xintercept = known_lambda),
#              linetype = "dashed") +
#   scale_fill_manual(values = c(c("#FF1984",
#                                  "#FF914A"))) +
#   facet_grid(known_lambda~bias_level,
#              labeller = label_value,
#              scales = "free") +
#   theme_bw() +
#   labs(x = "",
#        y = "") +
#   scale_x_continuous(n.breaks = 3)
# 
# lambda_fixed_strong_panel <- results |>
#   filter(bias_level == "strong",
#          cutoff == 3.63,
#          known_lambda == -1.25 |
#            known_lambda == -1.5 |
#            known_lambda == -1.75) |>
#   select(cutoff,
#          bias_level,
#          known_lambda,
#          lambda_under,
#          lambda_trimmed) |>
#   rename(`Biased` = lambda_under,
#          `Censored` = lambda_trimmed) |>
#   pivot_longer(`Biased`:`Censored`) |>
#   ggplot(aes(x = value,
#              y = name,
#              fill = name)) +
#   stat_halfeye(alpha = 0.6,
#                normalize = "groups") +
#   geom_vline(aes(xintercept = known_lambda),
#              linetype = "dashed") +
#   scale_fill_manual(values = c(c("#FF1984",
#                                  "#FF914A"))) +
#   facet_grid(known_lambda~bias_level,
#              labeller = label_value,
#              scales = "free") +
#   theme_bw() +
#   labs(x = "",
#        y = "") +
#   scale_x_continuous(n.breaks = 3)
# 
# lambda_fixed_extreme_panel <- results |>
#   filter(bias_level == "extreme",
#          cutoff == 10.1,
#          known_lambda == -1.25 |
#            known_lambda == -1.5 |
#            known_lambda == -1.75) |>
#   select(cutoff,
#          bias_level,
#          known_lambda,
#          lambda_under,
#          lambda_trimmed) |>
#   rename(`Biased` = lambda_under,
#          `Censored` = lambda_trimmed) |>
#   pivot_longer(`Biased`:`Censored`) |>
#   ggplot(aes(x = value,
#              y = name,
#              fill = name)) +
#   stat_halfeye(alpha = 0.6,
#                normalize = "groups") +
#   geom_vline(aes(xintercept = known_lambda),
#              linetype = "dashed") +
#   scale_fill_manual(values = c(c("#FF1984",
#                                  "#FF914A"))) +
#   facet_grid(known_lambda~bias_level,
#              labeller = label_value,
#              scales = "free") +
#   theme_bw() +
#   labs(x = "\u03bb estimate",
#        y = "") +
#   scale_x_continuous(n.breaks = 3)
# 
# ggpubr::ggarrange(lambda_fixed_minimal_panel,
#                   lambda_fixed_moderate_panel,
#                   lambda_fixed_strong_panel,
#                   lambda_fixed_extreme_panel,
#                   common.legend = TRUE,
#                   legend = "right",
#                   nrow = 1)
# ggsave("plots/lambda_morin_fixed.png",
#        units = "in",
#        height = 6,
#        width = 10)

# just 95% ####
lambda_plot_dat <- results |>
  mutate(keep = case_when(
    cutoff == 0.584 & M == 0.125 ~ TRUE,
    cutoff == 1.410 & M == 0.25 ~ TRUE,
    cutoff == 3.63 & M == 0.5 ~ TRUE,
    cutoff == 10.1 & M == 1 ~ TRUE,
    .default = FALSE
  )) |> 
  filter(keep == TRUE)


lambda_fixed_125_panel <- lambda_plot_dat |>
  filter(known_lambda == -1.25) |>
  select(cutoff,
         bias_level,
         known_lambda,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(known_lambda~bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "") +
  scale_x_continuous(n.breaks = 3)

lambda_fixed_15_panel <- lambda_plot_dat |>
  filter(known_lambda == -1.5) |>
  select(cutoff,
         bias_level,
         known_lambda,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(known_lambda~bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "") +
  scale_x_continuous(n.breaks = 3)

lambda_fixed_175_panel <- lambda_plot_dat |>
  filter(known_lambda == -1.75) |>
  select(cutoff,
         bias_level,
         known_lambda,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(known_lambda~bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "") +
  scale_x_continuous(n.breaks = 3)
ggpubr::ggarrange(lambda_fixed_125_panel,
                  lambda_fixed_15_panel,
                  lambda_fixed_175_panel,
                  common.legend = TRUE,
                  legend = "right",
                  ncol = 1)
ggsave("plots/lambda_morin_fixed.png",
       units = "in",
       height = 6,
       width = 10)

# gradient plots ----------------------------------------------------------
results |>
  filter(bias_level == "minimal",
         cutoff == 0.584) |>
  ggplot(aes(x = env_gradient,
             y = lambda_under,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.25,
              linewidth = 0.25,
              color = "#FF1984") +
  facet_grid(. ~ known_beta,
             labeller = label_value) +
  theme_bw() +
  labs(title = "Minimally Biased Data")

results |>
  filter(bias_level == "extreme",
         cutoff == 1.41) |>
  ggplot(aes(x = env_gradient,
             y = lambda_under,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.1,
              linewidth = 0.1,
              color = "#FF1984") +
  facet_grid(. ~ known_beta,
             labeller = label_value) +
  theme_bw() +
  labs(title = "Extreme Biased data")

results |>
  filter(bias_level == "minimal") |>
  ggplot(aes(x = env_gradient,
             y = lambda_trimmed,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.25,
              linewidth = 0.25,
              color = "#FF914A") +
  facet_grid(cutoff ~ known_beta,
             labeller = label_value) +
  theme_bw() +
  labs(title = "Censored data, Minimal Bias")

results |>
  filter(bias_level == "extreme") |>
  ggplot(aes(x = env_gradient,
             y = lambda_trimmed,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.25,
              linewidth = 0.25,
              color = "#FF914A") +
  facet_grid(cutoff ~ known_beta,
             labeller = label_value) +
  theme_bw() +
  labs(title = "Censored data, Extreme Bias")


results |>
  filter(bias_level == "moderate") |>
  ggplot(aes(x = env_gradient,
             y = lambda_trimmed,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.25,
              linewidth = 0.25,
              color = "#FF914A") +
  facet_grid(cutoff ~ known_beta,
             labeller = label_value) +
  theme_bw() +
  labs(title = "Censored data, Moderate Bias")

results |>
  filter(bias_level == "strong") |>
  ggplot(aes(x = env_gradient,
             y = lambda_trimmed,
             group = rep)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE,
              alpha = 0.25,
              linewidth = 0.25,
              color = "#FF914A") +
  facet_grid(cutoff ~ known_beta,
             labeller = label_value) +
  theme_bw() +
  labs(title = "Censored data, Strong Bias")

# plot rep lines ----------------------------------------------------------

results |>
  filter(rep %in% 1:100,
         bias_level == "extreme") |>
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
  scale_color_manual(values = c(c("#FF914A",
                                  "#FF1984"))) +
  theme_bw() +
  geom_line(inherit.aes = FALSE,
            aes(x = env_gradient,
                y = known_lambda),
            color = "black",
            linetype = "dashed",
            linewidth = 1) +
  labs(title = "Extreme Bias")
ggsave("plots/rep_lines_morin_fixed_extreme_MS.png",
       units = "in",
       height = 6,
       width = 10)

results |>
  filter(rep %in% 1:100,
         bias_level == "minimal") |>
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
  scale_color_manual(values = c(c("#FF914A",
                                  "#FF1984"))) +
  theme_bw() +
  geom_line(inherit.aes = FALSE,
            aes(x = env_gradient,
                y = known_lambda),
            color = "black",
            linetype = "dashed",
            linewidth = 1)+
  labs(title = "Minimal Bias")
ggsave("plots/rep_lines_morin_fixed_minimal_MS.png",
       units = "in",
       height = 6,
       width = 10)

results |>
  filter(rep %in% 1:100,
         bias_level == "moderate") |>
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
  scale_color_manual(values = c(c("#FF914A",
                                  "#FF1984"))) +
  theme_bw() +
  geom_line(inherit.aes = FALSE,
            aes(x = env_gradient,
                y = known_lambda),
            color = "black",
            linetype = "dashed",
            linewidth = 1)+
  labs(title = "Moderate Bias")
ggsave("plots/rep_lines_morin_fixed_moderate_MS.png",
       units = "in",
       height = 6,
       width = 10)

results |>
  filter(rep %in% 1:100,
         bias_level == "strong") |>
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
  scale_color_manual(values = c(c("#FF914A",
                                  "#FF1984"))) +
  theme_bw() +
  geom_line(inherit.aes = FALSE,
            aes(x = env_gradient,
                y = known_lambda),
            color = "black",
            linetype = "dashed",
            linewidth = 1)+
  labs(title = "Strong Bias")
ggsave("plots/rep_lines_morin_fixed_strong_MS.png",
       units = "in",
       height = 6,
       width = 10)

# lm models ---------------------------------------------------------------

lm_models <- results |>
  select(group,
         lambda_under,
         lambda_trimmed,
         bias_level,
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
           bias_level,
           cutoff,
           rep,
           known_beta,
           data_model) |>
  add_count() |>
  filter(n ==5) |>
  nest() |>
  mutate(lm_model = map(
    data, ~lm(
      lambda_est ~ env_gradient, .))) |>
  mutate(coefs = map(
    lm_model, tidy)) |>
  unnest(coefs) 


# AUC calculations --------------------------------------------------------


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
  filter(term == "env_gradient",
         bias_level == "minimal") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Minimal Bias",
       x = "beta estimate",
       subtitle = "Fixed cutoffs")
ggsave("plots/beta_morin_fixed_minimal_MS.png",
       units = "in",
       height = 6,
       width = 10)


lm_models |>
  filter(term == "env_gradient",
         #cutoff != 0.001,
         bias_level == "extreme") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Extreme Bias",
       x = "beta estimate",
       subtitle = "Fixed cutoffs")
ggsave("plots/beta_morin_fixed_extreme_MS.png",
       units = "in",
       height = 6,
       width = 10)


lm_models |>
  filter(term == "env_gradient",
         #cutoff != 0.001,
         bias_level == "moderate") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~cutoff,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Moderate Bias",
       x = "beta estimate",
       subtitle = "Fixed cutoffs")
ggsave("plots/beta_morin_fixed_moderate_MS.png",
       units = "in",
       height = 6,
       width = 10)

lm_models |>
  filter(term == "env_gradient",
         #cutoff != 0.001,
         bias_level == "strong") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(cutoff~known_beta,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Strong Bias",
       x = "beta estimate",
       subtitle = "Fixed cutoffs")
ggsave("plots/beta_morin_fixed_strong_MS.png",
       units = "in",
       height = 6,
       width = 10)

# appropriate cutoffs -----------------------------------------------------


lm_models |>
  filter(term == "env_gradient",
         cutoff == 3.63,
         bias_level == "strong") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(cutoff~known_beta,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Strong Bias",
       x = "beta estimate",
       subtitle = "95% retention probability cutoff")
ggsave("plots/beta_morin_fixed_strong_95_MS.png",
       units = "in",
       height = 6,
       width = 10)

lm_models |>
  filter(term == "env_gradient",
         cutoff == 0.584,
         bias_level == "minimal") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(cutoff~known_beta,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Minimal Bias",
       x = "beta estimate",
       subtitle = "95% retention probability cutoff")
ggsave("plots/beta_morin_fixed_minimal_95_MS.png",
       units = "in",
       height = 6,
       width = 10)
lm_models |>
  filter(term == "env_gradient",
         cutoff == 1.41,
         bias_level == "moderate") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(cutoff~known_beta,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Moderate Bias",
       x = "beta estimate",
       subtitle = "95% retention probability cutoff")
ggsave("plots/beta_morin_fixed_moderate_95_MS.png",
       units = "in",
       height = 6,
       width = 10)

lm_models |>
  filter(term == "env_gradient",
         cutoff == 10.1,
         bias_level == "extreme") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(cutoff~known_beta,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Extreme Bias",
       x = "beta estimate",
       subtitle = "95% retention probability cutoff")
ggsave("plots/beta_morin_fixed_extreme_95_MS.png",
       units = "in",
       height = 6,
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
