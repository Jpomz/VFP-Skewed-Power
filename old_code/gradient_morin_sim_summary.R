#gradient sim summary

library(tidyverse)
library(broom)
library(tidybayes)
source("custom_functions.R")

results <- readRDS("simulation_results/gradient_sim_run_morin.rds")

results <- as_tibble(results) 
results

# how many parameter sets?
distinct(results, original_n, group, bias_level)

results |>
  select(M, original_n, under_n) |>
  mutate(prop = under_n / original_n) |>
  group_by(M) |>
  summarize(med_prop = median(prop),
            mean_prop = mean(prop), 
            sd_prop = sd(prop))


# N comparisons -----------------------------------------------------------

results |>
  select(known_lambda, 
         bias_level,
         original_n,
         under_n, 
         trimmed_n) |>
  filter(original_n == 5000,
         known_lambda == -1.25 |
           known_lambda == -1.5 |
           known_lambda == -1.75) |>
  pivot_longer(under_n:trimmed_n) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.5, 
               normalize = "groups")+
  scale_x_log10(guide = "axis_logticks") +
  geom_vline(aes(xintercept = original_n)) +
  facet_grid(bias_level~known_lambda,
             scales = "free")

results |>
  select(known_lambda, 
         bias_level,
         original_n,
         under_n, 
         trimmed_n) |>
  pivot_longer(under_n:trimmed_n) |>
  group_by(known_lambda, 
           bias_level,
           original_n,
           name) |>
  filter(name == "trimmed_n") |>
  summarize(q10 = 
              quantile(value,
                       probs = c(0.1)),
            q50 = 
              quantile(value,
                       probs = c(0.5)),
            q90 = 
              quantile(value,
                       probs = c(0.9)))


results |>
  filter(original_n == 5000,
         known_lambda == -1.25 |
           known_lambda == -1.5 |
           known_lambda == -1.75) |>
  select(known_lambda, 
         bias_level,
         original_n,
         under_n, 
         trimmed_n) |>
  #pivot_longer(under_n:trimmed_n) |>
  group_by(bias_level) |>
  summarize(med_under = median(under_n),
            med_xmin = median(trimmed_n))

# estimated xmins ---------------------------------------------------------

results |>
  filter(original_n == 5000, 
         known_lambda == -1.25 |
           known_lambda == -1.5 |
           known_lambda == -1.75) |>
  ggplot(aes(x = est_xmin,
             fill = as.factor(known_lambda))) +
  stat_halfeye(alpha = 0.5, 
               normalize = "groups")+
  facet_grid(bias_level~known_lambda,
             scales = "free") +
  scale_x_log10(guide = "axis_logticks") +
  scale_fill_viridis_d(option = "mako",
                       end = 0.8) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = expression(Estimated~x[min]),
       y = "")
ggsave("plots/xmin_dist_n5000.png",
       units = "in",
       height = 6,
       width = 10)
# xmin gets bigger with more original_n
# except in shallow lambdas - approx. same center, but wider distribution with smaller n
results |>
  filter(#original_n == 1000, 
    known_lambda == -1 |
      known_lambda == -1.5 |
      known_lambda == -2) |>
  ggplot(aes(x = est_xmin,
             fill = as.factor(original_n))) +
  stat_halfeye(alpha = 0.5, 
               normalize = "groups")+
  facet_grid(bias_level~known_lambda,
             scales = "free") +
  scale_x_log10(guide = "axis_logticks")
results |>
  filter(#original_n == 1000, 
    known_lambda == -1 |
      known_lambda == -1.5 |
      known_lambda == -2) |>
  ggplot(aes(x = est_xmin,
             fill = as.factor(original_n))) +
  stat_halfeye(alpha = 0.5, 
               normalize = "groups")+
  facet_grid(known_lambda~bias_level,
             scales = "free") +
  scale_x_log10(guide = "axis_logticks")

results |>
  filter(#original_n == 1000, 
    known_lambda == -1.25 |
      known_lambda == -1.5 |
      known_lambda == -1.75) |>
  ggplot(aes(x = est_xmin,
             fill = as.factor(original_n))) +
  stat_halfeye(alpha = 0.5, 
               normalize = "groups")+
  facet_grid(known_lambda~bias_level,
             scales = "free") +
  scale_x_log10(guide = "axis_logticks")


# why are some of the xmins bimodal?
results |>
  filter(known_lambda == -1.5,
    original_n == 10000, 
      bias_level == "extreme" ) |>
  ggplot(aes(x = est_xmin,
             fill = as.factor(original_n))) +
  stat_halfeye(alpha = 0.5, 
               normalize = "groups")+
  facet_grid(bias_level~known_lambda,
             scales = "free") +
  scale_x_log10(guide = "axis_logticks")

# "missing" values from ~0.6 to 0.7
results |>
  filter(known_lambda == -1.5,
         original_n == 10000, 
         bias_level == "extreme",
         est_xmin>0.49, est_xmin < 0.89) |>
  ggplot(aes(y = est_xmin,
             x = trimmed_n)) +
  geom_point()

# missing masses
sizeSpectra::massToLength(c(0.5, 0.6, 0.7, 0.8),
                          LWa = 0.0064,
                          LWb = 2.788)
# probability of sampling those masses with Mesh = 1 ("extreme")
plogis(morin_ln_p(L = c(4.77, 5.10, 5.39, 5.65),
                  M = 1))
# 75-82% chance of being sampled...
plot_morin_bias(n = 10000,
                lambda = -1.5,
                M = 0.25,
                binwidth = 0.05) +
  geom_vline(aes(xintercept = 0.55))

# trimmed_n vs est_xmin
results |>
  filter(bias_level == "extreme") |>
  ggplot(aes(y = est_xmin,
             x = trimmed_n,
             color = known_lambda)) +
  geom_point() +
  facet_grid(known_lambda~original_n,
             scales = "free")

results |>
  filter(bias_level == "extreme",
         est_xmin < 2.5,
         known_lambda == -1.5) |>
  ggplot(aes(y = est_xmin,
             x = trimmed_n,
             color = known_lambda)) +
  geom_point() +
  facet_grid(known_lambda~original_n,
             scales = "free") +
  theme_classic() +
  theme(legend.position = "none")
ggsave("plots/xmin_gaps_SI.png",
       units = "in",
       height = 6,
       width = 10)


results |>
  filter(known_lambda == -1.5, 
         bias_level == "extreme",
         original_n == 5000) |>
  ggplot(aes(y = est_xmin,
             x = trimmed_n,
             color = original_n)) +
  geom_point() +
  facet_grid(.~original_n,
             scales = "free")

# lambdas -----------------------------------------------------------------
results |>
  filter(original_n == 5000,
         known_lambda == -1.25 |
           known_lambda == -1.5 |
           known_lambda == -1.75) |>
  select(known_lambda, 
         bias_level,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(alpha = 0.6,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#019AFF"))) +
  facet_grid(bias_level ~ known_lambda,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "")


lambda_xmin_125_panel <- results |>
  filter(original_n == 5000,
         known_lambda == -1.25) |>
  select(known_lambda, 
         bias_level,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(alpha = 0.6,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#019AFF"))) +
  facet_grid(known_lambda ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "") +
  scale_x_continuous(n.breaks = 3)

lambda_xmin_15_panel <- results |>
  filter(original_n == 5000,
         known_lambda == -1.5) |>
  select(known_lambda, 
         bias_level,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(alpha = 0.6,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#019AFF"))) +
  facet_grid(known_lambda ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "") +
  scale_x_continuous(n.breaks = 3)

lambda_xmin_175_panel <- results |>
  filter(original_n == 5000,
         known_lambda == -1.75) |>
  select(known_lambda, 
         bias_level,
         lambda_under,
         lambda_trimmed) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(alpha = 0.6,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#019AFF"))) +
  facet_grid(known_lambda ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "")+
  scale_x_continuous(n.breaks = 3)

ggpubr::ggarrange(lambda_xmin_125_panel,
                  lambda_xmin_15_panel,
                  lambda_xmin_175_panel,
                  common.legend = TRUE,
                  legend = "right",
                  ncol = 1)
ggsave("plots/lambda_morin_xmin.png",
       units = "in",
       height = 6,
       width = 10)

# plot rep lines ----------------------------------------------------------

results |>
  filter(rep %in% 1:100,
         original_n == 5000) |>
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
  facet_grid(known_beta~bias_level) +
  scale_color_manual(values = c(c("#019AFF",
                                  "#FF1984"))) +
  theme_bw() +
  geom_line(inherit.aes = FALSE,
            aes(x = env_gradient,
                y = known_lambda),
            color = "black",
            linetype = "dashed",
            linewidth = 1) +
  labs(title = "N = 5000")
ggsave("plots/rep_lines_morin_clauset_MS.png",
       units = "in",
       height = 6,
       width = 10)




lm_models <- results |>
  select(group,
         lambda_under,
         lambda_trimmed,
         bias_level,
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
           bias_level,
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

# lm_models <- lm_models |>
#   mutate(`Bias level` = case_when(
#     bias_level == "01" ~ "Minimal",
#     bias_level == "02" ~ "Moderate",
#     bias_level == "03" ~ "Strong",
#     bias_level == "04" ~ "Extreme"
#   ),
#   `Bias level` = factor(`Bias level`,
#                         levels = c("Minimal",
#                                    "Moderate", 
#                                    "Strong",
#                                    "Extreme")))

# beta by bias_level ---------------------------------------------------------
beta_plot_dat <- lm_models |>
  filter(original_n == 5000)
saveRDS(beta_plot_dat,
        "simulation_summaries/beta_xmin_plot_data.RDS")

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
                                 "#019AFF"))) +
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
                                 "#019AFF"))) +
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
                                 "#019AFF"))) +
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
                                 "#019AFF"))) +
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
ggsave("plots/beta_clauset.png",
       units = "in",
       height = 6,
       width = 10)

# AUC ---------------------------------------------------------------------


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
272 / 12000
228 / 12000
# ~ 2-2.5% had false positive relationships

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
5032 / 36000 # 14% - trimmed
453 / 36000 # 1% - under

lm_models |>
  filter(term == "env_gradient",
         p.value < 0.05,
         known_beta == 0,
         original_n == 5000) |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(bias_level~group) +
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
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(bias_level~known_beta,
             labeller = 
               labeller(bias_level = label_value,
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
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(bias_level~known_beta,
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
  facet_grid(bias_level~known_beta,
             labeller = 
               labeller(bias_level = label_value,
                        known_beta = label_both)) +
  labs(title = "True Negative Relationship") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984")))

lm_models |>
  filter(original_n == 5000,
         term == "env_gradient") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~bias_level,
             scales = "free") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  theme_bw()
ggsave("plots/beta_morin_clauset_MS.png",
       units = "in",
       height = 6,
       width = 10)

lm_models |>
  filter(original_n == 5000,
         term == "env_gradient",
         bias_level == "minimal") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(.~known_beta,
             scales = "free") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Minimal Bias",
       subtitle = "Clauset x_min",
       x = expression(beta ~ "estimate"))

ggsave("plots/beta_morin_clauset_minimal_MS.png",
       units = "in",
       height = 6,
       width = 10)

lm_models |>
  filter(original_n == 5000,
         term == "env_gradient",
         bias_level == "moderate") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(.~known_beta,
             scales = "free") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Moderate Bias",
       subtitle = "Clauset x_min",
       x = expression(beta ~ "estimate"))

ggsave("plots/beta_morin_clauset_moderate_MS.png",
       units = "in",
       height = 6,
       width = 10)

lm_models |>
  filter(original_n == 5000,
         term == "env_gradient",
         bias_level == "strong") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(.~known_beta,
             scales = "free") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Strong Bias",
       subtitle = "Clauset x_min",
       x = expression(beta ~ "estimate"))

ggsave("plots/beta_morin_clauset_strong_MS.png",
       units = "in",
       height = 6,
       width = 10)

lm_models |>
  filter(original_n == 5000,
         term == "env_gradient",
         bias_level == "minimal") |>
  ggplot(aes(x = estimate,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(.~known_beta,
             scales = "free") +
  scale_fill_manual(values = c(c("#019AFF",
                                 "#FF1984"))) +
  theme_bw() +
  labs(title = "Minimal Bias",
       subtitle = "Clauset x_min",
       x = expression(beta ~ "estimate"))

ggsave("plots/beta_morin_clauset_minimal_MS.png",
       units = "in",
       height = 6,
       width = 10)

lm_models |>
  filter(original_n == 5000,
         term == "env_gradient",
         bias_level == "extreme") |>
  select(data_model, estimate, bias_level, known_beta) |>
  pivot_wider(values_from = estimate, names_from = data_model) |>
  ggplot(aes(x = lambda_trimmed,
             y = lambda_under)) +
  geom_point() +
  facet_grid(bias_level~known_beta,
             scales = "free")
lm_models |>
  filter(term == "env_gradient",
         bias_level == "extreme") |>
  select(data_model, estimate, bias_level, known_beta) |>
  pivot_wider(values_from = estimate, names_from = data_model) |>
  ggplot(aes(x = lambda_trimmed,
             y = lambda_under)) +
  geom_point() +
  facet_grid(original_n~known_beta,
             scales = "free")

lm_models |>
  filter(term == "env_gradient",
         original_n == "5000",
         bias_level == "extreme") |>
  select(data_model, estimate, bias_level, known_beta) |>
  pivot_wider(values_from = estimate, names_from = data_model) |>
  ggplot(aes(x = lambda_trimmed,
             y = lambda_under)) +
  geom_point() +
  facet_grid(original_n~known_beta,
             scales = "free")

lm_models |>
  filter(term == "env_gradient",
         original_n == "5000",
         bias_level == "strong") |>
  select(data_model, estimate, bias_level, known_beta) |>
  pivot_wider(values_from = estimate, names_from = data_model) |>
  ggplot(aes(x = lambda_trimmed,
             y = lambda_under)) +
  geom_point() +
  facet_grid(original_n~known_beta,
             scales = "free") +
  stat_smooth(method = "lm")


lm_models |>
  filter(term == "env_gradient",
         original_n == "5000",
         bias_level == "strong") |>
  select(data_model, estimate, bias_level, known_beta) |>
  pivot_wider(values_from = estimate, names_from = data_model) |>
  group_by(known_beta) |>
  summarize(cor(lambda_under, lambda_trimmed))

# bias tables -------------------------------------------------------------

lm_table_dat <- lm_models |>
  ungroup() |>
  select(group, bias_level, rep, original_n,
         known_beta, data_model, term, 
         estimate, std.error) |>
  filter(term == "env_gradient")
  
lm_table_dat |>
  distinct(original_n, known_beta, bias_level)

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
           bias_level, 
           data_model) |>
  add_count() |>
  #group_by(target_name, name, n) %>%
  summarize(
    median_ci_range = median(conf_width,
                             na.rm = TRUE),
    median_abs_bias = median(abs_bias,
                             na.rm = TRUE),
    sd_abs_bias = sd(abs_bias,
                     na.rm = TRUE),
    n = n()) |>
  arrange(original_n,
          known_beta,
          bias_level,
          data_model) 
write_csv(beta_bias, 
          "simulation_summaries/beta_bias_morin_clauset_table.csv")
