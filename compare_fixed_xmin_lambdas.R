# comparing estimates for fixed and xmin cutoffs

library(tidyverse)
library(broom)
library(tidybayes)


# Read in data ------------------------------------------------------------


fixed_cut <- readRDS("simulation_results/fixed_cut_morin_sim_run.rds")
xmin_cut <- readRDS("simulation_results/gradient_sim_run_morin.rds")

fixed_cut <- as_tibble(fixed_cut)
xmin_cut <- as_tibble(xmin_cut)
dim(fixed_cut)  
dim(xmin_cut)

names(fixed_cut)
names(xmin_cut)


# wrangle data ------------------------------------------------------------


# make a variable to filter on
# only want to keep the ~95% cutoff for each Mesh size
cut_mass <- tibble(cutoff = unique(fixed_cut$cutoff),
                   M = c(0.125, 0.25, 0.5, 1),
                   cutoff_mass = sizeSpectra::lengthToMass(cutoff, 
                                                           LWa = 0.0064,
                                                           LWb = 2.788))
fixed_cols <- fixed_cut |>
  select(known_lambda, 
         lambda_trimmed,
         cutoff, 
         M,
         rep,
         group,
         bias_level,
         env_gradient) |>
  # mutate(keep = case_when(
  #   cutoff == 0.584 & M == 0.125 ~ TRUE,
  #   cutoff == 1.410 & M == 0.25 ~ TRUE,
  #   cutoff == 3.63 & M == 0.5 ~ TRUE,
  #   cutoff == 10.1 & M == 1 ~ TRUE,
  #   .default = FALSE
  # )) |> 
  # filter(keep == TRUE) |>
  as_tibble() |>
  rename(fixed_lambda = lambda_trimmed) |>
  select(known_lambda, fixed_lambda, rep,
         group,
         bias_level,
         env_gradient)

xmin_cols <- xmin_cut |>
  as_tibble() |>
  # filter(original_n == 10000) |>
  select(known_lambda, 
         lambda_trimmed,
         rep,
         group,
         bias_level,
         env_gradient) |>
  rename(xmin_lambda = lambda_trimmed) |>
  select(known_lambda, xmin_lambda, rep,
         group,
         bias_level,
         env_gradient)

distinct(fixed_cols, rep)
distinct(xmin_cols, rep)
distinct(fixed_cols, known_lambda)
distinct(xmin_cols, known_lambda)

# xmin_cols |>
#   filter(known_lambda == -1.5, 
#          group == "A", 
#          rep == 1)
# fixed_cols |>
#   filter(known_lambda == -1.5, 
#          group == "A", 
#          rep == 1)


# xmin distribution -------------------------------------------------------


xmin_cut |>
  filter(known_lambda == -1.9 |
           known_lambda == -2 |
           known_lambda == -2.1) |>
  ggplot(aes(x = est_xmin,
             fill = as.factor(known_lambda))) +
  stat_halfeye(alpha = 0.5, 
               normalize = "groups")+
  geom_vline(aes(xintercept = 0.005),
             linetype = "dashed") +
  facet_grid(known_lambda~bias_level,
             scales = "free") +
  scale_x_log10(guide = "axis_logticks",
                n.breaks = 3) +
  scale_fill_viridis_d(option = "mako",
                       end = 0.8) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = expression(Estimated~x[min]),
       y = "") 
ggsave("plots/xmin_dist_n10000.png",
       units = "in",
       height = 6,
       width = 10)

# Lambdas -----------------------------------------------------------------


both_lambdas <- left_join(xmin_cols, fixed_cols)

both_lambdas |>
  filter(known_lambda == -1.9 | known_lambda == -2 | known_lambda == -2.1) |>
  ggplot(aes(x = xmin_lambda,
             y = fixed_lambda,
             color = known_lambda)) +
  geom_point() +
  facet_grid(bias_level~known_lambda,
             scales = "free")

both_lambdas |>
  filter(!is.na(fixed_lambda),
         !is.na(xmin_lambda)) |>
  group_by(known_lambda) |>
  summarise(cor(xmin_lambda, fixed_lambda))

# all three lambda estimates ----------------------------------------------

fixed_under_lambda <- fixed_cut |>
  select(known_lambda, 
         lambda_trimmed,
         lambda_under,
         cutoff, 
         M,
         rep,
         group,
         bias_level,
         env_gradient) |>
  # mutate(keep = case_when(
  #   cutoff == 0.584 & M == 0.125 ~ TRUE,
  #   cutoff == 1.410 & M == 0.25 ~ TRUE,
  #   cutoff == 3.63 & M == 0.5 ~ TRUE,
  #   cutoff == 10.1 & M == 1 ~ TRUE,
  #   .default = FALSE
  # )) |> 
  # filter(keep == TRUE) |>
  as_tibble() |>
  rename(fixed_lambda = lambda_trimmed) |>
  select(known_lambda,
         fixed_lambda,
         lambda_under,
         bias_level) |>
  rename(`Biased` = lambda_under,
         `Fixed` = fixed_lambda) |>
  pivot_longer(`Biased`:`Fixed`)

xmin_lambda <- xmin_cols |>
  select(known_lambda, xmin_lambda, bias_level) |>
  rename(xmin = xmin_lambda) |>
  pivot_longer(xmin)

three_lambdas <- bind_rows(xmin_lambda, fixed_under_lambda) |>
  mutate(name = factor(name, 
                             levels = c(
                               "Fixed",
                               "xmin", 
                               "Biased")))


lambda_19_panel <- three_lambdas |>
  filter(known_lambda == -1.9) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF914A",
                                 
                                 "#019AFF",
                                 "#FF1984"))) +
  facet_grid(known_lambda ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "") +
  scale_x_continuous(n.breaks = 4)

lambda_20_panel <- three_lambdas |>
  filter(known_lambda == -2) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF914A",
                                 
                                 "#019AFF",
                                 "#FF1984"))) +
  facet_grid(known_lambda ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "") +
  scale_x_continuous(n.breaks = 4)

lambda_21_panel <- three_lambdas |>
  filter(known_lambda == -2.1) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF914A",
                                 
                                 "#019AFF",
                                 "#FF1984"))) +
  facet_grid(known_lambda ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "") +
  scale_x_continuous(n.breaks = 4)

# lambda_2_panel <- three_lambdas |>
#   filter(known_lambda == -2) |>
#   ggplot(aes(x = value,
#              y = name,
#              fill = name)) +
#   stat_halfeye(normalize = "groups") +
#   geom_vline(aes(xintercept = known_lambda),
#              linetype = "dashed") +
#   scale_fill_manual(values = c(c("#FF914A",
#                                  
#                                  "#019AFF",
#                                  "#FF1984"))) +
#   facet_grid(known_lambda ~ bias_level,
#              labeller = label_value,
#              scales = "free") +
#   theme_bw() +
#   labs(x = "\u03bb estimate",
#        y = "") +
#   scale_x_continuous(n.breaks = 4)

ggpubr::ggarrange(lambda_19_panel,
                  lambda_20_panel,
                  lambda_21_panel,
                  common.legend = TRUE,
                  legend = "right",
                  ncol = 1)
ggsave("plots/lambdas_all_three.png",
       units = "in",
       height = 10,
       width = 10)

three_lambdas |>
  filter(known_lambda == -1.9 |
           known_lambda == -2 |
           known_lambda == -2.1) |>
  mutate(
    name = factor(name, levels = c("Biased", "xmin", "Fixed")), 
    abs_delta = abs(known_lambda - value)) |>
  group_by(known_lambda, bias_level, name)|>
  summarise(median_abs_delta = median(abs_delta),
            sd_abs_delta = sd(abs_delta)) |>
  write_csv("simulation_summaries/three_lambdas.csv")


three_lambdas |>
  filter(known_lambda == -1.9 |
           known_lambda == -2 |
           known_lambda == -2.1) |>
  ungroup() |>
  mutate(abs_delta = abs(known_lambda - value)) |>
  ggplot(aes(x = abs_delta,
             y = name,
             fill = name)) +
  stat_halfeye(normalize = "groups") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#019AFF",
                                 "#FF1984"))) +
  facet_grid(known_lambda ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw()+
  labs(x = "\u0394 in \u03bb estimates",
       y = "") +
  scale_x_continuous(n.breaks = 4)
ggsave("plots/lambda_deltas_all_three.png",
       units = "in",
       height = 10,
       width = 10)



# Sample sizes ------------------------------------------------------------

n_fixed <- fixed_cut |>
  select(known_lambda, 
         under_n,
         trimmed_n,
         cutoff, 
         M,
         rep,
         group,
         bias_level,
         env_gradient) |>
  as_tibble() |>
  rename(Fixed = trimmed_n,
         Biased = under_n)|>
  select(known_lambda:Fixed, 
         bias_level) |>
  pivot_longer(Biased:Fixed)

n_xmin <- xmin_cut |>
  as_tibble() |>
  select(known_lambda, 
         trimmed_n,
         rep,
         group,
         bias_level,
         env_gradient) |>
  rename(value = trimmed_n) |>
  select(known_lambda,
         bias_level, 
         value) |>
  mutate(name = "xmin")

three_ns <- bind_rows(n_fixed, n_xmin) |>
  mutate(name = factor(name, 
                       levels = c("Fixed", "xmin", "Biased"))) |>
  filter(known_lambda == -1.9 |
           known_lambda == -2 |
           known_lambda == -2.1) |>
  mutate(known_lambda = factor(known_lambda,
                               levels = c("-1.9", 
                                          "-2", 
                                          "-2.1")))

three_ns  |>
  ggplot(aes(x = value, 
             y = name,
             fill = name)) +
  stat_halfeye(normalize = "groups") +
  facet_grid(known_lambda~bias_level,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#019AFF",
                                 "#FF1984"))) +
  theme_bw() +
  labs(x = "Sample Size, N",
       y = "") +
  scale_x_continuous(n.breaks = 4)
  
ggsave("plots/N_all_three.png",
       units = "in",
       height = 10,
       width = 10)

three_ns |>
  group_by(known_lambda, bias_level, name) |>
  summarize(med_n = median(value),
            sd_n = sd(value)) |>
  write_csv("simulation_summaries/three_ns.csv")

three_ns |>
  filter(bias_level == "moderate",
         name == "Fixed") |>
  group_by(known_lambda, bias_level, name) |>
  summarize(median(value))

# plot comparing xmin and fixed cutoffs beta ------------------------------

# fixed_cut lm models?
fixed_lm <- fixed_cut |>
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

fixed_beta <- fixed_lm

# xmin lm-Models
xmin_lm <- xmin_cut |>
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

xmin_beta <- xmin_lm

fixed <- fixed_beta |>
  filter(term == "env_gradient") |>
  mutate(data_model = case_when(data_model == "lambda_trimmed"~"Fixed",
                                .default = "Biased"))
xmin <- xmin_beta |>
  filter(term == "env_gradient",
         data_model == "lambda_trimmed") |>
  mutate(data_model = case_when(data_model == "lambda_trimmed"~"xmin",
                                .default = data_model))

betas <- bind_rows(fixed, xmin) |>
  mutate(data_model = factor(data_model, 
                             levels = c(
                                        "Fixed",
                                        "xmin", 
                                        "Biased")))
# 0, -0.05, -0.1, -0.2
beta_00_panel <- betas |>
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
  scale_fill_manual(values = c(c("#FF914A",
                                 "#019AFF",
                                 "#FF1984"))) +
  facet_grid(known_beta ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "")+
  scale_x_continuous(n.breaks = 4)

beta_05_panel <- betas |>
  filter(known_beta == -0.05,
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
  scale_fill_manual(values = c(c("#FF914A",
                                 "#019AFF",
                                 "#FF1984"))) +
  facet_grid(known_beta ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "")+
  scale_x_continuous(n.breaks = 4)

beta_10_panel <- betas |>
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
  scale_fill_manual(values = c(c("#FF914A",
                                 "#019AFF",
                                 "#FF1984"))) +
  facet_grid(known_beta ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "",
       y = "")+
  scale_x_continuous(n.breaks = 4)

beta_20_panel <- betas |>
  filter(known_beta == -0.2,
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
  scale_fill_manual(values = c(c("#FF914A",
                                 "#019AFF",
                                 "#FF1984"))) +
  facet_grid(known_beta ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03b2 estimate",
       y = "")+
  scale_x_continuous(n.breaks = 4)

ggpubr::ggarrange(beta_00_panel,
                  beta_05_panel,
                  beta_10_panel,
                  beta_20_panel,
                  common.legend = TRUE,
                  legend = "right",
                  ncol = 1)
ggsave("plots/betas_all_three.png",
       units = "in",
       height = 10,
       width = 10)

# beta bias table ####
betas |>
  ungroup() |>
  mutate(abs_delta = abs(known_beta - estimate)) |>
  group_by(known_beta, bias_level, data_model) |>
  summarise(median_abs_delta = median(abs_delta),
            sd_abs_delta = sd(abs_delta))|>
  write_csv("simulation_summaries/three_betas.csv")


betas |>
  ungroup() |>
  mutate(abs_delta = abs(known_beta - estimate)) |>
  ggplot(aes(x = abs_delta,
             y = data_model,
             fill = data_model)) +
  stat_halfeye(normalize = "groups") +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#019AFF",
                                 "#FF1984"))) +
  facet_grid(known_beta ~ bias_level,
             labeller = label_value,
             scales = "free") +
  theme_bw()+
  labs(x = "\u0394 in \u03b2 estimates",
       y = "") +
  scale_x_continuous(n.breaks = 4)
ggsave("plots/beta_deltas_all_three.png",
       units = "in",
       height = 10,
       width = 10)


# rep lines ---------------------------------------------------------------

fixed_cut |>
  filter(rep %in% c(1:100, 111, 314, 330, 345, 395, 453, 477, 487)) |>
  # numbers are reps which actually have lambda estimates for all 5 "sites"
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
  scale_color_manual(values = c(c("#FF914A",
                                  "#FF1984"))) +
  theme_bw() +
  geom_line(inherit.aes = FALSE,
            aes(x = env_gradient,
                y = known_lambda),
            color = "black",
            linetype = "dashed",
            linewidth = 1) +
  labs(title = "Fixed Cutoffs")
ggsave("plots/rep_lines_morin_fixed_SI.png",
       units = "in",
       height = 6,
       width = 10)

xmin_cut |>
  filter(rep %in% c(1:100)) |>
  # numbers are reps which actually have lambda estimates for all 5 "sites"
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
  scale_color_manual(values = c(c("#019AFF",
                                  "#FF1984"))) +
  theme_bw() +
  geom_line(inherit.aes = FALSE,
            aes(x = env_gradient,
                y = known_lambda),
            color = "black",
            linetype = "dashed",
            linewidth = 1) +
  labs(title = "x_min Cutoffs")
ggsave("plots/rep_lines_morin_xmin_SI.png",
       units = "in",
       height = 6,
       width = 10)

