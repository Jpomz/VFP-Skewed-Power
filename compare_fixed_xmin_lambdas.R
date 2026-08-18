# comparing estimates for fixed and xmin cutoffs

library(tidyverse)
library(broom)
library(tidybayes)


# Read in data ------------------------------------------------------------

# ~97.5% cutoff
# fixed_cut <- readRDS("simulation_results/fixed_cut_morin_sim_run.rds")

# 10xM cutoff
fixed_cut <- readRDS("simulation_results/fixed_cut_morin_10xM_sim_run.rds")
cutoff_df_10x <- readRDS("simulation_results/Morin_10x_cutoff.RDS")
xmin_cut <- readRDS("simulation_results/gradient_sim_run_morin.rds")

fixed_cut <- as_tibble(fixed_cut)
xmin_cut <- as_tibble(xmin_cut)
dim(fixed_cut)  
dim(xmin_cut)

names(fixed_cut)
names(xmin_cut)


# wrangle data ------------------------------------------------------------


# make a variable to filter on
# only want to keep the ~XX% cutoff for each Mesh size
# cut_mass <- tibble(
#   cutoff = unique(fixed_cut$cutoff),
#   M = c(0.125, 0.25, 0.5, 1),
#   cutoff_mass = sizeSpectra::lengthToMass(cutoff, 
#                                           LWa = 0.0064,
#                                           LWb = 2.788))
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

# xmin distribution -------------------------------------------------------
xmin_cut |>
  filter(known_lambda == -2.0 |
           known_lambda == -1.9 |
           known_lambda == -2.1) |>
  left_join(cutoff_df_10x) |>
  ggplot(aes(x = est_xmin,
             y = bias_level,
             fill = bias_level,
             group = bias_level)) +
  ggridges::geom_density_ridges(
    alpha = 0.5
  ) +
  scale_x_log10(guide = "axis_logticks",
                n.breaks = 3) +
  scale_fill_viridis_d(option = "mako",
                       end = 0.8) +
  theme_classic() +
  facet_grid(known_lambda~.) +
  theme(legend.position = "none") +
  labs(x = expression(Estimated~x[min]),
       y = "") 
## old version #
# xmin_cut |>
#   filter(known_lambda == -1.9 |
#            known_lambda == -2 |
#            known_lambda == -2.1) |>
#   left_join(cutoff_df_10x) |>
#   ggplot(aes(x = est_xmin,
#              fill = as.factor(known_lambda))) +
#   stat_halfeye(alpha = 0.5, 
#                normalize = "groups")+
#   geom_vline(aes(xintercept = cutoff_mass),
#              linetype = "dashed") +
#   facet_grid(known_lambda~bias_level) +
#   scale_x_log10(guide = "axis_logticks",
#                 n.breaks = 3) +
#   scale_fill_viridis_d(option = "mako",
#                        end = 0.8) +
#   theme_classic() +
#   theme(legend.position = "none") +
#   labs(x = expression(Estimated~x[min]),
#        y = "") 
ggsave("plots/xmin_dist_n10000_SI.png",
       units = "in",
       height = 6,
       width = 10)

# xmin_cut |>
#   filter(known_lambda == -2.0) |>
#   left_join(cutoff_df_10x) |>
#   ggplot(aes(x = est_xmin,
#              fill = bias_level)) +
#   stat_halfeye(alpha = 0.5, 
#                normalize = "groups")+
#   geom_vline(aes(xintercept = cutoff_mass),
#              linetype = "dashed") +
#   facet_grid(known_lambda~bias_level) +
#   scale_x_log10(guide = "axis_logticks",
#                 n.breaks = 3) +
#   scale_fill_viridis_d(option = "mako",
#                        end = 0.8) +
#   theme_classic() +
#   theme(legend.position = "none") +
#   labs(x = expression(Estimated~x[min]),
#        y = "") 

xmin_cut |>
  filter(known_lambda == -2.0) |>
  left_join(cutoff_df_10x) |>
  ggplot(aes(x = est_xmin,
             y = bias_level,
             fill = bias_level,
             group = bias_level)) +
  ggridges::geom_density_ridges(
    alpha = 0.5
  ) +
  scale_x_log10(guide = "axis_logticks",
                n.breaks = 3) +
  scale_fill_viridis_d(option = "mako",
                       end = 0.8) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = expression(Estimated~x[min]),
       y = "") 
ggsave("plots/xmin_dist_n10000_MS.png",
       units = "in",
       height = 6,
       width = 10)

# Lambdas -----------------------------------------------------------------


both_lambdas <- left_join(xmin_cols, fixed_cols)

both_lambdas |>
  filter(known_lambda == -1.9 | known_lambda == -2 | known_lambda == -2.1) |>
  filter(group == "D") |>
  ggplot(aes(x = xmin_lambda,
             y = fixed_lambda,
             color = known_lambda)) +
  geom_point() +
  facet_grid(bias_level~known_lambda,
             scales = "free")

both_lambdas |>
  filter(!is.na(fixed_lambda),
         !is.na(xmin_lambda),
         known_lambda == -1.9 | 
           known_lambda == -2 | 
           known_lambda == -2.1) |>
  group_by(known_lambda, bias_level) |>
  summarise(cor(xmin_lambda, fixed_lambda)) |>
  arrange(bias_level, known_lambda)

# all three lambda estimates ----------------------------------------------

fixed_under_lambda <- fixed_cut |>
  filter(group == "D") |>
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
  filter(group == "D") |>
  select(known_lambda, xmin_lambda, bias_level) |>
  rename(xmin = xmin_lambda) |>
  pivot_longer(xmin)

three_lambdas <- bind_rows(xmin_lambda, fixed_under_lambda) |>
  mutate(name = factor(name, 
                             levels = c(
                               "Fixed",
                               "xmin", 
                               "Biased")))

three_lambdas |>
  filter(known_lambda == -2.1 |
           known_lambda == -1.9 |
           known_lambda == -2) |>
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
             scale = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "") +
  scale_x_continuous(n.breaks = 4)

ggsave("plots/lambdas_all_three_SI.png",
       units = "in",
       height = 10,
       width = 10)

three_lambdas |>
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
             scale = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "") +
  scale_x_continuous(n.breaks = 4)
ggsave("plots/lambdas_one_panel_MS.png",
       units = "in",
       height = 10,
       width = 10)

lambdas_deviation <- three_lambdas |>
  filter(known_lambda == -1.9 |
           known_lambda == -2 |
           known_lambda == -2.1) |>
  mutate(
    data = factor(name, levels = c("Biased", "xmin", "Fixed")), 
    abs_delta = abs(known_lambda - value)) |>
  group_by(known_lambda, bias_level, data)|>
  summarise(
    n_reps = n(),
    median_abs_delta = median(abs_delta),
            sd_abs_delta = sd(abs_delta))

# plot of lambda delta
# silencing for now per JSW comment
# three_lambdas |>
#   filter(known_lambda == -1.9 |
#            known_lambda == -2 |
#            known_lambda == -2.1) |>
#   ungroup() |>
#   mutate(abs_delta = abs(known_lambda - value)) |>
#   ggplot(aes(x = abs_delta,
#              y = name,
#              fill = name)) +
#   stat_halfeye(normalize = "groups") +
#   scale_fill_manual(values = c(c("#FF914A",
#                                  "#019AFF",
#                                  "#FF1984"))) +
#   facet_grid(known_lambda ~ bias_level,
#              labeller = label_value,
#              scales = "free") +
#   theme_bw()+
#   labs(x = "\u0394 in \u03bb estimates",
#        y = "") +
#   scale_x_continuous(n.breaks = 4)
# ggsave("plots/lambda_deltas_all_three.png",
#        units = "in",
#        height = 10,
#        width = 10)

# lambda CIs ####
fixed_lambda_ci <- fixed_cut |>
  filter(group == "D") |>
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
fixed_lambda_ci |>
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
  labs(title = "Fixed Cutoff: 10x Mesh",
       y = "",
       x = "95% CI",
       color = "CI contains known value") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(n.breaks = 3)
ggsave("plots/fixed_10xM_lambda_CIs.png",
       units = "in",
       height = 10,
       width = 10)

fixed_lambda_ci_percent <- fixed_lambda_ci |>
  mutate(
    ci_width = lambda_trimmed_hi - lambda_trimmed_lo) |>
  group_by(known_lambda, M) |>
  summarize(n = n(), 
            prop_in = sum(in_fixed_ci),
            mean_ci_width = mean(ci_width),
            sd_ci_width = sd(ci_width)) |>
  mutate(percent_ci = prop_in / n,
         data = "Fixed") |>
  select(known_lambda,
         M,
         data,
         percent_ci,
         mean_ci_width, 
         sd_ci_width) 


# xmin and biased lambda CIs ####
xmin_lambda_ci <- xmin_cut |>
  filter(group == "D") |>
  select(known_lambda, lambda_trimmed:lambda_trimmed_hi, M) |>
  filter(known_lambda == -1.9|
           known_lambda == -2 |
           known_lambda == -2.1) |>
  mutate(in_xmin_ci = 
           known_lambda > lambda_trimmed_lo &
           known_lambda < lambda_trimmed_hi)

# xmin ci plot ####
xmin_lambda_ci |>
  group_by(M, known_lambda) |>
  arrange(lambda_trimmed_lo) |>
  mutate(y = (1:n())) |>
  ggplot(aes(x = lambda_trimmed, 
             xmin = lambda_trimmed_lo,
             xmax = lambda_trimmed_hi,
             y = y,
             color = in_xmin_ci)) +
  geom_pointrange(fatten = 0,
                  linewidth = 0.25,
                  alpha = 0.75) +
  geom_vline(aes(xintercept = known_lambda)) +
  facet_grid(known_lambda~M,
             scales = "free") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(title = "xmin Cutoff",
       y = "",
       x = "95% CI",
       color = "CI contains known value") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_x_continuous(n.breaks = 3)
ggsave("plots/xmin_lambda_CIs.png",
       units = "in",
       height = 10,
       width = 10)
# xmin percent CI table ####
xmin_lambda_ci_percent <- xmin_lambda_ci |>
  mutate(
    ci_width = lambda_trimmed_hi - lambda_trimmed_lo) |>
  group_by(known_lambda, M) |>
  summarize(n = n(), 
            prop_in = sum(in_xmin_ci),
            mean_ci_width = mean(ci_width),
            sd_ci_width = sd(ci_width)) |>
  mutate(percent_ci = prop_in / n,
         data = "xmin") |>
  select(known_lambda,
         M,
         data,
         percent_ci,
         mean_ci_width, 
         sd_ci_width) 

# biased lambda CI ####
biased_lambda_ci <- xmin_cut |>
  filter(group == "D") |>
  select(known_lambda, lambda_under:lambda_under_hi, M) |>
  filter(known_lambda == -1.9|
           known_lambda == -2 |
           known_lambda == -2.1) |>
  mutate(in_biased_ci = 
           known_lambda > lambda_under_lo &
           known_lambda < lambda_under_hi)
# biased ci plot ####
biased_lambda_ci |>
  group_by(M, known_lambda) |>
  arrange(lambda_under_lo) |>
  mutate(y = (1:n())) |>
  ggplot(aes(x = lambda_under, 
             xmin = lambda_under_lo,
             xmax = lambda_under_hi,
             y = y,
             color = in_biased_ci)) +
  geom_pointrange(fatten = 0,
                  linewidth = 0.25,
                  alpha = 0.75) +
  geom_vline(aes(xintercept = known_lambda)) +
  facet_grid(known_lambda~M,
             scales = "free") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(title = "Biased Data",
       y = "",
       x = "95% CI",
       color = "CI contains known value") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_x_continuous(n.breaks = 3)
ggsave("plots/biased_lambda_CIs.png",
       units = "in",
       height = 10,
       width = 10)
# biased percent CI table ####
biased_lambda_ci_percent <- biased_lambda_ci |>
  mutate(
    ci_width = lambda_under_hi - lambda_under_lo) |>
  group_by(known_lambda, M) |>
  summarize(n = n(), 
            prop_in = sum(in_biased_ci),
            mean_ci_width = mean(ci_width),
            sd_ci_width = sd(ci_width)) |>
  mutate(percent_ci = prop_in / n,
         data = "Biased") |>
  select(known_lambda,
         M,
         data,
         percent_ci,
         mean_ci_width, 
         sd_ci_width) 

lambda_summary_table <- xmin_lambda_ci_percent |>
  bind_rows(biased_lambda_ci_percent) |>
  bind_rows(fixed_lambda_ci_percent) |>
  mutate(bias_level = case_when(
    M == 0.125 ~ "minimal",
    M == 0.25 ~ "moderate",
    M == 0.5 ~ "strong",
    M == 1 ~ "extreme"
  )) |>
  left_join(lambdas_deviation) |>
  arrange(-known_lambda, M, data) |>
  select(known_lambda,
         bias_level,
         data:mean_ci_width,
         median_abs_delta) 

# Sample sizes ------------------------------------------------------------

n_fixed <- fixed_cut |>
  filter(group == "D") |>
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
  filter(group == "D") |>
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

three_ns <- bind_rows(n_fixed, n_xmin) 
three_ns_summary <- three_ns |>
  filter(known_lambda == -1.9 |
           known_lambda == -2 |
           known_lambda == -2.1) |>
  group_by(known_lambda, bias_level, name) |>
  summarise(mean_n_obs = mean(value),
            N_reps = n())

# save summary table
lambda_summary_table |>
  left_join(three_ns_summary |>
              rename(data = name)) |>
  write_csv("simulation_summaries/lambda_CI_percent.csv")

three_ns <- three_ns |>
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

three_ns  |>
  filter(known_lambda == -2) |>
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

ggsave("plots/N_just_one.png",
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
# multi panle manually stitched ####
# dropping for now
# beta_00_panel <- betas |>
#   filter(known_beta == 0,
#          term == "env_gradient") |>
#   select(known_beta, 
#          bias_level,
#          data_model,
#          estimate) |>
#   ggplot(aes(x = estimate,
#              y = data_model,
#              fill = data_model)) +
#   stat_halfeye(normalize = "groups") +
#   geom_vline(aes(xintercept = known_beta),
#              linetype = "dashed") +
#   scale_fill_manual(values = c(c("#FF914A",
#                                  "#019AFF",
#                                  "#FF1984"))) +
#   facet_grid(known_beta ~ bias_level,
#              labeller = label_value,
#              scales = "free") +
#   theme_bw() +
#   labs(x = "",
#        y = "")+
#   scale_x_continuous(n.breaks = 4)
# 
# beta_05_panel <- betas |>
#   filter(known_beta == -0.05,
#          term == "env_gradient") |>
#   select(known_beta, 
#          bias_level,
#          data_model,
#          estimate) |>
#   ggplot(aes(x = estimate,
#              y = data_model,
#              fill = data_model)) +
#   stat_halfeye(normalize = "groups") +
#   geom_vline(aes(xintercept = known_beta),
#              linetype = "dashed") +
#   scale_fill_manual(values = c(c("#FF914A",
#                                  "#019AFF",
#                                  "#FF1984"))) +
#   facet_grid(known_beta ~ bias_level,
#              labeller = label_value,
#              scales = "free") +
#   theme_bw() +
#   labs(x = "",
#        y = "")+
#   scale_x_continuous(n.breaks = 4)
# 
# beta_10_panel <- betas |>
#   filter(known_beta == -0.1,
#          term == "env_gradient") |>
#   select(known_beta, 
#          bias_level,
#          data_model,
#          estimate) |>
#   ggplot(aes(x = estimate,
#              y = data_model,
#              fill = data_model)) +
#   stat_halfeye(normalize = "groups") +
#   geom_vline(aes(xintercept = known_beta),
#              linetype = "dashed") +
#   scale_fill_manual(values = c(c("#FF914A",
#                                  "#019AFF",
#                                  "#FF1984"))) +
#   facet_grid(known_beta ~ bias_level,
#              labeller = label_value,
#              scales = "free") +
#   theme_bw() +
#   labs(x = "",
#        y = "")+
#   scale_x_continuous(n.breaks = 4)
# 
# beta_20_panel <- betas |>
#   filter(known_beta == -0.2,
#          term == "env_gradient") |>
#   select(known_beta, 
#          bias_level,
#          data_model,
#          estimate) |>
#   ggplot(aes(x = estimate,
#              y = data_model,
#              fill = data_model)) +
#   stat_halfeye(normalize = "groups") +
#   geom_vline(aes(xintercept = known_beta),
#              linetype = "dashed") +
#   scale_fill_manual(values = c(c("#FF914A",
#                                  "#019AFF",
#                                  "#FF1984"))) +
#   facet_grid(known_beta ~ bias_level,
#              labeller = label_value,
#              scales = "free") +
#   theme_bw() +
#   labs(x = "\u03b2 estimate",
#        y = "")+
#   scale_x_continuous(n.breaks = 4)
# 
# ggpubr::ggarrange(beta_00_panel,
#                   beta_05_panel,
#                   beta_10_panel,
#                   beta_20_panel,
#                   common.legend = TRUE,
#                   legend = "right",
#                   ncol = 1)
# one beta figure ####
betas |>
  filter(
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
             labeller = label_value) +
  theme_bw() +
  labs(x = "\u03b2 estimate",
       y = "")+
  scale_x_continuous(n.breaks = 4)
ggsave("plots/betas_all_three.png",
       units = "in",
       height = 10,
       width = 10)

# beta summary table ####
# beta deviation ####
# known value; Bias Level; Data; Percent Cis; Mean CI width; Median abs deviation; N replicates

beta_deviation <- betas |>
  ungroup() |>
  mutate(abs_delta = abs(known_beta - estimate)) |>
  group_by(known_beta, bias_level, data_model) |>
  summarise(median_abs_delta = median(abs_delta),
            n_reps = n())

# |>
#   write_csv("simulation_summaries/three_betas.csv")


# beta CIs ####
beta_ci <- betas |>
  mutate(ci_low = estimate - 1.96*std.error,
         ci_hi = estimate + 1.96*std.error,
         in_ci = known_beta > ci_low & known_beta < ci_hi) 

# 95% ci colored chart beta ####
# biased ####
beta_ci |>
  filter(data_model == "Biased") |>
  group_by(bias_level, known_beta) |>
  #slice_sample(prop = 0.5) |>
  arrange(ci_low) |>
  mutate(y = (1:n()*100)) |>
  ggplot(aes(x = estimate, 
             xmin = ci_low,
             xmax = ci_hi,
             y = y,
             color = in_ci)) +
  geom_pointrange(fatten = 0,
                  linewidth = 0.25,
                  alpha = 0.75) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~bias_level,
             scales = "free") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(title = "Biased Data",
       y = "",
       x = "95% CI") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("plots/biased_beta_CIs.png",
       units = "in",
       height = 10,
       width = 10)

# xmin ####
beta_ci |>
  filter(data_model == "xmin") |>
  group_by(bias_level, known_beta) |>
  #slice_sample(prop = 0.5) |>
  arrange(ci_low) |>
  mutate(y = (1:n()*100)) |>
  ggplot(aes(x = estimate, 
             xmin = ci_low,
             xmax = ci_hi,
             y = y,
             color = in_ci)) +
  geom_pointrange(fatten = 0,
                  linewidth = 0.25,
                  alpha = 0.75) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~bias_level,
             scales = "free") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(title = "Xmin Data",
       y = "",
       x = "95% CI") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("plots/xmin_beta_CIs.png",
       units = "in",
       height = 10,
       width = 10)
# fixed ####
beta_ci |>
  filter(data_model == "Fixed") |>
  group_by(bias_level, known_beta) |>
  #slice_sample(prop = 0.5) |>
  arrange(ci_low) |>
  mutate(y = (1:n()*100)) |>
  ggplot(aes(x = estimate, 
             xmin = ci_low,
             xmax = ci_hi,
             y = y,
             color = in_ci)) +
  geom_pointrange(fatten = 0,
                  linewidth = 0.25,
                  alpha = 0.75) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_grid(known_beta~bias_level,
             scales = "free") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(title = "Fixed Cutoff",
       y = "",
       x = "95% CI") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("plots/fixed_beta_CIs.png",
       units = "in",
       height = 10,
       width = 10)

# beta ci percent table ####
beta_ci_percent <- beta_ci |>
  ungroup() |>
  mutate(ci_width = ci_low - ci_hi) |>
  group_by(data_model, bias_level, known_beta, group) |>
  summarise(n = n(),
            prop_in = sum(in_ci),
            mean_ci_width = mean(ci_width)) |>
  select(n, prop_in, mean_ci_width) |>
  mutate(percent_ci = prop_in / n)

beta_ci_percent |>
  select(known_beta, bias_level, data_model, percent_ci, mean_ci_width) |>
  left_join(beta_deviation) |>
  arrange(known_beta, bias_level, data_model) |>
  write_csv("simulation_summaries/beta_CI_percent.csv")


# bar chart of percent beta CI ####
ggplot(beta_ci_percent,
       aes(x = known_beta, 
           y = percent_ci,
           fill = data_model)) +
  geom_col(position = "dodge") +
  facet_grid(.~bias_level) +
  scale_fill_manual(values = c(c("#FF914A",
                                 "#019AFF",
                                 "#FF1984"))) +
  theme_bw() 

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

