# summarize undersampling simulation

library(tidyverse)
library(tidybayes)

sim <- as_tibble(readRDS("simulation_results/undersampling_sim_run.rds"))

sim |>
  group_by(rep) |>
  count() |>
  ungroup() |>
  group_by(n) |>
  count() |>
  mutate(runs = n*nn) |>
  ungroup() |>
  summarize(sum(runs)/18000)
# 36 original parameter sets and 500 reps = 18,000
# 
sim |>
  summarize(min_xmax_obs = min(xmax_obs),
            mean_xmax_obs = mean(xmax_obs),
            med_xmax_obs = median(xmax_obs),
            max_xmax_obs = max(xmax_obs))
# xmax_obs has a lot of variation

sim |>
  summarize(min_xmin_obs = min(xmin_obs),
            mean_xmin_obs = mean(xmin_obs),
            med_xmin_obs = median(xmin_obs),
            max_xmin_obs = max(xmin_obs))
# xmin_obs is very consistent at the set boundary of 0.001
  
ggplot(sim, 
       aes(x = original_n,
           y = under_n)) +
  geom_point(position = position_jitter(width = 500)) +
  facet_wrap(known_lambda~pr_scenario)

# distribution of lambda estimates (not including CIs)
sim |>
  as_tibble() |>
  filter(original_n == 5000,
         known_lambda == -2) |>
  select(b, h, known_lambda, lambda_under, lambda_trimmed, original_n) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.6) +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#019AFF"))) +
  facet_wrap(b ~ h,
             #scales = "free_x",
             labeller = label_both) +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "density"#,
       # caption = "Distribution of \u03bb estimates for under sampled (blue) and data which has been trimmed at the estimated x_min (pink). \nFacet titles show variables for undersampling function and the number of data points in the under sampled data \nincreases from left to right and top to bottom. Dashed line shows the known value of \u03bb"
       )


# distribution of lambda estimates (not including CIs)
sim |>
  as_tibble() |>
  filter(original_n == 5000,
         known_lambda == -2) |>
  select(b, h, known_lambda, lambda_under, lambda_trimmed, original_n) |>
  rename(`Biased` = lambda_under,
         `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.6) +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#019AFF"))) +
  facet_wrap(b ~ h,
             #scales = "free_x",
             labeller = label_both) +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "density"#,
       # caption = "Distribution of \u03bb estimates for under sampled (blue) and data which has been trimmed at the estimated x_min (pink). \nFacet titles show variables for undersampling function and the number of data points in the under sampled data \nincreases from left to right and top to bottom. Dashed line shows the known value of \u03bb"
       )
ggsave("plots/under_trimmed_estimates_lambda_2.png", 
       scale = 2)
  

# proportion of samples
sim |>
  select(h, b, original_n, under_n) |>
  mutate(prop = under_n / original_n) |>
  group_by(h, b) |>
  summarize(med_prop = median(prop),
            sd_prop = sd(prop))
# h == 0.0001 
  # b == 1.5 ~ 89% +- 2.7
  # b == 2 ~ 40% +- 12.2
# h == 0.001 
  # b == 1.5 ~53% +- 9.6
  # b == 2 ~15% +- 11.5
sim |>
  select(h, b, original_n, under_n) |>
  mutate(prop = under_n / original_n) |>
  group_by(h, b, original_n) |>
  summarize(med_prop = median(prop),
            mean_prop = mean(prop), 
            sd_prop = sd(prop))|>
  write_csv("simulation_summaries/proportion_samples_pr_scenarios.csv")

set.seed(112)
sim |>
  filter(original_n == 5000,
         known_lambda == -2) |>
  group_by(pr_scenario) |>
  sample_n(100) |>
  arrange(lambda_under_lo) |>
  group_by(known_lambda, lambda_under_lo) |>
  mutate(id = cur_group_id(),
         color = lambda_under_lo < known_lambda & lambda_under_hi > known_lambda) |>
  ggplot(aes(y = id,
             x = lambda_under,
             xmin = lambda_under_lo,
             xmax = lambda_under_hi,
             color = color)) +
  geom_pointrange(
    position = position_jitter(height = 0.05)) +
  #geom_vline(aes(xintercept = known_lambda)) +
  facet_wrap(b~h,
             scales = "free",
             labeller = label_both) +
  theme_bw() +
  geom_vline(aes(xintercept = -2),
             linetype = "dashed") +
  scale_color_viridis_d(option = "plasma") +
  guides(
    color = guide_legend(
      title= "\u03bb in 95% CI?")) +
  labs(x = "Biased \u03bb estimate",
       y = "") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave("plots/ci_under_lambda_2.png", width = 9, height = 4.5)

# CI plot color if lambda is in or out
set.seed(112)
sim |>
  filter(original_n == 5000,
         known_lambda == -2) |>
  group_by(pr_scenario) |>
  sample_n(100) |>
  arrange(lambda_trimmed_lo) |>
  group_by(known_lambda, h, lambda_trimmed_lo) |>
  mutate(id = cur_group_id(),
         color = lambda_trimmed_lo < known_lambda & lambda_trimmed_hi > known_lambda) |>
  ggplot(aes(y = id,
             x = lambda_trimmed,
             xmin = lambda_trimmed_lo,
             xmax = lambda_trimmed_hi,
             color = color)) +
  geom_pointrange(position = position_jitter(height = 0.05)) +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  facet_wrap(b~h,
             scales = "free",
             labeller = label_both)+
  theme_bw() +
  scale_color_viridis_d(option = "plasma", end = 0.75) +
  guides(
    color = guide_legend(
      title= "\u03bb in 95% CI?")) +
  labs(x = "Censored \u03bb estimate") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave("plots/ci_trimmed_lambda_2.png",
       width = 9, height = 4.5)

# proportion of trimmed CIs with known lambda
sim |>
  filter(original_n == 5000) |>
  arrange(lambda_trimmed_lo) |>
  group_by(known_lambda, h, lambda_trimmed_lo) |>
  mutate(id = cur_group_id(),
         color = lambda_trimmed_lo < known_lambda & lambda_trimmed_hi > known_lambda) |>
  ungroup() |>
  group_by(known_lambda, h, b, original_n) |>
  summarize(ci_true = sum(color) / n()) |>
  print(n = 36)

# plot of proportion of ci with lambda
sim |>
  #filter(original_n == 5000) |>
  arrange(lambda_trimmed_lo) |>
  group_by(known_lambda, h, lambda_trimmed_lo) |>
  mutate(id = cur_group_id(),
         color = lambda_trimmed_lo < known_lambda & lambda_trimmed_hi > known_lambda) |>
  ungroup() |>
  group_by(known_lambda,
           pr_scenario,
           original_n) |>
  summarize(ci_true = mean(color),
            ci_sd = sd(color)) |>
  ggplot(aes(x = original_n,
             y = ci_true, 
             ymin = ci_true - ci_sd, 
             ymax = ci_true + ci_sd,
             color = pr_scenario)) +
  geom_pointrange(
    position = position_dodge(width = 1000)
  )+
  facet_grid(pr_scenario~known_lambda) +
  scale_color_viridis_d(option = "plasma")


sim |>
  # filter(known_lambda == -2, 
  #        original_n == 5000) |>
  ggplot(aes(x = est_xmin,
             fill = as.factor(known_lambda))) +
  stat_halfeye(alpha = 0.5, 
               normalize = "panels") +
  facet_wrap(original_n~pr_scenario,
             #scales = "free",
             ncol = 4) +
  scale_x_log10() +
  theme_bw() +
  guides(
    fill = guide_legend(
      title= "Known \u03bb"))
ggsave("plots/dist_xmin.png",
       width = 9, height = 4.5)

sim |>
  # filter(known_lambda == -2, 
  #        original_n == 5000) |>
  ggplot(aes(x = est_xmin,
             fill = as.factor(known_lambda))) +
  stat_halfeye(alpha = 0.5, 
               normalize = "panels") +
  facet_wrap(pr_scenario~original_n,
             #scales = "free",
             ncol = 3) +
  scale_x_log10() +
  theme_bw() +
  guides(
    fill = guide_legend(
      title= "Known \u03bb"))
ggsave("plots/dist_xmin_increasing_bias.png",
       width = 9, height = 4.5)

sim |>
  filter(known_lambda == -2) |>
  group_by(pr_scenario) |>
  reframe(q05 = quantile(est_xmin, probs = 0.05),
          q95 = quantile(est_xmin, probs = 0.95)) |>
  arrange(pr_scenario)

sim |>
  group_by(pr_scenario, known_lambda, original_n) |>
  summarize(median_xmin = round(median(est_xmin), 3),
            mean_xmin = round(mean(est_xmin), 3),
            sd_xmin = round(sd(est_xmin),4)) |>
    write_csv("simulation_summaries/est_xmin.csv")  

# bias tables
sim |>
  filter(original_n == 5000) |>
  mutate(deviation_under = known_lambda - lambda_under,
         deviation_trimmed = known_lambda - lambda_trimmed) |>
  group_by(pr_scenario, known_lambda) |>
  summarize(med_dev_und = median(deviation_under), 
            med_dev_tri = median(deviation_trimmed))

under <- sim |>
  select(known_lambda:lambda_under_hi, 
         original_n,
         pr_scenario,
         rep)|>
  rename(est = lambda_under,
         minCI = lambda_under_lo,
         maxCI = lambda_under_hi) |>
  mutate(data_source = "bias")
trimmed <- sim |>
  select(known_lambda,
         lambda_trimmed:lambda_trimmed_hi, 
         original_n,
         pr_scenario,
         rep) |>
  rename(est = lambda_trimmed,
         minCI = lambda_trimmed_lo,
         maxCI = lambda_trimmed_hi)|>
  mutate(data_source = "censored")

long_dat <- bind_rows(
  under, 
  trimmed
)

lambda_bias <- long_dat|>
  mutate(conf_width = maxCI - minCI,
         diff = est - known_lambda,
         abs_bias = abs(diff)) |>
  group_by(known_lambda, 
           original_n, 
           pr_scenario, 
           data_source) |>
  add_count() |>
  #group_by(target_name, name, n) %>%
  summarize(median_ci_range = median(conf_width),
            median_abs_bias = median(abs_bias),
            sd_abs_bias = sd(abs_bias)) |>
  arrange(original_n,
          known_lambda,
          pr_scenario,
          data_source) 
write_csv(lambda_bias, 
          "simulation_summaries/bias_table.csv")

lambda_ci_prop <- long_dat %>%
  mutate(in_ci = known_lambda > minCI & known_lambda < maxCI) %>%
  na.omit() %>%
  group_by(known_lambda,
           original_n,
           pr_scenario,
           data_source) %>%
  summarize(count = n(),
            proportion = sum(in_ci, na.rm = TRUE) / count) 
write_csv(lambda_ci_prop, 
          "simulation_summaries/lambda_ci_table.csv")

lambda_ci_prop |>
  filter(proportion<0.1) |>
  group_by(data_source, pr_scenario, original_n) |>
  count()

lambda_ci_prop |>
  filter(proportion>0.1) |>
  group_by(data_source, pr_scenario, original_n) |>
  count()

lambda_ci_prop |>
  group_by(pr_scenario) |>
  summarise(mean(proportion))
lambda_ci_prop |>
  group_by(known_lambda) |>
  summarise(mean(proportion))
lambda_ci_prop |>
  group_by(original_n) |>
  summarise(mean(proportion))
lambda_ci_prop |>
  group_by(data_source) |>
  summarise(mean(proportion))

lambda_ci_prop |>
  group_by(data_source, 
           pr_scenario) |>
  summarise(mean(proportion))

lambda_ci_prop |>
  group_by(data_source, 
           pr_scenario,
           known_lambda) |>
  summarise(mean(proportion))


lambda_bias |>
  group_by(pr_scenario) |>
  summarise(mean(median_ci_range))
lambda_bias |>
  group_by(known_lambda) |>
  summarise(mean(median_ci_range))
lambda_bias |>
  group_by(original_n) |>
  summarise(mean(median_ci_range))
lambda_bias |>
  group_by(data_source) |>
  summarise(mean(median_ci_range))


lambda_bias |>
  group_by(pr_scenario) |>
  summarise(mean(median_abs_bias))
lambda_bias |>
  group_by(known_lambda) |>
  summarise(mean(median_abs_bias))
lambda_bias |>
  group_by(original_n) |>
  summarise(mean(median_abs_bias))
lambda_bias |>
  group_by(data_source) |>
  summarise(mean(median_abs_bias))

