# morin search probabilities summary

library(tidyverse)
library(tidybayes)

results <- readRDS("simulation_results/fixed_cut_morin_search_sim_run.rds") |> 
  as_tibble() |>
  mutate(cutoff_probabilities = 
           factor(cutoff_probabilities, 
                  levels = c(c(
                    "10%",
                    "50%",
                    "90%",
                    "92.5%",
                    "93%",
                    "94%",
                    "95%",
                    "97.5%",
                    "99.0%",
                    "99.5%",
                    "10xM"))))

cutoff_df <- readRDS(
  "simulation_results/morin_fixed_cutoff_search.rds") |>
  as_tibble() |>
  select(M, cutoff_probabilities, 
         cutoff_length = cutoff, 
         cutoff_mass) |>
  mutate(cutoff_probabilities = 
           factor(cutoff_probabilities, 
                  levels = c(c(
                    "10%",
                    "50%",
                    "90%",
                    "92.5%",
                    "93%",
                    "94%",
                    "95%",
                    "97.5%",
                    "99.0%",
                    "99.5%",
                    "10xM"))))

# lambda ~ retention probability, M {0.125, 0.25, 0.5, 1}
results |>
  left_join(cutoff_df) |>
  select(known_lambda, M, cutoff_probabilities, lambda_under, lambda_trimmed) |>
  rename(`Biased` = lambda_under,
                 `Censored` = lambda_trimmed) |>
  pivot_longer(`Biased`:`Censored`) |>
  ggplot(aes(x = value,
             y = name,
             fill = name)) +
  stat_halfeye(alpha = 1,
               normalize = "groups") +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  facet_grid(cutoff_probabilities ~ M,
             labeller = label_value,
             scales = "free") +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "")+
  scale_x_continuous(n.breaks = 3)
ggsave("plots/lambda_morin_search_probabilities.png",
       units = "in",
       height = 10,
       width = 8)

# N ~ retention probability
results |>
  pivot_longer(c(under_n, trimmed_n)) |>
  ggplot(aes(x = value, 
             y = name, 
             fill = name))+
  stat_halfeye(normalize = "groups") +
  facet_grid(cutoff_probabilities~M,
             scales = "free") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#FF914A"))) +
  labs(x = "N",
       y = "")
ggsave("plots/lambda_morin_N_ret-prob.png",
       units = "in",
       height = 10,
       width = 8)



# summary table ------------------------------------

lambda_ci <- results |>
  select(known_lambda:lambda_trimmed_hi,
         M,
         cutoff_probabilities,
         trimmed_n) |>
  mutate(
    under_abs_delta = abs(known_lambda - lambda_under),
    fixed_abs_delta = abs(known_lambda - lambda_trimmed),
    in_fixed_ci = 
      known_lambda > lambda_trimmed_lo &
      known_lambda < lambda_trimmed_hi,
    fixed_ci_width = lambda_trimmed_hi - lambda_trimmed_lo) 

lambda_ci |>
  group_by(M, cutoff_probabilities) |>
  count()

lambda_ci |>
  group_by(M, cutoff_probabilities) |>
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
  facet_grid(cutoff_probabilities~M,
             scales = "free") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(title = "Fixed Cutoff",
       y = "",
       x = "95% CI") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("plots/lambda_morin_cutoff_probs_CIs.png",
       units = "in",
       height = 10,
       width = 10)

lambda_ci_percent <- lambda_ci |>
  group_by(M, cutoff_probabilities) |>
  summarize(n = n(), 
            mean_obs = mean(trimmed_n),
            sd_obs = sd(trimmed_n),
            prop_in = sum(in_fixed_ci),
            mean_fixed_ci_width = mean(fixed_ci_width)) |>
  mutate(percent_ci = prop_in / n) |>
  select(M, cutoff_probabilities, n, 
         mean_obs, sd_obs, 
         percent_ci, mean_fixed_ci_width) |>
  left_join(cutoff_df) |>
  select(M:cutoff_probabilities, n, 
         cutoff_length:cutoff_mass,
         mean_obs:mean_fixed_ci_width)

write_csv(lambda_ci_percent, 
          "simulation_summaries/lambda_ci_percent_morin_cutoff_probs.csv")

plogis(morin_ln_p(L = 2.5, M = 0.25))
plogis(morin_ln_p(L = 5, M = 0.5))
plogis(morin_ln_p(L = 1.25, M = 0.125))
