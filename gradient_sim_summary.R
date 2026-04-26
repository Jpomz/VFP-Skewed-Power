#gradient sim summary

library(tidyverse)
library(broom)
library(tidybayes)

results <- readRDS("simulation_results/gradient_sim_run.rds")

lm_models <- results |>
  select(group,
         lambda_under,
         lambda_trimmed,
         pr_scenario,
         rep,
         n,
         known_beta,
         env_gradient) |>
  pivot_longer(lambda_under:lambda_trimmed,
               names_to = "data_model", 
               values_to = "lambda_est") |>
  group_by(group, n, pr_scenario, rep, known_beta, data_model) |>
  nest() |>
  mutate(lm_model = map(data, ~lm(lambda_est ~ env_gradient, .))) |>
  mutate(coefs = map(lm_model, tidy)) |>
  unnest(coefs) 

# how many models incorrectly assumed a relationship when there was not one?
lm_models |>
  filter(term == "env_gradient", 
         known_beta == 0,
         p.value <0.05)

# how many models failed to detect a relationship when there was one?
lm_models |>
  filter(term == "env_gradient", 
         known_beta != 0,
         p.value >0.05)

lm_models |>
  filter(term == "env_gradient") |>
  ggplot(aes(x = estimate,
             fill = data_model)) +
  stat_halfeye(alpha = 0.5) +
  geom_vline(aes(xintercept = known_beta)) +
  facet_wrap(pr_scenario~.)
  

