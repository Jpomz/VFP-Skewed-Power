# morin search probabilities summary

library(tidyverse)
library(tidybayes)

results <- readRDS("simulation_results/fixed_cut_morin_search_sim_run.rds")

results |>
  group_by(known_lambda, M, cutoff_probabilities) |> 
  summarize(med_lambda_fixed = median(lambda_trimmed),
            sd_lambda_fixed = sd(lambda_trimmed),
            min_lambda = min(lambda_trimmed),
            max_lambda = max(lambda_trimmed)) |>
  mutate(range_delta = max_lambda - min_lambda) |>
  select(M, cutoff_probabilities, range_delta)



med_sd_summary <- results |>
  group_by(known_lambda, M, cutoff_probabilities) |> 
  summarize(med_lambda_under = median(lambda_under),
            med_lambda_fixed = median(lambda_trimmed),
            sd_lambda_fixed = sd(lambda_trimmed))
med_sd_summary
med_sd_summary |>
  mutate(delta = known_lambda - med_lambda_fixed) |>
  select(known_lambda, M, cutoff_probabilities, delta, sd_lambda_fixed) |>
  arrange(abs(delta))
  
write_csv(med_sd_summary,
  "simulation_summaries/search_prob_summary.csv")

results |>
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

results |>
  pivot_longer(c(under_n, trimmed_n)) |>
  ggplot(aes(x = value, 
             y = name, 
             fill = name))+
  stat_halfeye(normalize = "groups") +
  facet_grid(cutoff_probabilities~M,
             scales = "free")
