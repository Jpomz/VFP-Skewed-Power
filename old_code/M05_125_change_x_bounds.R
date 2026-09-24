# changing x_bounds on M 0.25 and 0.5

library(tidyverse)
library(tidybayes)

m05 <- readRDS("simulation_results/M0.5_x_bounds_sim_run.rds")
m125 <- readRDS("simulation_results/M0.125_x_bounds_sim_run.rds")


# helper functions --------------------------------------------------------

# functions to quickly calculate summaries
under_fixed_lambda_halfeye <- function(dat,
                                       m_size, 
                                       subtitle,
                                       my_colors = c("#FF1984",                                        "#FF914A")){
  dat |>
    select(known_lambda, lambda_under, lambda_trimmed) |>
    rename(`Biased` = lambda_under,
           `Fixed` = lambda_trimmed) |>
    pivot_longer(`Biased`:`Fixed`) |>
    ggplot(aes(x = value,
               y = name,
               fill = name)) +
    stat_halfeye(normalize = "groups") +
    geom_vline(aes(xintercept = known_lambda),
               linetype = "dashed") +
    scale_fill_manual(values = c(my_colors)) +
    facet_grid(known_lambda ~ .,
               labeller = label_value,
               scales = "free") +
    theme_bw() +
    labs(title = m_size,
         subtitle = subtitle,
         x = "\u03bb estimate",
         y = "") +
    scale_x_continuous(n.breaks = 4)
}
under_delta <- function(x, y){
  x |>
    mutate(
      abs_delta = abs(known_lambda - {{y}})) |>
    group_by(known_lambda)|>
    summarise(median_abs_delta = median(abs_delta),
              sd_abs_delta = sd(abs_delta))
}
ci_percent <- function(x, y_lo, y_hi){
  x |>
    mutate(in_ci = known_lambda > {{y_lo}} &
             known_lambda < {{y_hi}},
           ci_width = {{y_hi}} - {{y_lo}}) |>
    group_by(known_lambda) |>
    summarise(n = n(), 
              sum_ci = sum(in_ci),
              mean_ci_width = mean(ci_width)) |>
    mutate(prop_ci = sum_ci / n)
}
# ci_percent(m05,
#            lambda_trimmed_lo,
#            lambda_trimmed_hi)
delta_ci_table <- function(x, y1, y1_lo, y1_hi,
                           y2, y2_lo, y2_hi,
                           y1_name, y2_name,
                           out_name){
  
  delta1 <- under_delta(x, {{y1}})
  delta2 <- under_delta(x, {{y2}})
  
  ci_percent1 <- ci_percent(x,
                            {{y1_lo}},
                            {{y1_hi}}) 
  ci_percent2 <- ci_percent(x,
                            {{y2_lo}},
                            {{y2_hi}}) 
  
  
  table1 <- ci_percent1 |> 
    left_join(delta1) |>
    mutate(data = y1_name)
  
  table2 <- ci_percent2 |> 
    left_join(delta2) |>
    mutate(data = y2_name)
  
  bind_rows(table1, 
            table2) |>
    select(known_lambda, 
           data,
           prop_ci, 
           mean_ci_width,
           median_abs_delta,
           sd_abs_delta) |>
    write_csv(out_name)
  
}

# M125 --------------------------------------------------------------------

m125 <- as_tibble(m125)

under_fixed_lambda_halfeye(m125,
                           m_size = "M = 0.125",
                           subtitle = "Smaller body size bounds")
ggsave("plots/lambdas_M0125_x_bounds.png",
       units = "in",
       height = 10,
       width = 10)
delta_ci_table(x = m125, 
               y1 = lambda_under,
               y1_lo = lambda_under_lo,
               y1_hi = lambda_under_hi,
               y1_name = "Biased",
               y2 = lambda_trimmed,
               y2_lo = lambda_trimmed_lo,
               y2_hi = lambda_trimmed_hi,
               y2_name = "Fixed",
               out_name = "simulation_summaries/m125_xbounds.csv"
)

m125 |>
  mutate(
    in_fixed_ci = known_lambda > lambda_trimmed_lo &
           known_lambda < lambda_trimmed_hi) |>
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
  labs(title = "M=0.125",
       subtitle = "smaller x bounds",
       y = "",
       x = "95% CI") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("plots/m125_smaller_x_bounds_fixed_lambda_CIs.png",
       units = "in",
       height = 10,
       width = 10)
  
m125 |>
  mutate(
    in_under_ci = known_lambda > lambda_under_lo &
      known_lambda < lambda_under_hi) |>
  arrange(lambda_under_lo) |>
  mutate(y = (1:n())) |>
  ggplot(aes(x = lambda_under, 
             xmin = lambda_under_lo,
             xmax = lambda_under_hi,
             y = y,
             color = in_under_ci)) +
  geom_pointrange(fatten = 0,
                  linewidth = 0.25,
                  alpha = 0.75) +
  geom_vline(aes(xintercept = known_lambda)) +
  facet_grid(known_lambda~M,
             scales = "free") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(title = "M=0.125",
       subtitle = "smaller x bounds",
       y = "",
       x = "95% CI") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("plots/m125_smaller_x_bounds_under_lambda_CIs.png",
       units = "in",
       height = 10,
       width = 10)
# M500 --------------------------------------------------------------------

m05 <- as_tibble(m05)

under_fixed_lambda_halfeye(
  m05,
  m_size = "M = 0.5 mm",
  subtitle = "Larger body size bounds",
  my_colors = c("#FF1984",
                "springgreen"))

ggsave("plots/lambdas_M05_x_bounds.png",
       units = "in",
       height = 10,
       width = 10)
delta_ci_table(x = m05, 
               y1 = lambda_under,
               y1_lo = lambda_under_lo,
               y1_hi = lambda_under_hi,
               y1_name = "Biased",
               y2 = lambda_trimmed,
               y2_lo = lambda_trimmed_lo,
               y2_hi = lambda_trimmed_hi,
               y2_name = "Corrected",
               out_name = "simulation_summaries/m05_xbounds.csv"
)

m05 |>
  mutate(
    in_fixed_ci = known_lambda > lambda_trimmed_lo &
      known_lambda < lambda_trimmed_hi) |>
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
  labs(title = "M = 0.5; Corrected",
       subtitle = "larger x bounds",
       y = "",
       x = "95% CI") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("plots/m5_larger_x_bounds_corrected_lambda_CIs.png",
       units = "in",
       height = 10,
       width = 10)

m05 |>
  mutate(
    in_under_ci = known_lambda > lambda_under_lo &
      known_lambda < lambda_under_hi) |>
  arrange(lambda_under_lo) |>
  mutate(y = (1:n())) |>
  ggplot(aes(x = lambda_under, 
             xmin = lambda_under_lo,
             xmax = lambda_under_hi,
             y = y,
             color = in_under_ci)) +
  geom_pointrange(fatten = 0,
                  linewidth = 0.25,
                  alpha = 0.75) +
  geom_vline(aes(xintercept = known_lambda)) +
  facet_grid(known_lambda~M,
             scales = "free") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(title = "M = 0.5; Biased",
       subtitle = "smaller x bounds",
       y = "",
       x = "95% CI") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("plots/mo5_larger_x_bounds_under_lambda_CIs.png",
       units = "in",
       height = 10,
       width = 10)
