#figure out minimal to extreme bias

library(tidyverse)
library(sizeSpectra)

# function ####
sample_pr <- function(x, h, b){
  pr = 1 / (1 + (h / x**b))
  return(pr)
}

# SI figures for sampling probability with simulation values

set.seed(2112)
dat <- expand_grid(
  h = c(0.00001, 0.0001, 0.001, 0.01, 0.1),
  b = c(1.5, 1.75, 2),
  x = c(exp(seq(log(0.001), log(10), length.out = 1000)))) |>
mutate(pr = sample_pr(x = x, h = h, b = b))


dat <- dat |>
  group_by(h, b) |>
  mutate(scenario = cur_group_id())


# figure out new h and b values
# all combos on one plot:
dat |>
  ggplot(aes(x = x, 
             y = pr, 
             color = as.factor(scenario))) +
  geom_line(linewidth = 2) +
  geom_vline(aes(xintercept = 0.01), linetype = "dashed") +
  geom_hline(aes(yintercept = 0.9), linetype = "dashed") +
  scale_x_log10() +
  theme_bw() 
dat |>
  distinct(h, b, scenario)

dat |>
  filter(scenario == 1 |
           scenario == 4 |
           scenario == 7 |
           scenario == 10) |>
  distinct(h, b)

dat |>
  filter(scenario == 1 |
           scenario == 4 |
           scenario == 7 |
           scenario == 10) |>
  mutate(scenario = case_when(
    h == 0.00001 & b == 1.5 ~ "minimal",
    h == 0.0001 & b == 1.5 ~ "moderate",
    h == 0.001 & b == 1.5 ~ "strong",
    h == 0.01 & b == 1.5 ~ "extreme",
  ),
  scenario = factor(scenario, 
                    levels = c("minimal", "moderate", "strong", "extreme"))) |>
  ggplot(aes(x = x, 
             y = pr, 
             color = as.factor(scenario))) +
  geom_line(linewidth = 2) +
  geom_vline(aes(xintercept = 0.01), linetype = "dashed") +
  geom_hline(aes(yintercept = 0.9), linetype = "dashed") +
  scale_x_log10() +
  theme_bw() +
  scale_colour_manual(values = c("darkorchid1",
                                 "darkorchid2",
                                 "darkorchid3",
                                 "darkorchid4"),
                      name = "Bias") +
  labs(title = "Sampling probability as a function of body mass",
       x =expression(Log[10]~dry~mass),
       y = "Sampling probability") 

ggsave("plots/bias_scenarios_MS.png",
       units = "px",
       height = 1182,
       width = 2228)

## 4 scenarios on one plot
# dat |>
#   filter(scenario == 3 |
#            scenario == 6 |
#            scenario == 9 |
#            scenario == 12)|>
#   mutate(scenario = case_when(
#     h == 0.00001 & b == 2 ~ "minimal",
#     h == 0.0001 & b == 2 ~ "moderate",
#     h == 0.001 & b == 2 ~ "strong",
#     h == 0.01 & b == 2 ~ "extreme",
#   ),
#   scenario = factor(scenario, 
#                        levels = c("minimal", "moderate", "strong", "extreme"))) |>
#   ggplot(aes(x = x, 
#              y = pr, 
#              color = scenario)) +
#   geom_line(linewidth = 1.5) +
#   geom_vline(aes(xintercept = 0.01), linetype = "dashed") +
#   geom_hline(aes(yintercept = 0.9), linetype = "dashed") +
#   scale_x_log10() +
#   theme_bw() +
#   labs(title = "Sampling probability as a function of body mass",
#        x =expression(Log[10]~dry~mass),
#        y = "Sampling probability") +
#   scale_colour_manual(values = c("darkorchid1",
#                                  "darkorchid2",
#                                  "darkorchid3",
#                                  "darkorchid4")) +
#   NULL
# ggsave("plots/sample_pr_4_scenarios.png",
#        units = "px",
#        height = 1182,
#        width = 2228)

## 2 scenarios on one plot
dat |>
  filter(scenario == 3 |
           #scenario == 6 |
           scenario == 9 #|
           #scenario == 12
           )|>
  mutate(scenario = case_when(
    h == 0.00001 & b == 2 ~ "minimal",
    #h == 0.0001 & b == 2 ~ "moderate",
    h == 0.001 & b == 2 ~ "strong",
    #h == 0.01 & b == 2 ~ "extreme",
  ),
  scenario = factor(scenario, 
                    levels = c(
                      "minimal", 
                      #"moderate",
                      "strong"#,
                      #"extreme"
                      ))) |>
  ggplot(aes(x = x, 
             y = pr, 
             color = scenario)) +
  geom_line(linewidth = 1.5) +
  geom_vline(aes(xintercept = 0.01), linetype = "dashed") +
  geom_hline(aes(yintercept = 0.9), linetype = "dashed") +
  scale_x_log10() +
  theme_bw() +
  labs(title = "Sampling probability as a function of body mass",
       x =expression(Log[10]~dry~mass),
       y = "Sampling probability") +
  scale_colour_manual(values = c("darkorchid1",
                                 #"darkorchid2",
                                 "darkorchid3"#,
                                 #"darkorchid4"
                                 )) +
  NULL
ggsave("plots/sample_pr_2_scenarios.png",
       units = "px",
       height = 1182,
       width = 2228)
