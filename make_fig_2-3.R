# make figures 2 and 3 for conceptual presentation

# load custom functions
source("custom_functions.R")

library(tidyverse)

# make figure 2 ####
# plot retention probabilities for different lengths
# assumes a length-weight relationship of
# M = aL^b
# where L is length, a = 0.0064, b = 2.788 
# Using the values for "all insects" from table 2 in Benke et al. 1999

# lengths and sampling probabilities
p_length <- tibble(
  xl = seq(.1, 13, length.out = 1000)) |>
  mutate(p_125 = plogis(morin_ln_p(L = xl,
                                   M = .125)),
         p_250 = plogis(morin_ln_p(L = xl,
                                   M = .250)),
         p_500 = plogis(morin_ln_p(L = xl,
                                   M = .50)),
         p_1000 = plogis(morin_ln_p(L = xl,
                                    M = 1))) |>
  pivot_longer(p_125:p_1000, 
               names_to = "mesh", 
               values_to = "probability") 
p_length
mesh_factor <- tibble(mesh = c("p_125", 
                               "p_250",
                               "p_500", 
                               "p_1000"),
                      micron = factor(
                        c("125um",
                          "250um",
                          "500um",
                          "1000um"),
                        levels = c("125um",
                                   "250um",
                                   "500um",
                                   "1000um")))
p_length <- p_length |>
  left_join(mesh_factor)

p_length |>
  ggplot(aes(x = xl, 
             y = probability, 
             color = micron, 
             group = micron)) +
  geom_line(linewidth = 2) +
  scale_x_log10(guide = "axis_logticks") +
  theme_classic() +
  labs(x = "Length") +
  scale_color_manual(values = c("darkorchid1",
                                "darkorchid2",
                                "darkorchid3",
                                "darkorchid4")) +
  geom_vline(aes(xintercept = 1.0),
             linetype = "dashed",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = 0.99),
             linetype = "dashed",
             linewidth = 1.25)
ggsave("plots/morin_length_bias.png",
       units = "in",
       height = 6,
       width = 10)


# Make figure 3 ####
set.seed(1152)
n = 10000
p1 <- plot_morin_bias(n = n,
                      lambda = -2,
                      M = 0.125,
                      binwidth = 0.1) +
  labs(title = "Minimal Bias; 125 micron") +
  coord_cartesian(xlim = c(0.001, 100))
p2 <- plot_morin_bias(n = n,
                      lambda = -2,
                      M = 0.25,
                      binwidth = 0.1) +
  labs(title = "Moderate Bias; 250 micron") +
  coord_cartesian(xlim= c(0.001, 100))
p3 <- plot_morin_bias(n = n,
                      lambda = -2,
                      M = 0.5,
                      binwidth = 0.1) +
  labs(title = "Strong Bias; 500 micron") +
  coord_cartesian(xlim= c(0.001, 100))
p4 <- plot_morin_bias(n = n,
                      lambda = -2,
                      M = 1,
                      binwidth = 0.1) +
  labs(title = "Extreme Bias; 1000 micron") +
  coord_cartesian(xlim= c(0.001, 100))

ggpubr::ggarrange(p1, p2, p3, p4, 
                  common.legend = TRUE,
                  legend = "right")
ggsave("plots/morin_bias_levels.png",
       units = "in",
       height = 6,
       width = 10)