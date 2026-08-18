# Morin's equations
# morin coefficients
library(tidyverse)
library(sizeSpectra)

ln_p  <- function(a = -2.84,
                  b = 5.8, 
                  c = -3.18, 
                  L,
                  M = 0.25){
  RL = L / M
  a + b * log10(RL) + c * log10(RL)*log10(M)
}

what_length = function(p,
                       a = -2.84,
                       b = 5.8, 
                       c = -3.18, 
                       M = 0.25 ){
  logit_p <- log(p/(1-p))
  L <- M * 10^((logit_p - a) / (b + c*log10(M)))
  return(L)
}
plogis(ln_p(L = .5))

what_length(p = 0.3733928)
what_length(p = 0.99, M = c(0.125, 0.25, 0.5, 1))

plogis(ln_p(L = 1.712277))



a = -2.84
b = 5.8
c = -3.18 
M = 0.25

L_test <- .5   # arbitrary test value

# Forward: L -> logit -> p
logit_val <- a*b*log10(L_test/M) + c*log10(L_test/M)*log10(M)
p <- plogis(logit_val)

# Backward: p -> L
logit_p <- log(p/(1-p))
L_recovered <- M * 10^(logit_p / (a*b + c*log10(M)))

L_test - L_recovered   # should be ~0 if algebra and code are consistent



















plogis(ln_p(L = 1))

log(0.9 / (1 - 0.9))

l_dat <- tibble(
  L = seq(0.1, 3, length.out = 100))

l_dat |>
  mutate(ln_p = ln_p(L = L),
         p = plogis(ln_p)) |>
  ggplot(aes(x = L, 
             y = p)) +
  geom_line() +
  scale_x_log10(guide = "axis_logticks")


l_dat |>
  mutate(ln_p = ln_p(L = L),
         p = plogis(ln_p)) |>
  filter(p >=0.99)
plogis(ln_p(L = c(0.96, 4.47, 6.57)))

# ~99.00% probabilities
plogis(ln_p(L = 0.901, M = 0.125))
plogis(ln_p(L = 2.30, M = 0.25))
plogis(ln_p(L = 6.31, M = 0.5))
plogis(ln_p(L = 19.15, M = 1))

sizeSpectra::lengthToMass(c(0.901,
                          2.30,
                          6.31,
                          19.15),
                          LWa = 0.0064,
                          LWb = 2.788)

# 95% probabilities
plogis(ln_p(L = 0.901, M = 0.125))
plogis(ln_p(L = 2.30, M = 0.25))
plogis(ln_p(L = 6.31, M = 0.5))
plogis(ln_p(L = 19.15, M = 1))


# why is there bimodality in the estimated xmins?
# known lambda = -1.5
# bias_level = extreme
sizeSpectra::lengthToMass(19.15,
                          LWa = 0.0064,
                          LWb = 2.788)
sizeSpectra::massToLength(0.6, 
                          LWa = 0.0064,
                          LWb = 2.788)

sizeSpectra::massToLength(0.001, 
                          LWa = 0.0064,
                          LWb = 2.788)
sizeSpectra::massToLength(100, 
                          LWa = 0.0064,
                          LWb = 2.788)
# bias levels -------------------------------------------------------------

morin_ln_p  <- function(
    a = -2.84,
    b = 5.8, 
    c = -3.18, 
    L,
    M = 0.25){
  RL = L / M
  a + b * log10(RL) + c * log10(RL)*log10(M)
}

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

p_length_mass <- p_length |>
  mutate(xm = lengthToMass(xl,
                           LWa = 0.005, 
                           LWb = 2.78))
p_length_mass
p_length_mass |>
  ggplot(aes(x = xm, 
             y = probability , 
             color = micron, 
             group = micron)) +
  geom_line(linewidth = 2) +
  scale_x_log10(guide = "axis_logticks") +
  theme_classic() +
  labs(x = "Mass") +
  scale_color_manual(values = c("darkorchid1",
                                "darkorchid2",
                                "darkorchid3",
                                "darkorchid4")) +
  geom_vline(aes(xintercept = 0.01),
             linetype = "dashed",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = 0.99),
             linetype = "dashed",
             linewidth = 1.25) +
  coord_cartesian(xlim = c(0.001, 4.5))
ggsave("plots/morin_mass_bias.png",
       units = "in",
       height = 6,
       width = 10)

p_length_mass|>
  rename(mass = xm, length = xl) |>
  pivot_longer(c(mass,length)) |>
  ggplot(aes(x = value, 
             y = probability , 
             color = micron, 
             group = micron)) +
  geom_line(linewidth = 2) +
  scale_x_log10(guide = "axis_logticks") +
  theme_classic() +
  labs(x = "Length or Mass") +
  facet_grid(name~.) +
  scale_color_manual(values = c("darkorchid1",
                                "darkorchid2",
                                "darkorchid3",
                                "darkorchid4")) +
  #coord_cartesian(xlim = c(0.001, 13))
  NULL
ggsave("plots/morin_length_mass_bias.png",
       units = "in",
       height = 6,
       width = 10)

# 95% for each mesh size
p_length_mass |>
  filter(probability >=0.940,
         probability<=0.96) |>
  group_by(mesh) |>
  summarize(mean(xl))

sizeSpectra::lengthToMass(c(0.584, 1.41, 3.63, 10.1), 0.0064, 2.778)

# sampling probability of smallest body size for each mesh
p_length_mass |>
  group_by(mesh) |>
  filter(
    xm >= 0.001)

# example of bias data ----------------------------------------------------
source("custom_functions.R")
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


# search probability ------------------------------------------------------
# when M = mesh size = 500 microns = 0.5 mm
strong_probabilities <- tibble(
  xl = seq(.1, 20, length.out = 1000)) |>
  mutate(p_500 = plogis(morin_ln_p(L = xl,
                                   M = .50)))
strong_probabilities |>
  filter(p_500 >= 0.09,
         p_500 <= 0.11)
# 10% = 0.625
strong_probabilities |>
  filter(p_500 >= 0.49,
         p_500 <= 0.51)
# 50% = 1.32
strong_probabilities |>
  filter(p_500 >= 0.895,
         p_500 <= 0.901)
# 90% = 2.79
strong_probabilities |>
  filter(p_500 >= 0.92,
         p_500 <= 0.93)
# 92.5% = 3.1
strong_probabilities |>
  filter(p_500 >= 0.929,
         p_500 <= 0.93)
#93% = 3.17
strong_probabilities |>
  filter(p_500 >= 0.939,
         p_500 <= 0.94)
#94% = 3.35
strong_probabilities |>
  filter(p_500 >= 0.949,
         p_500 <= 0.95)
#95% = 3.59
strong_probabilities |>
  filter(p_500 >= 0.974,
         p_500 <= 0.975)
# 97.5% = 4.56
strong_probabilities |>
  filter(p_500 >= 0.9900,
         p_500 <0.9901)
#99.0% = 6.31
strong_probabilities |>
  filter(p_500 >= 0.9989,
         p_500 <0.9990)
#99.9% = 13.5

sizeSpectra::lengthToMass(c(0.625,# 10% 
                            1.32,# 50% 
                            2.79,# 90% 
                            3.1, #92.5%
                            3.17, 
                            3.35,
                            3.59,# 95% 
                            4.56, #97.5
                            6.31,# 99.0% 
                            13.5), # 99.9%
                          LWa = 0.0064,
                          LWb = 2.788)
# body masses of:
# 0.001726209 0.013878422 0.111821289 0.225830864 1.088082505 9.068814540

ggplot(strong_probabilities,
       aes(x = xl, 
           y = p_500)) +
  geom_line() +
  scale_x_log10(guide = "axis_logticks")
max(strong_probabilities$p_500)
plogis(morin_ln_p(L = 100,
                  M = c(0.125, 0.25, .50, 1)))

# moderate probabilities --------------------------------------------------

mod_probabilities <- tibble(
  xl = seq(.1, 20, length.out = 1000)) |>
  mutate(p_250 = plogis(morin_ln_p(L = xl,
                                   M = .250)))
mod_probabilities |>
  filter(p_250 >= 0.09,
         p_250 <= 0.11)
# 10% = 0.299
mod_probabilities |>
  filter(p_250 >= 0.49,
         p_250 <= 0.51)
# 50% = 0.578
mod_probabilities |>
  filter(p_250 >= 0.895,
         p_250 <= 0.901)
# 90% = 0.898
mod_probabilities |>
  filter(p_250 >= 0.92,
         p_250 <= 0.93)
# 92.5% = 1.24
mod_probabilities |>
  filter(p_250 >= 0.929,
         p_250 <= 0.935)
# 93% = 1.28
mod_probabilities |>
  filter(p_250 >= 0.939,
         p_250 <= 0.945)
# 94% = 1.34
mod_probabilities |>
  filter(p_250 >= 0.94,
         p_250 <= 0.955)
#95% = 1.41
mod_probabilities |>
  filter(p_250 >= 0.974,
         p_250 <= 0.975)
# 97.5% = 1.73
mod_probabilities |>
  filter(p_250 >= 0.9900,
         p_250 <0.991)
#99.0% = 2.33
mod_probabilities |>
  filter(p_250 >= 0.99899,
         p_250 <0.999)
#99.9% = 4.58

sizeSpectra::lengthToMass(c(0.299,# 10% 
                            0.578,# 50% 
                            0.898,# 90% 
                            1.24, #92.5%
                            1.28, 
                            1.34,
                            1.41,# 95% 
                            1.73, #97.5
                            2.33,# 99.0% 
                            4.58), # 99.9%
                          LWa = 0.0064,
                          LWb = 2.788)


# 93% cutof probabilities -------------------------------------------------
plogis(morin_ln_p(L = 0.530, M = 0.125))
plogis(morin_ln_p(L = 1.265, M = 0.25))
plogis(morin_ln_p(L = 3.19, M = 0.5))
plogis(morin_ln_p(L = 8.625, M = 1))

sizeSpectra::lengthToMass(c(0.530, 
                            1.265, 
                            3.19, 
                            8.625), # 99.9%
                          LWa = 0.0064,
                          LWb = 2.788)
# mass 0.001090087, 0.012325608, 0.162460729, 2.600616749

# xmin bound > xmin empirical ####
# keep data but set xmin inside calcLike()?
x <- rPLB(5000, xmin = 0.001)

calcLike(negLL.fn = negLL.PLB,
         x = x,
         xmin = 0.1, #min(x), 
         xmax = max(x), 
         n = length(x), 
         sumlogx = sum(log(x)), 
         p = -1.5,
         suppress.warnings = TRUE,
         vecDiff = 2)
# doesn't work
# estimate of ~ -309