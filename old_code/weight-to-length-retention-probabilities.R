# weight to length

library(tidyverse)
library(sizeSpectra)
# library(poweRlaw)

# morin's equation for 250um mesh
morin_ln_p  <- function(
    a = -2.84,
    b = 5.8, 
    c = -3.18, 
    L,
    M = 0.25){
  RL = L / M
  a + b * log10(RL) + c * log10(RL)*log10(M)
}

massToLength(masses = 1, 
            LWa = 0.005,
            LWb = 2.778)
#L = (M/a)^(1/b)
#lwa = 0.005
#lwb = 2.778

xm <- exp(seq(log(1e-5),
              log(1e2),
              length.out = 100))

xl <- massToLength(masses = xm,
                   LWa = 0.005,
                   LWb = 2.778)

sample_p <- plogis(morin_ln_p(L = xl))

tibble(x = c(xm, xl),
       variable = rep(c("mass", "length"), 
                      each = length(xm)),
       p = c(sample_p, sample_p)) |>
  ggplot() +
  geom_line(aes(x = x,
                y = p,
                color = variable,
                group = variable)) +
  scale_x_log10(guide = "axis_logticks")

# check probabilities for different meshes and one Length

tibble(xl = 1#exp(seq(log(1e-3),
       # log(1e2),
       # length.out = 100))
) |>
  mutate(p_125 = plogis(morin_ln_p(L = xl,
                                   M = .125)),
         p_250 = plogis(morin_ln_p(L = xl,
                                   M = .250)),
         p_500 = plogis(morin_ln_p(L = xl,
                                   M = .50)))
# This roughly matches what is in Morin

# do it in a tibble
tibble(xm = exp(seq(log(1e-3),
                    log(1e1),
                    length.out = 1000))
       ) |>
  mutate(xl = massToLength(xm, 
                           0.005, 
                           2.778),
         p_125 = plogis(morin_ln_p(L = xl,
                                   M = .125)),
         p_250 = plogis(morin_ln_p(L = xl,
                                   M = .250)),
         p_500 = plogis(morin_ln_p(L = xl,
                                   M = .50)),
         p_1000 = plogis(morin_ln_p(L = xl,
                                   M = 1))) |>
  pivot_longer(xm:xl) |>
  pivot_longer(p_125:p_1000, 
               names_to = "mesh", 
               values_to = "probablity") |>
  ggplot(aes(x = value, 
             y = probablity, 
             color = mesh, 
             group = mesh)) +
  geom_line() +
  scale_x_log10(guide = "axis_logticks") +
  theme_bw() +
  facet_grid(name~.)
  
       



# Length to Mass NEON Data ------------------------------------------------

# this section uses lengths from NEON macroinvertebrates (1-65 mm) and calculates sample probability based on length, and what the generic mass would look like
p_length_mass <- tibble(xl = seq(.1, 100, length.out = 1000)) |>
  mutate(xm = lengthToMass(xl, 
                           0.005, 
                           2.778),
         p_125 = plogis(morin_ln_p(L = xl,
                                   M = .125)),
         p_250 = plogis(morin_ln_p(L = xl,
                                   M = .250)),
         p_500 = plogis(morin_ln_p(L = xl,
                                   M = .50)),
         p_1000 = plogis(morin_ln_p(L = xl,
                                    M = 1))) |>
  pivot_longer(xm:xl) |>
  pivot_longer(p_125:p_1000, 
               names_to = "mesh", 
               values_to = "probablity") 

# both sets
p_length_mass |>
  ggplot(aes(x = value, 
             y = probablity, 
             color = mesh, 
             group = mesh)) +
  geom_line() +
  scale_x_log10(guide = "axis_logticks") +
  theme_bw() +
  labs(x = "Length or Mass") +
  facet_grid(name~.)

# just length
p_length_mass |>
  filter(name == "xl") |>
  ggplot(aes(x = value, 
             y = probablity, 
             color = mesh, 
             group = mesh)) +
  geom_line() +
  scale_x_log10(guide = "axis_logticks") +
  theme_bw() +
  geom_vline(aes(xintercept = 1),
             color = "black")+
  geom_vline(aes(xintercept = 65),
             color = "black")+
  labs(x = "Length") 

# just mass
# just length
p_length_mass |>
  filter(name == "xm") |>
  ggplot(aes(x = value, 
             y = probablity, 
             color = mesh, 
             group = mesh)) +
  geom_line() +
  scale_x_log10(guide = "axis_logticks") +
  theme_bw() +
  labs(x = "Mass") 

lengthToMass(1, 0.005, 2.778)
