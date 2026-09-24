# compare lengths from morin to weights in here

lw <- tibble(a = c(0.0082, 0.0077, 0.0082, 0.0028),
             b = c(2.168, 2.588, 2.813, 2.179))
lengths = c(1, 1.2, 1.3, 1.6)

lw_dat <- expand_grid(lw, lengths)
lw_dat |>
  mutate(dw = a*lengths**b) |>
  arrange(lengths)


expand_grid(lw, lengths = 2.5) |>
  mutate(dw = a*lengths**b) |>
  arrange(dw)

sample_pr <- function(x, h, b){
  pr = 1 / (1 + (h / x**b))
  return(pr)
}

dw_v <- c(0.001, 0.01, 0.1, 1)
# minimal - 0.01 ~ 99..%
sample_pr(x = 0.01,
          h = 0.00001,
          b = 1.5)
# moderate - 0.047 ~ 99.02%
sample_pr(x = 0.047,
          h = 0.0001,
          b = 1.5)
# strong - 0.22 ~ 99.04%
sample_pr(x = 0.22,
          h = 0.001,
          b = 1.5)

# extreme - 1.00 ~99.00%
sample_pr(x = 1.00,
          h = 0.01,
          b = 1.5)

# sampling probabilities for 4 cutoffs at different scenarios
expand_grid(x = c(0.01, 0.047, 0.22, 1.0),
            h = c(0.00001, 0.0001,0.001,0.01),
            b = 1.5) |>
  mutate(p = sample_pr(x = x, h = h, b = b))


# Calculate body size for given Pr ----------------------------------------
body_pr <- function(h, p, b){
  x = (h / ((1 / p) - 1))^(1/b)
  return(x)
}

body_pr(h = 0.00001,
        b = 1.5,
        p = c(0.9, 0.95, 0.99005))

sample_pr(x = c(0.002008299, 
                0.003304982,
                0.009933222),
          h = 0.00001, 
          b = 1.5)

body_pr(h = c(0.00001, 0.0001, 0.001, 0.01),
        b = 1.5,
        p = c(0.99))
sample_pr(x = c(.993322173, 1), 
          h = 0.01, 
          b = 1.5)
