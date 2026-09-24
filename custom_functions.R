# custom functions for simulating effects of under sampling small body sizes


# custom functions ####
# this function calculates the retention probability based on individual length (L, in mm), mesh size (M, in mm) and coefficients presented in Morin et al. 2004. 
morin_ln_p  <- function(
    L,
    M = 0.25, # 250 micron mesh is the default
    # a, b, and C are constants from Morin et al. 2004
    a = -2.84,
    b = 5.8, 
    c = -3.18 
    ){
  RL = L / M
  a + b * log10(RL) + c * log10(RL)*log10(M)
}

# this is a helper function which is used to make figure 3
# which shows retention probabilities based on indiviudal length and Mesh size
plot_morin_bias <- function(
    n = 5000, 
    lambda = -2,
    xmin = 0.001, 
    xmax = 100,
    M = .25,
    vecDiff = 2,
    LWa = 0.0064,
    LWb = 2.788,
    binwidth = 0.1,
    alpha = 0.75){
  x <- sizeSpectra::rPLB(n = n, b = lambda, xmin = xmin, xmax = xmax)
  
  # put masses into a df and convert to length
  x_df <- tibble(x = x)
  # convert mass to length
  # Using the values for "all insects" from table 2 in Benke et al. 1999
  # M = aL^b
  # where L is length, a = 0.0064, b = 2.788
  # solved for L = (M/LWa)^(1/b)
  x_df <- x_df |>
    mutate(xl = sizeSpectra::massToLength(x, LWa = LWa, LWb = LWb))
  
  # make a df with the sample probability and not/sampled columns
  x_df <- x_df |>
    mutate(probability = plogis(morin_ln_p(L = xl,
                                           M = M)), 
           sampled = rbinom(n(), 1, prob = probability),
           fill = case_when(sampled == 1 ~ "Biased", 
                            .default = "not sampled"))
  
  # under sampled ####
  # filter out the "sampled" data
  # this represents empirical data which has fewer little things than expected
  x_under <- x_df |>
    filter(fill == "Biased")
  
  plot_x_df <- tibble(x = x, 
                      fill = "Original")
  
  plot_df <- bind_rows(
    x_under |>
      select(x, fill),
    plot_x_df)
  
  p <- ggplot(plot_df,
              aes(x = x, 
                  fill = fill)) +
    geom_histogram(binwidth = binwidth, 
                   position = "dodge",
                   alpha = alpha) +
    scale_fill_manual(values = c("#FF1984", "black"))+
    scale_x_log10(guide = "axis_logticks") +
    # scale_y_log10() +
    theme_classic() 
  print(p)
}

# this is one of the main functions in the simulation
# this estimates lambda from biased data as well as estimates the x_min value, censors the data, and estimates lambda from x_min censored data
morin_bias_lambda <- function(
    n = 5000, 
    lambda = -2,
    xmin = 0.01, 
    xmax = 100,
    M = .25,
    vecDiff = 2,
    LWa = 0.0064,
    LWb = 2.788){
  # sample masses from a bounded power law
  x <- rPLB(n = n, b = lambda, xmin = xmin, xmax = xmax)
  xmin_obs = min(x)
  xmax_obs = max(x)
  
  # put masses into a df and convert to length
  x_df <- tibble(x = x)
  # convert mass to length
  # Using the values for "all insects" from table 2 in Benke et al. 1999
  # M = aL^b
  # where L is length, a = 0.0064, b = 2.788
  # solved for L = (M/LWa)^(1/b)
  x_df <- x_df |>
    mutate(xl = massToLength(x, LWa = LWa, LWb = LWb))
  
  # make a df with the sample probability and not/sampled columns
  x_df <- x_df |>
    mutate(probability = plogis(morin_ln_p(L = xl,
                                           M = M)), 
    sampled = rbinom(n(), 1, prob = probability),
    fill = case_when(sampled == 1 ~ "sampled", 
                     .default = "not sampled"))
  
  # under sampled ####
  # filter out the "sampled" data
  # this represents empirical data which has fewer little things than expected
  x_under <- x_df |>
    filter(fill == "sampled")
  # how many body sizes were sampled?
  under_n <- nrow(x_under)
  
  # estimate x_min ####
  # use poweRlaw to estimate where x_min is, i.e., x < x_min is under sampled
  x_power <- conpl$new(x_under$x)
  x_xmin <- estimate_xmin(x_power)$xmin
  
  # add estimated x_min to data frames
  x_df$est_xmin <- x_xmin
  x_under$est_xmin <- x_xmin
  
  # lambdas ####
  # lambda under ####
  # estimate lambda from undersampled data
  x_under_vector <- x_under$x
  lambda_under <- calcLike(negLL.fn = negLL.PLB,
                           x = x_under_vector,
                           xmin = min(x_under_vector), 
                           xmax = max(x_under_vector), 
                           n = length(x_under_vector), 
                           sumlogx = sum(log(x_under_vector)), 
                           p = -1.5,
                           suppress.warnings = TRUE,
                           vecDiff = vecDiff)
  
  # lambda trimmed ####
  # estimate lambda from trimmed data
  x_trimmed_vector <- x_under |>
    filter(x >= est_xmin) |>
    pull(x)
  trimmed_n <- length(x_trimmed_vector)
  lambda_trimmed <- calcLike(negLL.fn = negLL.PLB,
                             x = x_trimmed_vector,
                             xmin = min(x_trimmed_vector), 
                             xmax = max(x_trimmed_vector), 
                             n = length(x_trimmed_vector), 
                             sumlogx = sum(log(x_trimmed_vector)), 
                             p = -1.5,
                             suppress.warnings = TRUE,
                             vecDiff = 2)
  # end/return ####
  
  
  # return data frame ####
  out_df <- data.frame(
    known_lambda = lambda, 
    lambda_under = lambda_under$MLE,
    lambda_under_lo = lambda_under$conf[1],
    lambda_under_hi = lambda_under$conf[2],
    lambda_trimmed = lambda_trimmed$MLE,
    lambda_trimmed_lo = lambda_trimmed$conf[1],
    lambda_trimmed_hi = lambda_trimmed$conf[2],
    original_n = n,
    under_n = under_n,
    trimmed_n = trimmed_n,
    est_xmin = x_xmin, 
    M = M,
    xmin_obs = xmin_obs,
    xmin_under = min(x_under_vector),
    xmin_trimmed = min(x_trimmed_vector),
    xmax_obs = xmax_obs)
  return(out_df)
}

# this is the other main simulation function
# it biases data and censors it based on a fixed cutoff value
# it then estimates lambda from both the biased and censored data
morin_fixed_cut_lambda <- function(
    n = 5000, 
    lambda = -2,
    xmin = 0.01, 
    xmax = 100,
    M = .25,
    vecDiff = 2,
    LWa = 0.0064,
    LWb = 2.788,
    cutoff = 0.1 # cutoff body length
    ){
  # Mass ####
  # sample masses from a bounded power law
  x <- rPLB(n = n, b = lambda, xmin = xmin, xmax = xmax)
  xmin_obs = min(x)
  xmax_obs = max(x)
  
  # put masses into a df and convert to length
  x_df <- tibble(x = x)
  # lengths ####
  # convert mass to length
  # Using the values for "all insects" from table 2 in Benke et al. 1999
  # M = aL^b
  # where L is length, a = 0.0064, b = 2.788
  # solved for L = (M/LWa)^(1/b)
  x_df <- x_df |>
    mutate(xl = massToLength(x, LWa = LWa, LWb = LWb))
  
  # bias ####
  # make a df with the sample probability and not/sampled columns
  x_df <- x_df |>
    mutate(probability = plogis(morin_ln_p(L = xl,
                                           M = M)), 
           sampled = rbinom(n(), 1, prob = probability),
           fill = case_when(sampled == 1 ~ "sampled", 
                            .default = "not sampled"))
  
  # under sampled ####
  # filter out the "sampled" data
  # this represents empirical data which has fewer little things than expected
  x_under <- x_df |>
    filter(fill == "sampled")
  # how many body sizes were sampled?
  under_n <- nrow(x_under)
  
  # cutoff body length ####
  
  # add cutoff length to data frames
  x_df$cutoff <- cutoff
  x_under$cutoff <- cutoff
  
  # lambdas ####
  # lambda under ####
  # estimate lambda from undersampled data
  x_under_vector <- x_under$x
  lambda_under <- calcLike(negLL.fn = negLL.PLB,
                           x = x_under_vector,
                           xmin = min(x_under_vector), 
                           xmax = max(x_under_vector), 
                           n = length(x_under_vector), 
                           sumlogx = sum(log(x_under_vector)), 
                           p = -1.5,
                           suppress.warnings = TRUE,
                           vecDiff = vecDiff)
  
  # lambda trimmed ####
  # estimate lambda from trimmed data
  x_trimmed_vector <- x_under |>
    filter(xl >= cutoff) |>
    pull(x)
  trimmed_n <- length(x_trimmed_vector)
  lambda_trimmed <- calcLike(negLL.fn = negLL.PLB,
                             x = x_trimmed_vector,
                             xmin = min(x_trimmed_vector), 
                             xmax = max(x_trimmed_vector), 
                             n = length(x_trimmed_vector), 
                             sumlogx = sum(log(x_trimmed_vector)), 
                             p = -1.5,
                             suppress.warnings = TRUE,
                             vecDiff = 2)
  # end/return ####
  # return data frame ####
  out_df <- data.frame(
    known_lambda = lambda, 
    lambda_under = lambda_under$MLE,
    lambda_under_lo = lambda_under$conf[1],
    lambda_under_hi = lambda_under$conf[2],
    lambda_trimmed = lambda_trimmed$MLE,
    lambda_trimmed_lo = lambda_trimmed$conf[1],
    lambda_trimmed_hi = lambda_trimmed$conf[2],
    original_n = n,
    under_n = under_n,
    trimmed_n = trimmed_n,
    M = M,
    cutoff = cutoff, 
    xmin_obs = xmin_obs,
    xmin_under = min(x_under_vector),
    xmin_trimmed = min(x_trimmed_vector),
    xmax_obs = xmax_obs)
  return(out_df)
}




