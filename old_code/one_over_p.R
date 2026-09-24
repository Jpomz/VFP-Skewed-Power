# one_over_p
# exploring the Morin correction

library(sizeSpectra)
library(tidyverse)

n = 10000
lambda = -2
xmin = 0.01
xmax = 100
M = .5
vecDiff = 2
LWa = 0.0064
LWb = 2.788
  
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
  filter(fill == "sampled") |>
  mutate(weight = sampled)
# how many body sizes were sampled?
under_n <- nrow(x_under)

plot_x_df <- tibble(x = x, 
                    fill = "Original") |>
  mutate(weight = 1)

inv_df <- x_under |>
  mutate(inv_n = (1 / probability),
         fill = "inverse",
         weight = inv_n)

inv_df_plot <- inv_df |>
  select(x, fill, weight)

plot_df <- bind_rows(
  x_under |>
    select(x, fill, weight),
  plot_x_df) |>
  bind_rows(inv_df_plot)

plot_df |> 
  group_by(fill) |>
  count()

ggplot(plot_df,
       aes(x = x, 
           weight = weight,
           fill = fill)) +
  geom_histogram(binwidth = 0.1, 
                 position = "dodge",
                 alpha = 1) +
  scale_fill_manual(values = c("springgreen",
                               "black",
                               "#FF1984"))+
  scale_x_log10(guide = "axis_logticks") +
  # scale_y_log10() +
  theme_classic()

(inv_lambda <- calcLike(
  negLL.fn = negLL.PLB.counts,
  x = inv_df$x,
  c = inv_df$weight,
  p = -1.5,
  suppress.warnings = TRUE,
  vecDiff = vecDiff))

(orig_lambda <- calcLike(
  negLL.fn = negLL.PLB.counts,
  x = plot_x_df$x,
  c = plot_x_df$weight,
  p = -1.5,
  suppress.warnings = TRUE,
  vecDiff = vecDiff))

(under_lambda <- calcLike(
  negLL.fn = negLL.PLB.counts,
  x = x_under$x,
  c = x_under$weight,
  p = -1.5,
  suppress.warnings = TRUE,
  vecDiff = vecDiff))

# ggplot(inv_df,
#        aes(x = x,
#            weight = inv_n)) +
#   geom_histogram(binwidth = 0.1, 
#                  position = "dodge",
#                  alpha = 1) +
#   scale_fill_manual(values = c("springgreen", "black", "#FF1984"))+
#   scale_x_log10(guide = "axis_logticks") +
#   # scale_y_log10() +
#   theme_classic()
under_bin <- binData(x = x_under$x, binWidth = "2k")
LBN_bin_plot(under_bin$binVals,
             b.MLE = under_lambda$MLE, 
             b.confMin = under_lambda$conf[1], 
             b.confMax = under_lambda$conf[2])


num_bins <- nrow(under_bin$binVals)
bin_breaks <- c(dplyr::pull(under_bin$binVals, binMin),
                dplyr::pull(under_bin$binVals, binMax)[num_bins])

bin_counts <- dplyr::pull(under_bin$binVals, binCount)

calcLike(negLL.PLB.binned,
                          p = -1.5,
                          w = bin_breaks,
                          d = bin_counts,
                          J = length(bin_counts),   # = num.bins
                          vecDiff = 1,
                          suppress.warnings = TRUE)     



# correcting tvf data? ####
tvf <- read_csv("empirical_data_examples/tvf_dw.csv") |>
  filter(!is.na(dw))

tvf <- tvf |> 
  mutate(prob = plogis(morin_ln_p(L = Value, M = 0.5)))

tvf <- tvf |>
  mutate(weight = 1/prob,
         weight_1 = 1)

tvf_01 <- tvf |>
  filter(site_number == "03")
calcLike(
  negLL.fn = negLL.PLB.counts,
  x = tvf_01$dw,
  c = tvf_01$weight_1,
  p = -1.5,
  suppress.warnings = TRUE,
  vecDiff = vecDiff)
calcLike(
  negLL.fn = negLL.PLB.counts,
  x = tvf_01$dw,
  c = tvf_01$weight,
  p = -2,
  suppress.warnings = TRUE,
  vecDiff = vecDiff)
