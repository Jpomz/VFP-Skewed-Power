# back and forward estimating counts of body sizes

estimate_pareto_N = function(n, lambda, xmin, xmin2, xmax){
  lambdaPlus = lambda + 1
  n * (xmax^(lambdaPlus) - xmin2^(lambdaPlus)) /
    (xmax^(lambdaPlus) - xmin^(lambdaPlus))
}

N <- 1:10
xmin <- seq(0.001, 100, length.out = 10)
estimate_pareto_N(n=N, -2, xmin = xmin, 0.006, 100)

pareto_expectation = function(lambda, xmin, xmax){
  if(lambda == -1) {
    # Special case λ = -1
    return((xmax - xmin) / log(xmax / xmin))
  } else if(lambda == -2) {
    # Special case λ = -2
    return((xmin * xmax / (xmin - xmax)) * log(xmax / xmin))
  } else {
    # General case
    numerator <- (lambda + 1) * (xmax^(lambda + 2) - xmin^(lambda + 2))
    denominator <- (lambda + 2) * (xmax^(lambda + 1) - xmin^(lambda + 1))
    return(numerator / denominator)
  }
}

lambda <- seq(-1.99, -1.01, length.out = 10)
p_exp <- sapply(X = lambda, FUN = pareto_expectation, xmin = 0.006, xmax = 100)

plot(lambda, p_exp)  
