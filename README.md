# VFP-Skewed-Power

###   

Investigating how under sampling small body sizes affects estimates of power law exponents (lambda). The repository name (Very Fast Picket Skewed Power) is inspired by ship names in the Culture Series by Iain M. Banks. 

## Introduction  

Power laws describe many phenomena in the natural and biological sciences. In ecology, a common pattern in biological communities is the rapid decline in abundance with increasing organismal body size. This pattern has received much attention in the literature and is known as size spectra, mass-abundance relationships, or community biomass distributions. Here, we are specifically interested in constructing individual size distributions (ISD) as defined by [White et al. 2007.](https://www.cell.com/ajhg/abstract/S0169-5347(07)00098-5). ISD are a frequency distribution in the form of $f(m) \propto M^\lambda$, where $\lambda$ controls the rate of decline in abundance of large individuals. 

Empirical observations which are thought to follow a power law distribution, such as ISD, often display under sampling of the smallest values, often attributed to artefacts from the collection techniques. For example, when sampling benthic macroinvertebrates, the mesh size (typically 250-500 um) used has a range of sizes which it is optimized for. Small individuals are able to get through the mesh without being captured. Technically any individual which is larger than the net openings should be captured with high efficiency. Practically, though, large individuals are unlikely to occur in the relatively small areas typically sampled. For example, benthic macroinvertebrates are usually collected in areas < 1 m^2^, where as fish collections typically occur on the scale of 100-1000's m^2^. For small individuals, there is generally a sigmoidal response, where very small individuals have a very low probability of being retained in the net. This retention probability increases rapidly as individuals get larger and asymptotically approaches 100% probability. This functional response of retention probability leads to lower densities of small individuals than their true population densities would suggest.

Although the effect of under sampling the smallest individuals has been recognized widely in the ISD literature, it remains unclear how much estimates of power law exponents are affected by under sampling of the smallest body sizes. 

Here, we use a simulation framework with repeated sampling to explore how under sampling may affect estimates of power law exponents. We sample body size values from a known distribution and then artificially under sample some range of small body sizes to investigate the deviance of estimates from the known values. Furthermore, we explore methods for correcting bias in estimates through censoring the data in two different ways (Figure 1). Finally, we also simulate a known relationship between $\lambda$ and a hypothetical gradient and investigate the deviance between the known and estimated relationship using biased and censored data (Figure 2). 

![Figure 1. Conceptual figure showing the simulation framework. ](plots/conceptual-bias-correction.png)
![Figure 2. Conceptual figure showing the relationship across a hypothetical gradient. ](plots/gradient_concept.png)

## Basic workflow  

1. Sample body mass observations, $m$, from a bounded power law with a known value of $\lambda$. Body masses were bound from 0.001 to 10 000 mg dry mass. This size range was based on benthic macroinvertebrate samples collected at a continental scale [(Pomeranz et al., 2022)](https://doi.org/10.1111/gcb.15862). 
2. Convert body masses ($m$) to body lengths, ($l$). 
3. Calculate retention probability as a function of body lengths. Retention probability is calculated with the logistic regression model of [Morin et al. (2004)](https://www.journals.uchicago.edu/doi/10.1899/0887-3593(2004)023%3C0383:SRPOSB%3E2.0.CO;2). Briefly, this model is a function of relative length (RL) and mesh size (φ) in millimeters. Relative length is the ratio of body length to the mesh opening size. Here, we used mesh sizes of 0.125, 0.25, 0.5 and 1 mm to represent realistic sampling protocols used for freshwater benthic macroinvertebrates. The logistic model has a sigmoid function such that small body sizes have an extremely low probability of being sampled which rises rapidly with increasing body size and saturates around a probability of ~0.999. This results in small body sizes having low retention probabilities, while large body sizes have high retention probabilities.
4. Create a biased data set by performing Bernoulli trials according to the retention probabilities calculated in step 3. 
5. Create the censored data by removing all individual observations which are smaller than the cutoff size. The minimum body size is determined in two ways: 1) estimated using the method described by [Clauset et al., (2009)](https://doi.org/10.1137/070710111) and implemented as the `estimate_xmin()` algorithm in the [*poweRlaw* package](https://cran.r-project.org/web/packages/poweRlaw/index.html) and 2) using fixed cutoffs where xmin is equal to 10x the mesh size opening e.g., 1.25, 2.5, 5 and 10 mm body length, respectively.  
6. Estimate λ from the biased data and from the two censored data sets and compare them with known values of λ.
 
## Specific workflow  
*Update this with names/sequence of scripts to run*
