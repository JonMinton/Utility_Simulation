###################################################################################################
###################################################################################################

# Bespoke functions: Load these first

#####
# convert from RR with CIs to simulations

# This function takes the central estimate plus lower and upper 95% CIs of a log normal distribution and 
# produces 10000 simulated values from this distribution.
rr_simulate <- function(central, lower, upper, sims=10000){
  mu <- log(central)
  sigma_low <- (1/ 1.96) * (mu - log(lower))
  sigma_high <- (1/ 1.96) * (log(upper) - mu )
  sigma <- (sigma_low + sigma_high) / 2
  
  out <- list(
    sims = exp(rnorm(n=sims, mean=mu, sd=sigm)),
    mu=mu, sigma=sigma
    )
  return(out)
}




# This function produces bootstrapped CIs of means of a vector
bootstrap <- function(inputs, sims = 10000){
  out <- vector("numeric", sims)
  n_inputs <- length(inputs)
  for (i in 1:sims) {
    out[i] <- mean(inputs[sample(1:n_inputs, replace=T)])
  }
  return(out)
}

