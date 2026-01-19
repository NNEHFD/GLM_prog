library(tidymodels)
library(magrittr)
library(tidyverse)
library(furrr)

dRR_1 <- function(psi_0) {
  1/psi_0
}

dRR_0 <- function(psi_1, psi_0) {
  -psi_1/psi_0^2
}

RR_variance_bound <- function(
  psi0, # E[Y(0)]
  sigma2, # V[Y(0)]
  kappa, # E[ (Y(0) - mu0(X))^2 ] ^ 1/2
  target_effect, # risk ratio E[Y(1)] / E[Y(0)] of interest
  pi = 1/2 # E[A]
) {
  psi1 = target_effect*psi0 

  # Find the derivatives 
  d1 <- dRR_1(psi0) #rate ratio specific 
  d0 <- dRR_0(psi0, psi1) # rate ratio specific
  
  # Calculate the variance bound
  d0^2 * sigma2 + d1^2 * sigma2 + pi * (1 - pi) * ((abs(d0) * kappa / (1 - pi) + abs(d1) * kappa / pi)^2)
}


# function to determine sample size from variance bound ----------------------------------
calc_power <- function(target_effect, variance, n, alpha=0.05){
  sd <- sqrt(variance / n)
  # Critical values under H0 (mean=1)
  c_upper <- qnorm(1 - alpha/2, mean = 1, sd = sd)
  c_lower <- qnorm(alpha/2, mean = 1, sd = sd)
  
  # Power under H1 (mean=target_effect)
  p_upper <- 1 - pnorm(c_upper, mean = target_effect, sd = sd)
  p_lower <- pnorm(c_lower, mean = target_effect, sd = sd)
  
  p_upper + p_lower
}

sample_size_glm <- function(target_effect, 
                            variance,
                            alpha = 0.05,
                            power = 0.8) {
  # Analytic solution for normal approximation
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)
  n <- variance * (z_alpha + z_beta)^2 / (target_effect - 1)^2
  return(n)
}

sample_sizes =function(true_rr, power_stats, ...) {
  power_stats %$% list(
    'oracle-adj. GLM' = kappa_oracle,
    'prog-adj. GLM' = kappa_ml,
    'GLM' = kappa_glm,
    unadjusted = sqrt(sigma2) 
  ) |> imap(\(kappa, name) {
    tibble(
      estimator = name,
      n = sample_size_glm(
        target_effect=true_rr, 
        variance = power_stats %$% RR_variance_bound(psi0, sigma2, kappa, target_effect=true_rr),
        ...
      )
    )
  }) |> bind_rows()
}

