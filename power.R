# -------------------------------------------------------------------------------------------
# Purpose: Create function for determining the variance bound
#
# Dependencies: 
#
# Output: Function for determining the variance bound
#
# -------------------------------------------------------------------------------------------

# Load libraries ------------------------------------------------------------------
db <- NNaccess::nnaccess()

library(tidymodels)
library(magrittr)
library(tidyverse)
library(furrr)


#derivatives

dRR_1 <- function(psi_0) {
  1/psi_0
}

dRR_0 <- function(psi_1, psi_0) {
  -psi_1/psi_0^2
}


# function to determine the variance bound  ----------------------------------------------

variance_bound <- function(data, 
                           target_effect, 
                           pi, #probability of treatment
                           hatmu,
                           V,
                           pull = FALSE) {
  
  
  
  # Counterfactual mean determination
  psi_0 = data %$% {
    mean(Y) #dont filter for A=0 since all of the historical data have A=0
  }
  psi_1 = target_effect*psi_0 #specific for rate ratio
  
  # Marginal outcome variance 
  var_Y_0 = data %$% {
    var(Y) #dont filter for A=0 since all of the historical data have A=0
  }
  var_Y_1 = var_Y_0 #since we don't know anything about the structure as we would if data was binary
  
  # Find the derivatives 
  d1 <- dRR_1(psi_0) #rate ratio specific 
  d0 <- dRR_0(psi_0, psi_1) # rate ratio specific
  
  # CV expected MSE from best possible GLM fit 
  cv <- data$Y %>%
    createFolds(k = V, list = FALSE)
  
  data$predictions <- NA
  
  for (i in 1:V) {
    train.data  <- data[cv != i, ]
    test.data <- data[cv == i, ]
    # Build the model
    fit <- hatmu(train.data)
    
    if (pull == TRUE) {
      data$predictions[cv == i] <- predict(fit, test.data) %>% 
        pull(.pred)
    } else {
      data$predictions[cv == i] <- predict(fit, test.data, type = "response")
    }
    
    }
  
  
  
  # Calculate the estimate of kappa0^2
  n_tilde <- nrow(data)
  kappa0_squared <- 1/n_tilde * sum((data$Y - data$predictions)^2)
  kappa1_squared <- kappa0_squared
  
  # Calculate the variance bound
  v_bound <- d0^2 * var_Y_0 + d1^2 * var_Y_1 +
    pi * (1 - pi) * ((abs(d0) * sqrt(kappa0_squared) / (1 - pi) + abs(d1) * sqrt(kappa1_squared) / pi)^2)
  
  return(v_bound)
}


# function to determine sample size from variance bound ----------------------------------
calc_power <- function(alpha, target_effect, variance, n){
   f0 <- qnorm(1 - alpha/2, mean = 1, sd = sqrt(variance / n))
   f1 <- pnorm(f0, mean = target_effect, sd = sqrt(variance / n))
   1 - f1
}


sample_size_glm <- function(target_effect, 
                        variance,
                        alpha = 0.05,
                        gamma = 0.8,
                        initial_n = 1,
                        increment = 2) {
  power <- 0
  n <- initial_n
  
  while (power < gamma) {
    power <- calc_power(alpha, target_effect, variance, n)
    n <- n + increment
  }
  
  tibble(n = n)
}




