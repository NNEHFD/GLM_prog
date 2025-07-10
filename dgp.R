# -------------------------------------------------------------------------------------------
# Purpose: Script to generate data sets for different DPGs
#
# Dependencies: none
#
# Output: Functions for generating datasets for different DPGs, i.e. outcome_dpg and dgp
#
# -------------------------------------------------------------------------------------------

# Load libraries ---------------------------------------------------------
library(data.table)
library(tidyverse)
library(magrittr)
library(MASS)
library(postcard)


# Define link function  ------------------------------------------------------------
# define parameters
zeta <- 0.2425194

# define link function
inverse_canonical_link <- function(x) {
  exp(x)
}

link <- function(x) {
  log(x)
}

#help function
add_bernoulli <- function(x) {
  integer_part <- floor(x)  # Extract the integer part
  fractional_part <- x - integer_part  # Extract the fractional part
  result <- integer_part + rbinom(length(x), size = 1, prob = fractional_part)  # Perform the Bernoulli trial
  return(result)
}



# Outcome specification for DGP ---------------------------
outcome_dgp <- function(dgp_string = 'constant') {
  if (dgp_string == 'constant') {
    mu0 <<- function(X) {
      X %$%
        {
          4.1 * sin((abs(W2))) +
            2 * W4 +
            as.numeric(abs(W3) > 2.5) * 1.4 +
            as.numeric(abs(W4) > 0.25) * 1.5 +
            1.5 * sin(abs(W5)) -
            as.numeric(W1 < -4.1) * 5 * sin(abs(W2)) -
            as.numeric(W1 < -6.1) * 3 * sin(abs(W2)) -
            as.numeric(U > 1.1) * 5 * sin(abs(W2)) -
            as.numeric(U > 1.55) * 3 * sin(abs(W2))
        } %>% 
        abs()
        }
    
    mu1 <<- function(X) {
      X %$%
        {
          inverse_canonical_link(zeta + link(mu0(.)))
        }
    }
    
    
  } else if (dgp_string == 'heterogeneous') {
    mu0 <<- function(X) {
      X %$%
        {
          4.1 * sin((abs(W2))) +
            2 * W4 +
            as.numeric(abs(W3) > 2.5) * 1.4 +
            as.numeric(abs(W4) > 0.25) * 1.5 +
            1.5 * sin(abs(W5)) -
            as.numeric(W1 < -4.1) * 5 * sin(abs(W2)) -
            as.numeric(W1 < -6.1) * 3 * sin(abs(W2)) -
            as.numeric(U > 1.1) * 5 * sin(abs(W2)) -
            as.numeric(U > 1.55) * 3 * sin(abs(W2))
        } %>% 
        abs()
    }
    
    mu1 <<- function(X) {
      X %$%
        {
          5.3 * sin((abs(W2)))^2 +
            2.6 * W4 +
            as.numeric(abs(W3) > 2.5) * 1.4 +
            as.numeric(abs(W4) > 0.25) * 1.3 +
            4.1 * as.numeric(W2 > 0) * sin(abs(W5)) +
            1.6 * sin(abs(W6)) -
            as.numeric(W1 < -4.1) * 5 * sin(abs(W2)) -
            as.numeric(W1 < -6.1) * 3 * sin(abs(W2)) -
            as.numeric(U > 1.1) * 5 * sin(abs(W2)) -
            as.numeric(U > 1.55) * 3 * sin(abs(W2))
          
        } %>% 
        abs()
    }
    
  }
  
  invisible()
}


# data generating process -------------------------------------------------
dgp = function(n,
               p = 7,
               noise = 0.5,
               shift_W1 = 0,
               shift_U = 0) {
  p = pmax(p, 6)
  
  # Generate  observed covariates
  W1_lbound = -2 + shift_W1
  W1_ubound = 1 + shift_W1
  U_lbound = 0 + shift_U
  U_ubound = 1 + shift_U
  
  W1 <- runif(n, min = W1_lbound, max = W1_ubound)
  W2 <- runif(n, min = -2, max = 1)
  W3 <- rnorm(n, 0, 3)
  W4 <- rexp(n, .8)
  W5 <- rgamma(n, 5, 10)
  
  # Generate unobserved covariates
  U <- runif(n, min = U_lbound, max = U_ubound)
  
  # Generate the full dataset
  matrix(runif(n * (p - 5), 1, 2), ncol = (p - 5)) %>% # Generate covariates W6-Wp
    as.data.frame(check.names = FALSE) %>%
    set_names(str_c("W", 6:(p))) %>%
    tibble(U, W1, W2, W3, W4, W5) %>%
    inset("m0", value = mu0(.)) %>%
    inset("m1", value = mu1(.)) %>%
    mutate(
      A = rbinom(n, 1, prob = 1 / 2),
      Y0 = add_bernoulli(m0) + rbinom(n, size = floor(m0), prob = noise) - rbinom(n, size = floor(m0), prob = noise),
      Y1 = add_bernoulli(m1) + rbinom(n, size = floor(m1), prob = noise) - rbinom(n, size = floor(m1), prob = noise),
      Y = ifelse(A, Y1, Y0)
    )
}

# Generate data test ---------------------------------------------------------

outcome_dgp()

dat <- dgp(n = 100)



outcome_dgp(dgp_string = 'heterogeneous')

dat1 <- dgp(n = 10000)

dat1 %>% filter(Y == 0)
dat1 %>% filter(A == 1 & Y == 0)
dat1 %>% filter(A == 0 & Y == 0)

# Check the skewness ----------------------------------------------------------

ggplot(dat1 %>% filter(A == 1), aes(x = Y)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Histogram of Variable Y", x = "Y", y = "Frequency")


ggplot(dat1 %>% filter(A == 0), aes(x = Y)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Histogram of Variable Y", x = "Y", y = "Frequency")

ggplot(dat1, aes(x = Y)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Histogram of Variable Y", x = "Y", y = "Frequency")


library(e1071)
skewness(dat$Y)

# 
# n <- 6800000
# outcome_dgp(dgp_string = 'heterogeneous')
# dat <- dgp(n,
#           p = 7,
#           shift_W1 = 0,
#           shift_U = 0) %>% as.data.frame()
# 
# 
# RR <- rctglm(formula = Y ~ A,
#             exposure_indicator = A,
#             exposure_prob = 1/2,
#             data = dat,
#             family = negative.binomial(3, link = "log"),
#             estimand_fun = "rate_ratio")
# 
# 
# estimand(RR)
