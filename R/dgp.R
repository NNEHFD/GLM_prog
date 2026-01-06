# -------------------------------------------------------------------------------------------
# Purpose: Script to generate data sets for different DPGs
#
# Dependencies: none
#
# Output: Functions for generating datasets for different DPGs, i.e. outcome_dpg and dgp
#
# -------------------------------------------------------------------------------------------

# Load libraries ---------------------------------------------------------
library(tidyverse)
library(magrittr)

# Define link functions --------------------------------------------------
inverse_canonical_link <- function(x) exp(x)
link <- function(x) log(pmax(x, 1e-6))
rem <- function(x) x - floor(x)

E_abs_norm <- function(mu) {
  sqrt(2/pi) * exp(-0.5 * mu^2) + abs(mu) * (2 * pnorm(abs(mu)) - 1)
}

# Data generating process -------------------------------------------------
dgp <- function(
                shift_W1 = 0,
                shift_U = 0,
                zeta = 0.1,
                const_effect = F) {
  structure(
    list(
      shift_W1 = shift_W1,
      shift_U = shift_U,
      zeta = zeta,
      const_effect = const_effect
    ),
    class = "dgp"
  )
}

draw <- function(x, ...) UseMethod("draw")

hinge = \(x, x0=0) (x-x0)*(x > x0)

draw.dgp <- function(dgp, n, ...) {
  
  # Generate covariates
  data = tibble(
    U = rnorm(n, mean = dgp$shift_U),
    W1 = rnorm(n, mean = dgp$shift_W1),
    W2 = rnorm(n),
    W3 = rnorm(n),
    W4 = rnorm(n),
    W5 = rnorm(n)
  )
 
  # Generate outcomes
  data |>
    mutate(
      # main = 0.1 + 3/2*hinge(W1) + W2^2 + hinge(W1*W4),
      # inter = 3/2*hinge(W3,-2),
      # het = 2*W4^2,
      main = 0.1 + 3/2*hinge(W1,-1) + W2^2 + hinge(W1*W4),
      inter = hinge(W3,-2),
      het = 2*hinge(W4),
      m0 = main + abs(U) * inter,
      mu0 =  main + E_abs_norm(dgp$shift_U) * inter,
      m1 = if (dgp$const_effect) {
        inverse_canonical_link(dgp$zeta + link(m0))
      } else {
        inverse_canonical_link(dgp$zeta + link(m0 + het))
      },
      A = rbinom(n, 1, prob = 0.5),
      Y0 = rbinom(n, size = 2*floor(m0), prob = 1/2) + rbinom(n, 1, prob = rem(m0)),
      Y1 = rbinom(n, size = 2*floor(m1), prob = 1/2) + rbinom(n, 1, prob = rem(m1)),
      Y = if_else(A == 1, Y1, Y0)
    ) |>
    dplyr::select(-main, -inter, -het)
}

params <- function(x, ...) UseMethod("params")
params.dgp = function(dgp) {
  tibble(
    shift_W1 = dgp$shift_W1,
    shift_U = dgp$shift_U,
    zeta = dgp$zeta,
    const_effect = dgp$const_effect
  )
}

RR = function(x, ...) UseMethod("RR")
RR.dgp = function(dgp, n=1e6) {
  dgp |> draw(n) %$% {
    mean(m1) / mean(m0)
  }
}

r2 = function(x, ...) UseMethod("r2")
r2 = function(dgp, n=1e6) {
  data = dgp |> 
    draw(n) |>
    mutate(m = ifelse(A==1, m1, m0)) |>
    dplyr::select(m, Y, A, starts_with('W'))
  m_nonlin = earth::earth(m ~ . - Y, degree=3, nprune=50, data) |> 
    predict(data)
  m_lin = lm(m ~ . - Y, data) |> 
    predict(data)
  data %$% list(
    r2 = {var(m) / var(Y)},
    lin_r2 = {var(m_lin) / var(Y)},
    nonlin_r2 = {var(m_nonlin) / var(Y)},
    nonlin_rmse = sqrt(mean((Y - m_nonlin)^2)),
    bayes_rmse = sqrt(mean((Y - m)^2)),
    lin_rmse = sqrt(mean((Y - m_lin)^2))
  )
}

# trial = dgp(const_effect = T, zeta=0.1)
# c(lrnr_results, power_stats, lrnr) %<-% fit_prog(trial, 5000, mu_wf=MU_WF)
# as_tibble(power_stats) |>
#   mutate(rr= RR(trial))

# trial = dgp(zeta=-0.05)
# c(lrnr_results, power_stats, lrnr) %<-% fit_prog(trial, 5000, mu_wf=MU_WF)
# as_tibble(power_stats) |>
#   mutate(rr= RR(trial))