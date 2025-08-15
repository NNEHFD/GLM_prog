# -------------------------------------------------------------------------------------------
# Purpose: Plots for results 
#
# Dependencies: 
#
# Output: 
#
# -------------------------------------------------------------------------------------------

# load libraries and dependencies -----------------------------------------------------
set.seed(34264524)

library(tidymodels)
library(magrittr)
library(tidyverse)
library(furrr)
library(caret)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggpubr)

source("statprog/GLM/dgp.R") # to define zeta as in the dpg.R file
source("statprog/GLM/experiment_power.R")
source("statprog/GLM/power.R")

# Load datasets -----------------------------------------------------------------------
db <- NNaccess::nnaccess(project = "students", trial = "ehfd_phd", instance = "current")

data_constant <- db$output_datasets("data_constant", ext = "rds")
data_het <- db$output_datasets("data_het", ext = "rds")
data_obs_shift <- db$output_datasets("data_obs_shift", ext = "rds")
data_unobs_shift <- db$output_datasets("data_unobs_shift", ext = "rds")

# Vary both loading
inc <- c(20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 225, 250, 275, 300)

data_vary_both <- db$output_datasets(paste0("data_vary_both_", inc[1]), ext = "rds")

for (i in 2:length(inc)) {
  data_vary_both <- rbind(data_vary_both,
                          db$output_datasets(paste0("data_vary_both_", inc[i]), ext = "rds") )
}


# -------------------------------------------------------------------------------------------
# Determining the true rate ratio in the two scenarios
# -------------------------------------------------------------------------------------------

# Constant treatment effect means log(mu1) = zeta + log(mu0) = zeta + eta0 implying
# RR = E[mu1]/E[mu0] = E[exp(zeta + eta0)]/E[exp(eta0)] = E[exp(zeta)exp(eta0)]/E[exp(eta0)] = exp(zeta)

truth_constant <- exp(zeta)

# For the heterogeneous treatment effect scenario the true RR should be found
# by using a large sample and using the law of large numbers since RR = E[mu1]/E[mu0]
n <- 6820000
outcome_dgp(dgp_string = 'heterogeneous')
dat <- dgp(n,
           p = 7,
           shift_W1 = 0,
           shift_U = 0)

RR_old <- rctglm(formula = Y ~ A,
                 exposure_indicator = A,
                 exposure_prob = 1/2,
                 data = dat,
                 family = "gaussian",
                 estimand_fun = "rate_ratio",
                 verbose = 0,
                 cv_variance = FALSE)$estimand$Estimate
norm_diff <- 1

# Set up loop to increase data points
while (norm_diff > 0.0001) {
  n <- n + 10000
  dat_new <- dgp(n,
                 p = 7,
                 shift_W1 = 0,
                 shift_U = 0)
  
  # Fit the unadjusted estimator
  RR_new <- rctglm(formula = Y ~ A,
                   exposure_indicator = A,
                   exposure_prob = 1/2,
                   data = dat_new,
                   family = "gaussian",
                   estimand_fun = "rate_ratio",
                   verbose = 0,
                   cv_variance = FALSE)$estimand$Estimate
  
  # Compute norm difference between old and new coefficients
  norm_diff <- sqrt(sum((RR_new - RR_old)^2))
  
  # update old coefficients
  RR_old <- RR_new
}

list <- list(coefficients = RR_old, n = n)

truth_het <- list$coefficients

# -------------------------------------------------------------------------------------------
# Determining the oracle variance in the two scenarios
# -------------------------------------------------------------------------------------------

# Determining the oracle variance for the rate ratio in the two scenarios by using the same plug-in estimate method 
# Constant treatment effect using n from the determination of the true het effect (n=180) 

future::plan(multisession, workers = 90)

# marginal effect as the rate ratio
RR <- function(psi_1, psi_0) {
  psi_1/psi_0
}

# derivatives of marginal effect
dRR_1 <- function(psi_0) {
  1/psi_0
}

dRR_0 <- function(psi_1, psi_0) {
  -psi_1/psi_0^2
}


mean_IF_1 <- function(A, Y, mu1, pi, psi_1) {
  (A/pi * (Y - mu1) + (mu1 - psi_1))
}

mean_IF_0 <- function(A, Y, mu0, pi, psi_0) {
  ((1 - A)/(1 - pi) * (Y - mu0) + (mu0 - psi_0))
}


marginaleffect_IF <- function(data, pi, psi_1, psi_0) {
  data %$% {
    dRR_1(psi_0)*mean_IF_1(A, Y, m1, pi, psi_1) + dRR_0(psi_1, psi_0)*mean_IF_0(A, Y, m0, pi, psi_0)
  }
}


oracle_var <- function(n, dgp_string) {
  outcome_dgp(dgp_string = dgp_string)
  dat <- dgp(n = n,
             p = 7,
             shift_W1 = 0,
             shift_U = 0)

  psi_1 = dat %$% {mean(m1)}
  psi_0 = dat %$% {mean(m0)}

  IC <- marginaleffect_IF(dat, pi = 1/2, psi_1, psi_0)
  variance = var(IC)
  se = sqrt(variance/nrow(dat)) %>% as.data.frame()

}


N <- 1000

params = expand_grid(
  n = c(180),
  dgp_string = c("constant"),
)
results = 1:N %>% future_map( ~ params %>% pmap_df(oracle_var), .options = furrr_options(seed = 3021377))
results %<>% purrr::list_rbind()

oracle_se_constant <- results %>% summarize(est_eff = mean(.))

# Heterogeneous treatment
params = expand_grid(
  n = c(180),
  dgp_string = c("heterogeneous")
)
results = 1:N %>% future_map( ~ params %>% pmap_df(oracle_var), .options = furrr_options(seed = 3021377))
results %<>% purrr::list_rbind()

oracle_se_het <- results %>% summarize(est_eff = mean(.))



# -------------------------------------------------------------------------------------------
# Sample size from variance bound for the two scenarios (power = 0.8)
# -------------------------------------------------------------------------------------------
N <- 100

#Heterogeneous treatment effect
params = expand_grid(
  n_hist = c(4000),
  p = c(7),
  shift_W1 = c(0), # shift in historical data only
  shift_U = c(0), # shift in historical data only
  dgp_string = c("heterogeneous"),
  target_effect = truth_het,
  pi = 1/2,
  r = 3,
  alpha = 0.025,
  gamma = 0.9,
  initial_n = 20,
  increment = 2
)
results = 1:N %>% future_map_dfr( ~ params %>% pmap_df(experiment_power), .options = furrr_options(seed = 3021377))


db$output_datasets("data_ss_het", results, ext = "rds")

ss <- db$output_datasets("data_ss_het", ext = "rds")
ss_het <- ss %>% group_by(estr, prog) %>% summarize(est_eff = ceiling(mean(n))) 

future::plan(sequential)
# -------------------------------------------------------------------------------------------
# Defining a function for finding the results
# -------------------------------------------------------------------------------------------


res <- function(results, truth = truth_het, alpha = 0.05) {
  results %>%
    summarize(
      est_eff = mean(Estimate),
      emp_se = sd(Estimate),
      bias = mean(Estimate - truth),
      #bias_thing = mean((Estimate - truth)/emp_se),
      mean_est_se = mean(`Std. Error`),
      rmse = sqrt(mean((Estimate - truth) ^ 2)),
      power = mean(2*(1 - pnorm(abs((Estimate - 1) / `Std. Error`))) < alpha),
      coverage = mean(((qnorm(alpha/2) * `Std. Error` + Estimate) < truth) &
                        (truth < (-qnorm(alpha/2) * `Std. Error` + Estimate)))
      
    )
}
# above .1 then bias affecting the coverage  
# Define the number of simulations 
N <- 500

# -------------------------------------------------------------------------------------------
# Table with all four scenarios including:
# - Scenario
# - method for analysing data
# - True rate ratio
# - Empirical standard error
# - Mean estimated standard error
# - RMSE
# - Power
# - Coverage
# -------------------------------------------------------------------------------------------
results <- rbind(data_het %>% mutate(scenario = 'het'), 
                 data_obs_shift %>% mutate(scenario = 'obs_shift'),
                 data_unobs_shift %>% mutate(scenario = 'unobs_shift'))

tab <- results %>% group_by(scenario, estr, prog, shift_W1, shift_U) %>% res()

tab %<>% rbind(data_constant %>%
                 mutate(scenario = 'constant') %>% 
                 group_by(scenario, estr, prog, shift_W1, shift_U) %>% 
                 res(truth = truth_constant))

tab <- tab %>%
  mutate(scenario = case_when(
    scenario == "het" ~ "Heterogeneous",
    scenario == "constant" ~ "Constant",
    scenario == "obs_shift" ~ "Observable shift",
    scenario == "unobs_shift" ~ "Unobservable shift"
  )) %>%
  mutate(estr = case_when(
    estr == "glm" & prog == "fit" ~ "GLM with Super Learner prognostic score",
    estr == "glm" & prog == "fit random" ~ "GLM with non-informative prognostic score",
    estr == "glm" & prog == "none" ~ "GLM",
    estr == "glm" & prog == "oracle" ~ "GLM with oracle prognostic score",
    estr == "unadjusted" & prog == "none" ~ "Unadjusted",
    TRUE ~ NA_character_
  )) %>%
  mutate(shift = shift_W1 + shift_U) %>%
  group_by(scenario) %>%
  dplyr::select(-prog, -shift_W1, -shift_U) %>%  # Removing the 'estr' and 'prog' columns
  dplyr::select(scenario, estr, shift, est_eff, bias, everything()) %>% 
  arrange(scenario, shift, match(estr, c("Unadjusted", "GLM", "GLM with non-informative prognostic score", "GLM with Super Learner prognostic score", "GLM with oracle prognostic score")))  # Reordering the dataframe based on 'estr_prog

print(xtable(tab), include.rownames = FALSE)


# -------------------------------------------------------------------------------------------
# Plot of observed and unobserved change both small and large
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# Custom labels and colors
# -------------------------------------------------------------------------------------------

custom_labels <- c(
  "unadjusted" = "Unadjusted",
  "GLM" = "GLM",
  "GLM_w/prog" = "GLM with Super Learner prognostic score",
  "GLM_w/prog\nsmall obs. shift" = "GLM with Super Learner prognostic score\nsmall observed shift",
  "GLM_w/prog\nlarge obs. shift" = "GLM with Super Learner prognostic score\nlarge observed shift",
  "GLM_w/prog\nsmall unobs. shift" = "GLM with Super Learner prognostic score\nsmall shift",
  "GLM_w/prog\nlarge unobs. shift" = "GLM with Super Learner prognostic score\nlarge shift",
  "GLM_w/oracle" = "GLM with oracle prognostic score"
)

custom_colors <- c(
  "unadjusted" = "#FD8D3C",
  "GLM" = "#C51B8A",
  "GLM_w/prog" = "Light Sky Blue",
  "GLM_w/oracle" = "#31A354",
  "GLM_w/prog\nsmall obs. shift" = "#6BAED6",
  "GLM_w/prog\nlarge obs. shift" = "#08519C",  # New blue color
  "GLM_w/prog\nsmall unobs. shift" = "#6BAED6",
  "GLM_w/prog\nlarge unobs. shift" = "#08519C"  # New blue color
)


# -------------------------------------------------------------------------------------------
# Plot of observed change both small and large
# -------------------------------------------------------------------------------------------

results_obs <- rbind(
  data_het %>% filter(estr == "unadjusted"), 
  data_het %>% filter(estr == "glm" & prog == "none"),
  data_het %>% filter(estr == "glm" & prog == "fit"),
  data_obs_shift %>% filter(estr == "glm" & prog == "fit" & shift_W1 == -3),
  data_obs_shift %>% filter(estr == "glm" & prog == "fit" & shift_W1 == -5),
  data_het %>% filter(estr == "glm" & prog == "oracle")
)

results_obs <- data.frame(results_obs)

results_obs$estimator <- factor(c(
  rep('unadjusted', nrow(data_het %>% filter(estr == "unadjusted"))),
  rep('GLM', nrow(data_het %>% filter(estr == "glm" & prog == "none"))),
  rep('GLM_w/prog', nrow(data_het %>% filter(estr == "glm" & prog == "fit"))),
  rep('GLM_w/prog\nsmall obs. shift', nrow(data_obs_shift %>% filter(estr == "glm" & prog == "fit" & shift_W1 == -3))),
  rep('GLM_w/prog\nlarge obs. shift', nrow(data_obs_shift %>% filter(estr == "glm" & prog == "fit" & shift_W1 == -5))),
  rep('GLM_w/oracle', nrow(data_het %>% filter(estr == "glm" & prog == "oracle")))  # Place "GLM_w/oracle" last
), levels = c(
  'unadjusted',
  'GLM',
  'GLM_w/prog',
  'GLM_w/prog\nsmall obs. shift',
  'GLM_w/prog\nlarge obs. shift',
  'GLM_w/oracle'  # Ensure this is the last level
))

p1 <- ggplot(results_obs, aes(x = estimator , y = `Std..Error`, fill = estimator)) +
  geom_violin(trim = TRUE) +
  scale_fill_manual(values = custom_colors, labels = custom_labels) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black"),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  labs(
    title = "Observed shift",
    y = "Estimated standard error",
    x = "Estimator"
  )

# -------------------------------------------------------------------------------------------
# Plot of unobserved change both small and large
# -------------------------------------------------------------------------------------------

results_unobs <- rbind(
  data_het %>% filter(estr == "unadjusted"), 
  data_het %>% filter(estr == "glm" & prog == "none"),
  data_het %>% filter(estr == "glm" & prog == "fit"),
  data_unobs_shift %>% filter(estr == "glm" & prog == "fit" & shift_U == 0.5),
  data_unobs_shift %>% filter(estr == "glm" & prog == "fit" & shift_U == 1),
  data_het %>% filter(estr == "glm" & prog == "oracle")
)

results_unobs <- data.frame(results_unobs)

results_unobs$estimator <- factor(c(
  rep('unadjusted', nrow(data_het %>% filter(estr == "unadjusted"))),
  rep('GLM', nrow(data_het %>% filter(estr == "glm" & prog == "none"))),
  rep('GLM_w/prog', nrow(data_het %>% filter(estr == "glm" & prog == "fit"))),
  rep('GLM_w/prog\nsmall unobs. shift', nrow(data_unobs_shift %>% filter(estr == "glm" & prog == "fit" & shift_U == 0.5))),
  rep('GLM_w/prog\nlarge unobs. shift', nrow(data_unobs_shift %>% filter(estr == "glm" & prog == "fit" & shift_U == 1))),
  rep('GLM_w/oracle', nrow(data_het %>% filter(estr == "glm" & prog == "oracle")))  # Place "GLM_w/oracle" last
), levels = c(
  'unadjusted',
  'GLM',
  'GLM_w/prog',
  'GLM_w/prog\nsmall unobs. shift',
  'GLM_w/prog\nlarge unobs. shift',
  'GLM_w/oracle'  # Ensure this is the last level
))

p2 <- ggplot(results_unobs, aes(x = estimator, y = `Std..Error`, fill = estimator)) +
  geom_violin(trim = TRUE) +
  scale_fill_manual(values = custom_colors, labels = custom_labels) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black"),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  labs(
    title = "Unobserved shift",
    y = "Estimated standard error",
    x = "Estimator"
  )


# -------------------------------------------------------------------------------------------
# Combined plot
# -------------------------------------------------------------------------------------------
legend <- get_legend(p2)

# Combine the plots without legends
plot_row <- plot_grid(
  p1 + theme(legend.position = "none"), 
  p2 + theme(legend.position = "none"),
  labels = c('A', 'B'), label_size = 20, nrow = 1
)

# Add the title
title <- ggdraw() + 
  draw_label(
    "Standard error comparison with shifted covariates",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )

# Combine title, plot row, and legend
combined_plot <- plot_grid(
  title, plot_row, legend,
  ncol = 1,
  rel_heights = c(0.1, 1, 0.1)
)

# Display the combined plot
print(combined_plot)

db$exportOutput(combined_plot, "GLM_shift", Format = "pdf", FgHeight = 14, FgWidth = 25)
db$exportOutput(combined_plot, "GLM_shift", Format = "jpeg", FgHeight = 14, FgWidth = 25)



# -------------------------------------------------------------------------------------------
# performance diff scenarios 
# -------------------------------------------------------------------------------------------

results <- rbind(
  data_constant %>% mutate(scenario = 'constant'),
  data_het %>% mutate(scenario = 'het'), 
  data_obs_shift %>% filter(shift_W1 == -3) %>% mutate(scenario = 'obs_shift_small'),
  data_obs_shift %>% filter(shift_W1 == -5) %>% mutate(scenario = 'obs_shift_large'),
  data_unobs_shift %>% filter(shift_U == 0.5) %>% mutate(scenario = 'unobs_shift_small'),
  data_unobs_shift %>% filter(shift_U == 1) %>% mutate(scenario = 'unobs_shift_large')
) %>% mutate(estimator = case_when(
  estr == "unadjusted" ~ "unadjusted",
  estr == "glm" & prog == "none" ~ "GLM",
  estr == "glm" & prog == "fit random" ~ "GLM_w/rand_prog",
  estr == "glm" & prog == "fit" ~ "GLM_w/prog",
  estr == "glm" & prog == "oracle" ~ "GLM_w/oracle"
))

tab <- results %>% group_by(scenario, estimator, estr, prog, shift_W1, shift_U) %>% res()

tab$estimator <- factor(tab$estimator, levels = c(
  'unadjusted',
  'GLM',
  'GLM_w/rand_prog',
  'GLM_w/prog',
  'GLM_w/oracle'
))

tab$scenario <- factor(tab$scenario, levels = c(
  'constant',
  'het',
  'obs_shift_small',
  'obs_shift_large',
  'unobs_shift_small',
  'unobs_shift_large'
))

custom_labels <- c(
  "unadjusted" = "Unadjusted",
  "GLM" = "GLM",
  "GLM_w/rand_prog" = "GLM with non-informative prognostic score",
  "GLM_w/prog" = "GLM with Super Learner prognostic score",
  "GLM_w/oracle" = "GLM with oracle prognostic score"
)

custom_colors <- c(
  "unadjusted" = "#FD8D3C",
  "GLM" = "#C51B8A",
  "GLM_w/rand_prog" = "#6BAED6",
  "GLM_w/prog" = "#08519C",
  "GLM_w/oracle" = "#31A354"
)

# Create the plot
p <- ggplot(tab, aes(x = estimator, y = mean_est_se, fill = estimator)) +
  geom_point(aes(color = estimator), size = 3, position = position_dodge(width = .75)) +
  geom_point(aes(y = emp_se, color = estimator), size = 2, shape = 8, position = position_dodge(width = .75)) +
  facet_grid(. ~ factor(scenario), labeller = as_labeller(c(
    `constant` = "Additive",
    `het` = "Heterogeneous",
    `obs_shift_small` = "Observed shift small",
    `obs_shift_large` = "Observed shift large",
    `unobs_shift_small` = "Unobserved shift small",
    `unobs_shift_large` = "Unobserved shift large"
  ))) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    panel.grid.major = element_line(color = "grey95"),  # Adjust grid lines
    panel.background = element_rect(fill = "white"),  # Adjust panel background
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "bottom",  # Adjust legend position
    legend.title = element_blank(),  # Adjust legend title size
    legend.background = element_rect(color = NA, fill = NA),
    legend.text = element_text(size = 8),  # Adjust legend text size
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Title formatting
  ) +
  labs(
    x = "Estimator", 
    y = "Standard error estimates", 
    title = "Standard Error Estimates by Scenario"
  ) +
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Custom colors for points
  scale_fill_manual(values = custom_colors, labels = custom_labels) + 
  geom_hline(yintercept = oracle_se_het[[1]], linetype = "dashed", color = "black") 

# Display the plot
print(p)

db$exportOutput(p, "GLM_perf_dif_scen", Format = "pdf", FgHeight = 14, FgWidth = 25)
db$exportOutput(p, "GLM_perf_dif_scen", Format = "jpeg", FgHeight = 14, FgWidth = 25)
# -------------------------------------------------------------------------------------------
# Plot n_both change
# -------------------------------------------------------------------------------------------

# Custom labels and colors
custom_labels <- c(
  "unadjusted" = "Unadjusted",
  "GLM" = "GLM",
  "GLM_w/rand_prog" = "GLM with non-informative prognostic score",
  "GLM_w/prog" = "GLM with Super Learner prognostic score",
  "GLM_w/oracle" = "GLM with oracle prognostic score"
)

custom_colors <- c(
  "unadjusted" = "#FD8D3C",
  "GLM" = "#DD3497",
  "GLM_w/rand_prog" = "Light Sky Blue",
  "GLM_w/prog" = "#4292C6",
  "GLM_w/oracle" = "#31A354"
)


tab <- data_vary_both %>% 
  mutate(estimator = case_when(
    estr == "unadjusted" ~ "unadjusted",
    estr == "glm" & prog == "none" ~ "GLM",
    estr == "glm" & prog == "fit random" ~ "GLM_w/rand_prog",
    estr == "glm" & prog == "fit" ~ "GLM_w/prog",
    estr == "glm" & prog == "oracle" ~ "GLM_w/oracle"
  )) %>% 
  group_by(n_trial, estimator, estr, prog) %>% 
  res()

# Ensure the correct ordering of the factor levels
tab$estimator <- factor(tab$estimator, levels = c(
  'unadjusted',
  'GLM',
  'GLM_w/rand_prog',
  'GLM_w/prog',
  'GLM_w/oracle'  # Ensure this is the last level
))


# Add vertical lines for sample size at 90% power
ss_het <- ss_het %>% 
  mutate(estimator = case_when(
    estr == "unadjusted" & prog == "none" ~ "unadjusted",
    estr == "glm" & prog == "none" ~ "GLM",
    estr == "glm" & prog == "fit random" ~ "GLM_w/rand_prog",
    estr == "glm" & prog == "fit" ~ "GLM_w/prog",
    estr == "glm" & prog == "oracle" ~ "GLM_w/oracle"
  ))

# Ceiling the est_eff values
ss_het$est_eff <- ceiling(ss_het$est_eff)

# Adjust the positions for the annotations
ss_het$y_pos <- ifelse(ss_het$estimator == "GLM_w/prog", 0.25, 
                       ifelse(ss_het$estimator == "GLM_w/oracle", 0.22, 0.22))

ss_het$x_pos <- ss_het$est_eff - 2  # Move to the left

# -------------------------------------------------------------------------------------------
# Power Plot
# -------------------------------------------------------------------------------------------

p.pwr <- tab %>%
  ggplot(aes(x = n_trial, y = power, color = estimator)) +
  geom_line() +  # Removed linetype aesthetic
  xlab("n") +
  ylab("Empirically estimated power") +
  labs(color = "") +
  theme_minimal() +
  scale_y_continuous(breaks = c(.2, .4, .6, .8, 1),
                     labels = function(x){paste0(x*100, "%")}) +
  coord_cartesian(xlim = c(50, 300), ylim = c(.2, 1)) +
  scale_color_manual(values = custom_colors, labels = custom_labels) +
  theme(
    legend.position = "bottom",  # Place legend at the bottom
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black") +  # Add horizontal dashed line at 90% power
  geom_segment(data = ss_het, aes(x = est_eff, xend = est_eff, y = 0, yend = 0.9), 
               linetype = "solid", color = custom_colors[as.character(ss_het$estimator)], show.legend = FALSE) +
  geom_text(data = ss_het, aes(x = x_pos, y = y_pos, label = est_eff), 
            color = custom_colors[as.character(ss_het$estimator)], vjust = 1.5, hjust = 1, size = 3, show.legend = FALSE)

# -------------------------------------------------------------------------------------------
# Coverage Plot
# -------------------------------------------------------------------------------------------

p.cov <- tab %>%
  ggplot(aes(x = n_trial, y = coverage, color = estimator)) +
  geom_line() +
  xlab("n") +
  ylab("Empirically estimated coverage") +
  labs(color = "") +
  theme_minimal() +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_y_continuous(breaks = c(.7, .75, .8, .85, .9, .95, 1), labels = function(x){paste0(x*100, "%")}) +
  scale_color_manual(values = custom_colors, labels = custom_labels) +
  theme(
    legend.position = "bottom",  # Place legend at the bottom
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

# -------------------------------------------------------------------------------------------
# Combined Power and Coverage Plot
# -------------------------------------------------------------------------------------------

# Extract the legend from one of the plots
legend <- get_legend(p.pwr)

# Combine the plots without legends
plot_row <- plot_grid(
  p.pwr + theme(legend.position = "none"), 
  p.cov + theme(legend.position = "none"),
  labels = c('A', 'B'), label_size = 20, nrow = 1,
  label_x = c(0, 0.01), label_y = c(1.05, 1.05)
)

# Add the title
title <- ggdraw() + 
  draw_label(
    "Empirical power and coverage with increasing sample size",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )

# Combine title, plot row, and legend
combined_plot <- plot_grid(
  title, plot_row, legend,
  ncol = 1,
  rel_heights = c(0.2, 1, 0.1)
)

# Display the combined plot
print(combined_plot)



db$exportOutput(combined_plot, "GLM_power_cov", Format = "pdf", FgHeight = 18, FgWidth = 29)
db$exportOutput(combined_plot, "GLM_power_cov", Format = "jpeg", FgHeight = 18, FgWidth = 29)















