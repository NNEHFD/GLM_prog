
# Common definitions for plotting
library(ggtext)
library(glue)

# Colors
plot_colors = c(
  "Prog. + Covar." = "#02d769ff",
  "Prog." = "#0b9099ff",
  "Oracle + Covar." = "#2137c6ff",
  "Noise + Covar." = "#6f3637ff",
  "Covar." = "#333333ff",
  "None" = "#797979ff"
)

# Estimator formatting
format_estimator <- function(x) {
  factor(x, levels = c(
    "prog-adj. GLM", 
    "prog-adj. GLM (no covar.)", 
    "oracle-adj. GLM", 
    "noise-adj. GLM", 
    "GLM", 
    "unadjusted"
  ), labels = c(
    "Prog. + Covar.", 
    "Prog.", 
    "Oracle + Covar.", 
    "Noise + Covar.", 
    "Covar.", 
    "None"
  ))
}

# Scenario formatting
format_scenario_data <- function(data) {
  data |>
    tidyr::separate(scenario, into = c("trial_type", "hist_type"), sep = "-", remove = FALSE) |>
    mutate(
      # Define the top-level grouping
      `trial DGP` = case_when(
        trial_type == "null" ~ "Null",
        trial_type == "additive" ~ "Additive",
        trial_type == "hetero" ~ "Heterogeneous"
      ) |> factor(levels = c("Null", "Additive", "Heterogeneous")),
      # Define the subgrouping (shift type)
      `modification for historical DGP` = case_when(
        hist_type == "no_shift" ~ "No Shift",
        hist_type == "shift_obs_small" ~ "Small Observed Shift",
        hist_type == "shift_unobs_small" ~ "Small Unobserved Shift",
        hist_type == "shift_obs_large" ~ "Large Observed Shift",
        hist_type == "shift_unobs_large" ~ "Large Unobserved Shift"
      ) |> factor(levels = c(
        "No Shift", 
        "Small Observed Shift", 
        "Small Unobserved Shift", 
        "Large Observed Shift", 
        "Large Unobserved Shift"
      ))
    )
}
