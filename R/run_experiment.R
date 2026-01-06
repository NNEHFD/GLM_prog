# -------------------------------------------------------------------------------------------
# Purpose: Running the experiment in different scenarios
#
# Dependencies:
#
# Output:
#
# -------------------------------------------------------------------------------------------

# loading library and dependencies --------------------------
library(future)
library(tidymodels)
library(magrittr)
library(tidyverse)
library(furrr)
library(postcard)

source("R/dgp.R")
source("R/experiment.R")

# Define output directory
output_dir <- file.path("..", "outputs", "datasets")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# setting up the parallel processing ----------------------------------------
# future::plan(multisession, workers = 90)
# furrr_options(seed = TRUE)
N <- 25
set.seed(455321)

# print("init")
# # -------------------------------------------------------------------------------------------
# # Performance in different scenarios
# # -------------------------------------------------------------------------------------------
#
# Constant treatment effect - no shift
params = expand_grid(
  n_hist = c(4000),
  n_trial = c(200),
  shift_W1 = c(0), # shift in historical data only
  shift_U = c(0), # shift in historical data only
  dgp_string = c("constant")
)
results = 1:N %>% map( ~ params %>% pmap(experiment), .options = furrr_options(seed = 3021377))
results %<>% purrr::list_flatten() %>% purrr::list_rbind()

saveRDS(results, file.path(output_dir, "data_constant_cv.rds"))
print("constant")


# Heterogeneous - no shift
params = expand_grid(
  n_hist = c(4000),
  n_trial = c(200),
  shift_W1 = c(0),
  shift_U = c(0),
  dgp_string = c("heterogeneous")
)
results = 1:N %>% map( ~ params %>% pmap(experiment), .options = furrr_options(seed = 3021377))
results %<>% purrr::list_flatten() %>% purrr::list_rbind()

saveRDS(results, file.path(output_dir, "data_het_cv.rds"))
print("het")

# Heterogeneous - observable shift in W1 (small and large)
params = expand_grid(
  n_hist = c(4000),
  n_trial = c(200),
  shift_W1 = c(-3, -5),
  shift_U = c(0),
  dgp_string = c("heterogeneous")
)
results = 1:N %>% map( ~ params %>% pmap(experiment), .options = furrr_options(seed = 3021377))
results %<>% purrr::list_flatten() %>% purrr::list_rbind()

saveRDS(results, file.path(output_dir, "data_obs_shift_cv.rds"))

print("obs_shift")

# Heterogeneous - unobservable shift in U (small and large)
params = expand_grid(
  n_hist = c(4000),
  n_trial = c(200),
  shift_W1 = c(0),
  shift_U = c(.5, 1),
  dgp_string = c("heterogeneous")
)
results = 1:N %>% map( ~ params %>% pmap(experiment), .options = furrr_options(seed = 3021377))
results %<>% purrr::list_flatten() %>% purrr::list_rbind()

saveRDS(results, file.path(output_dir, "data_unobs_shift_cv.rds"))

print("unobs_shift")

# # -------------------------------------------------------------------------------------------
# # Performance varying sample size in heterogeneous treatment effect scenario
# # -------------------------------------------------------------------------------------------
# 
# # varying historical data size
# params = expand_grid(
#   n_hist = c(250, 500, 750, 1000, 2000, 4000, 6000, 8000),
#   n_trial = c(300),
#   p = c(7),
#   shift_W1 = c(0),
#   shift_U = c(0),
#   dgp_string = c("heterogeneous")
# )
# results = 1:N %>% future_map_dfr( ~ params %>% pmap_df(experiment), .options = furrr_options(seed = 3021377))
# 
# db$output_datasets("data_vary_hist", results, ext = "rds")
# 
# print("vary_hist")
# #
# # varying trial data size
# params = expand_grid(
#   n_hist = c(4000),
#   n_trial = c(50, 100, 150, 300, 400, 500, 600, 700, 800, 1000),
#   p = c(7),
#   shift_W1 = c(0),
#   shift_U = c(0),
#   dgp_string = c("heterogeneous")
# )
# results = 1:N %>% future_map_dfr( ~ params %>% pmap_df(experiment), .options = furrr_options(seed = 3021377))
# 
# db$output_datasets("data_vary_trial", results, ext = "rds")
# 
# print("vary_trial")
# 
# varying both historical and trial data size
inc <- c(20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 225, 250, 275, 300)


for (i in 1:length(inc)) {
  params = expand_grid(
    n_hist = inc[i]*10,
    n_trial = inc[i],
    p = c(7),
    shift_W1 = c(0),
    shift_U = c(0),
    dgp_string = c("heterogeneous")
  )
  results = 1:N %>% future_map( ~ params %>% pmap(experiment), .options = furrr_options(seed = 3021377))
  results %<>% purrr::list_flatten() %>% purrr::list_rbind()

  saveRDS(results, file.path(output_dir, paste0("data_vary_both_cv_", inc[i], ".rds")))
}


print("vary_both")

#future::plan(sequential)

