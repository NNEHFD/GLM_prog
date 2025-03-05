# -------------------------------------------------------------------------------------------
# Purpose: Create function for running the sample size calculation 
#
# Dependencies: dgp.R
#
# Output: Function for doing the sample size estimation 
#
# -------------------------------------------------------------------------------------------

# Load libraries ------------------------------------------------------------------

source("statprog/GLM/dgp.R")

# Help functions ---------------------------------------------------------------------

learners <- list(gbt = list(
  model = boost_tree(
    mode = "regression",
    trees = tune("trees"),
    tree_depth = tune("tree_depth"),
    learn_rate = 0.1
  ) %>%
    set_engine("xgboost"),
  grid = expand_grid(trees = seq.int(25, 500, by = 25), tree_depth = c(3))
))

lrnr <- function(data) {
  n_hist <- nrow(data)
  if(n_hist >= 5000) {V = 3} else if(n_hist >= 1000) {V = 5} else {V = 10}
  
  PostCard::fit_best_learner(data = data %>% 
                               dplyr::select(Y, starts_with('W')), 
                             formula = Y ~ .,
                             cv_folds = V,
                             learners = learners,
                             verbose = 0)
}


# Experiment ---------------------------------------------------------------------

experiment_power = function(n_hist,
                      p=p,
                      shift_W1=shift_W1, # historical only
                      shift_U=shift_U, # historical only
                      dgp_string =c("constant"),
                      target_effect,
                      pi = 1/2,
                      psi_1_func,
                      dmarginaleffect_1,
                      dmarginaleffect_0,
                      r = 3, 
                      alpha = 0.05,
                      gamma = 0.8,
                      initial_n = 100,
                      increment = 2) {
  
  outcome_dgp(dgp_string)
  data_hist = dgp(n = n_hist,
                  p = p, 
                  shift_W1 = shift_W1,
                  shift_U = shift_U) %>% 
    mutate(Y = Y0) %>%
    dplyr::select(Y, starts_with('W'), m0)
  
  if(n_hist >= 5000) {V = 3} else if(n_hist >= 1000) {V = 5} else {V = 10}
  
  
  
  
  var_unadjusted <- variance_bound(
    data_hist,
    target_effect,
    pi,
    hatmu = function(data) {
      glm(Y ~ 1, family = negative.binomial(1/r, link = "log"), data = data)
    },
    V
  )
  
  
  var_glm <- variance_bound(
    data_hist,
    target_effect,
    pi,
    hatmu = function(data) {
      glm(Y ~ ., family = negative.binomial(1/r, link = "log"), data = data %>% dplyr::select(-m0, -predictions))
    },
    V
  )
  
  var_prog <- variance_bound(
    data_hist,
    target_effect,
    pi,
    hatmu = lrnr,
    V,
    pull = TRUE
  )
  
  
  var_oracle <- variance_bound(
    data_hist,
    target_effect,
    pi,
    hatmu = function(data) {
      glm(Y ~ ., family = negative.binomial(1/r, link = "log"), data = data %>% dplyr::select(-predictions))
    },
    V
  )
  
  bind_rows(
    #unadjusted 
    sample_size_glm(target_effect, var_unadjusted, alpha, gamma, initial_n, increment) %>% 
      mutate(prog="none", estr="unadjusted"),
    #No prognostic score but glm adjusted with W
    sample_size_glm(target_effect, var_glm, alpha, gamma, initial_n, increment) %>%  
      mutate(prog="none", estr="glm"),
    #Prognostic score
    sample_size_glm(target_effect, var_prog, alpha, gamma, initial_n, increment) %>%
      mutate(prog="fit", estr="glm"),
    # Oracle prognostic score
    sample_size_glm(target_effect, var_oracle, alpha, gamma, initial_n, increment) %>%
      mutate(prog="oracle", estr="glm")
  ) %>%
    mutate(
      n_hist = n_hist,
      p = p,
      shift_W1 = shift_W1,
      shift_U = shift_U
    )
}










