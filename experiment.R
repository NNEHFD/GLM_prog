# -------------------------------------------------------------------------------------------
# Purpose: Function for running the experiment 
#
# Dependencies: 
#
# Output: 
#
# -------------------------------------------------------------------------------------------
# default_learners list defined by the model, grid for tuning if needed
default_learners = list( 
  mars = list(
    model = mars(mode = "regression", prod_degree = 3) %>%
      set_engine("earth"),
    grid = NULL
  ),
  
  lm_pois = list(
    model = linear_reg() %>%
      set_engine("glm", family = "poisson"),
    grid = NULL
  ),
  
  lm_neg = list(
    model = linear_reg() %>%
      set_engine("glm", family = negative.binomial(3, link = "log")),
    grid = NULL
  ),
  
  gbt = list(
    model = boost_tree(
      mode = "regression",
      trees = tune("trees"),
      tree_depth = tune("tree_depth"),
      learn_rate = 0.1
    ) %>%
      set_engine("xgboost"),
    grid = expand_grid(trees = seq.int(25, 500, by = 25),
                       tree_depth = c(3)
    )
  )
)

experiment = function(n_hist, n_trial,
                      p=p,
                      shift_W1=shift_W1, # historical only
                      shift_U=shift_U, # historical only
                      dgp_string =c("constant")) {
  
  outcome_dgp(dgp_string)
  data_hist = dgp(n = n_hist,
                  p = p, 
                  shift_W1 = shift_W1,
                  shift_U = shift_U) %>% 
    mutate(Y = Y0) %>%
    dplyr::select(Y, starts_with('W'))
  
  if(n_hist >= 5000) {V = 3} else if(n_hist >= 1000) {V = 5} else {V = 10}
  
  lrnr = PostCard::fit_best_learner(data = data_hist %>% 
                                    dplyr::select(Y, starts_with('W')), 
                                    formula = Y ~ .,
                                    cv_folds = V,
                                    learners = default_learners,
                                    verbose = 0)
  
  data_trial = dgp(n = n_trial,
                   p = p, 
                   shift_W1 = 0,
                   shift_U = 0) %>%
    dplyr::select(Y, A, starts_with('W'), m0) #m0 is the oracle prognostic score 
  
  data_trial %<>% mutate(prog = predict(lrnr, data_trial %>% 
                                          dplyr::select(Y, starts_with('W'))) %>% 
                           pull(.pred) %>% 
                           pmax(1e-10, .))  # add the estimated prognostic score to the data
  
  data_trial %<>% mutate(m0 = link(m0), prog = link(prog)) # add the true prognostic score to the data

  
  bind_rows(
    #unadjusted 
    data_trial %>% 
      dplyr::select(-prog, -m0) %>%
      rctglm(formula = Y ~ A,
             group_indicator = A,
             family = "gaussian",
             estimand_fun = "rate_ratio",
             verbose = 0) %>% 
      estimand() %>% 
      mutate(prog="none", estr="unadjusted"),
    #No prognostic score but glm adjusted with W
    data_trial %>% 
      dplyr::select(-prog, -m0) %>%
      rctglm(formula = Y ~ .,
             group_indicator = A,
             family = negative.binomial(3, link = "log"),
             estimand_fun = "rate_ratio",
             verbose = 0) %>% 
      estimand() %>% 
      mutate(prog="none", estr="glm"),
    # Random/non-informative prognostic score 
    data_trial %>% 
      dplyr::select(-prog, -m0) %>%
      mutate(prog = runif(nrow(data_trial), min = min(data_hist$Y), max = max(data_hist$Y))) %>%
      rctglm(formula = Y ~ .,
             group_indicator = A,
             family = negative.binomial(3, link = "log"),
             estimand_fun = "rate_ratio",
             verbose = 0) %>% 
      estimand() %>% 
      mutate(prog="fit random", estr="glm"),
    #Prognostic score
    data_trial %>% 
      dplyr::select(-m0) %>%
      rctglm(formula = Y ~ .,
             group_indicator = A,
             family = negative.binomial(3, link = "log"),
             estimand_fun = "rate_ratio",
             verbose = 0) %>% 
      estimand() %>% 
      mutate(prog="fit", estr="glm"),
    # Oracle prognostic score
    data_trial %>% 
      dplyr::select(-prog) %>%
      rctglm(formula = Y ~ .,
             group_indicator = A,
             family = negative.binomial(3, link = "log"),
             estimand_fun = "rate_ratio",
             verbose = 0) %>% 
      estimand() %>% 
      mutate(prog="oracle", estr="glm")
) %>%
    mutate(
      n_hist = n_hist,
      n_trial = n_trial,
      p = p,
      shift_W1 = shift_W1,
      shift_U = shift_U
    )
}
