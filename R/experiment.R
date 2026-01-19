library(furrr)
future::plan(future::multisession, workers = future::availableCores())
`%<-%` = zeallot::`%<-%`

library(FastKRR)
library(parsnip)
library(zeallot)
library(postcard)
library(MASS)
library(dplyr)
library(magrittr)
source('R/learn.R')

MU_WF = list(
    # gbt = list(
    #   model = boost_tree(trees=tune(), tree_depth=3, learn_rate=0.05) |>
    #     set_engine("xgboost", nthread = 1) |>
    #     set_mode("regression"),
    #   grid = expand_grid(
    #     trees = c(1:20 * 10, 201)
    #   )
    # ),
    # glm = list(
    #   model = linear_reg(penalty = tune(), mixture = tune()) |> 
    #     set_engine("glmnet") |>
    #     set_mode("regression"),
    #   grid = expand_grid(
    #     penalty = c(0.001, 0.1),
    #     mixture = c(0, 1)
    #   )
    # ),
    # krr = list(
    #   model = krr_reg(kernel="gaussian", penalty=tune(), rho=0.1) |>
    #     set_engine("fastkrr", verbose=F) |>
    #     set_mode("regression"),
    #   grid = expand_grid(
    #     penalty = c(0.001, 0.1)
    #   )
    # ),
    mars = list(
      model = mars(prod_degree = 3, num_terms=50) |> 
        set_engine("earth") |>
        set_mode("regression"),
      grid = NULL
    )
  ) |> 
  make_wfs(Y ~ .)

fit_prog = function(hist, n_hist, mu_wf=MU_WF) {
  data_hist_oracle = hist %>%
    draw(n_hist) 
  
  data_hist = data_hist_oracle %>%
    mutate(Y = Y0) %>%
    dplyr::select(Y, starts_with('W'))

  suppressMessages(suppressWarnings({
      c(lrnr_results, lrnr) %<-% train_wfset(mu_wf, data_hist)
  }))

  power_stats = list(
    psi0 = mean(data_hist$Y),
    sigma2 = var(data_hist$Y),
    kappa_oracle = data_hist_oracle %$% sqrt(mean((Y0 - mu0)^2)),
    kappa_ml = collect_all_metrics(lrnr_results) |>
      slice_min(mean, n=1) |>
      pull(mean),
    kappa_glm = glm(Y ~ ., data = data_hist, family = MASS::negative.binomial(theta = 3, link = "log")) |>
      residuals(type = "response") |>
      raise_to_power(2) |>
      mean() |>
      sqrt()
  )

  list(lrnr_results, power_stats, lrnr)
}

estimate_glm = function(df, cv_variance=T) {
  df |>
  postcard::rctglm(formula = Y ~ .,
        exposure_indicator = A,
        exposure_prob = 1/2, 
        family = MASS::negative.binomial(3, link = "log"), # theta=3 is an arbitrary choice, result is insensitive
        estimand_fun = "rate_ratio",
        verbose = 0,
        cv_variance = cv_variance,
        cv_variance_folds = 10) %>% 
  postcard::estimand()
}

analyze_trial = function(trial, n_trial, lrnr, cv=T) {
  data_trial = trial |>
    draw(n_trial) %>%
    dplyr::select(Y, A, starts_with('W'), mu0)
  
  processed_data = data_trial %>% 
    mutate(
      prog = data_trial |> 
        dplyr::select(starts_with('W')) |>
        lrnr() |>
        link(),
      mu0 = link(mu0)
    )
    
  bind_rows(
      #unadjusted  
      dplyr::select(processed_data, Y, A) |>
        estimate_glm(cv_variance = F) |>
        mutate(estimator = 'unadjusted'),
      #No prognostic score but glm adjusted with W 
      dplyr::select(processed_data, Y, A, starts_with('W')) |>
        estimate_glm(cv_variance=cv) |>
        mutate(estimator="GLM"),
      # Random/non-informative prognostic score  
      dplyr::select(processed_data, Y, A, starts_with('W'), prog) |>
        mutate(prog = sample(prog)) |> # shuffle
        estimate_glm(cv_variance=cv) |>
        mutate(estimator = "noise-adj. GLM"),
      #Prognostic score
      dplyr::select(processed_data, Y, A, starts_with('W'), prog) |>
        estimate_glm(cv_variance=cv) |>
        mutate(estimator = "prog-adj. GLM"),
      #Prognostic score (no covar)
      dplyr::select(processed_data, Y, A, prog) |>
        estimate_glm(cv_variance=cv) |>
        mutate(estimator = "prog-adj. GLM (no covar.)"),
      # Oracle prognostic score 
      dplyr::select(processed_data, Y, A, starts_with('W'), mu0) |>
        estimate_glm(cv_variance=cv) |>
        mutate(estimator = "oracle-adj. GLM")
    )
}

experiment = function(hist, trial, n_hist, n_trial, mu_wf = MU_WF, cv=T) {
  
  c(lrnr_results, power_stats, lrnr) %<-% fit_prog(hist, n_hist, mu_wf)
  estimates = analyze_trial(trial, n_trial, lrnr, cv=cv) |>
    mutate(
      n_hist = n_hist,
      n_trial = n_trial,
      hist = list(hist),
      trial = list(trial)
    )
  
  list(estimates, collect_all_metrics(lrnr_results))
}

rep_experiment = function(
  hist, 
  trial, 
  n_hist = 5000,
  n_trial = 250,
  reps = 100,
  cv = T
) {
  1:reps |> furrr::future_map(
    ~{
      draw.dgp # Force export of S3 method
      experiment(hist, trial, n_hist, n_trial, cv=cv)
    }, 
    .options = furrr::furrr_options(seed = 3021377)
  ) |> 
    list_transpose() |>
    map(bind_rows, .id = 'rep')
}

report = function(estimates) {

  estimates |>
    rename(estimate=Estimate, est_se = `Std. Error`) |>
    mutate(ci_lower= estimate - 2*est_se , ci_upper= estimate + 2*est_se) |>
    group_by(estimator, scenario, n_hist, n_trial) |>
    summarise(
      bias = mean(estimate - true_rr),
      se = sd(estimate),
      mean_est_se = mean(est_se),
      rmse = sqrt(mean((estimate - true_rr)^2)),
      coverage = mean(ci_lower <= true_rr & true_rr <= ci_upper),
      pct_null = mean(ci_lower <= 1 & 1 <= ci_upper),
      pct_significant = 1 - pct_null,
      .groups = "drop"
    ) |>
    arrange(rmse)
}