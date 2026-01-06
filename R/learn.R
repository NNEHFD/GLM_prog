# functions to learn nuisance functions needed in estimation of representable estimands

library(tidyverse)
library(magrittr)
library(zeallot)
library(tidymodels)
library(workflowsets)

lgl_fct = function(A) factor(A, levels=c(0,1))

cond_exp = function(model, df, ...){
  if (model$fit$actions$model$spec$mode == "classification") {
    predict(model, df, type='prob')$.pred_1
  } else {
    predict(model, df)$.pred
  }
}

train_wfset = function(wfset, df, v = 3) {
  format = \(df) df
  metrics = metric_set(yardstick::rmse)
  if (wfset_mode(wfset) == "classification") {   # this handles the case where outcome is binary
    outcome = wfset_outcome(wfset) |> rlang::sym() 
    format = \(df) mutate(df, {{ outcome }} := lgl_fct( {{ outcome }} ))
    metrics = metric_set(yardstick::mn_log_loss)
  }

  results = wfset |>
    workflow_map(
      fn        = "tune_grid",
      resamples = df |> format() |> vfold_cv(v = v),
      metrics   = metrics,
      control   = control_grid(save_pred = TRUE, save_workflow = T)
    )
  model_fit = fit_best(results)
  pred_fn = \(df) cond_exp(model_fit, format(df))
  list(results, pred_fn)
}

train <- function(x, ...) UseMethod("train")
train.ATE <- function(estimand, df, mu_wf, pi_wf) {

  c(mu_results, pred_mu) %<-% train_wfset(mu_wf, df)
  c(pi_results, pred_pi) %<-% train_wfset(pi_wf, df)

  list(
    list(
      pi = pred_pi,
      alpha = \(df) tibble(pi=pred_pi(df)) |> estimand$make_alpha() |> exec(df),
      m_alpha = \(df) tibble(pi=pred_pi(df)) |> estimand$make_alpha() |> estimand$m() |> exec(df),
      mu = pred_mu,
      m_mu = estimand$m(pred_mu)
    ),
    list(pi = pi_results, mu = mu_results)
  )
}
train.TSM1 = train.ATE

make_wfs <- function(specs, rec, tailor = NULL) {
  # 1) build the workflow_set from the models
  wfs <- workflow_set(
    preproc = list(base = rec),
    models  = purrr::map(specs, "model")
  )

  # 2) attach per-workflow grids
  for (nm in names(specs)) {
    grid_i <- specs[[nm]]$grid
    if (!is.null(grid_i)) {
      wfs <- option_add(
        wfs,
        grid = grid_i,
        id   = paste0("base_", nm)
      )
    }
  }
  wfs
}

wfset_mode = function(wfset) {
  wfset %>%
    mutate(
      wf   = map(wflow_id, ~ workflowsets::extract_workflow(wfset, id = .x)),
      spec = map(wf, workflows::extract_spec_parsnip),
      mode   = map_chr(spec, ~ purrr::pluck(.x, "mode",   .default = "unknown"))
    ) %>%
    pull(mode) |>
    unique()
}

wf_outcome = function(wf) {
  pre <- workflows::extract_preprocessor(wf)

  # recipe case
  if (inherits(pre, "recipe")) {
    return(
      summary(pre) |>
        dplyr::filter(role == "outcome") |>
        dplyr::pull(variable)
    )
  }

  # plain formula case
  if (inherits(pre, "formula")) {
    return(all.vars(pre[[2]]))  # LHS
  }
}

wfset_outcome = function(wfset) {
  wfset %>%
  mutate(
    wf = map(wflow_id, ~ workflowsets::extract_workflow(wfset, id = .x)),
    outcome = map(wf, wf_outcome),
    outcome = map_chr(outcome, ~ paste(.x, collapse = ", "))
  ) %>%
  pull(outcome) |>
  unique()
}

add_tailor_wfs <- function(wfs, tailor) {
  new_info <- wfs$info

  for (i in seq_len(nrow(wfs))) {
    id_i <- wfs$wflow_id[[i]]
    wf_i <- workflowsets::extract_workflow(wfs, id = id_i) |>
      workflows::add_tailor(tailor)
    new_info[[i]]$workflow <- list(wf_i)
  }
  wfs$info <- new_info
  wfs
}

adjust_probability_range <- function(x, lower = 0.01, upper = 0.99) {
  stopifnot(lower >= 0, upper <= 1, lower < upper)

  tailor::adjust_predictions_custom(
    x,
    # 1) tidymodels-style: .pred_<class>, but not .pred_class
    dplyr::across(
      dplyr::starts_with(".pred_") & !dplyr::any_of(".pred_class"),
      ~ pmin(pmax(.x, lower), upper)
    ),
    # 2) example-data-style: Class1, Class2, ...
    dplyr::across(
      dplyr::matches("^Class[0-9]+$"),
      ~ pmin(pmax(.x, lower), upper)
    ),
    .pkgs = c("dplyr")
  )
}

collect_all_metrics <- function(wfs_result, summarize = TRUE) {
  purrr::map2_dfr(
    wfs_result$wflow_id,
    wfs_result$result,
    ~ collect_metrics(.y, summarize = summarize) |>
      mutate(wflow_id = .x),
    .id = NULL
  ) |>
  arrange(mean)
}