test_that("pvEBayes", {
  #' @srrstats {BS5.3, BS5.4, BS5.5}
  #' The empirical Bayes methods implemented in \pkg{pvEBayes} do not rely on
  #' stochastic sampling, and therefore do not produce the types of
  #' convergence diagnostics typically associated with full Bayesian modeling.
  #' Convergence in the ECM algorithm is reached (at least to a sub-optimal).
  #' This is ensured by monotonically increased log joint marginal likelihood,
  #' as proved by Tan et al. (*Stat. in Med*, 2025). In the implementation
  #' the convergency of ECM algorithm is checked by comparing the absolute
  #' difference in log likelihood between two consecutive iterations to a
  #' tolerance argument "tol_ecm".
  #' @srrstats {G5.8, G5.8a, G5.8b, G5.8c, G5.8d}
  #' Edge condition tests are provided below.

  valid_matrix <- matrix(c(70, 40, 44, 50, 52, 60, 70, 80), nrow = 2)
  rownames(valid_matrix) <- c("AE_1", "AE_2")
  colnames(valid_matrix) <- c("drug_1", "drug_2", "drug_3", "drug_4")
  expect_equal(.is_valid_contin_table(valid_matrix), TRUE)
  result <- estimate_null_expected_count(valid_matrix)

  expect_equal(rownames(result), rownames(valid_matrix))
  expect_equal(colnames(result), colnames(valid_matrix))
  E <- estimate_null_expected_count(valid_matrix)

  # test invalid input
  invalid_matrix <- matrix(c(-70, 40, 44, 50, 52, 60, 70, 80), nrow = 2)
  invalid_E <- E
  invalid_E[1, 1] <- -10

  expect_warning(
    .is_valid_contin_table(contin_table = invalid_matrix)
  )
  expect_error(
    pvEBayes(contin_table = integer(0), model = "GPS")
  )
  expect_warning(expect_error(
    pvEBayes(contin_table = matrix(NA, 2, 2), model = "GPS")
  ))

  ## valid contin_table + invalid E
  expect_error(
    pvEBayes(contin_table = valid_matrix, model = "GPS", E = invalid_E)
  )


  ## general-gamma without specifying alpha
  expect_error(
    pvEBayes(contin_table = valid_matrix, model = "general-gamma")
  )

  ## K-gamma without specifying K
  expect_error(
    pvEBayes(contin_table = valid_matrix, model = "K-gamma")
  )

  expect_error(
    pvEBayes(contin_table = valid_matrix, model = "sdfa")
  )

  expect_message(
    pvEBayes(contin_table = valid_matrix, model = "GPS", alpha = 0.5)
  )

  expect_message(
    pvEBayes(
      contin_table = valid_matrix, model = "general-gamma",
      alpha = 0.5, K = 4
    )
  )




  ## check if errors can be caught

  expect_true(
    tryCatch(
      {
        pvEBayes(contin_table = valid_matrix, model = "GPS", E = invalid_E)
        FALSE
      },
      error = function(e) TRUE
    )
  )

  # pseudo sample generation
  grid <- .grid_based_on_hist_log_scale_sobol(valid_matrix,
    result,
    max_draws = 200
  )
  expect_equal(length(grid), 200)

  # check the main function
  ## GPS
  fit_gps <- pvEBayes(contin_table = valid_matrix, model = "GPS")
  expect_equal(is.pvEBayes(fit_gps), TRUE)
  qn_tmp <- .calc_qn_mix_gamma(fit_gps, valid_matrix, E, log = TRUE)
  print_tmp <- print(fit_gps)
  summary_tmp <- summary(fit_gps)
  return_tmp <- summary(fit_gps, return = "prior parameters")
  return_tmp <- summary(fit_gps, return = "likelihood")
  return_tmp <- summary(fit_gps, return = "detected signal")
  return_tmp <- summary(fit_gps, return = "posterior draws")
  ## K-gamma
  fit_4g <- pvEBayes(contin_table = valid_matrix, model = "K-gamma", K = 4)
  expect_equal(is.pvEBayes(fit_4g), TRUE)

  ## general-gamma
  fit_gg <- pvEBayes(
    contin_table = valid_matrix,
    model = "general-gamma", alpha = 0.5
  )
  expect_equal(is.pvEBayes(fit_gg), TRUE)

  ## check data.table input

  fit_gg_dt <- pvEBayes(
    contin_table = data.table::data.table(valid_matrix),
    model = "general-gamma", alpha = 0.5
  )
  expect_equal(is.pvEBayes(fit_gg_dt), TRUE)

  ## KM
  grid <- .grid_based_on_hist_log_scale_sobol(valid_matrix,
    E,
    max_draws = FALSE
  )

  fit_km <- pvEBayes(
    contin_table = valid_matrix,
    model = "KM"
  )
  print_tmp <- print(fit_km)
  summary_tmp <- summary(fit_km)
  expect_equal(is.pvEBayes(fit_km), TRUE)


  ## efron
  fit_e <- pvEBayes(
    contin_table = valid_matrix,
    model = "efron", p = 40, c0 = 0.01
  )
  print_tmp <- print(fit_e)
  summary_tmp <- summary(fit_e)
  expect_equal(is.pvEBayes(fit_e), TRUE)

  # check hyperparameter alpha selection
  gg_selection <- pvEBayes_tune(valid_matrix, model = "general-gamma")
  expect_equal(is.pvEBayes(gg_selection), TRUE)

  # check if the message can be suppressed
  expect_no_message(
    suppressMessages(
      gg_selection <- pvEBayes_tune(valid_matrix, model = "general-gamma")
    )
  )

  # check hyperparameter alpha selection
  e_selection <- pvEBayes_tune(valid_matrix, model = "efron")
  expect_equal(is.pvEBayes(e_selection), TRUE)

  all_fit_e <- extract_all_fitted_models(e_selection)
  # check for NA, NaN,Inf for log likelihood.
  check_failure_fit <- function(loglik) {
    res <- is.finite(loglik)
    res
  }

  expect_true(is.finite(fit_gps$loglik))
  expect_true(is.finite(fit_e$loglik))
  expect_true(is.finite(fit_km$loglik))
  expect_true(is.finite(fit_gg$loglik))
  expect_true(is.finite(fit_4g$loglik))
})
