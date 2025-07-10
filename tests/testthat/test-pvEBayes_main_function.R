

test_that("pvEBayes", {
  valid_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2)
  rownames(valid_matrix) <- c("AE_1", "AE_2")
  colnames(valid_matrix) <- c("drug_1", "drug_2", "drug_3", "drug_4")
  expect_equal(.is_valid_contin_table(valid_matrix), TRUE)
  result <- calculate_tilde_e(valid_matrix)

  expect_equal(rownames(result), rownames(valid_matrix))
  expect_equal(colnames(result), colnames(valid_matrix))

  # pseudo sample generation
  grid <- .grid_based_on_hist_log_scale_sobol(valid_matrix,
                                              result, max_draws = 200)
  expect_equal(length(grid), 200)

  # check the main function
  ## GPS
  fit_gps <- pvEBayes(contin_table = valid_matrix, model = "GPS")
  expect_equal(is.pvEBayes(fit_gps), TRUE)

  ## K-gamma
  fit_4g <- pvEBayes(contin_table = valid_matrix, model = "K-gamma", K = 4)
  expect_equal(is.pvEBayes(fit_4g), TRUE)

  ## general-gamma
  fit_gg <- pvEBayes(
    contin_table = valid_matrix,
    model = "general-gamma", alpha = 0.5
  )
  expect_equal(is.pvEBayes(fit_gg), TRUE)


  ## efron
  fit_e <- pvEBayes(
    contin_table = valid_matrix,
    model = "efron", p = 40, c0 = 0.01
  )
  expect_equal(is.pvEBayes(fit_e), TRUE)

  # check hyperparameter alpha selection
  gg_selection <- pvEBayes_tune(valid_matrix, model = "general-gamma")
  expect_equal(is.pvEBayes(gg_selection), TRUE )

  # check hyperparameter alpha selection
  e_selection <- pvEBayes_tune(valid_matrix, model = "efron")
  expect_equal(is.pvEBayes(e_selection), TRUE )

  #
} )

