test_that("correctness", {
  #' @srrstats {G5.4, G5.4a, G5.5} Correctness tests are provided here
  #' @srrstats {G5.6, G5.6a, G5.6b, G5.7, BS7.0, BS7.1, BS7.2, BS7.4, BS7.4a}
  #' recovery tests are provided here
  #' @srrstats {G5.9, G5.9a, G5.9b} Noise susceptibility test is provided here
  set.seed(1)
  ref_table <- statin42[-c(10:30), ] # limit the dimension to save time

  sig_mat <- matrix(1, nrow(ref_table), ncol(ref_table))
  sig_mat[1, 1] <- 3 # set signal strength

  simu_table <- generate_contin_table(
    ref_table = ref_table,
    signal_mat = sig_mat,
    n_table = 1
  )[[1]]

  {
    fit <- pvEBayes_tune(
      contin_table = simu_table,
      model = "general-gamma"
    )
  } %>%
    suppressMessages()

  # the number of observation is nrow(ref_table) * ncol(ref_table) = 154
  # the true signal strength distribution is
  # lambda = 3 with probability 1/154 and lambda = 1 with probability 153/154
  # the general-gamma model is a nonparametric empirical Bayes model where the
  # prior distribution of signal strength is estimated via a gamma mixture model
  # with parameters (r,h, omega). see vignette for detail.
  # Here we borrow the function that take posterior draws to generate draws from
  # the estimated prior distribution

  estimated_prior_draws <- post_draw_gmix_cpp(
    fit$r,
    1 / fit$h,
    fit$omega,
    matrix(0, 1, 1), matrix(0, 1, 1), 1000
  ) %>% c()

  # evaluate the scaled square distance between the estimated prior and the true
  # distribution via draws from the estimated prior distribution for testing the
  # prior recovery

  prior_error <- (((estimated_prior_draws - 1) / 1)^2 * (153 / 154) +
    ((estimated_prior_draws - 3) / 3)^2 * (1 / 154)) %>%
    mean()

  expect_true(prior_error <= 1e-1)

  post_squared_error_signal <- fit %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean(. - 3)^2
    }
  expect_true(abs(post_squared_error_signal - 0) <= 1e-3)




  # run the above test with different seed
  set.seed(2)
  ref_table <- statin42[-c(10:30), ] # limit the dimension to save time

  sig_mat <- matrix(1, nrow(ref_table), ncol(ref_table))
  sig_mat[1, 1] <- 3 # set signal strength

  simu_table <- generate_contin_table(
    ref_table = ref_table,
    signal_mat = sig_mat,
    n_table = 1
  )[[1]]

  {
    fit <- pvEBayes_tune(
      contin_table = simu_table,
      model = "general-gamma"
    )
  } %>%
    suppressMessages()


  estimated_prior_draws <- post_draw_gmix_cpp(
    fit$r,
    1 / fit$h,
    fit$omega,
    matrix(0, 1, 1), matrix(0, 1, 1), 1000
  ) %>% c()

  prior_error <- ((estimated_prior_draws - 1)^2 * (153 / 154) +
    (estimated_prior_draws - 1)^2 * (1 / 154)) %>%
    mean()

  expect_true(prior_error <= 1e-1)

  post_squared_error_signal <- fit %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean(. - 3)^2
    }
  expect_true(abs(post_squared_error_signal - 0) <= 1e-2)


  # test if a trivial noise added to signal strength lead to meaningfully change
  # in results
  set.seed(2)
  ref_table <- statin42[-c(10:30), ] # limit the dimension to save time

  sig_mat <- matrix(1, nrow(ref_table), ncol(ref_table))
  sig_mat[1, 1] <- 3 + runif(1, 1e-8, 2e-8) # set signal strength

  simu_table <- generate_contin_table(
    ref_table = ref_table,
    signal_mat = sig_mat,
    n_table = 1
  )[[1]]

  {
    fit <- pvEBayes_tune(
      contin_table = simu_table,
      model = "general-gamma"
    )
  } %>%
    suppressMessages()





  post_squared_error_signal_noise <- fit %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean(. - 3)^2
    }
  expect_true(abs(post_squared_error_signal - post_squared_error_signal_noise) <=
    1e-2)


  # run the test with different reference table
  set.seed(2)
  ref_table <- statin2025_44[-c(10:30), ] # limit the dimension to save time

  sig_mat <- matrix(1, nrow(ref_table), ncol(ref_table))
  sig_mat[1, 1] <- 3 # set signal strength

  simu_table <- generate_contin_table(
    ref_table = ref_table,
    signal_mat = sig_mat,
    n_table = 1
  )[[1]]

  {
    fit <- pvEBayes_tune(
      contin_table = simu_table,
      model = "general-gamma"
    )
  } %>%
    suppressMessages()


  estimated_prior_draws <- post_draw_gmix_cpp(
    fit$r,
    1 / fit$h,
    fit$omega,
    matrix(0, 1, 1), matrix(0, 1, 1), 1000
  ) %>% c()

  prior_error <- ((estimated_prior_draws - 1)^2 * (153 / 154) +
    (estimated_prior_draws - 1)^2 * (1 / 154)) %>%
    mean()

  expect_true(prior_error <= 1e-1)

  post_squared_error_signal <- fit %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean(. - 3)^2
    }
  expect_true(abs(post_squared_error_signal - 0) <= 1e-2)
})
