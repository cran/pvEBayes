test_that("correctness", {
  #' @srrstats {G5.4, G5.4a, G5.5} Correctness tests are provided here
  #' @srrstats {G5.6, G5.6a, G5.6b, BS7.4, BS7.4a}
  #' recovery tests are provided here
  #' @srrstats {G5.9, G5.9a, G5.9b} Noise susceptibility test is provided here
  set.seed(1)
  ref_table <- statin42[-c(2:42), ] # limit the dimension to save time

  sig_mat <- matrix(1, nrow(ref_table), ncol(ref_table))
  sig_mat[1, 1] <- 5 # set signal strength

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
    fit_gps <- pvEBayes(
      contin_table = simu_table,
      model = "GPS"
    )
    fit_km <- pvEBayes(
      contin_table = simu_table,
      model = "KM"
    )
    fit_efron <- pvEBayes_tune(
      contin_table = simu_table,
      model = "efron"
    )
  } %>%
    suppressMessages()

  # the number of observation is nrow(ref_table) * ncol(ref_table) = 14
  # the true signal strength distribution is
  # lambda = 3 with probability 1/14 and lambda = 1 with probability 13/14
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

  # Do the samething for GPS, KM, efron

  estimated_prior_draws_gps <- post_draw_gmix_cpp(
    fit_gps$r,
    1 / fit_gps$h,
    fit_gps$omega,
    matrix(0, 1, 1), matrix(0, 1, 1), 1000
  ) %>% c()

  estimated_prior_draws_gps <- post_draw_gmix_cpp(
    fit_gps$r,
    1 / fit_gps$h,
    fit_gps$omega,
    matrix(0, 1, 1), matrix(0, 1, 1), 1000
  ) %>% c()

  estimated_prior_draws_km <- fit_km %>%
    {
      sample(.$grid,
        size = 1000, replace = TRUE,
        prob = .$g
      )
    }

  estimated_prior_draws_efron <- fit_efron %>%
    {
      sample(.$grid,
        size = 1000, replace = TRUE,
        prob = .$g
      )
    }

  # evaluate the scaled Wasserstein-2 distance between the estimated prior and
  # the true distribution via equal size draws from the estimated prior
  # distribution and the true distribution
  # (lambda = 1 with probability 13/14; lambda = 3 with probability 1/14) for
  # testing the prior recovery.

  prior_error_general_gamma <- estimated_prior_draws %>%
    sort() %>%
    {
      (. - c(rep(1, 929), rep(5, 71)))^2
    } %>%
    mean()

  prior_error_gps <- estimated_prior_draws_gps %>%
    sort() %>%
    {
      (. - c(rep(1, 929), rep(5, 71)))^2
    } %>%
    mean()

  prior_error_km <- estimated_prior_draws_km %>%
    sort() %>%
    {
      (. - c(rep(1, 929), rep(5, 71)))^2
    } %>%
    mean()

  prior_error_efron <- estimated_prior_draws_efron %>%
    sort() %>%
    {
      (. - c(rep(1, 929), rep(5, 71)))^2
    } %>%
    mean()

  expect_true(prior_error_general_gamma <= 5e-1)
  expect_true(prior_error_gps <= 5e-1)
  expect_true(prior_error_km <= 5e-1)
  expect_true(prior_error_efron <= 5e-1)


  # test correctness of posterior estimates
  post_squared_error_general_gamma <- fit %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean((. - 5)^2)
    }
  post_squared_error_gps <- fit_gps %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean((. - 5)^2)
    }
  post_squared_error_km <- fit_km %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean((. - 5)^2)
    }
  post_squared_error_efron <- fit_efron %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean((. - 5)^2)
    }
  expect_true(abs(post_squared_error_general_gamma - 0) <= 5e-1)
  expect_true(abs(post_squared_error_gps - 0) <= 5e-1)
  expect_true(abs(post_squared_error_km - 0) <= 5e-1)
  expect_true(abs(post_squared_error_efron - 0) <= 5e-1)




  # run the above test with different seed
  set.seed(2)
  ref_table <- statin42[-c(2:42), ] # limit the dimension to save time

  sig_mat <- matrix(1, nrow(ref_table), ncol(ref_table))
  sig_mat[1, 1] <- 5 # set signal strength

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
    fit_gps <- pvEBayes(
      contin_table = simu_table,
      model = "GPS"
    )
    fit_km <- pvEBayes(
      contin_table = simu_table,
      model = "KM"
    )
    fit_efron <- pvEBayes_tune(
      contin_table = simu_table,
      model = "efron"
    )
  } %>%
    suppressMessages()

  # the number of observation is nrow(ref_table) * ncol(ref_table) = 14
  # the true signal strength distribution is
  # lambda = 3 with probability 1/14 and lambda = 1 with probability 13/14
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

  # Do the samething for GPS, KM, efron

  estimated_prior_draws_gps <- post_draw_gmix_cpp(
    fit_gps$r,
    1 / fit_gps$h,
    fit_gps$omega,
    matrix(0, 1, 1), matrix(0, 1, 1), 1000
  ) %>% c()

  estimated_prior_draws_gps <- post_draw_gmix_cpp(
    fit_gps$r,
    1 / fit_gps$h,
    fit_gps$omega,
    matrix(0, 1, 1), matrix(0, 1, 1), 1000
  ) %>% c()

  estimated_prior_draws_km <- fit_km %>%
    {
      sample(.$grid,
        size = 1000, replace = TRUE,
        prob = .$g
      )
    }

  estimated_prior_draws_efron <- fit_efron %>%
    {
      sample(.$grid,
        size = 1000, replace = TRUE,
        prob = .$g
      )
    }

  # evaluate the scaled Wasserstein-2 distance between the estimated prior and
  # the true distribution via equal size draws from the estimated prior
  # distribution and the true distribution
  # (lambda = 1 with probability 13/14; lambda = 3 with probability 1/14) for
  # testing the prior recovery.

  prior_error_general_gamma <- estimated_prior_draws %>%
    sort() %>%
    {
      (. - c(rep(1, 929), rep(5, 71)))^2
    } %>%
    mean()

  prior_error_gps <- estimated_prior_draws_gps %>%
    sort() %>%
    {
      (. - c(rep(1, 929), rep(5, 71)))^2
    } %>%
    mean()

  prior_error_km <- estimated_prior_draws_km %>%
    sort() %>%
    {
      (. - c(rep(1, 929), rep(5, 71)))^2
    } %>%
    mean()

  prior_error_efron <- estimated_prior_draws_efron %>%
    sort() %>%
    {
      (. - c(rep(1, 929), rep(5, 71)))^2
    } %>%
    mean()

  expect_true(prior_error_general_gamma <= 5e-1)
  expect_true(prior_error_gps <= 5e-1)
  expect_true(prior_error_km <= 5e-1)
  expect_true(prior_error_efron <= 5e-1)


  # test correctness of posterior estimates
  post_squared_error_general_gamma <- fit %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean((. - 5)^2)
    }
  post_squared_error_gps <- fit_gps %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean((. - 5)^2)
    }
  post_squared_error_km <- fit_km %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean((. - 5)^2)
    }
  post_squared_error_efron <- fit_efron %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean((. - 5)^2)
    }
  expect_true(abs(post_squared_error_general_gamma - 0) <= 5e-1)
  expect_true(abs(post_squared_error_gps - 0) <= 5e-1)
  expect_true(abs(post_squared_error_km - 0) <= 5e-1)
  expect_true(abs(post_squared_error_efron - 0) <= 5e-1)


  # test if a trivial noise added to signal strength lead to meaningfully change
  # in results
  set.seed(2)
  ref_table <- statin42[-c(2:42), ] # limit the dimension to save time

  sig_mat <- matrix(1, nrow(ref_table), ncol(ref_table))
  sig_mat[1, 1] <- 5 + runif(1, 1e-8, 2e-8) # set signal strength

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
      mean(. - 5)^2
    }
  expect_true(abs(post_squared_error_general_gamma -
    post_squared_error_signal_noise) <=
    1e-2)


  # run the test with different reference table
  set.seed(2)
  ref_table <- statin2025_44[-c(3:42), ] # limit the dimension to save time

  sig_mat <- matrix(1, nrow(ref_table), ncol(ref_table))
  sig_mat[1, 1] <- 5 # set signal strength

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

  prior_error_general_gamma <- estimated_prior_draws %>%
    sort() %>%
    {
      (. - c(rep(1, 952), rep(5, 48)))^2
    } %>%
    mean()

  expect_true(prior_error_general_gamma <= 5e-1)

  post_squared_error_signal <- fit %>%
    summary(return = "posterior draws") %>%
    .[, 1, 1] %>%
    {
      mean(. - 5)^2
    }
  expect_true(abs(post_squared_error_signal - 0) <= 1e-2)

  ################################################################################
})
