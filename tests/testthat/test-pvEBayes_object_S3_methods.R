test_that("pvEBayes object associated functions", {
  valid_matrix <- matrix(c(70, 40, 44, 50, 52, 60, 70, 80), nrow = 2)

  obj <- pvEBayes(
    contin_table = valid_matrix,
    model = "general-gamma", alpha = 0.1
  )
  obj_redraw <- posterior_draws(obj)
  plot1 <- eyeplot_pvEBayes(obj, text_shift = 0.1)
  plot1 <- eyeplot_pvEBayes(obj, text_shift = 0.1, log_scale = TRUE)
  expect_equal(inherits(plot1, "ggplot"), TRUE)

  plot2 <- heatmap_pvEBayes(obj)
  expect_equal(inherits(plot2, "ggplot"), TRUE)
  print_tmp <- print(obj)
  expect_equal(is.pvEBayes(print_tmp), TRUE)

  summary_res <- summary(obj)
  expect_equal(is.list(summary_res), TRUE)

  likelihood_obj <- logLik(obj)
  expect_equal(is.numeric(likelihood_obj), TRUE)

  plot3 <- plot(obj, type = "eyeplot")
  plot4 <- plot(obj, type = "heatmap")
  expect_equal(inherits(plot3, "ggplot"), TRUE)
  expect_equal(inherits(plot4, "ggplot"), TRUE)



  # check for errors
  expect_error(
    extract_all_fitted_models(1)
  )

  expect_error(
    posterior_draws(obj, n_posterior_draws = -1)
  )

  expect_error(
    plot(obj, type = "eyeplot", num_top_AEs = -1)
  )
  expect_error(
    plot(obj, type = "eyeplot", num_top_drugs = -1)
  )
  expect_error(
    plot(obj, type = "eyeplot", N_threshold = -1)
  )
  expect_error(
    plot(obj, type = "heatmap", num_top_AEs = -1)
  )
  expect_error(
    plot(obj, type = "eyeplot", text_size = "a")
  )
  expect_error(
    plot(obj, type = "eyeplot", text_shift = "a")
  )
  expect_error(
    plot(obj, type = "eyeplot", x_lim_scalar = "a")
  )
  expect_error(
    plot(obj, type = "eyeplot", specified_AEs = 1)
  )

  expect_error(
    plot(obj, type = "eyeplot", specified_drugs = 1)
  )
  expect_error(
    plot(obj, type = "heatmap", num_top_drugs = -1)
  )

  expect_error(
    plot(obj, type = "heatmap", specified_AEs = 1)
  )

  expect_error(
    plot(obj, type = "heatmap", specified_drugs = 1)
  )

  expect_error(
    plot(obj, type = "heatmap", cutoff_signal = -1)
  )

  expect_error(
    summary(obj, return = 123)
  )
})
