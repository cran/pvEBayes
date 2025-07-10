test_that("pvEBayes object associated functions", {
  valid_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2)

  obj <- pvEBayes(
    contin_table = valid_matrix,
    model = "general-gamma", alpha = 0.5
  )
  plot1 <- eyeplot_pvEBayes(obj, text_shift = 0.1)
  expect_equal(ggplot2::is_ggplot(plot1), TRUE)

  plot2 <- heatmap_pvEBayes(obj)
  expect_equal(ggplot2::is_ggplot(plot2), TRUE)
  print_tmp <- print(obj)
  expect_equal(is.pvEBayes(print_tmp), TRUE)

  summary_res <- summary(obj)
  expect_equal(is.list(summary_res), TRUE)

  likelihood_obj <- logLik(obj)
  expect_equal(is.numeric(likelihood_obj), TRUE)

  plot3 <- plot(obj, type = "eyeplot")
  plot4 <- plot(obj, type = "heatmap")
  expect_equal(ggplot2::is_ggplot(plot3), TRUE)
  expect_equal(ggplot2::is_ggplot(plot4), TRUE)
})
