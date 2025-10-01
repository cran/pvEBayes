test_that("pvEBayes object associated functions", {
  valid_matrix <- matrix(c(70, 40, 44, 50, 52, 60, 70, 80), nrow = 2)

  obj <- pvEBayes(
    contin_table = valid_matrix,
    model = "general-gamma", alpha = 0.1
  )
  plot1 <- eyeplot_pvEBayes(obj, text_shift = 0.1)
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
})
