test_that("multiplication works", {
  valid_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2) %>%
    .set_default_names()
  random_table <- generate_contin_table(ref_table = valid_matrix)[[1]]
  expect_equal(dim(random_table), dim(valid_matrix))

  random_table2 <- generate_contin_table(ref_table = valid_matrix,
    Variation = TRUE
  )[[1]][[1]]
  expect_equal(dim(random_table2), dim(valid_matrix))
})
