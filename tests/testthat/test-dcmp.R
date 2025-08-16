test_that("dcmp basic functionality", {
  # Very basic test - just check it returns a valid probability
  result <- dcmp(5, mu = 10, phi = 1)
  expect_true(is.numeric(result))
  expect_true(result >= 0)
  expect_true(result <= 1)
  expect_length(result, 1)
})
