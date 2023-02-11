test_that("generate.data() generates a data frame", {
  set.seed(12345)
  expect_equal(class(generate.data(N = 1000, b = param.causal.model())), "data.frame")
})
