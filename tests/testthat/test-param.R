test_that("param.causal.model() generates a list with the default values of parameters", {
  expect_equal(param.causal.model(), list(c(p_L1 = 0.5, p_L2 = 0.2, p_L3 = 0.7),
                                          c(b_A1 = 0.10, b_L1_A1 = 0.15, b_L2_A1 = 0.25),
                                          c(b_A2 = 0.15, b_L1_A2 = 0.20, b_L3_A2 = 0.20),
                                          c(b_Y = 0.10, b_L1_Y = 0.02, b_L2_Y = 0.02, b_L3_Y = -0.02,
                                            b_A1_Y = 0.30, b_A2_Y = 0.10, b_A1A2_Y = 0.40),
                                          NULL))
})
