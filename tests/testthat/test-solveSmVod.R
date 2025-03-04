test_that("Correct VOD values chosen", {
  # test values
  sm <- seq(0.2, 0.8, by = 0.2)
  vod_test <- seq(0.5, 1, by = 0.2)
  air <- 300
  soil <- 275
  omega <- 0.03
  angle <- 40 # degrees
  clay <- 0.232
  h <- 0.16
  tbH <- 280
  tbV <- 285

  # run function
  sol1 <- solveSmVod(
    smc = sm, vod = vod_test,
    tbH = tbH, tbV = tbV,
    Tair = air, Tsoil = soil,
    omega = omega, clay_frac = clay,
    roughness = h, inc_angle = angle,
    mat = T
  )

  # check that the correct values are chosen
  expect_equal(as.numeric(row.names(sol1$min_cf_index)), sol1$sm_est)
  expect_true(sol1$vod_est %in% vod_test)
  expect_equal(min(sol1$cf_mat), sol1$cf_tb)

  # check that the H and V cf sum is the same as what is reported for the total cf
  sum_cfHV <- sol1$cf_tbV + sol1$cf_tbH
  expect_equal(sum_cfHV, sol1$cf_tb)
})
