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

test_that("Same values backward and forward", {
  v <- c(
    268.4722473, 268.003972, 267.1787049,
    266.446253, 266.4795018, 266.1335705,
    265.6720042
  )
  h <- c(
    258.957763, 258.4377463, 257.8380787,
    257.2369798, 257.1168902, 256.7121893,
    256.5970527
  )

  air <- c(
    291.32, 291.3, 290.2,
    289.24, 288.76, 288.12,
    287.88
  )

  soil <- air * 0.90

  set.seed(2)
  sm <- retvod:::gen_sin(1000, rangeL = 0.2, rangeH = 0.35) |> sample(size = 7)
  vod <- sm * 2
  inc_angle <- 40

  sol2 <- solveSmVod(
    smc = sm[4], vod = vod, tbH = h[4], tbV = v[4],
    Tair = air[4], Tsoil = soil[4],
    omega = 0.5, clay_frac = 0.232,
    roughness = 0.1, inc_angle = 40,
    mat = T
  )

  expect_equal(min(sol2$cf_mat), sol2$cf_tb)
  expect_equal(sol2$sm_est, sm[4])

  solved_vod <- sol2$vod_est
  input_sm <- sm[4]

  sol2_back <- solveSmVod(
    smc = sm, vod = solved_vod, tbH = h[4], tbV = v[4],
    Tair = air[4], Tsoil = soil[4],
    omega = 0.5, clay_frac = 0.232,
    roughness = 0.1, inc_angle = 40,
    mat = T
  )
  expect_identical(sol2[-c(1,10)], sol2_back[-c(1,10)])
})
