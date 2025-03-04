test_that("Check errors", {
  # check that the lengths are the same.
  expect_error(
    retVOD(
      tbH = seq(100, 200, length.out = 10),
      tbV = seq(100, 200, length.out = 10),
      smc = seq(0.5, 0.9, length.out = 9), vod = 1, Tair = 1, Tsoil = 1,
      omega = 1, roughness = 1, inc_angle = 1
    ),
    "tbH, tbV, and smc lengths differ."
  )
})

test_that("Output matches manual single measurement solution", {
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
  sm <- retvod:::gen_sin(1000, rangeL = 0.2, rangeH = 0.45) |> sample(size = 7)
  vod <- seq(0,3, length.out = 30)
  inc_angle <- 40

  sol2 <- solveSmVod(
    smc = sm[4], vod = vod, tbH = h[4], tbV = v[4],
    Tair = air[4], Tsoil = soil[4],
    omega = 0.5, clay_frac = 0.232,
    roughness = 0.1, inc_angle = 40,
    mat = T
  )

  retrieval <- retVOD(h,v,
                      smc = sm, vod = vod,
                      Tair = air, Tsoil = soil,
                      omega= 0.5, roughness = 0.1, inc_angle = 40)

  expect_identical(sm, retrieval$smEst) # soil moisture input and output are same
  expect_equal(retrieval$vodEst[4],sol2$vod_est) # output from solution is same as retrieved values
  expect_equal(retrieval$cfEst[4],sol2$cf_tb)
  expect_equal(retrieval$tbVpred[4], sol2$pred_tbV)
  expect_equal(retrieval$tbHpred[4], sol2$pred_tbH)
  expect_equal(retrieval$tbHcost[4], sol2$cf_tbH)
  expect_equal(retrieval$tbVcost[4], sol2$cf_tbV)
})
