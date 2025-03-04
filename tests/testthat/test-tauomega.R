test_that("Output is the correct",{
  thetar = 40*pi/180
  h = 0.1
  rhfac = exp(-(h*cos(thetar)))
  Tair = 300
  Tsoil = 290
  vod = 0.8
  gamma = exp(-1 * vod / cos(thetar))
  omega = 0.08

  #calculate test reflectivity
  test_eps <- mironov(1.4e9, 0.4, 0.232)
  reflecs <- fresnelr(test_eps$dielectric, 40)
  rH <- reflecs$fH
# est brightness temps
tb <- tau_omega(
  r = rH,
  rhfac = rhfac,
  gamma = gamma,
  Tair = Tair,
  Tsoil = Tsoil,
  omega = omega,
  emiss = F
)

emiss <- tau_omega(
  r = rH,
  rhfac = rhfac,
  gamma = gamma,
  omega = omega,

  emiss = T
)
# check class
  expect_s3_class(tb, "tb")

# check that the emissivity and brightness temps are equivalent
  tbs_est <- tb$tbs/Tsoil
  tbc_est <- tb$tbc/Tair
  test_emiss <- tbs_est + tbc_est

  expect_equal(test_emiss, emiss$emissivity)
})


