test_that("Check output", {

  theta = 40*pi/180

  Tair = 300
  Tsoil = 290
  #r = 1
  roughness = 0.1
  vod=0.6
  rhfac = exp(-roughness*cos(theta))
  gamma = exp(-1 * vod / cos(theta))
  omega = 0.08

  #Check that you get the same value forward and back
  test_eps <- mironov(1.4e9, 0.4, 0.232)

  reflecs <- fresnelr(test_eps$dielectric, 40)

  tb <- tau_omega(reflecs$fH, rhfac=rhfac, gamma = gamma, Tair, Tsoil, 0.5, noT = T)

})
