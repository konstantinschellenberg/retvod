test_that("Check input values", {

theta = 40*pi/180

Tair = 300
Tsoil = 290
#r = 1
roughness = 0.1
vod=0.6
rhfac = exp(-roughness*cos(theta))
gamma = exp(-1 * vod / cos(theta))
omega = 0.08
noT = F

tau_omega(r = 0.5, rhfac = 0.1, gamma = 0.1, Tair = 0.1, Tsoil = 0.1, omega = 0.1, noT = F)


})



