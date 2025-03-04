test_that("All residuals are returned", {

H <- seq(200,250, length.out = 10)
V <- seq(260,280, length.out = 10)

pH <- seq(200,250, length.out = 10) + 1
pV <- seq(260,280, length.out = 10) + 2

res_test <- tbResiduals(tbH = H, tbV = V, tbHpred = pH, tbVpred = pV)

expect_true(all(res_test$tbH==1))
expect_true(all(res_test$tbV==4))
expect_true(all(res_test$totaltb==5))
expect_true(all(res_test$rootSE==sqrt(res_test$totaltb)))

# # Create data with oscillations
# time <- seq(0, 100, by = 0.1)  # Time index
# freq <- runif(1, 0.05, 0.2)  # Random frequency
# amp <- (300 - 250) / 2  # Amplitude to oscillate between 250 and 300
# noise_level <- 0.5  # Noise level
#
# signal <- 250 + amp + amp * sin(2 * pi * freq * time) + rnorm(length(time), mean = 0, sd = noise_level)
})


test_that("Error handling works", {

  expect_error(tbResiduals(tbH = H, tbV = V, tbHpred = c(pH,1), tbVpred = pV),"tbH and tbHpred have different lengths.")
  expect_error(tbResiduals(tbH = H, tbV = V, tbHpred = pH, tbVpred = c(pV,1)),"tbV and tbVpred have different lengths.")
  expect_warning(tbResiduals(tbH = H, tbV = V, tbHpred = H, tbVpred = pV),"Input vectors are identical.")

})
